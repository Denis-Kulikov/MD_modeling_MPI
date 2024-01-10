#include "../include/MD_modeling.hpp"
#include "../include/MPI.hpp"

#define FIRST(x, s)  (static_cast<int>(nMol / s) * x + std::min(x, nMol % s))
#define LAST(x, s)   (FIRST(x, s) + static_cast<int>(nMol / s) + (x < nMol % s ? 1 : 0))
#define NELEMS(x, s) (LAST(x, s) - FIRST(x, s) + (nMol % s) * (x == (s - 1)))

extern DataMol Mol;
extern Vector3f vSum, Total_vSum;
extern double uSum, virSum, vvSum, Total_uSum, Total_virSum, Total_vvSum;
extern int nMol, stepCount, stepLimit, stepAvg;

MPI_Datatype vector3f_type;
MPI_Datatype OneMol_type;

void Exchange(int rank, int commsize, int &nPros, int neighbours[6], MPI_Comm &cartcomm, Escapees escapees)
{
    MPI_Request reqs[12];
    MPI_Status status[12];
    MPI_Status statusProbe[6];
    OneMol* recv_buffer[6] = {nullptr};

    for (int i = RIGHT; i <= BACK; i++) { 
        MPI_Isend(escapees.esc[i]->data(), escapees.esc[i]->size(), OneMol_type, neighbours[i], 0, cartcomm, &reqs[i]); 
    }

    int recv_count[6];
    for (int i = RIGHT; i <= BACK; i++) {
        MPI_Probe(neighbours[i], 0, cartcomm, &statusProbe[i]); // Узнаём размер принимаемых данных
        MPI_Get_count(&statusProbe[i], OneMol_type, &recv_count[i]);
        if (recv_count[i] != 0) {
            TRY(((recv_buffer[i] = (OneMol*)malloc(sizeof(OneMol) * recv_count[i])) == nullptr), "Memory allocation error (recv_buffer[i]).");
            MPI_Irecv(recv_buffer[i], recv_count[i], OneMol_type, neighbours[i], 0, cartcomm, &reqs[6 + i]); // Получение в буфер
        } else {
            MPI_Irecv(nullptr, 0, OneMol_type, neighbours[i], 0, cartcomm, &reqs[6 + i]); // Повторное получение пустого сообщения
        }
    }

    if ((stepCount % stepAvg) == 0) {
        Total_uSum = 0, Total_virSum = 0, Total_vvSum = 0;
        Total_vSum.VZero();

        MPI_Reduce(&uSum,     &Total_uSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&virSum,   &Total_virSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&vvSum,    &Total_vvSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&vSum.x,   &Total_vSum.x, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    MPI_Waitall(12, reqs, status);

    int j;
    for (int i = RIGHT; i <= BACK; i++) {
        for (j = 0; (j < static_cast<int>(escapees.n[i]->size())) && (j < recv_count[i]); j++) { // замена
                Mol.p[(*escapees.n[i])[j]] = (recv_buffer[i])[j].p;
                Mol.f[(*escapees.n[i])[j]] = (recv_buffer[i])[j].f;
                Mol.v[(*escapees.n[i])[j]] = (recv_buffer[i])[j].v;
                Mol.m[(*escapees.n[i])[j]] = (recv_buffer[i])[j].m;
        }
        for (; j < recv_count[i]; j++) { // добавление в конец
            Mol.p[nPros] = (recv_buffer[i])[j].p;
            Mol.f[nPros] = (recv_buffer[i])[j].f;
            Mol.v[nPros] = (recv_buffer[i])[j].v;
            Mol.m[nPros] = (recv_buffer[i])[j].m;
            nPros++;
            if (nPros > (nMol / 2)) {
                printf("|| %d ||\n", nPros);
            }
        }
        for (int k = static_cast<int>(escapees.n[i]->size()) - j - 1; k >= 0; k--) { // удаление
                Mol.p[(*escapees.n[i])[j + k]] = Mol.p[nPros - 1];
                Mol.f[(*escapees.n[i])[j + k]] = Mol.f[nPros - 1];
                Mol.v[(*escapees.n[i])[j + k]] = Mol.v[nPros - 1];
                Mol.m[(*escapees.n[i])[j + k]] = Mol.m[nPros - 1];
                nPros--;
        }
    }

    for (int i = RIGHT; i <= BACK; i++) {
        if (recv_buffer[i] != nullptr) free(recv_buffer[i]);
        free(escapees.esc[i]);
        free(escapees.n[i]);
    }
}

int main(int argc, char *argv[])
{
    int commsize, rank;
    MPI_Init(&argc, &argv);
    double ttotal = -MPI_Wtime();
    double tLeapfrog = 0;
    double tExchange = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    MPI_Type_contiguous(3, MPI_DOUBLE, &vector3f_type);
    MPI_Type_commit(&vector3f_type);
    MPI_Type_contiguous(10, MPI_DOUBLE, &OneMol_type);
    MPI_Type_commit(&OneMol_type);

    MPI_Comm cartcomm;
    Vector3i dims(0, 0, 0), periodic(1, 1, 1);
    MPI_Dims_create(commsize, 3, dims.GetData());
    if (dims.x < 2 || dims.y < 2 || dims.z < 2) {
        fprintf(stderr, "Invalid number of processes %d: px %d py %d pz %d must be greater than 1\n", commsize, dims.x, dims.y, dims.z);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims.GetData(), periodic.GetData(), 0, &cartcomm);
    int coords[3];
    MPI_Cart_coords(cartcomm, rank, 3, coords);
    Vector3i crank(coords[0], coords[1], coords[2]);

    SetParams(dims);
    int first = FIRST(rank, commsize);
    int last = LAST(rank, commsize);
    int nPros = last - first;
    int neighbours[6];
    MPI_Cart_shift(cartcomm, 0, -1, &neighbours[RIGHT], &neighbours[LEFT]);
    MPI_Cart_shift(cartcomm, 1,  1, &neighbours[TOP], &neighbours[BOTTOM]);
    MPI_Cart_shift(cartcomm, 2,  1, &neighbours[FRONT], &neighbours[BACK]);

    TRY((!AllocArrays(commsize)), "Memory allocation error (AllocArrays).");
    Vector3f center = GetCenter(crank, dims);
    SetupJob (rank, nPros, center);

    if (rank == 0) {
        printf("NBody: %d\n", nMol);
        printf("Processes %d: px %d py %d pz %d\n", commsize, dims.x, dims.y, dims.z);
        WriteParams(commsize);
        OpenSystemFile();
    }

    printf("[%d]: NBody: %d, coords(%d, %d, %d)\n", rank, nPros, crank.x, crank.y, crank.z);

    OpenFile(rank);
    while (true) {
        tLeapfrog -= MPI_Wtime();
        SingleStep (nPros);
        tLeapfrog += MPI_Wtime();

        if ((stepCount % stepAvg) == 0) {
            GetvSum (nPros);
        }

        tExchange -= MPI_Wtime();
        Exchange(rank, commsize, nPros, neighbours, cartcomm, FindEscapees(nPros, crank, dims, center));
        tExchange += MPI_Wtime();

        if ((rank == 0) && ((stepCount % stepAvg) == 0)) {
            EvalProps ();
            AccumProps (1);
            AccumProps (2);
            WriteSystem();
            AccumProps (0);
        }

        if (stepCount >= stepLimit) break;
    }
    CloseFile(rank);
    if (rank == 0) CloseSystemFile();

    ttotal += MPI_Wtime();
    printf("[%d] ttotal: %f\n", rank, ttotal);

    if (rank == 0) printf("tLeapfrog: %f | tExchange: %f\n", tLeapfrog, tExchange);

    MPI_Type_free(&vector3f_type);
    MPI_Type_free(&OneMol_type);

    free(Mol.m);
    free(Mol.v);
    free(Mol.f);
    free(Mol.p);

    MPI_Finalize();

    return 0;
}
