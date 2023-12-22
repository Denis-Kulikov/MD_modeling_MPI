#include "../include/MD_modeling.hpp"
#include "../include/MPI.hpp"

#define FIRST(x, s)  (static_cast<int>(nMol / s) * x + std::min(x, nMol % s))
#define LAST(x, s)   (FIRST(x, s) + static_cast<int>(nMol / s) + (x < nMol % s ? 1 : 0))
#define NELEMS(x, s) (LAST(x, s) - FIRST(x, s) + (nMol % s) * (x == (s - 1)))

extern DataMol Mol;
extern Vector3f vSum, Total_vSum;
extern double uSum, virSum, vvSum, Total_uSum, Total_virSum, Total_vvSum;
extern int nMol, moreCycles, stepCount, stepLimit, stepAvg;

MPI_Datatype vector3f_type;
MPI_Datatype OneMol_type;

void ExchangeAndReduce(int rank, int commsize, int &nPros, int interface[6], MPI_Comm &cartcomm, Escapees escapees)
{
    // printf("[%d] Start ExchangeAndReduce\n", rank);
    MPI_Request reqs[12];
    MPI_Status status[12];
    MPI_Status statusProbe[6];
    OneMol* recv_buffer[6] = {nullptr};

    // printf("[%d] Start MPI_Isend\n", rank);
    for (int i = RIGHT; i <= BACK; i++) { 
        MPI_Isend(escapees.esc[i]->data(), escapees.esc[i]->size(), OneMol_type, interface[i], 0, cartcomm, &reqs[i]); 
    }

    int recv_count[6];
    // printf("[%d] Start MPI_Irecv\n", rank);
    for (int i = RIGHT; i <= BACK; i++) {
        MPI_Probe(interface[i], 0, cartcomm, &statusProbe[i]); // 
        MPI_Get_count(&statusProbe[i], OneMol_type, &recv_count[i]);
        if (recv_count[i] != 0) {
            TRY(((recv_buffer[i] = (OneMol*)malloc(sizeof(OneMol) * recv_count[i])) == nullptr), "Memory allocation error (recv_buffer[i]).");
            MPI_Irecv(recv_buffer[i], recv_count[i], OneMol_type, interface[i], 0, cartcomm, &reqs[6 + i]); // Получение в буфер
        } else {
            MPI_Irecv(nullptr, 0, OneMol_type, interface[i], 0, cartcomm, &reqs[6 + i]); // Повторное получение пустого сообщения
        }
    }

    MPI_Waitall(12, reqs, status);
    // printf("[%d] End MPI_Waitall\n", rank);

    int j;
    for (int i = RIGHT; i <= BACK; i++) {
        for (j = 0; (j < static_cast<int>(escapees.n[i]->size())) && (j < recv_count[i]); j++) {
            // if (j < static_cast<int>(escapees.n[i]->size())) { // замена
                Mol.p[(*escapees.n[i])[j]] = (recv_buffer[i])[j].p;
                Mol.f[(*escapees.n[i])[j]] = (recv_buffer[i])[j].f;
                Mol.v[(*escapees.n[i])[j]] = (recv_buffer[i])[j].v;
                Mol.m[(*escapees.n[i])[j]] = (recv_buffer[i])[j].m;
            // } else { // добавление в конец
            //     Mol.p[nPros] = (recv_buffer[i])[j].p;
            //     Mol.f[nPros] = (recv_buffer[i])[j].f;
            //     Mol.v[nPros] = (recv_buffer[i])[j].v;
            //     Mol.m[nPros] = (recv_buffer[i])[j].m;
            //     nPros++;
            //     if (nPros > (nMol / 2)) {
            //         printf("|| %d ||\n", nPros);
            //     }
            // }
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

        // for (; j < static_cast<int>(escapees.n[i]->size()); j++); // удаление
        
        for (int k = static_cast<int>(escapees.n[i]->size()) - j - 1; k >= 0; k--) { // удаление
                if (nPros < 5) printf("*****| %d\n", nPros);
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
    // printf("[%d] End ExchangeAndReduce\n", rank);
}

int main(int argc, char *argv[])
{
    int commsize, rank;
    MPI_Init(&argc, &argv);
    double ttotal = -MPI_Wtime();
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
        fprintf(stderr, "Invalid number of processes %d: px %d py %d pz %d must be greater than 1\n", commsize, dims.x, dims.y < 2,dims.z);
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
    int interface[6];
    MPI_Cart_shift(cartcomm, 0, 1, &interface[RIGHT], &interface[LEFT]);
    MPI_Cart_shift(cartcomm, 1, 1, &interface[TOP], &interface[BOTTOM]);
    MPI_Cart_shift(cartcomm, 2, 1, &interface[FRONT], &interface[BACK]);
    TRY((!AllocArrays(commsize)), "Memory allocation error (AllocArrays).");
    Vector3f center = GetCenter(crank, dims);
    SetupJob (rank, nPros, center);

    if (rank == 0) {
        printf("NBody: %d\n", nMol);
        WriteParams(commsize);
    }

    printf("[%d]: NBody: %d, coords(%d, %d, %d) center(%f, %f, %f)\n", rank, nPros, crank.x, crank.y, crank.z, center.x, center.y, center.x);

    OpenFile(rank);
    while (moreCycles) {
        SingleStep (nPros);
        ExchangeAndReduce(rank, commsize, nPros, interface, cartcomm, FindEscapees(nPros, crank, dims, center));
        if (stepCount >= stepLimit) moreCycles = 0;
    }
    CloseFile(rank);

    ttotal += MPI_Wtime();
    printf("[%d] ttotal: %f\n", rank, ttotal);

    int TotalnMol = 0;
    MPI_Reduce(&nPros, &TotalnMol, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) printf("Total nMol: %d\n", TotalnMol);

    MPI_Type_free(&vector3f_type);
    MPI_Type_free(&OneMol_type);

    free(Mol.m);
    free(Mol.v);
    free(Mol.f);
    free(Mol.p);

    MPI_Finalize();

    return 0;
}
