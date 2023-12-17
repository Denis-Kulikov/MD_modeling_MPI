#include "../include/n_body.hpp"

extern struct distance_by_index *distances;
extern Pipeline pipeline;

int iter = 0;

int nMol;
DataMol Mol;
const int size = 2;
Vector3i initUcell(size, size, size);
Vector3f region(size, size, size);
Vector3f vSum;
Prop kinEnergy, pressure, totEnergy;

double deltaT = 0.005;
double density = 0.8;

double rCut, temperature, timeNow, uSum, velMag, virSum, vvSum;
int moreCycles, stepAvg, stepCount, stepEquil, stepLimit;

void VRand(Vector3f &p) {
    double s, x, y;
    s = 2.;
    while (s > 1.) {
        x = 2. * RandR() - 1.;
        y = 2. * RandR() - 1.;
        s = Sqr(x) + Sqr(y);
    }
    p.z = 1. - 2. * s;
    s = 2. * sqrt(1. - s);
    p.x = s * x;
    p.y = s * y;
}

void move_particles(double dt)
{
    DO_MOL(nMol) {
        Vector3f dv(
            Mol.f[i].x / Mol.m[i] * dt,
            Mol.f[i].y / Mol.m[i] * dt,
            Mol.f[i].z / Mol.m[i] * dt
        );
        
        Vector3f dp(
            (Mol.v[i].x + dv.x * 0.5) * dt,
            (Mol.v[i].y + dv.y * 0.5) * dt,
            (Mol.v[i].z + dv.z * 0.5) * dt
        );

        Mol.v[i] = Mol.v[i].VAdd(dv);
        Mol.p[i] = Mol.p[i].VAdd(dp);

if (unlikely(Mol.p[i].x < -region.x)) Mol.p[i].x = fmod(region.x * 2 - Mol.p[i].x, 2 * region.x) * 0.5 + 0.5 * region.x;
if (unlikely(Mol.p[i].x > region.x)) Mol.p[i].x = fmod(-region.x * 2 + Mol.p[i].x, 2 * region.x) * 0.5 - 0.5 * region.x;

if (unlikely(Mol.p[i].y < -region.y)) Mol.p[i].y = fmod(region.y * 2 - Mol.p[i].y, 2 * region.y) * 0.5 + 0.5 * region.y;
if (unlikely(Mol.p[i].y > region.y)) Mol.p[i].y = fmod(-region.y * 2 + Mol.p[i].y, 2 * region.y) * 0.5 - 0.5 * region.y;

if (unlikely(Mol.p[i].z < -region.z)) Mol.p[i].z = fmod(region.z * 2 - Mol.p[i].z, 2 * region.z) * 0.5 + 0.5 * region.z;
if (unlikely(Mol.p[i].z > region.z)) Mol.p[i].z = fmod(-region.z * 2 + Mol.p[i].z, 2 * region.z) * 0.5 - 0.5 * region.z;


        // Mol.p[i].x = (Mol.p[i].x >= -region.x) * Mol.p[i].x + (Mol.p[i].x < -region.x) * (2 * region.x - Mol.p[i].x);  
        // Mol.p[i].x = (Mol.p[i].x <= region.x)  * Mol.p[i].x + (Mol.p[i].x > region.x) * (-2 * region.x + Mol.p[i].x); 
        // Mol.p[i].y = (Mol.p[i].y >= -region.y) * Mol.p[i].y + (Mol.p[i].y < -region.y) * (2 * region.y - Mol.p[i].y); 
        // Mol.p[i].y = (Mol.p[i].y <= region.y)  * Mol.p[i].y + (Mol.p[i].y > region.y) * (-2 * region.y + Mol.p[i].y); 
        // Mol.p[i].z = (Mol.p[i].z >= -region.z) * Mol.p[i].z + (Mol.p[i].z < -region.z) * (2 * region.z - Mol.p[i].z); 
        // Mol.p[i].z = (Mol.p[i].z <= region.z)  * Mol.p[i].z + (Mol.p[i].z > region.z) * (-2 * region.z + Mol.p[i].z); 

        // Mol.f[i].x = Mol.f[i].y = Mol.f[i].z = 0;

        distances[i].index = i;
        distances[i].dist = Mol.p[i].Distance(pipeline.camera.Params.WorldPos);
    }
}

void calculate_forces ()
{
    Vector3f dr;
    double fcVal, rr, rrCut, rri, rri3;
    int j1, j2;
    rrCut = rCut * rCut;
    DO_MOL (nMol) {
        Mol.f[i].x = 0.0;
        Mol.f[i].y = 0.0; 
        Mol.f[i].z = 0.0;
    }
    
    double uSum = 0.0;
    for (j1 = 0; j1 < nMol - 1; j1 ++) {
        for (j2 = j1 + 1; j2 < nMol; j2 ++) {
            dr.x = Mol.p[j1].x - Mol.p[j2].x; 
            dr.y = Mol.p[j1].y - Mol.p[j2].y;
            dr.z = Mol.p[j1].z - Mol.p[j2].z;

            dr.x += region.x * (1 - 2 * (dr.x >= 0.5 * region.x));
            dr.y += region.x * (1 - 2 * (dr.y >= 0.5 * region.y));
            dr.z += region.z * (1 - 2 * (dr.z >= 0.5 * region.z));

            rr = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
            if (rr < rrCut) {
                rri = 1. / rr;
                rri3 = rri * rri * rri;
                fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
                Mol.f[j1].x += fcVal * dr.x;
                Mol.f[j1].y += fcVal * dr.y;
                Mol.f[j1].z += fcVal * dr.z;
                Mol.f[j2].x -= fcVal * dr.x;
                Mol.f[j2].y -= fcVal * dr.y;
                Mol.f[j2].z -= fcVal * dr.z;
                uSum += 4. * rri3 * (rri3 - 1.) + 1.0;
                virSum += fcVal * rr;
            }
        }
    }

    // if ((iter % 500) == 0) {
        // DO_MOL(nMol) printf("C: %0.2f %0.2f %0.2f\n", Mol.p[i].x, Mol.p[i].y, Mol.p[i].z);
        // printf("\n");
    // }
    iter++;
}

void AccumProps (int icode)
{
    if (icode == 0) {
        totEnergy.PropZero();
        kinEnergy.PropZero();
        pressure .PropZero();
    } else if (icode == 1) {
        totEnergy.PropAccum();
        kinEnergy.PropAccum();
        pressure .PropAccum();
    } else if (icode == 2) {
        totEnergy.PropAvg(stepAvg);
        kinEnergy.PropAvg(stepAvg);
        pressure .PropAvg(stepAvg);
    }
}

void EvalVelDist ()
{
    // double deltaV, histSum;
    // int j, n;

    // if (countVel == 0) {
    // for (j = 0; j < sizeHistVel; j ++) histVel[j] = 0.;
    // }
    // deltaV = rangeVel / sizeHistVel;
    // DO_MOL { 10
    // j = VLen (mol[n].rv) / deltaV;
    // ++ histVel[Min (j, sizeHistVel - 1)];
    // }
    // ++ countVel;
    // if (countVel == limitVel) { 
    // histSum = 0.;
    // for (j = 0; j < sizeHistVel; j ++) histSum += histVel[j];
    // for (j = 0; j < sizeHistVel; j ++) histVel[j] /= histSum;
    // PrintVelDist (stdout);
    // countVel = 0;
    // }
}

void EvalProps ()
{
    double vv;
    vSum.VZero();
    vvSum = 0.0;
    DO_MOL(nMol) {
        vSum = vSum.VAdd(Mol.v[i]);
        vv = Mol.v[i].VLenSq();
        vvSum += vv;
    }
    kinEnergy.val = 0.5 * vvSum / nMol;
    totEnergy.val = kinEnergy.val + uSum / nMol;
    pressure.val = density * (vvSum + virSum) / (nMol * NDIM);
}

void InitMass ()
{    
    DO_MOL(nMol) {
        Mol.m[i] = rand() / (double)RAND_MAX * 0.5 + 0.5;
        Mol.m[i] = rand() / (double)RAND_MAX * 0.5 + 0.5;
        Mol.m[i] = rand() / (double)RAND_MAX * 0.5 + 0.5;
    }
}

void InitVels ()
{
    srand(time(NULL));
    vSum.VZero();
    DO_MOL(nMol) {
        VRand (Mol.v[i]);
        Mol.v[i].VScale(velMag);
        vSum = vSum.VAdd(Mol.v[i]);
    }
    DO_MOL(nMol) {
        Mol.v[i] = Mol.v[i].VAdd(vSum.VScale(-1.0f / nMol));
    }
}

void InitCoords ()
{
    Vector3f c, gap;
    int n = 0;
    gap = region.VDiv(initUcell);
    for (int ny = 0; ny < initUcell.y; ny++) {
        for (int nx = 0; nx < initUcell.x; nx++) {
            for (int nz = 0; nz < initUcell.z; nz++) {
                c.VSet (nx - 0.5f * region.x, ny - 0.5f * region.y, nz - 0.5f * region.z); 
                c = c.VMul(gap);
                c = c.VScale(2);
                Mol.p[n] = c;
                n++;
            }
        } 
    }
}

bool AllocArrays ()
{
    if (nMol < 2) {
        return false;
    }

    Mol.p = (Vector3f*)malloc(sizeof(Vector3f) * nMol);
    Mol.f = (Vector3f*)calloc(sizeof(Vector3f),  nMol);
    Mol.v = (Vector3f*)malloc(sizeof(Vector3f) * nMol);
    Mol.m = (double*)malloc(sizeof(double) * nMol);

    return (Mol.p != nullptr) && (Mol.f != nullptr) && (Mol.v != nullptr) && (Mol.m != nullptr);
}

void SetupJob ()
{
    TRY(!AllocArrays(), "Memory alloc error (SetupJob).");
    stepCount = 0;
    srand(time(NULL));
    InitCoords ();
    InitVels ();
    InitMass ();
    // AccumProps (0);
}

void SingleStep ()
{
    ++ stepCount;
    timeNow = stepCount * deltaT;
    // LeapfrogStep (1);
    // ApplyBoundaryCond ();
    // ComputeForces ();
    // LeapfrogStep (2);
    // EvalProps ();
    // AccumProps (1); 10
    // if (stepCount % stepAvg == 0) {
    //     AccumProps (2);
    //     PrintSummary (stdout);
    //     AccumProps (0);
    // }
}

void SetParams ()
{
    region.VSet(2.0f, 2.0f, 4.0f);
    nMol = static_cast<int>(region.x * region.y * region.z * density);
    Vector3f Relations(1, region.y / region.x, region.z / region.x);
    double SumRelations = Relations.x + Relations.y + Relations.z;
    printf("%f %f %f\n", Relations.x, Relations.y, Relations.z);
    printf("%f %f %f\n", nMol * Relations.x / SumRelations, nMol * Relations.y / SumRelations, nMol * Relations.z / SumRelations);
    initUcell.x = static_cast<int>(std::round(std::pow(nMol * region.x / (region.x * region.y * region.z), 1.0f / 3.0f)) + 1);
    initUcell.y = static_cast<int>(std::round(initUcell.x * region.y / region.x) + 1);
    initUcell.z = static_cast<int>(std::round(initUcell.x * region.z / region.x) + 1);
    printf("%d %d %d\n", initUcell.x, initUcell.y, initUcell.z);
    rCut = pow(2.0f, 1.0f / 6.0f);
    velMag = sqrt (NDIM * (1.0f - 1.0f / nMol) * temperature);
}

void MoveBody()
{
    double dt = 1e-6;
    calculate_forces(); 
    move_particles(dt);
}

void init()
{
    density = 10;
    SetParams();
    SetupJob ();
    printf("%d\n", nMol);
}
