#include "../include/n_body.hpp"

extern struct distance_by_index *distances;
extern Pipeline pipeline;

int nMol;
DataMol Mol;
Vector3i initUcell;
Vector3f region;
Vector3f vSum;
Prop kinEnergy, pressure, totEnergy;

double deltaT, density, rCut, temperature, timeNow, uSum, velMag, virSum, vvSum;
int moreCycles, stepAvg, stepCount, stepEquil, stepLimit;

void PrintSummary (const char *FileName)
{
    FILE *fp = fopen(FileName, "w");
    TRY((fp == nullptr), "Fopen error.");

    fprintf (fp, "%5d %8.4f %7.4f \n%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
             stepCount, timeNow, vSum.VCSum() / nMol,
             PropEst(totEnergy),PropEst (kinEnergy), PropEst (pressure));

    printf ("%5d %8.4f %7.4f \n%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
            stepCount, timeNow, vSum.VCSum() / nMol,
            PropEst(totEnergy),PropEst (kinEnergy), PropEst (pressure));

    fclose(fp);
}

void VRand(Vector3f &p)
{
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

void CalculateDistance()
{
    DO_MOL(nMol) {
        distances[i].index = i;
        distances[i].dist = Mol.p[i].Distance(pipeline.camera.Params.WorldPos);
    }
}

void CalculateForces ()
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
}

void LeapfrogStep (int part)
{
    if (part == 1) {
        DO_MOL(nMol) {
            Mol.v[i] = Mol.v[i].VAdd(Mol.f[i].VScale(0.5 * deltaT));  
            Mol.p[i] = Mol.p[i].VAdd(Mol.v[i].VScale(deltaT)); 

            // if (unlikely(Mol.p[i].x < -region.x)) Mol.p[i].x = fmod(region.x * 2 - Mol.p[i].x, 2 * region.x) * 0.5 + 0.5 * region.x;
            // if (unlikely(Mol.p[i].x > region.x)) Mol.p[i].x = fmod(-region.x * 2 + Mol.p[i].x, 2 * region.x) * 0.5 - 0.5 * region.x;

            // if (unlikely(Mol.p[i].y < -region.y)) Mol.p[i].y = fmod(region.y * 2 - Mol.p[i].y, 2 * region.y) * 0.5 + 0.5 * region.y;
            // if (unlikely(Mol.p[i].y > region.y)) Mol.p[i].y = fmod(-region.y * 2 + Mol.p[i].y, 2 * region.y) * 0.5 - 0.5 * region.y;

            // if (unlikely(Mol.p[i].z < -region.z)) Mol.p[i].z = fmod(region.z * 2 - Mol.p[i].z, 2 * region.z) * 0.5 + 0.5 * region.z;
            // if (unlikely(Mol.p[i].z > region.z)) Mol.p[i].z = fmod(-region.z * 2 + Mol.p[i].z, 2 * region.z) * 0.5 - 0.5 * region.z;


            if (unlikely(Mol.p[i].x < -region.x)) Mol.p[i].x = (region.x + fmod(Mol.p[i].x, region.x));
            if (unlikely(Mol.p[i].x > region.x)) Mol.p[i].x = (-region.x + fmod(Mol.p[i].x, region.x));

            if (unlikely(Mol.p[i].y < -region.y)) Mol.p[i].y = (region.y + fmod(Mol.p[i].y, region.y));
            if (unlikely(Mol.p[i].y > region.y)) Mol.p[i].y = (-region.y + fmod(Mol.p[i].y, region.y));

            if (unlikely(Mol.p[i].z < -region.z)) Mol.p[i].z = (region.z + fmod(Mol.p[i].z, region.z));
            if (unlikely(Mol.p[i].z > region.z)) Mol.p[i].z = (-region.z + fmod(Mol.p[i].z, region.z));
        }
    } else {
        DO_MOL(nMol) Mol.v[i] = Mol.v[i].VAdd(Mol.f[i].VScale(0.5 * deltaT));  
    }
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
    gap = region.VDiv(Vector3i(initUcell.x - 1, initUcell.x - 1, initUcell.x - 1));
    for (int ny = 0; ny < initUcell.y; ny++) {
        for (int nx = 0; nx < initUcell.x; nx++) {
            for (int nz = 0; nz < initUcell.z; nz++) {
                c.VSet (nx, ny, nz); 
                c = c.VMul(gap);
                c = c.VSub(Vector3f(region.x / 2, region.y / 2, region.z / 2));
                c = c.VScale(2 - 1e-16);
                Mol.p[n] = c;
                n++;
                if (n == nMol) return;
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
    AccumProps (0);
}

void SingleStep ()
{
    ++ stepCount;
    timeNow = stepCount * deltaT;
    LeapfrogStep (1);
    CalculateForces ();
    LeapfrogStep (2);
    EvalProps ();
    AccumProps (1);
    if (stepCount % stepAvg == 0) {
        AccumProps (2);
        PrintSummary ("log.out");
        AccumProps (0);
    }
    CalculateDistance();
}

void SetParams ()
{
    const int size = 3;
    moreCycles = 1;
    stepCount = 0;
    stepAvg = 500;
    stepEquil = 0;
    stepLimit = 1000 * 1000;
    deltaT = 1e-5;
    // density = 13;
    rCut = pow(size, 1.0f / 6.0f);
    density = pow(size / rCut, 3);
    // rCut = size;
    region.VSet(size, size, size);
    nMol = static_cast<int>(region.x * region.y * region.z * density);
    initUcell.x = static_cast<int>(pow(nMol, 1.0 / 3.0) + 1);
    initUcell.y = initUcell.x;
    initUcell.z = static_cast<int>(nMol / initUcell.x / initUcell.y) + 1;
    printf("%d %d %d\n", initUcell.x, initUcell.y, initUcell.z);
    velMag = sqrt (NDIM * (1.0f - 1.0f / nMol) * temperature);
}

void init()
{
    SetParams();
    SetupJob ();
    printf("%d\n", nMol);
}
