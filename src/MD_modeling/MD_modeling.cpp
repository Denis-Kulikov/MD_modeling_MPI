#include "../include/math_3d.h"
#include "../include/MD_modeling.hpp"

extern Pipeline pipeline;

int nMol;
DataMol Mol;
Vector3i initUcell;
Vector3f region;
Vector3f vSum, Total_vSum;
Prop kinEnergy, pressure, totEnergy;
FILE *result;
FILE *logFile;
double deltaT, alpha, density, rCut, temperature, size;
double timeNow, uSum, velMag, virSum, vvSum, Total_uSum, Total_virSum, Total_vvSum;
int moreCycles, stepAvg, stepCount, stepEquil, stepLimit, stepWrite;

void WritePosition()
{
    DO_MOL(nMol) fwrite(&Mol.p[i], sizeof(Vector3f), 1, result);
}

void WriteParams()
{
    fwrite(&nMol, sizeof(int), 1, result);
    fwrite(&size, sizeof(double), 1, result);
    DO_MOL(nMol) fwrite(&Mol.m[i], sizeof(double), 1, result);
    DO_MOL(nMol) fwrite(&Mol.p[i], sizeof(Vector3f), 1, result);
}

void PrintSummary ()
{
    fprintf (logFile, "%5d %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
             stepCount, timeNow, vSum.VCSum() / nMol,
             PropEst(totEnergy),PropEst (kinEnergy), PropEst (pressure));
}

void VRand(Vector3f &p)
{
    double s, x, y, z;
    s = 2.0f;
    while (s > 1.0f) {
        x = 2.0f * RandR() - 1.0f;
        y = 2.0f * RandR() - 1.0f;
        z = 2.0f * RandR() - 1.0f;
        s = Sqr(x) + Sqr(y) + Sqr(z);
    }
    p.z = 1.0f - 2.0f * s;
    s = 2.0f * sqrt(1.0f - s);
    p.x = s * x;
    p.y = s * y;
    p.z = s * z;
}

void CalculateForces (int first, int last)
{
    Vector3f dr;
    double fcVal, rr, rrCut, rri, rri3;
    rrCut = rCut * rCut;
    for (int i = first; i < last; i++) {
        Mol.f[i].x = 0.0;
        Mol.f[i].y = 0.0; 
        Mol.f[i].z = 0.0;
    }
    
    uSum = 0.0;
    virSum = 0.0;
    for (int j1 = first; j1 < last; j1++) {
        for (int j2 = j1 + 1; j2 != j1; j2 = (j2 + 1) % nMol) {
            dr.x = Mol.p[j1].x - Mol.p[j2].x; 
            dr.y = Mol.p[j1].y - Mol.p[j2].y;
            dr.z = Mol.p[j1].z - Mol.p[j2].z;

            dr.x += region.x * (1 - 2 * (dr.x >= 0.5 * region.x));
            dr.y += region.x * (1 - 2 * (dr.y >= 0.5 * region.y));
            dr.z += region.z * (1 - 2 * (dr.z >= 0.5 * region.z));

            rr = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
            if ((rr < rrCut) && (rr > 1e-2)) {
                double x = rr / rrCut;
                double smoothing = 0.5 * (1 + cos(M_PI * x)) * (1 - tanh(alpha * (1 - x)));

                rri = 1. / rr;
                rri3 = rri * rri * rri;
                fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
                Mol.f[j1].x += fcVal * dr.x * smoothing;
                Mol.f[j1].y += fcVal * dr.y * smoothing;
                Mol.f[j1].z += fcVal * dr.z * smoothing;
                uSum += 4. * rri3 * (rri3 - 1.) + 1.0;
                virSum += fcVal * rr;
            }
        }
    }
}

void LeapfrogStep (int part, int first, int last)
{
    if (part == 1) {
        for (int i = first; i < last; i++) {
            Mol.v[i] = Mol.v[i].VAdd(Mol.f[i].VScale(0.5 * deltaT));  
            Mol.p[i] = Mol.p[i].VAdd(Mol.v[i].VScale(deltaT)); 

            if (unlikely(Mol.p[i].x < -region.x)) Mol.p[i].x = (region.x + fmod(Mol.p[i].x, region.x));
            if (unlikely(Mol.p[i].x > region.x)) Mol.p[i].x = (-region.x + fmod(Mol.p[i].x, region.x));

            if (unlikely(Mol.p[i].y < -region.y)) Mol.p[i].y = (region.y + fmod(Mol.p[i].y, region.y));
            if (unlikely(Mol.p[i].y > region.y)) Mol.p[i].y = (-region.y + fmod(Mol.p[i].y, region.y));

            if (unlikely(Mol.p[i].z < -region.z)) Mol.p[i].z = (region.z + fmod(Mol.p[i].z, region.z));
            if (unlikely(Mol.p[i].z > region.z)) Mol.p[i].z = (-region.z + fmod(Mol.p[i].z, region.z));
        }
    } else {
        for (int i = first; i < last; i++) {
            Mol.v[i] = Mol.v[i].VAdd(Mol.f[i].VScale(0.5 * deltaT));  
        }
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

void GetvSum (int first, int last)
{
    vSum.VZero();
    vvSum = 0.0;
    for (int i = first; i < last; i++) {
        vSum = vSum.VAdd(Mol.v[i]);
        vvSum += Mol.v[i].VLenSq();
    }
}

void EvalProps ()
{
    kinEnergy.val = 0.5 * Total_vvSum / nMol;
    totEnergy.val = kinEnergy.val + Total_uSum / nMol;
    pressure.val = density * (Total_vvSum + Total_virSum) / (nMol * NDIM);
}

void SingleStep (int first, int last)
{
    ++stepCount;
    timeNow = stepCount * deltaT;
    if (first == 0) {
        if ((stepCount % stepWrite) == 0) {
            WritePosition();
        }
        if ((stepCount % stepAvg) == 0) {
            printf("|Step: %d\n", stepCount);
            EvalProps ();
            AccumProps (1);
            AccumProps (2);
            PrintSummary ();
            AccumProps (0);
            DO_MOL(nMol) {
                if ( Mol.v[i].x > 1000 || Mol.v[i].y > 1000 || Mol.v[i].z > 1000)
                    printf("%f %f %f\n", Mol.v[i].x, Mol.v[i].y, Mol.v[i].z);
            }
        }
    }

    LeapfrogStep(1, first, last);
    CalculateForces(first, last);
    LeapfrogStep(2, first, last);

    if (((stepCount + 1) % stepAvg) == 0) {
        GetvSum(first, last);
    }

    if (stepCount == stepLimit) {
        if (first == 0) {
            fclose(result);
            fclose(logFile);
        }
        moreCycles = 0;
    }
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
    DO_MOL(nMol) {
        Mol.p[i].x = (i % 2 == 0 ? -1 : 1) * (i / 2 % 2 == 0 ? 1 : -1) * region.x / 2.0 + 
                     (rand() / (double)RAND_MAX - 0.5) * region.x * (1 - 1e-4);
        Mol.p[i].y = (i / 4 % 2 == 0 ? 1 : -1) * region.x / 2.0 + 
                     (rand() / (double)RAND_MAX - 0.5) * region.x * (1 - 1e-4);
        Mol.p[i].z = (i / 8 % 2 == 0 ? 1 : -1) * region.z / 2.0 + 
                     (rand() / (double)RAND_MAX - 0.5) * region.z * (1 - 1e-4);
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
    result = fopen("data/result.bin", "wb");
    logFile = fopen("data/log.txt", "w");
    fprintf(logFile, "Step\tTime\tv\tE()\t\tE_kin\t\tPressure\n");
    srand(time(NULL));
    InitCoords ();
    InitVels ();
    InitMass ();
    AccumProps (0);
    WriteParams();
}

void SetParams ()
{
    moreCycles = 1;
    stepCount = 0;
    stepAvg = 200;
    stepEquil = 0;
    stepLimit = 10 * 1000;
    stepWrite = 50;
    temperature = 1;
    size = 20;
    density = 0.5;
    deltaT = 1e-4;
    alpha = 10;
    rCut = pow(size, 1.0f / 3.0f);
    region.VSet(size, size, size);
    nMol = static_cast<int>(region.x * region.y * region.z * density);
    initUcell.x = static_cast<int>(pow(nMol, 1.0 / 3.0) + 1);
    initUcell.y = initUcell.x;
    initUcell.z = static_cast<int>(nMol / initUcell.x / initUcell.y) + 1;
    velMag = sqrt (NDIM * (1.0f - 1.0f / nMol) * temperature);
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