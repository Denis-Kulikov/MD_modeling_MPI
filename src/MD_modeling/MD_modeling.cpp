#include "../include/math_3d.h"
#include "../include/MD_modeling.hpp"
#include "../include/try.hpp"

int nMol;
DataMol Mol;
Vector3f vSum, Total_vSum, region, cregion;
Prop kinEnergy, pressure, totEnergy;
FILE *result;
FILE *dataFile;
double deltaT, alpha, density, rCut, temperature, size;
double timeNow, uSum, velMag, virSum, vvSum, Total_uSum, Total_virSum, Total_vvSum;
int stepAvg, stepCount, stepEquil, stepLimit, stepWrite;


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

void WritePosition(int n)
{
    fwrite(&n, sizeof(int), 1, result);
    DO_MOL(n) fwrite(&Mol.m[i], sizeof(double), 1, result);
    DO_MOL(n) fwrite(&Mol.p[i], sizeof(Vector3f), 1, result);
}

void CloseFile(int rank)
{
    fclose(result);
}

void OpenFile(int rank)
{
    std::string file_name = "data/result" + std::to_string(rank) + ".bin";
    result = fopen(file_name.c_str(), "wb");
}

void WriteParams(int commsize)
{
    dataFile = fopen("data/params.bin", "wb");
    fwrite(&commsize, sizeof(int), 1, dataFile);
    fwrite(&nMol, sizeof(int), 1, dataFile);
    fwrite(&size, sizeof(double), 1, dataFile);
    fclose(dataFile);
}


void WriteSystem ()
{
    fprintf (dataFile, "%5d\t%8.4f\t%7.6f %7.6f %7.6f %7.6f %7.6f %7.6f %7.6f\n",
             stepCount, stepCount * deltaT, Total_vSum.VCSum() / nMol,
             PropEst(totEnergy), PropEst (kinEnergy), PropEst (pressure));
}

void CloseSystemFile()
{
    fclose(dataFile);
}

void OpenSystemFile()
{
    dataFile = fopen("data/system.txt", "w");
    fprintf (dataFile, "step, time, v, Energy, sqrEnergy, kinEnergy, sqrkinEnergy, Pressure, sqrPressure\n");
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

void EvalProps ()
{
    kinEnergy.val = 0.5 * Total_vvSum / nMol;
    totEnergy.val = kinEnergy.val + fabs(Total_uSum) / nMol;
    pressure.val = density * (Total_vvSum + Total_virSum) / (nMol * NDIM);
}

void GetvSum (int n)
{
    vSum.VZero();
    vvSum = 0.0;
    DO_MOL(n) {
        vSum = vSum.VAdd(Mol.v[i]);
        vvSum += Mol.v[i].VLenSq();
    }
}


Escapees FindEscapees(int n, const Vector3i &crank, const Vector3i &dims, const Vector3f &center)
{
    Escapees escapee;
    for (int i = RIGHT; i <= BACK; i++) {
        TRY(((escapee.esc[i] = new(std::vector<OneMol>)) == nullptr), "Memory allocation error (escapee->escapee[]).")
        TRY(((escapee.n[i] = new(std::vector<int>)) == nullptr), "Memory allocation error (escapee->escapee[]).")
    }

    DO_MOL(n) {
        if (Mol.p[i].x > (center.x + cregion.x / 2)) {
            if (crank.x == (dims.x - 1))
                Mol.p[i].x = (-region.x + fmod(Mol.p[i].x, cregion.x));
            escapee.esc[RIGHT]->push_back(OneMol(Mol.p[i], Mol.f[i], Mol.v[i], Mol.m[i])); 
            escapee.n[RIGHT]->push_back(i);
        } else if (Mol.p[i].x < (center.x - cregion.x / 2)) {
            if (crank.x == 0)
                Mol.p[i].x = (region.x + fmod(Mol.p[i].x, cregion.x));
            escapee.esc[LEFT]->push_back(OneMol(Mol.p[i], Mol.f[i], Mol.v[i], Mol.m[i]));
            escapee.n[LEFT]->push_back(i);
        }
        else if (Mol.p[i].y > (center.y + cregion.y / 2)) {
            if (crank.y == (dims.y - 1))
                Mol.p[i].y = (-region.y + fmod(Mol.p[i].y, cregion.y));
            escapee.esc[TOP]->push_back(OneMol(Mol.p[i], Mol.f[i], Mol.v[i], Mol.m[i]));
            escapee.n[TOP]->push_back(i);
        } else if (Mol.p[i].y < (center.y - cregion.y / 2)) {
            if (crank.y == 0)
                Mol.p[i].y = (region.y + fmod(Mol.p[i].y, cregion.y));
            escapee.esc[BOTTOM]->push_back(OneMol(Mol.p[i], Mol.f[i], Mol.v[i], Mol.m[i]));
            escapee.n[BOTTOM]->push_back(i);
        }
        else if (Mol.p[i].z > (center.z + cregion.z / 2)) {
            if (crank.z == (dims.z - 1))
                Mol.p[i].z = (-region.z + fmod(Mol.p[i].z, cregion.z));
            escapee.esc[FRONT]->push_back(OneMol(Mol.p[i], Mol.f[i], Mol.v[i], Mol.m[i]));
            escapee.n[FRONT]->push_back(i);
        } else if (Mol.p[i].z < (center.z - cregion.z / 2)) {
            if (crank.z == 0)
                Mol.p[i].z = (region.z + fmod(Mol.p[i].z, cregion.z));
            escapee.esc[BACK]->push_back(OneMol(Mol.p[i], Mol.f[i], Mol.v[i], Mol.m[i]));
            escapee.n[BACK]->push_back(i);
        }
    }
    return escapee;
}

void CalculateForces (int n)
{
    Vector3f dr;
    double fcVal, rr, rrCut, rri, rri3;
    rrCut = rCut * rCut;
    DO_MOL(n) {
        Mol.f[i].x = 0.0;
        Mol.f[i].y = 0.0; 
        Mol.f[i].z = 0.0;
    }
    
    uSum = 0.0;
    virSum = 0.0;
    for (int j1 = 0; j1 < n; j1++) {
        for (int j2 = j1 + 1; j2 < n - 1; j2++) {
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

void LeapfrogStep (int part, int n)
{
    if (part == 1) {
        DO_MOL(n) {
            Mol.v[i] = Mol.v[i].VAdd(Mol.f[i].VScale(0.5 * deltaT).VDiv(Vector3f(Mol.m[i], Mol.m[i], Mol.m[i])));  
            Mol.p[i] = Mol.p[i].VAdd(Mol.v[i].VScale(deltaT).VDiv(Vector3f(Mol.m[i], Mol.m[i], Mol.m[i]))); 
        }
    } else {
        DO_MOL(n) 
            Mol.v[i] = Mol.v[i].VAdd(Mol.f[i].VScale(0.5 * deltaT)); 
    }
}

void SingleStep (int n)
{
    stepCount++;

    LeapfrogStep(1, n);
    CalculateForces(n);
    LeapfrogStep(2, n);
    
    if ((stepCount % stepWrite) == 0) WritePosition(n);
}

Vector3f GetCenter(const Vector3i &crank, const Vector3i &dims)
{
    Vector3f center(
        region.x / dims.x * (crank.x * 2 + 1) - region.x,
        region.y / dims.y * (crank.y * 2 + 1) - region.y,
        region.z / dims.z * (crank.z * 2 + 1) - region.z
    );
    return center;
}

void InitMass (int n)
{    
    DO_MOL(n) {
        Mol.m[i] = rand() / (double)RAND_MAX * 0.5 + 0.5;
        Mol.m[i] = rand() / (double)RAND_MAX * 0.5 + 0.5;
        Mol.m[i] = rand() / (double)RAND_MAX * 0.5 + 0.5;
    }
}

void InitVels (int n)
{
    srand(time(NULL));
    vSum.VZero();
    DO_MOL(n) {
        VRand (Mol.v[i]);
        Mol.v[i].VScale(velMag);
        vSum = vSum.VAdd(Mol.v[i]);
    }
    DO_MOL(n) {
        Mol.v[i] = Mol.v[i].VAdd(vSum.VScale(-1.0f / nMol));
    }
}

void InitCoords (int n, Vector3f center)
{
    DO_MOL(n) {
        Mol.p[i].x = center.x + (i % 2 == 0 ? -1 : 1) * (i / 2 % 2 == 0 ? 1 : -1) * cregion.x / 4.0 + (rand() / (double)RAND_MAX / 2 - 0.25) * cregion.x * (1 - 1e-6);
        Mol.p[i].y = center.y + (i / 4 % 2 == 0 ? 1 : -1) * cregion.y / 4.0 + (rand() / (double)RAND_MAX / 2 - 0.25) * cregion.y * (1 - 1e-6);
        Mol.p[i].z = center.z + (i / 8 % 2 == 0 ? 1 : -1) * cregion.z / 4.0 + (rand() / (double)RAND_MAX / 2 - 0.25) * cregion.z * (1 - 1e-6);
    }
}

bool AllocArrays (int commsize)
{
    if (nMol < commsize * 2)
        return false;

    Mol.p = (Vector3f*)malloc(sizeof(Vector3f) * nMol / commsize * 2);
    Mol.f = (Vector3f*)calloc(sizeof(Vector3f),  nMol / commsize * 2);
    Mol.v = (Vector3f*)malloc(sizeof(Vector3f) * nMol / commsize * 2);
    Mol.m = (double*)malloc(sizeof(double) * nMol / commsize * 2);

    return (Mol.p != nullptr) && (Mol.f != nullptr) && (Mol.v != nullptr) && (Mol.m != nullptr);
}

void SetupJob (int rank, int n, Vector3f center)
{
    srand(rank);
    InitCoords (n, center);
    InitVels (n);
    InitMass (n);
    AccumProps (0);
}

void SetParams (const Vector3i &dims)
{
    stepCount = 0;
    stepAvg = 5000;
    stepEquil = 0;
    stepLimit = 300 * 1000;
    stepWrite = 100;
    temperature = 50;
    size = 15;
    density = 0.3;
    deltaT = 1e-4;
    alpha = 10;
    rCut = pow(size, 1.0f / 3.0f);
    region.VSet(size, size, size);
    cregion.VSet(
        region.x * 2 / dims.x,
        region.y * 2 / dims.y,
        region.z * 2 / dims.z
    );
    nMol = static_cast<int>(region.x * region.y * region.z * density);
    velMag = sqrt (NDIM * (1.0f - 1.0f / nMol) * temperature);
}
