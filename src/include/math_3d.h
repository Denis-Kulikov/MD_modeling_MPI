#ifndef MATH_3D_H
#define	MATH_3D_H

#include <stdio.h>
#include <math.h>

#define ToRadian(x) ((x) * M_PI / 180.0f)
#define ToDegree(x) ((x) * 180.0f / M_PI)

struct Vector2f
{
    double x;
    double y;

    Vector2f()
    {
    }

    Vector2f(double _x, double _y)
    {
        x = _x;
        y = _y;
    }
};

class Vector3i
{
public:
    int x, y, z;
    Vector3i() {};
    Vector3i(int _x, int _y, int _z)
    {
        x = _x;
        y = _y;
        z = _z;
    }

    Vector3i VMul(Vector3i v) {return Vector3i(x * v.x, y * v.y, z * v.z);};
    Vector3i VScale(double s) {return VMul(Vector3i(s, s, s));};
    int* GetData() { return &x; }
};

class Vector3f
{
public:
    double x;
    double y;
    double z;

    Vector3f()
    {
    }

    Vector3f(double _x, double _y, double _z)
    {
        x = _x;
        y = _y;
        z = _z;
    }

    Vector3f& Normalize();
    double VDot(const Vector3f &v);
    double VLenSq();
    double VLen();
    double VCSum();
    double Distance(const Vector3f &v);
    Vector3f VAdd(const Vector3f &v);
    Vector3f VSub(const Vector3f &v);
    Vector3f VMul(const Vector3f &v);
    Vector3f VDiv(const Vector3f &v);
    Vector3f VDiv(const Vector3i &v);
    Vector3f VScale(double s);
    void VSet(const Vector3f &v);
    void VSet(double s);
    void VSet(double sx, double sy, double sz);
    void VSetI(const Vector3i &v);
    void VZero();
    Vector3f VVV();

    Vector3f Cross(const Vector3f& v) const;

    void Print() const;
};

class Matrix4f
{
public:
    double m[4][4];

    Matrix4f()
    {        
    }


    inline void InitIdentity()
    {
        m[0][0] = 1.0f; m[0][1] = 0.0f; m[0][2] = 0.0f; m[0][3] = 0.0f;
        m[1][0] = 0.0f; m[1][1] = 1.0f; m[1][2] = 0.0f; m[1][3] = 0.0f;
        m[2][0] = 0.0f; m[2][1] = 0.0f; m[2][2] = 1.0f; m[2][3] = 0.0f;
        m[3][0] = 0.0f; m[3][1] = 0.0f; m[3][2] = 0.0f; m[3][3] = 1.0f;
    }

    inline Matrix4f operator*(const Matrix4f& Right) const
    {
        Matrix4f Ret;

        for (unsigned int i = 0 ; i < 4 ; i++) {
            for (unsigned int j = 0 ; j < 4 ; j++) {
                Ret.m[i][j] = m[i][0] * Right.m[0][j] +
                              m[i][1] * Right.m[1][j] +
                              m[i][2] * Right.m[2][j] +
                              m[i][3] * Right.m[3][j];
            }
        }

        return Ret;
    }

    void InitScaleTransform(double ScaleX, double ScaleY, double ScaleZ);
    void InitRotateTransform(double RotateX, double RotateY, double RotateZ);
    void InitTranslationTransform(double x, double y, double z);
    void InitCameraTransform(const Vector3f& Target, const Vector3f& Up);
    void InitPersProjTransform(double FOV, double Width, double Height, double zNear, double zFar);
};


#endif	/* MATH_3D_H */
