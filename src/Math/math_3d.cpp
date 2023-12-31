#include "../include/math_3d.h"

Vector3f Vector3f::Cross(const Vector3f& v) const
{
    const double _x = y * v.z - z * v.y;
    const double _y = z * v.x - x * v.z;
    const double _z = x * v.y - y * v.x;

    return Vector3f(_x, _y, _z);
}

Vector3f& Vector3f::Normalize()
{
    const double Length = sqrtf(x * x + y * y + z * z);

    x /= Length;
    y /= Length;
    z /= Length;

    return *this;
}

double Vector3f::VDot(const Vector3f &v)
{
    return x * v.x + y * v.y + z * v.z;
}

double Vector3f::VLenSq()
{
    return VDot(*this);
}

double Vector3f::VLen()
{
    return sqrt(VLenSq());
}

double Vector3f::VCSum()
{
    return x + y + z;
}

double Vector3f::Distance(const Vector3f &v)
{
    return sqrt(pow(v.x - x, 2) + pow(v.y - y, 2) + pow(v.z - z, 2));
}

Vector3f Vector3f::VAdd(const Vector3f &v)
{
    return Vector3f(x + v.x, y + v.y, z + v.z);
}

Vector3f Vector3f::VSub(const Vector3f &v)
{
    return Vector3f(x - v.x, y - v.y, z - v.z);
}

Vector3f Vector3f::VMul(const Vector3f &v)
{

    return Vector3f(x * v.x, y * v.y, z * v.z);
}

Vector3f Vector3f::VDiv(const Vector3f &v)
{
    return Vector3f(x / v.x, y / v.y, z / v.z);
}

Vector3f Vector3f::VDiv(const Vector3i &v)
{
    return Vector3f(x / v.x, y / v.y, z / v.z);
}

Vector3f Vector3f::VScale(double s)
{
    return VMul(Vector3f(s, s, s));
}

void Vector3f::VSet(const Vector3f &v)
{
    x = v.x;
    y = v.y;
    z = v.z;
}

void Vector3f::VSet(double s)
{
    x = s;
    y = s;
    z = s;
}

void Vector3f::VSet(double sx, double sy, double sz)
{
    x = sx;
    y = sy;
    z = sz;
}

void Vector3f::VSetI(const Vector3i &v)
{
    x = v.x;
    y = v.y;
    z = v.z;
}

void Vector3f::VZero()
{
    VSet(0);
}

Vector3f Vector3f::VVV()
{
    return Vector3f(0, 0, 0);
}

double VectorDot(Vector3f a, Vector3f b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector3f VectorNormalize(Vector3f v)
{
    double len = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    if (len > 0.0f) {
        double invLen = 1.0f / len;
        v.x *= invLen;
        v.y *= invLen;
        v.z *= invLen;
    }
    return v;
}

Vector3f VectorSubtract(Vector3f a, Vector3f b)
{
    Vector3f result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    return result;
}

Vector3f VectorAdd(Vector3f a, Vector3f b)
{
    Vector3f result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    return result;
}

Vector3f VectorScale(const Vector3f &v, double s)
{
    Vector3f result;
    result.x = v.x * s;
    result.y = v.y * s;
    result.z = v.z * s;
    return result;
}

void Vector3f::Print() const
{
    printf("%.02f, %.02f, %.02f\n", x, y, z);
}


void Matrix4f::InitScaleTransform(double ScaleX, double ScaleY, double ScaleZ)
{
    m[0][0] = ScaleX; m[0][1] = 0.0f;   m[0][2] = 0.0f;   m[0][3] = 0.0f;
    m[1][0] = 0.0f;   m[1][1] = ScaleY; m[1][2] = 0.0f;   m[1][3] = 0.0f;
    m[2][0] = 0.0f;   m[2][1] = 0.0f;   m[2][2] = ScaleZ; m[2][3] = 0.0f;
    m[3][0] = 0.0f;   m[3][1] = 0.0f;   m[3][2] = 0.0f;   m[3][3] = 1.0f;
}

void Matrix4f::InitRotateTransform(double RotateX, double RotateY, double RotateZ)
{
    Matrix4f rx, ry, rz;

    const double x = ToRadian(RotateX);
    const double y = ToRadian(RotateY);
    const double z = ToRadian(RotateZ);

    rx.m[0][0] = 1.0f; rx.m[0][1] = 0.0f   ; rx.m[0][2] = 0.0f    ; rx.m[0][3] = 0.0f;
    rx.m[1][0] = 0.0f; rx.m[1][1] = cosf(x); rx.m[1][2] = -sinf(x); rx.m[1][3] = 0.0f;
    rx.m[2][0] = 0.0f; rx.m[2][1] = sinf(x); rx.m[2][2] = cosf(x) ; rx.m[2][3] = 0.0f;
    rx.m[3][0] = 0.0f; rx.m[3][1] = 0.0f   ; rx.m[3][2] = 0.0f    ; rx.m[3][3] = 1.0f;

    ry.m[0][0] = cosf(y); ry.m[0][1] = 0.0f; ry.m[0][2] = -sinf(y); ry.m[0][3] = 0.0f;
    ry.m[1][0] = 0.0f   ; ry.m[1][1] = 1.0f; ry.m[1][2] = 0.0f    ; ry.m[1][3] = 0.0f;
    ry.m[2][0] = sinf(y); ry.m[2][1] = 0.0f; ry.m[2][2] = cosf(y) ; ry.m[2][3] = 0.0f;
    ry.m[3][0] = 0.0f   ; ry.m[3][1] = 0.0f; ry.m[3][2] = 0.0f    ; ry.m[3][3] = 1.0f;

    rz.m[0][0] = cosf(z); rz.m[0][1] = -sinf(z); rz.m[0][2] = 0.0f; rz.m[0][3] = 0.0f;
    rz.m[1][0] = sinf(z); rz.m[1][1] = cosf(z) ; rz.m[1][2] = 0.0f; rz.m[1][3] = 0.0f;
    rz.m[2][0] = 0.0f   ; rz.m[2][1] = 0.0f    ; rz.m[2][2] = 1.0f; rz.m[2][3] = 0.0f;
    rz.m[3][0] = 0.0f   ; rz.m[3][1] = 0.0f    ; rz.m[3][2] = 0.0f; rz.m[3][3] = 1.0f;

    *this = rz * ry * rx;
}

void Matrix4f::InitTranslationTransform(double x, double y, double z)
{
    m[0][0] = 1.0f; m[0][1] = 0.0f; m[0][2] = 0.0f; m[0][3] = x;
    m[1][0] = 0.0f; m[1][1] = 1.0f; m[1][2] = 0.0f; m[1][3] = y;
    m[2][0] = 0.0f; m[2][1] = 0.0f; m[2][2] = 1.0f; m[2][3] = z;
    m[3][0] = 0.0f; m[3][1] = 0.0f; m[3][2] = 0.0f; m[3][3] = 1.0f;
}


void Matrix4f::InitCameraTransform(const Vector3f& Target, const Vector3f& Up)
{
    Vector3f N = Target;
    N.Normalize();
    Vector3f U = Up;
    U.Normalize();
    U = U.Cross(N);
    Vector3f V = N.Cross(U);
       
    m[0][0] = U.x;   m[0][1] = U.y;   m[0][2] = U.z;   m[0][3] = 0.0f;
    m[1][0] = V.x;   m[1][1] = V.y;   m[1][2] = V.z;   m[1][3] = 0.0f;
    m[2][0] = N.x;   m[2][1] = N.y;   m[2][2] = N.z;   m[2][3] = 0.0f;
    m[3][0] = 0.0f;  m[3][1] = 0.0f;  m[3][2] = 0.0f;  m[3][3] = 1.0f;        
}

void Matrix4f::InitPersProjTransform(double FOV, double Width, double Height, double zNear, double zFar)
{
    const double ar         = Width / Height;
    const double zRange     = zNear - zFar;
    const double tanHalfFOV = tanf(ToRadian(FOV / 2.0f));

    m[0][0] = 1.0f/(tanHalfFOV * ar); m[0][1] = 0.0f;            m[0][2] = 0.0f;          m[0][3] = 0.0;
    m[1][0] = 0.0f;                   m[1][1] = 1.0f/tanHalfFOV; m[1][2] = 0.0f;          m[1][3] = 0.0;
    m[2][0] = 0.0f;                   m[2][1] = 0.0f;            m[2][2] = (-zNear -zFar)/zRange ; m[2][3] = 2.0f * zFar*zNear/zRange;
    m[3][0] = 0.0f;                   m[3][1] = 0.0f;            m[3][2] = 1.0f;          m[3][3] = 0.0;
}
