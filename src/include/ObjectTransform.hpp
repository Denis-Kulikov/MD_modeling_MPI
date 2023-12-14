#pragma once
#include "math_3d.h"

class ObjectTransform 
{
public:
    ObjectTransform()
    {
        Scale = Vector3f(1.0f, 1.0f, 1.0f);
        WorldPos = Vector3f(0.0f, 0.0f, 0.0f);
        Rotate = Vector3f(0.0f, 0.0f, 0.0f);
    }

    void SetScale(double ScaleX, double ScaleY, double ScaleZ)
    {
        Scale.x = ScaleX;
        Scale.y = ScaleY;
        Scale.z = ScaleZ;
    }

    void SetWorldPos(double x, double y, double z)
    {
        WorldPos.x = x;
        WorldPos.y = y;
        WorldPos.z = z;
    }

    void SetRotate(double RotateX, double RotateY, double RotateZ)
    {
        Rotate.x = RotateX;
        Rotate.y = RotateY;
        Rotate.z = RotateZ;
    }

    Vector3f Scale;
    Vector3f WorldPos;
    Vector3f Rotate;
};
