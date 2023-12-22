#version 330

layout (location = 0) in vec3 Position;

uniform mat4 gWorld;
uniform float gMass;

out vec4 Color;

const float gLight = 0.125f;

void main()
{
    gl_Position = gWorld * vec4(Position, 1.0);
    Color = vec4(gMass * (1 - gLight)       + Position.x * gLight,
                (1 - gMass * 0.95) * (1 - gLight)  + Position.x * gLight,
                Position.x * gLight * 2, 1.0);
}