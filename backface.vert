// for raycasting
#version 400 compatibility

layout(location = 0) in vec3 VerPos;
layout(location = 1) in vec3 VerClr;

out vec3 Color;

uniform mat4 MVP;


void main()
{
    Color = VerClr;
    //gl_Position = MVP * vec4(VerPos, 1.0);
	gl_Position = gl_ModelViewProjectionMatrix * vec4(VerPos, 1.0);
}
