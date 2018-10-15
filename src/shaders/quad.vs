
#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec2 UV;

out vec3 FragPos;
out vec2 texCoords;

void main()
{
    FragPos = position;
    texCoords = UV;

    gl_Position = vec4(position, 1.f);
}