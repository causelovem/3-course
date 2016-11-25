#version 330

in vec4 point;
in vec2 position;
in vec4 variance;
in vec2 params; // x - rot, y - height

out vec2 TexCoord;

uniform mat4 camera;
uniform vec4 variance_wind;

void main()
{
    mat4 scaleMatrix = mat4(1.0);
    scaleMatrix[0][0] = 0.01;
    //scaleMatrix[1][1] = 0.1;
    scaleMatrix[1][1] = params.y;
    mat4 positionMatrix = mat4(1.0);
    positionMatrix[3][0] = position.x;
    positionMatrix[3][2] = position.y;

    mat4 rotMatrix = mat4(0.0);
    rotMatrix[1][1] = 1.0;
    rotMatrix[3][3] = 1.0;
    rotMatrix[0][0] = cos(params.x);
    rotMatrix[0][2] = sin(params.x);
    rotMatrix[2][0] = -sin(params.x);
    rotMatrix[2][2] = cos(params.x);

	gl_Position = camera * (positionMatrix * rotMatrix * scaleMatrix * point + variance * point.y * point.y);
	//gl_Position = camera * (positionMatrix * rotMatrix * scaleMatrix * point + variance_wind * point.y * point.y);
    TexCoord = vec2(point.x, point.y);
}
