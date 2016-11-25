#version 330

out vec4 outColor;
in vec2 TexCoord;

uniform sampler2D inTexture;

void main()
{
    // outColor = vec4(0.5, 0.3, 0, 0);
    outColor = texture(inTexture, TexCoord);
}
