// Single-colour shader used for the wireframe overlay.

#ifndef SHADERFLATCOLOR_H
#define SHADERFLATCOLOR_H

#ifndef SHADERFLATCOLOR_VERT
#define SHADERFLATCOLOR_VERT \
  "#version 100\n \
uniform mat4 view;\n \
uniform mat4 proj;\n \
attribute vec3 position;\n \
void main() {\n \
  gl_Position = proj * view * vec4(position, 1);\n \
}"
#endif

#ifndef SHADERFLATCOLOR_FRAG
#define SHADERFLATCOLOR_FRAG \
  "#version 100\n \
precision mediump float;\n \
uniform vec4 color;\n \
void main() {\n \
  gl_FragColor = color;\n \
}"
#endif

#endif
