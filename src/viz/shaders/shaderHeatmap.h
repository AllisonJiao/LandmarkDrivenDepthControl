// Per-vertex colour shader for the heat map: the vertex colour is modulated
// by the brightness of the matcap so the surface keeps its shading cues.

#ifndef SHADERHEATMAP_H
#define SHADERHEATMAP_H

#ifndef SHADERHEATMAP_VERT
#define SHADERHEATMAP_VERT \
  "#version 100\n \
uniform mat4 view;\n \
uniform mat4 proj;\n \
uniform mat4 normalMatrix;\n \
attribute vec3 position;\n \
attribute vec3 normal;\n \
attribute vec3 color;\n \
varying vec3 e;\n \
varying vec3 n;\n \
varying vec3 c;\n \
void main() {\n \
  vec4 p = vec4(position, 1);\n \
  e = normalize( vec3( view * p ) );\n \
  n = normalize( mat3(normalMatrix) * normal );\n \
  c = color;\n \
  gl_Position = proj * view * p;\n \
}"
#endif

#ifndef SHADERHEATMAP_FRAG
#define SHADERHEATMAP_FRAG \
  "#version 100\n \
precision mediump float;\n \
uniform sampler2D tMatCap;\n \
varying vec3 e;\n \
varying vec3 n;\n \
varying vec3 c;\n \
void main() {\n \
  vec3 r = reflect( e, n );\n \
  float m = 2. * sqrt( pow( r.x, 2. ) + pow( r.y, 2. ) + pow( r.z + 1., 2. ) );\n \
  vec2 vN = r.xy / m + .5;\n \
  vec3 base = texture2D( tMatCap, vN ).rgb;\n \
  float shade = dot(base, vec3(0.299, 0.587, 0.114));\n \
  gl_FragColor = vec4(c * clamp(0.3 + shade, 0.0, 1.0), 1.0);\n \
}"
#endif

#endif
