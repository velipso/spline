// (c) Copyright 2017, Sean Connelly (@voidqk), http://syntheti.cc
// MIT License
// Project Home: https://github.com/voidqk/spline

// Catmull-Rom spline implementation (2d, 3d, 4d)

#ifndef SPLINE__H
#define SPLINE__H

#include <stdbool.h>

// vector definitions can come from NVQM (https://github.com/voidqk/nvqm), or just added here
#ifndef NVQM__H
typedef struct { float v[2]; } vec2;
typedef struct { float v[3]; } vec3;
typedef struct { float v[4]; } vec4;
#endif // NVQM__H

// `p0`, `p1`, `p2`, `p3` are points
// `tension` is number between 0 and 1 (0.5f is good default)
// `segs` is the number of segments to return via `out`
// outputs a list of points of length `segs`, where `out[0]` is `p1`, and `out[segs - 1]` is
// approaching `p2`
void spline2d_segment(vec2 p0, vec2 p1, vec2 p2, vec2 p3, float tension, int segs, vec2 *out);
void spline3d_segment(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float tension, int segs, vec3 *out);
void spline4d_segment(vec4 p0, vec4 p1, vec4 p2, vec4 p3, float tension, int segs, vec4 *out);

// calculate the size of an output spline with the given settings
int splinesize(int pointslen, int segs, bool close);

// `pointslen` is the length of `points`
// `points` is an array of points
// `tension` is a number between 0 and 1 (0.5f is good default)
// `segs` is the number of segments generated between two adjacent points
// `close` is a boolean that indicates whether the returned curve should be a closed loop
// outputs the result in `out` which must be big enough to hold the result -- use `splinesize(...)`
// to calculate the minimum size needed for a given set of parameters
void spline2d(int pointslen, const vec2 *points, float tension, int segs, bool close, vec2 *out);
void spline3d(int pointslen, const vec3 *points, float tension, int segs, bool close, vec3 *out);
void spline4d(int pointslen, const vec4 *points, float tension, int segs, bool close, vec4 *out);

#endif // SPLINE__H
