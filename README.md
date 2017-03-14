spline
======

Spline functions in different languages (Catmull-Rom spline).

C and JavaScript implementations.

```c
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
```
