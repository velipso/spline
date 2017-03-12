// (c) Copyright 2017, Sean Connelly (@voidqk), http://syntheti.cc
// MIT License
// Project Home: https://github.com/voidqk/spline

#include "spline.h"
#include <math.h>

static inline float tb2d(float ax, float ay, float bx, float by, float tension){
	float dx = bx - ax;
	float dy = by - ay;
	float res = powf(dx * dx + dy * dy, tension * 0.5f);
	return res < 0.0001f ? 0.0001f : res;
}

void spline2d_segment(vec2 p0, vec2 p1, vec2 p2, vec2 p3, float tension, int segs, vec2 *out){
	float
		p0x = p0.v[0], p0y = p0.v[1],
		p1x = p1.v[0], p1y = p1.v[1],
		p2x = p2.v[0], p2y = p2.v[1],
		p3x = p3.v[0], p3y = p3.v[1];

	float
		t1 = tb2d(p0x, p0y, p1x, p1y, tension),
		t2 = tb2d(p1x, p1y, p2x, p2y, tension) + t1,
		t3 = tb2d(p2x, p2y, p3x, p3y, tension) + t2;

	float
		dt21 = t2 - t1,
		dt31 = t3 - t1,
		dt32 = t3 - t2;

	for (int s = 0; s < segs; s++){
		float t = dt21 * s / segs + t1;

		float
			dt1 = t1 - t,
			dt2 = t2 - t,
			dt3 = t3 - t;

		float
			A1x = (dt1 * p0x + t * p1x) / t1,
			A1y = (dt1 * p0y + t * p1y) / t1;

		float
			A2x = (dt2 * p1x - dt1 * p2x) / dt21,
			A2y = (dt2 * p1y - dt1 * p2y) / dt21;

		float
			A3x = (dt3 * p2x - dt2 * p3x) / dt32,
			A3y = (dt3 * p2y - dt2 * p3y) / dt32;

		float
			B1x = (dt2 * A1x + t * A2x) / t2,
			B1y = (dt2 * A1y + t * A2y) / t2;

		float
			B2x = (dt3 * A2x - dt1 * A3x) / dt31,
			B2y = (dt3 * A2y - dt1 * A3y) / dt31;

		out[s] = (vec2){
			(dt2 * B1x - dt1 * B2x) / dt21,
			(dt2 * B1y - dt1 * B2y) / dt21
		};
	}
}

static inline float tb3d(float ax, float ay, float az, float bx, float by, float bz, float tension){
	float dx = bx - ax;
	float dy = by - ay;
	float dz = bz - az;
	float res = powf(dx * dx + dy * dy + dz * dz, tension * 0.5f);
	return res < 0.0001f ? 0.0001f : res;
}

void spline3d_segment(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float tension, int segs, vec3 *out){
	float
		p0x = p0.v[0], p0y = p0.v[1], p0z = p0.v[2],
		p1x = p1.v[0], p1y = p1.v[1], p1z = p1.v[2],
		p2x = p2.v[0], p2y = p2.v[1], p2z = p2.v[2],
		p3x = p3.v[0], p3y = p3.v[1], p3z = p3.v[2];

	float
		t1 = tb3d(p0x, p0y, p0z, p1x, p1y, p1z, tension),
		t2 = tb3d(p1x, p1y, p1z, p2x, p2y, p2z, tension) + t1,
		t3 = tb3d(p2x, p2y, p2z, p3x, p3y, p3z, tension) + t2;

	float
		dt21 = t2 - t1,
		dt31 = t3 - t1,
		dt32 = t3 - t2;

	for (int s = 0; s < segs; s++){
		float t = dt21 * s / segs + t1;

		float
			dt1 = t1 - t,
			dt2 = t2 - t,
			dt3 = t3 - t;

		float
			A1x = (dt1 * p0x + t * p1x) / t1,
			A1y = (dt1 * p0y + t * p1y) / t1,
			A1z = (dt1 * p0z + t * p1z) / t1;

		float
			A2x = (dt2 * p1x - dt1 * p2x) / dt21,
			A2y = (dt2 * p1y - dt1 * p2y) / dt21,
			A2z = (dt2 * p1z - dt1 * p2z) / dt21;

		float
			A3x = (dt3 * p2x - dt2 * p3x) / dt32,
			A3y = (dt3 * p2y - dt2 * p3y) / dt32,
			A3z = (dt3 * p2z - dt2 * p3z) / dt32;

		float
			B1x = (dt2 * A1x + t * A2x) / t2,
			B1y = (dt2 * A1y + t * A2y) / t2,
			B1z = (dt2 * A1z + t * A2z) / t2;

		float
			B2x = (dt3 * A2x - dt1 * A3x) / dt31,
			B2y = (dt3 * A2y - dt1 * A3y) / dt31,
			B2z = (dt3 * A2z - dt1 * A3z) / dt31;

		out[s] = (vec3){
			(dt2 * B1x - dt1 * B2x) / dt21,
			(dt2 * B1y - dt1 * B2y) / dt21,
			(dt2 * B1z - dt1 * B2z) / dt21
		};
	}
}

static inline float tb4d(float ax, float ay, float az, float aw, float bx, float by, float bz, float bw,
	float tension){
	float dx = bx - ax;
	float dy = by - ay;
	float dz = bz - az;
	float dw = bw - aw;
	float res = powf(dx * dx + dy * dy + dz * dz + dw * dw, tension * 0.5f);
	return res < 0.0001f ? 0.0001f : res;
}

void spline4d_segment(vec4 p0, vec4 p1, vec4 p2, vec4 p3, float tension, int segs, vec4 *out){
	float
		p0x = p0.v[0], p0y = p0.v[1], p0z = p0.v[2], p0w = p0.v[3],
		p1x = p1.v[0], p1y = p1.v[1], p1z = p1.v[2], p1w = p1.v[3],
		p2x = p2.v[0], p2y = p2.v[1], p2z = p2.v[2], p2w = p2.v[3],
		p3x = p3.v[0], p3y = p3.v[1], p3z = p3.v[2], p3w = p3.v[3];

	float
		t1 = tb4d(p0x, p0y, p0z, p0w, p1x, p1y, p1z, p1w, tension),
		t2 = tb4d(p1x, p1y, p1z, p1w, p2x, p2y, p2z, p2w, tension) + t1,
		t3 = tb4d(p2x, p2y, p2z, p2w, p3x, p3y, p3z, p3w, tension) + t2;

	float
		dt21 = t2 - t1,
		dt31 = t3 - t1,
		dt32 = t3 - t2;

	for (int s = 0; s < segs; s++){
		float t = dt21 * s / segs + t1;

		float
			dt1 = t1 - t,
			dt2 = t2 - t,
			dt3 = t3 - t;

		float
			A1x = (dt1 * p0x + t * p1x) / t1,
			A1y = (dt1 * p0y + t * p1y) / t1,
			A1z = (dt1 * p0z + t * p1z) / t1,
			A1w = (dt1 * p0w + t * p1w) / t1;

		float
			A2x = (dt2 * p1x - dt1 * p2x) / dt21,
			A2y = (dt2 * p1y - dt1 * p2y) / dt21,
			A2z = (dt2 * p1z - dt1 * p2z) / dt21,
			A2w = (dt2 * p1w - dt1 * p2w) / dt21;

		float
			A3x = (dt3 * p2x - dt2 * p3x) / dt32,
			A3y = (dt3 * p2y - dt2 * p3y) / dt32,
			A3z = (dt3 * p2z - dt2 * p3z) / dt32,
			A3w = (dt3 * p2w - dt2 * p3w) / dt32;

		float
			B1x = (dt2 * A1x + t * A2x) / t2,
			B1y = (dt2 * A1y + t * A2y) / t2,
			B1z = (dt2 * A1z + t * A2z) / t2,
			B1w = (dt2 * A1w + t * A2w) / t2;

		float
			B2x = (dt3 * A2x - dt1 * A3x) / dt31,
			B2y = (dt3 * A2y - dt1 * A3y) / dt31,
			B2z = (dt3 * A2z - dt1 * A3z) / dt31,
			B2w = (dt3 * A2w - dt1 * A3w) / dt31;

		out[s] = (vec4){
			(dt2 * B1x - dt1 * B2x) / dt21,
			(dt2 * B1y - dt1 * B2y) / dt21,
			(dt2 * B1z - dt1 * B2z) / dt21,
			(dt2 * B1w - dt1 * B2w) / dt21
		};
	}
}

int splinesize(int pointslen, int segs, bool close){
	if (pointslen <= 0)
		return 0;
	else if (pointslen == 1)
		return 1;
	if (close)
		return segs * pointslen;
	return segs * (pointslen - 1) + 1;
}

#define MAKEFUNC(splineXd, vecX, spline)                                                           \
	void splineXd(int ptslen, const vecX *points, float tension, int segs, bool close, vecX *out){ \
		const vecX *P = points;                                                                    \
		int L = ptslen;                                                                            \
		/* */                                                                                      \
		if (L <= 0)                                                                                \
			return;                                                                                \
		else if (L == 1){                                                                          \
			out[0] = P[0];                                                                         \
			return;                                                                                \
		}                                                                                          \
		else if (L == 2){                                                                          \
			if (close){                                                                            \
				spline(P[1], P[0], P[1], P[0], tension, segs, &out[0]);                            \
				spline(P[0], P[1], P[0], P[1], tension, segs, &out[segs]);                         \
			}                                                                                      \
			else{                                                                                  \
				spline(P[0], P[0], P[1], P[1], tension, segs, &out[0]);                            \
				out[segs] = P[1];                                                                  \
			}                                                                                      \
			return;                                                                                \
		}                                                                                          \
		/* */                                                                                      \
		if (close){                                                                                \
			spline(P[L - 1], P[0], P[1], P[2], tension, segs, &out[0]);                            \
			if (L == 3){                                                                           \
				spline(P[0], P[1], P[2], P[0], tension, segs, &out[segs]);                         \
				spline(P[1], P[2], P[0], P[1], tension, segs, &out[segs * 2]);                     \
				return;                                                                            \
			}                                                                                      \
		}                                                                                          \
		else{                                                                                      \
			spline(P[0], P[0], P[1], P[2], tension, segs, &out[0]);                                \
			if (L == 3){                                                                           \
				spline(P[0], P[1], P[2], P[2], tension, segs, &out[segs]);                         \
				out[segs * 2] = P[2];                                                              \
				return;                                                                            \
			}                                                                                      \
		}                                                                                          \
		/* */                                                                                      \
		for (int i = 0; i < L - 3; i++)                                                            \
			spline(P[i], P[i + 1], P[i + 2], P[i + 3], tension, segs, &out[(i + 1) * segs]);       \
		/* */                                                                                      \
		if (close){                                                                                \
			spline(P[L - 3], P[L - 2], P[L - 1], P[0], tension, segs, &out[(L - 2) * segs]);       \
			spline(P[L - 2], P[L - 1], P[0], P[1], tension, segs, &out[(L - 1) * segs]);           \
			return;                                                                                \
		}                                                                                          \
		spline(P[L - 3], P[L - 2], P[L - 1], P[L - 1], tension, segs, &out[(L - 2) * segs]);       \
		out[(L - 1) * segs] = P[L - 1];                                                            \
	}

MAKEFUNC(spline2d, vec2, spline2d_segment)
MAKEFUNC(spline3d, vec3, spline3d_segment)
MAKEFUNC(spline4d, vec4, spline4d_segment)

#undef MAKEFUNC
