// (c) Copyright 2017, Sean Connelly (@voidqk), http://syntheti.cc
// MIT License
// Project Home: https://github.com/voidqk/spline

// Catmull-Rom spline implementation (2d, 3d, 4d)

// `p0`, `p1`, `p2`, `p3` are points in format of [x, y]
// `tension` is number between 0 and 1 (0.5 is good default)
// `segs` is the number of segments to return
// returns a list of points of length `segs`, where `result[0]` is `p1`, and `result[segs - 1]` is
// approaching `p2`
function spline2d_segment(p0, p1, p2, p3, tension, segs){
	var p0x = p0[0], p0y = p0[1];
	var p1x = p1[0], p1y = p1[1];
	var p2x = p2[0], p2y = p2[1];
	var p3x = p3[0], p3y = p3[1];

	function tb(ax, ay, bx, by){
		var dx = bx - ax;
		var dy = by - ay;
		return Math.max(Math.pow(dx * dx + dy * dy, tension * 0.5), 0.0001);
	}
	var t1 = tb(p0x, p0y, p1x, p1y);
	var t2 = tb(p1x, p1y, p2x, p2y) + t1;
	var t3 = tb(p2x, p2y, p3x, p3y) + t2;

	var dt21 = t2 - t1;
	var dt31 = t3 - t1;
	var dt32 = t3 - t2;

	var res = [];
	for (var s = 0; s < segs; s++){
		var t = dt21 * s / segs + t1;

		var dt1 = t1 - t;
		var dt2 = t2 - t;
		var dt3 = t3 - t;

		var A1x = (dt1 * p0x + t * p1x) / t1;
		var A1y = (dt1 * p0y + t * p1y) / t1;

		var A2x = (dt2 * p1x - dt1 * p2x) / dt21;
		var A2y = (dt2 * p1y - dt1 * p2y) / dt21;

		var A3x = (dt3 * p2x - dt2 * p3x) / dt32;
		var A3y = (dt3 * p2y - dt2 * p3y) / dt32;

		var B1x = (dt2 * A1x + t * A2x) / t2;
		var B1y = (dt2 * A1y + t * A2y) / t2;

		var B2x = (dt3 * A2x - dt1 * A3x) / dt31;
		var B2y = (dt3 * A2y - dt1 * A3y) / dt31;

		res.push([
			(dt2 * B1x - dt1 * B2x) / dt21,
			(dt2 * B1y - dt1 * B2y) / dt21
		]);
	}
	return res;
}

// `p0`, `p1`, `p2`, `p3` are points in format of [x, y, z]
// `tension` is number between 0 and 1 (0.5 is good default)
// `segs` is the number of segments to return
// returns a list of points of length `segs`, where `result[0]` is `p1`, and `result[segs - 1]` is
// approaching `p2`
function spline3d_segment(p0, p1, p2, p3, tension, segs){
	var p0x = p0[0], p0y = p0[1], p0z = p0[2];
	var p1x = p1[0], p1y = p1[1], p1z = p1[2];
	var p2x = p2[0], p2y = p2[1], p2z = p2[2];
	var p3x = p3[0], p3y = p3[1], p3z = p3[2];

	function tb(ax, ay, az, bx, by, bz){
		var dx = bx - ax;
		var dy = by - ay;
		var dz = bz - az;
		return Math.max(Math.pow(dx * dx + dy * dy + dz * dz, tension * 0.5), 0.0001);
	}
	var t1 = tb(p0x, p0y, p0z, p1x, p1y, p1z);
	var t2 = tb(p1x, p1y, p1z, p2x, p2y, p2z) + t1;
	var t3 = tb(p2x, p2y, p2z, p3x, p3y, p3z) + t2;

	var dt21 = t2 - t1;
	var dt31 = t3 - t1;
	var dt32 = t3 - t2;

	var res = [];
	for (var s = 0; s < segs; s++){
		var t = dt21 * s / segs + t1;

		var dt1 = t1 - t;
		var dt2 = t2 - t;
		var dt3 = t3 - t;

		var A1x = (dt1 * p0x + t * p1x) / t1;
		var A1y = (dt1 * p0y + t * p1y) / t1;
		var A1z = (dt1 * p0z + t * p1z) / t1;

		var A2x = (dt2 * p1x - dt1 * p2x) / dt21;
		var A2y = (dt2 * p1y - dt1 * p2y) / dt21;
		var A2z = (dt2 * p1z - dt1 * p2z) / dt21;

		var A3x = (dt3 * p2x - dt2 * p3x) / dt32;
		var A3y = (dt3 * p2y - dt2 * p3y) / dt32;
		var A3z = (dt3 * p2z - dt2 * p3z) / dt32;

		var B1x = (dt2 * A1x + t * A2x) / t2;
		var B1y = (dt2 * A1y + t * A2y) / t2;
		var B1z = (dt2 * A1z + t * A2z) / t2;

		var B2x = (dt3 * A2x - dt1 * A3x) / dt31;
		var B2y = (dt3 * A2y - dt1 * A3y) / dt31;
		var B2z = (dt3 * A2z - dt1 * A3z) / dt31;

		res.push([
			(dt2 * B1x - dt1 * B2x) / dt21,
			(dt2 * B1y - dt1 * B2y) / dt21,
			(dt2 * B1z - dt1 * B2z) / dt21
		]);
	}
	return res;
}

// `p0`, `p1`, `p2`, `p3` are points in format of [x, y, z, w]
// `tension` is number between 0 and 1 (0.5 is good default)
// `segs` is the number of segments to return
// returns a list of points of length `segs`, where `result[0]` is `p1`, and `result[segs - 1]` is
// approaching `p2`
function spline4d_segment(p0, p1, p2, p3, tension, segs){
	var p0x = p0[0], p0y = p0[1], p0z = p0[2], p0w = p0[3];
	var p1x = p1[0], p1y = p1[1], p1z = p1[2], p1w = p1[3];
	var p2x = p2[0], p2y = p2[1], p2z = p2[2], p2w = p2[3];
	var p3x = p3[0], p3y = p3[1], p3z = p3[2], p3w = p3[3];

	function tb(ax, ay, az, aw, bx, by, bz, bw){
		var dx = bx - ax;
		var dy = by - ay;
		var dz = bz - az;
		var dw = bw - aw;
		return Math.max(Math.pow(dx * dx + dy * dy + dz * dz + dw * dw, tension * 0.5), 0.0001);
	}
	var t1 = tb(p0x, p0y, p0z, p0w, p1x, p1y, p1z, p1w);
	var t2 = tb(p1x, p1y, p1z, p1w, p2x, p2y, p2z, p2w) + t1;
	var t3 = tb(p2x, p2y, p2z, p2w, p3x, p3y, p3z, p3w) + t2;

	var dt21 = t2 - t1;
	var dt31 = t3 - t1;
	var dt32 = t3 - t2;

	var res = [];
	for (var s = 0; s < segs; s++){
		var t = dt21 * s / segs + t1;

		var dt1 = t1 - t;
		var dt2 = t2 - t;
		var dt3 = t3 - t;

		var A1x = (dt1 * p0x + t * p1x) / t1;
		var A1y = (dt1 * p0y + t * p1y) / t1;
		var A1z = (dt1 * p0z + t * p1z) / t1;
		var A1w = (dt1 * p0w + t * p1w) / t1;

		var A2x = (dt2 * p1x - dt1 * p2x) / dt21;
		var A2y = (dt2 * p1y - dt1 * p2y) / dt21;
		var A2z = (dt2 * p1z - dt1 * p2z) / dt21;
		var A2w = (dt2 * p1w - dt1 * p2w) / dt21;

		var A3x = (dt3 * p2x - dt2 * p3x) / dt32;
		var A3y = (dt3 * p2y - dt2 * p3y) / dt32;
		var A3z = (dt3 * p2z - dt2 * p3z) / dt32;
		var A3w = (dt3 * p2w - dt2 * p3w) / dt32;

		var B1x = (dt2 * A1x + t * A2x) / t2;
		var B1y = (dt2 * A1y + t * A2y) / t2;
		var B1z = (dt2 * A1z + t * A2z) / t2;
		var B1w = (dt2 * A1w + t * A2w) / t2;

		var B2x = (dt3 * A2x - dt1 * A3x) / dt31;
		var B2y = (dt3 * A2y - dt1 * A3y) / dt31;
		var B2z = (dt3 * A2z - dt1 * A3z) / dt31;
		var B2w = (dt3 * A2w - dt1 * A3w) / dt31;

		res.push([
			(dt2 * B1x - dt1 * B2x) / dt21,
			(dt2 * B1y - dt1 * B2y) / dt21,
			(dt2 * B1z - dt1 * B2z) / dt21,
			(dt2 * B1w - dt1 * B2w) / dt21
		]);
	}
	return res;
}

//
// convenience functions for generating splines using a list of points
//

var spline2d, spline3d, spline4d;
(function(){
	function makefunc(spline){
		// `points` is a list of points with the proper dimensions [[x, y], ...] or
		// [[x, y, z], ...] or [[x, y, z, w], ...]
		// `tension` is a number between 0 and 1 (0.5 is good default)
		// `segs` is the number of segments generated between two adjacent points
		// `close` is a boolean that indicates whether the returned curve should be a closed loop
		return function(points, tension, segs, close){
			var P = points;
			var L = P.length;

			if (L <= 0)
				return [];
			else if (L == 1)
				return [P[0].concat()];
			else if (L == 2){
				if (close){
					return spline(P[1], P[0], P[1], P[0], tension, segs).concat(
						spline(P[0], P[1], P[0], P[1], tension, segs));
				}
				return spline(P[0], P[0], P[1], P[1], tension, segs).concat([P[1].concat()]);
			}

			var res;
			if (close){
				res = spline(P[L - 1], P[0], P[1], P[2], tension, segs);
				if (L == 3){
					return res.concat(
						spline(P[0], P[1], P[2], P[0], tension, segs)).concat(
						spline(P[1], P[2], P[0], P[1], tension, segs));
				}
			}
			else{
				res = spline(P[0], P[0], P[1], P[2], tension, segs);
				if (L == 3){
					return res.concat(spline(P[0], P[1], P[2], P[2], tension, segs)).concat(
						[P[0].concat()]);
				}
			}

			for (var i = 0; i < L - 3; i++)
				res = res.concat(spline(P[i], P[i + 1], P[i + 2], P[i + 3], tension, segs));

			if (close){
				return res.concat(
					spline(P[L - 3], P[L - 2], P[L - 1], P[0], tension, segs)).concat(
					spline(P[L - 2], P[L - 1], P[0], P[1], tension, segs));
			}
			return res.concat(
				spline(P[L - 3], P[L - 2], P[L - 1], P[L - 1], tension, segs)).concat(
				P[L - 1].concat());
		};
	}
	spline2d = makefunc(spline2d_segment);
	spline3d = makefunc(spline3d_segment);
	spline4d = makefunc(spline4d_segment);
})();

if (typeof module !== 'undefined' && module.exports){
	// inside node.js, so export functions into global namespace
	global.spline2d_segment = spline2d_segment;
	global.spline3d_segment = spline3d_segment;
	global.spline4d_segment = spline4d_segment;
	global.spline2d = spline2d;
	global.spline3d = spline3d;
	global.spline4d = spline4d;
}
