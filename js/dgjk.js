// input: two convex sets in the form of point clouds
// (requirements of input sets):
// 1. l2w() is a function returning the local to world matrix 
// transformation. It must be a column major 3x3 matrix!
//
// 2. l2w_iv() is the world to local matrix transformation. It
// must also be a column major 3x3 matrix. It MUST correspond
// to the inverse of whatever is returned in l2w().

// (optional) the last simplex iteration if frame
// coherence is high enough to warrant it.
// eps -> epsilon, keep default if you don't know what you're
// doing.
// bias -> play around with this if you'd like, I don't believe
// it would help with minimization but early termination may
// be a problem in certain configurations.
const DGJK=(sa,sb,splx=null,eps=0.00001,bias=1)=> {
	if(!splx) splx = new GJKSimplex2D();
	else {
// recompute differences (assuming they have not been reassigned to)
// if you manually reassign points in your polygon this simplex will not be
// updated! In place modifications will preserve this attachment.
		splx.rebind();
	}
	const m_a = sa.l2w();     // l2w matrix for A
	const m_b = sb.l2w(); 	  // l2w matrix for B 
	const mi_a = sa.l2w_iv(); // w2l matrix for A
	const mi_b = sb.l2w_iv(); // w2l matrix for B
// helper function returning the translational component of the last column
// in a 3x3
	const origin=(mat)=> {
		const w = mat[8];
		return new vec2(w*mat[6],w*mat[7]);
	}
// starting point should exist in convex difference
	const nv = sub2(origin(m_a), origin(m_b));
	const d = copy2(nv);
	let i = 0;
	do {
// flip next vertex to maximize in opposite direction
		d._x = -nv._x; d._y = -nv._y;
		const a = FAST_SUPPORT(d, sa, m_a, mi_a);
		const b = FAST_SUPPORT(nv, sb, m_b, mi_b);
		const w = new MinkowskiVertex2D(a,b);
// floating point heuristic from the help of Gino Van's analysis
		const v1 = nv._x*nv._x + nv._y*nv._y;
		if(v1 - dot2(w.ab(),nv) < eps * v1) {
			return { splx, nv };
		}else {
			splx.push(w);
			MUTATE_SIMPLEX(splx,nv);
		}
	} while(i++ < bias*(sa.hull.length + sb.hull.length) && splx.dim() < 3);
	return { splx, nv };
}
// designed to run on the finishing simplex of the DGJK. Computes the
// barycentric decomposition of a point in standard coordinates.
// this is super easy in 2D as we only need to deal with vector projections.
const DGJK_CLOSEST=(splx,nv)=> {
	switch(splx.dim()) {
		default:
			return {a:zero2(),b:zero2()};
		case 1:
			const p1 = splx.peek(0);
			return {a:p1.a(),b:p1.b()};
		case 2:
			return BARY_2(splx,nv);
	}
}
// barycentric coordinates for one dimensional
// simplices (its as easy as a simple vector projection)
const BARY_2=(splx,nv)=> {
	const p1 = splx.peek(0);
	const p2 = splx.peek(1);
	const ab = sub2(p2.ab(),p1.ab());
	const t = dot2(ab,sub2(nv,p1.ab())) / dot2(ab,ab);
	return {
		a:lerp2(p1.a(),p2.a(),t),
		b:lerp2(p1.b(),p2.b(),t)
	};
}
// mutate the simplex to eliminate any facets not
// contributing to a convex combination of nv.
const MUTATE_SIMPLEX=(splx,nv)=> {
	switch(splx.dim()) {
		case 1: // convex, simple point (0D)
			MUTATE_0D(nv,splx);
		break;
		case 2: // convex, simple segment (1D)
			MUTATE_1D(nv,splx);
		break;
		case 3: // convex, simple area (2D)
			MUTATE_2D(nv,splx);
		break;
	}
}
// responsible for handling simplex state for 
// a zero dimensional simplex (point)
const MUTATE_0D=(nv,splx)=> {
	const a = splx.peek(0).ab();
	nv._x = a._x;
	nv._y = a._y;
}
// responsible for handling simplex state for
// a one dimensional simplex (convex segment)
const MUTATE_1D=(nv,splx)=> {
	const a = splx.peek(0);
	const b = splx.peek(1);
	const rep = VOROCENTRIC_1D({a:a.ab(), b:b.ab()}, {a:0x1, b:0x2});
	nv._x = rep.v._x;
	nv._y = rep.v._y;
	switch(rep.r) {
		case 0x1: // A
			const a = splx.peek(0);
			splx.clear();
			splx.push(a);
		break;
	}
}
// responsible for handling simplex state for
// a two dimensional simplex (convex triangle)
const MUTATE_2D=(nv,splx)=> {
	const a = splx.peek(0);
	const b = splx.peek(1);
	const c = splx.peek(2);
	const rep = VOROCENTRIC_2D({a:a.ab(),b:b.ab(),c:c.ab()},{a:0x1,b:0x2,c:0x4});
	nv._x = rep.v._x; nv._y = rep.v._y;
	switch(rep.r) {
		case 0x1: splx.clear(); splx.push(a);
		break;
		case 0x3: splx.clear(); splx.push(b); splx.push(a);
		break;
		case 0x5: splx.clear(); splx.push(c); splx.push(a);
		break;
	}
}
// Barycentric/Voronoi decomposition for a 
// line segment in R^2.
const VOROCENTRIC_1D=(edge,bits)=> {
	const ao = mul2(-1, edge.a);
	const bo = mul2(-1, edge.b);
	const ab = sub2(edge.b, edge.a);

	const abm = norm2(ab);
	if(abm > 0) {
		ab._x /= abm;
		ab._y /= abm;
	}

	const v = dot2(ao, ab);
	if(v > abm) {
		return {r:bits.b, v:edge.b };
	}else if(v < 0) {
		return {r:bits.a, v:edge.a};
	}else {
		return {r:bits.a | bits.b,
			v:new vec2(
				edge.a._x + ab._x*v,
				edge.a._y + ab._y*v
			),
		}
	}
}
// Barycentric/Voronoi decomposition for a
// triangle in R^2.
const VOROCENTRIC_2D=(tri,bits)=> {
	const o = zero2();
	const a = tri.a;
	const b = tri.b;
	const c = tri.c; 
	const same=(v,w)=> {
		return dot2(v,w) > 0;
	}
// computes the signed orthogonal projections onto each
// combinational segment on the simplex
	const det = DET_2D(sub2(a,c),sub2(a,b));
	const signed_areas=(det)=> {
		let nflags = 0;
		nflags |= ORIENT_2D(a,b,o,det) ? 0x1 : 0; // BA
		nflags |= ORIENT_2D(b,c,o,det) ? 0x2 : 0; // CB
		nflags |= ORIENT_2D(c,a,o,det) ? 0x4 : 0; // AC 
		return nflags;
	}
	const dual_edges=(pm,bits)=> {
		const xy = sub2(pm.y, pm.x);
		const zy = sub2(pm.y, pm.z);

		const xp = mul2(-1, pm.x);
		const yp = mul2(-1, pm.y);
		const zp = mul2(-1, pm.z);
// REGION Y
		if(same(yp,zy) && same(yp, xy)) {
			return {r:bits.y, v:pm.y};
		}
		if(!same(yp,xy)) {
			if(!same(xp,xy)) {
// REGION X
				return {r:bits.x, v:pm.x};
			}else {
// REGION XY
				return {r:bits.x|bits.y, 
					v:add2(onto2c(xp, xy), pm.x),
				};
			}
		}
		if(!same(yp, zy)) {
// REGION Z
			if(!same(zp,zy)) {
				return {r:bits.z, v:pm.z};
			}else {
// REGION ZY
				return {r:bits.y|bits.z,
					v:add2(onto2c(zp, zy), pm.z)
				};
			}
		}
		return {r:bits.x|bits.y|bits.z, v:zero2()};
	}
	switch(signed_areas(det)) {
		case 1: // AB-edge
			return VOROCENTRIC_1D(
				{a:tri.a, b:tri.b}, {a:bits.a,b:bits.b});
		case 2: // CB-edge
			return VOROCENTRIC_1D(
				{a:tri.c,b:tri.b}, {a:bits.c,b:bits.b});
		case 3: // ABC-chain
			return dual_edges(
					{x:tri.a,y:tri.b,z:tri.c},{x:bits.a,y:bits.b,z:bits.c}
				);
		case 4: // CA-chain
			return VOROCENTRIC_1D(
				{a:tri.c,b:tri.a},{a:bits.c,b:bits.a}
			);
		case 5: // CAB-chain
			return dual_edges(
				{x:tri.c,y:tri.a,z:tri.b}, {x:bits.c,y:bits.a,z:bits.b}
			);
		case 6: // BCA-chain
			return dual_edges(
				{x:tri.b,y:tri.c,z:tri.a}, {x:bits.b,y:bits.c,z:bits.a}
			);
		default: // inscribed => intersection. will recurse @ top
			return {r:bits.a|bits.b|bits.c, v:zero2()};
	}
}
// determinant of 2x2
const DET_2D=(v,w)=> { return v.x()*w.y() - v.y()*w.x(); }
// sign test for edge (a,b) and point c, with sign s used
// if you want consistent reports for unknown winding orders.
// leave as is if you do not know what this does.
const ORIENT_2D=(a,b,c,s=1)=> {
	const ac   = sub2(c,a);
	const p_ab = perp2(sub2(b,a));
	return s*dot2(p_ab, ac) > 0;
}
// optimize later, focus on getting the support function
// to run, work on greedy adjacency approach using DCEL
// using simulation coherency, this can be a near constant
// time operation if correctly implemented. As of now, we
// naively search the whole set for its maximum.
const FAST_SUPPORT=(v,set,mat,inv)=> {
	const hull = set.hull;
	let mv = hull[0];
	let max = Number.NEGATIVE_INFINITY;
// invert direction w.r.t coordinate frame
	tv = mTransform(inv, [v._x,v._y,0]);
	for(let i=0;i<hull.length;i++) {
		const dot = hull[i]._x*tv[0] + hull[i]._y*tv[1];
		if(dot >= max) {
			max = dot;
			mv = hull[i];
		}
	}
	mp = mTransform(mat, [mv._x, mv._y, 1]);
	return new vec2(mp[0],mp[1]);
}