// returns randomized DCEL and cloud representation for convex
// polygon in R^2.
const RANDOM_CONVEX_POLYGON=(n=7,w=64,h=64)=> {
// randomized point in R^2 scaled w.r.t AABB
	const rpoint=(w,h)=> { 
		return new vec2(w*random(),h*random());
	}
	const pts=[];
	for(let i=0;i<n;i++) {
		pts.push(rpoint(w,h));
	}
	return CONVEX_POLYGON(pts);
}

const DRAW_POLYGON=(gon,name="",mt=null)=> {
	const m = mt == null ? gon.l2w() : mt;
	ITERATE_DCEL(gon.dcel, (vert)=> {
		const next = vert.next();
		if(next) {
			arrow2(MT(m,vert.point()), MT(m,next.point()),10,30);
		}
	});
	text(name,...gon.origin().f2d());
}

// two convex sets subtracted
const CSO_POLYGON=(s1,s2)=> {
	const cso = [];
	const pts1 = s1.hull;
	const pts2 = s2.hull;
	const m1 = s1.l2w();
	const m2 = s2.l2w();
	let p1 = [0,0];
	let p2 = [0,0];
	for(let i=0;i<pts1.length;i++) {
		p1[0] = pts1[i]._x;
		p1[1] = pts1[i]._y;
		p1 = mTransform(m1, p1);
		for(let j=0;j<pts2.length;j++) {
			p2[0] = pts2[j]._x;
			p2[1] = pts2[j]._y;
			p2 = mTransform(m2, p2);
			cso.push(new vec2(p1[0]-p2[0],p1[1]-p2[1]));
		}	
	}
	return CONVEX_POLYGON(cso);
}

// This is an example generating function you could use
// to have objects compatible with DGJK:
const CONVEX_POLYGON=(pts)=>{
// get the convex hull of an arbitary point set in R^2
	const gon = QUICKHULL(pts);
	const hull = gon.hull;

// compute centroid to then use as origin for object space
// representation
	const centroid = zero2();
	for(let i=0;i<hull.length;i++) {
		centroid._x += hull[i]._x;
		centroid._y += hull[i]._y;
	}
// centroid is a uniform distribution of weight along
// all contributing points
	centroid._x /= hull.length;
	centroid._y /= hull.length;		
// reconfigure polygon to be w.r.t the centroid of its 
// vertices. This essentially allows us to apply linear transformations
// from the viewpoint of the polygon's centroid.
	for(let i=0;i<hull.length;i++) {
		hull[i]._x = hull[i]._x - centroid._x;
		hull[i]._y = hull[i]._y - centroid._y;
	}

	gon.mat = new Matrix3x3();
	gon.mat.translate(centroid._x, centroid._y);		
// set up the matrix dependencies for the polygon
// l2w() and l2w_iv() are REQUIRED for DGJK to function. These must
// be part of your input set to function.
	gon.l2w =()=> { return gon.mat.get(); }				// local to world matrix accessor
	gon.l2w_iv=()=> { return mInverse(gon.l2w()); }		// world to local matrix accessor
	gon.pts =()=> { return gon.hull; }					// point set accessor
// helper method for center calculation (not necessarily needed for DGJK to run)
	gon.origin=()=> {
		const m = gon.l2w();
		const w = m[8];
		return new vec2(w*m[6],w*m[7]);
	}
	return gon;
}

// linear least squares for a point set in R^2.
const LEAST_SQR=(pts)=> {
	const n = pts.length;
	let spw=0;
	for(let i=0;i<n;i++) { spw = spw + pts[i].x()*pts[i].y(); }
	let sx=0; let sy=0;
	for(let i=0;i<n;i++) { sx += pts[i].x(); sy += pts[i].y(); }
	let ssum=0;
	for(let i=0;i<n;i++) { ssum += pts[i].x()*pts[i].x(); }
	let sums=0;
	for(let i=0;i<n;i++) { sums += pts[i].x(); }
	sums *= sums;
// y = ax + b
	const a = (n*spw - sx*sy) / (n*ssum - sums);
	const b = (sy - a*sx) / n;
	return {a,b};
}

// assumes general position coordinate frame for input
const SUPPORT=(p,v,pts)=> {
	const n = pts.length;
	let idx = -1; let mdp = Number.NEGATIVE_INFINITY;
	for(let i=0;i<n;i++) {
		const dif = sub2(pts[i], p);
		const dot = dot2(dif, v);
		if(dot > mdp) {
			idx = i;
			mdp = dot;
		}
	} return {i:idx,d:mdp};
}

// assumes general position coordinate frame for input
const USUPPORT=(p,v,pts)=> {
	const n = pts.length;
	let idx = -1; let mdp = 0;
	for(let i=0;i<n;i++) {
		const dif = sub2(pts[i], p);
		let dot = dot2(dif, v);
		dot = dot > 0 ? dot : -dot;
		if(dot > mdp) {
			idx = i;
			mdp = dot;
		}
	} return {i:idx,d:mdp};
}

// construct a convex hull out of a general point set in R^2. input is a two dimensional
// array of vec2s.
const QUICKHULL=(set)=> {
	if(!set) return { error:true, msg:'set was null', hull:null, dcel:null }; // set was null
	const lst = LEAST_SQR(set);
	const org = new vec2(0, lst.b);
	const dif = new vec2(1, lst.a);
// maximize along linear least squares differential
	const s1 = SUPPORT(org, dif, set);
// maximize along the negative of the differential
	const s2 = SUPPORT(org, mul2(-1,dif), set);
	if(s1.i < 0 || s2.i < 0) { 
		return {  error:true, msg:'least squares failed to maximize differential', hull:null, dcel:null };
	}
	const p1 = set[s1.i];
	const p2 = set[s2.i];

	const perp = perp2(sub2(p2,p1));
	const mid = lerp2(p1,p2,0.5);
// maximize along the differential's orthogonal complement
	const s3 = USUPPORT(mid, perp, set);
	if(s3.i < 0) {
		return { error:true, msg:'least squares failed to maximize orthogonal', hull:[p1,p2], dcel:null };
	}
	const p3 = set[s3.i];
	const A = new Vertex2D(p1);
	const B = new Vertex2D(p2);
	const C = new Vertex2D(p3);
	if(!ORIENT_2D(A.point(), B.point(), C.point())) {
		A.insert_next(B); B.insert_next(C); C.insert_next(A);
	}else {
		C.insert_prev(A); B.insert_prev(C); A.insert_prev(B);
	}
// remove all interior elements
	for(let i=set.length-1;i>=0;i--) {
		const p = set[i];
		let dummy = A;
		do {
			const sp = dummy.point();
			const next = dummy.next();
			if(next) {
				const spx = next.point();
				if(ORIENT_2D(sp, spx, p)) {
					dummy=null; break;
				}
			}
			dummy = next;
		}while(dummy != null && dummy != A);
// reached end of hull, must be contained!
		if(dummy == A) {
			set[i] = set[set.length-1];
			set.pop();
		}
	}
	const queue = new Queue();
	queue.push(A); queue.push(B); queue.push(C);

	while(!queue.empty() && set.length > 0) {
		const vert = queue.pop();
		const next = vert.next();
		if(!next) { 
			return { error:true, 'msg':'DCEL not connected', hull:[p1,p2,p3], dcel:A }
		}

		const perp    = perp2(sub2(next.point(), vert.point()));
		const maximal = SUPPORT(vert.point(), perp, set);
		if(maximal.i >= 0 && maximal.d > 0) {
			const sup = set[maximal.i];
			const nv = new Vertex2D(sup);
			let cw = vert;
			do {
				if(ORIENT_2D(cw.point(), cw.prev().point(), sup)) break;
				cw = cw.prev();
			} while(cw != vert);
			let ccw = vert;
			do {
				if(ORIENT_2D(ccw.next().point(), ccw.point(), sup)) break;
				ccw = ccw.next();
			} while(cw != next);

			for(let i=set.length-1;i>=0;i--) {
				const pt = set[i];
				if(ORIENT_2D(ccw.point(), cw.point(), pt)) continue;
				if(ORIENT_2D(cw.point(), sup, pt)) continue;
				if(ORIENT_2D(sup, ccw.point(), pt)) continue;
				set[i]=set[set.length-1]; set.pop();
			}
// clockwise and counterclockwise horizon forms a diagonal along the set.
// We'll connect our new vertex to these guys, the GC will most likely
// recognize that these vertices in the other chains no longer have any
// references.
			cw.bind_next(nv);
			nv.bind_prev(cw);
			ccw.bind_prev(nv);
			nv.bind_next(ccw);
// notify hull that new edges must be detected. We'll check these
// last as we may have heavily converged in this direction. This way,
// less points will be checked for each incoming iteration.
			queue.push(cw);
			queue.push(nv);
		}
	}
// construct hull from DCEL
	const hull = [];
	let dummy = A;
	do {
		hull.push(dummy.point());
		dummy = dummy.next();
	}while(dummy != A);
	return { error:false, msg:'', dcel:A, hull:hull };
}