// GEOMETRIC PRIMITIVES: {
// 		class Vertex2D 	 := DCEL node wrapper for polygonal/edge structures in R^2.
// 		class MinkowskiVertex2D := a three R^2 tuple containing, { A, B, A - B } for use in DGJK(...)
// 		class GJKSimplex2D := iteratively constructed simplex in R^2 for minimization in DGJK(...)
// }
class Vertex2D {
	#_point; #_next; #_prev;
	constructor(point) { this.#_point = point; this.#_next = null; this.#_prev = null; };
	bind_prev=(prev)=> { this.#_prev = prev; } // set back pointer
	bind_next=(next)=> { this.#_next = next; } // set fwds pointer
// inserts a new vertex in between us and our next door neighbor
	insert_next=(middle)=> {
		const next = this.#_next;
// connect ourselves to middle
		this.bind_next(middle);
		middle.bind_prev(this);
// if next is not null
		if(next) next.bind_prev(middle);
	}
// inserts a new vertex in between us and our previous door neighbor
	insert_prev=(middle)=> {
		const prev = this.#_prev;
		this.bind_prev(middle);
		middle.bind_next(this);

		if(prev) prev.bind_next(middle);
	}
// removes itself and stitches its' neighbors together
	dissolve=()=> {
		const prev = this.prev();
		const next = this.next();
		if(prev) prev.bind_next(next);
		if(next) next.bind_prev(prev);
	}
	next=()=>{ return this.#_next; }
	prev=()=>{ return this.#_prev; }
	point=()=>{ return this.#_point; } 
}

// A Minkowski Vertex is a vertex containing the vertex
// of convex sets A, and B and their difference. These will all
// be used in calculating barycentric coordinates for closest
// point computations.
class MinkowskiVertex2D {
	#_a; #_b; #_ab;
	constructor(a,b) {
		this.#_a = a;
		this.#_b = b;
		this.#_ab = sub2(a,b);
	}
	a=()=>this.#_a;
	b=()=>this.#_b;
	ab=()=>this.#_ab;
	rebind=()=>{this.#_ab = sub2(this.#_a,this.#_b);}
}

// d-simplex restricted to two dimensions for the GJK algorithm. This will
// function similar to that of a stack as the elements will be reordered.
class GJKSimplex2D {
	#_pts; #_dim;
	constructor() {
// dimensionality of simplex
		this.#_dim = 0;
// pseudo-stack arraylist
		this.#_pts = [];
	}
// requires support points a and b
	push=(w)=> {
		if(this.#_dim < 3) {
			this.#_pts.push(w);
			this.#_dim++;
// swap
			for(let i=this.#_dim-1;i>0;i--) {
				const c = this.#_pts[i];
				const p = this.#_pts[i-1];
				this.#_pts[i-1]=c;
				this.#_pts[i]=p;
			}
		}else { return; }
	}
	peek=(i)=> {
		if(i > this.#_dim) { return null; }
		else { return this.#_pts[i]; }
	}
	pop=()=> {
		if(this.#_dim > 0) { return this.#_pts[this.#_dim--]; }
		else return null;
	}
	clear=()=> {
		this.#_dim = 0;
		while(this.#_pts.length > 0) this.#_pts.pop();
	}
	dim=()=> { return this.#_dim; }
	rebind=()=> { for(let i=0;i<this.#_pts.length;i++) this.#_pts[i].rebind(); }
}
