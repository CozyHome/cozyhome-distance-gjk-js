// A Minkowski Vertex is a vertex containing the vertex
// of convex sets A, and B and their difference. These will all
// be used in calculating barycentric coordinates for closest
// point computations.
class MinkowskiVertex2D {
	#_a; #_b; #_ab;
	constructor(a,b,ab) {
		this.#_a = a;
		this.#_b = b;
		this.#_ab = ab;
	}
	a=()=>this.#_a;
	b=()=>this.#_b;
	ab=()=>this.#_ab;
	rebind=(ab)=>{this.#_ab = ab; }
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
	rebind=(m_a,m_b)=> { 
		for(let i=0;i<this.#_pts.length;i++) { 
			this.#_pts[i].rebind(
				sub2(MT(m_a,this.#_pts[i].a()),
					 MT(m_b,this.#_pts[i].b())
			));
		}
	}
}

// iterator for a DCEL (double connected edge list)
const ITERATE_DCEL=(vertex, yank, MAX=1000)=> {
	let dummy = vertex;
	let count = 0;
	if(!vertex) return;
	do {
		yank(dummy);
		dummy = dummy.next();
		count++;
		if(count > MAX) { 
			console.log("threshold exceeded");
			break;
		}
	}while(dummy != null && dummy != vertex);
}

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

const ITERATE_QUEUE=(queue, yoink)=> {
	if(!queue.head()) return; // is queue null?
	let dummy=queue.head().get_next();
	if(!dummy) return; // is next entry null?
	do{
		yoink(dummy.data());
		dummy=dummy.get_next();
	}while(dummy != null);
}

class QNode {
	#_prev; #_next; #_obj;
	constructor(obj) { this.#_obj = obj; }
	set_prev=(prev)=> { this.#_prev = prev; }
	set_next=(next)=> { this.#_next = next; }
	get_prev=()=> { return this.#_prev; }
	get_next=()=> { return this.#_next; }
	data=()=>{ return this.#_obj; }
}

class Queue {
	#_head; #_tail; #_count;
	constructor() { this.#_count = 0; }
	head=()=>{ return this.#_head; }
	push=(obj)=> {
		if(this.#_count > 0) {
			const next = new QNode(obj);
			next.set_prev(this.#_tail);
			this.#_tail.set_next(next);
			this.#_tail = next;
		}else {
			this.#_head = new QNode();
			this.#_tail = new QNode(obj);
			this.#_head.set_next(this.#_tail);
		}
		this.#_count++;
	}
	skip=(obj)=> {
		const next = new QNode(obj);
		const hn = this.#_head.get_next();
		this.#_head.set_next(next);
		next.set_prev(this.#_head);
		next.set_next(hn);
		if(hn) hn.set_prev(next);
		this.#_count++;
	}
	pop=()=> {
		if(this.#_count > 0) {
			const hn = this.#_head.get_next();
			this.#_head = hn;
			hn.set_prev(null);
			this.#_count--;
			return hn.data();
		}else {
			return null;
		}
	}
	peek=()=> {
		if(this.#_count > 0) return this.#_head.get_next().data();
		return null;
	}
	count=()=> { return this.#_count; }
	empty=()=> { return this.#_count <= 0; }
}