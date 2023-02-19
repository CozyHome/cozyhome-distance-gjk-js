let GJK_E;
let SPLX;

const GJK_FSM = new FSM([{
	key:'init',
	setup:function(fsm,man) {},
	enter:function(prev,fsm,man) {
		fsm.cswitch(man, 'gjk');
	},
	exit:function(next,fsm,man) {},
	pulse:function(fsm,man) {}
},
{
	key:'gjk',
	setup:function(fsm,man) {},
	enter:function(prev,fsm,man) {
		man.cset = [RANDOM_CONVEX_POLYGON(5,192,192),
					RANDOM_CONVEX_POLYGON(5,192,192)];
	},
	exit:function(next,fsm,man) {},
	pulse:function(fsm,man) {
		noFill();
		translate(width/2,height/2);
		scale(1.5,1.5);

		strokeWeight(1);
		const A = man.cset[0];
		stroke(255); DRAW_POLYGON(A, "A");
		const B = man.cset[1];
		stroke(0,255,120); DRAW_POLYGON(B, "B");
		const cso = CSO_POLYGON(A,B);
		stroke(255,0,0);
		circle(0,0,16);

// Here is how you would typically go about executing the DGJK:
// run DGJK by passing the two point sets in. You can also choose
// to pass in a simplex from a previous frame to converge quicker for 
// this iteration. This is known as temporal coherency. Check js/dgjk.js for more details.
		const query = DGJK(A,B,SPLX);
// the final simplex after termination. This will be used with 'nv' to
// get a non-unique closest pair of points on each convex hull of the point sets.
		const splx 	= query.splx; 	/* type: GJKSimplex2D */
		SPLX = splx;
// the minimized vector in A - B returned by the DGJK.
		const nv 	= query.nv; 	/* type: vec2 */
// the norm of nv is the metric (euclidean) distance between the two
// convex sets.
		const dist 	= norm2(nv);	/* type: Number */
// the pairwise tuple of (non-unique) closest vertices that contribute
// to the minimized norm of A - B.
		const pair 	= DGJK_CLOSEST(splx,nv,A.l2w(),B.l2w()); /* type: {a:vec2, b:vec2 } */
// dereferencing the two points computed in DGJK_CLOSEST:
		const pa = pair.a; 
		const pb = pair.b;
		
		if(dist < 0.1) stroke(255,0,0);
		else stroke(255,255,0);
		
		DRAW_POLYGON(cso, "A - B");

		arrow2(pa,pb);

		stroke(0,255,255);
		circle(nv.x(),nv.y(),16);
		arrow2(nv,zero2());
		stroke(255,255,0);
		for(let i=0;i<splx.dim();i++) {
			const w = splx.peek(i);
			circle(w.ab().x(),w.ab().y(),16);
			const wv = MT(A.l2w(), w.a());
			const wv2 = MT(B.l2w(), w.b());
		}
		stroke(255);
		for(let j=0;j<splx.dim();j++) {
			const w = splx.peek(j);
			const v = splx.peek((1+j)%splx.dim());
			line2(w.ab(),v.ab());
		}
		const sp = 128;
		if(keyIsDown(65)) man.cset[0].mat.translate(-sp*deltaTime/1000,0);
		if(keyIsDown(68)) man.cset[0].mat.translate(sp*deltaTime/1000,0);
		if(keyIsDown(87)) man.cset[0].mat.translate(0,-sp*deltaTime/1000);
		if(keyIsDown(83)) man.cset[0].mat.translate(0,sp*deltaTime/1000);
		if(keyIsDown(81)) man.cset[0].mat.rot(0.5*deltaTime/1000);
		if(keyIsDown(69)) man.cset[0].mat.rot(-0.5*deltaTime/1000);
	}
}]);

function preload() {}
function setup() {
	const canvas = createCanvas(800,600);
	canvas.id("p5canvas");
	canvas.parent("#center_flexbox");

	GJK_E = CONSTRUCTOR_FSE(GJK_FSM);
}
function mousePressed() {
	const mv = new vec2(mouseX, mouseY);
//	GJK_E.man.onClick(mv);
}
function keyTyped() {
//	GJK_E.man.onKey(keyCode, new vec2(mouseX, mouseY));
}
function draw() {
	background(0);
	const mv = new vec2(mouseX, mouseY);
	noFill(); stroke(255); circle(0,0,32);
	GJK_E.fsm.pulse(GJK_E.man);
}
