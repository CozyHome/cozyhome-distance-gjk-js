// CREDITS: Daniel J. Cucuzza
// DATE: February 19th, 2023
// You can contact me at gaunletgames@gmail.com if you have
// any questions about the implementation or if you notice
// any errors.

// entity for housing state & data
let GJK_E;
// state machine for program
const GJK_FSM=new FSM([{
	key:'init',
	setup:function(fsm,man) {},
	enter:function(prev,fsm,man) {
		fsm.cswitch(man,'gjk');
	},
	exit:function(next,fsm,man) {},
	pulse:function(fsm,man) {}

},
{
	key:'gjk',
	setup:function(fsm,man) {},
	enter:function(prev,fsm,man) {
		man.cset = [RANDOM_CONVEX_POLYGON(12,192,192),
			RANDOM_CONVEX_POLYGON(12,192,192)
		];
	},
	exit:function(next,fsm,man) {},
	pulse:function(fsm,man) {
		noFill();
		translate(width/2,height/2);
		const mv = new vec2(mouseX - width/2, mouseY - height/2);

		const A = man.cset[0];
		stroke(255); DRAW_POLYGON(A, "A");
		const B = man.cset[1];
		stroke(0,255,120); DRAW_POLYGON(B, "B");
		const cso = CSO_POLYGON(A,B);
		stroke(255,255,0);
		circle(0,0,16);
		
		const dv = sub2(mv, A.origin());
		draw2p(A.origin(),dv);

		DRAW_POLYGON(cso, "A - B");

		const cast = CAST_DGJK(A,dv,B);
		const t = cast.toi;

		const n_mt = A.l2w().slice();
		if(t > 0) {
			n_mt[6] += t*dv._x;
			n_mt[7] += t*dv._y;

			DRAW_POLYGON(A,"A",n_mt);
			const p = DGJK_CLOSEST(cast.query.splx,cast.query.nv,n_mt,B.l2w());
			circle(p.b.x(),p.b.y(),32);
			const n = DGJK_NORMAL(cast.query.splx,dv,A.l2w(),B.l2w());
			stroke(255,0,0);
			draw2p(p.a, mul2(100,n));
		}else {
			n_mt[6] += dv._x;
			n_mt[7] += dv._y;

			DRAW_POLYGON(A,"A",n_mt);
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
//	DEL_E.man.onClick(mv);
}
function keyTyped() {
//	DEL_E.man.onKey(keyCode, new vec2(mouseX, mouseY));
}
function draw() {
	background(0);
//	const mv = new vec2(mouseX, mouseY);
	noFill(); stroke(255); circle(0,0,32);
	GJK_E.fsm.pulse(GJK_E.man);
}
