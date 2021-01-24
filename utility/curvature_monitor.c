/**************************
 * Monitor.c
 * Generated by KeYmaera X
 **************************/

#include <math.h>
#include <stdbool.h>

typedef struct plantparameters {
  long double A;
  long double B;
  long double T;
  long double a;
  long double eps;
  long double k;
  long double vh;
  long double vl;
} plantparameters;

typedef struct plantstate {
  long double t;
  long double t_0;
  long double v;
  long double v_0;
  long double xg;
  long double xg_0;
  long double yg;
  long double yg_0;
} plantstate;

typedef struct parameters {
  long double A;
  long double B;
  long double T;
  long double eps;
} parameters;

typedef struct state {
  long double a;
  long double k;
  long double t;
  long double v;
  long double vh;
  long double vl;
  long double xg;
  long double yg;
} state;


typedef struct input input;

typedef struct verdict { int id; long double val; } verdict;

/* Computes distance to safety boundary on prior and current state (>=0 is safe, <0 is unsafe) */
verdict plantBoundaryDist(plantstate pre, plantstate curr, const plantparameters* const params) {
  if (pre.v >= 0.0L) {
if (pre.t <= params->T) {
if (((pre.t >= 0.0L) && (((((params->k)*((params->eps)*(params->eps)))-((200.0L)*(params->eps)))*(100.0L) < ((params->k)*(((pre.xg)*(pre.xg))+((pre.yg)*(pre.yg))))-((((2.0L)*(pre.xg))*(100.0L))*(10.0L))) && (((params->k)*(((pre.xg)*(pre.xg))+((pre.yg)*(pre.yg))))-((((2.0L)*(pre.xg))*(100.0L))*(10.0L)) < (((params->k)*((params->eps)*(params->eps)))+((200.0L)*(params->eps)))*(100.0L)))) || ((((10.0L)*(pre.v))+((params->a)*((params->T)-(pre.t))) >= 0.0L) && ((((params->a >= 0.0L) && (((10.0L)*(pre.v))+((params->a)*((params->T)-(pre.t))) <= (10.0L)*(params->vh))) || (((params->a <= 0.0L) && (pre.v <= params->vh)) || (((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->B)*(((((2.0L)*(pre.v))*((params->T)-(pre.t)))*(10.0L))+((params->a)*(((params->T)-(pre.t))*((params->T)-(pre.t))))))+(((((pre.v)*(10.0L))+((params->a)*((params->T)-(pre.t))))*(((pre.v)*(10.0L))+((params->a)*((params->T)-(pre.t)))))-(((10.0L)*(params->vh))*((10.0L)*(params->vh))))) <= ((((2.0L)*(params->B))*((pre.yg)-((10.0L)*(params->eps))))*(10000.0L))*(100.0L)) || ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->B)*(((((2.0L)*(pre.v))*((params->T)-(pre.t)))*(10.0L))+((params->a)*(((params->T)-(pre.t))*((params->T)-(pre.t))))))+(((((pre.v)*(10.0L))+((params->a)*((params->T)-(pre.t))))*(((pre.v)*(10.0L))+((params->a)*((params->T)-(pre.t)))))-(((10.0L)*(params->vh))*((10.0L)*(params->vh))))) <= ((((2.0L)*(params->B))*((fabsl(pre.xg))-((10.0L)*(params->eps))))*(10000.0L))*(100.0L))))) && (((params->a >= 0.0L) && (pre.v >= params->vl)) || (((params->a <= 0.0L) && (((10.0L)*(pre.v))+((params->a)*((params->T)-(pre.t))) >= (10.0L)*(params->vl))) || (((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->A)*(((((2.0L)*(pre.v))*((params->T)-(pre.t)))*(10.0L))+((params->a)*(((params->T)-(pre.t))*((params->T)-(pre.t))))))+((((params->vl)*(10.0L))*((params->vl)*(10.0L)))-((((pre.v)*(10.0L))+((params->a)*((params->T)-(pre.t))))*(((pre.v)*(10.0L))+((params->a)*((params->T)-(pre.t))))))) <= ((((2.0L)*(params->A))*((pre.yg)-((10.0L)*(params->eps))))*(10000.0L))*(100.0L)) || ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->A)*(((((2.0L)*(pre.v))*((params->T)-(pre.t)))*(10.0L))+((params->a)*(((params->T)-(pre.t))*((params->T)-(pre.t))))))+((((params->vl)*(10.0L))*((params->vl)*(10.0L)))-((((pre.v)*(10.0L))+((params->a)*((params->T)-(pre.t))))*(((pre.v)*(10.0L))+((params->a)*((params->T)-(pre.t))))))) <= ((((2.0L)*(params->A))*((fabsl(pre.xg))-((10.0L)*(params->eps))))*(10000.0L))*(100.0L)))))))) {
if (curr.v >= 0.0L) {
if (curr.t <= params->T) {
if (((curr.t >= 0.0L) && (((((params->k)*((params->eps)*(params->eps)))-((200.0L)*(params->eps)))*(100.0L) < ((params->k)*(((curr.xg)*(curr.xg))+((curr.yg)*(curr.yg))))-((((2.0L)*(curr.xg))*(100.0L))*(10.0L))) && (((params->k)*(((curr.xg)*(curr.xg))+((curr.yg)*(curr.yg))))-((((2.0L)*(curr.xg))*(100.0L))*(10.0L)) < (((params->k)*((params->eps)*(params->eps)))+((200.0L)*(params->eps)))*(100.0L)))) || ((((10.0L)*(curr.v))+((params->a)*((params->T)-(curr.t))) >= 0.0L) && ((((params->a >= 0.0L) && (((10.0L)*(curr.v))+((params->a)*((params->T)-(curr.t))) <= (10.0L)*(params->vh))) || (((params->a <= 0.0L) && (curr.v <= params->vh)) || (((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->B)*(((((2.0L)*(curr.v))*((params->T)-(curr.t)))*(10.0L))+((params->a)*(((params->T)-(curr.t))*((params->T)-(curr.t))))))+(((((curr.v)*(10.0L))+((params->a)*((params->T)-(curr.t))))*(((curr.v)*(10.0L))+((params->a)*((params->T)-(curr.t)))))-(((10.0L)*(params->vh))*((10.0L)*(params->vh))))) <= ((((2.0L)*(params->B))*((curr.yg)-((10.0L)*(params->eps))))*(10000.0L))*(100.0L)) || ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->B)*(((((2.0L)*(curr.v))*((params->T)-(curr.t)))*(10.0L))+((params->a)*(((params->T)-(curr.t))*((params->T)-(curr.t))))))+(((((curr.v)*(10.0L))+((params->a)*((params->T)-(curr.t))))*(((curr.v)*(10.0L))+((params->a)*((params->T)-(curr.t)))))-(((10.0L)*(params->vh))*((10.0L)*(params->vh))))) <= ((((2.0L)*(params->B))*((fabsl(curr.xg))-((10.0L)*(params->eps))))*(10000.0L))*(100.0L))))) && (((params->a >= 0.0L) && (curr.v >= params->vl)) || (((params->a <= 0.0L) && (((10.0L)*(curr.v))+((params->a)*((params->T)-(curr.t))) >= (10.0L)*(params->vl))) || (((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->A)*(((((2.0L)*(curr.v))*((params->T)-(curr.t)))*(10.0L))+((params->a)*(((params->T)-(curr.t))*((params->T)-(curr.t))))))+((((params->vl)*(10.0L))*((params->vl)*(10.0L)))-((((curr.v)*(10.0L))+((params->a)*((params->T)-(curr.t))))*(((curr.v)*(10.0L))+((params->a)*((params->T)-(curr.t))))))) <= ((((2.0L)*(params->A))*((curr.yg)-((10.0L)*(params->eps))))*(10000.0L))*(100.0L)) || ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->A)*(((((2.0L)*(curr.v))*((params->T)-(curr.t)))*(10.0L))+((params->a)*(((params->T)-(curr.t))*((params->T)-(curr.t))))))+((((params->vl)*(10.0L))*((params->vl)*(10.0L)))-((((curr.v)*(10.0L))+((params->a)*((params->T)-(curr.t))))*(((curr.v)*(10.0L))+((params->a)*((params->T)-(curr.t))))))) <= ((((2.0L)*(params->A))*((fabsl(curr.xg))-((10.0L)*(params->eps))))*(10000.0L))*(100.0L)))))))) {
if (curr.t_0 == pre.t) {
if (curr.v_0 == pre.v) {
if (curr.xg_0 == pre.xg) {
if (curr.yg_0 == pre.yg) {
verdict result = { .id=1, .val=((((((0.0L)+(-((0.0L)-(pre.v))))+(-((pre.t)-(params->T))))+(-(fminl(fmaxl((0.0L)-(pre.t), fmaxl(((((params->k)*((params->eps)*(params->eps)))+(-((200.0L)*(params->eps))))*(100.0L))-(((params->k)*(((pre.xg)*(pre.xg))+((pre.yg)*(pre.yg))))+(-((((2.0L)*(pre.xg))*(100.0L))*(10.0L)))), (((params->k)*(((pre.xg)*(pre.xg))+((pre.yg)*(pre.yg))))+(-((((2.0L)*(pre.xg))*(100.0L))*(10.0L))))-((((params->k)*((params->eps)*(params->eps)))+((200.0L)*(params->eps)))*(100.0L)))), fmaxl((0.0L)-(((10.0L)*(pre.v))+((params->a)*((params->T)+(-(pre.t))))), fmaxl(fminl(fmaxl((0.0L)-(params->a), (((10.0L)*(pre.v))+((params->a)*((params->T)+(-(pre.t)))))-((10.0L)*(params->vh))), fminl(fmaxl(params->a, (pre.v)-(params->vh)), fminl(((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->B)*(((((2.0L)*(pre.v))*((params->T)+(-(pre.t))))*(10.0L))+((params->a)*(((params->T)+(-(pre.t)))*((params->T)+(-(pre.t)))))))+(((((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t)))))*(((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t))))))+(-(((10.0L)*(params->vh))*((10.0L)*(params->vh)))))))-(((((2.0L)*(params->B))*((pre.yg)+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)), ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->B)*(((((2.0L)*(pre.v))*((params->T)+(-(pre.t))))*(10.0L))+((params->a)*(((params->T)+(-(pre.t)))*((params->T)+(-(pre.t)))))))+(((((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t)))))*(((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t))))))+(-(((10.0L)*(params->vh))*((10.0L)*(params->vh)))))))-(((((2.0L)*(params->B))*((fabsl(pre.xg))+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L))))), fminl(fmaxl((0.0L)-(params->a), (params->vl)-(pre.v)), fminl(fmaxl(params->a, ((10.0L)*(params->vl))-(((10.0L)*(pre.v))+((params->a)*((params->T)+(-(pre.t)))))), fminl(((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->A)*(((((2.0L)*(pre.v))*((params->T)+(-(pre.t))))*(10.0L))+((params->a)*(((params->T)+(-(pre.t)))*((params->T)+(-(pre.t)))))))+((((params->vl)*(10.0L))*((params->vl)*(10.0L)))+(-((((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t)))))*(((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t))))))))))-(((((2.0L)*(params->A))*((pre.yg)+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)), ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->A)*(((((2.0L)*(pre.v))*((params->T)+(-(pre.t))))*(10.0L))+((params->a)*(((params->T)+(-(pre.t)))*((params->T)+(-(pre.t)))))))+((((params->vl)*(10.0L))*((params->vl)*(10.0L)))+(-((((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t)))))*(((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t))))))))))-(((((2.0L)*(params->A))*((fabsl(pre.xg))+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)))))))))))+(-((0.0L)-(curr.v))))+(-((curr.t)-(params->T))))+(-(fminl(fmaxl((0.0L)-(curr.t), fmaxl(((((params->k)*((params->eps)*(params->eps)))+(-((200.0L)*(params->eps))))*(100.0L))-(((params->k)*(((curr.xg)*(curr.xg))+((curr.yg)*(curr.yg))))+(-((((2.0L)*(curr.xg))*(100.0L))*(10.0L)))), (((params->k)*(((curr.xg)*(curr.xg))+((curr.yg)*(curr.yg))))+(-((((2.0L)*(curr.xg))*(100.0L))*(10.0L))))-((((params->k)*((params->eps)*(params->eps)))+((200.0L)*(params->eps)))*(100.0L)))), fmaxl((0.0L)-(((10.0L)*(curr.v))+((params->a)*((params->T)+(-(curr.t))))), fmaxl(fminl(fmaxl((0.0L)-(params->a), (((10.0L)*(curr.v))+((params->a)*((params->T)+(-(curr.t)))))-((10.0L)*(params->vh))), fminl(fmaxl(params->a, (curr.v)-(params->vh)), fminl(((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->B)*(((((2.0L)*(curr.v))*((params->T)+(-(curr.t))))*(10.0L))+((params->a)*(((params->T)+(-(curr.t)))*((params->T)+(-(curr.t)))))))+(((((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t)))))*(((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t))))))+(-(((10.0L)*(params->vh))*((10.0L)*(params->vh)))))))-(((((2.0L)*(params->B))*((curr.yg)+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)), ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->B)*(((((2.0L)*(curr.v))*((params->T)+(-(curr.t))))*(10.0L))+((params->a)*(((params->T)+(-(curr.t)))*((params->T)+(-(curr.t)))))))+(((((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t)))))*(((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t))))))+(-(((10.0L)*(params->vh))*((10.0L)*(params->vh)))))))-(((((2.0L)*(params->B))*((fabsl(curr.xg))+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L))))), fminl(fmaxl((0.0L)-(params->a), (params->vl)-(curr.v)), fminl(fmaxl(params->a, ((10.0L)*(params->vl))-(((10.0L)*(curr.v))+((params->a)*((params->T)+(-(curr.t)))))), fminl(((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->A)*(((((2.0L)*(curr.v))*((params->T)+(-(curr.t))))*(10.0L))+((params->a)*(((params->T)+(-(curr.t)))*((params->T)+(-(curr.t)))))))+((((params->vl)*(10.0L))*((params->vl)*(10.0L)))+(-((((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t)))))*(((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t))))))))))-(((((2.0L)*(params->A))*((curr.yg)+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)), ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->A)*(((((2.0L)*(curr.v))*((params->T)+(-(curr.t))))*(10.0L))+((params->a)*(((params->T)+(-(curr.t)))*((params->T)+(-(curr.t)))))))+((((params->vl)*(10.0L))*((params->vl)*(10.0L)))+(-((((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t)))))*(((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t))))))))))-(((((2.0L)*(params->A))*((fabsl(curr.xg))+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)))))))))) }; return result;
} else {
verdict result = { .id=-1, .val=-1.0L }; return result;
}
} else {
verdict result = { .id=-2, .val=-1.0L }; return result;
}
} else {
verdict result = { .id=-3, .val=-1.0L }; return result;
}
} else {
verdict result = { .id=-4, .val=-1.0L }; return result;
}
} else {
verdict result = { .id=-5, .val=((-1.0L))+(-(fminl(fmaxl((0.0L)-(curr.t), fmaxl(((((params->k)*((params->eps)*(params->eps)))+(-((200.0L)*(params->eps))))*(100.0L))-(((params->k)*(((curr.xg)*(curr.xg))+((curr.yg)*(curr.yg))))+(-((((2.0L)*(curr.xg))*(100.0L))*(10.0L)))), (((params->k)*(((curr.xg)*(curr.xg))+((curr.yg)*(curr.yg))))+(-((((2.0L)*(curr.xg))*(100.0L))*(10.0L))))-((((params->k)*((params->eps)*(params->eps)))+((200.0L)*(params->eps)))*(100.0L)))), fmaxl((0.0L)-(((10.0L)*(curr.v))+((params->a)*((params->T)+(-(curr.t))))), fmaxl(fminl(fmaxl((0.0L)-(params->a), (((10.0L)*(curr.v))+((params->a)*((params->T)+(-(curr.t)))))-((10.0L)*(params->vh))), fminl(fmaxl(params->a, (curr.v)-(params->vh)), fminl(((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->B)*(((((2.0L)*(curr.v))*((params->T)+(-(curr.t))))*(10.0L))+((params->a)*(((params->T)+(-(curr.t)))*((params->T)+(-(curr.t)))))))+(((((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t)))))*(((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t))))))+(-(((10.0L)*(params->vh))*((10.0L)*(params->vh)))))))-(((((2.0L)*(params->B))*((curr.yg)+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)), ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->B)*(((((2.0L)*(curr.v))*((params->T)+(-(curr.t))))*(10.0L))+((params->a)*(((params->T)+(-(curr.t)))*((params->T)+(-(curr.t)))))))+(((((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t)))))*(((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t))))))+(-(((10.0L)*(params->vh))*((10.0L)*(params->vh)))))))-(((((2.0L)*(params->B))*((fabsl(curr.xg))+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L))))), fminl(fmaxl((0.0L)-(params->a), (params->vl)-(curr.v)), fminl(fmaxl(params->a, ((10.0L)*(params->vl))-(((10.0L)*(curr.v))+((params->a)*((params->T)+(-(curr.t)))))), fminl(((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->A)*(((((2.0L)*(curr.v))*((params->T)+(-(curr.t))))*(10.0L))+((params->a)*(((params->T)+(-(curr.t)))*((params->T)+(-(curr.t)))))))+((((params->vl)*(10.0L))*((params->vl)*(10.0L)))+(-((((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t)))))*(((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t))))))))))-(((((2.0L)*(params->A))*((curr.yg)+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)), ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->A)*(((((2.0L)*(curr.v))*((params->T)+(-(curr.t))))*(10.0L))+((params->a)*(((params->T)+(-(curr.t)))*((params->T)+(-(curr.t)))))))+((((params->vl)*(10.0L))*((params->vl)*(10.0L)))+(-((((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t)))))*(((curr.v)*(10.0L))+((params->a)*((params->T)+(-(curr.t))))))))))-(((((2.0L)*(params->A))*((fabsl(curr.xg))+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)))))))))) }; return result;
}
} else {
verdict result = { .id=-6, .val=((-1.0L))+(-((curr.t)-(params->T))) }; return result;
}
} else {
verdict result = { .id=-7, .val=((-1.0L))+(-((0.0L)-(curr.v))) }; return result;
}
} else {
verdict result = { .id=-8, .val=((-1.0L))+(-(fminl(fmaxl((0.0L)-(pre.t), fmaxl(((((params->k)*((params->eps)*(params->eps)))+(-((200.0L)*(params->eps))))*(100.0L))-(((params->k)*(((pre.xg)*(pre.xg))+((pre.yg)*(pre.yg))))+(-((((2.0L)*(pre.xg))*(100.0L))*(10.0L)))), (((params->k)*(((pre.xg)*(pre.xg))+((pre.yg)*(pre.yg))))+(-((((2.0L)*(pre.xg))*(100.0L))*(10.0L))))-((((params->k)*((params->eps)*(params->eps)))+((200.0L)*(params->eps)))*(100.0L)))), fmaxl((0.0L)-(((10.0L)*(pre.v))+((params->a)*((params->T)+(-(pre.t))))), fmaxl(fminl(fmaxl((0.0L)-(params->a), (((10.0L)*(pre.v))+((params->a)*((params->T)+(-(pre.t)))))-((10.0L)*(params->vh))), fminl(fmaxl(params->a, (pre.v)-(params->vh)), fminl(((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->B)*(((((2.0L)*(pre.v))*((params->T)+(-(pre.t))))*(10.0L))+((params->a)*(((params->T)+(-(pre.t)))*((params->T)+(-(pre.t)))))))+(((((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t)))))*(((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t))))))+(-(((10.0L)*(params->vh))*((10.0L)*(params->vh)))))))-(((((2.0L)*(params->B))*((pre.yg)+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)), ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->B)*(((((2.0L)*(pre.v))*((params->T)+(-(pre.t))))*(10.0L))+((params->a)*(((params->T)+(-(pre.t)))*((params->T)+(-(pre.t)))))))+(((((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t)))))*(((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t))))))+(-(((10.0L)*(params->vh))*((10.0L)*(params->vh)))))))-(((((2.0L)*(params->B))*((fabsl(pre.xg))+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L))))), fminl(fmaxl((0.0L)-(params->a), (params->vl)-(pre.v)), fminl(fmaxl(params->a, ((10.0L)*(params->vl))-(((10.0L)*(pre.v))+((params->a)*((params->T)+(-(pre.t)))))), fminl(((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->A)*(((((2.0L)*(pre.v))*((params->T)+(-(pre.t))))*(10.0L))+((params->a)*(((params->T)+(-(pre.t)))*((params->T)+(-(pre.t)))))))+((((params->vl)*(10.0L))*((params->vl)*(10.0L)))+(-((((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t)))))*(((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t))))))))))-(((((2.0L)*(params->A))*((pre.yg)+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)), ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(params->k)))*(100.0L)))+(((params->eps)*(params->eps))*((params->k)*(params->k))))*(((params->A)*(((((2.0L)*(pre.v))*((params->T)+(-(pre.t))))*(10.0L))+((params->a)*(((params->T)+(-(pre.t)))*((params->T)+(-(pre.t)))))))+((((params->vl)*(10.0L))*((params->vl)*(10.0L)))+(-((((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t)))))*(((pre.v)*(10.0L))+((params->a)*((params->T)+(-(pre.t))))))))))-(((((2.0L)*(params->A))*((fabsl(pre.xg))+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)))))))))) }; return result;
}
} else {
verdict result = { .id=-9, .val=((-1.0L))+(-((pre.t)-(params->T))) }; return result;
}
} else {
verdict result = { .id=-10, .val=((-1.0L))+(-((0.0L)-(pre.v))) }; return result;
};
}

verdict boundaryDist(state pre, state curr, const parameters* const params) {
  if (((curr.xg >= 0.0L) && (curr.k >= 0.0L)) || ((curr.xg <= 0.0L) && (curr.k <= 0.0L))) {
if (curr.yg > 0.0L) {
if ((fabsl(curr.k))*(params->eps) <= 100.0L) {
if ((((curr.k)*((params->eps)*(params->eps)))-((200.0L)*(params->eps)))*(100.0L) < ((curr.k)*(((curr.xg)*(curr.xg))+((curr.yg)*(curr.yg))))-((((2.0L)*(curr.xg))*(100.0L))*(10.0L))) {
if (((curr.k)*(((curr.xg)*(curr.xg))+((curr.yg)*(curr.yg))))-((((2.0L)*(curr.xg))*(100.0L))*(10.0L)) < (((curr.k)*((params->eps)*(params->eps)))+((200.0L)*(params->eps)))*(100.0L)) {
if (0.0L <= curr.vl) {
if (curr.vl < curr.vh) {
if ((params->A)*(params->T) <= (10.0L)*((curr.vh)-(curr.vl))) {
if ((params->B)*(params->T) <= (10.0L)*((curr.vh)-(curr.vl))) {
if (-(params->B) <= curr.a) {
if (curr.a <= params->A) {
if (((10.0L)*(pre.v))+((curr.a)*(params->T)) >= 0.0L) {
if (((pre.v <= curr.vh) && (((10.0L)*(pre.v))+((curr.a)*(params->T)) <= (10.0L)*(curr.vh))) || (((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(curr.k)))*(100.0L)))+(((params->eps)*(params->eps))*((curr.k)*(curr.k))))*(((params->B)*(((((2.0L)*(pre.v))*(params->T))*(10.0L))+((curr.a)*((params->T)*(params->T)))))+(((((pre.v)*(10.0L))+((curr.a)*(params->T)))*(((pre.v)*(10.0L))+((curr.a)*(params->T))))-(((10.0L)*(curr.vh))*((10.0L)*(curr.vh))))) <= ((((2.0L)*(params->B))*((fabsl(curr.xg))-((10.0L)*(params->eps))))*(10000.0L))*(100.0L)) || ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(curr.k)))*(100.0L)))+(((params->eps)*(params->eps))*((curr.k)*(curr.k))))*(((params->B)*(((((2.0L)*(pre.v))*(params->T))*(10.0L))+((curr.a)*((params->T)*(params->T)))))+(((((pre.v)*(10.0L))+((curr.a)*(params->T)))*(((pre.v)*(10.0L))+((curr.a)*(params->T))))-(((10.0L)*(curr.vh))*((10.0L)*(curr.vh))))) <= ((((2.0L)*(params->B))*((curr.yg)-((10.0L)*(params->eps))))*(10000.0L))*(100.0L)))) {
if (((curr.vl <= pre.v) && (((10.0L)*(pre.v))+((curr.a)*(params->T)) >= (10.0L)*(curr.vl))) || (((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(curr.k)))*(100.0L)))+(((params->eps)*(params->eps))*((curr.k)*(curr.k))))*(((params->A)*(((((2.0L)*(pre.v))*(params->T))*(10.0L))+((curr.a)*((params->T)*(params->T)))))+((((curr.vl)*(10.0L))*((curr.vl)*(10.0L)))-((((pre.v)*(10.0L))+((curr.a)*(params->T)))*(((pre.v)*(10.0L))+((curr.a)*(params->T)))))) <= ((((2.0L)*(params->A))*((fabsl(curr.xg))-((10.0L)*(params->eps))))*(10000.0L))*(100.0L)) || ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(curr.k)))*(100.0L)))+(((params->eps)*(params->eps))*((curr.k)*(curr.k))))*(((params->A)*(((((2.0L)*(pre.v))*(params->T))*(10.0L))+((curr.a)*((params->T)*(params->T)))))+((((curr.vl)*(10.0L))*((curr.vl)*(10.0L)))-((((pre.v)*(10.0L))+((curr.a)*(params->T)))*(((pre.v)*(10.0L))+((curr.a)*(params->T)))))) <= ((((2.0L)*(params->A))*((curr.yg)-((10.0L)*(params->eps))))*(10000.0L))*(100.0L)))) {
if (curr.v == pre.v) {
if (curr.t == 0.0L) {
verdict result = { .id=1, .val=((((((((((((((0.0L)+(-(fminl(fmaxl((0.0L)-(curr.xg), (0.0L)-(curr.k)), fmaxl(curr.xg, curr.k)))))+(-((0.0L)-(curr.yg))))+(-(((fabsl(curr.k))*(params->eps))-(100.0L))))+(-(((((curr.k)*((params->eps)*(params->eps)))+(-((200.0L)*(params->eps))))*(100.0L))-(((curr.k)*(((curr.xg)*(curr.xg))+((curr.yg)*(curr.yg))))+(-((((2.0L)*(curr.xg))*(100.0L))*(10.0L)))))))+(-((((curr.k)*(((curr.xg)*(curr.xg))+((curr.yg)*(curr.yg))))+(-((((2.0L)*(curr.xg))*(100.0L))*(10.0L))))-((((curr.k)*((params->eps)*(params->eps)))+((200.0L)*(params->eps)))*(100.0L)))))+(-((0.0L)-(curr.vl))))+(-((curr.vl)-(curr.vh))))+(-(((params->A)*(params->T))-((10.0L)*((curr.vh)+(-(curr.vl)))))))+(-(((params->B)*(params->T))-((10.0L)*((curr.vh)+(-(curr.vl)))))))+(-((-(params->B))-(curr.a))))+(-((curr.a)-(params->A))))+(-((0.0L)-(((10.0L)*(pre.v))+((curr.a)*(params->T))))))+(-(fminl(fmaxl((pre.v)-(curr.vh), (((10.0L)*(pre.v))+((curr.a)*(params->T)))-((10.0L)*(curr.vh))), fminl(((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(curr.k)))*(100.0L)))+(((params->eps)*(params->eps))*((curr.k)*(curr.k))))*(((params->B)*(((((2.0L)*(pre.v))*(params->T))*(10.0L))+((curr.a)*((params->T)*(params->T)))))+(((((pre.v)*(10.0L))+((curr.a)*(params->T)))*(((pre.v)*(10.0L))+((curr.a)*(params->T))))+(-(((10.0L)*(curr.vh))*((10.0L)*(curr.vh)))))))-(((((2.0L)*(params->B))*((fabsl(curr.xg))+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)), ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(curr.k)))*(100.0L)))+(((params->eps)*(params->eps))*((curr.k)*(curr.k))))*(((params->B)*(((((2.0L)*(pre.v))*(params->T))*(10.0L))+((curr.a)*((params->T)*(params->T)))))+(((((pre.v)*(10.0L))+((curr.a)*(params->T)))*(((pre.v)*(10.0L))+((curr.a)*(params->T))))+(-(((10.0L)*(curr.vh))*((10.0L)*(curr.vh)))))))-(((((2.0L)*(params->B))*((curr.yg)+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)))))))+(-(fminl(fmaxl((curr.vl)-(pre.v), ((10.0L)*(curr.vl))-(((10.0L)*(pre.v))+((curr.a)*(params->T)))), fminl(((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(curr.k)))*(100.0L)))+(((params->eps)*(params->eps))*((curr.k)*(curr.k))))*(((params->A)*(((((2.0L)*(pre.v))*(params->T))*(10.0L))+((curr.a)*((params->T)*(params->T)))))+((((curr.vl)*(10.0L))*((curr.vl)*(10.0L)))+(-((((pre.v)*(10.0L))+((curr.a)*(params->T)))*(((pre.v)*(10.0L))+((curr.a)*(params->T))))))))-(((((2.0L)*(params->A))*((fabsl(curr.xg))+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)), ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(curr.k)))*(100.0L)))+(((params->eps)*(params->eps))*((curr.k)*(curr.k))))*(((params->A)*(((((2.0L)*(pre.v))*(params->T))*(10.0L))+((curr.a)*((params->T)*(params->T)))))+((((curr.vl)*(10.0L))*((curr.vl)*(10.0L)))+(-((((pre.v)*(10.0L))+((curr.a)*(params->T)))*(((pre.v)*(10.0L))+((curr.a)*(params->T))))))))-(((((2.0L)*(params->A))*((curr.yg)+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)))))) }; return result;
} else {
verdict result = { .id=-1, .val=-fabs(curr.t) }; return result;
}
} else {
verdict result = { .id=-2, .val=-fabs(curr.v-pre.v) }; return result;
}
} else {
verdict result = { .id=-3, .val=((-1.0L))+(-(fminl(fmaxl((curr.vl)-(pre.v), ((10.0L)*(curr.vl))-(((10.0L)*(pre.v))+((curr.a)*(params->T)))), fminl(((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(curr.k)))*(100.0L)))+(((params->eps)*(params->eps))*((curr.k)*(curr.k))))*(((params->A)*(((((2.0L)*(pre.v))*(params->T))*(10.0L))+((curr.a)*((params->T)*(params->T)))))+((((curr.vl)*(10.0L))*((curr.vl)*(10.0L)))+(-((((pre.v)*(10.0L))+((curr.a)*(params->T)))*(((pre.v)*(10.0L))+((curr.a)*(params->T))))))))-(((((2.0L)*(params->A))*((fabsl(curr.xg))+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)), ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(curr.k)))*(100.0L)))+(((params->eps)*(params->eps))*((curr.k)*(curr.k))))*(((params->A)*(((((2.0L)*(pre.v))*(params->T))*(10.0L))+((curr.a)*((params->T)*(params->T)))))+((((curr.vl)*(10.0L))*((curr.vl)*(10.0L)))+(-((((pre.v)*(10.0L))+((curr.a)*(params->T)))*(((pre.v)*(10.0L))+((curr.a)*(params->T))))))))-(((((2.0L)*(params->A))*((curr.yg)+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)))))) }; return result;
}
} else {
verdict result = { .id=-4, .val=((-1.0L))+(-(fminl(fmaxl((pre.v)-(curr.vh), (((10.0L)*(pre.v))+((curr.a)*(params->T)))-((10.0L)*(curr.vh))), fminl(((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(curr.k)))*(100.0L)))+(((params->eps)*(params->eps))*((curr.k)*(curr.k))))*(((params->B)*(((((2.0L)*(pre.v))*(params->T))*(10.0L))+((curr.a)*((params->T)*(params->T)))))+(((((pre.v)*(10.0L))+((curr.a)*(params->T)))*(((pre.v)*(10.0L))+((curr.a)*(params->T))))+(-(((10.0L)*(curr.vh))*((10.0L)*(curr.vh)))))))-(((((2.0L)*(params->B))*((fabsl(curr.xg))+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)), ((((10000.0L)+((((2.0L)*(params->eps))*(fabsl(curr.k)))*(100.0L)))+(((params->eps)*(params->eps))*((curr.k)*(curr.k))))*(((params->B)*(((((2.0L)*(pre.v))*(params->T))*(10.0L))+((curr.a)*((params->T)*(params->T)))))+(((((pre.v)*(10.0L))+((curr.a)*(params->T)))*(((pre.v)*(10.0L))+((curr.a)*(params->T))))+(-(((10.0L)*(curr.vh))*((10.0L)*(curr.vh)))))))-(((((2.0L)*(params->B))*((curr.yg)+(-((10.0L)*(params->eps)))))*(10000.0L))*(100.0L)))))) }; return result;
}
} else {
verdict result = { .id=-5, .val=((-1.0L))+(-((0.0L)-(((10.0L)*(pre.v))+((curr.a)*(params->T))))) }; return result;
}
} else {
verdict result = { .id=-6, .val=((-1.0L))+(-((curr.a)-(params->A))) }; return result;
}
} else {
verdict result = { .id=-7, .val=((-1.0L))+(-((-(params->B))-(curr.a))) }; return result;
}
} else {
verdict result = { .id=-8, .val=((-1.0L))+(-(((params->B)*(params->T))-((10.0L)*((curr.vh)+(-(curr.vl)))))) }; return result;
}
} else {
verdict result = { .id=-9, .val=((-1.0L))+(-(((params->A)*(params->T))-((10.0L)*((curr.vh)+(-(curr.vl)))))) }; return result;
}
} else {
verdict result = { .id=-10, .val=((-1.0L))+(-((curr.vl)-(curr.vh))) }; return result;
}
} else {
verdict result = { .id=-11, .val=((-1.0L))+(-((0.0L)-(curr.vl))) }; return result;
}
} else {
verdict result = { .id=-12, .val=((-1.0L))+(-((((curr.k)*(((curr.xg)*(curr.xg))+((curr.yg)*(curr.yg))))+(-((((2.0L)*(curr.xg))*(100.0L))*(10.0L))))-((((curr.k)*((params->eps)*(params->eps)))+((200.0L)*(params->eps)))*(100.0L)))) }; return result;
}
} else {
verdict result = { .id=-13, .val=((-1.0L))+(-(((((curr.k)*((params->eps)*(params->eps)))+(-((200.0L)*(params->eps))))*(100.0L))-(((curr.k)*(((curr.xg)*(curr.xg))+((curr.yg)*(curr.yg))))+(-((((2.0L)*(curr.xg))*(100.0L))*(10.0L)))))) }; return result;
}
} else {
verdict result = { .id=-14, .val=((-1.0L))+(-(((fabsl(curr.k))*(params->eps))-(100.0L))) }; return result;
}
} else {
verdict result = { .id=-15, .val=((-1.0L))+(-((0.0L)-(curr.yg))) }; return result;
}
} else {
verdict result = { .id=-16, .val=((-1.0L))+(-(fminl(fmaxl((0.0L)-(curr.xg), (0.0L)-(curr.k)), fmaxl(curr.xg, curr.k)))) }; return result;
};
}