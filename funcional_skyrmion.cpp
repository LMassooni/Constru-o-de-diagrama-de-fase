/* PARA MAIS DETALHES DA IMPLEMENTAÇÃO, CONSULTAR O ARQUIVO "funcional_helicoidal.cpp" CUJA IMPLEMENTAÇÃO É ANALOGA Á DESTE ARQUIVO*/


#include "comum.hpp"
#include "entropia.hpp"

double funcional_skyrmion(const real_1d_array &y,double H, double alpha, double B, double T){
double f= 0.0;
	for(int j=1; j<Nm-2;j++){
	f+= 1./2.*dx*(4.*dx*j*(alpha*y[Nm-2]*y[Nm-2]/(y[Nm-1]*y[Nm-1])*(pow((y[j]-y[j-1])/dx,2)+sin(y[j])*sin(y[j])/(dx*dx*j*j))
	+y[Nm-2]*y[Nm-2]*alpha/(y[Nm-1])*2.*(y[j]-y[j-1])/dx
	+B*y[Nm-2]*y[Nm-2]*cos(y[j])*cos(y[j])
	-H*y[Nm-2]*cos(y[j]))
	+2*y[Nm-2]*y[Nm-2]*alpha/y[Nm-1]*sin(2.*y[j]));
	}

	f+= 1./2.*dx*(2.*dx*n*(alpha*y[Nm-2]*y[Nm-2]/(y[Nm-1]*y[Nm-1])*(pow((0-y[n-1])/dx,2)))
	+y[Nm-2]*y[Nm-2]*alpha/(y[Nm-1])*2.*(0-y[n-1])/dx
	+B*y[Nm-2]*y[Nm-2]
	-H*y[Nm-2])+2./3.*T*entropia(y[Nm-2]);
return f;
}

struct Parametros {
    double H;
    double B;
    double alpha;
    double T;
};


void nlcfunc2_fvec_skyrmion(const real_1d_array &y, real_1d_array &fi,void *ptr){

Parametros *p = static_cast<Parametros*>(ptr);
double H = p->H;
double B = p->B;
double alpha = p->alpha;
double T = p->T;

fi[0] = funcional_skyrmion(y,H,alpha,B,T);

}

real_1d_array skyrmion(double h, double alpha, double b, double t,const real_1d_array &y0){

try{

Parametros params;
params.H = h;
params.B = b;
params.alpha = alpha;
params.T = t;

real_1d_array s;
real_1d_array bndl;
real_1d_array bndu;
s.setlength(Nm);
bndl.setlength(Nm);
bndu.setlength(Nm);

//y0 é a condição inicial, s é a escala das variaveis, bndl e bndu são os limites menores e maiores para as variaveis
for(int j=0;j<Nm-2;j++){
	s[j] = 1;
	bndl[j] = -PI;
	bndu[j] = PI;
}

//escala s
s[Nm-2] = 1;
s[Nm-1] = 1;

//Condição de contorno do skyrmion
bndl[0] = PI;
bndu[0] = PI;

bndl[Nm-3] = 0;
bndu[Nm-3] = 0;

//m

if(t==0){
bndl[Nm-2] = 0;
bndu[Nm-2] = 1;
}
else{
bndl[Nm-2] = 0;
bndu[Nm-2] = 0.998;
}

//lambda
bndl[Nm-1] = 0.001;
bndu[Nm-1] = INFINITY;

double epsx = 0.000000001;
double diffstep = 0.000001;
ae_int_t maxits = 5000;

minnlcstate state;
real_1d_array x1;

//Parametros internos do alglib
minnlccreatef(Nm, y0, diffstep, state);
minnlcsetcond(state, epsx, maxits);
minnlcsetscale(state, s);
minnlcreport rep;
minnlcsetstpmax(state, 10.0);

//Definição do tipo de método utilizado para a minimização (escolhi SQP-BFGS)
minnlcsetalgosqp(state);
// Condições de contorno. Coloquei uma folga pra ficarem entre [pi-1,pi+1].
minnlcsetbc(state, bndl, bndu);

// Condições para o constrain. nl é o inferior e nu o superior. Nesse caso, temos 0<m<1.
/*
real_1d_array nl = "[0,-inf]";
real_1d_array nu = "[0,0]";
minnlcsetnlc2(state, nl, nu);
*/

// Passando as condições pro constrain
//minnlcsetnlc2(state,nl,nu);

//O resolvedor
alglib::minnlcoptimize(state, nlcfunc2_fvec_skyrmion, NULL, &params);
minnlcresults(state, x1, rep);
double energia = funcional_skyrmion(x1,h,alpha,b,t);
printf("Energia = %lf\n m = %lf\n lambda = %lf\n",energia,x1[Nm-2],x1[Nm-1]);
return x1;
}
catch(alglib::ap_error alglib_exception)
{
	real_1d_array x2 = "[1]";
    printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
    return x2;
}
	
}