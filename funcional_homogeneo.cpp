/* PARA MAIS DETALHES DA IMPLEMENTAÇÃO, CONSULTAR O ARQUIVO "funcional_helicoidal.cpp" CUJA IMPLEMENTAÇÃO É ANALOGA Á DESTE ARQUIVO*/

#include "comum.hpp"
#include "entropia.hpp"

//Diferente da fase helicoidal e de skyrmions, para o homogeneo temos apenas 2 parâmetros: theta, aqui como y[0], o valor final de todos os spins
// e a magnetização m, y[1]
double funcional_homogeneo(const real_1d_array &y,double H, double alpha ,double B, double T){
double f=0.0;
f = B*cos(y[0])*cos(y[0]) - H*cos(y[0]) + 2./3.*T*entropia(y[1]);
return f;
}


struct Parametros {
    double H;
    double B;
    double alpha;
    double T;
};

void nlcfunc2_fvec_homogeneo(const real_1d_array &y, real_1d_array &fi,void *ptr){
Parametros *p = static_cast<Parametros*>(ptr);
double H = p-> H;
double B = p-> B;
double T = p -> T;
double alpha = p -> alpha;
fi[0] = funcional_homogeneo(y, H, alpha ,B, T);
}


//Diferente dos outros dois casos, a condição inicial usada será sempre a mesma: theta = 0.1, m = 1.
real_1d_array homogeneo(double h, double alpha ,double b ,double t,const real_1d_array &y0){
try{
Parametros params;
params.H = h;
params.B = b;
params.T = t;
params.alpha = alpha;
real_1d_array s;
real_1d_array bndl;
real_1d_array bndu;
s.setlength(2);
bndl.setlength(2);
bndu.setlength(2);

s[0] = 1;
s[1] = 1;

//theta
bndl[0]=-PI;
bndu[0]=PI;

//m
if(t==0){
bndl[1] = 0;
bndu[1] = 1;
}
else{
bndl[1] = 0;
bndu[1] = 0.998;
}
double epsx = 0.000000001;
double diffstep = 0.000001;
ae_int_t maxits = 5000;
minnlcstate state;
real_1d_array x1;

//Parametros internos do alglib
minnlccreatef(2, y0, diffstep, state);
minnlcsetcond(state, epsx, maxits);
minnlcsetscale(state, s);
minnlcreport rep;
minnlcsetstpmax(state, 10.0);

//Definição do tipo de método utilizado para a minimização (escolhi SQP-BFGS)
minnlcsetalgosqpbfgs(state);

// Condições de contorno. Coloquei uma folga pra ficarem entre [pi-1,pi+1].
minnlcsetbc(state, bndl, bndu);

// Condições para o constrain. nl é o inferior e nu o superior. Nesse caso, temos 0<m<1.


//O resolvedor
alglib::minnlcoptimize(state, nlcfunc2_fvec_homogeneo,NULL ,&params);
minnlcresults(state, x1, rep);
double energia = funcional_homogeneo(x1, h, alpha ,b, t);
printf("Energia = %lf\n m = %lf\n theta = %lf\n",energia,x1[1],x1[0]);

	return x1;
	
}
catch(alglib::ap_error alglib_exception)
{
	real_1d_array x2 = "[1]";
    printf("ALGLIB exception with message homogeneo'%s'\n", alglib_exception.msg.c_str());
    return x2;
}
}