#include "comum.hpp" //Informações em comum de todos os arquivos (#define e etc)
#include "entropia.hpp" //Traz a função interpolada da entropia

//y é  o valor de theta da solução
//y[Nm-1] é o comprimento de modulação enquanto y[Nm-2] é o modulo do spin
//Usei a lei dos trapézios para o funcional, com derivada backward
double funcional_helicoidal(const real_1d_array &y,double H, double alpha, double B, double T){
double f=0.0;
	for(int j=1; j<Nm-3;j++){
		f += dx*((y[Nm-2]*y[Nm-2]/(y[Nm-1]*y[Nm-1])*alpha*pow((y[j]-y[j-1])/dx,2)
		-2.*y[Nm-2]*y[Nm-2]*alpha/y[Nm-1]*(y[j]-y[j-1])/dx
		-H*y[Nm-2]*cos(y[j])
		+B*y[Nm-2]*y[Nm-2]*pow(cos(y[j]),2)));
	}

f += dx/2.*(-H*y[Nm-2] + B*y[Nm-2]*y[Nm-2]  //termo de j = 0
    +(y[Nm-2]*y[Nm-2]*alpha*pow((y[n]-y[n-1])/(dx*y[Nm-1]),2)// termo de j = Nm-3
				-2.*y[Nm-2]*y[Nm-2]*alpha/y[Nm-1]*(y[n]	-y[n-1])/dx
				-H*y[Nm-2]*cos(y[n])
				+B*y[Nm-2]*y[Nm-2]*cos(y[n])*cos(y[n])))+2./3.*T*entropia(y[Nm-2]);

return f;
}

//struct para passar os parametros pra função que minimiza o funcional
struct Parametros {
    double H;
    double B;
    double alpha;
    double T;
};

// Função do alglib que minimiza o funcional. fi[0] é a função objetivo, enquanto constrains podem ser adicionados com fi[1], fi[2], etc.
// Se é um constrain do tipo a<f<b ou f = a, é uma informação passada através dos boundaries de constrain com a função minnlcsetnlc2(state, nl, nu);
// void *ptr é a função para passar os parâmetros extras.
void nlcfunc2_fvec_helicoidal(const real_1d_array &y, real_1d_array &fi,void *ptr){
Parametros *p = static_cast<Parametros*>(ptr);
double H = p->H;
double B = p->B;
double alpha = p->alpha;
double T = p->T;

fi[0] = funcional_helicoidal(y,H,alpha,B,T);
}

// "main" do arquivo. está definido como real_1d_array pois será retornado o perfil final para o arquivo main.cpp onde chama este arquivo. y0 é o perfil inicial para começar a minimização
// Irei usar o perfil da minimização anterior como condição inicial para o próximo. 
real_1d_array helicoidal(double h, double alpha, double b, double t,const real_1d_array &y0){

//Inicializa os parâmetros
try{
Parametros params;
params.H = h;
params.B = b;
params.alpha = alpha;
params.T =t;

//Inicializa os arrays de escala (s), e os limites para todas as variaveis (bndl e bndu)
real_1d_array s;
real_1d_array bndl;
real_1d_array bndu;
s.setlength(Nm);
bndl.setlength(Nm);
bndu.setlength(Nm);


for(int j=0;j<Nm-2;j++){
	s[j] = 1;
	bndl[j] = -PI-1;
	bndu[j] = PI+1;
}

//m e lambda
s[Nm-2] = 1;
s[Nm-1] = 1;

//Condição de contorno do helicoidal
bndl[0]=0;
bndu[0]=0;

bndl[n]=PI;
bndu[n]=PI;

//limites da magnetização
if(t==0){
bndl[Nm-2] = 0;
bndu[Nm-2] = 1;
}
else{
bndl[Nm-2] = 0;
bndu[Nm-2] = 0.998;
}
//lambda
bndl[Nm-1] = 0;
bndu[Nm-1] = INFINITY;

//Parâmetros de minimização
double epsx = 0.000000001; //Precisão da minimização
double diffstep = 0.000001; // Passo da diferença finita que o alglib faz
ae_int_t maxits = 5000; //Número máximo de iteraçoes

//x1 será o array com o perfil minimizado final
minnlcstate state;
real_1d_array x1;

// Passando os parametros internos para o alglib
minnlccreatef(Nm, y0, diffstep, state);
minnlcsetcond(state, epsx, maxits);
minnlcsetscale(state, s);
minnlcreport rep;
minnlcsetstpmax(state, 10.0);

//Definição do tipo de método utilizado para a minimização (escolhi SQP-BFGS) - centenas de variáveis - com boa precisão
minnlcsetalgosqpbfgs(state);

// Condições de contorno. Coloquei uma folga pra ficarem entre [pi-1,pi+1].
minnlcsetbc(state, bndl, bndu);

// Condições para o constrain. nl é o inferior e nu o superior. Nesse caso, temos 0<m<1.
/*
real_1d_array nl = "[0,-inf]";
real_1d_array nu = "[0,0]";
minnlcsetnlc2(state, nl, nu);
*/

//O resolvedor
alglib::minnlcoptimize(state, nlcfunc2_fvec_helicoidal,NULL ,&params);
minnlcresults(state, x1, rep);

// energia guarda o valor final do funcional.
double energia = funcional_helicoidal(x1, h, alpha, b, t);
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