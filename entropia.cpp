#include "comum.hpp" // Informações comuns a todos os arquivos

//Variável global contendo a função entropia interpolada
static spline1dinterpolant s;

// < m > (ver dissertação do Matheus para mais informações)
static double mh(double h){
	return 1./tanh(h)-1./h;
}

// Função partição
static double z(double ni){
	return 4.*PI*sinh(ni)/ni;	
}

//Inicialização da interpolação. Será chamada apenas no main.cpp
void init_entropia(){
    try
    {

//Número de pontos usados para interpolar
int m = 10001;
double v[m];
	for(int j=0;j<m;j++){
		v[j] =0.000001+ 0.07*(double)j;
	}

//x são os pontos referentes há magnetização, igualmente espaçados
//y são os valores da entropia propriamente dita
real_1d_array x,y;
x.setlength(m);
y.setlength(m);

//Tipo das condições de contorno. Nesse caso, 2 refere-se há condições normais
ae_int_t natural_bound_type = 2;


for(int j = 0; j<m;j++){
	x[j] = mh(v[j]);
	y[j] = v[j]*mh(v[j])-log(z(v[j]));
}

// Atribuição à variavel s a função entropia
spline1dbuildcubic(x, y, m, natural_bound_type, 0.0, natural_bound_type, 0.0, s);


}
    catch(alglib::ap_error alglib_exception)
    {
        printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
    }

}


//Função que os arquivos que contem os funcionais irão chamar pra calcular a entropia durante a minimização
double entropia(double m){
	return  spline1dcalc(s, m);	
}