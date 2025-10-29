#include "comum.hpp" //Informações em comum (#define)
#include "funcional_helicoidal.hpp" //minimiza o helicoidal
#include "funcional_skyrmion.hpp"  //minimiza o skyrmion
#include "funcional_homogeneo.hpp" //minimiza o homogeneo
#include "entropia.hpp" //Constrói a interpolação da entropia

int main(){

FILE *diagrama;
FILE *perfil_skyrmion;
FILE *perfil_helicoidal;
diagrama = fopen("sky-hel-0.txt","w+");
perfil_helicoidal = fopen("perfil_skyrmion.txt","w+");
perfil_skyrmion = fopen("perfil_helicoidal.txt","w+");

//Inicialização da entropia
init_entropia();

//Parâmetros do diagrama
double H = 0.0;
double B = 0.0;
double T = 0.1;

//Força relativa das interações de troca e DM
double alpha = 0.2;

//Quantidade de pontos de H calculados
int intervalo = 1; 

// Guarda as energias minimizadas
double ene_hel;
double ene_sky;
double ene_hom;

// Inicialização das condições iniciais
real_1d_array y0hel,y0sky,y0hom,hel,sky,hom,sol_sky,sol_hel;

y0hel.setlength(Nm);
y0sky.setlength(Nm);
sol_hel.setlength(Nm);
sol_sky.setlength(Nm);
y0hom.setlength(2);
hom.setlength(2);
hel.setlength(Nm);
sky.setlength(Nm);


for(int l = 0;l<Nm-2;l++){
	y0hel[l] = PI*l/n; //Condição inicial helicoidal -> reta crescente
	y0sky[l] = PI - PI*l/n; //Condição inicial skyrmion -> reta decrescente
}
if(T==0){
y0hel[Nm-2] = 1; //magnetização m helicoidal
y0hel[Nm-1] = 3.; // lambda helicoidal

y0sky[Nm-2] = 1; //magnetização m
y0sky[Nm-1] = 4.; //lambda skyrmion

y0hom[0] = 0.1; //theta do homogeneo
y0hom[1] = 1; //magnetização m homogeneo
}
else{
y0hel[Nm-2] = 0.998;
y0hel[Nm-1] = 3.;

y0sky[Nm-2] = 0.998; 
y0sky[Nm-1] = 4.;

y0hom[0] = 0.1; 
y0hom[1] = 0.998; 
}

//Minimização dos funcionais
for(int h = 0;h<intervalo;h++){
hom = homogeneo(H, alpha, B, T, y0hom);
hel = helicoidal(H, alpha, B, T, y0hel);
sky = skyrmion(H, alpha, B, T, y0sky);
ene_hel = funcional_helicoidal(hel, H, alpha, B, T);
ene_sky = funcional_skyrmion(sky, H, alpha, B, T);
ene_hom = funcional_homogeneo(hom, H, alpha, B, T);

fprintf(diagrama, "%lf %lf %lf %lf\n", H, ene_hel, ene_sky, ene_hom);
y0hel = hel;
y0sky = sky;

//Usar o if abaixo para salvar o perfil das soluções para algum valor dos parâmetros
if(h==0){
sol_sky = sky;
sol_hel = hel;
}
H += 0.1;
}

for(int v = 0;v<Nm-2 ;v++){
	fprintf(perfil_skyrmion,"%lf %lf\n",v*dx, sol_sky[v]);
	fprintf(perfil_helicoidal,"%lf %lf\n",v*dx, sol_hel[v]);
	
}

fclose(diagrama);
fclose(perfil_helicoidal);
fclose(perfil_skyrmion);
	return 0;
}