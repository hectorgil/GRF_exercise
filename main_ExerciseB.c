//Generate a GRF from a toy model power spectrum or an cosmological-like input power spectrum
//The step to generate a set of particles given a delta(x) field is not implemented. Therefore some lines are commented.
//Hector Gil Marin, all rights reserved.
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>//complex.h always BEFORE fftw3.h
#include <fftw3.h>
#include <math.h>
#include <time.h>

double Norm(double mu,double var);

double P_power_law(double A, int n, double k);

int countlines(char *filename);

double P_interpol(double k0, double *k, double *P, int N);

int main()
{
int mode;
//mode=0;//write grid
mode=1;//write particles

double Pi=(4.*atan(1.));
double Lbox=1000;//size of box
double nbar=0.01;//mean density
double A0=100;//400000;
double P0;
int n0=-3;
int particles=(int)(nbar*pow(Lbox,3));
char output[200],output2[200],output3[200];
char input_pk[200];
double *pcamb,*kcamb;
int Ncamb;
double s8=1;//0.5;
double mu2;
double growth;
sprintf(output,"particles_nbar10-2.txt");
sprintf(output2,"deltax_field.dat");
sprintf(output3,"deltax_threshold.dat");
sprintf(input_pk,"Pkfiducial.txt");

Ncamb=countlines(input_pk);


FILE *f,*g,*h;
int npower=8;
int N=pow(2,npower);
int i,j,k,l,l2,i2,j2,k2;
int Numpart;
double Numpart2;
double *delta_x;
int Ntot=pow(N,3);
double *delta_k_re,*delta_k_im;
double ki,kj,kk,ktot;
double random;
double *kx;
double P,var;
double sx,sy,sz;
double cx,cy,cz;

double delta_mean;
double Normalization;//No volume normalizaation for delta(k) -> delta(x)
double delta_max;
double delta_min;

pcamb=(double *)calloc(Ncamb,sizeof(double));
kcamb=(double *)calloc(Ncamb,sizeof(double));

f=fopen(input_pk,"r");
for(i=0;i<Ncamb;i++)
{
fscanf(f,"%lf %lf\n",&kcamb[i],&pcamb[i]);
}
fclose(f);

P0=P_interpol(0.1,kcamb, pcamb, Ncamb);
printf("nP_powerlaw=%lf\n",nbar*A0);
printf("nP_powerinputfile=%lf\n",nbar*P0);
printf("Expected number of particles %d \n",particles);

delta_k_re=(double *)calloc(Ntot,sizeof(double));
delta_k_im=(double *)calloc(Ntot,sizeof(double));
delta_x=(double *)calloc(Ntot,sizeof(double));


        kx=malloc(sizeof(double)*(N));
        for(i=0;i<N;i++)
        {
                        if(i<N/2+1)
                        {
                                kx[i]=i*1.0*(2.0*Pi/Lbox);
                        }
                        else
                        {
                                kx[i]=-(N-i)*1.0*(2.0*Pi/Lbox);
                        }
        }


srand48(time(NULL));

for(i=0;i<N;i++)
{

ki=kx[i];
if(i!=0 && i!=N/2){i2=N-i;}
else{i2=i;}

for(j=0;j<N;j++)
{

kj=kx[j];
if(j!=0 && j!=N/2){j2=N-j;}
else{j2=j;}

for(k=0;k<N;k++)
{

kk=kx[k];
if(k!=0 && k!=N/2){k2=N-k;}
else{k2=k;}


ktot=sqrt(ki*ki+kj*kj+kk*kk);
mu2=kk*kk/(ktot*ktot);
if(ktot==0){//Integral constrain in real space (no effect in k-space)
delta_k_re[0]=0;
delta_k_im[0]=0;
}
else
{

l=i*N*N+j*N+k;
l2=i2*N*N+j2*N+k2;

growth=0;
//P=(1+growth*mu2)*P_power_law(A0,n0,ktot);
P=(1+growth*mu2)*s8*P_interpol(ktot,kcamb,pcamb,Ncamb);


var=sqrt(P/2.*pow(Lbox,3));//Units of (Mpc/h)^(3). Note the relation between P_continious * L^3 P-discrete, where P-continious has L^3 units, but P-discrete L^6
delta_k_re[l]=Norm(0,1)*var;//Units of (Mpc/h)^3 This is Hk
delta_k_im[l]=Norm(0,1)*var;//Units of (Mpc/h)^3


if(l2!=l)
{
delta_k_re[l2]=delta_k_re[l];
delta_k_im[l2]=-delta_k_im[l];
}
else
{
delta_k_im[l2]=0;
}

}

}}}


//compute P(k) at this point to check. 
double kbin=0.01;
double kmin=0;
double kmax=0.3;
double *Pk,*K,*RePk,*ImPk;
int *modes;
int Nbins=(int)((kmax-kmin)/kbin)+1;
Pk=(double *)calloc(Nbins,sizeof(double));
RePk=(double *)calloc(Nbins,sizeof(double));
ImPk=(double *)calloc(Nbins,sizeof(double));

K=(double *)calloc(Nbins,sizeof(double));
modes=(int *)calloc(Nbins,sizeof(double));

for(l=0;l<Ntot;l++)
{
i=(int)(l*1./(N*N*1.));
j=(int)( (l-N*N*i)*1./N*1.);
k=l-N*N*i-N*j;

ktot=sqrt(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k]);

l2=(int)((ktot-kmin)/kbin);

if(l2<0){l2=0;}

if(l2<Nbins)
{
K[l2]=K[l2]+ktot;
Pk[l2]=Pk[l2]+(delta_k_re[l]*delta_k_re[l]+delta_k_im[l]*delta_k_im[l]);
RePk[l2]=RePk[l2]+delta_k_re[l]*delta_k_re[l];
ImPk[l2]=ImPk[l2]+delta_k_im[l]*delta_k_im[l];
modes[l2]++;
}
}

f=fopen("Pk_in.txt","w");
for(l2=0;l2<Nbins;l2++)
{
K[l2]=K[l2]/modes[l2];
Pk[l2]=Pk[l2]/modes[l2]*pow(Lbox,-3);
RePk[l2]=RePk[l2]/modes[l2]*pow(Lbox,-3);
ImPk[l2]=ImPk[l2]/modes[l2]*pow(Lbox,-3);

if(modes[l2]>0){
fprintf(f,"%d %lf %lf %lf %lf %d\n",l2,K[l2],Pk[l2],RePk[l2],ImPk[l2],modes[l2]);
}
}
fclose(f);


  fftw_complex *in,*out;
  fftw_plan p;

    in =(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Ntot));
    out =(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Ntot));

//assigne real and imag
for(l=0;l<Ntot;l++)
{
i=(int)(l*1./(N*N*1.));
j=(int)( (l-i*N*N)*1./N*1.);
k=l-N*N*i-N*j;

//Grid (inverse) correction for particles (NGC)
if(mode==1)//particles
{
           cx=sin( kx[i]*Pi/(2.*kx[N/2]) )/( kx[i]*Pi/(2.*kx[N/2]) );
           cy=sin( kx[j]*Pi/(2.*kx[N/2]) )/( kx[j]*Pi/(2.*kx[N/2]) );
           cz=sin( kx[k]*Pi/(2.*kx[N/2]) )/( kx[k]*Pi/(2.*kx[N/2]) );
           if(kx[i]==0){cx=1.;}
           if(kx[j]==0){cy=1.;}
           if(kx[k]==0){cz=1.;}
}

if(mode==0)//grid
{
cx=1;
cy=1;
cz=1;
}
      __real__ in[l]=delta_k_re[l]*pow(cx*cy*cz,+1);
      __imag__ in[l]=delta_k_im[l]*pow(cx*cy*cz,+1);
}

        p =  fftw_plan_dft_3d(N,N,N,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);

    fftw_execute(p);//FFT
    fftw_destroy_plan(p);


//Delta_x needs to be normalized! Now has units of (Mpc/h)^3 coming from deltak(k)

if(mode==1){f=fopen(output,"w");}//particles
if(mode==0){
g=fopen(output2,"w");
h=fopen(output3,"w");}

delta_mean=0;
Normalization=pow(N,-3);//1/N^3normalizaation for delta(k) -> delta(x)
delta_max=-9999999999;
delta_min=9999999999;
//int num_part_min,num_part_max;
for(l=0;l<Ntot;l++)
{

delta_x[l]=creal(out[l])*Normalization*pow(Lbox/N,-3.);// Now delta_x unitless Normalization factor for FT

if(mode==0){
fprintf(g,"%e\n",delta_x[l]);
if(delta_x[l]>-1){fprintf(h,"%e\n",delta_x[l]);}
else{fprintf(h,"-1.00000\n");}
}

if(delta_x[l]<delta_min){delta_min=delta_x[l];}
if(delta_x[l]>delta_max){delta_max=delta_x[l];}


i=(int)(l*1./(N*N*1.));
j=(int)( (l-N*N*i)*1./N*1.);
k=l-N*N*i-N*j;


delta_mean=delta_mean+delta_x[l];

if(mode==1){
Numpart2=((delta_x[l]+1)*nbar*pow(Lbox/N*1.,3));//real value
Numpart=(int)(Numpart2);//parte entera
//Numpart is at max -1 from Numpart2
random=drand48();

if(random<Numpart2-Numpart){Numpart=Numpart+1;}//we randomly add this missing particle

if(Numpart>0){

for(l2=0;l2<Numpart;l2++)
{

//random assignment within the grid!
random=drand48()*Lbox/N*1.;
sx=i*Lbox/N*1.+random;

random=drand48()*Lbox/N*1.;
sy=j*Lbox/N*1.+random;

random=drand48()*Lbox/N*1.;
sz=k*Lbox/N*1.+random;
//write particles
fprintf(f,"%lf %lf %lf\n",sx,sy,sz);

}

}


}

}

if(mode==1){fclose(f);}
if(mode==0){fclose(g);fclose(h);}


printf("Sum delta=%e \t",delta_mean);
delta_mean=delta_mean*pow(N,-3);
printf("<delta>=%e, N=%d\n",delta_mean,N);
//Numpart=countlines(output);
//printf("%d particles, density %e (%e)\n",Numpart,Numpart*pow(Lbox,-3),nbar);
printf("Delta max=%lf, Delta min=%lf\n",delta_max,delta_min);


  fftw_complex *in2,*out2;
  fftw_plan p2;

    in2 =(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Ntot));
    out2 =(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Ntot));

//corss check
for(l=0;l<Ntot;l++)
{

      __real__ in2[l]=delta_x[l]*pow(Lbox/N,3);// no normalization and give units of deltak (aka Hk)
      __imag__ in2[l]=0;

}

        p2 =  fftw_plan_dft_3d(N,N,N,in2,out2,FFTW_FORWARD,FFTW_ESTIMATE);

    fftw_execute(p2);//FFT
    fftw_destroy_plan(p2);

free(Pk);
free(K);
free(RePk);
free(ImPk);
free(modes);

Pk=(double *)calloc(Nbins,sizeof(double));
RePk=(double *)calloc(Nbins,sizeof(double));
ImPk=(double *)calloc(Nbins,sizeof(double));
K=(double *)calloc(Nbins,sizeof(double));
modes=(int *)calloc(Nbins,sizeof(double));

for(l=0;l<Ntot;l++)
{
i=(int)(l*1./(N*N*1.));
j=(int)( (l-N*N*i)*1./N*1.);
k=l-N*N*i-N*j;

ktot=sqrt(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k]);

l2=(int)((ktot-kmin)/kbin);
if(l2<0){l2=0;}

if(l2<Nbins)
{
K[l2]=K[l2]+ktot;
Pk[l2]=Pk[l2]+(creal(out2[l])*creal(out2[l])+cimag(out2[l])*cimag(out2[l]));
RePk[l2]=RePk[l2]+creal(out2[l])*creal(out2[l]);
ImPk[l2]=ImPk[l2]+cimag(out2[l])*cimag(out2[l]);
modes[l2]++;
}
}
f=fopen("Pk_out.txt","w");
for(l2=0;l2<Nbins;l2++)
{
K[l2]=K[l2]/modes[l2];
Pk[l2]=Pk[l2]/modes[l2]*pow(Lbox,-3);
RePk[l2]=RePk[l2]/modes[l2]*pow(Lbox,-3);
ImPk[l2]=ImPk[l2]/modes[l2]*pow(Lbox,-3);

if(modes[l2]>0){
fprintf(f,"%d %lf %lf %lf %lf %d\n",l2,K[l2],Pk[l2],RePk[l2],ImPk[l2],modes[l2]);
}
}
fclose(f);

return 0;
}

double P_power_law(double A, int n, double k)
{
double p;
double kp=0.1;
p=A*pow(k/kp,n);//Units of (Mpc/h)^3

return p;
}


double Norm(double mu,double var)
{
double gauss;
double x1,x2;
double Pi=(4.*atan(1.));
x1=drand48();
x2=drand48();

gauss=sqrt(-2.0*log(x1))*cos(2.*Pi*x2);
//gauss=sqrt(-2.0*log(x1))*sin(2.*Pi*x2);

gauss=gauss*sqrt(var)+mu;
return gauss;
}

int countlines(char *filename)
{
FILE* myfile = fopen(filename, "r");
int ch, number_of_lines = -1;

do
{
            ch = fgetc(myfile);
                    if(ch == '\n')
                                        number_of_lines++;
} while (ch != EOF);

if(ch != '\n' && number_of_lines != 0)
    number_of_lines++;

        fclose(myfile);
        return number_of_lines;
}

double P_interpol(double k0, double *k, double *P, int N)
{
double P0,m,n;
int i;
i=-1;
do
{
i=i+1;
}while(k0>k[i] && i<N);
if(i==0)
{
m=( log10(P[i]) - log10(P[i+1]) )/( log10(k[i]) - log10(k[i+1]) );
n=log10(P[i])-m*log10(k[i]);
P0=m*log10(k0)+n;
P0=pow(10.,P0);
}
else
{
m=( log10(P[i]) - log10(P[i-1]) )/( log10(k[i]) - log10(k[i-1]) );
n=log10(P[i])-m*log10(k[i]);
P0=m*log10(k0)+n;
P0=pow(10.,P0);
}
return P0;
}

