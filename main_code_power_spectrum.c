#include <stdio.h>
#include <stdlib.h>
#include <complex.h>//complex.h always BEFORE fftw3.h
#include <fftw3.h>
#include <math.h>

long int countlines(char *filename);

double Leg2(double x);

int main()
{
long int N,Ncell,i,j,k,l,l2;
FILE *f;
double *delta;
double *Pk,*K,*Pk2;
long int *modes;
double mu2;
double *ki;
double kmin,kmax,Deltak;
long int Nbins;
char name_file[200];
char name_out[200];
double Pi=(4.*atan(1.));
double Lbox=1000;//size of box
double ktot;

//Read the delta field
sprintf(name_file,"deltaxfield.dat");
//sprintf(name_file,"deltaxfieldaniso.dat");
N=countlines(name_file);
Ncell=(long int)( round(pow(N*1.,1./3.)) );
delta=(double *)calloc(N,sizeof(double));
f=fopen(name_file,"r");
for(i=0;i<N;i++)
{
fscanf(f,"%lf\n",&delta[i]);

}
fclose(f);


//Fourier Transform the delta field
  fftw_complex *in,*out;
  fftw_plan p;

    in =(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(N));
    out =(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(N));


for(i=0;i<N;i++)
{
      __real__ in[i]=delta[i]*pow(Lbox,3)/N*1.;// no normalization and give units of deltak (aka Hk)
      __imag__ in[i]=0;
}

     p =  fftw_plan_dft_3d(Ncell,Ncell,Ncell,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

    fftw_execute(p);
    fftw_destroy_plan(p);

free(delta);

//printf("%lf %lf\n",creal(out[0])/(pow(Lbox,3)/N*1.),cimag(out[0])/(pow(Lbox,3)/N*1.));
//printf("%lf %lf\n",creal(out[256*256])/(pow(Lbox,3)/N*1.),cimag(out[256*256])/(pow(Lbox,3)/N*1.));
//printf("%lf %lf\n",creal(out[256*256*2])/(pow(Lbox,3)/N*1.),cimag(out[256*256*2])/(pow(Lbox,3)/N*1.));
//return 0;

        ki=malloc(sizeof(double)*Ncell);
        for(i=0;i<Ncell;i++)
        {
                        if(i<Ncell/2+1)
                        {
                                ki[i]=i*1.0*(2.0*Pi/Lbox);
                        }
                        else
                        {
                                ki[i]=-(N-i)*1.0*(2.0*Pi/Lbox);
                        }
        }

kmin=2.*Pi/Lbox;
kmax=kmin*Ncell/2.;
Deltak=0.001;
Nbins=(long int)((kmax-kmin)/Deltak)+1;
Pk=(double *)calloc(Nbins,sizeof(double));
Pk2=(double *)calloc(Nbins,sizeof(double));
K=(double *)calloc(Nbins,sizeof(double));
modes=(long int *)calloc(Nbins,sizeof(long int));

for(l=0;l<N;l++)
{
i=(long int)(l*1./(Ncell*Ncell*1.));
j=(long int)( (l-Ncell*Ncell*i)*1./Ncell*1.);
k=l-Ncell*Ncell*i-Ncell*j;

ktot=sqrt(ki[i]*ki[i]+ki[j]*ki[j]+ki[k]*ki[k]);
mu2=ki[k]*ki[k]/(ktot*ktot);
mu2=sqrt(mu2);

l2=(long int)((ktot-kmin)/Deltak);
if(l2<0){l2=0;}

if(l2<Nbins)
{
K[l2]=K[l2]+ktot;
Pk[l2]=Pk[l2]+(creal(out[l])*creal(out[l])+cimag(out[l])*cimag(out[l]));
Pk2[l2]=Pk2[l2]+(creal(out[l])*creal(out[l])+cimag(out[l])*cimag(out[l]))*Leg2(mu2);
modes[l2]++;
}
}
sprintf(name_out,"PowerSpectrum.txt");
//sprintf(name_out,"PowerSpectrum_aniso.txt");
f=fopen(name_out,"w");
for(l2=0;l2<Nbins;l2++)
{
K[l2]=K[l2]/modes[l2];
Pk[l2]=Pk[l2]/modes[l2]*pow(Lbox,-3);
Pk2[l2]=5.*Pk2[l2]/modes[l2]*pow(Lbox,-3);

if(modes[l2]>0){
fprintf(f,"%lf %lf %lf %ld\n",K[l2],Pk[l2],Pk2[l2],modes[l2]);
}
}
fclose(f);
free(K);
free(Pk);
free(Pk2);
free(modes);
free(ki);
return 0;
}


long int countlines(char *filename)
{
FILE* myfile = fopen(filename, "r"); if(myfile==NULL){printf("Cannot open %s file. Exiting now...\n",filename);exit(0);}
long int ch, number_of_lines = -1;

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

double Leg2(double x)
{
	double f=0.5*(3.*x*x-1.);

	return f;
}

