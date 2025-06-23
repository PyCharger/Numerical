#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define M_PI 3.14159265358979323846 


double f_1(double x, int N)
{
    return M_PI*x*cos(M_PI*x/2.0) - 1.0/4.0*(-16*x+M_PI*M_PI*(-2+pow(x,2)))*sin(M_PI*x/2.0);
}


double u_1(double x, int N)
{
    return sin(M_PI*x/2);
}



double residual(int N, double* c, double u(double,int))
{
    double h = 1.0/N;
    
    double res = h*(u(0, N)-c[0])*(u(0, N)-c[0])/2.0;
    for(int i = 1; i<N; i++)
    {
        //printf("%d\n", i);
	res += h*(u(i*h, N)-c[i])*(u(i*h, N)-c[i]);
    }

    return sqrt(res);
}



void init_Ritz(int N, double* a, double* b, double* c, double* d, double f(double,int), double u(double,int))
{

    double h = 1.0/N;
    
    //vector b

    d[0] = f(0, N)*h/3 + f(h, N)*h/6;
    d[N] = f(1, N)*h/3 + f((N-1)*h, N)*h/6;
    for(int i=1; i<N; i++)
    {
        d[i] = f(i*h, N)*2*h/3 + f((i-1)*h, N)*h/6 + f((i+1)*h, N)*h/6;    
    }
    //главная диагональ
    c[0] = 2.0/h - h/3.0 + h*h/3.0;
    for (int i=1; i<N; i++)
    {
     	c[i] = 4.0/h - (2.0/3.0+2*i*i)*h + 8.0/3.0*i*h*h;
    }
    c[N] = (3+3*h+3*h*h-h*h*h)/(3*h);
    //боковые диагонали
    for(int i = 0; i<N-1; i++)
    {
	a[i]=-(-6+pow(h,2)-pow(h,3)-3*pow(h,2)*(i+1)+2*pow(h,3)*(i+1)+3*pow(h,2)*(i+1)*(i+1))/(3*h);
	//b[i]=a[i];
	b[i] =-(-6+pow(h,2)+pow(h,3)+3*pow(h,2)*i+2*pow(h,3)*i+3*pow(h,2)*i*i)/(3*h);
    }
}



void right_progonka(int N, double* a, double* b, double* c, double* f, double* y)
{
    double* alpha;
    double* beta;

    alpha = new double[N+1];
    beta = new double[N+1];

    alpha[1]=b[0]/c[0];
    beta[1]=f[0]/c[0];
    
    //прямой ход 
    for (int i=1; i<=N-1; i++)
    { 
         alpha[i+1]=b[i]/(c[i]-a[i-1]*alpha[i]);
         beta[i+1]=(f[i]+a[i-1]*beta[i])/(c[i]-a[i-1]*alpha[i]);
    }

    //обратный ход
    y[N] = (f[N]+a[N-1]*beta[N])/(c[N]-a[N-1]*alpha[N]);
    //printf("%0.3e  %0.3e\n", (f[N]+a[N]*beta[N]), (c[N]-a[N]*alpha[N]));
    for (int i=N-1; i>=0; i-=1)
    {
        y[i]=alpha[i+1]*y[i+1]+beta[i+1];
    }
    //for (int i=0; i<=N; i++) printf("%0.3e  %0.3e  %0.3e\n", y[i], alpha[i], beta[i]);
    delete[] alpha;
    delete[] beta;    
}


void test(int N, double eps, double f(double,int), double u(double,int))
{
    //tridiagonal matrix
    double* a;
    double* b;
    double* c;

    a = new double[N];
    b = new double[N];
    c = new double[N+1];
    
    //test matrix
    /*
    for (int i=0; i<N; i++)
    {
        a[i] = 1;
        b[i] = 2;
        c[i] = 3;
    }
    */
    //c[N]=4;

    //right vector
    double* d;
    d = new double[N+1];
    
    //test d
    /* 

    for (int i=0; i<N+1; i++)
    {
        d[i] = 1;
    }
    */
    //matrix and d define
    init_Ritz(N, a, b, c, d, f, u);
    //for (int i=0; i<N; i++) printf("%0.3e  %0.3e  %0.3e\n", c[i], a[i], b[i]);

    double* y;
    y = new double[N];//вектор решений
    
    right_progonka(N, a, b, c, d, y);
    
    for (int i=0; i<N; i++) printf("%0.3e\n", y[i]);
    printf("%0.3e\n", residual(N, y, u));
    delete[] a;
    delete[] b;
    delete[] c;  
    delete[] y;
}


int main()
{ 
    
    int N = 32;
    double eps = 1e-4;
    
    for(int i=0;i<3;i++)
    {
        printf("n = %d ",N);
        test(N,eps,f_1, u_1);
        N *= 2;
    }
    return 0;
}