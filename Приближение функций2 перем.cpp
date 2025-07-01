#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double
f1 (double x, double y)
{
  return x*x;
}

double
f2 (double x, double y)
{
  return x*y;
}

double
f3 (double x, double y)
{
  return y*y;
}
void mat_print(int m,int n,double* val, int* ind,int nonzero,double h)
{
    int bol;
    for  (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            bol = 0;
            for(int k=0;k< nonzero;k++)
            {
                if (ind[k]==(i+1) && ind[k+nonzero] == (j+1))
                {
                    printf("%.2f ", val[k]/(h*h)*24);
                    bol = 1;
                    break;
                }
            }
            if (bol == 0)
            {
                printf("0.00 ");
            }
        }
        printf("\n");
    }
}
int sum(int k,int n)
{
    int res = 0, i = n;
    while (k>0)
    {
        res += i;
        k -= 1;
        i -= 1;
    }
    return res;
}
///corretto
void vector_b(int n, double f(double,double),double* b)
{
    int xk,yk;
    int t;
    double h = (double) 1/n;
    for(xk = 0; xk < n+1; xk++)
        {
        for(yk = 0; yk < n+1+n-xk; yk++)
            {
                if (yk<=n)
                    t = (n+1)*yk+xk;
                else 
                    t = (n+1)*(n+1)+sum(yk-(n+1),n)+xk;
                if (xk == 0 && yk == 0)
                    {
                        b[t] = ((h*h)/12.0)*(f(0,h/2.0)+f(h/2.0,0));
                    }
                else if(xk == n && yk == 0)
                    {
                        b[t] = ((h*h)/12.0)*(f(1-h/2.0,0)+2*f(1-h/2.0,h/2.0)+f(1,h/2.0));
                    }
                else if(xk == n && yk == n)
                    {
                        b[t] = ((h*h)/12.0)*(f(1,1-h/2.0)+2*f(1-h/2.0,1)+f(1-h/2.0,1+h/2.0));
                    }
                else if(xk == 0 && yk == 2*n)
                    {
                        b[t] = ((h*h)/12.0)*(f(0,2-h/2.0)+f(h/2.0,2-h/2.0));
                    }
                else if(xk == 0)
                    {
                        b[t] = ((h*h)/12.0)*(f(0,yk*h-h/2.0)+2*f(h/2.0,yk*h-h/2.0)+2*f(h/2.0,yk*h)+f(0,yk*h+h/2.0));
                    }
                else if(yk == 0)
                    {
                        b[t] = ((h*h)/12.0)*(f(xk*h-h/2.0,0)+2*f(xk*h-h/2.0,h/2.0)+2*f(xk*h,h/2.0)+f(xk*h+h/2.0,0));
                    }
                else if(xk == n)
                    {
                        b[t] = ((h*h)/12.0)*(f(1,yk*h-h/2.0)+2*f(1-h/2.0,yk*h)+2*f(1-h/2.0,yk*h+h/2.0)+f(1,yk*h+h/2.0));
                    }
                else if((yk>n) && (xk == 2*n-yk))
                    {
                        b[t] = ((h*h)/12.0)*(f(xk*h+h/2.0,yk*h-h/2.0)+2*f(xk*h,yk*h-h/2.0)+2*f(xk*h-h/2.0,yk*h)+f(xk*h-h/2.0,yk*h+h/2.0));
                    }
                else
                    {
                        b[t] = ((h*h)/6.0)*(f(xk*h+h/2.0,yk*h)+f(xk*h,yk*h+h/2.0)+f(xk*h-h/2.0,yk*h+h/2.0)+f(xk*h-h/2.0,yk*h)+f(xk*h,yk*h-h/2.0)+f(xk*h+h/2.0,yk*h-h/2.0));
                    }
            }
        }
}
void gram_matrix(int n,double* val,int* ind)
{
    int counter = 0;
    int ne = 2*n*(2*n+1)+(n+1)*(3*n+1)+(2+3*(n-1))*n/2+(4+2*(n-1))*n+1;
    double hh = (double) 1.0/n/n/24.0;
    val[counter] = 2*hh;
    ind[counter] = 1;
    ind[counter+ne] = 1;
    counter ++;
    for(int i = 2; i < n+1 ; i++)
    {
        val[counter] = 6*hh;
        ind[counter] = i;
        ind[counter+ne] = i;
        counter++;
    }
    val[counter] = 4*hh;
    ind[counter] = n+1;
    ind[counter+ne] = n+1;
    counter++;
    for(int i = 1; i < n+1 ; i++)
    {
        val[counter] = hh;
        ind[counter] = i;
        ind[counter+ne] = i+1;
        counter++;
        val[counter] = hh;
        ind[counter] = i+1;
        ind[counter+ne] = i; 
        counter++;
    } 
    /////K is ready!Lets go M
    for (int j = 1;j<n;j++)
    {
        val[counter] = 6*hh;
        ind[counter] = j*(n+1)+1;
        ind[counter+ne] = j*(n+1)+1;
        counter ++;
        for(int i = 2; i < n+1 ; i++)
        {
            val[counter] = 12*hh;
            ind[counter] = j*(n+1)+i;
            ind[counter+ne] = j*(n+1)+i;
            counter++;
        }
        val[counter] = 6*hh;
        ind[counter] = j*(n+1)+n+1;
        ind[counter+ne] = j*(n+1)+n+1;
        counter++;
        for(int i = 1; i < n+1 ; i++)
        {
            val[counter] = 2*hh;
            ind[counter] = j*(n+1)+i;
            ind[counter+ne] = j*(n+1)+i+1;
            counter++;
            val[counter] = 2*hh;
            ind[counter] = j*(n+1)+i+1;
            ind[counter+ne] = j*(n+1)+i; 
            counter++;
        } 
    }
    //////M is ready!Lets go N 
    val[counter] = 6*hh;
    ind[counter] = n*(n+1)+1;
    ind[counter+ne] = n*(n+1)+1;
    counter ++;
    for(int i = 2; i < n+1 ; i++)
    {
        val[counter] = 12*hh;
        ind[counter] = n*(n+1)+i;
        ind[counter+ne] = n*(n+1)+i;
        counter++;
    }
    val[counter] = 4*hh;
    ind[counter] = n*(n+1)+n+1;
    ind[counter+ne] = n*(n+1)+n+1;
    counter++;
    for(int i = 1; i < n+1 ; i++)
    {
        val[counter] = 2*hh;
        ind[counter] = n*(n+1)+i;
        ind[counter+ne] = n*(n+1)+i+1;
        counter++;
        val[counter] = 2*hh;
        ind[counter] = n*(n+1)+i+1;
        ind[counter+ne] = n*(n+1)+i; 
        counter++;
    }
    ///////N is ready! Lets go L and L.T
    for (int j = 0;j<n;j++)
    {
        val[counter] = hh;
        ind[counter] = j*(n+1)+1;
        ind[counter+ne] = (j+1)*(n+1)+1;
        counter ++;
        val[counter] = hh;
        ind[counter] = (j+1)*(n+1)+1;
        ind[counter+ne] = j*(n+1)+1;
        counter ++;
        for(int i = 2; i < n+1 ; i++)
        {
            val[counter] = 2*hh;
            ind[counter] = j*(n+1)+i;
            ind[counter+ne] = (j+1)*(n+1)+i;
            counter++;
            val[counter] = 2*hh;
            ind[counter] = (j+1)*(n+1)+i;
            ind[counter+ne] = j*(n+1)+i;
            counter++;
        }
        val[counter] = hh;
        ind[counter] = j*(n+1)+n+1;
        ind[counter+ne] = (j+1)*(n+1)+n+1;
        counter++;
        val[counter] = hh;
        ind[counter] = (j+1)*(n+1)+n+1;
        ind[counter+ne] = j*(n+1)+n+1;
        counter++;
        for(int i = 2; i < n+2 ; i++)
        {
            val[counter] = 2*hh;
            ind[counter] = j*(n+1)+i;
            ind[counter+ne] = (j+1)*(n+1)+i-1;
            counter++;
            val[counter] = 2*hh;
            ind[counter] = (j+1)*(n+1)+i-1;
            ind[counter+ne] = j*(n+1)+i; 
            counter++;
        } 
    }
    int t = 0;
    ////L and L.T are ready!!! Lets go M of unconstant size
    for (int k = n; k>1;k--)
    {
        val[counter] = 6*hh;
        ind[counter] = (n+1)*(n+1)+1+t;
        ind[counter+ne] = (n+1)*(n+1)+1+t;
        counter ++;
        for(int i = 2; i < k ; i++)
        {
            val[counter] = 12*hh;
            ind[counter] = (n+1)*(n+1)+i+t;
            ind[counter+ne] = (n+1)*(n+1)+i+t;
            counter++;
        }
        val[counter] = 6*hh;
        ind[counter] = (n+1)*(n+1)+k+t;
        ind[counter+ne] = (n+1)*(n+1)+k+t;
        counter ++;
        for(int i = 1; i < k ; i++)
        {
            val[counter] = 2*hh;
            ind[counter] = (n+1)*(n+1)+i+t;
            ind[counter+ne] = (n+1)*(n+1)+i+1+t;
            counter++;
            val[counter] = 2*hh;
            ind[counter] = (n+1)*(n+1)+i+1+t;
            ind[counter+ne] = (n+1)*(n+1)+i+t; 
            counter++;
        } 
        t += k;
    }
    val[counter] = 2*hh;
    ind[counter] = (n+1)*(n+1)+1+t;
    ind[counter+ne] = (n+1)*(n+1)+1+t;
    counter ++;
    int tx = 0,ty = 0;
    ///M of unconstant size are ready!!! Lets go L-star and L-star.T 
    for (int k = n; k>0;k--)
    {
        val[counter] = hh;
        ind[counter] = (n+1)*(n+1)+1+tx;
        ind[counter+ne] = (n)*(n+1)+1+ty;
        counter ++;
        val[counter] = hh;
        ind[counter] = (n)*(n+1)+1+ty;
        ind[counter+ne] = (n+1)*(n+1)+1+tx;
        counter ++;
        for(int i = 2; i < k+1 ; i++)
        {
            val[counter] = 2*hh;
            ind[counter] = (n+1)*(n+1)+i+tx;
            ind[counter+ne] = (n)*(n+1)+i+ty;
            counter++;
            val[counter] = 2*hh;
            ind[counter] = (n)*(n+1)+i+ty;
            ind[counter+ne] = (n+1)*(n+1)+i+tx;
            counter++;
        }
        val[counter] = hh;
        ind[counter] = (n+1)*(n+1)+k+tx;
        ind[counter+ne] = (n)*(n+1)+k+1+ty;
        counter++;
        val[counter] = hh;
        ind[counter] = (n)*(n+1)+k+1+ty;
        ind[counter+ne] = (n+1)*(n+1)+k+tx;
        counter++;
        for(int i = 2; i < k+1 ; i++)
        {
            val[counter] = 2*hh;
            ind[counter] = (n+1)*(n+1)+i-1+tx;
            ind[counter+ne] = (n)*(n+1)+i+ty;
            counter++;
            val[counter] = 2*hh;
            ind[counter] = (n)*(n+1)+i+ty;
            ind[counter+ne] = (n+1)*(n+1)+i-1+tx; 
            counter++;
        }
        tx += k;
        ty += k+1;
    }
    ///that's all 
}
double sp (double *a, double *b, int n)
{
  double res = 0;
  for (int i = 0; i < n; i++)
    res += a[i] * b[i];
  return res;
}
void multiply (int n, int ne, double* val, int* ind, double* v, double* Av)
{
  for (int i = 0; i < n; i++)
    {
      Av[i] = 0.0;
    }
    int i, j;
    for(int k = 0; k < ne; k++)
    {
        i = ind[k]-1;
        j = ind[k+ne]-1;
        Av[i] += val[k]*v[j];
    }
}
double norm(double *b,int n)
{
    double s = 0;
    for (int i=0;i<n;++i)
    {
        s += b[i]*b[i];
    }
    return sqrt(s);
}
int DirectLanczos(double* val, int* ind, double* b, double* c,int n,double eps)
{
    double psi1,eta,lambda = 0,beta = 0,*temp,*v0,*v1,*w,alpha,*p,*r,*rnew;
    int iteration = 0;
    int N = (n+1)*(n+1)+n*(n+1)/2;
    int ne = 2*n*(2*n+1)+(n+1)*(3*n+1)+(2+3*(n-1))*n/2+(4+2*(n-1))*n+1;
    p = (double *)malloc(N * sizeof (double));
    temp = (double *)malloc(N * sizeof (double));
    v0 = (double *)malloc(N * sizeof (double));
    v1 = (double *)malloc(N * sizeof (double));
    w = (double *)malloc(N * sizeof (double));
    r = (double *)malloc(N * sizeof (double));
    rnew = (double *)malloc(N * sizeof (double));
    multiply (N, ne, val,ind ,c,temp);
    for (int i = 0; i < N; i++)
    {
        r[i] = b[i] - temp[i];
    }
    psi1 = norm(r,N);
    for (int i = 0; i < N; i++)
    {
        v1[i]= r[i]/psi1;
    }
    for (int j = 1; j < 100; j++)
    {
        multiply(N,ne,val,ind,v1, temp);
        for (int i = 0; i < N; i++)
        {
            w[i] = temp[i] - beta*v0[i];
        } 
        alpha = sp(w,v1,N);
        if (j>1)
        {
            lambda = beta/eta;
            psi1 = psi1*(-lambda);
        }
        eta = alpha - lambda*beta;
        for (int i = 0; i < N; i++)
        {
            p[i] = (v1[i]-beta*p[i])/eta;
        }
        for (int i = 0; i < N; i++)
        {
            c[i] = c[i]+psi1*p[i];
        }
        for (int i = 0; i < N; i++)
        {
            w[i] = w[i]-alpha*v1[i];
        }
        beta = norm(w,N);
        for (int i = 0; i < N; i++)
        {
            v0[i] = v1[i];
            v1[i] = w[i]/beta;
        }
        iteration += 1;
        multiply(N,ne,val,ind, c, temp);
        for (int i = 0; i < N; i++)
        {
            rnew[i] = b[i] - temp[i];
        }
        if (norm(rnew,N)/norm(r,N)<eps)
        {
            break;
        }
    }
    return iteration;
    free(p);
    free(temp);
    free(v0);
    free(v1);
    free(w);
    free(r);
    free(rnew);
}
void test(int n, double eps, double f(double,double),double ff)
{

    double* val;
    int* ind;
    double* b;
    double* c;
    double* temp;
    int iter_num = 0;
    double delta;
    int N = (n+1)*(n+1)+n*(n+1)/2;
    int ne = 2*n*(2*n+1)+(n+1)*(3*n+1)+(2+3*(n-1))*n/2+(4+2*(n-1))*n+1;
    val = (double *)malloc(ne*sizeof(double));
    ind = (int *)malloc(2*ne*sizeof(int));
    b = (double *)malloc(N*sizeof(double));
    c = (double *)malloc(N*sizeof(double));
    temp = (double *)malloc(N*sizeof(double));
    double h = 1/(double)n;
    gram_matrix(n,val,ind);
    vector_b(n,f,b);
    for(int i=0;i<N;i++)
        c[i] = 0.0;
    //mat_print((n+1)*(n+1)+n*(n+1)/2,(n+1)*(n+1)+n*(n+1)/2,val, ind,ne,h);
    printf("%d iterations ",DirectLanczos(val,ind,b,c,n,eps));
    multiply(N,ne,val,ind,c,temp);
    delta = sqrt(fabs(ff-2*sp(b,c,N)+sp(temp,c,N)));
    printf("delta = %0.3e\n",delta);
    free(val);
    free(ind);
    free(b);
    free(c);
    free(temp);
}
int main()
{
    int n = 8;
    double eps = 1e-6;
    //double ff = 7.0/30.0;
    //double ff = 7.0/30.0;
    double ff = 21.0/10.0;
    for(int i=0;i<4;i++)
    {
        printf("n = %d ",n);
        test(n,eps,f3,ff);
        n *= 2;
    }
    return 1;
}
