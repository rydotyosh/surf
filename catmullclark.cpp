
/*
   --- simple program that prints out the Catmull-Clark eigenstructures that 
		 were precomputed.

       Author : Jos Stam (jstam@aw.sgi.com), Alias|wavefront
*/

/*
 * Modify : Ryogo Yoshimura (ryogo.yoshimura@gmail.com)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <limits>

using namespace std;

#define IX(i,j,n) ((i)+(n)*(j))

typedef vector<double> vd;

struct EigenDat
{
	vd val;
	vd vecI;
	vd Phi[3];
};

bool read_eval (vector<EigenDat> &res)
{
	res.clear();
	FILE * f;
	int Nmax, i, N, K;

#if defined(_WIN32) || defined(__CYGWIN__) || defined(__APPLE__)
	if ( !(f = fopen ( "ccdata50NT.dat", "rb" )) ) return false;
#else
	if ( !(f = fopen ( "ccdata50.dat", "r" )) ) return false;
#endif

	fread ( &Nmax, sizeof(int), 1, f );
	res.resize(Nmax-2);
	for ( i=0 ; i<Nmax-2 ; i++ )
	{
		N = i+3;
		K = 2*N+8;

		EigenDat &e = res[i];
		e.val.resize(K);
		e.vecI.resize(K*K);
		e.Phi[0].resize(K*16);
		e.Phi[1].resize(K*16);
		e.Phi[2].resize(K*16);

		fread ( &e.val[0], sizeof(double), K, f );
		fread ( &e.vecI[0], sizeof(double), K*K, f );
		fread ( &e.Phi[0][0], sizeof(double), K*16, f );
		fread ( &e.Phi[1][0], sizeof(double), K*16, f );
		fread ( &e.Phi[2][0], sizeof(double), K*16, f );
	}
	fclose ( f );
	return true;
}

static void print_eval (vector<EigenDat> &ev)
{
	int i, j, k, l, N, K;

	int sz = ev.size();
	printf ( "Nmax = %d\n\n", sz );

	for ( i=0 ; i<sz ; i++ )
	{
		N = i+3;
		K = 2*N+8;

		printf ( "N=%d\n\n", N );

		printf ( "Eigenvalues:\n\n" );
		for ( j=0 ; j<K ; j++ ) printf ( "%6.3f ", ev[i].val[j] );
		printf ( "\n\n" );

		printf ( "Inverse of Eigenvectors:\n\n" );

		for ( j=0 ; j<K ; j++ )
		{
			for ( k=0 ; k<K ; k++ ) printf ( "%6.3f ", ev[i].vecI[IX(j,k,K)] );
			printf ( "\n" );
		}
		printf ( "\n\n" );

		printf ( "Coefficients of the Eigenbasis functions:\n" );

		for ( k=0 ; k<3 ; k++ )
		{
			printf ( "k=%d:\n", j );
			for ( j=0 ; j<K ; j++ )
			{
				for ( l=0 ; l<16 ; l++ ) printf("%6.3f ",ev[i].Phi[k][IX(j,l,K)]);
				printf ( "\n" );
			}
			printf ( "\n" );
		}
	}

	return;
}

class point
{
public:
	double x,y,z;
	point(){}
	point(double x_,double y_,double z_):x(x_),y(y_),z(z_){}
	point &operator+=(const point &p){x+=p.x;y+=p.y;z+=p.z;return *this;}
};
point operator*(double m, const point &p){return point(p.x*m,p.y*m,p.z*m);}

void ProjectPoints(vector<EigenDat> &es, point *Cp, const point *C, int N)
{
	const int K = 2*N+8;
	const int ni = N-3;
	for(int i=0;i<K;++i)
	{
		Cp[i] = point(0,0,0);
		for(int j=0;j<K;++j)
		{
			Cp[i] += es[ni].vecI[IX(i,j,K)] * C[j];
		}
	}
}

double log2(double x)
{
	if(x==0) return -300;
	return log(x)/log(2);
}

double NFunc(int i, double t)
{
	if(i==0) return (1+t*(-3+t*(3+t*(-1))))/6.0;
	if(i==1) return (4+t*(0+t*(-6+t*3)))/6.0;
	if(i==2) return (1+t*(3+t*(3+t*(-3))))/6.0;
	else return (t*t*t)/6.0;
}

double EvalSpline(const int K, const double *coef, double u, double v)
{
	double r=0;
	for(int i=0;i<16;++i)
	{
		r += coef[K*i] * NFunc(i%4,u) * NFunc(i/4,v);
	}
	return r;
}

point EvalSurf (vector<EigenDat> &es, double u, double v, const point *Cp, int N)
{
	const double n = floor(min(-log2(u), -log2(v)));
	const double p2 = pow(2, n);

	const double orgu = u,orgv = v;

	u *= p2;
	v *= p2;
	int k;
	if(v<0.5)
	{
		k=0;
		u=2*u-1;
		v=2*v;
	}
	else if(u<0.5)
	{
		k=2;
		u=2*u;
		v=2*v-1;
	}
	else
	{
		k=1;
		u=2*u-1;
		v=2*v-1;
	}
	point P(0,0,0);

	const int K = 2*N+8;
	const int ni = N-3;
	for(int i=0;i<K;++i)
	{
		P += pow(es[ni].val[i],n) * EvalSpline(K,&es[ni].Phi[k][0]+i,u,v) * Cp[i];
	}
	return P;
}

main ( int argc, char ** argv )
{
	vector<EigenDat> ev;
	int Nmax;

	if(!read_eval ( ev ))
		exit ( 1 );

	//print_eval ( ev, Nmax );
	int N=3;
	int K=14;
	vector<point> pts(K);
	pts[0]=point(.6,.6,.6);
	pts[1]=point(.7,0,.7);
	pts[2]=point(0,0,1);
	pts[3]=point(0,.7,.7);
	pts[4]=point(0,1,0);
	pts[5]=point(.7,.7,0);
	pts[6]=point(1,0,0);
	pts[7]=point(-.6,.6,-.6);
	pts[8]=point(0,.7,-.7);
	pts[9]=point(.6,.6,-.6);
	pts[10]=point(.7,0,-.7);
	pts[11]=point(-.7,.7,0);
	pts[12]=point(-.6,.6,.6);
	pts[13]=point(-.7,0,.7);

	vector<point> cp(K);
	ProjectPoints(ev, &cp[0], &pts[0], N);

	int div=10;
	for(int i=0;i<=div;++i)
	{
		for(int j=0;j<=div;++j)
		{
			double u=i/(double)div;
			double v=j/(double)div;
			point p=EvalSurf(ev, u, v, &cp[0], N);
			printf("%f %f %f\n", p.x,p.y,p.z);
		}
		printf("\n");
	}
	/*
	   for(int i=0;i<pts.size();++i)
	   {
	   point p = pts[i];
	   printf("%f %f %f\n", p.x,p.y,p.z);
	   }
	   */

	exit ( 0 );
}

