#include <R.h>

void mylikC(int *nin1, int *nin2, double *h1y1, double *h2y1, double *h1y2, double *h2y2, double *mypar, double *result)
{
	int i, j;
	result[0] =0;
	int n1=nin1[0];
	int n2=nin2[0];
	double R=0;
	for (i=0; i<n1; i++)
	for (j=0; j<n2; j++)
	{
		R=(1-mypar[3]+mypar[3]*exp(mypar[0]+mypar[1]*h1y1[i]+mypar[2]*h2y1[i]))/(1-mypar[3]+mypar[3]*exp(mypar[0]+mypar[1]*h1y2[j]+mypar[2]*h2y2[j]));
		result[0] +=-log(1+R);
	}
}
