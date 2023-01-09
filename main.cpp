/************************************************************************
Main program to call greens
Version 2.0, May 1, 2010.
Version 3.0, May 17, 2011.
Version 4.0, March 1, 2018.
See greens.cpp for description of changes.
***********************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void input(void);
void analyzenet(void);
void picturenetwork(float *nodvar, float *segvar, const char fname[]);
void greens(void);
void contour(const char fname[]);
void histogram(const char fname[]);
void setuparrays0();
void setuparrays1(int nseg, int nnod);
void setuparrays2(int nnv, int nnt);
void cmgui(float *segvar);
void postgreens(void);

int max = 200, nmaxvessel, nmaxtissue, nmax, rungreens, initgreens, g0method, linmethod, is2d;
int mxx, myy, mzz, nnt, nseg, nnod, nnodfl, nnv, nsp, nnodbc, nodsegm, nsegfl, kmain, imain;
int slsegdiv, nsl1, nsl2;
int nvaryparams, nruns, ntissparams, npostgreensparams, npostgreensout;	//needed for varying parameters, postgreens
int *mainseg, *permsolute, *nodrank, *nodtyp, *nodout, *bcnodname, *bcnod, *bctyp, *lowflow;
int *nodname, *segname, *segtyp, *nspoint, *istart, *nl, *nk, *indx, *ista, *iend;
int *errvesselcount, *errtissuecount;
int *imaxerrvessel, *imaxerrtissue, *nresis, *oxygen, *diffsolute;
int **segnodname, **nodseg, **tisspoints, **nodnod, **ivaryparams;
int ***nbou;

float gtt, fn, c, alphab, p50, cs, cext, hext, req, q0fac, totalq, flowfac = 1.e6 / 60.;
float plow, phigh, clowfac, chighfac, pphighfac;
float pi1 = atan(1.)*4., fac = 1. / 4. / pi1;
float lb, maxl, v, vol, vdom, errfac, tlength, alx, aly, alz, lowflowcrit;
float tlengthq, tlengthqhd, xmax, ymax, scalefac, w2d, r2d;
float *axt, *ayt, *azt, *ds, *diff, *pmin, *pmax, *pmeant, *pmeanv, *psdt, *psdv, *pref, *g0, *g0fac, *g0facnew, *sumal, *dtmin;
float *diam, *rseg, *q, *qdata, *qq, *hd, *oxflux, *segc, *bcprfl, *bchd, *nodvar, *segvar, *qvtemp, *qvfac;
float *x, *y, *lseg, *ss, *cbar, *mtiss, *mptiss, *dqvsumdg0, *dqtsumdg0;
float *epsvessel, *epstissue, *eps, *errvessel, *errtissue, *pinit, *p;
float *rhs,*rhstest,*g0old,*ptt,*ptpt,*qtsum,*qvsum, *xsl0,*xsl1,*xsl2,*clmin,*clint,*cl;
float **start,**scos,**ax,**cnode,**resisdiam,**resis,**bcp, **qv,**qt,**pv,**pev,**pt, **qvseg,**pvseg,**pevseg;
float **paramvalue, *solutefac, *intravascfac, *postgreensparams, *postgreensout;
float **pvt, **pvprev, **qvprev, **cv, **dcdp, **tissparam;
float **ptprev, **ptv, **gamma1, **cv0, **conv0, **gvv,**end,**al,**zv;
float ***rsta,***rend,***dtt,***psl;
double **mat, **rhsg, *rhsl, *matx;

int main(int argc, char *argv[])
{
	int iseg, inod, j, isp;
	char fname[80];
	FILE *ofp;

	input();

	is2d = 0; //set to 1 for 2d version, 0 otherwise
	if (mzz == 1) is2d = 1; //assumes 2d version if all tissue points lie in one z-plane

	setuparrays0();

	setuparrays1(nseg, nnod);

	analyzenet();

	setuparrays2(nnv, nnt);

	for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = segname[iseg];
	for (inod = 1; inod <= nnod; inod++) nodvar[inod] = nodname[inod];

	for (iseg = 1; iseg <= nseg; iseg++) {
		if (segtyp[iseg] == 4 || segtyp[iseg] == 5) segvar[iseg] = log(fabs(qdata[iseg]));
		else segvar[iseg] = 0.;
	}

		//*************************************
		greens();		//run greens
		//*************************************

	return 0;
}