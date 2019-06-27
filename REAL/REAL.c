/****************************************************************************
 *REAL - Rapid Earthquake Association and Location
 *
 *What you need:
 *  1. Traveltime table for P or/and S waves (dist,dep,P arrival,S arrival ...) 
 *  2. Station information (stlo,stla,net,sta,chan,elev)
 *  3. Picks at each station and their weighting factors and amplitudes
 *  4. Control parameters (see usage)
 *      a. searching range and grid size 
 *      b. average velocities for both P and S waves
 *      c. date of the day
 *      d. thresholds
 *
 *Output:
 *  1. Associated and located earthquakes with origin time, magnitude, and location
 *  2. Associated picks for each earthquake
 *  (local magnitude is preliminarily estimated based on HUTTON and BOORE, BSSA, 1987)
 *
 *Usage:
 *  See usage as below
 *
 *Author: 
 *  Miao Zhang, Stanford University
 *  Now at Dalhousie University (miao.zhang@dal.ca)
 *
 *Reference:
 *  Miao Zhang, William Ellsworth and Greg Beroza, Rapid Earthquake Association and Location, 2019 (Submitted)
 *
 *Revision history:
 *  June     2018       M. Zhang    Initial version in C 
 *  June     2019		M. Zhang	Release version 1.0
 ************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define PHASESEL "phase_sel.txt"
#define CATALOGSEL "catalog_sel.txt"
#define RESOLUTION "resolution.txt"

#define D2R  .017453292519943295769237
#define R2D    57.2957795130823208768
/*default to process one day continuous data*/
//#define MAXTIME 86400.00 //one day
#define MAXTIME 2700000.00 //one month

typedef struct ttable{
    double gdist;
    double dep;
    double ptime;
    double stime;
    double prayp;
    double srayp;
    double phslow;
    double shslow;
    char pphase[10];
    char sphase[10];
} TTT;


typedef struct reselect{
	int num1;
	char otime1[50];
	double atime1;
	double std1;
	double lat1;
	double lon1;
	double dep1;
	int nofp1;
	int nofs1;
	int ntotal1;
} SELECT;

typedef struct trigg{
	double trig;
	double weight;
	double amp;
} TRIG;

typedef struct stationinfo{
    double stlo;
    double stla;
    char   net[5];
    char   sta[8];
    char   comp[4];
    double elev;
} STATION;

void GCinit(const double*,const double*,const double*,const double*,double*);
double CalculateMedian(double *,int);
double CalculateMean(double *,int);
double CalculateStd(double *,double,int);
double Find_min(double **,int,int);
double Find_max(double **,int,int);
void Find_min_loc(double **,int,int,double *,int *,int *);
int Readttime(char *, TTT *,int);
int Readstation(char *, STATION *,int);
int DetermineNp(double **,int, int);
int DetermineNg(TRIG **,TRIG **,int, int);
void SortTriggers(TRIG **,TRIG **,double **,double **,int,int);
void SortTriggers0(TRIG **,TRIG **,double **,double **,double**,double**,int,int);
void DeleteOne(double **,int,int,int);
int DetermineNprange(double **,double,int,int);
void DetermineNps0range(double **,double **,double,double,double,double,int,int);
int ReselectFinal(SELECT *,int);
void Accounttriggers_homo(double,double,double,double,double,double,int);
void Accounttriggers_layer(double,double,double,double,double,double,int);
void Sortpscounts(double **,int);


//global 
double tint;
double **pscounts;
STATION *ST;
TRIG **TGP,**TGS;
TTT *TB;
double ptw,stw,nrt;
double  **ptrig,**strig,**ptrig0,**strig0;
double vp0,vs0,s_vp0,s_vs0;
int NNps,Nps2;
int np0,ns0,nps0;
int *np0_start,*np0_end,*ns0_start,*ns0_end;
double tpmin0,tdx,tdh,trx,trh;
int ppp;
double dtps;

int Nst = 500; //maximum number of stations
int Nps = 20000; //maximum number of P/S picks recorded at one station
int Ntb = 20000; //maximum number of lines in traveltime table
int ispeed = 1; //default setting ispeed = 1

//main function
int main(int argc, char **argv){
	int i,j,k,l,m,n;
	FILE *fp,*fpr,*fp1,*fp2;
	char output1[256],output2[256],dir[256],input[256];
	int test,error,pcount,scount,nnn,ps,nppp,nselect;
	double dx,dh,rx,rh;
	double tp0_cal,ts0_cal,tp_cal,ts_cal,tp_pre,ts_pre,tp_pre_b,ts_pre_b,tp_pre_e,ts_pre_e;
	double median,std,GCarc,distmax;
	double told,lonref,latref,elevref,latref0,lonref0;
	double lonmin,lonmax,latmin,latmax,lat0,lon0,dep;
	double stlamin,stlamax,stlomin,stlomax;
	int ttd,tth,tts,mmm;
	int nlat,nlon,ndep;
	double ttm,std0,ptemp,rsel;
	char otime[50];
	int igrid,ires,ielev,ig,ih,im,iremove;
	SELECT *RELC;
	double **pamp0,**samp0,*mag;
	double mag_median,mag_std,p_mag,s_mag,pamp,samp;
	int nyear,nmon,nday;
	double tpmin,tpmax,tsmin,tsmax,Maxt0;
	
	//initiating parameters
	error = 0;
	igrid = 0;
	ielev = 0;
	ires = 0;
	rsel = 5;
	latref0 = -10000;
	lonref0 = -10000;
    s_vp0 = 10000;
    s_vs0 = 10000;

	for(i=1;!error && i<argc; i++){
		if(argv[i][0] == '-'){
			switch(argv[i][1]){
				case 'R':
					sscanf(&argv[i][2],"%lf/%lf/%lf/%lf/%lf/%lf/%lf",&rx,&rh,&dx,&dh,&tint,&latref0,&lonref0);
					break;
				case 'S':
					sscanf(&argv[i][2],"%d/%d/%d/%d/%lf/%lf/%lf/%lf/%d",&np0,&ns0,&nps0,&nppp,&std0,&dtps,&nrt,&rsel,&ires);
					break;
				case 'V':
					sscanf(&argv[i][2],"%lf/%lf/%lf/%lf/%d",&vp0,&vs0,&s_vp0,&s_vs0,&ielev);
					break;
				case 'G':
					sscanf(&argv[i][2],"%lf/%lf/%lf/%lf",&trx,&trh,&tdx,&tdh);
					igrid = 1; 
					break;
				case 'D':
					sscanf(&argv[i][2],"%d/%d/%d",&nyear,&nmon,&nday);
					break;
				default:
					error = 1;
					break;
			}
		}
	}
	
	//Usage
    if(argc < 3 || error == 1){
        fprintf(stderr, "Usage:  Rapid Earthquake Association and Location (REAL)\n");
        fprintf(stderr, "   -D(nyear/nmon/nday) -R(rx/rh/tdx/tdh/tint[/latref0/lonref0]]) -V(vp0/vs0/[s_vp0/s_vs0/ielev])\n");
        fprintf(stderr, "   -S(np0/ns0/nps0/nppp/std0/dtps/nrt[/rsel/ires]) [-G(trx/trh/tdx/tdh)] station pickdir [ttime]\n");
        fprintf(stderr, "   ------------------------------------explanation--------------------------------------------\n");
        fprintf(stderr, "   -D: date of the day (year/month/day)\n");
        fprintf(stderr, "   -R: search ranges and grids around the nearest station in horizontal direction and depth,\n");
        fprintf(stderr, "       phase and event interval, reference location (degree/km/degree/km/sec[/degree/degree])\n");
        fprintf(stderr, "   -V: average velocities and near-surface velocities of P and S waves, station elevation_or_not\n");
		fprintf(stderr, "       (km/s|km/s|[km/s|km/s|int])\n");
        fprintf(stderr, "   -S: thresholds: number of picks (P,S,P+S),effective number of grids,residual threshold, \n");
		fprintf(stderr, "       S-P interval,nrt*length of time window,only keep picks with residuals < rsel*residual,\n");
		fprintf(stderr,	"       resolution_or_not (int/int/int/int/double/double/double/[double/int])\n");
        fprintf(stderr, "   -G: ranges and grids for traveltime table in horizontal direction and depth (degree/km/degree/km)\n");
        fprintf(stderr, "   station: station information; pickdir: directory of picks; ttime: [traveltime table]\n");
        exit(-1);
    }    

	fprintf(stderr,"Max Setting: Nst %-5d Nps %-5d Ntb %-5d\n",Nst,Nps,Ntb);
	
	/* read station information */
	if(igrid == 0){strcpy(input,argv[5]);}else{strcpy(input,argv[6]);}
	ST = (STATION *)malloc(sizeof(STATION)*Nst);
	Nst = Readstation(input,ST,Nst);

	if(ielev == 0){
		for(i=0;i<Nst;i++)ST[i].elev = 0.0;
	}
	
	stlamin = 1.0e8; stlomin = 1.0e8;
	stlamax = -1.0e8; stlomax = -1.0e8;
	for(i=0;i<Nst;i++){
		if(ST[i].stla > stlamax) stlamax = ST[i].stla;
		if(ST[i].stlo > stlomax) stlomax = ST[i].stlo;
		if(ST[i].stla < stlamin) stlamin = ST[i].stla;
		if(ST[i].stlo < stlomin) stlomin = ST[i].stlo;
	}
	GCinit(&stlamin,&stlomin,&stlamax,&stlomax,&distmax);


	/* read triggers */
	if(igrid == 0){strcpy(dir,argv[6]);}else{strcpy(dir,argv[7]);}

	TGP = (TRIG **)malloc(sizeof(TRIG*)*Nst);
	TGS = (TRIG **)malloc(sizeof(TRIG*)*Nst);
	for(i=0;i<Nst;i++){
		TGP[i] = (TRIG *)malloc(sizeof(TRIG)*Nps);
		TGS[i] = (TRIG *)malloc(sizeof(TRIG)*Nps);
	}

	for(i=0;i<Nst;i++){
		for(j=0;j<Nps;j++){
			TGP[i][j].trig = 1.0e8;
			TGP[i][j].weight = 0.0;
			TGP[i][j].amp = 0.0;
			TGS[i][j].trig = 1.0e8;
			TGS[i][j].weight = 0.0;
			TGS[i][j].amp = 0.0;
		}
	}

	for(i=0;i<Nst;i++){
		sprintf(input,"%s/%s.%s.P.txt",dir,ST[i].net,ST[i].sta);
        if((fp=fopen(input,"r"))==NULL){
            fprintf(stderr,"Can not open file in ReadFile %s\n", input);
        }else{
			test = 0;
			for (j = 0; j < Nps; j++){
				if(fscanf(fp, "%lf %lf %lf",&TGP[i][j].trig,&TGP[i][j].weight,&TGP[i][j].amp)==EOF)test=1;
				if(test==1)break;
			}
		}
        fclose(fp);

		sprintf(input,"%s/%s.%s.S.txt",dir,ST[i].net,ST[i].sta);
        if((fp=fopen(input,"r"))==NULL){
            fprintf(stderr,"Can not open file in ReadFile %s\n", input);
        }else{
			test = 0;
			for (j = 0; j < Nps; j++){
				if(fscanf(fp, "%lf %lf %lf",&TGS[i][j].trig,&TGS[i][j].weight,&TGS[i][j].amp)==EOF)test=1;
				if(test==1)break;
			}
		}
        fclose(fp);
	}

	/* read travel time table */
	if(igrid == 1){
		strcpy(input,argv[8]);
		if((TB = malloc(sizeof(TTT)*Ntb)) == NULL){
			fprintf(stderr,"malloc memory error for TTT\n");
			exit(-1);
		}
		Ntb=Readttime(input,TB,Ntb);
	}	
	
	Nps = DetermineNg(TGP,TGS,Nst,Nps);
	NNps = Nps;
	
	fprintf(stderr,"Actual     : Nst %-5d Nps %-5d Ntb %-5d\n",Nst,Nps-1,Ntb);
    if(latref0 > -999 && lonref0 > -999){
        fprintf(stderr,"searching range: %lf %lf %lf %lf\n",latref0-rx,latref0+rx,lonref0-rx,lonref0+rx);
    }else{
        fprintf(stderr,"searching range: %lf %lf %lf %lf\n",stlamin,stlamax,stlomin,stlomax);
    }

	
	ptrig = (double **) malloc(sizeof(double*)*Nst);
	strig = (double **) malloc(sizeof(double*)*Nst);
	ptrig0 = (double **) malloc(sizeof(double*)*Nst);
	strig0 = (double **) malloc(sizeof(double*)*Nst);
	pamp0  = (double **) malloc(sizeof(double*)*Nst);
	samp0  = (double **) malloc(sizeof(double*)*Nst);
	for(i=0;i<Nst;i++){
		ptrig[i] = (double *)malloc(sizeof(double)*Nps);
		strig[i] = (double *)malloc(sizeof(double)*Nps);
		ptrig0[i] = (double *)malloc(sizeof(double)*Nps);
		strig0[i] = (double *)malloc(sizeof(double)*Nps);
		pamp0[i] = (double *)malloc(sizeof(double)*Nps);
		samp0[i] = (double *)malloc(sizeof(double)*Nps);
	}
	
	//default number of events (picks*Nst)
    RELC = (SELECT *)malloc(sizeof(SELECT)*Nst*Nps);
	/* determine traveltime within one grid*/
	ptw = sqrt(2*(dx*111.19)*(dx*111.19)+dh*dh)/vp0;
	stw = sqrt(2*(dx*111.19)*(dx*111.19)+dh*dh)/vs0;
	if(tint < stw/2.0)tint=stw/2.0;
	fprintf(stderr,"p-window= %.2f sec; s-window= %.2f sec; phase&event-window= %.2f sec\n",nrt*ptw,nrt*stw,tint);
    SortTriggers(TGP,TGS,ptrig,strig,Nst,Nps);
	//sort triggers and their amplitudes
	SortTriggers0(TGP,TGS,ptrig0,strig0,pamp0,samp0,Nst,Nps);
	
	nlat = (int)(2*rx/dx + 1);
    nlon = (int)(2*rx/dx + 1);
    ndep = (int)(rh/dh + 1);
	nnn = nlat*nlon*ndep;
    printf("Nlat= %d Nlon= %d Ndep= %d\n",nlat,nlon,ndep);

	pscounts = (double **)malloc(nnn*sizeof(double*));
	for(k=0;k<nnn;k++){
		pscounts[k] = (double *)malloc(8*sizeof(double));
	}	

	np0_start = (int *)malloc(sizeof(int)*Nst);
	np0_end = (int *)malloc(sizeof(int)*Nst);
	ns0_start = (int *)malloc(sizeof(int)*Nst);
	ns0_end = (int *)malloc(sizeof(int)*Nst);

	told = 0.0;
	mmm = 0;
	m = 0;
	
	Maxt0 = Find_max(strig,Nst,Nps);
	//search each initiating P picks
	while(Find_min(ptrig,Nst,Nps) < Maxt0){
		//jump:
		Nps = DetermineNp(ptrig,Nst,Nps);
		Find_min_loc(ptrig,Nst,1,&tpmin0,&m,&n);
		if(fabs(tpmin0 - 1.0e8) < 1) break;

		if(latref0 > -999 && lonref0 > -999){
			latref = latref0;
			lonref = lonref0;
		}else{
			lonref = ST[m].stlo;
			latref = ST[m].stla;
		}
		elevref = ST[m].elev;
		
		//Make sure you know what you are doing!
		//if(fabs(tpmin0 - told) < tint/2){
		//	DeleteOne(ptrig,m,Nps,n);
		//	goto jump;
		//}
		
		tpmin = tpmin0 - 1.2*(distmax*111.19/vp0);
		tpmax = tpmin0 + 1.2*(distmax*111.19/vp0);
		tsmin = tpmin0 - 1.2*(distmax*111.19/vs0);
		tsmax = tpmin0 + 1.2*(distmax*111.19/vs0);
		
		Nps2 = DetermineNprange(ptrig,tpmax,Nst,Nps);
		//printf("%d %lf %lf\n",Nps,told,tpmin0);
		
		if(tpmin < 0.0)tpmin = 0.0;
		if(tsmin < 0.0)tsmin = 0.0;
		if(tpmax > MAXTIME)tpmax = MAXTIME;
		if(tsmax > MAXTIME)tsmax = MAXTIME;

		DetermineNps0range(ptrig0,strig0,tpmin,tpmax,tsmin,tsmax,Nst,NNps);
		
		for(k=0;k<nnn;k++){
			for(l=0;l<8;l++){
				pscounts[k][l] = 0.0; 
			}
		}
		ppp = 0;
		
		//homo model
		if(igrid == 0){
			#pragma omp parallel for shared(pscounts) firstprivate(latref,lonref,elevref,nlon,ndep,dx,dh) private(lat0,lon0,dep,l,i,j,k)
			for(l=0; l<nnn; ++l){
				i = (int)(l/(nlon*ndep));
				j = (int)((l - i*nlon*ndep)/ndep);
				k = l - i*nlon*ndep - j*ndep;
				//In case that searched location is co-located with the station position (gcarc == 0).
				lat0 = latref - rx + i*dx + 0.1234*dx;
				lon0 = lonref - rx + j*dx + 0.1234*dx;
				dep =  k*dh;
				Accounttriggers_homo(lat0,lon0,dep,latref,lonref,elevref,l);
			}
			#pragma omp barrier
		//layer model
		}else{
			#pragma omp parallel for shared(pscounts) firstprivate(latref,lonref,elevref,nlon,ndep,dx,dh) private(lat0,lon0,dep,l,i,j,k)
			for(l=0; l<nnn; ++l){
				i = (int)(l/(nlon*ndep));
				j = (int)((l - i*nlon*ndep)/ndep);
				k = l - i*nlon*ndep - j*ndep;
				//In case that searched location is co-located with the station position (gcarc == 0).
				lat0 = latref - rx + i*dx + 0.1234*dx;
				lon0 = lonref - rx + j*dx + 0.1234*dx;
				dep =  k*dh;
				Accounttriggers_layer(lat0,lon0,dep,latref,lonref,elevref,l);
			}
			#pragma omp barrier
		}
		//only output the resolution file for the first effective event (the first pick should be true)
		if(ires == 1){
			fpr = fopen(RESOLUTION,"w");
			for(k=0;k<nnn;k++){
				fprintf(fpr,"%12.4lf %12.4lf %12.4lf %12.4lf %4d %4d %4d %8.4lf\n",pscounts[k][3],pscounts[k][0],pscounts[k][1],pscounts[k][2],(int)pscounts[k][4],(int)pscounts[k][5],(int)pscounts[k][7],pscounts[k][6]);
			}
			fclose(fpr);
			exit(-1);
		}
		
		//sort pscounts
		Sortpscounts(pscounts,nnn);

		if(ppp >= nppp && pscounts[nnn-1][7] >= nps0  && pscounts[nnn-1][6] <= std0 && pscounts[nnn-1][0] >= stlamin && pscounts[nnn-1][0] <= stlamax && pscounts[nnn-1][1] >= stlomin && pscounts[nnn-1][1] <= stlomax){
			told = pscounts[nnn-1][3];
			ttd = (int)(pscounts[nnn-1][3]/86400);
			tth = (int)((pscounts[nnn-1][3] - ttd*86400)/3600);
			tts = (int)((pscounts[nnn-1][3] - ttd*86400 - tth*3600)/60);
			ttm = pscounts[nnn-1][3] - ttd*86400 - tth*3600 - tts*60;
			sprintf(otime,"%04d %02d %02d %02d:%02d:%06.3f",nyear,nmon,ttd+nday,tth,tts,ttm);

			RELC[mmm].num1 = mmm+1;
			strcpy(RELC[mmm].otime1,otime);
			RELC[mmm].atime1 = pscounts[nnn-1][3];
			RELC[mmm].std1 = pscounts[nnn-1][6];
			RELC[mmm].lat1 = pscounts[nnn-1][0];
			RELC[mmm].lon1 = pscounts[nnn-1][1];
			RELC[mmm].dep1 = pscounts[nnn-1][2];
			RELC[mmm].nofp1 = pscounts[nnn-1][4];
			RELC[mmm].nofs1 = pscounts[nnn-1][5];
			RELC[mmm].ntotal1 = pscounts[nnn-1][7];

			fprintf(stderr,"%5d %10s %12.4lf %8.4lf %12.4lf %12.4lf %12.4lf %4d %4d %4d\n",mmm+1,otime,pscounts[nnn-1][3],pscounts[nnn-1][6],pscounts[nnn-1][0],pscounts[nnn-1][1],pscounts[nnn-1][2],(int)(pscounts[nnn-1][4]),(int)(pscounts[nnn-1][5]),(int)(pscounts[nnn-1][7]));
			mmm++;
			
            iremove = 0;
			//Ispeed is recommended to save time (use a strict threshold)	
			//Otherwise, one event would be associated and located by many initiating P picks. That's huge!
			if(ispeed > 1.0e-5){
				for(k=0;k<Nst;k++){
					lat0 = pscounts[nnn-1][0];
					lon0 = pscounts[nnn-1][1];
					dep  = pscounts[nnn-1][2];
				
					GCinit(&lat0,&lon0,&ST[k].stla,&ST[k].stlo,&GCarc);
					if(igrid == 0){
						tp_cal = sqrt((GCarc*111.19)*(GCarc*111.19) + dep*dep)/vp0 + ST[k].elev/s_vp0;
					}else{
						ih = rint(dep/tdh);
						ig = ih*rint(trx/tdx) + rint(GCarc/tdx);
						tp_cal = TB[ig].ptime + (GCarc - TB[ig].gdist)*TB[ig].prayp + (dep - TB[ig].dep)*TB[ig].phslow + ST[k].elev/s_vp0;
					}

					tp_pre = pscounts[nnn-1][3] + tp_cal;
					tp_pre_b = tp_pre - nrt*ptw/2.0;
					tp_pre_e = tp_pre + nrt*ptw/2.0;


					if(tp_pre_b < 0.0)tp_pre_b = 0.0;
					if(tp_pre_e > MAXTIME)tp_pre_e = MAXTIME;
					
					//To speed up, remove those associated P picks
					for(j=0;j<Nps2;j++){
						if(ptrig[k][j] > tp_pre_b && ptrig[k][j] < tp_pre_e){
							DeleteOne(ptrig,k,Nps,j);
							iremove++;
							break;
						}	
					}
				}
          }
          //make sure the current initiating P is removed
		 if(iremove < 1.0e-5){
		    DeleteOne(ptrig,m,Nps,n);
		 }
        }else{DeleteOne(ptrig,m,Nps,n);}
	}

	fp1 = fopen(CATALOGSEL,"w");
	fp2 = fopen(PHASESEL,"w");
	/*select again to keep the most reliable event within a time window and output the selected events only*/
	nselect = ReselectFinal(RELC,mmm);
	fprintf(stderr,"before selection: %d\n after selection: %d\n",mmm,nselect);

	mag = (double *)malloc(Nst*sizeof(double));
	for(i=0;i<nselect;i++){
		pcount = 0;scount = 0;ps = 0;
		fprintf(fp2,"%5d %20s %12.4lf %8.4lf %12.4lf %12.4lf %12.4lf\n",i+1,RELC[i].otime1,RELC[i].atime1,RELC[i].std1,RELC[i].lat1,RELC[i].lon1,RELC[i].dep1);
		im = 0;
		for(k=0;k<Nst;k++){
			mag[k] = -100;
			lat0 = RELC[i].lat1;
			lon0 = RELC[i].lon1;
			dep  = RELC[i].dep1;
				
			GCinit(&lat0,&lon0,&ST[k].stla,&ST[k].stlo,&GCarc);
			if(igrid == 0){
				tp_cal = sqrt((GCarc*111.19)*(GCarc*111.19) + dep*dep)/vp0 + ST[k].elev/s_vp0;
				ts_cal = sqrt((GCarc*111.19)*(GCarc*111.19) + dep*dep)/vs0 + ST[k].elev/s_vs0;
			}else{
				ih = rint(dep/tdh);
				ig = ih*rint(trx/tdx) + rint(GCarc/tdx);
				tp_cal = TB[ig].ptime + (GCarc - TB[ig].gdist)*TB[ig].prayp + (dep - TB[ig].dep)*TB[ig].phslow + ST[k].elev/s_vp0;
				ts_cal = TB[ig].stime + (GCarc - TB[ig].gdist)*TB[ig].srayp + (dep - TB[ig].dep)*TB[ig].shslow + ST[k].elev/s_vs0;
			}

			tp_pre = RELC[i].atime1 + tp_cal;
			ts_pre = RELC[i].atime1 + ts_cal;
				
			tp_pre_b = tp_pre - nrt*ptw/2.0;
			tp_pre_e = tp_pre + nrt*ptw/2.0;
			ts_pre_b = ts_pre - nrt*stw/2.0;
			ts_pre_e = ts_pre + nrt*stw/2.0;
			if(tp_pre_b < 0.0)tp_pre_b = 0.0;
			if(ts_pre_b < 0.0)ts_pre_b = 0.0;
			if(tp_pre_e > MAXTIME)tp_pre_e = MAXTIME;
			if(ts_pre_e > MAXTIME)ts_pre_e = MAXTIME;
			
			p_mag = -100;
			ptemp = -100;
			for(j=0;j<NNps;j++){
				//rsel*std to remove some picks with large residuals
				if(ptrig0[k][j] > tp_pre_b && ptrig0[k][j] < tp_pre_e && fabs(ptrig0[k][j]-tp_pre) < rsel*RELC[i].std1){
					fprintf(fp2,"%s %s P %12.4lf %12.4lf %12.4e %12.4lf\n",ST[k].net,ST[k].sta,ptrig0[k][j],ptrig0[k][j]-RELC[i].atime1,pamp0[k][j],ptrig0[k][j]-tp_pre);
					pamp = pamp0[k][j];
					p_mag = log(pamp)/log(10) + 1.110*log(GCarc*111.19/100)/log(10) + 0.00189*(GCarc*111.19-100)+3.0;	
					pcount++; ps++;
					ptemp = ptrig0[k][j];
					break;
				}	
			}
			
			s_mag = -100;
			//dtps: to remove some false S picks (they may be P picks but wrongly identified as S picks, it happens!)
			//rsel*std to remove some picks with large residuals
			for(j=0;j<NNps;j++){
				if((ts_pre-tp_pre) > dtps && fabs(ptemp - strig0[k][j]) > dtps && strig0[k][j] > ts_pre_b && strig0[k][j] < ts_pre_e && fabs(strig0[k][j]-ts_pre) < rsel*RELC[i].std1){
					fprintf(fp2,"%s %s S %12.4lf %12.4lf %12.4e %12.4f\n",ST[k].net,ST[k].sta,strig0[k][j],strig0[k][j]-RELC[i].atime1,samp0[k][j],strig0[k][j]-ts_pre);
					samp = samp0[k][j];
					s_mag = log(samp)/log(10) + 1.110*log(GCarc*111.19/100)/log(10) + 0.00189*(GCarc*111.19-100)+3.0;	
					scount++; ps++;
					break;
				}
			}
			
			//if(GCarc*111.19 > 10 && (p_mag > -90 || s_mag > -90)){
			if(p_mag > -90 || s_mag > -90){
				if( p_mag > s_mag){
					mag[im] = p_mag;
					im++;
				}else{
					mag[im] = s_mag;
					im++;
				}
			}
		}
		
		if(im < 2){
			mag_median = -100.0;
			mag_std = -100.0;
		}else{
			mag_median = CalculateMedian(mag,im);
			mag_std = CalculateStd(mag,mag_median,im);
		}
		
		fprintf(fp1,"%5d %20s %12.4lf %8.4lf %12.4lf %12.4lf %12.4lf %8.3lf %8.3lf %4d %4d %4d\n",i+1,RELC[i].otime1,RELC[i].atime1,RELC[i].std1,RELC[i].lat1,RELC[i].lon1,RELC[i].dep1,mag_median,mag_std,pcount,scount,ps);
	}

	fclose(fp1);
	fclose(fp2);
	free(np0_start);
	free(np0_end);
	free(ns0_start);
	free(ns0_end);
	for(i=0;i<Nst;i++){
		free(ptrig[i]);free(strig[i]);
		free(ptrig0[i]);free(strig0[i]);
		free(pamp0[i]);free(samp0[i]);
		free(TGP[i]);free(TGS[i]);
	}
	for(i=0;i<nnn;i++)free(pscounts[i]);
	free(pscounts);
	free(TGP); free(TGS);
	free(ptrig);free(strig);
	free(ptrig0);free(strig0);
	free(pamp0);free(samp0);
	free(ST);free(RELC);
	free(TB);free(mag);
	return 0;
}

int ReselectFinal(SELECT *RELC,int m){
	int i,k,nps;
	char b[50];
	double a,c,d,e,f,g,h,o,p;
	extern int np0,ns0,nps0;

	for(i=0;i<m;i++){
		for(k=(i+1);k<m;k++){
			if(RELC[i].atime1 > RELC[k].atime1){
				a = RELC[i].num1;
				strcpy(b,RELC[i].otime1);
				c = RELC[i].atime1;
				d = RELC[i].std1;
				e = RELC[i].lat1;
				f = RELC[i].lon1;
				g = RELC[i].dep1;
				h = RELC[i].nofp1;
				o = RELC[i].nofs1;
				p = RELC[i].ntotal1;
				
				RELC[i].num1 = RELC[k].num1;
				strcpy(RELC[i].otime1,RELC[k].otime1);
				RELC[i].atime1 = RELC[k].atime1;
				RELC[i].std1 = RELC[k].std1;
				RELC[i].lat1 = RELC[k].lat1;
				RELC[i].lon1 = RELC[k].lon1;
				RELC[i].dep1 = RELC[k].dep1;
				RELC[i].nofp1 = RELC[k].nofp1;
				RELC[i].nofs1 = RELC[k].nofs1;
				RELC[i].ntotal1 = RELC[k].ntotal1;

				RELC[k].num1 = a;
				strcpy(RELC[k].otime1,b);
				RELC[k].atime1 = c;
				RELC[k].std1 = d;
				RELC[k].lat1 = e;
				RELC[k].lon1 = f;
				RELC[k].dep1 = g;
				RELC[k].nofp1 = h;
				RELC[k].nofs1 = o;
				RELC[k].ntotal1 = p;
			}
		}
    }
	
	//exclude the case – one event is associated twice
	for(i=1;i<m;i++){
        for(k=0;k<m,k!= i;k++){
                if(fabs(RELC[i].atime1 - RELC[k].atime1) < 1.05*tint ){
			        if(RELC[i].ntotal1 > RELC[k].ntotal1 || (RELC[i].ntotal1 == RELC[k].ntotal1 && RELC[i].std1 < RELC[k].std1)){
				        RELC[k].atime1 = 1.0e8;
			        }else{
				        RELC[i].atime1 = 1.0e8;
			        }
		    }
        }
	}

	for(i=0;i<m;i++){
		if(RELC[i].nofp1 < np0 || RELC[i].nofs1 < ns0 || RELC[i].ntotal1 < nps0) RELC[i].atime1 = 1.0e8;
	}

	for(i=0;i<m;i++){
		for(k=(i+1);k<m;k++){
			if(RELC[i].atime1 > RELC[k].atime1){
				a = RELC[i].num1;
				strcpy(b,RELC[i].otime1);
				c = RELC[i].atime1;
				d = RELC[i].std1;
				e = RELC[i].lat1;
				f = RELC[i].lon1;
				g = RELC[i].dep1;
				h = RELC[i].nofp1;
				o = RELC[i].nofs1;
				p = RELC[i].ntotal1;
				
				RELC[i].num1 = RELC[k].num1;
				strcpy(RELC[i].otime1,RELC[k].otime1);
				RELC[i].atime1 = RELC[k].atime1;
				RELC[i].std1 = RELC[k].std1;
				RELC[i].lat1 = RELC[k].lat1;
				RELC[i].lon1 = RELC[k].lon1;
				RELC[i].dep1 = RELC[k].dep1;
				RELC[i].nofp1 = RELC[k].nofp1;
				RELC[i].nofs1 = RELC[k].nofs1;
				RELC[i].ntotal1 = RELC[k].ntotal1;

				RELC[k].num1 = a;
				strcpy(RELC[k].otime1,b);
				RELC[k].atime1 = c;
				RELC[k].std1 = d;
				RELC[k].lat1 = e;
				RELC[k].lon1 = f;
				RELC[k].dep1 = g;
				RELC[k].nofp1 = h;
				RELC[k].nofs1 = o;
				RELC[k].ntotal1 = p;
			}
		}
	}
	
	nps = m;
	for(i=0;i<m;i++){
		if(fabs(RELC[i].atime1 - 1.0e8) < 1 && RELC[i-1].atime1 < MAXTIME){
			nps = i;
			break;
		}
	}
	return nps;
}

double CalculateStd(double *arrValue,double median, int max){
    int i;
    double std,temp;

	temp = 0.0;
    for(i = 0; i < max; i++){
		temp += (arrValue[i] - median)*(arrValue[i] - median);
	}

	std = sqrt(temp/max);
    return std;
}

double CalculateMean(double *arrValue, int max){
	double mean = 0.0;
	int i;
	for(i = 0; i < max; i++) mean = mean + arrValue[i];
	return mean/max;
}


double CalculateMedian(double *arrValue, int max){
    double median = 0;
    double *value;
    int i, j;
    double temp;
    value = (double *)malloc(max*sizeof(double));
    for(i = 0; i < max; i++) value[i] = arrValue[i];

    for(i = 0; i < max; i++){
        for(j = 0; j < max - i - 1; j++){
            if(value[j] > value[j + 1]){
                temp = value[j];
                value[j] = value[j + 1];
                value[j + 1] = temp;
            }
        }
    }
    if( (max % 2) == 1){
        median =  value[ (max + 1) / 2 - 1];
    }else{
        median = (value[max / 2] + value[max / 2 - 1]) / 2;
    }
    free(value);
    return median;
}

//Compute great circle (gcarc)
void GCinit(const double *lat1,const double *lon1,const double *lat2,const double *lon2,double *GCarc){
    double x1,yp1,z1,x2,y2,z2;
    double the1,phe1,the2,phe2;
    the1=(90.0-*lat1)*D2R;         /* convert to radius */
    phe1=(*lon1)*D2R;
    the2=(90.0-*lat2)*D2R;
    phe2=(*lon2)*D2R;

    x1=sin(the1)*cos(phe1);
    yp1=sin(the1)*sin(phe1);
    z1=cos(the1);
    x2=sin(the2)*cos(phe2);
    y2=sin(the2)*sin(phe2);
    z2=cos(the2);
    *GCarc=acos((x1)*(x2)+(yp1)*(y2)+(z1)*(z2));
	if( fabs(*GCarc-M_PI) <=1.e-16) {
        fprintf(stderr," The great circle is not determined! (GCinit)\n");
        exit(1);
    }
    *GCarc *= R2D;
}


int Readttime(char *name, TTT *TB, int nmax){
        int i, test;
        FILE *infile;

        test=0;

        while((infile=fopen(name,"r"))==NULL){
            fprintf(stdout,"Can not open file in ReadFile %s\n", name);
            exit(-1);
        }

        for (i = 0; i <= nmax; i++){
            if(fscanf(infile, "%lf %lf %lf %lf %lf %lf %lf %lf %s %s\n",&TB[i].gdist,&TB[i].dep,&TB[i].ptime,&TB[i].stime,&TB[i].prayp,&TB[i].srayp,&TB[i].phslow,&TB[i].shslow,TB[i].pphase,TB[i].sphase)==EOF)test=1;
            if(test==1)break;
        }
        fclose(infile);
		return i;
}

int Readstation(char *name, STATION *ST, int nmax){
        int i, test;
        FILE *infile;

        test=0;

        while((infile=fopen(name,"r"))==NULL){
            fprintf(stdout,"Can not open file in ReadFile %s\n", name);
            exit(-1);
        }

        for (i = 0; i <= nmax; i++){
            if(fscanf(infile, "%lf %lf %s %s %s %lf\n",&ST[i].stlo,&ST[i].stla,ST[i].net,ST[i].sta,ST[i].comp,&ST[i].elev)==EOF)test=1;
            if(test==1)break;
        }
        fclose(infile);
        return i;
}

double Find_min(double **array,int n1, int n2){
	int i,j;
	double amin;
	
	amin = 1.0e8;
	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			if(array[i][j] < amin){
				amin = array[i][j];
			}
		}
	}
	return amin;
}

double Find_max(double **array,int n1, int n2){
	int i,j;
	double amin;
	
	amin = -1.0e8;
	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			if(array[i][j] > amin && array[i][j] <1.0e8){
				amin = array[i][j];
			}
		}
	}
	return amin;
}


void Find_min_loc(double **array, int n1, int n2, double *amin, int *m, int *n){
	int i,j;
	
	*amin = 1.0e8;
	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			if(array[i][j] < *amin){
				*amin = array[i][j];
				*m = i;
				*n = j;
			}
		}
	}
}

// find largest Nps with effective triggers
int DetermineNg(TRIG **ar1, TRIG **ar2, int n1, int n2){
	int i,j,Nps1,Nps0;	
	Nps1 = 0; Nps0 = 0;
	for(i=0;i<n1;i++){
		for(j=1;j<n2;j++){
			if(fabs(ar1[i][j].trig - 1.0e8) < 1 && ar1[i][j-1].trig <= MAXTIME){
				Nps0 = j;
				break;
			}
		}
		if(Nps0 > Nps1){Nps1 = Nps0;}
	}

	for(i=0;i<n1;i++){
		for(j=1;j<n2;j++){
			if(fabs(ar2[i][j].trig - 1.0e8) < 1 && ar2[i][j-1].trig <= MAXTIME){
				Nps0 = j;
				break;
			}
		}
		if(Nps0 > Nps1){Nps1 = Nps0;}
	}
	return Nps1+1;
}

// find largest Np with effective triggers
int DetermineNp(double **ar1,int n1, int n2){
    int i,j,Nps1,Nps0;  
    Nps1 = 0; Nps0 = 0;
    for(i=0;i<n1;i++){
        for(j=1;j<n2;j++){
            if(fabs(ar1[i][j] - 1.0e8) < 1 && ar1[i][j-1] <= MAXTIME){
                Nps0 = j;
                break;
            }   
        }   
        if(Nps0 >= Nps1){Nps1 = Nps0;}
    }   
    return Nps1+1;
}


// find Np range with effective time window
int DetermineNprange(double **ar1,double tpmax,int Nst,int Nps){
	int i,j,Nps0,Nps00;	
	Nps00 = 0; Nps0 = 0;

	// determine the upper bound for tpmax
	for(i=0;i<Nst;i++){
		for(j=1;j<Nps;j++){
			if(ar1[i][j] > tpmax && ar1[i][j-1] < tpmax){
				Nps0 = j;
				break;
			}
		}
		if(Nps0 >= Nps00){Nps00 = Nps0;}
	}
	return Nps00+1;
}


void DetermineNps0range(double **ar1, double **ar2, double tpmin, double tpmax, double tsmin, double tsmax, int Nst, int Nps){
	int i,j;	
	extern int *np0_start, *np0_end, *ns0_start, *ns0_end;

	// determine the lower bound for tpmin and upper bound for tpmax
	for(i=0;i<Nst;i++){
				np0_start[i] = 0;
		for(j=1;j<Nps;j++){
			if(ar1[i][j] > tpmin && ar1[i][j-1] < tpmin){
				np0_start[i] = j-1;
				break;
			}
		}
	}
	for(i=0;i<Nst;i++){
				np0_end[i] = 0;
		for(j=1;j<Nps;j++){
			if(ar1[i][j] > tpmax && ar1[i][j-1] < tpmax){
				np0_end[i] = j;
				break;
			}
		}
	}
	// determine the lower bound for tsmin and upper bound for tsmax
	for(i=0;i<Nst;i++){
				ns0_start[i] = 0;
		for(j=1;j<Nps;j++){
			if(ar2[i][j] > tsmin && ar2[i][j-1] < tsmin){
				ns0_start[i] = j-1;
				break;
			}
		}
	}
	for(i=0;i<Nst;i++){
				ns0_end[i] = 0;
		for(j=1;j<Nps;j++){
			if(ar2[i][j] > tsmax && ar2[i][j-1] < tsmax){
				ns0_end[i] = j;
				break;
			}
		}
	}
}

void SortTriggers(TRIG **tgp,TRIG **tgs,double **array1,double **array2,int m, int n){
	int i,j,k,l;
	double a,b;
	
	for (i = 0; i < m; ++i){
		for (j = 0; j < n; ++j){
			for (k =(j + 1); k < n; ++k){
				if (tgp[i][j].trig > tgp[i][k].trig){
					a = tgp[i][j].trig;
					b = tgp[i][j].weight;
					tgp[i][j].trig = tgp[i][k].trig;
					tgp[i][j].weight = tgp[i][k].weight;
					tgp[i][k].trig = a;
					tgp[i][k].weight = b;
				}
				if (tgs[i][j].trig > tgs[i][k].trig){
					a = tgs[i][j].trig;
					b = tgs[i][j].weight;
					tgs[i][j].trig = tgs[i][k].trig;
					tgs[i][j].weight = tgs[i][k].weight;
					tgs[i][k].trig = a;
					tgs[i][k].weight = b;
				}
            }
        }
    }
	
	for(i=0;i<m;i++){
		array1[i][0] = tgp[i][0].trig;
		array2[i][0] = tgs[i][0].trig;
		for(j=1;j<n;j++){
			if(tgp[i][j].trig - tgp[i][j-1].trig < tint){
			   if(tgp[i][j].weight > tgp[i][j-1].weight){
					array1[i][j] = tgp[i][j].trig;
					array1[i][j-1] = 1.0e8;
			   }else{
					array1[i][j] = 1.0e8;
			   }
			}else{
				array1[i][j] = tgp[i][j].trig;
			}
			
			if(tgs[i][j].trig - tgs[i][j-1].trig < tint){
			   if(tgs[i][j].weight > tgs[i][j-1].weight){
					array2[i][j] = tgs[i][j].trig;
					array2[i][j-1] = 1.0e8;
			   }else{
					array2[i][j] = 1.0e8;
			   }
			}else{
				array2[i][j] = tgs[i][j].trig;
			}

		}
	}
	
	for (i = 0; i < m; ++i){
        for (j = 0; j < n; ++j){
            for (k =(j + 1); k < n; ++k){
                if (array1[i][j] > array1[i][k]){
                    a = array1[i][j];
                    array1[i][j] = array1[i][k];
                    array1[i][k] = a;
                }
                if (array2[i][j] > array2[i][k]){
                    a = array2[i][j];
                    array2[i][j] = array2[i][k];
                    array2[i][k] = a;
                }
            }
        }
    }
}

void SortTriggers0(TRIG **tgp,TRIG **tgs,double **array1,double **array2,double **pamp,double **samp,int m, int n){
	int i,j,k,l;
	double a,b,c;
	
	for (i = 0; i < m; ++i){
		for (j = 0; j < n; ++j){
			for (k =(j + 1); k < n; ++k){
				if (tgp[i][j].trig > tgp[i][k].trig){
					a = tgp[i][j].trig;
					b = tgp[i][j].weight;
					c = tgp[i][j].amp;
					tgp[i][j].trig = tgp[i][k].trig;
					tgp[i][j].weight = tgp[i][k].weight;
					tgp[i][j].amp = tgp[i][k].amp;
					tgp[i][k].trig = a;
					tgp[i][k].weight = b;
					tgp[i][k].amp = c;
				}
				if (tgs[i][j].trig > tgs[i][k].trig){
					a = tgs[i][j].trig;
					b = tgs[i][j].weight;
					c = tgs[i][j].amp;
					tgs[i][j].trig = tgs[i][k].trig;
					tgs[i][j].weight = tgs[i][k].weight;
					tgs[i][j].amp = tgs[i][k].amp;
					tgs[i][k].trig = a;
					tgs[i][k].weight = b;
					tgs[i][k].amp = c;
				}
            }
        }
    }
	
	for(i=0;i<m;i++){
		array1[i][0] = tgp[i][0].trig;
		array2[i][0] = tgs[i][0].trig;
		pamp[i][0] = tgp[i][0].amp;
		samp[i][0] = tgs[i][0].amp;
		for(j=1;j<n;j++){
			if(tgp[i][j].trig - tgp[i][j-1].trig < tint){
			   if(tgp[i][j].weight > tgp[i][j-1].weight){
					array1[i][j] = tgp[i][j].trig;
					pamp[i][j] = tgp[i][j].amp;
					array1[i][j-1] = 1.0e8;
					pamp[i][j-1] = 0.0;
			   }else{
					array1[i][j] = 1.0e8;
					pamp[i][j] = 0.0;
			   }
			}else{
				array1[i][j] = tgp[i][j].trig;
				pamp[i][j] = tgp[i][j].amp;
			}
			
			if(tgs[i][j].trig - tgs[i][j-1].trig < tint){
			   if(tgs[i][j].weight > tgs[i][j-1].weight){
					array2[i][j] = tgs[i][j].trig;
					samp[i][j] = tgs[i][j].amp;
					array2[i][j-1] = 1.0e8;
					samp[i][j-1] = 0.0;
			   }else{
					array2[i][j] = 1.0e8;
					samp[i][j] = 0.0;
			   }
			}else{
				array2[i][j] = tgs[i][j].trig;
				samp[i][j] = tgs[i][j].amp;
			}
		}
	}
	
	for (i = 0; i < m; ++i){
        for (j = 0; j < n; ++j){
            for (k =(j + 1); k < n; ++k){
                if (array1[i][j] > array1[i][k]){
                    a = array1[i][j];
					b = pamp[i][j];
                    array1[i][j] = array1[i][k];
					pamp[i][j] = pamp[i][k];
                    array1[i][k] = a;
					pamp[i][k] = b;
                }
                if (array2[i][j] > array2[i][k]){
                    a = array2[i][j];
					b = samp[i][j];
                    array2[i][j] = array2[i][k];
					samp[i][j] = samp[i][k];
                    array2[i][k] = a;
					samp[i][k] = b;
                }
            }
        }
    }
}

void DeleteOne(double **array,int Nst0,int Nps0,int Nloc){
	int i;
	for (i = Nloc; i < Nps0-1; i++){array[Nst0][i] = array[Nst0][i+1];}
	array[Nst0][Nps0-1] = 1.0e8;
}

void Sortpscounts(double **pscounts0,int np){
	int i,j,k;
	double a,b,c,d,e,f,g,h;
	
	for(i=0;i<np;i++){
		for (j =(i + 1); j < np; j++){
                if (pscounts0[i][7] > pscounts0[j][7] || (pscounts0[i][7] == pscounts0[j][7] && pscounts0[i][6] < pscounts0[j][6])){

					a = pscounts0[i][0];
					b = pscounts0[i][1];
					c = pscounts0[i][2];
					d = pscounts0[i][3];
					e = pscounts0[i][4];
					f = pscounts0[i][5];
					g = pscounts0[i][6];
					h = pscounts0[i][7];
					
					for(k=0;k<8;k++){
						pscounts0[i][k] = pscounts0[j][k];
					}

					pscounts0[j][0] = a;
					pscounts0[j][1] = b;
					pscounts0[j][2] = c;
					pscounts0[j][3] = d;
					pscounts0[j][4] = e;
					pscounts0[j][5] = f;
					pscounts0[j][6] = g;
					pscounts0[j][7] = h;

                }
		}
	}
}


void Accounttriggers_homo(double lat0,double lon0,double dep,double latref,double lonref,double elevref, int l){
	int pcount,scount,ps;
	int i,j,k;
	double GCarc,median,std,ptemp;
	double tp0_cal,tp_cal,ts_cal,tp_pre,ts_pre,tp_pre_b,tp_pre_e,ts_pre_b,ts_pre_e;
	extern double vp0,vs0,s_vp0,s_vs0;
	extern double nrt,ptw,stw,tpmin0;
	extern int np0,ns0,nps0,Nst,NNps,ppp;
	extern double **ptrig0,**strig0;
	extern int *np0_start, *np0_end, *ns0_start, *ns0_end;
	extern STATION *ST;
	extern double **pscounts;
	double *torg;
	extern double dtps;
	

	pcount = 0;scount = 0;ps = 0;
		
	torg = (double *)malloc(2*Nst*sizeof(double));   
	for(k=0;k<2*Nst;k++)torg[k] = 0.0;
	
	GCinit(&lat0,&lon0,&latref,&lonref,&GCarc);
	tp0_cal = sqrt((GCarc*111.19)*(GCarc*111.19) + dep*dep)/vp0 + elevref/s_vp0;


	for(i=0;i<Nst;i++){
		GCinit(&lat0,&lon0,&ST[i].stla,&ST[i].stlo,&GCarc);
		tp_cal = sqrt((GCarc*111.19)*(GCarc*111.19) + dep*dep)/vp0 + ST[i].elev/s_vp0;
		ts_cal = sqrt((GCarc*111.19)*(GCarc*111.19) + dep*dep)/vs0 + ST[i].elev/s_vs0;
						
		tp_pre = tpmin0 - tp0_cal + tp_cal;
		ts_pre = tpmin0 - tp0_cal + ts_cal;

		tp_pre_b = tp_pre - nrt*ptw/2.0;
		tp_pre_e = tp_pre + nrt*ptw/2.0;
		ts_pre_b = ts_pre - nrt*stw/2.0;
		ts_pre_e = ts_pre + nrt*stw/2.0;
		
		if(tp_pre_b < 0.0)tp_pre_b = 0.0;
		if(ts_pre_b < 0.0)ts_pre_b = 0.0;
		if(tp_pre_e > MAXTIME)tp_pre_e = MAXTIME;
		if(ts_pre_e > MAXTIME)ts_pre_e = MAXTIME;
		
		ptemp = -100;
		for(j=np0_start[i];j<np0_end[i];j++){
			if(ptrig0[i][j] > tp_pre_b && ptrig0[i][j] < tp_pre_e){
				torg[ps] = ptrig0[i][j] - tp_cal;
				pcount = pcount + 1; ps = ps + 1;
				ptemp = ptrig0[i][j];
				break;
			}
		}
		
		//dtps: to remove some false S picks (they may be P picks but wrongly identified as S picks, it happens!)
		for(j=ns0_start[i];j<ns0_end[i];j++){
			if((ts_pre - tp_pre) > dtps && fabs(ptemp - strig0[i][j]) > dtps && strig0[i][j] > ts_pre_b && strig0[i][j] < ts_pre_e){
				torg[ps] = strig0[i][j] - ts_cal;
				scount = scount + 1; ps = ps + 1;
				break;
			}
		}
	}
	
	if(pcount >= np0 && scount >= ns0 && ps >= nps0){
		median = CalculateMedian(torg,ps);
		//median = CalculateMean(torg,ps);
		std = CalculateStd(torg,median,ps);
		pscounts[l][0] = lat0;
		pscounts[l][1] = lon0;
		pscounts[l][2] = dep;
		pscounts[l][3] = median;
		pscounts[l][4] = pcount;
		pscounts[l][5] = scount;
		pscounts[l][6] = std;
		pscounts[l][7] = ps;
		ppp++;
	}else{
		for(j=0;j<8;++j)pscounts[l][j] = 0.0;
	}
	free(torg);
}


void Accounttriggers_layer(double lat0,double lon0,double dep,double latref,double lonref,double elevref, int l){
	int pcount,scount,ps;
	int i,j,k,ig,ih;
	double GCarc,median,std,ptemp;
	double tp0_cal,tp_cal,ts_cal,tp_pre,ts_pre,tp_pre_b,tp_pre_e,ts_pre_b,ts_pre_e;
	extern double vp0,vs0,s_vp0,s_vs0;
	extern double nrt,ptw,stw,tpmin0;
	extern int np0,ns0,nps0,Nst,NNps,ppp;
	extern double **ptrig0,**strig0;
	extern int *np0_start, *np0_end, *ns0_start, *ns0_end;
	extern STATION *ST;
	extern double **pscounts;
	extern double trx,tdx,tdh,dtps;
	double *torg;
	

	pcount = 0;scount = 0;ps = 0;
		
	torg = (double *)malloc(2*Nst*sizeof(double));   
	for(k=0;k<2*Nst;k++)torg[k] = 0.0;
	
	GCinit(&lat0,&lon0,&latref,&lonref,&GCarc);
	ih = round(dep/tdh);
	ig = ih*rint(trx/tdx) + rint(GCarc/tdx);
	tp0_cal = TB[ig].ptime + (GCarc - TB[ig].gdist)*TB[ig].prayp + (dep - TB[ig].dep)*TB[ig].phslow + elevref/s_vp0;
	
	for(i=0;i<Nst;i++){
		GCinit(&lat0,&lon0,&ST[i].stla,&ST[i].stlo,&GCarc);
		ih = rint(dep/tdh);
		ig = ih*rint(trx/tdx) + rint(GCarc/tdx);
		tp_cal = TB[ig].ptime + (GCarc - TB[ig].gdist)*TB[ig].prayp + (dep - TB[ig].dep)*TB[ig].phslow + ST[i].elev/s_vp0;
		ts_cal = TB[ig].stime + (GCarc - TB[ig].gdist)*TB[ig].srayp + (dep - TB[ig].dep)*TB[ig].shslow + ST[i].elev/s_vs0;
						
		tp_pre = tpmin0 - tp0_cal + tp_cal;
		ts_pre = tpmin0 - tp0_cal + ts_cal;

		tp_pre_b = tp_pre - nrt*ptw/2.0;
		tp_pre_e = tp_pre + nrt*ptw/2.0;
		ts_pre_b = ts_pre - nrt*stw/2.0;
		ts_pre_e = ts_pre + nrt*stw/2.0;


		if(tp_pre_b < 0.0)tp_pre_b = 0.0;
		if(ts_pre_b < 0.0)ts_pre_b = 0.0;
		if(tp_pre_e > MAXTIME)tp_pre_e = MAXTIME;
		if(ts_pre_e > MAXTIME)ts_pre_e = MAXTIME;
		
		ptemp = -100;
		for(j=np0_start[i];j<np0_end[i];j++){
			if(ptrig0[i][j] > tp_pre_b && ptrig0[i][j] < tp_pre_e){
				torg[ps] = ptrig0[i][j] - tp_cal;
				pcount = pcount + 1; ps = ps + 1;
				ptemp = ptrig0[i][j];
				break;
			}
		}
		
		//dtps: to remove some false S picks (they may be P picks but wrongly identified as S picks, it happens!)
		for(j=ns0_start[i];j<ns0_end[i];j++){
			if((ts_pre - tp_pre) > dtps && fabs(ptemp - strig0[i][j]) > dtps && strig0[i][j] > ts_pre_b && strig0[i][j] < ts_pre_e){
				torg[ps] = strig0[i][j] - ts_cal;
				scount = scount + 1; ps = ps + 1;
				break;
			}
		}
	}
	
	if(pcount >= np0 && scount >= ns0 && ps >= nps0){
		//median = CalculateMean(torg,ps);
		median = CalculateMedian(torg,ps);
		std = CalculateStd(torg,median,ps);
		pscounts[l][0] = lat0;
		pscounts[l][1] = lon0;
		pscounts[l][2] = dep;
		pscounts[l][3] = median;
		pscounts[l][4] = pcount;
		pscounts[l][5] = scount;
		pscounts[l][6] = std;
		pscounts[l][7] = ps;
		ppp++;
	}else{
		for(j=0;j<8;++j)pscounts[l][j] = 0.0;
	}
	free(torg);
}
