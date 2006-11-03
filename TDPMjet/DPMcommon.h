#ifndef ROOT_DPMCommon
#define ROOT_DPMCommon

#ifndef __CFORTRAN_LOADED
//*KEEP,cfortran.
#include "cfortran.h"
//*KEND.
#endif

extern "C" {
/*========================================================*/
/* COMMON/DTEVNO/NEVENT,ICASCA				  */
/*--------------------------------------------------------*/
typedef struct {
   Int_t    nevent;
   Int_t    icasca;
} DtevnoCommon;

#define DTEVNO COMMON_BLOCK(DTEVNO,dtevno)
COMMON_BLOCK_DEF(DtevnoCommon,DTEVNO);

/**********************************************************/
/*           D E S C R I P T I O N :                      */
/*--------------------------------------------------------*/
/* 			Event flag		          */
/*========================================================*/

/*========================================================*/
/* COMMON/DTEVT1/NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),*/
/*	JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),		  */
/*	PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)	  */
/*--------------------------------------------------------*/
typedef struct {
   Int_t    nhkk;
   Int_t    nevhkk;
   Int_t    isthkk[200000];
   Int_t    idhkk[200000];
   Int_t    jmohkk[200000][2];
   Int_t    jdahkk[200000][2];
   Double_t  phkk[200000][5];
   Double_t  vhkk[200000][4];
   Double_t  whkk[200000][4];
} Dtevt1Common;

#define DTEVT1 COMMON_BLOCK(DTEVT1,dtevt1)
COMMON_BLOCK_DEF(Dtevt1Common,DTEVT1);

/**********************************************************/
/*           D E S C R I P T I O N :                      */
/*--------------------------------------------------------*/
/* 			Event history		          */
/* 						          */
/* NHKK - number of entries in common block		  */
/* NEVHKK - number of the event				  */
/* ISTHKK(i) - status code for entry i with following     */
/*	       meanings:			  	  */
/*	     = 1 final state particle produced in         */
/*	     	photon-/hadron-/nucleon-nucleon collisions*/
/*	     	 or in intranuclear cascade proc.         */
/*	     =-1 nucleons, deuterons, H3, He3, He4        */
/*	     	 evaporated from excited nucleus and      */
/*	     	 photons produced in nuclear deexcitation */
/*	     	 processes;			          */
/*	     = 1001  residual nucleus (ground state).     */
/* IDHKK(i) - particle identity according to PDG code;    */
/*	      for nuclei (evaporation products and     	  */
/* 	      residual nucleus): IDHKK(IHKK)=80000     	  */
/* JMOHKK(1,i) - pointer to the position where the mother */
/*		 is stored; the initial value is 0  	  */
/* JMOHKK(2,i) - pointer to the position of the 2nd mother*/
/* 		 Normally only 1 mother exists, in which  */
/* 		case the value 0 is used. In cluster      */
/*		fragmentation models, the 2 mothers would */
/*		correspond to the q and qbar which join to*/
/*		form a cluster. In string fragmentation,  */
/*		the two mothers of a particle produced in */
/*		the fragmentation would be the 2 endpoints*/
/*		of the string.				  */
/* JDAHKK(1,i) - pointer to the position of the 1st       */
/*		daughter; if an entry has not decayed =0. */
/* JDAHKK(1,i) - pointer to the position of the last      */
/*		daughter; if an entry has not decayed =0. */
/* PHKK(1,i) - momentum in x direction in GeV/c           */
/* PHKK(2,i) - momentum in y direction in GeV/c           */
/* PHKK(3,i) - momentum in z direction in GeV/c           */
/* PHKK(4,i) - energy in GeV			          */
/* PHKK(5,i) - mass in GeV/c^2; for space-like partons    */
/*	    it is allowed to use a mass<0, according      */
/*	    PHKK(5,IHKK) = -sqrt(-m^2)  	          */ 
/* VHKK(1,i) - production vertex in x position in mm      */
/* VHKK(2,i) - production vertex in y position in mm      */
/* VHKK(3,i) - production vertex in z position in mm      */
/* VHKK(4,i) - production time in mm/c (=3.33*10^(-12)s   */
/* WHKK(I,i) - gives positions and times in projectile    */ 
/*		frame, the chains are created on the posi-*/
/*		tions of the projectile nucleons in the   */
/*		projectile frame (target nucleons in target*/
/*		frame); both positions are therefore not  */
/*		completely consistent. The times in the   */
/*		projectile frame are obtained by a Lorentz*/
/*		transformation from the LAB system.	  */
/*========================================================*/

/*========================================================*/
/* COMMON/DTEVT2/IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),*/
/*               IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10),   */
/*               IHIST(2,NMXHKK)			  */
/*--------------------------------------------------------*/
typedef struct {
   Int_t    idres[200000];
   Int_t    idxres[200000];
   Int_t    nobam[200000];
   Int_t    idbam[200000];
   Int_t    idch[200000];
   Int_t    npoint[10];
   Int_t    ihist[200000][2];
} Dtevt2Common;

#define DTEVT2 COMMON_BLOCK(DTEVT2,dtevt2)
COMMON_BLOCK_DEF(Dtevt2Common,DTEVT2);

/**********************************************************/
/*           D E S C R I P T I O N :                      */
/*--------------------------------------------------------*/
/* 		Extended event history		          */
/* 						          */
/* NMXHKK - max. num. of entries (partons/particles) that */
/* 	  	can be stored in the common block	  */
/* IDRES(IHKK) - mass num. A in case of nuclear fragments */
/*		or residual nuclei (IDHKK(IHKK)=80000).   */
/* IDXRES(IHKK) - charge Zin case of nuclear fragments    */
/*		or residual nuclei (IDHKK(IHKK)=80000).   */
/* NOBAM(IHKK) =1 for particles from proj. fragmentation  */
/* 	       =2 for particles from target fragmentation.*/
/* IDBAM(IHKK) - internal dpmjet particle code(BAMJET code)*/
/*========================================================*/

}

/*========================================================*/
/* COMMON/DTPRTA/IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG*/
/*--------------------------------------------------------*/
typedef struct {
   Int_t    it;
   Int_t    itz;
   Int_t    ip;
   Int_t    ipz;
   Int_t    ijproj;
   Int_t    ibproj;
   Int_t    ijtarg;
   Int_t    ibtarg;
} DtprtaCommon;

#define DTPRTA COMMON_BLOCK(DTPRTA,dtprta)
COMMON_BLOCK_DEF(DtprtaCommon,DTPRTA);
/**********************************************************/
/*	     D E S C R I P T I O N :			  */
/*--------------------------------------------------------*/
/* IT, ITZ - nucleon/atomic number of target nucleus	  */
/* IP, IPZ - nucleon/atomic number of projectile nucleus  */
/*	   for incident hadrons IP=IPZ=1		  */
/*========================================================*/

/*========================================================*/
/* COMMON /DTIMPA/ BIMIN,BIMAX,XSFRAC,ICENTR		  */
/*--------------------------------------------------------*/
typedef struct {
   Double_t    bimin;
   Double_t    bimax;
   Double_t    xsfrac;
   Double_t    icent;
} DtimpaCommon;

#define DTIMPA COMMON_BLOCK(DTIMPA,dtimpa)
COMMON_BLOCK_DEF(DtimpaCommon,DTIMPA);
/**********************************************************/
/*	     D E S C R I P T I O N :			  */
/*--------------------------------------------------------*/
/* BIMIN, BIMAX - min., max. b values (default bmin = 0)  */
/* XSFRAC - fraction of x-section (default: 1)  	  */
/* ICENTR =1. central production forced (default: 0)	  */
/*	<0 && >-100 -> bmin = BIMIN, bmax = BIMAX	  */
/*	<-99 -> fraction of x-sec. = XSFRAC		  */
/*	=-1. -> evaporation/fzc suppressed		  */
/*	<-1. -> evaporation/fzc suppressed		  */
/*========================================================*/


/*========================================================*/
/* COMMON /DTGLCP/RPROJ,RTARG,BIMPAC,NWTSAM,NWASAM,NWBSAM,*/
/*		  NWTACC,NWAACC,NWBACC		  	  */
/*--------------------------------------------------------*/
typedef struct {
   Double_t	rproj;
   Double_t	rtarg;
   Double_t	bimpac;
   Int_t	nwtsam;
   Int_t	nwasam;
   Int_t	nwbsam;
   Int_t	nwtacc;
   Int_t	nwtaac;
   Int_t	nwtbac;
} DtglcpCommon;

#define DTGLCP COMMON_BLOCK(DTGLCP,dtglcp)
COMMON_BLOCK_DEF(DtglcpCommon,DTGLCP);
/**********************************************************/
/*	     D E S C R I P T I O N :			  */
/*--------------------------------------------------------*/
/* RPROJ  = radius of projectile nucleus		  */
/* RPROJ  = radius of target nucleus			  */
/* BIMPAC = impact parameter of the collision		  */
/* NWTSAM = total number of wounded nucleons		  */
/* NWASAM = number of wounded nucleons in projectile	  */
/* NWBSAM = number of wounded nucleons in target	  */
/* NWTACC = total number of interacting nucleons	  */
/* NWTAAC = total number of interacting nucleons in proj. */
/* NWTBAC = total number of interacting nucleons in target*/
/*========================================================*/
/*========================================================*/
/* COMMON /POPRCS/IPROCE,IDNODF,IDIFR1,IDIFR2,IDDPOM,     */
/*		  IPRON(15,4)   		  	          */
/*--------------------------------------------------------*/
typedef struct {
   Int_t	iproce;  
   Int_t        idnodf;
   Int_t        idifr1;
   Int_t        idifr2;
   Int_t        iddpom;
   Int_t        ipron[15][4];
} PoprcsCommon;

#define POPRCS COMMON_BLOCK(POPRCS,poprcs)
COMMON_BLOCK_DEF(PoprcsCommon,POPRCS);


/**********************************************************/
/*	     D E S C R I P T I O N :			  */
/*--------------------------------------------------------*/

#endif
