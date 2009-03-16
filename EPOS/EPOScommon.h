/*
 *###################################################################
 *#        EPOS 1.67     K. WERNER, T. PIEROG, S. PORTEBOEUF.       #
 *#                      Contact: werner@subatech.in2p3.fr          #
 *###################################################################
 *
 * EPOScommon.h
 * 
 * Definitions of common blocks
 *
 *      Author: Piotr Ostrowski, postrow@if.pw.edu.pl
 */

#ifndef EPOSCOMMON_H_
#define EPOSCOMMON_H_

#ifndef __CFORTRAN_LOADED
//*KEEP,cfortran.
#include <cfortran.h>
//*KEND.
#endif

#define MXJERR 10
#define MXPTL 200000

/* Common blocks */
extern "C" {
	//--------------------------------------------------------------------------
	//                   epos event common block
	//---------------------------------------------------------------------------
	//
	//      common/cevt/phievt,nevt,bimevt,kolevt,koievt,pmxevt,egyevt,npjevt
	//     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
	//     *,xbjevt,qsqevt,nglevt,zppevt,zptevt,minfra,maxfra
	//
	//     nevt .......... error code. 1=valid event, 0=invalid event
	//     bimevt ........ absolute value of impact parameter
	//     phievt ........ angle of impact parameter
	//     kolevt ........ number of collisions
	//     koievt ........ number of inelastic collisions
	//     pmxevt ........ reference momentum
	//     egyevt ........ pp cm energy (hadron) or string energy (lepton)
	//     npjevt ........ number of primary projectile participants
	//     ntgevt ........ number of primary target participants
	//     npnevt ........ number of primary projectile neutron spectators
	//     nppevt ........ number of primary projectile proton spectators
	//     ntnevt ........ number of primary target neutron spectators
	//     ntpevt ........ number of primary target proton spectators
	//     jpnevt ........ number of absolute projectile neutron spectators
	//     jppevt ........ number of absolute projectile proton spectators
	//     jtnevt ........ number of absolute target neutron spectators
	//     jtpevt ........ number of absolute target proton spectators
	//     xbjevt ........ bjorken x for dis
	//     qsqevt ........ q**2 for dis
	//     sigtot ........ total cross section
	//     nglevt ........ number of collisions acc to  Glauber
	//     zppevt ........ average Z-parton-proj
	//     zptevt ........ average Z-parton-targ

typedef struct {
	Float_t phievt;
	Int_t nevt;
	Float_t bimevt;
	Int_t kolevt,koievt;
	Float_t pmxevt,egyevt;
	Int_t npjevt,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt;
	Float_t xbjevt,qsqevt;
	Int_t nglevt;
	Float_t zppevt,zptevt;
	Int_t minfra,maxfra;
} CevtCommon;

#define cevt COMMON_BLOCK(cevt,cevt)
COMMON_BLOCK_DEF(CevtCommon,cevt);

//---------------------------------------------------------------------------
//                   epos particle list common block
//---------------------------------------------------------------------------
//
//     common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
//    *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
//    *,xorptl(4,mxptl),ibptl(4,mxptl),ityptl(mxptl)
//     common/c1ptl/itsptl(mxptl)
//
//     nptl .......... current particle index (=number of ptls stored)
//     idptl(i) ...... particle id
//     pptl(1,i) ..... x-component of particle momentum
//     pptl(2,i) ..... y-component of particle momentum
//     pptl(3,i) ..... z-component of particle momentum
//     pptl(4,i) ..... particle energy
//     pptl(5,i) ..... particle mass
//     iorptl(i) ..... particle number of father (if .le. 0 : no father)
//     jorptl(i) ..... particle number of mother (if .le. 0 : no mother)
//     istptl(i) ..... status: 40 and 41 : Remnant
//                             30 and 31 : Pomeron
//                             20 and 21 : Parton
//                             10 and 11 : Droplet
//                             00 and 01 : Particle
//                            last digit = 0 : last generation
//                            last digit = 1 : not last generation
//     xorptl(1,i) ... x-component of formation point
//     xorptl(2,i) ... y-component of formation point
//     xorptl(3,i) ... z-component of formation point
//     xorptl(4,i) ... formation time
//     tivptl(1,i) ... formation time (always in the pp-cms!)
//     tivptl(2,i) ... destruction time (always in the pp-cms!)
//     ityptl(i)  .... type of particles origin:
//                         10-19: target
//                         20-29: soft Pom
//                         30-39: hard Pom
//                         40-49: projectile
//                         50: string, droplet
//     itsptl(i) ..... string type of particles origin (if string)
//     common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
//    *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
//    *,xorptl(4,mxptl),ibptl(4,mxptl),ityptl(mxptl)

typedef struct {
	Int_t nptl;
	Float_t pptl[MXPTL][5];
	Int_t iorptl[MXPTL];
	Int_t idptl[MXPTL];
	Int_t istptl[MXPTL];
	Float_t tivptl[MXPTL][2];
	Int_t ifrptl[MXPTL][2];
	Int_t jorptl[MXPTL];
	Float_t xorptl[MXPTL][4];
	Int_t ibptl[MXPTL][4];
	Int_t ityptl[MXPTL];
} CptlCommon;

#define cptl COMMON_BLOCK(cptl,cptl)
COMMON_BLOCK_DEF(CptlCommon,cptl);

typedef struct {
	Int_t itsptl[MXPTL];
} ClptCommon;

#define clpt COMMON_BLOCK(clpt,clpt)
COMMON_BLOCK_DEF(ClptCommon,clpt);

/*=========================*/
/*common/copen/nopen,nopenr*/
/*-------------------------*/
typedef struct {
	Int_t	nopen;
	Int_t	nopenr;
} CopenCommon;

#define copen COMMON_BLOCK(copen,copen)
COMMON_BLOCK_DEF(CopenCommon,copen);

/*================================================================*/
/*common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi*/
/*----------------------------------------------------------------*/
typedef struct {
	Int_t iprmpt, ish, ishsub, irandm, irewch, iecho, modsho, idensi;
} Prnt1Common;

#define prnt1 COMMON_BLOCK(prnt1,prnt1)
COMMON_BLOCK_DEF(Prnt1Common,prnt1);

/*========================*/
/*common/appli/iappl,model*/
/*------------------------*/
typedef struct {
	Int_t iappl, model;
} AppliCommon;

#define appli COMMON_BLOCK(appli,appli)
COMMON_BLOCK_DEF(AppliCommon,appli);

/*=========================================*/
/*common/enrgy/egymin,egymax,elab,ecms,ekin*/
/*-----------------------------------------*/
typedef struct {
	Float_t egymin, egymax, elab, ecms, ekin;
} EnrgyCommon;

#define enrgy COMMON_BLOCK(enrgy,enrgy)
COMMON_BLOCK_DEF(EnrgyCommon,enrgy);

/*=============================================*/
/*common/lept1/engy,elepti,elepto,angmue,icinpu*/
/*---------------------------------------------*/
typedef struct {
	Float_t engy,elepti,elepto,angmue;
	Int_t icinpu;
} Lept1Common;

#define lept1 COMMON_BLOCK(lept1,lept1)
COMMON_BLOCK_DEF(Lept1Common,lept1);

/*=====================================================*/
/*common/ebin/noebin,engmin,engmax,nrebin,iologe,iologl*/
/*-----------------------------------------------------*/
typedef struct {
	Int_t noebin;
	Float_t engmin,engmax;
	Int_t nrebin,iologe,iologl;
} EbinCommon;

#define ebin COMMON_BLOCK(ebin,ebin)
COMMON_BLOCK_DEF(EbinCommon,ebin);

/*==========================================*/
/*common/events/nevent,nfull,nfreeze,ninicon*/
/*------------------------------------------*/
typedef struct {
	Int_t nevent,nfull,nfreeze,ninicon;
} EventsCommon;

#define events COMMON_BLOCK(events,events)
COMMON_BLOCK_DEF(EventsCommon,events);

/*=============================================*/
/*common/sprio/ispherio,icotabm,icotabr,icocore*/
/*---------------------------------------------*/
typedef struct {
	Int_t ispherio,icotabm,icotabr,icocore;
} SprioCommon;

#define sprio COMMON_BLOCK(sprio,sprio)
COMMON_BLOCK_DEF(SprioCommon,sprio);

/*============================================================*/
/*common/accum/imsg,jerr(mxjerr),ntevt,nrevt,naevt,nrstr,nrptl*/
/*------------------------------------------------------------*/
typedef struct {
	Int_t imsg,jerr[MXJERR],ntevt,nrevt,naevt,nrstr,nrptl;
} AccumCommon;

#define accum COMMON_BLOCK(accum,accum)
COMMON_BLOCK_DEF(AccumCommon,accum);

/*==================================================================*/
/*common/othe1/istore,istmax,gaumx,irescl,ntrymx,nclean,iopdg,ioidch*/
/*------------------------------------------------------------------*/
typedef struct {
	Int_t istore,istmax;
	Double_t gaumx;
	Int_t irescl,ntrymx,nclean,iopdg,ioidch;
} Othe1Common;

#define othe1 COMMON_BLOCK(othe1,othe1)
COMMON_BLOCK_DEF(Othe1Common,othe1);

/*========================================*/
/*common/othe2/ifrade,iframe,idecay,jdecay*/
/*----------------------------------------*/
typedef struct {
	Int_t ifrade,iframe,idecay,jdecay;
} Othe2Common;

#define othe2 COMMON_BLOCK(othe2,othe2)
COMMON_BLOCK_DEF(Othe2Common,othe2);





/*====================================================*/
/*common/nucl1/laproj,maproj,latarg,matarg,core,fctrmx*/
/*----------------------------------------------------*/
typedef struct {
	Int_t laproj,maproj,latarg,matarg;
	Float_t core,fctrmx;
} Nucl1Common;

#define nucl1 COMMON_BLOCK(nucl1,nucl1)
COMMON_BLOCK_DEF(Nucl1Common,nucl1);

/*========================================*/
/*common/nucl2/bmaxim,bminim,phimax,phimin*/
/*----------------------------------------*/
typedef struct {
	Float_t bmaxim,bminim,phimax,phimin;
} Nucl2Common;

#define nucl2 COMMON_BLOCK(nucl2,nucl2)
COMMON_BLOCK_DEF(Nucl2Common,nucl2);

/*==============================================================*/
/* character*80 fnch,fnhi,fndt,fnii,fnid,fnie,fnrj,fnmt
  * ,fngrv,fncp,fnnx,fncs,fndr,fnhy
   common/fname/  fnch, fnhi, fndt, fnii, fnid, fnie, fnrj, fnmt*/
/*--------------------------------------------------------------*/
typedef struct {
	 char fnch[80], fnhi[80], fndt[80], fnii[80], fnid[80], fnie[80], fnrj[80];
	 char fnmt[80], fngrv[80], fncp[80], fnnx[80], fncs[80], fndr[80], fnhy[80];
} FnameCommon;

#define fname COMMON_BLOCK(fname,fname)
COMMON_BLOCK_DEF(FnameCommon,fname);

/*=============================================================*/
/*common/nfname/nfnch,nfnhi,nfndt,nfnii,nfnid,nfnie,nfnrj,nfnmt*/
/**,nfngrv,nfncp,nfnnx,nfncs,nfndr,nfnhy                       */
/*-------------------------------------------------------------*/
typedef struct {
    Int_t nfnch,nfnhi,nfndt,nfnii,nfnid,nfnie,nfnrj,nfnmt;
    Int_t nfngrv,nfncp,nfnnx,nfncs,nfndr,nfnhy;
} NfnameCommon;

#define nfname COMMON_BLOCK(nfname,nfname)
COMMON_BLOCK_DEF(NfnameCommon,nfname);

/*character*80 fndat,fnncs,fnIIdat,fnIIncs*/
/*============================================*/
/*common/qgsfname/  fndat, fnncs, ifdat, ifncs*/
/*--------------------------------------------*/
typedef struct {
	char fndat[80], fnncs[80];
	Int_t ifdat, ifncs;
} QgsfnameCommon;

#define qgsfname COMMON_BLOCK(qgsfname,qgsfname)
COMMON_BLOCK_DEF(QgsfnameCommon,qgsfname);

/*====================================================*/
/*common/qgsIIfname/fnIIdat, fnIIncs, ifIIdat, ifIIncs*/
/*----------------------------------------------------*/
typedef struct {
	char fnIIdat[80], fnIIncs[80];
	Int_t ifIIdat, ifIIncs;
} QgsIIfnameCommon;

#define qgsiifname COMMON_BLOCK(qgsiifname,qgsiifname)
COMMON_BLOCK_DEF(QgsIIfnameCommon,qgsiifname);

/*================================*/
/*common/qgsnfname/ nfndat, nfnncs*/
/*--------------------------------*/
typedef struct {
	Int_t nfndat, nfnncs;
} QgsnfnameCommon;

#define qgsnfname COMMON_BLOCK(qgsnfname,qgsnfname)
COMMON_BLOCK_DEF(QgsnfnameCommon,qgsnfname);

/*======================================*/
/*common/qgsIInfname/ nfnIIdat, nfnIIncs*/
/*--------------------------------------*/
typedef struct {
	Int_t nfnIIdat, nfnIIncs;
} QgsIInfnameCommon;

#define qgsiinfname COMMON_BLOCK(qgsiinfname,qgsiinfname)
COMMON_BLOCK_DEF(QgsIInfnameCommon,qgsiinfname);

/*================================*/
/*parameter(mxnody=200)           */
/*common/nodcy/nrnody,nody(mxnody)*/
/*--------------------------------*/
#define MXNODY 200
typedef struct {
	Int_t nrnody, nody[MXNODY];
} NodcyCommon;

#define nodcy COMMON_BLOCK(nodcy,nodcy)
COMMON_BLOCK_DEF(NodcyCommon,nodcy);

/*===============================================================*/
/*common/hadr3/iregge,isopom,ishpom,iscreen,nprmax,inueff,irmdrop*/
/*---------------------------------------------------------------*/
typedef struct {
	Int_t iregge,isopom,ishpom,iscreen,nprmax,inueff,irmdrop;
} Hadr3Common;

#define hadr3 COMMON_BLOCK(hadr3,hadr3)
COMMON_BLOCK_DEF(Hadr3Common,hadr3);

/*===============================================================*/
/*common/hadr4/alppom,slopom,gamhad(4),r2had(4),chad(4),wdiff(4) */
/*&            ,gamtil,facdif,r2hads(4),gamhads(4),slopoms,isplit*/
/*---------------------------------------------------------------*/
typedef struct {
	Float_t alppom,slopom,gamhad[4],r2had[4],chad[4],wdiff[4];
	Float_t gamtil,facdif,r2hads[4],gamhads[4],slopoms;
	Int_t isplit;
} Hadr4Common;

#define hadr4 COMMON_BLOCK(hadr4,hadr4)
COMMON_BLOCK_DEF(Hadr4Common,hadr4);

/*================================================*/
/*common /psar5/  qcdlam,q2min,q2ini,q2fin,pt2cut,*/
/**betpom,glusea,naflav,alfe                      */
/*------------------------------------------------*/
typedef struct {
	Float_t qcdlam,q2min,q2ini,q2fin,pt2cut,betpom,glusea;
	Int_t naflav;
	Float_t alfe;
} Psar5Common;

#define psar5 COMMON_BLOCK(psar5,psar5)
COMMON_BLOCK_DEF(Psar5Common,psar5);

/*===========================*/
/*common/hadr17/edmaxi,epmaxi*/
/*---------------------------*/
typedef struct {
	Float_t edmaxi,epmaxi;
} Hadr17Common;

#define hadr17 COMMON_BLOCK(hadr17,hadr17)
COMMON_BLOCK_DEF(Hadr17Common,hadr17);

/*======================================================*/
/*common/frag1/ndecay,maxres,pud,pudk,pudr,strcut,diqcut*/
/*------------------------------------------------------*/
typedef struct {
	Int_t ndecay, maxres;
	Float_t pud,pudk,pudr,strcut,diqcut;
} Frag1Common;

#define frag1 COMMON_BLOCK(frag1,frag1)
COMMON_BLOCK_DEF(Frag1Common,frag1);

/*==================================*/
/*double precision seedi,seedj,seedc*/
/*common/cseed/seedi,seedj,seedc    */
/*----------------------------------*/
typedef struct {
	Double_t seedi,seedj,seedc;
} CseedCommon;

#define cseed COMMON_BLOCK(cseed,cseed)
COMMON_BLOCK_DEF(CseedCommon,cseed);

}
#endif /* EPOSCOMMON_H_ */
