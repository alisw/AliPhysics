#ifndef ALIITSALIGNMENTTRACK_H
#define ALIITSALIGNMENTTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */
/* $Author$ */
/* $Date$ */
/* $Name$ */
/* $Header$ */
/*
   $Log$
   Revision 1.1.2.4  2000/06/04 16:35:24  nilsen
   One more try to fix log comments.

   Revision 1.1.2.3  2000/03/12 16:05:55  nilsen
   Fixed but in $Log$
   Fixed but in Revision 1.1.2.4  2000/06/04 16:35:24  nilsen
   Fixed but in One more try to fix log comments.
   Fixed but in, hopefully.
  
   Revision 1.1.2.2  2000/03/04 23:40:19  nilsen
   Fixed the logs???
  
   Revision 1.1.2.1  2000/03/02 20:14:25  nilsen
   A new class useful for ITS detector alignment studdies. 
*/
/*   $Revision$ */

// Standard C & C++ libraries

// Standard Root Libraries
#include <TObject.h>
#include <TParticle.h>

// ITS libraries
#include "AliITSgeom.h"

struct TClustAl_sl{
    Int_t   findex;
    Float_t fxg[3],fExg[3][3];
    Float_t fxl[3],fExl[3][3];
};



class AliITSAlignmentTrack : public TObject{
///////////////////////////////////////////////////////////////////////////
//      A track class define exclusively for the use in doing ITS detector
// alignment studdies. Not intended for general use.
// Author: B. S. Nilsen
// Date:   January 17 2000
///////////////////////////////////////////////////////////////////////////

 protected:

    Int_t        ftrack,fnclust,fnclustMax;
    TClustAl_sl  *fclust;
    Int_t        ffunc;
    Double_t     fpar[10];
    Float_t      fChi2;
    Float_t      fpx,fpy,fpz,fp,fpt;
    // x=fp[0]+fp[1]*z and y=fp[2]+fp[3]*z                 : ffunc = 0
    //
    // x=fp[0]+fp[1]*y and z=fp[2]+fp[3]*y                 : ffunc = 1
    //
    // x=fp[0]+fp[2]*cos(th), y=fp[1]+fp[2]*sin(th),
    // and z=fp[3]+f[4]*th   th[i]=atan2(y[i]-y0,x[i]-x0)  : ffunc = 2

 public:

    AliITSAlignmentTrack();
    virtual ~AliITSAlignmentTrack();

    void CreatePoints(Int_t n){fclust = new TClustAl_sl[n];fnclustMax=n;
                  fnclust=0;for(Int_t j=0;j<fnclustMax;j++)
	          for(Int_t i=0;i<3;i++)fclust[j].fxl[i]=fclust[j].fxg[i]=0.0;}
    Int_t AddPointL(Int_t n,Int_t indx,Float_t *lp,Float_t **lep){
	if(n>=fnclustMax)return 0;fclust[n].findex=indx;
	for(Int_t i=0;i<3;i++){fclust[n].fxl[i]=lp[i];for(Int_t j=0;j<3;j++)
	    fclust[n].fExl[i][j]=lep[i][j];}return -1;}
    Int_t AddPointLastL(Int_t indx,Float_t *lp,Float_t **lep){
	if(++fnclust>=fnclustMax)return 0;//+first
	fclust[fnclust].findex=indx;
	for(Int_t i=0;i<3;i++){fclust[fnclust].fxl[i]=lp[i];
	   for(Int_t j=0;j<3;j++)fclust[fnclust].fExl[i][j]=lep[i][j];
	}return -1;}
    Int_t GetPointL(Int_t n,Float_t *lp){if(n>=fnclustMax)return 0;
                       for(Int_t i=0;i<3;i++)
			   lp[i]=fclust[n].fxl[i];return -1;}
    Int_t GetPointL(Int_t n,Double_t *lp){if(n>=fnclustMax)return 0;
                       for(Int_t i=0;i<3;i++)
			   lp[i]=(Double_t)fclust[n].fxl[i];return -1;}
    Int_t GetPointG(Int_t n,Float_t *gp){if(n>=fnclustMax)return 0;
                       for(Int_t i=0;i<3;i++)
			   gp[i]=fclust[n].fxg[i];return -1;}
    Int_t GetPointG(Int_t n,Double_t *gp){if(n>=fnclustMax)return 0;
                       for(Int_t i=0;i<3;i++)
			   gp[i]=(Double_t)fclust[n].fxg[i];return -1;}
    Int_t GetErrorG(Int_t n,Float_t **gp){if(n>=fnclustMax)return 0;
                       for(Int_t i=0;i<3;i++)for(Int_t j=0;j<3;j++)
			   gp[i][j]=fclust[n].fExg[i][j];return -1;}
    Int_t GetErrorG(Int_t n,Double_t **gp){if(n>=fnclustMax)return 0;
                       for(Int_t i=0;i<3;i++)for(Int_t j=0;j<3;j++)
			   gp[i][j]=(Double_t)fclust[n].fExg[i][j];return -1;}
    Int_t GetErrorL(Int_t n,Float_t **gp){if(n>=fnclustMax)return 0;
                       for(Int_t i=0;i<3;i++)for(Int_t j=0;j<3;j++)
			   gp[i][j]=fclust[n].fExl[i][j];return -1;}
    Int_t GetErrorL(Int_t n,Double_t **gp){if(n>=fnclustMax)return 0;
                       for(Int_t i=0;i<3;i++)for(Int_t j=0;j<3;j++)
			   gp[i][j]=(Double_t)fclust[n].fExl[i][j];return -1;}
    Int_t GetIndex(Int_t n,Int_t &indx){if(n>=fnclustMax)return 0;
                                              indx=fclust[n].findex;return -1;}

    void SetGlobalPosition(AliITSgeom *gm){
	for(Int_t i=0;i<fnclust;i++){
	    gm->LtoG(fclust[i].findex,fclust[i].fxl,fclust[i].fxg);
         gm->LtoGErrorMatrix(fclust[i].findex,(Double_t **) fclust[i].fExl,
			     (Double_t **) fclust[i].fExg);}}
    void SetTrackNumber(Int_t trk) {ftrack=trk;}
    void SetTParticle(TParticle *prt){fpx=prt->Px();fpy=prt->Py();
                                      fpz=prt->Pz();fp=prt->P();fpt=prt->Pt();}
    void SetPx(Float_t px) {fpx=px;}
    void SetPy(Float_t py) {fpy=py;}
    void SetPz(Float_t pz) {fpz=pz;}
    void SetP(Float_t p) {fp=p;}
    void SetPt(Float_t pt) {fpt=pt;}
    void SetParameter(Int_t n,Double_t a) {if(n>=0&&n<10);fpar[n]=a;}
    void SetChi2(Float_t chi2) {fChi2=chi2;}
    void SetFunctionNumber(Int_t fn) {ffunc = fn;}

    Int_t GetTrackNumber() {return ftrack;}
    Int_t GetNumberOfClustersSl() {return fnclust;}
    Float_t GetPx() {return fpx;}
    Float_t GetPy() {return fpy;}
    Float_t GetPz() {return fpz;}
    Float_t GetP() {return fp;}
    Float_t GetPt() {return fpt;}
    Double_t GetParameter(Int_t n) {if(n>=0&&n<10) return fpar[n]; return -1.;}
    Float_t GetChi2() {return fChi2;}
    Int_t   GetfunID() {return ffunc;}

    void FitToFunction(Int_t n,AliITSgeom *gm);
    void func(Double_t *go,Double_t *gi);
    void func0(Double_t *go,Double_t *gi);
    void func1(Double_t *go,Double_t *gi);
    void func2(Double_t *go,Double_t *gi);
    Int_t FitTrackToLineG();
    Int_t FitTrackToLineL(AliITSgeom *gm);
    Int_t FindCircleCenter(Double_t *xc,Double_t *x1,
			   Double_t *x2,Double_t *x3);

    private:
    Double_t ComputeChi2();

    ClassDef(AliITSAlignmentTrack,1) // Track class for ITS Alignment

};
#endif
