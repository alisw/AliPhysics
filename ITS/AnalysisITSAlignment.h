#ifndef ALIGNITS_H
#define ALIGNITS_H
//////////////////////////////////////////////////////////////////////////
//  Alice ITS first detector alignment program.                         //
//                                                                      //
// version: 0.0.0 Draft.                                                //
// Date: April 18 1999                                                  //
// By: Bjorn S. Nilsen                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "AliITS.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "AliRun.h"
#include "AliITSAlignmentTrack.h"

// Data structure to hold averaged clusts.
struct ClustAl_sl{
    Int_t lay,lad,det;
    Float_t xg,yg,zg,xl,yl,zl;
};
struct ClustAl_tl{
    Int_t    track,nclust;  // track number and number of data points.
    ClustAl_sl *clust;        // data points to fit.
    Float_t  a,b,c,d,a0,b0,c0,d0,qual;  // fit parameters and fit quality.
    Float_t  px,py,pz,p,pt;
    // x=a+b*z and y=c+d*z;
    // x=a0+b0*y and z=c0+d0*y in coordinate system of clust[0].lay,lad,det
};

void HitsTo(ClustAl_tl *trk,Int_t &ntrk,Int_t,TTree *TH,AliITS *ITS,
	    Float_t nsigmaT1,Float_t nsigmaT2,Float_t nsigmaT3,
	    Float_t nsigmaR1,Float_t nsigmaR2,Float_t nsigmaR3);
void HitsToClustAl(ClustAl_tl *trk,Int_t &ntrk,Int_t nt,TTree *TH,
		   AliITS *ITS,Float_t fraction);
void PlotGeomChanges(AliITSgeom *gt,AliITSgeom *gc,TFile *Hfile,Float_t *Rdta);
void FillGlobalPositions(ClustAl_tl *trk,Int_t ntrk,AliITSgeom *g);
void FitAllTracks(ClustAl_tl *trk,Int_t ntrk,Float_t *v0,AliITSgeom *gm,
		  const char *sfile,TFile *Hfile,Float_t *Fdta,Int_t *Ndta);
void FitVertexAll(ClustAl_tl *trk,Int_t ntrk,const char *sfile,TFile *Hfile);
void OnlyOneGeometry(char *filename,AliITSgeom *gm,AliITSgeom &gm2,
		     Float_t trans[],Float_t rot[]);
void deleteClustAl(ClustAl_tl *trk,Int_t ntrk);
void FillAliITSAlignmentTrack(AliITSAlignmentTrack *trk,Int_t &ntrk,Int_t nt,
			      TTree *TH,AliITS *ITS,Float_t fraction);
void RunAlignment(Int_t evnt,Float_t fraction);

#endif
