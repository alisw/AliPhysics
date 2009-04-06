#ifndef ALIVERTEXERTRACKS_H
#define ALIVERTEXERTRACKS_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------
// Class for vertex determination with tracks
//
//   Origin: AliITSVertexerTracks  
//           A.Dainese, Padova, andrea.dainese@pd.infn.it
//           M.Masera,  Torino, massimo.masera@to.infn.it 
//   Moved to STEER and adapted to ESD tracks: 
//           F.Prino, Torino, prino@to.infn.it 
//-------------------------------------------------------

/*****************************************************************************
 *                                                                           *
 * This class determines the vertex of a set of tracks.                      *
 * Different algorithms are implemented, see data member fAlgo.              *
 *                                                                           *
 *****************************************************************************/

#include <TObjArray.h>
#include <TMatrixD.h>

#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"

class AliVEvent;
class AliESDEvent;
class AliStrLine;

class AliVertexerTracks : public TObject {
  
 public:
  AliVertexerTracks(); 
  AliVertexerTracks(Double_t fieldkG); 
  virtual ~AliVertexerTracks();

  AliESDVertex* FindPrimaryVertex(AliVEvent *vEvent);
  AliESDVertex* FindPrimaryVertex(TObjArray *trkArrayOrig,UShort_t *idOrig);
  AliESDVertex* VertexForSelectedTracks(TObjArray *trkArray,UShort_t *id,
					Bool_t optUseFitter=kTRUE,
					Bool_t optPropagate=kTRUE);
  AliESDVertex* VertexForSelectedESDTracks(TObjArray *trkArray,
					Bool_t optUseFitter=kTRUE,
					Bool_t optPropagate=kTRUE);
  AliESDVertex* RemoveTracksFromVertex(AliESDVertex *inVtx,
				       TObjArray *trkArray,UShort_t *id,
				       Float_t *diamondxy); 
  void  SetITSMode(Double_t dcacut=0.1,
		   Double_t dcacutIter0=0.1,
		   Double_t maxd0z0=0.5,
		   Int_t minCls=5,
		   Int_t mintrks=1,
		   Double_t nsigma=3.,
		   Double_t mindetfitter=100.,
		   Double_t maxtgl=1000.,
		   Double_t fidR=3.,
		   Double_t fidZ=30.,
		   Int_t finderAlgo=1,
		   Int_t finderAlgoIter0=4); 
  void  SetTPCMode(Double_t dcacut=0.1,
		   Double_t dcacutIter0=1.0,
		   Double_t maxd0z0=5.0,
		   Int_t minCls=10,
		   Int_t mintrks=1,
		   Double_t nsigma=3.,
		   Double_t mindetfitter=0.1,
		   Double_t maxtgl=1.5, 
		   Double_t fidR=3.,
		   Double_t fidZ=30.,
		   Int_t finderAlgo=1,
		   Int_t finderAlgoIter0=4); 
  void  SetCuts(Double_t *cuts);
  void  SetConstraintOff() { fConstraint=kFALSE; return; }
  void  SetConstraintOn() { fConstraint=kTRUE; return; }
  void  SetDCAcut(Double_t maxdca) { fDCAcut=maxdca; return; }
  void  SetDCAcutIter0(Double_t maxdca) { fDCAcutIter0=maxdca; return; }
  void  SetFinderAlgorithm(Int_t opt=1) { fAlgo=opt; return; }
  void  SetITSrefitRequired() { fITSrefit=kTRUE; return; }
  Bool_t GetITSrefitRequired() const { return fITSrefit; }
  void  SetITSrefitNotRequired() { fITSrefit=kFALSE; return; }
  void  SetFiducialRZ(Double_t r=3,Double_t z=30) { fFiducialR=r; fFiducialZ=z; return; }
  void  SetMaxd0z0(Double_t maxd0z0=0.5) { fMaxd0z0=maxd0z0; return; }
  void  SetMinClusters(Int_t n=5) { fMinClusters=n; return; }
  Int_t GetMinClusters() const { return fMinClusters; }
  void  SetMinTracks(Int_t n=1) { fMinTracks=n; return; }
  void  SetNSigmad0(Double_t n=3) { fNSigma=n; return; }
  Double_t GetNSigmad0() const { return fNSigma; }
  void  SetMinDetFitter(Double_t mindet=100.) { fMinDetFitter=mindet; return; }
  void  SetMaxTgl(Double_t maxtgl=1.) { fMaxTgl=maxtgl; return; }
  void  SetOnlyFitter() { if(!fConstraint) AliFatal("Set constraint first!"); 
     fOnlyFitter=kTRUE; return; }
  void  SetSkipTracks(Int_t n,Int_t *skipped);
  void  SetVtxStart(Double_t x=0,Double_t y=0,Double_t z=0) 
    { fNominalPos[0]=x; fNominalPos[1]=y; fNominalPos[2]=z; return; }
  void  SetVtxStartSigma(Double_t sx=3.,Double_t sy=3.,Double_t sz=15.) 
    { fNominalCov[0]=sx*sx; fNominalCov[2]=sy*sy; fNominalCov[5]=sz*sz;
      fNominalCov[1]=0.; fNominalCov[3]=0.; fNominalCov[4]=0.; return; }
  void  SetVtxStart(AliESDVertex *vtx);
  static Double_t GetStrLinMinDist(Double_t *p0,Double_t *p1,Double_t *x0);
  static Double_t GetDeterminant3X3(Double_t matr[][3]);
  static void GetStrLinDerivMatrix(Double_t *p0,Double_t *p1,Double_t (*m)[3],Double_t *d);
  static void GetStrLinDerivMatrix(Double_t *p0,Double_t *p1,Double_t *sigmasq,Double_t (*m)[3],Double_t *d);
  static AliESDVertex TrackletVertexFinder(TClonesArray *lines, Int_t optUseWeights=0);
  static AliESDVertex TrackletVertexFinder(AliStrLine **lines, const Int_t knacc, Int_t optUseWeights=0);
  void     SetFieldkG(Double_t field=-999.) { fFieldkG=field; return; }
  Double_t GetFieldkG() const { 
    if(fFieldkG<-99.) AliFatal("Field value not set");
    return fFieldkG; } 
  void SetNSigmaForUi00(Double_t n=1.5) { fnSigmaForUi00=n; return; }
  Double_t GetNSigmaForUi00() const { return fnSigmaForUi00; }

 protected:
  void     HelixVertexFinder();
  void     OneTrackVertFinder();
  Int_t    PrepareTracks(TObjArray &trkArrayOrig,UShort_t *idOrig,
			 Int_t optImpParCut);
  Bool_t   PropagateTrackTo(AliExternalTrackParam *track,
			    Double_t xToGo);
  Bool_t   TrackToPoint(AliExternalTrackParam *t,
		        TMatrixD &ri,TMatrixD &wWi,
			Bool_t uUi3by3=kFALSE) const;     
  void     VertexFinder(Int_t optUseWeights=0);
  void     VertexFitter();
  void     StrLinVertexFinderMinDist(Int_t optUseWeights=0);
  void     TooFewTracks();

  AliESDVertex fVert;         // vertex after vertex finder
  AliESDVertex *fCurrentVertex;  // ESD vertex after fitter
  UShort_t  fMode;            // 0 ITS+TPC; 1 TPC
  Double_t  fFieldkG;         // z component of field (kGauss) 
  Double_t  fNominalPos[3];   // initial knowledge on vertex position
  Double_t  fNominalCov[6];   // initial knowledge on vertex position
  TObjArray fTrkArraySel;     // array with tracks to be processed
  UShort_t  *fIdSel;          // IDs of the tracks (AliESDtrack::GetID())
  Int_t     *fTrksToSkip;     // track IDs to be skipped for find and fit 
  Int_t     fNTrksToSkip;     // number of tracks to be skipped 
  Bool_t    fConstraint;      // true when "mean vertex" was set in 
                              // fNominal ... and must be used in the fit
  Bool_t    fOnlyFitter;      // primary with one fitter shot only
                              // (use only with beam constraint)
  Int_t     fMinTracks;       // minimum number of tracks
  Int_t     fMinClusters;     // minimum number of ITS or TPC clusters per track
  Double_t  fDCAcut;          // maximum DCA between 2 tracks used for vertex
  Double_t  fDCAcutIter0;     // maximum DCA between 2 tracks used for vertex
  Double_t  fNSigma;          // number of sigmas for d0 cut in PrepareTracks()
  Double_t  fMaxd0z0;         // value for sqrt(d0d0+z0z0) cut 
                              // in PrepareTracks(1) if fConstraint=kFALSE
  Double_t  fMinDetFitter;    // minimum determinant to try to invertex matrix
  Double_t  fMaxTgl;          // maximum tgl of tracks
  Bool_t    fITSrefit;        // if kTRUE (default), use only kITSrefit tracks
                              // if kFALSE, use all tracks (also TPC only)
  Double_t  fFiducialR;       // radius of fiducial cylinder for tracks 
  Double_t  fFiducialZ;       // length of fiducial cylinder for tracks
  Double_t  fnSigmaForUi00;   // n. sigmas from finder in TrackToPoint
  Int_t     fAlgo;            // option for vertex finding algorythm
  Int_t     fAlgoIter0;       // this is for iteration 0
  // fAlgo=1 (default) finds minimum-distance point among all selected tracks
  //         approximated as straight lines 
  //         and uses errors on track parameters as weights
  // fAlgo=2 finds minimum-distance point among all the selected tracks
  //         approximated as straight lines 
  // fAlgo=3 finds the average point among DCA points of all pairs of tracks
  //         treated as helices
  // fAlgo=4 finds the average point among DCA points of all pairs of tracks
  //         approximated as straight lines 
  //         and uses errors on track parameters as weights
  // fAlgo=5 finds the average point among DCA points of all pairs of tracks
  //         approximated as straight lines 


 private:
  AliVertexerTracks(const AliVertexerTracks & source);
  AliVertexerTracks & operator=(const AliVertexerTracks & source);

  ClassDef(AliVertexerTracks,12) // 3D Vertexing with tracks 
};

#endif

