#ifndef ALIVERTEXERTRACKS_H
#define ALIVERTEXERTRACKS_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------
// Class for vertex determination with ESD tracks
//
//   Origin: AliITSVertexerTracks  
//           A.Dainese, Padova, andrea.dainese@pd.infn.it
//           M.Masera,  Torino, massimo.masera@to.infn.it 
//   Moved to STEER and adapted to ESD tracks: 
//           F.Prino, Torino, prino@to.infn.it 
//-------------------------------------------------------

/*****************************************************************************
 *                                                                           *
 * This class determines the vertex of a set of ESD tracks.                  *
 * Different algorithms are implemented, see data member fAlgo.              *
 *                                                                           *
 *****************************************************************************/

#include "AliESDVertex.h"
#include "AliTracker.h"
#include "AliLog.h"

#include <TObjArray.h>

class TTree; 
class AliESD;

class AliVertexerTracks : public TObject {
  
 public:
  AliVertexerTracks(); 
  AliVertexerTracks(Double_t xStart, Double_t yStart); 
  virtual ~AliVertexerTracks();


  // computes the vertex from the set of tracks in the tree
  AliVertex* VertexForSelectedTracks(TTree *trkTree);
  AliVertex* VertexForSelectedTracks(TObjArray *trkArray);
  AliESDVertex* FindPrimaryVertex(const AliESD *esdEvent);
  AliESDVertex* FindPrimaryVertexOld(const AliESD *esdEvent);
  void  SetMinTracks(Int_t n=2) { fMinTracks = n; return; }
  void  SetMinITSClusters(Int_t n=5) { fMinITSClusters = n; return; }
  void  SetSkipTracks(Int_t n,Int_t *skipped);
  void SetDebug(Int_t optdebug=0) {fDebug=optdebug;}
  void  SetVtxStart(Double_t x=0,Double_t y=0,Double_t z=0) 
    { fNominalPos[0]=x; fNominalPos[1]=y; fNominalPos[2]=z; return; }
  void  SetVtxStartSigma(Double_t sx=3,Double_t sy=3,Double_t sz=6) 
    { fNominalSigma[0]=sx; fNominalSigma[1]=sy; fNominalSigma[2]=sz; return; }
  void  SetVtxStart(AliESDVertex *vtx) 
    { SetVtxStart(vtx->GetXv(),vtx->GetYv(),vtx->GetZv());
      SetVtxStartSigma(vtx->GetXRes(),vtx->GetYRes(),vtx->GetZRes()); return; }
  void  SetDCAcut(Double_t maxdca)
    { fDCAcut=maxdca; return;}
  void SetFinderAlgorithm(Int_t opt=1) 
    { fAlgo=opt; return;}
  void  SetNSigmad0(Double_t n=3) 
    { fNSigma=n; return; }
  static Double_t GetStrLinMinDist(Double_t *p0,Double_t *p1,Double_t *x0);
  static Double_t GetDeterminant3X3(Double_t matr[][3]);
  static void GetStrLinDerivMatrix(Double_t *p0,Double_t *p1,Double_t (*m)[3],Double_t *d);
  static void GetStrLinDerivMatrix(Double_t *p0,Double_t *p1,Double_t *sigmasq,Double_t (*m)[3],Double_t *d);

 protected:
  Double_t GetField() const { 
    if(!AliTracker::GetFieldMap())
      AliFatal("Field map not set; use AliTracker::SetFieldMap()!");
    return AliTracker::GetBz(); } 
  Int_t    PrepareTracks(TTree &trkTree, Int_t OptImpParCut);
  Double_t Sigmad0rphi(Double_t pt) const;
  void     VertexFinder(Int_t optUseWeights=0);
  void     HelixVertexFinder();
  void     StrLinVertexFinderMinDist(Int_t OptUseWeights=0);
  void     VertexFitter(Bool_t useNominaVtx=kFALSE);
  void     TooFewTracks(const AliESD *esdEvent);

   
  AliVertex fVert;         // vertex after vertex finder
  AliESDVertex *fCurrentVertex;  // ESD vertex after fitter
  Double_t  fNominalPos[3];   // initial knowledge on vertex position
  Double_t  fNominalSigma[3]; // initial knowledge on vertex position
  Int_t     fMinTracks;       // minimum number of tracks
  Int_t     fMinITSClusters;  // minimum number of ITS clusters per track
  TObjArray fTrkArray;        // array with tracks to be processed
  Int_t     *fTrksToSkip;     // tracks to be skipped for find and fit 
  Int_t     fNTrksToSkip;     // number of tracks to be skipped 
  Double_t  fDCAcut;          // maximum DCA between 2 tracks used for vertex
  Int_t     fAlgo;            // option for vertex finding algorythm
  Double_t  fNSigma;          // number of sigmas for d0 cut in PrepareTracks()
  Int_t fDebug;               //! debug flag - verbose printing if >0
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

  ClassDef(AliVertexerTracks,4) // 3D Vertexing with ESD tracks 
};

#endif
