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
  void  SetMinTracks(Int_t n=2) { fMinTracks = n; return; }
  void  SetSkipTracks(Int_t n,Int_t *skipped);
  void SetDebug(Int_t optdebug=0) {fDebug=optdebug;}
  void  SetVtxStart(Double_t x=0,Double_t y=0) 
    { fNominalPos[0]=x; fNominalPos[1]=y; return; }
  void  SetDCAcut(Double_t maxdca)
    { fDCAcut=maxdca; return;}
  void SetFinderAlgorithm(Int_t opt=1) 
    { fAlgo=opt; return;}
  
 protected:
  Double_t   GetField() const { return AliTracker::GetBz();} 
  void     ComputeMaxChi2PerTrack(Int_t nTracks);
  Int_t PrepareTracks(TTree &trkTree, Int_t OptImpParCut);
  void     VertexFinder(Int_t optUseWeights=0);
  void     HelixVertexFinder();
  void     StrLinVertexFinderMinDist(Int_t OptUseWeights=0);
  static void GetStrLinDerivMatrix(Double_t *p0,Double_t *p1,Double_t m[][3],Double_t *d);
  static void GetStrLinDerivMatrix(Double_t *p0,Double_t *p1,Double_t *sigmasq,Double_t m[][3],Double_t *d);
  static Double_t GetStrLinMinDist(Double_t *p0,Double_t *p1,Double_t *x0);
  static Double_t GetDeterminant3X3(Double_t matr[][3]);
  void     VertexFitter();

  AliVertex fVert;         // vertex after vertex finder
  AliESDVertex *fCurrentVertex;  // ESD vertex after fitter
  Double_t  fNominalPos[2];   // initial knowledge on vertex position
  Int_t     fMinTracks;       // minimum number of tracks
  Double_t  fMaxChi2PerTrack; // maximum contribition to the chi2 
  TObjArray fTrkArray;        // array with tracks to be processed
  Int_t     *fTrksToSkip;     // tracks to be skipped for find and fit 
  Int_t     fNTrksToSkip;     // number of tracks to be skipped 
  Double_t  fDCAcut;          // maximum DCA between 2 tracks used for vertex
  Int_t     fAlgo;            // option for vertex finding algorythm
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


  ClassDef(AliVertexerTracks,3) // 3D Vertexing with ESD tracks 
};

#endif



