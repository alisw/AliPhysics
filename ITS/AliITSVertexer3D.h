#ifndef ALIITSVERTEXER3D_H
#define ALIITSVERTEXER3D_H

#include<AliITSVertexer.h>

///////////////////////////////////////////////////////////////////
//                                                               //
// Class for primary vertex finding  (3D reconstruction)         //
//                                                               //
///////////////////////////////////////////////////////////////////

/* $Id$ */

#include <AliESDVertex.h>

class TClonesArray;
class AliITSVertexer3D : public AliITSVertexer {

 public:

  AliITSVertexer3D();
  AliITSVertexer3D(TString filename);
  virtual ~AliITSVertexer3D();
  virtual AliESDVertex* FindVertexForCurrentEvent(Int_t evnumb);
  virtual void FindVertices();
  AliESDVertex GetVertex3D() const {return fVert3D;}
//   virtual void MakeTracklet(Double_t *pA, Double_t *pB, Int_t &nolines); */
//   virtual void MakeTracklet(Float_t *pA, Float_t *pB, Int_t &nolines); */
  virtual void PrintStatus() const;
  void SetCoarseDiffPhiCut(Float_t dphi = 0.5){fCoarseDiffPhiCut=dphi;}
  void SetCoarseMaxRCut(Float_t rad = 2.5){fCoarseMaxRCut=rad;}
  void SetMaxRCut(Float_t rad = 0.5){fMaxRCut=rad;}
  void SetZCutDiamond(Float_t zcut = 14.0){fZCutDiamond=zcut;}
  void SetMaxZCut(Float_t dz = 0.5){fMaxZCut=dz;}
  void SetDCAcut(Float_t dca=0.1){fDCAcut=dca;} 
  void SetDiffPhiMax(Float_t pm = 0.01){fDiffPhiMax = pm;}
  void SetMeanPSelTracks(Float_t pGeV=0.875){fMeanPSelTrk = pGeV;}
  void SetMeanPtSelTracks(Float_t ptGeV=0.630){fMeanPtSelTrk = ptGeV;}
  void SetMeanPPtSelTracks(Float_t fieldTesla);

protected:
  AliITSVertexer3D(const AliITSVertexer3D& vtxr);
  AliITSVertexer3D& operator=(const AliITSVertexer3D& /* vtxr */);
  Int_t FindTracklets(Int_t evnumber, Int_t optCuts);
  Int_t Prepare3DVertex(Int_t optCuts);
  void ResetVert3D();


  TClonesArray *fLines;      //! array of tracklets
  AliESDVertex fVert3D;        // 3D Vertex
  Float_t fCoarseDiffPhiCut; // loose cut on DeltaPhi for RecPoint matching 
  Float_t fCoarseMaxRCut; // cut on tracklet DCA to Z axis
  Float_t fMaxRCut; // cut on tracklet DCA to beam axis
  Float_t fZCutDiamond;   // cut on +-Z of the diamond
  Float_t fMaxZCut;   // cut on Z distance from estimated vertex
  Float_t fDCAcut; // cut on tracklet to tracklet and tracklet to vertex DCA
  Float_t fDiffPhiMax;     // Maximum delta phi allowed among corr. pixels
  Float_t fMeanPSelTrk; // GeV, mean P for tracks with dphi<0.01 rad
  Float_t fMeanPtSelTrk; // GeV, mean Pt for tracks with dphi<0.01 rad

  ClassDef(AliITSVertexer3D,3);

};

#endif
