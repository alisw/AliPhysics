#ifndef ALICASCADEVERTEXER_H
#define ALICASCADEVERTEXER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//------------------------------------------------------------------
//                    Cascade Vertexer Class
//          Reads V0s and tracks, writes out cascade vertices
//    Origin: Christian Kuhn, IReS, Strasbourg, christian.kuhn@ires.in2p3.fr
//------------------------------------------------------------------

#include "TObject.h"

class AliESD;
class TTree;
class AliESDv0;
class AliExternalTrackParam;

//_____________________________________________________________________________
class AliCascadeVertexer : public TObject {
public:
  AliCascadeVertexer();
  AliCascadeVertexer(const Double_t cuts[8]);
  void SetCuts(const Double_t cuts[8]);
  void SetVertex(Double_t *vtx) { fX=vtx[0]; fY=vtx[1]; fZ=vtx[2]; }

  Int_t V0sTracks2CascadeVertices(AliESD *event);
  Double_t PropagateToDCA(AliESDv0 *vtx,AliExternalTrackParam *trk,Double_t b);

  void GetCuts(Double_t cuts[8]) const;
  void GetVertex(Double_t *vtx) const { vtx[0]=fX; vtx[1]=fY; vtx[2]=fZ; }

private:
  Double_t fChi2max;    // maximal allowed chi2 
  Double_t fDV0min;     // min. allowed V0 impact parameter
  Double_t fMassWin;    // window around the Lambda mass
  Double_t fDBachMin;   // min. allowed bachelor impact parameter
  Double_t fDCAmax;     // maximal allowed DCA between the V0 and the track 
  Double_t fCPAmax;     // maximal allowed cosine of the cascade pointing angle
  Double_t fRmin, fRmax;// max & min radii of the fiducial volume
  
  Double_t fX;            // X-coordinate of the primary vertex
  Double_t fY;            // Y-coordinate of the primary vertex
  Double_t fZ;            // Z-coordinate of the primary vertex

  ClassDef(AliCascadeVertexer,1)  // cascade verterxer 
};

inline AliCascadeVertexer::AliCascadeVertexer() {
 fChi2max=33.; 
 fDV0min=0.015; fMassWin=0.05; fDBachMin=0.015;
 fDCAmax=0.01;  fCPAmax=0.025; 
 fRmin=0.5;     fRmax=2.5; 
}

inline AliCascadeVertexer::AliCascadeVertexer(const Double_t cuts[8]) {
  fChi2max=cuts[0]; 
  fDV0min=cuts[1];   fMassWin=cuts[2]; fDBachMin=cuts[3];
  fDCAmax=cuts[4];   fCPAmax=cuts[5];
  fRmin=cuts[6];     fRmax=cuts[7]; 
  fX=fY=fZ=0.;
}

inline void AliCascadeVertexer::SetCuts(const Double_t cuts[8]) {
  fChi2max=cuts[0]; 
  fDV0min=cuts[1];   fMassWin=cuts[2]; fDBachMin=cuts[3];
  fDCAmax=cuts[4];   fCPAmax=cuts[5];
  fRmin=cuts[6];     fRmax=cuts[7]; 
}

inline void AliCascadeVertexer::GetCuts(Double_t cuts[8]) const {
  cuts[0]=fChi2max; 
  cuts[1]=fDV0min;   cuts[2]=fMassWin;  cuts[3]=fDBachMin;
  cuts[4]=fDCAmax;   cuts[5]=fCPAmax;
  cuts[6]=fRmin;     cuts[7]=fRmax; 
}

#endif

