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
#include "TMath.h"

class AliESDEvent;
class AliESDv0;
class AliExternalTrackParam;

//_____________________________________________________________________________
class AliCascadeVertexer : public TObject {
public:
  AliCascadeVertexer();
  void SetCuts(const Double_t cuts[8]);
  static void SetDefaultCuts(const Double_t cuts[8]);

  Int_t V0sTracks2CascadeVertices(AliESDEvent *event);
  Double_t Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const;
  Double_t Det(Double_t a00,Double_t a01,Double_t a02,
	       Double_t a10,Double_t a11,Double_t a12,
	       Double_t a20,Double_t a21,Double_t a22) const;

  Double_t PropagateToDCA(AliESDv0 *vtx,AliExternalTrackParam *trk,Double_t b);

  void GetCuts(Double_t cuts[8]) const;
  static void GetDefaultCuts(Double_t cuts[8]);

private:
  static
  Double_t fgChi2max;   // maximal allowed chi2 
  static
  Double_t fgDV0min;    // min. allowed V0 impact parameter
  static
  Double_t fgMassWin;   // window around the Lambda mass
  static
  Double_t fgDBachMin;  // min. allowed bachelor impact parameter
  static
  Double_t fgDCAmax;    // maximal allowed DCA between the V0 and the track 
  static
  Double_t fgCPAmin;    // minimal allowed cosine of the cascade pointing angle
  static
  Double_t fgRmin, fgRmax;// max & min radii of the fiducial volume
  
  Double_t fChi2max;    // maximal allowed chi2 
  Double_t fDV0min;     // min. allowed V0 impact parameter
  Double_t fMassWin;    // window around the Lambda mass
  Double_t fDBachMin;   // min. allowed bachelor impact parameter
  Double_t fDCAmax;     // maximal allowed DCA between the V0 and the track 
  Double_t fCPAmin;     // minimal allowed cosine of the cascade pointing angle
  Double_t fRmin2, fRmax2;// max & min radii^2 of the fiducial volume
  Double_t fRmaxMargin2;// (rmax+margin)^2  for fast checks
  
  ClassDef(AliCascadeVertexer,4)  // cascade verterxer 
};

inline AliCascadeVertexer::AliCascadeVertexer() 
:TObject(),
  fChi2max(fgChi2max), 
  fDV0min(fgDV0min),
  fMassWin(fgMassWin),
  fDBachMin(fgDBachMin),
  fDCAmax(fgDCAmax),
  fCPAmin(fgCPAmin), 
  fRmin2(fgRmin*fgRmin),
  fRmax2(fgRmax*fgRmax),
  fRmaxMargin2((fgRmax+5)*(fgRmax+5))
{
}

inline void AliCascadeVertexer::SetCuts(const Double_t cuts[8]) {
  fChi2max=cuts[0]; 
  fDV0min=cuts[1];   fMassWin=cuts[2]; fDBachMin=cuts[3];
  fDCAmax=cuts[4];   fCPAmin=cuts[5];
  fRmin2=cuts[6]*cuts[6]; 
  fRmax2=cuts[7]*cuts[7]; 
  fRmaxMargin2 = (cuts[7]+5)*(cuts[7]+5);
}

inline void AliCascadeVertexer::SetDefaultCuts(const Double_t cuts[8]) {
  fgChi2max=cuts[0]; 
  fgDV0min=cuts[1];   fgMassWin=cuts[2]; fgDBachMin=cuts[3];
  fgDCAmax=cuts[4];   fgCPAmin=cuts[5];
  fgRmin=cuts[6];     fgRmax=cuts[7]; 
}

inline void AliCascadeVertexer::GetCuts(Double_t cuts[8]) const {
  cuts[0]=fChi2max; 
  cuts[1]=fDV0min;   cuts[2]=fMassWin;  cuts[3]=fDBachMin;
  cuts[4]=fDCAmax;   cuts[5]=fCPAmin;
  cuts[6]=TMath::Sqrt(fRmin2);     cuts[7]=TMath::Sqrt(fRmax2); 
}

inline void AliCascadeVertexer::GetDefaultCuts(Double_t cuts[8]) {
  cuts[0]=fgChi2max; 
  cuts[1]=fgDV0min;   cuts[2]=fgMassWin;  cuts[3]=fgDBachMin;
  cuts[4]=fgDCAmax;   cuts[5]=fgCPAmin;
  cuts[6]=fgRmin;     cuts[7]=fgRmax; 
}

#endif

