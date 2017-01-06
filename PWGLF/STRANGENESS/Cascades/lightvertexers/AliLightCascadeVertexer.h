#ifndef AliLightCascadeVertexer_H
#define AliLightCascadeVertexer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//------------------------------------------------------------------
//                    Cascade Vertexer Class
//          Reads V0s and tracks, writes out cascade vertices
//    Origin: Christian Kuhn, IReS, Strasbourg, christian.kuhn@ires.in2p3.fr
//------------------------------------------------------------------

#include "TObject.h"

class AliESDEvent;
class AliESDv0;
class AliExternalTrackParam;

//_____________________________________________________________________________
class AliLightCascadeVertexer : public TObject {
public:
  AliLightCascadeVertexer();
  void SetCuts(const Double_t cuts[8]);
  static void SetDefaultCuts(const Double_t cuts[8]);

  Int_t V0sTracks2CascadeVertices(AliESDEvent *event);
  Double_t Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const;
  Double_t Det(Double_t a00,Double_t a01,Double_t a02,
	       Double_t a10,Double_t a11,Double_t a12,
	       Double_t a20,Double_t a21,Double_t a22) const;

  Double_t PropagateToDCA(AliESDv0 *vtx,AliExternalTrackParam *trk,Double_t b);
    void CheckChargeV0(AliESDv0 *v0);

  void GetCuts(Double_t cuts[8]) const;
    static void GetDefaultCuts(Double_t cuts[8]);
    static void SetDefaultMaxEta(Double_t lMaxEta);
    static void SetDefaultMinClusters(Int_t lMinClusters);
    static void SetDefaultUseOnTheFlyV0 (Bool_t lOption);
    void SetMaxEta(Double_t lMaxEta);
    void SetMinClusters(Int_t lMinClusters);
    void SetSwitchCharges(Bool_t lOption);
    void SetUseOnTheFlyV0 (Bool_t lOption);
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
    static Double_t fgMaxEta;       // maximum eta value for track pre-selection
    static Int_t fgMinClusters;  // minimum single-track clusters value (>=)
  static Bool_t fgSwitchCharges;  // minimum single-track clusters value (>=)
  static Bool_t fgUseOnTheFlyV0;  // minimum single-track clusters value (>=)
    
  Double_t fChi2max;    // maximal allowed chi2 
  Double_t fDV0min;     // min. allowed V0 impact parameter
  Double_t fMassWin;    // window around the Lambda mass
  Double_t fDBachMin;   // min. allowed bachelor impact parameter
  Double_t fDCAmax;     // maximal allowed DCA between the V0 and the track 
  Double_t fCPAmin;     // minimal allowed cosine of the cascade pointing angle
  Double_t fRmin, fRmax;// max & min radii of the fiducial volume
    Double_t fMaxEta;       // maximum eta value for track pre-selection
    Int_t fMinClusters;  // minimum single-track clusters value (>=)
    Bool_t fSwitchCharges; //switch to change bachelor charge
    Bool_t fUseOnTheFlyV0; //switch to use on-the-fly V0s (HIGHLY EXPERIMENTAL)
  
  ClassDef(AliLightCascadeVertexer,3)  // cascade verterxer 
};

inline AliLightCascadeVertexer::AliLightCascadeVertexer() :
  TObject(),
  fChi2max(fgChi2max), 
  fDV0min(fgDV0min),
  fMassWin(fgMassWin),
  fDBachMin(fgDBachMin),
  fDCAmax(fgDCAmax),
  fCPAmin(fgCPAmin), 
  fRmin(fgRmin),
fRmax(fgRmax),
fMaxEta(fgMaxEta),
fMinClusters(fgMinClusters),
fSwitchCharges(fgSwitchCharges),
fUseOnTheFlyV0(fgUseOnTheFlyV0)
{
}

inline void AliLightCascadeVertexer::SetCuts(const Double_t cuts[8]) {
  fChi2max=cuts[0]; 
  fDV0min=cuts[1];   fMassWin=cuts[2]; fDBachMin=cuts[3];
  fDCAmax=cuts[4];   fCPAmin=cuts[5];
  fRmin=cuts[6];     fRmax=cuts[7]; 
}

inline void AliLightCascadeVertexer::SetDefaultCuts(const Double_t cuts[8]) {
  fgChi2max=cuts[0]; 
  fgDV0min=cuts[1];   fgMassWin=cuts[2]; fgDBachMin=cuts[3];
  fgDCAmax=cuts[4];   fgCPAmin=cuts[5];
  fgRmin=cuts[6];     fgRmax=cuts[7]; 
}

inline void AliLightCascadeVertexer::GetCuts(Double_t cuts[8]) const {
  cuts[0]=fChi2max; 
  cuts[1]=fDV0min;   cuts[2]=fMassWin;  cuts[3]=fDBachMin;
  cuts[4]=fDCAmax;   cuts[5]=fCPAmin;
  cuts[6]=fRmin;     cuts[7]=fRmax; 
}

inline void AliLightCascadeVertexer::GetDefaultCuts(Double_t cuts[8]) {
  cuts[0]=fgChi2max; 
  cuts[1]=fgDV0min;   cuts[2]=fgMassWin;  cuts[3]=fgDBachMin;
  cuts[4]=fgDCAmax;   cuts[5]=fgCPAmin;
  cuts[6]=fgRmin;     cuts[7]=fgRmax; 
}

inline void AliLightCascadeVertexer::SetDefaultMaxEta(Double_t lMaxEta) {
    fgMaxEta = lMaxEta;
}
inline void AliLightCascadeVertexer::SetDefaultMinClusters(Int_t lMinClusters) {
    fgMinClusters = lMinClusters;
}
inline void AliLightCascadeVertexer::SetDefaultUseOnTheFlyV0(Bool_t lOption) {
    fgUseOnTheFlyV0 = lOption;
}
inline void AliLightCascadeVertexer::SetMaxEta(Double_t lMaxEta) {
    fMaxEta = lMaxEta;
}
inline void AliLightCascadeVertexer::SetMinClusters(Int_t lMinClusters) {
    fMinClusters = lMinClusters;
}
inline void AliLightCascadeVertexer::SetSwitchCharges(Bool_t lOption) {
    fSwitchCharges = lOption;
}
inline void AliLightCascadeVertexer::SetUseOnTheFlyV0(Bool_t lOption) {
    fUseOnTheFlyV0 = lOption;
}

#endif

