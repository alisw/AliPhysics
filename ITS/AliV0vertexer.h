#ifndef ALIV0VERTEXER_H
#define ALIV0VERTEXER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//------------------------------------------------------------------
//                    V0 Vertexer Class
//
//   Origin: Iouri Belikov, IReS, Strasbourg, Jouri.Belikov@cern.ch 
//------------------------------------------------------------------

#include "TObject.h"

class TFile;
class AliITStrackV2;

//_____________________________________________________________________________
class AliV0vertexer : public TObject {
public:
  AliV0vertexer();
  AliV0vertexer(const Double_t cuts[7]);
  void SetCuts(const Double_t cuts[7]);

  Int_t Tracks2V0vertices(const TFile *in, TFile *out);
  Double_t PropagateToDCA(AliITStrackV2 *nt, AliITStrackV2 *pt);

  void GetCuts(Double_t cuts[7]) const;

private:
  Double_t fChi2max;      // maximal allowed chi2 
  Double_t fDNmin;        // min. allowed negative daughter's impact parameter
  Double_t fDPmin;        // min. allowed positive daughter's impact parameter
  Double_t fDCAmax;       // maximal allowed DCA between the daughter tracks 
  Double_t fCPAmax;       // maximal allowed cosine of V0's pointing angle
  Double_t fRmin, fRmax;  // max & min radii of the fiducial volume
  
  ClassDef(AliV0vertexer,1)  // V0 verterxer 
};

inline AliV0vertexer::AliV0vertexer() {
 fChi2max=33.; 
 fDNmin=0.015; fDPmin=0.015;
 fDCAmax=0.01; fCPAmax=0.025; 
 fRmin=0.5;    fRmax=2.5; 
}

inline AliV0vertexer::AliV0vertexer(const Double_t cuts[7]) {
  fChi2max=cuts[0]; 
  fDNmin=cuts[1];   fDPmin=cuts[2];
  fDCAmax=cuts[3];  fCPAmax=cuts[4];
  fRmin=cuts[5];    fRmax=cuts[6]; 
}

inline void AliV0vertexer::SetCuts(const Double_t cuts[7]) {
  fChi2max=cuts[0]; 
  fDNmin=cuts[1];   fDPmin=cuts[2];
  fDCAmax=cuts[3];  fCPAmax=cuts[4];
  fRmin=cuts[5];    fRmax=cuts[6]; 
}

inline void AliV0vertexer::GetCuts(Double_t cuts[7]) const {
  cuts[0]=fChi2max; 
  cuts[1]=fDNmin;   cuts[2]=fDPmin; 
  cuts[3]=fDCAmax;  cuts[4]=fCPAmax;
  cuts[5]=fRmin;    cuts[6]=fRmax; 
}

#endif


