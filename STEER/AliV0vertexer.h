#ifndef ALIV0VERTEXER_H
#define ALIV0VERTEXER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//------------------------------------------------------------------
//                    V0 Vertexer Class
//            reads tracks writes out V0 vertices
//   Origin: Iouri Belikov, IReS, Strasbourg, Jouri.Belikov@cern.ch 
//------------------------------------------------------------------

#include "TObject.h"

class TTree;
class AliESD;

//_____________________________________________________________________________
class AliV0vertexer : public TObject {
public:
  AliV0vertexer();
  void SetCuts(const Double_t cuts[7]);
  static void SetDefaultCuts(const Double_t cuts[7]);

  Int_t Tracks2V0vertices(AliESD *event);

  void GetCuts(Double_t cuts[7]) const;
  static void GetDefaultCuts(Double_t cuts[7]);

private:
  static
  Double_t fgChi2max;      // maximal allowed chi2 
  static
  Double_t fgDNmin;        // min allowed impact parameter for the 1st daughter
  static
  Double_t fgDPmin;        // min allowed impact parameter for the 2nd daughter
  static
  Double_t fgDCAmax;       // maximal allowed DCA between the daughter tracks 
  static
  Double_t fgCPAmax;       // maximal allowed cosine of V0's pointing angle
  static
  Double_t fgRmin, fgRmax; // max & min radii of the fiducial volume
  
  Double_t fChi2max;      // maximal allowed chi2 
  Double_t fDNmin;        // min allowed impact parameter for the 1st daughter
  Double_t fDPmin;        // min allowed impact parameter for the 2nd daughter
  Double_t fDCAmax;       // maximal allowed DCA between the daughter tracks 
  Double_t fCPAmax;       // maximal allowed cosine of V0's pointing angle
  Double_t fRmin, fRmax;  // max & min radii of the fiducial volume
  
  ClassDef(AliV0vertexer,2)  // V0 verterxer 
};

inline AliV0vertexer::AliV0vertexer() :
  TObject(),
  fChi2max(fgChi2max), 
  fDNmin(fgDNmin),
  fDPmin(fgDPmin),
  fDCAmax(fgDCAmax),
  fCPAmax(fgCPAmax), 
  fRmin(fgRmin),
  fRmax(fgRmax) 
{
}

inline void AliV0vertexer::SetCuts(const Double_t cuts[7]) {
  fChi2max=cuts[0]; 
  fDNmin=cuts[1];   fDPmin=cuts[2];
  fDCAmax=cuts[3];  fCPAmax=cuts[4];
  fRmin=cuts[5];    fRmax=cuts[6]; 
}

inline void AliV0vertexer::SetDefaultCuts(const Double_t cuts[7]) {
  fgChi2max=cuts[0]; 
  fgDNmin=cuts[1];   fgDPmin=cuts[2];
  fgDCAmax=cuts[3];  fgCPAmax=cuts[4];
  fgRmin=cuts[5];    fgRmax=cuts[6]; 
}

inline void AliV0vertexer::GetCuts(Double_t cuts[7]) const {
  cuts[0]=fChi2max; 
  cuts[1]=fDNmin;   cuts[2]=fDPmin; 
  cuts[3]=fDCAmax;  cuts[4]=fCPAmax;
  cuts[5]=fRmin;    cuts[6]=fRmax; 
}

inline void AliV0vertexer::GetDefaultCuts(Double_t cuts[7]) {
  cuts[0]=fgChi2max; 
  cuts[1]=fgDNmin;   cuts[2]=fgDPmin; 
  cuts[3]=fgDCAmax;  cuts[4]=fgCPAmax;
  cuts[5]=fgRmin;    cuts[6]=fgRmax; 
}

#endif


