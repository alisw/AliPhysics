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
#include "TObjArray.h"

class AliV0HypSel;
class TTree;
class AliESDEvent;
class AliExternalTrackParam;

//_____________________________________________________________________________
class AliV0vertexer : public TObject {
public:
  AliV0vertexer();
  void SetCuts(const Double_t *cuts);
  static void SetDefaultCuts(const Double_t *cuts);
  
  Int_t Tracks2V0vertices(AliESDEvent *event);
  
  void GetCuts(Double_t *cuts) const;
  static void GetDefaultCuts(Double_t *cuts);
  
  void SetV0HypSel(const TObjArray* selArr);
  const TObjArray* GetV0HypSelArray() const {return fV0HypSelArray;}
  
  void SetUseImprovedFinding(const Bool_t lInput);
  
  //For V0 finding improvements
  static void GetHelixCenter(const AliExternalTrackParam *track,Double_t center[2], const Double_t b);
  static Bool_t Preoptimize(const AliExternalTrackParam *nt, AliExternalTrackParam *pt, Double_t *lPreprocessxn, Double_t *lPreprocessxp, const Double_t b);
  
private:
  AliV0vertexer(const AliV0vertexer&) : fEtaMax(), fChi2max(), fDNmin(), fDPmin(), fDCAmax(), fCPAmin(), fRmin(), fRmax(),
    fUseImprovedFinding(kTRUE), fV0HypSelArray(NULL) {}
  const AliV0vertexer& operator=(const AliV0vertexer&) {return *this;}

  static
  Double_t fgEtaMax;      // maximal allowed Eta
  static
  Double_t fgChi2max;      // maximal allowed chi2
  static
  Double_t fgDNmin;        // min allowed impact parameter for the 1st daughter
  static
  Double_t fgDPmin;        // min allowed impact parameter for the 2nd daughter
  static
  Double_t fgDCAmax;       // maximal allowed DCA between the daughter tracks
  static
  Double_t fgCPAmin;       // minimal allowed cosine of V0's pointing angle
  static
  Double_t fgRmin, fgRmax; // max & min radii of the fiducial volume
  
  Double_t fEtaMax;       // max eta
  Double_t fChi2max;      // maximal allowed chi2
  Double_t fDNmin;        // min allowed impact parameter for the 1st daughter
  Double_t fDPmin;        // min allowed impact parameter for the 2nd daughter
  Double_t fDCAmax;       // maximal allowed DCA between the daughter tracks
  Double_t fCPAmin;       // minimal allowed cosine of V0's pointing angle
  Double_t fRmin, fRmax;  // max & min radii of the fiducial volume
  
  Bool_t fUseImprovedFinding; //use DCA pre-optimization + V0 refit
  
  const TObjArray* fV0HypSelArray; // array of V0 hypothesis to select
  
  ClassDef(AliV0vertexer,5)  // V0 verterxer
};

inline AliV0vertexer::AliV0vertexer() :
TObject(),
fEtaMax(fgEtaMax),
fChi2max(fgChi2max),
fDNmin(fgDNmin),
fDPmin(fgDPmin),
fDCAmax(fgDCAmax),
fCPAmin(fgCPAmin),
fRmin(fgRmin),
fRmax(fgRmax),
fUseImprovedFinding(kTRUE),
fV0HypSelArray(0)
{
}

inline void AliV0vertexer::SetCuts(const Double_t *cuts) {
  fChi2max=cuts[0]; 
  fDNmin=cuts[1];   fDPmin=cuts[2];
  fDCAmax=cuts[3];  fCPAmin=cuts[4];
  fRmin=cuts[5];    fRmax=cuts[6];
  fEtaMax = cuts[7];
}

inline void AliV0vertexer::SetDefaultCuts(const Double_t *cuts) {
  fgChi2max=cuts[0];
  fgDNmin=cuts[1];   fgDPmin=cuts[2];
  fgDCAmax=cuts[3];  fgCPAmin=cuts[4];
  fgRmin=cuts[5];    fgRmax=cuts[6];
  fgEtaMax = cuts[7];
}

inline void AliV0vertexer::GetCuts(Double_t *cuts) const {
  cuts[0]=fChi2max;
  cuts[1]=fDNmin;   cuts[2]=fDPmin;
  cuts[3]=fDCAmax;  cuts[4]=fCPAmin;
  cuts[5]=fRmin;    cuts[6]=fRmax;
  cuts[7] = fEtaMax;
}

inline void AliV0vertexer::GetDefaultCuts(Double_t *cuts) {
  cuts[0]=fgChi2max;
  cuts[1]=fgDNmin;   cuts[2]=fgDPmin;
  cuts[3]=fgDCAmax;  cuts[4]=fgCPAmin;
  cuts[5]=fgRmin;    cuts[6]=fgRmax;
  cuts[7] = fgEtaMax;
}

#endif


