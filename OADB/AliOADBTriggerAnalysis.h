#ifndef AliOADBTriggerAnalysis_H
#define AliOADBTriggerAnalysis_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
// OADB container for filling scheme information (BX ids, name ...)
// Author: Michele Floris, CERN
// Current support and development: Evgeny Kryshen, PNPI
//-------------------------------------------------------------------------

#include "TNamed.h"

class AliOADBTriggerAnalysis : public TNamed {

 public :
  AliOADBTriggerAnalysis(TString name="default");
  virtual ~AliOADBTriggerAnalysis();
  
  // Getters
  Float_t GetZDCCutRefSumCorr()     { return fZDCCutRefSumCorr;     }
  Float_t GetZDCCutRefDeltaCorr()   { return fZDCCutRefDeltaCorr;   }
  Float_t GetZDCCutSigmaSumCorr()   { return fZDCCutSigmaSumCorr;   }
  Float_t GetZDCCutSigmaDeltaCorr() { return fZDCCutSigmaDeltaCorr; }
  Float_t GetZDCCutZNATimeCorrMax() { return fZDCCutZNATimeCorrMax; }
  Float_t GetZDCCutZNATimeCorrMin() { return fZDCCutZNATimeCorrMin; }
  Float_t GetZDCCutZNCTimeCorrMax() { return fZDCCutZNCTimeCorrMax; }
  Float_t GetZDCCutZNCTimeCorrMin() { return fZDCCutZNCTimeCorrMin; }
  Float_t GetSPDClsVsTklA()         { return fSPDClsVsTklA;         }
  Float_t GetSPDClsVsTklB()         { return fSPDClsVsTklB;         }
  Float_t GetV0C012vsTklA()         { return fV0C012vsTklA;         }
  Float_t GetV0C012vsTklB()         { return fV0C012vsTklB;         }
  Float_t GetV0MOnVsOfA()           { return fV0MOnVsOfA;           }
  Float_t GetV0MOnVsOfB()           { return fV0MOnVsOfB;           }
  Float_t GetSPDOnVsOfA()           { return fSPDOnVsOfA;           }
  Float_t GetSPDOnVsOfB()           { return fSPDOnVsOfB;           }
  Int_t   GetVtxMinContributors()   { return fVtxMinContributors;   }
  Float_t GetVtxMinZdist()          { return fVtxMinZdist;          }
  Float_t GetVtxNSigmaZdist()       { return fVtxNSigmaZdist;       }
  Float_t GetVtxNSigmaDiamXY()      { return fVtxNSigmaDiamXY;      }
  Float_t GetVtxNSigmaDiamZ()       { return fVtxNSigmaDiamZ;       }
  Float_t GetV0CasymA()             { return fV0CasymA;             }
  Float_t GetV0CasymB()             { return fV0CasymB;             }
  Int_t   GetNBCsPast()             { return fNBCsPast;             }
  Int_t   GetNBCsFuture()           { return fNBCsFuture;           }
  Int_t   GetVIRBBAflags()          { return fVIRBBAflags;          }
  Int_t   GetVIRBBCflags()          { return fVIRBBCflags;          }
  Int_t   GetVIRBGAflags()          { return fVIRBGAflags;          }
  Int_t   GetVIRBGCflags()          { return fVIRBGCflags;          }
  Int_t   GetVHMBBAflags()          { return fVHMBBAflags;          }
  Int_t   GetVHMBBCflags()          { return fVHMBBCflags;          }
  Int_t   GetVHMBGAflags()          { return fVHMBGAflags;          }
  Int_t   GetVHMBGCflags()          { return fVHMBGCflags;          }
  Int_t   GetV0MOnThreshold()       { return fV0MOnThreshold;       }
  Float_t GetV0MOfThreshold()       { return fV0MOfThreshold;       }
  Int_t   GetSPDGFOThreshhold()     { return fSPDGFOThreshold;      }
  Int_t   GetSH1OuterThreshold()    { return fSH1OuterThreshold;    }
  Int_t   GetSH2OuterThreshold()    { return fSH2OuterThreshold;    }
  Int_t   GetTklThreshold()         { return fTklThreshold;         }
  Float_t GetFMDLowThreshold()      { return fFMDLowCut;            }
  Float_t GetFMDHitThreshold()      { return fFMDHitCut;            }
  Float_t GetTRDptHSE()             { return fTRDptHSE;             }
  UChar_t GetTRDpidHSE()            { return fTRDpidHSE;            }
  Float_t GetTRDptHQU()             { return fTRDptHQU;             }
  UChar_t GetTRDpidHQU()            { return fTRDpidHQU;            }
  Float_t GetTRDptHEE()             { return fTRDptHEE;             }
  UChar_t GetTRDpidHEE()            { return fTRDpidHEE;            }
  UChar_t GetTRDminSectorHEE()      { return fTRDminSectorHEE;      }
  UChar_t GetTRDmaxSectorHEE()      { return fTRDmaxSectorHEE;      }
  Float_t GetTRDptHJT()             { return fTRDptHJT;             }
  UChar_t GetTRDnHJT()              { return fTRDnHJT;              }
  
  // Setters
  void SetSPDClsVsTklA(Float_t val)     { fSPDClsVsTklA       = val; }
  void SetSPDClsVsTklB(Float_t val)     { fSPDClsVsTklB       = val; }
  void SetV0C012vsTklA(Float_t val)     { fV0C012vsTklA       = val; }
  void SetV0C012vsTklB(Float_t val)     { fV0C012vsTklB       = val; }
  void SetV0MOnVsOfA(Float_t val)       { fV0MOnVsOfA         = val; }
  void SetV0MOnVsOfB(Float_t val)       { fV0MOnVsOfB         = val; }
  void SetSPDOnVsOfA(Float_t val)       { fSPDOnVsOfA         = val; }
  void SetSPDOnVsOfB(Float_t val)       { fSPDOnVsOfB         = val; }
  void SetVtxMinContributors(Int_t val) { fVtxMinContributors = val; }
  void SetVtxMinZdist(Float_t val)      { fVtxMinZdist        = val; }
  void SetVtxNSigmaZdist(Float_t val)   { fVtxNSigmaZdist     = val; }
  void SetVtxNSigmaDiamXY(Float_t val)  { fVtxNSigmaDiamXY    = val; }
  void SetVtxNSigmaDiamZ(Float_t val)   { fVtxNSigmaDiamZ     = val; }
  void SetV0CasymA(Float_t val)         { fV0CasymA           = val; }
  void SetV0CasymB(Float_t val)         { fV0CasymB           = val; }
  void SetNBCsPast(Int_t val)           { fNBCsPast           = val; }
  void SetNBCsFuture(Int_t val)         { fNBCsFuture         = val; }
  void SetVIRBBAflags(Int_t val)        { fVIRBBAflags        = val; }
  void SetVIRBBCflags(Int_t val)        { fVIRBBCflags        = val; }
  void SetVIRBGAflags(Int_t val)        { fVIRBGAflags        = val; }
  void SetVIRBGCflags(Int_t val)        { fVIRBGCflags        = val; }
  void SetVHMBBAflags(Int_t val)        { fVHMBBAflags        = val; }
  void SetVHMBBCflags(Int_t val)        { fVHMBBCflags        = val; }
  void SetVHMBGAflags(Int_t val)        { fVHMBGAflags        = val; }
  void SetVHMBGCflags(Int_t val)        { fVHMBGCflags        = val; }
  void SetV0MOnThreshold(Int_t val)     { fV0MOnThreshold     = val; }
  void SetV0MOfThreshold(Float_t val)   { fV0MOfThreshold     = val; }
  void SetSPDGFOThreshhold(Int_t val)   { fSPDGFOThreshold    = val; }
  void SetSH1OuterThreshold(Int_t val)  { fSH1OuterThreshold  = val; }
  void SetSH2OuterThreshold(Int_t val)  { fSH2OuterThreshold  = val; }
  void SetTklThreshold(Int_t val)       { fTklThreshold       = val; }
  
  void SetZDCCorrParameters(Float_t sumCorr, Float_t deltaCorr, Float_t sigmaSumCorr, Float_t sigmaDeltaCorr){ 
    fZDCCutRefSumCorr = sumCorr; 
    fZDCCutRefDeltaCorr = deltaCorr; 
    fZDCCutSigmaSumCorr = sigmaSumCorr; 
    fZDCCutSigmaDeltaCorr = sigmaDeltaCorr;
  }
  void SetZNCorrParameters(Float_t znaTimeCorrMin, Float_t znaTimeCorrMax, Float_t zncTimeCorrMin, Float_t zncTimeCorrMax){ 
    fZDCCutZNATimeCorrMin = znaTimeCorrMin; 
    fZDCCutZNATimeCorrMax = znaTimeCorrMax; 
    fZDCCutZNCTimeCorrMin = zncTimeCorrMin; 
    fZDCCutZNCTimeCorrMax = zncTimeCorrMax;
  }
  void SetTRDTriggerParameters(Float_t ptHSE, UChar_t pidHSE, Float_t ptHQU, UChar_t pidHQU, Float_t ptHEE, UChar_t pidHEE, UChar_t minSectorHEE, UChar_t maxSectorHEE, Float_t ptHJT, UChar_t nHJT) {
    fTRDptHSE = ptHSE; fTRDpidHSE = pidHSE;
    fTRDptHQU = ptHQU; fTRDpidHQU = pidHQU;
    fTRDptHEE = ptHEE; fTRDpidHEE = pidHEE;
    fTRDminSectorHEE = minSectorHEE; fTRDmaxSectorHEE = maxSectorHEE;
    fTRDptHJT = ptHJT; fTRDnHJT = nHJT;
  }
  void SetFMDThreshold(Float_t low, Float_t hit) { fFMDLowCut = low; fFMDHitCut = hit; }

  virtual Bool_t IsFolder() const { return kTRUE; }
  void Browse(TBrowser *b);
  virtual void Print(Option_t* option = "") const;

 protected:
  Float_t fFMDLowCut;             // 
  Float_t fFMDHitCut;             // 
  Float_t fZDCCutRefSum;          // ZDC time cut configuration
  Float_t fZDCCutRefDelta;        // ZDC time cut configuration
  Float_t fZDCCutSigmaSum;        // ZDC time cut configuration
  Float_t fZDCCutSigmaDelta;      // ZDC time cut configuration
  Float_t fZDCCutRefSumCorr;      // Corrected ZDC time cut configuration
  Float_t fZDCCutRefDeltaCorr;    // Corrected ZDC time cut configuration
  Float_t fZDCCutSigmaSumCorr;    // Corrected ZDC time cut configuration
  Float_t fZDCCutSigmaDeltaCorr;  // Corrected ZDC time cut configuration  
  Float_t fZDCCutZNATimeCorrMin;  // Corrected ZNA minimum time cut configuration
  Float_t fZDCCutZNATimeCorrMax;  // Corrected ZNA maximum time cut configuration
  Float_t fZDCCutZNCTimeCorrMin;  // Corrected ZNC minimum time cut configuration
  Float_t fZDCCutZNCTimeCorrMax;  // Corrected ZNC maximum time cut configuration
  Float_t fSPDClsVsTklA;          // constant for the linear cut in SPD clusters vs tracklets
  Float_t fSPDClsVsTklB;          // slope for the linear cut in SPD clusters vs tracklets
  Float_t fV0C012vsTklA;          // constant for the linear cut in V0C012 vs tracklets
  Float_t fV0C012vsTklB;          // slope for the linear cut in V0C012 vs tracklets
  Float_t fV0MOnVsOfA;            // constant for the linear pileup cut in Online vs Offline V0M
  Float_t fV0MOnVsOfB;            // slope for the linear pileup cut in Online vs Offline V0M
  Float_t fSPDOnVsOfA;            // constant for the linear pileup cut in Online vs Offline SPD
  Float_t fSPDOnVsOfB;            // slope for the linear pileup cut in Online vs Offline SPD
  Int_t   fVtxMinContributors;    // SPD vertex pileup cut: minimum number of contributors
  Float_t fVtxMinZdist;           // SPD vertex pileup cut: minimum z-vertex distance
  Float_t fVtxNSigmaZdist;        // SPD vertex pileup cut: n sigma distrance
  Float_t fVtxNSigmaDiamXY;       // SPD vertex pileup cut: n sigma xy diam
  Float_t fVtxNSigmaDiamZ;        // SPD vertex pileup cut: n sigma z diam
  Float_t fV0CasymA;              // constant for the linear cut on V0C012 vs V0C3 asymmetry
  Float_t fV0CasymB;              // slope for the linear cut on V0C012 vs V0C3 asymmetry
  Int_t fNBCsPast;                // VIR past-future protection: number of past BCs (BC%4=0)
  Int_t fNBCsFuture;              // VIR past-future protection: number of future BCs (BC%4=0)
  Int_t fVIRBBAflags;             // VIR past-future protection: min number of BBA flags in VIR definition
  Int_t fVIRBBCflags;             // VIR past-future protection: min number of BBC flags in VIR definition
  Int_t fVIRBGAflags;             // VIR past-future protection: min number of BGA flags in VIR definition
  Int_t fVIRBGCflags;             // VIR past-future protection: min number of BGC flags in VIR definition
  Int_t fVHMBBAflags;             // VHM trigger: min number of BBA flags (read out from OCDB)
  Int_t fVHMBBCflags;             // VHM trigger: min number of BBC flags (read out from OCDB)
  Int_t fVHMBGAflags;             // VHM trigger: min number of BGA flags (read out from OCDB)
  Int_t fVHMBGCflags;             // VHM trigger: min number of BGC flags (read out from OCDB)
  Int_t fV0MOnThreshold;          // V0M HM trigger: min V0M threshold (read out from OCDB)
  Float_t fV0MOfThreshold;        // V0M HM offline: min V0M threshold
  Int_t fSPDGFOThreshold;         // SPD GFO trigger: min number of outer chips
  Int_t fSH1OuterThreshold;       // SPD 0SH1 trigger: min number of outer chips
  Int_t fSH2OuterThreshold;       // SPD 0SH2 trigger: min number of outer chips
  Int_t fTklThreshold;            // Offline cut on number of tracklets (for high-multiplicity SPD trigger)
  Float_t fTRDptHSE;              // pt threshold for HSE trigger
  UChar_t fTRDpidHSE;             // PID threshold for HSE trigger
  Float_t fTRDptHQU;              // pt threshold for HQU trigger
  UChar_t fTRDpidHQU;             // PID threshold for HQU trigger
  Float_t fTRDptHEE;              // pt threshold for HEE trigger
  UChar_t fTRDpidHEE;             // PID threshold for HEE trigger
  UChar_t fTRDminSectorHEE;       // min sector for HEE trigger
  UChar_t fTRDmaxSectorHEE;       // max sector for HEE trigger
  Float_t fTRDptHJT;              // pt threshold for HJT trigger
  UChar_t fTRDnHJT;               // no of track threshold for HJT trigger

  AliOADBTriggerAnalysis(const AliOADBTriggerAnalysis& cont);  // not implemented
  AliOADBTriggerAnalysis& operator=(const AliOADBTriggerAnalysis& cont); // not implemented

  ClassDef(AliOADBTriggerAnalysis, 6);
};

#endif
