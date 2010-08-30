#ifndef ALIESDPIDCUTS_H
#define ALIESDPIDCUTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//
// Class for PID cuts
// Cuts the track based on numbers of sigmas in the detectors TPC and TOF
// The sigma cuts can be applied symmetrically or assymetrically
//

#ifndef ALIANALYSISCUTS_H
#include "AliAnalysisCuts.h"
#endif

#ifndef ALIPID_H
#include "AliPID.h"
#endif

class TCollection;
class TH1F;
class TH1I;
class TH2I;
class TList;
class AliESDtrack;
class AliESDEvent;
class AliESDpid;

class AliESDpidCuts : public AliAnalysisCuts{
  enum{
    kHasHistograms = BIT(17)
  };
  public: 
    AliESDpidCuts(const Char_t *name = "AliESDpidCuts", const Char_t *title = "");
    AliESDpidCuts(const AliESDpidCuts &ref);  // Copy constructor
    AliESDpidCuts &operator=(const AliESDpidCuts &ref);
    virtual ~AliESDpidCuts();

    virtual void Copy(TObject &c) const;
    virtual Long64_t Merge(TCollection *coll);

    Bool_t HasHistograms() const { return TestBit(kHasHistograms); }
    void DefineHistograms(Color_t color = kRed);
    void DrawHistograms();
    void SaveHistograms(const Char_t *location = NULL);
    virtual Bool_t IsSelected(TObject *);
    virtual Bool_t IsSelected(TList * /*lst*/) {return kTRUE; }
    virtual Bool_t AcceptTrack(const AliESDtrack *track, const AliESDEvent *event);

    AliESDpid *GetESDpid() { return fESDpid; };
    
    void SetTPCclusterRatioCut(Float_t clr) { fCutTPCclusterRatio = clr; }
    inline void SetTPCnSigmaCut(AliPID::EParticleType itype, Float_t nSigma);
    inline void SetTPCnSigmaCut(AliPID::EParticleType itype, Float_t negSigma, Float_t posSigma);
    inline void SetTOFnSigmaCut(AliPID::EParticleType itype, Float_t nSigma);
    inline void SetTOFnSigmaCut(AliPID::EParticleType itype, Float_t negSigma, Float_t posSigma);
    void SetMinMomentumTOF(Float_t mom) { fMinMomentumTOF = mom; }
  
  protected:
    static const Int_t kNcuts;                      // Number of Cuts
    AliESDpid *fESDpid;                             //! PID helper (n-sigma-cut)
    Char_t  fTPCsigmaCutRequired;                   // Sigma cut Requirement for TPC and Particle Species
    Char_t  fTOFsigmaCutRequired;                   // Sigma cut Requirement for TOF and Particle Species
    Float_t fCutTPCnSigma[AliPID::kSPECIES * 2];    // Species dependent cut on the distance to the TPC dE/dx line
    Float_t fCutTOFnSigma[AliPID::kSPECIES * 2];    // Species dependent cut on the distance to the TOF calculated time of flight line
    Float_t fCutTPCclusterRatio;                    // Cut on Ratio of found clusters with repect to findable clusters in the TPC
    Float_t fMinMomentumTOF;                        // Apply TOF PID only above a certain momentum

    //------------------------------------------
    // QA histograms
    TH1I *fHcutStatistics;                       // Cut Statistics
    TH2I *fHcutCorrelation;                      // Cut Correlation
    TH1F *fHclusterRatio[2];                     // TPC cluster Ratio
    TH1F *fHnSigmaTPC[AliPID::kSPECIES][2];      // TPC n-sigma cut
    TH1F *fHnSigmaTOF[AliPID::kSPECIES][2];      // TOF n-sigma cut
    //------------------------------------------
    
    ClassDef(AliESDpidCuts, 3)
};

//_____________________________________________________________________
void AliESDpidCuts::SetTPCnSigmaCut(AliPID::EParticleType itype, Float_t nSigma){ 
  //
  // symmetric sigma cut for TPC PID
  //
  fCutTPCnSigma[itype * 2]      = -nSigma;
  fCutTPCnSigma[itype * 2 + 1]  = nSigma; 
  fTPCsigmaCutRequired |= 1 << static_cast<Int_t >(itype);
}    

//_____________________________________________________________________
void AliESDpidCuts::SetTPCnSigmaCut(AliPID::EParticleType itype, Float_t negSigma, Float_t posSigma){
  //
  // assymetric sigma cut for TPC PID
  //
  fCutTPCnSigma[itype * 2]      = negSigma;
  fCutTPCnSigma[itype * 2 + 1]  = posSigma;
  fTPCsigmaCutRequired |= 1 << static_cast<Int_t >(itype);
} 

//_____________________________________________________________________
void AliESDpidCuts::SetTOFnSigmaCut(AliPID::EParticleType itype, Float_t nSigma){ 
  //
  // symmetric sigma cut for TOF PID
  //
  fCutTOFnSigma[itype * 2]      = -nSigma;
  fCutTOFnSigma[itype * 2 + 1]  = nSigma; 
  fTOFsigmaCutRequired |= 1 << static_cast<Int_t >(itype);
}    

//_____________________________________________________________________
void AliESDpidCuts::SetTOFnSigmaCut(AliPID::EParticleType itype, Float_t negSigma, Float_t posSigma){
  //
  // assymetric sigma cut for TOF PID
  //
  fCutTOFnSigma[itype * 2]      = negSigma;
  fCutTOFnSigma[itype * 2 + 1]  = posSigma;
  fTOFsigmaCutRequired |= 1 << static_cast<Int_t >(itype);
}
#endif
