#ifndef ALIESDPIDCUTS_H
#define ALIESDPIDCUTS_H

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
class AliTPCpidESD;
class AliTOFpidESD;

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
    virtual Bool_t IsSelected(TObject *o);
    virtual Bool_t IsSelected(TList *) { return kTRUE; }
    virtual Bool_t AcceptTrack(const AliESDtrack *track);
    
    void SetTPCclusterRatioCut(Float_t clr) { fCutTPCclusterRatio = clr; }
    void SetTPCnSigmaCut(AliPID::EParticleType itype, Float_t nSigma){ fCutTPCnSigma[itype] = nSigma; }
    void SetTOFnSigmaCut(AliPID::EParticleType itype, Float_t nSigma){ fCutTOFnSigma[itype] = nSigma; }
    void SetMinMomentumTOF(Float_t mom) { fMinMomentumTOF = mom; }
  
  protected:
    static const Int_t kNcuts;                      // Number of Cuts
    AliTPCpidESD *fTPCpid;                          //! TPC PID (n-sigma cut)
    AliTOFpidESD *fTOFpid;                          //! TOF PID (n-sigma-cut)
    Float_t fCutTPCnSigma[AliPID::kSPECIES];        // Species dependent cut on the distance to the TPC dE/dx line
    Float_t fCutTOFnSigma[AliPID::kSPECIES];        // Species dependent cut on the distance to the TOF calculated time of flight line
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
    
    ClassDef(AliESDpidCuts, 1)
};
#endif
