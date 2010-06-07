#ifndef ALIRSNANALYSISPHI900GEV
#define ALIRSNANALYSISPHI900GEV

#include "AliAnalysisTaskSE.h"
#include "AliRsnTOFT0maker.h"

class TTree;
class AliESDEvent;
class AliESDVertex;
class AliStack;

class AliRsnAnalysisPhi900GeV : public AliAnalysisTaskSE
{
  public:

    AliRsnAnalysisPhi900GeV(const char *name = "Phi900GeV");
    AliRsnAnalysisPhi900GeV(const AliRsnAnalysisPhi900GeV& copy);
    AliRsnAnalysisPhi900GeV& operator=(const AliRsnAnalysisPhi900GeV& copy);
    virtual ~AliRsnAnalysisPhi900GeV();

    void            SetUseMC(Bool_t yn = kTRUE) {fUseMC = yn;}
    void            SetMaxDCAr(Double_t v) {fDCAr = v;}
    void            SetMaxDCAz(Double_t v) {fDCAz = v;}
    void            SetMaxChi2(Double_t v) {fChi2 = v;}
    void            SetMinNTPC(Int_t    n) {fNTPC = n;}
    void            SetTPCrange(Double_t min, Double_t max) {fMinTPC = min; fMaxTPC = max;}
    void            SetTPCpar(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4)
                    {fTPCpar[0]=p0;fTPCpar[1]=p1;fTPCpar[2]=p2;fTPCpar[3]=p3;fTPCpar[4]=p4;}

    virtual void    UserCreateOutputObjects();
    virtual void    UserExec(Option_t *option = "");
    virtual void    Terminate(Option_t *option = "");

  private:

    void     ProcessESD(AliESDEvent *esd, const AliESDVertex *v, Double_t time0, AliStack *stack);
    void     ProcessMC(AliStack *stack);
    Double_t AlephBB(Double_t p, Double_t mass = 0.493677);
    Double_t RemakeTOFtimeMC(AliESDEvent *& esd);

    Bool_t   fUseMC;
    
    Float_t  fPDG;
    Float_t  fIM;
    Float_t  fPt;
    Float_t  fY;
    Float_t  fEta;
    
    Double_t fDCAr;
    Double_t fDCAz;
    Double_t fChi2;
    Int_t    fNTPC;

    Double_t fTPCpar[5];
    Double_t fMinTPC;
    Double_t fMaxTPC;

    TTree   *fOutTree;

    Bool_t                       fTOFESD;              //  TOF flag to check if ESD data should be used
    Double_t                     fTOFSigma;            //  TOF default resolution
    AliRsnTOFT0maker            *fTOFmaker;            //! TOF time0 computator
    AliRsnTOFT0maker::ESettings  fTOFSettings;         //  TOF settings

    // ROOT dictionary
    ClassDef(AliRsnAnalysisPhi900GeV,1)
};

#endif
