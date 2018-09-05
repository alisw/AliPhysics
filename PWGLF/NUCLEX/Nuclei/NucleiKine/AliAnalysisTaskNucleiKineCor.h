#ifndef ALIANALYSISTASKNUCLEIKINECOR_H
#define ALIANALYSISTASKNUCLEIKINECOR_H

#include <AliAnalysisTaskSE.h>
#include <AliGenLightNuclei.h>

class TH1D;
class TH3D;
class TList;
#include <vector>

class AliAnalysisTaskNucleiKineCor : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskNucleiKineCor(const char* name = "AliAnalysisTaskNucleiKineCor");
    virtual ~AliAnalysisTaskNucleiKineCor();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *opt) {}
    void SetPt(Double_t pt) {fPt = pt;}
    void SetP0(Double_t p0) {fP0 = p0;}

    std::vector<int> fPdgCodes;
    enum Species {
      kPiPlus, kPiMinus, kKplus, kKminus, kProton, kAntiProton, kNeutron, kAntiNeutron, kDeuteron, kAntiDeuteron,
      kLambda, kAntiLambda, kXiMinus, kXiPlus, kOmegaMinus, kOmegaPlus
    };
    std::vector<std::string> fParticleNames;

  protected:
    Double_t DeltaPhi(Double_t phia, Double_t phib,
		      Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2) const;
    Double_t          fPt;          //  pt of trigger hadron
    Double_t          fP0;          //  coalescence momentum
    TList*            fOutputList;  //! output list for histograms
    TH1D*             fEvents;      //! number of events
    TH1D*             fPtHist;      //! pt dist
    TH1D*             fPtLead;      //! leading pt dist
    TH1D*             fPtPro;       //! proton pt dist
    TH1D*             fPtDeu;       //! deuteron pt dist
    TH1*              fHists[999];  //! hists
    AliGenLightNuclei fAfterBurner; // afterburner

    AliAnalysisTaskNucleiKineCor(const AliAnalysisTaskNucleiKineCor& other);
    AliAnalysisTaskNucleiKineCor& operator=(const AliAnalysisTaskNucleiKineCor& other);

  ClassDef(AliAnalysisTaskNucleiKineCor, 2)
};
#endif
