#ifndef ALIANALYSISTASKNUCLEIKINE_H
#define ALIANALYSISTASKNUCLEIKINE_H

#include "AliAnalysisTaskSE.h"

#include <TMath.h>

#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/LorentzVector.h"
using ROOT::Math::PxPyPzE4D;
using ROOT::Math::LorentzVector;

#include <string>
using std::string;

#include <vector>
using std::vector;

class TH1D;
class TH2D;
class TList;


class AliAnalysisTaskNucleiKine: public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskNucleiKine(const string name = "AliAnalysisTaskNucleiKine");
    virtual ~AliAnalysisTaskNucleiKine();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *opt) {}

    template<typename F> static F GetPcm(const LorentzVector<PxPyPzE4D<F> > &a, const LorentzVector<PxPyPzE4D<F> > &b);

    /// Coalescence parameters
    float fPmax;
    float fRmax;
    float fSpinProb;

    // POIs
    vector<int> fPdgCodes;
    enum Species {
      kPiPlus, kPiMinus, kKplus, kKminus, kProton, kAntiProton, kNeutron, kAntiNeutron, kDeuteron, kAntiDeuteron
    };

    struct Particle {
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > pos;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > mom;
    };

  protected:
    AliAnalysisTaskNucleiKine(const AliAnalysisTaskNucleiKine& other);
    AliAnalysisTaskNucleiKine& operator=(const AliAnalysisTaskNucleiKine& other);

    TList* fOutputList;    //! output list for histograms

    TH1D*  fEventCounter;
    TH2D*  fPtSpectra;

    ClassDef(AliAnalysisTaskNucleiKine, 1)
};

template<typename F> F AliAnalysisTaskNucleiKine::GetPcm(const ROOT::Math::LorentzVector<PxPyPzE4D<F> > &a, const ROOT::Math::LorentzVector<PxPyPzE4D<F> > &b) {
  LorentzVector<PxPyPzE4D<F> > c = a + b;
  const F s = c.mass2();
  const F dm2 = (a.mass() - b.mass()) * (a.mass() - b.mass());
  return TMath::Abs(s - dm2) / (2.f * c.mass());
}

#endif
