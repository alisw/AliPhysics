#ifndef ALIANALYSISTASKNUCLEIKINE_H
#define ALIANALYSISTASKNUCLEIKINE_H

#include <AliAnalysisTaskSE.h>
#include <AliGenLightNuclei.h>

class TH1D;
class TH3D;
class TList;
#include <vector>

class AliAnalysisTaskNucleiKine : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskNucleiKine(const char* name = "AliAnalysisTaskNucleiKine");
    virtual ~AliAnalysisTaskNucleiKine();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *opt) {}

    // POIs
    std::vector<int> fPdgCodes;
    enum Species {
      kPiPlus, kPiMinus, kKplus, kKminus, kProton, kAntiProton, kNeutron, kAntiNeutron, kDeuteron, kAntiDeuteron,
      kLambda, kAntiLambda, kXiMinus, kXiPlus, kOmegaMinus, kOmegaPlus
    };
    std::vector<std::string> fParticleNames;

    bool   fIgnoreCentrality;
    bool   fUseAfterburner;
    AliGenLightNuclei fAfterburner; // Afterburner

  protected:
    AliAnalysisTaskNucleiKine(const AliAnalysisTaskNucleiKine& other);
    AliAnalysisTaskNucleiKine& operator=(const AliAnalysisTaskNucleiKine& other);


    TList* fOutputList;    //! output list for histograms

    TH1D*  fReactionPlane;       //!
    TH1D*  fEventCounter;        //!
    TH3D*  fPtSpectra;           //!
    TH3D*  fPtSpectraNoRapidity; //!

    ClassDef(AliAnalysisTaskNucleiKine, 1)
};



#endif
