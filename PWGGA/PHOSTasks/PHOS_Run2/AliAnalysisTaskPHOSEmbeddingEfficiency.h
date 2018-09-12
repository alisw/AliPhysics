#ifndef AliAnalysisTaskPHOSEmbeddingEfficiency_cxx
#define AliAnalysisTaskPHOSEmbeddingEfficiency_cxx

// Author: Daiki Sekihata (Hiroshima University)

class TList;
class AliVEvent;
class AliESDEvent;
class AliAODEvent;
class AliPHOSGeometry;
class TClonesArray;
class AliESDtrackCuts;
class AliMultSelection;

#include "AliAnalysisTaskPHOSPi0EtaToGammaGamma.h"

class AliAnalysisTaskPHOSEmbeddingEfficiency : public AliAnalysisTaskPHOSPi0EtaToGammaGamma {
  public:

    AliAnalysisTaskPHOSEmbeddingEfficiency(const char *name = "EmbeddingEfficiency");
    virtual ~AliAnalysisTaskPHOSEmbeddingEfficiency(); 
    void SetEmbeddedParticle(TString parname) {fParticleName = parname;}

  protected:
    void UserCreateOutputObjects();
    void UserExec(Option_t *option);
    void Terminate(Option_t *option);
    void GetEmbeddedMCInfo();
    virtual void EstimatePIDCutEfficiency();

    virtual void ProcessMC();
    void SetWeightToClusters();

    void FillPhoton();
    void FillMgg();
    void FillMixMgg();
    void MCPhotonPurity();

  protected:
    TString fParticleName;

  private:
    AliAnalysisTaskPHOSEmbeddingEfficiency(const AliAnalysisTaskPHOSEmbeddingEfficiency&); // not implemented
    AliAnalysisTaskPHOSEmbeddingEfficiency& operator=(const AliAnalysisTaskPHOSEmbeddingEfficiency&); // not implemented

    ClassDef(AliAnalysisTaskPHOSEmbeddingEfficiency, 13); // example of analysis
};

#endif
