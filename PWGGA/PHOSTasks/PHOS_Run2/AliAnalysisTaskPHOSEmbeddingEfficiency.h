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

    virtual void ProcessMC();
    void SetWeightToClusters();
    Int_t FindCommonParent(Int_t iPart, Int_t jPart);

    void FillPhoton();
    void FillMgg();
    void FillMixMgg();
    Double_t REMB(AliAODMCParticle *p);//in cylindrical system
    Double_t RhoEMB(AliAODMCParticle *p);//in sperical system

  protected:
    TString fParticleName;
    TClonesArray *fMCArray;
    TF1 *fWeightCen0005;

  private:
    AliAnalysisTaskPHOSEmbeddingEfficiency(const AliAnalysisTaskPHOSEmbeddingEfficiency&); // not implemented
    AliAnalysisTaskPHOSEmbeddingEfficiency& operator=(const AliAnalysisTaskPHOSEmbeddingEfficiency&); // not implemented

    ClassDef(AliAnalysisTaskPHOSEmbeddingEfficiency, 8); // example of analysis
};

#endif
