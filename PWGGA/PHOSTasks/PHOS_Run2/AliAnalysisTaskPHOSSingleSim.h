#ifndef AliAnalysisTaskPHOSSingleSim_cxx
#define AliAnalysisTaskPHOSSingleSim_cxx

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

class AliAnalysisTaskPHOSSingleSim : public AliAnalysisTaskPHOSPi0EtaToGammaGamma {
  public:

    AliAnalysisTaskPHOSSingleSim(const char *name = "SingleSim");
    virtual ~AliAnalysisTaskPHOSSingleSim(); 
    void SetParticle(TString parname) {fParticleName = parname;}

  protected:
    void UserCreateOutputObjects();
    void UserExec(Option_t *option);
    void Terminate(Option_t *option);
    void GetMCInfo();
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
    AliAnalysisTaskPHOSSingleSim(const AliAnalysisTaskPHOSSingleSim&); // not implemented
    AliAnalysisTaskPHOSSingleSim& operator=(const AliAnalysisTaskPHOSSingleSim&); // not implemented

    ClassDef(AliAnalysisTaskPHOSSingleSim, 5);
};

#endif
