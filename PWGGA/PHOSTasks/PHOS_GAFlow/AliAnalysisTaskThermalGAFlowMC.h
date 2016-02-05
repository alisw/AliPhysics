#ifndef AliAnalysisTaskThermalGAFlowMC_h
#define AliAnalysisTaskThermalGAFlowMC_h

class TParticle;
class AliStack;

#include "AliAnalysisTaskThermalGAFlow.h"

class AliAnalysisTaskThermalGAFlowMC : public AliAnalysisTaskThermalGAFlow {

  public:
    AliAnalysisTaskThermalGAFlowMC( const char *name = "AliAnalysisTaskThermalGAFlow");
    virtual ~AliAnalysisTaskThermalGAFlowMC();
    virtual void UserCreateOutputObjects();
    virtual void UserExec( Option_t *);
    virtual void Terminate( Option_t *);

  private:
    //Global Variables
  
    AliStack *fStack;

    //Functions

    AliStack* GetMCStack();
    Bool_t IsPhotonShower(AliESDCaloCluster *cluster);
    void ScanClustersMC();
    void ExtractRealPions();
    void GuessPure(AliESDCaloCluster *cluster);
    void TasteFlavor();

    void CaptainsLog(Int_t s);
    void FillHist(const char * key, Double_t x, Double_t y) const;
    void FillHist(const char * key, Double_t x) const;
    void FillHist(const char * key, Double_t x, Double_t y, Double_t z) const;


//*******************More CERN Voodoo to make the class def work***********//
AliAnalysisTaskThermalGAFlowMC(const AliAnalysisTaskThermalGAFlowMC & g) ; //Copy Constructor
AliAnalysisTaskThermalGAFlowMC& operator=(const AliAnalysisTaskThermalGAFlowMC&); //Mystery Line
ClassDef(AliAnalysisTaskThermalGAFlowMC, 4);

};

#endif
