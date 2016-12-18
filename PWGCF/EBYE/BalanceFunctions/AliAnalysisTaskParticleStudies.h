#ifndef ALIANALYSISTASKPARTICLESTUDIES_H
#define ALIANALYSISTASKPARTICLESTUDIES_H

// Analysis task for an analysis on rho mesons and other particles
// mainly to understand if BF behaviour could come from boosted rho mesons
// Authors: m.weber@cern.ch

class TList;
class TH2D;
class TH3D;


class AliAnalysisTaskParticleStudies : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskParticleStudies(const char *name = "AliAnalysisTaskParticleStudies");
  virtual ~AliAnalysisTaskParticleStudies(); 
   
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   FinishTaskOutput();
  virtual void   Terminate(Option_t *);


  void SetPtRange(Double_t min, Double_t max){
    fPtMin = min;
    fPtMax = max;

    Printf("Thresholds Set");
    Printf("pT = %f - %f",fPtMin,fPtMax);
  }

  void SetEtaRange(Double_t min, Double_t max){
    fEtaMin = min;
    fEtaMax = max;

    Printf("Thresholds Set");
    Printf("eta = %f - %f",fEtaMin,fEtaMax);
  }

  void SetPdgValue(Int_t pdg){ 
    fPdgCode = pdg; 
    Printf("PDG set = %d",fPdgCode);
  }

  void SetMotherPdgValue(Int_t pdg){ 
    fMotherPdgCode = pdg; 
    Printf("Mother PDG set = %d",fPdgCode);
  }



 private:
  Double_t IsEventAccepted(AliVEvent *event);
  Int_t GetAcceptedTracks(AliVEvent *event, Double_t centrality);
  
  AliAnalysisTaskParticleStudies(const AliAnalysisTaskParticleStudies&); // not implemented
  AliAnalysisTaskParticleStudies& operator=(const AliAnalysisTaskParticleStudies&); // not implemented

  TList *fListQA;//!output list for QA histograms

  Int_t    fPdgCode;//PDG value under study
  Int_t    fMotherPdgCode;//mother PDG value under study
  Double_t fPtMin;// minimum pT threshold (default = 0)
  Double_t fPtMax;// maximum pT threshold (default = 1000)
  Double_t fEtaMin;// minimum eta threshold (default = -10)
  Double_t fEtaMax;// maximum eta threshold (default = 10)

  TH2D* fHistTrackStats;//!QA histogram for track statistics
  TH2D* fHistCentStatsUsed;//!QA histogram for centralities
  TH2D* fHistIPPt;//!QA histogram for impact parameter/pt distribution
  TH3D* fHistEtaPhiPt;//!QA histogram for eta/phi/pt distribution
   
  ClassDef(AliAnalysisTaskParticleStudies, 1); //
};



#endif
