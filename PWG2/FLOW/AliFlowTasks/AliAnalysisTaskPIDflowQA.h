//Class to extract data to do ITS+TPC global Spectra
//Autor Marek Chojnacki, Marek.Chojnacki@cern.ch
//modified Mikolaj Krzewicki, Mikolaj.Krzewicki@cern.ch

#ifndef ALIANALYSISFLOWPIDTASK_H
#define  ALIANALYSISFLOWPIDTASK_H
class TH2F;
class TH1F;
class TProfile;
class AliESDEvent;
class AliESDtrack;
class AliESDpid;
class AliFlowEventCuts;
class AliFlowTrackCuts;
#include "AliAnalysisTaskSE.h"
#include "AliESDpid.h"

class  AliAnalysisTaskPIDflowQA : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskPIDflowQA();
  AliAnalysisTaskPIDflowQA(const char *name);
  virtual ~AliAnalysisTaskPIDflowQA() {}
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *); 

  AliFlowTrackCuts* GetAliESDtrackCuts() const {return fCuts;}
  void SetTrackCuts(AliFlowTrackCuts* cuts ){fCuts=cuts;}
  void SetEventCuts(AliFlowEventCuts* c) {fEventCuts=c;}
  AliESDpid* GetESDpid() const {return fESDpid;}
  void SetUseDebugFile(Bool_t b=kTRUE) {fUseDebugFile=b;}
  Bool_t TPCTOFagree(const AliESDtrack *track);  

  Float_t Beta(Float_t m, Float_t p);

private:
  AliESDEvent *fESD;            //!ESD object    
  AliFlowTrackCuts *fCuts;       //cuts 
  AliFlowEventCuts *fEventCuts; //event cuts
  AliESDpid *fESDpid;           //pid object
  Bool_t fUseDebugFile; //write debug file
  FILE* fFile; //debug output file

  TH2F* fTPCsignal; //!TPC signal as function of p
  TH2F* fTPCsignalPi;//!TPC PID signal as function of p for pi+
  TH2F* fTPCsignalK;//!TPC PID signal as function of p for K+
  TH2F* fTPCsignalP;//!TPC PID signal as function of p for p

  TH2F* fTPCsignalPimc;//!TPC PID signal as function of p for pi+
  TH2F* fTPCsignalKmc;//!TPC PID signal as function of p for K+
  TH2F* fTPCsignalPmc;//!TPC PID signal as function of p for p

  TH2F* fTOFtime;//!vs time
  TH2F* fTOFtimeE;//!vs time
  TH2F* fTOFtimePi;//!vs time
  TH2F* fTOFtimeK;//!vs time
  TH2F* fTOFtimeP;//!vs time

  TH2F* fTOFbeta;//!vs beta
  TH2F* fTOFbetaE;//!vs beta
  TH2F* fTOFbetaPi;//!vs beta
  TH2F* fTOFbetaK;//!vs beta
  TH2F* fTOFbetaP;//!vs beta

  TH2F* fTOFinvbeta;//!vs beta
  TH2F* fTOFinvbetaE;//!vs beta
  TH2F* fTOFinvbetaPi;//!vs beta
  TH2F* fTOFinvbetaK;//!vs beta
  TH2F* fTOFinvbetaP;//!vs beta

  TH2F* fTOFrawtime;//!vs time
  TH2F* fTOFrawtimeE;//!vs time
  TH2F* fTOFrawtimePi;//!vs time
  TH2F* fTOFrawtimeK;//!vs time
  TH2F* fTOFrawtimeP;//!vs time

  TH2F* fTOFrawbeta;//!vs beta
  TH2F* fTOFrawbetaE;//!vs beta
  TH2F* fTOFrawbetaPi;//!vs beta
  TH2F* fTOFrawbetaK;//!vs beta
  TH2F* fTOFrawbetaP;//!vs beta

  TH2F* fTOFrawinvbeta;//!vs beta
  TH2F* fTOFrawinvbetaE;//!vs beta
  TH2F* fTOFrawinvbetaPi;//!vs beta
  TH2F* fTOFrawinvbetaK;//!vs beta
  TH2F* fTOFrawinvbetaP;//!vs beta

  TH2F* fPvsPt; //!P vs Pt yield
  TProfile* fMeanPvsP; //!mean p per bin
  TH2F* fTPCvsGlobalMult; //! correlation tpc only tracks vs global tracks
  AliFlowTrackCuts* fStandardGlobalCuts; //! cuts
  AliFlowTrackCuts* fStandardTPCCuts; //! cuts

  AliFlowTrackCuts* fCutsTOFbetaElectrons; //!
  AliFlowTrackCuts* fCutsTOFbetaPions; //!
  AliFlowTrackCuts* fCutsTOFbetaKaons; //!
  AliFlowTrackCuts* fCutsTOFbetaProtons; //!

  AliFlowTrackCuts* fCutsTOFbetaSimpleElectrons; //!
  AliFlowTrackCuts* fCutsTOFbetaSimplePions; //!
  AliFlowTrackCuts* fCutsTOFbetaSimpleKaons; //!
  AliFlowTrackCuts* fCutsTOFbetaSimpleProtons; //!

  AliFlowTrackCuts* fCutsTOFbayesianElectrons; //!
  AliFlowTrackCuts* fCutsTOFbayesianPions; //!
  AliFlowTrackCuts* fCutsTOFbayesianKaons; //!
  AliFlowTrackCuts* fCutsTOFbayesianProtons; //!

  AliFlowTrackCuts* fCutsTPCdedxElectrons; //!
  AliFlowTrackCuts* fCutsTPCdedxPions; //!
  AliFlowTrackCuts* fCutsTPCdedxKaons; //!
  AliFlowTrackCuts* fCutsTPCdedxProtons; //!

  AliFlowTrackCuts* fCutsTPCpidElectrons; //!
  AliFlowTrackCuts* fCutsTPCpidPions; //!
  AliFlowTrackCuts* fCutsTPCpidKaons; //!
  AliFlowTrackCuts* fCutsTPCpidProtons; //!

  AliFlowTrackCuts* fCutsTPCbayesianElectrons; //!
  AliFlowTrackCuts* fCutsTPCbayesianPions; //!
  AliFlowTrackCuts* fCutsTPCbayesianKaons; //!
  AliFlowTrackCuts* fCutsTPCbayesianProtons; //!

  AliFlowTrackCuts* fCutsMCelectrons;
  AliFlowTrackCuts* fCutsMCpions;
  AliFlowTrackCuts* fCutsMCkaons;
  AliFlowTrackCuts* fCutsMCprotons;
  AliFlowTrackCuts* fCutsMCprimaryelectrons;
  AliFlowTrackCuts* fCutsMCprimarypions;
  AliFlowTrackCuts* fCutsMCprimarykaons;
  AliFlowTrackCuts* fCutsMCprimaryprotons;

  TList* fOutputList;//!output list
	
  AliAnalysisTaskPIDflowQA(const  AliAnalysisTaskPIDflowQA&); // not implemented
  AliAnalysisTaskPIDflowQA& operator=(const  AliAnalysisTaskPIDflowQA&); // not implemented
 
  void pidITS(AliESDtrack* t, Int_t pdgcode);
  void pidTPC(AliESDtrack* t, Int_t pdgcode);
  void pidTOF(AliESDtrack* t, Int_t pdgcode);
    
  ClassDef( AliAnalysisTaskPIDflowQA, 4); // example of analysis
};

#endif
