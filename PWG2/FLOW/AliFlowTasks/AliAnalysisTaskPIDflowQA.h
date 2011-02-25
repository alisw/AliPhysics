//Class to extract data to do ITS+TPC global Spectra
//Autor Marek Chojnacki, Marek.Chojnacki@cern.ch
//modified Mikolaj Krzewicki, Mikolaj.Krzewicki@cern.ch

#ifndef ALIANALYSISFLOWPIDTASK_H
#define  ALIANALYSISFLOWPIDTASK_H
class TH2F;
class TProfile;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDpidCuts;
class AliESDpid;
class TGraph;
class AliStack;
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

  AliESDtrackCuts* GetAliESDtrackCuts() const {return fCuts;}
  void SetAliESDtrackCuts(AliESDtrackCuts* const cuts ){fCuts=cuts;}
  void SetEventCuts(AliFlowEventCuts* c) {fEventCuts=c;}
  void SetMCOn(){fMC=kTRUE;}
  AliESDpid* GetESDpid() const {return fESDpid;}

  Float_t Beta(Float_t m, Float_t p);

private:
  AliESDEvent *fESD;            //!ESD object    
  AliESDtrackCuts *fCuts;       //cuts 
  AliFlowEventCuts *fEventCuts; //event cuts
  AliESDpid *fESDpid;           //pid object
  Bool_t fMC;                   //if TRUE use MC 

  TH2F* fITSsignal; //!ITS signal as function of p
  TH2F* fTPCsignal; //!TPC signal as function of p
  TH2F* fTOFsignal; //!TOF signal as function of p 
  TH2F* fITSsignalpip;//!ITS PID signal as function of p for pi+
  TH2F* fTPCsignalpip;//!TPC PID signal as function of p for pi+
  TH2F* fTOFsignalpip;//!TOF PID signal as function of p for pi+
  TH2F* fITSsignalKp;//!ITS PID signal as function of p for K+
  TH2F* fTPCsignalKp;//!TPC PID signal as function of p for K+
  TH2F* fTOFsignalKp;//!TOF PID signal as function of p for K+
  TH2F* fITSsignalpp;//!ITS PID signal as function of p for p
  TH2F* fTPCsignalpp;//!TPC PID signal as function of p for p
  TH2F* fTOFsignalpp;//!TOF PID signal as function of p for p
  TH2F* fITSsignalpiMCp;//!ITS PID signal as function of p for pi+
  TH2F* fTPCsignalpiMCp;//!TPC PID signal as function of p for pi+
  TH2F* fTOFsignalpiMCp;//!TOF PID signal as function of p for pi+
  TH2F* fITSsignalKMCp;//!ITS PID signal as function of p for K+
  TH2F* fTPCsignalKMCp;//!TPC PID signal as function of p for K+
  TH2F* fTOFsignalKMCp;//!TOF PID signal as function of p for K+
  TH2F* fITSsignalpMCp;//!ITS PID signal as function of p for p
  TH2F* fTPCsignalpMCp;//!TPC PID signal as function of p for p
  TH2F* fTOFsignalpMCp;//!TOF PID signal as function of p for p
  TH2F* fTOFsignalBeta;//!vs beta
  TH2F* fTOFsignalPiBeta;//!vs beta
  TH2F* fTOFsignalKBeta;//!vs beta
  TH2F* fTOFsignalPBeta;//!vs beta
  TH2F* fPvsPt; //!P vs Pt yield
  TProfile* fMeanPvsP; //!mean p per bin
  TH2F* fTPCvsGlobalMult; //! correlation tpc only tracks vs global tracks
  AliFlowTrackCuts* fStandardGlobalCuts; //! cuts
  AliFlowTrackCuts* fStandardTPCCuts; //! cuts
  TList* fOutputList;//!output list
	
  AliAnalysisTaskPIDflowQA(const  AliAnalysisTaskPIDflowQA&); // not implemented
  AliAnalysisTaskPIDflowQA& operator=(const  AliAnalysisTaskPIDflowQA&); // not implemented
 
  void pidITS(AliESDtrack* t, Int_t pdgcode);
  void pidTPC(AliESDtrack* t, Int_t pdgcode);
  void pidTOF(AliESDtrack* t, Int_t pdgcode);
    
  ClassDef( AliAnalysisTaskPIDflowQA, 2); // example of analysis
};

#endif
