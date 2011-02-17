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
  TH2F* fTOFsignal; //!TOF signal as function of  p 
  TH2F* fITSsignalpi;//!ITS PID signal as function of pt for pi+
  TH2F* fTPCsignalpi;//!TPC PID signal as function of pt for pi+
  TH2F* fTOFsignalpi;//!TOF PID signal as function of pt for pi+
  TH2F* fITSsignalK;//!ITS PID signal as function of pt for K+
  TH2F* fTPCsignalK;//!TPC PID signal as function of pt for K+
  TH2F* fTOFsignalK;//!TOF PID signal as function of pt for K+
  TH2F* fITSsignalp;//!ITS PID signal as function of pt for p
  TH2F* fTPCsignalp;//!TPC PID signal as function of pt for p
  TH2F* fTOFsignalp;//!TOF PID signal as function of pt for p
  TH2F* fITSsignalpiMC;//!ITS PID signal as function of pt for pi+
  TH2F* fTPCsignalpiMC;//!TPC PID signal as function of pt for pi+
  TH2F* fTOFsignalpiMC;//!TOF PID signal as function of pt for pi+
  TH2F* fITSsignalKMC;//!ITS PID signal as function of pt for K+
  TH2F* fTPCsignalKMC;//!TPC PID signal as function of pt for K+
  TH2F* fTOFsignalKMC;//!TOF PID signal as function of pt for K+
  TH2F* fITSsignalpMC;//!ITS PID signal as function of pt for p
  TH2F* fTPCsignalpMC;//!TPC PID signal as function of pt for p
  TH2F* fTOFsignalpMC;//!TOF PID signal as function of pt for p
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
  TH2F* fTOFsignalPiExpKvsPt;//!TOF expected signal
  TH2F* fTOFsignalPiExpPvsPt;//!TOF expected signal
  TH2F* fTOFsignalKExpPivsPt;//!TOF expected signal
  TH2F* fTOFsignalKExpPvsPt;//!TOF expected signal
  TH2F* fTOFsignalPExpPivsPt;//!TOF expected signal
  TH2F* fTOFsignalPExpKvsPt;//!TOF expected signal
  TH2F* fTOFsignalPiExpKvsP;//!TOF expected signal
  TH2F* fTOFsignalPiExpPvsP;//!TOF expected signal
  TH2F* fTOFsignalKExpPivsP;//!TOF expected signal
  TH2F* fTOFsignalKExpPvsP;//!TOF expected signal
  TH2F* fTOFsignalPExpPivsP;//!TOF expected signal
  TH2F* fTOFsignalPExpKvsP;//!TOF expected signal
  TH2F* fTOFsignalBeta;//!vs beta
  TH2F* fTOFsignalPiBeta;//!vs beta
  TH2F* fTOFsignalKBeta;//!vs beta
  TH2F* fTOFsignalPBeta;//!vs beta
  TH2F* fPvsPt; //!P vs Pt yield
  TProfile* fMeanPvsP; //!mean p per bin
  TProfile* fMeanPtvsPt; //!mean pt per bin
  TList* fOutputList;//!output list
	
  AliAnalysisTaskPIDflowQA(const  AliAnalysisTaskPIDflowQA&); // not implemented
  AliAnalysisTaskPIDflowQA& operator=(const  AliAnalysisTaskPIDflowQA&); // not implemented
 
  void pidITS(AliESDtrack* t, Int_t pdgcode);
  void pidTPC(AliESDtrack* t, Int_t pdgcode);
  void pidTOF(AliESDtrack* t, Int_t pdgcode);
    
  ClassDef( AliAnalysisTaskPIDflowQA, 2); // example of analysis
};

#endif
