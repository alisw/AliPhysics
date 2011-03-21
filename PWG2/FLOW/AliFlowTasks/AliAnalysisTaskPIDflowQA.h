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
  void SetAliESDtrackCuts(AliFlowTrackCuts* const cuts ){fCuts=cuts;}
  void SetEventCuts(AliFlowEventCuts* c) {fEventCuts=c;}
  void SetMCOn(){fMC=kTRUE;}
  AliESDpid* GetESDpid() const {return fESDpid;}
  void SetUseDebugFile(Bool_t b=kTRUE) {fUseDebugFile=b;}

  Float_t Beta(Float_t m, Float_t p);

private:
  AliESDEvent *fESD;            //!ESD object    
  AliFlowTrackCuts *fCuts;       //cuts 
  AliFlowEventCuts *fEventCuts; //event cuts
  AliESDpid *fESDpid;           //pid object
  Bool_t fMC;                   //if TRUE use MC 
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

  TH2F* fTOFbetaAfterElectronsCuts; //!
  TH2F* fTOFbetaAfterPionCuts; //!
  TH2F* fTOFbetaAfterKaonCuts; //!
  TH2F* fTOFbetaAfterProtonCuts; //!
  TH2F* fTPCsignalAfterPionCuts; //!
  TH2F* fTPCsignalAfterKaonCuts; //!
  TH2F* fTPCsignalAfterProtonCuts; //!

  TH2F* fTOFbetaAfterElectronsCuts1; //!
  TH2F* fTOFbetaAfterPionCuts1; //!
  TH2F* fTOFbetaAfterKaonCuts1; //!
  TH2F* fTOFbetaAfterProtonCuts1; //!
  TH2F* fTPCsignalAfterPionCuts1; //!
  TH2F* fTPCsignalAfterKaonCuts1; //!
  TH2F* fTPCsignalAfterProtonCuts1; //!

  TH2F* fTOFbetaEafter; //!
  TH2F* fTOFbetaPiafter; //!
  TH2F* fTOFbetaKafter; //!
  TH2F* fTOFbetaPafter; //!

  TH2F* fTPCsignalPiafter; //!
  TH2F* fTPCsignalKafter; //!
  TH2F* fTPCsignalPafter; //!

  //MC
  TH1F* fTOFyieldSelEmcE;//!
  TH1F* fTOFyieldSelPimcE;//!
  TH1F* fTOFyieldSelKmcE;//!
  TH1F* fTOFyieldSelPmcE;//!
  TH1F* fTOFyieldSelEmcM;//!
  TH1F* fTOFyieldSelPimcM;//!
  TH1F* fTOFyieldSelKmcM;//!
  TH1F* fTOFyieldSelPmcM;//!
  TH1F* fTOFyieldSelEmcPi;//!
  TH1F* fTOFyieldSelPimcPi;//!
  TH1F* fTOFyieldSelKmcPi;//!
  TH1F* fTOFyieldSelPmcPi;//!
  TH1F* fTOFyieldSelEmcK;//!
  TH1F* fTOFyieldSelPimcK;//!
  TH1F* fTOFyieldSelKmcK;//!
  TH1F* fTOFyieldSelPmcK;//!
  TH1F* fTOFyieldSelEmcP;//!
  TH1F* fTOFyieldSelPimcP;//!
  TH1F* fTOFyieldSelKmcP;//!
  TH1F* fTOFyieldSelPmcP;//!
  TH1F* fTPCyieldSelEmcE;//!
  TH1F* fTPCyieldSelPimcE;//!
  TH1F* fTPCyieldSelKmcE;//!
  TH1F* fTPCyieldSelPmcE;//!
  TH1F* fTPCyieldSelEmcM;//!
  TH1F* fTPCyieldSelPimcM;//!
  TH1F* fTPCyieldSelKmcM;//!
  TH1F* fTPCyieldSelPmcM;//!
  TH1F* fTPCyieldSelEmcPi;//!
  TH1F* fTPCyieldSelPimcPi;//!
  TH1F* fTPCyieldSelKmcPi;//!
  TH1F* fTPCyieldSelPmcPi;//!
  TH1F* fTPCyieldSelEmcK;//!
  TH1F* fTPCyieldSelPimcK;//!
  TH1F* fTPCyieldSelKmcK;//!
  TH1F* fTPCyieldSelPmcK;//!
  TH1F* fTPCyieldSelEmcP;//!
  TH1F* fTPCyieldSelPimcP;//!
  TH1F* fTPCyieldSelKmcP;//!
  TH1F* fTPCyieldSelPmcP;//!
  
  TH2F* fTPCdedxAfterTOFpidPions;
  TH2F* fTPCdedxAfterTOFpidKaons;
  TH2F* fTPCdedxAfterTOFpidProtons;

  TH2F* fPvsPt; //!P vs Pt yield
  TProfile* fMeanPvsP; //!mean p per bin
  TH2F* fTPCvsGlobalMult; //! correlation tpc only tracks vs global tracks
  AliFlowTrackCuts* fStandardGlobalCuts; //! cuts
  AliFlowTrackCuts* fStandardTPCCuts; //! cuts
  AliFlowTrackCuts* fCutsTOFElectrons; //!
  AliFlowTrackCuts* fCutsTOFPions; //!
  AliFlowTrackCuts* fCutsTOFKaons; //!
  AliFlowTrackCuts* fCutsTOFProtons; //!
  AliFlowTrackCuts* fCutsTPCElectrons; //!
  AliFlowTrackCuts* fCutsTPCPions; //!
  AliFlowTrackCuts* fCutsTPCKaons; //!
  AliFlowTrackCuts* fCutsTPCProtons; //!
  TList* fOutputList;//!output list
	
  AliAnalysisTaskPIDflowQA(const  AliAnalysisTaskPIDflowQA&); // not implemented
  AliAnalysisTaskPIDflowQA& operator=(const  AliAnalysisTaskPIDflowQA&); // not implemented
 
  void pidITS(AliESDtrack* t, Int_t pdgcode);
  void pidTPC(AliESDtrack* t, Int_t pdgcode);
  void pidTOF(AliESDtrack* t, Int_t pdgcode);
    
  ClassDef( AliAnalysisTaskPIDflowQA, 3); // example of analysis
};

#endif
