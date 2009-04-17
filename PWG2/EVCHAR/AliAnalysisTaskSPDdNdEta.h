#ifndef AliAnalysisTaskSPDdNdEta_cxx
#define AliAnalysisTaskSPDdNdEta_cxx

class TH1F;
class TH2F;
class TH3F;
class AliESDEvent;

#include "AliAnalysisTask.h"

class AliAnalysisTaskSPDdNdEta : public AliAnalysisTask {
 public:
  AliAnalysisTaskSPDdNdEta(const char *name = "AliAnalysisTaskSPDdNdEta");
  virtual ~AliAnalysisTaskSPDdNdEta(); 
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void SetReadMC(Bool_t readmc = kFALSE) { fCorr = readmc; }
  void SetTrigger(Int_t MBtrigger) { fTrigger = MBtrigger; } 
 
 protected:
  AliESDEvent *fESD;               //ESD object
  TList* fOutput;                  //! list send on output slot 0
  
  Bool_t fCorr;
  Int_t fTrigger;
 
  //Data to be corrected 
  TH2F        *fHistSPDRAWMultvsZ;
  TH2F        *fHistSPDRAWMultvsZTriggEvts;
  TH2F        *fHistSPDRAWEtavsZ;

  //Clusters inner layer and tracklets
  TH1F        *fHistSPDmultEtacut;
  TH1F        *fHistSPDmult;
  TH1F        *fHistSPDeta;
  TH1F        *fHistSPDcl1multEtacutLay1;
  TH1F        *fHistSPDcl1mult;
  TH1F        *fHistSPDcl1eta;
  TH1F        *fHistSPDphi;
  TH1F        *fHistSPDcl1phi;
  TH1F        *fHistSPDtheta;
  TH1F        *fHistSPDcl1theta;
  TH1F        *fHistSPDdePhi;
  TH1F        *fHistSPDdePhiZ;
  TH1F        *fHistSPDdePhi3D;
  TH2F        *fHistSPDphivsSPDeta;
  TH2F        *fHistSPDcl1phivsSPDcl1eta;

  //SPD vertex distributions
  TH1F        *fHistSPDvtx;    
  TH3F        *fHistSPDvtx3D;
  TH3F        *fHistSPDvtxZ;
  TH2F        *fHistNcontribSPDvtxvsSPDvtx;
  TH1F        *fHistNcontribSPDvtx3D;
  TH1F        *fHistNcontribSPDvtxZ;
  TH1F        *fHistNcontribSPDvtxall;
  TH2F        *fHistSPDmultvsSPDvtx;

  //SPD fired chips distributions            
  TH2F        *fHistSPDcl1multvsnFiredChipsLay1;
  TH2F        *fHistSPDmultvsnFiredChipsLay1;
  TH2F        *fHistSPDmultvsnFiredChipsLay2;
  TH2F        *fHistnFiredChipsLay2vsnFiredChipsLay1;
  TH2F        *fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec;

 
  //Track level correction histograms
  TH2F* fHistBkgCorrNum;
  TH2F* fHistBkgCorrDen;
  TH2F* fHistAlgEffNum;
  TH2F* fHistNonDetectableCorrNum;
  TH2F* fHistNonDetectableCorrDen;
  TH2F* fHistTrackTrigVtxCorrNum;
  TH2F* fHistTrackTrigCorrDen;
   
  //Event level correction histograms
  TH2F* fHistTrigVtxCorrNum;
  TH2F* fHistTrigVtxCorrDen;

  TH2F* fHistTrigCorrDen;

  //MC distributions
  TH2F* fHistMCEtavsZTriggMCvtxEvts;
  TH2F* fHistMCEtavsZTriggESDvtxEvts;
  TH2F* fHistMCEtavsZ;

  //Check histos
  TH1F* fHistTRradius;
  TH2F* fHistContributorsvsDeVtx;
  TH3F* fHistoDetectableNotr;
  TH2F* fHistoDetectabletr;
  TH2F* fHistoNonStoppingTracks;
  TH2F* fHistoDetectedLay1;
  TH2F* fHistoDetectedLay2;
  TH1F* fHistoPt;
  TH2F* fHistoDetectableTRm1;
  TH2F* fHistoDetectableTR0;
  TH2F* fHistoDetectableTR1;
  TH2F* fHistoDetectableTR2;
  TH2F* fHistoDetectableTR3;
  TH2F* fHistoDetectableTR4;
  TH2F* fHistoDetectableTR5;
  TH2F* fHistoDetectableTR6;
  TH1F* fHistoRTRm1;
 private:    
  AliAnalysisTaskSPDdNdEta(const AliAnalysisTaskSPDdNdEta&); 
  AliAnalysisTaskSPDdNdEta& operator=(const AliAnalysisTaskSPDdNdEta&); 
  
  ClassDef(AliAnalysisTaskSPDdNdEta, 1); 
};

#endif
