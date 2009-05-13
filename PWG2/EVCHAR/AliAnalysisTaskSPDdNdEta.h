#ifndef ALIANALYSISTASKSPDDNDETA_H
#define ALIANALYSISTASKSPDDNDETA_H

///////////////////////////////////////////////////////////////////////////
// Class AliAnalysisTaskSPDdNdEta                                        //
// Analysis task for dN/dEta reconstruction with the SPD                 //
//                                                                       //
// Author:  M. Nicassio (INFN Bari)                                      //
// Contact: Maria.Nicassio@ba.infn.it, Domenico.Elia@ba.infn.it          //
///////////////////////////////////////////////////////////////////////////

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
  void SetppAnalysis(Bool_t ppAna) { fppAna = ppAna; }
 
 protected:
  AliESDEvent *fESD;               //ESD object
  TList* fOutput;                  //! list send on output slot 0
  
  Bool_t fCorr;                    // flag to anable the correction histo calculation
  Int_t fTrigger;                  // to set the MBtrigger selection
  Bool_t fppAna;                   // flag to set the proper multiplicity binning (p-p or Pb-Pb)
 
  TH2F        *fHistSPDRAWMultvsZ;          // data to be corrected 
  TH2F        *fHistSPDRAWMultvsZTriggEvts; // data to be corrected
  TH2F        *fHistSPDRAWEtavsZ;           // data to be corrected

  TH1F        *fHistSPDmultEtacut;          // cluster inner layer and tracklet check histos
  TH1F        *fHistSPDmult;                // cluster inner layer and tracklet check histos
  TH1F        *fHistSPDeta;                 // cluster inner layer and tracklet check histos
  TH1F        *fHistSPDcl1multEtacutLay1;   // cluster inner layer and tracklet check histos
  TH1F        *fHistSPDcl1mult;             // cluster inner layer and tracklet check histos
  TH1F        *fHistSPDcl1eta;              // cluster inner layer and tracklet check histos
  TH1F        *fHistSPDphi;                 // cluster inner layer and tracklet check histos
  TH1F        *fHistSPDcl1phi;              // cluster inner layer and tracklet check histos
  TH1F        *fHistSPDtheta;               // cluster inner layer and tracklet check histos
  TH1F        *fHistSPDcl1theta;            // cluster inner layer and tracklet check histos
  TH1F        *fHistSPDdePhi;               // cluster inner layer and tracklet check histos
  TH1F        *fHistSPDdePhiZ;              // cluster inner layer and tracklet check histos
  TH1F        *fHistSPDdePhi3D;             // cluster inner layer and tracklet check histos
  TH2F        *fHistSPDphivsSPDeta;         // cluster inner layer and tracklet check histos
  TH2F        *fHistSPDcl1phivsSPDcl1eta;   // cluster inner layer and tracklet check histos

  TH1F        *fHistSPDvtx;                 // SPD vertex distributions
  TH3F        *fHistSPDvtx3D;               // SPD vertex distributions
  TH3F        *fHistSPDvtxZ;                // SPD vertex distributions
  TH2F        *fHistNcontribSPDvtxvsSPDvtx; // SPD vertex distributions
  TH1F        *fHistNcontribSPDvtx3D;       // SPD vertex distributions
  TH1F        *fHistNcontribSPDvtxZ;        // SPD vertex distributions
  TH1F        *fHistNcontribSPDvtxall;      // SPD vertex distributions
  TH2F        *fHistSPDmultvsSPDvtx;        // SPD vertex distributions

  TH2F        *fHistSPDcl1multvsnFiredChipsLay1;              // SPD fired chips distributions
  TH2F        *fHistSPDmultvsnFiredChipsLay1;                 // SPD fired chips distributions
  TH2F        *fHistSPDmultvsnFiredChipsLay2;                 // SPD fired chips distributions
  TH2F        *fHistnFiredChipsLay2vsnFiredChipsLay1;         // SPD fired chips distributions
  TH2F        *fHistnFiredChipsLay2vsnFiredChipsLay1novtxrec; // SPD fired chips distributions

  TH2F* fHistBkgCorrNum;           // track level correction histograms
  TH2F* fHistBkgCorrDen;           // track level correction histograms
  TH2F* fHistAlgEffNum;            // track level correction histograms
  TH2F* fHistNonDetectableCorrNum; // track level correction histograms
  TH2F* fHistNonDetectableCorrDen; // track level correction histograms
  TH2F* fHistTrackTrigVtxCorrNum;  // track level correction histograms
  TH2F* fHistTrackTrigCorrDen;     // track level correction histograms
   
  TH2F* fHistTrigVtxCorrNum;       // event level correction histograms
  TH2F* fHistTrigVtxCorrDen;       // event level correction histograms

  TH2F* fHistTrigCorrDen;          // event level correction histograms

  TH2F* fHistMCEtavsZTriggMCvtxEvts;  // MC distributions
  TH2F* fHistMCEtavsZTriggESDvtxEvts; // MC distributions
  TH2F* fHistMCEtavsZ;                // MC distributions 

  TH1F* fHistTRradius;                // additional check histos
  TH2F* fHistContributorsvsDeVtx;     // additional check histos
  TH3F* fHistoDetectableNotr;         // additional check histos
  TH2F* fHistoDetectabletr;           // additional check histos
  TH2F* fHistoNonStoppingTracks;      // additional check histos
  TH2F* fHistoDetectedLay1;           // additional check histos
  TH2F* fHistoDetectedLay2;           // additional check histos
  TH1F* fHistoPt;                     // additional check histos 
  TH2F* fHistoDetectableTRm1;         // additional check histos
  TH2F* fHistoDetectableTrITS;        // additional check histos
  TH2F* fHistoDetectableTrTPC;        // additional check histos
  TH2F* fHistoDetectableTrFRAME;      // additional check histos
  TH2F* fHistoDetectableTrTRD;        // additional check histos
  TH2F* fHistoDetectableTrTOF;        // additional check histos
  TH2F* fHistoDetectableTrMUON;       // additional check histos
  TH2F* fHistoDetectableTrHMPID;      // additional check histos
  TH2F* fHistoDetectableTrT0;         // additional check histos
  TH2F* fHistoDetectableTrEMCAL;      // additional check histos
  TH2F* fHistoDetectableTrFMD;        // additional check histos
  TH2F* fHistoDetectableTrVZERO;      // additional check histos
  TH1F* fHistoRTRm1;                  // additional check histos

 private:    
  AliAnalysisTaskSPDdNdEta(const AliAnalysisTaskSPDdNdEta&); 
  AliAnalysisTaskSPDdNdEta& operator=(const AliAnalysisTaskSPDdNdEta&); 
  
  ClassDef(AliAnalysisTaskSPDdNdEta, 1); 
};

#endif
