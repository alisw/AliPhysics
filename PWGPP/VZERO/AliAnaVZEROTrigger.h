#ifndef AliAnaVZEROTrigger_cxx
#define AliAnaVZEROTrigger_cxx

class TH1F;
class TH2F;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"

class AliAnaVZEROTrigger : public AliAnalysisTaskSE {
 public:
  AliAnaVZEROTrigger();
  AliAnaVZEROTrigger(const char *name);
  virtual ~AliAnaVZEROTrigger() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  virtual void Init();

  void Setup(const char *filename);
  void SetMBTrigName(const char *name = "CPBI") {fMBTrigName = name;}
  void SetUsePhysSel(Bool_t usePhysSel) {fUsePhysSel = usePhysSel;}

  Float_t GetMinThr() const {return fMinThr;}
  Float_t GetMaxThr() const {return fMaxThr;}
  Float_t GetRatio() const {return fRatio;}
  Int_t GetNThr() const {return fNThr;}

  Float_t GetCentCuts(Int_t i) const {return fCentCuts[i];}
  Float_t GetSemiCentCuts(Int_t i) const {return fSemiCentCuts[i];}
  Float_t GetThrA(Int_t j) const;
  Float_t GetThrC(Int_t j) const;

 private:
  AliESDEvent *fESD;    //! ESD object
  TList       *fOutputList; //! Output list

  Float_t fMinThr; // Minimum threshold for the scan
  Float_t fMaxThr; // Maximum threshold for the scan
  Float_t fRatio;  // Ratio between C and A side trigger charge
  Int_t   fNThr;   // Number of thresholds in the scan

  Float_t fCentCuts[2]; // Central cut extracted from the eff fit
  Float_t fSemiCentCuts[2]; // Semi-central cut extracted from the eff fit

  TString fMBTrigName; // MB trigger name (for evt sel)
  Bool_t  fUsePhysSel; // Use or not phys sel

  TH1F *fV0Percent; //! V0 centrality percentile
  TH1F *fV0PercentAll; //! V0 centrality percentile (no evt sel)
  TH1F *fZvtx; //! Z vertex distribution after phys sel
  TH2F *fXYvtx; //! XY vertex distribution after phys sel
  TH1F *fV0Mult1d; //! V0 multiplicity distribution
  TH2F *fV0Charge2d; //! V0 trigger charges distribution
  TH2F *fV0Charge2dPercent; //! V0 trigger charges distribution weighted by centrality
  TH2F *fV0Charge2dAll; //! V0 trigger charges distribution (no evt sel)

  TH1F **fV0PercentBins; //! V0 centrality percentile with sequential thr
  TH1F **fV0PercentBinsAll; //! V0 centrality percentile with sequential thr (no evt sel)

  TH1F *fV0Cent; //! centrality percentile with central trigger using custom thresholds
  TH1F *fV0CentAll; //! centrality percentile with central trigger using custom threshold (no evt sel)
  TH1F *fV0SemiCent; //! centrality percentile with semi-central trigger using custom thresholds
  TH1F *fV0SemiCentAll; //! centrality percentile with semi-central trigger using custom thresholds (no evt sel)

  TH1F *fV0CentHw; //! centrality percentile with central trigger using hardware thresholds
  TH1F *fV0CentHwAll; //! centrality percentile with central trigger using hardware threshold (no evt sel)
  TH1F *fV0SemiCentHw; //! centrality percentile with semi-central trigger using hardware thresholds
  TH1F *fV0SemiCentHwAll; //! centrality percentile with semi-central trigger using hardware thresholds (no evt sel)

  TH1F *fV0CentTr; //! centrality percentile with central trigger using hardware thresholds
  TH1F *fV0CentTrAll; //! centrality percentile with central trigger using hardware threshold (no evt sel)
  TH1F *fV0SemiCentTr; //! centrality percentile with semi-central trigger using hardware thresholds
  TH1F *fV0SemiCentTrAll; //! centrality percentile with semi-central trigger using hardware thresholds (no evt sel)

  TH1F *fV0Percent63; //! centrality percentile when >=63 cells are fired
  TH1F *fV0Percent63All; //! centrality percentile when >=63 cells are fired (no evt sel)

  TH2F *fV0MultAll; //! total V0 multiplicity for all events
  TH2F *fV0Mult63; //! total V0 multiplicity for events with at least 63 bb flags

  TH2F *fV0CentVsMult; //! V0 centrality vs total V0 multiplicity
  TH2F *fV0CentVsCharge; //! V0 centrality vs total V0 offline charge
  TH2F *fV0CentVsTrCharge; //! V0 centrality vs total V0 online charge
  TH2F *fV0TrChargeVsChargeA; //! V0 nline charge (A side) vs V0 offline charge (A side)
  TH2F *fV0TrChargeVsChargeC; //! V0 nline charge (C side) vs V0 offline charge (C side)

  TH2F *fAdcPmt; //! Charge per channel
  TH2F *fAdcPmt0; //! Charge per channel
  TH2F *fAdcPmt1; //! Charge per channel
  TH2F *fAdcPmt2; //! Charge per channel
  TH2F *fAdcPmt3; //! Charge per channel

  AliAnaVZEROTrigger(const AliAnaVZEROTrigger&); // not implemented
  AliAnaVZEROTrigger& operator=(const AliAnaVZEROTrigger&); // not implemented
  
  ClassDef(AliAnaVZEROTrigger, 1); // VZERO analysis task for setting up of centrality trigger
};

#endif
