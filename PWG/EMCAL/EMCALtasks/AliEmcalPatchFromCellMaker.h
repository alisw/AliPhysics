#ifndef ALIEMCALPATCHFROMCELLMAKER_H
#define ALIEMCALPATCHFROMCELLMAKER_H

class TClonesArray;
class TH3F;
class AliEMCALTriggerBitConfig;

#include "AliAnalysisTaskEmcal.h"

class AliEmcalPatchFromCellMaker : public AliAnalysisTaskEmcal {
 public:
  AliEmcalPatchFromCellMaker();
  AliEmcalPatchFromCellMaker(const char *name);
  virtual ~AliEmcalPatchFromCellMaker();

  void               SetCaloTriggersOutName(const char *name)          { fCaloTriggersOutName = name; }
  void               SetPatchDimension(Int_t i)                        { fPatchDim            = i;    }
  void               SetMinCellE(Double_t e)                           { fMinCellE            = e;    }
  void               SetCellTimeCuts(Double_t min, Double_t max)       { fCellTimeMin = min; fCellTimeMax = max; }
  void               ActivateSlidingPatch(Bool_t b)                    { fL1Slide             = b;    }

 protected:
  enum{
    kPatchCols = 48,
    kPatchRows = 64
  };

  void               ExecOnce();
  Bool_t             Run();
  // Bool_t             FillHistograms();
  void               UserCreateOutputObjects();

  Bool_t             FillPatchADCSimple();
  void               RunSimpleOfflineTrigger();

  //Getters
  Int_t              GetPatchDimension() const                         { return fPatchDim;  }
  Double_t           GetPatchArea() const                              { return (Double_t)(fPatchDim*fPatchDim)*0.014*0.014; }
  Int_t              GetDimFastor() const;
  Int_t              GetSlidingStepSizeFastor() const;

  TString            fCaloTriggersOutName;  // name of output patch array
  TClonesArray      *fCaloTriggersOut;      //!trigger array out
      
  Double_t           fPatchADCSimple[kPatchCols][kPatchRows];   // patch map for simple offline trigger
  Double_t           fPatchESimple[kPatchCols][kPatchRows];     // patch map for simple offline trigger

  Int_t              fPatchDim;             // dimension of patch in #cells
  Double_t           fMinCellE;             // minimum cell energy
  Double_t           fCellTimeMin;          // minimum time cell
  Double_t           fCellTimeMax;          // maximum time cell
  Bool_t             fL1Slide;              // sliding window on
  AliEMCALTriggerBitConfig *fTriggerBitConfig; // dummy trigger bit config

 private:
  TH3F     *fh3EEtaPhiCell;                    //! cell E, eta, phi
  TH2F     *fh2CellEnergyVsTime;               //! emcal cell energy vs time
  TH1F     *fh1CellEnergySum;                  //! sum of energy in all emcal cells

  AliEmcalPatchFromCellMaker(const AliEmcalPatchFromCellMaker&);            // not implemented
  AliEmcalPatchFromCellMaker &operator=(const AliEmcalPatchFromCellMaker&); // not implemented

  ClassDef(AliEmcalPatchFromCellMaker, 1); // Task to make PicoTracks in a grid corresponding to EMCAL/DCAL acceptance
};
#endif
