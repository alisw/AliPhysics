#ifndef ALIANALYSISTASKEMCALBADCELLS_H
#define ALIANALYSISTASKEMCALBADCELLS_H

class TH1;
class TH2;

#include "AliAnalysisTaskEmcal.h"

class AliAnalysisTaskEmcalBadCells : public AliAnalysisTaskEmcal {
 public:

  AliAnalysisTaskEmcalBadCells();
  AliAnalysisTaskEmcalBadCells(const char *name, Bool_t histo = kFALSE);
  virtual ~AliAnalysisTaskEmcalBadCells();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

 protected:
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;

  // General histograms
  TH2F                       *fh2AmplitudeCellNumber; //! amplitude vs cell number

 private:
  AliAnalysisTaskEmcalBadCells(const AliAnalysisTaskEmcalBadCells&);            // not implemented
  AliAnalysisTaskEmcalBadCells &operator=(const AliAnalysisTaskEmcalBadCells&); // not implemented

  ClassDef(AliAnalysisTaskEmcalBadCells, 2) // jet sample analysis task
};
#endif
