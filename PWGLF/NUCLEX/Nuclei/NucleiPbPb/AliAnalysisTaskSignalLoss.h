/// \author Luca Barioglio <luca.barioglio@cern.ch>, University and INFN Torino
/// \date Mar 30, 2018

#ifndef __AliAnalysisTaskSignalLoss__
#define __AliAnalysisTaskSignalLoss__

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

#include <Rtypes.h>
#include <TList.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TArrayF.h>
#include <TArrayI.h>

class AliAnalysisTaskSignalLoss : public AliAnalysisTaskSE{

public:

  AliAnalysisTaskSignalLoss(const char* task_name = "SignalLossTask");
  virtual ~ AliAnalysisTaskSignalLoss();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);

  AliEventCuts  fEventCuts;

  void SetCentBins (int nbins, float *bins);
  void SetPtBins   (int nbins, float *bins);
  void SetPDGcodes (int nbins, int *bins);

private:

  TList*    fOutputList;

  TH1F*     fHistAccEvents;
  TH1F*     fHistTrueINELgt0Events;
  TH3F*     fHistGenPartsAcc[2];
  TH3F*     fHistGenPartsINELgt0[2];

  TArrayF   fPtBins;
  TArrayF   fCentBins;
  TArrayI   fPDGcodes;

  int       fNspecies;
  float     fRequireYmin;
  float     fRequireYmax;

  ClassDef(AliAnalysisTaskSignalLoss,1);

};

#endif
