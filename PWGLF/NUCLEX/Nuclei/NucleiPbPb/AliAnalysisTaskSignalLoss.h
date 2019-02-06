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
  void SetYRange (float ymin, float ymax) { fRequireYmin = ymin; fRequireYmax = ymax; }

private:
  AliAnalysisTaskSignalLoss(const AliAnalysisTaskSignalLoss &source);
  AliAnalysisTaskSignalLoss &operator=(const AliAnalysisTaskSignalLoss &source);

  TList*    fOutputList; //! Output List

  TH1F*     fHistAccEvents; //! Number of accepted events
  TH1F*     fHistTrueINELgt0Events; //! Number of true INEL>0 events
  TH3F*     fHistGenPartsAcc[2]; //! Number of particles in accepted events
  TH3F*     fHistGenPartsINELgt0[2]; //! Number of particles in true INEL>0 events

  TArrayF   fPtBins;
  TArrayF   fCentBins;
  TArrayI   fPDGcodes;

  int       fNspecies;
  float     fRequireYmin;
  float     fRequireYmax;

  ClassDef(AliAnalysisTaskSignalLoss,1);

};

#endif
