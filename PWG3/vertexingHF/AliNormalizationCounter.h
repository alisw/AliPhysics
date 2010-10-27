#ifndef ALINORMALIZATIONCOUNTER_H
#define ALINORMALIZATIONCOUNTER_H

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH1D.h>
#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include <AliAODEvent.h>
#include <AliVParticle.h>
#include "AliAnalysisTaskSE.h"
#include "AliCounterCollection.h"
//#include "AliAnalysisVertexingHF.h"

class AliNormalizationCounter : public AliCounterCollection
{
 public:

  AliNormalizationCounter();
  AliNormalizationCounter(const char *name);
  virtual ~AliNormalizationCounter();

  void SetESD(Bool_t flag){fESD=flag;}
  void StoreEvent(AliVEvent*,Bool_t mc=kFALSE);
  void StoreCandidates(AliVEvent*, Int_t nCand=0,Bool_t flagFilter=kTRUE);
  TH1D* DrawAgainstRuns(TString candle);
  TH1D* DrawRatio(TString candle1,TString candle2);
 private:
  AliNormalizationCounter(const AliNormalizationCounter &source);
  AliNormalizationCounter& operator=(const AliNormalizationCounter& source);
  Bool_t fESD; //flag for ESD vs AOD

  ClassDef(AliNormalizationCounter,1);

};
#endif
