#ifndef ALIANALYSISMUMUTRIGGERRESPONSE_H
#define ALIANALYSISMUMUTRIGGERRESPONSE_H

/**
 *
 * \class AliAnalysisMuMuTriggerResponse
 * \brief Trigger response syst. analysis
 * \author B. Audurier (Subatech)
 */

#include "AliAnalysisMuMuBase.h"
#include "AliAnalysisMuMuBinning.h"
#include "AliMergeableCollection.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "THnSparse.h"
#include "TH1.h"

class TH1;
class AliVParticle;
class TLorentzVector;
class AliMergeableCollectionProxy;

class AliAnalysisMuMuTriggerResponse : public AliAnalysisMuMuBase
{
public:

  AliAnalysisMuMuTriggerResponse(AliMergeableCollection* oc=0x0);
  virtual ~AliAnalysisMuMuTriggerResponse();

  void SetLptAptEtaRange(Double_t* x, int xsize );

protected:

  void DefineHistogramCollection(const char* eventSelection,
                                 const char* triggerClassName,
                                 const char* centrality,
                                 Bool_t mix =kFALSE);

  void FillHistosForTrack(
                          const char* eventSelection,
                          const char* triggerClassName,
                          const char* centrality,
                          const char* trackCutName,
                          const AliVParticle& track);

  virtual void FillHistosForPair(const char* eventSelection,const char* triggerClassName,
                                 const char* centrality,
                                 const char* pairCutName,
                                 const AliVParticle& part,
                                 const AliVParticle& part2,
                                 const Bool_t IsMixedHisto =kFALSE);

  Double_t WeightFromEta(const AliVParticle& track,
                                AliMergeableCollectionProxy& proxylptaptHisto,
                                int value);

  Double_t WeightFromLocalBoard(const AliVParticle& track,
                                Double_t LocalBoardNumber,
                                AliMergeableCollectionProxy& proxylptaptHisto,
                                int value);

  Double_t WeightFromLocalBoardGroup(const AliVParticle& track,
                                Double_t LocalBoardNumber,
                                AliMergeableCollectionProxy& proxylptaptHisto,
                                int value);


private:

  AliMergeableCollection* fOC; // mergeable collection from first run of the task
  Int_t fNEtaBin; // Number of Pt bin
  Double_t* fLptAptEtaRanges; // [fNEtaBin]


  ClassDef(AliAnalysisMuMuTriggerResponse,3) // implementation of AliAnalysisMuMuBase for trigger response studies
};

#endif
