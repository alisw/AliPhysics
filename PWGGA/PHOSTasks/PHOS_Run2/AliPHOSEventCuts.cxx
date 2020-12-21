#include "TMath.h"
#include "TString.h"

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVHeader.h"
#include "AliAnalysisUtils.h"

#include "AliPHOSEventCuts.h"

using namespace std;

// Author: Daiki Sekihata (Hiroshima University)

ClassImp(AliPHOSEventCuts)

//________________________________________________________________________
AliPHOSEventCuts::AliPHOSEventCuts(const char *name):
	fIsMC(kFALSE),
  fMaxAbsZvtx(10.),
  fRejectPileup(kTRUE),
  fRejectDAQIncomplete(kTRUE),
  fPF(AliPHOSEventCuts::kSPDInMultBins)
{
  // Constructor

}
//________________________________________________________________________
AliPHOSEventCuts::~AliPHOSEventCuts()
{



}
//________________________________________________________________________
Bool_t AliPHOSEventCuts::AcceptEvent(AliVEvent *event)
{

//  //select event which PHOS was readout from trigger cluster point of view.
//  //for example, PHOS was not in MUFAST cluster.
//  TString trigClasses = event->GetFiredTriggerClasses();
//
//  if(!fIsMC 
//      && !trigClasses.Contains("-CENT") //accept CENT, CENTNO[TRD|PMD]
//      && !trigClasses.Contains("-FAST") //accept FAST
//      && !trigClasses.Contains("-CALO") //accept CALO, CALOFAST
//    ){
//    //At least, PHOS must be in CENT[|NOTRD|NOPMD] or [CALO|FAST] or as a readout cluster. INT7 or PHI7 do not matter.
//    AliWarning(Form("Skip event with triggers %s",trigClasses.Data()));
//    return kFALSE;
//  }

  Int_t run = event->GetRunNumber();

  Bool_t IsZvtxOut       = kFALSE;
  Bool_t IsDAQIncomplete = kFALSE;
  Bool_t eventPileup     = kFALSE;

  if(event->IsIncompleteDAQ()){
    AliWarning("This is IncompleteDAQ.");
    IsDAQIncomplete = kTRUE;
  }

  const AliVVertex *vVertex = event->GetPrimaryVertex();

  if(vVertex->GetNContributors() < 1){
    AliInfo("vertex N contributors is less than 1. reject.");
    return kFALSE;
  }

  Double_t vertex[3] = {};
  vertex[0] = vVertex->GetX();
  vertex[1] = vVertex->GetY();
  vertex[2] = vVertex->GetZ();

  if(TMath::Abs(vertex[2]) > fMaxAbsZvtx){
    AliInfo(Form("Zvtx %f cm is out of threshold %f cm.", vertex[2], fMaxAbsZvtx));
    IsZvtxOut = kTRUE;
  }

  if(244917 <= run && run <= 246994){
    const AliVVertex* vtTrc = event->GetPrimaryVertex();
    const AliVVertex* vtSPD = event->GetPrimaryVertexSPD();
    double covTrc[6],covSPD[6];
    vtTrc->GetCovarianceMatrix(covTrc);
    vtSPD->GetCovarianceMatrix(covSPD);
    double dz = vtTrc->GetZ()-vtSPD->GetZ();
    double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
    double errTrc = TMath::Sqrt(covTrc[5]);
    double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
    if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20){
      // reject, bad reconstructed track vertex
      AliInfo("reject, bad reconstructed track vertex");
      return kFALSE;
    }
  }

  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(event);
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(event);

//  if(fESDEvent){
//    if(fESDEvent->IsPileupFromSPD()) {
//      eventPileup = kTRUE;
//      AliInfo("This is pile up event.");
//    }
//  }//end of ESD pile up
//  else if(fAODEvent){
//    if(fAODEvent->IsPileupFromSPD()){
//      eventPileup = kTRUE;
//      AliInfo("This is pile up event.");
//    }
//  }//end of AOD pileup

  const Int_t minContributors=5;
  const Float_t minChi2=5.;
  const Float_t minWeiZDiff=15;
  const Bool_t checkPlpFromDifferentBC=kFALSE;

  AliAnalysisUtils utils;
  utils.SetMinPlpContribMV(minContributors);
  utils.SetMaxPlpChi2MV(minChi2);
  utils.SetMinWDistMV(minWeiZDiff);
  utils.SetCheckPlpFromDifferentBCMV(checkPlpFromDifferentBC);
  eventPileup = utils.IsPileUpMV(event);

  switch(fPF){
    case AliPHOSEventCuts::kSPD:
      if(esd)      eventPileup = esd->IsPileupFromSPD();
      else if(aod) eventPileup = aod->IsPileupFromSPD();
      break;

    case AliPHOSEventCuts::kSPDInMultBins:
      if(esd)      eventPileup = esd->IsPileupFromSPDInMultBins();
      else if(aod) eventPileup = aod->IsPileupFromSPDInMultBins();
      break;

    case AliPHOSEventCuts::kMultiVertexer:
      eventPileup = utils.IsPileUpMV(event);
      break;
    default:
      eventPileup = kFALSE;
      break;
  }

  if(IsZvtxOut)                               return kFALSE; //reject event with Zvtx > threshold
  if(fRejectPileup && eventPileup)            return kFALSE; //reject pile up event
  if(fRejectDAQIncomplete && IsDAQIncomplete) return kFALSE; //reject DAQ imcopmelete event

  return kTRUE;
}
//________________________________________________________________________

