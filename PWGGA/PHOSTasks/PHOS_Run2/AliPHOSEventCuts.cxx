#include "stdio.h"
#include "iostream"
#include "TH2.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TString.h"

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliPHOSGeometry.h"
#include "AliOADBContainer.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliVHeader.h"
#include "AliVCaloTrigger.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliAnalysisUtils.h"

#include "AliPHOSTriggerHelper.h"
#include "AliPHOSEventCuts.h"

using namespace std;

// Author: Daiki Sekihata (Hiroshima University)

ClassImp(AliPHOSEventCuts)

//________________________________________________________________________
AliPHOSEventCuts::AliPHOSEventCuts(const char *name):
	fIsMC(kFALSE),
  fUsePHOSTender(kTRUE),
  fMaxAbsZvtx(10.),
  fRejectPileup(kTRUE),
  fRejectDAQIncomplete(kTRUE),
  fIsPHOSTriggerAnalysis(kFALSE),
  fTriggerHelper(0x0),
  fPHOSGeo(0x0)
{
  // Constructor

  for(Int_t i=0;i<6;i++){
    fPHOSTRUBadMap[i] = 0x0;
  }

  AliInfo("event selection constructor");
}
//________________________________________________________________________
AliPHOSEventCuts::~AliPHOSEventCuts()
{
  AliInfo("event selection destructor");



}
//________________________________________________________________________
Bool_t AliPHOSEventCuts::AcceptEvent(AliVEvent *event)
{
  AliAODEvent *fAODEvent = dynamic_cast<AliAODEvent*>(event);
  AliESDEvent *fESDEvent = dynamic_cast<AliESDEvent*>(event);

  //select event which PHOS was readout from trigger cluster point of view.
  //for example, PHOS was not in MUFAST cluster.
  TString trigClasses = event->GetFiredTriggerClasses();

  if(!fIsMC && !trigClasses.Contains("CENT") && !trigClasses.Contains("FAST")){
    //At least, PHOS must be in CENT or FAST as a readout cluster. INT7 or PHI7 do not matter.
    AliWarning(Form("Skip event with triggers %s",trigClasses.Data()));
    return kFALSE;
  }

  Int_t run = event->GetRunNumber();

  if(fUsePHOSTender){
    if(run<209122) //Run1
      fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP");
    else
      fPHOSGeo = AliPHOSGeometry::GetInstance("Run2");
  }
  else{
    if(run<209122) //Run1
      fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP");
    else
      fPHOSGeo = AliPHOSGeometry::GetInstance("Run2");

    AliOADBContainer geomContainer("phosGeo");
    geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
    TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(run,"PHOSRotationMatrixes");

    for(Int_t mod=0; mod<6; mod++) {
      if(!matrixes->At(mod)) {
        AliError(Form("No PHOS Matrix for mod:%d, geo=%p\n", mod, fPHOSGeo));
        continue;
      }
      else {
        fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
        AliInfo(Form("Adding PHOS Matrix for mod:%d, geo=%p\n", mod, fPHOSGeo));
      }
    }//end of module loop

  }

  Bool_t IsZvtxOut       = kFALSE;
  Bool_t IsDAQIncomplete = kFALSE;
  Bool_t eventPileup     = kFALSE;

  if(event->IsIncompleteDAQ()){
    AliWarning("This is IncompleteDAQ.");
    IsDAQIncomplete = kTRUE;
  }

  const AliVVertex *vVertex    = event->GetPrimaryVertex();
  //const AliVVertex *vVertexSPD = event->GetPrimaryVertexSPD();

  if(vVertex->GetNContributors()<1){
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

  if(IsZvtxOut)                               return kFALSE; //reject event with Zvtx > threshold
  if(fRejectPileup && eventPileup)            return kFALSE; //reject pile up event
  if(fRejectDAQIncomplete && IsDAQIncomplete) return kFALSE; //reject DAQ imcopmelete event

  Bool_t IsPHI7fired = kFALSE;
  //additional criteriat for only PHOS trigger analysis
  if(fIsPHOSTriggerAnalysis){
    IsPHI7fired = fTriggerHelper->IsPHI7(event);
    if(!IsPHI7fired) return kFALSE;
  }//end of PHOS trigger decision.


  return kTRUE;
}
//________________________________________________________________________

