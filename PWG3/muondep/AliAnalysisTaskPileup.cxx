/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-------------------------------------------------------------------------
//             Implementation of class AliAnalysisTaskPileup
//
// The class re-weights the number of analysed events with the correction
// factor to account for the event pileup.
//
// The final idea will be to include all of the methods for the pileup
// calculation.
//
// For the moment, the implemented method is based on a work 
// by L. Aphecetche.
// The probability of having at least 1 collision is given by:
// MB L0 trigger / max collision rate
// The max collision rate is given by:
//    number of bunches X revolution frequency ~ CBEAMB L0 trigger
// Assuming a poissonian distribution of the collision, with:
//    mu = mean number of collisions
// the probability of having 0 collisions is:
//    exp{-mu} = 1 - L0b_MB / max collision rate
// Since the MB trigger is measuring the probability of having at least
// 1 collisions, the number should be corrected by the factor:
//    CF = mu / ( 1 - exp{-mu} )
// The class weights the number of events with this correction factor
//
// CAVEATS:
// - so far the class needs to access the OCDB
//   (hopefully it will be changed in the future)
//   Hence the following libraries need to be included (in addition to
//   the standard analysis ones):
//   libRAWDatabase.so libCDB.so libSTEER.so libPWG3base.so
//-------------------------------------------------------------------------

#include <Riostream.h>

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TMath.h"

// STEER includes
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliTriggerScalersESD.h"
#include "AliLog.h"

// ANALYSIS includes
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"

#include "AliCounterCollection.h"

#include "AliAnalysisTaskPileup.h"

#if defined(READOCDB)
#include "AliTriggerRunScalers.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerClass.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBPath.h"
#endif


ClassImp(AliAnalysisTaskPileup)


//________________________________________________________________________
AliAnalysisTaskPileup::AliAnalysisTaskPileup(const char *name) :
  AliAnalysisTaskSE(name),
  fEventCounters(0x0),
  fHistoEventsList(0x0),
  fTriggerClasses(0x0),
  fTriggerClassIndex(new TArrayI(50))
#if defined(READOCDB)
  , fTriggerRunScalers(0x0),
  fDefaultStorage(new TString())
#endif

{
  /// Constructor

  DefineOutput(1,AliCounterCollection::Class());
  DefineOutput(2,TObjArray::Class());
}

//________________________________________________________________________
AliAnalysisTaskPileup::~AliAnalysisTaskPileup()
{
  /// Destructor
  delete fEventCounters;
  delete fHistoEventsList;
  delete fTriggerClasses;
  delete fTriggerClassIndex;

#if defined(READOCDB)
  delete fTriggerRunScalers;
  delete fDefaultStorage;
#endif

}


//___________________________________________________________________________
void AliAnalysisTaskPileup::NotifyRun()
{
  /// Notify run

#if defined(READOCDB)
  if ( ! AliCDBManager::Instance()->GetDefaultStorage() ) {
    AliCDBManager::Instance()->SetDefaultStorage(fDefaultStorage->Data());
  }

  AliCDBManager::Instance()->SetRun(InputEvent()->GetRunNumber());

  AliCDBEntry *entry = 0x0;
  entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
  if ( entry ) {
    AliTriggerConfiguration* trigConf = (AliTriggerConfiguration*)entry->GetObject();
    const TObjArray& classesArray = trigConf->GetClasses();
    Int_t nclasses = classesArray.GetEntriesFast();

    if ( fTriggerClasses ) delete fTriggerClasses;

    fTriggerClasses = new TObjArray(nclasses);
    fTriggerClasses->SetOwner();

    Int_t currActive = -1;
    for( Int_t iclass=0; iclass < nclasses; iclass++ ) {
      AliTriggerClass* trclass = (AliTriggerClass*)classesArray.At(iclass);
      if (trclass && trclass->GetMask()>0) {
	currActive++;
	Int_t currPos = currActive;

	// Store the CBEAMB class at the first position
	TString trigName = trclass->GetName();
	if ( trigName.Contains("CBEAMB") && ! trigName.Contains("WU") )
	  currPos = 0;
	else if ( ! fTriggerClasses->At(0) )
	  currPos++;

	Int_t trindex = TMath::Nint(TMath::Log2(trclass->GetMask()));
	TObjString* objString = new TObjString(trigName.Data());
	fTriggerClasses->AddAtAndExpand(objString, currPos);
	(*fTriggerClassIndex)[currPos] = trindex;
	AliDebug(3, Form("Current class %s  index %i  position %i", trigName.Data(), trindex, currPos));
      } // is class active
    } // loop on trigger classes
  } // if entry

  entry = AliCDBManager::Instance()->Get("GRP/CTP/Scalers");
  if (entry) {    
       AliInfo("Found an AliTriggerRunScalers in GRP/CTP/Scalers, reading it");
       fTriggerRunScalers = dynamic_cast<AliTriggerRunScalers*> (entry->GetObject());
       entry->SetOwner(0);
       if (fTriggerRunScalers->CorrectScalersOverflow() == 0) AliInfo("32bit Trigger counters corrected for overflow");
  }
#endif

}


//___________________________________________________________________________
void AliAnalysisTaskPileup::UserCreateOutputObjects()
{
  /// Create histograms and counters
  
  // The framework has problems if the name of the object
  // and the one of container differ
  // To be complaint, get the name from container and set it
  TString containerName = GetOutputSlot(1)->GetContainer()->GetName();

  // initialize event counters
  fEventCounters = new AliCounterCollection(containerName.Data());
  fEventCounters->AddRubric("event", "any/correctedL0");
  fEventCounters->AddRubric("trigger", 1000000);
  fEventCounters->AddRubric("selection", "any/hasVtxContrib/nonPileupSPD");
  fEventCounters->Init(kTRUE);

  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fEventCounters);
}

//________________________________________________________________________
void AliAnalysisTaskPileup::UserExec(Option_t *)
{
  /// Called for each event

  AliAODEvent* aodEvent = 0x0;

  AliESDEvent* esdEvent = dynamic_cast<AliESDEvent*> (InputEvent());
  if ( ! esdEvent ) {
    aodEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  }

  if ( ! aodEvent && ! esdEvent ) {
    AliError ("AOD or ESD event not found. Nothing done!");
    return;
  }
 
  // check physics selection
  //Bool_t isPhysicsSelected = (fInputHandler && fInputHandler->IsEventSelected());

#if defined(READOCDB)

  TString firedTrigClasses = ( esdEvent ) ? esdEvent->GetFiredTriggerClasses() : aodEvent->GetFiredTriggerClasses();

  if ( firedTrigClasses.IsNull() ) return; // Reject non-physics events

  AliDebug(2, Form("\nEvent %lli\n", Entry()));

  Int_t nPoints = fTriggerRunScalers->GetScalersRecordsESD()->GetEntriesFast();

  AliTimeStamp timeStamp(InputEvent()->GetOrbitNumber(), InputEvent()->GetPeriodNumber(), InputEvent()->GetBunchCrossNumber());
  Int_t position = fTriggerRunScalers->FindNearestScalersRecord(&timeStamp);
  if ( position < 0 ) {
    AliWarning("Position out of range: put to 1");
    position = 1;
  } 
  if ( position == 0 ) position++; // Don't trust the first one
  else if ( position + 1 >= nPoints ) position--;
  AliDebug(2, Form("position %i\n", position));
  AliTriggerScalersRecordESD* trigScalerRecords1 = (AliTriggerScalersRecordESD*)fTriggerRunScalers->GetScalersRecordsESD()->At(position);
  AliTriggerScalersRecordESD* trigScalerRecords2 = 0x0;

  // Sometimes scalers are filled very close to each others
  // in this case skip to the next entry
  for ( Int_t ipos=position+1; ipos<nPoints; ipos++ ) {
    trigScalerRecords2 = (AliTriggerScalersRecordESD*)fTriggerRunScalers->GetScalersRecordsESD()->At(ipos);
    Double_t deltaTime = (Double_t)( trigScalerRecords2->GetTimeStamp()->GetSeconds() - trigScalerRecords1->GetTimeStamp()->GetSeconds() );
    AliDebug(2, Form("Pos %i  TimeStamp %u - %u = %.0f\n", ipos, trigScalerRecords2->GetTimeStamp()->GetSeconds(), trigScalerRecords1->GetTimeStamp()->GetSeconds(), deltaTime));

    if ( deltaTime > 1 )
      break;
  }

#endif

  ULong64_t trigMask = 0;
  Double_t correctFactor = 1.;

  Int_t nVtxContrib = ( esdEvent ) ? esdEvent->GetPrimaryVertex()->GetNContributors() : aodEvent->GetPrimaryVertex()->GetNContributors();
  Bool_t isPileupSPD = ( esdEvent ) ? esdEvent->IsPileupFromSPD(3,0.8) : aodEvent->IsPileupFromSPD(3,0.8);

  TString selKey[3] = {"any","hasVtxContrib","nonPileupSPD"};
  Bool_t fillSel[3] = {kTRUE, ( nVtxContrib > 0 ), ( ( nVtxContrib > 0 ) && ( ! isPileupSPD ) )};

  //const AliTriggerScalersRecordESD* trigScalerRecords = esdEvent->GetHeader()->GetTriggerScalersRecord(); // REMEMBER TO CUT

  TString trigName = "";
  TString eventType = "";

  Int_t nTriggerClasses = fTriggerClasses->GetEntries();
  Int_t classIndex = -1;
  Double_t deltaScalersBeam = 0., deltaScalers = 0.;
  Bool_t isFiredOnce = kFALSE;
  for (Int_t itrig=0; itrig<nTriggerClasses+1; itrig++) {

    Double_t correctFactorL0 = 1.;

    Bool_t isClassFired = kFALSE;

    if ( itrig < nTriggerClasses ) {

      // Check if current mask contains trigger
      trigName = ((TObjString*)fTriggerClasses->At(itrig))->GetString();
      classIndex = (*fTriggerClassIndex)[itrig];
      trigMask = ( 1ull << classIndex );
      isClassFired = ( trigMask & InputEvent()->GetTriggerMask() );

      if ( isClassFired || itrig == 0 ) {
	// Get scalers
	const AliTriggerScalersESD* scaler1 = trigScalerRecords1->GetTriggerScalersForClass(classIndex+1);
	const AliTriggerScalersESD* scaler2 = trigScalerRecords2->GetTriggerScalersForClass(classIndex+1);
	deltaScalers = scaler2->GetLOCB() - scaler1->GetLOCB();
	
	if ( itrig == 0 )
	  deltaScalersBeam = deltaScalers;
	else if ( isClassFired ) {
	  correctFactorL0 = GetL0Correction(deltaScalers, deltaScalersBeam);
	  AliDebug(2, Form("Scalers: %s %.0f  %s %.0f -> CF %f\n", fTriggerClasses->At(itrig)->GetName(), deltaScalers, fTriggerClasses->At(0)->GetName(), deltaScalersBeam, correctFactorL0));
	}
      }
    }
    else {
      trigName = "any";
      classIndex = -1;
      isClassFired = isFiredOnce;
    }

    if ( ! isClassFired ) continue;
    isFiredOnce = kTRUE;

    //const AliTriggerScalersESD* trigScaler = trigScalerRecords->GetTriggerScalersForClass(classIndex+1); // REMEMBER TO CUT
    //if ( classIndex > 1 ) printf("Index: trigger %i  scaler %i\n", classIndex+1, trigScaler->GetClassIndex()); // REMEMBER TO CUT

    AliDebug(2, Form("Fired trig %s\n", trigName.Data()));

    for ( Int_t ievType=0; ievType<2; ievType++ ){
      switch ( ievType ) {
#if defined(READOCDB)
      case kHeventsCorrectL0:
	correctFactor = correctFactorL0;
	eventType = "correctedL0";
	break;
#endif
      default:
	correctFactor = 1.;
	eventType = "any";
      }

      for ( Int_t isel=0; isel<3; isel++ ) {
	if ( ! fillSel[isel] ) continue;
	fEventCounters->Count(Form("event:%s/trigger:%s/selection:%s",eventType.Data(),trigName.Data(), selKey[isel].Data()),correctFactor);
      } // loop on vertex selection
    } // loop on event type
  } // loop on trigger classes

  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fEventCounters);
}


//________________________________________________________________________
void AliAnalysisTaskPileup::Terminate(Option_t *)
{
  //
  /// Save final histograms
  /// and draw result to the screen
  //

  fEventCounters = dynamic_cast<AliCounterCollection*>(GetOutputData(1));
  if ( ! fEventCounters ) return;

  if ( ! fHistoEventsList ) fHistoEventsList = new TObjArray(2);
  fHistoEventsList->SetOwner();

  TH2D* histo = 0x0;
  histo = fEventCounters->Draw("trigger","selection","event:any");
  if ( histo ) {
    histo->SetName("hEvents");
    histo->SetTitle("Events per trigger");
    fHistoEventsList->AddAtAndExpand(histo, kHevents);
  }

  histo = fEventCounters->Draw("trigger","selection","event:correctedL0");
  if ( histo ) {
    histo->SetName("hEventsCorrectL0");
    histo->SetTitle("L0 corrected events");
    fHistoEventsList->AddAtAndExpand(histo, kHeventsCorrectL0);
  }

  PostData(2, fHistoEventsList);

  if ( gROOT->IsBatch() )
    return;

  TH2D* histoPileupL0 = (TH2D*)fHistoEventsList->At(kHeventsCorrectL0)->Clone("hPileupL0");
  histoPileupL0->Divide((TH2D*)fHistoEventsList->At(kHevents));
  
  TCanvas *can = new TCanvas("can1_Pileup","Pileup",10,10,310,310);
  can->SetFillColor(10); can->SetHighLightColor(10);
  can->SetLeftMargin(0.15); can->SetBottomMargin(0.15);
  histoPileupL0->DrawCopy("text");
}


//________________________________________________________________________
Double_t AliAnalysisTaskPileup::GetL0Correction(Double_t nCINT1B, Double_t nCBEAMB)
{
  /// Get the correction factor for L0 calculation

  if ( nCBEAMB == 0. )
    return 1.;

  Double_t ratio = nCINT1B / nCBEAMB;

  if ( ratio >= 1. || ratio == 0. )
    return 1.;

  Double_t mu = -TMath::Log(1-ratio);

  return mu / ( 1. - TMath::Exp(-mu) );

}
 
