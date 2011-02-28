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

/* $Id$ */

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
#include "AliCentrality.h"

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
AliAnalysisTaskPileup::AliAnalysisTaskPileup() :
AliAnalysisTaskSE(),
  fEventCounters(0x0),
  fHistoEventsList(0x0),
  fTriggerClasses(0x0),
  fTriggerClassIndex(0x0),
  fIsInitCDB(0),
  fCentralityClasses(0x0)
#if defined(READOCDB)
  , fTriggerRunScalers(0x0),
  fStorageList(0x0)
#endif

{
  /// Default Constructor
}


//________________________________________________________________________
AliAnalysisTaskPileup::AliAnalysisTaskPileup(const char *name) :
  AliAnalysisTaskSE(name),
  fEventCounters(0x0),
  fHistoEventsList(0x0),
  fTriggerClasses(0x0),
  fTriggerClassIndex(0x0),
  fIsInitCDB(0),
  fCentralityClasses(0x0)
#if defined(READOCDB)
  , fTriggerRunScalers(0x0),
  fStorageList(0x0)
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

  // For proof: do not delete output containers
  if ( ! AliAnalysisManager::GetAnalysisManager()->IsProofMode() ) {
    delete fEventCounters;
  }

  delete fHistoEventsList;
  delete fTriggerClasses;
  delete fTriggerClassIndex;

#if defined(READOCDB)
  //delete fTriggerRunScalers; // Not owner -> Owned by OCDB
  delete fStorageList;
#endif

}


//___________________________________________________________________________
void AliAnalysisTaskPileup::NotifyRun()
{
  /// Notify run

#if defined(READOCDB)
  fStorageList->Compress();
  if ( ! AliCDBManager::Instance()->GetDefaultStorage() ) {
    for ( Int_t ientry=0; ientry<fStorageList->GetEntries(); ientry++ ) {
      TObjString* calibStr = (TObjString*)fStorageList->At(ientry);
      ientry++;
      TObjString* dbStr = (TObjString*)fStorageList->At(ientry);
      TString calibName = calibStr->GetString();
      if ( ! calibName.CompareTo("default") ) {
	AliCDBManager::Instance()->SetDefaultStorage(dbStr->GetName());
      }
      else {
	AliCDBManager::Instance()->SetSpecificStorage(calibStr->GetName(), dbStr->GetName());
      }
    }
  }

  // Default storage was not correclty set: nothing done
  if ( ! AliCDBManager::Instance()->GetDefaultStorage() ) return;

  AliCDBManager::Instance()->SetRun(InputEvent()->GetRunNumber());

  AliCDBEntry *entry = 0x0;
  entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
  if ( ! entry ) return;
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
      if ( ( trigName.Contains("CBEAMB") ||  trigName.Contains("CTRUE") ) && ! trigName.Contains("WU") )
	currPos = 0;
      else if ( ! fTriggerClasses->At(0) )
	currPos++;

      Int_t trindex = TMath::Nint(TMath::Log2(trclass->GetMask()));
      TObjString* objString = new TObjString(trigName.Data());
      fTriggerClasses->AddAtAndExpand(objString, currPos);
      (*fTriggerClassIndex)[currPos] = trindex;
      if ( fDebug >= 3 ) printf("AliAnalysisTaskPileup: Current class %s  index %i  position %i\n", trigName.Data(), trindex, currPos);
    } // is class active
  } // loop on trigger classes

  entry = AliCDBManager::Instance()->Get("GRP/CTP/Scalers");
  if ( ! entry ) return;
  AliInfo("Found an AliTriggerRunScalers in GRP/CTP/Scalers, reading it");
  fTriggerRunScalers = static_cast<AliTriggerRunScalers*> (entry->GetObject());
  if (fTriggerRunScalers->CorrectScalersOverflow() == 0) AliInfo("32bit Trigger counters corrected for overflow");

  fIsInitCDB = kTRUE;

#endif

}


//___________________________________________________________________________
void AliAnalysisTaskPileup::UserCreateOutputObjects()
{
  /// Create histograms and counters

  fTriggerClassIndex = new TArrayI(50);
  fTriggerClassIndex->Reset(-1);

  fCentralityClasses = new TAxis(20, 0., 100.);

  TString centralityClassesStr = "", currClass = "";
  for ( Int_t ibin=1; ibin<=fCentralityClasses->GetNbins(); ibin++ ){
    if ( ! centralityClassesStr.IsNull() )
      centralityClassesStr.Append("/");
    currClass = Form("%.0f-%.0f",fCentralityClasses->GetBinLowEdge(ibin),fCentralityClasses->GetBinUpEdge(ibin));
    centralityClassesStr += currClass;
    fCentralityClasses->SetBinLabel(ibin, currClass.Data());
  }
  
  // The framework has problems if the name of the object
  // and the one of container differ
  // To be complaint, get the name from container and set it
  TString containerName = GetOutputSlot(1)->GetContainer()->GetName();

  // initialize event counters
  fEventCounters = new AliCounterCollection(containerName.Data());
  fEventCounters->AddRubric("event", "any/correctedL0");
  fEventCounters->AddRubric("trigger", 1000000);
  fEventCounters->AddRubric("vtxSelection", "any/hasVtxContrib/nonPileupSPD");
  fEventCounters->AddRubric("selected", "yes/no");
  fEventCounters->AddRubric("centrality", centralityClassesStr.Data());
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
  Bool_t isPhysicsSelected = (fInputHandler && fInputHandler->IsEventSelected());
  TString selected = ( isPhysicsSelected ) ? "yes" : "no";

  TString firedTrigClasses = ( esdEvent ) ? esdEvent->GetFiredTriggerClasses() : aodEvent->GetFiredTriggerClasses();

  if ( ! fIsInitCDB ) {
    delete fTriggerClasses;
    fTriggerClasses = firedTrigClasses.Tokenize(" ");
    fTriggerClasses->SetOwner();
  }

#if defined(READOCDB)

  //if ( firedTrigClasses.IsNull() ) return; // Reject non-physics events

  if ( fDebug >= 2 ) printf("\nAliAnalysisTaskPileup: Event %lli\n", Entry());

  // If scalers is not correctly loaded, the task becomes a mere event counter
  Int_t nPoints = ( fTriggerRunScalers ) ? fTriggerRunScalers->GetScalersRecordsESD()->GetEntriesFast() : 0;

  AliTriggerScalersRecordESD* trigScalerRecords1 = 0x0;
  AliTriggerScalersRecordESD* trigScalerRecords2 = 0x0;
  if ( nPoints > 1 ) {
    // Add protection for MC (no scalers there!)

    AliTimeStamp timeStamp(InputEvent()->GetOrbitNumber(), InputEvent()->GetPeriodNumber(), InputEvent()->GetBunchCrossNumber());
    Int_t position = fTriggerRunScalers->FindNearestScalersRecord(&timeStamp);
    if ( position < 0 ) {
      AliWarning("Position out of range: put to 1");
      position = 1;
    } 
    if ( position == 0 ) position++; // Don't trust the first one
    else if ( position + 1 >= nPoints ) position--;
    if ( fDebug >= 2 ) printf("AliAnalysisTaskPileup: position %i\n", position);
    trigScalerRecords1 = (AliTriggerScalersRecordESD*)fTriggerRunScalers->GetScalersRecordsESD()->At(position);

    // Sometimes scalers are filled very close to each others
    // in this case skip to the next entry
    for ( Int_t ipos=position+1; ipos<nPoints; ipos++ ) {
      trigScalerRecords2 = (AliTriggerScalersRecordESD*)fTriggerRunScalers->GetScalersRecordsESD()->At(ipos);
      Double_t deltaTime = (Double_t)( trigScalerRecords2->GetTimeStamp()->GetSeconds() - trigScalerRecords1->GetTimeStamp()->GetSeconds() );
      if ( fDebug >= 2 ) printf("AliAnalysisTaskPileup: Pos %i  TimeStamp %u - %u = %.0f\n", ipos, trigScalerRecords2->GetTimeStamp()->GetSeconds(), trigScalerRecords1->GetTimeStamp()->GetSeconds(), deltaTime);

      if ( deltaTime > 1 )
	break;
    } // loop on position
  } // nPoins > 0

#endif

  ULong64_t trigMask = 0;
  Double_t correctFactor = 1.;

  Int_t nVtxContrib = ( esdEvent ) ? esdEvent->GetPrimaryVertex()->GetNContributors() : aodEvent->GetPrimaryVertex()->GetNContributors();
  Bool_t isPileupSPD = ( esdEvent ) ? esdEvent->IsPileupFromSPD(3,0.8) : aodEvent->IsPileupFromSPD(3,0.8);

  Double_t centralityClass = InputEvent()->GetCentrality()->GetCentralityPercentile("V0M");
  // If percentile is exactly 100. the bin chosen is in overflow
  // fix it by setting 99.999 instead of 100.
  if ( centralityClass >= 100. ) centralityClass = 99.999;
  Int_t centralityBin = fCentralityClasses->FindBin(centralityClass);

  TString vtxSelKey[3] = {"any","hasVtxContrib","nonPileupSPD"};
  Bool_t fillSel[3] = {kTRUE, ( nVtxContrib > 0 ), ( ( nVtxContrib > 0 ) && ( ! isPileupSPD ) )};

  //const AliTriggerScalersRecordESD* trigScalerRecords = esdEvent->GetHeader()->GetTriggerScalersRecord(); // REMEMBER TO CUT

  TString trigName = "";
  TString eventType = "";

  // Get entries counts the number of entries != 0
  // However the position 0 can be empty in MC (it is reserved to CBEAMB)
  // In this case with GetEntries the last entry is lost
  // Use GetEntriesFast instead
  Int_t nTriggerClasses = fTriggerClasses->GetEntriesFast();
  Int_t classIndex = -1;
#if defined(READOCDB)
  Double_t deltaScalersBeam = 0., deltaScalers = 0.;
#endif
  //Bool_t isFiredOnce = kFALSE;
  for (Int_t itrig=0; itrig<nTriggerClasses+1; itrig++) {

    Double_t correctFactorL0 = 1.;

    Bool_t isClassFired = kFALSE;

    if ( itrig < nTriggerClasses ) {

      if ( fIsInitCDB ) {
	// Check if current mask contains trigger
	classIndex = (*fTriggerClassIndex)[itrig];
	if ( classIndex < 0 ) continue; // Protection for MC (where BEAMB not present
	trigMask = ( 1ull << classIndex );
	isClassFired = ( trigMask & InputEvent()->GetTriggerMask() );
      }
      else
	isClassFired = kTRUE;

      trigName = ((TObjString*)fTriggerClasses->At(itrig))->GetString();

#if defined(READOCDB)
      if ( trigScalerRecords2 && ( isClassFired || itrig == 0 ) ) {
	// Get scalers
	const AliTriggerScalersESD* scaler1 = trigScalerRecords1->GetTriggerScalersForClass(classIndex+1);
	const AliTriggerScalersESD* scaler2 = trigScalerRecords2->GetTriggerScalersForClass(classIndex+1);
	deltaScalers = scaler2->GetLOCB() - scaler1->GetLOCB();
 
	if ( itrig == 0 )
	  deltaScalersBeam = deltaScalers;
	else if ( isClassFired ) {
	  correctFactorL0 = GetL0Correction(deltaScalers, deltaScalersBeam);
	  if ( fDebug >= 2 ) printf("AliAnalysisTaskPileup: Scalers: %s %.0f  %s %.0f -> CF %f\n", fTriggerClasses->At(itrig)->GetName(), deltaScalers, fTriggerClasses->At(0)->GetName(), deltaScalersBeam, correctFactorL0);
	}
      }
#endif
    } // if ( itrig < nTriggerClasses )
    else {
      classIndex = -1;
      trigName = "any";
      isClassFired = kTRUE; // isFiredOnce;
    }

    if ( ! isClassFired ) continue;
    //isFiredOnce = kTRUE;

    //const AliTriggerScalersESD* trigScaler = trigScalerRecords->GetTriggerScalersForClass(classIndex+1); // REMEMBER TO CUT
    //if ( classIndex > 1 ) printf("Index: trigger %i  scaler %i\n", classIndex+1, trigScaler->GetClassIndex()); // REMEMBER TO CUT

    if ( fDebug >= 2 ) printf("AliAnalysisTaskPileup: Fired trig %s\n", trigName.Data());

    for ( Int_t ievType=0; ievType<2; ievType++ ){
      switch ( ievType ) {
      case kHeventsCorrectL0:
	correctFactor = correctFactorL0;
	eventType = "correctedL0";
	break;
      default:
	correctFactor = 1.;
	eventType = "any";
	break;
      }

      for ( Int_t isel=0; isel<3; isel++ ) {
	if ( ! fillSel[isel] ) continue;
	fEventCounters->Count(Form("event:%s/trigger:%s/vtxSelection:%s/selected:%s/centrality:%s",eventType.Data(),trigName.Data(), vtxSelKey[isel].Data(), selected.Data(),fCentralityClasses->GetBinLabel(centralityBin)),correctFactor);
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

  // Fill the container only at the very last step
  // i.e. in local, when done interactively
  if ( gROOT->IsBatch() )
    return;

  fEventCounters = dynamic_cast<AliCounterCollection*>(GetOutputData(1));
  if ( ! fEventCounters ) return;

  if ( ! fHistoEventsList ) fHistoEventsList = new TObjArray(0);
  fHistoEventsList->SetOwner();

  TH2D* histo = 0x0;
  const Int_t kNevtTypes = 2;
  TString evtSel[kNevtTypes] = {"event:any", "event:correctedL0"};
  TString evtName[kNevtTypes] = {"", "CorrectL0"};
  TString evtTitle[kNevtTypes] = {"Events", "L0 corrected events"};
  const Int_t kNphysSel = 2;
  TString physSel[kNphysSel] = {"selected:any","selected:yes"};
  TString physName[kNphysSel] = {"", "PhysSel"};
  TString physTitle[kNphysSel] = {"", "w/ physics selection"};

  TString typeName[2] = {"", "Centrality"};

  Int_t currHisto = -1;
  TString currName = "";
  for ( Int_t itype=0; itype<2; itype++ ) {
    for ( Int_t isel=0; isel<kNphysSel; isel++ ) {
      for ( Int_t iev=0; iev<kNevtTypes; iev++ ) {
	currName = Form("%s/%s", evtSel[iev].Data(), physSel[isel].Data());
	if ( itype == 0 )
	  histo = fEventCounters->Get("trigger","vtxSelection",currName.Data());
	else {
	  currName.Append("/vtxSelection:any");
	  histo = fEventCounters->Get("trigger","centrality",currName.Data());
	}
	if ( ! histo ) continue;
	currHisto++;
	currName = Form("hEvents%s%s%s", typeName[itype].Data(), evtName[iev].Data(), physName[isel].Data());
	histo->SetName(currName.Data());
	currName = Form("%s %s", evtTitle[iev].Data(), physTitle[isel].Data());
	histo->SetTitle(currName.Data());
	fHistoEventsList->AddAtAndExpand(histo, currHisto);
      } // loop on event type
      TH2D* num = (TH2D*)fHistoEventsList->At(currHisto);
      TH2D* den = (TH2D*)fHistoEventsList->At(currHisto-1);
      if ( ! num || ! den ) continue;
      currName = Form("hPileupL0%s%s_%s", typeName[itype].Data(), physName[isel].Data(), GetName());
      TH2D* histoPileupL0 = (TH2D*)num->Clone(currName.Data());
      histoPileupL0->Divide(den);
      currName.ReplaceAll("hPileupL0","canPileupL0");
      TCanvas *can = new TCanvas(currName.Data(),"Pileup",10,10,310,310);
      can->SetFillColor(10); can->SetHighLightColor(10);
      can->SetLeftMargin(0.15); can->SetBottomMargin(0.15);
      histoPileupL0->DrawCopy("text");
      delete histoPileupL0;
    } // loop on physics selection
  } // loop on histo type

  PostData(2, fHistoEventsList);
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

//________________________________________________________________________
void AliAnalysisTaskPileup::SetDefaultStorage(TString dbString)
{
  /// Set default storage
  SetSpecificStorage("default", dbString);
}


//________________________________________________________________________
void AliAnalysisTaskPileup::SetSpecificStorage(TString calibType, TString dbString)
{
  /// Set specific storage
#if defined(READOCDB)
  if ( ! fStorageList ) {
    fStorageList = new TObjArray(5);
    fStorageList->SetOwner();
  }
  fStorageList->AddLast(new TObjString(calibType));
  fStorageList->AddLast(new TObjString(dbString));
#else
  calibType = "";
  dbString  = "";
  AliWarning(Form("Class was not compiled to run on OCDB. Command will not have effect"));
#endif
}

