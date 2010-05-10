/* $Id: AliPhysicsSelection.cxx 35782 2009-10-22 11:54:31Z jgrosseo $ */

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
//                      Implementation of   Class AliPhysicsSelection
// This class selects collision candidates from data runs, applying selection cuts on triggers 
// and background rejection based on the content of the ESD
//
// Usage:
//
// Create the object:
//   fPhysicsSelection = new AliPhysicsSelection;
//
// For MC data, call
//   fPhysicsSelection->SetAnalyzeMC()
//
// To check if an event is a collision candidate, use:
//   fPhysicsSelection->IsCollisionCandidate(fESD)
//
// After processing save the resulting histograms to a file with (a folder physics_selection 
//   will be created that contains the histograms):
//   fPhysicsSelection->SaveHistograms("physics_selection")
//
// To print statistics after processing use:
//   fPhysicsSelection->Print();
//
// The BX ids corresponding to real bunches crossings p2 are
// automatically selected. You cannot process runs with different
// filling schemes if this option is set. If you want to disable this,
// use: 
//   fPhysicsSelection->SetUseBXNumbers(0);
//
//
// If you are analizing muons and you want to keep the muon triggers
// besides the CINT1B you can set:
//   fPhysicsSelection->SetUseMuonTriggers();
//
//
// To compute the Background automatically using the control triggers
// use: 
//   fPhysicsSelection->SetComputeBG();
// this will show the value of the Beam Gas, accidentals and good
// events as additional rows in the statistic tables, but it will NOT
// subtract the background automatically.
// This option will only work for runs taken with the CINT1
// suite. This options enables automatically also the usage of BX
// numbers. You can only process one run at a time if you require this
// options, because it uses the bunch intensity estimated run by run.
// 
// The BG will usually be more important in the so-called "bin 0": the
// class can also compute the statistics table for events in this
// bin. Since the definition of bin 0 may in general change from
// analysis to analysis, the user needs to provide a callback
// implementing the definition of bin zero. The callback should be
// implemented as a method in the analysis task and should override
// the IsEventInBinZero method of AliAnalysisTaskSE, and should thus
// have the the following prototype:
//   Bool_t IsEventInBinZero(); 
// It should return true if the event is in the bin 0 and it is set by
// passing to the physics selection the NAME of the task where the
// callback is implemented: 
//   fPhysicsSelection->SetBin0Callback("MyTask").
//
//
// Usually the class selects the trigger scheme by itself depending on the run number.
// Nevertheless, you can do that manually by calling AddCollisionTriggerClass() and AddBGTriggerClass()
// Example:
// To define the class CINT1B-ABCE-NOPF-ALL as collision trigger (those will be accepted as  
// collision candidates when they pass the selection):
//   AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL #769 #3119");
// To select on bunch crossing IDs in addition, use:
//   AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL #769 #3119");
// To define the class CINT1A-ABCE-NOPF-ALL as a background trigger (those will only be counted
// for the control histograms):
//   AddBGTriggerClass("+CINT1A-ABCE-NOPF-ALL");
// You can also specify more than one trigger class in a string or you can require that some are *not*
// present. The following line would require CSMBA-ABCE-NOPF-ALL, but CSMBB-ABCE-NOPF-ALL is not allowed
// to be present:
//   AddBGTriggerClass("+CSMBA-ABCE-NOPF-ALL -CSMBB-ABCE-NOPF-ALL");
//
//   Origin: Jan Fiete Grosse-Oetringhaus, CERN 
//           Michele Floris, CERN
//-------------------------------------------------------------------------

#include <Riostream.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TIterator.h>
#include <TDirectory.h>
#include <TObjArray.h>

#include <AliPhysicsSelection.h>

#include <AliTriggerAnalysis.h>
#include <AliLog.h>

#include <AliESDEvent.h>
#include <AliAnalysisTaskSE.h>
#include "AliAnalysisManager.h"

ClassImp(AliPhysicsSelection)

AliPhysicsSelection::AliPhysicsSelection() :
  AliAnalysisCuts("AliPhysicsSelection", "AliPhysicsSelection"),
  fCurrentRun(-1),
  fMC(kFALSE),
  fCollTrigClasses(),
  fBGTrigClasses(),
  fTriggerAnalysis(),
  fBackgroundIdentification(0),
  fHistBunchCrossing(0),
  fHistTriggerPattern(0),
  fSkipTriggerClassSelection(0),
  fUsingCustomClasses(0),
  fSkipV0(0),
  fBIFactorA(1),
  fBIFactorC(1),
  fRatioBEEE(2),
  fComputeBG(0),
  fUseBXNumbers(1),
  fUseMuonTriggers(0),
  fFillingScheme(""),
  fBin0CallBack(""),
  fBin0CallBackPointer(0)
{
  // constructor
  
  fCollTrigClasses.SetOwner(1);
  fBGTrigClasses.SetOwner(1);
  fTriggerAnalysis.SetOwner(1);
  fHistStatistics[0] = 0;
  fHistStatistics[1] = 0;
  
  AliLog::SetClassDebugLevel("AliPhysicsSelection", AliLog::kWarning);
}
    
AliPhysicsSelection::~AliPhysicsSelection()
{
  // destructor

  fCollTrigClasses.Delete();
  fBGTrigClasses.Delete();
  fTriggerAnalysis.Delete();

  if (fHistStatistics[0])
  {
    delete fHistStatistics[0];
    fHistStatistics[0] = 0;
  }
  if (fHistStatistics[1])
  {
    delete fHistStatistics[1];
    fHistStatistics[1] = 0;
  }
  
  if (fHistBunchCrossing)
  {
    delete fHistBunchCrossing;
    fHistBunchCrossing = 0;
  }
  if (fHistTriggerPattern)
  {
    delete fHistTriggerPattern;
    fHistTriggerPattern = 0;
  }

}

Bool_t AliPhysicsSelection::CheckTriggerClass(const AliESDEvent* aEsd, const char* trigger) const
{
  // checks if the given trigger class(es) are found for the current event
  // format of trigger: +TRIGGER1 -TRIGGER2
  //   requires TRIGGER1 and rejects TRIGGER2
  
  Bool_t foundBCRequirement = kFALSE;
  Bool_t foundCorrectBC = kFALSE;
  
  TString str(trigger);
  TObjArray* tokens = str.Tokenize(" ");
  
  for (Int_t i=0; i < tokens->GetEntries(); i++)
  {
    TString str2(((TObjString*) tokens->At(i))->String());
    
    if (str2[0] == '+' || str2[0] == '-')
    {
      Bool_t flag = (str2[0] == '+');
      
      str2.Remove(0, 1);
      
      if (flag && !aEsd->IsTriggerClassFired(str2))
      {
        AliDebug(AliLog::kDebug, Form("Rejecting event because trigger class %s is not present", str2.Data()));
        delete tokens;
        return kFALSE;
      }
      if (!flag && aEsd->IsTriggerClassFired(str2))
      {
        AliDebug(AliLog::kDebug, Form("Rejecting event because trigger class %s is present", str2.Data()));
        delete tokens;
        return kFALSE;
      }
    }
    else if (str2[0] == '#')
    {
      foundBCRequirement = kTRUE;
    
      str2.Remove(0, 1);
      
      Int_t bcNumber = str2.Atoi();
      AliDebug(AliLog::kDebug, Form("Checking for bunch crossing number %d", bcNumber));
      
      if (aEsd->GetBunchCrossNumber() == bcNumber)
      {
        foundCorrectBC = kTRUE;
        AliDebug(AliLog::kDebug, Form("Found correct bunch crossing %d", bcNumber));
      }
    }
    else
      AliFatal(Form("Invalid trigger syntax: %s", trigger));
  }
  
  delete tokens;
  
  if (foundBCRequirement && !foundCorrectBC)
    return kFALSE;
  
  return kTRUE;
}
    
Bool_t AliPhysicsSelection::IsCollisionCandidate(const AliESDEvent* aEsd)
{
  // checks if the given event is a collision candidate
  
  if (fCurrentRun != aEsd->GetRunNumber())
    if (!Initialize(aEsd->GetRunNumber()))
      AliFatal(Form("Could not initialize for run %d", aEsd->GetRunNumber()));
    
  const AliESDHeader* esdHeader = aEsd->GetHeader();
  if (!esdHeader)
  {
    AliError("ESD Header could not be retrieved");
    return kFALSE;
  }
  
  // check event type; should be PHYSICS = 7 for data and 0 for MC
  if (!fMC)
  {
    if (esdHeader->GetEventType() != 7)
      return kFALSE;
  }
  else
  {
    if (esdHeader->GetEventType() != 0)
      AliFatal(Form("Invalid event type for MC: %d", esdHeader->GetEventType()));
  }
  
  Bool_t accept = kFALSE;
    
  Int_t count = fCollTrigClasses.GetEntries() + fBGTrigClasses.GetEntries();
  for (Int_t i=0; i < count; i++)
  {
    const char* triggerClass = 0;
    if (i < fCollTrigClasses.GetEntries())
      triggerClass = ((TObjString*) fCollTrigClasses.At(i))->String();
    else
      triggerClass = ((TObjString*) fBGTrigClasses.At(i - fCollTrigClasses.GetEntries()))->String();
  
    AliDebug(AliLog::kDebug, Form("Processing trigger class %s", triggerClass));
  
    AliTriggerAnalysis* triggerAnalysis = static_cast<AliTriggerAnalysis*> (fTriggerAnalysis.At(i));
  
    triggerAnalysis->FillTriggerClasses(aEsd);
    
    if (CheckTriggerClass(aEsd, triggerClass))
    {
      triggerAnalysis->FillHistograms(aEsd);
  
      Bool_t isBin0 = kFALSE;
      if (fBin0CallBack != "") {
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
	  AliError("Cannot get the analysis manager");
	}  
	else {
	  isBin0 = ((AliAnalysisTaskSE*)mgr->GetTask(fBin0CallBack.Data()))->IsEventInBinZero();
	}
      } else if (fBin0CallBackPointer) {
	  isBin0 = (*fBin0CallBackPointer)(aEsd);
	
      }
      

      
      // hardware trigger (should only remove events for MC)
      // replay CINT1B hardware trigger
      // TODO this has to depend on the actual hardware trigger (and that depends on the run...)
      Int_t fastORHW = triggerAnalysis->SPDFiredChips(aEsd, 1); // SPD number of chips from trigger bits (!)
      Bool_t v0A       = fSkipV0 ? 0 :triggerAnalysis->IsOfflineTriggerFired(aEsd, AliTriggerAnalysis::kV0A);
      Bool_t v0C       = fSkipV0 ? 0 :triggerAnalysis->IsOfflineTriggerFired(aEsd, AliTriggerAnalysis::kV0C);
      Bool_t v0AHW     = fSkipV0 ? 0 :(triggerAnalysis->V0Trigger(aEsd, AliTriggerAnalysis::kASide, kTRUE) == AliTriggerAnalysis::kV0BB);// should replay hw trigger
      Bool_t v0CHW     = fSkipV0 ? 0 :(triggerAnalysis->V0Trigger(aEsd, AliTriggerAnalysis::kCSide, kTRUE) == AliTriggerAnalysis::kV0BB);// should replay hw trigger
      // offline trigger
      Int_t fastOROffline = triggerAnalysis->SPDFiredChips(aEsd, 0); // SPD number of chips from clusters (!)
      Bool_t v0ABG = fSkipV0 ? 0 :triggerAnalysis->IsOfflineTriggerFired(aEsd, AliTriggerAnalysis::kV0ABG);
      Bool_t v0CBG = fSkipV0 ? 0 :triggerAnalysis->IsOfflineTriggerFired(aEsd, AliTriggerAnalysis::kV0CBG);
      Bool_t v0BG = v0ABG || v0CBG;

      // fmd
      Bool_t fmdA =  triggerAnalysis->IsOfflineTriggerFired(aEsd, AliTriggerAnalysis::kFMDA);
      Bool_t fmdC =  triggerAnalysis->IsOfflineTriggerFired(aEsd, AliTriggerAnalysis::kFMDC);
      Bool_t fmd  = fmdA || fmdC;
    
      // SSD
      Int_t ssdClusters = triggerAnalysis->SSDClusters(aEsd);

      // Some "macros"
      Bool_t mb1 = (fastOROffline > 0 || v0A || v0C) && (!v0BG);
      Bool_t mb1prime = (fastOROffline > 1 || (fastOROffline > 0 && (v0A || v0C)) || (v0A && v0C) ) && (!v0BG);

      // Background rejection
      Bool_t bgID = ! fBackgroundIdentification->IsSelected(const_cast<AliESDEvent*> (aEsd));

      Int_t ntrig = fastOROffline; // any 2 hits
      if(v0A)              ntrig += 1;
      if(v0C)              ntrig += 1; //v0C alone is enough
      if(fmd)              ntrig += 1;
      if(ssdClusters>1)    ntrig += 1;


      Bool_t hwTrig = fastORHW > 0 || v0AHW || v0CHW;

      // Fill trigger pattern histo
      Int_t tpatt = 0;
      if (fastORHW>0) tpatt+=1;
      if (v0AHW)      tpatt+=2;
      if (v0CHW)      tpatt+=4;
      fHistTriggerPattern->Fill( tpatt );

      // fill statistics and return decision
      const Int_t nHistStat = 2;
      for(Int_t iHistStat = 0; iHistStat < nHistStat; iHistStat++){
	if (iHistStat == kStatIdxBin0 && !isBin0) continue; // skip the filling of bin0 stats if the event is not in the bin0
      
	fHistStatistics[iHistStat]->Fill(kStatTriggerClass, i);



	// We fill the rest only if hw trigger is ok
	if (!hwTrig)
	  {
	    AliDebug(AliLog::kDebug, "Rejecting event because hardware trigger is not fired");
	    continue;
	  } else {       
	  fHistStatistics[iHistStat]->Fill(kStatHWTrig, i);
	}
      

	// v0 BG stats
	if (v0ABG)
	  fHistStatistics[iHistStat]->Fill(kStatV0ABG, i);
	if (v0CBG)
	  fHistStatistics[iHistStat]->Fill(kStatV0CBG, i);

	// We fill the rest only if mb1 && ! v0BG
	if (mb1)
	  fHistStatistics[iHistStat]->Fill(kStatMB1, i);
	else continue;

	if (mb1prime)
	  fHistStatistics[iHistStat]->Fill(kStatMB1Prime, i);

	if (fmd)
	  fHistStatistics[iHistStat]->Fill(kStatFMD, i);

	if(ssdClusters>1)
	  fHistStatistics[iHistStat]->Fill(kStatSSD1, i);

	if(ntrig >= 2 && !v0BG) 
	  fHistStatistics[iHistStat]->Fill(kStatAny2Hits, i);

	if (fastOROffline > 0)
	  fHistStatistics[iHistStat]->Fill(kStatFO1, i);
	if (fastOROffline > 1)
	  fHistStatistics[iHistStat]->Fill(kStatFO2, i);
        
	if (v0A)
	  fHistStatistics[iHistStat]->Fill(kStatV0A, i);
	if (v0C)
	  fHistStatistics[iHistStat]->Fill(kStatV0C, i);
        
	//       if (fastOROffline > 1 && !v0BG)
	//         fHistStatistics[iHistStat]->Fill(kStatFO2NoBG, i);
            
	if (fastOROffline > 0 && (v0A || v0C) && !v0BG)
	  fHistStatistics[iHistStat]->Fill(kStatFO1AndV0, i);
  
	if (v0A && v0C && !v0BG && !bgID)
	  fHistStatistics[iHistStat]->Fill(kStatV0, i);


	if ( mb1 )
	
	  {
	    if (!v0BG || fSkipV0)
	      {
		if (!v0BG) fHistStatistics[iHistStat]->Fill(kStatOffline, i);
      
		if (fBackgroundIdentification && bgID)
		  {
		    AliDebug(AliLog::kDebug, "Rejecting event because of background identification");
		    fHistStatistics[iHistStat]->Fill(kStatBG, i);
		  }
		else
		  {
		    AliDebug(AliLog::kDebug, "Accepted event for histograms");
            
		    fHistStatistics[iHistStat]->Fill(kStatAccepted, i);
		    if(iHistStat == kStatIdxAll) fHistBunchCrossing->Fill(aEsd->GetBunchCrossNumber(), i); // Fill only for all (avoid double counting)
		    if((i < fCollTrigClasses.GetEntries() || fSkipTriggerClassSelection) && (iHistStat==kStatIdxAll))
		      accept = kTRUE; // only set for "all" (should not really matter)
		  }
	      }
	    else
	      AliDebug(AliLog::kDebug, "Rejecting event because of V0 BG flag");
	  }
	else
	  AliDebug(AliLog::kDebug, "Rejecting event because trigger condition is not fulfilled");
      }
    }
  }
 
  if (accept)
    AliDebug(AliLog::kDebug, "Accepted event as collision candidate");
  
  return accept;
}

Int_t AliPhysicsSelection::GetTriggerScheme(UInt_t runNumber) const
{
  // returns the current trigger scheme (classes that are accepted/rejected)
  
  if (fMC)
    return 0;
    
  // TODO dependent on run number
  
  switch (runNumber)
  {
    // CSMBB triggers
    case 104044:
    case 105054:
    case 105057:
      return 2;
  }
  
  // default: CINT1 suite
  return 1;
}  

const char * AliPhysicsSelection::GetFillingScheme(UInt_t runNumber)  {

  if(fMC) return "MC";

  if      (runNumber >= 104065 && runNumber <= 104160) {
    return "4x4a";
  } 
  else if (runNumber >= 104315 && runNumber <= 104321) {
    return "4x4a*";
  }
  else if (runNumber >= 104792 && runNumber <= 104803) {
    return "4x4b";
  }
  else if (runNumber >= 104824 && runNumber <= 104892) {
    return "4x4c";
  }
  else if (runNumber == 105143 || runNumber == 105160) {
    return "16x16a";
  }
  else if (runNumber >= 105256 && runNumber <= 105268) {
    return "4x4c";
  } else if (runNumber >= 114786 && runNumber <= 116684) {
    return "Single_2b_1_1_1";
  }
  else {
    AliError(Form("Unknown filling scheme (run %d)", runNumber));
  }

  return "Unknown";
}

Int_t AliPhysicsSelection::GetRatioBBBE(Int_t runNumber) {


  if(fMC) return 1;

  if (runNumber == 105143 || runNumber == 105160) {
    return 8;
  }else if (runNumber == 114786 || runNumber == 114798 ) {
    return 1;
  } else if (runNumber >= 114783 && runNumber <= 116684){
    return 1;
  }
  else if (fComputeBG &&
	   !(runNumber >= 105256 && runNumber <= 105268) &&
	   !(runNumber >= 104065 && runNumber <= 104160) &&
	   !(runNumber >= 104315 && runNumber <= 104321) &&
	   !(runNumber >= 104792 && runNumber <= 104803) &&
	   !(runNumber >= 104824 && runNumber <= 104892)
	   ){     

    AliError(Form("Unknown run %d, assuming ratio BE/EE = 2",runNumber));

  }

  return 2;
}


const char * AliPhysicsSelection::GetBXIDs(UInt_t runNumber, const char * trigger)  {

  if (!fUseBXNumbers || fMC) return "";

  if      (runNumber >= 104065 && runNumber <= 104160) {
    if     (!strcmp("CINT1B-ABCE-NOPF-ALL",trigger)) return " #2128 #3019";
    else if(!strcmp("CINT1A-ABCE-NOPF-ALL",trigger)) return " #346 #3465";
    else if(!strcmp("CINT1C-ABCE-NOPF-ALL",trigger)) return " #1234 #1680";
    else if(!strcmp("CINT1-E-NOPF-ALL",trigger))     return " #790";
    else AliError(Form("Unknown trigger: %s", trigger));
  } 
  else if (runNumber >= 104315 && runNumber <= 104321) {
    if     (!strcmp("CINT1B-ABCE-NOPF-ALL",trigger)) return " #2000 #2891";
    else if(!strcmp("CINT1A-ABCE-NOPF-ALL",trigger)) return " #218 #3337";
    else if(!strcmp("CINT1C-ABCE-NOPF-ALL",trigger)) return " #1106 #1552";
    else if(!strcmp("CINT1-E-NOPF-ALL",trigger))     return " #790";
    else AliError(Form("Unknown trigger: %s", trigger));
  }
  else if (runNumber >= 104792 && runNumber <= 104803) {
    if     (!strcmp("CINT1B-ABCE-NOPF-ALL",trigger)) return " #2228 #3119";
    else if(!strcmp("CINT1A-ABCE-NOPF-ALL",trigger)) return " #2554 #446";
    else if(!strcmp("CINT1C-ABCE-NOPF-ALL",trigger)) return " #1334 #769";
    else if(!strcmp("CINT1-E-NOPF-ALL",trigger))     return " #790";
    else AliError(Form("Unknown trigger: %s", trigger));
  }
  else if (runNumber >= 104824 && runNumber <= 104892) {
    if     (!strcmp("CINT1B-ABCE-NOPF-ALL",trigger)) return " #3119 #769";
    else if(!strcmp("CINT1A-ABCE-NOPF-ALL",trigger)) return " #2554 #446";
    else if(!strcmp("CINT1C-ABCE-NOPF-ALL",trigger)) return " #1334 #2228";
    else if(!strcmp("CINT1-E-NOPF-ALL",trigger))     return " #790";
    else AliError(Form("Unknown trigger: %s", trigger));
  }
  else if (runNumber == 105143 || runNumber == 105160) {
    fRatioBEEE = 8;
    if     (!strcmp("CINT1B-ABCE-NOPF-ALL",trigger)) return " #1337 #1418 #2228 #2309 #3119 #3200 #446 #527";
    else if(!strcmp("CINT1A-ABCE-NOPF-ALL",trigger)) return " #1580  #1742  #1904  #2066  #2630  #2792  #2954  #3362";
    else if(!strcmp("CINT1C-ABCE-NOPF-ALL",trigger)) return "  #845  #1007  #1169   #1577 #3359 #3521 #119  #281 ";
    else if(!strcmp("CINT1-E-NOPF-ALL",trigger))     return " #790";
    else AliError(Form("Unknown trigger: %s", trigger));
  }
  else if (runNumber >= 105256 && runNumber <= 105268) {
    if     (!strcmp("CINT1B-ABCE-NOPF-ALL",trigger)) return " #3019 #669";
    else if(!strcmp("CINT1A-ABCE-NOPF-ALL",trigger)) return " #2454 #346";
    else if(!strcmp("CINT1C-ABCE-NOPF-ALL",trigger)) return " #1234 #2128";
    else if(!strcmp("CINT1-E-NOPF-ALL",trigger))     return " #790";
    else AliError(Form("Unknown trigger: %s", trigger));
  } else if (runNumber >= 114786 && runNumber <= 116684) { // 7 TeV 2010, assume always the same filling scheme
    if     (!strcmp("CINT1B-ABCE-NOPF-ALL",trigger)) return " #346";
    else if(!strcmp("CINT1A-ABCE-NOPF-ALL",trigger)) return " #2131";
    else if(!strcmp("CINT1C-ABCE-NOPF-ALL",trigger)) return " #3019";
    else if(!strcmp("CINT1-E-NOPF-ALL",trigger))     return " #1238";
    else AliError(Form("Unknown trigger: %s", trigger));
  }
  else {
    AliError(Form("Unknown run %d, using all BXs!",runNumber));
  }

  return "";
}
    
Bool_t AliPhysicsSelection::Initialize(Int_t runNumber)
{
  // initializes the object for the given run
  // TODO having the run number here and parameters hardcoded is clearly temporary, a way needs to be found to have a CDB-like configuration also for analysis
  
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  if(!fBin0CallBack) 
    AliError("Bin0 Callback not set: will not fill the statistics for the bin 0");

  if (fMC) {
    // ovverride BX and bg options in case of MC
    fComputeBG    = kFALSE;
    fUseBXNumbers = kFALSE;
  }

  Int_t triggerScheme = GetTriggerScheme(runNumber);
  if (!fUsingCustomClasses && fCurrentRun != -1 && triggerScheme != GetTriggerScheme(fCurrentRun))
    AliFatal("Processing several runs with different trigger schemes is not supported");
  
  if(fComputeBG && fCurrentRun != -1 && fCurrentRun != runNumber) 
    AliFatal("Cannot process several runs because BG computation is requested");

  if(fComputeBG && !fUseBXNumbers) 
    AliFatal("Cannot compute BG id BX numbers are not used");
  
  if(fUseBXNumbers && fFillingScheme != "" && fFillingScheme != GetFillingScheme(runNumber))
    AliFatal("Cannot process runs with different filling scheme if usage of BX numbers is requested");

  fRatioBEEE          = GetRatioBBBE(runNumber);
  fFillingScheme      = GetFillingScheme(runNumber);

  if(fComputeBG) SetBIFactors(runNumber);

  AliInfo(Form("Initializing for run %d", runNumber));
  
  // initialize first time?
  if (fCurrentRun == -1)
  {
    if (fUsingCustomClasses) {
      AliInfo("Using user-provided trigger classes");
    } else {
      switch (triggerScheme)
      {
      case 0:
        fCollTrigClasses.Add(new TObjString(""));
        break;
        
      case 1:
	{ // need a new scope to avoid cross-initialization errors

	  if (fUseMuonTriggers) {
	    // Muon trigger have the same BXIDs of the corresponding CINT triggers
	    fCollTrigClasses.Add(new TObjString(Form("%s%s ","+CMUS1B-ABCE-NOPF-MUON",  GetBXIDs(runNumber,"CINT1B-ABCE-NOPF-ALL"))));
	    fBGTrigClasses.Add  (new TObjString(Form("%s%s ","+CMUS1A-ABCE-NOPF-MUON",  GetBXIDs(runNumber,"CINT1A-ABCE-NOPF-ALL"))));
	    fBGTrigClasses.Add  (new TObjString(Form("%s%s ","+CMUS1C-ABCE-NOPF-MUON",  GetBXIDs(runNumber,"CINT1C-ABCE-NOPF-ALL"))));	    
	    fBGTrigClasses.Add  (new TObjString(Form("%s%s ","+CMUS1-E-NOPF-MUON"    ,  GetBXIDs(runNumber,"CINT1-E-NOPF-ALL"))));
	  }
	  TObjString * cint1b = new TObjString(Form("%s%s","+CINT1B-ABCE-NOPF-ALL",  GetBXIDs(runNumber,"CINT1B-ABCE-NOPF-ALL")));
	  TObjString * cint1a = new TObjString(Form("%s%s","+CINT1A-ABCE-NOPF-ALL",  GetBXIDs(runNumber,"CINT1A-ABCE-NOPF-ALL")));
	  TObjString * cint1c = new TObjString(Form("%s%s","+CINT1C-ABCE-NOPF-ALL",  GetBXIDs(runNumber,"CINT1C-ABCE-NOPF-ALL")));
	  TObjString * cint1e = new TObjString(Form("%s%s","+CINT1-E-NOPF-ALL",      GetBXIDs(runNumber,"CINT1-E-NOPF-ALL"))    );
	  //
	  fCollTrigClasses.Add(cint1b);
	  fBGTrigClasses.Add(cint1a);
	  fBGTrigClasses.Add(cint1c);
	  fBGTrigClasses.Add(cint1e);

	}
        break;
        
      case 2:
        fCollTrigClasses.Add(new TObjString("+CSMBB-ABCE-NOPF-ALL"));
        fBGTrigClasses.Add(new TObjString("+CSMBA-ABCE-NOPF-ALL -CSMBB-ABCE-NOPF-ALL"));
        fBGTrigClasses.Add(new TObjString("+CSMBC-ABCE-NOPF-ALL -CSMBB-ABCE-NOPF-ALL"));
        break;
        
      default:
        AliFatal(Form("Unsupported trigger scheme %d", triggerScheme));
      }
    }
    
    Int_t count = fCollTrigClasses.GetEntries() + fBGTrigClasses.GetEntries();
    
    for (Int_t i=0; i<count; i++)
    {
      AliTriggerAnalysis* triggerAnalysis = new AliTriggerAnalysis;
      triggerAnalysis->SetAnalyzeMC(fMC);
      triggerAnalysis->EnableHistograms();
      triggerAnalysis->SetSPDGFOThreshhold(1);
      fTriggerAnalysis.Add(triggerAnalysis);
    }
      
    // TODO: shall I really delete this?
    if (fHistStatistics[0])
      delete fHistStatistics[0];
    if (fHistStatistics[1])
      delete fHistStatistics[1];
  
    fHistStatistics[kStatIdxBin0] = BookHistStatistics("_Bin0");
    fHistStatistics[kStatIdxAll]  = BookHistStatistics("");
    
    if (fHistBunchCrossing)
      delete fHistBunchCrossing;
  
    fHistBunchCrossing = new TH2F("fHistBunchCrossing", "fHistBunchCrossing;bunch crossing number;", 4000, -0.5, 3999.5,  count, -0.5, -0.5 + count);

    if (fHistTriggerPattern)
      delete fHistTriggerPattern;
    
    const int ntrig=3;
    Int_t n = 1;
    const Int_t nbinTrig = TMath::Nint(TMath::Power(2,ntrig));

    fHistTriggerPattern = new TH1F("fHistTriggerPattern", "Trigger pattern: FO + 2*v0A + 4*v0C", 
				   nbinTrig, -0.5, nbinTrig-0.5);    
    fHistTriggerPattern->GetXaxis()->SetBinLabel(1,"NO TRIG");
    fHistTriggerPattern->GetXaxis()->SetBinLabel(2,"FO");
    fHistTriggerPattern->GetXaxis()->SetBinLabel(3,"v0A");
    fHistTriggerPattern->GetXaxis()->SetBinLabel(4,"FO & v0A");
    fHistTriggerPattern->GetXaxis()->SetBinLabel(5,"v0C");
    fHistTriggerPattern->GetXaxis()->SetBinLabel(6,"FO & v0C");
    fHistTriggerPattern->GetXaxis()->SetBinLabel(7,"v0A & v0C");
    fHistTriggerPattern->GetXaxis()->SetBinLabel(8,"FO & v0A & v0C");

  
    n = 1;
    for (Int_t i=0; i < fCollTrigClasses.GetEntries(); i++)
    {
      fHistBunchCrossing->GetYaxis()->SetBinLabel(n, ((TObjString*) fCollTrigClasses.At(i))->String());
      n++;
    }
    for (Int_t i=0; i < fBGTrigClasses.GetEntries(); i++)
    {
      fHistBunchCrossing->GetYaxis()->SetBinLabel(n, ((TObjString*) fBGTrigClasses.At(i))->String());
      n++;
    }

    

  }
    
  Int_t count = fCollTrigClasses.GetEntries() + fBGTrigClasses.GetEntries();
  for (Int_t i=0; i<count; i++)
  {
    AliTriggerAnalysis* triggerAnalysis = static_cast<AliTriggerAnalysis*> (fTriggerAnalysis.At(i));
  
    switch (runNumber)
    {
      case 104315:
      case 104316:
      case 104320:
      case 104321:
      case 104439:
        triggerAnalysis->SetV0TimeOffset(7.5);
        break;
      default:
        triggerAnalysis->SetV0TimeOffset(0);
    }
  }
    
  fCurrentRun = runNumber;
  
  TH1::AddDirectory(oldStatus);
  
  return kTRUE;
}

TH2F * AliPhysicsSelection::BookHistStatistics(const char * tag) {

    // add 6 rows to count for the estimate of good, accidentals and
    // BG and the ratio of BG and accidentals to total +ratio goot to
    // first col + 2 for error on good.

  Int_t count = fCollTrigClasses.GetEntries() + fBGTrigClasses.GetEntries();
#ifdef VERBOSE_STAT
  Int_t extrarows = fComputeBG ? 8 : 0;
#else
  Int_t extrarows = fComputeBG ? 3 : 0;
#endif
  TH2F * h = new TH2F(Form("fHistStatistics%s",tag), Form("fHistStatistics - %s ;;",tag), kStatAccepted, 0.5, kStatAccepted+0.5, count+extrarows, -0.5, -0.5 + count+extrarows);

  h->GetXaxis()->SetBinLabel(kStatTriggerClass,  "Trigger class");
  h->GetXaxis()->SetBinLabel(kStatHWTrig,	 "Hardware trigger");
  h->GetXaxis()->SetBinLabel(kStatFO1,	         "FO >= 1");
  h->GetXaxis()->SetBinLabel(kStatFO2,	         "FO >= 2");
  h->GetXaxis()->SetBinLabel(kStatV0A,	         "V0A");
  h->GetXaxis()->SetBinLabel(kStatV0C,	         "V0C");
  h->GetXaxis()->SetBinLabel(kStatFMD,	         "FMD");
  h->GetXaxis()->SetBinLabel(kStatSSD1,	         "SSD >= 2");
  h->GetXaxis()->SetBinLabel(kStatV0ABG,	 "V0A BG");
  h->GetXaxis()->SetBinLabel(kStatV0CBG,	 "V0C BG");
  h->GetXaxis()->SetBinLabel(kStatMB1,	         "(FO >= 1 | V0A | V0C) & !V0 BG");
  h->GetXaxis()->SetBinLabel(kStatMB1Prime,      "(FO >= 2 | (FO >= 1 & (V0A | V0C)) | (V0A &v0C) ) & !V0 BG");
  h->GetXaxis()->SetBinLabel(kStatFO1AndV0,	 "FO >= 1 & (V0A | V0C) & !V0 BG");
  h->GetXaxis()->SetBinLabel(kStatV0,	         "V0A & V0C & !V0 BG & !BG ID");
  h->GetXaxis()->SetBinLabel(kStatOffline,	 "Offline Trigger");
  h->GetXaxis()->SetBinLabel(kStatAny2Hits,	 "2 Hits & !V0 BG");
  h->GetXaxis()->SetBinLabel(kStatBG,	         "Background identification");
  h->GetXaxis()->SetBinLabel(kStatAccepted,      "Accepted");

  Int_t n = 1;
  for (Int_t i=0; i < fCollTrigClasses.GetEntries(); i++)
    {
      h->GetYaxis()->SetBinLabel(n, ((TObjString*) fCollTrigClasses.At(i))->String());
      n++;
    }
  for (Int_t i=0; i < fBGTrigClasses.GetEntries(); i++)
    {
      h->GetYaxis()->SetBinLabel(n, ((TObjString*) fBGTrigClasses.At(i))->String());
      n++;
    }

  if(fComputeBG) {
    h->GetYaxis()->SetBinLabel(n++, "BG (A+C)");
    h->GetYaxis()->SetBinLabel(n++, "ACC");
#ifdef VERBOSE_STAT
    h->GetYaxis()->SetBinLabel(n++, "BG (A+C) %  (rel. to CINT1B)");
    h->GetYaxis()->SetBinLabel(n++, "ACC % (rel. to CINT1B)");
    h->GetYaxis()->SetBinLabel(n++, "ERR GOOD %");
    h->GetYaxis()->SetBinLabel(n++, "GOOD % (rel. to 1st col)");
    h->GetYaxis()->SetBinLabel(n++, "ERR GOOD");
#endif
    h->GetYaxis()->SetBinLabel(n++, "GOOD");
  }

  return h;
}

void AliPhysicsSelection::Print(Option_t* /* option */) const
{
  // print the configuration
  
  Printf("Configuration initialized for run %d (MC: %d):", fCurrentRun, fMC);
  
  Printf("Collision trigger classes:");
  for (Int_t i=0; i < fCollTrigClasses.GetEntries(); i++)
    Printf("%s", ((TObjString*) fCollTrigClasses.At(i))->String().Data());
  
  Printf("Background trigger classes:");
  for (Int_t i=0; i < fBGTrigClasses.GetEntries(); i++)
    Printf("%s", ((TObjString*) fBGTrigClasses.At(i))->String().Data());

  AliTriggerAnalysis* triggerAnalysis = dynamic_cast<AliTriggerAnalysis*> (fTriggerAnalysis.At(0));
  
  if (triggerAnalysis)
  {
    if (triggerAnalysis->GetV0TimeOffset() > 0)
      Printf("V0 time offset active: %.2f ns", triggerAnalysis->GetV0TimeOffset());
    
    Printf("\nTotal available events:");
    
    triggerAnalysis->PrintTriggerClasses();
  }
  
  if (fHistStatistics[kStatIdxAll] && fCollTrigClasses.GetEntries() > 0)
  {
    Printf("\nSelection statistics for first collision trigger (%s):", ((TObjString*) fCollTrigClasses.First())->String().Data());
    
    Printf("Total events with correct trigger class: %d", (Int_t) fHistStatistics[kStatIdxAll]->GetBinContent(1, 1));
    Printf("Selected collision candidates: %d", (Int_t) fHistStatistics[kStatIdxAll]->GetBinContent(fHistStatistics[kStatIdxAll]->GetXaxis()->FindBin("Accepted"), 1));
  }
  
  if (fHistBunchCrossing)
  {
    Printf("\nBunch crossing statistics:");
    
    for (Int_t i=1; i<=fHistBunchCrossing->GetNbinsY(); i++)
    {
      TString str;
      str.Form("Trigger %s has accepted events in the bunch crossings: ", fHistBunchCrossing->GetYaxis()->GetBinLabel(i));
      
      for (Int_t j=1; j<=fHistBunchCrossing->GetNbinsX(); j++)
        if (fHistBunchCrossing->GetBinContent(j, i) > 0)
          str += Form("%d, ", (Int_t) fHistBunchCrossing->GetXaxis()->GetBinCenter(j));
       
      Printf("%s", str.Data());
    }
    
    for (Int_t j=1; j<=fHistBunchCrossing->GetNbinsX(); j++)
    {
      Int_t count = 0;
      for (Int_t i=1; i<=fHistBunchCrossing->GetNbinsY(); i++)
      {
        if (fHistBunchCrossing->GetBinContent(j, i) > 0)
          count++;
      }
      if (count > 1)
        Printf("WARNING: Bunch crossing %d has more than one trigger class active. Check BPTX functioning for this run!", (Int_t) fHistBunchCrossing->GetXaxis()->GetBinCenter(j));
    }
  }

  if (fUsingCustomClasses)        
    Printf("WARNING: Using custom trigger classes!");
  if (fSkipTriggerClassSelection) 
    Printf("WARNING: Skipping trigger class selection!");
  if (fSkipV0) 
    Printf("WARNING: Ignoring V0 information in selection");
  if(!fBin0CallBack) 
    Printf("WARNING: Callback not set: will not fill the statistics for the bin 0");

}

Long64_t AliPhysicsSelection::Merge(TCollection* list)
{
  // Merge a list of AliMultiplicityCorrection objects with this (needed for
  // PROOF).
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;
  
  // collections of all histograms
  const Int_t nHists = 9;
  TList collections[nHists];

  Int_t count = 0;
  while ((obj = iter->Next())) {

    AliPhysicsSelection* entry = dynamic_cast<AliPhysicsSelection*> (obj);
    if (entry == 0) 
      continue;
      
    collections[0].Add(&(entry->fTriggerAnalysis));
    if (entry->fHistStatistics[0])
      collections[1].Add(entry->fHistStatistics[0]);
    if (entry->fHistStatistics[1])
      collections[2].Add(entry->fHistStatistics[1]);
    if (entry->fHistBunchCrossing)
      collections[3].Add(entry->fHistBunchCrossing);
    if (entry->fHistTriggerPattern)
      collections[4].Add(entry->fHistTriggerPattern);
    if (entry->fBackgroundIdentification)
      collections[5].Add(entry->fBackgroundIdentification);

    count++;
  }

  fTriggerAnalysis.Merge(&collections[0]);
  if (fHistStatistics[0])
    fHistStatistics[0]->Merge(&collections[1]);
  if (fHistStatistics[1])
    fHistStatistics[1]->Merge(&collections[2]);
  if (fHistBunchCrossing)
    fHistBunchCrossing->Merge(&collections[3]);
  if (fHistTriggerPattern)
    fHistTriggerPattern->Merge(&collections[4]);
  if (fBackgroundIdentification)
    fBackgroundIdentification->Merge(&collections[5]);
  
  delete iter;

  return count+1;
}

void AliPhysicsSelection::SaveHistograms(const char* folder) const
{
  // write histograms to current directory
  
  if (!fHistStatistics[0] || !fHistStatistics[1])
    return;
    
  if (folder)
  {
    gDirectory->mkdir(folder);
    gDirectory->cd(folder);
  }
  

  // Fill the last rows of fHistStatistics before saving
  if (fComputeBG) {
    Int_t triggerScheme = GetTriggerScheme(UInt_t(fCurrentRun));
    if(triggerScheme != 1){
      AliWarning("BG estimate only supported for trigger scheme \"1\" (CINT1 suite)");
    } else if (fUseMuonTriggers) {
      AliWarning("BG estimate with muon triggers to be implemented");
    } else {
      Int_t nHistStat = 2;
      // TODO: get number of rows in a more flexible way
      // 1. loop over all cols

      for(Int_t iHistStat = 0; iHistStat < nHistStat; iHistStat++){
    
	Int_t ncol = fHistStatistics[iHistStat]->GetNbinsX();
	Float_t good1 = 0;
	for(Int_t icol = 1; icol <= ncol; icol++) {
	  Int_t cint1B = (Int_t) fHistStatistics[iHistStat]->GetBinContent(icol,1);	
	  Int_t cint1A = (Int_t) fHistStatistics[iHistStat]->GetBinContent(icol,2);	
	  Int_t cint1C = (Int_t) fHistStatistics[iHistStat]->GetBinContent(icol,3);	
	  Int_t cint1E = (Int_t) fHistStatistics[iHistStat]->GetBinContent(icol,4);      
      
	  if (cint1B>0) {
	    Int_t acc  = fRatioBEEE*cint1E; 
	    Double_t acc_err = TMath::Sqrt(fRatioBEEE*fRatioBEEE*cint1E);
	    //      Int_t bg   = cint1A + cint1C - 2*acc;
	    Float_t bg   = fBIFactorA*(cint1A-acc) + fBIFactorC*(cint1C-acc) ;
	    Float_t good = Float_t(cint1B) - bg - acc;
	    if (icol ==1) good1 = good;
	    //      Float_t errGood     = TMath::Sqrt(2*(cint1A+cint1C+cint1E));// Error on the number of goods assuming only bg fluctuates
	    //      DeltaG^2 = B + FA^2 A + FC^2 C + Ratio^2 (FA+FC-1)^2 E.
	    Float_t errGood     = TMath::Sqrt( cint1B + 
					       fBIFactorA*fBIFactorA*cint1A +
					       fBIFactorC*fBIFactorC*cint1C +
					       fRatioBEEE * fRatioBEEE * 
					       (fBIFactorA + fBIFactorC - 1)*(fBIFactorA + fBIFactorC - 1)*cint1E);
	    Float_t errBG = TMath::Sqrt(fBIFactorA*fBIFactorA*cint1A+
					fBIFactorC*fBIFactorC*cint1C+
					fRatioBEEE*fRatioBEEE*(fBIFactorA+fBIFactorC)*(fBIFactorA+fBIFactorC)*cint1E);
	
	
	    fHistStatistics[iHistStat]->SetBinContent(icol,kStatRowBG,bg);	
	    fHistStatistics[iHistStat]->SetBinError  (icol,kStatRowBG,errBG);	
	    fHistStatistics[iHistStat]->SetBinContent(icol,kStatRowAcc,acc);	
	    fHistStatistics[iHistStat]->SetBinError  (icol,kStatRowAcc,acc_err);	
	    fHistStatistics[iHistStat]->SetBinContent(icol,kStatRowGood,good);    
	    fHistStatistics[iHistStat]->SetBinError  (icol,kStatRowGood,errGood);    

#ifdef VERBOSE_STAT
	    //kStatRowBG=5,kStatRowAcc,kStatRowBGFrac,kStatRowAccFrac,kStatRowErrGoodFrac,kStatRowGoodFrac,kStatRowGood,kStatRowErrGood
	    Float_t accFrac   = Float_t(acc) / cint1B  *100;
	    Float_t errAccFrac= Float_t(acc_err) / cint1B  *100;
	    Float_t bgFrac    = Float_t(bg)  / cint1B  *100;
	    Float_t goodFrac  = Float_t(good)  / good1 *100;
	    Float_t errGoodFrac = errGood/good1 * 100;
	    Float_t errFracBG = bg > 0 ? TMath::Sqrt((errBG/bg)*(errBG/bg) + 1/cint1B)*bgFrac : 0;
	    fHistStatistics[iHistStat]->SetBinContent(icol,kStatRowBGFrac,bgFrac);	
	    fHistStatistics[iHistStat]->SetBinError  (icol,kStatRowBGFrac,errFracBG);	
	    fHistStatistics[iHistStat]->SetBinContent(icol,kStatRowAccFrac,accFrac);    
	    fHistStatistics[iHistStat]->SetBinError  (icol,kStatRowAccFrac,errAccFrac);    
	    fHistStatistics[iHistStat]->SetBinContent(icol,kStatRowGoodFrac,goodFrac);    
	    fHistStatistics[iHistStat]->SetBinContent(icol,kStatRowErrGoodFrac,errGoodFrac);    
	    fHistStatistics[iHistStat]->SetBinContent(icol,kStatRowErrGood,errGood);    
#endif
	  }
	}
      }
    }  
  }

  fHistStatistics[0]->Write();
  fHistStatistics[1]->Write();
  if(fHistBunchCrossing ) fHistBunchCrossing ->Write();
  if(fHistTriggerPattern) fHistTriggerPattern->Write();
  
  Int_t count = fCollTrigClasses.GetEntries() + fBGTrigClasses.GetEntries();
  for (Int_t i=0; i < count; i++)
  {
    TString triggerClass = "trigger_histograms_";
    if (i < fCollTrigClasses.GetEntries())
      triggerClass += ((TObjString*) fCollTrigClasses.At(i))->String();
    else
      triggerClass += ((TObjString*) fBGTrigClasses.At(i - fCollTrigClasses.GetEntries()))->String();
  
    gDirectory->mkdir(triggerClass);
    gDirectory->cd(triggerClass);
  
    static_cast<AliTriggerAnalysis*> (fTriggerAnalysis.At(i))->SaveHistograms();
    
    gDirectory->cd("..");
  }
 
  if (fBackgroundIdentification)
  {
    gDirectory->mkdir("background_identification");
    gDirectory->cd("background_identification");
      
    fBackgroundIdentification->GetOutput()->Write();
      
    gDirectory->cd("..");
  }
  
  if (folder)
    gDirectory->cd("..");
}

void AliPhysicsSelection::SetBIFactors(Int_t run) {

  switch(run) {
  case 104155:
    fBIFactorA = 0.961912722908;
    fBIFactorC = 1.04992336081;
    break;
  case 104157:
    fBIFactorA = 0.947312854998;
    fBIFactorC = 1.01599706417;
    break;
  case 104159:
    fBIFactorA = 0.93659320151;
    fBIFactorC = 0.98580804207;
    break;
  case 104160:
    fBIFactorA = 0.929664189926;
    fBIFactorC = 0.963467679851;
    break;
  case 104315:
    fBIFactorA = 1.08939104979;
    fBIFactorC = 0.931113921925;
    break;
  case 104316:
    fBIFactorA = 1.08351880974;
    fBIFactorC = 0.916068345845;
    break;
  case 104320:
    fBIFactorA = 1.07669281245;
    fBIFactorC = 0.876818744763;
    break;
  case 104321:
    fBIFactorA = 1.00971079602;
    fBIFactorC = 0.773781299076;
    break;
  case 104792:
    fBIFactorA = 0.787215863962;
    fBIFactorC = 0.778253173071;
    break;
  case 104793:
    fBIFactorA = 0.692211363661;
    fBIFactorC = 0.733152456667;
    break;
  case 104799:
    fBIFactorA = 1.04027825161;
    fBIFactorC = 1.00530825942;
    break;
  case 104800:
    fBIFactorA = 1.05309910671;
    fBIFactorC = 1.00376801855;
    break;
  case 104801:
    fBIFactorA = 1.0531231922;
    fBIFactorC = 0.992439666758;
    break;
  case 104802:
    fBIFactorA = 1.04191478134;
    fBIFactorC = 0.979368585208;
    break;
  case 104803:
    fBIFactorA = 1.03121314094;
    fBIFactorC = 0.973379962609;
    break;
  case 104824:
    fBIFactorA = 0.969945926722;
    fBIFactorC = 0.39549745806;
    break;
  case 104825:
    fBIFactorA = 0.968627213937;
    fBIFactorC = 0.310100412205;
    break;
  case 104841:
    fBIFactorA = 0.991601393212;
    fBIFactorC = 0.83762204722;
    break;
  case 104845:
    fBIFactorA = 0.98040863886;
    fBIFactorC = 0.694824205793;
    break;
  case 104867:
    fBIFactorA = 1.10646173412;
    fBIFactorC = 0.841407246916;
    break;
  case 104876:
    fBIFactorA = 1.12063452421;
    fBIFactorC = 0.78726542895;
    break;
  case 104890:
    fBIFactorA = 1.02346137453;
    fBIFactorC = 1.03355663595;
    break;
  case 104892:
    fBIFactorA = 1.05406025913;
    fBIFactorC = 1.00029166135;
    break;
  case 105143:
    fBIFactorA = 0.947343384349;
    fBIFactorC = 0.972637444408;
    break;
  case 105160:
    fBIFactorA = 0.908854622177;
    fBIFactorC = 0.958851103977;
    break; 
  case 105256:
    fBIFactorA = 0.810076150206;
    fBIFactorC = 0.884663561883;
    break;
  case 105257:
    fBIFactorA = 0.80974912303;
    fBIFactorC = 0.878859123479;
    break;
  case 105268:
    fBIFactorA = 0.809052110679;
    fBIFactorC = 0.87233890989;
    break;
  default:
    fBIFactorA = 1;
    fBIFactorC = 1;
  }


}
