/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved.      *
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
// --- ROOT system ---
#include <iostream>
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH1I.h> 

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliFMDQADataMakerSim.h"
#include "AliFMDDigit.h"
#include "AliFMDHit.h"
#include "AliQAChecker.h"
#include "AliFMDParameters.h"

//_____________________________________________________________________
// This is the class that collects the QA data for the FMD during simulation.
// The following data types are picked up:
// - hits
// - digits
// The following data types are not supported (yet):
// - raws
// - sdigits
// Author : Hans Hjersing Dalsgaard, Niels Bohr Institute, hans.dalsgaard@cern.ch
//_____________________________________________________________________

ClassImp(AliFMDQADataMakerSim)
#if 0
; // This line is for Emacs - do not delete!
#endif
//_____________________________________________________________________
AliFMDQADataMakerSim::AliFMDQADataMakerSim() 
  :  AliQADataMakerSim(AliQA::GetDetName(AliQA::kFMD),
		       "FMD Quality Assurance Data Maker"),
     fDigitsArray(0),
     fHitsArray(0)
{
  // ctor
  fDigitsArray = new TClonesArray("AliFMDDigit", 1000) ; 
  fHitsArray   = new TClonesArray("AliFMDHit", 10);
}

//_____________________________________________________________________
AliFMDQADataMakerSim::AliFMDQADataMakerSim(const AliFMDQADataMakerSim& qadm) 
  : AliQADataMakerSim(),
    fDigitsArray(qadm.fDigitsArray),
    fHitsArray(qadm.fHitsArray)
{
  //copy ctor 
  
  // Parameters: 
  //    qadm    Object to copy from
  
}
//_____________________________________________________________________
AliFMDQADataMakerSim& AliFMDQADataMakerSim::operator = (const AliFMDQADataMakerSim& qadm ) 
{
  fDigitsArray = qadm.fDigitsArray;
  fHitsArray = qadm.fHitsArray;
  
  return *this;
}
//_____________________________________________________________________
AliFMDQADataMakerSim::~AliFMDQADataMakerSim()
{
  delete fDigitsArray;
  delete fHitsArray;
}

//_____________________________________________________________________
void AliFMDQADataMakerSim::EndOfDetectorCycle(AliQA::TASKINDEX_t task, 
					      TObjArray * list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliLog::Message(5,"FMD: end of detector cycle",
		  "AliFMDQADataMakerSim","AliFMDQADataMakerSim",
		  "AliFMDQADataMakerSim::EndOfDetectorCycle",
		  "AliFMDQADataMakerSim.cxx",83);
  AliQAChecker::Instance()->Run(AliQA::kFMD, task, list) ;  
  
}

//____________________________________________________________________ 
void AliFMDQADataMakerSim::InitHits()
{
  // create Digits histograms in Digits subdir
  TH1F* hEnergyOfHits = new TH1F("hEnergyOfHits","Energy distribution",100,0,3);
  hEnergyOfHits->SetXTitle("Edep");
  hEnergyOfHits->SetYTitle("Counts");
  Add2HitsList(hEnergyOfHits, 0);
}

//_____________________________________________________________________
void AliFMDQADataMakerSim::InitDigits()
{
  // create Digits histograms in Digits subdir
  TH1I* hADCCounts = new TH1I("hADCCounts","Dist of ADC counts",1024,0,1024);
  hADCCounts->SetXTitle("ADC counts");
  Add2DigitsList(hADCCounts, 0);
}

//_____________________________________________________________________
void AliFMDQADataMakerSim::MakeHits(TClonesArray * hits)
{
  TIter next(hits);
  AliFMDHit * hit;
  while ((hit = static_cast<AliFMDHit *>(next()))) 
    GetHitsData(0)->Fill(hit->Edep()/hit->Length()*0.032);
}

//_____________________________________________________________________
void AliFMDQADataMakerSim::MakeHits(TTree * hitTree)
{
  // make QA data from Hit Tree
  
  TBranch * branch = hitTree->GetBranch("FMD") ;
  if (!branch) {
    AliWarning("FMD branch in Hit Tree not found") ; 
    return;
  }
  fHitsArray->Clear();
  branch->SetAddress(&fHitsArray) ;
  
  for (Int_t ientry = 0 ; ientry < branch->GetEntries() ; ientry++) {
    branch->GetEntry(ientry);
    MakeHits(fHitsArray);   //tmp); 
  } 	
  
}

//_____________________________________________________________________
void AliFMDQADataMakerSim::MakeDigits(TClonesArray * digits)
{
  // makes data from Digits
  if(!digits) return;

  for(Int_t i = 0 ; i < fDigitsArray->GetEntriesFast() ; i++) {
    //Raw ADC counts
    AliFMDDigit* digit = static_cast<AliFMDDigit*>(digits->At(i));
    GetDigitsData(0)->Fill(digit->Counts());
  }
}

//_____________________________________________________________________
void AliFMDQADataMakerSim::MakeDigits(TTree * digitTree)
{
  
  fDigitsArray->Clear();
  TBranch * branch = digitTree->GetBranch("FMD") ;
  if (!branch)    {
      AliWarning("FMD branch in Digit Tree not found") ; 
      return;
  } 
  branch->SetAddress(&fDigitsArray) ;
  branch->GetEntry(0) ; 
  MakeDigits(fDigitsArray) ; 
}

//_____________________________________________________________________ 
void AliFMDQADataMakerSim::StartOfDetectorCycle()
{
   
}
//_____________________________________________________________________ 
//
// EOF
//
