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



//---
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.

//  Authors:
//
//  Luciano Diaz Gonzalez <luciano.diaz@nucleares.unam.mx> (ICN-UNAM)
//  Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch> (FCFM-BUAP)
//  Arturo Fernandez Tellez <afernan@mail.cern.ch (FCFM-BUAP)
//
//  Created: June 13th 2008
//---
// Last Update: Aug. 27th 2008 ---> Implementation to declare QA expert histogram 


// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH2F.h>
#include <TDirectory.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliACORDEdigit.h" 
#include "AliACORDEhit.h"
#include "AliACORDEQADataMakerSim.h"
#include "AliQAChecker.h"
#include "AliACORDERawReader.h"
ClassImp(AliACORDEQADataMakerSim)
           
//____________________________________________________________________________ 
AliACORDEQADataMakerSim::AliACORDEQADataMakerSim():AliQADataMakerSim(AliQAv1::GetDetName(AliQAv1::kACORDE), "ACORDE Quality Assurance Data Maker")
{
}
//____________________________________________________________________________ 
AliACORDEQADataMakerSim::AliACORDEQADataMakerSim(const AliACORDEQADataMakerSim& qadm) :
  AliQADataMakerSim() 
{
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}
//__________________________________________________________________
AliACORDEQADataMakerSim& AliACORDEQADataMakerSim::operator = (const AliACORDEQADataMakerSim& qadm )
{
  // Equal operator.
  this->~AliACORDEQADataMakerSim();
  new(this) AliACORDEQADataMakerSim(qadm);
  return *this;
}
//____________________________________________________________________________
void AliACORDEQADataMakerSim::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  ResetEventTrigClasses(); // reset triggers list to select all histos
  AliDebug(AliQAv1::GetQADebugLevel(), "ACORDE---->Detector specific actions at END of cycle\n................\n");

  AliQAChecker::Instance()->Run(AliQAv1::kACORDE, task, list) ;
}
//____________________________________________________________________________
void AliACORDEQADataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  AliDebug(AliQAv1::GetQADebugLevel(), "ACORDE---->Detector specific actions at START of cycle\n................\n");
}
//____________________________________________________________________________ 
void AliACORDEQADataMakerSim::InitHits()
{
  // create Hits histograms in Hits subdir
  
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F * fHitsACORDE;  
  
  fHitsACORDE = new TH1F("ACORDEBitPatternfromHits","Distribution of ACORDE fired modules from HITS",60,1,60);
  Add2HitsList(fHitsACORDE,0,!expert,image);
  
  const char *acoModule[60]={"0_0","0_1","0_2","0_3","0_4","0_5","0_6","0_7","0_8","0_9",
			     "1_0","1_1","1_2","1_3","1_4","1_5","1_6","1_7","1_8","1_9",
			     "2_0","2_1","2_2","2_3","2_4","2_5","2_6","2_7","2_8","2_9",
			     "3_0","3_1","3_2","3_3","3_4","3_5","3_6","3_7","3_8","3_9",
			     "4_0","4_1","4_2","4_3","4_4","4_5","4_6","4_7","4_8","4_9",
			     "5_0","5_1","5_2","5_3","5_4","5_5","5_6","5_7","5_8","5_9"};
  
  
  fHitsACORDE->SetXTitle("Modules");
  fHitsACORDE->SetYTitle("Counts");
  for (Int_t i=0;i<60;i++)
    {
      fHitsACORDE->GetXaxis()->SetBinLabel(i+1,acoModule[i]);
    }
  //
  //
  ClonePerTrigClass(AliQAv1::kHITS); // this should be the last line
}

//____________________________________________________________________________ 
void AliACORDEQADataMakerSim::InitDigits()
{
  // create Digits histograms in Digits subdir

  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F *    fhDigitsModule;
  fhDigitsModule = new TH1F("ACORDEBitPatternfromDigits","Distribution of ACORDE fired modules from DIGITS",60,1,60);
  Add2DigitsList(fhDigitsModule,0,!expert,image);
  const char *acoModule[60]={"0_0","0_1","0_2","0_3","0_4","0_5","0_6","0_7","0_8","0_9",
                        "1_0","1_1","1_2","1_3","1_4","1_5","1_6","1_7","1_8","1_9",
                        "2_0","2_1","2_2","2_3","2_4","2_5","2_6","2_7","2_8","2_9",
                        "3_0","3_1","3_2","3_3","3_4","3_5","3_6","3_7","3_8","3_9",
                        "4_0","4_1","4_2","4_3","4_4","4_5","4_6","4_7","4_8","4_9",
                        "5_0","5_1","5_2","5_3","5_4","5_5","5_6","5_7","5_8","5_9"};
  
  
  fhDigitsModule->SetXTitle("Modules");
  fhDigitsModule->SetYTitle("Counts");
  for (Int_t i=0;i<60;i++)
    {
      fhDigitsModule->GetXaxis()->SetBinLabel(i+1,acoModule[i]);
    }
    //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}

//____________________________________________________________________________
void AliACORDEQADataMakerSim::MakeHits(TTree *hitTree)
{
  // Here we fill the QA histos for Hits declared above
  IncEvCountCycleHits();
  IncEvCountTotalHits();

	if (fHitsArray) 
    	fHitsArray->Clear() ; 
  	else
    	fHitsArray = new TClonesArray("AliACORDEhit",1000);
	TBranch * branch = hitTree->GetBranch("ACORDE");
	if (!branch) 
	{
		AliWarning("ACORDE branch in Hit Tree not found");
	} else 
	{
    		branch->SetAddress(&fHitsArray);
		for(Int_t track = 0 ; track < branch->GetEntries() ; track++) 
		{
			branch->GetEntry(track);
			for(Int_t ihit=0 ; ihit < fHitsArray->GetEntriesFast() ; ihit++) 
			{
				AliACORDEhit *AcoHit = (AliACORDEhit*) fHitsArray->UncheckedAt(ihit);
				if(!AcoHit) 
				{
					AliError("The unchecked hit doesn't exist");
					continue ;
				}
				FillHitsData(0,AcoHit->GetModule());
      			}			
    		}
  	}

}
//____________________________________________________________________________
void AliACORDEQADataMakerSim::MakeDigits( TTree *digitsTree)
{
  //fills QA histos for Digits
  IncEvCountCycleDigits();
  IncEvCountTotalDigits();

  if (fDigitsArray) 
    fDigitsArray->Clear() ; 
  else
    fDigitsArray = new TClonesArray("AliACORDEdigit",1000);
  
  TBranch * branch = digitsTree->GetBranch("ACORDEdigit");
  if (!branch) {
    AliWarning("ACORDE branch in Digits Tree not found");
  } else {
   branch->SetAddress(&fDigitsArray);
    for(Int_t track = 0 ; track < branch->GetEntries() ; track++) {
      branch->GetEntry(track);
      for(Int_t idigit = 0 ; idigit < fDigitsArray->GetEntriesFast() ; idigit++) {
        AliACORDEdigit *AcoDigit = (AliACORDEdigit*) fDigitsArray->UncheckedAt(idigit);
        if (!AcoDigit) {
          AliError("The unchecked digit doesn't exist");
          continue ;
        }
        FillDigitsData(0,AcoDigit->GetModule());
      }
    }
  }
}
