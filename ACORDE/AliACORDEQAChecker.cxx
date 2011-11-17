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
//...
//  Checks the quality assurance for ACORDE. 
//  Default implementation
//  Authors:
//	Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch> (FCFM-BUAP) 
//	Luciano Diaz Gonzalez <luciano.diaz@nucleares.unam.mx> (ICN-UNAM)
//	Arturo Fernandez <afernan@mail.cern.ch> (FCFM-BUAP)
//  Last update: Nov. 14t 2009 --> MRC <mrodrigu@mail.cern.ch> (FCFM-BUAP) 
//...

// --- ROOT system ---
#include <TClass.h>
#include <TH1F.h> 
#include <TH1I.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 
#include <TPaveText.h>
#include <TLine.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliACORDEQAChecker.h"
#include "AliCDBEntry.h"
#include "AliQAManager.h"

/*************************************************************************
   Last update: Oct. 10th 2011 from Mario RC, <mrodrigu@mail.cern.ch>
	|-> Adding the checker class for raw and esd index
	|-> Setting the threshold lines and box for DQM shifter

*************************************************************************/

ClassImp(AliACORDEQAChecker)

//____________________________________________________________________________

AliACORDEQAChecker::AliACORDEQAChecker()  :
AliQACheckerBase("ACORDE","ACORDE Quality Assurance Data Maker"),
fTextDQMShifterInfo(new TPaveText(35,0.5,55,0.7,"T")),
fMax(new TLine(0,0.9,60,0.9)),
fMin(new TLine(0,0.3,60,0.3))
{
	// default constructor
	fMax->SetLineColor(kGreen);
	fMax->SetLineWidth(3);
	fMin->SetLineColor(kGreen);
	fMax->SetLineWidth(3);
}
//____________________________________________________________________________
AliACORDEQAChecker::~AliACORDEQAChecker()
{
	// destructor
	delete fTextDQMShifterInfo;
	delete fMax;
	delete fMin;
}
//____________________________________________________________________________
AliACORDEQAChecker::AliACORDEQAChecker(const AliACORDEQAChecker& qac) :
AliQACheckerBase(qac.GetName(), qac.GetTitle()),
fTextDQMShifterInfo(new TPaveText(35,0.5,55,0.7,"T")),
fMax(static_cast<TLine*>(qac.fMax->Clone())),
fMin(static_cast<TLine*>(qac.fMin->Clone()))
{
	//
}
//____________________________________________________________________________
AliACORDEQAChecker& AliACORDEQAChecker::operator = (const AliACORDEQAChecker &qac)
{
	
	fTextDQMShifterInfo = new TPaveText(35,0.5,55,0.7,"T");
	fMax = static_cast<TLine*>(qac.fMax->Clone());
	fMin = static_cast<TLine*>(qac.fMin->Clone());
	return *this;
}
//____________________________________________________________________________
void AliACORDEQAChecker::Check(Double_t * test, AliQAv1::ALITASK_t /*index*/)
{
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    test[specie] = 0.0 ; 
}
//____________________________________________________________________________
void AliACORDEQAChecker::Check(Double_t * test, AliQAv1::ALITASK_t /*index*/, TObjArray ** list, const AliDetectorRecoParam * /*recoParam*/)
{
// Close version to the final one for the ACORDE QA Checker
  
// Loop over the run species (for specie!= cosmic by now we set QA to INFO) 
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
  {
	if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) continue ; 
	if (list[specie]->GetEntries() == 0) test[specie] = 1.; // Nothing to check
	else 
	{
		TIter next(list[specie]) ; 
		TH1 * hdata ; // Data created by the AliACORDEQADataMakerXXX (Sim/Rec)
		while ( (hdata = dynamic_cast<TH1 *>(next())) ) 
		{
			if (hdata) 
			{ 
				Double_t rv = 0.0 ; 
				if(hdata->GetEntries()>0) rv=1; 
				AliDebug(AliQAv1::GetQADebugLevel(), Form("%s -> %f", hdata->GetName(), rv)) ; 
				TString hdataName = hdata->GetName();
				if (hdata->GetListOfFunctions()->GetEntries() == 0  && (hdataName.Contains("HitsSL0_DQM")||(hdataName.Contains("HitsAMU_DQM")))){
					hdata->GetListOfFunctions()->Add(fTextDQMShifterInfo);
                                        hdata->GetListOfFunctions()->Add(fMax);
                                        hdata->GetListOfFunctions()->Add(fMin);

				}
				// Here we use the QAref ACORDE data from fRefOCDBSubDir
				
				if ( (fRefOCDBSubDir[specie]) && (hdataName.Contains("ACORDEBitPattern")) ) 
				{
					TH1 * href = NULL;
					if (fRefSubDir) href = static_cast<TH1*>(fRefSubDir->Get(hdata->GetName()));
					else if (fRefOCDBSubDir[specie]) href = static_cast<TH1*>(fRefOCDBSubDir[specie]->FindObject(hdata->GetName()));
					test[specie] = CheckAcordeRefHits(href,hdata);
				}else if (hdataName.Contains("ACORDEBitPattern")||hdataName.Contains("HitsSL0_DQM")||hdataName.Contains("HitsAMU_DQM"))
				// Here we use an inner QA Checher without the QAref data
				{
					Float_t acoDataMax = hdata->GetMaximum();
					Int_t flagAcoQAChecker = 0;
					Int_t flagAcoQA = 0;
					if (acoDataMax == 0) continue;
					for(Int_t i=0;i<60;i++)
					{
						if ((hdata->GetBinContent(i)/acoDataMax) < 0.75) flagAcoQAChecker++; 
						if (hdataName.Contains("HitsSL0_DQM")||hdataName.Contains("HitsAMU_DQM")){
							if (hdata->GetBinContent(i)<0.29) flagAcoQA++;
						}
					}
					if (hdataName.Contains("ACORDEBitPattern")){
						Double_t simpleFlag = 1.-flagAcoQAChecker/60.;
						if ( (simpleFlag >= 0.90) && (simpleFlag <= 1.0) ) test[specie] = 0.75; // INFO
						if ( (simpleFlag >= 0.70) && (simpleFlag < 0.90) ) test[specie] = 0.50; // WARNING
						if ( (simpleFlag >= 0.25) && (simpleFlag < 0.70) ) test[specie] = 0.25; // ERROR
						if ( (simpleFlag >= 0.0) && (simpleFlag < 0.25) )  test[specie] = -1.0; // FATAL
					} else if (hdataName.Contains("HitsSL0_DQM")||hdataName.Contains("HitsAMU_DQM")){
						if (flagAcoQA < 100) test[specie] = 0.75;
						else test[specie] = 0.3;
					}


				}	

				// Setting Warning message for possible Empty Events with the ACORDE-Trigger
					
				if (hdataName.Contains("ACORDEMultiplicity") && (hdata->GetBinContent(0)!=0)) AliWarning("Empty event registered with ACORDE Trigger !!!");
					
				
			}else AliError("Data type cannot be processed") ;
			if (fTextDQMShifterInfo)
			{
				fTextDQMShifterInfo->Clear();
				if (test[specie] > 0.3){
					fTextDQMShifterInfo->SetFillColor(kGreen);
					fTextDQMShifterInfo->AddText("ACORDE: O.K.");
				}else{
					fTextDQMShifterInfo->SetFillColor(kRed);
					fTextDQMShifterInfo->AddText("ACORDE: Not, O.K.");
					fTextDQMShifterInfo->AddText("CALL THE EXPERTS");
				}
			}
        			
		}
	}
    if ( (specie == AliRecoParam::kHighMult) || (specie == AliRecoParam::kLowMult) || (specie == AliRecoParam::kCalib) ) test[specie] = 0.75;
  }
}
//____________________________________________________________________________
Double_t AliACORDEQAChecker::CheckAcordeRefHits(const TH1 * href, const TH1 * hdata) const
{
	Double_t test = 0.;
	Int_t flag=0;
	for (Int_t i=0;i<60;i++)
	{
		if (TMath::Abs(href->GetBinContent(i)-hdata->GetBinContent(i))>10) flag++;
	}
	if ((flag>50)&&(flag<=60)) test = -1.;
	if ((flag>30)&&(flag<=50)) test = 0.25;
	if ((flag>10)&&(flag<=30)) test = 0.5;
	if ((flag>0)&&(flag<=10)) test = 0.75;
	return test;
}
