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
//...

// --- ROOT system ---
#include <TClass.h>
#include <TH1F.h> 
#include <TH1I.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliACORDEQAChecker.h"
#include "AliCDBEntry.h"
#include "AliQAManager.h"

ClassImp(AliACORDEQAChecker)

//__________________________________________________________________
Double_t * AliACORDEQAChecker::Check(AliQAv1::ALITASK_t /*index*/, TObjArray ** list, AliDetectorRecoParam * /*recoParam*/)
{

	Double_t * test = new Double_t[AliRecoParam::kNSpecies] ; 
  	Int_t * count   = new Int_t[AliRecoParam::kNSpecies] ; 
  	Double_t * acoTest = new Double_t[AliRecoParam::kNSpecies];
 

	// Look at the QAref data for ACORDE

	char * acoOCDBDir = Form("ACORDE/%s/%s",AliQAv1::GetRefOCDBDirName(),AliQAv1::GetRefDataDirName());
	AliCDBEntry *acoQARefDir = AliQAManager::QAManager()->Get(acoOCDBDir);

	// Check variables set to 0

  	for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
	{
    		test[specie] = 0.0 ; 
    		count[specie] = 0 ; 
		acoTest[specie] = 0.0;
  	}
  
  	for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
	{
    		if (list[specie]->GetEntries() == 0)
		{  
      			test[specie] = 1. ; // nothing to check
			acoTest[specie] = 1.;
    		}
    		else 
		{
      			TIter next(list[specie]) ; 
      			TH1 * hdata ;
      			while ( (hdata = dynamic_cast<TH1 *>(next())) ) 
			{
        			if (hdata) 
				{ 
          				Double_t rv = 0.0 ; 
          				if(hdata->GetEntries()>0) rv=1; 
          				AliDebug(AliQAv1::GetQADebugLevel(), Form("%s -> %f", hdata->GetName(), rv)) ; 
          				count[specie]++ ; 
          				test[specie] += rv ; 
					Double_t acoHitsNorm;
					if (hdata->GetMaximum()==1) acoHitsNorm = 1;
					else acoHitsNorm = (hdata->GetMaximum() - 0.50)/hdata->GetMaximum();
					// here we implement the second version for ACORDEQAChecker
					// by the moment we only compare the Mean between the QA histograms and the Reference data 
					if  (acoQARefDir)
					{
						//AliWarning("Using the QA Reference data for ACORDE !!!");
						Double_t acoHistChecked = CheckAcordeRefHits(list[specie],(TObjArray *)acoQARefDir->GetObject());
						if ( (acoHistChecked>0.75) && (acoHistChecked<=1) ) acoTest[specie] = 0.86;
						if ( (acoHistChecked>0.0020) && (acoHistChecked<=0.75) ) acoTest[specie] = 0.251;
						if ( (acoHistChecked>0.0) && (acoHistChecked<=0.0020) ) acoTest[specie] = 0.0010;
						if ( (acoHistChecked>-1.0) && (acoHistChecked<=0.0) ) acoTest[specie] = -0.5;
	
					}else
					{
						//AliWarning("Using the inner ACORDE QA Checker !!!");
						if ( (acoHitsNorm>0.40) && (acoHitsNorm<=1) ) acoTest[specie] = 0.86;
						if ( (acoHitsNorm>0.0020) && (acoHitsNorm<=0.40) ) acoTest[specie] = 0.251;
						if ( (acoHitsNorm>0.0) && (acoHitsNorm<=0.0020) ) acoTest[specie] = 0.0010;
						if ( (acoHitsNorm>-1.0) && (acoHitsNorm<=0.0) ) acoTest[specie] = -0.5;
					}
        			}
        			else
				{
          				AliError("Data type cannot be processed") ;
        			}
      			}
      			if (count[specie] != 0) 
			{ 
        			if (test[specie]==0) 
				{
          				test[specie] = 0.5;  //upper limit value to set kWARNING flag for a task
        			}
        			else 
				{
					if (acoQARefDir) test[specie] = acoTest[specie];
					else
					{
						test[specie] = acoTest[specie];
					}
        			}
      			}
    		}
  	}
  	return test ; 
}
Double_t AliACORDEQAChecker::CheckAcordeRefHits(TObjArray *HistAcordeList, TObjArray *AcordeRef) const
{
	Double_t acordeTest = 0;
	TIter next(AcordeRef);
	TIter next1(HistAcordeList);
	TH1 *histoAcordeRef;
	TH1 *histoAcorde;
	Float_t acordeHistoQAMaker=0;
	Float_t meanACOQAReference=0;
	Float_t meanACOQAMaker=0;
	Float_t test1ACORDE = 0;
	while((histoAcordeRef=(TH1*)next()) && (histoAcorde=(TH1*)next1())) 
	{
		for(Int_t i=0;i<60;i++) acordeHistoQAMaker=acordeHistoQAMaker + histoAcorde->GetBinContent(i)/histoAcorde->GetMaximum();
		meanACOQAReference = histoAcordeRef->GetMean();
		meanACOQAMaker = acordeHistoQAMaker/60;
		test1ACORDE = TMath::Abs(meanACOQAReference-meanACOQAMaker);
		if (test1ACORDE<0.45) acordeTest = 0.86;
		if (test1ACORDE > 0.45) acordeTest = 0.50;
		if (test1ACORDE > 0.70) acordeTest = 0.25;
	}
	return acordeTest;
}
