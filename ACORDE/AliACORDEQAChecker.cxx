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
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliACORDEQAChecker.h"
#include "AliCDBEntry.h"
#include "AliQAManager.h"

ClassImp(AliACORDEQAChecker)

//____________________________________________________________________________
Double_t * AliACORDEQAChecker::Check(AliQA::ALITASK_t /*index*/)
{
  Double_t * rv = new Double_t[AliRecoParam::kNSpecies] ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    rv[specie] = 0.0 ; 
  return rv ;  
}

//__________________________________________________________________
Double_t * AliACORDEQAChecker::Check(AliQA::ALITASK_t /*index*/, TObjArray ** list)
{

	// We added one check to the ACORDE's QA histograms:
	// 1.- We check if they are empty
	// we check for the reference histogram to start the QAChecker. If not QAref object
	// is found, we check that the number of hits per channel is not so far from
	// the maximum number of hits.
  Double_t * test = new Double_t[AliRecoParam::kNSpecies] ; 
  Int_t * count   = new Int_t[AliRecoParam::kNSpecies] ; 
  Double_t * acoTest = new Double_t[AliRecoParam::kNSpecies];
  Double_t acoHitsNorm = 0;
 Double_t * acoRefTest = new Double_t[AliRecoParam::kNSpecies];

	// Look at the QAref data for ACORDE

	char * acoOCDBDir = Form("ACORDE/%s/%s",AliQA::GetRefOCDBDirName(),AliQA::GetRefDataDirName());
	AliCDBEntry *acoQARefDir = AliQAManager::QAManager()->Get(acoOCDBDir);


  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    test[specie]    = 0.0 ; 
    count[specie] = 0 ; 
	acoTest[specie] = 0.0;
  }
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (list[specie]->GetEntries() == 0){  
      test[specie] = 1. ; // nothing to check
	acoTest[specie] = 1.;
    }
    else {
      TIter next(list[specie]) ; 
      TH1 * hdata ;
      while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
        if (hdata) { 
          Double_t rv = 0.0 ; 
          if(hdata->GetEntries()>0)rv=1; 
          AliInfo(Form("%s -> %f", hdata->GetName(), rv)) ; 
          count[specie]++ ; 
          test[specie] += rv ; 
	

	// here we implement the second version for ACORDEQAChecker
	// by the moment we only compare that the hits in every ACORDE's channel
	// are close and > 0 
	for (Int_t i=0;i<60;i++)
	{
		acoHitsNorm =  hdata->GetBinContent(i)/hdata->GetMaximum();
		if  (acoQARefDir)
		{
	//		AliWarning("Using the QA Reference data for ACORDE !!!");
			test[specie] = CheckAcordeRefHits(list[specie],(TObjArray *)acoQARefDir->GetObject());
			if ((test[specie] = 0.86) || (acoHitsNorm>0.50)) 
			{
				acoRefTest[specie]=0.78;//printf("testMario: %f\n",acoRefTest[specie]);printf("histo:%f\n",hdata->GetMaximum());
			}
		}else{
	//	AliWarning("Using the inner ACORDE QA Checker !!!");
		if ( (acoHitsNorm>0.40) && (acoHitsNorm<=1) ) acoTest[specie] = 0.75;
		if ( (acoHitsNorm>0.0020) && (acoHitsNorm<=0.40) ) acoTest[specie] = 0.251;
		if ( (acoHitsNorm>0.0) && (acoHitsNorm<=0.0020) ) acoTest[specie] = 0.0010;
		if ( (acoHitsNorm>-1.0) && (acoHitsNorm<=0.0) ) acoTest[specie] = -0.5;
		}
	}
        }
        else{
          AliError("Data type cannot be processed") ;
        }
      }
      if (count[specie] != 0) { 
        if (test[specie]==0) {
         // AliWarning("Histograms are there, but they are all empty: setting flag to kWARNING");
          test[specie] = 0.5;  //upper limit value to set kWARNING flag for a task
        }
        else {
	if (acoQARefDir) test[specie] = acoRefTest[specie];
	else{
	test[specie] = acoTest[specie];//printf("testDyMa: %f\n",test[specie]);
	}
        }
      }
    }
   // AliInfo(Form("Test Result = %f", test[specie])) ; 
  }
  return test ; 
}
Double_t AliACORDEQAChecker::CheckAcordeRefHits(TObjArray *AcordeList, TObjArray *AcordeRef) const
{
	Double_t acoTest = 0;
	TIter next(AcordeList);
	TH1 *histo;
	for (Int_t i=0;i<60;i++)
	{
		while ( (histo = dynamic_cast<TH1 *>(next())) )
		{	
			if ( (histo->GetBinContent(i)/histo->GetMaximum())<1.0 ) acoTest = 0.86;
//		if( histo->KolmogorovTest((TH1F *)AcordeRef->At(0))<0.8)  acoTest = 0.86;
			//printf("href:%f\n",histo->GetMaximum());
		}
	}	
	return acoTest;
}
