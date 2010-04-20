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

// $Id$

#include "AliMUONTriggerQAChecker.h"

/// \class AliMUONTriggerQAChecker
///
/// Implementation of QAChecker for MTR
///
/// For the moment we only implement the checking of raw data QA for the trigger
/// by looking at the local structure and trigger response errors.
///
/// \author Diego Stocco, Subatech


#include "AliRecoParam.h"
#include "AliMUONQAIndices.h"
#include "AliQAv1.h"
#include "TH1.h"
#include "TPaveText.h"
#include "TString.h"
#include "TObjArray.h"
#include "TList.h"

/// \cond CLASSIMP
ClassImp(AliMUONTriggerQAChecker)
/// \endcond

//__________________________________________________________________
AliMUONTriggerQAChecker::AliMUONTriggerQAChecker() : AliMUONVQAChecker()
{
	/// ctor
}          

//__________________________________________________________________
AliMUONTriggerQAChecker::~AliMUONTriggerQAChecker() 
{
	/// dtor
}

//______________________________________________________________________________
AliMUONVQAChecker::ECheckCode 
AliMUONTriggerQAChecker::MarkHisto(TH1& histo, AliMUONVQAChecker::ECheckCode value) const
{
  /// Mark histo as originator of some QA error/warning
  
  if ( value != AliMUONVQAChecker::kInfo )
  {
    histo.SetBit(AliQAv1::GetQABit());
  }
  
  return value;
}

//__________________________________________________________________
AliMUONVQAChecker::ECheckCode* 
AliMUONTriggerQAChecker::CheckRaws(TObjArray** list, const AliMUONRecoParam* )
{
  /// Check raw data

  AliMUONVQAChecker::ECheckCode * rv = new AliMUONVQAChecker::ECheckCode[AliRecoParam::kNSpecies] ; 

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    rv[specie] = AliMUONVQAChecker::kInfo; 
  }

  Int_t histoRawsPercentIndex[] = {
    AliMUONQAIndices::kTriggerErrorSummaryNorm, 
    AliMUONQAIndices::kTriggerCalibSummaryNorm,
    AliMUONQAIndices::kTriggerReadOutErrorsNorm
  };
  const Int_t kNrawsHistos = sizeof(histoRawsPercentIndex)/sizeof(histoRawsPercentIndex[0]);

  // MOVE THESE TO REFERENCE HISTOS
// START WITH THOSE COMMENTED OUT UNTIL WE GAIN CONFIDENCE...
//  Float_t safeFactor = 5.;
//  Float_t alarmPercentTrigAlgo[AliMUONTriggerQADataMakerRec::kNtrigAlgoErrorBins] = {safeFactor*1., safeFactor*1., safeFactor*1., 100., 100., 100., 100., safeFactor*1., safeFactor*1., safeFactor*1.};
//  Float_t alarmPercentCalib[AliMUONTriggerQADataMakerRec::kNtrigCalibSummaryBins] = {safeFactor*0.4, safeFactor*1., 6.2, 0.0001, safeFactor*0.4};
//  Float_t alarmPercentReadout[AliMUONTriggerQADataMakerRec::kNtrigStructErrorBins] = {safeFactor*1., safeFactor*1., safeFactor*1., safeFactor*1.};
//
//  Float_t* alarmPercent[kNrawsHistos] = {alarmPercentTrigAlgo, alarmPercentCalib, alarmPercentReadout};
// END OF COWARD COMMENTING...
  
  TH1* currHisto = 0x0;
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    TH1* hAnalyzedEvents = AliQAv1::GetData(list,AliMUONQAIndices::kTriggerRawNAnalyzedEvents,AliRecoParam::ConvertIndex(specie));
    Int_t nAnalyzedEvents = 0;
    if ( hAnalyzedEvents ) 
      nAnalyzedEvents = hAnalyzedEvents->GetBinContent(1);
    for(Int_t ihisto = 0; ihisto<kNrawsHistos; ihisto++){
      currHisto = AliQAv1::GetData(list,histoRawsPercentIndex[ihisto],AliRecoParam::ConvertIndex(specie));
      if ( currHisto ){
	currHisto->SetBarWidth(0.5);
	currHisto->SetBarOffset(0.25);
	TPaveText* text = new TPaveText(0.65,0.65,0.99,0.99,"NDC");
	TString binName;
	Bool_t isOk = kTRUE;
	Int_t nbins = currHisto->GetXaxis()->GetNbins();
	for (Int_t ibin = 1; ibin<=nbins; ibin++){
	  binName = currHisto->GetXaxis()->GetBinLabel(ibin);
	  binName.ReplaceAll("#splitline","");
	  binName.ReplaceAll("{","");
	  binName.ReplaceAll("}","");
	  Float_t binContent = currHisto->GetBinContent(ibin);
//	  if (binContent > alarmPercent[ihisto][ibin-1]) isOk = kFALSE;
	  text->AddText(Form("%5.2f %% in %s", binContent, binName.Data()));
	  //text->AddText(Form("%5.2f %% in %s (limit %5.2f %%)", binContent, binName.Data(), alarmPercent[ihisto][ibin-1]));
	}
	text->AddText(Form("Total events %i", nAnalyzedEvents));
	if ( ! isOk || nAnalyzedEvents == 0 ) {
	  text->SetFillColor(kRed);
	  rv[specie] = MarkHisto(*currHisto, AliMUONVQAChecker::kError);
	}
	else text->SetFillColor(kGreen);
	currHisto->GetListOfFunctions()->Add(text);
	currHisto->SetStats(kFALSE);
	currHisto->GetYaxis()->SetRangeUser(0., 110.);
      }
    } // loop on histos
  } // loop on species

  return rv;
}

//__________________________________________________________________
AliMUONVQAChecker::ECheckCode* 
AliMUONTriggerQAChecker::CheckRecPoints(TObjArray** , const AliMUONRecoParam* )
{
  /// Check rec points
  return 0x0;
}

//__________________________________________________________________
AliMUONVQAChecker::ECheckCode* 
AliMUONTriggerQAChecker::CheckESD(TObjArray** , const AliMUONRecoParam* )
{
  /// Check esd
  return 0x0;
}
