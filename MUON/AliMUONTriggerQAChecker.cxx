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

  AliMUONVQAChecker::ECheckCode * rv = new AliMUONVQAChecker::ECheckCode[AliRecoParam::kNSpecies];

  Int_t histoRawsPercentIndex[] = {
    AliMUONQAIndices::kTriggerErrorSummaryNorm, 
    AliMUONQAIndices::kTriggerCalibSummaryNorm,
    AliMUONQAIndices::kTriggerReadOutErrorsNorm
  };
  const Int_t kNrawsHistos = sizeof(histoRawsPercentIndex)/sizeof(histoRawsPercentIndex[0]);

  // BEGIN OF LIMITS
  // Fixme: Move me to reference histos
  Float_t safeFactor = 5.;
  Float_t warningPercentTrigAlgo[AliMUONQAIndices::kNtrigAlgoErrorBins] = {safeFactor*1., safeFactor*1., safeFactor*1., 100., 100., 100., 100., safeFactor*1., safeFactor*1., safeFactor*1.};
  Float_t warningPercentCalib[AliMUONQAIndices::kNtrigCalibSummaryBins] = {safeFactor*0.4, safeFactor*1., 6.2, 0.0001, safeFactor*0.4};
  Float_t warningPercentReadout[AliMUONQAIndices::kNtrigStructErrorBins] = {safeFactor*1., safeFactor*1., safeFactor*1., safeFactor*1.};

  Float_t* warningPercent[kNrawsHistos] = {warningPercentTrigAlgo, warningPercentCalib, warningPercentReadout};
  
  Float_t errorFactor = 30.;
  Float_t errorPercentTrigAlgo[AliMUONQAIndices::kNtrigAlgoErrorBins] = {errorFactor*1., errorFactor*1., errorFactor*1., 100., 100., 100., 100., errorFactor*1., errorFactor*1., errorFactor*1.};
  Float_t errorPercentCalib[AliMUONQAIndices::kNtrigCalibSummaryBins] = {errorFactor*0.4, errorFactor*1., 3.*6.2, 3.*0.0001, errorFactor*0.4};
  Float_t errorPercentReadout[AliMUONQAIndices::kNtrigStructErrorBins] = {errorFactor*1., errorFactor*1., errorFactor*1., errorFactor*1.};
  // END OF LIMTS

  Float_t* errorPercent[kNrawsHistos] = {errorPercentTrigAlgo, errorPercentCalib, errorPercentReadout};
  
  TObjArray messages;
  messages.SetOwner(kTRUE);

  TH1* currHisto = 0x0;
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    rv[specie] = AliMUONVQAChecker::kInfo;

    TH1* hAnalyzedEvents = AliQAv1::GetData(list,AliMUONQAIndices::kTriggerRawNAnalyzedEvents,AliRecoParam::ConvertIndex(specie));
    Int_t nAnalyzedEvents = 0;
    if ( hAnalyzedEvents ) 
      nAnalyzedEvents = TMath::Nint(hAnalyzedEvents->GetBinContent(1));

    //if ( nAnalyzedEvents == 0 ) rv[specie] = AliMUONVQAChecker::kFatal;

    for(Int_t ihisto = 0; ihisto<kNrawsHistos; ihisto++){
      AliMUONVQAChecker::ECheckCode currRv = AliMUONVQAChecker::kInfo;
      messages.Clear();
      currHisto = AliQAv1::GetData(list,histoRawsPercentIndex[ihisto],AliRecoParam::ConvertIndex(specie));
      if ( ! currHisto ) continue;

      Int_t nbins = currHisto->GetXaxis()->GetNbins();
      for (Int_t ibin = 1; ibin<=nbins; ibin++){
        Double_t binContent = currHisto->GetBinContent(ibin);
        if ( binContent > errorPercent[ihisto][ibin-1] )
          currRv = AliMUONVQAChecker::kError;
        else if ( binContent > warningPercent[ihisto][ibin-1] )
          currRv = AliMUONVQAChecker::kWarning;
        else if ( ibin == 4 && binContent > 50. && AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCalib) {
          messages.Add(new TObjString("Do not panic:")); 
          messages.Add(new TObjString("copy errors do not affect data"));
        }
      } // loop on bins
      if ( currRv != AliMUONVQAChecker::kInfo ) {
	switch ( histoRawsPercentIndex[ihisto] ) {
	case AliMUONQAIndices::kTriggerErrorSummaryNorm:
	case AliMUONQAIndices::kTriggerCalibSummaryNorm:
	  messages.Add(new TObjString("Trigger algo errors"));
	  break;
	case AliMUONQAIndices::kTriggerReadOutErrorsNorm:
	  messages.Add(new TObjString("Readout errors"));
	  break;
	}
	if ( currRv == AliMUONVQAChecker::kWarning )
	  messages.Add(new TObjString("are a little bit high"));
	else if ( currRv == AliMUONVQAChecker::kError || 
		  currRv == AliMUONVQAChecker::kFatal )
	  messages.Add(new TObjString("are too high"));
      }
      else if ( nAnalyzedEvents != 0 ) {
	switch ( histoRawsPercentIndex[ihisto] ) {
	case AliMUONQAIndices::kTriggerErrorSummaryNorm:
	case AliMUONQAIndices::kTriggerCalibSummaryNorm:
	case AliMUONQAIndices::kTriggerReadOutErrorsNorm:
	  messages.Add(new TObjString("Values within limits"));
	  break;
	}
      }
      if ( MarkHisto(*currHisto, currRv) < rv[specie] )
	rv[specie] = currRv;
      currHisto->GetYaxis()->SetRangeUser(0., 110.);
      SetupHisto(nAnalyzedEvents, messages, *currHisto, currRv);
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


//___________________________________________________________________ 
void AliMUONTriggerQAChecker::SetupHisto(Int_t nevents, const TObjArray& messages, TH1& histo, AliMUONVQAChecker::ECheckCode code)
{
  //
  /// Add text to histos
  //

  Double_t y1 = 0.87 - (messages.GetLast()+2)*0.075;
  TPaveText* text = new TPaveText(0.25,y1,0.75,0.89,"NDC");
    
  text->AddText(Form("MTR - Total events %i", nevents));
    
  TIter next(&messages);
  TObjString* str;
    
  while ( ( str = static_cast<TObjString*>(next()) ) ){
    text->AddText(str->String());
  }

  TString defaultText = "";  
  Int_t color = 0;
  
  if ( nevents == 0 ) {
    color = AliMUONVQAChecker::kFatalColor;
    defaultText = "No events analyzed!";
  }
  else if ( nevents <= 20 ) {
    color = AliMUONVQAChecker::kWarningColor;
    text->AddText("Not enough events to judge");
    text->AddText("Please wait for more statistics");
  }
  else {
    switch ( code ) {
      case AliMUONVQAChecker::kInfo:
        color = AliMUONVQAChecker::kInfoColor;
        defaultText = "All is fine!";
        break;
      case AliMUONVQAChecker::kWarning:
        color = AliMUONVQAChecker::kWarningColor;
        defaultText = "Please keep an eye on it!";
        break;
      case AliMUONVQAChecker::kFatal:
        color = AliMUONVQAChecker::kFatalColor;
        defaultText = "PLEASE CALL MUON TRIGGER EXPERT!!!";
        break;
      default:
        color = AliMUONVQAChecker::kErrorColor;
        defaultText = "PLEASE CALL MUON TRIGGER EXPERT!";
        break;
    }
  }

  text->AddText(defaultText.Data());
  text->SetFillColor(color);
                      
  histo.SetFillStyle(1001);
  histo.SetFillColor(color);

  histo.SetStats(kFALSE);
    
  histo.GetListOfFunctions()->Delete();
  histo.GetListOfFunctions()->Add(text);
}
