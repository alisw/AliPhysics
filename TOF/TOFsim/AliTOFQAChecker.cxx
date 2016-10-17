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

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  Checks the quality assurance.                                  //
//  By analysis of the histograms & comparing with reference data  //
//  F. Bellini last modified on 17.10.2016                         //
//                                                                 //
/////////////////////////////////////////////////////////////////////


#include "TList.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TString.h" 
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliQADataMakerRec.h"
#include "AliTOFQAChecker.h"
#include <TPaveText.h>


ClassImp(AliTOFQAChecker)
//____________________________________________________________________________
void AliTOFQAChecker::Check(Double_t * test, AliQAv1::ALITASK_t /*index*/,
				  TObjArray ** list,
				  const AliDetectorRecoParam * /*recoParam*/) 
{
  // Super-basic check on the QA histograms on the input list: 
  // look whether they are empty!
  Float_t rvarr[2][4] = {1.0, 0.75, 0.5, 0.25,
			 1.0, 0.5, 0.002, 0.0};
  AliInfo(Form("Checking only golden plots? %s", (fCheckGolden? "YES":"NO")));
  Int_t count[AliRecoParam::kNSpecies] = { 0 }; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (! AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::ConvertIndex(specie)))
      continue ;
    test[specie] = 1.0 ;
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue ;
    if (!list || (list[specie]->GetEntries() == 0)){
      // nothing to check, setting to FATAL
      test[specie] = -1.0 ;
    } else {
      TIter next(list[specie]) ;
      TH1 * hdata ;
      count[specie] = 0 ;
      while ( (hdata = static_cast<TH1 *>(next())) ) {
	if (!hdata || !hdata->InheritsFrom("TH1")) {
	  AliError("Data type cannot be processed");
	} else {
	  TString histname = hdata->GetName();
	  Bool_t isGoldenHisto = (histname.Contains("TOFRawsMulti") || histname.Contains("TOFRawsTime") || histname.Contains("hTOFRawsToT"));
	  if (fCheckGolden && !isGoldenHisto) continue;
	  Double_t rv = 0.;
	  Int_t checkOutput = CheckRaws(hdata,specie);
	  if (checkOutput>AliQAv1::kNULLBit && checkOutput<AliQAv1::kNBIT) rv = rvarr[fCheckGolden][checkOutput];
	  else rv = AliQAv1::kNULLBit;
	  AliDebug(AliQAv1::GetQADebugLevel(), Form("%s -> %f", hdata->GetName(), rv));
	  count[specie]++;
	  test[specie] += rv;
	}   
      } //end while loop on the list
      if (count[specie]<=0) {
	//the only case in which this condition is verified at this stage is
	//when the golden plots are missing from the list and were requested,
	//therefore this generates a FATAL
	AliError("Error: cannot find requested golden histograms");
	test[specie] = -1.0;
      } else {
	test[specie] /= count[specie] ;
	AliDebug(AliQAv1::GetQADebugLevel(), Form("Test Result = %f", test[specie]));
      }
    } //END ELSE: is there anything to check?
  }  //END FOR: cycle over all possible species
}  

//------------------------------------------------------
AliTOFQAChecker& AliTOFQAChecker::operator = (const AliTOFQAChecker& qac)
{
  
  if (this==&qac) return *this;
  return *this;
  
}

//____________________________________________________________________________
Int_t AliTOFQAChecker::CheckRaws(TH1* histo, Int_t specie)
{
  /*
  checker for RAWS
  */
  Int_t nTrgCl=AliQADataMaker::GetNTrigClasses();
  Int_t nToDrawTrgCl;

  /*customize the summary image*/
  Bool_t drawRawsSumImage=kTRUE; //hTOF_Raws: kTRUE shows it in the summary image, kFALSE does not
  Bool_t drawRawsTimeSumImage=kTRUE; //hTOF_RawsTime: kTRUE shows it in the summary image, kFALSE does not
  Bool_t drawRawsToTSumImage=kTRUE; //hTOF_RawsToT: kTRUE shows it in the summary image, kFALSE does not
  TString ClToDraw[]={"kINT7", "kCalibBarell", "0"}; //trigger classes shown in the summary image (it MUST end with "0")
  for(nToDrawTrgCl=0; !ClToDraw[nToDrawTrgCl].EqualTo("0"); nToDrawTrgCl++) {}

  Int_t flag = AliQAv1::kNULLBit;
  Int_t trgCl;
  Int_t trigId=0;
  Bool_t suffixTrgCl=kFALSE;
  if(histo->GetEntries()>0) flag = AliQAv1::kINFO; 
  
  Double_t binWidthTOFrawTime = 2.44;
  Float_t minTOFrawTime, maxTOFrawTime;
  if (AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCosmic){
    minTOFrawTime=150.;//ns
    maxTOFrawTime=250.;//ns
  } else {
    minTOFrawTime=150.;//ns
    maxTOFrawTime=225.;//ns
  } 
  Float_t minTOFrawTot = 10.;
  Double_t maxTOFrawTot = 15.;

  Float_t minTOFrawINT7hits = 10;
  Float_t maxTOFrawINT7hits = 70;
 
  TString histname = histo->GetName();
  TPaveText * text = 0x0;
  for (trgCl=0; trgCl<nTrgCl; trgCl++) 
    if (histname.EndsWith(AliQADataMaker::GetTrigClassName(trgCl))) suffixTrgCl=kTRUE;
  if ((histname.EndsWith("TOFRawsMulti")) || (histname.Contains("TOFRawsMulti") && suffixTrgCl)) {
    if (!suffixTrgCl) histo->SetBit(AliQAv1::GetImageBit(), drawRawsSumImage);
    if (suffixTrgCl) {
      histo->SetBit(AliQAv1::GetImageBit(), kFALSE); //clones not shown by default
      for(int i=0; i<nToDrawTrgCl; i++) {
        if(histname.EndsWith(ClToDraw[i]))
          histo->SetBit(AliQAv1::GetImageBit(), kTRUE);
       }
    }
    text = (TPaveText *) histo->GetListOfFunctions()->FindObject("hitsMsg");
    if (histo->GetEntries()==0) {
        if (text){
            text->Clear();
            text->AddText("No entries. IF TOF IN RUN");
            text->AddText("check the TOF TWiki");
            text->SetFillColor(kYellow);
        }
      flag = AliQAv1::kWARNING;
    } else {
      Float_t multiMean = histo->GetMean();
      Float_t zeroBinIntegral = histo->Integral(1,1);
      Float_t lowMIntegral = histo->Integral(1,10);
      Float_t totIntegral = histo->Integral(2, histo->GetNbinsX());
      
      if (totIntegral==0){ //if only "0 hits per event" bin is filled -> error
	if (histo->GetBinContent(1)>0) {
        if (text){
	  text->Clear();
	  text->AddText("No TOF hits for all events."); 
	  text->AddText("Call TOF on-call."); 
	  text->SetFillColor(kRed);
        }
        flag = AliQAv1::kERROR;
	}
      } else { 
	if (AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCosmic) {
	  if (multiMean<10.){
          if (text){
	    text->Clear();
	    text->AddText(Form("Multiplicity within limits"));
	    text->AddText("for COSMICS: OK!!!");
	    text->SetFillColor(kGreen);
          }
          flag = AliQAv1::kINFO;
	  } else {
          if (text){
	    text->Clear();
	    text->AddText(Form("Multiplicity too high"));
	    text->AddText("for COSMICS: email TOF on-call");
	    text->SetFillColor(kYellow);
          }
          flag = AliQAv1::kWARNING;
	  }
	} else {
	  Bool_t isZeroBinContentHigh = kFALSE;
	  Bool_t isLowMultContentHigh = kFALSE;
	  Bool_t isINT7AverageLow = kFALSE;
	  Bool_t isINT7AverageHigh = kFALSE;
	  
	  if (zeroBinIntegral>0.75*totIntegral) isZeroBinContentHigh = kTRUE;
	  if (lowMIntegral>0.75*totIntegral) isLowMultContentHigh = kTRUE;
	  if (multiMean<minTOFrawINT7hits) isINT7AverageLow = kTRUE;
	  if (multiMean>maxTOFrawINT7hits) isINT7AverageHigh = kTRUE;
	  
	  if (AliRecoParam::ConvertIndex(specie) == AliRecoParam::kLowMult) {
	      if (isZeroBinContentHigh && (multiMean>10.)) {
		if (text){
		  text->Clear();
		  text->AddText(Form("Events with 0 hits = %5.2f%%", zeroBinIntegral*100./totIntegral));
		  text->AddText("OK for COSMICS and technical.");
		  text->AddText("For pp check trigger class and TWiki.");
		  text->SetFillColor(kYellow);
		}
	      } else {
		if (!histname.Contains("INT7") && (multiMean>100.)){
		  if (text){
		    text->Clear();
		    text->AddText(Form("Unexpected mean value = %5.2f",multiMean));
		    text->AddText("Check trigger class and TOF TWiki.");
		    text->SetFillColor(kYellow);
		  }
		  flag = AliQAv1::kWARNING;
		} else {
		  if ( histname.Contains("INT7") && (isINT7AverageLow || isINT7AverageHigh)){
		    if (text){
		      text->Clear();
		      text->AddText(Form("Unexpected mean value = %5.2f",multiMean));
		      text->AddText("Check trigger class and TOF TWiki.");
		      text->SetFillColor(kYellow);
		    }
		    flag = AliQAv1::kWARNING;
		  } else {
		    if (text){
		      text->Clear();
		      text->AddText(Form("Mean value = %5.2f",multiMean));
		      text->AddText(Form("Events with 0 hits = %5.2f%%",zeroBinIntegral*100./totIntegral));
		      text->AddText("OK!");
		      text->SetFillColor(kGreen);
		    }
		    flag = AliQAv1::kINFO;
		  }
		}
	      }
	  }
	  else if ( (AliRecoParam::ConvertIndex(specie) == AliRecoParam::kHighMult)
		    &&(isLowMultContentHigh || (multiMean>500.))){
	    //assume that good range of multi in PbPb goes from 20 to 500 tracks
	    if (text){
              text->Clear();
	      text->AddText(Form("Unexpected mean value for Pb-Pb = %5.2f",multiMean));
	      text->AddText("OK for technical.");
	      text->AddText("Check TOF TWiki for PbPb.");
	      text->SetFillColor(kYellow);
	    }
	    flag = AliQAv1::kWARNING;
	  } else {
	    if (text){
              text->Clear();
	      text->AddText(Form("Multiplicity within limits"));
	      text->AddText("    OK!!!    ");
	      text->SetFillColor(kGreen);
	    }
	    flag = AliQAv1::kINFO;
	  }
	}
      }
    }
  }
  if ((histname.EndsWith("RawsTime")) || (histname.Contains("RawsTime") && suffixTrgCl)) {
    if (!suffixTrgCl) histo->SetBit(AliQAv1::GetImageBit(), drawRawsTimeSumImage);
    if (suffixTrgCl) {
      histo->SetBit(AliQAv1::GetImageBit(), kFALSE);
      for(int i=0; i<nToDrawTrgCl; i++) {
        if(histname.EndsWith(ClToDraw[i]))
          histo->SetBit(AliQAv1::GetImageBit(), kTRUE);
       }
    }
    text = (TPaveText *) histo->GetListOfFunctions()->FindObject("timeMsg");
    if (histo->GetEntries()==0) {
        if (text){
      text->Clear();
      text->AddText("No entries. If TOF in the run"); 
      text->AddText("check TOF TWiki"); 
      text->SetFillColor(kYellow);
        }
        flag = AliQAv1::kWARNING;
    } else {
      Float_t timeMean = histo->GetMean();
      Int_t lowBinId = TMath::Nint(minTOFrawTime/binWidthTOFrawTime);
      Int_t highBinId = TMath::Nint(maxTOFrawTime/binWidthTOFrawTime);
      Float_t peakIntegral = histo->Integral(lowBinId,highBinId);
      Float_t totIntegral = histo->Integral(1, histo->GetNbinsX());      
      if ( (timeMean > minTOFrawTime) && (timeMean < maxTOFrawTime)) {
	flag = AliQAv1::kINFO;
	if (text){
	  text->Clear();
	  text->AddText("Mean inside limits: OK!!!"); 
	  text->SetFillColor(kGreen);
	}
      } else {
	if (peakIntegral/totIntegral > 0.20) { 
	  AliWarning(Form("Raw time: peak/total integral = %5.2f, mean = %5.2f ns -> Check filling scheme...",peakIntegral/totIntegral,timeMean));
	  if (text){
	    text->Clear();
	    text->AddText(Form("Raw time peak/total integral = %5.2f%%", peakIntegral*100./totIntegral));
	    text->AddText(Form("Mean = %5.2f ns",timeMean));
	    text->AddText("If multiple peaks,"); 
	    text->AddText("check filling scheme."); 
	    text->AddText("See TOF TWiki."); 
	    text->SetFillColor(kYellow);
	  }
	  flag = AliQAv1::kWARNING;
	} else {
	  AliWarning(Form("Raw time peak/total integral = %5.2f, mean = %5.2f ns", peakIntegral/totIntegral,timeMean));
	  if (text){
	    text->Clear();
	    text->AddText("Call TOF on-call."); 
	    text->AddText("Signal mean outside limits."); 
	    text->AddText(Form("Raw time peak/total integral = %5.2f%%", peakIntegral*100./totIntegral));
	    text->AddText(Form("Mean = %5.2f ns",timeMean));
	    
	    text->SetFillColor(kRed);
	  }
	  flag = AliQAv1::kERROR;
	} 	
      }
    }
  }
  
  if ((histname.EndsWith("RawsToT")) || (histname.Contains("RawsToT") && suffixTrgCl)) {
    if (!suffixTrgCl) histo->SetBit(AliQAv1::GetImageBit(), drawRawsToTSumImage);
    if (suffixTrgCl) {
      histo->SetBit(AliQAv1::GetImageBit(), kFALSE);
      for(int i=0; i<nToDrawTrgCl; i++) {
        if(histname.EndsWith(ClToDraw[i]))
          histo->SetBit(AliQAv1::GetImageBit(), kTRUE);
       }
    }
    text = (TPaveText *) histo->GetListOfFunctions()->FindObject("totMsg");
    if (histo->GetEntries()==0) {
      if (text){
	  text->Clear();
	  text->AddText("No entries. Check TOF TWiki"); 
	  text->SetFillColor(kYellow);
      }
      flag = AliQAv1::kWARNING;
    } else {
      Float_t timeMean = histo->GetMean();
      if ( (timeMean > minTOFrawTot) && (timeMean < maxTOFrawTot) ) {
	flag = AliQAv1::kINFO;
	if (text){
	  text->Clear();
	  text->AddText("Mean inside limits: OK!!!"); 
	  text->SetFillColor(kGreen);
	}
      } else {
	flag = AliQAv1::kERROR;
	AliWarning(Form("ToT mean = %5.2f ns", timeMean));
          if (text){
	    text->Clear();
	    text->AddText("Mean outside limits."); 
	text->AddText("If NOT a technical run,"); 
	text->AddText("call TOF on-call."); 
	text->SetFillColor(kRed);
          }
      }
    }
  }
  
  // if (histname.EndsWith("HitMap24")) {
  //   TPaveText * stats = (TPaveText *) histo->GetListOfFunctions()->FindObject("stats");
  //   if (stats) stats->Delete();
  // }
  return flag;
}
