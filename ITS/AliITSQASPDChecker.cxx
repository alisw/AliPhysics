/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id $ */

// *****************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  P. Cerello Apr 2008
//  INFN Torino

// --- ROOT system ---
#include "TH1.h"
#include "TString.h"

// --- AliRoot header files ---
#include "AliITSQASPDChecker.h"
#include "AliITSQADataMakerRec.h"
#include "AliLog.h"

ClassImp(AliITSQASPDChecker)
//__________________________________________________________________
AliITSQASPDChecker& AliITSQASPDChecker::operator = (const AliITSQASPDChecker& qac ) 
{
  // Equal operator.
  this->~AliITSQASPDChecker();
  new(this) AliITSQASPDChecker(qac);
  return *this;
}


//__________________________________________________________________
Double_t AliITSQASPDChecker::Check(AliQAv1::ALITASK_t index, TObjArray * list, const AliDetectorRecoParam * /*recoParam*/)
{
//
// General methods for SPD Cheks to be used in RAWS and REC ALITASK_t
//

  AliDebug(2, Form("AliITSQASPDChecker called with offset: %d\n", fSubDetOffset));

  Double_t test = 0.0;
  Int_t count = 0;
  // Checks for ALITASK_t AliQAv1::kRAW
  if(index == AliQAv1::kRAW) {
  return CheckRawData(list);
  } else {
  if (list->GetEntries() == 0) {
    test = 1.; // nothing to check
  }
  else {
    TIter next(list);
    TH1 * hdata;
    count = 0;
    while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
      if (hdata) {
        TString histName = hdata->GetName();
        if (!histName.Contains("_SPD")) continue;
        Double_t rv = 0.;
        if (hdata->GetEntries()>0) rv = 1;
        if (histName.Contains("LayPattern")) {
         if (hdata->GetBinContent(1)) {
           Double_t ratio=hdata->GetBinContent(2)/hdata->GetBinContent(1);
           AliDebug(2, Form("%s: ratio RecPoints lay2 / lay1 = %f", hdata->GetName(), ratio));
         }
         else
           AliDebug(AliQAv1::GetQADebugLevel(), "No RecPoints in lay1");
       }
        else if(histName.Contains("ModPattern")) {
           Int_t ndead=0;
           for(Int_t ibin=0;ibin<hdata->GetNbinsX();ibin++) {
             if(histName.Contains("SPD1") && ibin<80 && hdata->GetBinContent(ibin+1)>0) ndead++;
             if(histName.Contains("SPD2") && ibin>79 && hdata->GetBinContent(ibin+1)>0) ndead++;
           }
           AliDebug(2, Form("%s: Entries = %d  number of empty modules = %d", 
                        hdata->GetName(),(Int_t)hdata->GetEntries(),ndead));
        }
        else if(histName.Contains("SizeYvsZ")) {
           Double_t meanz=hdata->GetMean(1);
           Double_t meany=hdata->GetMean(2);
           Double_t rmsz=hdata->GetRMS(1);
           Double_t rmsy=hdata->GetRMS(2);
           AliDebug(AliQAv1::GetQADebugLevel(), Form("%s: Cluster sizeY mean = %f  rms = %f", hdata->GetName(),meany,rmsy));
           AliDebug(AliQAv1::GetQADebugLevel(), Form("%s: Cluster sizeZ mean = %f  rms = %f", hdata->GetName(),meanz,rmsz));
        }
        else if(histName.Contains("SPDMultiplicity")) {
           AliDebug(2, Form("%s: Events = %d  mean = %f  rms = %f",
                        hdata->GetName(),(Int_t)hdata->GetEntries(),hdata->GetMean(),hdata->GetRMS()));}

       // else AliDebug(AliQAv1::GetQADebugLevel(), Form("%s -> %f", hdata->GetName(), rv));
        count++;
        test += rv;
      }
      else {
        AliError("Data type cannot be processed") ;
      }
    }

    if (count != 0) {
      if (AliITSQADataMakerRec::AreEqual(test,0)) {
        AliWarning("Histograms are there, but they are all empty: setting flag to kWARNING");
        test = fHighSPDValue[AliQAv1::kWARNING];  //upper limit value to set kWARNING flag for a task
      }
      else {
        test /= count;
      }
    }
  }
}
  AliDebug(AliQAv1::GetQADebugLevel(), Form("Test Result = %f", test));
  return test ;

}
//__________________________________________________________________
Double_t AliITSQASPDChecker::CheckRawData(const TObjArray * list) {
//
// Checks on the raw data histograms [ preliminary version ]
// The output of this method is the fraction of SPD histograms which are processed by the checker. 
// The methods returns fHighSPDValue[AliQAv1::kFATAL] in case of data format errors or MEB errors
// 
// A. Mastroserio

Double_t test =0;

// basic checks on input data
if(!list) {
 AliError("NO histogram list for RAWS");
 return test;
 }

 if(list->GetEntries() == 0) {
 AliWarning("No histograms in RAW list \n");
 return test;
 }

// loop over the raw data histograms
TIter next(list);
TH1 * hdata;
Double_t totalHistos = 0;
Double_t goodHistos = 0; // number of histograms which passed the checks
Double_t response =0;
Bool_t fatalProblem = kFALSE;

while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
if (hdata) {
        TString histName = hdata->GetName();
        if(!histName.Contains("SPD")) continue;
        totalHistos++;
        // data format error
        if(histName.Contains("SPDErrorsAll")){
        if(hdata->Integral(0,hdata->GetNbinsX())>0){
           response = fHighSPDValue[AliQAv1::kFATAL];
	   fatalProblem=kTRUE;
           break;
         }
        }
        // MEB error
        if(histName.Contains("MEB")){
          if(hdata->GetEntries()>0){
	   AliWarning("************* MEB ERROR!!!! ****************\n");
           response = fHighSPDValue[AliQAv1::kFATAL];
	   fatalProblem=kTRUE;
           break;
          }
         }
        goodHistos++;
        }
}
    if(!fatalProblem) response = goodHistos/totalHistos;
   // printf("n histos %f - good ones %f ----> ratio %f , fatal response %i\n",totalHistos,goodHistos,goodHistos/totalHistos,(Int_t)fatalProblem);
    return response;
}

//__________________________________________________________________
void AliITSQASPDChecker::SetTaskOffset(Int_t TaskOffset)
{
// Offset for SPD within ITS QA
  fSubDetOffset = TaskOffset;
}

//__________________________________________________________________
void AliITSQASPDChecker::SetStepBit(const Double_t *steprange) 
{
// Step bit for SPD within ITS QA
  fStepBitSPD = new Double_t[AliQAv1::kNBIT];
  for(Int_t bit=0;bit<AliQAv1::kNBIT;bit++)
    {
      fStepBitSPD[bit]=steprange[bit];
    }
}


//__________________________________________________________________
void  AliITSQASPDChecker::SetSPDLimits(const Float_t *lowvalue, const Float_t * highvalue)
{
// SPD limints for QA bit within general ITS QA
  fLowSPDValue = new Float_t[AliQAv1::kNBIT];
  fHighSPDValue= new Float_t[AliQAv1::kNBIT];

  for(Int_t bit=0;bit<AliQAv1::kNBIT;bit++)
    {
      fLowSPDValue[bit]=lowvalue[bit];
      fHighSPDValue[bit]= highvalue[bit];
    }

}

