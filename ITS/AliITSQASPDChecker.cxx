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
Double_t AliITSQASPDChecker::Check(AliQAv1::ALITASK_t /*index*/, TObjArray * list)
{
  AliDebug(2, Form("AliITSQASPDChecker called with offset: %d\n", fSubDetOffset));

  Double_t test = 0.0;
  Int_t count = 0;

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
      if (test==0) {
        AliWarning("Histograms are there, but they are all empty: setting flag to kWARNING");
        test = 0.5;  //upper limit value to set kWARNING flag for a task
      }
      else {
        test /= count;
      }
    }
  }

  AliDebug(AliQAv1::GetQADebugLevel(), Form("Test Result = %f", test));
  return test ;

}

//__________________________________________________________________
void AliITSQASPDChecker::SetTaskOffset(Int_t TaskOffset)
{
  fSubDetOffset = TaskOffset;
}
