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

/* $Id$ */

// *****************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  W.Ferrarese  P.Cerello  Mag 2008
//  INFN Torino

// --- ROOT system ---
#include "TH1.h"
#include <Riostream.h>

// --- AliRoot header files ---
#include "AliITSQAChecker.h"
#include "AliITSQASPDChecker.h"
#include "AliITSQASDDChecker.h"
#include "AliITSQASSDChecker.h"

ClassImp(AliITSQAChecker)

//____________________________________________________________________________
AliITSQAChecker::AliITSQAChecker(Bool_t kMode, Short_t subDet, Short_t ldc) :
AliQACheckerBase("ITS","SDD Quality Assurance Checker"),
fkOnline(0),
fDet(0),  
fLDC(0),
fSPDOffset(0), 
fSDDOffset(0), 
fSSDOffset(0),
fSPDChecker(0),  // SPD Checker
fSDDChecker(0),  // SDD Checker
fSSDChecker(0)  // SSD Checker
{
  // Standard constructor
  fkOnline = kMode; fDet = subDet; fLDC = ldc;
  if(fDet == 0 || fDet == 1) {
    AliDebug(1,"AliITSQAChecker::Create SPD Checker\n");
  }
  if(fDet == 0 || fDet == 2) {
    AliDebug(1,"AliITSQAChecker::Create SDD Checker\n");
  }
  if(fDet == 0 || fDet == 3) {
    AliDebug(1,"AliITSQAChecker::Create SSD Checker\n");
  }

}

//____________________________________________________________________________
AliITSQAChecker::AliITSQAChecker(const AliITSQAChecker& qac):
AliQACheckerBase(qac.GetName(), qac.GetTitle()), 
fkOnline(qac.fkOnline), 
fDet(qac.fDet), 
fLDC(qac.fLDC), 
fSPDOffset(qac.fSPDOffset), 
fSDDOffset(qac.fSDDOffset), 
fSSDOffset(qac.fSSDOffset), 
fSPDChecker(0), 
fSDDChecker(0), 
fSSDChecker(0) {
  // copy constructor
  AliError("Copy should not be used with this class\n");
}
//____________________________________________________________________________
AliITSQAChecker& AliITSQAChecker::operator=(const AliITSQAChecker& qac){
  // assignment operator
  this->~AliITSQAChecker();
  new(this)AliITSQAChecker(qac);
  return *this;
}

//____________________________________________________________________________
Double_t AliITSQAChecker::Check(AliQA::ALITASK_t index, TObjArray * list)
{
  
  // Super-basic check on the QA histograms on the input list:
  // look whether they are empty!
  if(index == AliQA::kESD){
    AliDebug(1,"Checker for ESD");
    Int_t tested = 0;
    Int_t empty = 0;
    Double_t rv = 0.;
    // The following flags are set to kTRUE if the corresponding
    // QA histograms exceed a given quality threshold
    Bool_t cluMapSA = kFALSE;
    Bool_t cluMapMI = kFALSE;
    Bool_t cluMI = kFALSE;
    Bool_t cluSA = kFALSE;
    Bool_t verSPDZ = kFALSE;
    if (list->GetEntries() == 0) {
      rv = 0.; // nothing to check
    }
    else {
      TIter next1(list);
      TH1 * hdata;
      Int_t nskipped=0;
      Bool_t skipped[6]={kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      // look for layers that we wanted to skip
      while ( (hdata = dynamic_cast<TH1 *>(next1())) ) {
	if(!hdata) continue;
	TString hname = hdata->GetName();
	if(!hname.Contains("hESDSkippedLayers")) continue;
	for(Int_t k=1; k<7; k++) {
	  if(hdata->GetBinContent(k)>0) { 
	    nskipped++; 
	    skipped[k-1]=kTRUE; 
	  } 
	} 
      }
      TIter next(list);
      while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
	if(hdata){
	  TString hname = hdata->GetName();
	  Double_t entries = hdata->GetEntries();
	  ++tested;
	  if(!(entries>0.))++empty;
	  AliDebug(1,Form("ESD hist name %s - entries %12.1g",hname.Data(),entries));

	  if(hname.Contains("hESDClusterMapSA") && entries>0.){
	    cluMapSA = kTRUE;
	    AliDebug(1,Form("Processing histogram %s",hname.Data()));
	    // Check if there are layers with anomalously low 
	    // contributing points to SA reconstructed tracks
	    for(Int_t k=1;k<7;k++){
	      // check if the layer was skipped
	      if(skipped[k-1]) continue;
	      if(hdata->GetBinContent(k)<0.5*(entries/6.)){
		cluMapSA = kFALSE;
		AliInfo(Form("SA tracks have few points on layer %d - look at histogram hESDClustersSA",k));
	      }
	    }  
	  }

	  else if(hname.Contains("hESDClusterMapMI") && entries>0.){
	    // Check if there are layers with anomalously low 
	    // contributing points to MI reconstructed tracks
	    AliDebug(1,Form("Processing histogram %s",hname.Data()));
	    cluMapMI = kTRUE;
	    for(Int_t k=1;k<7;k++){
	      // check if the layer was skipped
	      if(skipped[k-1]) continue;
	      if(hdata->GetBinContent(k)<0.5*(entries/6.)){
		cluMapMI = kFALSE;
		AliInfo(Form("MI tracks have few points on layer %d - look at histogram hESDClustersMI",k));
	      }
	    }  
	  }

	  else if(hname.Contains("hESDClustersMI") && entries>0.){
	    // Check if 6 clusters MI tracks are the majority
	    AliDebug(1,Form("Processing histogram %s",hname.Data()));
	    cluMI = kTRUE;
	    Double_t maxlaytracks = hdata->GetBinContent(7-nskipped);
	    for(Int_t k=2; k<7-nskipped; k++){
	      if(hdata->GetBinContent(k)>maxlaytracks){
		cluMI = kFALSE;
		AliInfo(Form("MI Tracks with %d clusters are more than tracks with %d clusters. Look at histogram hESDClustersMI",k-1,6-nskipped));
	      }
	    }
	  }

	  else if(hname.Contains("hESDClustersSA") && entries>0.){
	    // Check if 6 clusters SA tracks are the majority
	    AliDebug(1,Form("Processing histogram %s",hname.Data()));
	    cluSA = kTRUE;
	    Double_t maxlaytracks = hdata->GetBinContent(7-nskipped);
	    for(Int_t k=2; k<7-nskipped; k++){
	      if(hdata->GetBinContent(k)>maxlaytracks){
		cluSA = kFALSE;
		AliInfo(Form("SA Tracks with %d clusters are more than tracks with %d clusters. Look at histogram hESDClustersSA",k-1,6-nskipped));
	      }
	    }
	  }

	  else if(hname.Contains("hSPDVertexZ") && entries>0.){
	    // Check if average Z vertex coordinate is -5 < z < 5 cm
	    AliDebug(1,Form("Processing histogram %s",hname.Data()));
	    verSPDZ = kTRUE;
	    if(hdata->GetMean()<-5. && hdata->GetMean()>5.){
	      verSPDZ = kFALSE;
	      AliInfo(Form("Average z vertex coordinate is at z= %10.4g cm",hdata->GetMean()));
	    }
	  }

	}
	else{
	  AliError("ESD Checker - invalid data type");
	}
	  
	rv = 0.;
	if(tested>0){
	  if(tested == empty){
	    rv = 0.1;
	    AliWarning("All ESD histograms are empty");
	  }
	  else {
	    rv = 0.1+0.4*(static_cast<Double_t>(tested-empty)/static_cast<Double_t>(tested));
	    if(cluMapSA)rv+=0.1;
	    if(cluMapMI)rv+=0.1;
	    if(cluMI)rv+=0.1;
	    if(cluSA)rv+=0.1;
	    if(verSPDZ)rv+=0.1;
	  }
	}
  
      }
    }  
    AliInfo(Form("ESD - Tested %d histograms, Return value %f \n",tested,rv));
    return rv;
  }  // end of ESD QA

  //____________________________________________________________________________

  Double_t spdCheck, sddCheck, ssdCheck;
  Double_t retval = 1.;
  //pixel
  if(fDet == 0 || fDet == 1) {
    AliDebug(1,"AliITSQAChecker::Create SPD Checker\n");
    if(!fSPDChecker) {
      fSPDChecker = new AliITSQASPDChecker();
    }
    fSPDChecker->SetTaskOffset(fSPDOffset);
    spdCheck = fSPDChecker->Check(index, list);
    if(spdCheck<retval)retval = spdCheck;
  }
  //drift
  if(fDet == 0 || fDet == 2) {
    AliDebug(1,"AliITSQAChecker::Create SDD Checker\n");
    if(!fSDDChecker) {
      fSDDChecker = new AliITSQASDDChecker();
    }
    fSDDChecker->SetTaskOffset(fSDDOffset);
    sddCheck = fSDDChecker->Check(index, list);
    if(sddCheck<retval)retval = sddCheck;
  }
  //strip
  if(fDet == 0 || fDet == 3) {
    AliDebug(1,"AliITSQAChecker::Create SSD Checker\n");
    if(!fSSDChecker) {
      fSSDChecker = new AliITSQASSDChecker();
      AliInfo(Form("Number of monitored objects SSD: %d", list->GetEntries()));
    }
    fSSDChecker->SetTaskOffset(fSSDOffset);
    ssdCheck = fSSDChecker->Check(index, list);
    if(ssdCheck<retval)retval = ssdCheck;  }
  // here merging part for common ITS QA result
  // 
  AliDebug(1,Form("AliITSQAChecker::QAChecker returned value is %f \n",retval));
  return retval;
  
}


//____________________________________________________________________________
void AliITSQAChecker::SetTaskOffset(Int_t SPDOffset, Int_t SDDOffset, Int_t SSDOffset)
{
  //Setting the 3 offsets for each task called
  fSPDOffset = SPDOffset;
  fSDDOffset = SDDOffset;
  fSSDOffset = SSDOffset;
}

 //____________________________________________________________________________
 void AliITSQAChecker::SetDetTaskOffset(Int_t subdet,Int_t offset)
 {
   switch(subdet){
   case 1:
     SetSPDTaskOffset(offset);
     break;
   case 2:
     SetSDDTaskOffset(offset);
     break;
   case 3:
     SetSSDTaskOffset(offset);
     break;
   default:
     AliWarning("No specific (SPD,SDD or SSD) subdetector correspond to to this number!!! all offsets set to zero for all the detectors\n");
     SetTaskOffset(0, 0, 0);
     break;
   }
 }
