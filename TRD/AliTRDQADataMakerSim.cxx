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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Produces the data needed to calculate the quality assurance.          //
//  All data must be mergeable objects.                                   //
//                                                                        //
//  Author:                                                               //
//    Sylwester Radomski (radomski@physi.uni-heidelberg.de)               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1D.h> 
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TTree.h>

// --- AliRoot header files ---
//#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliTRDdigit.h"
#include "AliTRDhit.h"
//#include "AliTRDcluster.h"
#include "AliTRDQADataMakerSim.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDgeometry.h"
#include "AliTRDarrayADC.h"
#include "AliTRDarraySignal.h"
//#include "AliTRDrawStream.h"

#include "AliQAChecker.h"

ClassImp(AliTRDQADataMakerSim)

//____________________________________________________________________________ 
  AliTRDQADataMakerSim::AliTRDQADataMakerSim() : 
  AliQADataMakerSim(AliQAv1::GetDetName(AliQAv1::kTRD), "TRD Quality Assurance Data Maker")
{
  //
  // Default constructor
}

//____________________________________________________________________________ 
AliTRDQADataMakerSim::AliTRDQADataMakerSim(const AliTRDQADataMakerSim& qadm) :
  AliQADataMakerSim()
{
  //
  // Copy constructor 
  //

  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 

}

//__________________________________________________________________
AliTRDQADataMakerSim& AliTRDQADataMakerSim::operator=(const AliTRDQADataMakerSim& qadm)
{
  //
  // Equal operator.
  //

  this->~AliTRDQADataMakerSim();
  new(this) AliTRDQADataMakerSim(qadm);
  return *this;

}

//____________________________________________________________________________ 
void AliTRDQADataMakerSim::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //
  // Detector specific actions at end of cycle
  //

  //AliInfo(Form("EndOfCycle", "Fitting RecPoints %d", task));

  // call the checker
  AliQAChecker::Instance()->Run(AliQAv1::kTRD, task, list) ;    


}

//____________________________________________________________________________ 
void AliTRDQADataMakerSim::InitHits()
{
  //
  // Create Hits histograms in Hits subdir
  //

  const Int_t kNhist = 4;
  TH1D *hist[kNhist];

  hist[0] = new TH1D("qaTRD_hits_det", ";Detector Id of the hit", 540, -0.5, 539.5) ; 

  hist[1] = new TH1D("qaTRD_hist_Qdrift", ";Charge from tracks", 100, 0, 100);
  hist[2] = new TH1D("qaTRD_hist_Qamp", ";Charge from TRD photon", 100, 0, 100);
  hist[3] = new TH1D("qaTRD_hist_Qphoton", ";Charge from TRD photon", 100, 0, 100);

  for(Int_t i=0; i<kNhist; i++) {
    //hist[i]->Sumw2();
    Add2HitsList(hist[i], i);
  }

}

//____________________________________________________________________________ 
void AliTRDQADataMakerSim::InitDigits()
{
  //
  // Create Digits histograms in Digits subdir
  //

  const Int_t kNhist = 3;
  TH1D *hist[kNhist];

  hist[0] = new TH1D("qaTRD_digits_det", ";Detector Id of the digit", 540, -0.5, 539.5);
  hist[1] = new TH1D("qaTRD_digits_time", ";Time bin", 40, -0.5, 39.5);
  hist[2] = new TH1D("qaTRD_digits_amp", ";Amplitude", 100, -5.5, 94.5);

  for(Int_t i=0; i<kNhist; i++) {
    hist[i]->Sumw2();
    Add2DigitsList(hist[i], i);
  }

}

//____________________________________________________________________________ 
void AliTRDQADataMakerSim::InitSDigits()
{
  //
  // Create SDigits histograms in SDigits subdir
  //

  const Int_t kNhist = 3;
  TH1D *hist[kNhist];

  hist[0] = new TH1D("qaTRD_sdigits_det", ";Detector Id of the digit", 540, -0.5, 539.5);
  hist[1] = new TH1D("qaTRD_sdigits_time", ";Time bin", 40, -0.5, 39.5);
  hist[2] = new TH1D("qaTRD_sdigits_amp", ";Amplitude", 100, 0, 1e7);

  for(Int_t i=0; i<kNhist; i++) {
    hist[i]->Sumw2();
    Add2SDigitsList(hist[i], i);
  }

}

//____________________________________________________________________________
void AliTRDQADataMakerSim::MakeHits(TClonesArray * hits)
{
  //
  // Make QA data from Hits
  //

  TIter next(hits); 
  AliTRDhit * hit; 

  while ( (hit = dynamic_cast<AliTRDhit *>(next())) ) {
    GetHitsData(0)->Fill(hit->GetDetector());
    Double_t q = TMath::Abs(hit->GetCharge());

    if (hit->FromDrift()) GetHitsData(1)->Fill(q);
    if (hit->FromAmplification()) GetHitsData(2)->Fill(q);
    if (hit->FromTRphoton()) GetHitsData(3)->Fill(q);
  }

}

//____________________________________________________________________________
void AliTRDQADataMakerSim::MakeHits(TTree * hitTree)
{
  //
  // Make QA data from Hits
  //

  if (!CheckPointer(hitTree, "TRD hits tree")) return;

  TBranch *branch = hitTree->GetBranch("TRD");
  if (!CheckPointer(branch, "TRD hits branch")) return;

  Int_t nhits = (Int_t)(hitTree->GetTotBytes()/sizeof(AliTRDhit));
  TClonesArray *hits = new TClonesArray("AliTRDhit", nhits+1000);
  TClonesArray *tmp = new TClonesArray("AliTRDhit", 1000);
  branch->SetAddress(&tmp);

  Int_t index = 0;
  Int_t nEntries = (Int_t)branch->GetEntries();
  for(Int_t i = 0; i < nEntries; i++) {
    branch->GetEntry(i);
    Int_t nHits = (Int_t)tmp->GetEntries();
    for(Int_t j=0; j<nHits; j++) {
      AliTRDhit *hit = (AliTRDhit*)tmp->At(j);
      new((*hits)[index++]) AliTRDhit(*hit);
    }
  }

  tmp->Delete();
  delete tmp;
  MakeHits(hits);
  hits->Delete();
  delete hits;

}

//____________________________________________________________________________
void AliTRDQADataMakerSim::MakeDigits(TClonesArray * digits)
{
  //
  // Makes data from Digits
  //

  TIter next(digits) ; 
  AliTRDdigit * digit ; 
  
  // Info("Make digits", "From the arrya");
  
  while ( (digit = dynamic_cast<AliTRDdigit *>(next())) ) {
    if (digit->GetAmp() < 1) continue;
    GetDigitsData(0)->Fill(digit->GetDetector());
    GetDigitsData(1)->Fill(digit->GetTime());
    GetDigitsData(2)->Fill(digit->GetAmp());
  }

}

//____________________________________________________________________________
void AliTRDQADataMakerSim::MakeDigits(TTree * digits)
{
  //
  // Makes data from digits tree
  //

  // Info("Make digits", "From a tree");

  AliTRDdigitsManager *digitsManager = new AliTRDdigitsManager();
  digitsManager->CreateArrays();
  digitsManager->ReadDigits(digits);

  TH1D *histDet = (TH1D*)GetDigitsData(0);
  TH1D *histTime = (TH1D*)GetDigitsData(1);
  TH1D *histSignal = (TH1D*)GetDigitsData(2);

  for (Int_t i = 0; i < AliTRDgeometry::kNdet; i++) 
    {
      AliTRDarrayADC *digitsIn = (AliTRDarrayADC *) digitsManager->GetDigits(i);      

      // This is to take care of switched off super modules
      if (digitsIn->GetNtime() == 0) continue;
      
      digitsIn->Expand();
      
      //AliTRDSignalIndex* indexes = digitsManager->GetIndexes(i);
      //if (indexes->IsAllocated() == kFALSE) digitsManager->BuildIndexes(i);
      
      Int_t nRows = digitsIn->GetNrow();
      Int_t nCols = digitsIn->GetNcol();
      Int_t nTbins = digitsIn->GetNtime();
      
      for(Int_t row = 0; row < nRows; row++) 
	for(Int_t col = 0; col < nCols; col++) 
	  for(Int_t time = 0; time < nTbins; time++) 
	    {   
	      Float_t signal = digitsIn->GetData(row,col,time);
	      if (signal < 1) continue;
	      histDet->Fill(i);
	      histTime->Fill(time);
	      histSignal->Fill(signal);
	    }

    //delete digitsIn;
  }
  delete digitsManager;
}

//____________________________________________________________________________
void AliTRDQADataMakerSim::MakeSDigits(TClonesArray * sdigits)
{
  //
  // Makes data from Digits
  //

  TIter next(sdigits) ; 
  AliTRDdigit * digit ; 
  while ( (digit = dynamic_cast<AliTRDdigit *>(next())) ) {
    GetDigitsData(0)->Fill(digit->GetDetector());
    GetDigitsData(1)->Fill(digit->GetTime());
    GetDigitsData(2)->Fill(digit->GetAmp());
  }

}

//____________________________________________________________________________
void AliTRDQADataMakerSim::MakeSDigits(TTree * digits)
{
  //
  // Makes data from SDigits
  //

  AliTRDdigitsManager *digitsManager = new AliTRDdigitsManager();
  digitsManager->SetSDigits(kTRUE);
  digitsManager->CreateArrays();
  digitsManager->ReadDigits(digits);

  TH1D *histDet = (TH1D*)GetSDigitsData(0);
  TH1D *histTime = (TH1D*)GetSDigitsData(1);
  TH1D *histSignal = (TH1D*)GetSDigitsData(2);

  for (Int_t i = 0; i < AliTRDgeometry::kNdet; i++) 
    {
      AliTRDarraySignal *digitsIn = (AliTRDarraySignal *) digitsManager->GetSDigits(i);
  
    // This is to take care of switched off super modules
      if (digitsIn->GetNtime() == 0) continue;
      
    digitsIn->Expand(); 
      
      //AliTRDSignalIndex* indexes = digitsManager->GetIndexes(i);
      //if (indexes->IsAllocated() == kFALSE) digitsManager->BuildIndexes(i);
      Int_t nRows = digitsIn->GetNrow();
      Int_t nCols = digitsIn->GetNcol();
      Int_t nTbins = digitsIn->GetNtime();

      for(Int_t row = 0; row < nRows; row++) 
	{
	  for(Int_t col = 0; col < nCols; col++) 
	    {
	      for(Int_t time = 0; time < nTbins; time++) 
		{	    
		  Float_t signal = digitsIn->GetData(row,col,time);
		  if (signal < 1) continue;
		  histDet->Fill(i);
		  histTime->Fill(time);
		  histSignal->Fill(signal);
		}
	    }
	}
      // delete digitsIn;
    }
  delete digitsManager;
}

//____________________________________________________________________________ 
void AliTRDQADataMakerSim::StartOfDetectorCycle()
{
  //
  // Detector specific actions at start of cycle
  //

}

//__________________________________________________________________________
Int_t AliTRDQADataMakerSim::CheckPointer(TObject *obj, const char *name) 
{
  //
  // Checks initialization of pointers
  //

  if (!obj) AliWarning(Form("null pointer: %s", name));
  return !!obj;

}
