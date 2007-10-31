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

/*
Produces the data needed to calculate the quality assurance. 
All data must be mergeable objects.
S.Radomski Uni-Heidelberg October 2007
*/

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1D.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliTRDdigit.h"
#include "AliTRDhit.h"
#include "AliTRDcluster.h"
#include "AliTRDQADataMaker.h"

#include "AliTRDRawStreamV2.h"

ClassImp(AliTRDQADataMaker)
           
//____________________________________________________________________________ 
  AliTRDQADataMaker::AliTRDQADataMaker() : 
  AliQADataMaker(AliQA::GetDetName(AliQA::kTRD), "TRD Quality Assurance Data Maker")
{
  // ctor
}

//____________________________________________________________________________ 
AliTRDQADataMaker::AliTRDQADataMaker(const AliTRDQADataMaker& qadm) :
  AliQADataMaker()
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliTRDQADataMaker& AliTRDQADataMaker::operator = (const AliTRDQADataMaker& qadm )
{
  // Equal operator.
  this->~AliTRDQADataMaker();
  new(this) AliTRDQADataMaker(qadm);
  return *this;
}
 
//____________________________________________________________________________ 
void AliTRDQADataMaker::EndOfDetectorCycle()
{
  //Detector specific actions at end of cycle
}

//____________________________________________________________________________ 
void AliTRDQADataMaker::InitESDs()
{
  //create ESDs histograms in ESDs subdir
  const Int_t nhist = 1;
  TH1D *hist[nhist];
 
  hist[0] = new TH1D("qaTRD_esd_bits", ";Bits", 64, -0.5, 63.5);

  for(Int_t i=0; i<nhist; i++) {
    //hist[i]->Sumw2();
    Add2ESDsList(hist[i], i);
  }
}

//____________________________________________________________________________ 
void AliTRDQADataMaker::InitHits()
{
  // create Hits histograms in Hits subdir
  const Int_t nhist = 4;
  TH1D *hist[nhist];
  
  hist[0] = new TH1D("qaTRD_hits_det", ";Detector Id of the hit", 540, -0.5, 539.5) ; 
  
  hist[1] = new TH1D("qaTRD_hist_Qdrift", ";Charge from tracks", 100, 0, 100);
  hist[2] = new TH1D("qaTRD_hist_Qamp", ";Charge from TRD photon", 100, 0, 100);
  hist[3] = new TH1D("qaTRD_hist_Qphoton", ";Charge from TRD photon", 100, 0, 100);

  for(Int_t i=0; i<nhist; i++) {
    //hist[i]->Sumw2();
    Add2HitsList(hist[i], i);
  }
}

//____________________________________________________________________________ 
void AliTRDQADataMaker::InitDigits()
{
  // create Digits histograms in Digits subdir
  
  const Int_t nhist = 3;
  TH1D *hist[nhist];
  
  hist[0] = new TH1D("qaTRD_digits_det", ";Detector Id of the digit", 540, -0.5, 539.5);
  hist[1] = new TH1D("qaTRD_digits_time", ";Time bin", 40, -0.5, 39.5);
  hist[2] = new TH1D("qaTRD_digits_amp", ";Amplitude", 100, 0, 100.);

  for(Int_t i=0; i<nhist; i++) {
    hist[i]->Sumw2();
    Add2DigitsList(hist[i], i);
  }

}

//____________________________________________________________________________ 
void AliTRDQADataMaker::InitRecPoints()
{
  // create Reconstructed Points histograms in RecPoints subdir
  const Int_t nhist = 7;
  TH1D *hist[nhist];

  hist[0] = new TH1D("qaTRD_recPoints_det", ";Detector ID of the cluster", 540, -0.5, 539.5);
  hist[1] = new TH1D("qaTRD_recPoints_amp", ";Amplitude", 200, -0.5, 199.5);
  hist[2] = new TH1D("qaTRD_recPoints_npad", ";Number of Pads", 12, -0.5, 11.5);
  hist[3] = new TH1D("qaTRD_recPoints_dist2", ";residuals [2pad]", 100, -1, 1);
  hist[4] = new TH1D("qaTRD_recPoints_dist3", ";residuals [3pad]", 100, -1, 1);
  hist[5] = new TH1D("qaTRD_recPoints_dist4", ";residuals [4pad]", 100, -1, 1);
  hist[6] = new TH1D("qaTRD_recPoints_dist5", ";residuals [5pad]", 100, -1, 1);

  for(Int_t i=0; i<nhist; i++) {
    //hist[i]->Sumw2();
    Add2RecPointsList(hist[i], i);
  }
}

//____________________________________________________________________________ 
void AliTRDQADataMaker::InitRaws()
{
  // create Raws histograms in Raws subdir
  const Int_t kSM = 18;
  const Int_t nhist = 6+kSM;
  TH1D *hist[nhist];
 
  hist[0] = new TH1D("qaTRD_raws_det", ";detector", 540, -0.5, 539.5);
  hist[1] = new TH1D("qaTRD_raws_sig", ";signal", 100, -0.5, 99.5);
  hist[2] = new TH1D("qaTRD_raws_sigCentral", "; signal central bin", 100, -0.5, 99.5);
  hist[3] = new TH1D("qaTRD_raws_sigTail", ";signal cluster", 100, -0.5, 99.5);
  hist[4] = new TH1D("qaTRD_raws_tiemBin", ";time bin", 40, -0.5, 39.5); 
  hist[5] = new TH1D("qaTRD_rows_smId", ";supermodule", 18, -0.5, 17.5);

  // one char per ADC chanell
  const Int_t nADC = 30 * 8 * 16 * 22;
  for(Int_t i=0; i<kSM; i++)
    hist[6+i] = new TH1D(Form("qaTRD_raws_sm%d",i),"",nADC, -0.5, nADC-0.5); 

  for(Int_t i=0; i<nhist; i++) {
    //hist[i]->Sumw2();
    Add2RawsList(hist[i], i);
  }
}

//____________________________________________________________________________ 
void AliTRDQADataMaker::InitSDigits()
{
  // create SDigits histograms in SDigits subdir
  
  const Int_t nhist = 2;
  TH1D *hist[nhist];
  
  hist[0] = new TH1D("qaTRD_digits_det", ";Detector Id of the digit", 540, -0.5, 539.5);
  hist[1] = new TH1D("qaTRD_digits_amp", ";Amplitude", 100, -0.5, 99.5);

  for(Int_t i=0; i<nhist; i++) {
    hist[i]->Sumw2();
    Add2SDigitsList(hist[i], i);
  }
}

//____________________________________________________________________________
void AliTRDQADataMaker::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs
  
   Int_t nTracks = esd->GetNumberOfTracks();

   for(Int_t i=0; i<nTracks; i++) {

     AliESDtrack *track = esd->GetTrack(i);
     UInt_t status = track->GetStatus();     
     UInt_t u = 1;
     for(Int_t bit=0; bit<64; bit++) 
       if (u<<bit & status) GetESDsData(0)->Fill(bit);
   }
}

//____________________________________________________________________________
void AliTRDQADataMaker::MakeHits(TClonesArray * hits)
{
  //make QA data from Hits
  //printf("making QA for TRD hits from an array %d\n", hits->GetEntriesFast());

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
void AliTRDQADataMaker::MakeHits(TTree * hitTree)
{
  //make QA data from Hits
  //printf("making QA for TRD hits from a tree\n");
  
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
}
//____________________________________________________________________________
void AliTRDQADataMaker::MakeDigits(TClonesArray * digits)
{
  // makes data from Digits
  
  TIter next(digits) ; 
  AliTRDdigit * digit ; 
  while ( (digit = dynamic_cast<AliTRDdigit *>(next())) ) {
    GetDigitsData(0)->Fill(digit->GetDetector());
    GetDigitsData(1)->Fill(digit->GetTime());
    GetDigitsData(2)->Fill(digit->GetAmp());
  }  
}

//____________________________________________________________________________
void AliTRDQADataMaker::MakeSDigits(TClonesArray * sdigits)
{
  // makes data from Digits
  
  TIter next(sdigits) ; 
  AliTRDdigit * digit ; 
  while ( (digit = dynamic_cast<AliTRDdigit *>(next())) ) {
    GetDigitsData(0)->Fill(digit->GetDetector());
    GetDigitsData(1)->Fill(digit->GetAmp());
  }  
}

//____________________________________________________________________________
void AliTRDQADataMaker::MakeRaws(AliRawReader* rawReader)
{
  // 157
  // T9 -- T10

  //const Int_t kSM = 18;
  //const Int_t kROC = 30;
  const Int_t kROB = 8;
  //const Int_t kLayer = 6;
  //const Int_t kStack = 5;
  const Int_t kMCM = 16;
  const Int_t kADC = 22;

  AliTRDRawStreamV2 *raw = new AliTRDRawStreamV2(rawReader);
  
  raw->SetRawVersion(3);
  raw->Init();
  
  while (raw->Next()) {
    
    GetRawsData(0)->Fill(raw->GetDet());
    
    Int_t *sig = raw->GetSignals();
    for(Int_t i=0; i<3; i++) GetRawsData(1)->Fill(sig[i]);
    
    GetRawsData(2)->Fill(sig[1]);
    GetRawsData(3)->Fill(sig[0]);
    GetRawsData(3)->Fill(sig[2]);
    
    GetRawsData(4)->Fill(raw->GetTimeBin());
    
    // calculate the index;
    Int_t sm = raw->GetSM();
    Int_t roc = raw->GetROC();
    Int_t rob = raw->GetROB();
    Int_t mcm = raw->GetMCM();
    Int_t adc = raw->GetADC();
    
    Int_t index = roc * (kROB*kMCM*kADC) + rob * (kMCM*kADC) + mcm * kADC + adc;
    GetRawsData(5)->Fill(sm);
    //printf("SM = %d\n", sm);
    GetRawsData(6+sm)->Fill(index);
  }
}

//____________________________________________________________________________
void AliTRDQADataMaker::MakeRecPoints(TTree * clustersTree)
{
  
  // makes data from RecPoints
  // 
  
  Int_t nsize = Int_t(clustersTree->GetTotBytes() / (sizeof(AliTRDcluster))); 
  TObjArray *clusterArray = new TObjArray(nsize+1000); 
  
  TBranch *branch = clustersTree->GetBranch("TRDcluster");
  if (!branch) {
    AliError("Can't get the branch !");
    return;
  }
  branch->SetAddress(&clusterArray); 
  
  // Loop through all entries in the tree
  Int_t nEntries   = (Int_t) clustersTree->GetEntries();
  Int_t nbytes     = 0;
  AliTRDcluster *c = 0;
  
  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {    
    
    // Import the tree
    nbytes += clustersTree->GetEvent(iEntry);  
    
    // Get the number of points in the detector
    Int_t nCluster = clusterArray->GetEntriesFast();  
    
    // Loop through all TRD digits
    for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) { 
      c = (AliTRDcluster *) clusterArray->UncheckedAt(iCluster);
      
      GetRecPointsData(0)->Fill(c->GetDetector());
      GetRecPointsData(1)->Fill(c->GetQ());
      GetRecPointsData(2)->Fill(c->GetNPads());
      if (c->GetNPads() < 6)
	GetRecPointsData(1+c->GetNPads())->Fill(c->GetCenter());
    }
  }
  
  delete clusterArray;
}

//____________________________________________________________________________ 
void AliTRDQADataMaker::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle

}
//__________________________________________________________________________
Int_t AliTRDQADataMaker::CheckPointer(TObject *obj, const char *name) {

  if (!obj) AliWarning(Form("null pointer: %s", name));
  return !!obj;
}
//__________________________________________________________________________
