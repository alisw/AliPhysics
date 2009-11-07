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
//
//
// 
//Task for analysis if TPC base QA 

// QA histogram
// Dimensions
// 0 - chi2
// 1 - number of clusters
// 2 - number of findable clusters
// 3 - number of clusters/ findable clusters  
// 4 - pt          - at the entrance of the TPC
// 5 - eta         - at the entrance of the TPC
// 6 - phi         - at the entrance of the TPC

//-----------------------------------------------------------------------
// Author : M.Ivanov  marian.ivanov@cern.ch - 
//-----------------------------------------------------------------------



//
// ROOT includes
#include <TChain.h>
#include <TMath.h>
#include <TVectorD.h>
#include <TSystem.h>
#include <TFile.h>
#include <TParticle.h>

// ALIROOT includes
#include <AliAnalysisManager.h>
#include <AliESDInputHandler.h>
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMathBase.h"

#include <AliESD.h>
#include "AliExternalTrackParam.h"
#include "AliTracker.h"
#include "AliTPCseed.h"
//
#include "AliTPCtaskQA.h"
//
#include <THnSparse.h>

//

// STL includes
#include <iostream>

using namespace std;


ClassImp(AliTPCtaskQA)

//________________________________________________________________________
AliTPCtaskQA::AliTPCtaskQA() : 
  AliAnalysisTask(), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  fList(0),
  fTPCqa(0)
{
  //
  // Default constructor (should not be used)
  //
}

AliTPCtaskQA::AliTPCtaskQA(const AliTPCtaskQA& info) : 
  AliAnalysisTask(info), 
  fMCinfo(info.fMCinfo),     //! MC event handler
  fESD(info.fESD),        //!
  fList(0),
  fTPCqa(0)
{
  //
  // Dummy Copy  constructor - no copy constructor for THnSparse 
  //
  fList = (TObjArray*)(info.fList->Clone());
}



//________________________________________________________________________
AliTPCtaskQA::AliTPCtaskQA(const char *name) : 
  AliAnalysisTask(name, "AliTPCtaskQA"), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  fList(0),
  fTPCqa(0)
{
  //
  // Normal constructor
  //
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList
  DefineOutput(0, TObjArray::Class());
  //
  //make histos
  Init(); 
}

void AliTPCtaskQA::Init(){
  //
  // Init qa histogram
  // Dimensions
  //
  // 0 - chi2
  // 1 - number of clusters
  // 2 - number of findable clusters
  // 3 - number of clusters/ findable clusters  
  // 4 - pt          - at the entrance of the TPC
  // 5 - eta         - at the entrance of the TPC
  // 6 - phi         - at the entrance of the TPC
  


  Double_t xmin[7],  xmax[7];
  Int_t    nbins[7];
  // 
  nbins[0]=100;                // chi2
  xmin[0]=0; xmax[0]=10;
  //
  nbins[1]=80;                 // ncls
  xmin[1]=0; xmax[1]=160;
  //
  nbins[2]=80;                 // nclsf
  xmin[2]=0; xmax[2]=160;
  //
  nbins[3]=40;                 // ncls/nclsf
  xmin[3] =-0.1; xmax[3]=1.1;
  //
  nbins[4]=50;                 // pt
  xmin[4] =0.1; xmax[4]=100;

  nbins[5]=40;                 // eta
  xmin[5] =-2; xmax[5]=2;

  nbins[6]= 360;                 // phi - 10 bins per sector
  xmin[6] = -TMath::Pi(); xmax[6]=TMath::Pi();



  fTPCqa = new THnSparseS("TPC qa","TPC qa",7,nbins,xmin,xmax);
  //
  //
  BinLogX(fTPCqa->GetAxis(4));
  
  char *hisAxisName[7] ={"chi2/N_{cl}","N_{cl}","N_{clF}","N_{clR}","p_{t}","#eta","#phi"};
  //  
  for (Int_t i=0;i<7;i++) {
    fTPCqa->GetAxis(i)->SetTitle(hisAxisName[i]);
    fTPCqa->GetAxis(i)->SetName(hisAxisName[i]);
  }
  fList = new TObjArray(3);
  fList->AddAt(fTPCqa,0);
}




AliTPCtaskQA::~AliTPCtaskQA(){
  //
  //
  //
  delete fTPCqa;
}


//________________________________________________________________________
void AliTPCtaskQA::ConnectInputData(Option_t *) 
{
  //
  // Connect the input data
  //

  TTree* tree=dynamic_cast<TTree*>(GetInputData(0));
  if (!tree) {
    //Printf("ERROR: Could not read chain from input slot 0");
  }
  else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!esdH) {
      //Printf("ERROR: Could not get ESDInputHandler");
    }
    else {
      fESD = esdH->GetEvent();
      //Printf("*** CONNECTED NEW EVENT ****");
    }  
  }
  AliMCEventHandler* mcinfo = (AliMCEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());  
  mcinfo->SetReadTR(kTRUE);
  
  fMCinfo = mcinfo->MCEvent();


}






//________________________________________________________________________
void AliTPCtaskQA::CreateOutputObjects() 
{
  //
  // Connect the output objects
  //

}


//________________________________________________________________________
void AliTPCtaskQA::Exec(Option_t *) {
  //
  // Execute analysis for current event 
  //


  // If MC has been connected   

  if (!fMCinfo){
    cout << "Not MC info\n" << endl;
  }else{
    ProcessMCInfo();
  }
  //
  PostData(0, fList);
}      










void  AliTPCtaskQA::ProcessMCInfo(){
  //
  //
  //
  //

  if (!fTPCqa) Init();
  Int_t npart   = fMCinfo->GetNumberOfTracks();
  Int_t ntracks = fESD->GetNumberOfTracks(); 
  if (npart<=0) return;
  if (ntracks<=0) return;
  //
  //
  TParticle * particle= new TParticle;
  TClonesArray * trefs = new TClonesArray("AliTrackReference");
  
  for (Int_t itrack=0;itrack<ntracks;itrack++){
    AliESDtrack *track = fESD->GetTrack(itrack);
    //
    if ((track->GetStatus()&AliESDtrack::kTPCrefit)==0) continue;   // only refited tracks
    Int_t ipart = TMath::Abs(track->GetLabel());
    //
    Int_t status = fMCinfo->GetParticleAndTR(ipart, particle, trefs);
    if (status<0) continue;
    if (!particle) continue;
    if (!trefs) continue;
    //
    //
    AliTrackReference *tpcRef=0;
    for (Int_t iref=0; iref<trefs->GetEntries(); iref++){
      AliTrackReference *ref = (AliTrackReference*)trefs->At(iref);
      if (ref->DetectorId()== AliTrackReference::kTPC){
	tpcRef=ref;
	break;
      }
    }
    if (!tpcRef) continue;
    
    //
    // Fill histos
    //
    Double_t x[7];
    x[0]= track->GetTPCchi2()/track->GetTPCNcls();
    x[1]= track->GetTPCNcls();
    x[2]= track->GetTPCNclsF();
    x[3]= Float_t(track->GetTPCNcls())/Float_t(track->GetTPCNclsF());
    x[4]= tpcRef->Pt();
    x[5]= -0.5*TMath::Log((tpcRef->P()+tpcRef->Pz())/(tpcRef->P()-tpcRef->Pz()));
    x[6]= TMath::ATan2(tpcRef->Y(),tpcRef->X());
    //
    fTPCqa->Fill(x);
  }
}






void AliTPCtaskQA::BinLogX(TAxis *axis) {
  //
  //
  //
  Int_t bins = axis->GetNbins();
  
  Double_t from = axis->GetXmin();
  Double_t to   = axis->GetXmax();
  Double_t *new_bins = new Double_t[bins + 1];
  
  new_bins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (Int_t i = 1; i <= bins; i++) {
    new_bins[i] = factor * new_bins[i-1];
  }
  axis->Set(bins, new_bins);
  delete new_bins;
}


AliTPCtaskQA* AliTPCtaskQA::ReadFromFile(const char *fname){
  //
  //
  //
  AliTPCtaskQA *qa = new AliTPCtaskQA;
  TFile *f = new TFile(fname);
  TObjArray *array= (TObjArray*)f->Get("tpcTaskQA");
  qa->fTPCqa =  (THnSparse*)array->At(0);  
  delete f;
  return qa;
}

