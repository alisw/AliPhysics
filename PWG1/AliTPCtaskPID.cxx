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
//Task for analysis if TPC dEdx and PID information 
//
//-----------------------------------------------------------------------
// Author : M.Ivanov  marian.ivanov@cern.ch - 
//-----------------------------------------------------------------------


// 3 6D histograms  - THnSparse created in the task:
// TPC raw dEdx
// TPC normalized dEdx (dEdx_rec/dNdx_mc)
// TPC PID probabilities
//
// The values are binned in following variables:
// Some of them are correlated - but THnSpase handle it  
//                               ~ 14 MBy per object needed
//
// 0 - MC particle species as defined in the AliPID - negatives+5 <0,9>
// 1 - momenta - at the entrance of the TPC
// 2 - tan lambda- fP[3]
// 3 - betagamma
// 4 - npoints
// 5 - measurement - dEdx, dEdx/BB resp.  PID probability
// 6 - BB resp. rec particle


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
#include "AliTPCtaskPID.h"
//
#include <THnSparse.h>

//

// STL includes
#include <iostream>

using namespace std;


ClassImp(AliTPCtaskPID)

//________________________________________________________________________
AliTPCtaskPID::AliTPCtaskPID() : 
  AliAnalysisTask(), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  fList(0),
  fTPCsignal(0),
  fTPCsignalNorm(0),
  fTPCr(0)
{
  //
  // Default constructor (should not be used)
  //
}

AliTPCtaskPID::AliTPCtaskPID(const AliTPCtaskPID& info) : 
  AliAnalysisTask(info), 
  fMCinfo(info.fMCinfo),     //! MC event handler
  fESD(info.fESD),        //!
  fList(0),
  fTPCsignal(0),
  fTPCsignalNorm(0),
  fTPCr(0)
{
  //
  // Dummy Copy  constructor - no copy constructor for THnSparse 
  //
  fList = (TObjArray*)(info.fList->Clone());
}



//________________________________________________________________________
AliTPCtaskPID::AliTPCtaskPID(const char *name) : 
  AliAnalysisTask(name, "AliTPCtaskPID"), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  fList(0),
  fTPCsignal(0),
  fTPCsignalNorm(0),
  fTPCr(0)
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

void AliTPCtaskPID::Init(){
  //
  // Init dEdx histogram
  // Dimensions
  // 0 - particle specie as defined in the AliPID - negatives+5 <0,9>
  // 1 - momenta - at the entrance of the TPC
  // 2 - tan lambda- fP[3]
  // 3 - betagamma
  // 4 - npoints
  // 5 - measurement - dEdx or dEdx/BB
  // 6 - BB
  //
  Double_t xmin[7],  xmax[7];
  Int_t    nbins[7];
  // pid
  nbins[0]=10;
  xmin[0]=0; xmax[0]=10;
  // momenta
  nbins[1]=50;
  xmin[1]=0.1; xmax[1]=100;
  //pseudorapidity
  nbins[2]=40;
  xmin[2]=-1.4; xmax[2]=1.4;
  //
  // betagamma
  //
  nbins[3]=100;
  xmin[3]=0.1; xmax[3]=1000;
  //
  nbins[4]=11;
  xmin[4] =50; xmax[4]=160;
  //
  // 
  nbins[6]=40;
  xmin[6]=1; xmax[6]=4;

  nbins[5]=400;
  xmin[5]=20; xmax[5]=400;
  fTPCsignal = new THnSparseF("TPC signal","TPC signal",7,nbins,xmin,xmax);
  nbins[5]=100;
  xmin[5]=25; xmax[5]=75;
  fTPCsignalNorm = new THnSparseF("TPC signal Norm","TPC signal Norm",7,nbins,xmin,xmax);
  //
  nbins[5]=256;
  xmin[5]=-0.001; xmax[5]=1.001;
  nbins[6]=10;
  xmin[6]=0; xmax[6]=10;
  fTPCr = new THnSparseF("TPC pid probability ","TPC pid probability",7,nbins,xmin,xmax);
  //
  //
  BinLogX(fTPCsignal->GetAxis(1));
  BinLogX(fTPCsignal->GetAxis(3));
  BinLogX(fTPCsignal->GetAxis(6));
  BinLogX(fTPCsignalNorm->GetAxis(1));
  BinLogX(fTPCsignalNorm->GetAxis(3));
  BinLogX(fTPCsignalNorm->GetAxis(6));
  BinLogX(fTPCr->GetAxis(1));
  BinLogX(fTPCr->GetAxis(3));
  
  char *hisAxisName[7] ={"pid","p (GeV/c)","#eta","#beta#gamma","Number of cluster","dEdx_{rec},dedx_{mc}"};
  char *hisAxisNameNorm[7] ={"pid","p (GeV/c)","#eta","#beta#gamma","Number of cluster","dEdx_{rec}/dEdx_{mc},dedx_{mc}"};
  char *hisAxisNameR[7] ={"pid","p (GeV/c)","#eta","#beta#gamma","Number of cluster","TPCr","pid2"};
  //  
  for (Int_t i=0;i<7;i++) {
    fTPCsignal->GetAxis(i)->SetTitle(hisAxisName[i]);
    fTPCsignal->GetAxis(i)->SetName(hisAxisName[i]);
    fTPCsignalNorm->GetAxis(i)->SetTitle(hisAxisNameNorm[i]);
    fTPCsignalNorm->GetAxis(i)->SetName(hisAxisNameNorm[i]);
    fTPCr->GetAxis(i)->SetTitle(hisAxisNameR[i]);
    fTPCr->GetAxis(i)->SetName(hisAxisNameR[i]);
  }

  fList = new TObjArray(3);
  fList->AddAt(fTPCsignal,0);
  fList->AddAt(fTPCsignalNorm,1);
  fList->AddAt(fTPCr,2);
}




AliTPCtaskPID::~AliTPCtaskPID(){
  //
  //
  //
  delete fTPCsignal;
  delete fTPCsignalNorm;
  delete fTPCr;
}


//________________________________________________________________________
void AliTPCtaskPID::ConnectInputData(Option_t *) 
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
void AliTPCtaskPID::CreateOutputObjects() 
{
  //
  // Connect the output objects
  //

}


//________________________________________________________________________
void AliTPCtaskPID::Exec(Option_t *) {
  //
  // Execute analysis for current event 
  //


  // If MC has been connected   

  if (!fMCinfo){
    cout << "Not MC info\n" << endl;
  }else{
    ProcessMCInfo();
    //mcinfo->Print();
    //DumpInfo();
  }
  //
  PostData(0, fList);
}      




//________________________________________________________________________
void AliTPCtaskPID::Terminate(Option_t *) {
    //
    // Terminate loop
    //
  //
  return;
}







void  AliTPCtaskPID::ProcessMCInfo(){
  //
  //
  //
  //
  if (!fTPCsignal) Init();
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
    const AliExternalTrackParam *in=track->GetInnerParam();
    const AliExternalTrackParam *out=track->GetOuterParam();
    if (!in) continue;
    Int_t ipart = TMath::Abs(track->GetLabel());
    //
    Int_t status = fMCinfo->GetParticleAndTR(ipart, particle, trefs);
    if (status<0) continue;
    if (!particle) continue;
    if (!trefs) continue;
    //
    //
    Double_t mom = in->GetP();
    if (track->GetP()>5)  mom= track->GetP();
    if (out&&out->GetX()<300)  mom= (in->GetP()+out->GetP())*0.5;
    Double_t dedx=track->GetTPCsignal();
    Double_t mass = particle->GetMass();
    Double_t bg  =mom/mass;
    Double_t betheBloch = AliMathBase::BetheBlochAleph(bg);
    //
    // Fill histos
    //
    Double_t x[5];
    //PID
    Int_t pdg = particle->GetPdgCode();
    for (Int_t iType=0;iType<5;iType++) {
      if (AliPID::ParticleCode(iType)==TMath::Abs(pdg)){
	x[0]=iType;
	if (pdg<0) x[0]+=5;
      }
    }
    x[1]= mom;
    x[2]= -0.5*TMath::Log((track->P()+track->Pz())/(track->P()-track->Pz()));
    x[3]= bg;
    x[4]=track->GetTPCsignalN();
    x[5]= dedx;
    x[6]= betheBloch;
    fTPCsignal->Fill(x);
    x[5]=dedx/betheBloch;
    fTPCsignalNorm->Fill(x);
    for (Int_t ipart2=0;ipart2<5;ipart2++){
      x[6]=ipart2;
      Double_t tpcpid[AliPID::kSPECIES];
      track->GetTPCpid(tpcpid);
      x[5]=tpcpid[ipart2];
      fTPCr->Fill(x);
    }
  }
}




void AliTPCtaskPID::FinishTaskOutput()
{
  //
  // According description in AliAnalisysTask this method is call
  // on the slaves before sending data
  //
  Terminate("slave");
  gSystem->Exec("pwd");
}




void AliTPCtaskPID::BinLogX(TAxis *axis) {
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
