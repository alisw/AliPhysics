/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------
//         AliAnalysisTaskSpectraAllChNanoAOD class
// NanoAOD edition: this is only meant to test the developement
// version of NanoAODs
//
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliVParticle.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskSpectraAllChNanoAOD.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliSpectraAODTrackCuts.h"
#include "AliSpectraAODEventCuts.h"
#include "AliHelperPID.h"
#include "AliPIDCombined.h"
#include "AliCentrality.h"
#include "TProof.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include <TMCProcess.h>

#include <iostream>
#include "AliNanoAODHeader.h"
#include "AliNanoAODTrack.h"

using namespace AliHelperPIDNameSpace;
using namespace std;

ClassImp(AliAnalysisTaskSpectraAllChNanoAOD)

//________________________________________________________________________
AliAnalysisTaskSpectraAllChNanoAOD::AliAnalysisTaskSpectraAllChNanoAOD(const char *name) : AliAnalysisTaskSE(name),  
  fAOD(0x0),
  fTrackCuts(0x0),
  fEventCuts(0x0),
  fHelperPID(0x0),
  fIsMC(0),
  fDoDoubleCounting(0),
  fFillOnlyEvents(0),
  fCharge(0),
  fVZEROside(0),
  fOutput(0x0),
  fnCentBins(20),
  fnQvecBins(40),
  fnNchBins(200)
{
  // Default constructor
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliSpectraAODEventCuts::Class());
  DefineOutput(3, AliSpectraAODTrackCuts::Class());
  DefineOutput(4, AliHelperPID::Class());
}

//________________________________________________________________________
void AliAnalysisTaskSpectraAllChNanoAOD::UserCreateOutputObjects()
{
  // create output objects
  fOutput=new TList();
  fOutput->SetOwner();
  fOutput->SetName("fOutput");
  
  if (!fTrackCuts) AliFatal("Track Cuts should be set in the steering macro");
  if (!fEventCuts) AliFatal("Event Cuts should be set in the steering macro");
  if (!fHelperPID)  AliFatal("HelperPID should be set in the steering macro");
  
  // binning common to all the THn
  const Double_t ptBins[] = {0.20,0.30,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.2,3.6,4.0,5.0,6.0,7.0,8.0,9.0,10.,12.,15.};
  const Int_t nptBins=26;
  
  //dimensions of THnSparse for tracks
  const Int_t nvartrk=8;
  //                                             pt          cent        Q vec     IDrec     IDgen       isph           iswd      y
  Int_t    binsHistRealTrk[nvartrk] = {      nptBins, fnCentBins,   fnQvecBins,        4,        3,         2,          2,       2};
  Double_t xminHistRealTrk[nvartrk] = {         0.,          0.,            0.,      -.5,      -0.5,      -0.5,        -0.5,   -0.5};
  Double_t xmaxHistRealTrk[nvartrk] = {       10.,       100.,            8.,      3.5,      2.5,       1.5,         1.5,     0.5};    
  THnSparseF* NSparseHistTrk = new THnSparseF("NSparseHistTrk","NSparseHistTrk",nvartrk,binsHistRealTrk,xminHistRealTrk,xmaxHistRealTrk);
  NSparseHistTrk->GetAxis(0)->SetTitle("#it{p}_{T,rec}");
  NSparseHistTrk->GetAxis(0)->SetName("pT_rec");
  NSparseHistTrk->SetBinEdges(0,ptBins);
  NSparseHistTrk->GetAxis(1)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  NSparseHistTrk->GetAxis(1)->SetName(Form("%s_cent",fEventCuts->GetCentralityMethod().Data()));
  NSparseHistTrk->GetAxis(2)->SetTitle("Q vec");
  NSparseHistTrk->GetAxis(2)->SetName("Q_vec");
  NSparseHistTrk->GetAxis(3)->SetTitle("ID rec");
  NSparseHistTrk->GetAxis(3)->SetName("ID_rec");
  NSparseHistTrk->GetAxis(4)->SetTitle("ID gen");
  NSparseHistTrk->GetAxis(4)->SetName("ID_gen");
  NSparseHistTrk->GetAxis(5)->SetTitle("isph");
  NSparseHistTrk->GetAxis(5)->SetName("isph");
  NSparseHistTrk->GetAxis(6)->SetTitle("iswd");
  NSparseHistTrk->GetAxis(6)->SetName("iswd");
  NSparseHistTrk->GetAxis(7)->SetTitle("y");
  NSparseHistTrk->GetAxis(7)->SetName("y");
  fOutput->Add(NSparseHistTrk);
  
  //dimensions of THnSparse for stack
  const Int_t nvarst=5;
  //                                             pt          cent    IDgen        isph        y
  Int_t    binsHistRealSt[nvarst] = {      nptBins,  fnCentBins,        3,         2,        2};
  Double_t xminHistRealSt[nvarst] = {         0.,           0.,      -0.5,      -0.5,    -0.5};
  Double_t xmaxHistRealSt[nvarst] = {       10.,        100.,      2.5,       1.5,      0.5};
  THnSparseF* NSparseHistSt = new THnSparseF("NSparseHistSt","NSparseHistSt",nvarst,binsHistRealSt,xminHistRealSt,xmaxHistRealSt);
  NSparseHistSt->GetAxis(0)->SetTitle("#it{p}_{T,gen}");
  NSparseHistSt->SetBinEdges(0,ptBins);
  NSparseHistSt->GetAxis(0)->SetName("pT_rec");
  NSparseHistSt->GetAxis(1)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  NSparseHistSt->GetAxis(1)->SetName(Form("%s_cent",fEventCuts->GetCentralityMethod().Data()));
  NSparseHistSt->GetAxis(2)->SetTitle("ID gen");
  NSparseHistSt->GetAxis(2)->SetName("ID_gen");
  NSparseHistSt->GetAxis(3)->SetTitle("isph");
  NSparseHistSt->GetAxis(3)->SetName("isph");
  NSparseHistSt->GetAxis(4)->SetTitle("y");
  NSparseHistSt->GetAxis(4)->SetName("y");
  fOutput->Add(NSparseHistSt);
  
  //dimensions of THnSparse for the normalization
  const Int_t nvarev=3;
  //                                             cent             Q vec                Nch
  Int_t    binsHistRealEv[nvarev] = {    fnCentBins,      fnQvecBins,           fnNchBins};
  Double_t xminHistRealEv[nvarev] = {           0.,               0.,                   0.};
  Double_t xmaxHistRealEv[nvarev] = {       100.,               8.,               2000.};
  THnSparseF* NSparseHistEv = new THnSparseF("NSparseHistEv","NSparseHistEv",nvarev,binsHistRealEv,xminHistRealEv,xmaxHistRealEv);
  NSparseHistEv->GetAxis(0)->SetTitle(Form("%s cent",fEventCuts->GetCentralityMethod().Data()));
  NSparseHistEv->GetAxis(0)->SetName(Form("%s_cent",fEventCuts->GetCentralityMethod().Data()));
  NSparseHistEv->GetAxis(1)->SetTitle("Q vec");
  NSparseHistEv->GetAxis(1)->SetName("Q_vec");
  NSparseHistEv->GetAxis(2)->SetTitle("N charged");
  NSparseHistEv->GetAxis(2)->SetName("N_ch");
  fOutput->Add(NSparseHistEv);
  
  PostData(1, fOutput  );
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);
  PostData(4, fHelperPID);
}

//________________________________________________________________________

void AliAnalysisTaskSpectraAllChNanoAOD::UserExec(Option_t *)
{
  //Printf("An event");
  // main event loop
  fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
  if (!fAOD) {
    AliWarning("ERROR: AliAODEvent not available \n");
    return;
  }
  
  if (strcmp(fAOD->ClassName(), "AliAODEvent"))
    {
      AliFatal("Not processing AODs");
    } 

  AliNanoAODHeader * headNano = dynamic_cast<AliNanoAODHeader*>((TObject*)fAOD->GetHeader());
  
  Bool_t isNano = (headNano != 0);
 
  if(!isNano) {
    if(!fEventCuts->IsSelected(fAOD,fTrackCuts))return;//event selection
  }
  //Default TPC priors
  if(fHelperPID->GetPIDType()==kBayes)fHelperPID->GetPIDCombined()->SetDefaultTPCPriors();//FIXME maybe this can go in the UserCreateOutputObject?
  
  Double_t Qvec=0.;//in case of MC we save space in the memory
  if(!fIsMC){
    if(!isNano) {
      if(fVZEROside==0)Qvec=fEventCuts->GetqV0A();
      else if (fVZEROside==1)Qvec=fEventCuts->GetqV0C();
    } else {
      
      if(fVZEROside==0)      Qvec=headNano->GetVar(0); // FIXME: we need proper getters here
      else if (fVZEROside==1)Qvec=headNano->GetVar(1);
      
    }
  }

  Double_t Cent=isNano ? headNano->GetVar(2) : fEventCuts->GetCent();
    
  // First do MC to fill up the MC particle array
  TClonesArray *arrayMC = 0;
  if (fIsMC)
    {
      arrayMC = (TClonesArray*) fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
      if (!arrayMC) {
	AliFatal("Error: MC particles branch not found!\n");
      }
      Int_t nMC = arrayMC->GetEntries();
      for (Int_t iMC = 0; iMC < nMC; iMC++)
	{
	  AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(iMC);
	  if(!partMC->Charge()) continue;//Skip neutrals
	  if(fCharge != 0 && partMC->Charge()*fCharge < 0.) continue;//if fCharge != 0 only select fCharge
	  
	  if(partMC->Eta() < fTrackCuts->GetEtaMin() || partMC->Eta() > fTrackCuts->GetEtaMax())continue;//ETA CUT ON GENERATED!!!!!!!!!!!!!!!!!!!!!!!!!!
	  
	  //pt     cent    IDgen        isph        y
 	  Double_t varSt[5];
	  varSt[0]=partMC->Pt();
	  varSt[1]=Cent;
	  varSt[2]=(Double_t)fHelperPID->GetParticleSpecies(partMC);
	  varSt[3]=(Double_t)partMC->IsPhysicalPrimary();
	  varSt[4]=partMC->Y();
	  ((THnSparseF*)fOutput->FindObject("NSparseHistSt"))->Fill(varSt);//stack loop
	  
	  //Printf("a particle");
	  
	}
    }
  
  //main loop on tracks
  
  Int_t Nch = 0.;
  
  for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
    AliVTrack* track = (AliVTrack*) fAOD->GetTrack(iTracks);
    if(fCharge != 0 && track->Charge() != fCharge) continue;//if fCharge != 0 only select fCharge 
    if(!isNano) {
      if (!fTrackCuts->IsSelected((AliAODTrack*)track,kTRUE)) continue; //track selection (rapidity selection NOT in the standard cuts)
    }
    
    if(!fFillOnlyEvents){
      Int_t IDrec=isNano ? GetNanoTrackID (track) : fHelperPID->GetParticleSpecies(track,kTRUE);//id from detector      
      Double_t y= 0;
      if(isNano) y = ((AliNanoAODTrack*)track)->Y(fHelperPID->GetMass((AliHelperParticleSpecies_t)IDrec));
      else y = ((AliAODTrack*)track)->Y(fHelperPID->GetMass((AliHelperParticleSpecies_t)IDrec));
      Int_t IDgen=kSpUndefined;//set if MC
      Int_t isph=-999;
      Int_t iswd=-999;
      
      if (arrayMC) {
	AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(TMath::Abs(track->GetLabel()));
	if (!partMC) { 
	  AliError("Cannot get MC particle");
	  continue; 
	}
	IDgen=fHelperPID->GetParticleSpecies(partMC);
	isph=partMC->IsPhysicalPrimary();
	iswd=partMC->IsSecondaryFromWeakDecay();//FIXME not working on old productions
      }
      
    //pt     cent    Q vec     IDrec     IDgen       isph           iswd      y
      Double_t varTrk[8];
      varTrk[0]=track->Pt();
      varTrk[1]=Cent;
      varTrk[2]=Qvec;
      varTrk[3]=(Double_t)IDrec;
      varTrk[4]=(Double_t)IDgen;
      varTrk[5]=(Double_t)isph;
      varTrk[6]=(Double_t)iswd;
      varTrk[7]=y;
      ((THnSparseF*)fOutput->FindObject("NSparseHistTrk"))->Fill(varTrk);//track loop
      
      //for nsigma PID fill double counting of ID
      if(fHelperPID->GetPIDType()<kBayes && fDoDoubleCounting){//only nsigma
	Bool_t *HasDC;
	HasDC=fHelperPID->GetDoubleCounting(track,kTRUE);//get the array with double counting
	for(Int_t ipart=0;ipart<kNSpecies;ipart++){
	  if(HasDC[ipart]==kTRUE){
	    varTrk[3]=(Double_t)ipart;
	    ((THnSparseF*)fOutput->FindObject("NSparseHistTrk"))->Fill(varTrk);//track loop
	  }
	}
      }
      
      //fill all charged (3)
      varTrk[3]=3.;
      varTrk[4]=3.;
      ((THnSparseF*)fOutput->FindObject("NSparseHistTrk"))->Fill(varTrk);//track loop
    }//end if fFillOnlyEvents
    
    //Printf("a track");
      
    Nch++;
  } // end loop on tracks
  
  Double_t varEv[3];
  varEv[0]=Cent;
  varEv[1]=Qvec;
  varEv[2]=Nch;
  ((THnSparseF*)fOutput->FindObject("NSparseHistEv"))->Fill(varEv);//event loop
  
  PostData(1, fOutput  );
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);
  PostData(4, fHelperPID);
}

//_________________________________________________________________
void   AliAnalysisTaskSpectraAllChNanoAOD::Terminate(Option_t *)
{
  // Terminate
}


Int_t AliAnalysisTaskSpectraAllChNanoAOD::GetNanoTrackID(AliVTrack * track) {
  // Applies nsigma PID to nano tracks
  AliNanoAODTrack * nanoTrack = dynamic_cast<AliNanoAODTrack*>(track);
  if(!nanoTrack) AliFatal("Not a nano AOD track");

  // Cache indexes
  static const Int_t kcstNSigmaTPCPi  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstNSigmaTPCPi");
  static const Int_t kcstNSigmaTPCKa  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstNSigmaTPCKa");
  static const Int_t kcstNSigmaTPCPr  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstNSigmaTPCPr");
  static const Int_t kcstNSigmaTOFPi  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstNSigmaTOFPi");
  static const Int_t kcstNSigmaTOFKa  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstNSigmaTOFKa");
  static const Int_t kcstNSigmaTOFPr  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstNSigmaTOFPr");

  Double_t nSigmaPID = 3.0;



  //get the identity of the particle with the minimum Nsigma
  Double_t nsigmaPion=999., nsigmaKaon=999., nsigmaProton=999.;
  if(nanoTrack->Pt() > fTrackCuts->GetPtTOFMatching()) {
    nsigmaProton =  TMath::Sqrt(nanoTrack->GetVar(kcstNSigmaTPCPr)*nanoTrack->GetVar(kcstNSigmaTPCPr)+nanoTrack->GetVar(kcstNSigmaTOFPr)*nanoTrack->GetVar(kcstNSigmaTOFPr));
    nsigmaKaon   =  TMath::Sqrt(nanoTrack->GetVar(kcstNSigmaTPCKa)*nanoTrack->GetVar(kcstNSigmaTPCKa)+nanoTrack->GetVar(kcstNSigmaTOFKa)*nanoTrack->GetVar(kcstNSigmaTOFKa));
    nsigmaPion   =  TMath::Sqrt(nanoTrack->GetVar(kcstNSigmaTPCPi)*nanoTrack->GetVar(kcstNSigmaTPCPi)+nanoTrack->GetVar(kcstNSigmaTOFPi)*nanoTrack->GetVar(kcstNSigmaTOFPi));
  }
  else {
    nsigmaProton =  TMath::Abs(nanoTrack->GetVar(kcstNSigmaTPCPr));
    nsigmaKaon   =  TMath::Abs(nanoTrack->GetVar(kcstNSigmaTPCKa));  
    nsigmaPion   =  TMath::Abs(nanoTrack->GetVar(kcstNSigmaTPCPi));  
  }

  

  // guess the particle based on the smaller nsigma (within nSigmaPID)
  if( ( nsigmaKaon==nsigmaPion ) && ( nsigmaKaon==nsigmaProton )) return kSpUndefined;//if is the default value for the three
  
  if( ( nsigmaKaon   < nsigmaPion ) && ( nsigmaKaon < nsigmaProton ) && (nsigmaKaon   < nSigmaPID)){
    return kSpKaon;
  }
  if( ( nsigmaPion   < nsigmaKaon ) && ( nsigmaPion < nsigmaProton ) && (nsigmaPion   < nSigmaPID)){    
    return kSpPion;
  }
  if( ( nsigmaProton < nsigmaKaon ) && ( nsigmaProton < nsigmaPion ) && (nsigmaProton < nSigmaPID)){
    return kSpProton;
  }
  // else, return undefined
  return kSpUndefined;
  

}
