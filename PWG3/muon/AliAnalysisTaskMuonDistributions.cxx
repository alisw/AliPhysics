/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------------
// Analysis task to compute muon/dimuon kinematic distributions
// The output is a list of histograms.
// The macro class can run on AOD or in the train with the ESD filter.
// R. Arnaldi
//
//-----------------------------------------------------------------------------

//#ifndef ALIANALYSISTASKMUONDISTRIBUTIONS_CXX
//#define ALIANALYSISTASKMUONDISTRIBUTIONS_CXX

#include <TChain.h>
#include <TTree.h>
#include <TList.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TH1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TPaveText.h>

#include "AliAnalysisTaskMuonDistributions.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliVEvent.h"
#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliHeader.h"
#include "AliESDHeader.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "AliESDMuonTrack.h"
#include "AliESDtrack.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"

ClassImp(AliAnalysisTaskMuonDistributions)

//__________________________________________________________________________
AliAnalysisTaskMuonDistributions::AliAnalysisTaskMuonDistributions() :
  fBeamEnergy(0.),
  fInvMassFitLimitMin(2.),
  fInvMassFitLimitMax(5.),
  fPsiFitLimitMin(2.9),
  fPsiFitLimitMax(3.3),
  fBckFitLimitMin(2.2),
  fBckFitLimitMax(2.85),
  fInvariantMassFit(kFALSE),
  fAnalysisType(0x0),
  fOutput(0x0)
{
}
//___________________________________________________________________________
AliAnalysisTaskMuonDistributions::AliAnalysisTaskMuonDistributions(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fBeamEnergy(0.),
  fInvMassFitLimitMin(2.),
  fInvMassFitLimitMax(5.),
  fPsiFitLimitMin(2.9),
  fPsiFitLimitMax(3.3),
  fBckFitLimitMin(2.2),
  fBckFitLimitMax(2.85),
  fInvariantMassFit(kFALSE),
  fAnalysisType(0x0),
  fOutput(0x0)
{
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskMuonDistributions","Calling Constructor");
  
  DefineOutput(1,TList::Class());

}

//___________________________________________________________________________
AliAnalysisTaskMuonDistributions& AliAnalysisTaskMuonDistributions::operator=(const AliAnalysisTaskMuonDistributions& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskMuonDistributions::AliAnalysisTaskMuonDistributions(const AliAnalysisTaskMuonDistributions& c) :
  AliAnalysisTaskSE(c),
  fBeamEnergy(c.fBeamEnergy),
  fInvMassFitLimitMin(c.fInvMassFitLimitMin),
  fInvMassFitLimitMax(c.fInvMassFitLimitMax),
  fPsiFitLimitMin(c.fPsiFitLimitMin),
  fPsiFitLimitMax(c.fPsiFitLimitMax),
  fBckFitLimitMin(c.fBckFitLimitMin),
  fBckFitLimitMax(c.fBckFitLimitMax),
  fInvariantMassFit(c.fInvariantMassFit),
  fAnalysisType(c.fAnalysisType),
  fOutput(c.fOutput)
 {
  //
  // Copy Constructor										
  //
}

//___________________________________________________________________________
AliAnalysisTaskMuonDistributions::~AliAnalysisTaskMuonDistributions() {
  //
  //destructor
  //
  Info("~AliAnalysisTaskMuonDistributions","Calling Destructor");
}

//___________________________________________________________________________
void AliAnalysisTaskMuonDistributions::UserCreateOutputObjects(){
	 
 fOutput = new TList();
 fOutput->SetOwner(); 
 //
 // various histos
 //
 TH1D *hNumberMuonTracks = new TH1D("hNumberMuonTracks","hNumberMuonTracks;N_{#mu tracks}",10,0.,10.);
 //
 // dimuon histos
 //
 TH1D *hMass_Dimu   = new TH1D("hMass_Dimu","hMass_Dimu;M_{#mu#mu} (GeV/c^{2})",180,0,9.);	
 TH1D *hPt_Dimu  = new TH1D("hPt_Dimu","hPt_Dimu;p_{T} (GeV/c)",100,0,20);	
 TH1D *hRapidity_Dimu  = new TH1D("hRapidity_Dimu","hRapidity_Dimu;y",100,-5,-2);	
 TH1D *hCosThetaCS_Dimu  = new TH1D("hCosThetaCS_Dimu","hCosThetaCS_Dimu;cos#theta_{CS}",100,-1.,1.);	
 TH1D *hCosThetaHE_Dimu  = new TH1D("hCosThetaHE_Dimu","hCosThetaHE_Dimu;cos#theta_{HE}",100,-1.,1.);	
 //
 // muon histos
 //
 TH1D *hP  = new TH1D("hP","hP;p (GeV/c)",100,0,500);	
 TH1D *hPt  = new TH1D("hPt","hPt;p_{T} (GeV/c)",100,0,20);	
 TH1D *hRapidity  = new TH1D("hRapidity","hRapidity;y",100,-5,-2);	
	
 fOutput->Add(hNumberMuonTracks); 	
 fOutput->Add(hMass_Dimu); 	
 fOutput->Add(hPt_Dimu); 	
 fOutput->Add(hRapidity_Dimu); 	
 fOutput->Add(hCosThetaCS_Dimu); 	
 fOutput->Add(hCosThetaHE_Dimu); 	
 fOutput->Add(hP); 	
 fOutput->Add(hPt); 	
 fOutput->Add(hRapidity); 	
 fOutput->ls(); 
} 

//_________________________________________________
void AliAnalysisTaskMuonDistributions::UserExec(Option_t *)
{
  AliESDEvent *esd=0x0;
  AliAODEvent *aod=0x0;
  
  if(strcmp(fAnalysisType,"ESD")==0){
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>
        (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    esd = esdH->GetEvent();
  } else if(strcmp(fAnalysisType,"AOD")==0){
    aod = dynamic_cast<AliAODEvent*> (InputEvent());
  }
  
  Int_t ntracks=-999;
  if(strcmp(fAnalysisType,"ESD")==0) ntracks=esd->GetNumberOfMuonTracks();
  else if(strcmp(fAnalysisType,"AOD")==0) ntracks=aod->GetNumberOfTracks();
  Int_t nmuontracks=0;
  
  for (Int_t j = 0; j<ntracks; j++) {
    Float_t px_mu1=-999; Float_t py_mu1=-999; Float_t pz_mu1=-999; Float_t pt_mu1=-999; Float_t p_mu1=-999;
    Float_t e_mu1=-999; Float_t rapidity_mu1=-999;  Float_t charge_mu1=-999; 
    if(strcmp(fAnalysisType,"ESD")==0){ 
      AliESDMuonTrack* mu1 = new AliESDMuonTrack(*(esd->GetMuonTrack(j)));
      if (!mu1->ContainTrackerData()) continue;
      charge_mu1 = mu1->Charge();
      px_mu1 = mu1->Px();
      py_mu1 = mu1->Py();
      pz_mu1 = mu1->Pz();
      e_mu1 = mu1->E();
      p_mu1 = mu1->P();
      pt_mu1 = mu1->Pt();
      rapidity_mu1 = Rapidity(e_mu1,pz_mu1);
    } else if(strcmp(fAnalysisType,"AOD")==0){
      AliAODTrack *mu1 = aod->GetTrack(j);
      if(!mu1->IsMuonTrack()) continue;
      charge_mu1 = mu1->Charge();
      px_mu1 = mu1->Px();
      py_mu1 = mu1->Py();
      pz_mu1 = mu1->Pz();
      e_mu1 = mu1->E();
      p_mu1 = mu1->P();
      pt_mu1 = mu1->Pt();
      rapidity_mu1 = Rapidity(e_mu1,pz_mu1);
    }
    ((TH1D*)(fOutput->FindObject("hP")))->Fill(p_mu1);    
    ((TH1D*)(fOutput->FindObject("hPt")))->Fill(pt_mu1);    
    ((TH1D*)(fOutput->FindObject("hRapidity")))->Fill(rapidity_mu1);	  
    nmuontracks++;
    if(charge_mu1<0){
      for (Int_t jj = 0; jj<ntracks; jj++) {
        Float_t px_mu2=-999; Float_t py_mu2=-999; Float_t pz_mu2=-999;
        Float_t e_mu2=-999;Float_t charge_mu2=-999; 
        if(strcmp(fAnalysisType,"ESD")==0){ 
          AliESDMuonTrack* mu2 = new AliESDMuonTrack(*(esd->GetMuonTrack(jj)));
          if (!mu2->ContainTrackerData()) continue;
	  charge_mu2 = mu2->Charge();
          px_mu2 = mu2->Px();
          py_mu2 = mu2->Py();
          pz_mu2 = mu2->Pz();
	  e_mu2 = mu2->E();
        } else if(strcmp(fAnalysisType,"AOD")==0){
          AliAODTrack *mu2 = aod->GetTrack(jj);
          if(!mu2->IsMuonTrack()) continue; 
	  charge_mu2 = mu2->Charge();
          px_mu2 = mu2->Px();
          py_mu2 = mu2->Py();
          pz_mu2 = mu2->Pz();
	  e_mu2 = mu2->E();
        }
        if(charge_mu2>0){
	  Float_t pt_dimu = TMath::Sqrt((px_mu1+px_mu2)*(px_mu1+px_mu2)+(py_mu1+py_mu2)*(py_mu1+py_mu2));
 	  Float_t mass_dimu = InvMass(e_mu1,px_mu1,py_mu1,pz_mu1,e_mu2,px_mu2,py_mu2,pz_mu2);
 	  Float_t rapidity_dimu = Rapidity((e_mu1+e_mu2),(pz_mu1+pz_mu2));
	  Double_t costhetaCS_dimu = CostCS((Double_t) px_mu1,(Double_t) py_mu1,(Double_t)pz_mu1,(Double_t) e_mu1,(Double_t)charge_mu1,(Double_t) px_mu2,(Double_t) py_mu2,(Double_t)pz_mu2,(Double_t) e_mu2);
	  Double_t costhetaHE_dimu = CostHE((Double_t) px_mu1,(Double_t) py_mu1,(Double_t)pz_mu1,(Double_t) e_mu1,(Double_t)charge_mu1,(Double_t) px_mu2,(Double_t) py_mu2,(Double_t)pz_mu2,(Double_t) e_mu2);
	  ((TH1D*)(fOutput->FindObject("hMass_Dimu")))->Fill(mass_dimu);
	  ((TH1D*)(fOutput->FindObject("hPt_Dimu")))->Fill(pt_dimu);	
	  ((TH1D*)(fOutput->FindObject("hRapidity_Dimu")))->Fill(rapidity_dimu);	
	  ((TH1D*)(fOutput->FindObject("hCosThetaCS_Dimu")))->Fill(costhetaCS_dimu);	
	  ((TH1D*)(fOutput->FindObject("hCosThetaHE_Dimu")))->Fill(costhetaHE_dimu);	
        }
        //delete mu2;
      }      // second mu Loop
    }          // mu- Selection
    //delete mu1;
  }        
  ((TH1D*)(fOutput->FindObject("hNumberMuonTracks")))->Fill(nmuontracks); 
  
  PostData(1,fOutput);
  }


//________________________________________________________________________
void AliAnalysisTaskMuonDistributions::Terminate(Option_t *) 
{
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  Int_t xmin=20; 
  Int_t ymin=20;
  
  printf("Using beam Energy=%f \n",fBeamEnergy);

  TH1D *h_NumberMuonTracks = dynamic_cast<TH1D*> (fOutput->FindObject("hNumberMuonTracks"));  
  TH1D *h_Mass_Dimu = dynamic_cast<TH1D*> (fOutput->FindObject("hMass_Dimu"));  
  TH1D *h_Pt_Dimu = dynamic_cast<TH1D*> (fOutput->FindObject("hPt_Dimu"));  
  TH1D *h_Rapidity_Dimu = dynamic_cast<TH1D*> (fOutput->FindObject("hRapidity_Dimu"));  
  TH1D *h_CostCS_Dimu = dynamic_cast<TH1D*> (fOutput->FindObject("hCosThetaCS_Dimu"));  
  TH1D *h_CostHE_Dimu = dynamic_cast<TH1D*> (fOutput->FindObject("hCosThetaHE_Dimu"));  
  TH1D *h_P = dynamic_cast<TH1D*> (fOutput->FindObject("hP"));  
  TH1D *h_Pt = dynamic_cast<TH1D*> (fOutput->FindObject("hPt"));  
  TH1D *h_Rapidity = dynamic_cast<TH1D*> (fOutput->FindObject("hRapidity"));  

  TCanvas *c0 = new TCanvas("c0","Plots",xmin,ymin,600,600);
  c0->Divide(2,2);
  c0->cd(1);
  h_NumberMuonTracks->Draw();
  
  xmin+=20; ymin+=20;
  TCanvas *c1 = new TCanvas("c1","Muon kinematic distributions Plots",xmin,ymin,600,600);
  c1->Divide(2,2);  
  c1->cd(1);
  gPad->SetLogy(1);
  h_P->Draw();
  c1->cd(2);
  gPad->SetLogy(1);
  h_Pt->Draw();
  c1->cd(3);
  h_Rapidity->Draw();
  
  xmin+=20; ymin+=20;
  TCanvas *c2 = new TCanvas("c2","Dimuon kinematic distributions Plots",xmin,ymin,600,600);
  c2->Divide(2,2);  
  c2->cd(1);
  gPad->SetLogy(1);
  h_Pt_Dimu->Draw();
  c2->cd(2);
  h_Rapidity_Dimu->Draw();
  c2->cd(3);
  h_CostCS_Dimu->Draw();
  c2->cd(4);
  h_CostHE_Dimu->Draw();
  
  xmin+=20; ymin+=20;
  TCanvas *c3 = new TCanvas("c3","Invariant Mass Plots",xmin,ymin,600,600);
  gPad->SetLogy(1);
  h_Mass_Dimu->Draw();  
  if(fInvariantMassFit) FitInvMass(h_Mass_Dimu);
  c3->Update();

}

//________________________________________________________________________
Float_t AliAnalysisTaskMuonDistributions::InvMass(Float_t e1, Float_t px1, Float_t py1, Float_t pz1,
				   Float_t e2, Float_t px2, Float_t py2, Float_t pz2) 
{
// invariant mass calculation
    Float_t imassrec = TMath::Sqrt((e1+e2)*(e1+e2)-((px1+px2)*(px1+px2)+
                                    (py1+py2)*(py1+py2)+(pz1+pz2)*(pz1+pz2)));
    return imassrec;
}
//________________________________________________________________________
Float_t AliAnalysisTaskMuonDistributions::Rapidity(Float_t e, Float_t pz) 
{
// calculate rapidity
    Float_t rap;
    if(e!=pz){
      rap = 0.5*TMath::Log((e+pz)/(e-pz));
      return rap;
    } else {
      rap = -200;
      return rap;
    }
}
//________________________________________________________________________
Double_t AliAnalysisTaskMuonDistributions::CostCS(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2)
{
  TLorentzVector PMu1CM, PMu2CM, PProjCM, PTargCM, PDimuCM; // In the CM. frame
  TLorentzVector PMu1Dimu, PMu2Dimu, PProjDimu, PTargDimu; // In the dimuon rest frame
  TVector3 beta,zaxisCS;
  Double_t mp=0.93827231;
  //
  // --- Fill the Lorentz vector for projectile and target in the CM frame
  //
  PProjCM.SetPxPyPzE(0.,0.,-fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+mp*mp)); 
  PTargCM.SetPxPyPzE(0.,0.,fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+mp*mp)); 
  //
  // --- Get the muons parameters in the CM frame 
  //
  PMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
  PMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  PDimuCM=PMu1CM+PMu2CM;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  beta=(-1./PDimuCM.E())*PDimuCM.Vect();
  PMu1Dimu=PMu1CM;
  PMu2Dimu=PMu2CM;
  PProjDimu=PProjCM;
  PTargDimu=PTargCM;
  PMu1Dimu.Boost(beta);
  PMu2Dimu.Boost(beta);
  PProjDimu.Boost(beta);
  PTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the CS angle 
  //
  zaxisCS=(((PProjDimu.Vect()).Unit())-((PTargDimu.Vect()).Unit())).Unit();
  //				     
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  Double_t cost;
  
  if(charge1>0) {cost = zaxisCS.Dot((PMu1Dimu.Vect()).Unit());}
  else {cost = zaxisCS.Dot((PMu2Dimu.Vect()).Unit());}
  return cost;
}

//________________________________________________________________________
Double_t AliAnalysisTaskMuonDistributions::CostHE(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2)
{
  TLorentzVector PMu1CM, PMu2CM, PDimuCM; // In the CM frame 
  TLorentzVector PMu1Dimu, PMu2Dimu; // In the dimuon rest frame
  TVector3 beta,zaxisCS;
  //
  // --- Get the muons parameters in the CM frame
  //
  PMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
  PMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  PDimuCM=PMu1CM+PMu2CM;
  //
  // --- Translate the muon parameters in the dimuon rest frame
  //
  beta=(-1./PDimuCM.E())*PDimuCM.Vect();
  PMu1Dimu=PMu1CM;
  PMu2Dimu=PMu2CM;
  PMu1Dimu.Boost(beta);
  PMu2Dimu.Boost(beta);
  //
  // --- Determine the z axis for the calculation of the polarization angle (i.e. the direction of the dimuon in the CM system)
  //
  TVector3 zaxis;
  zaxis=(PDimuCM.Vect()).Unit();
  
  // --- Calculation of the polarization angle (angle between mu+ and the z axis defined above)
  Double_t cost;
  if(charge1>0) {cost = zaxis.Dot((PMu1Dimu.Vect()).Unit());} 
  else {cost = zaxis.Dot((PMu2Dimu.Vect()).Unit());} 
  return cost;
}  

//________________________________________________________________________
void AliAnalysisTaskMuonDistributions::FitInvMass(TH1D *histo)
{
  TF1 *gau = new TF1("gau","gaus",fPsiFitLimitMin,fPsiFitLimitMax);
  TF1 *ex = new TF1("ex","expo",fBckFitLimitMin,fBckFitLimitMax);    
  TF1 *tot = new TF1("mtot","gaus(0)+expo(3)",fInvMassFitLimitMin,fInvMassFitLimitMax);
  Double_t par[5];
  Double_t BinWidth= histo->GetBinWidth(1);
  printf("BW=%f\n",BinWidth);
  gau->SetLineColor(3);
  gau->SetLineWidth(2);
  histo->Fit(gau,"RlQ"); 
  ex->SetLineColor(4);
  ex->SetLineWidth(2);
  histo->Fit(ex,"RlQ+");
  gau->GetParameters(&par[0]);
  ex->GetParameters(&par[3]);
  tot->SetParameters(par);   
  tot->SetLineColor(2);
  tot->SetLineWidth(2);
  histo->Fit(tot,"Rl+");
  histo->Draw("e");
  Double_t Chi2 = tot->GetChisquare();
  Double_t NDF = tot->GetNDF();
  Float_t MeanPsi= tot->GetParameter(1);
  Float_t SigPsi= tot->GetParameter(2)*1000.;
  Double_t NPsiFit = TMath::Sqrt(2*3.1415)*tot->GetParameter(0)*tot->GetParameter(2)/BinWidth;
  TF1 *exfix = new TF1("exfix","expo",2.,5.);  
  exfix->SetParameter(0,tot->GetParameter(3));  
  exfix->SetParameter(1,tot->GetParameter(4));  
  Double_t NBck = exfix->Integral(2.9,3.3)/BinWidth;
  
  printf("\n\n****************************************************************************\n");
  char psi_text[100];
  sprintf(psi_text,"N. J/#psi (2.9-3.3)=%10.0f",NPsiFit);
  printf("\nN. J/psi (2.9-3.3)=%10.0f\n",NPsiFit);
  TLatex *psi_latex = new TLatex(4.5,0.85*histo->GetMaximum(),psi_text);
  psi_latex->SetTextColor(2);
  psi_latex->SetTextSize(0.03);
  psi_latex->SetTextAlign(2);
  psi_latex->Draw();
  
  char psi2_text[100];
  sprintf(psi2_text,"J/#psi m=%4.3f GeV #sigma=%4.2f MeV",MeanPsi,SigPsi);
  printf("J/psi m=%4.3f GeV sigma=%4.2f MeV\n",MeanPsi,SigPsi);
  TLatex *psi2_latex = new TLatex(4.5,0.425*histo->GetMaximum(),psi2_text);
  psi2_latex->SetTextColor(2);
  psi2_latex->SetTextSize(0.03);
  psi2_latex->SetTextAlign(2);
  psi2_latex->Draw();
  
  char sb_text[100];
  sprintf(sb_text,"S/B (2.9-3.3)=%4.2f ",NPsiFit/NBck);
  printf("S/B (2.9-3.3)=%4.2f\n",NPsiFit/NBck);
  TLatex *sb_latex = new TLatex(4.5,0.212*histo->GetMaximum(),sb_text);
  sb_latex->SetTextColor(2);
  sb_latex->SetTextSize(0.03);
  sb_latex->SetTextAlign(2);
  sb_latex->Draw();
  
  char chi2_text[100];
  sprintf(chi2_text,"#chi^2/ndf =%4.2f ",Chi2/NDF);
  printf("chi^2/ndf =%4.2f\n",Chi2/NDF);
  TLatex *chi2_latex = new TLatex(4.5,0.106*histo->GetMaximum(),chi2_text);
  chi2_latex->SetTextColor(2);
  chi2_latex->SetTextSize(0.03);
  chi2_latex->SetTextAlign(2);
  chi2_latex->Draw();
  printf("\n****************************************************************************\n");
  
}


