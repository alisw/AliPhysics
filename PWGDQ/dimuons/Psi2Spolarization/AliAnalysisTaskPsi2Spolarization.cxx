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

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TLatex.h>
#include <TList.h>
#include <TStyle.h>

#include "AliAnalysisTaskMuonDistributions.h"
#include "AliAnalysisTaskSE.h"

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
#include "AliESDInputHandler.h"
#include "AliMuonTrackCuts.h"
#include "AliAODHeader.h"
#include "AliAODHandler.h"
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskPsi2Spolarization.h"
#include "AliAnalysisMuonUtility.h"
class AliAnalysisTaskPsi2Spolarization;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskPsi2Spolarization) // classimp: necessary for root
Float_t invm,invpt,invy,zvertexf, costhetacs, costhetahe, phihe, phics;
Int_t ntracks;
Int_t ntracklet;


void AliAnalysisTaskPsi2Spolarization::NotifyRun()
{
  fMuonTrackCuts->SetRun(fInputHandler);
}


AliAnalysisTaskPsi2Spolarization::AliAnalysisTaskPsi2Spolarization() : AliAnalysisTaskSE(), 
						 fAOD(0), fOutputList(0), fVertexZ(0), fHistInvMass(0), fHistPtDimu(0), fHistRap(0), fHistCosThetaCSDimu(0),fHistPhiHEDimu(0), fHistCosThetaHEDimu(0),fHistPhiCSDimu(0), tree(0), tree1(0), fMuonTrackCuts(new AliMuonTrackCuts("stdMuonCuts","stdMuonCuts"))
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskPsi2Spolarization::AliAnalysisTaskPsi2Spolarization(const char* name) : AliAnalysisTaskSE(name),
								 fAOD(0), fOutputList(0), fVertexZ(0), fHistInvMass(0), fHistPtDimu(0), fHistRap(0), fHistCosThetaCSDimu(0), fHistCosThetaHEDimu(0),tree(0),fHistPhiCSDimu(0),fHistPhiHEDimu(0), tree1(0), fMuonTrackCuts(new AliMuonTrackCuts("stdMuonCuts","stdMuonCuts"))
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
    
    fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuMatchHpt |AliMuonTrackCuts::kMuPdca);
    fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
}
//_____________________________________________________________________________
AliAnalysisTaskPsi2Spolarization::~AliAnalysisTaskPsi2Spolarization()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskPsi2Spolarization::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file

    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    
    fVertexZ = new TH1F("fVertexZ", "Event Counter",160,-40,40);
    
    fHistPtDimu  = new TH1F("fHistPtDimu","PtDimu;p_{T} (GeV/c)",10000,0.0,20.0);
    fHistPtDimu->GetXaxis()->SetTitle("P (GeV/c)");
    fHistPtDimu->GetYaxis()->SetTitle("dN/dP (c/GeV)");
    fHistPtDimu->SetMarkerStyle(kFullCircle);
    
    fHistInvMass = new TH1F("fHistInvMass", "M_{#mu#mu} dist",10000,0.0,10.0);
    fHistInvMass->GetXaxis()->SetTitle("M_{#mu#mu} (GeV/c^{2})");
    fHistInvMass->GetYaxis()->SetTitle("dN/dM_{#mu#mu} (GeV/c^{2})^{-1}");
    fHistInvMass->SetMarkerStyle(kFullCircle);
    
    fHistRap = new TH1F("fHistRap", "fHistRap", 1000, -4.0, -2.5);
    fHistCosThetaCSDimu  = new TH1F("fHistCosThetaCSDimu","CosThetaCSDimu;cos#theta_{CS}",100,-1.,1.);
    fHistCosThetaHEDimu  = new TH1F("fHistCosThetaHEDimu","CosThetaHEDimu;cos#theta_{HE}",100,-1.,1.);
    fHistPhiCSDimu  = new TH1F("fHistPhiCSDimu","PhiCSDimu;#phi_{CS}",100,-5.,5.);
    fHistPhiHEDimu  = new TH1F("fHistPhiHEDimu","PhiHEDimu;#phi_{HE}",100,-5.,5.);
    
    
    TTree *tree= new TTree("tree","Invmumu");
    tree->Branch("invm",&invm,"invm/F");
    tree->Branch("invpt",&invpt,"invpt/F");
    tree->Branch("invy",&invy,"invy/F");
    tree->Branch("zvertexf",&zvertexf,"zvertexf/F");
    tree->Branch("costhetacs",&costhetacs,"costhetacs/F");
    tree->Branch("costhetahe",&costhetahe,"costhetahe/F");
    tree->Branch("phics",&phics,"phics/F");
    tree->Branch("phihe",&phihe,"phihe/F");
    tree->Branch("ntracks",&ntracks,"ntracks/I");


       
    TTree *tree1= new TTree("tree1","trackletinfo");
    tree1->Branch("ntracklet",&ntracklet,"ntracklet/I");
    
    fOutputList->Add(fVertexZ);
    fOutputList->Add(fHistPtDimu);
    fOutputList->Add(fHistInvMass);
    fOutputList->Add(fHistRap);
    fOutputList->Add(fHistCosThetaCSDimu);
    fOutputList->Add(fHistCosThetaHEDimu);
    fOutputList->Add(fHistPhiCSDimu);
    fOutputList->Add(fHistPhiHEDimu);
    fOutputList->Add(tree);
    fOutputList->Add(tree1);
    fOutputList->ls();
    
    
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskPsi2Spolarization::UserExec(Option_t *)
{
    
 

  UInt_t fSelectMask= fInputHandler->IsEventSelected();
  Bool_t isINT7selected = fSelectMask& AliVEvent::kINT7;

    if (!isINT7selected) return;
    
  AliVVertex* vertex = NULL;
  AliAODVertex* vertexSPD = NULL;
    
  Bool_t physSel = kTRUE;
    
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  Double_t px1,px2,py1,py2,pz1,pz2,E1,invM,E2,M1,M2;
  
  vertexSPD = const_cast<AliAODVertex *>(fAOD->GetPrimaryVertexSPD());
  double zVertex = vertexSPD->GetZ();
    
  if(!fAOD) return;
    
  TString trigger = fAOD->GetFiredTriggerClasses();
  if (trigger.Contains("CMUL7-B-NOPF-MUFAST"))
    {
      Int_t iTracks(fAOD->GetNumberOfTracks());                                    // see how many tracks there are in the event

    
    for(Int_t i(0); i < iTracks; i++)
    {                                                                            // loop over all these tracks

        AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
	if (!fMuonTrackCuts->IsSelected(track)) continue;
        
	    if(track->Charge()==1)
	      {
		Float_t M1 = track->M();
		Float_t pxr1 = track->Px();
		Float_t pyr1 = track->Py();
		Float_t pzr1 = track->Pz();
		Float_t er1 = track->E();
		Float_t eta1 = track->Eta();
		Float_t charger1 = track->Charge();
		Float_t matchmu1 = AliAnalysisMuonUtility::GetMatchTrigger(track);
		Float_t P=pxr1*pxr1+pyr1*pyr1;
		Float_t Pt1=TMath::Sqrt(pxr1*pxr1+pyr1*pyr1);
		Float_t rap1=Rap(er1,pzr1);
		Float_t Rabs1=AliAnalysisMuonUtility::GetRabs(track);
        
   		 
        for(Int_t j(0); j < iTracks; j++)
        {                                                                                // loop ove rall these tracks
            AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(j));         // get a track (type AliAODTrack) from the event
                                                                                         // if(!track2 || !track2->TestFilterBit(1)) continue;
		if(track2->Charge()==-1)
		  {
		    Float_t M2 = track2->M();
		    Float_t pxr2 = track2->Px();
		    Float_t pyr2 = track2->Py();
		    Float_t pzr2 = track2->Pz();
		    Float_t er2 = track2->E();
		    Float_t eta2 = track2->Eta();
		    Float_t charger2 = track2->Charge();
		    Float_t matchmu2 = AliAnalysisMuonUtility::GetMatchTrigger(track2);
		    Float_t Pt2=TMath::Sqrt(pxr2*pxr2+pyr2*pyr2);
		    Float_t rap2=Rap(er2,pzr2);
		    Float_t Rabs2=AliAnalysisMuonUtility::GetRabs(track2);

		     if(Rabs1 <= 17.6 || Rabs1 >= 89.5) continue;
		     if(Rabs2 <= 17.6 || Rabs2 >= 89.5) continue;
		     if(matchmu1<1 || matchmu2<1) continue;
		     if(eta1 > -2.5 || eta1 < -4.0) continue;
		     if(eta2 > -2.5 || eta2 < -4.0) continue;
		        if( TMath::Abs( zVertex ) > 10 ) continue;

            Float_t ptdimu = TMath::Sqrt((pxr1+pxr2)*(pxr1+pxr2)+(pyr1+pyr2)*(pyr1+pyr2));
            Float_t raprec = Rap((er1+er2),(pzr1+pzr2));
            Float_t imassrec = TMath::Sqrt((er1+er2)*(er1+er2)-((pxr1+pxr2)*(pxr1+pxr2)+(pyr1+pyr2)*(pyr1+pyr2)+(pzr1+pzr2)*(pzr1+pzr2)));
            Double_t costhetaCSdimu = CostCS(track, track2);
            Double_t costhetaHEdimu = CostHE(track, track2);
            Double_t phiCSdimu = PhiCS(track, track2);
            Double_t phiHEdimu = PhiHE(track, track2);
              
	     if(raprec> -2.5 || raprec < -4.0) continue;

            fVertexZ->Fill(zVertex);
            fHistPtDimu->Fill(ptdimu);
            fHistRap->Fill(raprec);
            fHistInvMass->Fill(imassrec);
            fHistCosThetaCSDimu->Fill(costhetaCSdimu);
            fHistCosThetaHEDimu->Fill(costhetaHEdimu);
            fHistPhiCSDimu->Fill(phiCSdimu);
            fHistPhiHEDimu->Fill(phiHEdimu);
              
            invm=imassrec;
            invpt=ptdimu;
            invy=raprec;
            zvertexf=zVertex;
            costhetacs = costhetaCSdimu;
            costhetahe = costhetaHEdimu;
            phics = phiCSdimu;
            phihe = phiHEdimu;

	    


              
            ((TTree*)(fOutputList->FindObject("tree")))->Fill();

            
		  }
	      }
	    }
    }

    AliAODTracklets* mult = fAOD->GetTracklets();
    Int_t nTracklets = mult->GetNumberOfTracklets();
    ntracklet=nTracklets;

    ((TTree*)(fOutputList->FindObject("tree1")))->Fill();
    
    } // trigger loop
    
    PostData(1, fOutputList);                           


}


//_____________________________________________________________________________

Float_t AliAnalysisTaskPsi2Spolarization::Rap(Float_t e, Float_t pz) const
{
// calculate rapidity
     Float_t rap;
         if(TMath::Abs(e-pz)>1e-7){
         	rap = 0.5*TMath::Log((e+pz)/(e-pz));
       		return rap;
	    }
	        else{
     	       	rap = -200;
     		return rap;
  	    }
}
//
//_____________________________________________________________________________
//________________________________________________________________________
Double_t AliAnalysisTaskPsi2Spolarization::CostHE(AliAODTrack* track, AliAODTrack* track2)
{
    Double_t EBeam = 6500;
    Double_t mp = 0.93827231;
    Double_t pbeam = TMath::Sqrt(EBeam*EBeam - mp*mp);
    Double_t pla10 = track -> Px();
    Double_t pla11 = track -> Py();
    Double_t pla12 = track -> Pz();
    Double_t e1 = track -> E();
    Double_t mu1Charge = track -> Charge();
    Double_t pla20 = track2 -> Px();
    Double_t pla21 = track2 -> Py();
    Double_t pla22 = track2 -> Pz();
    Double_t e2 = track2 -> E();
    if(pla10==0 && pla11==0 && pla12==0 && e1==0 && mu1Charge==0 && pla20==0 && pla21==0 && pla22==0 && e2==0.){return -666.;}

    // Fill the Lorentz vector for projectile and target
    // For the moment we consider no crossing angle
    // Projectile runs towards the MUON arm
    TLorentzVector pProjLab(0.,0.,-pbeam,EBeam); // projectile
    TLorentzVector pTargLab(0.,0., pbeam,EBeam); // target
    //
    // --- Get the muons parameters in the LAB frame
    //
    TLorentzVector pMu1Lab(pla10,pla11,pla12,e1);
    TLorentzVector pMu2Lab(pla20,pla21,pla22,e2);
    //
    // --- Obtain the dimuon parameters in the LAB frame
    //
    TLorentzVector pDimuLab = pMu1Lab + pMu2Lab;
    //
    // --- Translate the dimuon parameters in the dimuon rest frame
    //
    TVector3 beta = (-1./pDimuLab.E())*pDimuLab.Vect();
    TLorentzVector pMu1Dimu = pMu1Lab;
    TLorentzVector pMu2Dimu = pMu2Lab;
    pMu1Dimu.Boost(beta);
    pMu2Dimu.Boost(beta);

    TVector3 zaxis;
    zaxis=(pDimuLab.Vect()).Unit();
    //
    // --- Calculation of the polarization angle (Helicity)
    // (angle between mu+ and the z axis defined above)
    //
    Double_t cost;
    if(mu1Charge > 0) {
      cost = zaxis.Dot((pMu1Dimu.Vect()).Unit());
    } else {
      cost = zaxis.Dot((pMu2Dimu.Vect()).Unit());
    }
    return cost;
  }
//________________________________________________________________________
Double_t AliAnalysisTaskPsi2Spolarization::CostCS(AliAODTrack* track, AliAODTrack* track2){
  Double_t EBeam = 6500.;
  Double_t mp = 0.93827231;
  Double_t pbeam = TMath::Sqrt(EBeam*EBeam - mp*mp);
  Double_t pla10 = track -> Px();
  Double_t pla11 = track -> Py();
  Double_t pla12 = track -> Pz();
  Double_t e1 = track -> E();
  Double_t mu1Charge = track -> Charge();
  Double_t pla20 = track2 -> Px();
  Double_t pla21 = track2 -> Py();
  Double_t pla22 = track2 -> Pz();
  Double_t e2 = track2 -> E();
  Double_t mu2Charge = track2 -> Charge();
  if(pla10==0 && pla11==0 && pla12==0 && e1==0 && mu1Charge==0 && pla20==0 && pla21==0 && pla22==0 && e2==0.){return -666.;}

  // Fill the Lorentz vector for projectile and target
  // For the moment we do not consider the crossing angle
  // Projectile runs towards the MUON arm
  TLorentzVector pProjCM(0.,0.,-pbeam,EBeam); // projectile
  TLorentzVector pTargCM(0.,0., pbeam,EBeam); // target
  //
  // --- Get the muons parameters in the CM frame
  //
  TLorentzVector pMu1CM(pla10,pla11,pla12,e1);
  TLorentzVector pMu2CM(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  TLorentzVector pDimuCM = pMu1CM + pMu2CM;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta = (-1./pDimuCM.E())*pDimuCM.Vect();
  TLorentzVector pMu1Dimu = pMu1CM;
  TLorentzVector pMu2Dimu = pMu2CM;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the CS angle
  //
  TVector3 zaxisCS = (((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  //
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  //
  Double_t cost;
  if(mu1Charge > 0) {
    cost = zaxisCS.Dot((pMu1Dimu.Vect()).Unit());
    // Theta CS is not properly defined for Like-Sign muons
    if(mu2Charge > 0 && cost<0) cost = -cost;
  } else {
    // Theta CS is not properly defined for Like-Sign muons
    cost = zaxisCS.Dot((pMu2Dimu.Vect()).Unit());
    if(mu2Charge < 0 && cost<0) cost = -cost;
  }
  return cost;
}

Double_t AliAnalysisTaskPsi2Spolarization::PhiCS(AliAODTrack* track, AliAODTrack* track2){
  // Cosinus of the Collins-Soper polar decay angle
  Double_t EBeam = 6500.;
  if(EBeam <= 0){
    printf("Can not compute phiCS with EBeam=%f\n",EBeam);
    return -999999999;
  }
  Double_t mp = 0.93827231;
  Double_t pbeam = TMath::Sqrt(EBeam*EBeam - mp*mp);
  Double_t pla10 = track -> Px();
  Double_t pla11 = track->Py();
  Double_t pla12 = track->Pz();
  Double_t e1 = track->E();
  Double_t mu1Charge = track -> Charge();
  Double_t pla20 = track2 -> Px();
  Double_t pla21 = track2 -> Py();
  Double_t pla22 = track2 -> Pz();
  Double_t e2 = track2 -> E();
  //Double_t mu2Charge=Mu1->Charge();
  if(pla10==0 && pla11==0 && pla12==0 && e1==0 && mu1Charge==0 && pla20==0 && pla21==0 && pla22==0 && e2==0.){return -666.;}

  // Fill the Lorentz vector for projectile and target
  // For the moment we do not consider the crossing angle
  // Projectile runs towards the MUON arm
  TLorentzVector pProjCM(0.,0.,-pbeam,EBeam); // projectile
  TLorentzVector pTargCM(0.,0., pbeam,EBeam); // target
  //
  // --- Get the muons parameters in the CM frame
  //
  TLorentzVector pMu1CM(pla10,pla11,pla12,e1);
  TLorentzVector pMu2CM(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  TLorentzVector pDimuCM = pMu1CM + pMu2CM;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta = (-1./pDimuCM.E())*pDimuCM.Vect();
  TLorentzVector pMu1Dimu = pMu1CM;
  TLorentzVector pMu2Dimu = pMu2CM;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the CS angle
  //
  TVector3 zaxisCS = (((pProjDimu.Vect()).Unit()) - ((pTargDimu.Vect()).Unit())).Unit();
  //
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  //
   TVector3 yaxisCS = (((pProjDimu.Vect()).Unit()).Cross((pTargDimu.Vect()).Unit())).Unit();
   TVector3 xaxisCS = (yaxisCS.Cross(zaxisCS)).Unit();

   Double_t phi;
   if(mu1Charge>0) phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxisCS),((pMu1Dimu.Vect()).Dot(xaxisCS)));
   else phi = TMath::ATan2((pMu2Dimu.Vect()).Dot(yaxisCS),((pMu2Dimu.Vect()).Dot(xaxisCS)));

   return phi;
}

Double_t AliAnalysisTaskPsi2Spolarization::PhiHE(AliAODTrack* track, AliAODTrack* track2){
  // Calculation the Helicity aimuthal angle (adapted from code by R. Arnaldi)
  Double_t EBeam = 6500.;
  if(EBeam <= 0){
    printf("Can not compute phiHE with EBeam=%f\n",EBeam);
    return -999999999;
  }
  Double_t mp = 0.93827231;
  Double_t pbeam = TMath::Sqrt(EBeam*EBeam - mp*mp);
  Double_t pla10 = track -> Px();
  Double_t pla11 = track -> Py();
  Double_t pla12 = track -> Pz();
  Double_t e1 = track -> E();
  Double_t mu1Charge = track -> Charge();
  Double_t pla20 = track2 -> Px();
  Double_t pla21 = track2 -> Py();
  Double_t pla22 = track2 -> Pz();
  Double_t e2 = track2 -> E();
  //Double_t mu2Charge=Mu1->Charge();
  if(pla10==0 && pla11==0 && pla12==0 && e1==0 && mu1Charge==0 && pla20==0 && pla21==0 && pla22==0 && e2==0.){return -666.;}

  // Fill the Lorentz vector for projectile and target
  // For the moment we consider no crossing angle
  // Projectile runs towards the MUON arm
  TLorentzVector pProjCM(0.,0.,-pbeam,EBeam); // projectile
  TLorentzVector pTargCM(0.,0., pbeam,EBeam); // target
  //
  // --- Get the muons parameters in the CM frame
  //
  TLorentzVector pMu1CM(pla10,pla11,pla12,e1);
  TLorentzVector pMu2CM(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  TLorentzVector pDimuCM = pMu1CM + pMu2CM;
  //
  // --- Translate the muon parameters in the dimuon rest frame
  //
  TVector3 zaxis = (pDimuCM.Vect()).Unit();
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta = (-1./pDimuCM.E())*pDimuCM.Vect();
  TLorentzVector pMu1Dimu = pMu1CM;
  TLorentzVector pMu2Dimu = pMu2CM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);

  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);

  TVector3 yaxis = ((pProjDimu.Vect()).Cross(pTargDimu.Vect())).Unit();
  TVector3 xaxis = (yaxis.Cross(zaxis)).Unit();
  //
  // --- Calculation of the azimuthal angle (Helicity)
  //
   Double_t phi;
   if(mu1Charge>0) phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxis),(pMu1Dimu.Vect()).Dot(xaxis));
   else phi = TMath::ATan2((pMu2Dimu.Vect()).Dot(yaxis),(pMu2Dimu.Vect()).Dot(xaxis));

   return phi;
}

void AliAnalysisTaskPsi2Spolarization::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
