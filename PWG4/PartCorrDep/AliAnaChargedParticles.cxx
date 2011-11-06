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
/* $Id: $ */

//_________________________________________________________________________
//
// Class for track selection and identification (not done now)
// Tracks from the CTS are kept in the AOD.
// Few histograms produced.
//
//-- Author: Gustavo Conesa (INFN-LNF)
//_________________________________________________________________________


// --- ROOT system ---
#include "TParticle.h"
#include "TH2F.h"

//---- AliRoot system ----
#include "AliAnaChargedParticles.h"
#include "AliCaloTrackReader.h"
#include "AliAODPWG4Particle.h"
#include "AliStack.h"
#include "AliFiducialCut.h"
#include "AliVTrack.h"
#include "AliAODMCParticle.h"

ClassImp(AliAnaChargedParticles)
  
//__________________________________________________
  AliAnaChargedParticles::AliAnaChargedParticles() : 
    AliAnaPartCorrBaseClass(),
    fPdg(0), 
    fhNtracks(0),   fhPt(0),
    fhPhiNeg(0),    fhEtaNeg(0), 
    fhPhiPos(0),    fhEtaPos(0), 
    fhEtaPhiPos(0), fhEtaPhiNeg(0),
    //MC
    fhPtPion(0),    fhPhiPion(0),     fhEtaPion(0),
    fhPtProton(0),  fhPhiProton(0),   fhEtaProton(0),
    fhPtElectron(0),fhPhiElectron(0), fhEtaElectron(0),
    fhPtKaon(0),    fhPhiKaon(0),     fhEtaKaon(0),
    fhPtUnknown(0), fhPhiUnknown(0),  fhEtaUnknown(0)
{
  //Default Ctor

  //Initialize parameters
  InitParameters();

}

//_______________________________________________________
TList *  AliAnaChargedParticles::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in fOutputContainer
  
  
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("ExampleHistos") ; 
  
  Int_t nptbins  = GetHistoPtBins(); Int_t nphibins = GetHistoPhiBins(); Int_t netabins = GetHistoEtaBins();
  Float_t ptmax  = GetHistoPtMax();  Float_t phimax = GetHistoPhiMax();  Float_t etamax = GetHistoEtaMax();
  Float_t ptmin  = GetHistoPtMin();  Float_t phimin = GetHistoPhiMin();  Float_t etamin = GetHistoEtaMin();	

  fhNtracks  = new TH1F ("hNtracks","# of tracks", 1000,0,1000); 
  fhNtracks->SetXTitle("# of tracks");
  outputContainer->Add(fhNtracks);
    
  fhPt  = new TH1F ("hPt","p_T distribution", nptbins,ptmin,ptmax); 
  fhPt->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPt);
  
  fhPhiNeg  = new TH2F ("hPhiNegative","#phi of negative charges distribution",
                        nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
  fhPhiNeg->SetYTitle("#phi (rad)");
  fhPhiNeg->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPhiNeg);
  
  fhEtaNeg  = new TH2F ("hEtaNegative","#eta of negative charges distribution",
                        nptbins,ptmin,ptmax, netabins,etamin,etamax); 
  fhEtaNeg->SetYTitle("#eta ");
  fhEtaNeg->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhEtaNeg);
  
  fhPhiPos  = new TH2F ("hPhiPositive","#phi of negative charges distribution",
                        nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
  fhPhiPos->SetYTitle("#phi (rad)");
  fhPhiPos->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPhiPos);
  
  fhEtaPos  = new TH2F ("hEtaPositive","#eta of negative charges distribution",
                        nptbins,ptmin,ptmax, netabins,etamin,etamax); 
  fhEtaPos->SetYTitle("#eta ");
  fhEtaPos->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhEtaPos);
  
  fhEtaPhiPos  = new TH2F ("hPtEtaPhiPositive","pt/eta/phi of positive charge",netabins,etamin,etamax, nphibins,phimin,phimax); 
  fhEtaPhiPos->SetXTitle("p_{T}^{h^{+}} (GeV/c)");
  fhEtaPhiPos->SetYTitle("#eta ");
  fhEtaPhiPos->SetZTitle("#phi (rad)");  
  outputContainer->Add(fhEtaPhiPos);
  
  fhEtaPhiNeg  = new TH2F ("hEtaPhiNegative","eta vs phi of negative charge",netabins,etamin,etamax, nphibins,phimin,phimax); 
  fhEtaPhiNeg->SetXTitle("p_{T}^{h^{-}} (GeV/c)");
  fhEtaPhiNeg->SetYTitle("#eta ");
  fhEtaPhiNeg->SetZTitle("#phi (rad)");  
  outputContainer->Add(fhEtaPhiNeg);
  
  if(IsDataMC()){
    
    fhPtPion  = new TH1F ("hPtMCPion","p_T distribution from #pi", nptbins,ptmin,ptmax); 
    fhPtPion->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtPion);
    
    fhPhiPion  = new TH2F ("hPhiMCPion","#phi distribution from #pi",nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
    fhPhiPion->SetXTitle("#phi (rad)");
    outputContainer->Add(fhPhiPion);
    
    fhEtaPion  = new TH2F ("hEtaMCPion","#eta distribution from #pi",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
    fhEtaPion->SetXTitle("#eta ");
    outputContainer->Add(fhEtaPion);
    
    fhPtProton  = new TH1F ("hPtMCProton","p_T distribution from proton", nptbins,ptmin,ptmax); 
    fhPtProton->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtProton);
    
    fhPhiProton  = new TH2F ("hPhiMCProton","#phi distribution from proton",nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
    fhPhiProton->SetXTitle("#phi (rad)");
    outputContainer->Add(fhPhiProton);
    
    fhEtaProton  = new TH2F ("hEtaMCProton","#eta distribution from proton",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
    fhEtaProton->SetXTitle("#eta ");
    outputContainer->Add(fhEtaProton);
    
    fhPtKaon  = new TH1F ("hPtMCKaon","p_T distribution from kaon", nptbins,ptmin,ptmax); 
    fhPtKaon->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtKaon);
    
    fhPhiKaon  = new TH2F ("hPhiMCKaon","#phi distribution from kaon",nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
    fhPhiKaon->SetXTitle("#phi (rad)");
    outputContainer->Add(fhPhiKaon);
    
    fhEtaKaon  = new TH2F ("hEtaMCKaon","#eta distribution from kaon",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
    fhEtaKaon->SetXTitle("#eta ");
    outputContainer->Add(fhEtaKaon);
    
    fhPtElectron  = new TH1F ("hPtMCElectron","p_T distribution from electron", nptbins,ptmin,ptmax); 
    fhPtElectron->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtElectron);
    
    fhPhiElectron  = new TH2F ("hPhiMCElectron","#phi distribution from electron",nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
    fhPhiElectron->SetXTitle("#phi (rad)");
    outputContainer->Add(fhPhiElectron);
    
    fhEtaElectron  = new TH2F ("hEtaMCElectron","#eta distribution from electron",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
    fhEtaElectron->SetXTitle("#eta ");
    outputContainer->Add(fhEtaElectron);
    
    fhPtUnknown  = new TH1F ("hPtMCUnknown","p_T distribution from unknown", nptbins,ptmin,ptmax); 
    fhPtUnknown->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtUnknown);
    
    fhPhiUnknown  = new TH2F ("hPhiMCUnknown","#phi distribution from unknown",nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
    fhPhiUnknown->SetXTitle("#phi (rad)");
    outputContainer->Add(fhPhiUnknown);
    
    fhEtaUnknown  = new TH2F ("hEtaMCUnknown","#eta distribution from unknown",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
    fhEtaUnknown->SetXTitle("#eta ");
    outputContainer->Add(fhEtaUnknown);
    
  }
  
  return outputContainer;

}

//___________________________________________
void AliAnaChargedParticles::InitParameters()
{ 
  //Initialize the parameters of the analysis.
  SetOutputAODClassName("AliAODPWG4Particle");
  SetOutputAODName("PWG4Particle");

  AddToHistogramsName("AnaCharged_");

  fPdg = -1; //Select all tracks 
  
}

//____________________________________________________________
void AliAnaChargedParticles::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");	
	
  printf("Min Pt = %3.2f\n", GetMinPt());
  printf("Max Pt = %3.2f\n", GetMaxPt());
  printf("Select clusters with pdg %d \n",fPdg);
  
} 

//_________________________________
void AliAnaChargedParticles::Init()
{  
  //Init
  //Do some checks
  if(!GetReader()->IsCTSSwitchedOn()){
    printf("AliAnaChargedParticles::Init() - STOP!: You want to use CTS tracks in analysis but not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  
}

//_________________________________________________
void  AliAnaChargedParticles::MakeAnalysisFillAOD() 
{
  //Do analysis and fill aods
  if(!GetCTSTracks() || GetCTSTracks()->GetEntriesFast() == 0) return ;
  Int_t ntracks = GetCTSTracks()->GetEntriesFast();
  Double_t vert[3] = {0,0,0}; //vertex ;
  //Some prints
  if(GetDebug() > 0)
    printf("AliAnaChargedParticles::MakeAnalysisFillAOD() - In CTS aod entries %d\n", ntracks);
  
  //Fill AODParticle with CTS aods
  TVector3 p3;
  Int_t evtIndex = 0;
  for(Int_t i = 0; i < ntracks; i++){
    
    AliVTrack * track =  (AliVTrack*) (GetCTSTracks()->At(i));
    
    //Fill AODParticle after some selection       
    Double_t mom[3] = {track->Px(),track->Py(),track->Pz()};
    p3.SetXYZ(mom[0],mom[1],mom[2]);
    
    //Acceptance selection
    Bool_t in =  GetFiducialCut()->IsInFiducialCut(mom,"CTS") ;
    if(GetDebug() > 1) 
      printf("AliAnaChargedParticles::MakeAnalysisFillAOD() - Track pt %2.2f, phi %2.2f, in fiducial cut %d\n",p3.Pt(), p3.Phi(), in);

    if(p3.Pt() > GetMinPt() && in) {
      //Keep only particles identified with fPdg
      //Selection not done for the moment
      //Should be done here.
      
      // Mixed event
      if (GetMixedEvent()){
        evtIndex = GetMixedEvent()->EventIndex(track->GetID()) ;
      } 
      GetVertex(vert,evtIndex); 
      if(TMath::Abs(vert[2])> GetZvertexCut()) return; 
      
      AliAODPWG4Particle tr = AliAODPWG4Particle(mom[0],mom[1],mom[2],0);
      tr.SetDetector("CTS");
      tr.SetLabel(track->GetLabel());
      tr.SetTrackLabel(track->GetID(),-1);
      tr.SetChargedBit(track->Charge()>0);
		
      AddAODParticle(tr);
      
    }//selection
  }//loop
  
  if(GetDebug() > 0) 	
    printf("AliAnaChargedParticles::MakeAnalysisFillAOD() - Final aod branch entries %d\n", GetOutputAODBranch()->GetEntriesFast());   
} 

//__________________________________________________________________
void  AliAnaChargedParticles::MakeAnalysisFillHistograms() 
{
  //Do analysis and fill histograms
  
  //Loop on stored AODParticles
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  
  fhNtracks->Fill(GetReader()->GetTrackMultiplicity()) ;
  
  if(GetDebug() > 0) 
    printf("AliAnaChargedParticles::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
  
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* tr =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
        
    fhPt->Fill(tr->Pt());
    
    if(tr->GetChargedBit()){
      fhPhiPos   ->Fill(tr->Pt(), tr->Phi());
      fhEtaPhiPos->Fill(tr->Eta(),tr->Phi());
    }
    else{
      fhEtaNeg   ->Fill(tr->Pt(), tr->Eta());
      fhEtaPhiNeg->Fill(tr->Eta(),tr->Phi());
    }
    
    if(IsDataMC()){
      //Play with the MC stack if available		
      Int_t mompdg = -1;
      Int_t label  = TMath::Abs(tr->GetLabel());
      if(GetReader()->ReadStack()){
        TParticle * mom = GetMCStack()->Particle(label);
        mompdg =TMath::Abs(mom->GetPdgCode());
      }
      else if(GetReader()->ReadAODMCParticles()){
        AliAODMCParticle * aodmom = 0;
        //Get the list of MC particles
        aodmom = (AliAODMCParticle*) (GetReader()->GetAODMCParticles(tr->GetInputFileIndex()))->At(label);
        mompdg =TMath::Abs(aodmom->GetPdgCode());
      }
      
      if(mompdg==211){
        fhPtPion ->Fill(tr->Pt());
        fhPhiPion->Fill(tr->Pt(), tr->Phi());
        fhEtaPion->Fill(tr->Pt(), tr->Eta());
      }
      else if(mompdg==2212){
        fhPtProton ->Fill(tr->Pt());
        fhPhiProton->Fill(tr->Pt(), tr->Phi());
        fhEtaProton->Fill(tr->Pt(), tr->Eta());
      }
      else if(mompdg==321){
        fhPtKaon ->Fill(tr->Pt());
        fhPhiKaon->Fill(tr->Pt(), tr->Phi());
        fhEtaKaon->Fill(tr->Pt(), tr->Eta());
      }
      else if(mompdg==11){
        fhPtElectron ->Fill(tr->Pt());
        fhPhiElectron->Fill(tr->Pt(), tr->Phi());
        fhEtaElectron->Fill(tr->Pt(), tr->Eta());
      }
      else {
        //printf("unknown pdg %d\n",mompdg);
        fhPtUnknown ->Fill(tr->Pt());
        fhPhiUnknown->Fill(tr->Pt(), tr->Phi());
        fhEtaUnknown->Fill(tr->Pt(), tr->Eta());
      }
      
    }//Work with stack also
    
  }// aod branch loop
  
}
