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

/* History of cvs commits:
 *
 * $Log$
 *
 */

//_________________________________________________________________________
// Class for the analysis of gamma-jet correlations 
//  Basically it seaches for a prompt photon in the Calorimeters acceptance, 
//  if so we construct a jet around the highest pt particle in the opposite 
//  side in azimuth, inside the Central Tracking System (ITS+TPC) and 
//  EMCAL acceptances (only when PHOS detects the gamma). First the leading 
//  particle and then the jet have to fullfill several conditions 
//  (energy, direction ..) to be accepted. Then the fragmentation function 
//  of this jet is constructed   
//  Class created from old AliPHOSGammaPion 
//  (see AliRoot versions previous Release 4-09)
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

#include <TFile.h>
#include <TParticle.h>
#include <TH2.h>

#include "AliAnaGammaHadron.h" 
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDCaloCluster.h"
#include "Riostream.h"
#include "AliLog.h"

ClassImp(AliAnaGammaHadron)

//____________________________________________________________________________
AliAnaGammaHadron::AliAnaGammaHadron(const char *name) : 
  AliAnaGammaDirect(name), 
  fPhiMaxCut(0.), fPhiMinCut(0.), 
  fInvMassMaxCut(0.), fInvMassMinCut(0.),
  fMinPtPion(0),
  fAngleMaxParam()
{

  // ctor
  fAngleMaxParam.Set(4) ;
  fAngleMaxParam.Reset(0.);
        
  TList * list = gDirectory->GetListOfKeys() ; 
  TIter next(list) ; 
  TH2F * h = 0 ;
  Int_t index ; 
  for (index = 0 ; index < list->GetSize()-1 ; index++) { 
    //-1 to avoid GammaPion Task
    h = dynamic_cast<TH2F*>(gDirectory->Get(list->At(index)->GetName())) ; 
    fOutputContainer->Add(h) ; 
  }

  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
  
}


//____________________________________________________________________________
AliAnaGammaHadron::AliAnaGammaHadron(const AliAnaGammaHadron & gj) : 
  AliAnaGammaDirect(gj), 
  fPhiMaxCut(gj.fPhiMaxCut), fPhiMinCut(gj.fPhiMinCut), 
  fInvMassMaxCut(gj.fInvMassMaxCut), fInvMassMinCut(gj.fInvMassMinCut),
  fMinPtPion(gj.fMinPtPion),
  fOutputContainer(0), fAngleMaxParam(gj.fAngleMaxParam)

{
  // cpy ctor
  SetName (gj.GetName()) ; 
  SetTitle(gj.GetTitle()) ; 

}

//____________________________________________________________________________
AliAnaGammaHadron::~AliAnaGammaHadron() 
{
  // Remove all pointers
  fOutputContainer->Clear() ; 
  delete fOutputContainer ;
  
  delete fhPhiCharged  ;  
  delete fhPhiNeutral   ; 
  delete fhEtaCharged  ; 
  delete fhEtaNeutral  ; 
  delete fhDeltaPhiGammaCharged  ;  
  delete fhDeltaPhiGammaNeutral   ; 
  delete fhDeltaEtaGammaCharged  ; 
  delete fhDeltaEtaGammaNeutral  ; 

  delete fhCorrelationGammaNeutral  ; 
  delete fhCorrelationGammaCharged  ; 

  delete fhAnglePairNoCut  ; 
  delete fhAnglePairAzimuthCut  ; 
  delete fhAnglePairOpeningAngleCut   ; 
  delete fhAnglePairAllCut   ;  
  delete fhInvMassPairNoCut    ; 
  delete fhInvMassPairAzimuthCut  ; 
  delete fhInvMassPairOpeningAngleCut  ; 
  delete fhInvMassPairAllCut   ;    

}



//____________________________________________________________________________
void AliAnaGammaHadron::Exec(Option_t *) 
{
  
  // Processing of one event
    
  //Get ESDs
  Long64_t entry = GetChain()->GetReadEntry() ;
  
  if (!GetESD()) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  if (GetPrintInfo()) 
    AliInfo(Form("%s ----> Processing event # %lld",  (dynamic_cast<TChain *>(GetChain()))->GetFile()->GetName(), entry)) ; 

  //CreateTLists with arrays of TParticles. Filled with particles only relevant for the analysis.

  TClonesArray * particleList = new TClonesArray("TParticle",1000); // All particles refitted in CTS and detected in EMCAL (jet)
  TClonesArray * plCTS         = new TClonesArray("TParticle",1000); // All particles refitted in Central Tracking System (ITS+TPC)
  TClonesArray * plNe          = new TClonesArray("TParticle",1000);   // All particles measured in Jet Calorimeter
  TClonesArray * plCalo     = new TClonesArray("TParticle",1000);  // All particles measured in Prompt Gamma calorimeter


  TParticle *pGamma = new TParticle(); //It will contain the kinematics of the found prompt gamma

 
  Bool_t iIsInPHOSorEMCAL = kFALSE ; //To check if Gamma was in any calorimeter
  
  //Fill lists with photons, neutral particles and charged particles
  //look for the highest energy photon in the event inside fCalorimeter
  //Fill particle lists 
  AliDebug(2, "Fill particle lists, get prompt gamma");

  //Fill particle lists 
  if(GetCalorimeter() == "PHOS")
    CreateParticleList(particleList, plCTS,plNe,plCalo); 
  else if(GetCalorimeter() == "EMCAL")
    CreateParticleList(particleList, plCTS,plCalo,plNe); 
  else
    AliError("No calorimeter selected");
 
  //Search highest energy prompt gamma in calorimeter
  GetPromptGamma(plCalo,  plCTS, pGamma, iIsInPHOSorEMCAL) ; 


  AliDebug(1, Form("Is Gamma in %s? %d",GetCalorimeter().Data(),iIsInPHOSorEMCAL));
    AliDebug(3,Form("Charged list entries %d, Neutral list entries %d, %s list entries %d",
		    plCTS->GetEntries(),plNe->GetEntries(), GetCalorimeter().Data(),plCalo->GetEntries()));
    
  //If there is any prompt photon  in fCalorimeter, 
  //search jet leading particle

  if(iIsInPHOSorEMCAL){

    if (GetPrintInfo())
      AliInfo(Form("Prompt Gamma: pt %f, phi %f, eta %f", pGamma->Pt(),pGamma->Phi(),pGamma->Eta())) ;
    
    AliDebug(2, "Make correlation");
    
    //Search correlation 
    MakeGammaChargedCorrelation(plCTS, pGamma);
    MakeGammaNeutralCorrelation(plNe, pGamma);

  }//Gamma in Calo
     
  AliDebug(2, "End of analysis, delete pointers");

  particleList->Delete() ; 
  plCTS->Delete() ;
  plNe->Delete() ;
  plCalo->Delete() ;
  pGamma->Delete();

  delete plNe ;
  delete plCalo ;
  delete plCTS ;
  delete particleList ;
  //  delete pGamma;

  PostData(0, fOutputContainer);
}    


//____________________________________________________________________________
void  AliAnaGammaHadron::MakeGammaChargedCorrelation(TClonesArray * pl, TParticle * pGamma) const 
{  
  //Search for the charged particle with highest with 
  //Phi=Phi_gamma-Pi and pT=0.1E_gamma 
  Double_t ptg  = pGamma->Pt();
  Double_t phig = pGamma->Phi();
  Double_t pt    = -100.;
  Double_t rat   = -100.; 
  Double_t phi   = -100. ;

  for(Int_t ipr = 0;ipr < pl->GetEntries() ; ipr ++ ){
    
    TParticle * particle = dynamic_cast<TParticle *>(pl->At(ipr)) ;

    pt    = particle->Pt();
    rat   = pt/ptg ;
    phi   = particle->Phi() ;
    
    AliDebug(3,Form("pt %f, phi %f, phi gamma %f. Cuts:  delta phi min %f,  max%f, pT min %f",pt,phi,phig,fPhiMinCut,fPhiMaxCut,fMinPtPion));
    
    fhEtaCharged->Fill(ptg,particle->Eta());
    fhPhiCharged->Fill(ptg,phi);
    fhDeltaEtaGammaCharged->Fill(ptg,pGamma->Eta()-particle->Eta());
    fhDeltaPhiGammaCharged->Fill(ptg,phig-phi);
    //Selection within angular and energy limits
    if(((phig-phi)> fPhiMinCut) && ((phig-phi)<fPhiMaxCut) && pt > fMinPtPion){
      AliDebug(2,Form("Selected: pt %f, phi %f",pt,phi));
      fhCorrelationGammaCharged->Fill(ptg,rat);
    } 
  }//particle loop
}

//____________________________________________________________________________
void  AliAnaGammaHadron::MakeGammaNeutralCorrelation(TClonesArray * pl, TParticle * pGamma)  
{  

  //Search for the neutral pion with highest with 
  //Phi=Phi_gamma-Pi and pT=0.1E_gamma 
  Double_t pt = -100.;
  Double_t rat = -100.; 
  Double_t phi = -100. ;
  Double_t ptg  = pGamma->Pt();
  Double_t phig = pGamma->Phi();

  TIter next(pl);
  TParticle * particle = 0;
  
  Int_t iPrimary = -1;
  TLorentzVector gammai,gammaj;
  Double_t angle = 0., e = 0., invmass = 0.;
  Int_t ksPdg = 0;
  Int_t jPrimary=-1;

  while ( (particle = (TParticle*)next()) ) {
    iPrimary++;	  
    ksPdg = particle->GetPdgCode();

    //2 gamma overlapped, found with PID
    if(ksPdg == 111){ 
      pt  = particle->Pt();
      rat = pt/ptg ;
      phi = particle->Phi() ;
      fhEtaCharged->Fill(ptg,particle->Eta());
      fhPhiCharged->Fill(ptg,phi);
      fhDeltaEtaGammaCharged->Fill(ptg,pGamma->Eta()-particle->Eta());
      fhDeltaPhiGammaCharged->Fill(ptg,phig-phi);
      //AliDebug(3,Form("pt %f, phi %f",pt,phi));
      if (GetPrintInfo())
	AliInfo(Form("pt %f, phi %f",pt,phi));
      //Selection within angular and energy limits
      if((pt> ptg)&& ((phig-phi)>fPhiMinCut)&&((phig-phi)<fPhiMaxCut)){
	fhCorrelationGammaNeutral ->Fill(ptg,rat);
	//AliDebug(2,Form("Selected: pt %f, phi %f",pt,phi));
	if (GetPrintInfo())
	  AliInfo(Form("Selected: pt %f, phi %f",pt,phi));
      }// cuts
    }// pdg = 111

    //Make invariant mass analysis
    else if(ksPdg == 22){
      //Search the photon companion in case it comes from  a Pi0 decay
      //Apply several cuts to select the good pair
      particle->Momentum(gammai);
      jPrimary=-1;
      TIter next2(pl);
      while ( (particle = (TParticle*)next2()) ) {
	jPrimary++;
	if(jPrimary>iPrimary){
	  ksPdg = particle->GetPdgCode();

	  if(ksPdg == 22){
	    particle->Momentum(gammaj);
	    //Info("GetLeadingPi0","Egammai %f, Egammaj %f", 
	    //gammai.Pt(),gammaj.Pt());
	    pt  = (gammai+gammaj).Pt();
	    phi = (gammai+gammaj).Phi();
	    if(phi < 0)
	      phi+=TMath::TwoPi();
	    rat          = pt/ptg ;
	    invmass = (gammai+gammaj).M();
	    angle      = gammaj.Angle(gammai.Vect());
	    e             = (gammai+gammaj).E();
	    fhEtaNeutral->Fill(ptg,(gammai+gammaj).Eta());
	    fhPhiNeutral->Fill(ptg,phi);
	    fhDeltaEtaGammaNeutral->Fill(ptg,pGamma->Eta()-(gammai+gammaj).Eta());
	    fhDeltaPhiGammaNeutral->Fill(ptg,phig-phi);
	    // AliDebug(3,Form("pt %f, phi %f",pt,phi));
	    if (GetPrintInfo())
	      AliInfo(Form("pt %f, phi %f",pt,phi));

	    //Fill histograms with no cuts applied.
	    fhAnglePairNoCut->Fill(e,angle);
	    fhInvMassPairNoCut->Fill(ptg,invmass);

	    //First cut on the energy and azimuth of the pair
	    if((phig-phi) > fPhiMinCut && (phig-phi) < fPhiMaxCut 
	       && pt > fMinPtPion){
	      fhAnglePairAzimuthCut     ->Fill(e,angle);
	      fhInvMassPairAzimuthCut->Fill(ptg,invmass);
	      AliDebug(3,Form("1st cut: pt %f, phi %f",pt,phi));

	      //Second cut on the aperture of the pair
	      if(IsAngleInWindow(angle,e)){
		fhAnglePairOpeningAngleCut     ->Fill(e,angle);
		fhInvMassPairOpeningAngleCut->Fill(ptg,invmass);
		 AliDebug(3,Form("2nd cut: pt %f, phi %f",pt,phi));

		//Third cut on the invariant mass of the pair
		if((invmass>fInvMassMinCut) && (invmass<fInvMassMaxCut)){ 
		  fhInvMassPairAllCut  ->Fill(ptg,invmass);
		  fhAnglePairAllCut       ->Fill(e,angle);
		  //Fill correlation histogram
		  fhCorrelationGammaNeutral ->Fill(ptg,rat);
		  //AliDebug(2,Form("Selected: pt %f, phi %f",pt,phi));
		  if (GetPrintInfo())
		    AliInfo(Form("Selected: pt %f, phi %f",pt,phi));
		}//(invmass>0.125) && (invmass<0.145)
	      }//Opening angle cut
	    }//Azimuth and pt cut.
	  }//if pdg = 22
	}//iprimary<jprimary
      }//while
    }// if pdg = 22
  }//while
}

  //____________________________________________________________________________
void AliAnaGammaHadron::Init(const Option_t * )
{
  // Initialisation of branch container 
  AliAnaGammaDirect::Init();
 
  //Initialize the parameters of the analysis.
  //fCalorimeter="PHOS";
  fAngleMaxParam.Set(4) ;
  fAngleMaxParam.AddAt(0.4,0);//={0.4,-0.25,0.025,-2e-4};
  fAngleMaxParam.AddAt(-0.25,1) ;
  fAngleMaxParam.AddAt(0.025,2) ;
  fAngleMaxParam.AddAt(-2e-4,3) ;

  //fPrintInfo           = kTRUE;
  fInvMassMaxCut  = 0.16 ;
  fInvMassMinCut  = 0.11 ;
  fPhiMaxCut      = 4.5;
  fPhiMinCut      = 1.5 ;

  fMinPtPion = 0.   ;

  //Fill particle lists when PID is ok
  // fEMCALPID = kFALSE;
  // fPHOSPID = kFALSE;

  //Initialization of histograms 

  MakeHistos() ; 

}

//__________________________________________________________________________-
Bool_t AliAnaGammaHadron::IsAngleInWindow(const Float_t angle,const Float_t e) {
  //Check if the opening angle of the candidate pairs is inside 
  //our selection windowd

  Bool_t result = kFALSE;
  Double_t mpi0 = 0.1349766;
  Double_t max =  fAngleMaxParam.At(0)*TMath::Exp(fAngleMaxParam.At(1)*e)
    +fAngleMaxParam.At(2)+fAngleMaxParam.At(3)*e;
  Double_t arg = (e*e-2*mpi0*mpi0)/(e*e);
  Double_t min = 100. ;
  if(arg>0.)
    min = TMath::ACos(arg);

  if((angle<max)&&(angle>=min))
    result = kTRUE;
 
  return result;
}

//____________________________________________________________________________
void AliAnaGammaHadron::MakeHistos()
{
  // Create histograms to be saved in output file and 
  // stores them in fOutputContainer
  
  fOutputContainer = new TObjArray(10000) ;

  //Use histograms in AliAnaGammaDirect
  TObjArray  * outputContainer =GetOutputContainer();
  for(Int_t i = 0; i < outputContainer->GetEntries(); i++ )
    fOutputContainer->Add(outputContainer->At(i)) ;
    
  fhPhiCharged  = new TH2F
    ("PhiCharged","#phi_{#pi^{#pm}}  vs p_{T #gamma}",
     120,0,120,120,0,7); 
  fhPhiCharged->SetYTitle("#phi_{#pi^{#pm}} (rad)");
  fhPhiCharged->SetXTitle("p_{T #gamma} (GeV/c)");
  fOutputContainer->Add(fhPhiCharged) ;
  
  fhPhiNeutral  = new TH2F
    ("PhiNeutral","#phi_{#pi^{0}}  vs p_{T #gamma}",
     120,0,120,120,0,7); 
  fhPhiNeutral->SetYTitle("#phi_{#pi^{0}} (rad)");
  fhPhiNeutral->SetXTitle("p_{T #gamma} (GeV/c)");
  fOutputContainer->Add(fhPhiNeutral) ;  
  
  fhEtaCharged  = new TH2F
    ("EtaCharged","#eta_{#pi^{#pm}}  vs p_{T #gamma}",
     120,0,120,120,-1,1); 
  fhEtaCharged->SetYTitle("#eta_{#pi^{#pm}} (rad)");
  fhEtaCharged->SetXTitle("p_{T #gamma} (GeV/c)");
  fOutputContainer->Add(fhEtaCharged) ;

  fhEtaNeutral  = new TH2F
    ("EtaNeutral","#eta_{#pi^{0}}  vs p_{T #gamma}",
     120,0,120,120,-1,1); 
  fhEtaNeutral->SetYTitle("#eta_{#pi^{0}} (rad)");
  fhEtaNeutral->SetXTitle("p_{T #gamma} (GeV/c)");
  fOutputContainer->Add(fhEtaNeutral) ;  

  fhDeltaPhiGammaCharged  = new TH2F
    ("DeltaPhiGammaCharged","#phi_{#gamma} - #phi_{charged #pi} vs p_{T #gamma}",
     200,0,120,200,0,6.4); 
  fhDeltaPhiGammaCharged->SetYTitle("#Delta #phi");
  fhDeltaPhiGammaCharged->SetXTitle("p_{T #gamma} (GeV/c)");
  fOutputContainer->Add(fhDeltaPhiGammaCharged) ; 
  
  fhDeltaEtaGammaCharged  = new TH2F
    ("DeltaEtaGammaCharged","#eta_{#gamma} - #eta_{#pi^{#pm}} vs p_{T #gamma}",
     200,0,120,200,-2,2); 
  fhDeltaEtaGammaCharged->SetYTitle("#Delta #eta");
  fhDeltaEtaGammaCharged->SetXTitle("p_{T #gamma} (GeV/c)");
  fOutputContainer->Add(fhDeltaEtaGammaCharged) ; 

  fhDeltaPhiGammaNeutral  = new TH2F
    ("DeltaPhiGammaNeutral","#phi_{#gamma} - #phi_{#pi^{0}} vs p_{T #gamma}",
     200,0,120,200,0,6.4); 
  fhDeltaPhiGammaNeutral->SetYTitle("#Delta #phi");
  fhDeltaPhiGammaNeutral->SetXTitle("p_{T #gamma} (GeV/c)");
  fOutputContainer->Add(fhDeltaPhiGammaNeutral) ; 
  
  fhDeltaEtaGammaNeutral  = new TH2F
    ("DeltaEtaGammaNeutral","#eta_{#gamma} - #eta_{#pi^{#pm}} vs p_{T #gamma}",
     200,0,120,200,-2,2); 
  fhDeltaEtaGammaNeutral->SetYTitle("#Delta #eta");
  fhDeltaEtaGammaNeutral->SetXTitle("p_{T #gamma} (GeV/c)");
  fOutputContainer->Add(fhDeltaEtaGammaNeutral) ; 
  
  //
  fhAnglePairAccepted  = new TH2F
    ("AnglePairAccepted",
     "Angle between #pi^{0} #gamma pair vs p_{T  #pi^{0}}, both #gamma in eta<0.7, inside window",
     200,0,50,200,0,0.2); 
  fhAnglePairAccepted->SetYTitle("Angle (rad)");
  fhAnglePairAccepted->SetXTitle("E_{ #pi^{0}} (GeV/c)");
  fOutputContainer->Add(fhAnglePairAccepted) ; 
  
  fhAnglePairNoCut  = new TH2F
    ("AnglePairNoCut",
     "Angle between all #gamma pair vs p_{T  #pi^{0}}",200,0,50,200,0,0.2); 
  fhAnglePairNoCut->SetYTitle("Angle (rad)");
  fhAnglePairNoCut->SetXTitle("E_{ #pi^{0}} (GeV/c)");
  fOutputContainer->Add(fhAnglePairNoCut) ; 
  
  fhAnglePairAzimuthCut  = new TH2F
    ("AnglePairAzimuthCut",
     "Angle between all #gamma pair that have a good phi and pt vs p_{T  #pi^{0}}",
     200,0,50,200,0,0.2); 
  fhAnglePairAzimuthCut->SetYTitle("Angle (rad)");
  fhAnglePairAzimuthCut->SetXTitle("E_{ #pi^{0}} (GeV/c)");
  fOutputContainer->Add(fhAnglePairAzimuthCut) ; 
  
    fhAnglePairOpeningAngleCut  = new TH2F
      ("AnglePairOpeningAngleCut",
       "Angle between all #gamma pair (opening angle + azimuth cut) vs p_{T  #pi^{0}}"
       ,200,0,50,200,0,0.2); 
    fhAnglePairOpeningAngleCut->SetYTitle("Angle (rad)");
    fhAnglePairOpeningAngleCut->SetXTitle("E_{ #pi^{0}} (GeV/c)");
    fOutputContainer->Add(fhAnglePairOpeningAngleCut) ;
    
    fhAnglePairAllCut  = new TH2F
      ("AnglePairAllCut",
       "Angle between all #gamma pair (opening angle + inv mass cut+azimuth) vs p_{T  #pi^{0}}"
       ,200,0,50,200,0,0.2); 
    fhAnglePairAllCut->SetYTitle("Angle (rad)");
    fhAnglePairAllCut->SetXTitle("E_{ #pi^{0}} (GeV/c)");
    fOutputContainer->Add(fhAnglePairAllCut) ; 
    
    
    //
    fhInvMassPairNoCut  = new TH2F
      ("InvMassPairNoCut","Invariant Mass of all #gamma pair vs p_{T #gamma}",
       120,0,120,360,0,0.5); 
    fhInvMassPairNoCut->SetYTitle("Invariant Mass (GeV/c^{2})");
    fhInvMassPairNoCut->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhInvMassPairNoCut) ; 
    
    fhInvMassPairAzimuthCut  = new TH2F
      ("InvMassPairAzimuthCut",
       "Invariant Mass of #gamma pair (azimuth cuts) vs p_{T #gamma}",
       120,0,120,360,0,0.5); 
    fhInvMassPairAzimuthCut->SetYTitle("Invariant Mass (GeV/c^{2})");
    fhInvMassPairAzimuthCut->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhInvMassPairAzimuthCut) ; 
    
    fhInvMassPairOpeningAngleCut  = new TH2F
      ("InvMassPairOpeningAngleCut",
       "Invariant Mass of #gamma pair (angle cut) vs p_{T #gamma}",
       120,0,120,360,0,0.5); 
    fhInvMassPairOpeningAngleCut->SetYTitle("Invariant Mass (GeV/c^{2})");
    fhInvMassPairOpeningAngleCut->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhInvMassPairOpeningAngleCut) ; 
    
    fhInvMassPairAllCut  = new TH2F
      ("InvMassPairAllCut",
       "Invariant Mass of #gamma pair (opening angle+invmass cut+azimuth) vs p_{T #gamma}",
       120,0,120,360,0,0.5); 
    fhInvMassPairAllCut->SetYTitle("Invariant Mass (GeV/c^{2})");
    fhInvMassPairAllCut->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhInvMassPairAllCut) ; 
 
    //   
    fhCorrelationGammaCharged  = 
      new TH2F("CorrelationGammaCharged","z_{#gamma #pi} = p_{T #pi^{#pm}} / p_{T #gamma}",
	       240,0.,120.,1000,0.,1.2); 
    fhCorrelationGammaCharged->SetYTitle("z_{#gamma #pi}");
    fhCorrelationGammaCharged->SetXTitle("p_{T #gamma}");
    fOutputContainer->Add(fhCorrelationGammaCharged) ;

    fhCorrelationGammaNeutral  = 
      new TH2F("CorrelationGammaNeutral","z_{#gamma #pi} = p_{T #pi^{0}} / p_{T #gamma}",
	       240,0.,120.,1000,0.,1.2); 
    fhCorrelationGammaNeutral->SetYTitle("z_{#gamma #pi}");
    fhCorrelationGammaNeutral->SetXTitle("p_{T #gamma}");
    fOutputContainer->Add(fhCorrelationGammaNeutral) ;

}

//____________________________________________________________________________
void AliAnaGammaHadron::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;

  Info("Print", "%s %s", GetName(), GetTitle() ) ;
  printf("pT Pion       >    %f\n", fMinPtPion) ; 
  printf("Phi Pion      <     %f\n", fPhiMaxCut) ; 
  printf("Phi Pion      >     %f\n", fPhiMinCut) ;
  printf("M_pair        <     %f\n", fInvMassMaxCut) ; 
  printf("M_pair        >     %f\n", fInvMassMinCut) ; 
 
} 

//__________________________________________
void AliAnaGammaHadron::Terminate(Option_t *)
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
    

}
