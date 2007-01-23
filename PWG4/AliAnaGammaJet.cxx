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
//  Class created from old AliPHOSGammaJet 
//  (see AliRoot versions previous Release 4-09)
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

#include <TFile.h>
#include <TParticle.h>
#include <TH2.h>

#include "AliAnaGammaJet.h" 
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDCaloCluster.h"
#include "Riostream.h"
#include "AliLog.h"

ClassImp(AliAnaGammaJet)

//____________________________________________________________________________
AliAnaGammaJet::AliAnaGammaJet(const char *name) : 
  AliAnaGammaDirect(name), 
  fSeveralConeAndPtCuts(0),
  fPbPb(kFALSE), 
  fJetsOnlyInCTS(0),
  fEtaEMCALCut(0.),fPhiMaxCut(0.),
  fPhiMinCut(0.), 
  fInvMassMaxCut(0.), fInvMassMinCut(0.),
  fJetCTSRatioMaxCut(0.),
  fJetCTSRatioMinCut(0.), fJetRatioMaxCut(0.),
  fJetRatioMinCut(0.), fNCone(0),
  fNPt(0), fCone(0), fPtThreshold(0),
  fPtJetSelectionCut(0.0),
  fAngleMaxParam(), fSelect(0)
{

  // ctor
  fAngleMaxParam.Set(4) ;
  fAngleMaxParam.Reset(0.);

  for(Int_t i = 0; i<10; i++){
    fCones[i]         = 0.0 ;
    fNameCones[i]     = ""  ;
    fPtThres[i]      = 0.0 ;
    fNamePtThres[i]  = ""  ;
    if( i < 6 ){
      fJetXMin1[i]     = 0.0 ;
      fJetXMin2[i]     = 0.0 ;
      fJetXMax1[i]     = 0.0 ;
      fJetXMax2[i]     = 0.0 ;
      fBkgMean[i]      = 0.0 ;
      fBkgRMS[i]       = 0.0 ;
      if( i < 2 ){
	fJetE1[i]        = 0.0 ;
	fJetE2[i]        = 0.0 ;
	fJetSigma1[i]    = 0.0 ;
	fJetSigma2[i]    = 0.0 ;
	fPhiEMCALCut[i]  = 0.0 ;
      }
    }
  }
        
  TList * list = gDirectory->GetListOfKeys() ; 
  TIter next(list) ; 
  TH2F * h = 0 ;
  Int_t index ; 
  for (index = 0 ; index < list->GetSize()-1 ; index++) { 
    //-1 to avoid GammaJet Task
    h = dynamic_cast<TH2F*>(gDirectory->Get(list->At(index)->GetName())) ; 
    fOutputContainer->Add(h) ; 
  }

  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
  
}


//____________________________________________________________________________
AliAnaGammaJet::AliAnaGammaJet(const AliAnaGammaJet & gj) : 
  AliAnaGammaDirect(gj), 
  fSeveralConeAndPtCuts(gj.fSeveralConeAndPtCuts), 
  fPbPb(gj.fPbPb), fJetsOnlyInCTS(gj.fJetsOnlyInCTS),
  fEtaEMCALCut(gj.fEtaEMCALCut),
  fPhiMaxCut(gj.fPhiMaxCut), fPhiMinCut(gj.fPhiMinCut), 
  fInvMassMaxCut(gj.fInvMassMaxCut), fInvMassMinCut(gj.fInvMassMinCut),
  fRatioMinCut(gj.fRatioMinCut), 
  fJetCTSRatioMaxCut(gj.fJetCTSRatioMaxCut),
  fJetCTSRatioMinCut(gj.fJetCTSRatioMinCut), fJetRatioMaxCut(gj.fJetRatioMaxCut),
  fJetRatioMinCut(gj.fJetRatioMinCut),  fNCone(gj.fNCone),
  fNPt(gj.fNPt), fCone(gj.fCone), fPtThreshold(gj.fPtThreshold),
  fPtJetSelectionCut(gj.fPtJetSelectionCut),
  fOutputContainer(0), fAngleMaxParam(gj.fAngleMaxParam), 
  fSelect(gj.fSelect)
{
  // cpy ctor
  SetName (gj.GetName()) ; 
  SetTitle(gj.GetTitle()) ; 

  for(Int_t i = 0; i<10; i++){
    fCones[i]        = gj.fCones[i] ;
    fNameCones[i]    = gj.fNameCones[i] ;
    fPtThres[i]      = gj.fPtThres[i] ;
    fNamePtThres[i]  = gj.fNamePtThres[i] ;
    if( i < 6 ){
      fJetXMin1[i]       = gj.fJetXMin1[i] ;
      fJetXMin2[i]       = gj.fJetXMin2[i] ;
      fJetXMax1[i]       = gj.fJetXMax1[i] ;
      fJetXMax2[i]       = gj.fJetXMax2[i] ;
      fBkgMean[i]        = gj.fBkgMean[i] ;
      fBkgRMS[i]         = gj.fBkgRMS[i] ;
      if( i < 2 ){
	fJetE1[i]        = gj.fJetE1[i] ;
	fJetE2[i]        = gj.fJetE2[i] ;
	fJetSigma1[i]    = gj.fJetSigma1[i] ;
	fJetSigma2[i]    = gj.fJetSigma2[i] ;
	fPhiEMCALCut[i]  = gj.fPhiEMCALCut[i] ;
      }
    }          
  } 
}

//____________________________________________________________________________
AliAnaGammaJet::~AliAnaGammaJet() 
{
  // Remove all pointers
  fOutputContainer->Clear() ; 
  delete fOutputContainer ;
  
  delete fhChargeRatio  ; 
  delete fhPi0Ratio   ; 
  delete fhDeltaPhiCharge  ;  
  delete fhDeltaPhiPi0   ; 
  delete fhDeltaEtaCharge  ; 
  delete fhDeltaEtaPi0  ; 
  delete fhAnglePair  ; 
  delete fhAnglePairAccepted  ; 
  delete fhAnglePairNoCut  ; 
  delete fhAnglePairLeadingCut  ; 
  delete fhAnglePairAngleCut   ; 
  delete fhAnglePairAllCut   ; 
  delete fhAnglePairLeading  ; 
  delete fhInvMassPairNoCut    ; 
  delete fhInvMassPairLeadingCut  ; 
  delete fhInvMassPairAngleCut  ; 
  delete fhInvMassPairAllCut   ; 
  delete fhInvMassPairLeading  ; 
  delete fhNBkg   ; 
  delete fhNLeading  ; 
  delete fhNJet  ; 
  delete fhJetRatio  ; 
  delete fhJetPt   ; 
  delete fhBkgRatio   ; 
  delete fhBkgPt  ; 
  delete fhJetFragment  ; 
  delete fhBkgFragment  ; 
  delete fhJetPtDist  ; 
  delete fhBkgPtDist  ; 
  
  delete [] fhJetRatios;  
  delete [] fhJetPts;  
  delete [] fhBkgRatios;
  delete [] fhBkgPts;  

  delete [] fhNLeadings;
  delete [] fhNJets;  
  delete [] fhNBkgs;
  
  delete [] fhJetFragments;
  delete [] fhBkgFragments;
  delete [] fhJetPtDists;
  delete [] fhBkgPtDists;
  
}

//____________________________________________________________________________
Double_t AliAnaGammaJet::CalculateJetRatioLimit(const Double_t ptg, 
						 const Double_t *par, 
						 const Double_t *x) {

  //Parametrized cut for the energy of the jet.

  Double_t epp = par[0] + par[1] * ptg ;
  Double_t spp = par[2] + par[3] * ptg ;
  Double_t f   = x[0]   + x[1]   * ptg ;
  Double_t epb = epp + par[4] ;
  Double_t spb = TMath::Sqrt(spp*spp+ par[5]*par[5]) ;
  Double_t rat = (epb - spb * f) / ptg ;
  //Info("CalculateLimit","epp %f, spp %f, f %f", epp, spp, f);
  //Info("CalculateLimit","epb %f, spb %f, rat %f", epb, spb, rat);
  return rat ;
}


//____________________________________________________________________________
void AliAnaGammaJet::Exec(Option_t *) 
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
  TClonesArray * plNe         = new TClonesArray("TParticle",1000);   // All particles measured in Jet Calorimeter
  TClonesArray * plCalo     = new TClonesArray("TParticle",1000);  // All particles measured in Prompt Gamma calorimeter 
  
  
  TParticle *pGamma = new TParticle(); //It will contain the kinematics of the found prompt gamma
  TParticle *pLeading = new TParticle(); //It will contain the kinematics of the found leading particle
  
  Bool_t iIsInPHOSorEMCAL = kFALSE ; //To check if Gamma was in any calorimeter
  
  //Fill lists with photons, neutral particles and charged particles
  //look for the highest energy photon in the event inside fCalorimeter
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
		  plCTS->GetEntries(),plNe->GetEntries(), GetCalorimeter().Data(), plCalo->GetEntries()));

  //If there is any prompt photon in fCalorimeter, 
  //search jet leading particle
  
  if(iIsInPHOSorEMCAL){
    if (GetPrintInfo())
      AliInfo(Form("Prompt Gamma: pt %f, phi %f, eta %f", pGamma->Pt(),pGamma->Phi(),pGamma->Eta())) ;
 
    AliDebug(2, "Get Leading Particles, Make jet");

    //Search leading particles in CTS and EMCAL 
    if(GetLeadingParticle(plCTS, plNe, pGamma, pLeading)){
      if (GetPrintInfo())
	AliInfo(Form("Leading: pt %f, phi %f, eta %f", pLeading->Pt(),pLeading->Phi(),pLeading->Eta())) ;

      //Search Jet
      if(!fSeveralConeAndPtCuts)
	MakeJet(particleList,pGamma,pLeading,"");
      else{
	for(Int_t icone = 0; icone<fNCone; icone++) {
	  for(Int_t ipt = 0; ipt<fNPt;ipt++) {  
	    TString lastname ="Cone"+ fNameCones[icone]+"Pt"+ fNamePtThres[ipt];
	    MakeJet(particleList, pGamma, pLeading,lastname);
	  }//icone
	}//ipt
      }//fSeveralConeAndPtCuts
    }//Leading
  }//Gamma in Calo
     
  AliDebug(2, "End of analysis, delete pointers");

  particleList->Delete() ; 
  plCTS->Delete() ;
  plNe->Delete() ;
  plCalo->Delete() ;
  pLeading->Delete();
  pGamma->Delete();

  delete plNe ;
  delete plCalo ;
  delete plCTS ;
  delete particleList ;
  //  delete pLeading;
  //  delete pGamma;

  PostData(0, fOutputContainer);
}    

//____________________________________________________________________________
void AliAnaGammaJet::FillJetHistos(TClonesArray * pl, Double_t ptg, Double_t ptl, TString type, TString lastname)
{
  //Fill histograms wth jet fragmentation 
  //and number of selected jets and leading particles
  //and the background multiplicity
  TParticle * particle = 0 ;
  Int_t ipr = 0;
  Float_t  charge = 0;

  TIter next(pl) ; 
  while ( (particle = dynamic_cast<TParticle*>(next())) ) {
    ipr++ ;
    Double_t pt = particle->Pt();

    charge = TDatabasePDG::Instance()
      ->GetParticle(particle->GetPdgCode())->Charge();
    if(charge != 0){//Only jet Charged particles 
      dynamic_cast<TH2F*>
	(fOutputContainer->FindObject(type+"Fragment"+lastname))
	->Fill(ptg,pt/ptg);
      dynamic_cast<TH2F*>
	(fOutputContainer->FindObject(type+"PtDist"+lastname))
	->Fill(ptg,pt);
    }//charged

  }//while

  if(type == "Bkg")
    dynamic_cast<TH1F*>
      (fOutputContainer->FindObject("NBkg"+lastname))
      ->Fill(ipr);
  else{
    dynamic_cast<TH1F*>
      (fOutputContainer->FindObject("NJet"+lastname))->
      Fill(ptg);
    dynamic_cast<TH2F*>
      (fOutputContainer->FindObject("NLeading"+lastname))
      ->Fill(ptg,ptl);
  }
  
}

//____________________________________________________________________________
void  AliAnaGammaJet::GetLeadingCharge(TClonesArray * pl, TParticle * pGamma, TParticle * pLeading) const 
{  
  //Search for the charged particle with highest with 
  //Phi=Phi_gamma-Pi and pT=0.1E_gamma 
  Double_t pt  = -100.;
  Double_t phi = -100.;

  for(Int_t ipr = 0;ipr < pl->GetEntries() ; ipr ++ ){

    TParticle * particle = dynamic_cast<TParticle *>(pl->At(ipr)) ;

    Double_t ptl  = particle->Pt();
    Double_t rat  = ptl/pGamma->Pt() ;
    Double_t phil = particle->Phi() ;
    Double_t phig = pGamma->Phi();

    //Selection within angular and energy limits
    if(((phig-phil)> fPhiMinCut) && ((phig-phil)<fPhiMaxCut) &&
        (rat > fRatioMinCut) && (rat < fRatioMaxCut)  && (ptl  > pt)) {
       phi = phil ;
       pt  = ptl ;
       pLeading->SetMomentum(particle->Px(),particle->Py(),particle->Pz(),particle->Energy());
       AliDebug(4,Form("Charge in CTS: pt %f eta %f phi %f pt/Eg %f \n", pt, particle->Eta(), phi,rat)) ;
     }
  }
  
  AliDebug(3,Form("Leading in CTS: pt %f eta %f phi %f pt/Eg %f \n", pt, pLeading->Eta(), phi,pt/pGamma->Pt())) ;

}


//____________________________________________________________________________
void  AliAnaGammaJet::GetLeadingPi0(TClonesArray * pl, TParticle * pGamma, TParticle * pLeading)  
{  

  //Search for the neutral pion with highest with 
  //Phi=Phi_gamma-Pi and pT=0.1E_gamma 
  Double_t pt  = -100.;
  Double_t phi = -100.;
  Double_t ptl = -100.;
  Double_t rat = -100.; 
  Double_t phil = -100. ;
  Double_t ptg  = pGamma->Pt();
  Double_t phig = pGamma->Phi();

  TIter next(pl);
  TParticle * particle = 0;
  
  Int_t iPrimary = -1;
  TLorentzVector gammai,gammaj;
  Double_t angle = 0., e = 0., invmass = 0.;
  Double_t anglef = 0., ef = 0., invmassf = 0.;
  Int_t ksPdg = 0;
  Int_t jPrimary=-1;

  while ( (particle = (TParticle*)next()) ) {
    iPrimary++;	  
    
    ksPdg = particle->GetPdgCode();
    ptl  = particle->Pt();
    if(ksPdg == 111){ //2 gamma overlapped, found with PID
      rat = ptl/ptg ;
      phil = particle->Phi() ;
      //Selection within angular and energy limits
      if((ptl> pt)&& (rat > fRatioMinCut) && (rat < fRatioMaxCut) && 
	 ((phig-phil)>fPhiMinCut)&&((phig-phil)<fPhiMaxCut)){
	phi = phil ;
	pt  = ptl ;
	pLeading->SetMomentum(particle->Px(),particle->Py(),particle->Pz(),particle->Energy());
	AliDebug(4,Form("Pi0 candidate: pt %f eta %f phi %f pt/Eg %f \n",  pLeading->Pt(), pLeading->Eta(),  pLeading->Phi(),  pLeading->Pt()/pGamma->Pt())) ;
      }// cuts
    }// pdg = 111
    else if(ksPdg == 22){//1 gamma
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
	    ptl  = (gammai+gammaj).Pt();
	    phil = (gammai+gammaj).Phi();
	    if(phil < 0)
	      phil+=TMath::TwoPi();
	    rat = ptl/ptg ;
	    invmass = (gammai+gammaj).M();
	    angle = gammaj.Angle(gammai.Vect());
	    //Info("GetLeadingPi0","Angle %f", angle);
	    e = (gammai+gammaj).E();
	    //Fill histograms with no cuts applied.
	    fhAnglePairNoCut->Fill(e,angle);
	    fhInvMassPairNoCut->Fill(ptg,invmass);
	    //First cut on the energy and azimuth of the pair
	    if((rat > fRatioMinCut) && (rat < fRatioMaxCut) && 
	       ((phig-phil)>fPhiMinCut)&&((phig-phil)<fPhiMaxCut)){
	      
	      fhAnglePairLeadingCut->Fill(e,angle);
	      fhInvMassPairLeadingCut->Fill(ptg,invmass);
	      //Second cut on the aperture of the pair
	      if(IsAngleInWindow(angle,e)){
		fhAnglePairAngleCut->Fill(e,angle);
		fhInvMassPairAngleCut->Fill(ptg,invmass);
		
		//Info("GetLeadingPi0","InvMass %f", invmass);
		//Third cut on the invariant mass of the pair
		if((invmass>fInvMassMinCut) && (invmass<fInvMassMaxCut)){ 
		  fhInvMassPairAllCut->Fill(ptg,invmass);
		  fhAnglePairAllCut->Fill(e,angle);
		  if(ptl > pt ){
		    pt       = ptl;
		    phi      = phil ;
		    ef       = e ;
		    anglef   = angle ;
		    invmassf = invmass ;
		    pLeading->SetMomentum(particle->Px(),particle->Py(),particle->Pz(),particle->Energy());
		    AliDebug(4,Form("Pi0 candidate: pt %f eta %f phi %f pt/Eg %f \n",  pLeading->Pt(), pLeading->Eta(),  pLeading->Phi(),  pLeading->Pt()/pGamma->Pt())) ;
		  }
		}//cuts
	      }//(invmass>0.125) && (invmass<0.145)
	    }//gammaj.Angle(gammai.Vect())<0.04
	  }//if pdg = 22
	}//iprimary<jprimary
      }//while
    }// if pdg = 22
    //     }
  }//while
  
  if(ef > 0.){//Final pi0 found, highest pair energy, fill histograms
    fhInvMassPairLeading->Fill(ptg,invmassf);
    fhAnglePairLeading->Fill(ef,anglef);
  }
  
  AliDebug(3,Form("Leading EMCAL: pt %f eta %f phi %f pt/Eg %f \n",  pLeading->Pt(), pLeading->Eta(),  pLeading->Phi(),  pLeading->Pt()/pGamma->Pt())) ;
}

//____________________________________________________________________________
Bool_t  AliAnaGammaJet::GetLeadingParticle(TClonesArray * plCTS, TClonesArray * plNe,  
					 TParticle * pGamma, TParticle * pLeading) 
{
  //Search Charged or Neutral leading particle, select the highest one.
  
  TParticle * pLeadingCh = new TParticle();
  TParticle * pLeadingPi0 = new TParticle();
  
  Double_t ptg  =  pGamma->Pt(); 
  Double_t phig = pGamma->Phi(); 
  Double_t etag = pGamma->Eta(); 
  
  if(GetCalorimeter() == "PHOS" && !fJetsOnlyInCTS)
    {
      AliDebug(3, "GetLeadingPi0");
      GetLeadingPi0   (plNe, pGamma, pLeadingPi0) ;
      AliDebug(3, "GetLeadingCharge");
      GetLeadingCharge(plCTS, pGamma, pLeadingCh) ;
      
      Double_t ptch = pLeadingCh->Pt(); 
      Double_t phich = pLeadingCh->Phi(); 
      Double_t etach = pLeadingCh->Eta(); 
      Double_t ptpi = pLeadingPi0->Pt(); 
      Double_t phipi = pLeadingPi0->Phi(); 
      Double_t etapi = pLeadingPi0->Eta(); 

      //Is leading cone inside EMCAL acceptance?
      
      Bool_t insidech = kFALSE ;
      if((phich - fCone) >  fPhiEMCALCut[0] && 
	 (phich + fCone) <  fPhiEMCALCut[1] && 
	(etach-fCone) < fEtaEMCALCut )
	insidech = kTRUE ;
      
      Bool_t insidepi = kFALSE ;
      if((phipi - fCone) >  fPhiEMCALCut[0] && 
	 (phipi + fCone) <  fPhiEMCALCut[1] &&
	(etapi-fCone) < fEtaEMCALCut )
	insidepi = kTRUE ;

      AliDebug(2,Form("Leading:  charged pt %f, pi0 pt  %f",ptch,ptpi)) ;
      
      if (ptch > 0 || ptpi > 0)
	{
	  if(insidech && (ptch > ptpi))
	    {
	      if (GetPrintInfo())
		AliInfo("Leading found in CTS");
	      pLeading->SetMomentum(pLeadingCh->Px(),pLeadingCh->Py(),pLeadingCh->Pz(),pLeadingCh->Energy());
	      AliDebug(3,Form("Final leading found in CTS, pt %f, phi %f, eta %f",ptch,phich,etach)) ;
	      fhChargeRatio->Fill(ptg,ptch/pGamma->Pt());
	      fhDeltaPhiCharge->Fill(ptg,pGamma->Phi()-phich);
	      fhDeltaEtaCharge->Fill(ptg,pGamma->Eta()-etach);
	      return 1;
	    }
	  
	  else if((ptpi > ptch) && insidepi)
	    {
	      if (GetPrintInfo())
		AliInfo("Leading found in EMCAL");
	      pLeading->SetMomentum(pLeadingPi0->Px(),pLeadingPi0->Py(),pLeadingPi0->Pz(),pLeadingPi0->Energy());
	      AliDebug(3,Form("Final leading found in EMCAL, pt %f, phi %f, eta %f",ptpi,phipi,etapi)) ;
	      fhPi0Ratio     ->Fill(ptg,ptpi/ptg);
	      fhDeltaPhiPi0->Fill(ptg,phig-phipi);
	      fhDeltaEtaPi0->Fill(ptg,etag-etapi);
	      return 1;	    
	    }

	  else{
	    AliDebug(3,"NO LEADING PARTICLE FOUND");}
	  return 0; 
	}
      else{
	AliDebug(3,"NO LEADING PARTICLE FOUND");
	return 0;
      }
    }

  else
    {
      //No calorimeter present for Leading particle detection
      GetLeadingCharge(plCTS, pGamma, pLeading) ;
      Double_t ptch = pLeading->Pt(); 
      Double_t phich = pLeading->Phi(); 
      Double_t etach = pLeading->Eta(); 
      if(ptch > 0){
	fhChargeRatio->Fill(ptg,ptch/ptg);
	fhDeltaPhiCharge->Fill(ptg,phig-phich);
	fhDeltaEtaCharge->Fill(ptg,etag-etach);
	AliDebug(3,Form("Leading found :  pt %f, phi %f, eta %f",ptch,phich,etach)) ;
	return 1;
      }
      else
	{
	  AliDebug(3,"NO LEADING PARTICLE FOUND");	
	  return 0;
	}
    }
}

  //____________________________________________________________________________
void AliAnaGammaJet::Init(const Option_t * )
{
//   // Initialisation of branch container 
  AliAnaGammaDirect::Init();
 
  //Initialize the parameters of the analysis.
  //fCalorimeter="PHOS";
  fAngleMaxParam.Set(4) ;
  fAngleMaxParam.AddAt(0.4,0);//={0.4,-0.25,0.025,-2e-4};
  fAngleMaxParam.AddAt(-0.25,1) ;
  fAngleMaxParam.AddAt(0.025,2) ;
  fAngleMaxParam.AddAt(-2e-4,3) ;
  fSeveralConeAndPtCuts   = kFALSE ;
  //fPrintInfo           = kTRUE;
  fPbPb                = kFALSE ;
  fInvMassMaxCut  = 0.16 ;
  fInvMassMinCut  = 0.11 ;
  //fJetsOnlyInCTS    = kTRUE ;
  fEtaEMCALCut     = 0.7 ;
  fPhiEMCALCut[0] = 80. *TMath::Pi()/180.;
  fPhiEMCALCut[1] = 190.*TMath::Pi()/180.;
  fPhiMaxCut      = 3.4 ;
  fPhiMinCut      = 2.9 ;

  //Jet selection parameters
  //Fixed cut (old)
  fRatioMaxCut    = 1.0 ;
  fRatioMinCut    = 0.1 ; 
  fJetRatioMaxCut = 1.2 ; 
  fJetRatioMinCut = 0.3 ; 
  fJetsOnlyInCTS = kFALSE ;
  fJetCTSRatioMaxCut = 1.2 ;
  fJetCTSRatioMinCut = 0.3 ;
  fSelect         = 0  ;

  //Cut depending on gamma energy

  fPtJetSelectionCut = 20.; //For Low pt jets+BKG, another limits applied
  //Reconstructed jet energy dependence parameters 
  //e_jet = a1+e_gamma b2. 
  //Index 0-> Pt>2 GeV r = 0.3; Index 1-> Pt>0.5 GeV r = 0.3
  fJetE1[0] = -5.75; fJetE1[1] = -4.1;
  fJetE2[0] = 1.005; fJetE2[1] = 1.05;

  //Reconstructed sigma of jet energy dependence parameters 
  //s_jet = a1+e_gamma b2. 
  //Index 0-> Pt>2 GeV r = 0.3; Index 1-> Pt>0.5 GeV r = 0.3
  fJetSigma1[0] = 2.65;   fJetSigma1[1] = 2.75;
  fJetSigma2[0] = 0.0018; fJetSigma2[1] = 0.033;

  //Background mean energy and RMS
  //Index 0-> No BKG; Index 1-> BKG > 2 GeV; 
  //Index 2-> (low pt jets)BKG > 0.5 GeV;
  //Index > 2, same for CTS conf
  fBkgMean[0] = 0.; fBkgMean[1] = 8.8 ; fBkgMean[2] = 69.5;
  fBkgMean[3] = 0.; fBkgMean[4] = 6.4;  fBkgMean[5] = 48.6;
  fBkgRMS[0]  = 0.; fBkgRMS[1]  = 7.5;  fBkgRMS[2]  = 22.0; 
  fBkgRMS[3]  = 0.; fBkgRMS[4]  = 5.4;  fBkgRMS[5]  = 13.2; 

  //Factor x of min/max = E -+ x * sigma. Obtained after selecting the
  //limits for monoenergetic jets.
  //Index 0-> No BKG; Index 1-> BKG > 2 GeV; 
  //Index 2-> (low pt jets) BKG > 0.5 GeV;
  //Index > 2, same for CTS conf

  fJetXMin1[0] =-0.69 ; fJetXMin1[1] = 0.39 ; fJetXMin1[2] =-0.88 ; 
  fJetXMin1[3] =-2.0  ; fJetXMin1[4] =-0.442 ; fJetXMin1[5] =-1.1  ;
  fJetXMin2[0] = 0.066; fJetXMin2[1] = 0.038; fJetXMin2[2] = 0.034; 
  fJetXMin2[3] = 0.25 ; fJetXMin2[4] = 0.113; fJetXMin2[5] = 0.077 ;
  fJetXMax1[0] =-3.8  ; fJetXMax1[1] =-0.76 ; fJetXMax1[2] =-3.6  ; 
  fJetXMax1[3] =-2.7  ; fJetXMax1[4] =-1.21 ; fJetXMax1[5] =-3.7  ;
  fJetXMax2[0] =-0.012; fJetXMax2[1] =-0.022; fJetXMax2[2] = 0.016; 
  fJetXMax2[3] =-0.024; fJetXMax2[4] =-0.008; fJetXMax2[5] = 0.027;


  //Different cones and pt thresholds to construct the jet

  fCone        = 0.3  ;
  fPtThreshold = 0.   ;
  fNCone       = 3    ;
  fNPt         = 3    ;
  fCones[1]    = 0.2  ; fNameCones[1]   = "02" ;
  fPtThres[0]  = 0.0  ; fNamePtThres[0] = "00" ;
  fCones[0]    = 0.3  ; fNameCones[0]   = "03" ;
  fPtThres[1]  = 0.5  ; fNamePtThres[1] = "05" ;
  fCones[2]    = 0.4  ; fNameCones[2]   = "04" ;
  fPtThres[2]  = 1.0  ; fNamePtThres[2] = "10" ;
  //Fill particle lists when PID is ok
  // fEMCALPID = kFALSE;
  // fPHOSPID = kFALSE;

  //Initialization of histograms 

  MakeHistos() ; 

}

//__________________________________________________________________________-
Bool_t AliAnaGammaJet::IsAngleInWindow(const Float_t angle,const Float_t e) {
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

//__________________________________________________________________________-
Bool_t AliAnaGammaJet::IsJetSelected(const Double_t ptg, const Double_t ptj){
  //Check if the energy of the reconstructed jet is within an energy window

  Double_t par[6];
  Double_t xmax[2];
  Double_t xmin[2];

  Int_t iCTS = 0;
  if(fJetsOnlyInCTS)
    iCTS = 3 ;

  if(!fPbPb){
    //Phythia alone, jets with pt_th > 0.2, r = 0.3 
    par[0] = fJetE1[0]; par[1] = fJetE2[0]; 
    //Energy of the jet peak
    //e_jet = fJetE1[0]+fJetE2[0]*e_gamma, simulation fit
    par[2] = fJetSigma1[0]; par[3] = fJetSigma2[0];
    //Sigma  of the jet peak
    //sigma_jet = fJetSigma1[0]+fJetSigma2[0]*e_gamma, simulation fit
    par[4] = fBkgMean[0 + iCTS]; par[5] = fBkgRMS[0 + iCTS];
    //Parameters reserved for PbPb bkg.
    xmax[0] = fJetXMax1[0 + iCTS]; xmax[1] = fJetXMax2[0 + iCTS];
    xmin[0] = fJetXMin1[0 + iCTS]; xmin[1] = fJetXMin2[0 + iCTS];
    //Factor that multiplies sigma to obtain the best limits, 
    //by observation, of mono jet ratios (ptjet/ptg)
    //X_jet = fJetX1[0]+fJetX2[0]*e_gamma
   
  }
  else{
    if(ptg > fPtJetSelectionCut){
      //Phythia +PbPb with  pt_th > 2 GeV/c, r = 0.3 
      par[0] = fJetE1[0]; par[1] = fJetE2[0]; 
      //Energy of the jet peak, same as in pp
      //e_jet = fJetE1[0]+fJetE2[0]*e_gamma, simulation fit
      par[2] = fJetSigma1[0]; par[3] = fJetSigma2[0];
      //Sigma  of the jet peak, same as in pp
      //sigma_jet = fJetSigma1[0]+fJetSigma2[0]*e_gamma, simulation fit
      par[4] = fBkgMean[1 + iCTS]; par[5] = fBkgRMS[1 + iCTS];
      //Mean value and RMS of PbPb Bkg 
      xmax[0] = fJetXMax1[1 + iCTS]; xmax[1] = fJetXMax2[1 + iCTS];
      xmin[0] = fJetXMin1[1 + iCTS]; xmin[1] = fJetXMin2[1 + iCTS];
      //Factor that multiplies sigma to obtain the best limits, 
      //by observation, of mono jet ratios (ptjet/ptg) mixed with PbPb Bkg, 
      //pt_th > 2 GeV, r = 0.3
      //X_jet = fJetX1[0]+fJetX2[0]*e_gamma
     
    }
    else{
      //Phythia + PbPb with  pt_th > 0.5 GeV/c, r = 0.3
      par[0] = fJetE1[1]; par[1] = fJetE2[1]; 
      //Energy of the jet peak, pt_th > 2 GeV/c, r = 0.3 
      //e_jet = fJetE1[0]+fJetE2[0]*e_gamma, simulation fit
      par[2] = fJetSigma1[1]; par[3] = fJetSigma2[1];
      //Sigma  of the jet peak, pt_th > 2 GeV/c, r = 0.3
      //sigma_jet = fJetSigma1[0]+fJetSigma2[0]*e_gamma, simulation fit
      par[4] = fBkgMean[2 + iCTS]; par[5] = fBkgRMS[2 + iCTS];
      //Mean value and RMS of PbPb Bkg in a 0.3 cone, pt > 2 GeV.
      xmax[0] = fJetXMax1[2 + iCTS]; xmax[1] = fJetXMax2[2 + iCTS];
      xmin[0] = fJetXMin1[2 + iCTS]; xmin[1] = fJetXMin2[2 + iCTS];
      //Factor that multiplies sigma to obtain the best limits, 
      //by observation, of mono jet ratios (ptjet/ptg) mixed with PbPb Bkg, 
      //pt_th > 2 GeV, r = 0.3
      //X_jet = fJetX1[0]+fJetX2[0]*e_gamma
     
    }//If low pt jet in bkg
  }//if Bkg

 //Calculate minimum and maximum limits of the jet ratio.
  Double_t min = CalculateJetRatioLimit(ptg, par, xmin);
  Double_t max = CalculateJetRatioLimit(ptg, par, xmax);
  
  AliDebug(3,Form("Jet selection?  : Limits min %f, max %f,  pt_jet %f,  pt_gamma %f, pt_jet / pt_gamma %f",min,max,ptj,ptg,ptj/ptg));

  if(( min < ptj/ptg ) && ( max > ptj/ptg))
    return kTRUE;
  else
    return kFALSE;

}

//____________________________________________________________________________
void AliAnaGammaJet::MakeHistos()
{
  // Create histograms to be saved in output file and 
  // stores them in fOutputContainer
  
  fOutputContainer = new TObjArray(10000) ;
 
  TObjArray  * outputContainer =GetOutputContainer();
  for(Int_t i = 0; i < outputContainer->GetEntries(); i++ )
    fOutputContainer->Add(outputContainer->At(i)) ;

  //
  fhChargeRatio  = new TH2F
    ("ChargeRatio","p_{T leading charge} /p_{T #gamma} vs p_{T #gamma}",
     120,0,120,120,0,1); 
  fhChargeRatio->SetYTitle("p_{T lead charge} /p_{T #gamma}");
  fhChargeRatio->SetXTitle("p_{T #gamma} (GeV/c)");
  fOutputContainer->Add(fhChargeRatio) ;
  
  fhDeltaPhiCharge  = new TH2F
    ("DeltaPhiCharge","#phi_{#gamma} - #phi_{charge} vs p_{T #gamma}",
     200,0,120,200,0,6.4); 
  fhDeltaPhiCharge->SetYTitle("#Delta #phi");
  fhDeltaPhiCharge->SetXTitle("p_{T #gamma} (GeV/c)");
  fOutputContainer->Add(fhDeltaPhiCharge) ; 
  
  fhDeltaEtaCharge  = new TH2F
    ("DeltaEtaCharge","#eta_{#gamma} - #eta_{charge} vs p_{T #gamma}",
     200,0,120,200,-2,2); 
  fhDeltaEtaCharge->SetYTitle("#Delta #eta");
  fhDeltaEtaCharge->SetXTitle("p_{T #gamma} (GeV/c)");
  fOutputContainer->Add(fhDeltaEtaCharge) ; 
  
  //
  if(!fJetsOnlyInCTS){
    fhPi0Ratio  = new TH2F
      ("Pi0Ratio","p_{T leading  #pi^{0}} /p_{T #gamma} vs p_{T #gamma}",
       120,0,120,120,0,1); 
    fhPi0Ratio->SetYTitle("p_{T lead  #pi^{0}} /p_{T #gamma}");
    fhPi0Ratio->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhPi0Ratio) ; 
    
    fhDeltaPhiPi0  = new TH2F
      ("DeltaPhiPi0","#phi_{#gamma} - #phi_{ #pi^{0}} vs p_{T #gamma}",
       200,0,120,200,0,6.4); 
    fhDeltaPhiPi0->SetYTitle("#Delta #phi");
    fhDeltaPhiPi0->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhDeltaPhiPi0) ; 
    
    fhDeltaEtaPi0  = new TH2F
      ("DeltaEtaPi0","#eta_{#gamma} - #eta_{ #pi^{0}} vs p_{T #gamma}",
       200,0,120,200,-2,2); 
    fhDeltaEtaPi0->SetYTitle("#Delta #eta");
    fhDeltaEtaPi0->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhDeltaEtaPi0) ; 
 
    //
    fhAnglePair  = new TH2F
      ("AnglePair",
       "Angle between #pi^{0} #gamma pair vs p_{T  #pi^{0}}",
       200,0,50,200,0,0.2); 
    fhAnglePair->SetYTitle("Angle (rad)");
    fhAnglePair->SetXTitle("E_{ #pi^{0}} (GeV/c)");
    fOutputContainer->Add(fhAnglePair) ; 
    
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
    
    fhAnglePairLeadingCut  = new TH2F
      ("AnglePairLeadingCut",
       "Angle between all #gamma pair that have a good phi and pt vs p_{T  #pi^{0}}",
       200,0,50,200,0,0.2); 
    fhAnglePairLeadingCut->SetYTitle("Angle (rad)");
    fhAnglePairLeadingCut->SetXTitle("E_{ #pi^{0}} (GeV/c)");
    fOutputContainer->Add(fhAnglePairLeadingCut) ; 
    
    fhAnglePairAngleCut  = new TH2F
      ("AnglePairAngleCut",
       "Angle between all #gamma pair (angle + leading cut) vs p_{T  #pi^{0}}"
       ,200,0,50,200,0,0.2); 
    fhAnglePairAngleCut->SetYTitle("Angle (rad)");
    fhAnglePairAngleCut->SetXTitle("E_{ #pi^{0}} (GeV/c)");
    fOutputContainer->Add(fhAnglePairAngleCut) ;
    
    fhAnglePairAllCut  = new TH2F
      ("AnglePairAllCut",
       "Angle between all #gamma pair (angle + inv mass cut+leading) vs p_{T  #pi^{0}}"
       ,200,0,50,200,0,0.2); 
    fhAnglePairAllCut->SetYTitle("Angle (rad)");
    fhAnglePairAllCut->SetXTitle("E_{ #pi^{0}} (GeV/c)");
    fOutputContainer->Add(fhAnglePairAllCut) ; 
    
    fhAnglePairLeading  = new TH2F
      ("AnglePairLeading",
       "Angle between all #gamma pair finally selected vs p_{T  #pi^{0}}",
       200,0,50,200,0,0.2); 
    fhAnglePairLeading->SetYTitle("Angle (rad)");
    fhAnglePairLeading->SetXTitle("E_{ #pi^{0}} (GeV/c)");
    fOutputContainer->Add(fhAnglePairLeading) ; 
    
    //
    fhInvMassPairNoCut  = new TH2F
      ("InvMassPairNoCut","Invariant Mass of all #gamma pair vs p_{T #gamma}",
       120,0,120,360,0,0.5); 
    fhInvMassPairNoCut->SetYTitle("Invariant Mass (GeV/c^{2})");
    fhInvMassPairNoCut->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhInvMassPairNoCut) ; 
    
    fhInvMassPairLeadingCut  = new TH2F
      ("InvMassPairLeadingCut",
       "Invariant Mass of #gamma pair (leading cuts) vs p_{T #gamma}",
       120,0,120,360,0,0.5); 
    fhInvMassPairLeadingCut->SetYTitle("Invariant Mass (GeV/c^{2})");
    fhInvMassPairLeadingCut->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhInvMassPairLeadingCut) ; 
    
    fhInvMassPairAngleCut  = new TH2F
      ("InvMassPairAngleCut",
       "Invariant Mass of #gamma pair (angle cut) vs p_{T #gamma}",
       120,0,120,360,0,0.5); 
    fhInvMassPairAngleCut->SetYTitle("Invariant Mass (GeV/c^{2})");
    fhInvMassPairAngleCut->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhInvMassPairAngleCut) ; 
    
    fhInvMassPairAllCut  = new TH2F
      ("InvMassPairAllCut",
       "Invariant Mass of #gamma pair (angle+invmass cut+leading) vs p_{T #gamma}",
       120,0,120,360,0,0.5); 
    fhInvMassPairAllCut->SetYTitle("Invariant Mass (GeV/c^{2})");
    fhInvMassPairAllCut->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhInvMassPairAllCut) ; 
    
    fhInvMassPairLeading  = new TH2F
      ("InvMassPairLeading",
       "Invariant Mass of #gamma pair selected vs p_{T #gamma}",
       120,0,120,360,0,0.5); 
    fhInvMassPairLeading->SetYTitle("Invariant Mass (GeV/c^{2})");
    fhInvMassPairLeading->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhInvMassPairLeading) ; 
  }
  
  
  if(!fSeveralConeAndPtCuts){
    
    //Count
    fhNBkg = new TH1F("NBkg","bkg multiplicity",9000,0,9000); 
    fhNBkg->SetYTitle("counts");
    fhNBkg->SetXTitle("N");
    fOutputContainer->Add(fhNBkg) ; 
    
    fhNLeading  = new TH2F
      ("NLeading","Accepted Jet Leading", 240,0,120,240,0,120); 
    fhNLeading->SetYTitle("p_{T charge} (GeV/c)");
    fhNLeading->SetXTitle("p_{T #gamma}(GeV/c)");
    fOutputContainer->Add(fhNLeading) ; 
    
    fhNJet  = new TH1F("NJet","Accepted jets",240,0,120); 
    fhNJet->SetYTitle("N");
    fhNJet->SetXTitle("p_{T #gamma}(GeV/c)");
    fOutputContainer->Add(fhNJet) ; 
    
    //Ratios and Pt dist of reconstructed (not selected) jets
    //Jet
    fhJetRatio  = new TH2F
      ("JetRatio","p_{T jet lead}/p_{T #gamma} vs p_{T #gamma}",
       240,0,120,200,0,10);
    fhJetRatio->SetYTitle("p_{T jet lead #pi^{0}}/p_{T #gamma}");
    fhJetRatio->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhJetRatio) ; 
    
    fhJetPt  = new TH2F
      ("JetPt", "p_{T jet lead} vs p_{T #gamma}",240,0,120,400,0,200);
    fhJetPt->SetYTitle("p_{T jet}");
    fhJetPt->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhJetPt) ; 
    
    //Bkg
    
    fhBkgRatio  = new TH2F
      ("BkgRatio","p_{T bkg lead}/p_{T #gamma} vs p_{T #gamma}",
       240,0,120,200,0,10);
    fhBkgRatio->SetYTitle("p_{T bkg lead charge}/p_{T #gamma}");
    fhBkgRatio->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhBkgRatio) ;
    
    fhBkgPt  = new TH2F
      ("BkgPt","p_{T jet lead} vs p_{T #gamma}",240,0,120,400,0,200);
    fhBkgPt->SetYTitle("p_{T jet lead charge}/p_{T #gamma}");
    fhBkgPt->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhBkgPt) ;
    
    //Jet Distributions
    
    fhJetFragment  = 
      new TH2F("JetFragment","x = p_{T i charged}/p_{T #gamma}",
	       240,0.,120.,1000,0.,1.2); 
    fhJetFragment->SetYTitle("x_{T}");
    fhJetFragment->SetXTitle("p_{T #gamma}");
    fOutputContainer->Add(fhJetFragment) ;
    
    fhBkgFragment  = new TH2F
      ("BkgFragment","x = p_{T i charged}/p_{T #gamma}",
       240,0.,120.,1000,0.,1.2);
    fhBkgFragment->SetYTitle("x_{T}");
    fhBkgFragment->SetXTitle("p_{T #gamma}");
    fOutputContainer->Add(fhBkgFragment) ;
    
    fhJetPtDist  = 
      new TH2F("JetPtDist","x = p_{T i charged}",240,0.,120.,400,0.,200.); 
    fhJetPtDist->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhJetPtDist) ;
    
    fhBkgPtDist  = new TH2F
      ("BkgPtDist","x = p_{T i charged}",240,0.,120.,400,0.,200.); 
    fhBkgPtDist->SetXTitle("p_{T #gamma} (GeV/c)");
    fOutputContainer->Add(fhBkgPtDist) ;

  }
  else{
    //If we want to study the jet for different cones and pt
    
    for(Int_t icone = 0; icone<fNCone; icone++){
      for(Int_t ipt = 0; ipt<fNPt;ipt++){ 
	
	//Jet
	
	fhJetRatios[icone][ipt]  = new TH2F
	  ("JetRatioCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt], 
	   "p_{T jet lead #pi^{0}}/p_{T #gamma} vs p_{T #gamma}, cone ="
	   +fNameCones[icone]+", pt>" +fNamePtThres[ipt]+" GeV/c",
	   240,0,120,200,0,10);
	fhJetRatios[icone][ipt]->
	  SetYTitle("p_{T jet lead #pi^{0}}/p_{T #gamma}");
	fhJetRatios[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	fOutputContainer->Add(fhJetRatios[icone][ipt]) ; 
	
	
	fhJetPts[icone][ipt]  = new TH2F
	  ("JetPtCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt], 
	   "p_{T jet lead #pi^{0}}/p_{T #gamma} vs p_{T #gamma}, cone ="
	   +fNameCones[icone]+", pt>" +fNamePtThres[ipt]+" GeV/c",
	   240,0,120,400,0,200);
	fhJetPts[icone][ipt]->
	  SetYTitle("p_{T jet lead #pi^{0}}/p_{T #gamma}");
	fhJetPts[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	fOutputContainer->Add(fhJetPts[icone][ipt]) ; 
	
	//Bkg
	fhBkgRatios[icone][ipt]  = new TH2F
	  ("BkgRatioCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt], 
	   "p_{T bkg lead #pi^{0}}/p_{T #gamma} vs p_{T #gamma}, cone ="
	   +fNameCones[icone]+", pt>" +fNamePtThres[ipt]+" GeV/c",
	   240,0,120,200,0,10);
	fhBkgRatios[icone][ipt]->
	  SetYTitle("p_{T bkg lead #pi^{0}}/p_{T #gamma}");
	fhBkgRatios[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	fOutputContainer->Add(fhBkgRatios[icone][ipt]) ; 
	
	fhBkgPts[icone][ipt]  = new TH2F
	  ("BkgPtCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt], 
	   "p_{T jet lead #pi^{0}}/p_{T #gamma} vs p_{T #gamma}, cone ="
	   +fNameCones[icone]+", pt>" +fNamePtThres[ipt]+" GeV/c",
	   240,0,120,400,0,200);
	fhBkgPts[icone][ipt]->
	  SetYTitle("p_{T jet lead #pi^{0}}/p_{T #gamma}");
	fhBkgPts[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	fOutputContainer->Add(fhBkgPts[icone][ipt]) ; 
	
	//Counts
	fhNBkgs[icone][ipt]  = new TH1F
	  ("NBkgCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	   "bkg multiplicity cone ="+fNameCones[icone]+", pt>" 
	   +fNamePtThres[ipt]+" GeV/c",9000,0,9000); 
	fhNBkgs[icone][ipt]->SetYTitle("counts");
	fhNBkgs[icone][ipt]->SetXTitle("N");
	fOutputContainer->Add(fhNBkgs[icone][ipt]) ; 
	
	fhNLeadings[icone][ipt]  = new TH2F
	  ("NLeadingCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	   "p_{T #gamma} vs p_{T #pi^{0}} cone ="+fNameCones[icone]+", pt>" 
	   +fNamePtThres[ipt]+" GeV/c",120,0,120,120,0,120); 
	fhNLeadings[icone][ipt]->SetYTitle("p_{T #pi^{0}}(GeV/c)");
	fhNLeadings[icone][ipt]->SetXTitle("p_{T #gamma}(GeV/c)");
	fOutputContainer->Add(fhNLeadings[icone][ipt]) ; 
	
	fhNJets[icone][ipt]  = new TH1F
	  ("NJetCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	   "Number of neutral jets, cone ="+fNameCones[icone]+", pt>" 
	   +fNamePtThres[ipt]+" GeV/c",120,0,120); 
	fhNJets[icone][ipt]->SetYTitle("N");
	fhNJets[icone][ipt]->SetXTitle("p_{T #gamma}(GeV/c)");
	fOutputContainer->Add(fhNJets[icone][ipt]) ; 
	
	//Fragmentation Function
	fhJetFragments[icone][ipt]  = new TH2F
	  ("JetFragmentCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	   "x_{T} = p_{T i}/p_{T #gamma}, cone ="+fNameCones[icone]+", pt>" 
	   +fNamePtThres[ipt]+" GeV/c",120,0.,120.,240,0.,1.2); 
	fhJetFragments[icone][ipt]->SetYTitle("x_{T}");
	fhJetFragments[icone][ipt]->SetXTitle("p_{T #gamma}");
	fOutputContainer->Add(fhJetFragments[icone][ipt]) ; 
	
	fhBkgFragments[icone][ipt]  = new TH2F
	  ("BkgFragmentCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	   "x_{T} = p_{T i}/p_{T #gamma}, cone ="+fNameCones[icone]+", pt>" 
	   +fNamePtThres[ipt]+" GeV/c",120,0.,120.,240,0.,1.2); 
	fhBkgFragments[icone][ipt]->SetYTitle("x_{T}");
	fhBkgFragments[icone][ipt]->SetXTitle("p_{T #gamma}");
	fOutputContainer->Add(fhBkgFragments[icone][ipt]) ; 
	
	//Jet particle distribution
	
	fhJetPtDists[icone][ipt]  = new TH2F
	  ("JetPtDistCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	   "p_{T i}, cone ="+fNameCones[icone]+", pt>" +fNamePtThres[ipt]+
	   " GeV/c",120,0.,120.,120,0.,120.); 
	fhJetPtDists[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	fOutputContainer->Add(fhJetPtDists[icone][ipt]) ; 
	
	fhBkgPtDists[icone][ipt]  = new TH2F
	  ("BkgPtDistCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	   "p_{T i}, cone ="+fNameCones[icone]+", pt>" +fNamePtThres[ipt]+
	   " GeV/c",120,0.,120.,120,0.,120.); 
	fhBkgPtDists[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	fOutputContainer->Add(fhBkgPtDists[icone][ipt]) ; 
	
      }//ipt
    } //icone
  }//If we want to study any cone or pt threshold
}

//____________________________________________________________________________
void AliAnaGammaJet::MakeJet(TClonesArray * pl, TParticle * pGamma, TParticle* pLeading,TString lastname)
{
  //Fill the jet with the particles around the leading particle with 
  //R=fCone and pt_th = fPtThres. Caclulate the energy of the jet and 
  //check if we select it. Fill jet histograms
  
  TClonesArray * jetList = new TClonesArray("TParticle",1000);
  TClonesArray * bkgList = new TClonesArray("TParticle",1000);

  TLorentzVector jet   (0,0,0,0);  
  TLorentzVector bkg(0,0,0,0);
  TLorentzVector lv (0,0,0,0);

  Double_t ptjet = 0.0;
  Double_t ptbkg = 0.0;
  Int_t n0 = 0;
  Int_t n1 = 0;  
  Bool_t b1 = kFALSE;
  Bool_t b0 = kFALSE;
  
  Double_t ptg  = pGamma->Pt();
  Double_t phig = pGamma->Phi();
  Double_t ptl  = pLeading->Pt();
  Double_t phil = pLeading->Phi();
  Double_t etal = pLeading->Eta();

  Float_t ptcut = 0. ;
  if(fPbPb){
    if(ptg > fPtJetSelectionCut)  ptcut = 2. ;
    else                          ptcut = 0.5;
  }

  TIter next(pl) ; 
  TParticle * particle = 0 ; 
  while ( (particle = dynamic_cast<TParticle*>(next())) ) {
    
    b0 = kFALSE;
    b1 = kFALSE;

    //Particles in jet 
    SetJet(particle, b0, fCone, etal, phil) ;  

    if(b0){
      new((*jetList)[n0++]) TParticle(*particle) ;
      particle->Momentum(lv);
      if(particle->Pt() > ptcut ){
	jet+=lv;
	ptjet+=particle->Pt();
      }
    }

    //Background around (phi_gamma-pi, eta_leading)
    SetJet(particle, b1, fCone,etal, phig) ;

    if(b1) { 
      new((*bkgList)[n1++]) TParticle(*particle) ;
      particle->Momentum(lv);
      if(particle->Pt() > ptcut ){
	bkg+=lv;
	ptbkg+=particle->Pt();    
      }  
    }
  }
  
  ptjet = jet.Pt();
  ptbkg = bkg.Pt();

  if(ptjet > 0.) {

    AliDebug(2,Form("Gamma   pt %f, Jet pt %f, Bkg pt %f",ptg,ptjet,ptbkg));
    
    //Fill histograms
    
    Double_t ratjet   = ptjet/ptg ;
    Double_t ratbkg  = ptbkg/ptg ;
    
    dynamic_cast<TH2F*>
      (fOutputContainer->FindObject("JetRatio"+lastname))
      ->Fill(ptg,ratjet);	 
    dynamic_cast<TH2F*>
      (fOutputContainer->FindObject("JetPt"+lastname))
      ->Fill(ptg,ptjet);
    
    dynamic_cast<TH2F*>
      (fOutputContainer->FindObject("BkgRatio"+lastname))
      ->Fill(ptg,ratbkg);
    
    dynamic_cast<TH2F*>
      (fOutputContainer->FindObject("BkgPt"+lastname))
      ->Fill(ptg,ptbkg);


    //Jet selection
    Bool_t kSelect = kFALSE;
    if(fSelect == 0)
      kSelect = kTRUE; //Accept all jets, no restriction
    else if(fSelect == 1){
      //Selection with parametrized cuts
      if(IsJetSelected(ptg,ptjet))   kSelect = kTRUE;
    }
    else if(fSelect == 2){
      //Simple selection
      if(!fJetsOnlyInCTS){
	if((ratjet <  fJetRatioMaxCut) && (ratjet > fJetRatioMinCut )) kSelect = kTRUE;
      }
      else{
	if((ratjet <  fJetCTSRatioMaxCut) && (ratjet > fJetCTSRatioMinCut )) kSelect = kTRUE;
      }
    }
    else
      AliError("Jet selection option larger than 2, DONT SELECT JETS");
    
    
    if(kSelect){
      if (GetPrintInfo())
	AliInfo(Form("Jet Selected: pt %f ", ptjet)) ;
      
      FillJetHistos(jetList, ptg, ptl,"Jet",lastname);
      FillJetHistos(bkgList, ptg, ptl, "Bkg",lastname);
    }
  } //ptjet > 0
  
  jetList ->Delete();
  bkgList ->Delete();
  
}

//____________________________________________________________________________
void AliAnaGammaJet::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;

  Info("Print", "%s %s", GetName(), GetTitle() ) ;
  if(!fJetsOnlyInCTS)
    printf("EMCAL Acceptance cut  : phi min %f, phi max %f, eta %f \n", fPhiEMCALCut[0],  fPhiEMCALCut[1], fEtaEMCALCut) ;
  
  printf("Phi_Leading                                <     %f\n", fPhiMaxCut) ; 
  printf("Phi_Leading                                >     %f\n", fPhiMinCut) ;
  printf("pT Leading / pT Gamma             <     %f\n", fRatioMaxCut) ; 
  printf("pT Leading / pT Gamma             >     %f\n", fRatioMinCut) ;
  printf("M_pair                                        <     %f\n", fInvMassMaxCut) ; 
  printf("M_pair                                        >     %f\n", fInvMassMinCut) ; 
  if(fSelect == 2){
    printf("pT Jet / pT Gamma                     <    %f\n", fJetRatioMaxCut) ; 
    printf("pT Jet / pT Gamma                     >    %f\n", fJetRatioMinCut) ;
    printf("pT Jet (Only CTS)/ pT Gamma   <    %f\n", fJetCTSRatioMaxCut) ; 
    printf("pT Jet (Only CTS)/ pT Gamma   >    %f\n", fJetCTSRatioMinCut) ;
  }

 
} 

//___________________________________________________________________
void AliAnaGammaJet::SetJet(TParticle * part, Bool_t & b, Float_t cone, 
			     Double_t eta, Double_t phi)
{

  //Check if the particle is inside the cone defined by the leading particle
  b = kFALSE;
  
  if(phi > TMath::TwoPi())
    phi-=TMath::TwoPi();
  if(phi < 0.)
    phi+=TMath::TwoPi();
  
  Double_t  rad = 10000 + cone;
  
  if(TMath::Abs(part->Phi()-phi) <= (TMath::TwoPi() - cone))
    rad = TMath::Sqrt(TMath::Power(part->Eta()-eta,2)+
		      TMath::Power(part->Phi()-phi,2));
  else{
    if(part->Phi()-phi > TMath::TwoPi() - cone)
      rad = TMath::Sqrt(TMath::Power(part->Eta()-eta,2)+
			TMath::Power((part->Phi()-TMath::TwoPi())-phi,2));
    if(part->Phi()-phi < -(TMath::TwoPi() - cone))
      rad = TMath::Sqrt(TMath::Power(part->Eta()-eta,2)+
			TMath::Power((part->Phi()+TMath::TwoPi())-phi,2));
  }

  if(rad < cone )
    b = kTRUE;
  
}


void AliAnaGammaJet::Terminate(Option_t *)
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
    

}
