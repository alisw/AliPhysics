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
 * Revision 1.13  2006/04/26 07:32:37  hristov
 * Coding conventions, clean-up and related changes
 *
 * Revision 1.12  2006/03/10 13:23:36  hristov
 * Using AliESDCaloCluster instead of AliESDtrack
 *
 * Revision 1.11  2006/01/31 20:30:52  hristov
 * Including TFile.h
 *
 * Revision 1.10  2006/01/23 18:04:08  hristov
 * Removing meaningless const
 *
 * Revision 1.9  2006/01/12 16:23:26  schutz
 * ESD is properly read with methods of macros/AliReadESD.C copied in it
 *
 * Revision 1.8  2005/12/20 07:08:32  schutz
 * corrected error in call AliReadESD
 *
 * Revision 1.6  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
// Class for the analysis of gamma-jet correlations 
//  Basically it seaches for a prompt photon in the PHOS acceptance, 
//  if so we construct a jet around the highest pt particle in the opposite 
//  side in azimuth, inside the TPC and EMCAL acceptances. First the leading 
//  particle and then the jet have to fullfill several conditions 
//  (energy, direction ..) to be accepted. Then the fragmentation function 
//  of this jet is constructed   
// 
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, CERN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

#include <TFile.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TPaveLabel.h>
#include <TPad.h>
#include <TH2.h>

#include "AliPHOSGammaJet.h" 
#include "AliPHOSGetter.h" 
#include "AliPHOSGeometry.h"
#include "AliPHOSFastGlobalReconstruction.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDCaloCluster.h"
#include "Riostream.h"


ClassImp(AliPHOSGammaJet)

//____________________________________________________________________________
AliPHOSGammaJet::AliPHOSGammaJet() : 
  TTask(), fAnyConeOrPt(0), fOptionGJ(""),
  fOutputFile(new TFile(gDirectory->GetName())),
  fOutputFileName(gDirectory->GetName()),
  fInputFileName(gDirectory->GetName()),
  fHIJINGFileName(gDirectory->GetName()),
  fHIJING(0), fESDdata(0), fEtaCut(0.),
  fOnlyCharged(0), fPhiMaxCut(0.),
  fPhiMinCut(0.), fPtCut(0.),
  fNeutralPtCut(0.), fChargedPtCut(0.),
  fInvMassMaxCut(0.), fInvMassMinCut(0.),
  fMinDistance(0.), fRatioMaxCut(0.), fRatioMinCut(0.),
  fTPCCutsLikeEMCAL(0), fDirName(""), fESDTree(""),
  fPattern(""), fJetTPCRatioMaxCut(0.),
  fJetTPCRatioMinCut(0.), fJetRatioMaxCut(0.),
  fJetRatioMinCut(0.), fNEvent(0), fNCone(0),
  fNPt(0), fCone(0), fPtThreshold(0),
  fPtJetSelectionCut(0.0),
  fListHistos(new TObjArray(100)),
  fFastRec(0), fOptFast(0),
  fRan(0), fResPara1(0.), fResPara2(0.), fResPara3(0.),  
  fPosParaA(0.), fPosParaB(0.), fAngleMaxParam(), fSelect(0)
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
    fListHistos->Add(h) ; 
  }
  List() ; 
}

//____________________________________________________________________________
AliPHOSGammaJet::AliPHOSGammaJet(const TString inputfilename) : 
  TTask("GammaJet","Analysis of gamma-jet correlations"),
  fAnyConeOrPt(0), fOptionGJ(),
  fOutputFile(0),
  fOutputFileName(),
  fInputFileName(),
  fHIJINGFileName(),
  fHIJING(0), fESDdata(0), fEtaCut(0.),
  fOnlyCharged(0), fPhiMaxCut(0.),
  fPhiMinCut(0.), fPtCut(0.),
  fNeutralPtCut(0.), fChargedPtCut(0.),
  fInvMassMaxCut(0.), fInvMassMinCut(0.),
  fMinDistance(0.), fRatioMaxCut(0.), fRatioMinCut(0.),
  fTPCCutsLikeEMCAL(0), fDirName(), fESDTree(),
  fPattern(), fJetTPCRatioMaxCut(0.),
  fJetTPCRatioMinCut(0.), fJetRatioMaxCut(0.),
  fJetRatioMinCut(0.), fNEvent(0), fNCone(0),
  fNPt(0), fCone(0), fPtThreshold(0),
  fPtJetSelectionCut(0.0),
  fListHistos(0),
  fFastRec(0), fOptFast(0),
  fRan(0), fResPara1(0.), fResPara2(0.), fResPara3(0.),  
  fPosParaA(0.), fPosParaB(0.), fAngleMaxParam(), fSelect(0)

{
  // ctor
  fInputFileName = inputfilename;
  fFastRec = new AliPHOSFastGlobalReconstruction(fInputFileName);  
  AliPHOSGetter *  gime = AliPHOSGetter::Instance(fInputFileName) ;
  fNEvent = gime->MaxEvent();
  InitParameters();
}

//____________________________________________________________________________
AliPHOSGammaJet::AliPHOSGammaJet(const AliPHOSGammaJet & gj) : 
  TTask(gj),
  fAnyConeOrPt(gj.fAnyConeOrPt), fOptionGJ(gj.fOptionGJ),
  fOutputFile(gj.fOutputFile),
  fOutputFileName(gj.fOutputFileName),
  fInputFileName(gj.fInputFileName),
  fHIJINGFileName(gj.fHIJINGFileName),
  fHIJING(gj.fHIJING), fESDdata(gj.fESDdata), fEtaCut(gj.fEtaCut),
  fOnlyCharged(gj.fOnlyCharged), fPhiMaxCut(gj.fPhiMaxCut),
  fPhiMinCut(gj.fPhiMinCut), fPtCut(gj.fPtCut),
  fNeutralPtCut(gj.fNeutralPtCut), fChargedPtCut(gj.fChargedPtCut),
  fInvMassMaxCut(gj.fInvMassMaxCut), fInvMassMinCut(gj.fInvMassMinCut),
  fMinDistance(gj.fMinDistance), fRatioMaxCut(gj.fRatioMaxCut), 
  fRatioMinCut(gj.fRatioMinCut), fTPCCutsLikeEMCAL(gj.fTPCCutsLikeEMCAL), 
  fDirName(gj.fDirName), fESDTree(gj.fESDTree),
  fPattern(gj.fPattern), fJetTPCRatioMaxCut(gj.fJetTPCRatioMaxCut),
  fJetTPCRatioMinCut(gj.fJetTPCRatioMinCut), fJetRatioMaxCut(gj.fJetRatioMaxCut),
  fJetRatioMinCut(gj.fJetRatioMinCut), fNEvent(gj.fNEvent), fNCone(gj.fNCone),
  fNPt(gj.fNPt), fCone(gj.fCone), fPtThreshold(gj.fPtThreshold),
  fPtJetSelectionCut(gj.fPtJetSelectionCut),
  fListHistos(0),//?????
  fFastRec(gj.fFastRec), fOptFast(gj.fOptFast),
  fRan(0), //???
  fResPara1(gj.fResPara1), fResPara2(gj.fResPara2), fResPara3(gj.fResPara3),  
  fPosParaA(gj.fPosParaA), fPosParaB(gj.fPosParaB), 
  fAngleMaxParam(gj.fAngleMaxParam), fSelect(gj.fSelect)
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
AliPHOSGammaJet::~AliPHOSGammaJet() 
{
  fOutputFile->Close() ;
}

//____________________________________________________________________________
void AliPHOSGammaJet::AddHIJINGToList(Int_t iEvent, TClonesArray * particleList, 
				      TClonesArray * plCh, 
				      TClonesArray * plNe,  
				      TClonesArray * plNePHOS)
{

  // List of particles copied to a root file.

//   Char_t sb[] = "bgrd/";
//   //  cout<<sb<<endl;
//   Char_t si[10];
//   Char_t sf[] = "/list.root";
//   //  cout<<sf<<endl;
//   sprintf(si,"%d",iEvent);
//   strcat(sb,si) ;
//   //  cout<<si<<endl;
//   strcat(sb,sf) ;
//   //  cout<<si<<endl;
//   TFile * f = TFile::Open(sb) ;
//   //cout<<f->GetName()<<endl;

  Char_t fi[100];
  sprintf(fi,"bgrd/%d/list.root",iEvent);
  TFile * f = TFile::Open(fi) ;
  //cout<<f->GetName()<<endl;

  TParticle *particle = new TParticle();
  TTree *t = (TTree*) f->Get("T");
  TBranch *branch = t->GetBranch("primaries");
  branch->SetAddress(&particle);

  Int_t index   = particleList->GetEntries() ; 
  Int_t indexNe = plNe->GetEntries() ;
  Int_t indexCh = plCh->GetEntries() ;
  Int_t indexNePHOS = plNePHOS->GetEntries() ;
  Double_t charge = 0.;
  Int_t iParticle = 0 ; 
  Int_t m = 0;
  Double_t x = 0., z = 0.;
  //  cout<<"bkg entries "<<t->GetEntries()<<endl;

  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  const AliPHOSGeometry * geom = gime->PHOSGeometry() ;

  if(!fOptFast){
    for (iParticle=0 ; iParticle < t->GetEntries() ; iParticle++) {
      t->GetEvent(iParticle) ;
 
      m = 0 ;
      x = 0. ;
      z = 0. ;

      charge = TDatabasePDG::Instance()
	->GetParticle(particle->GetPdgCode())->Charge();
     
      if((charge != 0) && (particle->Pt() > fChargedPtCut)){
	if(TMath::Abs(particle->Eta())<fEtaCut){  
	  new((*particleList)[index]) TParticle(*particle) ;
	  (dynamic_cast<TParticle*>(particleList->At(index)))
	    ->SetStatusCode(0) ;
	  index++ ; 
	  
	  new((*plCh)[indexCh])       TParticle(*particle) ;
	  (dynamic_cast<TParticle*>(plCh->At(indexCh)))->SetStatusCode(0) ;
	  indexCh++ ;
	  
	  if(strstr(fOptionGJ,"deb all")||strstr(fOptionGJ,"deb")){
	    dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
	      Fill(particle->Pt());
	  }
	}
	else if((charge == 0) && (particle->Pt() > fNeutralPtCut) ){
	  geom->ImpactOnEmc(particle->Theta(),particle->Phi(), m,z,x);
	 
	  if(m != 0)
	    {//Is in PHOS
	      if(strstr(fOptionGJ,"deb all")|| strstr(fOptionGJ,"deb"))
		dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
		  Fill(particle->Pt());
	
	      new((*plNePHOS)[indexNePHOS])       TParticle(*particle) ;
	      (dynamic_cast<TParticle*>(plNePHOS->At(indexNePHOS)))->SetStatusCode(0) ;
	      indexNePHOS++ ; 
	    }
	  
	  if((particle->Phi()>fPhiEMCALCut[0]) && 
	     (particle->Phi()<fPhiEMCALCut[1]) && m == 0)
	    {//Is in EMCAL 
	      if(strstr(fOptionGJ,"deb all")|| strstr(fOptionGJ,"deb"))
		dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
		  Fill(particle->Pt());

	      new((*particleList)[index]) TParticle(*particle) ;
	      (dynamic_cast<TParticle*>(particleList->At(index)))
		->SetStatusCode(0) ;
	      index++ ; 
	      new((*plNe)[indexNe])       TParticle(*particle) ;
	      (dynamic_cast<TParticle*>(plNe->At(indexNe)))->SetStatusCode(0) ;
	      indexNe++ ; 
	    }
	}
      }
    }
  } //No OptFast
  else{
    
    Double_t mass = TDatabasePDG::Instance()->GetParticle(111)->Mass(); 
    TLorentzVector pPi0, pGamma1, pGamma2 ;
    Double_t angle = 0, cellDistance = 0.;
    Bool_t p1 = kFALSE;

    //     fFastRec = new AliPHOSFastGlobalReconstruction(fHIJINGFileName);
    //     fFastRec->FastReconstruction(iiEvent);
    
    for (iParticle=0 ; iParticle <   t->GetEntries() ; iParticle++) {
      t->GetEvent(iParticle) ;
      m = 0 ;
      x = 0. ;
      z = 0. ;
      charge = TDatabasePDG::Instance()
	->GetParticle(particle->GetPdgCode())->Charge();
      
      if((charge != 0) && (particle->Pt() > fChargedPtCut)
	 && (TMath::Abs(particle->Eta())<fEtaCut)){
	
	new((*particleList)[index]) TParticle(*particle) ;
	(dynamic_cast<TParticle*>(particleList->At(index)))
	  ->SetStatusCode(0) ;
	index++ ; 
	
	new((*plCh)[indexCh])       TParticle(*particle) ;
	(dynamic_cast<TParticle*>(plCh->At(indexCh)))->SetStatusCode(0) ;
	indexCh++ ;
      }
      else if(charge == 0){
	geom->ImpactOnEmc(particle->Theta(),particle->Phi(), m,z,x);
	if((particle->GetPdgCode() != 111) && (particle->Pt() > fNeutralPtCut) &&
	   (TMath::Abs(particle->Eta())<fEtaCut) ){
	  TLorentzVector part(particle->Px(),particle->Py(),
			      particle->Pz(),particle->Energy());
	  MakePhoton(part);
	  if(part.Pt() > fNeutralPtCut){

	    if(particle->Phi()>fPhiEMCALCut[0] && 
	       particle->Phi()<fPhiEMCALCut[1] && m == 0)
	      {
		new((*particleList)[index]) TParticle(*particle) ;
		(dynamic_cast<TParticle*>(particleList->At(index)))
		  ->SetStatusCode(0) ;
		(dynamic_cast<TParticle*>(particleList->At(index)))
		  ->SetMomentum(part.Px(),part.Py(),part.Pz(),part.E());
		index++ ; 	
    

		new((*plNe)[indexNe])       TParticle(*particle) ;
		(dynamic_cast<TParticle*>(plNe->At(indexNe)))
		  ->SetMomentum(part.Px(),part.Py(),part.Pz(),part.E());
		(dynamic_cast<TParticle*>(plNe->At(indexNe)))->SetStatusCode(0) ;
		indexNe++ ; 
	      }
	    if(m != 0)
	      {
		new((*plNePHOS)[indexNePHOS])       TParticle(*particle) ;
		(dynamic_cast<TParticle*>(plNePHOS->At(indexNePHOS)))
		  ->SetMomentum(part.Px(),part.Py(),part.Pz(),part.E());
		(dynamic_cast<TParticle*>(plNePHOS->At(indexNePHOS)))->SetStatusCode(0) ;
		indexNePHOS++ ; 
	      }
	  }
	}
	if((particle->GetPdgCode() == 111) && (particle->Pt() > fNeutralPtCut) && 
	   (TMath::Abs(particle->Eta())<fEtaCut+1))
	  {
	    
	    pPi0.SetPxPyPzE(particle->Px(),particle->Py(),particle->Pz(),
			    particle->Energy());
	    
	    //Decay
	    
	    Pi0Decay(mass,pPi0,pGamma1,pGamma2,angle);
	    //Check if decay photons are too close for PHOS
	    cellDistance = angle*460; //cm
	    if (cellDistance < fMinDistance) {
	      if(strstr(fOptionGJ,"deb all")|| strstr(fOptionGJ,"deb"))
		dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
		  Fill(particle->Pt());
	      
	      //Pi0 inside phi EMCAL acceptance
	      
	      TLorentzVector part(particle->Px(),particle->Py(),
				  particle->Pz(),particle->Energy());
	      MakePhoton(part);
	      if(part.Pt() > fNeutralPtCut){
		if(particle->Phi()>fPhiEMCALCut[0] && 
		   particle->Phi()<fPhiEMCALCut[1] && m == 0){
		  
		  new((*particleList)[index]) TParticle(*particle) ;
		  (dynamic_cast<TParticle*>(particleList->At(index)))->SetStatusCode(0) ;
		  (dynamic_cast<TParticle*>(particleList->At(index)))
		    ->SetMomentum(part.Px(),part.Py(),part.Pz(),part.E());
		  index++ ;
		  
		  new((*plNe)[indexNe])       TParticle(*particle) ;
		  (dynamic_cast<TParticle*>(plNe->At(indexNe)))      ->SetStatusCode(0) ;
		  (dynamic_cast<TParticle*>(plNe->At(indexNe)))
		    ->SetMomentum(part.Px(),part.Py(),part.Pz(),part.E());
		  indexNe++ ; 
		}
		if(m != 0){
		  if(strstr(fOptionGJ,"deb all")|| strstr(fOptionGJ,"deb"))
		    dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
		      Fill(particle->Pt());
		  new((*plNePHOS)[indexNePHOS])       TParticle(*particle) ;
		  (dynamic_cast<TParticle*>(plNePHOS->At(indexNePHOS))) ->SetStatusCode(0) ;
		  (dynamic_cast<TParticle*>(plNePHOS->At(indexNePHOS)))
		    ->SetMomentum(part.Px(),part.Py(),part.Pz(),part.E());
		  indexNePHOS++;
		}//In PHOS
	      }
	    }// if cell<distance	
	    else {
	      
	      dynamic_cast<TH2F*>(fListHistos->FindObject("AnglePair"))
		->Fill(pPi0.E(),angle);
	      
	      p1 = kFALSE;
	      if(pGamma1.Pt() > 0. && TMath::Abs(pGamma1.Eta())<fEtaCut){
		
		MakePhoton(pGamma1);
		
		if(pGamma1.Pt() > fNeutralPtCut ){
		  
		  TParticle * photon1 = 
		    new TParticle(22,1,0,0,0,0,pGamma1.Px(),pGamma1.Py(),
				  pGamma1.Pz(),pGamma1.E(),0,0,0,0);
		  geom->ImpactOnEmc(photon1->Theta(),photon1->Phi(), m,z,x);
		  if( photon1->Phi()>fPhiEMCALCut[0] && photon1->Phi()<fPhiEMCALCut[1]
		      && m == 0){
		    if(strstr(fOptionGJ,"deb all") || strstr(fOptionGJ,"deb"))
		      dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
			Fill(photon1->Pt());
		    new((*particleList)[index]) TParticle(*photon1) ;
		    (dynamic_cast<TParticle*>(particleList->At(index)))->SetStatusCode(0) ;
		    index++ ; 

		    new((*plNe)[indexNe])       TParticle(*photon1) ;
		    (dynamic_cast<TParticle*>(plNe->At(indexNe)))      ->SetStatusCode(0) ;
		    indexNe++ ; 
		    p1 = kTRUE;
		  }
		  if(m != 0){
		    if(strstr(fOptionGJ,"deb all") || strstr(fOptionGJ,"deb"))
		      dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
			Fill(photon1->Pt());
		    new((*plNePHOS)[indexNePHOS])       TParticle(*photon1) ;
		    (dynamic_cast<TParticle*>(plNePHOS->At(indexNePHOS)))->SetStatusCode(0) ;
		    indexNePHOS++;
		    p1 = kTRUE;
		  }//Is in PHOS
		}
	      }
	      if(pGamma2.Pt() > 0. && TMath::Abs(pGamma2.Eta())<fEtaCut){
		
		MakePhoton(pGamma2);
		if(pGamma2.Pt() > fNeutralPtCut){
		  
		  TParticle * photon2 =
		    new TParticle(22,1,0,0,0,0,pGamma2.Px(), pGamma2.Py(),
				  pGamma2.Pz(),pGamma2.E(),0,0,0,0);
		  geom->ImpactOnEmc(photon2->Theta(),photon2->Phi(), m,z,x);
		  if(photon2->Phi()>fPhiEMCALCut[0] && 
		     photon2->Phi()<fPhiEMCALCut[1] && m == 0){
		    if(strstr(fOptionGJ,"deb all") || strstr(fOptionGJ,"deb"))
		      dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
			Fill(photon2->Pt());
		    new((*particleList)[index]) TParticle(*photon2) ;
		    (dynamic_cast<TParticle*>(particleList->At(index)))->SetStatusCode(0) ;
		    index++ ; 
		    
		    new((*plNe)[indexNe])       TParticle(*photon2) ;
		    (dynamic_cast<TParticle*>(plNe->At(indexNe)))      ->SetStatusCode(0) ;
		    indexNe++ ; 
		  }
		  if(m != 0){
		    if(strstr(fOptionGJ,"deb all") || strstr(fOptionGJ,"deb"))
		      dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
			Fill(photon2->Pt());
		    new((*plNePHOS)[indexNePHOS])       TParticle(*photon2) ;
		    (dynamic_cast<TParticle*>(plNePHOS->At(indexNePHOS))) ->SetStatusCode(0) ;
		    indexNePHOS++;
		  }
		  if(p1){
		    //		    e = (pGamma1+pGamma2).E();
		    //		    if(IsAngleInWindow(angle,e))  	    
		      dynamic_cast<TH2F*>
			(fListHistos->FindObject("AnglePairAccepted"))->
			Fill(pPi0.E(),angle);
		  }
		}
	      }//photon2 in acceptance
	    }//if angle > mindist
	  }//if pi0
      }
    }//for (iParticle<nParticle)
  }
  
  //Info("AddHIJINGToList","End HIJING");
}  

//____________________________________________________________________________
Double_t AliPHOSGammaJet::CalculateJetRatioLimit(const Double_t ptg, 
						 const Double_t *par, 
						 const Double_t *x) {

  //Info("CalculateLimit","x1 %f, x2%f",x[0],x[1]);
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
void AliPHOSGammaJet::CreateParticleList(Int_t iEvent, 
					 TClonesArray * particleList, 
					 TClonesArray * plCh, 
					 TClonesArray * plNe,  
					 TClonesArray * plNePHOS)
{
  //Info("CreateParticleList","Inside");
  AliPHOSGetter * gime = AliPHOSGetter::Instance(fInputFileName) ;
  const AliPHOSGeometry * geom = gime->PHOSGeometry() ; 
  gime->Event(iEvent, "X") ;


  Int_t index = particleList->GetEntries() ; 
  Int_t indexCh     = plCh->GetEntries() ;
  Int_t indexNe     = plNe->GetEntries() ;
  Int_t indexNePHOS = plNePHOS->GetEntries() ;
  Int_t iParticle = 0 ;
  Double_t charge = 0.;
  Int_t m = 0;
  Double_t x = 0., z = 0.;
  if(!fOptFast){
    
    for (iParticle=0 ; iParticle <  gime->NPrimaries() ; iParticle++) {
      const TParticle * particle = gime->Primary(iParticle) ;
      
      m = 0 ;
      x = 0. ;
      z = 0. ;
      
      //Keep Stable particles within eta range 
      if((particle->GetStatusCode() == 1) &&
	 (particle->Pt() > 0)){
	if(TMath::Abs(particle->Eta())<fEtaCut){  
	  //Fill lists
	  
	  charge = TDatabasePDG::Instance()
	    ->GetParticle(particle->GetPdgCode())->Charge();
	  if((charge != 0) && (particle->Pt() > fChargedPtCut)){
	    
	    if(strstr(fOptionGJ,"deb all")|| strstr(fOptionGJ,"deb"))
	      dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
		Fill(particle->Pt());
	    new((*plCh)[indexCh++])       TParticle(*particle) ;
	    new((*particleList)[index++]) TParticle(*particle) ;
	  }
	  else if((charge == 0) && (particle->Pt() > fNeutralPtCut)){
	    geom->ImpactOnEmc(particle->Theta(),particle->Phi(), m,z,x);
	    if(m != 0)
	      {//Is in PHOS
		if(strstr(fOptionGJ,"deb all")|| strstr(fOptionGJ,"deb"))
		  dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
		    Fill(particle->Pt());
		
		new((*plNePHOS)[indexNePHOS++])       TParticle(*particle) ;
	      }
	    if((particle->Phi()>fPhiEMCALCut[0]) && 
	       (particle->Phi()<fPhiEMCALCut[1]) && m == 0)
	      {//Is in EMCAL 
		if(strstr(fOptionGJ,"deb all")|| strstr(fOptionGJ,"deb"))
		  dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
		    Fill(particle->Pt());
		new((*plNe)[indexNe++])       TParticle(*particle) ;
		new((*particleList)[index++]) TParticle(*particle) ;
	      }
	  }
	}
      }//final particle, etacut
    }//for (iParticle<nParticle)
  }// No Op
  else
    {
      Double_t mass = TDatabasePDG::Instance()->GetParticle(111)->Mass(); 
      TLorentzVector pPi0, pGamma1, pGamma2 ;
      Double_t angle = 0, cellDistance = 0.;
      
      fFastRec = new AliPHOSFastGlobalReconstruction(fInputFileName);
      fFastRec->FastReconstruction(iEvent);
      
      for (iParticle=0 ; iParticle <  gime->NPrimaries() ; iParticle++) {
	const TParticle * particle = gime->Primary(iParticle) ;
	m = 0 ;
	x = 0. ;
	z = 0. ;
	//Keep Stable particles within eta range 
	if((particle->GetStatusCode() == 1) && (particle->Pt() > 0)){
	  
	  //Fill lists
	  
	  charge = TDatabasePDG::Instance()
	    ->GetParticle(particle->GetPdgCode())->Charge();
	  if((charge != 0) && (particle->Pt() > fChargedPtCut) && (TMath::Abs(particle->Eta())<fEtaCut)){
	    new((*plCh)[indexCh++])       TParticle(*particle) ;
	    new((*particleList)[index++]) TParticle(*particle) ;
	  }
	  else if(charge == 0) {
	    geom->ImpactOnEmc(particle->Theta(),particle->Phi(), m,z,x);
	    if((particle->GetPdgCode() != 111) && particle->Pt() > 0 &&
	       (TMath::Abs(particle->Eta())<fEtaCut))
	    {                
	      
	      TLorentzVector part(particle->Px(),particle->Py(),
				  particle->Pz(),particle->Energy());
	      
	      MakePhoton(part);
	      
	      if(part.Pt() > fNeutralPtCut){
		if(particle->Phi()>fPhiEMCALCut[0] && 
		   particle->Phi()<fPhiEMCALCut[1] && m == 0)
		  {
		    new((*particleList)[index]) TParticle(*particle) ;
		    (dynamic_cast<TParticle*>(particleList->At(index)))
		      ->SetMomentum(part.Px(),part.Py(),part.Pz(),part.E());
		    index++ ; 

		    new((*plNe)[indexNe])       TParticle(*particle) ;
		    (dynamic_cast<TParticle*>(plNe->At(indexNe)))
		      ->SetMomentum(part.Px(),part.Py(),part.Pz(),part.E());
		    indexNe++ ; 
		  }
		if(m != 0)
		  {
		    new((*plNePHOS)[indexNePHOS])       TParticle(*particle) ;
		    (dynamic_cast<TParticle*>(plNePHOS->At(indexNePHOS)))
		      ->SetMomentum(part.Px(),part.Py(),part.Pz(),part.E());
		    indexNePHOS++ ; 
		  }
 	      }// Small pt
	    } //No Pi0
	    if((particle->GetPdgCode() == 111) && (particle->Pt() > fNeutralPtCut) && 
	       (TMath::Abs(particle->Eta())<fEtaCut+1))
	      {
	      
		pPi0.SetPxPyPzE(particle->Px(),particle->Py(),particle->Pz(),
				particle->Energy());
		
		//Decay
		
		Pi0Decay(mass,pPi0,pGamma1,pGamma2,angle);
		//Check if decay photons are too close for PHOS
		cellDistance = angle*460; //cm
		
		if (cellDistance < fMinDistance) {
		  
		  //Pi0 inside phi EMCAL acceptance
	
		    
		    TLorentzVector part(particle->Px(),particle->Py(),
					particle->Pz(),particle->Energy());
		    MakePhoton(part);
		    
		    if(part.Pt() > fNeutralPtCut){
		      if(particle->Phi()>fPhiEMCALCut[0] && 
			 particle->Phi()<fPhiEMCALCut[1] && m == 0){
			if(strstr(fOptionGJ,"deb all")|| strstr(fOptionGJ,"deb"))
			  dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
			    Fill(particle->Pt());
			
			new((*plNe)[indexNe])       TParticle(*particle) ;
			(dynamic_cast<TParticle*>(plNe->At(indexNe)))
			  ->SetMomentum(part.Px(),part.Py(),part.Pz(),part.E());
			new((*particleList)[index]) TParticle(*particle) ;
			(dynamic_cast<TParticle*>(particleList->At(index)))
			  ->SetMomentum(part.Px(),part.Py(),part.Pz(),part.E());
			index++;
			indexNe++;
		      }//InEMCAL
		      if(m != 0){
			if(strstr(fOptionGJ,"deb all")|| strstr(fOptionGJ,"deb"))
			  dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
			    Fill(particle->Pt());
			new((*plNePHOS)[indexNePHOS])       TParticle(*particle) ;
			(dynamic_cast<TParticle*>(plNePHOS->At(indexNePHOS)))
			  ->SetMomentum(part.Px(),part.Py(),part.Pz(),part.E());
			indexNePHOS++;
		      }//In PHOS
		    }//Small Pt
		}// if cell<distance	
		else {
		  
		  dynamic_cast<TH2F*>(fListHistos->FindObject("AnglePair"))
		    ->Fill(pPi0.E(),angle);
		  
		  Bool_t p1 = kFALSE;
		  
		  if(pGamma1.Pt() > 0 && TMath::Abs(pGamma1.Eta())<fEtaCut){
		    MakePhoton(pGamma1);
		    
		    if(pGamma1.Pt() > fNeutralPtCut){
		      TParticle * photon1 = 
			new TParticle(22,1,0,0,0,0,pGamma1.Px(),pGamma1.Py(),
				      pGamma1.Pz(),pGamma1.E(),0,0,0,0);
		      geom->ImpactOnEmc(photon1->Theta(),photon1->Phi(), m,z,x);
		      if( photon1->Phi()>fPhiEMCALCut[0] && photon1->Phi()<fPhiEMCALCut[1]
			  && m == 0){
		      if(strstr(fOptionGJ,"deb all") || strstr(fOptionGJ,"deb"))
			dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
			  Fill(photon1->Pt());
		      new((*plNe)[indexNe++])       TParticle(*photon1) ;
		      new((*particleList)[index++]) TParticle(*photon1) ;
		      //photon1->Print();
		      p1 = kTRUE;
		      }
		      if(m != 0){
			if(strstr(fOptionGJ,"deb all") || strstr(fOptionGJ,"deb"))
			  dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
			    Fill(photon1->Pt());
			new((*plNePHOS)[indexNePHOS++])       TParticle(*photon1) ;
		      p1 = kTRUE;
		      }
		    }
		  }
		  if(pGamma2.Pt() > 0 && TMath::Abs(pGamma2.Eta())<fEtaCut){
		    MakePhoton(pGamma2);
		    
		    if(pGamma2.Pt() > fNeutralPtCut ){
		      TParticle * photon2 =
			new TParticle(22,1,0,0,0,0,pGamma2.Px(), pGamma2.Py(),
				      pGamma2.Pz(),pGamma2.E(),0,0,0,0);
		      geom->ImpactOnEmc(photon2->Theta(),photon2->Phi(), m,z,x);
		      if(photon2->Phi()>fPhiEMCALCut[0] && 
			 photon2->Phi()<fPhiEMCALCut[1] && m == 0){
			if(strstr(fOptionGJ,"deb all") || strstr(fOptionGJ,"deb"))
			  dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
			    Fill(photon2->Pt());
			new((*plNe)[indexNe++])       TParticle(*photon2) ;
			new((*particleList)[index++]) TParticle(*photon2) ;
		      }
		      if(m != 0){
			if(strstr(fOptionGJ,"deb all") || strstr(fOptionGJ,"deb"))
			  dynamic_cast<TH1F*>(fListHistos->FindObject("PtSpectra"))->
			    Fill(photon2->Pt());
			new((*plNePHOS)[indexNePHOS++])       TParticle(*photon2) ;
		      }
		      
		      if(p1){
			// Float_t e = (pGamma1+pGamma2).E();
			// if(IsAngleInWindow(angle,e))		    
			  dynamic_cast<TH2F*>
			    (fListHistos->FindObject("AnglePairAccepted"))->
			    Fill(pPi0.E(),angle);
		      }
		    }
		  }//photon2 in acceptance
		}//if angle > mindist
	      }//if pi0
	  }//If neutral
	}//final particle, etacut
      }//for (iParticle<nParticle)
    }//OptFast
  //gime->Delete() ; 
}
//____________________________________________________________________________
void AliPHOSGammaJet::CreateParticleListFromESD(TClonesArray * pl, 
						TClonesArray * plCh, 
						TClonesArray * plNe,  
						TClonesArray * plNePHOS,
						const AliESD * esd){
  
  //Create a list of particles from the ESD. These particles have been measured 
  //by the Central Tracking system (TPC), PHOS and EMCAL 
  //(EMCAL not available for the moment). 
  //Info("CreateParticleListFromESD","Inside");

  Int_t index = pl->GetEntries() ; 
  Int_t npar  = 0 ;
  Float_t *pid = new Float_t[AliPID::kSPECIESN];  

  //########### PHOS ##############
  //Info("CreateParticleListFromESD","Fill ESD PHOS list");
  Int_t begphos = esd->GetFirstPHOSCluster();  
  Int_t endphos = esd->GetFirstPHOSCluster() + 
    esd->GetNumberOfPHOSClusters() ;  
  Int_t indexNePHOS = plNePHOS->GetEntries() ;
  if(strstr(fOptionGJ,"deb all"))
    Info("CreateParticleListFromESD","PHOS: first particle %d, last particle %d",
	 begphos,endphos);

  for (npar =  begphos; npar <  endphos; npar++) {//////////////PHOS track loop
      AliESDCaloCluster * clus = esd->GetCaloCluster(npar) ; // retrieve track from esd
   
    //Create a TParticle to fill the particle list

    Float_t en = clus->GetClusterEnergy() ;
    Float_t *p = new Float_t();
    clus->GetGlobalPosition(p) ;
    TVector3 pos(p[0],p[1],p[2]) ; 
    Double_t phi  = pos.Phi();
    Double_t theta= pos.Theta();
    Double_t px = en*TMath::Cos(phi)*TMath::Sin(theta);;
    Double_t py = en*TMath::Sin(phi)*TMath::Sin(theta);
    Double_t pz = en*TMath::Cos(theta);

    TParticle * particle = new TParticle() ;
    particle->SetMomentum(px,py,pz,en) ;

    //Select only photons
    
    pid=clus->GetPid();
    //cout<<"pid "<<pid[AliPID::kPhoton]<<endl ;
    if( pid[AliPID::kPhoton] > 0.75)
      new((*plNePHOS)[indexNePHOS++])   TParticle(*particle) ;
  }

  //########### TPC #####################
  //Info("CreateParticleListFromESD","Fill ESD TPC list");
  Int_t begtpc = 0 ;  
  Int_t endtpc = esd->GetNumberOfTracks() ;
  Int_t indexCh     = plCh->GetEntries() ;
  if(strstr(fOptionGJ,"deb all"))
    Info("CreateParticleListFromESD","TPC: first particle %d, last particle %d",
	 begtpc,endtpc);

  for (npar =  begtpc; npar <  endtpc; npar++) {////////////// track loop
    AliESDtrack * track = esd->GetTrack(npar) ; // retrieve track from esd
    
    Double_t en = track ->GetTPCsignal() ;
    Double_t mom[3];
    track->GetPxPyPz(mom) ;
    Double_t px = mom[0];
    Double_t py = mom[1];
    Double_t pz = mom[2]; //Check with TPC people if this is correct.

    //cout<<"TPC signal "<<en<<endl;
    //cout<<"px "<<px<<"; py "<<py<<"; pz "<<pz<<endl;
    TParticle * particle = new TParticle() ;
    particle->SetMomentum(px,py,pz,en) ;

    new((*plCh)[indexCh++])       TParticle(*particle) ;    
    new((*pl)[index++])           TParticle(*particle) ;

  }

  //################ EMCAL ##############
  Double_t v[3] ; //vertex ;
  esd->GetVertex()->GetXYZ(v) ; 
  //##########Uncomment when ESD for EMCAL works ##########  
  //Info("CreateParticleListFromESD","Fill ESD EMCAL list");
 
  Int_t begem = esd->GetFirstEMCALCluster();  
  Int_t endem = esd->GetFirstEMCALCluster() + 
    esd->GetNumberOfEMCALClusters() ;  
  Int_t indexNe  = plNe->GetEntries() ; 
  if(strstr(fOptionGJ,"deb all"))
    Info("CreateParticleListFromESD","EMCAL: first particle %d, last particle %d",
	 begem,endem);
   
  for (npar =  begem; npar <  endem; npar++) {//////////////EMCAL track loop
     AliESDCaloCluster * clus = esd->GetCaloCluster(npar) ; // retrieve track from esd
  
    Float_t en = clus->GetClusterEnergy() ;
    Float_t *p = new Float_t();
    clus->GetGlobalPosition(p) ;
    TVector3 pos(p[0],p[1],p[2]) ;
    Double_t phi  = pos.Phi();
    Double_t theta= pos.Theta();
    Double_t px = en*TMath::Cos(phi)*TMath::Sin(theta);;
    Double_t py = en*TMath::Sin(phi)*TMath::Sin(theta);
    Double_t pz = en*TMath::Cos(theta);
    //cout<<"EMCAL signal "<<en<<endl;
    //cout<<"px "<<px<<"; py "<<py<<"; pz "<<pz<<endl;
    //TParticle * particle = new TParticle() ;
    //particle->SetMomentum(px,py,pz,en) ;
  
    Int_t pdg = 0;
    //  //Uncomment if PID IS WORKING, photon and pi0 idenitification.
    //  //if( pid[AliPID::kPhoton] > 0.75) //This has to be fixen.
    //  //pdg = 22;
    //  //else if( pid[AliPID::kPi0] > 0.75)
    //  //pdg = 111;
    pdg = 22; //No PID, assume all photons
    TParticle * particle = new TParticle(pdg, 1, -1, -1, -1, -1, 
					 px, py, pz, en, v[0], v[1], v[2], 0);

    new((*plNe)[indexNe++])       TParticle(*particle) ; 
    new((*pl)[index++])           TParticle(*particle) ;
  }
  
  //  Info("CreateParticleListFromESD","End Inside");
  
}



//____________________________________________________________________________
void AliPHOSGammaJet::Exec(Option_t *option) 
{
  // does the job
  fOptionGJ = option;
  MakeHistos() ; 
  
  AliESD * esd = 0;
  TChain * t = 0 ;

  if(fESDdata){
    // Create chain of esd trees
    const UInt_t kNevent = static_cast<UInt_t>(GetNEvent()) ; 
    t = ReadESD(kNevent, fDirName, fESDTree, fPattern) ; 
    if(!t) {
      AliError("Could not create the TChain") ; 
      //break ;
    }
    
    // ESD  
    t->SetBranchAddress("ESD",&esd); // point to the container esd where to put the event from the esdTree  

  }
 

//   AliGenPythia* pyth = (AliGenPythia*) gAlice->Generator();
//   pyth->Init();

  TClonesArray * particleList = new TClonesArray("TParticle",1000);
  TClonesArray * plCh         = new TClonesArray("TParticle",1000);
  TClonesArray * plNe         = new TClonesArray("TParticle",1000);
  TClonesArray * plNePHOS     = new TClonesArray("TParticle",1000);

  for (Int_t iEvent = 0 ; iEvent < fNEvent ; iEvent++) {
    if(strstr(fOptionGJ,"deb")||strstr(fOptionGJ,"deb all"))
      Info("Exec", "Event %d", iEvent) ;

    fRan.SetSeed(0);

    Double_t phig = 0., phil = 0., phich = 0 , phipi = 0;
    Double_t etag = 0., etal = 0., etach = 0., etapi = 0. ;  
    Double_t ptg  = 0., ptl  = 0., ptch  = 0., ptpi  = 0.;
 
    TLorentzVector jet   (0,0,0,0);
    TLorentzVector jettpc(0,0,0,0);

    if(fESDdata){

      Int_t iNbytes = t->GetEntry(iEvent); // store event in esd
      //cout<<"nbytes "<<iNbytes<<endl;
      if ( iNbytes == 0 ) {
	AliError("Empty TChain") ; 
	break ;
      }      
      CreateParticleListFromESD(particleList, plCh,plNe,plNePHOS, esd); //,iEvent);
    }
    else{   
      CreateParticleList(iEvent, particleList, plCh,plNe,plNePHOS); 
      
      //    TLorentzVector pyjet(0,0,0,0);
      
      //     Int_t nJ, nJT;
      //     Float_t jets[4][10];
      //     pyth->SetJetReconstructionMode(1);
      //     pyth->LoadEvent();
      //     pyth->GetJets(nJ, nJT, jets);
      
      //     Float_t pxJ = jets[0][0];
      //     Float_t pyJ = jets[1][0];
      //     Float_t pzJ = jets[2][0];
      //     Float_t eJ  = jets[3][0];
      //     pyjet.SetPxPyPzE(pxJ,pyJ,pzJ,eJ ) ;
      
      //     if(nJT > 1){
      //       //Info("Exec",">>>>>>>>>>Number of jets !!!!   %d",nJT);
      //       for (Int_t iJ = 1; iJ < nJT; iJ++) {
      // 	Float_t pxJ = jets[0][iJ];
      // 	Float_t pyJ = jets[1][iJ];
      // 	Float_t pzJ = jets[2][iJ];
      // 	Float_t eJ  = jets[3][iJ];
      // 	pyjet.SetPxPyPzE(pxJ,pyJ,pzJ,eJ ) ;
      // 	//Info("Exec",">>>>>Pythia Jet: %d, Phi %f, Eta %f, Pt %f",
      // 	//	     iJ,pyjet.Phi(),pyjet.Eta(),pyjet.Pt());
      //       }
      
      //     }
      
      if(fHIJING)
	AddHIJINGToList(iEvent, particleList, plCh,plNe, plNePHOS);
    }
    
    Bool_t iIsInPHOS = kFALSE ;
    GetGammaJet(plNePHOS, ptg, phig, etag, iIsInPHOS) ; 

    if(iIsInPHOS){

      //Info("Exec"," In PHOS") ;
      dynamic_cast<TH1F*>(fListHistos->FindObject("NGamma"))->Fill(ptg);
      dynamic_cast<TH2F*>(fListHistos->FindObject("PhiGamma"))
 	->Fill(ptg,phig);
      dynamic_cast<TH2F*>(fListHistos->FindObject("EtaGamma"))
 	->Fill(ptg,etag);
      if(strstr(fOptionGJ,"deb")||strstr(fOptionGJ,"deb all"))
 	Info("Exec", "Gamma: pt %f, phi %f, eta %f", ptg,
 	     phig,etag) ;

//       cout<<"n charged "<<plCh->GetEntries()<<endl;
//       cout<<"n neutral "<<plNe->GetEntries()<<endl;
//       cout<<"n All     "<<particleList->GetEntries()<<endl;
      
      GetLeadingCharge(plCh, ptg, phig, ptch, etach, phich) ;
      GetLeadingPi0   (plNe, ptg, phig, ptpi, etapi, phipi) ;

//       cout<<"n2 charged "<<plCh->GetEntries()<<endl;
//       cout<<"n2 neutral "<<plNe->GetEntries()<<endl;
//       cout<<"n2 All     "<<particleList->GetEntries()<<endl;


      //TPC+EMCAL

      //Is the leading cone inside EMCAL?
      Bool_t insidech = kFALSE ;
      if((phich - fCone) >  fPhiEMCALCut[0] && 
	 (phich + fCone) <  fPhiEMCALCut[1]){
	insidech = kTRUE ;
      }
      Bool_t insidepi = kFALSE ;
      if((phipi - fCone) >  fPhiEMCALCut[0] && 
	 (phipi + fCone) <  fPhiEMCALCut[1]){
	insidepi = kTRUE ;
      }

      if ((ptch > 0 || ptpi > 0)){
	if((ptch > ptpi) && insidech){
	  phil = phich ;
	  etal = etach ;
	  ptl  = ptch ;
	  dynamic_cast<TH2F*>(fListHistos->FindObject("ChargeRatio"))
	    ->Fill(ptg,ptch/ptg);
	  dynamic_cast<TH2F*>(fListHistos->FindObject("DeltaPhiCharge"))
	    ->Fill(ptg,phig-phich);
	  dynamic_cast<TH2F*>(fListHistos->FindObject("DeltaEtaCharge"))
	    ->Fill(ptg,etag-etach);
	  if(strstr(fOptionGJ,"deb"))
	  Info("Exec"," Charged Leading") ;
	}
	if((ptpi > ptch) && insidepi){
	  phil = phipi ;
	  etal = etapi ;
	  ptl  = ptpi ;
	  
	  dynamic_cast<TH2F*>(fListHistos->FindObject("Pi0Ratio"))
	    ->Fill(ptg,ptpi/ptg);
	  dynamic_cast<TH2F*>(fListHistos->FindObject("DeltaPhiPi0"))
	    ->Fill(ptg,phig-phipi);
	  dynamic_cast<TH2F*>(fListHistos->FindObject("DeltaEtaPi0"))
	    ->Fill(ptg,etag-etapi);
	  
	  if(ptpi > 0. && strstr(fOptionGJ,"deb"))
	    Info("Exec"," Pi0 Leading") ;
	}
		
	if(strstr(fOptionGJ,"deb"))
	  Info("Exec","Leading pt %f, phi %f",ptl,phil);
	if(insidech || insidepi){
	  if(!fAnyConeOrPt){
	   
	    MakeJet(particleList, ptg, phig, ptl, phil, etal, "", jet);

	    if(strstr(fOptionGJ,"deb")){
// 	      Info("Exec","Pythia Jet: Phi %f, Eta %f, Pt %f",
// 		   pyjet.Phi(),pyjet.Eta(),pyjet.Pt());
	      Info("Exec","TPC+EMCAL Jet: Phi %f, Eta %f, Pt %f",
		   jet.Phi(),jet.Eta(),jet.Pt());
	    }
// 	    dynamic_cast<TH2F*>(fListHistos->FindObject("DeltaPhiJet"))
// 	      ->Fill(ptg,pyjet.Phi()-jet.Phi());
// 	    dynamic_cast<TH2F*>(fListHistos->FindObject("DeltaEtaJet"))
// 	      ->Fill(ptg,pyjet.Eta()-jet.Eta());
// 	    dynamic_cast<TH2F*>(fListHistos->FindObject("DeltaPtJet"))
// 	      ->Fill(ptg,pyjet.Pt()-jet.Pt());
	  }
	  else
	    MakeJetAnyConeOrPt(particleList, ptg, phig, ptl, phil, etal, "");
	}

	//TPC
	if(fOnlyCharged && ptch > 0.)
	  {
	    if(strstr(fOptionGJ,"deb"))
	      Info("Exec","Leading TPC pt %f, phi %f",ptch,phich);

	    dynamic_cast<TH2F*>(fListHistos->FindObject("TPCRatio"))
	      ->Fill(ptg,ptch/ptg);
	    dynamic_cast<TH2F*>(fListHistos->FindObject("DeltaPhiTPC"))
	      ->Fill(ptg,phig-phich);
	    dynamic_cast<TH2F*>(fListHistos->FindObject("DeltaEtaTPC"))
	      ->Fill(ptg,etag-etach);
	    
	    if(!fAnyConeOrPt){
	   
	      MakeJet(plCh, ptg, phig, ptch, phich, etach, "TPC",jettpc);

	      if(strstr(fOptionGJ,"deb")){
// 		Info("Exec","Pythia Jet: Phi %f, Eta %f, Pt %f",
// 		     pyjet.Phi(),pyjet.Eta(),pyjet.Pt());
		Info("Exec","TPC Jet: Phi %f, Eta %f, Pt %f",
		     jettpc.Phi(),jettpc.Eta(),jettpc.Pt());
	      }
// 	      dynamic_cast<TH2F*>(fListHistos->FindObject("DeltaPhiTPCJet"))
// 		->Fill(ptg,pyjet.Phi()-jettpc.Phi());
// 	      dynamic_cast<TH2F*>(fListHistos->FindObject("DeltaEtaTPCJet"))
// 		->Fill(ptg,pyjet.Eta()-jettpc.Eta());
// 	      dynamic_cast<TH2F*>(fListHistos->FindObject("DeltaPtTPCJet"))
// 		->Fill(ptg,pyjet.Pt()-jettpc.Pt());
	    }
	    else
	      MakeJetAnyConeOrPt(plCh, ptg, phig, ptch, phich, etach, "TPC");
	    
	  }
      }
    }
    
    particleList->Delete() ; 
    plCh->Delete() ;
    plNe->Delete() ;
    plNePHOS->Delete() ;
  }//loop: events
  
  delete plNe ;
  delete plCh ;
  delete particleList ;

  fOutputFile->Write() ; 
  fOutputFile->cd();
  this->Write();
}    

//____________________________________________________________________________
void AliPHOSGammaJet::FillJetHistos(TClonesArray * pl, Double_t ptg, 
				    TString conf, TString type)
{
  //Fill jet fragmentation histograms if !fAnyCone, 
  //only for fCone and fPtThres 
  TParticle * particle = 0 ;
  Int_t ipr = -1 ;
  Float_t  charge = 0;
  
  TIter next(pl) ; 
  while ( (particle = dynamic_cast<TParticle*>(next())) ) {
    ipr++ ;
    Double_t pt = particle->Pt();
    
    charge = TDatabasePDG::Instance()
      ->GetParticle(particle->GetPdgCode())->Charge();
    if(charge != 0){//Only jet Charged particles 
      dynamic_cast<TH2F*>
	(fListHistos->FindObject(type+conf+"Fragment"))
	->Fill(ptg,pt/ptg);
      dynamic_cast<TH2F*>
	(fListHistos->FindObject(type+conf+"PtDist"))
	->Fill(ptg,pt);
    }
  }
  if(type == "Bkg"){
    dynamic_cast<TH1F*>
      (fListHistos->FindObject("NBkg"+conf))->Fill(ipr);
  }
}
//____________________________________________________________________________
void AliPHOSGammaJet::FillJetHistosAnyConeOrPt(TClonesArray * pl, Double_t ptg, 
					       TString conf, TString type,
					       TString cone, TString ptcut)
{
   //Fill jet fragmentation histograms if fAnyCone, 
   //for several cones and pt thresholds  
   TParticle *particle = 0;
   Int_t ipr=-1;
   Float_t  charge = 0;
  
   TIter next(pl) ; 
   while ( (particle = dynamic_cast<TParticle*>(next())) ) {
     ipr++;  
     Double_t pt = particle->Pt();
     charge = TDatabasePDG::Instance()
       ->GetParticle(particle->GetPdgCode())->Charge();
     if(charge != 0){//Only jet Charged particles
       dynamic_cast<TH2F*>
	 (fListHistos->FindObject(type+conf+"FragmentCone"+cone+"Pt"+ptcut))
	 ->Fill(ptg,pt/ptg);
       dynamic_cast<TH2F*>
	 (fListHistos->FindObject(type+conf+"PtDistCone"+cone+"Pt"+ptcut))
	 ->Fill(ptg,pt);  
     } 
   }//while
   
   if(type == "Bkg"){
     dynamic_cast<TH1F*>
       (fListHistos->FindObject("NBkg"+conf+"Cone"+cone+"Pt"+ptcut))
       ->Fill(ipr);
   }  
}

//____________________________________________________________________________
void AliPHOSGammaJet::GetGammaJet(TClonesArray * pl, Double_t &pt, 
				  Double_t &phi, Double_t &eta, Bool_t &Is) const 
{
  //Search for the prompt photon in PHOS with pt > fPtCut
  pt  = -10.;
  eta = -10.;
  phi = -10.;

  for(Int_t ipr = 0;ipr < pl->GetEntries() ; ipr ++ ){
    TParticle * particle = dynamic_cast<TParticle *>(pl->At(ipr)) ;

    if((particle->Pt() > fPtCut) && (particle->Pt() > pt)){

      pt  = particle->Pt();          
      phi = particle->Phi() ;
      eta = particle->Eta() ;
      Is  = kTRUE;
    }
  }
}

//____________________________________________________________________________
void  AliPHOSGammaJet::GetLeadingCharge(TClonesArray * pl, 
					Double_t ptg, Double_t phig, 
					Double_t &pt, Double_t &eta, Double_t &phi) const 
{  
  //Search for the charged particle with highest with 
  //Phi=Phi_gamma-Pi and pT=0.1E_gamma 
  pt  = -100.;
  eta = -100;
  phi = -100;

  for(Int_t ipr = 0;ipr < pl->GetEntries() ; ipr ++ ){

    TParticle * particle = dynamic_cast<TParticle *>(pl->At(ipr)) ;

    Double_t ptl  = particle->Pt();
    Double_t rat  = ptl/ptg ;
    Double_t phil = particle->Phi() ;
   
     if(((phig-phil)> fPhiMinCut) && ((phig-phil)<fPhiMaxCut) &&
        (rat > fRatioMinCut) && (rat < fRatioMaxCut)  && (ptl  > pt)) {
      eta = particle->Eta() ;
      phi = phil ;
      pt  = ptl ;
      //printf("GetLeadingCharge: %f %f %f %f \n", pt, eta, phi,rat) ;
    }
  }
  //printf("GetLeadingCharge: %f %f %f \n", pt, eta, phi) ; 

}


//____________________________________________________________________________
void  AliPHOSGammaJet::GetLeadingPi0(TClonesArray * pl, 
				     Double_t ptg, Double_t phig, 
				     Double_t &pt,  Double_t &eta, Double_t &phi)  
{  

  //Search for the neutral pion with highest with 
  //Phi=Phi_gamma-Pi and pT=0.1E_gamma 
  pt  = -100.;
  eta = -100.;
  phi = -100.;
  Double_t ptl = -100.;
  Double_t rat = -100.; 
  Double_t phil = -100. ;

  TIter next(pl);
  TParticle * particle = 0;
  Float_t ef = 0;
  if(!fOptFast){
    Float_t e = 0;
    while ( (particle = (TParticle*)next()) ) {  
      if( particle->GetPdgCode() == 111){
	ptl  = particle->Pt();
	rat  = ptl/ptg ;
	phil = particle->Phi() ;
	e    = particle->Energy();
	dynamic_cast<TH2F*>
	  (fListHistos->FindObject("AnglePairNoCut"))->
	  Fill(e,0.1);
	dynamic_cast<TH2F*>
	  (fListHistos->FindObject("InvMassPairNoCut"))->
	  Fill(ptg,1.35);

	if(((phig-phil)> fPhiMinCut) && ((phig-phil)<fPhiMaxCut) &&
	   (rat > fRatioMinCut) && (rat < fRatioMaxCut)) {

	  dynamic_cast<TH2F*>
	    (fListHistos->FindObject("AnglePairLeadingCut"))->
	    Fill(e,0.1);
	  dynamic_cast<TH2F*>
	    (fListHistos->FindObject("InvMassPairLeadingCut"))->
	    Fill(ptg,1.35);

	  dynamic_cast<TH2F*>
	    (fListHistos->FindObject("AnglePairAngleCut"))->
	    Fill(e,0.15);
	  dynamic_cast<TH2F*>
	    (fListHistos->FindObject("InvMassPairAngleCut"))->
	    Fill(ptg,1.36);
	  
	  dynamic_cast<TH2F*>
	    (fListHistos->FindObject("InvMassPairAllCut"))->
	    Fill(ptg,0.27);
	  dynamic_cast<TH2F*>
	    (fListHistos->FindObject("AnglePairAllCut"))->
	    Fill(e,1.34);


	  if(ptl  > pt){
	    eta = particle->Eta() ; 
	    phi = phil ;
	    pt  = ptl ;
	    ef  = e;
	    //printf("GetLeadingPi0: %f %f %f %f %f \n", pt, eta, phi, rat, ptg) ;
	  } 	  
	}
      }

      dynamic_cast<TH2F*>
	(fListHistos->FindObject("InvMassPairLeading"))->
	Fill(ptg,1.35);
      dynamic_cast<TH2F*>
	(fListHistos->FindObject("AnglePairLeading"))->
	Fill(ef,0.1);
    }
  }//No fOptfast
  else{
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
      if(ksPdg == 111){ //2 gamma
	rat = ptl/ptg ;
        phil = particle->Phi() ;
	if((ptl> pt)&& (rat > fRatioMinCut) && (rat < fRatioMaxCut) && 
	   ((phig-phil)>fPhiMinCut)&&((phig-phil)<fPhiMaxCut)){
	  eta = particle->Eta() ; 
	  phi = phil ;
	  pt  = ptl ;
	}// cuts
      }// pdg = 111
      if(ksPdg == 22){//1 gamma
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
  
	      dynamic_cast<TH2F*>
		(fListHistos->FindObject("AnglePairNoCut"))->
		Fill(e,angle);
	      dynamic_cast<TH2F*>
		(fListHistos->FindObject("InvMassPairNoCut"))->
		Fill(ptg,invmass);

	      if((rat > fRatioMinCut) && (rat < fRatioMaxCut) && 
		 ((phig-phil)>fPhiMinCut)&&((phig-phil)<fPhiMaxCut)){
	    
		dynamic_cast<TH2F*>
		  (fListHistos->FindObject("AnglePairLeadingCut"))->
		  Fill(e,angle);
		dynamic_cast<TH2F*>
		  (fListHistos->FindObject("InvMassPairLeadingCut"))->
		  Fill(ptg,invmass);
		
		if(IsAngleInWindow(angle,e)){
		  dynamic_cast<TH2F*>
		    (fListHistos->FindObject("AnglePairAngleCut"))->
		    Fill(e,angle);
		  dynamic_cast<TH2F*>
		    (fListHistos->FindObject("InvMassPairAngleCut"))->
		    Fill(ptg,invmass);
		  
		  //Info("GetLeadingPi0","InvMass %f", invmass);
		  if((invmass>fInvMassMinCut) && (invmass<fInvMassMaxCut)){ 
		    dynamic_cast<TH2F*>
		      (fListHistos->FindObject("InvMassPairAllCut"))->
		      Fill(ptg,invmass);
		    dynamic_cast<TH2F*>
		      (fListHistos->FindObject("AnglePairAllCut"))->
		      Fill(e,angle);
		    if(ptl > pt ){
		      pt       = ptl;
		      eta      = particle->Eta() ; 
		      phi      = phil ;
		      ef       = e ;
		      anglef   = angle ;
		      invmassf = invmass ;

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

    if(ef > 0.){
      dynamic_cast<TH2F*>
	(fListHistos->FindObject("InvMassPairLeading"))->
	Fill(ptg,invmassf);
      dynamic_cast<TH2F*>
	(fListHistos->FindObject("AnglePairLeading"))->
	Fill(ef,anglef);
    }
  }//fOptFast
  //  printf("GetLeadingPi0: %f %f %f \n", pt, eta, phi) ;
}


//____________________________________________________________________________
void AliPHOSGammaJet::InitParameters()
{

  //Initialize the parameters of the analysis.

  fAngleMaxParam.Set(4) ;
  fAngleMaxParam.AddAt(0.4,0);//={0.4,-0.25,0.025,-2e-4};
  fAngleMaxParam.AddAt(-0.25,1) ;
  fAngleMaxParam.AddAt(0.025,2) ;
  fAngleMaxParam.AddAt(-2e-4,3) ;
  fAnyConeOrPt    = kFALSE ;
  fOutputFileName = "GammaJet.root" ;
  fOptionGJ         = "";
  fHIJINGFileName = "galice.root" ;
  fHIJING         = kFALSE ;
  fMinDistance    = 3.6 ;
  fEtaCut         = 0.7 ;
  fInvMassMaxCut  = 0.15 ;
  fInvMassMinCut  = 0.12 ;
  fOnlyCharged    = kFALSE ;
  fOptFast        = kFALSE ;
  fESDdata        = kTRUE ;
  fPhiEMCALCut[0] = 60. *TMath::Pi()/180.;
  fPhiEMCALCut[1] = 180.*TMath::Pi()/180.;
  fPhiMaxCut      = 3.4 ;
  fPhiMinCut      = 2.9 ;
  fPtCut          = 10. ;
  fNeutralPtCut   = 0.5 ;
  fChargedPtCut   = 0.5 ;
  fTPCCutsLikeEMCAL  = kFALSE ;
  //Jet selection parameters
  //Fixed cut (old)
  fRatioMaxCut    = 1.0 ;
  fRatioMinCut    = 0.1 ; 
  fJetRatioMaxCut = 1.2 ; 
  fJetRatioMinCut = 0.8 ; 
  fJetTPCRatioMaxCut = 1.2 ;
  fJetTPCRatioMinCut = 0.3 ;
  fSelect         = kFALSE  ;

  fDirName      = "./" ;
  fESDTree      = "esdTree" ;
  fPattern      = "." ;

  //Cut depending on gamma energy

  fPtJetSelectionCut = 20.; //For Low pt jets+BKG, another limits applyed
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
  //Index > 2, same for TPC conf
  fBkgMean[0] = 0.; fBkgMean[1] = 8.8 ; fBkgMean[2] = 69.5;
  fBkgMean[3] = 0.; fBkgMean[4] = 6.4;  fBkgMean[5] = 48.6;
  fBkgRMS[0]  = 0.; fBkgRMS[1]  = 7.5;  fBkgRMS[2]  = 22.0; 
  fBkgRMS[3]  = 0.; fBkgRMS[4]  = 5.4;  fBkgRMS[5]  = 13.2; 

  //Factor x of min/max = E -+ x * sigma. Obtained after selecting the
  //limits for monoenergetic jets.
  //Index 0-> No BKG; Index 1-> BKG > 2 GeV; 
  //Index 2-> (low pt jets) BKG > 0.5 GeV;
  //Index > 2, same for TPC conf

  fJetXMin1[0] =-0.69 ; fJetXMin1[1] = 0.39 ; fJetXMin1[2] =-0.88 ; 
  fJetXMin1[3] =-2.0  ; fJetXMin1[4] =-0.442 ; fJetXMin1[5] =-1.1  ;
  fJetXMin2[0] = 0.066; fJetXMin2[1] = 0.038; fJetXMin2[2] = 0.034; 
  fJetXMin2[3] = 0.25 ; fJetXMin2[4] = 0.113; fJetXMin2[5] = 0.077 ;
  fJetXMax1[0] =-3.8  ; fJetXMax1[1] =-0.76 ; fJetXMax1[2] =-3.6  ; 
  fJetXMax1[3] =-2.7  ; fJetXMax1[4] =-1.21 ; fJetXMax1[5] =-3.7  ;
  fJetXMax2[0] =-0.012; fJetXMax2[1] =-0.022; fJetXMax2[2] = 0.016; 
  fJetXMax2[3] =-0.024; fJetXMax2[4] =-0.008; fJetXMax2[5] = 0.027;


  //Photon fast reconstruction
  fResPara1       = 0.0255 ;    // GeV
  fResPara2       = 0.0272 ; 
  fResPara3       = 0.0129 ; 
  
  fPosParaA      = 0.096 ;    // cm
  fPosParaB      = 0.229 ;  
  
  //Different cones and pt thresholds to construct the jet

  fCone        = 0.3  ;
  fPtThreshold = 0.   ;
  fNCone       = 1    ;
  fNPt         = 1    ;
  fCones[0]    = 0.3  ; fNameCones[0]   = "03" ;
  fPtThres[0]  = 0.5  ; fNamePtThres[0] = "05" ;

}

//__________________________________________________________________________-
Bool_t AliPHOSGammaJet::IsAngleInWindow(const Float_t angle,const Float_t e) {
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
Bool_t AliPHOSGammaJet::IsJetSelected(const Double_t ptg, const Double_t ptj, 
				      const TString type ){
  //Check if the energy of the reconstructed jet is within an energy window

  Double_t par[6];
  Double_t xmax[2];
  Double_t xmin[2];

  Int_t iTPC = 0;

  if(type == "TPC" && !fTPCCutsLikeEMCAL){
    iTPC = 3 ;//If(fTPCCutsLikeEMCAL) take jet energy cuts like EMCAL
  }


  if(!fHIJING){
    //Phythia alone, jets with pt_th > 0.2, r = 0.3 
    par[0] = fJetE1[0]; par[1] = fJetE2[0]; 
    //Energy of the jet peak
    //e_jet = fJetE1[0]+fJetE2[0]*e_gamma, simulation fit
    par[2] = fJetSigma1[0]; par[3] = fJetSigma2[0];
    //Sigma  of the jet peak
    //sigma_jet = fJetSigma1[0]+fJetSigma2[0]*e_gamma, simulation fit
    par[4] = fBkgMean[0 + iTPC]; par[5] = fBkgRMS[0 + iTPC];
    //Parameters reserved for HIJING bkg.
    xmax[0] = fJetXMax1[0 + iTPC]; xmax[1] = fJetXMax2[0 + iTPC];
    xmin[0] = fJetXMin1[0 + iTPC]; xmin[1] = fJetXMin2[0 + iTPC];
    //Factor that multiplies sigma to obtain the best limits, 
    //by observation, of mono jet ratios (ptjet/ptg)
    //X_jet = fJetX1[0]+fJetX2[0]*e_gamma
   
  }
  else{
    if(ptg > fPtJetSelectionCut){
      //Phythia +HIJING with  pt_th > 2 GeV/c, r = 0.3 
      par[0] = fJetE1[0]; par[1] = fJetE2[0]; 
      //Energy of the jet peak, same as in pp
      //e_jet = fJetE1[0]+fJetE2[0]*e_gamma, simulation fit
      par[2] = fJetSigma1[0]; par[3] = fJetSigma2[0];
      //Sigma  of the jet peak, same as in pp
      //sigma_jet = fJetSigma1[0]+fJetSigma2[0]*e_gamma, simulation fit
      par[4] = fBkgMean[1 + iTPC]; par[5] = fBkgRMS[1 + iTPC];
      //Mean value and RMS of HIJING Bkg 
      xmax[0] = fJetXMax1[1 + iTPC]; xmax[1] = fJetXMax2[1 + iTPC];
      xmin[0] = fJetXMin1[1 + iTPC]; xmin[1] = fJetXMin2[1 + iTPC];
      //Factor that multiplies sigma to obtain the best limits, 
      //by observation, of mono jet ratios (ptjet/ptg) mixed with HIJING Bkg, 
      //pt_th > 2 GeV, r = 0.3
      //X_jet = fJetX1[0]+fJetX2[0]*e_gamma
     
    }
    else{
      //Phythia + HIJING with  pt_th > 0.5 GeV/c, r = 0.3
      par[0] = fJetE1[1]; par[1] = fJetE2[1]; 
      //Energy of the jet peak, pt_th > 2 GeV/c, r = 0.3 
      //e_jet = fJetE1[0]+fJetE2[0]*e_gamma, simulation fit
      par[2] = fJetSigma1[1]; par[3] = fJetSigma2[1];
      //Sigma  of the jet peak, pt_th > 2 GeV/c, r = 0.3
      //sigma_jet = fJetSigma1[0]+fJetSigma2[0]*e_gamma, simulation fit
      par[4] = fBkgMean[2 + iTPC]; par[5] = fBkgRMS[2 + iTPC];
      //Mean value and RMS of HIJING Bkg in a 0.3 cone, pt > 2 GeV.
      xmax[0] = fJetXMax1[2 + iTPC]; xmax[1] = fJetXMax2[2 + iTPC];
      xmin[0] = fJetXMin1[2 + iTPC]; xmin[1] = fJetXMin2[2 + iTPC];
      //Factor that multiplies sigma to obtain the best limits, 
      //by observation, of mono jet ratios (ptjet/ptg) mixed with HIJING Bkg, 
      //pt_th > 2 GeV, r = 0.3
      //X_jet = fJetX1[0]+fJetX2[0]*e_gamma
     
    }//If low pt jet in bkg
  }//if Bkg

 //Calculate minimum and maximum limits of the jet ratio.
  Double_t min = CalculateJetRatioLimit(ptg, par, xmin);
  Double_t max = CalculateJetRatioLimit(ptg, par, xmax);
  
  //Info("IsJetSeleted","%s : Limits min %f, max %f, ptg / ptj %f",
  //   type.Data(),min,max,ptj/ptg);
  if(( min < ptj/ptg ) && ( max > ptj/ptg))
    return kTRUE;
  else
    return kFALSE;

}

//____________________________________________________________________________
void AliPHOSGammaJet::List() const
{
  // List the histos

  Info("List", "%d histograms found", fListHistos->GetEntries() ) ; 
  TIter next(fListHistos) ; 
  TH2F * h = 0 ; 
  while ( (h = dynamic_cast<TH2F*>(next())) )
    Info("List", "%s", h->GetName()) ; 
}

//____________________________________________________________________________
Double_t AliPHOSGammaJet::MakeEnergy(const Double_t energy)
{  
  // Smears the energy according to the energy dependent energy resolution.
  // A gaussian distribution is assumed
  
  Double_t sigma  = SigmaE(energy) ; 
  return  fRan.Gaus(energy, sigma) ;   


}
//____________________________________________________________________________
void AliPHOSGammaJet::MakeHistos()
{
  // Create histograms to be saved in output file and 
  // stores them in a TObjectArray
  
  fOutputFile = new TFile(fOutputFileName, "recreate") ;
  
  fListHistos = new TObjArray(10000) ;
  
   // Histos gamma pt vs leading pt
  
  TH1F * hPtSpectra  = new TH1F
    ("PtSpectra","p_{T i} vs p_{T #gamma}",200,0,200); 
  hPtSpectra->SetXTitle("p_{T} (GeV/c)");
  fListHistos->Add(hPtSpectra) ;
  
  //Histos ratio charged leading pt / gamma pt vs pt 
  
  TH2F * hChargeRatio  = new TH2F
    ("ChargeRatio","p_{T leading charge} /p_{T #gamma} vs p_{T #gamma}",
     120,0,120,120,0,1); 
  hChargeRatio->SetYTitle("p_{T lead charge} /p_{T #gamma}");
  hChargeRatio->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hChargeRatio) ;
  
  TH2F * hTPCRatio  = new TH2F
    ("TPCRatio","p_{T leading charge} /p_{T #gamma} vs p_{T #gamma}",
     120,0,120,120,0,1); 
  hTPCRatio->SetYTitle("p_{T lead charge} /p_{T #gamma}");
  hTPCRatio->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hTPCRatio) ;
  
  
  TH2F * hPi0Ratio  = new TH2F
    ("Pi0Ratio","p_{T leading  #pi^{0}} /p_{T #gamma} vs p_{T #gamma}",
     120,0,120,120,0,1); 
  hPi0Ratio->SetYTitle("p_{T lead  #pi^{0}} /p_{T #gamma}");
  hPi0Ratio->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hPi0Ratio) ;
  
  TH2F * hPhiGamma  = new TH2F
    ("PhiGamma","#phi_{#gamma}",200,0,120,200,0,7); 
  hPhiGamma->SetYTitle("#phi");
  hPhiGamma->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hPhiGamma) ; 
  
  TH2F * hEtaGamma  = new TH2F
    ("EtaGamma","#phi_{#gamma}",200,0,120,200,-0.8,0.8); 
  hEtaGamma->SetYTitle("#eta");
  hEtaGamma->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hEtaGamma) ;
 
//   //Jet reconstruction check
//   TH2F * hDeltaPhiJet  = new TH2F
//     ("DeltaPhiJet","#phi_{jet} - #phi_{pyth jet} vs p_{T #gamma}",
//      200,0,120,200,-1.5,1.5); 
//   hDeltaPhiJet->SetYTitle("#Delta #phi");
//   hDeltaPhiJet->SetXTitle("p_{T #gamma} (GeV/c)");
//   fListHistos->Add(hDeltaPhiJet) ;   

//   TH2F * hDeltaPhiTPCJet  = new TH2F
//     ("DeltaPhiTPCJet","#phi_{jet TPC} - #phi_{pyth jet} vs p_{T #gamma}",
//      200,0,120,200,-1.5,1.5); 
//   hDeltaPhiTPCJet->SetYTitle("#Delta #phi");
//   hDeltaPhiTPCJet->SetXTitle("p_{T #gamma} (GeV/c)");
//   fListHistos->Add(hDeltaPhiTPCJet) ; 

//   TH2F * hDeltaEtaJet  = new TH2F
//     ("DeltaEtaJet","#phi_{jet} - #phi_{pyth jet} vs p_{T #gamma}",
//      200,0,120,200,-1.5,1.5); 
//   hDeltaEtaJet->SetYTitle("#Delta #phi");
//   hDeltaEtaJet->SetXTitle("p_{T #gamma} (GeV/c)");
//   fListHistos->Add(hDeltaEtaJet) ;   

//   TH2F * hDeltaEtaTPCJet  = new TH2F
//     ("DeltaEtaTPCJet","#phi_{jet TPC} - #phi_{pyth jet} vs p_{T #gamma}",
//      200,0,120,200,-1.5,1.5); 
//   hDeltaEtaTPCJet->SetYTitle("#Delta #phi");
//   hDeltaEtaTPCJet->SetXTitle("p_{T #gamma} (GeV/c)");
//   fListHistos->Add(hDeltaEtaTPCJet) ; 

//   TH2F * hDeltaPtJet  = new TH2F
//     ("DeltaPtJet","#phi_{jet} - #phi_{pyth jet} vs p_{T #gamma}",
//      200,0,120,200,0.,100.); 
//   hDeltaPtJet->SetYTitle("#Delta #phi");
//   hDeltaPtJet->SetXTitle("p_{T #gamma} (GeV/c)");
//   fListHistos->Add(hDeltaPtJet) ;   

//   TH2F * hDeltaPtTPCJet  = new TH2F
//     ("DeltaPtTPCJet","#phi_{jet TPC} - #phi_{pyth jet} vs p_{T #gamma}",
//      200,0,120,200,0.,100.); 
//   hDeltaPtTPCJet->SetYTitle("#Delta #phi");
//   hDeltaPtTPCJet->SetXTitle("p_{T #gamma} (GeV/c)");
//   fListHistos->Add(hDeltaPtTPCJet) ; 

  //
  TH2F * hDeltaPhiCharge  = new TH2F
    ("DeltaPhiCharge","#phi_{#gamma} - #phi_{charge} vs p_{T #gamma}",
     200,0,120,200,0,6.4); 
  hDeltaPhiCharge->SetYTitle("#Delta #phi");
  hDeltaPhiCharge->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hDeltaPhiCharge) ; 
  
  TH2F * hDeltaPhiTPC  = new TH2F
    ("DeltaPhiTPC","#phi_{#gamma} - #phi_{charge} vs p_{T #gamma}",
     200,0,120,200,0,6.4); 
  hDeltaPhiTPC->SetYTitle("#Delta #phi");
  hDeltaPhiTPC->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hDeltaPhiTPC) ; 	 
  TH2F * hDeltaPhiPi0  = new TH2F
    ("DeltaPhiPi0","#phi_{#gamma} - #phi_{ #pi^{0}} vs p_{T #gamma}",
     200,0,120,200,0,6.4); 
  hDeltaPhiPi0->SetYTitle("#Delta #phi");
  hDeltaPhiPi0->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hDeltaPhiPi0) ; 

  TH2F * hDeltaEtaCharge  = new TH2F
    ("DeltaEtaCharge","#eta_{#gamma} - #eta_{charge} vs p_{T #gamma}",
     200,0,120,200,-2,2); 
  hDeltaEtaCharge->SetYTitle("#Delta #eta");
  hDeltaEtaCharge->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hDeltaEtaCharge) ; 

  TH2F * hDeltaEtaTPC  = new TH2F
    ("DeltaEtaTPC","#eta_{#gamma} - #eta_{charge} vs p_{T #gamma}",
     200,0,120,200,-2,2); 
  hDeltaEtaTPC->SetYTitle("#Delta #eta");
  hDeltaEtaTPC->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hDeltaEtaTPC) ; 
	 
  TH2F * hDeltaEtaPi0  = new TH2F
    ("DeltaEtaPi0","#eta_{#gamma} - #eta_{ #pi^{0}} vs p_{T #gamma}",
     200,0,120,200,-2,2); 
  hDeltaEtaPi0->SetYTitle("#Delta #eta");
  hDeltaEtaPi0->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hDeltaEtaPi0) ; 

  if(fOptFast){
      
    TH2F * hAnglePair  = new TH2F
      ("AnglePair",
       "Angle between #pi^{0} #gamma pair vs p_{T  #pi^{0}}",
       200,0,50,200,0,0.2); 
    hAnglePair->SetYTitle("Angle (rad)");
    hAnglePair->SetXTitle("E_{ #pi^{0}} (GeV/c)");
    fListHistos->Add(hAnglePair) ; 
    
    TH2F * hAnglePairAccepted  = new TH2F
      ("AnglePairAccepted",
       "Angle between #pi^{0} #gamma pair vs p_{T  #pi^{0}}, both #gamma in eta<0.7, inside window",
       200,0,50,200,0,0.2); 
    hAnglePairAccepted->SetYTitle("Angle (rad)");
    hAnglePairAccepted->SetXTitle("E_{ #pi^{0}} (GeV/c)");
    fListHistos->Add(hAnglePairAccepted) ; 
    
    TH2F * hAnglePairNoCut  = new TH2F
      ("AnglePairNoCut",
       "Angle between all #gamma pair vs p_{T  #pi^{0}}",200,0,50,200,0,0.2); 
    hAnglePairNoCut->SetYTitle("Angle (rad)");
    hAnglePairNoCut->SetXTitle("E_{ #pi^{0}} (GeV/c)");
    fListHistos->Add(hAnglePairNoCut) ; 
    
    TH2F * hAnglePairLeadingCut  = new TH2F
      ("AnglePairLeadingCut",
       "Angle between all #gamma pair that have a good phi and pt vs p_{T  #pi^{0}}",
       200,0,50,200,0,0.2); 
    hAnglePairLeadingCut->SetYTitle("Angle (rad)");
    hAnglePairLeadingCut->SetXTitle("E_{ #pi^{0}} (GeV/c)");
    fListHistos->Add(hAnglePairLeadingCut) ; 
    
    TH2F * hAnglePairAngleCut  = new TH2F
      ("AnglePairAngleCut",
       "Angle between all #gamma pair (angle + leading cut) vs p_{T  #pi^{0}}"
       ,200,0,50,200,0,0.2); 
    hAnglePairAngleCut->SetYTitle("Angle (rad)");
    hAnglePairAngleCut->SetXTitle("E_{ #pi^{0}} (GeV/c)");
    fListHistos->Add(hAnglePairAngleCut) ;
    
    TH2F * hAnglePairAllCut  = new TH2F
      ("AnglePairAllCut",
       "Angle between all #gamma pair (angle + inv mass cut+leading) vs p_{T  #pi^{0}}"
       ,200,0,50,200,0,0.2); 
    hAnglePairAllCut->SetYTitle("Angle (rad)");
    hAnglePairAllCut->SetXTitle("E_{ #pi^{0}} (GeV/c)");
    fListHistos->Add(hAnglePairAllCut) ; 
      
    TH2F * hAnglePairLeading  = new TH2F
      ("AnglePairLeading",
       "Angle between all #gamma pair finally selected vs p_{T  #pi^{0}}",
       200,0,50,200,0,0.2); 
    hAnglePairLeading->SetYTitle("Angle (rad)");
    hAnglePairLeading->SetXTitle("E_{ #pi^{0}} (GeV/c)");
    fListHistos->Add(hAnglePairLeading) ; 
    
    
    TH2F * hInvMassPairNoCut  = new TH2F
      ("InvMassPairNoCut","Invariant Mass of all #gamma pair vs p_{T #gamma}",
       120,0,120,360,0,0.5); 
    hInvMassPairNoCut->SetYTitle("Invariant Mass (GeV/c^{2})");
    hInvMassPairNoCut->SetXTitle("p_{T #gamma} (GeV/c)");
    fListHistos->Add(hInvMassPairNoCut) ; 
    
    TH2F * hInvMassPairLeadingCut  = new TH2F
      ("InvMassPairLeadingCut",
       "Invariant Mass of #gamma pair (leading cuts) vs p_{T #gamma}",
       120,0,120,360,0,0.5); 
    hInvMassPairLeadingCut->SetYTitle("Invariant Mass (GeV/c^{2})");
    hInvMassPairLeadingCut->SetXTitle("p_{T #gamma} (GeV/c)");
    fListHistos->Add(hInvMassPairLeadingCut) ; 
    
    TH2F * hInvMassPairAngleCut  = new TH2F
      ("InvMassPairAngleCut",
       "Invariant Mass of #gamma pair (angle cut) vs p_{T #gamma}",
       120,0,120,360,0,0.5); 
    hInvMassPairAngleCut->SetYTitle("Invariant Mass (GeV/c^{2})");
    hInvMassPairAngleCut->SetXTitle("p_{T #gamma} (GeV/c)");
    fListHistos->Add(hInvMassPairAngleCut) ; 
    
    
    TH2F * hInvMassPairAllCut  = new TH2F
      ("InvMassPairAllCut",
       "Invariant Mass of #gamma pair (angle+invmass cut+leading) vs p_{T #gamma}",
       120,0,120,360,0,0.5); 
    hInvMassPairAllCut->SetYTitle("Invariant Mass (GeV/c^{2})");
    hInvMassPairAllCut->SetXTitle("p_{T #gamma} (GeV/c)");
    fListHistos->Add(hInvMassPairAllCut) ; 
    
    TH2F * hInvMassPairLeading  = new TH2F
      ("InvMassPairLeading",
       "Invariant Mass of #gamma pair selected vs p_{T #gamma}",
       120,0,120,360,0,0.5); 
    hInvMassPairLeading->SetYTitle("Invariant Mass (GeV/c^{2})");
    hInvMassPairLeading->SetXTitle("p_{T #gamma} (GeV/c)");
    fListHistos->Add(hInvMassPairLeading) ; 
  }
  
  //Count
  
  TH1F * hNGamma  = new TH1F("NGamma","Number of #gamma over PHOS",240,0,120); 
  hNGamma->SetYTitle("N");
  hNGamma->SetXTitle("p_{T #gamma}(GeV/c)");
  fListHistos->Add(hNGamma) ; 
  
  TH1F * hNBkg = new TH1F("NBkg","bkg multiplicity",9000,0,9000); 
  hNBkg->SetYTitle("counts");
  hNBkg->SetXTitle("N");
  fListHistos->Add(hNBkg) ; 
  
  TH2F * hNLeading  = new TH2F
    ("NLeading","Accepted Jet Leading", 240,0,120,240,0,120); 
  hNLeading->SetYTitle("p_{T charge} (GeV/c)");
  hNLeading->SetXTitle("p_{T #gamma}(GeV/c)");
  fListHistos->Add(hNLeading) ; 
  
  
  TH1F * hN  = new TH1F("NJet","Accepted jets",240,0,120); 
  hN->SetYTitle("N");
  hN->SetXTitle("p_{T #gamma}(GeV/c)");
  fListHistos->Add(hN) ; 
  
  
  //Ratios and Pt dist of reconstructed (not selected) jets
  //Jet
  TH2F * hJetRatio  = new TH2F
    ("JetRatio","p_{T jet lead}/p_{T #gamma} vs p_{T #gamma}",
     240,0,120,200,0,10);
  hJetRatio->SetYTitle("p_{T jet lead #pi^{0}}/p_{T #gamma}");
  hJetRatio->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hJetRatio) ; 
  

  TH2F * hJetPt  = new TH2F
    ("JetPt", "p_{T jet lead} vs p_{T #gamma}",240,0,120,400,0,200);
  hJetPt->SetYTitle("p_{T jet lead #pi^{0}}/p_{T #gamma}");
  hJetPt->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hJetPt) ; 
  

  //Bkg
  
  TH2F * hBkgRatio  = new TH2F
    ("BkgRatio","p_{T bkg lead}/p_{T #gamma} vs p_{T #gamma}",
     240,0,120,200,0,10);
  hBkgRatio->SetYTitle("p_{T bkg lead charge}/p_{T #gamma}");
  hBkgRatio->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hBkgRatio) ;
  
  
  TH2F * hBkgPt  = new TH2F
    ("BkgPt","p_{T jet lead} vs p_{T #gamma}",240,0,120,400,0,200);
  hBkgPt->SetYTitle("p_{T jet lead charge}/p_{T #gamma}");
  hBkgPt->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hBkgPt) ;
 

  //Jet Distributions
  
  TH2F * hJetFragment  = 
    new TH2F("JetFragment","x = p_{T i charged}/p_{T #gamma}",
	     240,0.,120.,1000,0.,1.2); 
  hJetFragment->SetYTitle("x_{T}");
  hJetFragment->SetXTitle("p_{T #gamma}");
  fListHistos->Add(hJetFragment) ;
  
  TH2F * hBkgFragment  = new TH2F
    ("BkgFragment","x = p_{T i charged}/p_{T #gamma}",
     240,0.,120.,1000,0.,1.2);
  hBkgFragment->SetYTitle("x_{T}");
  hBkgFragment->SetXTitle("p_{T #gamma}");
  fListHistos->Add(hBkgFragment) ;

  TH2F * hJetPtDist  = 
    new TH2F("JetPtDist","x = p_{T i charged}",240,0.,120.,400,0.,200.); 
  hJetPtDist->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hJetPtDist) ;
  
  TH2F * hBkgPtDist  = new TH2F
    ("BkgPtDist","x = p_{T i charged}",240,0.,120.,400,0.,200.); 
  hBkgPtDist->SetXTitle("p_{T #gamma} (GeV/c)");
  fListHistos->Add(hBkgPtDist) ;
  

  if(fOnlyCharged){
    //Counts 
    TH1F * hNBkgTPC  = new TH1F
      ("NBkgTPC","TPC bkg multiplicity ",9000,0,9000); 
    hNBkgTPC->SetYTitle("counts");
    hNBkgTPC->SetXTitle("N");
    fListHistos->Add(hNBkgTPC) ;
    
    TH2F * hNTPCLeading  = new TH2F
      ("NTPCLeading","Accepted TPC jet leading",240,0,120,240,0,120); 
    hNTPCLeading->SetYTitle("p_{T charge} (GeV/c)");
    hNTPCLeading->SetXTitle("p_{T #gamma}(GeV/c)");
    fListHistos->Add(hNTPCLeading) ; 

    TH1F * hNTPC  = new TH1F("NTPCJet","Number of TPC jets",240,0,120); 
    hNTPC->SetYTitle("N");
    hNTPC->SetXTitle("p_{T #gamma}(GeV/c)");
    fListHistos->Add(hNTPC) ;
 
    TH2F * hJetTPCRatio  = new TH2F
      ("JetTPCRatio", "p_{T jet lead TPC}/p_{T #gamma} vs p_{T #gamma}",
       240,0,120,200,0,10);
    hJetTPCRatio->SetYTitle("p_{T jet lead TPC}/p_{T #gamma}");
    hJetTPCRatio->SetXTitle("p_{T #gamma} (GeV/c)");
    fListHistos->Add(hJetTPCRatio) ; 
    
    TH2F * hBkgTPCRatio  = new TH2F
      ("BkgTPCRatio","p_{T bkg lead TPC}/p_{T #gamma} vs p_{T #gamma}",
       240,0,120,200,0,10);
    hBkgTPCRatio->SetYTitle("p_{T bkg lead TPC}/p_{T #gamma}");
    hBkgTPCRatio->SetXTitle("p_{T #gamma} (GeV/c)");
    fListHistos->Add(hBkgTPCRatio) ; 
    
    TH2F * hJetTPCPt  = new TH2F
      ("JetTPCPt", "p_{T jet lead TPC} vs p_{T #gamma}",240,0,120,400,0,200);
    hJetTPCPt->SetYTitle("p_{T jet lead TPC}/p_{T #gamma}");
    hJetTPCPt->SetXTitle("p_{T #gamma} (GeV/c)");
    fListHistos->Add(hJetTPCPt) ; 
    
    TH2F * hBkgTPCPt  = new TH2F
      ("BkgTPCPt", "p_{T bkg lead TPC} vs p_{T #gamma}",240,0,120,400,0,200);
    hBkgTPCPt->SetYTitle("p_{T bkg lead TPC}/p_{T #gamma}");
    hBkgTPCPt->SetXTitle("p_{T #gamma} (GeV/c)");
    fListHistos->Add(hBkgTPCPt) ; 

    //JetDistributions
    
    TH2F * hJetTPCFragment  = 
      new TH2F("JetTPCFragment","x = p_{T i charged}/p_{T #gamma}",
	       240,0.,120.,1000,0.,1.2); 
    hJetTPCFragment->SetYTitle("x_{T}");
    hJetTPCFragment->SetXTitle("p_{T #gamma}");
    fListHistos->Add(hJetTPCFragment) ;
    
    TH2F * hBkgTPCFragment  = new TH2F
      ("BkgTPCFragment","x = p_{T i charged}/p_{T #gamma}",
       240,0.,120.,1000,0.,1.2); 
    hBkgTPCFragment->SetYTitle("x_{T}");
    hBkgTPCFragment->SetXTitle("p_{T #gamma}");
    fListHistos->Add(hBkgTPCFragment) ;


    TH2F * hJetTPCPtDist  = new TH2F("JetTPCPtDist",
       "x = p_{T i charged}",240,0.,120.,400,0.,200.); 
    hJetTPCPtDist->SetXTitle("p_{T #gamma} (GeV/c)");
    fListHistos->Add(hJetTPCPtDist) ;
    
    TH2F * hBkgTPCPtDist  = new TH2F
      ("BkgTPCPtDist","x = p_{T i charged}",240,0.,120.,400,0.,200.); 
    hBkgTPCPtDist->SetXTitle("p_{T #gamma} (GeV/c)");
    fListHistos->Add(hBkgTPCPtDist) ;
    
  }
  

  if(fAnyConeOrPt){
    //If we want to study the jet for different cones and pt. Old version
 
    TH2F * hJetRatios[5][5];
    TH2F * hJetTPCRatios[5][5];
    
    TH2F * hJetPts[5][5];
    TH2F * hJetTPCPts[5][5];

    TH2F * hBkgRatios[5][5];
    TH2F * hBkgTPCRatios[5][5];
    
    TH2F * hBkgPts[5][5];
    TH2F * hBkgTPCPts[5][5];
        
    TH2F * hNLeadings[5][5];
    TH2F * hNTPCLeadings[5][5];
    
    TH1F * hNs[5][5];
    TH1F * hNTPCs[5][5];
    
    TH1F * hNBkgs[5][5];
    TH1F * hNBkgTPCs[5][5];

    TH2F * hJetFragments[5][5];
    TH2F * hBkgFragments[5][5];
    TH2F * hJetPtDists[5][5];
    TH2F * hBkgPtDists[5][5];
   
    TH2F * hJetTPCFragments[5][5];
    TH2F * hBkgTPCFragments[5][5];
    TH2F * hJetTPCPtDists[5][5];
    TH2F * hBkgTPCPtDists[5][5];

    
    for(Int_t icone = 0; icone<fNCone; icone++){
      for(Int_t ipt = 0; ipt<fNPt;ipt++){
	//Jet Pt / Gamma Pt 
	
	//Jet

	hJetRatios[icone][ipt]  = new TH2F
	  ("JetRatioCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt], 
	 "p_{T jet lead #pi^{0}}/p_{T #gamma} vs p_{T #gamma}, cone ="
	   +fNameCones[icone]+", pt>" +fNamePtThres[ipt]+" GeV/c",
	   240,0,120,200,0,10);
	hJetRatios[icone][ipt]->
	  SetYTitle("p_{T jet lead #pi^{0}}/p_{T #gamma}");
	hJetRatios[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	fListHistos->Add(hJetRatios[icone][ipt]) ; 
	

	hJetPts[icone][ipt]  = new TH2F
	  ("JetPtCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt], 
	   "p_{T jet lead #pi^{0}}/p_{T #gamma} vs p_{T #gamma}, cone ="
	   +fNameCones[icone]+", pt>" +fNamePtThres[ipt]+" GeV/c",
	   240,0,120,400,0,200);
	hJetPts[icone][ipt]->
	  SetYTitle("p_{T jet lead #pi^{0}}/p_{T #gamma}");
	hJetPts[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	fListHistos->Add(hJetPts[icone][ipt]) ; 
	
	
	//Bkg
	hBkgRatios[icone][ipt]  = new TH2F
	  ("BkgRatioCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt], 
	   "p_{T bkg lead #pi^{0}}/p_{T #gamma} vs p_{T #gamma}, cone ="
	   +fNameCones[icone]+", pt>" +fNamePtThres[ipt]+" GeV/c",
	   240,0,120,200,0,10);
	hBkgRatios[icone][ipt]->
	  SetYTitle("p_{T bkg lead #pi^{0}}/p_{T #gamma}");
	hBkgRatios[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	fListHistos->Add(hBkgRatios[icone][ipt]) ; 
	
	
	
	hBkgPts[icone][ipt]  = new TH2F
	  ("BkgPtCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt], 
	   "p_{T jet lead #pi^{0}}/p_{T #gamma} vs p_{T #gamma}, cone ="
	   +fNameCones[icone]+", pt>" +fNamePtThres[ipt]+" GeV/c",
	   240,0,120,400,0,200);
	hBkgPts[icone][ipt]->
	  SetYTitle("p_{T jet lead #pi^{0}}/p_{T #gamma}");
	hBkgPts[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	fListHistos->Add(hBkgPts[icone][ipt]) ; 
	
	
	if(fOnlyCharged){
	  
	  hJetTPCRatios[icone][ipt]  = new TH2F
	    ("JetTPCRatioCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt], 
	     "p_{T jet lead TPC}/p_{T #gamma} vs p_{T #gamma}, cone ="
	     +fNameCones[icone]+", pt>" +fNamePtThres[ipt]+" GeV/c",
	     240,0,120,200,0,10);
	  hJetTPCRatios[icone][ipt]->SetYTitle("p_{T jet lead TPC}/p_{T #gamma}");
	  hJetTPCRatios[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	  fListHistos->Add(hJetTPCRatios[icone][ipt]) ; 
	  
	  hBkgTPCRatios[icone][ipt]  = new TH2F
	    ("BkgTPCRatioCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt], 
	     "p_{T bkg lead TPC}/p_{T #gamma} vs p_{T #gamma}, cone ="
	     +fNameCones[icone]+", pt>" +fNamePtThres[ipt]+" GeV/c",
	     240,0,120,200,0,10);
	  hBkgTPCRatios[icone][ipt]->SetYTitle("p_{T bkg lead TPC}/p_{T #gamma}");
	  hBkgTPCRatios[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	  fListHistos->Add(hBkgTPCRatios[icone][ipt]) ; 
	  
	  hJetTPCPts[icone][ipt]  = new TH2F
	    ("JetTPCPtCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt], 
	     "p_{T jet lead TPC}/p_{T #gamma} vs p_{T #gamma}, cone ="
	     +fNameCones[icone]+", pt>" +fNamePtThres[ipt]+" GeV/c",
	     240,0,120,400,0,200);
	  hJetTPCPts[icone][ipt]->SetYTitle("p_{T jet lead TPC}/p_{T #gamma}");
	  hJetTPCPts[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	  fListHistos->Add(hJetTPCPts[icone][ipt]) ; 
	  
	  hBkgTPCPts[icone][ipt]  = new TH2F
	    ("BkgTPCPtCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt], 
	     "p_{T bkg lead TPC}/p_{T #gamma} vs p_{T #gamma}, cone ="
	     +fNameCones[icone]+", pt>" +fNamePtThres[ipt]+" GeV/c",
	     240,0,120,400,0,200);
	  hBkgTPCPts[icone][ipt]->SetYTitle("p_{T bkg lead TPC}/p_{T #gamma}");
	  hBkgTPCPts[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	  fListHistos->Add(hBkgTPCPts[icone][ipt]) ; 
	  
	}
		
 	//Counts
 	hNBkgs[icone][ipt]  = new TH1F
	  ("NBkgCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	   "bkg multiplicity cone ="+fNameCones[icone]+", pt>" 
	   +fNamePtThres[ipt]+" GeV/c",9000,0,9000); 
	hNBkgs[icone][ipt]->SetYTitle("counts");
	hNBkgs[icone][ipt]->SetXTitle("N");
	fListHistos->Add(hNBkgs[icone][ipt]) ; 
	
	hNBkgTPCs[icone][ipt]  = new TH1F
	  ("NBkgTPCCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	   "bkg multiplicity cone ="+fNameCones[icone]+", pt>" 
	   +fNamePtThres[ipt]+" GeV/c",9000,0,9000); 
	hNBkgTPCs[icone][ipt]->SetYTitle("counts");
	hNBkgTPCs[icone][ipt]->SetXTitle("N");
	fListHistos->Add(hNBkgTPCs[icone][ipt]) ;
	
	
	hNLeadings[icone][ipt]  = new TH2F
	  ("NLeadingCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	   "p_{T #gamma} vs p_{T #pi^{0}} cone ="+fNameCones[icone]+", pt>" 
	   +fNamePtThres[ipt]+" GeV/c",120,0,120,120,0,120); 
	hNLeadings[icone][ipt]->SetYTitle("p_{T #pi^{0}}(GeV/c)");
	hNLeadings[icone][ipt]->SetXTitle("p_{T #gamma}(GeV/c)");
	fListHistos->Add(hNLeadings[icone][ipt]) ; 
	
	hNs[icone][ipt]  = new TH1F
	  ("NJetCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	   "Number of neutral jets, cone ="+fNameCones[icone]+", pt>" 
	   +fNamePtThres[ipt]+" GeV/c",120,0,120); 
	hNs[icone][ipt]->SetYTitle("N");
	hNs[icone][ipt]->SetXTitle("p_{T #gamma}(GeV/c)");
	fListHistos->Add(hNs[icone][ipt]) ; 
	
	//Fragmentation Function
	hJetFragments[icone][ipt]  = new TH2F
	  ("JetFragmentCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	   "x_{T} = p_{T i}/p_{T #gamma}, cone ="+fNameCones[icone]+", pt>" 
	   +fNamePtThres[ipt]+" GeV/c",120,0.,120.,240,0.,1.2); 
	hJetFragments[icone][ipt]->SetYTitle("x_{T}");
	hJetFragments[icone][ipt]->SetXTitle("p_{T #gamma}");
	fListHistos->Add(hJetFragments[icone][ipt]) ; 
	
	hBkgFragments[icone][ipt]  = new TH2F
	  ("BkgFragmentCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	   "x_{T} = p_{T i}/p_{T #gamma}, cone ="+fNameCones[icone]+", pt>" 
	   +fNamePtThres[ipt]+" GeV/c",120,0.,120.,240,0.,1.2); 
	hBkgFragments[icone][ipt]->SetYTitle("x_{T}");
	hBkgFragments[icone][ipt]->SetXTitle("p_{T #gamma}");
	fListHistos->Add(hBkgFragments[icone][ipt]) ; 
	
	
	//Jet particle distribution
	
	hJetPtDists[icone][ipt]  = new TH2F
	  ("JetPtDistCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	   "p_{T i}, cone ="+fNameCones[icone]+", pt>" +fNamePtThres[ipt]+
	   " GeV/c",120,0.,120.,120,0.,120.); 
	hJetPtDists[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	fListHistos->Add(hJetPtDists[icone][ipt]) ; 
	
	hBkgPtDists[icone][ipt]  = new TH2F
	  ("BkgPtDistCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	   "p_{T i}, cone ="+fNameCones[icone]+", pt>" +fNamePtThres[ipt]+
	   " GeV/c",120,0.,120.,120,0.,120.); 
	hBkgPtDists[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	fListHistos->Add(hBkgPtDists[icone][ipt]) ; 
	

	if(fOnlyCharged){ 
	  hNTPCLeadings[icone][ipt]  = new TH2F
	    ("NTPCLeadingCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	     "p_{T #gamma} vs p_{T charge} cone ="+fNameCones[icone]+", pt>" 
	     +fNamePtThres[ipt]+" GeV/c",120,0,120,120,0,120); 
	  hNTPCLeadings[icone][ipt]->SetYTitle("p_{T charge} (GeV/c)");
	  hNTPCLeadings[icone][ipt]->SetXTitle("p_{T #gamma}(GeV/c)");
	  fListHistos->Add(hNTPCLeadings[icone][ipt]) ; 
	  
	  hNTPCs[icone][ipt]  = new TH1F
	    ("NTPCJetCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	     "Number of charged jets, cone ="+fNameCones[icone]+", pt>" 
	     +fNamePtThres[ipt]+" GeV/c",120,0,120); 
	  hNTPCs[icone][ipt]->SetYTitle("N");
	  hNTPCs[icone][ipt]->SetXTitle("p_{T #gamma}(GeV/c)");
	  fListHistos->Add(hNTPCs[icone][ipt]) ; 
	  
	  hJetTPCFragments[icone][ipt]  = new TH2F
	    ("JetTPCFragmentCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	     "x = p_{T i charged}/p_{T #gamma}, cone ="+fNameCones[icone]+", pt>" 
	     +fNamePtThres[ipt]+" GeV/c",120,0.,120.,240,0.,1.2); 
	  hJetTPCFragments[icone][ipt]->SetYTitle("x_{T}");
	  hJetTPCFragments[icone][ipt]->SetXTitle("p_{T #gamma}");
	  fListHistos->Add(hJetTPCFragments[icone][ipt]) ;
	  
	  hBkgTPCFragments[icone][ipt]  = new TH2F
	    ("BkgTPCFragmentCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	     "x = p_{T i charged}/p_{T #gamma}, cone ="+fNameCones[icone]+", pt>" 
	     +fNamePtThres[ipt]+" GeV/c",120,0.,120.,240,0.,1.2); 
	  hBkgTPCFragments[icone][ipt]->SetYTitle("x_{T}");
	  hBkgTPCFragments[icone][ipt]->SetXTitle("p_{T #gamma}");
	  fListHistos->Add(hBkgTPCFragments[icone][ipt]) ;
	  
	  hJetTPCPtDists[icone][ipt]  = new TH2F
	    ("JetTPCPtDistCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	     "x = p_{T i charged}, cone ="+fNameCones[icone]+", pt>" 
	     +fNamePtThres[ipt]+" GeV/c",120,0.,120.,120,0.,120.); 
	  hJetTPCPtDists[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	  fListHistos->Add(hJetTPCPtDists[icone][ipt]) ;
	  
	  hBkgTPCPtDists[icone][ipt]  = new TH2F
	    ("BkgTPCPtDistCone"+fNameCones[icone]+"Pt"+fNamePtThres[ipt],
	     "x = p_{T i charged}, cone ="+fNameCones[icone]+", pt>" +
	     fNamePtThres[ipt]+" GeV/c",120,0.,120.,120,0.,120.); 
	  hBkgTPCPtDists[icone][ipt]->SetXTitle("p_{T #gamma} (GeV/c)");
	  fListHistos->Add(hBkgTPCPtDists[icone][ipt]) ;
	}
      }//ipt
    } //icone
  }//If we want to study any cone or pt threshold
}


//____________________________________________________________________________
void AliPHOSGammaJet::MakeJet(TClonesArray * pl, 
			      Double_t ptg, Double_t phig, 
			      Double_t ptl, Double_t  phil, Double_t  etal,
			      TString conf, TLorentzVector & jet)
{
  //Fill the jet with the particles around the leading particle with 
  //R=fCone and pt_th = fPtThres. Calculate the energy of the jet and 
  //check if we select it. Fill jet histograms
  Float_t ptcut = 0. ;
  if(fHIJING){
    if(ptg > fPtJetSelectionCut)  ptcut = 2. ;
    else                          ptcut = 0.5;
  }
  
  TClonesArray * jetList = new TClonesArray("TParticle",1000);
  TClonesArray * bkgList = new TClonesArray("TParticle",1000);
  
  TLorentzVector bkg(0,0,0,0);
  TLorentzVector lv (0,0,0,0);

  Double_t ptjet = 0.0;
  Double_t ptbkg = 0.0;
  Int_t n0 = 0;
  Int_t n1 = 0;  
  Bool_t b1 = kFALSE;
  Bool_t b0 = kFALSE;

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
      if(particle->Pt() > ptcut ){
	bkg+=lv;
	ptbkg+=particle->Pt();    
      }  
    }
  }
  
  ptjet = jet.Pt();
  ptbkg = bkg.Pt();

  if(strstr(fOptionGJ,"deb") || strstr(fOptionGJ,"deb all"))
    Info("MakeJet","Gamma   pt %f, Jet pt %f, Bkg pt %f",ptg,ptjet,ptbkg);
 

  //Fill histograms

  Double_t rat     = ptjet/ptg ;
  Double_t ratbkg  = ptbkg/ptg ;

  dynamic_cast<TH2F*> 	
    (fListHistos->FindObject("Jet"+conf+"Ratio"))->Fill(ptg,rat);
  dynamic_cast<TH2F*>
    (fListHistos->FindObject("Jet"+conf+"Pt"))   ->Fill(ptg,ptjet);
  dynamic_cast<TH2F*>
    (fListHistos->FindObject("Bkg"+conf+"Ratio"))->Fill(ptg,ratbkg);
  dynamic_cast<TH2F*>
    (fListHistos->FindObject("Bkg"+conf+"Pt"))   ->Fill(ptg,ptbkg);


  if(IsJetSelected(ptg,ptjet,conf) || fSelect){
    if(strstr(fOptionGJ,"deb") || strstr(fOptionGJ,"deb all"))
      Info("MakeJet","JetSelected");
    dynamic_cast<TH1F*>(fListHistos->FindObject("N"+conf+"Jet"))->
      Fill(ptg);
    dynamic_cast<TH2F*>(fListHistos->FindObject("N"+conf+"Leading"))
      ->Fill(ptg,ptl);
    FillJetHistos(jetList, ptg, conf, "Jet");
    FillJetHistos(bkgList, ptg, conf, "Bkg");
  }
  
  jetList ->Delete();
  bkgList ->Delete();
}

//____________________________________________________________________________
void AliPHOSGammaJet::MakeJetAnyConeOrPt(TClonesArray * pl, Double_t ptg, 
					 Double_t  phig, Double_t ptl, 
					 Double_t  phil, Double_t  etal, 
					 TString conf){

  //Fill the jet with the particles around the leading particle with 
  //R=fCone(i) and pt_th = fPtThres(i). Calculate the energy of the jet and 
  //check if we select it. Fill jet i histograms
  
  TClonesArray * jetList = new TClonesArray("TParticle",1000);
  TClonesArray * bkgList = new TClonesArray("TParticle",1000);
  
  Double_t ptjet  = 0.0;
  Double_t ptbkg  = 0.0;

  Int_t n1  = 0;
  Int_t n0  = 0;  
  Bool_t b1 = kFALSE;
  Bool_t b0 = kFALSE;
  
  //Create as many jets as cones and pt thresholds are defined
  Double_t maxcut = fJetRatioMaxCut;
  Double_t mincut = fJetRatioMinCut;
  
  if(conf == "TPC"){
    maxcut = fJetTPCRatioMaxCut;
    mincut = fJetTPCRatioMinCut;
  }

  Double_t ratjet  = 0;
  Double_t ratbkg  = 0;
  
  for(Int_t icone = 0; icone<fNCone; icone++) {
    for(Int_t ipt = 0; ipt<fNPt;ipt++) {
      
      TString cone  = fNameCones[icone]  ;
      TString ptcut = fNamePtThres[ipt] ;
      
      TIter next(pl) ; 
      TParticle * particle = 0 ;

      ptjet  = 0 ;
      ptbkg = 0 ;
 
      while ( (particle = dynamic_cast<TParticle*>(next())) ) {
	b1 = kFALSE;
	b0 = kFALSE;

	SetJet(particle, b0, fCones[icone],etal, phil) ;  
	SetJet(particle, b1, fCones[icone],etal, phig) ;  

	if(b0){
	  new((*jetList)[n0++]) TParticle(*particle) ;
	  if(particle->Pt() > fPtThres[ipt] )
	    ptjet+=particle->Pt();
	}
	if(b1) { 
	  new((*bkgList)[n1++]) TParticle(*particle) ;
	  if(particle->Pt() > fPtThres[ipt] )
	    ptbkg+=particle->Pt();
	}

      }

      //Fill histograms
      if(ptjet > 0.) {

	if(strstr(fOptionGJ,"deb")){
	  Info("MakeJetAnyPt","cone %f, ptcut %f",fCones[icone],fPtThres[ipt]);
	  Info("MakeJetAnyPt","pT: Gamma %f, Jet %f, Bkg  %f",ptg,ptjet,ptbkg);
	}

	ratjet  = ptjet /ptg;
	ratbkg  = ptbkg/ptg;
  
	dynamic_cast<TH2F*>
	  (fListHistos->FindObject("Jet"+conf+"RatioCone"+cone+"Pt"+ptcut))
	  ->Fill(ptg,ratjet);	 
	dynamic_cast<TH2F*>
	  (fListHistos->FindObject("Jet"+conf+"PtCone"+cone+"Pt"+ptcut))
	  ->Fill(ptg,ptjet);

	dynamic_cast<TH2F*>
	  (fListHistos->FindObject("Bkg"+conf+"RatioCone"+cone+"Pt"+ptcut))
	  ->Fill(ptg,ratbkg);

	dynamic_cast<TH2F*>
	  (fListHistos->FindObject("Bkg"+conf+"PtCone"+cone+"Pt"+ptcut))
	  ->Fill(ptg,ptbkg);

	
	//Select Jet 
	if((ratjet < maxcut) && (ratjet > mincut)){
	  
	  dynamic_cast<TH1F*>
	    (fListHistos->FindObject("N"+conf+"JetCone"+cone+"Pt"+ptcut))->
	    Fill(ptg);
	  dynamic_cast<TH2F*>
	    (fListHistos->FindObject("N"+conf+"LeadingCone"+cone+"Pt"+ptcut))
	    ->Fill(ptg,ptl);
	  
 	    FillJetHistosAnyConeOrPt
 	      (jetList,ptg,conf,"Jet",fNameCones[icone],fNamePtThres[ipt]);
 	    FillJetHistosAnyConeOrPt
 	      (bkgList,ptg,conf,"Bkg",fNameCones[icone],fNamePtThres[ipt]);

	}
      } //ptjet > 0 
      jetList ->Delete();
      bkgList ->Delete();
    }//for pt threshold
  }// for cone
}

//____________________________________________________________________________
void  AliPHOSGammaJet::MakePhoton(TLorentzVector & particle)
{
  //Fast reconstruction for photons
  Double_t energy = particle.E()  ;
  Double_t modenergy = MakeEnergy(energy) ;
  //Info("MakePhoton","Energy %f, Modif %f",energy,modenergy);

  // get the detected direction
    TVector3 pos = particle.Vect(); 
    pos*=460./energy;
    TVector3 modpos = MakePosition(energy, pos) ;
    modpos *= modenergy / 460.;
    
    Float_t modtheta = modpos.Theta();
    Float_t modphi   = modpos.Phi(); 
    
    // Set the modified 4-momentum of the reconstructed particle
    Float_t py = modenergy*TMath::Sin(modphi)*TMath::Sin(modtheta);
    Float_t px = modenergy*TMath::Cos(modphi)*TMath::Sin(modtheta);
    Float_t pz = modenergy*TMath::Cos(modtheta); 
    
    particle.SetPxPyPzE(px,py,pz,modenergy);

}

//____________________________________________________________________________
TVector3 AliPHOSGammaJet::MakePosition(const Double_t energy, const TVector3 pos)
{
  // Smears the impact position according to the energy dependent position resolution
  // A gaussian position distribution is assumed

  TVector3 newpos ;

  Double_t sigma = SigmaP(energy) ;
  Double_t x = fRan.Gaus( pos.X(), sigma ) ;
  Double_t z = fRan.Gaus( pos.Z(), sigma ) ;
  Double_t y = pos.Y() ; 
  
  newpos.SetX(x) ; 
  newpos.SetY(y) ; 
  newpos.SetZ(z) ; 

  // Info("MakePosition","Theta dif %f",pos.Theta()-newpos.Theta());
//   Info("MakePosition","Phi   dif %f",pos.Phi()-newpos.Phi());	      
  return newpos ; 
}

//____________________________________________________________________________
void AliPHOSGammaJet::Pi0Decay(Double_t mPi0, TLorentzVector &p0, 
			       TLorentzVector &p1, TLorentzVector &p2, Double_t &angle) {
  // Perform isotropic decay pi0 -> 2 photons
  // p0 is pi0 4-momentum (inut)
  // p1 and p2 are photon 4-momenta (output)
  //  cout<<"Boost vector"<<endl;
  TVector3 b = p0.BoostVector();
  //cout<<"Parameters"<<endl;
  //Double_t mPi0   = p0.M();
  Double_t phi    = TMath::TwoPi() * gRandom->Rndm();
  Double_t cosThe = 2 * gRandom->Rndm() - 1;
  Double_t cosPhi = TMath::Cos(phi);
  Double_t sinPhi = TMath::Sin(phi);
  Double_t sinThe = TMath::Sqrt(1-cosThe*cosThe);
  Double_t ePi0   = mPi0/2.;
  //cout<<"ePi0 "<<ePi0<<endl;
  //cout<<"Components"<<endl;
  p1.SetPx(+ePi0*cosPhi*sinThe);
  p1.SetPy(+ePi0*sinPhi*sinThe);
  p1.SetPz(+ePi0*cosThe);
  p1.SetE(ePi0);
  //cout<<"p1: "<<p1.Px()<<" "<<p1.Py()<<" "<<p1.Pz()<<" "<<p1.E()<<endl;
  //cout<<"p1 Mass: "<<p1.Px()*p1.Px()+p1.Py()*p1.Py()+p1.Pz()*p1.Pz()-p1.E()*p1.E()<<endl;
  p2.SetPx(-ePi0*cosPhi*sinThe);
  p2.SetPy(-ePi0*sinPhi*sinThe);
  p2.SetPz(-ePi0*cosThe);
  p2.SetE(ePi0);
  //cout<<"p2: "<<p2.Px()<<" "<<p2.Py()<<" "<<p2.Pz()<<" "<<p2.E()<<endl;
  //cout<<"p2 Mass: "<<p2.Px()*p2.Px()+p2.Py()*p2.Py()+p2.Pz()*p2.Pz()-p2.E()*p2.E()<<endl;
  //cout<<"Boost "<<b.X()<<" "<<b.Y()<<" "<<b.Z()<<endl;
  p1.Boost(b);
  //cout<<"p1: "<<p1.Px()<<" "<<p1.Py()<<" "<<p1.Pz()<<" "<<p1.E()<<endl;
  p2.Boost(b);
  //cout<<"p2: "<<p2.Px()<<" "<<p2.Py()<<" "<<p2.Pz()<<" "<<p2.E()<<endl;
  //cout<<"angle"<<endl;
  angle = p1.Angle(p2.Vect());
  //cout<<angle<<endl;
}

//____________________________________________________________________________
void AliPHOSGammaJet::Plot(TString what, Option_t * option) const
{
  //Plot some relevant histograms of the analysis
  TH2F * h = dynamic_cast<TH2F*>(fOutputFile->Get(what));
  if(h){
    h->Draw();
  }
  else if (what == "all") {
  TCanvas * can = new TCanvas("GammaJet", "Gamma-Jet Study",10,40,1600,1200);
  can->cd() ; 
  can->Range(0,0,22,20);
  TPaveLabel *pl1 = new TPaveLabel(1,18,20,19.5,"Titre","br");
  pl1->SetFillColor(18);
  pl1->SetTextFont(32);
  pl1->SetTextColor(49);
  pl1->Draw();
  Int_t index ; 
  TPad * pad = 0; 
  Float_t begx = -0.29, begy = 0., endx = 0., endy = 0.30 ; 
  for (index = 0 ; index < fListHistos->GetEntries() ; index++) {
    TString name("pad") ; 
    name += index ; 
    begx += 0.30 ;
    endx += 0.30 ;
    if (begx >= 1.0 || endx >= 1.0) {
      begx = 0.01 ; 
      endx = 0.30 ; 
      begy += 0.30 ;
      endy += 0.30 ; 
    } 
    printf("%f %f %f %f \n", begx, begy, endx, endy) ; 
    pad = new TPad(name,"This is a pad",begx,begy,endx,endy,33);
    pad->Draw();
    pad->cd() ; 
    fListHistos->At(index)->Draw(option) ; 
    pad->Modified() ; 
    pad->Update() ; 
    can->cd() ; 
  }

  }
  else{
    Info("Draw", "Histogram %s does not exist or unknown option", what.Data()) ;
     fOutputFile->ls();
  }
}

//____________________________________________________________________________
void AliPHOSGammaJet::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;

  Info("Print", "%s %s", GetName(), GetTitle() ) ;
 
  printf("Eta cut           : %f\n", fEtaCut) ;
  printf("D phi max cut     : %f\n", fPhiMaxCut) ; 
  printf("D phi min cut     : %f\n", fPhiMinCut) ;
  printf("Leading Ratio max cut     : %f\n", fRatioMaxCut) ; 
  printf("Leading Ratio min cut     : %f\n", fRatioMinCut) ;
  printf("Jet Ratio max cut     : %f\n", fJetRatioMaxCut) ; 
  printf("Jet Ratio min cut     : %f\n", fJetRatioMinCut) ;
  printf("Jet TPC Ratio max cut     : %f\n", fJetTPCRatioMaxCut) ; 
  printf("Jet TPC Ratio min cut     : %f\n", fJetTPCRatioMinCut) ;
  printf("Fast recons       : %d\n", fOptFast);
  printf("Inv Mass max cut  : %f\n", fInvMassMaxCut) ; 
  printf("Inv Mass min cut  : %f\n", fInvMassMinCut) ; 
 
} 

//__________________________________________________________________________
TChain * AliPHOSGammaJet::ReadESDfromdisk(const UInt_t eventsToRead, 
			 const TString dirName, 
			 const TString esdTreeName, 
			 const char *  pattern) 
{
  // Reads ESDs from Disk
 TChain *  rv = 0; 

  AliInfo( Form("\nReading files in %s \nESD tree name is %s \nReading %d events", 
		dirName.Data(), esdTreeName.Data(), eventsToRead) ) ;
  
  // create a TChain of all the files 
  TChain * cESDTree = new TChain(esdTreeName) ; 

  // read from the directory file until the require number of events are collected
  void * from = gSystem->OpenDirectory(dirName) ;
  if (!from) {
    AliError( Form("Directory %s does not exist") ) ;
    rv = 0 ;
  }
  else{ // reading file names from directory
    const char * subdir ; 
    // search all subdirectories witch matching pattern
    while( (subdir = gSystem->GetDirEntry(from))  && 
	   (cESDTree->GetEntries() < eventsToRead)) {
      if ( strstr(subdir, pattern) != 0 ) { 
	char file[200] ; 
        sprintf(file, "%s%s/AliESDs.root", dirName.Data(), subdir); 	
	AliInfo( Form("Adding %s\n", file) );
	cESDTree->Add(file) ;
      }
    } // while file
  
    AliInfo( Form(" %d events read", cESDTree->GetEntriesFast()) ) ;
    rv = cESDTree ; 
    
  } // reading file names from directory
  return rv ; 
}

//__________________________________________________________________________
TChain * AliPHOSGammaJet::ReadESD(const UInt_t eventsToRead,
		  const TString dirName, 
		  const TString esdTreeName, 
		  const char *  pattern)  
{
  // Read AliESDs files and return a Chain of events
 
  if ( dirName == "" ) {
    AliError("Give the name of the DIR where to find files") ; 
    return 0 ; 
  }
  if ( esdTreeName == "" ) 
    return ReadESDfromdisk(eventsToRead, dirName) ;
  else if ( strcmp(pattern, "") == 0 )
    return ReadESDfromdisk(eventsToRead, dirName, esdTreeName) ;
  else 
    return ReadESDfromdisk(eventsToRead, dirName, esdTreeName, pattern) ;	    
}

//___________________________________________________________________
void AliPHOSGammaJet::SetJet(TParticle * part, Bool_t & b, Float_t cone, 
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

//____________________________________________________________________________
Double_t AliPHOSGammaJet::SigmaE(Double_t energy)
{
  // Calculates the energy dependent energy resolution
  
  Double_t rv = -1 ; 
  
  rv = TMath::Sqrt( TMath::Power(fResPara1/energy, 2) 
	       + TMath::Power(fResPara2/TMath::Sqrt(energy), 2) 
	       + TMath::Power(fResPara3, 2) ) ;  

  return rv * energy ; 
}

//____________________________________________________________________________
Double_t AliPHOSGammaJet::SigmaP(Double_t energy)
{
  // Calculates the energy dependent position resolution 

  Double_t sigma = TMath::Sqrt(TMath::Power(fPosParaA,2) + 
			       TMath::Power(fPosParaB,2) / energy) ; 


  return sigma   ; // in cm  
}

