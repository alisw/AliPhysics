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

//_________________________________________________________________________
// Algorythm class to analyze PHOS events. In this class we demostrate, 
// how to handle reconstructed objects with AliPHSOIndexToObject.
// As an example we propose sulotions for four most frequently used tasks:
//    DrawRecon(...) - to draw reconstructed objects in the PHOS plane,
//                     very usefull in the debuging
//    InvarianMass(...) - to calculate "REAL" and "MIXED" photon pairs 
//                        invariant mass distributions
//    EnergyResoluition(...) -\ Energy and position resolutions of the 
//    PositionResolution(...)-/ reconstructed photons
//    Contamination(...) - calculates contamination of the photon spectrum and 
//                         pobability of reconstruction of several primaries as 
//                         kGAMMA,kELECTRON etc.
//
//    User Case:
//    root [0] AliPHOSAnalyze * a = new AliPHOSAnalyze("galice.root")
//                    // set the file you want to analyse
//    root [1] a->DrawRecon(1,3)
//                    // plot RecObjects, made in event 1, PHOS module 3 
//    root [2] a->DrawRecon(1,3,"PHOSRP","another PID")
//                    // plot RecObjets made in the event 1, PHOS module 3,
//                    // produced in the another reconstruction pass,
//                    // which produced PHOS RecParticles ("PHOSRP") with 
//                    // title "another PID".
//    root [3] a->InvariantMass()
//                    // Calculates "REAL" and "MIXED" invariant mass 
//                    // distributions of kGAMMA and (kGAMMA+kNEUTRALEM)
//                    // and APPENDS this to the file "invmass.root"
//    root [4] a->PositionResolution()
//                    // calculates two dimentional histos: energy of the primary
//                    // photon vs distance betwin incedence point and reconstructed 
//                    // poisition. One can analyse the produced file position.root 
//                    // with macro PhotonPosition.C
//    root [5] a->EnergyResolution()
//                    // calculates two dimentional histos: energy of the primary
//                    // photon vs energy of the reconstructed particle. One can 
//                    // analyse the produced file energy.root 
//                    // with macro PhotonEnergy.C
//    root [6] a->Contamination()
//                    // fills spectra of primary photons and several kinds of 
//                    // reconstructed particles, so that analyzing them one can 
//                    // estimate conatmination, efficiency of registration etc.
//*--
//*-- Author: Dmitri Peressounko (SUBATECH & RRC Kurchatov Institute)
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TMath.h"
#include "TROOT.h"
#include "TFolder.h"

// --- Standard library ---

#include <iomanip.h>

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliPHOSv1.h"
#include "AliPHOSAnalyze.h"
#include "AliPHOSDigit.h"
#include "AliPHOSSDigitizer.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSRecParticle.h"
#include "AliPHOSIndexToObject.h"
#include "AliPHOSCpvRecPoint.h"
#include "AliPHOSPpsdRecPoint.h"


ClassImp(AliPHOSAnalyze)

//____________________________________________________________________________
  AliPHOSAnalyze::AliPHOSAnalyze()
{
  // default ctor (useless)
  fCorrection = 1.2 ;  //Value calculated for default parameters of reconstruction  
}

//____________________________________________________________________________
AliPHOSAnalyze::AliPHOSAnalyze(Text_t * fileName)
{
  // ctor: analyze events from root file "name"
  ffileName = fileName ;
  fCorrection = 1.05 ;  //Value calculated for default parameters of reconstruction  
  fObjGetter = 0 ;  // should be instantiated
 
}

//____________________________________________________________________________
AliPHOSAnalyze::AliPHOSAnalyze(const AliPHOSAnalyze & ana)
{
  // copy ctor
  ( (AliPHOSAnalyze &)ana ).Copy(*this) ;
}

//____________________________________________________________________________
AliPHOSAnalyze::~AliPHOSAnalyze()
{
  // dtor

}
//____________________________________________________________________________
void AliPHOSAnalyze::DrawRecon(Int_t Nevent,Int_t Nmod,const char * branchName,const char* branchTitle){
  //Draws pimary particles and reconstructed 
  //digits, RecPoints, RecPartices etc 
  //for event Nevent in the module Nmod.
  
  TH2F * digitOccupancy  = new TH2F("digitOccupancy","EMC digits", 64,-71.,71.,64,-71.,71.);
  TH2F * sdigitOccupancy = new TH2F("sdigitOccupancy","EMC sdigits", 64,-71.,71.,64,-71.,71.);
  TH2F * emcOccupancy    = new TH2F("emcOccupancy","EMC RecPoints",64,-71.,71.,64,-71.,71.);
  TH2F * ppsdUp          = new TH2F("ppsdUp","PPSD Up digits",     128,-71.,71.,128,-71.,71.) ;
  TH2F * ppsdUpCl        = new TH2F("ppsdUpCl","PPSD Up RecPoints",128,-71.,71.,128,-71.,71.) ;
  TH2F * ppsdLow         = new TH2F("ppsdLow","PPSD Low digits",     128,-71.,71.,128,-71.,71.) ;
  TH2F * ppsdLowCl       = new TH2F("ppsdLowCl","PPSD Low RecPoints",128,-71.,71.,128,-71.,71.) ;
  TH2F * nbar            = new TH2F("nbar","Primary nbar",    64,-71.,71.,64,-71.,71.);
  TH2F * phot            = new TH2F("phot","Primary Photon",  64,-71.,71.,64,-71.,71.);
  TH2F * charg           = new TH2F("charg","Primary charged",64,-71.,71.,64,-71.,71.);
  TH2F * recPhot         = new TH2F("recPhot","RecParticles with primary Photon",64,-71.,71.,64,-71.,71.);
  TH2F * recNbar         = new TH2F("recNbar","RecParticles with primary Nbar",  64,-71.,71.,64,-71.,71.);
  
  //========== Create IndexToObjectGetter
  fObjGetter = AliPHOSIndexToObject::GetInstance(ffileName.Data(),branchName,branchTitle) ;
  fObjGetter->GetEvent(Nevent);
  
  //Plot Primary Particles
  TParticle * primary ;
  Int_t iPrimary ;
  for ( iPrimary = 0 ; iPrimary < fObjGetter->GimeNPrimaries() ; iPrimary++)
    {
      primary = fObjGetter->GimePrimary(iPrimary) ;
      Int_t primaryType = primary->GetPdgCode() ;
      if( (primaryType == 211)||(primaryType == -211)||(primaryType == 2212)||(primaryType == -2212) ) {
        Int_t moduleNumber ;
        Double_t primX, primZ ;
        fGeom->ImpactOnEmc(primary->Theta(), primary->Phi(), moduleNumber, primX, primZ) ;
        if(moduleNumber==Nmod)
          charg->Fill(primZ,primX,primary->Energy()) ;
      }
      if( primaryType == 22 ) {
        Int_t moduleNumber ;
        Double_t primX, primZ ;
        fGeom->ImpactOnEmc(primary->Theta(), primary->Phi(), moduleNumber, primX, primZ) ;
        if(moduleNumber==Nmod)
          phot->Fill(primZ,primX,primary->Energy()) ;
      }
      else{
        if( primaryType == -2112 ) {
          Int_t moduleNumber ;
          Double_t primX, primZ ;
          fGeom->ImpactOnEmc(primary->Theta(), primary->Phi(), moduleNumber, primX, primZ) ;
          if(moduleNumber==Nmod)
            nbar->Fill(primZ,primX,primary->Energy()) ;
        }
      }
    }  
  
  Int_t iSDigit ;
  AliPHOSDigit * sdigit ;
  
  for(iSDigit = 0; iSDigit < fObjGetter->GimeNSDigits(); iSDigit++)
      {
	sdigit = fObjGetter->GimeSDigit(iSDigit) ;
	Int_t relid[4];
	fGeom->AbsToRelNumbering(sdigit->GetId(), relid) ;
	Float_t x,z ;
	fGeom->RelPosInModule(relid,x,z) ;
  	Float_t e = fObjGetter->GimeSDigitizer()->Calibrate(sdigit->GetAmp()) ;
	if(relid[0]==Nmod){
	  if(relid[1]==0)  //EMC
	    sdigitOccupancy->Fill(x,z,e) ;
	  if((relid[1]>0)&&(relid[1]<17))
	    ppsdUp->Fill(x,z,e) ;
	  if(relid[1]>16)
	    ppsdLow->Fill(x,z,e) ;
	}
      }
  
  //Plot digits
  Int_t iDigit ;
  AliPHOSDigit * digit ;
  for(iDigit = 0; iDigit < fObjGetter->GimeNDigits(); iDigit++)
    {
      digit = fObjGetter->GimeDigit(iDigit) ;
      Int_t relid[4];
      fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
      Float_t x,z ;
      fGeom->RelPosInModule(relid,x,z) ;
      Float_t e = fObjGetter->GimeSDigitizer()->Calibrate(digit->GetAmp()) ;
      if(relid[0]==Nmod){
	if(relid[1]==0)  //EMC
	  digitOccupancy->Fill(x,z,e) ;
	if((relid[1]>0)&&(relid[1]<17))
	  ppsdUp->Fill(x,z,e) ;
	if(relid[1]>16)
	  ppsdLow->Fill(x,z,e) ;
      }
    }


  //Plot RecPoints
  Int_t irecp ;
  TVector3 pos ;
  
  for(irecp = 0; irecp < fObjGetter->GimeNEmcRecPoints() ; irecp ++){
    AliPHOSEmcRecPoint * emc= fObjGetter->GimeEmcRecPoint(irecp) ;
    if(emc->GetPHOSMod()==Nmod){
      emc->GetLocalPosition(pos) ;
      emcOccupancy->Fill(pos.X(),pos.Z(),emc->GetEnergy());
    }
  }


  for(irecp = 0; irecp < fObjGetter->GimeNCpvRecPoints() ; irecp ++){
    AliPHOSRecPoint * cpv = fObjGetter->GimeCpvRecPoint(irecp) ;
    if((strcmp(cpv->ClassName(),"AliPHOSPpsdRecPoint" )) == 0){  // PPSD Rec Point
      AliPHOSPpsdRecPoint * ppsd = (AliPHOSPpsdRecPoint*) cpv ;
      ppsd->GetLocalPosition(pos) ;
      if(ppsd->GetPHOSMod()==Nmod){
	ppsd->GetLocalPosition(pos) ;
	if(ppsd->GetUp())
	  ppsdUpCl->Fill(pos.X(),pos.Z(),ppsd->GetEnergy());
	else
	  ppsdLowCl->Fill(pos.X(),pos.Z(),ppsd->GetEnergy());
      }
    }
    else{
      AliPHOSCpvRecPoint * cpv1 = (AliPHOSCpvRecPoint*) cpv ;
      if(cpv1->GetPHOSMod()==Nmod){
	cpv1->GetLocalPosition(pos) ;
	ppsdUpCl->Fill(pos.X(),pos.Z(),cpv1->GetEnergy());
      }
    }
  }
  
  
  //Plot RecParticles
  AliPHOSRecParticle * recParticle ;
  Int_t iRecParticle ;
  for(iRecParticle = 0; iRecParticle < fObjGetter->GimeNRecParticles() ;iRecParticle++ )
    {
      recParticle =  fObjGetter->GimeRecParticle(iRecParticle) ;
      Int_t moduleNumberRec ;
      Double_t recX, recZ ;
      fGeom->ImpactOnEmc(recParticle->Theta(), recParticle->Phi(), moduleNumberRec, recX, recZ) ;
      if(moduleNumberRec == Nmod){

	Double_t minDistance = 5. ;
	Int_t closestPrimary = -1 ;
	

	//extract list of primaries: it is stored at EMC RecPoints
	Int_t emcIndex = fObjGetter->GimeTrackSegment(recParticle->GetPHOSTSIndex())->GetEmcIndex() ;
	Int_t numberofprimaries ;
	Int_t * listofprimaries  = fObjGetter->GimeEmcRecPoint(emcIndex)->GetPrimaries(numberofprimaries)  ;
	Int_t index ;
	TParticle * primary ;
	Double_t distance = minDistance ;
	
	for ( index = 0 ; index < numberofprimaries ; index++){
	  primary = fObjGetter->GimePrimary(listofprimaries[index]) ;
	  Int_t moduleNumber ;
	  Double_t primX, primZ ;
	  fGeom->ImpactOnEmc(primary->Theta(), primary->Phi(), moduleNumber, primX, primZ) ;
	  if(moduleNumberRec == moduleNumber)
	    distance = TMath::Sqrt((recX-primX)*(recX-primX)+(recZ-primZ)*(recZ-primZ) ) ;
	  if(minDistance > distance)
	    {
	      minDistance = distance ;
	      closestPrimary = listofprimaries[index] ;
	    }
	}
	
	if(closestPrimary >=0 ){
	  
	  Int_t primaryType = fObjGetter->GimePrimary(closestPrimary)->GetPdgCode() ;
	  
	  if(primaryType==22)
	    recPhot->Fill(recZ,recX,recParticle->Energy()) ;
	  else
	    if(primaryType==-2112)
	      recNbar->Fill(recZ,recX,recParticle->Energy()) ; 
	}
      }
    }

  
  //Plot made histograms
  digitOccupancy->Draw("box") ;
  sdigitOccupancy->SetLineColor(5) ;
  sdigitOccupancy->Draw("box") ;
  emcOccupancy->SetLineColor(2) ;
  emcOccupancy->Draw("boxsame") ;
  ppsdUp->SetLineColor(3) ;
  ppsdUp->Draw("boxsame") ;
  ppsdLow->SetLineColor(4) ;
  ppsdLow->Draw("boxsame") ;
  phot->SetLineColor(8) ;
  phot->Draw("boxsame") ;
  nbar->SetLineColor(6) ;
  nbar->Draw("boxsame") ;
  
}
//____________________________________________________________________________
void AliPHOSAnalyze::Ls(){
  //lists branches and titles of PHOS-related branches of TreeR, TreeD, TreeS
  
  if(fObjGetter == 0)
    fObjGetter = AliPHOSIndexToObject::GetInstance(ffileName.Data()) ;

  Int_t ibranch;
  TObjArray * branches; 
  
  branches = gAlice->TreeS()->GetListOfBranches() ;
 
  cout << "TreeS: " << endl ;
  for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
    TBranch * branch=(TBranch *) branches->At(ibranch) ;
    if(strstr(branch->GetName(),"PHOS") )
      cout << "       " << branch->GetName() << "     " << branch->GetTitle() << endl ;
  }
  cout << endl ;
 
  branches = gAlice->TreeD()->GetListOfBranches() ;
  
  cout << "TreeD: " << endl ;
  for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
    TBranch * branch=(TBranch *) branches->At(ibranch) ;
    if(strstr(branch->GetName(),"PHOS") )
      cout << "       " << branch->GetName() << "     " << branch->GetTitle() << endl ;
  }
  cout << endl ;
  

  branches = gAlice->TreeR()->GetListOfBranches() ;
  
  cout << "TreeR: " << endl ;
  for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
    TBranch * branch=(TBranch *) branches->At(ibranch) ;
    if(strstr(branch->GetName(),"PHOS") )
      cout << "       " << branch->GetName() << "     " << branch->GetTitle() << endl ;
  }
  cout << endl ;


}
//____________________________________________________________________________
 void AliPHOSAnalyze::InvariantMass(const char* branchTitle)    
{
  // Calculates Real and Mixed invariant mass distributions
  fObjGetter = AliPHOSIndexToObject::GetInstance(ffileName.Data(),"PHOSRP",branchTitle) ;

  Int_t nMixedEvents = 4 ; //# of events used for calculation of 'mixed' distribution 

  //opening file
  TFile * mfile = new TFile("invmass.root","update");
  
  //========== Reading /Booking Histograms
  TH2D * hRealEM = 0 ;
  hRealEM = (TH2D*) mfile->Get("hRealEM") ;
  if(hRealEM == 0) 
    hRealEM = new TH2D("hRealEM",   "Real for EM particles",      250,0.,1.,40,0.,4.) ;
  TH2D * hRealPhot = 0 ;

  hRealPhot = (TH2D*)mfile->Get("hRealPhot");
  if(hRealPhot == 0)
    hRealPhot = new TH2D("hRealPhot", "Real for kPhoton particles", 250,0.,1.,40,0.,4.) ;

  TH2D * hMixedEM = 0 ;
  hMixedEM = (TH2D*) mfile->Get("hMixedEM") ;
  if(hMixedEM == 0)
    hMixedEM = new TH2D("hMixedEM",  "Mixed for EM particles",     250,0.,1.,40,0.,4.) ;

  TH2D * hMixedPhot = 0 ;
  hMixedPhot = (TH2D*) mfile->Get("hMixedPhot") ;
  if(hMixedPhot == 0)
    hMixedPhot = new TH2D("hMixedPhot","Mixed for kPhoton particles",250,0.,1.,40,0.,4.) ;
  

  //reading event and copyng it to TConesArray of all photons

  TClonesArray * allRecParticleList  = new TClonesArray("AliPHOSRecParticle", 1000) ;
  Int_t * nRecParticles = new Int_t[nMixedEvents] ; // to mark boundaries of each event in the total list
  for(Int_t index = 0; index < nMixedEvents; index ++)
    nRecParticles[index] = 0 ;
  Int_t iRecPhot = 0 ;                // number of EM particles in total list
  
  //scan over all events
  Int_t event ;
  for(event = 0; event < fObjGetter->GetMaxEvent(); event++  ){
    
    fObjGetter->GetEvent(event);
    
    //copy EM RecParticles to the "total" list        
    AliPHOSRecParticle * recParticle ;
    Int_t iRecParticle ;
    for(iRecParticle = 0; iRecParticle < fObjGetter->GimeNRecParticles()  ;iRecParticle++ )
      {
	recParticle = fObjGetter->GimeRecParticle(iRecParticle) ;
	if((recParticle->GetType() == AliPHOSFastRecParticle::kGAMMA)||
	   (recParticle->GetType() == AliPHOSFastRecParticle::kNEUTRALEM))
	  new( (*allRecParticleList)[iRecPhot++] ) AliPHOSRecParticle(*recParticle) ;
      }
    
    Int_t mevent = event%nMixedEvents ; //event number in the "mixed" cicle
    nRecParticles[mevent] = iRecPhot-1 ;  

    //check, if it is time to calculate invariant mass?
    if((mevent == 0) && (event +1 == fObjGetter->GetMaxEvent())){
      
      //calculate invariant mass:
      Int_t irp1,irp2 ;
      Int_t nCurEvent = 0 ;
      
      for(irp1 = 0; irp1 < allRecParticleList->GetEntries()-1; irp1++){
	AliPHOSRecParticle * rp1 = (AliPHOSRecParticle *)allRecParticleList->At(irp1) ;
	
	for(irp2 = irp1+1; irp2 < allRecParticleList->GetEntries(); irp2++){
	  AliPHOSRecParticle * rp2 = (AliPHOSRecParticle *)allRecParticleList->At(irp2) ;
	  
	  Double_t invMass ;
	  invMass = (rp1->Energy()+rp2->Energy())*(rp1->Energy()+rp2->Energy())-
	    (rp1->Px()+rp2->Px())*(rp1->Px()+rp2->Px())-
	    (rp1->Py()+rp2->Py())*(rp1->Py()+rp2->Py())-
	    (rp1->Pz()+rp2->Pz())*(rp1->Pz()+rp2->Pz()) ;
	  
	  if(invMass> 0)
	    invMass = TMath::Sqrt(invMass);
	  
	  Double_t pt ; 
	  pt = TMath::Sqrt((rp1->Px()+rp2->Px() )*( rp1->Px()+rp2->Px() ) + 
			   (rp1->Py()+rp2->Py() )*( rp1->Py()+rp2->Py() ) );
	  
	  if(irp1 > nRecParticles[nCurEvent])
	    nCurEvent++;
	  
	  if(irp2 <= nRecParticles[nCurEvent]){ //'Real' event
	    hRealEM->Fill(invMass,pt);
	    if((rp1->GetType() == AliPHOSFastRecParticle::kGAMMA)&&
	       (rp2->GetType() == AliPHOSFastRecParticle::kGAMMA) )
	      hRealPhot->Fill(invMass,pt);
	  }
	  else{
	    hMixedEM->Fill(invMass,pt);
	    if((rp1->GetType() == AliPHOSFastRecParticle::kGAMMA)&&
	       (rp2->GetType() == AliPHOSFastRecParticle::kGAMMA) )
	      hMixedPhot->Fill(invMass,pt);
	  } //real-mixed
	  
      } //loop over second rp
      }//loop over first rp
      
      //Make some cleanings 
      for(Int_t index = 0; index < nMixedEvents; index ++)
	nRecParticles[index] = 0 ;
      iRecPhot = 0 ;              
      allRecParticleList->Clear() ;
      
    }
  }
  delete allRecParticleList ;
  
  //writing output
  mfile->cd();
  
  hRealEM->Write(0,kOverwrite) ;
  hRealPhot->Write(0,kOverwrite) ;
  hMixedEM->Write(0,kOverwrite) ;
  hMixedPhot->Write(0,kOverwrite) ;
  
  mfile->Write();
  mfile->Close();
  delete mfile ;
  delete nRecParticles;

}

//____________________________________________________________________________
 void AliPHOSAnalyze::EnergyResolution(const char * branchTitle)    
{
  //fills two dimentional histo: energy of primary vs. energy of reconstructed

  TH2F * hAllEnergy = 0 ;  //all reconstructed with primary photon   
  TH2F * hPhotEnergy= 0 ;  //kGamma  with primary photon   
  TH2F * hEMEnergy  = 0 ;  //electromagnetic with primary photon   

  //opening file and reading histograms if any
  TFile * efile = new TFile("energy.root","update");

  hAllEnergy = (TH2F*)efile->Get("hAllEnergy") ;
  if(hAllEnergy == 0)
    hAllEnergy  = new TH2F("hAllEnergy",  "Energy of any RP with primary photon",100, 0., 5., 100, 0., 5.);

  hPhotEnergy =(TH2F*) efile->Get("hPhotEnergy") ;
  if(hPhotEnergy == 0)
    hPhotEnergy = new TH2F("hPhotEnergy", "Energy of kGAMMA with primary photon",100, 0., 5., 100, 0., 5.);

  hEMEnergy =(TH2F*) efile->Get("hEMEnergy");
  if(hEMEnergy == 0)
    hEMEnergy   = new TH2F("hEMEnergy",   "Energy of EM with primary photon",    100, 0., 5., 100, 0., 5.);


  fObjGetter = AliPHOSIndexToObject::GetInstance(ffileName.Data(),"PHOSRP",branchTitle) ;
  fGeom = fObjGetter->GetPHOSGeometry() ; 

  Int_t ievent;
  for ( ievent=0; ievent < fObjGetter->GetMaxEvent() ; ievent++){

    //read the current event
    fObjGetter->GetEvent(ievent) ;

    AliPHOSRecParticle * recParticle ;
    Int_t iRecParticle ;
    for(iRecParticle = 0; iRecParticle < fObjGetter->GimeNRecParticles() ;iRecParticle++ ){
      recParticle = fObjGetter->GimeRecParticle(iRecParticle) ;
      
      //find the closest primary
      Int_t moduleNumberRec ;
      Double_t recX, recZ ;
      fGeom->ImpactOnEmc(recParticle->Theta(), recParticle->Phi(), moduleNumberRec, recX, recZ) ;
      
      Double_t minDistance  = 100. ;
      Int_t closestPrimary = -1 ;
      
      //extract list of primaries: it is stored at EMC RecPoints
      Int_t emcIndex = fObjGetter->GimeTrackSegment(recParticle->GetPHOSTSIndex())->GetEmcIndex() ;
      Int_t numberofprimaries ;
      Int_t * listofprimaries  = fObjGetter->GimeEmcRecPoint(emcIndex)->GetPrimaries(numberofprimaries)  ;

      Int_t index ;
      TParticle * primary ;
      Double_t distance = minDistance ;
      Double_t dX, dZ; 
      Double_t dXmin = 0.; 
      Double_t dZmin = 0. ;
      for ( index = 0 ; index < numberofprimaries ; index++){
	primary = fObjGetter->GimePrimary(listofprimaries[index]) ;
	Int_t moduleNumber ;
	Double_t primX, primZ ;
	fGeom->ImpactOnEmc(primary->Theta(), primary->Phi(), moduleNumber, primX, primZ) ;
	if(moduleNumberRec == moduleNumber) {
	  dX = recX - primX;
	  dZ = recZ - primZ;
	  distance = TMath::Sqrt(dX*dX + dZ*dZ) ;
	  if(minDistance > distance) {
	    minDistance = distance ;
	    dXmin = dX;
	    dZmin = dZ;
	    closestPrimary = listofprimaries[index] ;
	  }
	}
      }

      //if found primary, fill histograms
      if(closestPrimary >=0 ){
	TParticle * primary = fObjGetter->GimePrimary(closestPrimary) ;
	if(primary->GetPdgCode() == 22){
	  hAllEnergy->Fill(primary->Energy(), recParticle->Energy()) ;
	  if(recParticle->GetType() == AliPHOSFastRecParticle::kGAMMA){
	    hPhotEnergy->Fill(primary->Energy(), recParticle->Energy() ) ; 
	    hEMEnergy->Fill(primary->Energy(), recParticle->Energy() ) ; 
	  }
	  else
	    if(recParticle->GetType() == AliPHOSFastRecParticle::kNEUTRALEM)
	      hEMEnergy->Fill(primary->Energy(), recParticle->Energy() ) ; 
	}
      }
    }
  }

  //write filled histograms
  efile->cd() ;
  hAllEnergy->Write(0,kOverwrite) ;
  hPhotEnergy->Write(0,kOverwrite) ;
  hEMEnergy->Write(0,kOverwrite)  ;
  //  efile->Write() ;
  efile->Close() ;
  delete efile ;

}
//____________________________________________________________________________
void AliPHOSAnalyze::PositionResolution(const char * branchTitle)    
{
  //fills two dimentional histo: energy vs. primary - reconstructed distance  



  TH2F * hAllPosition  = 0;     // Position of any RP with primary photon
  TH2F * hPhotPosition = 0;    // Position of kGAMMA with primary photon
  TH2F * hEMPosition   = 0;      // Position of EM with primary photon

  TH1F * hAllPositionX = 0;    // X-Position Resolution of photons with photon primary
  TH1F * hAllPositionZ = 0;    // Z-Position Resolution of photons with photon primary


  //opening file and reading histograms if any
  TFile * pfile = new TFile("position.root","update");

  hAllPosition = (TH2F*)pfile->Get("hAllPosition");
  if(hAllPosition == 0)
    hAllPosition  = new TH2F("hAllPosition",  
			     "Position of any RP with primary photon",100, 0., 5., 100, 0., 5.);
  hPhotPosition= (TH2F*)pfile->Get("hPhotPosition");
  if(hPhotPosition == 0)
    hPhotPosition = new TH2F("hPhotPosition", 
			     "Position of kGAMMA with primary photon",100, 0., 5., 100, 0., 5.);
  hEMPosition= (TH2F*)pfile->Get("hEMPosition") ;
  if(hEMPosition == 0)
    hEMPosition   = new TH2F("hEMPosition",   
			     "Position of EM with primary photon",    100, 0., 5., 100, 0., 5.);
  hAllPositionX = (TH1F*)pfile->Get("hAllPositionX") ;			   
  if(hAllPositionX == 0)
    hAllPositionX = new TH1F("hAllPositionX", 
			     "Delta X of any RP with primary photon",100, -2., 2.);
  hAllPositionZ =(TH1F*) pfile->Get("hAllPositionZ") ;
  if(hAllPositionZ == 0)
    hAllPositionZ = new TH1F("hAllPositionZ", 
			     "Delta X of any RP with primary photon",100, -2., 2.);


  fObjGetter = AliPHOSIndexToObject::GetInstance(ffileName.Data(),"PHOSRP",branchTitle) ;


  fGeom = AliPHOSGeometry::GetInstance(((AliPHOS*)gAlice->GetModule("PHOS"))->GetGeometry()->GetName(),
				       ((AliPHOS*)gAlice->GetModule("PHOS"))->GetGeometry()->GetTitle() );
  
  Int_t ievent;
  for ( ievent=0; ievent < fObjGetter->GetMaxEvent() ; ievent++){

    //read the current event
    fObjGetter->GetEvent(ievent) ;
    
    AliPHOSRecParticle * recParticle ;
    Int_t iRecParticle ;
    for(iRecParticle = 0; iRecParticle < fObjGetter->GimeNRecParticles() ;iRecParticle++ ){
      recParticle = fObjGetter->GimeRecParticle(iRecParticle) ;
      
      //find the closest primary
      Int_t moduleNumberRec ;
      Double_t recX, recZ ;
      fGeom->ImpactOnEmc(recParticle->Theta(), recParticle->Phi(), moduleNumberRec, recX, recZ) ;
      
      Double_t minDistance  = 100. ;
      Int_t closestPrimary = -1 ;
      
      //extract list of primaries: it is stored at EMC RecPoints
      Int_t emcIndex = fObjGetter->GimeTrackSegment(recParticle->GetPHOSTSIndex())->GetEmcIndex() ;
      Int_t numberofprimaries ;
      Int_t * listofprimaries  = fObjGetter->GimeEmcRecPoint(emcIndex)->GetPrimaries(numberofprimaries)  ;

      Int_t index ;
      TParticle * primary ;
      Double_t distance = minDistance ;
      Double_t dX = 1000; // incredible number
      Double_t dZ = 1000; // for the case if no primary will be found
      Double_t dXmin = 0.; 
      Double_t dZmin = 0. ;
      for ( index = 0 ; index < numberofprimaries ; index++){
	primary = fObjGetter->GimePrimary(listofprimaries[index]) ;
	Int_t moduleNumber ;
	Double_t primX, primZ ;
	fGeom->ImpactOnEmc(primary->Theta(), primary->Phi(), moduleNumber, primX, primZ) ;
	if(moduleNumberRec == moduleNumber) {
	  dX = recX - primX;
	  dZ = recZ - primZ;
	  distance = TMath::Sqrt(dX*dX + dZ*dZ) ;
	  if(minDistance > distance) {
	    minDistance = distance ;
	    dXmin = dX;
	    dZmin = dZ;
	    closestPrimary = listofprimaries[index] ;
	  }
	}
      }
      
      //if found primary, fill histograms
      if(closestPrimary >=0 ){
	TParticle * primary = fObjGetter->GimePrimary(closestPrimary) ;
	if(primary->GetPdgCode() == 22){
	  hAllPosition->Fill(primary->Energy(), minDistance) ;
	  hAllPositionX->Fill(primary->Energy(), dX) ;
	  hAllPositionZ->Fill(primary->Energy(), dZ) ;
	  if(recParticle->GetType() == AliPHOSFastRecParticle::kGAMMA){
	    hPhotPosition->Fill(primary->Energy(), minDistance ) ; 
	    hEMPosition->Fill(primary->Energy(), minDistance ) ; 
	  }
	  else
	    if(recParticle->GetType() == AliPHOSFastRecParticle::kNEUTRALEM)
	      hEMPosition->Fill(primary->Energy(), minDistance ) ; 
	}
      }
    }
  }
  
  //Write output histgrams
  pfile->cd() ;
  hAllPosition->Write(0,kOverwrite) ;
  hAllPositionX->Write(0,kOverwrite) ;
  hAllPositionZ->Write(0,kOverwrite) ;
  hPhotPosition->Write(0,kOverwrite) ;
  hEMPosition->Write(0,kOverwrite) ;
  pfile->Write() ;
  pfile->Close() ;
  delete pfile ;


}
//____________________________________________________________________________
void AliPHOSAnalyze::Contamination(const char* RecPointsTitle){
// fills spectra of primary photons and several kinds of 
// reconstructed particles, so that analyzing them one can 
// estimate conatmination, efficiency of registration etc.

  //define several general histograms
  TH1F * hPrimary = 0;   //spectrum (P_t distribution) of primary photons         
  TH1F * hAllRP   = 0;   //spectrum of all RecParticles in PHOS
  TH1F * hPhot    = 0;   //spectrum of kGAMMA RecParticles
  TH1F * hPPSD    = 0;   //spectrum of all RecParticles with PPSD signal
  TH1F * hShape   = 0;   //spectrum of all EM RecParticles
  TH1F * hVeto    = 0;   //spectrum of all neutral RecParticles

  //Now separate histograms in accoradance with primary
  //primary - photon
  TH1F * hPhotReg = 0;   //Registeres as photon
  TH1F * hPhotEM  = 0;   //Registered as EM       
  TH1F * hPhotPPSD= 0;   //Registered as RecParticle with PPSD signal        

  //primary - n
  TH1F * hNReg = 0;   //Registeres as photon          
  TH1F * hNEM  = 0;   //Registered as EM            
  TH1F * hNPPSD= 0;   //Registered as RecParticle with PPSD signal           

  //primary - nBar
  TH1F * hNBarReg = 0;   //Registeres as photon
  TH1F * hNBarEM  = 0;   //Registered as EM          
  TH1F * hNBarPPSD= 0;   //Registered as RecParticle with PPSD signal             

  //primary - charged hadron (pBar excluded)
  TH1F * hChargedReg = 0;   //Registeres as photon  
  TH1F * hChargedEM  = 0;   //Registered as EM           
  TH1F * hChargedPPSD= 0;   //Registered as RecParticle with PPSD signal           

  //primary - pBar
  TH1F * hPbarReg = 0;   //Registeres as photon  
  TH1F * hPbarEM  = 0;   //Registered as EM 
  TH1F * hPbarPPSD= 0;   //Registered as RecParticle with PPSD signal   


  //Reading histograms from the file
  TFile * cfile = new TFile("contamination.root","update") ;

  //read general histograms
  hPrimary = (TH1F*) cfile->Get("hPrimary") ;
  if(hPrimary == 0)
    hPrimary= new TH1F("hPrimary", "Primary photon spectrum", 100, 0., 5.);
  hAllRP = (TH1F*)cfile->Get("hAllRP") ;
  if(hAllRP == 0)
    hAllRP = new TH1F("hAllRP","All Reconstructed particles", 100, 0., 5.);
  hPhot  = (TH1F*)cfile->Get("hPhot") ;
  if(hPhot == 0)
    hPhot = new TH1F("hPhot","All kGAMMA RecParticles",100, 0., 5.);
  hPPSD = (TH1F*)cfile->Get("hPPSD") ;
  if(hPPSD == 0)
    hPPSD  = new TH1F("hPPSD", "All PPSD Recparticles",    100, 0., 5.);
  hShape = (TH1F*) cfile->Get("hShape") ;
  if(hShape == 0)
    hShape = new TH1F("hShape","All particles with EM shower",100, 0., 5.);
  hVeto= (TH1F*)cfile->Get("hVeto") ;
  if(hVeto == 0) 
    hVeto  = new TH1F("hVeto", "All uncharged particles",      100, 0., 5.);


  //primary - photon
  hPhotReg = (TH1F*)cfile->Get("hPhotReg");
  if(hPhotReg == 0)
    hPhotReg   = new TH1F("hPhotReg","Photon registered as photon",100, 0., 5.);
  hPhotEM  =(TH1F*)cfile->Get("hPhotEM");
  if(hPhotEM== 0)
    hPhotEM   = new TH1F("hPhotEM",  "Photon registered as EM", 100, 0., 5.);
  hPhotPPSD= (TH1F*)cfile->Get("hPhotPPSD");
  if(hPhotPPSD== 0)
    hPhotPPSD   = new TH1F("hPhotPPSD","Photon registered as PPSD", 100, 0., 5.);

  //primary - n
  hNReg = (TH1F*)cfile->Get("hNReg");
  if(hNReg== 0)
   hNReg      = new TH1F("hNReg",   "N registered as photon",              100, 0., 5.);
  hNEM  = (TH1F*)cfile->Get("hNEM"); 
  if(hNEM== 0)
    hNEM      = new TH1F("hNEM",    "N registered as EM",      100, 0., 5.);
  hNPPSD=(TH1F*)cfile->Get("hNPPSD"); 
  if(hNPPSD== 0)
    hNPPSD      = new TH1F("hNPPSD","N registered as PPSD",      100, 0., 5.);

  //primary - nBar
  hNBarReg =(TH1F*)cfile->Get("hNBarReg");
  if(hNBarReg== 0)
   hNBarReg   = new TH1F("hNBarReg", "NBar registered as photon",           100, 0., 5.);
  hNBarEM  =(TH1F*)cfile->Get("hNBarEM"); 
  if(hNBarEM== 0)
    hNBarEM   = new TH1F("hNBarEM",  "NBar registered as EM",   100, 0., 5.);
  hNBarPPSD=(TH1F*)cfile->Get("hNBarPPSD");
  if(hNBarPPSD== 0)
    hNBarPPSD   = new TH1F("hNBarPPSD","NBar registered as PPSD",   100, 0., 5.);

  //primary - charged hadron (pBar excluded)
  hChargedReg = (TH1F*)cfile->Get("hChargedReg");
  if(hChargedReg== 0)
    hChargedReg= new TH1F("hChargedReg", "Charged hadron registered as photon",100, 0., 5.);
  hChargedEM  = (TH1F*)cfile->Get("hChargedEM"); 
  if(hChargedEM== 0)
    hChargedEM= new TH1F("hChargedEM","Charged registered as EM",100, 0., 5.);
  hChargedPPSD= (TH1F*)cfile->Get("hChargedPPSD");
  if(hChargedPPSD== 0)
    hChargedPPSD= new TH1F("hChargedPPSD","Charged registered as PPSD",100, 0., 5.);
  
  //primary - pBar
  hPbarReg = (TH1F*)cfile->Get("hPbarReg");
  if(hPbarReg== 0)
    hPbarReg= new TH1F("hPbarReg", "pBar registered as photon",100, 0., 5.);
  hPbarEM  = (TH1F*)cfile->Get("hPbarEM");
  if(hPbarEM== 0)
    hPbarEM= new TH1F("hPbarEM","Pbar registered as EM",100, 0., 5.);
  hPbarPPSD= (TH1F*)cfile->Get("hPbarPPSD");
  if(hPbarPPSD== 0)
    hPbarPPSD= new TH1F("hPbarPPSD","Pbar as PPSD",100, 0., 5.);
  

  //Now make some initializations

  Int_t counter[8][5] ;      //# of registered particles 
  Int_t i1,i2 ;
  for(i1 = 0; i1<8; i1++)
    for(i2 = 0; i2<5; i2++)
      counter[i1][i2] = 0 ;



  fObjGetter = AliPHOSIndexToObject::GetInstance(ffileName.Data(),"PHOSRP",RecPointsTitle) ;
  fGeom = AliPHOSGeometry::GetInstance(((AliPHOS*)gAlice->GetModule("PHOS"))->GetGeometry()->GetName(),
				       ((AliPHOS*)gAlice->GetModule("PHOS"))->GetGeometry()->GetTitle() );
  
  Int_t ievent;
  for ( ievent=0; ievent < fObjGetter->GetMaxEvent() ; ievent++){
    
    fObjGetter->GetEvent(ievent) ;
    
    //=========== Make spectrum of the primary photons
    TParticle * primary ;
    Int_t iPrimary ;
    for( iPrimary = 0 ; iPrimary < fObjGetter->GimeNPrimaries() ; iPrimary++){
      primary = fObjGetter->GimePrimary(iPrimary) ;
      Int_t primaryType = primary->GetPdgCode() ;
      if( primaryType == 22 ) {
	//check, if photons folls onto PHOS
	Int_t moduleNumber ;
	Double_t primX, primZ ;
	fGeom->ImpactOnEmc(primary->Theta(), primary->Phi(), moduleNumber, primX, primZ) ;
	if(moduleNumber)
	  hPrimary->Fill(primary->Energy()) ;
	
      }
      
    }
    //========== Now scan over RecParticles            
    AliPHOSRecParticle * recParticle ;
    Int_t iRecParticle ;
    for(iRecParticle = 0; iRecParticle < fObjGetter->GimeNRecParticles(); iRecParticle++ ){
      recParticle = fObjGetter->GimeRecParticle(iRecParticle) ;
      //fill histo spectrum of all RecParticles
      hAllRP->Fill(CorrectedEnergy(recParticle->Energy())) ;
      
      //==========find the closest primary	
      Int_t moduleNumberRec ;
      Double_t recX, recZ ;
      fGeom->ImpactOnEmc(recParticle->Theta(), recParticle->Phi(), moduleNumberRec, recX, recZ) ;
      
      Double_t minDistance  = 100. ;
      Int_t closestPrimary = -1 ;
      
      //extract list of primaries: it is stored at EMC RecPoints
      Int_t emcIndex = fObjGetter->GimeTrackSegment(recParticle->GetPHOSTSIndex())->GetEmcIndex() ;
      Int_t numberofprimaries ;
      Int_t * listofprimaries  = fObjGetter->GimeEmcRecPoint(emcIndex)->GetPrimaries(numberofprimaries)  ;
      Int_t index ;
      TParticle * primary ;
      Double_t distance = minDistance ;
      Double_t dX, dZ; 
      Double_t dXmin = 0.; 
      Double_t dZmin = 0. ;
      for ( index = 0 ; index < numberofprimaries ; index++){
	primary = fObjGetter->GimePrimary(listofprimaries[index]) ;
	Int_t moduleNumber ;
	Double_t primX, primZ ;
	fGeom->ImpactOnEmc(primary->Theta(), primary->Phi(), moduleNumber, primX, primZ) ;
	if(moduleNumberRec == moduleNumber) {
	  dX = recX - primX;
	  dZ = recZ - primZ;
	  distance = TMath::Sqrt(dX*dX + dZ*dZ) ;
	  if(minDistance > distance) {
	    minDistance = distance ;
	    dXmin = dX;
	    dZmin = dZ;
	    closestPrimary = listofprimaries[index] ;
	  }
	}
      }
      
      //===========define the "type" of closest primary
      if(closestPrimary >=0 ){
	Int_t primaryCode = -1;
	TParticle * primary = fObjGetter->GimePrimary(closestPrimary) ;
	Int_t primaryType = primary->GetPdgCode() ;
	if(primaryType == 22) // photon ?
	  primaryCode = 0 ;
	else
	  if(primaryType == 2112) // neutron
	    primaryCode = 1 ; 
	  else
	    if(primaryType == -2112) // Anti neutron
	      primaryCode = 2 ;
	    else
	      if(primaryType == -2122) //Anti proton
		primaryCode = 4 ;
	      else
		if(fObjGetter->GimePrimary(closestPrimary)->GetPDG()->Charge())
		  primaryCode = 3 ;
	
	//==========Now look at the type of RecParticle
	Float_t energy = CorrectedEnergy(recParticle->Energy()) ;
	if(recParticle->GetType() == AliPHOSFastRecParticle::kGAMMA){
	  hPhot->Fill(energy ) ; 	
	  switch(primaryCode){
	  case 0:
	    hPhotReg->Fill(energy ) ; 
	    break ;
	  case 1:
	    hNReg->Fill(energy ) ; 
	    break ;
	  case 2:
	    hNBarReg->Fill(energy ) ; 
	    break ;
	  case 3:
	    hChargedReg->Fill(energy ) ;
	    break ;
	  case 4:
	    hPbarReg->Fill(energy ) ;
	    break ;
	  default:
	    break ;
	  }
	}
	if((recParticle->GetType() == AliPHOSFastRecParticle::kGAMMA)||
	   (recParticle->GetType() == AliPHOSFastRecParticle::kGAMMAHA)){ //with PPSD signal
	  hPPSD->Fill(energy ) ; 
	  switch(primaryCode){
	  case 0:
	    hPhotPPSD->Fill(energy ) ; 
	    break ;
	  case 1:
	    hNPPSD->Fill(energy ) ; 
	    break ;
	  case 2:
	    hNBarPPSD->Fill(energy ) ;
	    break ;
	  case 3:
	    hChargedPPSD->Fill(energy ) ;
	    break ;
	  case 4:
	    hPbarPPSD->Fill(energy ) ;
	    break ;
	  default:
	    break ;
	  }
	}
	if((recParticle->GetType() == AliPHOSFastRecParticle::kGAMMA)||
	   (recParticle->GetType() == AliPHOSFastRecParticle::kELECTRON)||
	   (recParticle->GetType() == AliPHOSFastRecParticle::kNEUTRALEM)||
	   (recParticle->GetType() == AliPHOSFastRecParticle::kABSURDEM) ){ //with EM shower
	  hShape->Fill(energy ) ;
	  switch(primaryCode){
	  case 0:
	    hPhotEM->Fill(energy ) ; 
	    break ;
	  case 1:
	    hNEM->Fill(energy ) ; 
	    break ;
	  case 2:
	    hNBarEM->Fill(energy ) ; 
	    break ;
	  case 3:
	    hChargedEM->Fill(energy ) ; 
	    break ;
	  case 4:
	    hPbarEM->Fill(energy ) ; 
	    break ;
	  default:
	    break ;
	  }
	}
	
	if((recParticle->GetType() == AliPHOSFastRecParticle::kGAMMA)||
	   (recParticle->GetType() == AliPHOSFastRecParticle::kNEUTRALHA) ||
	   (recParticle->GetType() == AliPHOSFastRecParticle::kGAMMAHA) ||
	   (recParticle->GetType() == AliPHOSFastRecParticle::kNEUTRALEM) ) //nuetral
	  hVeto->Fill(energy ) ;
	
	//fill number of primaries identified as ...
	if(primaryCode >= 0) // Primary code defined
	  counter[recParticle->GetType()][primaryCode]++ ; 
	
      }
      
    } // no closest primary found
  }     
  
  
  //===================  SaveHistograms
  cfile->cd() ;
  hPrimary->Write(0,kOverwrite); 
  hAllRP->Write(0,kOverwrite);  
  hPhot->Write(0,kOverwrite);  
  hPPSD->Write(0,kOverwrite);  
  hShape->Write(0,kOverwrite); 
  hVeto->Write(0,kOverwrite);  
  hPhotReg->Write(0,kOverwrite); 
  hPhotEM->Write(0,kOverwrite);   
  hPhotPPSD->Write(0,kOverwrite); 
  hNReg ->Write(0,kOverwrite);  
  hNEM  ->Write(0,kOverwrite); 
  hNPPSD->Write(0,kOverwrite);  
  hNBarReg ->Write(0,kOverwrite); 
  hNBarEM  ->Write(0,kOverwrite); 
  hNBarPPSD->Write(0,kOverwrite); 
  hChargedReg ->Write(0,kOverwrite); 
  hChargedEM  ->Write(0,kOverwrite); 
  hChargedPPSD->Write(0,kOverwrite); 
  hPbarReg ->Write(0,kOverwrite); 
  hPbarEM  ->Write(0,kOverwrite); 
  hPbarPPSD->Write(0,kOverwrite); 
  
  cfile->Write(0,kOverwrite); 
  cfile->Close();
  delete cfile ;
 
  
  //print Final Table

  cout << "Resolutions: Analyzed " << fObjGetter->GetMaxEvent() << " event(s)" << endl ;
  cout << endl ;
  
  cout << "        Primary:    Photon  Neutron  Antineutron  Charged hadron   AntiProton" << endl ; 
  cout << "--------------------------------------------------------------------------------" << endl ;
  cout << "         kGAMMA: " 
       << setw(8) << counter[2][0] << setw(9)  << counter[2][1] << setw(13) << counter[2][2] 
       << setw(15)<< counter[2][3] << setw(13) << counter[2][4] << endl ;
  cout << "       kGAMMAHA: " 
       << setw(8) << counter[3][0] << setw(9)  << counter[3][1] << setw(13) << counter[3][2] 
       << setw(15)<< counter[3][3] << setw(13) << counter[3][4] << endl ;
  cout << "     kNEUTRALEM: " 
       << setw(8) << counter[0][0] << setw(9)  << counter[0][1] << setw(13) << counter[0][2] 
       << setw(15)<< counter[0][3] << setw(13) << counter[0][4] << endl ;
  cout << "     kNEUTRALHA: " 
       << setw(8) << counter[1][0] << setw(9)  << counter[1][1] << setw(13) << counter[1][2] 
       << setw(15)<< counter[1][3] << setw(13) << counter[1][4] << endl ;
  cout << "      kABSURDEM: " 
       << setw(8) << counter[4][0] << setw(9)  << counter[4][1] << setw(13) << counter[4][2] 
       << setw(15)<< counter[4][3] << setw(13) << counter[4][4] << endl ;
  cout << "      kABSURDHA: " 
       << setw(8) << counter[5][0] << setw(9)  << counter[5][1] << setw(13) << counter[5][2] 
       << setw(15)<< counter[5][3] << setw(13) << counter[5][4] << endl ;
  cout << "      kELECTRON: " 
       << setw(8) << counter[6][0] << setw(9)  << counter[6][1] << setw(13) << counter[6][2] 
       << setw(15)<< counter[6][3] << setw(13) << counter[6][4] << endl ; 
  cout << "     kCHARGEDHA: " 
       << setw(8) << counter[7][0] << setw(9)  << counter[7][1] << setw(13) << counter[7][2] 
       << setw(15)<< counter[7][3] << setw(13) << counter[7][4] << endl ;
  cout << "--------------------------------------------------------------------------------" << endl ;
      
  Int_t totalInd = 0 ;
  for(i1 = 0; i1<8; i1++)
    for(i2 = 0; i2<5; i2++)
      totalInd+=counter[i1][i2] ;
  cout << "Indentified particles: " << totalInd << endl ;
  
}



