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

//_________________________________________________________________________
// Algorythm class to analyze PHOS events
//*-- Y. Schutz :   SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TFile.h"
#include "TH1.h"
#include "TPad.h"
#include "TH2.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h" 

// --- Standard library ---

#include <iostream>
#include <cstdio>

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliPHOSAnalyze.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSPIDv1.h"
#include "AliPHOSReconstructioner.h"
#include "AliPHOSDigit.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSRecParticle.h"

ClassImp(AliPHOSAnalyze)


//____________________________________________________________________________
  AliPHOSAnalyze::AliPHOSAnalyze()
{
  // ctor
  
  fRootFile = 0 ; 
}

//____________________________________________________________________________
AliPHOSAnalyze::AliPHOSAnalyze(Text_t * name)
{
  // ctor
  
  Bool_t ok = OpenRootFile(name)  ; 
  if ( !ok ) {
    cout << " AliPHOSAnalyze > Error opening " << name << endl ; 
  }
  else { 
    gAlice = (AliRun*) fRootFile->Get("gAlice");
    fPHOS  = (AliPHOSv0 *)gAlice->GetDetector("PHOS") ;
    fGeom  = AliPHOSGeometry::GetInstance( fPHOS->GetGeometry()->GetName(), fPHOS->GetGeometry()->GetTitle() ) ;
    fEvt = -999 ; 
  }
}

//____________________________________________________________________________
AliPHOSAnalyze::~AliPHOSAnalyze()
{
  // dtor

  fRootFile->Close() ; 
  delete fRootFile ; 
  fRootFile = 0 ; 

  delete fPHOS ; 
  fPHOS = 0 ; 

  delete fClu ; 
  fClu = 0 ;

  delete fPID ; 
  fPID = 0 ;

  delete fRec ; 
  fRec = 0 ;

  delete fTrs ; 
  fTrs = 0 ;

}

//____________________________________________________________________________
void AliPHOSAnalyze::AnalyzeOneEvent(Int_t evt)
{
  Bool_t ok = Init(evt) ; 
  
  if ( ok ) {
    //=========== Get the number of entries in the Digits array
    
    Int_t nId = fPHOS->Digits()->GetEntries();    
    printf("AnalyzeOneEvent > Number of entries in the Digit array is %d \n",nId);
    
    //=========== Do the reconstruction
    
    cout << "AnalyzeOneEvent > Found  " << nId << "  digits in PHOS"   << endl ;  
    
   
    fPHOS->Reconstruction(fRec);  
    
    // =========== End of reconstruction
    
    cout << "AnalyzeOneEvent > event # " << fEvt << " processed" << endl ;   
  } // ok
  else
    cout << "AnalyzeOneEvent > filed to process event # " << evt << endl ;   

}

//____________________________________________________________________________
 void AliPHOSAnalyze::AnalyzeManyEvents(Int_t Nevents, Int_t module)   // analyzes many events   
{

  if ( fRootFile == 0 ) 
    cout << "AnalyzeManyEvents > " << "Root File not openned" << endl ;  
  else
    {
      //========== Get AliRun object from file 
      gAlice = (AliRun*) fRootFile->Get("gAlice") ;
      //=========== Get the PHOS object and associated geometry from the file      
      fPHOS  = (AliPHOSv0 *)gAlice->GetDetector("PHOS") ;
      fGeom  = AliPHOSGeometry::GetInstance( fPHOS->GetGeometry()->GetName(), fPHOS->GetGeometry()->GetTitle() );
      //========== Booking Histograms
      cout << "AnalyzeManyEvents > " << "Booking Histograms" << endl ; 
      BookingHistograms();
      Int_t ievent;
      Int_t relid[4] ; 
      AliPHOSDigit * digit ;
      AliPHOSEmcRecPoint * emc ;
      AliPHOSPpsdRecPoint * ppsd ;
      //      AliPHOSTrackSegment * tracksegment ;
      AliPHOSRecParticle * recparticle;
      for ( ievent=0; ievent<Nevents; ievent++)
	{  
          if (ievent==0)  cout << "AnalyzeManyEvents > " << "Starting Analyzing " << endl ; 
	  //========== Create the Clusterizer
	  fClu = new AliPHOSClusterizerv1() ; 
	  fClu->SetEmcEnergyThreshold(0.025) ; 
	  fClu->SetEmcClusteringThreshold(0.50) ; 
	  fClu->SetPpsdEnergyThreshold    (0.0000002) ; 
	  fClu->SetPpsdClusteringThreshold(0.0000001) ; 
	  fClu->SetLocalMaxCut(0.03) ;
	  fClu->SetCalibrationParameters(0., 0.00000001) ;  
	  //========== Creates the track segment maker
	  fTrs = new AliPHOSTrackSegmentMakerv1()  ;
	  fTrs->UnsetUnfoldFlag() ; 
	  //========== Creates the particle identifier
	  fPID = new AliPHOSPIDv1() ;
	  fPID->SetShowerProfileCuts(0.3, 1.8, 0.3, 1.8 ) ; 
	  fPID->Print() ; 	    
	  //========== Creates the Reconstructioner  
	  fRec = new AliPHOSReconstructioner(fClu, fTrs, fPID) ; 
	  //========== Event Number
	  if ( ( log10(ievent+1) - (Int_t)(log10(ievent+1)) ) == 0. ) cout <<  "AnalyzeManyEvents > " << "Event is " << ievent << endl ;  
	  //=========== Connects the various Tree's for evt
	  gAlice->GetEvent(ievent);
	  //=========== Gets the Digit TTree
	  gAlice->TreeD()->GetEvent(0) ;
	  //=========== Gets the number of entries in the Digits array 
	  TIter nextdigit(fPHOS->Digits()) ;
	  while( ( digit = (AliPHOSDigit *)nextdigit() ) )
	    {
	      fGeom->AbsToRelNumbering(digit->GetId(), relid) ;         
	      if (fClu->IsInEmc(digit)) fhEmcDigit->Fill(fClu->Calibrate(digit->GetAmp())) ; 
	      else    
		{  
		  if (relid[1]<17) fhVetoDigit->Fill(fClu->Calibrate(digit->GetAmp())); 
		  if (relid[1]>16) fhConvertorDigit->Fill(fClu->Calibrate(digit->GetAmp()));
		}
	    }
	  //=========== Do the reconstruction
	  fPHOS->Reconstruction(fRec);
	  //=========== Cluster in module
	  TIter nextEmc(fPHOS->EmcClusters()  ) ;
	  while((emc = (AliPHOSEmcRecPoint *)nextEmc())) 
	    {
	      if ( emc->GetPHOSMod() == module )
		{  
		  fhEmcCluster->Fill(  emc->GetTotalEnergy()  ); 
		  TIter nextPpsd( fPHOS->PpsdClusters()) ;
		  while((ppsd = (AliPHOSPpsdRecPoint *)nextPpsd())) 
		    {
		      if ( ppsd->GetPHOSMod() == module )
			{ 			  
			  if (!ppsd->GetUp()) fhConvertorEmc->Fill(emc->GetTotalEnergy(),ppsd->GetTotalEnergy()) ;
			}
		    } 
		}
	    }
	  //=========== Cluster in module PPSD Down
	  TIter nextPpsd(fPHOS->PpsdClusters() ) ;
	  while((ppsd = (AliPHOSPpsdRecPoint *)nextPpsd())) 
	    {
	      if ( ppsd->GetPHOSMod() == module )
		{ 
		  if (!ppsd->GetUp()) fhConvertorCluster->Fill(ppsd->GetTotalEnergy()) ;
		  if (ppsd->GetUp())  fhVetoCluster     ->Fill(ppsd->GetTotalEnergy()) ;
		}
	    }
	  //========== TRackSegments in the event
	  TIter nextRecParticle(fPHOS->RecParticles() ) ; 
	  while((recparticle = (AliPHOSRecParticle *)nextRecParticle())) 
	    {
	      if ( recparticle->GetPHOSTrackSegment()->GetPHOSMod() == module )
		{ 
		  cout << "Particle type is " << recparticle->GetType() << endl ;  
		  switch(recparticle->GetType())
		    {
		    case kGAMMA:
		      fhPhotonEnergy->Fill(recparticle->Energy() ) ; 
		      //fhPhotonPositionX->Fill(recpart. ) ;
		      //fhPhotonPositionY->Fill(recpart. ) ;                 
		      cout << "PHOTON" << endl;
		      break;
		    case kELECTRON:
		      fhElectronEnergy->Fill(recparticle->Energy() ) ; 
		      //fhElectronPositionX->Fill(recpart. ) ;
		      //fhElectronPositionY->Fill(recpart. ) ; 
		      cout << "ELECTRON" << endl;
		      break;
		    case kNEUTRALHADRON:
		      fhNeutralHadronEnergy->Fill(recparticle->Energy() ) ; 
		      //fhNeutralHadronPositionX->Fill(recpart. ) ;
		      //fhNeutralHadronPositionY->Fill(recpart. ) ; 
		      cout << "NEUTRAl HADRON" << endl;
		      break ;
		    case kNEUTRALEM:
		      fhNeutralEMEnergy->Fill(recparticle->Energy() ) ; 
		      //fhNeutralEMPositionX->Fill(recpart. ) ;
		      //fhNeutralEMPositionY->Fill(recpart. ) ; 
		      //cout << "NEUTRAL EM" << endl;
		      break ;
		    case kCHARGEDHADRON:
		      fhChargedHadronEnergy->Fill(recparticle->Energy() ) ; 
		      //fhChargedHadronPositionX->Fill(recpart. ) ;
		      //fhChargedHadronPositionY->Fill(recpart. ) ; 
		      cout << "CHARGED HADRON" << endl;
		      break ;
		    case kGAMMAHADRON:
		      fhPhotonHadronEnergy->Fill(recparticle->Energy() ) ; 
		      //fhPhotonHadronPositionX->Fill(recpart. ) ;
		      //fhPhotonHadronPositionY->Fill(recpart. ) ; 
		      cout << "PHOTON HADRON" << endl;
		      break ;		      
		    }
		}
	    }
	  // Deleting fClu, fTrs, fPID et fRec
	  fClu->Delete();
	  fTrs->Delete();
	  fPID->Delete();
	  fRec->Delete();

	}   // endfor
      SavingHistograms();
    }       // endif
}           // endfunction


//____________________________________________________________________________
void  AliPHOSAnalyze::BookingHistograms()
{
  if (fhEmcDigit )         
    delete fhEmcDigit  ;
  if (fhVetoDigit )     
    delete fhVetoDigit  ;
  if (fhConvertorDigit ) 
    delete fhConvertorDigit   ;
  if (fhEmcCluster   )  
    delete  fhEmcCluster   ;
  if (fhVetoCluster )     
    delete fhVetoCluster   ;
  if (fhConvertorCluster )
    delete fhConvertorCluster  ;
  if (fhConvertorEmc )    
    delete fhConvertorEmc  ;
 
  fhEmcDigit                = new TH1F("hEmcDigit",      "hEmcDigit",         1000,  0. ,  25.);
  fhVetoDigit               = new TH1F("hVetoDigit",     "hVetoDigit",         500,  0. ,  3.e-5);
  fhConvertorDigit          = new TH1F("hConvertorDigit","hConvertorDigit",    500,  0. ,  3.e-5);
  fhEmcCluster              = new TH1F("hEmcCluster",    "hEmcCluster",       1000,  0. ,  30.);
  fhVetoCluster             = new TH1F("hVetoCluster",   "hVetoCluster",       500,  0. ,  3.e-5);
  fhConvertorCluster        = new TH1F("hConvertorCluster","hConvertorCluster",500,  0. ,  3.e-5);
  fhConvertorEmc            = new TH2F("hConvertorEmc",  "hConvertorEmc",      200,  1. ,  3., 200, 0., 3.e-5);
  fhPhotonEnergy            = new TH1F("hPhotonEnergy",  "hPhotonEnergy",     1000,  0. ,  30.);
  fhElectronEnergy          = new TH1F("hElectronEnergy","hElectronEnergy",   1000,  0. ,  30.);
  fhNeutralHadronEnergy     = new TH1F("hNeutralHadronEnergy", "hNeutralHadronEnergy",    1000,  0. ,  30.);
  fhNeutralEMEnergy         = new TH1F("hNeutralEMEnergy", "hNeutralEMEnergy",    1000,  0. ,  30.);
  fhChargedHadronEnergy     = new TH1F("hChargedHadronEnergy", "hChargedHadronEnergy",    1000,  0. ,  30.);
  fhPhotonHadronEnergy      = new TH1F("hPhotonHadronEnergy","hPhotonHadronEnergy",500,-80. , 80.);
  fhPhotonPositionX         = new TH1F("hPhotonPositionX","hPhotonPositionX",   500,-80. , 80.);
  fhElectronPositionX       = new TH1F("hElectronPositionX","hElectronPositionX",500,-80. , 80.);
  fhNeutralHadronPositionX  = new TH1F("hNeutralHadronPositionX","hNeutralHadronPositionX",500,-80. , 80.);
  fhNeutralEMPositionX      = new TH1F("hNeutralEMPositionX","hNeutralEMPositionX",500,-80. , 80.);
  fhChargedHadronPositionX  = new TH1F("hChargedHadronPositionX","hChargedHadronPositionX",500,-80. , 80.);
  fhPhotonHadronPositionX   = new TH1F("hPhotonHadronPositionX","hPhotonHadronPositionX",500,-80. , 80.);
  fhPhotonPositionY         = new TH1F("hPhotonPositionY","hPhotonPositionY",   500,-80. , 80.);
  fhElectronPositionY       = new TH1F("hElectronPositionY","hElectronPositionY",500,-80. , 80.);
  fhNeutralHadronPositionY  = new TH1F("hNeutralHadronPositionY","hNeutralHadronPositionY",500,-80. , 80.);
  fhNeutralEMPositionY      = new TH1F("hNeutralEMPositionY","hNeutralEMPositionY",500,-80. , 80.);
  fhChargedHadronPositionY  = new TH1F("hChargedHadronPositionY","hChargedHadronPositionY",500,-80. , 80.);
  fhPhotonHadronPositionY   = new TH1F("hPhotonHadronPositionY","hPhotonHadronPositionY",500,-80. , 80.);


}
//____________________________________________________________________________
Bool_t AliPHOSAnalyze::Init(Int_t evt)
{

  Bool_t ok = kTRUE ; 
  
   //========== Open galice root file  

  if ( fRootFile == 0 ) {
    Text_t * name  = new Text_t[80] ; 
    cout << "AnalyzeOneEvent > Enter file root file name : " ;  
    cin >> name ; 
    Bool_t ok = OpenRootFile(name) ; 
    if ( !ok )
      cout << " AliPHOSAnalyze > Error opening " << name << endl ; 
    else { 
      //========== Get AliRun object from file 
      
      gAlice = (AliRun*) fRootFile->Get("gAlice") ;
      
      //=========== Get the PHOS object and associated geometry from the file 
      
      fPHOS  = (AliPHOSv0 *)gAlice->GetDetector("PHOS") ;
      fGeom  = AliPHOSGeometry::GetInstance( fPHOS->GetGeometry()->GetName(), fPHOS->GetGeometry()->GetTitle() );
    } // else !ok
  } // if fRootFile
  
  if ( ok ) {
    
    //========== Create the Clusterizer

    fClu =  new AliPHOSClusterizerv1() ; 
    fClu->SetEmcEnergyThreshold(0.030) ; 
    fClu->SetEmcClusteringThreshold(0.50) ; 
    fClu->SetPpsdEnergyThreshold    (0.0000002) ; 
    fClu->SetPpsdClusteringThreshold(0.0000001) ; 
    fClu->SetLocalMaxCut(0.03) ;
    fClu->SetCalibrationParameters(0., 0.00000001) ;  
    cout <<  "AnalyzeOneEvent > using clusterizer " << fClu->GetName() << endl ; 
    fClu->PrintParameters() ; 
    
    //========== Creates the track segment maker
    
    fTrs = new AliPHOSTrackSegmentMakerv1() ;
    cout <<  "AnalyzeOneEvent > using tack segment maker " << fTrs->GetName() << endl ; 
    fTrs->UnsetUnfoldFlag() ;
    
    //========== Creates the particle identifier
    
    fPID = new AliPHOSPIDv1() ;
    cout <<  "AnalyzeOneEvent > using particle identifier " << fPID->GetName() << endl ; 
    
    //========== Creates the Reconstructioner  
    
    fRec = new AliPHOSReconstructioner(fClu, fTrs, fPID) ;     
    
    //=========== Connect the various Tree's for evt
    
    if ( evt == -999 ) {
      cout <<  "AnalyzeOneEvent > Enter event number : " ; 
      cin >> evt ;  
      cout << evt << endl ; 
    }
    fEvt = evt ; 
    gAlice->GetEvent(evt);
    
    //=========== Get the Digit TTree
    
    gAlice->TreeD()->GetEvent(0) ;     
    
  } // ok
  
  return ok ; 
}


//____________________________________________________________________________
void AliPHOSAnalyze::DisplayKineEvent(Int_t evt)
{
  if (evt == -999) 
    evt = fEvt ;

  Int_t module ; 
  cout <<  "DisplayKineEvent > which module (1-5,  -1: all) ? " ; 
  cin >> module ; cout << module << endl ; 

  Int_t testparticle ; 
  cout << " 22      : PHOTON " << endl 
       << " (-)11   : (POSITRON)ELECTRON " << endl 
       << " (-)2112 : (ANTI)NEUTRON " << endl  
       << " -999    : Everything else " << endl ; 
  cout  <<  "DisplayKineEvent > enter PDG particle code to display " ; 
  cin >> testparticle ; cout << testparticle << endl ; 

  Text_t histoname[80] ;
  sprintf(histoname,"Event %d: Incident particles in module %d", evt, module) ; 

  Double_t tm, tM, pm, pM ; // min and Max theta and phi covered by module   
  fGeom->EmcModuleCoverage(module, tm, tM, pm, pM, kDegre) ;

  Double_t theta, phi ; 
  fGeom->EmcXtalCoverage(theta, phi, kDegre) ;

  Int_t tdim = (Int_t)( (tM - tm) / theta ) ; 
  Int_t pdim = (Int_t)( (pM - pm) / phi ) ; 

  tm -= theta ; 
  tM += theta ; 
  pm -= phi ; 
  pM += phi ; 

  TH2F * histoparticle = new TH2F("histoparticle",  histoname, 
				 	  pdim, pm, pM, tdim, tm, tM) ; 
  histoparticle->SetStats(kFALSE) ;

  // Get pointers to Alice Particle TClonesArray

  TParticle * particle;
  TClonesArray * particlearray  = gAlice->Particles();    

  Text_t canvasname[80];
  sprintf(canvasname,"Particles incident in PHOS/EMC module # %d",module) ;
  TCanvas * kinecanvas = new TCanvas("kinecanvas", canvasname, 650, 500) ; 

  // get the KINE Tree

  TTree * kine =  gAlice->TreeK() ; 
  Stat_t nParticles =  kine->GetEntries() ; 
  cout << "DisplayKineEvent > events in kine " << nParticles << endl ; 
  
  // loop over particles

  Double_t kRADDEG = 180. / TMath::Pi() ; 
  Int_t index1 ; 
  Int_t nparticlein = 0 ; 
  for (index1 = 0 ; index1 < nParticles ; index1++){
    Int_t nparticle = particlearray->GetEntriesFast() ;
    Int_t index2 ; 
    for( index2 = 0 ; index2 < nparticle ; index2++) {         
      particle            = (TParticle*)particlearray->UncheckedAt(index2) ;
      Int_t  particletype = particle->GetPdgCode() ;
      if (testparticle == -999 || testparticle == particletype) { 
	Double_t phi        = particle->Phi() ;
	Double_t theta      = particle->Theta() ;
	Int_t mod ; 
	Double_t x, z ; 
	fGeom->ImpactOnEmc(theta, phi, mod, z, x) ;
	if ( mod == module ) {
	  nparticlein++ ; 
	  histoparticle->Fill(phi*kRADDEG, theta*kRADDEG, particle->Energy() ) ; 
	} 
      } 
    }
  }
  kinecanvas->Draw() ; 
  histoparticle->Draw("color") ; 
  TPaveText *  pavetext = new TPaveText(294, 100, 300, 101); 
  Text_t text[40] ; 
  sprintf(text, "Particles: %d ", nparticlein) ;
  pavetext->AddText(text) ; 
  pavetext->Draw() ; 
  kinecanvas->Update(); 

}
//____________________________________________________________________________
void AliPHOSAnalyze::DisplayRecParticles()
{
  if (fEvt == -999) {
    cout << "DisplayRecParticles > Analyze an event first ... (y/n) " ; 
    Text_t answer[1] ; 
    cin >> answer ; cout << answer ; 
    if ( answer == "y" ) 
      AnalyzeOneEvent() ;
  } 
    if (fEvt != -999) {
      
      Int_t module ; 
      cout <<  "DisplayRecParticles > which module (1-5,  -1: all) ? " ; 
      cin >> module ; cout << module << endl ;
      Text_t histoname[80] ; 
      sprintf(histoname,"Event %d: Reconstructed particles in module %d", fEvt, module) ; 
      Double_t tm, tM, pm, pM ; // min and Max theta and phi covered by module   
      fGeom->EmcModuleCoverage(module, tm, tM, pm, pM, kDegre) ;
      Double_t theta, phi ; 
      fGeom->EmcXtalCoverage(theta, phi, kDegre) ;
      Int_t tdim = (Int_t)( (tM - tm) / theta ) ; 
      Int_t pdim = (Int_t)( (pM - pm) / phi ) ; 
      tm -= theta ; 
      tM += theta ; 
      pm -= phi ; 
      TH2F * histoRparticle = new TH2F("histoRparticle",  histoname, 
				       pdim, pm, pM, tdim, tm, tM) ; 
      histoRparticle->SetStats(kFALSE) ;
      Text_t canvasname[80] ; 
      sprintf(canvasname, "Reconstructed particles in PHOSmodule # %d", module) ;
      TCanvas * rparticlecanvas = new TCanvas("RparticleCanvas", canvasname, 650, 500) ; 
      RecParticlesList * rpl = fPHOS->RecParticles() ; 
      Int_t nRecParticles = rpl->GetEntries() ; 
      Int_t nRecParticlesInModule = 0 ; 
      TIter nextRecPart(rpl) ; 
      AliPHOSRecParticle * rp ; 
      cout << "DisplayRecParticles > " << nRecParticles << " reconstructed particles " << endl ; 
      Double_t kRADDEG = 180. / TMath::Pi() ; 
      while ( (rp = (AliPHOSRecParticle *)nextRecPart() ) ) {
	AliPHOSTrackSegment * ts = rp->GetPHOSTrackSegment() ; 
	if ( ts->GetPHOSMod() == module ) {  
	  nRecParticlesInModule++ ; 
	  Double_t theta = rp->Theta() * kRADDEG ;
	  Double_t phi   = rp->Phi() * kRADDEG ;
	  Double_t energy = rp->Energy() ; 
	  histoRparticle->Fill(phi, theta, energy) ;
	}
      }
      histoRparticle->Draw("color") ; 

      nextRecPart.Reset() ; 
      while ( (rp = (AliPHOSRecParticle *)nextRecPart() ) ) {
	AliPHOSTrackSegment * ts = rp->GetPHOSTrackSegment() ; 
	if ( ts->GetPHOSMod() == module )  
	  rp->Draw("P") ; 
      }

      Text_t text[80] ; 
      sprintf(text, "reconstructed particles: %d", nRecParticlesInModule) ;
      TPaveText *  pavetext = new TPaveText(292, 100, 300, 101); 
      pavetext->AddText(text) ; 
      pavetext->Draw() ; 
      rparticlecanvas->Update() ; 
    }
}

//____________________________________________________________________________
void AliPHOSAnalyze::DisplayRecPoints()
{
  if (fEvt == -999) {
    cout << "DisplayRecPoints > Analyze an event first ... (y/n) " ; 
    Text_t answer[1] ; 
    cin >> answer ; cout << answer ; 
    if ( answer == "y" ) 
      AnalyzeOneEvent() ;
  } 
    if (fEvt != -999) {
      
      Int_t module ; 
      cout <<  "DisplayRecPoints > which module (1-5,  -1: all) ? " ; 
      cin >> module ; cout << module << endl ; 

      Text_t canvasname[80];
      sprintf(canvasname,"Digits in PHOS/EMC module # %d",module) ;
      TCanvas * modulecanvas = new TCanvas("module", canvasname, 650, 500) ; 
      modulecanvas->Draw() ;

      //=========== Creating 2d-histogram of the PHOS module
      // a little bit junkie but is used to test Geom functinalities

      Double_t tm, tM, pm, pM ; // min and Max theta and phi covered by module   
      
      fGeom->EmcModuleCoverage(module, tm, tM, pm, pM); 
      // convert angles into coordinates local to the EMC module of interest

      Int_t emcModuleNumber ;
      Double_t emcModulexm, emcModulezm ; // minimum local coordinate in a given EMCA module
      Double_t emcModulexM, emcModulezM ; // maximum local coordinate in a given EMCA module
      fGeom->ImpactOnEmc(tm, pm, emcModuleNumber, emcModulezm, emcModulexm) ;
      fGeom->ImpactOnEmc(tM, pM, emcModuleNumber, emcModulezM, emcModulexM) ;
      Int_t xdim = (Int_t)( ( emcModulexM - emcModulexm ) / fGeom->GetCrystalSize(0) ) ;  
      Int_t zdim = (Int_t)( ( emcModulezM - emcModulezm ) / fGeom->GetCrystalSize(2) ) ;
      Float_t xmin = emcModulexm - fGeom->GetCrystalSize(0) ; 
      Float_t xMax = emcModulexM + fGeom->GetCrystalSize(0) ; 
      Float_t zmin = emcModulezm - fGeom->GetCrystalSize(2) ; 
      Float_t zMax = emcModulezM + fGeom->GetCrystalSize(2) ;     
      Text_t histoname[80];
      sprintf(histoname,"Event %d: Digits and RecPoints in module %d", fEvt, module) ;
      TH2F * hModule = new TH2F("HistoReconstructed", histoname,
				xdim, xmin, xMax, zdim, zmin, zMax) ;  
      hModule->SetMaximum(2.0);
      hModule->SetMinimum(0.0);
      hModule->SetStats(kFALSE); 

      TIter next(fPHOS->Digits()) ;
      Float_t energy, y, z;
      Float_t etot=0.;
      Int_t relid[4]; Int_t nDigits = 0 ;
      AliPHOSDigit * digit ; 

      // Making 2D histogram of the EMC module
      while((digit = (AliPHOSDigit *)next())) 
	{  
	  fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
	  if (relid[0] == module && relid[1] == 0)  
	    {  
	      energy = fClu->Calibrate(digit->GetAmp()) ;
              cout << "Energy is " << energy << " and threshold is " << fClu->GetEmcEnergyThreshold() << endl; 
	      if (energy >  fClu->GetEmcEnergyThreshold()  ){
		nDigits++ ;
		etot += energy ; 
		fGeom->RelPosInModule(relid,y,z) ;   
		hModule->Fill(y, z, energy) ;
	      }
	    } 
	}
      cout <<"DrawRecPoints >  Found in module " 
	   << module << " " << nDigits << "  digits with total energy " << etot << endl ;
      hModule->Draw("col2") ;

      //=========== Cluster in module

      TClonesArray * emcRP = fPHOS->EmcClusters() ;
      etot = 0.; 
      Int_t totalnClusters = 0 ; 
      Int_t nClusters = 0 ;
      TIter nextemc(emcRP) ;
      AliPHOSEmcRecPoint * emc ;
      while((emc = (AliPHOSEmcRecPoint *)nextemc())) 
	{
	  //	  Int_t numberofprimaries ;
	  //	  Int_t * primariesarray = new Int_t[10] ;
	  //	  emc->GetPrimaries(numberofprimaries, primariesarray) ;
	  totalnClusters++ ;
	  if ( emc->GetPHOSMod() == module )
	    { 
	      nClusters++ ; 
	      energy = emc->GetTotalEnergy() ;   
	      etot+= energy ;  
	      emc->Draw("M") ;
	    }
	}
      cout << "DrawRecPoints > Found " << totalnClusters << " EMC Clusters in PHOS" << endl ; 
      cout << "DrawRecPoints > Found in module " << module << "  " << nClusters << " EMC Clusters " << endl ;
      cout << "DrawRecPoints > total energy  " << etot << endl ; 

      TPaveText *  pavetext = new TPaveText(22, 80, 83, 90); 
      Text_t text[40] ; 
      sprintf(text, "digits: %d;  clusters: %d", nDigits, nClusters) ;
      pavetext->AddText(text) ; 
      pavetext->Draw() ; 
      modulecanvas->Update(); 
 
      //=========== Cluster in module PPSD Down

      TClonesArray * ppsdRP = fPHOS->PpsdClusters() ;
      etot = 0.; 
      TIter nextPpsd(ppsdRP) ;
      AliPHOSPpsdRecPoint * ppsd ;
      while((ppsd = (AliPHOSPpsdRecPoint *)nextPpsd())) 
	{
	  totalnClusters++ ;
	  if ( ppsd->GetPHOSMod() == module )
	    { 
	      nClusters++ ; 
	      energy = ppsd->GetEnergy() ;   
	      etot+=energy ;  
	      if (!ppsd->GetUp()) ppsd->Draw("P") ;
	    }
	}
      cout << "DrawRecPoints > Found " << totalnClusters << " Ppsd Down Clusters in PHOS" << endl ; 
      cout << "DrawRecPoints > Found in module " << module << "  " << nClusters << " Ppsd Down Clusters " << endl ;
      cout << "DrawRecPoints > total energy  " << etot << endl ; 

      //=========== Cluster in module PPSD Up
  
      ppsdRP = fPHOS->PpsdClusters() ;
      etot = 0.; 
      TIter nextPpsdUp(ppsdRP) ;
      while((ppsd = (AliPHOSPpsdRecPoint *)nextPpsdUp())) 
	{
	  totalnClusters++ ;
	  if ( ppsd->GetPHOSMod() == module )
	    { 
	      nClusters++ ; 
	      energy = ppsd->GetEnergy() ;   
	      etot+=energy ;  
	      if (ppsd->GetUp()) ppsd->Draw("P") ;
	    }
	}
  cout << "DrawRecPoints > Found " << totalnClusters << " Ppsd Up Clusters in PHOS" << endl ; 
  cout << "DrawRecPoints > Found in module " << module << "  " << nClusters << " Ppsd Up Clusters " << endl ;
  cout << "DrawRecPoints > total energy  " << etot << endl ; 
    
    } // if !-999
}

//____________________________________________________________________________
void AliPHOSAnalyze::DisplayTrackSegments()
{
  if (fEvt == -999) {
    cout << "DisplayTrackSegments > Analyze an event first ... (y/n) " ; 
    Text_t answer[1] ; 
    cin >> answer ; cout << answer ; 
    if ( answer == "y" ) 
      AnalyzeOneEvent() ;
  } 
    if (fEvt != -999) {

      Int_t module ; 
      cout <<  "DisplayTrackSegments > which module (1-5,  -1: all) ? " ; 
      cin >> module ; cout << module << endl ; 
      //=========== Creating 2d-histogram of the PHOS module
      // a little bit junkie but is used to test Geom functinalities
      
      Double_t tm, tM, pm, pM ; // min and Max theta and phi covered by module   
      
      fGeom->EmcModuleCoverage(module, tm, tM, pm, pM); 
      // convert angles into coordinates local to the EMC module of interest
      
      Int_t emcModuleNumber ;
      Double_t emcModulexm, emcModulezm ; // minimum local coordinate in a given EMCA module
      Double_t emcModulexM, emcModulezM ; // maximum local coordinate in a given EMCA module
      fGeom->ImpactOnEmc(tm, pm, emcModuleNumber, emcModulezm, emcModulexm) ;
      fGeom->ImpactOnEmc(tM, pM, emcModuleNumber, emcModulezM, emcModulexM) ;
      Int_t xdim = (Int_t)( ( emcModulexM - emcModulexm ) / fGeom->GetCrystalSize(0) ) ;  
      Int_t zdim = (Int_t)( ( emcModulezM - emcModulezm ) / fGeom->GetCrystalSize(2) ) ;
      Float_t xmin = emcModulexm - fGeom->GetCrystalSize(0) ; 
      Float_t xMax = emcModulexM + fGeom->GetCrystalSize(0) ; 
      Float_t zmin = emcModulezm - fGeom->GetCrystalSize(2) ; 
      Float_t zMax = emcModulezM + fGeom->GetCrystalSize(2) ;     
      Text_t histoname[80];
      sprintf(histoname,"Event %d: Track Segments in module %d", fEvt, module) ; 
      TH2F * histotrack = new TH2F("histotrack",  histoname, 
				   xdim, xmin, xMax, zdim, zmin, zMax) ;  
      histotrack->SetStats(kFALSE); 
      Text_t canvasname[80];
      sprintf(canvasname,"Track segments in PHOS/EMC-PPSD module # %d", module) ;
      TCanvas * trackcanvas = new TCanvas("TrackSegmentCanvas", canvasname, 650, 500) ; 
      histotrack->Draw() ; 

      TrackSegmentsList * trsegl = fPHOS->TrackSegments() ;
      AliPHOSTrackSegment * trseg ;
 
      Int_t nTrackSegments = trsegl->GetEntries() ;
      Int_t index ;
      Float_t etot = 0 ;
      Int_t nTrackSegmentsInModule = 0 ; 
      for(index = 0; index < nTrackSegments ; index++){
	trseg = (AliPHOSTrackSegment * )trsegl->At(index) ;
	etot+= trseg->GetEnergy() ;
	if ( trseg->GetPHOSMod() == module ) { 
	  nTrackSegmentsInModule++ ; 
	  trseg->Draw("P");
	}
      } 
      Text_t text[80] ; 
      sprintf(text, "track segments: %d", nTrackSegmentsInModule) ;
      TPaveText *  pavetext = new TPaveText(22, 80, 83, 90); 
      pavetext->AddText(text) ; 
      pavetext->Draw() ; 
      trackcanvas->Update() ; 
      cout << "DisplayTrackSegments > Found " << trsegl->GetEntries() << " Track segments with total energy "<< etot << endl ;
    
   }
}
//____________________________________________________________________________
Bool_t AliPHOSAnalyze::OpenRootFile(Text_t * name)
{
  fRootFile   = new TFile(name) ;
  return  fRootFile->IsOpen() ; 
}
//____________________________________________________________________________
void AliPHOSAnalyze::SavingHistograms()
{
  Text_t outputname[80] ;// = fRootFile->GetName();
  sprintf(outputname,"%s.analyzed",fRootFile->GetName());
  TFile output(outputname,"RECREATE");
  output.cd();
  if (fhEmcDigit )         
    fhEmcDigit->Write()  ;
  if (fhVetoDigit )  
    fhVetoDigit->Write()  ;
  if (fhConvertorDigit ) 
    fhConvertorDigit->Write()   ;
  if (fhEmcCluster   )
    fhEmcCluster->Write()   ;
  if (fhVetoCluster ) 
    fhVetoCluster->Write()   ;
  if (fhConvertorCluster )
    fhConvertorCluster->Write()  ;
  if (fhConvertorEmc ) 
    fhConvertorEmc->Write()  ;
  if (fhPhotonEnergy)    
    fhPhotonEnergy->Write() ;
  if (fhPhotonPositionX)  
    fhPhotonPositionX->Write() ;
  if (fhPhotonPositionY)  
    fhPhotonPositionX->Write() ;
  if (fhElectronEnergy)  
    fhElectronEnergy->Write() ;
  if (fhElectronPositionX)
    fhElectronPositionX->Write() ;
  if (fhElectronPositionY) 
    fhElectronPositionX->Write() ;
  if (fhNeutralHadronEnergy) 
    fhNeutralHadronEnergy->Write() ;
  if (fhNeutralHadronPositionX)
    fhNeutralHadronPositionX->Write() ;
  if (fhNeutralHadronPositionY) 
    fhNeutralHadronPositionX->Write() ;
  if (fhNeutralEMEnergy)   
    fhNeutralEMEnergy->Write() ;
  if (fhNeutralEMPositionX)
    fhNeutralEMPositionX->Write() ;
  if (fhNeutralEMPositionY) 
    fhNeutralEMPositionX->Write() ;
  if (fhChargedHadronEnergy) 
    fhChargedHadronEnergy->Write() ;
  if (fhChargedHadronPositionX) 
    fhChargedHadronPositionX->Write() ;
  if (fhChargedHadronPositionY)
    fhChargedHadronPositionX->Write() ;
  if (fhPhotonHadronEnergy) 
    fhPhotonHadronEnergy->Write() ;
  if (fhPhotonHadronPositionX) 
    fhPhotonHadronPositionX->Write() ;
  if (fhPhotonHadronPositionY)
    fhPhotonHadronPositionX->Write() ;

  output.Write();
  output.Close();
}
