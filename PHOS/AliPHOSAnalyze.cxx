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

#include "TH2.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h" 

// --- Standard library ---

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliPHOSAnalyze.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSParticleGuesserv1.h"
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
  
  Bool_t OK = OpenRootFile(name)  ; 
  if ( !OK ) {
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

  delete fPag ; 
  fPag = 0 ;

  delete fRec ; 
  fRec = 0 ;

  delete fTrs ; 
  fTrs = 0 ;

}

//____________________________________________________________________________
void AliPHOSAnalyze::AnalyzeOneEvent(Int_t evt)
{
  Bool_t OK = Init(evt) ; 
  
  if ( OK ) {
    //=========== Get the number of entries in the Digits array
    
    Int_t nId = fPHOS->Digits()->GetEntries();    
    printf("AnalyzeOneEvent > Number of entries in the Digit array is %d \n",nId);
    
    //=========== Do the reconstruction
    
    cout << "AnalyzeOneEvent > Found  " << nId << "  digits in PHOS"   << endl ;  
    
    fPHOS->Reconstruction(fRec);  
    
    // =========== End of reconstruction
    
    cout << "AnalyzeOneEvent > event # " << fEvt << " processed" << endl ;   
  } // OK
  else
    cout << "AnalyzeOneEvent > filed to process event # " << evt << endl ;   

}

//____________________________________________________________________________
Bool_t AliPHOSAnalyze::Init(Int_t evt)
{

  Bool_t OK = kTRUE ; 
  
   //========== Open galice root file  

  if ( fRootFile == 0 ) {
    Text_t * name  = new Text_t[80] ; 
    cout << "AnalyzeOneEvent > Enter file root file name : " ;  
    cin >> name ; 
    Bool_t OK = OpenRootFile(name) ; 
    if ( !OK )
      cout << " AliPHOSAnalyze > Error opening " << name << endl ; 
    else { 
      //========== Get AliRun object from file 
      
      gAlice = (AliRun*) fRootFile->Get("gAlice") ;
      
      //=========== Get the PHOS object and associated geometry from the file 
      
      fPHOS  = (AliPHOSv0 *)gAlice->GetDetector("PHOS") ;
      fGeom  = AliPHOSGeometry::GetInstance( fPHOS->GetGeometry()->GetName(), fPHOS->GetGeometry()->GetTitle() );
    } // else !OK
  } // if fRootFile
  
  if ( OK ) {
    
    //========== Create the Clusterizer

    fClu = new AliPHOSClusterizerv1() ; 
    fClu->SetEmcEnergyThreshold(0.01) ; 
    fClu->SetEmcClusteringThreshold(0.1) ; 
    fClu->SetPpsdEnergyThreshold(0.0000001) ; 
    fClu->SetPpsdClusteringThreshold(0.0000002) ; 
    fClu->SetLocalMaxCut(0.03) ;
    fClu->SetCalibrationParameters(0., 0.0000001) ;  
    cout <<  "AnalyzeOneEvent > using clusterizer " << fClu->GetName() << endl ; 
    fClu->PrintParameters() ; 
    
    //========== Creates the track segment maker
    
    fTrs = new AliPHOSTrackSegmentMakerv1() ;
    cout <<  "AnalyzeOneEvent > using tack segment maker " << fTrs->GetName() << endl ; 
    
    //========== Creates the particle guesser
    
    fPag = new AliPHOSParticleGuesserv1() ;
    cout <<  "AnalyzeOneEvent > using particle guess " << fPag->GetName() << endl ; 
    
    //========== Creates the Reconstructioner  
    
    fRec = new AliPHOSReconstructioner(fClu, fTrs, fPag) ;     
    
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
    
  } // OK
  
  return OK ; 
}

//____________________________________________________________________________
void AliPHOSAnalyze:: DisplayKineEvent(Int_t evt)
{
  if (evt == -999) 
    evt = fEvt ;

  Int_t Module ; 
  cout <<  "DisplayKineEvent > which module (1-5,  -1: all) ? " ; 
  cin >> Module ; cout << Module << endl ; 

  Int_t TestParticle ; 
  cout << " 22      : PHOTON " << endl 
       << " (-)11   : (POSITRON)ELECTRON " << endl 
       << " (-)2112 : (ANTI)NEUTRON " << endl  
       << " -999    : Everything else " << endl ; 
  cout  <<  "DisplayKineEvent > enter PDG particle code to display " ; 
  cin >> TestParticle ; cout << TestParticle << endl ; 

  Text_t HistoName[80] ;
  sprintf(HistoName,"Event %d: Incident particles in module %d", evt, Module) ; 

  Double_t tm, tM, pm, pM ; // min and Max theta and phi covered by Module   
  fGeom->EmcModuleCoverage(Module, tm, tM, pm, pM, kDegre) ;

  Double_t theta, phi ; 
  fGeom->EmcXtalCoverage(theta, phi, kDegre) ;

  Int_t tdim = (Int_t)( (tM - tm) / theta ) ; 
  Int_t pdim = (Int_t)( (pM - pm) / phi ) ; 

  tm -= theta ; 
  tM += theta ; 
  pm -= phi ; 
  pM += phi ; 

  TH2F * HistoParticle = new TH2F("HistoParticle",  HistoName, 
				 	  pdim, pm, pM, tdim, tm, tM) ; 
  HistoParticle->SetStats(kFALSE) ;

  // Get pointers to Alice Particle TClonesArray

  TParticle * Particle;
  TClonesArray * ArrayOfParticles  = gAlice->Particles();    

  Text_t CanvasName[80];
  sprintf(CanvasName,"Particles incident in PHOS/EMC module # %d",Module) ;
  TCanvas * KineCanvas = new TCanvas("KineCanvas", CanvasName, 650, 500) ; 

  // get the KINE Tree

  TTree * Kine =  gAlice->TreeK() ; 
  Stat_t NumberOfParticles =  Kine->GetEntries() ; 
  cout << "DisplayKineEvent > events in Kine " << NumberOfParticles << endl ; 
  
  // loop over particles

  Double_t raddeg = 180. / TMath::Pi() ; 
  Int_t index1 ; 
  Int_t nparticlein = 0 ; 
  for (index1 = 0 ; index1 < NumberOfParticles ; index1++){
    Int_t nparticle = ArrayOfParticles->GetEntriesFast() ;
    Int_t index2 ; 
    for( index2 = 0 ; index2 < nparticle ; index2++) {         
      Particle            = (TParticle*)ArrayOfParticles->UncheckedAt(index2) ;
      Int_t  ParticleType = Particle->GetPdgCode() ;
      if (TestParticle == -999 || TestParticle == ParticleType) { 
	Double_t Phi        = Particle->Phi() ;
	Double_t Theta      = Particle->Theta() ;
	Int_t mod ; 
	Double_t x, z ; 
	fGeom->ImpactOnEmc(Theta, Phi, mod, z, x) ;
	if ( mod == Module ) {
	  nparticlein++ ; 
	  HistoParticle->Fill(Phi*raddeg, Theta*raddeg, Particle->Energy() ) ; 
	} 
      } 
    }
  }
  KineCanvas->Draw() ; 
  HistoParticle->Draw("color") ; 
  TPaveText *  PaveText = new TPaveText(294, 100, 300, 101); 
  Text_t text[40] ; 
  sprintf(text, "Particles: %d ", nparticlein) ;
  PaveText->AddText(text) ; 
  PaveText->Draw() ; 
  KineCanvas->Update(); 

}
//____________________________________________________________________________
void AliPHOSAnalyze::DisplayRecParticles()
{
  if (fEvt == -999) {
    cout << "DisplayRecPoints > Analyze an event first ... (y/n) " ; 
    Text_t answer[1] ; 
    cin >> answer ; cout << answer ; 
    if ( answer == "y" ) 
      AnalyzeOneEvent() ;
  } 
    if (fEvt != -999) {
      
      Int_t Module ; 
      cout <<  "DisplayRecPoints > which module (1-5,  -1: all) ? " ; 
      cin >> Module ; cout << Module << endl ;
      Text_t HistoName[80] ; 
      sprintf(HistoName,"Event %d: Reconstructed particles in module %d", fEvt, Module) ; 
      Double_t tm, tM, pm, pM ; // min and Max theta and phi covered by Module   
      fGeom->EmcModuleCoverage(Module, tm, tM, pm, pM, kDegre) ;
      Double_t theta, phi ; 
      fGeom->EmcXtalCoverage(theta, phi, kDegre) ;
      Int_t tdim = (Int_t)( (tM - tm) / theta ) ; 
      Int_t pdim = (Int_t)( (pM - pm) / phi ) ; 
      tm -= theta ; 
      tM += theta ; 
      pm -= phi ; 
      TH2F * HistoRParticle = new TH2F("HistoRParticle",  HistoName, 
				       pdim, pm, pM, tdim, tm, tM) ; 
      HistoRParticle->SetStats(kFALSE) ;
      Text_t CanvasName[80] ; 
      sprintf(CanvasName, "Reconstructed particles in PHOSmodule # %d", Module) ;
      TCanvas * RParticleCanvas = new TCanvas("RparticleCanvas", CanvasName, 650, 500) ; 
      RecParticlesList * rpl = fPHOS->RecParticles() ; 
      Int_t NRecParticles = rpl->GetEntries() ; 
      Int_t NRecParticlesInModule = 0 ; 
      TIter nextRecPart(rpl) ; 
      AliPHOSRecParticle * rp ; 
      cout << "DisplayRecParticles > " << NRecParticles << " reconstructed particles " << endl ; 
      Double_t raddeg = 180. / TMath::Pi() ; 
      while ( (rp = (AliPHOSRecParticle *)nextRecPart() ) ) {
	AliPHOSTrackSegment * ts = rp->GetPHOSTrackSegment() ; 
	if ( ts->GetPHOSMod() == Module ) {  
	  NRecParticlesInModule++ ; 
	  Double_t theta = rp->Theta() * raddeg ;
	  Double_t phi   = rp->Phi() * raddeg ;
	  Double_t energy = rp->Energy() ; 
	  HistoRParticle->Fill(phi, theta, energy) ;
	}
      }
      HistoRParticle->Draw("color") ; 
      Text_t text[80] ; 
      sprintf(text, "reconstructed particles: %d", NRecParticlesInModule) ;
      TPaveText *  PaveText = new TPaveText(292, 100, 300, 101); 
      PaveText->AddText(text) ; 
      PaveText->Draw() ; 
      RParticleCanvas->Update() ; 
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
      
      Int_t Module ; 
      cout <<  "DisplayRecPoints > which module (1-5,  -1: all) ? " ; 
      cin >> Module ; cout << Module << endl ; 

      Text_t CanvasName[80];
      sprintf(CanvasName,"Digits in PHOS/EMC module # %d",Module) ;
      TCanvas * ModuleCanvas = new TCanvas("Module", CanvasName, 650, 500) ; 
      ModuleCanvas->Draw() ;

      //=========== Creating 2d-histogram of the PHOS Module
      // a little bit junkie but is used to test Geom functinalities

      Double_t tm, tM, pm, pM ; // min and Max theta and phi covered by Module   
      
      fGeom->EmcModuleCoverage(Module, tm, tM, pm, pM); 
      // convert angles into coordinates local to the EMC module of interest

      Int_t EmcModuleNumber ;
      Double_t EmcModulexm, EmcModulezm ; // minimum local coordinate in a given EMCA module
      Double_t EmcModulexM, EmcModulezM ; // maximum local coordinate in a given EMCA module
      fGeom->ImpactOnEmc(tm, pm, EmcModuleNumber, EmcModulezm, EmcModulexm) ;
      fGeom->ImpactOnEmc(tM, pM, EmcModuleNumber, EmcModulezM, EmcModulexM) ;
      Int_t xdim = (Int_t)( ( EmcModulexM - EmcModulexm ) / fGeom->GetCrystalSize(0) ) ;  
      Int_t zdim = (Int_t)( ( EmcModulezM - EmcModulezm ) / fGeom->GetCrystalSize(2) ) ;
      Float_t xmin = EmcModulexm - fGeom->GetCrystalSize(0) ; 
      Float_t xMax = EmcModulexM + fGeom->GetCrystalSize(0) ; 
      Float_t zmin = EmcModulezm - fGeom->GetCrystalSize(2) ; 
      Float_t zMax = EmcModulezM + fGeom->GetCrystalSize(2) ;     
      Text_t HistoName[80];
      sprintf(HistoName,"Event %d: Digits and RecPoints in module %d", fEvt, Module) ;
      TH2F * hModule = new TH2F("HistoReconstructed", HistoName,
				xdim, xmin, xMax, zdim, zmin, zMax) ;  
      hModule->SetMaximum(2.0);
      hModule->SetMinimum(0.0);
      hModule->SetStats(kFALSE); 

      TIter next(fPHOS->Digits()) ;
      Float_t Energy, y, z;
      Int_t RelId[4]; Int_t NumberOfDigits = 0 ;
      AliPHOSDigit * digit ; 
      Float_t Etot ; 
      while((digit = (AliPHOSDigit *)next())) 
	{  
	  fGeom->AbsToRelNumbering(digit->GetId(), RelId) ;
	  if (RelId[0] == Module)  
	    {  
	      NumberOfDigits++ ;
	      Energy = fClu->Calibrate(digit->GetAmp()) ;
	      Etot += Energy ; 
	      fGeom->RelPosInModule(RelId,y,z) ; 
	      if (Energy > 0.01 )  
		hModule->Fill(y, z, Energy) ;
	    } 
	}
      cout <<"DrawRecPoints >  Found in Module " 
	   << Module << " " << NumberOfDigits << "  digits with total energy " << Etot << endl ;
      hModule->Draw("col2") ;

      //=========== Cluster in Module

      TClonesArray * EmcRP = fPHOS->EmcClusters() ;
      Etot = 0.; 
      Int_t TotalNumberOfClusters = 0 ; 
      Int_t NumberOfClusters = 0 ;
      TIter nextemc(EmcRP) ;
      AliPHOSEmcRecPoint * emc ;
      while((emc = (AliPHOSEmcRecPoint *)nextemc())) 
	{
	  TotalNumberOfClusters++ ;
	  if ( emc->GetPHOSMod() == Module )
	    { 
	      NumberOfClusters++ ; 
	      Energy = emc->GetTotalEnergy() ;   
	      Etot+= Energy ;  
	      emc->Draw("P") ;
	    }
	}
      cout << "DrawRecPoints > Found " << TotalNumberOfClusters << " EMC Clusters in PHOS" << endl ; 
      cout << "DrawRecPoints > Found in Module " << Module << "  " << NumberOfClusters << " EMC Clusters " << endl ;
      cout << "DrawRecPoints > Total energy  " << Etot << endl ; 

      TPaveText *  PaveText = new TPaveText(22, 80, 83, 90); 
      Text_t text[40] ; 
      sprintf(text, "digits: %d;  clusters: %d", NumberOfDigits, NumberOfClusters) ;
      PaveText->AddText(text) ; 
      PaveText->Draw() ; 
      ModuleCanvas->Update(); 
 
      //=========== Cluster in Module PPSD Down

      TClonesArray * PpsdRP = fPHOS->PpsdClusters() ;
      Etot = 0.; 
      TIter nextPpsd(PpsdRP) ;
      AliPHOSPpsdRecPoint * Ppsd ;
      while((Ppsd = (AliPHOSPpsdRecPoint *)nextPpsd())) 
	{
	  TotalNumberOfClusters++ ;
	  if ( Ppsd->GetPHOSMod() == Module )
	    { 
	      NumberOfClusters++ ; 
	      Energy = Ppsd->GetEnergy() ;   
	      Etot+=Energy ;  
	      if (!Ppsd->GetUp()) Ppsd->Draw("P") ;
	    }
	}
      cout << "DrawRecPoints > Found " << TotalNumberOfClusters << " Ppsd Down Clusters in PHOS" << endl ; 
      cout << "DrawRecPoints > Found in Module " << Module << "  " << NumberOfClusters << " Ppsd Down Clusters " << endl ;
      cout << "DrawRecPoints > Total energy  " << Etot << endl ; 

      //=========== Cluster in Module PPSD Up
  
      PpsdRP = fPHOS->PpsdClusters() ;
      Etot = 0.; 
      TIter nextPpsdUp(PpsdRP) ;
      while((Ppsd = (AliPHOSPpsdRecPoint *)nextPpsdUp())) 
	{
	  TotalNumberOfClusters++ ;
	  if ( Ppsd->GetPHOSMod() == Module )
	    { 
	      NumberOfClusters++ ; 
	      Energy = Ppsd->GetEnergy() ;   
	      Etot+=Energy ;  
	      if (Ppsd->GetUp()) Ppsd->Draw("P") ;
	    }
	}
  cout << "DrawRecPoints > Found " << TotalNumberOfClusters << " Ppsd Up Clusters in PHOS" << endl ; 
  cout << "DrawRecPoints > Found in Module " << Module << "  " << NumberOfClusters << " Ppsd Up Clusters " << endl ;
  cout << "DrawRecPoints > Total energy  " << Etot << endl ; 
    
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

      Int_t Module ; 
      cout <<  "DisplayTrackSegments > which module (1-5,  -1: all) ? " ; 
      cin >> Module ; cout << Module << endl ; 
      //=========== Creating 2d-histogram of the PHOS Module
      // a little bit junkie but is used to test Geom functinalities
      
      Double_t tm, tM, pm, pM ; // min and Max theta and phi covered by Module   
      
      fGeom->EmcModuleCoverage(Module, tm, tM, pm, pM); 
      // convert angles into coordinates local to the EMC module of interest
      
      Int_t EmcModuleNumber ;
      Double_t EmcModulexm, EmcModulezm ; // minimum local coordinate in a given EMCA module
      Double_t EmcModulexM, EmcModulezM ; // maximum local coordinate in a given EMCA module
      fGeom->ImpactOnEmc(tm, pm, EmcModuleNumber, EmcModulezm, EmcModulexm) ;
      fGeom->ImpactOnEmc(tM, pM, EmcModuleNumber, EmcModulezM, EmcModulexM) ;
      Int_t xdim = (Int_t)( ( EmcModulexM - EmcModulexm ) / fGeom->GetCrystalSize(0) ) ;  
      Int_t zdim = (Int_t)( ( EmcModulezM - EmcModulezm ) / fGeom->GetCrystalSize(2) ) ;
      Float_t xmin = EmcModulexm - fGeom->GetCrystalSize(0) ; 
      Float_t xMax = EmcModulexM + fGeom->GetCrystalSize(0) ; 
      Float_t zmin = EmcModulezm - fGeom->GetCrystalSize(2) ; 
      Float_t zMax = EmcModulezM + fGeom->GetCrystalSize(2) ;     
      Text_t HistoName[80];
      sprintf(HistoName,"Event %d: Track Segments in module %d", fEvt, Module) ; 
      TH2F * HistoTrack = new TH2F("HistoTrack",  HistoName, 
				   xdim, xmin, xMax, zdim, zmin, zMax) ;  
      HistoTrack->SetStats(kFALSE); 
      Text_t CanvasName[80];
      sprintf(CanvasName,"Track segments in PHOS/EMC-PPSD module # %d", Module) ;
      TCanvas * TrackCanvas = new TCanvas("TrackSegmentCanvas", CanvasName, 650, 500) ; 
      HistoTrack->Draw() ; 

      TrackSegmentsList * trsegl = fPHOS->TrackSegments() ;
      AliPHOSTrackSegment * trseg ;
 
      Int_t NTrackSegments = trsegl->GetEntries() ;
      Int_t index ;
      Float_t Etot = 0 ;
      Int_t NTrackSegmentsInModule = 0 ; 
      for(index = 0; index < NTrackSegments ; index++){
	trseg = (AliPHOSTrackSegment * )trsegl->At(index) ;
	Etot+= trseg->GetEnergy() ;
	if ( trseg->GetPHOSMod() == Module ) { 
	  NTrackSegmentsInModule++ ; 
	  trseg->Draw("P");
	}
      } 
      Text_t text[80] ; 
      sprintf(text, "track segments: %d", NTrackSegmentsInModule) ;
      TPaveText *  PaveText = new TPaveText(22, 80, 83, 90); 
      PaveText->AddText(text) ; 
      PaveText->Draw() ; 
      TrackCanvas->Update() ; 
      cout << "DisplayTrackSegments > Found " << trsegl->GetEntries() << " Track segments with total energy "<< Etot << endl ;
    
   }
}
//____________________________________________________________________________
Bool_t AliPHOSAnalyze::OpenRootFile(Text_t * name)
{
  fRootFile   = new TFile(name) ;
  return  fRootFile->IsOpen() ; 
}
