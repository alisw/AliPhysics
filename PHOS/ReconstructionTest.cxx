
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
// Cheking Reconstruction procedure of PHOS
//*-- Author : Gines MARTINEZ  SUBATECH january 2000
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TFile.h"
#include "TH1.h"
#include "TPad.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TParticle.h"

// --- Standard library ---

#include <iostream>
#include <cstdio>


// --- AliRoot header files ---
#include "AliRun.h"
#include "TFile.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSv0.h"
#include "AliPHOSDigit.h"
#include "AliPHOSRecPoint.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSPpsdRecPoint.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSReconstructioner.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "ReconstructionTest.h"


void ReconstructionTest(Text_t* infile,Int_t evt, Int_t Module)  
{
  //========== Opening galice.root file  
  TFile * file   = new TFile(infile); 
  //========== Get AliRun object from file 
  gAlice = (AliRun*) file->Get("gAlice");
  //=========== Gets the PHOS object and associated geometry from the file 
  AliPHOSv0 * PHOS  = (AliPHOSv0 *)gAlice->GetDetector("PHOS");
  AliPHOSGeometry * Geom = AliPHOSGeometry::GetInstance(PHOS->GetGeometry()->GetName(),PHOS->GetGeometry()->GetTitle());
  //========== Creates the Clusterizer
  AliPHOSClusterizerv1 clusterizer; 
  clusterizer.SetEmcEnergyThreshold(0.01) ; 
  clusterizer.SetEmcClusteringThreshold(0.1) ; 
  clusterizer.SetPpsdEnergyThreshold(0.0000001) ; 
  clusterizer.SetPpsdClusteringThreshold(0.0000002) ; 
  clusterizer.SetLocalMaxCut(0.03) ;
  clusterizer.SetCalibrationParameters(0., 0.0000001) ;
  //========== Creates the track segment maker
  AliPHOSTrackSegmentMakerv1 tracksegmentmaker ;
  //========== Creates the Reconstructioner  
  AliPHOSReconstructioner Reconstructioner(clusterizer,tracksegmentmaker);     
  //=========== Connects the various Tree's for evt
  gAlice->GetEvent(evt);
  //=========== Gets the Digit TTree
  gAlice->TreeD()->GetEvent(0) ;     
  //=========== Gets the number of entries in the Digits array
  Int_t nId = PHOS->Digits()->GetEntries();    
  printf("ReconstructionTest> Number of entries in the Digit array is %d \n",nId);
  //=========== Do the reconstruction
  AliPHOSDigit * digit ;
  TIter next(PHOS->Digits()) ;
  Float_t Etot=0 ;
  while( ( digit = (AliPHOSDigit *)next() ) )
    Etot+=clusterizer.Calibrate(digit->GetAmp()) ;
  cout <<"ReconstructionTest> Found  " << nId << "  digits in PHOS with total energy " << Etot << endl ;  
  PHOS->Reconstruction(Reconstructioner);  
  //================Make checks===========================
  
  //=========== Creating Canvas
  TCanvas * ModuleCanvas = new TCanvas("Module","Events in a single PHOS Module", 650, 500) ; 
  ModuleCanvas->Draw() ;
 
  //=========== Creating 2d-histogram of the PHOS Module
  // a little bit junkie but is used to test Geom functinalities

  Double_t tm, tM, pm, pM ; // min and Max theta and phi covered by Module 1  
 
  Geom->EmcModuleCoverage(1, tm, tM, pm, pM); 
  // convert angles into coordinates local to the EMC module of interest

  Int_t EmcModuleNumber ;
  Double_t EmcModulexm, EmcModulezm ; // minimum local coordinate in a given EMCA module
  Double_t EmcModulexM, EmcModulezM ; // maximum local coordinate in a given EMCA module
  Geom->ImpactOnEmc(tm, pm, EmcModuleNumber, EmcModulezm, EmcModulexm) ;
  Geom->ImpactOnEmc(tM, pM, EmcModuleNumber, EmcModulezM, EmcModulexM) ;
  Int_t xdim = (Int_t)( ( EmcModulexM - EmcModulexm ) / Geom->GetCrystalSize(0) ) ;  
  Int_t zdim = (Int_t)( ( EmcModulezM - EmcModulezm ) / Geom->GetCrystalSize(2) ) ;
  Float_t xmin = EmcModulexm - Geom->GetCrystalSize(0) ; 
  Float_t xMax = EmcModulexM + Geom->GetCrystalSize(0) ; 
  Float_t zmin = EmcModulezm - Geom->GetCrystalSize(2) ; 
  Float_t zMax = EmcModulezM + Geom->GetCrystalSize(2) ; 
  // histogram of reconstructed events
  Text_t HistoName[80];
  sprintf(HistoName,"Event %d: Reconstructed particles in module %d", evt, Module) ;
  TH2F * hModule = new TH2F("HistoReconstructed", HistoName,
			    xdim, xmin, xMax, zdim, zmin, zMax) ;  
  hModule->SetMaximum(2.0);
  hModule->SetMinimum(0.0);
  hModule->SetStats(kFALSE);
  // histogram of generated  particles
  sprintf(HistoName,"Event %d: Incident particles in module %d", evt, Module) ; 
  TH2F * HistoParticle = new TH2F("HistoParticle",  HistoName, 
				  xdim, xmin, xMax, zdim, zmin, zMax) ;  
  HistoParticle->SetStats(kFALSE) ;

  //=========== Digits of Module 
  TIter next2(PHOS->Digits()) ;
  Float_t Energy, y, z;
  Int_t RelId[4]; Int_t NumberOfDigits = 0 ;
  while((digit = (AliPHOSDigit *)next2())) 
    {  
      Geom->AbsToRelNumbering(digit->GetId(), RelId) ;
      if (RelId[0] == Module)  
	{  
	  NumberOfDigits++ ;
	  Energy = clusterizer.Calibrate(digit->GetAmp()) ;
	  Etot+=Energy ; 
	  Geom->RelPosInModule(RelId,y,z) ; 
	  if (Energy>0.01 )  
	    hModule->Fill(y,z,Energy) ;
	} 
  }
  cout <<"TestRec> Found in Module " << Module << " " << NumberOfDigits << "  digits with total energy " << Etot << endl ;
  hModule->Draw("col2") ;
  
  //=========== Cluster in Module
  TClonesArray * EmcRP = PHOS->EmcClusters() ;
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
	  Etot+=Energy ;  
	  emc->Draw("P") ;
    }
    }
  cout << "TestRec> Found " << TotalNumberOfClusters << " EMC Clusters in PHOS" << endl ; 
  cout << "TestRec> Found in Module " << Module << "  " << NumberOfClusters << " EMC Clusters " << endl ;
  cout << "TestRec> Total energy  " <<Etot << endl ; 

  TPaveText *  PaveText = new TPaveText(22, 80, 83, 90); 
  Text_t text[40] ; 
  sprintf(text, "digits: %d;  clusters: %d", NumberOfDigits, NumberOfClusters) ;
  PaveText->AddText(text) ; 
  PaveText->Draw() ; 
  ModuleCanvas->Update(); 
 
  //=========== Cluster in Module PPSD Down
  TClonesArray * PpsdRP = PHOS->PpsdClusters() ;
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
  cout << "TestRec> Found " << TotalNumberOfClusters << " Ppsd Down Clusters in PHOS" << endl ; 
  cout << "TestRec> Found in Module " << Module << "  " << NumberOfClusters << " Ppsd Down Clusters " << endl ;
  cout << "TestRec> Total energy  " <<Etot << endl ; 

//=========== Cluster in Module PPSD Up
  PpsdRP = PHOS->PpsdClusters() ;
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
  cout << "TestRec> Found " << TotalNumberOfClusters << " Ppsd Up Clusters in PHOS" << endl ; 
  cout << "TestRec> Found in Module " << Module << "  " << NumberOfClusters << " Ppsd Up Clusters " << endl ;
  cout << "TestRec> Total energy  " <<Etot << endl ; 
  
   // Get pointers to Alice Particle TClonesArray
  TParticle * Particle;
  TClonesArray * ArrayOfParticles  = gAlice->Particles();    
  TCanvas * KineCanvas = new TCanvas("KineCnvas", "Incident particles", 650, 500) ; 
  // get the KINE Tree
  TTree * Kine =  gAlice->TreeK() ; 
  Stat_t NumberOfParticles =  Kine->GetEntries() ; 
  cout << "events in Kine " << NumberOfParticles << endl ; 
  
  // loop over particles
  Int_t index1 ; 
  Int_t nparticlein = 0 ; 
  for (index1 = 0 ; index1 < NumberOfParticles ; index1++){
    Int_t nparticle = ArrayOfParticles->GetEntriesFast() ;
    cout << nparticle << endl ; 
    Int_t index2 ; 
    for( index2 = 0 ; index2 < nparticle; index2++) {         
      Particle            = (TParticle*)ArrayOfParticles->UncheckedAt(index2) ;
      Int_t  ParticleType = Particle->GetPdgCode() ;
      Double_t Phi        = Particle->Phi() ;
      Double_t Theta      = Particle->Theta() ;
      Int_t mod ; 
      Double_t x, z ; 
      Geom->ImpactOnEmc(Theta, Phi, mod, z, x) ;
      if ( mod == Module ) {
	nparticlein++ ; 
	HistoParticle->Fill( x, -z, Particle->Energy() ) ; //-z don't know why, but that is how it works
      } 
    } 
  }
  KineCanvas->Draw() ; 
  HistoParticle->Draw("color") ; 
  TPaveText *  PaveText2 = new TPaveText(22, 80, 83, 90); 
  sprintf(text, "Particles: %d ", nparticlein) ;
  PaveText2->AddText(text) ; 
  PaveText2->Draw() ; 
  KineCanvas->Update(); 

//  TObjArray * trsegl = PHOS->TrackSegments() ;
//    AliPHOSTrackSegment trseg ;
 
//    Int_t NTrackSegments = trsegl->GetEntries() ;
//    Int_t index ;
//    Etot = 0 ;
//    for(index = 0; index < NTrackSegments ; index++){
//      trseg = (AliPHOSTrackSegment * )trsegl->At(index) ;
//      Etot+= trseg->GetEnergy() ;
//      if ( trseg->GetPHOSMod() == Module )
//      { 
// //       trseg->Draw("P");
//        trseg->Print() ;
//      }
//    } 
//    cout << "Found " << trsegl->GetEntries() << " Track segments with total energy "<< Etot << endl ;
//  

 
  
}


