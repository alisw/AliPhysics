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
// Algorythm class to analyze PHOSv1 events:
// Construct histograms and displays them.
// Use the macro EditorBar.C for best access to the functionnalities
//*--
//*-- Author: Y. Schutz (SUBATECH) & Gines Martinez (SUBATECH)
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
#include "TStyle.h" 

// --- Standard library ---

#include <iostream.h>
#include <stdio.h>

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
#include "AliPHOSIndexToObject.h"
#include "AliPHOSCPVHit.h"
#include "AliPHOSCpvRecPoint.h"

ClassImp(AliPHOSAnalyze)

//____________________________________________________________________________
  AliPHOSAnalyze::AliPHOSAnalyze()
{
  // default ctor (useless)
  
  fRootFile = 0 ; 
}

//____________________________________________________________________________
AliPHOSAnalyze::AliPHOSAnalyze(Text_t * name)
{
  // ctor: analyze events from root file "name"
  
  Bool_t ok = OpenRootFile(name)  ; 
  if ( !ok ) {
    cout << " AliPHOSAnalyze > Error opening " << name << endl ; 
  }
  else { 
      //========== Get AliRun object from file 
      gAlice = (AliRun*) fRootFile->Get("gAlice") ;

      //=========== Get the PHOS object and associated geometry from the file      
      fPHOS  = (AliPHOSv1 *)gAlice->GetDetector("PHOS") ;
      fGeom  = AliPHOSGeometry::GetInstance( fPHOS->GetGeometry()->GetName(), fPHOS->GetGeometry()->GetTitle() );
 
      //========== Initializes the Index to Object converter
      fObjGetter = AliPHOSIndexToObject::GetInstance(fPHOS) ; 
      //========== Current event number 
      fEvt = -999 ; 

  }
  fDebugLevel = 0;
  fClu = 0 ;
  fPID = 0 ;
  fTrs = 0 ;
  fRec = 0 ;
  ResetHistograms() ;
}

//____________________________________________________________________________
AliPHOSAnalyze::AliPHOSAnalyze(const AliPHOSAnalyze & ana)
{
  // copy ctor
  ( (AliPHOSAnalyze &)ana ).Copy(*this) ;
}

//____________________________________________________________________________
void AliPHOSAnalyze::Copy(TObject & obj)
{
  // copy an analysis into an other one
  TObject::Copy(obj) ;
  // I do nothing more because the copy is silly but the Code checkers requires one
}

//____________________________________________________________________________
AliPHOSAnalyze::~AliPHOSAnalyze()
{
  // dtor

  if(fRootFile->IsOpen()) fRootFile->Close() ; 
  if(fRootFile)   {delete fRootFile ; fRootFile=0 ;}
  if(fPHOS)       {delete fPHOS     ; fPHOS    =0 ;}
  if(fClu)        {delete fClu      ; fClu     =0 ;}
  if(fPID)        {delete fPID      ; fPID     =0 ;}
  if(fRec)        {delete fRec      ; fRec     =0 ;}
  if(fTrs)        {delete fTrs      ; fTrs     =0 ;}

}

//____________________________________________________________________________
void AliPHOSAnalyze::ActivePPSD(Int_t Nevents=1){
  
  fhEnergyCorrelations  = new TH2F("hEnergyCorrelations","hEnergyCorrelations",40,  0., 0.15, 30, 0., 3.e-5);
  //========== Create the Clusterizer
  fClu = new AliPHOSClusterizerv1() ; 
  fClu->SetEmcEnergyThreshold(0.01) ; 
  fClu->SetEmcClusteringThreshold(0.20) ; 
  fClu->SetPpsdEnergyThreshold    (0.0000002) ; 
  fClu->SetPpsdClusteringThreshold(0.0000001) ; 
  fClu->SetLocalMaxCut(0.02) ;
  fClu->SetCalibrationParameters(0., 0.00000001) ; 

  Int_t ievent;
  
  for ( ievent=0; ievent<Nevents; ievent++)
    {  
      
      //========== Event Number>         
      if ( ( log10((Float_t)(ievent+1)) - (Int_t)(log10((Float_t)(ievent+1))) ) == 0. ) 
	cout <<  "AnalyzeResolutions > " << "Event is " << ievent << endl ;  
      
      //=========== Connects the various Tree's for evt
      gAlice->GetEvent(ievent);

      //=========== Gets the Kine TTree
      gAlice->TreeK()->GetEvent(0) ;
      
      //=========== Get the Digit Tree
      gAlice->TreeD()->GetEvent(0) ;
      
      //========== Creating branches ===================================       
      AliPHOSRecPoint::RecPointsList ** EmcRecPoints =  fPHOS->EmcRecPoints() ;
      gAlice->TreeR()->SetBranchAddress( "PHOSEmcRP", EmcRecPoints ) ;
      
      AliPHOSRecPoint::RecPointsList ** PpsdRecPoints = fPHOS->PpsdRecPoints() ;
      gAlice->TreeR()->SetBranchAddress( "PHOSPpsdRP", PpsdRecPoints ) ;
      
      AliPHOSTrackSegment::TrackSegmentsList **  TrackSegmentsList = fPHOS->TrackSegments() ;
      if( (*TrackSegmentsList) )
	(*TrackSegmentsList)->Clear() ;
      gAlice->TreeR()->SetBranchAddress( "PHOSTS", TrackSegmentsList ) ;
      
      AliPHOSRecParticle::RecParticlesList ** RecParticleList  = fPHOS->RecParticles() ;
      if( (*RecParticleList) )
	(*RecParticleList)->Clear() ;
      gAlice->TreeR()->SetBranchAddress( "PHOSRP", RecParticleList ) ;
      
      
      //=========== Gets the Reconstraction TTree
      gAlice->TreeR()->GetEvent(0) ;
            
      AliPHOSPpsdRecPoint * RecPoint ;
      Int_t relid[4] ; 
      TIter nextRP(*fPHOS->PpsdRecPoints() ) ;
      while( ( RecPoint = (AliPHOSPpsdRecPoint *)nextRP() ) )
	{
 	  if(!(RecPoint->GetUp()) ) {
	    AliPHOSDigit *digit ;
	    Int_t iDigit ;
	    for(iDigit = 0; iDigit < fPHOS->Digits()->GetEntries(); iDigit++) 
	      {
		digit = (AliPHOSDigit *) fPHOS->Digits()->At(iDigit) ;
		fGeom->AbsToRelNumbering(digit->GetId(), relid) ;    
                if((relid[2]==1)&&(relid[3]==1)&&(relid[0]==RecPoint->GetPHOSMod())){
		  Float_t ConvertorEnergy = fClu->Calibrate(digit->GetAmp()) ;
		  fhEnergyCorrelations->Fill(ConvertorEnergy,RecPoint->GetTotalEnergy() );  
		  break ; 
		} 
	      }
	    break ;
 	  }
 	}
    }
  SaveHistograms() ;
  fhEnergyCorrelations->Draw("BOX") ;
}


//____________________________________________________________________________
void AliPHOSAnalyze::AnalyzeManyEvents(Int_t Nevents, Int_t module)    
{
  // analyzes Nevents events in a single PHOS module  
  // Events should be reconstructed by Reconstruct()

  if ( fRootFile == 0 ) 
    cout << "AnalyzeManyEvents > " << "Root File not openned" << endl ;  
  else
    {
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
	  //========== Event Number>         
	  if ( ( log10((Float_t)(ievent+1)) - (Int_t)(log10((Float_t)(ievent+1))) ) == 0. ) 
	    cout <<  "AnalyzeManyEvents > " << "Event is " << ievent << endl ;  

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
	  

	  //=========== Cluster in module
	  TIter nextEmc(*fPHOS->EmcRecPoints()  ) ;
	  while((emc = (AliPHOSEmcRecPoint *)nextEmc())) 
	    {
	      if ( emc->GetPHOSMod() == module )
		{  
		  fhEmcCluster->Fill(  emc->GetTotalEnergy()  ); 
		  TIter nextPpsd( *fPHOS->PpsdRecPoints()) ;
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
	  TIter nextPpsd(*fPHOS->PpsdRecPoints() ) ;
	  while((ppsd = (AliPHOSPpsdRecPoint *)nextPpsd())) 
	    {
	      if ( ppsd->GetPHOSMod() == module )
		{ 
		  if (!ppsd->GetUp()) fhConvertorCluster->Fill(ppsd->GetTotalEnergy()) ;
		  if (ppsd->GetUp())  fhVetoCluster     ->Fill(ppsd->GetTotalEnergy()) ;
		}
	    }

	  //========== TRackSegments in the event
	  TIter nextRecParticle(*fPHOS->RecParticles() ) ; 
	  while((recparticle = (AliPHOSRecParticle *)nextRecParticle())) 
	    {
	      if ( recparticle->GetPHOSTrackSegment()->GetPHOSMod() == module )
		{ 
		  cout << "Particle type is " << recparticle->GetType() << endl ; 
		  Int_t numberofprimaries = 0 ;
		  Int_t * listofprimaries = recparticle->GetPrimaries(numberofprimaries) ;
		  cout << "Number of primaries = " << numberofprimaries << endl ; 
		  Int_t index ;
		  for ( index = 0 ; index < numberofprimaries ; index++)
		    cout << "    primary # " << index << " =  " << listofprimaries[index] << endl ;  
		}
	    }
	}   // endfor
      SaveHistograms();
    }       // endif
}           // endfunction

//____________________________________________________________________________
 void AliPHOSAnalyze::Reconstruct(Int_t Nevents,Int_t FirstEvent )    
{     

  // Performs reconstruction of EMC and CPV (GPS2 or IHEP)
  // for events from FirstEvent to Nevents

  Int_t ievent ;   
  for ( ievent=FirstEvent; ievent<Nevents; ievent++) {  
    if (ievent==FirstEvent) {
      cout << "Analyze > Starting Reconstructing " << endl ; 
      //========== Create the Clusterizer
      fClu = new AliPHOSClusterizerv1() ; 
      fClu->SetEmcEnergyThreshold(0.05) ; 
      fClu->SetEmcClusteringThreshold(0.20) ; 
      fClu->SetLocalMaxCut(0.03) ;
      if      (strcmp(fGeom->GetName(),"GPS2") == 0) {
	fClu->SetPpsdEnergyThreshold    (0.0000002) ; 
	fClu->SetPpsdClusteringThreshold(0.0000001) ; 
      }
      else if (strcmp(fGeom->GetName(),"IHEP") == 0) {
	fClu->SetLocalMaxCutCPV(0.03) ;
	fClu->SetLogWeightCutCPV(4.0) ;
	fClu->SetPpsdEnergyThreshold    (0.09) ;
      }
      fClu->SetCalibrationParameters(0., 0.00000001) ; 
      
      //========== Creates the track segment maker
      fTrs = new AliPHOSTrackSegmentMakerv1()  ;
 	  //	  fTrs->UnsetUnfoldFlag() ; 
     
      //========== Creates the particle identifier for GPS2 only
      if      (strcmp(fGeom->GetName(),"GPS2") == 0) {
	fPID = new AliPHOSPIDv1() ;
	fPID->SetShowerProfileCuts(0.3, 1.8, 0.3, 1.8 ) ; 
      }	  
      
      //========== Creates the Reconstructioner
      fRec = new AliPHOSReconstructioner(fClu, fTrs, fPID) ; 
      if (fDebugLevel != 0) fRec -> SetDebugReconstruction(kTRUE);     
    }
      
    if (fDebugLevel != 0 ||
	(ievent+1) % (Int_t)TMath::Power( 10, (Int_t)TMath::Log10(ievent+1) ) == 0)
      cout <<  "======= Analyze ======> Event " << ievent+1 << endl ;
    
    //=========== Connects the various Tree's for evt
    gAlice->GetEvent(ievent);
    
    //=========== Gets the Digit TTree
    gAlice->TreeD()->GetEvent(0) ;
    
    //=========== Do the reconstruction
    fPHOS->Reconstruction(fRec);
  }

  if(fClu)      {delete fClu      ; fClu     =0 ;}
  if(fPID)      {delete fPID      ; fPID     =0 ;}
  if(fRec)      {delete fRec      ; fRec     =0 ;}
  if(fTrs)      {delete fTrs      ; fTrs     =0 ;}
  
}

//-------------------------------------------------------------------------------------
void AliPHOSAnalyze::ReadAndPrintCPV(Int_t EvFirst, Int_t EvLast)
{
  //
  // Read and print generated and reconstructed hits in CPV
  // for events from EvFirst to Nevent.
  // If only EvFirst is defined, print only this one event.
  // Author: Yuri Kharlov
  // 12 October 2000
  //

  if (EvFirst!=0 && EvLast==0) EvLast=EvFirst;
  for ( Int_t ievent=EvFirst; ievent<=EvLast; ievent++) {  
    
    //========== Event Number>
    cout << endl <<  "==== ReadAndPrintCPV ====> Event is " << ievent+1 << endl ;
    
    //=========== Connects the various Tree's for evt
    Int_t ntracks = gAlice->GetEvent(ievent);

    //========== Creating branches ===================================
    AliPHOSRecPoint::RecPointsList ** emcRecPoints = fPHOS->EmcRecPoints() ;
    gAlice->TreeR()->SetBranchAddress( "PHOSEmcRP" , emcRecPoints  ) ;
    
    AliPHOSRecPoint::RecPointsList ** cpvRecPoints = fPHOS->PpsdRecPoints() ;
    gAlice->TreeR()->SetBranchAddress( "PHOSPpsdRP", cpvRecPoints ) ;

    // Read and print CPV hits
      
    AliPHOSCPVModule cpvModule;
    TClonesArray    *cpvHits;
    Int_t           nCPVhits;
    AliPHOSCPVHit   *cpvHit;
    TLorentzVector   p;
    Float_t          xgen, zgen;
    Int_t            ipart;
    Int_t            nGenHits = 0;
    for (Int_t itrack=0; itrack<ntracks; itrack++) {
      //=========== Get the Hits Tree for the Primary track itrack
      gAlice->ResetHits();
      gAlice->TreeH()->GetEvent(itrack);
      for (Int_t iModule=0; iModule < fGeom->GetNModules(); iModule++) {
	cpvModule = fPHOS->GetCPVModule(iModule);
	cpvHits   = cpvModule.Hits();
	nCPVhits  = cpvHits->GetEntriesFast();
	for (Int_t ihit=0; ihit<nCPVhits; ihit++) {
	  nGenHits++;
	  cpvHit = (AliPHOSCPVHit*)cpvHits->UncheckedAt(ihit);
	  p      = cpvHit->GetMomentum();
	  xgen   = cpvHit->GetX();
	  zgen   = cpvHit->GetY();
	  ipart  = cpvHit->GetIpart();
	  printf("CPV hit in module %d: ",iModule+1);
	  printf(" p = (%f, %f, %f, %f) GeV,\n",
		 p.Px(),p.Py(),p.Pz(),p.Energy());
	  printf("                  (X,Z) = (%8.4f, %8.4f) cm, ipart = %d\n",
		 xgen,zgen,ipart);
	}
      }
    }

    // Read and print CPV reconstructed points

    //=========== Gets the Reconstruction TTree
    gAlice->TreeR()->GetEvent(0) ;
    TIter nextRP(*fPHOS->PpsdRecPoints() ) ;
    AliPHOSPpsdRecPoint *cpvRecPoint ;
    Int_t nRecPoints = 0;
    while( ( cpvRecPoint = (AliPHOSPpsdRecPoint *)nextRP() ) ) {
      nRecPoints++;
      TVector3  locpos;
      cpvRecPoint->GetLocalPosition(locpos);
      Int_t phosModule = cpvRecPoint->GetPHOSMod();
      printf("CPV recpoint in module %d: (X,Z) = (%f,%f) cm\n",
	     phosModule,locpos.X(),locpos.Z());
    }
    printf("This event has %d generated hits and %d reconstructed points\n",
	   nGenHits,nRecPoints);
  }
}

//____________________________________________________________________________
void AliPHOSAnalyze::AnalyzeCPV(Int_t Nevents)
{
  //
  // Analyzes CPV characteristics
  // Author: Yuri Kharlov
  // 9 October 2000
  //

  // Book histograms

  TH1F *hDx   = new TH1F("hDx"  ,"CPV x-resolution@reconstruction",100,-5. , 5.);
  TH1F *hDz   = new TH1F("hDz"  ,"CPV z-resolution@reconstruction",100,-5. , 5.);
  TH1F *hDr   = new TH1F("hDr"  ,"CPV r-resolution@reconstruction",100, 0. , 5.);
  TH1S *hNrp  = new TH1S("hNrp" ,"CPV rec.point multiplicity",      21,-0.5,20.5);
  TH1S *hNrpX = new TH1S("hNrpX","CPV rec.point Phi-length"  ,      21,-0.5,20.5);
  TH1S *hNrpZ = new TH1S("hNrpZ","CPV rec.point Z-length"    ,      21,-0.5,20.5);

  cout << "Start CPV Analysis"<< endl ;
  for ( Int_t ievent=0; ievent<Nevents; ievent++) {  
      
    //========== Event Number>         
//      if ( (ievent+1) % (Int_t)TMath::Power( 10, (Int_t)TMath::Log10(ievent+1) ) == 0)
      cout << endl <<  "==== AnalyzeCPV ====> Event is " << ievent+1 << endl ;
    
    //=========== Connects the various Tree's for evt
    Int_t ntracks = gAlice->GetEvent(ievent);
    
    //========== Creating branches ===================================
    AliPHOSRecPoint::RecPointsList ** emcRecPoints = fPHOS->EmcRecPoints() ;
    gAlice->TreeR()->SetBranchAddress( "PHOSEmcRP" , emcRecPoints  ) ;
    
    AliPHOSRecPoint::RecPointsList ** cpvRecPoints = fPHOS->PpsdRecPoints() ;
    gAlice->TreeR()->SetBranchAddress( "PHOSPpsdRP", cpvRecPoints ) ;

    // Create and fill arrays of hits for each CPV module
      
    Int_t nOfModules = fGeom->GetNModules();
    TClonesArray *hitsPerModule[nOfModules];
    for (Int_t iModule=0; iModule < nOfModules; iModule++)
      hitsPerModule[iModule] = new TClonesArray("AliPHOSCPVHit",100);

    AliPHOSCPVModule cpvModule;
    TClonesArray    *cpvHits;
    Int_t           nCPVhits;
    AliPHOSCPVHit   *cpvHit;
    TLorentzVector   p;
    Float_t          xzgen[2];
    Int_t            ipart;

    // First go through all primary tracks and fill the arrays
    // of hits per each CPV module

    for (Int_t itrack=0; itrack<ntracks; itrack++) {
      // Get the Hits Tree for the Primary track itrack
      gAlice->ResetHits();
      gAlice->TreeH()->GetEvent(itrack);
      for (Int_t iModule=0; iModule < nOfModules; iModule++) {
	cpvModule = fPHOS->GetCPVModule(iModule);
	cpvHits   = cpvModule.Hits();
	nCPVhits  = cpvHits->GetEntriesFast();
	for (Int_t ihit=0; ihit<nCPVhits; ihit++) {
	  cpvHit   = (AliPHOSCPVHit*)cpvHits->UncheckedAt(ihit);
	  p        = cpvHit->GetMomentum();
	  xzgen[0] = cpvHit->GetX();
	  xzgen[1] = cpvHit->GetY();
	  ipart    = cpvHit->GetIpart();
	  TClonesArray &lhits = *(TClonesArray *)hitsPerModule[iModule];
//  	  new(lhits[hitsPerModule[iModule]->GetEntriesFast()]) AliPHOSCPVHit(p,xzgen,ipart);
	  new(lhits[hitsPerModule[iModule]->GetEntriesFast()]) AliPHOSCPVHit(*cpvHit);
	}
	cpvModule.Clear();
      }
    }
    for (Int_t iModule=0; iModule < nOfModules; iModule++) {
      Int_t nsum = hitsPerModule[iModule]->GetEntriesFast();
      printf("Module %d has %d hits\n",iModule,nsum);
    }

    // Then go through reconstructed points and for each find
    // the closeset hit
    // The distance from the rec.point to the closest hit
    // gives the coordinate resolution of the CPV

    // Get the Reconstruction Tree
    gAlice->TreeR()->GetEvent(0) ;
    TIter nextRP(*fPHOS->PpsdRecPoints() ) ;
    AliPHOSCpvRecPoint *cpvRecPoint ;
    Float_t xgen, zgen;
    while( ( cpvRecPoint = (AliPHOSCpvRecPoint *)nextRP() ) ) {
      TVector3  locpos;
      cpvRecPoint->GetLocalPosition(locpos);
      Int_t phosModule = cpvRecPoint->GetPHOSMod();
      Int_t rpMult     = cpvRecPoint->GetDigitsMultiplicity();
      Int_t rpMultX, rpMultZ;
      cpvRecPoint->GetClusterLengths(rpMultX,rpMultZ);
      Float_t xrec  = locpos.X();
      Float_t zrec  = locpos.Z();
      Float_t dxmin = 1.e+10;
      Float_t dzmin = 1.e+10;
      Float_t r2min = 1.e+10;
      Float_t r2;

      cpvHits = hitsPerModule[phosModule-1];
      Int_t nCPVhits  = cpvHits->GetEntriesFast();
      for (Int_t ihit=0; ihit<nCPVhits; ihit++) {
	cpvHit = (AliPHOSCPVHit*)cpvHits->UncheckedAt(ihit);
	xgen   = cpvHit->GetX();
	zgen   = cpvHit->GetY();
	r2 = TMath::Power((xgen-xrec),2) + TMath::Power((zgen-zrec),2);
	if ( r2 < r2min ) {
	  r2min = r2;
	  dxmin = xgen - xrec;
	  dzmin = zgen - zrec;
	}
      }
      hDx  ->Fill(dxmin);
      hDz  ->Fill(dzmin);
      hDr  ->Fill(TMath::Sqrt(r2min));
      hNrp ->Fill(rpMult);
      hNrpX->Fill(rpMultX);
      hNrpZ->Fill(rpMultZ);
    }
  }
  // Save histograms

  Text_t outputname[80] ;
  sprintf(outputname,"%s.analyzed",fRootFile->GetName());
  TFile output(outputname,"RECREATE");
  output.cd();

  hDx  ->Write() ;
  hDz  ->Write() ;
  hDr  ->Write() ;
  hNrp ->Write() ;
  hNrpX->Write() ;
  hNrpZ->Write() ;

  // Plot histograms

  TCanvas *cpvCanvas = new TCanvas("CPV","CPV analysis",20,20,800,400);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1);
  gStyle->SetOptDate(1);
  cpvCanvas->Divide(3,2);

  cpvCanvas->cd(1);
  gPad->SetFillColor(10);
  hNrp->SetFillColor(16);
  hNrp->Draw();

  cpvCanvas->cd(2);
  gPad->SetFillColor(10);
  hNrpX->SetFillColor(16);
  hNrpX->Draw();

  cpvCanvas->cd(3);
  gPad->SetFillColor(10);
  hNrpZ->SetFillColor(16);
  hNrpZ->Draw();

  cpvCanvas->cd(4);
  gPad->SetFillColor(10);
  hDx->SetFillColor(16);
  hDx->Fit("gaus");
  hDx->Draw();

  cpvCanvas->cd(5);
  gPad->SetFillColor(10);
  hDz->SetFillColor(16);
  hDz->Fit("gaus");
  hDz->Draw();

  cpvCanvas->cd(6);
  gPad->SetFillColor(10);
  hDr->SetFillColor(16);
  hDr->Draw();

  cpvCanvas->Print("CPV.ps");

}

//____________________________________________________________________________
 void AliPHOSAnalyze::InvariantMass(Int_t Nevents )    
{
  // Calculates Real and Mixed invariant mass distributions
  Int_t NMixedEvents = 4 ; //# of events used for calculation of 'mixed' distribution 
  Int_t MixedLoops = (Int_t )TMath::Ceil(Nevents/NMixedEvents) ;
  
  //========== Booking Histograms
  TH2D * hRealEM   = new TH2D("hRealEM",   "Real for EM particles",      250,0.,1.,40,0.,4.) ;
  TH2D * hRealPhot = new TH2D("hRealPhot", "Real for kPhoton particles", 250,0.,1.,40,0.,4.) ;
  TH2D * hMixedEM  = new TH2D("hMixedEM",  "Mixed for EM particles",     250,0.,1.,40,0.,4.) ;
  TH2D * hMixedPhot= new TH2D("hMixedPhot","Mixed for kPhoton particles",250,0.,1.,40,0.,4.) ;
  
  Int_t ievent;
  Int_t EventInMixedLoop ;
  
  Int_t NRecParticles[NMixedEvents] ;
  
  AliPHOSRecParticle::RecParticlesList * AllRecParticleList  = new TClonesArray("AliPHOSRecParticle", NMixedEvents*1000) ;
  
  for(EventInMixedLoop = 0; EventInMixedLoop < MixedLoops; EventInMixedLoop++  ){
    Int_t iRecPhot = 0 ;
    
    for ( ievent=0; ievent < NMixedEvents; ievent++){        
      
      Int_t AbsEventNumber = EventInMixedLoop*NMixedEvents + ievent ;
      
      //=========== Connects the various Tree's for evt
      gAlice->GetEvent(AbsEventNumber);
      
      //=========== Get the Digit Tree
      gAlice->TreeD()->GetEvent(0) ;
      
      //========== Creating branches ===================================       
      
      AliPHOSRecParticle::RecParticlesList ** RecParticleList  = fPHOS->RecParticles() ;
      if( (*RecParticleList) )
	(*RecParticleList)->Clear() ;
      gAlice->TreeR()->SetBranchAddress( "PHOSRP", RecParticleList ) ;
      
      //=========== Gets the Reconstraction TTree
      gAlice->TreeR()->GetEvent(0) ;
      
      AliPHOSRecParticle * RecParticle ;
      Int_t iRecParticle ;
      for(iRecParticle = 0; iRecParticle < (*RecParticleList)->GetEntries() ;iRecParticle++ )
 	{
 	  RecParticle = (AliPHOSRecParticle *) (*RecParticleList)->At(iRecParticle) ;
	  if((RecParticle->GetType() == AliPHOSFastRecParticle::kGAMMA)||
	     (RecParticle->GetType() == AliPHOSFastRecParticle::kNEUTRALEM)){ 
	    new( (*AllRecParticleList)[iRecPhot] ) AliPHOSRecParticle(*RecParticle) ;
	    iRecPhot++;
	  }
	}
      
	NRecParticles[ievent] = iRecPhot-1 ;  
    }
    
    //Now calculate invariant mass:
    Int_t irp1,irp2 ;
    Int_t NCurEvent = 0 ;

    for(irp1 = 0; irp1 < AllRecParticleList->GetEntries()-1; irp1++){
      AliPHOSRecParticle * rp1 = (AliPHOSRecParticle *)AllRecParticleList->At(irp1) ;

      for(irp2 = irp1+1; irp2 < AllRecParticleList->GetEntries(); irp2++){
	AliPHOSRecParticle * rp2 = (AliPHOSRecParticle *)AllRecParticleList->At(irp2) ;
	    
	Double_t InvMass ;
	InvMass = (rp1->Energy()+rp2->Energy())*(rp1->Energy()+rp2->Energy())-
	  (rp1->Px()+rp2->Px())*(rp1->Px()+rp2->Px())-
	  (rp1->Py()+rp2->Py())*(rp1->Py()+rp2->Py())-
	  (rp1->Pz()+rp2->Pz())*(rp1->Pz()+rp2->Pz()) ;
	
	if(InvMass> 0)
	  InvMass = TMath::Sqrt(InvMass);
	
	Double_t Pt ; 
	Pt = TMath::Sqrt((rp1->Px()+rp2->Px() )*( rp1->Px()+rp2->Px() ) +(rp1->Py()+rp2->Py())*(rp1->Py()+rp2->Py()));

	if(irp1 > NRecParticles[NCurEvent])
	  NCurEvent++;
	    
	if(irp2 <= NRecParticles[NCurEvent]){ //'Real' event
	  hRealEM->Fill(InvMass,Pt);
	  if((rp1->GetType() == AliPHOSFastRecParticle::kGAMMA)&&(rp2->GetType() == AliPHOSFastRecParticle::kGAMMA))
	    hRealPhot->Fill(InvMass,Pt);
	}
	else{
	  hMixedEM->Fill(InvMass,Pt);
	  if((rp1->GetType() == AliPHOSFastRecParticle::kGAMMA)&&(rp2->GetType() == AliPHOSFastRecParticle::kGAMMA))
	    hMixedPhot->Fill(InvMass,Pt);
	} //real-mixed
	    
      } //loop over second rp
    }//loop over first rp
    AllRecParticleList->Delete() ;
  } //Loop over events
  
  delete AllRecParticleList ;
  
  //writing output
  TFile output("invmass.root","RECREATE");
  output.cd();
  
  hRealEM->Write() ;
  hRealPhot->Write() ;
  hMixedEM->Write() ;
  hMixedPhot->Write() ;
  
  output.Write();
  output.Close();

}

//____________________________________________________________________________
 void AliPHOSAnalyze::AnalyzeResolutions(Int_t Nevents )    
{
  // analyzes Nevents events and calculate Energy and Position resolution as well as
  // probaility of correct indentifiing of the incident particle

  //========== Booking Histograms
  cout << "AnalyzeResolutions > " << "Booking Histograms" << endl ; 
  BookResolutionHistograms();

  Int_t Counter[9][5] ;     
  Int_t i1,i2,TotalInd = 0 ;
  for(i1 = 0; i1<9; i1++)
    for(i2 = 0; i2<5; i2++)
      Counter[i1][i2] = 0 ;
  
  Int_t TotalPrimary = 0 ;
  Int_t TotalRecPart = 0 ;
  Int_t TotalRPwithPrim = 0 ;
  Int_t ievent;

  cout << "Start Analysing"<< endl ;
  for ( ievent=0; ievent<Nevents; ievent++)
    {  
      
      //========== Event Number>         
      //      if ( ( log10((Float_t)(ievent+1)) - (Int_t)(log10((Float_t)(ievent+1))) ) == 0. ) 
	cout <<  "AnalyzeResolutions > " << "Event is " << ievent << endl ;  
      
      //=========== Connects the various Tree's for evt
      gAlice->GetEvent(ievent);

      //=========== Gets the Kine TTree
      gAlice->TreeK()->GetEvent(0) ;
      
      //=========== Gets the list of Primari Particles
      TClonesArray * PrimaryList  = gAlice->Particles();     

      TParticle * Primary ;
      Int_t iPrimary ;
      for ( iPrimary = 0 ; iPrimary < PrimaryList->GetEntries() ; iPrimary++)
	{
	  Primary = (TParticle*)PrimaryList->UncheckedAt(iPrimary) ;
	  Int_t PrimaryType = Primary->GetPdgCode() ;
	  if( PrimaryType == 22 ) {
	    Int_t ModuleNumber ;
	    Double_t PrimX, PrimZ ;
	    fGeom->ImpactOnEmc(Primary->Theta(), Primary->Phi(), ModuleNumber, PrimX, PrimZ) ;
	    if(ModuleNumber){
	      fhPrimary->Fill(Primary->Energy()) ;
	      if(Primary->Energy() > 0.3)
		TotalPrimary++ ;
	    }
	  } 
	}
      
      //=========== Get the Digit Tree
      gAlice->TreeD()->GetEvent(0) ;
      
      //========== Creating branches ===================================       
      AliPHOSRecPoint::RecPointsList ** EmcRecPoints =  fPHOS->EmcRecPoints() ;
      gAlice->TreeR()->SetBranchAddress( "PHOSEmcRP", EmcRecPoints ) ;
      
      AliPHOSRecPoint::RecPointsList ** PpsdRecPoints = fPHOS->PpsdRecPoints() ;
      gAlice->TreeR()->SetBranchAddress( "PHOSPpsdRP", PpsdRecPoints ) ;
      
      AliPHOSTrackSegment::TrackSegmentsList **  TrackSegmentsList = fPHOS->TrackSegments() ;
      if( (*TrackSegmentsList) )
	(*TrackSegmentsList)->Clear() ;
      gAlice->TreeR()->SetBranchAddress( "PHOSTS", TrackSegmentsList ) ;
      
      AliPHOSRecParticle::RecParticlesList ** RecParticleList  = fPHOS->RecParticles() ;
      if( (*RecParticleList) )
	(*RecParticleList)->Clear() ;
      gAlice->TreeR()->SetBranchAddress( "PHOSRP", RecParticleList ) ;
      
      //=========== Gets the Reconstraction TTree
      gAlice->TreeR()->GetEvent(0) ;
      
      AliPHOSRecParticle * RecParticle ;
      Int_t iRecParticle ;
      for(iRecParticle = 0; iRecParticle < (*RecParticleList)->GetEntries() ;iRecParticle++ )
	{
	  RecParticle = (AliPHOSRecParticle *) (*RecParticleList)->At(iRecParticle) ;
	  fhAllRP->Fill(CorrectEnergy(RecParticle->Energy())) ;
	  
	  Int_t ModuleNumberRec ;
	  Double_t RecX, RecZ ;
	  fGeom->ImpactOnEmc(RecParticle->Theta(), RecParticle->Phi(), ModuleNumberRec, RecX, RecZ) ;
	  
	  Double_t MinDistance = 2. ;
	  Int_t ClosestPrimary = -1 ;
	  
	  Int_t numberofprimaries ;
	  Int_t * listofprimaries  = RecParticle->GetPrimaries(numberofprimaries)  ;
	  Int_t index ;
	  TParticle * Primary ;
	  Double_t Distance = MinDistance ;
	  for ( index = 0 ; index < numberofprimaries ; index++){
	    Primary = (TParticle*)PrimaryList->UncheckedAt(listofprimaries[index]) ;
	    Int_t ModuleNumber ;
	    Double_t PrimX, PrimZ ;
	    fGeom->ImpactOnEmc(Primary->Theta(), Primary->Phi(), ModuleNumber, PrimX, PrimZ) ;
	    if(ModuleNumberRec == ModuleNumber)
	      Distance = TMath::Sqrt((RecX-PrimX)*(RecX-PrimX)+(RecZ-PrimZ)*(RecZ-PrimZ) ) ;
	    if(MinDistance > Distance)
	      {
		MinDistance = Distance ;
		ClosestPrimary = listofprimaries[index] ;
	      }
	  }
	  TotalRecPart++ ;

	  if(ClosestPrimary >=0 ){
	    TotalRPwithPrim++;
	    
	    Int_t PrimaryType = ((TParticle *)PrimaryList->At(ClosestPrimary))->GetPdgCode() ;
//  	    TParticlePDG* PDGparticle = ((TParticle *)PrimaryList->At(ClosestPrimary))->GetPDG();
//  	    Double_t charge =  PDGparticle->Charge() ;
// 	    if(charge)
// 	      cout <<"Primary " <<PrimaryType << " E " << ((TParticle *)PrimaryList->At(ClosestPrimary))->Energy() << endl ;
	    Int_t PrimaryCode ;
	    switch(PrimaryType)
	      {
	      case 22:
		PrimaryCode = 0;  //Photon
		fhAllEnergy->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(), RecParticle->Energy()) ;
		fhAllPosition->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(),MinDistance) ;
		break;
	      case 11 :
		PrimaryCode = 1;  //Electron
		break;
	      case -11 :
		PrimaryCode = 1;  //positron
		break;
	      case 321 :
		PrimaryCode = 4;  //K+
		break;
	      case -321 :
		PrimaryCode = 4;  //K-
		break;
	      case 310 :
		PrimaryCode = 4;  //K0s
		break;
	      case 130 :
		PrimaryCode = 4;  //K0l
		break;
	      case 211 :
		PrimaryCode = 2;  //K0l
		break;
	      case -211 :
		PrimaryCode = 2;  //K0l
		break;
	      case 2212 :
		PrimaryCode = 2;  //K0l
		break;
	      case -2212 :
		PrimaryCode = 2;  //K0l
		break;
	      default:
		PrimaryCode = 3; //ELSE
		break;
	      }
	    
	    switch(RecParticle->GetType())
	      {
	      case AliPHOSFastRecParticle::kGAMMA:
		if(PrimaryType == 22){
		  fhPhotEnergy->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(), RecParticle->Energy() ) ; 
		  fhEMEnergy->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(), RecParticle->Energy() ) ; 
		  fhPPSDEnergy->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(), RecParticle->Energy() ) ; 

		  fhPhotPosition->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(),MinDistance) ;
		  fhEMPosition->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(),MinDistance) ;
		  fhPPSDPosition->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(),MinDistance) ;

		  fhPhotReg->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhPhotEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhPhotPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;

		  fhPhotPhot->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		}
		if(PrimaryType == 2112){ //neutron
		  fhNReg->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhNEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhNPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		}
		
		if(PrimaryType == -2112){ //neutron ~
		  fhNBarReg->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhNBarEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhNBarPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  
		}
		if(PrimaryCode == 2){
		  fhChargedReg->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhChargedEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhChargedPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		}
		
		fhAllReg->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		fhAllEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		fhAllPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		fhShape->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		fhVeto->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		fhPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		Counter[0][PrimaryCode]++;
		break;
	      case  AliPHOSFastRecParticle::kELECTRON:
		if(PrimaryType == 22){ 
		  fhPhotElec->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhEMEnergy->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(), RecParticle->Energy() ) ; 
		  fhEMPosition->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(),MinDistance) ;
		  fhPhotEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhPhotPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		}	  
		if(PrimaryType == 2112){ //neutron
		  fhNEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhNPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		}
		
		if(PrimaryType == -2112){ //neutron ~
		  fhNBarEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhNBarPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  
		}
		if(PrimaryCode == 2){
		  fhChargedEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhChargedPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		}
		
		fhAllEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		fhAllPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		fhShape->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		fhPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		Counter[1][PrimaryCode]++;
		break;
	      case  AliPHOSFastRecParticle::kNEUTRALHA:
		if(PrimaryType == 22) 
		  fhPhotNeuH->Fill(CorrectEnergy(RecParticle->Energy()) ) ;

		fhVeto->Fill(CorrectEnergy(RecParticle->Energy()) ) ;		
		Counter[2][PrimaryCode]++;
		break ;
	      case  AliPHOSFastRecParticle::kNEUTRALEM:
		if(PrimaryType == 22){
		  fhEMEnergy->Fill(((TParticle *)PrimaryList->At(ClosestPrimary))->Energy(),RecParticle->Energy() ) ; 
		  fhEMPosition->Fill(((TParticle *)PrimaryList->At(ClosestPrimary))->Energy(),MinDistance ) ;
		
		  fhPhotNuEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhPhotEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		}
		if(PrimaryType == 2112) //neutron
		  fhNEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		
		if(PrimaryType == -2112) //neutron ~
		  fhNBarEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		
		if(PrimaryCode == 2)
		  fhChargedEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		
		fhAllEM->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		fhShape->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		fhVeto->Fill(CorrectEnergy(RecParticle->Energy()) ) ;

		Counter[3][PrimaryCode]++;
		break ;
	      case  AliPHOSFastRecParticle::kCHARGEDHA:
		if(PrimaryType == 22) //photon
		  fhPhotChHa->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		
		Counter[4][PrimaryCode]++ ;
		break ;
	      case  AliPHOSFastRecParticle::kGAMMAHA:
		  if(PrimaryType == 22){ //photon
		    fhPhotGaHa->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		    fhPPSDEnergy->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(), RecParticle->Energy() ) ; 
		    fhPPSDPosition->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(),MinDistance) ;
		    fhPhotPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  }
		  if(PrimaryType == 2112){ //neutron
		    fhNPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  }
		
		  if(PrimaryType == -2112){ //neutron ~
		    fhNBarPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ; 
		  }
		  if(PrimaryCode == 2){
		    fhChargedPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  }
		
		  fhAllPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhVeto->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  fhPPSD->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		  Counter[5][PrimaryCode]++ ;
		  break ;	
	      case  AliPHOSFastRecParticle::kABSURDEM:	      
		Counter[6][PrimaryCode]++ ;
		fhShape->Fill(CorrectEnergy(RecParticle->Energy()) ) ;
		break;
	      case  AliPHOSFastRecParticle::kABSURDHA:
		Counter[7][PrimaryCode]++ ;
		break;
	      default:
		Counter[8][PrimaryCode]++ ;
		break;
	      }
	  }
	}  
    }   // endfor
  SaveHistograms();
  cout << "Resolutions: Analyzed " << Nevents << " event(s)" << endl ;
  cout << "Resolutions: Total primary       " << TotalPrimary << endl ;
  cout << "Resoluitons: Total reconstracted " << TotalRecPart << endl ;
  cout << "TotalReconstructed with Primarie " << TotalRPwithPrim << endl ;
  cout << "                        Primary:   Photon   Electron   Ch. Hadr.  Neutr. Hadr  Kaons" << endl ; 
  cout << "             Detected as photon       " << Counter[0][0] << "          " << Counter[0][1] << "          " << Counter[0][2] << "          " <<Counter[0][3] << "          " << Counter[0][4] << endl ;
  cout << "           Detected as electron       " << Counter[1][0] << "          " << Counter[1][1] << "          " << Counter[1][2] << "          " <<Counter[1][3] << "          " << Counter[1][4] << endl ; 
  cout << "     Detected as neutral hadron       " << Counter[2][0] << "          " << Counter[2][1] << "          " << Counter[2][2] << "          " <<Counter[2][3] << "          " << Counter[2][4] << endl ;
  cout << "         Detected as neutral EM       " << Counter[3][0] << "          " << Counter[3][1] << "          " << Counter[3][2] << "          " <<Counter[3][3] << "          " << Counter[3][4] << endl ;
  cout << "     Detected as charged hadron       " << Counter[4][0] << "          " << Counter[4][1] << "          " << Counter[4][2] << "          " <<Counter[4][3] << "          " << Counter[4][4] << endl ;
  cout << "       Detected as gamma-hadron       " << Counter[5][0] << "          " << Counter[5][1] << "          " << Counter[5][2] << "          " <<Counter[5][3] << "          " << Counter[5][4] << endl ;
  cout << "          Detected as Absurd EM       " << Counter[6][0] << "          " << Counter[6][1] << "          " << Counter[6][2] << "          " <<Counter[6][3] << "          " << Counter[6][4] << endl ;
  cout << "      Detected as absurd hadron       " << Counter[7][0] << "          " << Counter[7][1] << "          " << Counter[7][2] << "          " <<Counter[7][3] << "          " << Counter[7][4] << endl ;
  cout << "          Detected as undefined       " << Counter[8][0] << "          " << Counter[8][1] << "          " << Counter[8][2] << "          " <<Counter[8][3] << "          " << Counter[8][4] << endl ;
      
      for(i1 = 0; i1<9; i1++)
	for(i2 = 0; i2<5; i2++)
	  TotalInd+=Counter[i1][i2] ;
      cout << "Indentified particles            " << TotalInd << endl ;
      
}           // endfunction


//____________________________________________________________________________
void  AliPHOSAnalyze::BookingHistograms()
{
  // Books the histograms where the results of the analysis are stored (to be changed)

  delete fhEmcDigit  ;
  delete fhVetoDigit  ;
  delete fhConvertorDigit   ;
  delete  fhEmcCluster   ;
  delete fhVetoCluster   ;
  delete fhConvertorCluster  ;
  delete fhConvertorEmc  ;
  
  fhEmcDigit                = new TH1F("hEmcDigit",      "hEmcDigit",         1000,  0. ,  25.);
  fhVetoDigit               = new TH1F("hVetoDigit",     "hVetoDigit",         500,  0. ,  3.e-5);
  fhConvertorDigit          = new TH1F("hConvertorDigit","hConvertorDigit",    500,  0. ,  3.e-5);
  fhEmcCluster              = new TH1F("hEmcCluster",    "hEmcCluster",       1000,  0. ,  30.);
  fhVetoCluster             = new TH1F("hVetoCluster",   "hVetoCluster",       500,  0. ,  3.e-5);
  fhConvertorCluster        = new TH1F("hConvertorCluster","hConvertorCluster",500,  0. ,  3.e-5);
  fhConvertorEmc            = new TH2F("hConvertorEmc",  "hConvertorEmc",      200,  1. ,  3., 200, 0., 3.e-5);

}
//____________________________________________________________________________
void  AliPHOSAnalyze::BookResolutionHistograms()
{
  // Books the histograms where the results of the Resolution analysis are stored

//   if(fhAllEnergy)
//     delete fhAllEnergy ;
//   if(fhPhotEnergy)
//     delete fhPhotEnergy ;
//   if(fhEMEnergy)
//     delete fhEMEnergy ;
//   if(fhPPSDEnergy)
//     delete fhPPSDEnergy ;


  fhAllEnergy  = new TH2F("hAllEnergy",  "Energy of any RP with primary photon",100, 0., 5., 100, 0., 5.);
  fhPhotEnergy = new TH2F("hPhotEnergy", "Energy of kGAMMA with primary photon",100, 0., 5., 100, 0., 5.);
  fhEMEnergy   = new TH2F("hEMEnergy",   "Energy of EM with primary photon",    100, 0., 5., 100, 0., 5.);
  fhPPSDEnergy = new TH2F("hPPSDEnergy", "Energy of PPSD with primary photon",  100, 0., 5., 100, 0., 5.);

//   if(fhAllPosition)
//     delete fhAllPosition ;
//   if(fhPhotPosition)
//     delete fhPhotPosition ;
//   if(fhEMPosition)
//     delete fhEMPosition ;
//   if(fhPPSDPosition)
//     delete fhPPSDPosition ;


  fhAllPosition  = new TH2F("hAllPosition",  "Position of any RP with primary photon",100, 0., 5., 100, 0., 5.);
  fhPhotPosition = new TH2F("hPhotPosition", "Position of kGAMMA with primary photon",100, 0., 5., 100, 0., 5.);
  fhEMPosition   = new TH2F("hEMPosition",   "Position of EM with primary photon",    100, 0., 5., 100, 0., 5.);
  fhPPSDPosition = new TH2F("hPPSDPosition", "Position of PPSD with primary photon",  100, 0., 5., 100, 0., 5.);

//   if(fhAllReg)
//     delete fhAllReg ;
//   if(fhPhotReg)
//     delete fhPhotReg ;
//   if(fhNReg)
//     delete fhNReg ;
//   if(fhNBarReg)
//     delete fhNBarReg ;
//   if(fhChargedReg)
//     delete fhChargedReg ;
  
  fhAllReg    = new TH1F("hAllReg",    "All primaries registered as photon",  100, 0., 5.);
  fhPhotReg   = new TH1F("hPhotReg",   "Photon registered as photon",         100, 0., 5.);
  fhNReg      = new TH1F("hNReg",      "N registered as photon",              100, 0., 5.);
  fhNBarReg   = new TH1F("hNBarReg",   "NBar registered as photon",           100, 0., 5.);
  fhChargedReg= new TH1F("hChargedReg", "Charged hadron registered as photon",100, 0., 5.);
  
//   if(fhAllEM)
//     delete fhAllEM ;
//   if(fhPhotEM)
//     delete fhPhotEM ;
//   if(fhNEM)
//     delete fhNEM ;
//   if(fhNBarEM)
//     delete fhNBarEM ;
//   if(fhChargedEM)
//     delete fhChargedEM ;
  
  fhAllEM    = new TH1F("hAllEM",    "All primary registered as EM",100, 0., 5.);
  fhPhotEM   = new TH1F("hPhotEM",   "Photon registered as EM", 100, 0., 5.);
  fhNEM      = new TH1F("hNEM",      "N registered as EM",      100, 0., 5.);
  fhNBarEM   = new TH1F("hNBarEM",   "NBar registered as EM",   100, 0., 5.);
  fhChargedEM= new TH1F("hChargedEM","Charged registered as EM",100, 0., 5.);

//   if(fhAllPPSD)
//     delete fhAllPPSD ;
//   if(fhPhotPPSD)
//     delete fhPhotPPSD ;
//   if(fhNPPSD)
//     delete fhNPPSD ;
//   if(fhNBarPPSD)
//     delete fhNBarPPSD ;
//   if(fhChargedPPSD)
//     delete fhChargedPPSD ;
  
  fhAllPPSD    = new TH1F("hAllPPSD",    "All primary registered as PPSD",100, 0., 5.);
  fhPhotPPSD   = new TH1F("hPhotPPSD",   "Photon registered as PPSD", 100, 0., 5.);
  fhNPPSD      = new TH1F("hNPPSD",      "N registered as PPSD",      100, 0., 5.);
  fhNBarPPSD   = new TH1F("hNBarPPSD",   "NBar registered as PPSD",   100, 0., 5.);
  fhChargedPPSD= new TH1F("hChargedPPSD","Charged registered as PPSD",100, 0., 5.);
  
//   if(fhPrimary)
//     delete fhPrimary ;
  fhPrimary= new TH1F("hPrimary", "hPrimary",  100, 0., 5.);

//   if(fhAllRP)
//     delete fhAllRP ;
//   if(fhVeto)
//     delete fhVeto ;
//   if(fhShape)
//     delete fhShape ;
//   if(fhPPSD)
//     delete fhPPSD ;

  fhAllRP = new TH1F("hAllRP","All Reconstructed particles",  100, 0., 5.);
  fhVeto  = new TH1F("hVeto", "All uncharged particles",      100, 0., 5.);
  fhShape = new TH1F("hShape","All particles with EM shaower",100, 0., 5.);
  fhPPSD  = new TH1F("hPPSD", "All PPSD photon particles",    100, 0., 5.);


//   if(fhPhotPhot)
//     delete fhPhotPhot ;
//   if(fhPhotElec)
//     delete fhPhotElec ;
//   if(fhPhotNeuH)
//     delete fhPhotNeuH ;
//   if(fhPhotNuEM)
//     delete fhPhotNuEM ;
//   if(fhPhotChHa)
//     delete fhPhotChHa ;
//   if(fhPhotGaHa)
//     delete fhPhotGaHa ;

  fhPhotPhot = new TH1F("hPhotPhot","hPhotPhot", 100, 0., 5.);   //Photon registered as photon
  fhPhotElec = new TH1F("hPhotElec","hPhotElec", 100, 0., 5.);   //Photon registered as Electron
  fhPhotNeuH = new TH1F("hPhotNeuH","hPhotNeuH", 100, 0., 5.);   //Photon registered as Neutral Hadron
  fhPhotNuEM = new TH1F("hPhotNuEM","hPhotNuEM", 100, 0., 5.);   //Photon registered as Neutral EM
  fhPhotChHa = new TH1F("hPhotChHa","hPhotChHa", 100, 0., 5.);   //Photon registered as Charged Hadron
  fhPhotGaHa = new TH1F("hPhotGaHa","hPhotGaHa", 100, 0., 5.);   //Photon registered as Gamma-Hadron


}
//____________________________________________________________________________
Bool_t AliPHOSAnalyze::Init(Int_t evt)
{
  // Do a few initializations: open the root file
  //                           get the AliRun object
  //                           defines the clusterizer, tracksegment maker and particle identifier
  //                           sets the associated parameters

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
      
      fPHOS  = (AliPHOSv1 *)gAlice->GetDetector("PHOS") ;
      fGeom = fPHOS->GetGeometry();
      //      fGeom  = AliPHOSGeometry::GetInstance( fPHOS->GetGeometry()->GetName(), fPHOS->GetGeometry()->GetTitle() );

    } // else !ok
  } // if fRootFile
  
  if ( ok ) {
    
    //========== Create the Clusterizer

    fClu =  new AliPHOSClusterizerv1() ; 
    fClu->SetEmcEnergyThreshold(0.030) ; 
    fClu->SetEmcClusteringThreshold(0.20) ; 
    fClu->SetPpsdEnergyThreshold    (0.0000002) ; 
    fClu->SetPpsdClusteringThreshold(0.0000001) ; 
    fClu->SetLocalMaxCut(0.03) ;
    fClu->SetCalibrationParameters(0., 0.00000001) ;  
    cout <<  "AnalyzeOneEvent > using clusterizer " << fClu->GetName() << endl ; 
    fClu->PrintParameters() ; 
    
    //========== Creates the track segment maker
    
    fTrs = new AliPHOSTrackSegmentMakerv1() ;
    cout <<  "AnalyzeOneEvent > using tack segment maker " << fTrs->GetName() << endl ; 
    //   fTrs->UnsetUnfoldFlag() ;
    
    //========== Creates the particle identifier
    
    fPID = new AliPHOSPIDv1() ;
    cout <<  "AnalyzeOneEvent > using particle identifier " << fPID->GetName() << endl ; 
    //fPID->SetShowerProfileCuts(Float_t l1m, Float_t l1M, Float_t l2m, Float_t l2M) ; 
    fPID->SetShowerProfileCuts(0.7, 2.0 , 0.6 , 1.5) ; 

    //========== Creates the Reconstructioner  
    
    fRec = new AliPHOSReconstructioner(fClu, fTrs, fPID) ;
//      fRec -> SetDebugReconstruction(kFALSE);     
    fRec -> SetDebugReconstruction(kTRUE);     
    
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
  // Display particles from the Kine Tree in global Alice (theta, phi) coordinates. 
  // One PHOS module at the time.
  // The particle type can be selected.
  
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
  fGeom->EmcModuleCoverage(module, tm, tM, pm, pM, AliPHOSGeometry::Degre() ) ;

  Double_t theta, phi ; 
  fGeom->EmcXtalCoverage(theta, phi, AliPHOSGeometry::Degre() ) ;

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
	  if (particle->Energy() >  fClu->GetEmcClusteringThreshold()  )
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
  // Display reconstructed particles in global Alice(theta, phi) coordinates. 
  // One PHOS module at the time.
  // Click on symbols indicate the reconstructed particle type. 

  if (fEvt == -999) {
    cout << "DisplayRecParticles > Analyze an event first ... (y/n) " ; 
    Text_t answer[1] ; 
    cin >> answer ; cout << answer ; 
//     if ( answer == "y" ) 
//       AnalyzeOneEvent() ;
  } 
    if (fEvt != -999) {
      
      Int_t module ; 
      cout <<  "DisplayRecParticles > which module (1-5,  -1: all) ? " ; 
      cin >> module ; cout << module << endl ;
      Text_t histoname[80] ; 
      sprintf(histoname,"Event %d: Reconstructed particles in module %d", fEvt, module) ; 
      Double_t tm, tM, pm, pM ; // min and Max theta and phi covered by module   
      fGeom->EmcModuleCoverage(module, tm, tM, pm, pM, AliPHOSGeometry::Degre() ) ;
      Double_t theta, phi ; 
      fGeom->EmcXtalCoverage(theta, phi, AliPHOSGeometry::Degre() ) ;
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
      AliPHOSRecParticle::RecParticlesList * rpl = *fPHOS->RecParticles() ; 
      Int_t nRecParticles = rpl->GetEntries() ; 
      Int_t nRecParticlesInModule = 0 ; 
      TIter nextRecPart(rpl) ; 
      AliPHOSRecParticle * rp ; 
      cout << "DisplayRecParticles > " << nRecParticles << " reconstructed particles " << endl ; 
      Double_t kRADDEG = 180. / TMath::Pi() ; 
      while ( (rp = (AliPHOSRecParticle *)nextRecPart() ) ) {
	AliPHOSTrackSegment * ts = rp->GetPHOSTrackSegment() ; 
	if ( ts->GetPHOSMod() == module ) {
	  Int_t numberofprimaries = 0 ;
	  Int_t * listofprimaries = 0;
	  rp->GetPrimaries(numberofprimaries) ;
	  cout << "Number of primaries = " << numberofprimaries << endl ; 
	  Int_t index ;
	  for ( index = 0 ; index < numberofprimaries ; index++)
	    cout << "    primary # " << index << " =  " << listofprimaries[index] << endl ;  
	  
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
  // Display reconstructed points in local PHOS-module (x, z) coordinates. 
  // One PHOS module at the time.
  // Click on symbols displays the EMC cluster, or PPSD information.

  if (fEvt == -999) {
    cout << "DisplayRecPoints > Analyze an event first ... (y/n) " ; 
    Text_t answer[1] ; 
    cin >> answer ; cout << answer ; 
//     if ( answer == "y" ) 
//       AnalyzeOneEvent() ;
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

      //      TClonesArray * emcRP = fPHOS->EmcClusters() ; 
      TObjArray * emcRP = *(fPHOS->EmcRecPoints()) ; 
      
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

      //      TClonesArray * ppsdRP = fPHOS->PpsdClusters() ;
      TObjArray * ppsdRP = *(fPHOS->PpsdRecPoints() );
 
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
  
      ppsdRP = *(fPHOS->PpsdRecPoints()) ;
     
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
  // Display track segments in local PHOS-module (x, z) coordinates. 
  // One PHOS module at the time.
  // One symbol per PHOS subsystem: EMC, upper PPSD, lower PPSD.

  if (fEvt == -999) {
    cout << "DisplayTrackSegments > Analyze an event first ... (y/n) " ; 
    Text_t answer[1] ; 
    cin >> answer ; cout << answer ; 
//     if ( answer == "y" ) 
//       AnalyzeOneEvent() ;
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

      AliPHOSTrackSegment::TrackSegmentsList * trsegl = *(fPHOS->TrackSegments()) ;
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
  // Open the root file named "name"
  
  fRootFile   = new TFile(name, "update") ;
  return  fRootFile->IsOpen() ; 
}
//____________________________________________________________________________
void AliPHOSAnalyze::SaveHistograms()
{
  // Saves the histograms in a root file named "name.analyzed" 

  Text_t outputname[80] ;
  sprintf(outputname,"%s.analyzed",fRootFile->GetName());
  TFile output(outputname,"RECREATE");
  output.cd();

  if (fhAllEnergy)    
    fhAllEnergy->Write() ;
  if (fhPhotEnergy)    
    fhPhotEnergy->Write() ;
  if(fhEMEnergy)
    fhEMEnergy->Write()  ;
  if(fhPPSDEnergy)
    fhPPSDEnergy->Write() ;
  if(fhAllPosition)
    fhAllPosition->Write() ;
  if(fhPhotPosition)
    fhPhotPosition->Write() ;
  if(fhEMPosition)
    fhEMPosition->Write() ;
  if(fhPPSDPosition)
    fhPPSDPosition->Write() ;
  if (fhAllReg) 
    fhAllReg->Write() ;
  if (fhPhotReg) 
    fhPhotReg->Write() ;
  if(fhNReg)
    fhNReg->Write() ;
  if(fhNBarReg)
    fhNBarReg->Write() ;
  if(fhChargedReg)
    fhChargedReg->Write() ;
  if (fhAllEM) 
    fhAllEM->Write() ;
  if (fhPhotEM) 
    fhPhotEM->Write() ;
  if(fhNEM)
    fhNEM->Write() ;
  if(fhNBarEM)
    fhNBarEM->Write() ;
  if(fhChargedEM)
    fhChargedEM->Write() ;
  if (fhAllPPSD) 
    fhAllPPSD->Write() ;
  if (fhPhotPPSD) 
    fhPhotPPSD->Write() ;
  if(fhNPPSD)
    fhNPPSD->Write() ;
  if(fhNBarPPSD)
    fhNBarPPSD->Write() ;
  if(fhChargedPPSD)
    fhChargedPPSD->Write() ;
  if(fhPrimary)
    fhPrimary->Write() ;
  if(fhAllRP)
    fhAllRP->Write()  ;
  if(fhVeto)
    fhVeto->Write()  ;
  if(fhShape)
    fhShape->Write()  ;
  if(fhPPSD)
    fhPPSD->Write()  ;
  if(fhPhotPhot)
    fhPhotPhot->Write() ;
  if(fhPhotElec)
    fhPhotElec->Write() ;
  if(fhPhotNeuH)
    fhPhotNeuH->Write() ;
  if(fhPhotNuEM)
    fhPhotNuEM->Write() ;
  if(fhPhotNuEM)
    fhPhotNuEM->Write() ;
  if(fhPhotChHa)
    fhPhotChHa->Write() ;
  if(fhPhotGaHa)
    fhPhotGaHa->Write() ;
  if(fhEnergyCorrelations)
    fhEnergyCorrelations->Write() ;
  
  output.Write();
  output.Close();
}
//____________________________________________________________________________
Float_t AliPHOSAnalyze::CorrectEnergy(Float_t ERecPart)
{
  return ERecPart/0.8783 ;
}

//____________________________________________________________________________
void AliPHOSAnalyze::ResetHistograms()
{
   fhEnergyCorrelations = 0 ;     //Energy correlations between Eloss in Convertor and PPSD(2)

   fhEmcDigit = 0 ;               // Histo of digit energies in the Emc 
   fhVetoDigit = 0 ;              // Histo of digit energies in the Veto 
   fhConvertorDigit = 0 ;         // Histo of digit energies in the Convertor
   fhEmcCluster = 0 ;             // Histo of Cluster energies in Emc
   fhVetoCluster = 0 ;            // Histo of Cluster energies in Veto
   fhConvertorCluster = 0 ;       // Histo of Cluster energies in Convertor
   fhConvertorEmc = 0 ;           // 2d Convertor versus Emc energies

   fhAllEnergy = 0 ;       
   fhPhotEnergy = 0 ;        // Total spectrum of detected photons
   fhEMEnergy = 0 ;         // Spectrum of detected electrons with electron primary
   fhPPSDEnergy = 0 ;
   fhAllPosition = 0 ; 
   fhPhotPosition = 0 ; 
   fhEMPosition = 0 ; 
   fhPPSDPosition = 0 ; 

   fhPhotReg = 0 ;          
   fhAllReg = 0 ;          
   fhNReg = 0 ;          
   fhNBarReg = 0 ;          
   fhChargedReg = 0 ;          
   fhPhotEM = 0 ;          
   fhAllEM = 0 ;          
   fhNEM = 0 ;          
   fhNBarEM = 0 ;          
   fhChargedEM = 0 ;          
   fhPhotPPSD = 0 ;          
   fhAllPPSD = 0 ;          
   fhNPPSD = 0 ;          
   fhNBarPPSD = 0 ;          
   fhChargedPPSD = 0 ;          

   fhPrimary = 0 ;          

   fhPhotPhot = 0 ;
   fhPhotElec = 0 ;
   fhPhotNeuH = 0 ;
   fhPhotNuEM = 0 ; 
   fhPhotChHa = 0 ;
   fhPhotGaHa = 0 ;


}
