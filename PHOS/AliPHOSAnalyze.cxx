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
#include "AliPHOSHit.h"
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
void AliPHOSAnalyze::DrawRecon(Int_t Nevent,Int_t Nmod){
  //Draws pimary particles and reconstructed 
  //digits, RecPoints, RecPartices etc 
  //for event Nevent in the module Nmod.

  TH2F * digitOccupancy  = new TH2F("digitOccupancy","EMC digits", 64,-71.,71.,64,-71.,71.);
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

  //========== Create the Clusterizer
  fClu = new AliPHOSClusterizerv1() ; 

  fClu->SetEmcEnergyThreshold(0.05) ;
  fClu->SetEmcClusteringThreshold(0.20) ;
  fClu->SetPpsdEnergyThreshold    (0.0000002) ;
  fClu->SetPpsdClusteringThreshold(0.0000001) ;
  fClu->SetLocalMaxCut(0.03) ;
  fClu->SetCalibrationParameters(0., 0.00000001) ;
  
  gAlice->GetEvent(Nevent);
  
  TClonesArray * primaryList  = gAlice->Particles();
  
  TParticle * primary ;
  Int_t iPrimary ;
  for ( iPrimary = 0 ; iPrimary < primaryList->GetEntries() ; iPrimary++)
    {
      primary = (TParticle*)primaryList->At(iPrimary) ;
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

  fPHOS->SetTreeAddress() ;

  gAlice->TreeD()->GetEvent(0) ;
  gAlice->TreeR()->GetEvent(0) ;
  
  TObjArray ** emcRecPoints =  fPHOS->EmcRecPoints() ;
  TObjArray ** ppsdRecPoints = fPHOS->PpsdRecPoints() ;
  TClonesArray ** recParticleList  = fPHOS->RecParticles() ;
  
  Int_t iDigit ;
  AliPHOSDigit * digit ;
  
  for(iDigit = 0; iDigit < fPHOS->Digits()->GetEntries(); iDigit++)
    {
      digit = (AliPHOSDigit *) fPHOS->Digits()->At(iDigit) ;
      Int_t relid[4];
      fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
      Float_t x,z ;
      fGeom->RelPosInModule(relid,x,z) ;
      Float_t e = fClu->Calibrate(digit->GetAmp()) ;
      if(relid[0]==Nmod){
        if(relid[1]==0)  //EMC
          digitOccupancy->Fill(x,z,e) ;
        if((relid[1]>0)&&(relid[1]<17))
          ppsdUp->Fill(x,z,e) ;
        if(relid[1]>16)
          ppsdLow->Fill(x,z,e) ;
      }
    }
  
  Int_t irecp ;
  TVector3 pos ;
  
  for(irecp = 0; irecp < (*emcRecPoints)->GetEntries() ; irecp ++){
    AliPHOSEmcRecPoint * emc= (AliPHOSEmcRecPoint*)(*emcRecPoints)->At(irecp) ;
    if(emc->GetPHOSMod()==Nmod){
      emc->GetLocalPosition(pos) ;
      emcOccupancy->Fill(pos.X(),pos.Z(),emc->GetEnergy());
    }
  }
  
  for(irecp = 0; irecp < (*ppsdRecPoints)->GetEntries() ; irecp ++){
    AliPHOSPpsdRecPoint * ppsd= (AliPHOSPpsdRecPoint *)(*ppsdRecPoints)->At(irecp) ;
    if(ppsd->GetPHOSMod()==Nmod){
      ppsd->GetLocalPosition(pos) ;
      if(ppsd->GetUp())
        ppsdUpCl->Fill(pos.X(),pos.Z(),ppsd->GetEnergy());
      else
        ppsdLowCl->Fill(pos.X(),pos.Z(),ppsd->GetEnergy());
    }
  }
  
  AliPHOSRecParticle * recParticle ;
  Int_t iRecParticle ;
  for(iRecParticle = 0; iRecParticle < (*recParticleList)->GetEntries() ;iRecParticle++ )
    {
      recParticle = (AliPHOSRecParticle *) (*recParticleList)->At(iRecParticle) ;
      
      Int_t moduleNumberRec ;
      Double_t recX, recZ ;
      fGeom->ImpactOnEmc(recParticle->Theta(), recParticle->Phi(), moduleNumberRec, recX, recZ) ;
      if(moduleNumberRec == Nmod){
	
        Double_t minDistance = 5. ;
        Int_t closestPrimary = -1 ;
	
        Int_t numberofprimaries ;
        Int_t * listofprimaries  = recParticle->GetPrimaries(numberofprimaries)  ;
        Int_t index ;
        TParticle * primary ;
        Double_t distance = minDistance ;
	  
	for ( index = 0 ; index < numberofprimaries ; index++){
	  primary = (TParticle*)primaryList->At(listofprimaries[index]) ;
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
	  
          Int_t primaryType = ((TParticle *)primaryList->At(closestPrimary))->GetPdgCode() ;
	  
          if(primaryType==22)
            recPhot->Fill(recZ,recX,recParticle->Energy()) ;
          else
            if(primaryType==-2112)
              recNbar->Fill(recZ,recX,recParticle->Energy()) ; 
        }
      }
    }
  
  
  digitOccupancy->Draw("box") ;
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
 void AliPHOSAnalyze::Reconstruct(Int_t nevents,Int_t firstEvent )    
{     

  // Performs reconstruction of EMC and CPV (GPS2, IHEP or MIXT)
  // for events from FirstEvent to Nevents

  Int_t ievent ;   
  for ( ievent=firstEvent; ievent<nevents; ievent++) {  
    if (ievent==firstEvent) {
      cout << "Analyze > Starting Reconstructing " << endl ; 
      //========== Create the Clusterizer
      fClu = new AliPHOSClusterizerv1() ; 
      fClu->SetEmcEnergyThreshold(0.05) ; 
      fClu->SetEmcClusteringThreshold(0.20) ; 
      fClu->SetLocalMaxCut(0.03) ;
      if      (strcmp(fGeom->GetName(),"GPS2") == 0 || strcmp(fGeom->GetName(),"MIXT") == 0) {
	fClu->SetPpsdEnergyThreshold    (0.0000002) ; 
	fClu->SetPpsdClusteringThreshold(0.0000001) ; 
      }
      else if (strcmp(fGeom->GetName(),"IHEP") == 0 || strcmp(fGeom->GetName(),"MIXT") == 0) {
	fClu->SetLocalMaxCutCPV(0.03) ;
	fClu->SetLogWeightCutCPV(4.0) ;
	fClu->SetCpvEnergyThreshold(0.2) ;
      }
      fClu->SetCalibrationParameters(0., 0.00000001) ; 
      
      //========== Creates the track segment maker
      fTrs = new AliPHOSTrackSegmentMakerv1()  ;
      //	  fTrs->UnsetUnfoldFlag() ; 
     
      //========== Creates the particle identifier
      fPID = new AliPHOSPIDv1() ;
      fPID->SetShowerProfileCuts(0.3, 1.8, 0.3, 1.8 ) ;       
      
      //========== Creates the Reconstructioner
      fRec = new AliPHOSReconstructioner(fClu, fTrs, fPID) ; 
      if (fDebugLevel != 0) fRec -> SetDebugReconstruction(kTRUE);     
    }
    
    if (fDebugLevel != 0 ||
	(ievent+1) % (Int_t)TMath::Power( 10, (Int_t)TMath::Log10(ievent+1) ) == 0)
      cout <<  "======= Analyze ======> Event " << ievent+1 << endl ;
    
    //=========== Connects the various Tree's for evt
    Int_t tracks = gAlice->GetEvent(ievent);
    
    fPHOS->Hit2Digit(tracks) ;
    
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
//   //
//   // Read and print generated and reconstructed hits in CPV
//   // for events from EvFirst to Nevent.
//   // If only EvFirst is defined, print only this one event.
//   // Author: Yuri Kharlov
//   // 12 October 2000
//   //

//   if (EvFirst!=0 && EvLast==0) EvLast=EvFirst;
//   for ( Int_t ievent=EvFirst; ievent<=EvLast; ievent++) {  
    
//     //========== Event Number>
//     cout << endl <<  "==== ReadAndPrintCPV ====> Event is " << ievent+1 << endl ;
    
//     //=========== Connects the various Tree's for evt
//     Int_t ntracks = gAlice->GetEvent(ievent);

//     //========== Creating branches ===================================
//     AliPHOSRecPoint::RecPointsList ** emcRecPoints = fPHOS->EmcRecPoints() ;
//     gAlice->TreeR()->SetBranchAddress( "PHOSEmcRP" , emcRecPoints  ) ;
    
//     AliPHOSRecPoint::RecPointsList ** cpvRecPoints = fPHOS->PpsdRecPoints() ;
//     gAlice->TreeR()->SetBranchAddress( "PHOSPpsdRP", cpvRecPoints ) ;

//     // Read and print CPV hits
      
//     AliPHOSCPVModule cpvModule;
//     TClonesArray    *cpvHits;
//     Int_t           nCPVhits;
//     AliPHOSCPVHit   *cpvHit;
//     TLorentzVector   p;
//     Float_t          xgen, zgen;
//     Int_t            ipart;
//     Int_t            nGenHits = 0;
//     for (Int_t itrack=0; itrack<ntracks; itrack++) {
//       //=========== Get the Hits Tree for the Primary track itrack
//       gAlice->ResetHits();
//       gAlice->TreeH()->GetEvent(itrack);
//       Int_t iModule = 0 ; 	
//       for (iModule=0; iModule < fGeom->GetNCPVModules(); iModule++) {
// 	cpvModule = fPHOS->GetCPVModule(iModule);
// 	cpvHits   = cpvModule.Hits();
// 	nCPVhits  = cpvHits->GetEntriesFast();
// 	for (Int_t ihit=0; ihit<nCPVhits; ihit++) {
// 	  nGenHits++;
// 	  cpvHit = (AliPHOSCPVHit*)cpvHits->UncheckedAt(ihit);
// 	  p      = cpvHit->GetMomentum();
// 	  xgen   = cpvHit->X();
// 	  zgen   = cpvHit->Y();
// 	  ipart  = cpvHit->GetIpart();
// 	  printf("CPV hit in module %d: ",iModule+1);
// 	  printf(" p = (%f, %f, %f, %f) GeV,\n",
// 		 p.Px(),p.Py(),p.Pz(),p.Energy());
// 	  printf("                  (X,Z) = (%8.4f, %8.4f) cm, ipart = %d\n",
// 		 xgen,zgen,ipart);
// 	}
//       }
//     }

//     // Read and print CPV reconstructed points

//     //=========== Gets the Reconstruction TTree
//     gAlice->TreeR()->GetEvent(0) ;
//     printf("Recpoints: %d\n",(*fPHOS->CpvRecPoints())->GetEntries());
//     TIter nextRP(*fPHOS->CpvRecPoints() ) ;
//     AliPHOSCpvRecPoint *cpvRecPoint ;
//     Int_t nRecPoints = 0;
//     while( ( cpvRecPoint = (AliPHOSCpvRecPoint *)nextRP() ) ) {
//       nRecPoints++;
//       TVector3  locpos;
//       cpvRecPoint->GetLocalPosition(locpos);
//       Int_t phosModule = cpvRecPoint->GetPHOSMod();
//       printf("CPV recpoint in module %d: (X,Z) = (%f,%f) cm\n",
// 	     phosModule,locpos.X(),locpos.Z());
//     }
//     printf("This event has %d generated hits and %d reconstructed points\n",
// 	   nGenHits,nRecPoints);
//   }
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
    TClonesArray **hitsPerModule = new TClonesArray *[nOfModules];
    Int_t iModule = 0; 	
    for (iModule=0; iModule < nOfModules; iModule++)
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
	  xzgen[0] = cpvHit->X();
	  xzgen[1] = cpvHit->Y();
	  ipart    = cpvHit->GetIpart();
	  TClonesArray &lhits = *(TClonesArray *)hitsPerModule[iModule];
	  new(lhits[hitsPerModule[iModule]->GetEntriesFast()]) AliPHOSCPVHit(*cpvHit);
	}
	cpvModule.Clear();
      }
    }
    for (iModule=0; iModule < nOfModules; iModule++) {
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
	xgen   = cpvHit->X();
	zgen   = cpvHit->Y();
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
    delete [] hitsPerModule;
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

  const Int_t nMixedEvents = 4 ; //# of events used for calculation of 'mixed' distribution 
  Int_t mixedLoops = (Int_t )TMath::Ceil(Nevents/nMixedEvents) ;
  
  //========== Booking Histograms
  TH2D * hRealEM   = new TH2D("hRealEM",   "Real for EM particles",      250,0.,1.,40,0.,4.) ;
  TH2D * hRealPhot = new TH2D("hRealPhot", "Real for kPhoton particles", 250,0.,1.,40,0.,4.) ;
  TH2D * hMixedEM  = new TH2D("hMixedEM",  "Mixed for EM particles",     250,0.,1.,40,0.,4.) ;
  TH2D * hMixedPhot= new TH2D("hMixedPhot","Mixed for kPhoton particles",250,0.,1.,40,0.,4.) ;
  
  Int_t ievent;
  Int_t eventInMixedLoop ;
  
  Int_t nRecParticles[4];//nMixedEvents] ;
  
  AliPHOSRecParticle::RecParticlesList * allRecParticleList  = new TClonesArray("AliPHOSRecParticle", nMixedEvents*1000) ;
  
  for(eventInMixedLoop = 0; eventInMixedLoop < mixedLoops; eventInMixedLoop++  ){
    Int_t iRecPhot = 0 ;
    
    for ( ievent=0; ievent < nMixedEvents; ievent++){        
      
      Int_t absEventNumber = eventInMixedLoop*nMixedEvents + ievent ;
      
      //=========== Connects the various Tree's for evt
      gAlice->GetEvent(absEventNumber);

      //========== Creating branches ===================================       
      fPHOS->SetTreeAddress() ;
      
      gAlice->TreeD()->GetEvent(0) ;
      gAlice->TreeR()->GetEvent(0) ;
      
      TClonesArray ** recParticleList  = fPHOS->RecParticles() ;
      
            
      AliPHOSRecParticle * recParticle ;
      Int_t iRecParticle ;
      for(iRecParticle = 0; iRecParticle < (*recParticleList)->GetEntries() ;iRecParticle++ )
 	{
 	  recParticle = (AliPHOSRecParticle *) (*recParticleList)->At(iRecParticle) ;
	  if((recParticle->GetType() == AliPHOSFastRecParticle::kGAMMA)||
	     (recParticle->GetType() == AliPHOSFastRecParticle::kNEUTRALEM)){ 
	    new( (*allRecParticleList)[iRecPhot] ) AliPHOSRecParticle(*recParticle) ;
	    iRecPhot++;
	  }
	}
      
	nRecParticles[ievent] = iRecPhot-1 ;  
    }
    
    //Now calculate invariant mass:
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
	pt = TMath::Sqrt((rp1->Px()+rp2->Px() )*( rp1->Px()+rp2->Px() ) +(rp1->Py()+rp2->Py())*(rp1->Py()+rp2->Py()));

	if(irp1 > nRecParticles[nCurEvent])
	  nCurEvent++;
	    
	if(irp2 <= nRecParticles[nCurEvent]){ //'Real' event
	  hRealEM->Fill(invMass,pt);
	  if((rp1->GetType() == AliPHOSFastRecParticle::kGAMMA)&&(rp2->GetType() == AliPHOSFastRecParticle::kGAMMA))
	    hRealPhot->Fill(invMass,pt);
	}
	else{
	  hMixedEM->Fill(invMass,pt);
	  if((rp1->GetType() == AliPHOSFastRecParticle::kGAMMA)&&(rp2->GetType() == AliPHOSFastRecParticle::kGAMMA))
	    hMixedPhot->Fill(invMass,pt);
	} //real-mixed
	    
      } //loop over second rp
    }//loop over first rp
    allRecParticleList->Delete() ;
  } //Loop over events
  
  delete allRecParticleList ;
  
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
 void AliPHOSAnalyze::ReadAndPrintEMC(Int_t EvFirst, Int_t EvLast)    
{
  //
  // Read and print generated and reconstructed hits in EMC
  // for events from EvFirst to Nevent.
  // If only EvFirst is defined, print only this one event.
  // Author: Yuri Kharlov
  // 24 November 2000
  //

  if (EvFirst!=0 && EvLast==0) EvLast=EvFirst;
  Int_t ievent;
  for (ievent=EvFirst; ievent<=EvLast; ievent++) {  
    
    //========== Event Number>
    cout << endl <<  "==== ReadAndPrintEMC ====> Event is " << ievent+1 << endl ;

    //=========== Connects the various Tree's for evt
    Int_t ntracks = gAlice->GetEvent(ievent);
    fPHOS->SetTreeAddress() ;
    
    gAlice->TreeD()->GetEvent(0) ;
    gAlice->TreeR()->GetEvent(0) ;

    // Loop over reconstructed particles
      
    TClonesArray ** recParticleList  = fPHOS->RecParticles() ;     
    AliPHOSRecParticle * recParticle ;
    Int_t iRecParticle ;
    Int_t *primList;
    Int_t nPrimary;
    for(iRecParticle = 0; iRecParticle < (*recParticleList)->GetEntries() ;iRecParticle++ ) {
      recParticle = (AliPHOSRecParticle *) (*recParticleList)->At(iRecParticle) ;
      Float_t recE = recParticle->Energy();
      primList     = recParticle->GetPrimaries(nPrimary);
      Int_t moduleNumberRec ;
      Double_t recX, recZ ;
      fGeom->ImpactOnEmc(recParticle->Theta(), recParticle->Phi(), moduleNumberRec, recX, recZ) ;
      printf("Rec point: module %d, (X,Z) = (%8.4f,%8.4f) cm, E = %.3f GeV, primary = %d\n",
	     moduleNumberRec,recX,recZ,recE,*primList);
    }

    // Read and print EMC hits from EMCn branches
      
    AliPHOSCPVModule emcModule;
    TClonesArray    *emcHits;
    Int_t           nEMChits;
    AliPHOSCPVHit   *emcHit;
    TLorentzVector   p;
    Float_t          xgen, zgen;
    Int_t            ipart, primary;
    Int_t            nGenHits = 0;
    for (Int_t itrack=0; itrack<ntracks; itrack++) {
      //=========== Get the Hits Tree for the Primary track itrack
      gAlice->ResetHits();
      gAlice->TreeH()->GetEvent(itrack);
      Int_t iModule = 0 ;
      for (iModule=0; iModule < fGeom->GetNModules(); iModule++) {
	emcModule = fPHOS->GetEMCModule(iModule);
	emcHits   = emcModule.Hits();
	nEMChits  = emcHits->GetEntriesFast();
	for (Int_t ihit=0; ihit<nEMChits; ihit++) {
	  nGenHits++;
	  emcHit = (AliPHOSCPVHit*)emcHits->UncheckedAt(ihit);
	  p      = emcHit->GetMomentum();
	  xgen   = emcHit->X();
	  zgen   = emcHit->Y();
	  ipart  = emcHit->GetIpart();
	  primary= emcHit->GetTrack();
	  printf("EMC hit A: module %d, ",iModule+1);
	  printf("    p = (%f .4, %f .4, %f .4, %f .4) GeV,\n",
		 p.Px(),p.Py(),p.Pz(),p.Energy());
	  printf("                     (X,Z) = (%8.4f, %8.4f) cm, ipart = %d, primary = %d\n",
		 xgen,zgen,ipart,primary);
	}
      }
    }

//      // Read and print EMC hits from PHOS branch

//      for (Int_t itrack=0; itrack<ntracks; itrack++) {
//        //=========== Get the Hits Tree for the Primary track itrack
//        gAlice->ResetHits();
//        gAlice->TreeH()->GetEvent(itrack);
//        TClonesArray *hits = fPHOS->Hits();
//        AliPHOSHit   *hit ;
//        Int_t ihit;
//        for ( ihit = 0 ; ihit < hits->GetEntries() ; ihit++ ) {
//  	hit = (AliPHOSHit*)hits->At(ihit) ;
//  	Float_t hitXYZ[3];
//  	hitXYZ[0]   = hit->X();
//  	hitXYZ[1]   = hit->Y();
//  	hitXYZ[2]   = hit->Z();
//  	ipart       = hit->GetPid();
//  	primary     = hit->GetPrimary();
//  	Int_t absId = hit->GetId();
//  	Int_t relId[4];
//  	fGeom->AbsToRelNumbering(absId, relId) ;
//  	Int_t module = relId[0];
//  	if (relId[1]==0 && !(hitXYZ[0]==0 && hitXYZ[2]==0))
//  	  printf("EMC hit B: module %d, (X,Z) = (%8.4f, %8.4f) cm, ipart = %d, primary = %d\n",
//  		 module,hitXYZ[0],hitXYZ[2],ipart,primary);
//        }
//      }

  }
}

//____________________________________________________________________________
 void AliPHOSAnalyze::AnalyzeEMC(Int_t Nevents)
{
  //
  // Read generated and reconstructed hits in EMC for Nevents events.
  // Plots the coordinate and energy resolution histograms.
  // Coordinate resolution is a difference between the reconstructed
  // coordinate and the exact coordinate on the face of the PHOS
  // Author: Yuri Kharlov
  // 27 November 2000
  //

  // Book histograms

  TH1F *hDx1   = new TH1F("hDx1"  ,"EMC x-resolution", 100,-5. , 5.);
  TH1F *hDz1   = new TH1F("hDz1"  ,"EMC z-resolution", 100,-5. , 5.);
  TH1F *hDE1   = new TH1F("hDE1"  ,"EMC E-resolution", 100,-2. , 2.);

  TH2F *hDx2   = new TH2F("hDx2"  ,"EMC x-resolution", 100, 0., 10., 100,-5. , 5.);
  TH2F *hDz2   = new TH2F("hDz2"  ,"EMC z-resolution", 100, 0., 10., 100,-5. , 5.);
  TH2F *hDE2   = new TH2F("hDE2"  ,"EMC E-resolution", 100, 0., 10., 100, 0. , 5.);

  cout << "Start EMC Analysis"<< endl ;
  for (Int_t ievent=0; ievent<Nevents; ievent++) {  
      
    //========== Event Number>         
    if ( (ievent+1) % (Int_t)TMath::Power( 10, (Int_t)TMath::Log10(ievent+1) ) == 0)
      cout << "==== AnalyzeEMC ====> Event is " << ievent+1 << endl ;
    
    //=========== Connects the various Tree's for evt
    Int_t ntracks = gAlice->GetEvent(ievent);

    fPHOS->SetTreeAddress() ;
    
    gAlice->TreeD()->GetEvent(0) ;
    gAlice->TreeR()->GetEvent(0) ;

    // Create and fill arrays of hits for each EMC module
      
    Int_t nOfModules = fGeom->GetNModules();
    TClonesArray **hitsPerModule = new TClonesArray *[nOfModules];
    Int_t iModule;
    for (iModule=0; iModule < nOfModules; iModule++)
      hitsPerModule[iModule] = new TClonesArray("AliPHOSCPVHit",100);

    AliPHOSCPVModule emcModule;
    TClonesArray    *emcHits;
    Int_t           nEMChits;
    AliPHOSCPVHit   *emcHit;

    // First go through all primary tracks and fill the arrays
    // of hits per each EMC module

    for (Int_t itrack=0; itrack<ntracks; itrack++) {
      // Get the Hits Tree for the Primary track itrack
      gAlice->ResetHits();
      gAlice->TreeH()->GetEvent(itrack);
      for (Int_t iModule=0; iModule < nOfModules; iModule++) {
	emcModule = fPHOS->GetEMCModule(iModule);
	emcHits   = emcModule.Hits();
	nEMChits  = emcHits->GetEntriesFast();
	for (Int_t ihit=0; ihit<nEMChits; ihit++) {
	  emcHit   = (AliPHOSCPVHit*)emcHits->UncheckedAt(ihit);
	  TClonesArray &lhits = *(TClonesArray *)hitsPerModule[iModule];
	  new(lhits[hitsPerModule[iModule]->GetEntriesFast()]) AliPHOSCPVHit(*emcHit);
	}
	emcModule.Clear();
      }
    }

    // Loop over reconstructed particles
      
    TClonesArray ** recParticleList  = fPHOS->RecParticles() ;     
    AliPHOSRecParticle * recParticle ;
    Int_t nEMCrecs = (*recParticleList)->GetEntries();
    if (nEMCrecs == 1) {
      recParticle = (AliPHOSRecParticle *) (*recParticleList)->At(0) ;
      Float_t recE = recParticle->Energy();
      Int_t phosModule;
      Double_t recX, recZ ;
      fGeom->ImpactOnEmc(recParticle->Theta(), recParticle->Phi(), phosModule, recX, recZ) ;

      // for this rec.point take the hit list in the same PHOS module

      emcHits = hitsPerModule[phosModule-1];
      Int_t nEMChits  = emcHits->GetEntriesFast();
      if (nEMChits == 1) {
	Float_t genX, genZ, genE;
	for (Int_t ihit=0; ihit<nEMChits; ihit++) {
	  emcHit = (AliPHOSCPVHit*)emcHits->UncheckedAt(ihit);
	  genX   = emcHit->X();
	  genZ   = emcHit->Y();
	  genE   = emcHit->GetMomentum().E();
	}
	Float_t dx = recX - genX;
	Float_t dz = recZ - genZ;
	Float_t de = recE - genE;
	hDx1  ->Fill(dx);
	hDz1  ->Fill(dz);
	hDE1  ->Fill(de);
	hDx2  ->Fill(genE,dx);
	hDz2  ->Fill(genE,dz);
	hDE2  ->Fill(genE,recE);
      }
    }
    delete [] hitsPerModule;
  }
  // Save histograms

  Text_t outputname[80] ;
  sprintf(outputname,"%s.analyzed",fRootFile->GetName());
  TFile output(outputname,"RECREATE");
  output.cd();

  hDx1  ->Write() ;
  hDz1  ->Write() ;
  hDE1  ->Write() ;
  hDx2  ->Write() ;
  hDz2  ->Write() ;
  hDE2  ->Write() ;

  // Plot histograms

  TCanvas *emcCanvas = new TCanvas("EMC","EMC analysis",20,20,700,300);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1);
  gStyle->SetOptDate(1);
  emcCanvas->Divide(3,1);

  emcCanvas->cd(1);
  gPad->SetFillColor(10);
  hDx1->SetFillColor(16);
  hDx1->Draw();

  emcCanvas->cd(2);
  gPad->SetFillColor(10);
  hDz1->SetFillColor(16);
  hDz1->Draw();

  emcCanvas->cd(3);
  gPad->SetFillColor(10);
  hDE1->SetFillColor(16);
  hDE1->Draw();

  emcCanvas->Print("EMC.ps");

}

//____________________________________________________________________________
 void AliPHOSAnalyze::AnalyzeResolutions(Int_t Nevents )    
{
  // analyzes Nevents events and calculate Energy and Position resolution as well as
  // probaility of correct indentifiing of the incident particle

  //========== Booking Histograms
  cout << "AnalyzeResolutions > " << "Booking Histograms" << endl ; 
  BookResolutionHistograms();

  Int_t counter[9][5] ;     
  Int_t i1,i2,totalInd = 0 ;
  for(i1 = 0; i1<9; i1++)
    for(i2 = 0; i2<5; i2++)
      counter[i1][i2] = 0 ;
  
  Int_t totalPrimary = 0 ;
  Int_t totalRecPart = 0 ;
  Int_t totalRPwithPrim = 0 ;
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
      TClonesArray * primaryList  = gAlice->Particles();     

      TParticle * primary ;
      Int_t iPrimary ;
      for ( iPrimary = 0 ; iPrimary < primaryList->GetEntries() ; iPrimary++)
	{
	  primary = (TParticle*)primaryList->UncheckedAt(iPrimary) ;
	  Int_t primaryType = primary->GetPdgCode() ;
	  if( primaryType == 22 ) {
	    Int_t moduleNumber ;
	    Double_t primX, primZ ;
	    fGeom->ImpactOnEmc(primary->Theta(), primary->Phi(), moduleNumber, primX, primZ) ;
	    if(moduleNumber){
	      fhPrimary->Fill(primary->Energy()) ;
	      if(primary->Energy() > 0.3)
		totalPrimary++ ;
	    }
	  } 
	}
      
      fPHOS->SetTreeAddress() ;
      
      gAlice->TreeD()->GetEvent(0) ;
      gAlice->TreeR()->GetEvent(0) ;
      
      TClonesArray ** recParticleList  = fPHOS->RecParticles() ;     
      
      AliPHOSRecParticle * recParticle ;
      Int_t iRecParticle ;
      for(iRecParticle = 0; iRecParticle < (*recParticleList)->GetEntries() ;iRecParticle++ )
	{
	  recParticle = (AliPHOSRecParticle *) (*recParticleList)->At(iRecParticle) ;
	  fhAllRP->Fill(CorrectEnergy(recParticle->Energy())) ;
	  
	  Int_t moduleNumberRec ;
	  Double_t recX, recZ ;
	  fGeom->ImpactOnEmc(recParticle->Theta(), recParticle->Phi(), moduleNumberRec, recX, recZ) ;
	  
	  Double_t minDistance  = 100. ;
	  Int_t closestPrimary = -1 ;
	  
	  Int_t numberofprimaries ;
	  Int_t * listofprimaries  = recParticle->GetPrimaries(numberofprimaries)  ;
	  Int_t index ;
	  TParticle * primary ;
	  Double_t distance = minDistance ;
	  Double_t dX, dZ; 
	  Double_t dXmin = 0.; 
	  Double_t dZmin = 0. ;
	  for ( index = 0 ; index < numberofprimaries ; index++){
	    primary = (TParticle*)primaryList->UncheckedAt(listofprimaries[index]) ;
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
	  totalRecPart++ ;

	  if(closestPrimary >=0 ){
	    totalRPwithPrim++;
	    
	    Int_t primaryType = ((TParticle *)primaryList->At(closestPrimary))->GetPdgCode() ;
//  	    TParticlePDG* pDGparticle = ((TParticle *)primaryList->At(closestPrimary))->GetPDG();
//  	    Double_t charge =  PDGparticle->Charge() ;
// 	    if(charge)
// 	      cout <<"Primary " <<primaryType << " E " << ((TParticle *)primaryList->At(closestPrimary))->Energy() << endl ;
	    Int_t primaryCode ;
	    switch(primaryType)
	      {
	      case 22:
		primaryCode = 0;  //Photon
		fhAllEnergy   ->Fill(((TParticle *) primaryList->At(closestPrimary))->Energy(), recParticle->Energy()) ;
		fhAllPosition ->Fill(((TParticle *) primaryList->At(closestPrimary))->Energy(), minDistance) ;
		fhAllPositionX->Fill(dXmin);
		fhAllPositionZ->Fill(dZmin);
		break;
	      case 11 :
		primaryCode = 1;  //Electron
		break;
	      case -11 :
		primaryCode = 1;  //positron
		break;
	      case 321 :
		primaryCode = 4;  //K+
		break;
	      case -321 :
		primaryCode = 4;  //K-
		break;
	      case 310 :
		primaryCode = 4;  //K0s
		break;
	      case 130 :
		primaryCode = 4;  //K0l
		break;
	      case 211 :
		primaryCode = 2;  //K0l
		break;
	      case -211 :
		primaryCode = 2;  //K0l
		break;
	      case 2212 :
		primaryCode = 2;  //K0l
		break;
	      case -2212 :
		primaryCode = 2;  //K0l
		break;
	      default:
		primaryCode = 3; //ELSE
		break;
	      }
	    
	    switch(recParticle->GetType())
	      {
	      case AliPHOSFastRecParticle::kGAMMA:
		if(primaryType == 22){
		  fhPhotEnergy->Fill(((TParticle *) primaryList->At(closestPrimary))->Energy(), recParticle->Energy() ) ; 
		  fhEMEnergy->Fill(((TParticle *) primaryList->At(closestPrimary))->Energy(), recParticle->Energy() ) ; 
		  fhPPSDEnergy->Fill(((TParticle *) primaryList->At(closestPrimary))->Energy(), recParticle->Energy() ) ; 

		  fhPhotPosition->Fill(((TParticle *) primaryList->At(closestPrimary))->Energy(),minDistance) ;
		  fhEMPosition->Fill(((TParticle *) primaryList->At(closestPrimary))->Energy(),minDistance) ;
		  fhPPSDPosition->Fill(((TParticle *) primaryList->At(closestPrimary))->Energy(),minDistance) ;

		  fhPhotReg->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhPhotEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhPhotPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;

		  fhPhotPhot->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		}
		if(primaryType == 2112){ //neutron
		  fhNReg->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhNEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhNPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		}
		
		if(primaryType == -2112){ //neutron ~
		  fhNBarReg->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhNBarEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhNBarPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  
		}
		if(primaryCode == 2){
		  fhChargedReg->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhChargedEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhChargedPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		}
		
		fhAllReg->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		fhAllEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		fhAllPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		fhShape->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		fhVeto->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		fhPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		counter[0][primaryCode]++;
		break;
	      case  AliPHOSFastRecParticle::kELECTRON:
		if(primaryType == 22){ 
		  fhPhotElec->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhEMEnergy->Fill(((TParticle *) primaryList->At(closestPrimary))->Energy(), recParticle->Energy() ) ; 
		  fhEMPosition->Fill(((TParticle *) primaryList->At(closestPrimary))->Energy(),minDistance) ;
		  fhPhotEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhPhotPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		}	  
		if(primaryType == 2112){ //neutron
		  fhNEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhNPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		}
		
		if(primaryType == -2112){ //neutron ~
		  fhNBarEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhNBarPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  
		}
		if(primaryCode == 2){
		  fhChargedEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhChargedPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		}
		
		fhAllEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		fhAllPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		fhShape->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		fhPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		counter[1][primaryCode]++;
		break;
	      case  AliPHOSFastRecParticle::kNEUTRALHA:
		if(primaryType == 22) 
		  fhPhotNeuH->Fill(CorrectEnergy(recParticle->Energy()) ) ;

		fhVeto->Fill(CorrectEnergy(recParticle->Energy()) ) ;		
		counter[2][primaryCode]++;
		break ;
	      case  AliPHOSFastRecParticle::kNEUTRALEM:
		if(primaryType == 22){
		  fhEMEnergy->Fill(((TParticle *)primaryList->At(closestPrimary))->Energy(),recParticle->Energy() ) ; 
		  fhEMPosition->Fill(((TParticle *)primaryList->At(closestPrimary))->Energy(),minDistance ) ;
		
		  fhPhotNuEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhPhotEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		}
		if(primaryType == 2112) //neutron
		  fhNEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		
		if(primaryType == -2112) //neutron ~
		  fhNBarEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		
		if(primaryCode == 2)
		  fhChargedEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		
		fhAllEM->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		fhShape->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		fhVeto->Fill(CorrectEnergy(recParticle->Energy()) ) ;

		counter[3][primaryCode]++;
		break ;
	      case  AliPHOSFastRecParticle::kCHARGEDHA:
		if(primaryType == 22) //photon
		  fhPhotChHa->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		
		counter[4][primaryCode]++ ;
		break ;
	      case  AliPHOSFastRecParticle::kGAMMAHA:
		  if(primaryType == 22){ //photon
		    fhPhotGaHa->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		    fhPPSDEnergy->Fill(((TParticle *) primaryList->At(closestPrimary))->Energy(), recParticle->Energy() ) ; 
		    fhPPSDPosition->Fill(((TParticle *) primaryList->At(closestPrimary))->Energy(),minDistance) ;
		    fhPhotPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  }
		  if(primaryType == 2112){ //neutron
		    fhNPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  }
		
		  if(primaryType == -2112){ //neutron ~
		    fhNBarPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ; 
		  }
		  if(primaryCode == 2){
		    fhChargedPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  }
		
		  fhAllPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhVeto->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  fhPPSD->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		  counter[5][primaryCode]++ ;
		  break ;	
	      case  AliPHOSFastRecParticle::kABSURDEM:	      
		counter[6][primaryCode]++ ;
		fhShape->Fill(CorrectEnergy(recParticle->Energy()) ) ;
		break;
	      case  AliPHOSFastRecParticle::kABSURDHA:
		counter[7][primaryCode]++ ;
		break;
	      default:
		counter[8][primaryCode]++ ;
		break;
	      }
	  }
	}  
    }   // endfor
  SaveHistograms();
  cout << "Resolutions: Analyzed " << Nevents << " event(s)" << endl ;
  cout << "Resolutions: Total primary       " << totalPrimary << endl ;
  cout << "Resoluitons: Total reconstracted " << totalRecPart << endl ;
  cout << "TotalReconstructed with Primarie " << totalRPwithPrim << endl ;
  cout << "                        Primary:   Photon   Electron   Ch. Hadr.  Neutr. Hadr  Kaons" << endl ; 
  cout << "             Detected as photon       " << counter[0][0] << "          " << counter[0][1] << "          " << counter[0][2] << "          " <<counter[0][3] << "          " << counter[0][4] << endl ;
  cout << "           Detected as electron       " << counter[1][0] << "          " << counter[1][1] << "          " << counter[1][2] << "          " <<counter[1][3] << "          " << counter[1][4] << endl ; 
  cout << "     Detected as neutral hadron       " << counter[2][0] << "          " << counter[2][1] << "          " << counter[2][2] << "          " <<counter[2][3] << "          " << counter[2][4] << endl ;
  cout << "         Detected as neutral EM       " << counter[3][0] << "          " << counter[3][1] << "          " << counter[3][2] << "          " <<counter[3][3] << "          " << counter[3][4] << endl ;
  cout << "     Detected as charged hadron       " << counter[4][0] << "          " << counter[4][1] << "          " << counter[4][2] << "          " <<counter[4][3] << "          " << counter[4][4] << endl ;
  cout << "       Detected as gamma-hadron       " << counter[5][0] << "          " << counter[5][1] << "          " << counter[5][2] << "          " <<counter[5][3] << "          " << counter[5][4] << endl ;
  cout << "          Detected as Absurd EM       " << counter[6][0] << "          " << counter[6][1] << "          " << counter[6][2] << "          " <<counter[6][3] << "          " << counter[6][4] << endl ;
  cout << "      Detected as absurd hadron       " << counter[7][0] << "          " << counter[7][1] << "          " << counter[7][2] << "          " <<counter[7][3] << "          " << counter[7][4] << endl ;
  cout << "          Detected as undefined       " << counter[8][0] << "          " << counter[8][1] << "          " << counter[8][2] << "          " <<counter[8][3] << "          " << counter[8][4] << endl ;
      
      for(i1 = 0; i1<9; i1++)
	for(i2 = 0; i2<5; i2++)
	  totalInd+=counter[i1][i2] ;
      cout << "Indentified particles            " << totalInd << endl ;
      
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

  fhAllPositionX = new TH1F("hAllPositionX", "#Delta X of any RP with primary photon",100, -2., 2.);
  fhAllPositionZ = new TH1F("hAllPositionZ", "#Delta X of any RP with primary photon",100, -2., 2.);

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
  if(fhAllPositionX)
    fhAllPositionX->Write() ;
  if(fhAllPositionZ)
    fhAllPositionZ->Write() ;
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
   fhAllPositionX = 0 ; 
   fhAllPositionZ = 0 ; 
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
