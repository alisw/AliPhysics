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
//
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
      fPHOS  = (AliPHOSv5 *)gAlice->GetDetector("PHOS") ;
      fGeom  = AliPHOSGeometry::GetInstance( fPHOS->GetGeometry()->GetName(), fPHOS->GetGeometry()->GetTitle() );
 
      //========== Initializes the Index to Object converter
      fObjGetter = AliPHOSIndexToObject::GetInstance(fPHOS) ; 
      //========== Current event number 
      fEvt = -999 ; 

  }
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

  if (fRootFile->IsOpen() ) 
    fRootFile->Close() ; 
  if(fRootFile)
    delete fRootFile ; 

  if(fPHOS)
    delete fPHOS ; 

  if(fClu)
    delete fClu ; 

  if(fPID)
    delete fPID ; 

  if(fRec)
    delete fRec ; 

  if(fTrs)
    delete fTrs ; 

}

//____________________________________________________________________________
void AliPHOSAnalyze::ActivePPSD(Int_t Nevents=1){
  
  fhEnergyCorrelations  = new TH2F("hEnergyCorrelations","hEnergyCorrelations",40,  0., 0.15, 30, 0., 3.e-5);
  //========== Create the Clusterizer
  fClu = new AliPHOSClusterizerv1() ; 
  fClu->SetEmcEnergyThreshold(0.05) ; 
  fClu->SetEmcClusteringThreshold(0.20) ; 
  fClu->SetPpsdEnergyThreshold    (0.0000002) ; 
  fClu->SetPpsdClusteringThreshold(0.0000001) ; 
  fClu->SetLocalMaxCut(0.03) ;
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
 void AliPHOSAnalyze::Reconstruct(Int_t Nevents )    
{     
  Int_t ievent ;   
  for ( ievent=0; ievent<Nevents; ievent++)
    {  
      if (ievent==0) 
	{
	  cout << "Analyze > Starting Reconstructing " << endl ; 
	  //========== Create the Clusterizer
	  fClu = new AliPHOSClusterizerv1() ; 
	  fClu->SetEmcEnergyThreshold(0.05) ; 
	  fClu->SetEmcClusteringThreshold(0.20) ; 
	  fClu->SetPpsdEnergyThreshold    (0.0000002) ; 
	  fClu->SetPpsdClusteringThreshold(0.0000001) ; 
	  fClu->SetLocalMaxCut(0.03) ;
	  fClu->SetCalibrationParameters(0., 0.00000001) ; 
	  
	  //========== Creates the track segment maker
	  fTrs = new AliPHOSTrackSegmentMakerv1()  ;
	  //	  fTrs->UnsetUnfoldFlag() ; 
	  
	  //========== Creates the particle identifier
	  fPID = new AliPHOSPIDv1() ;
	  fPID->SetShowerProfileCuts(0.3, 1.8, 0.3, 1.8 ) ; 
	  
	  //========== Creates the Reconstructioner  
	  fRec = new AliPHOSReconstructioner(fClu, fTrs, fPID) ; 
	  //  fRec -> SetDebugReconstruction(kTRUE);     

	}
      
      //========== Event Number>         
      if ( ( log10((Float_t)(ievent+1)) - (Int_t)(log10((Float_t)(ievent+1))) ) == 0. ) 
	cout <<  "Analyze > Event is " << ievent << endl ;  
      
      //=========== Connects the various Tree's for evt
      gAlice->GetEvent(ievent);

      //=========== Gets the Digit TTree
      gAlice->TreeD()->GetEvent(0) ;

      //=========== Do the reconstruction
      fPHOS->Reconstruction(fRec);
    }

    fClu->Delete();
    fTrs->Delete();
    fPID->Delete();
    fRec->Delete();

}
//-------------------------------------------------------------------------------------

//   TClonesArray AllDigitArray = TClonesArray("AliPHOSDigit",1000) ;
//   TClonesArray * PhotonsList ;
//   TClonesArray * FalsDigitsList ;
//   TClonesArray AllPrimary = TClonesArray("TParticle",5000) ;
//   TFile * file2 = new TFile("ph100.root") ; // file with added photons
//   gAlice = (AliRun*) file2->Get("gAlice") ;
//   Int_t ievent;
//   Int_t NDigits[Nevents+1] ;
//   NDigits[0]=0 ;
//   Int_t NAllDigits = 0;
//   Int_t NprimPerEvent = 20 ;
//   for (ievent=0; ievent <Nevents; ievent++)
//     {
//       PhotonsList  = gAlice->Particles();  //Primary
//       FalsDigitsList  = ((AliPHOSv1 *)gAlice->GetDetector("PHOS"))->Digits();  //Digits
//       gAlice->GetEvent(ievent) ;
//       gAlice->TreeD()->GetEvent(0) ;
//       gAlice->TreeK()->GetEvent(0) ;
//       //Copy Primary
//       Int_t Nprim ;
//       for(Nprim = 0 ;Nprim < NprimPerEvent ; Nprim++)
// 	new (AllPrimary[Nprim+ievent*NprimPerEvent])  TParticle(*((TParticle *) PhotonsList->At(Nprim))) ;

//       //Copy Digits
//       TIter nextDigit(FalsDigitsList) ;
//       AliPHOSDigit * FalseDigit ;
//       NDigits[ievent+1] = NDigits[ievent]+ FalsDigitsList->GetEntriesFast() ; 
//       while( (FalseDigit = (AliPHOSDigit *) nextDigit()))
// 	{	 
// 	  new (AllDigitArray[NAllDigits])  AliPHOSDigit(FalseDigit->GetPrimary(1),FalseDigit->GetId(),FalseDigit->GetAmp()) ;
// 	  NAllDigits++ ;
// 	}	  
//     }
//   file2->Close() ;



//       //Add primary particles
//       cout << "# of Primaries before add " << PrimaryList->GetEntriesFast() << endl;
//      Int_t NTruePrimary = 0 ;  //PrimaryList->GetEntriesFast() ;
//       Int_t Nprim ;
//       for(Nprim = 0; Nprim < NprimPerEvent; Nprim++)
// 	new ((*PrimaryList)[NTruePrimary+Nprim])  TParticle(*((TParticle *) AllPrimary.At(Nprim+ievent*NprimPerEvent))) ;
      
//       cout << "# of Primaries after add " << PrimaryList->GetEntriesFast() <<endl;

//       cout << "Digits before add " << DigitsList->GetEntries() << endl ;
//       cout << "Digits to add " <<  NDigits[ievent+1]-  NDigits[ievent]<< endl ;
      
      //=========== Add fals digits ==============================
//       TIter nextDigit(DigitsList) ;
//       AliPHOSDigit * FalseDigit ;
//       AliPHOSDigit * RealDigit ;
//       Int_t NTrueDigits = DigitsList->GetEntriesFast() ; 
//       Int_t Ndigit ;
//       for(Ndigit=NDigits[ievent];Ndigit<NDigits[ievent+1];Ndigit++)
// 	{	 
// 	  FalseDigit = (AliPHOSDigit*) AllDigitArray.At(Ndigit) ;
// 	  Bool_t Add = kTRUE ; 
// 	  AliPHOSDigit tmpDigit=AliPHOSDigit(FalseDigit->GetPrimary(1)+NTruePrimary,FalseDigit->GetId(),FalseDigit->GetAmp()) ;
	  
// 	  while( (RealDigit = (AliPHOSDigit *) nextDigit()) && Add)
// 	    {
// 	      if((*RealDigit) == (tmpDigit)) 
// 		{
// 		  *RealDigit=*RealDigit+tmpDigit ;
// 		  Add = kFALSE ;
// 		}
// 	    }
// 	  if(Add)
// 	    {
// 	      new ((*DigitsList)[NTrueDigits])  AliPHOSDigit(FalseDigit->GetPrimary(1)+NTruePrimary,FalseDigit->GetId(),FalseDigit->GetAmp()) ;
// 	      ((AliPHOSDigit *)DigitsList->At(NTrueDigits))->SetIndexInList(NTrueDigits) ;
// 	      NTrueDigits++ ;
// 	    }	  
// 	}
//       cout << "Digits after add " << DigitsList->GetEntries() << endl ;


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
      if ( ( log10((Float_t)(ievent+1)) - (Int_t)(log10((Float_t)(ievent+1))) ) == 0. ) 
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
	      if( PrimaryType == 22 ) 
		fhPrimary->Fill(Primary->Energy()) ;
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
 
      cout << ievent << " " << (*EmcRecPoints) << " " <<(*PpsdRecPoints) <<fPHOS->Digits()<< endl ; 
      cout << "    " << " " << (*EmcRecPoints)->GetEntries() << " " <<(*PpsdRecPoints)->GetEntries() <<fPHOS->Digits()->GetEntries()<< endl ; 

      AliPHOSRecParticle * RecParticle ;
      Int_t iRecParticle ;
      for(iRecParticle = 0; iRecParticle < (*RecParticleList)->GetEntries() ;iRecParticle++ )
	{
	  RecParticle = (AliPHOSRecParticle *) (*RecParticleList)->At(iRecParticle) ;
	  
	  Int_t ModuleNumberRec ;
	  Double_t RecX, RecZ ;
	  fGeom->ImpactOnEmc(RecParticle->Theta(), RecParticle->Phi(), ModuleNumberRec, RecX, RecZ) ;
	  
	  Double_t MinDistance = 10000 ;
	  Int_t ClosestPrimary = -1 ;
	  
	  Int_t numberofprimaries ;
	  Int_t * listofprimaries  = RecParticle->GetPrimaries(numberofprimaries)  ;
	  Int_t index ;
	  TParticle * Primary ;
	  Double_t Distance = MinDistance ;
	  for ( index = 0 ; index < numberofprimaries ; index++)
	    {
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
	  if(ClosestPrimary >=0 )
	    {
	      fhPhotonAllEnergy->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(), RecParticle->Energy() ) ; 
	      fhPhotonAllPosition->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(),Distance) ;
	      TotalRPwithPrim++;
	      Int_t PrimaryType = ((TParticle *)PrimaryList->At(ClosestPrimary))->GetPdgCode() ;
	      TParticlePDG* PDGparticle = ((TParticle *)PrimaryList->At(ClosestPrimary))->GetPDG();
	      Double_t charge =  PDGparticle->Charge() ;
	      Int_t PrimaryCode ;
	      switch(PrimaryType)
		{
		case 22:
		  PrimaryCode = 0;  //Photon
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
		default:
		  if(charge)
		    PrimaryCode = 2; //Charged hadron
		  else
		    PrimaryCode = 3; //Neutral hadron
		  break;
		}

	      switch(RecParticle->GetType())
		{
		case AliPHOSFastRecParticle::kGAMMA:
		  if(PrimaryType == 22){
		    fhPhotonEnergy->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(), RecParticle->Energy() ) ; 
		    fhPhotonPosition->Fill(((TParticle *) PrimaryList->At(ClosestPrimary))->Energy(),Distance) ;
		    fhPhotonReg->Fill(RecParticle->Energy() ) ;
		    fhPhotonEM->Fill(RecParticle->Energy() ) ;
		    fhPhotPhot->Fill(RecParticle->Energy() ) ;
		  }
		  if(PrimaryType == 2112){ //neutron
		    fhNReg->Fill(RecParticle->Energy() ) ;
		    fhNEM->Fill(RecParticle->Energy() ) ;
		  }
		  
		  if(PrimaryType == -2112){ //neutron ~
		    fhNBarReg->Fill(RecParticle->Energy() ) ;
		    fhNBarEM->Fill(RecParticle->Energy() ) ;

		  }
		  if(PrimaryCode == 2){
		    fhChargedReg->Fill(RecParticle->Energy() ) ;
		    fhChargedEM->Fill(RecParticle->Energy() ) ;
		  }

		  fhAllReg->Fill(RecParticle->Energy() ) ;
		  fhAllEM->Fill(RecParticle->Energy() ) ;
		  Counter[0][PrimaryCode]++;
		  break;
		case  AliPHOSFastRecParticle::kELECTRON:
		  if(PrimaryType == 11 || PrimaryType == -11){
		    fhElectronEnergy->Fill(((TParticle *)PrimaryList->At(ClosestPrimary))->Energy(), RecParticle->Energy() ) ; 
		    fhElectronPosition->Fill(((TParticle *)PrimaryList->At(ClosestPrimary))->Energy(),Distance ) ;
		  }
		  if(PrimaryType == 22) 
		    fhPhotElec->Fill(RecParticle->Energy() ) ;

		  Counter[1][PrimaryCode]++;
		  break;
		case  AliPHOSFastRecParticle::kNEUTRALHA:
		  if(PrimaryType == 22) 
		    fhPhotNeuH->Fill(RecParticle->Energy() ) ;

		  fhNeutralHadronEnergy->Fill( ((TParticle *)PrimaryList->At(ClosestPrimary))->Energy(), RecParticle->Energy()) ; 
		  fhNeutralHadronPosition->Fill(((TParticle *)PrimaryList->At(ClosestPrimary))->Energy() ,Distance  ) ;
		  Counter[2][PrimaryCode]++;
		  break ;
		case  AliPHOSFastRecParticle::kNEUTRALEM:
		  if(PrimaryType == 22 || PrimaryType == 11 || PrimaryType == -11){
		    fhNeutralEMEnergy->Fill(((TParticle *)PrimaryList->At(ClosestPrimary))->Energy(),RecParticle->Energy() ) ; 
		    fhNeutralEMPosition->Fill(((TParticle *)PrimaryList->At(ClosestPrimary))->Energy(),Distance ) ;
		  }

		  if(PrimaryType == 22){ //photon
		    fhPhotNuEM->Fill(RecParticle->Energy() ) ;
		    fhPhotonEM->Fill(RecParticle->Energy() ) ;
		  }
		  if(PrimaryType == 2112) //neutron
		    fhNEM->Fill(RecParticle->Energy() ) ;
		  
		  if(PrimaryType == -2112) //neutron ~
		    fhNBarEM->Fill(RecParticle->Energy() ) ;

		  if(PrimaryCode == 2)
		    fhChargedEM->Fill(RecParticle->Energy() ) ;

		  fhAllEM->Fill(RecParticle->Energy() ) ;

		  Counter[3][PrimaryCode]++;
		  break ;
		case  AliPHOSFastRecParticle::kCHARGEDHA:
		  if(PrimaryType == 22) //photon
		    fhPhotChHa->Fill(RecParticle->Energy() ) ;
		  
		  fhChargedHadronEnergy->Fill(((TParticle *)PrimaryList->At(ClosestPrimary))->Energy(),RecParticle->Energy() ) ; 
		  fhChargedHadronPosition->Fill(((TParticle *)PrimaryList->At(ClosestPrimary))->Energy(),Distance ) ;
		  Counter[4][PrimaryCode]++ ;
		  break ;
		case  AliPHOSFastRecParticle::kGAMMAHA:
		  if(PrimaryType == 22) //photon
		    fhPhotGaHa->Fill(RecParticle->Energy() ) ;
		  fhPhotonHadronEnergy->Fill(((TParticle *)PrimaryList->At(ClosestPrimary))->Energy(), RecParticle->Energy()) ; 
		  fhPhotonHadronPosition->Fill(((TParticle *)PrimaryList->At(ClosestPrimary))->Energy(),Distance ) ;
		  Counter[5][PrimaryCode]++ ;
		  break ;	
		case  AliPHOSFastRecParticle::kABSURDEM:	      
		  Counter[6][PrimaryCode]++ ;
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

  if(fhPhotonEnergy)
    delete fhPhotonEnergy ;
  if(fhPhotonAllEnergy)
    delete fhPhotonAllEnergy ;
  if(fhElectronEnergy)
    delete fhElectronEnergy ;
  if(fhElectronAllEnergy)
    delete fhElectronAllEnergy ;
  if(fhNeutralHadronEnergy)
    delete fhNeutralHadronEnergy ;
  if(fhNeutralEMEnergy)
    delete fhNeutralEMEnergy ;
  if(fhNeutralEMAllEnergy)
    delete fhNeutralEMAllEnergy ;
  if(fhChargedHadronEnergy)
    delete fhChargedHadronEnergy ;
  if(fhPhotonHadronEnergy)
    delete fhPhotonHadronEnergy ;
  if(fhPhotonPosition)
    delete fhPhotonPosition ;
  if(fhPhotonAllPosition)
    delete fhPhotonAllPosition ;
  if(fhElectronPosition)
    delete fhElectronPosition ;
  if(fhElectronAllPosition)
    delete fhElectronAllPosition ;
  if(fhNeutralHadronPosition)
    delete fhNeutralHadronPosition ;
  if(fhNeutralEMPosition)
    delete fhNeutralEMPosition ;
  if(fhNeutralEMAllPosition)
    delete fhNeutralEMAllPosition ;
  if(fhChargedHadronPosition)
    delete fhChargedHadronPosition ;
  if(fhPhotonHadronPosition)
    delete fhPhotonHadronPosition ;

  fhPhotonEnergy            = new TH2F("hPhotonEnergy",  "hPhotonEnergy",              100, 0., 5., 100, 0., 5.);
  fhPhotonAllEnergy         = new TH2F("hPhotonAllEnergy",  "hPhotonAllEnergy",        100, 0., 5., 100, 0., 5.);
  fhElectronEnergy          = new TH2F("hElectronEnergy","hElectronEnergy",            100, 0., 5., 100, 0., 5.);
  fhElectronAllEnergy       = new TH2F("hElectronAllEnergy","hElectronAllEnergy",      100, 0., 5., 100, 0., 5.);
  fhNeutralHadronEnergy     = new TH2F("hNeutralHadronEnergy", "hNeutralHadronEnergy", 100, 0., 5., 100, 0., 5.);
  fhNeutralEMEnergy         = new TH2F("hNeutralEMEnergy", "hNeutralEMEnergy",         100, 0., 5., 100, 0., 5.);
  fhNeutralEMAllEnergy      = new TH2F("hNeutralEMAllEnergy", "hNeutralEMAllEnergy",   100, 0., 5., 100, 0., 5.);
  fhChargedHadronEnergy     = new TH2F("hChargedHadronEnergy", "hChargedHadronEnergy", 100, 0., 5., 100, 0., 5.);
  fhPhotonHadronEnergy      = new TH2F("hPhotonHadronEnergy","hPhotonHadronEnergy",    100, 0., 5., 100, 0., 5.);
  fhPhotonPosition          = new TH2F("hPhotonPosition","hPhotonPosition",                20, 0., 5., 100, 0., 5.);
  fhPhotonAllPosition       = new TH2F("hPhotonAllPosition","hPhotonAllPosition",          20, 0., 5., 100, 0., 5.);
  fhElectronPosition        = new TH2F("hElectronPosition","hElectronPosition",            20, 0., 5., 100, 0., 5.);
  fhElectronAllPosition     = new TH2F("hElectronAllPosition","hElectronAllPosition",      20, 0., 5., 100, 0., 5.);
  fhNeutralHadronPosition   = new TH2F("hNeutralHadronPosition","hNeutralHadronPosition",  20, 0., 5., 100, 0., 5.);
  fhNeutralEMPosition       = new TH2F("hNeutralEMPosition","hNeutralEMPosition",          20, 0., 5., 100, 0., 5.);
  fhNeutralEMAllPosition    = new TH2F("hNeutralEMAllPosition","hNeutralEMAllPosition",    20, 0., 5., 100, 0., 5.);
  fhChargedHadronPosition   = new TH2F("hChargedHadronPosition","hChargedHadronPosition",  20, 0., 5., 100, 0., 5.);
  fhPhotonHadronPosition    = new TH2F("hPhotonHadronPosition","hPhotonHadronPosition",    20, 0., 5., 100, 0., 5.);

  if(fhPhotonReg)
    delete fhPhotonReg ;
  if(fhAllReg)
    delete fhAllReg ;
  if(fhNReg)
    delete fhNReg ;
  if(fhNReg)
    delete fhNReg ;
  if(fhNReg)
    delete fhNReg ;
  
  fhPhotonReg = new TH1F("hPhotonReg","hPhotonReg", 20, 0., 5.);
  fhAllReg    = new TH1F("hAllReg", "hAllReg",  20, 0., 5.);
  fhNReg      = new TH1F("hNReg", "hNReg",  20, 0., 5.);
  fhNBarReg   = new TH1F("hNBarReg", "hNBarReg",  20, 0., 5.);
  fhChargedReg= new TH1F("hChargedReg", "hChargedReg",  20, 0., 5.);
  
  if(fhPhotonEM)
    delete fhPhotonEM ;
  if(fhAllEM)
    delete fhAllEM ;
  if(fhNEM)
    delete fhNEM ;
  if(fhNBarEM)
    delete fhNBarEM ;
  if(fhChargedEM)
    delete fhChargedEM ;
  
  fhPhotonEM = new TH1F("hPhotonEM","hPhotonEM", 20, 0., 5.);
  fhAllEM    = new TH1F("hAllEM", "hAllEM",  20, 0., 5.);
  fhNEM      = new TH1F("hNEM", "hNEM",  20, 0., 5.);
  fhNBarEM   = new TH1F("hNBarEM", "hNBarEM",  20, 0., 5.);
  fhChargedEM= new TH1F("hChargedEM", "hChargedEM",  20, 0., 5.);
  
  if(fhPrimary)
    delete fhPrimary ;
  fhPrimary= new TH1F("hPrimary", "hPrimary",  20, 0., 5.);

  if(fhPhotPhot)
    delete fhPhotPhot ;
  if(fhPhotElec)
    delete fhPhotElec ;
  if(fhPhotNeuH)
    delete fhPhotNeuH ;
  if(fhPhotNuEM)
    delete fhPhotNuEM ;
  if(fhPhotChHa)
    delete fhPhotChHa ;
  if(fhPhotGaHa)
    delete fhPhotGaHa ;

  fhPhotPhot = new TH1F("hPhotPhot","hPhotPhot", 20, 0., 5.);   //Photon registered as photon
  fhPhotElec = new TH1F("hPhotElec","hPhotElec", 20, 0., 5.);   //Photon registered as Electron
  fhPhotNeuH = new TH1F("hPhotNeuH","hPhotNeuH", 20, 0., 5.);   //Photon registered as Neutral Hadron
  fhPhotNuEM = new TH1F("hPhotNuEM","hPhotNuEM", 20, 0., 5.);   //Photon registered as Neutral EM
  fhPhotChHa = new TH1F("hPhotChHa","hPhotChHa", 20, 0., 5.);   //Photon registered as Charged Hadron
  fhPhotGaHa = new TH1F("hPhotGaHa","hPhotGaHa", 20, 0., 5.);   //Photon registered as Gamma-Hadron


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
      
      fPHOS  = (AliPHOSv5 *)gAlice->GetDetector("PHOS") ;
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
    fRec -> SetDebugReconstruction(kFALSE);     
    
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
// void AliPHOSAnalyze::SavingHistograms()
// {
//   // Saves the histograms in a root file named "name.analyzed" 

//   Text_t outputname[80] ;
//   sprintf(outputname,"%s.analyzed",fRootFile->GetName());
//   TFile output(outputname,"RECREATE");
//   output.cd();
//   if (fhEmcDigit )         
//     fhEmcDigit->Write()  ;
//   if (fhVetoDigit )  
//     fhVetoDigit->Write()  ;
//   if (fhConvertorDigit ) 
//     fhConvertorDigit->Write()   ;
//   if (fhEmcCluster   )
//     fhEmcCluster->Write()   ;
//   if (fhVetoCluster ) 
//     fhVetoCluster->Write()   ;
//   if (fhConvertorCluster )
//     fhConvertorCluster->Write()  ;
//   if (fhConvertorEmc ) 
//     fhConvertorEmc->Write()  ;
//   if (fhPhotonEnergy)    
//     fhPhotonEnergy->Write() ;
//   if (fhPhotonPositionX)  
//     fhPhotonPositionX->Write() ;
//   if (fhPhotonPositionY)  
//     fhPhotonPositionX->Write() ;
//   if (fhElectronEnergy)  
//     fhElectronEnergy->Write() ;
//   if (fhElectronPositionX)
//     fhElectronPositionX->Write() ;
//   if (fhElectronPositionY) 
//     fhElectronPositionX->Write() ;
//   if (fhNeutralHadronEnergy) 
//     fhNeutralHadronEnergy->Write() ;
//   if (fhNeutralHadronPositionX)
//     fhNeutralHadronPositionX->Write() ;
//   if (fhNeutralHadronPositionY) 
//     fhNeutralHadronPositionX->Write() ;
//   if (fhNeutralEMEnergy)   
//     fhNeutralEMEnergy->Write() ;
//   if (fhNeutralEMPositionX)
//     fhNeutralEMPositionX->Write() ;
//   if (fhNeutralEMPositionY) 
//     fhNeutralEMPositionX->Write() ;
//   if (fhChargedHadronEnergy) 
//     fhChargedHadronEnergy->Write() ;
//   if (fhChargedHadronPositionX) 
//     fhChargedHadronPositionX->Write() ;
//   if (fhChargedHadronPositionY)
//     fhChargedHadronPositionX->Write() ;
//   if (fhPhotonHadronEnergy) 
//     fhPhotonHadronEnergy->Write() ;
//   if (fhPhotonHadronPositionX) 
//     fhPhotonHadronPositionX->Write() ;
//   if (fhPhotonHadronPositionY)
//     fhPhotonHadronPositionX->Write() ;

//   output.Write();
//   output.Close();
// }
//____________________________________________________________________________
void AliPHOSAnalyze::SaveHistograms()
{
  // Saves the histograms in a root file named "name.analyzed" 

  Text_t outputname[80] ;
  sprintf(outputname,"%s.analyzed",fRootFile->GetName());
  TFile output(outputname,"RECREATE");
  output.cd();

  if (fhPhotonEnergy)    
    fhPhotonEnergy->Write() ;
  if (fhPhotonAllEnergy)    
    fhPhotonAllEnergy->Write() ;
  if (fhPhotonPosition)  
    fhPhotonPosition->Write() ;
  if (fhPhotonAllPosition)  
    fhPhotonAllPosition->Write() ;
  if (fhElectronEnergy)  
    fhElectronEnergy->Write() ;
  if (fhElectronAllEnergy)  
    fhElectronAllEnergy->Write() ;
  if (fhElectronPosition)
    fhElectronPosition->Write() ;
  if (fhElectronAllPosition)
    fhElectronAllPosition->Write() ;
  if (fhNeutralHadronEnergy) 
    fhNeutralHadronEnergy->Write() ;
  if (fhNeutralHadronPosition)
    fhNeutralHadronPosition->Write() ;
  if (fhNeutralEMEnergy)   
    fhNeutralEMEnergy->Write() ;
  if (fhNeutralEMAllEnergy)   
    fhNeutralEMAllEnergy->Write() ;
  if (fhNeutralEMPosition)
    fhNeutralEMPosition->Write() ;
  if (fhNeutralEMAllPosition)
    fhNeutralEMAllPosition->Write() ;
  if (fhChargedHadronEnergy) 
    fhChargedHadronEnergy->Write() ;
  if (fhChargedHadronPosition) 
    fhChargedHadronPosition->Write() ;
  if (fhPhotonHadronEnergy) 
    fhPhotonHadronEnergy->Write() ;
  if (fhPhotonHadronPosition) 
    fhPhotonHadronPosition->Write() ;
  if (fhPhotonReg) 
    fhPhotonReg->Write() ;
  if (fhAllReg) 
    fhAllReg->Write() ;
  if(fhNReg)
    fhNReg->Write() ;
  if(fhNBarReg)
    fhNBarReg->Write() ;
  if(fhChargedReg)
    fhChargedReg->Write() ;
  if (fhPhotonEM) 
    fhPhotonEM->Write() ;
  if (fhAllEM) 
    fhAllEM->Write() ;
  if(fhNEM)
    fhNEM->Write() ;
  if(fhNBarEM)
    fhNBarEM->Write() ;
  if(fhChargedEM)
    fhChargedEM->Write() ;
  if(fhPrimary)
    fhPrimary->Write() ;
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

   fhPhotonEnergy = 0 ;           // Spectrum of detected photons with photon primary
   fhPhotonAllEnergy = 0 ;        // Total spectrum of detected photons
   fhElectronEnergy = 0 ;         // Spectrum of detected electrons with electron primary
   fhElectronAllEnergy = 0 ;      // Total spectrum of detected electrons
   fhNeutralHadronEnergy = 0 ;    // Spectrum of detected neutral hadron
   fhNeutralEMEnergy = 0 ;        // Spectrum of detected neutral EM with EM primary
   fhNeutralEMAllEnergy = 0 ;     // Spectrum of detected neutral EM
   fhChargedHadronEnergy = 0 ;    // Spectrum of detected charged
   fhPhotonHadronEnergy = 0 ;     // Spectrum of detected Photon-Hadron
   fhPhotonPosition = 0 ;        // Position Resolution of  photons with photon primary
   fhPhotonAllPosition = 0 ;     // Position Resolution of  photons
   fhElectronPosition = 0 ;      // Position Resolution of electrons with electron primary
   fhElectronAllPosition = 0 ;   // Position Resolution of electrons
   fhNeutralHadronPosition = 0 ; // Position Resolution of neutral hadron
   fhNeutralEMPosition = 0 ;     // Position Resolution of neutral EM with EM primary
   fhNeutralEMAllPosition = 0 ;  // Position Resolution of neutral EM
   fhChargedHadronPosition = 0 ; // Position Resolution of charged
   fhPhotonHadronPosition = 0 ;  // Position Resolution of Photon-Hadron
   fhPhotonPositionY = 0 ;        // Y distribution of detected photons
   fhElectronPositionY = 0 ;      // Y distribution of detected electrons
   fhNeutralHadronPositionY = 0 ; // Y distribution of detected neutral hadron
   fhNeutralEMPositionY = 0 ;     // Y distribution of detected neutral EM
   fhChargedHadronPositionY = 0 ; // Y distribution of detected charged
   fhPhotonHadronPositionY = 0 ;  // Y distribution of detected Photon-Hadron
   fhPhotonReg = 0 ;          
   fhAllReg = 0 ;          
   fhNReg = 0 ;          
   fhNBarReg = 0 ;          
   fhChargedReg = 0 ;          
   fhPhotonEM = 0 ;          
   fhAllEM = 0 ;          
   fhNEM = 0 ;          
   fhNBarEM = 0 ;          
   fhChargedEM = 0 ;          
   fhPrimary = 0 ;          

   fhPhotPhot = 0 ;
   fhPhotElec = 0 ;
   fhPhotNeuH = 0 ;
   fhPhotNuEM = 0 ; 
   fhPhotChHa = 0 ;
   fhPhotGaHa = 0 ;


}


