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
//*--
//*-- Yves Schutz (SUBATECH) 
// Reconstruction class. Redesigned from the old AliReconstructionner class and 
// derived from STEER/AliReconstructor. 
// 
// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliESD.h"
#include "AliPHOSReconstructor.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSPIDv1.h"
#include "AliPHOSGetter.h"
#include "AliPHOSTracker.h"
#include "AliRawReader.h"

 
ClassImp(AliPHOSReconstructor)

Bool_t AliPHOSReconstructor::fgDebug = kFALSE ; 

//____________________________________________________________________________
  AliPHOSReconstructor::AliPHOSReconstructor() 
{
  // ctor

} 

//____________________________________________________________________________
  AliPHOSReconstructor::~AliPHOSReconstructor()
{
  // dtor

} 

//____________________________________________________________________________
void AliPHOSReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
  // method called by AliReconstruction; 
  // Only the clusterization is performed,; the rest of the reconstruction is done in FillESD because the track
  // segment maker needs access to the AliESD object to retrieve the tracks reconstructed by 
  // the global tracking.
 
  TString headerFile(runLoader->GetFileName()) ; 
  TString branchName(runLoader->GetEventFolder()->GetName()) ;  
  
  AliPHOSClusterizerv1 clu(headerFile, branchName);
  clu.SetEventRange(0, -1) ; // do all the events
  if ( Debug() ) 
    clu.ExecuteTask("deb all") ; 
  else 
    clu.ExecuteTask("") ;  

}

//____________________________________________________________________________
void AliPHOSReconstructor::Reconstruct(AliRunLoader* runLoader, AliRawReader* rawreader) const
{
  // method called by AliReconstruction; 
  // Only the clusterization is performed,; the rest of the reconstruction is done in FillESD because the track
  // segment maker needs access to the AliESD object to retrieve the tracks reconstructed by 
  // the global tracking.
  // Here we reconstruct from Raw Data

  rawreader->Reset() ; 
  TString headerFile(runLoader->GetFileName()) ; 
  TString branchName(runLoader->GetEventFolder()->GetName()) ;  
  
  AliPHOSClusterizerv1 clu(headerFile, branchName);
  clu.SetEventRange(0, -1) ; // do all the events
  clu.SetRawReader(rawreader);

  TString option = GetOption();
  if (option.Contains("OldRCUFormat"))
    clu.SetOldRCUFormat(kTRUE);

  if ( Debug() ) 
    clu.ExecuteTask("deb all") ; 
  else 
    clu.ExecuteTask("") ;

}

//____________________________________________________________________________
void AliPHOSReconstructor::FillESD(AliRunLoader* runLoader, AliESD* esd) const
{
  // This function creates AliESDtracks from AliPHOSRecParticles
  //         and
  // writes them to the ESD

  Int_t eventNumber = runLoader->GetEventNumber() ;

  AliPHOSGetter *gime = AliPHOSGetter::Instance();
  gime->Event(eventNumber, "DRTP") ; 
  TClonesArray *recParticles  = gime->RecParticles();
  Int_t nOfRecParticles = recParticles->GetEntries();

  esd->SetNumberOfPHOSClusters(nOfRecParticles) ; 
  esd->SetFirstPHOSCluster(esd->GetNumberOfCaloClusters()) ;
  
  AliDebug(2,Form("%d digits and %d rec. particles in event %d, option %s",gime->Digits()->GetEntries(),nOfRecParticles,eventNumber,GetOption()));

  for (Int_t recpart = 0 ; recpart < nOfRecParticles ; recpart++) {
    AliPHOSRecParticle * rp = dynamic_cast<AliPHOSRecParticle*>(recParticles->At(recpart));
    if (Debug()) 
      rp->Print();
    // Get track segment and EMC rec.point associated with this rec.particle
    AliPHOSTrackSegment *ts    = gime->TrackSegment(rp->GetPHOSTSIndex());
    AliPHOSEmcRecPoint  *emcRP = gime->EmcRecPoint(ts->GetEmcIndex());
    AliESDCaloCluster   *ec    = new AliESDCaloCluster() ; 

    // fills the ESDCaloCluster
    Float_t xyz[3];
    for (Int_t ixyz=0; ixyz<3; ixyz++) 
      xyz[ixyz] = rp->GetPos()[ixyz];
    
    AliDebug(2,Form("Global position xyz=(%f,%f,%f)",xyz[0],xyz[1],xyz[2]));
    
    Int_t  digitMult  = emcRP->GetDigitsMultiplicity();
    Int_t *digitsList = emcRP->GetDigitsList();
    UShort_t *amplList  = new UShort_t[digitMult];
    UShort_t *timeList  = new UShort_t[digitMult];
    UShort_t *digiList  = new UShort_t[digitMult];

    // Convert Float_t* and Int_t* to UShort_t* to save memory
    for (Int_t iDigit=0; iDigit<digitMult; iDigit++) {
      AliPHOSDigit *digit = gime->Digit(digitsList[iDigit]);
      amplList[iDigit] = (UShort_t)(digit->GetEnergy()*500); // Energy in units of GeV/500
      timeList[iDigit] = (UShort_t)(digit->GetTime()*1e9*100); // time in units of 0.01 ns
      digiList[iDigit] = (UShort_t)(digit->GetId());
    }
    ec->SetGlobalPosition(xyz);                 //rec.point position in MARS
    ec->SetClusterEnergy(rp->Energy());         //total particle energy
    ec->SetClusterDisp(emcRP->GetDispersion()); //cluster dispersion
    ec->SetPid          (rp->GetPID()) ;        //array of particle identification
    ec->SetPrimaryIndex (rp->GetPrimaryIndex());//index of primary particle (for simulations)
    ec->SetM02(emcRP->GetM2x()) ;               //second moment M2x
    ec->SetM20(emcRP->GetM2z()) ;               //second moment M2z
    ec->SetNExMax(emcRP->GetNExMax());          //number of local maxima
    ec->SetNumberOfDigits(digitMult);           //digit multiplicity
    ec->SetDigitAmplitude(amplList);            //energies in 1/500 of GeV
    ec->SetDigitTime(timeList);                 //times in 1/100 on ns
    ec->SetDigitIndex(digiList);                //abs id of the cell
    ec->SetEmcCpvDistance(-1);                  //not yet implemented
    ec->SetClusterChi2(-1);                     //not yet implemented
    ec->SetM11(-1) ;                            //not yet implemented

    // add the track to the esd object
    esd->AddCaloCluster(ec);
    delete ec;    
  }  
}

void AliPHOSReconstructor::FillESD(AliRunLoader* runLoader,
				   AliRawReader* rawReader, AliESD* esd) const
{
  //This function creates AliESDtracks from AliPHOSRecParticles 
  //and writes them to the ESD in the case of raw data reconstruction.

  Int_t eventNumber = runLoader->GetEventNumber() ;

  AliPHOSGetter *gime = AliPHOSGetter::Instance();
  gime->Event(eventNumber, "DRTP") ; 

  TClonesArray *recParticles  = gime->RecParticles();
  Int_t nOfRecParticles = recParticles->GetEntries();

  esd->SetNumberOfPHOSClusters(nOfRecParticles) ; 
  esd->SetFirstPHOSCluster(esd->GetNumberOfCaloClusters()) ;
  
  AliDebug(2,Form("%d digits and %d rec. particles in event %d, option %s",gime->Digits()->GetEntries(),nOfRecParticles,eventNumber,GetOption()));

  for (Int_t recpart = 0 ; recpart < nOfRecParticles ; recpart++) {
    AliPHOSRecParticle * rp = dynamic_cast<AliPHOSRecParticle*>(recParticles->At(recpart));

    if(rp) {
    Float_t xyz[3];
    for (Int_t ixyz=0; ixyz<3; ixyz++) 
      xyz[ixyz] = rp->GetPos()[ixyz];

    AliDebug(2,Form("Global position xyz=(%f,%f,%f)",xyz[0],xyz[1],xyz[2]));
    
    AliPHOSTrackSegment *ts    = gime->TrackSegment(rp->GetPHOSTSIndex());
    AliPHOSEmcRecPoint  *emcRP = gime->EmcRecPoint(ts->GetEmcIndex());
    AliESDCaloCluster   *ec    = new AliESDCaloCluster() ; 

    Int_t  digitMult  = emcRP->GetDigitsMultiplicity();
    Int_t *digitsList = emcRP->GetDigitsList();
    UShort_t *amplList  = new UShort_t[digitMult];
    UShort_t *digiList  = new UShort_t[digitMult];

    // Convert Float_t* and Int_t* to UShort_t* to save memory
    for (Int_t iDigit=0; iDigit<digitMult; iDigit++) {
      AliPHOSDigit *digit = gime->Digit(digitsList[iDigit]);
      if(!digit) {
	AliFatal(Form("Digit not found at the expected position %d!",iDigit));
      }
      else {
	amplList[iDigit] = (UShort_t)digit->GetEnergy();
	digiList[iDigit] = (UShort_t)(digit->GetId());
      }
    }

    ec->SetGlobalPosition(xyz);                 //rec.point position in MARS
    ec->SetClusterEnergy(rp->Energy());         //total particle energy
    ec->SetClusterDisp(emcRP->GetDispersion()); //cluster dispersion
    ec->SetPid          (rp->GetPID()) ;        //array of particle identification
    ec->SetM02(emcRP->GetM2x()) ;               //second moment M2x
    ec->SetM20(emcRP->GetM2z()) ;               //second moment M2z
    ec->SetNExMax(emcRP->GetNExMax());          //number of local maxima
    ec->SetNumberOfDigits(digitMult);           //digit multiplicity
    ec->SetDigitAmplitude(amplList);            //digit energies
    ec->SetDigitIndex(digiList);                //abs id of the cell
    ec->SetEmcCpvDistance(-1);                  //not yet implemented
    ec->SetClusterChi2(-1);                     //not yet implemented
    ec->SetM11(-1) ;                            //not yet implemented

    // add the track to the esd object
    esd->AddCaloCluster(ec);
    delete ec;    

    }
  }


}

AliTracker* AliPHOSReconstructor::CreateTracker(AliRunLoader* runLoader) const
{
// creates the PHOS tracker
  if (!runLoader) return NULL; 
  return new AliPHOSTracker(runLoader);
}

