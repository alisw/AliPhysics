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
#include "AliPHOSTrigger.h"
#include "AliPHOSGeometry.h"

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


  //#########Calculate trigger and set trigger info###########
 
  AliPHOSTrigger tr ;
  //   tr.SetPatchSize(1);//create 4x4 patches
  tr.Trigger();
  
  Float_t maxAmp2x2  = tr.Get2x2MaxAmplitude();
  Float_t maxAmpnxn  = tr.GetnxnMaxAmplitude();
  Float_t ampOutOfPatch2x2  = tr.Get2x2AmpOutOfPatch() ;
  Float_t ampOutOfPatchnxn  = tr.GetnxnAmpOutOfPatch() ;

  AliPHOSGeometry * geom = gime->PHOSGeometry();

  Int_t iSM2x2      = tr.Get2x2SuperModule();
  Int_t iSMnxn      = tr.GetnxnSuperModule();
  Int_t iCrystalPhi2x2 = tr.Get2x2CrystalPhi();
  Int_t iCrystalPhinxn = tr.GetnxnCrystalPhi();
  Int_t iCrystalEta2x2 = tr.Get2x2CrystalEta();
  Int_t iCrystalEtanxn = tr.GetnxnCrystalEta();

  AliDebug(2, Form("Trigger 2x2 max amp %f, out amp %f, SM %d, iphi %d ieta %d",  maxAmp2x2, ampOutOfPatch2x2, iSM2x2,iCrystalPhi2x2, iCrystalEta2x2));
  AliDebug(2, Form("Trigger 4x4 max amp %f , out amp %f, SM %d, iphi %d, ieta %d",  maxAmpnxn, ampOutOfPatchnxn, iSMnxn,iCrystalPhinxn, iCrystalEtanxn));

  Int_t iRelId2x2 []= {iSM2x2+1,0,iCrystalPhi2x2,iCrystalEta2x2};// PHOS modules in order to calculate AbsId need to be 1-5 not 0-4 as returns trigger.
  Int_t iAbsId2x2 =-1;
  Int_t iRelIdnxn []= {iSMnxn+1,0,iCrystalPhinxn,iCrystalEtanxn};// PHOS modules in order to calculate AbsId need to be 1-5 not 0-4 as returns trigger.
  Int_t iAbsIdnxn =-1;
  TVector3    pos2x2(-1,-1,-1);
  TVector3    posnxn(-1,-1,-1);
  geom->RelToAbsNumbering(iRelId2x2, iAbsId2x2);
  geom->RelToAbsNumbering(iRelIdnxn, iAbsIdnxn);
  geom->RelPosInAlice(iAbsId2x2, pos2x2);
  geom->RelPosInAlice(iAbsIdnxn, posnxn);

  TArrayF triggerPosition(6);
  triggerPosition[0] = pos2x2(0) ;   
  triggerPosition[1] = pos2x2(1) ;   
  triggerPosition[2] = pos2x2(2) ;  
  triggerPosition[3] = posnxn(0) ;   
  triggerPosition[4] = posnxn(1) ;   
  triggerPosition[5] = posnxn(2) ;  

  TArrayF triggerAmplitudes(4);
  triggerAmplitudes[0] = maxAmp2x2 ;   
  triggerAmplitudes[1] = ampOutOfPatch2x2 ;    
  triggerAmplitudes[2] = maxAmpnxn ;   
  triggerAmplitudes[3] = ampOutOfPatchnxn ;   

  //esd->SetPHOSTriggerCells(triggerPosition);
  esd->AddPHOSTriggerPosition(triggerPosition);
  esd->AddPHOSTriggerAmplitudes(triggerAmplitudes);
  
  //######################################

  //Fill CaloClusters 
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

    ec->SetPHOS(kTRUE);
    ec->SetGlobalPosition(xyz);                 //rec.point position in MARS
    ec->SetClusterEnergy(rp->Energy());         //total particle energy
    ec->SetClusterDisp(emcRP->GetDispersion()); //cluster dispersion
    ec->SetPid          (rp->GetPID()) ;        //array of particle identification
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

    //Primaries
    ec->SetPrimaryIndex(rp->GetPrimaryIndex());
    Int_t  primMult  = 0;
    Int_t *primInts =  emcRP->GetPrimaries(primMult);
    ec->SetNumberOfPrimaries(primMult);           //primary multiplicity
    UShort_t *primList = new UShort_t[primMult];
    for (Int_t ipr=0; ipr<primMult; ipr++) 
      primList[ipr] = (UShort_t)(primInts[ipr]);	 
    ec->SetListOfPrimaries(primList);                  //primary List for a cluster
    
    // add the track to the esd object
    esd->AddCaloCluster(ec);
    delete ec;    
  }  
}

//____________________________________________________________________________
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

    ec->SetPHOS(kTRUE);
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

