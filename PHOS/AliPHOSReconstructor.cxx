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
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliPHOSReconstructor.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSPIDv1.h"
#include "AliPHOSGetter.h"
#include "AliPHOSTracker.h"
#include "AliRawReader.h"
#include "AliPHOSTrigger.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSRecoParamEmc.h"
#include "AliPHOSRecoParamCpv.h"

ClassImp(AliPHOSReconstructor)

Bool_t AliPHOSReconstructor::fgDebug = kFALSE ; 
AliPHOSRecoParam* AliPHOSReconstructor::fgkRecoParamEmc =0;  // EMC rec. parameters
AliPHOSRecoParam* AliPHOSReconstructor::fgkRecoParamCpv =0;  // CPV rec. parameters

//____________________________________________________________________________
  AliPHOSReconstructor::AliPHOSReconstructor() 
{
  // ctor

  if (!fgkRecoParamEmc) {
    AliWarning("The Reconstruction parameters for EMC nonitialized - Used default one");
    fgkRecoParamEmc = AliPHOSRecoParamEmc::GetEmcDefaultParameters();
  }

  if (!fgkRecoParamCpv) {
    AliWarning("The Reconstruction parameters for CPV nonitialized - Used default one");
    fgkRecoParamCpv = AliPHOSRecoParamCpv::GetCpvDefaultParameters();
  }

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
  // segment maker needs access to the AliESDEvent object to retrieve the tracks reconstructed by 
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
  // segment maker needs access to the AliESDEvent object to retrieve the tracks reconstructed by 
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
void AliPHOSReconstructor::FillESD(AliRunLoader* runLoader, AliESDEvent* esd) const
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

  Int_t maxClu = esd->GetNumberOfPHOSClusters() ; 

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
  const Float_t kBigShort = std::numeric_limits<short int>::max() - 1;
  const Float_t nsec100   = 1e9*100.; // units of 0.01 ns
  const Float_t gev500    = 500.;     // units of GeV/500

  for (Int_t recpart = 0 ; recpart < nOfRecParticles ; recpart++) {
    AliPHOSRecParticle * rp = dynamic_cast<AliPHOSRecParticle*>(recParticles->At(recpart));
    if (Debug()) 
      rp->Print();
    // Get track segment and EMC rec.point associated with this rec.particle
    AliPHOSTrackSegment *ts    = gime->TrackSegment(rp->GetPHOSTSIndex());
    AliPHOSEmcRecPoint  *emcRP = gime->EmcRecPoint(ts->GetEmcIndex());
    AliESDCaloCluster   *ec    = new AliESDCaloCluster() ; 
        
    Float_t xyz[3];
    for (Int_t ixyz=0; ixyz<3; ixyz++) 
      xyz[ixyz] = rp->GetPos()[ixyz];
    
    AliDebug(2,Form("Global position xyz=(%f,%f,%f)",xyz[0],xyz[1],xyz[2]));
    
    //Create digits lists
    Int_t  digitMult  = emcRP->GetDigitsMultiplicity();
    Int_t *digitsList = emcRP->GetDigitsList();
    Short_t *amplList  = new Short_t[digitMult];
    Short_t *timeList  = new Short_t[digitMult];
    Short_t *digiList  = new Short_t[digitMult];

    // Convert Float_t* and Int_t* to Short_t* to save memory
    for (Int_t iDigit=0; iDigit<digitMult; iDigit++) {
      AliPHOSDigit *digit = gime->Digit(digitsList[iDigit]);
      amplList[iDigit] =
	(Short_t)(TMath::Min(digit->GetEnergy()*gev500,kBigShort)); // Energy in units of GeV/500
      timeList[iDigit] =
	(Short_t)(TMath::Min(digit->GetTime()*nsec100,kBigShort)); // time in units of 0.01 ns
      digiList[iDigit] = (Short_t)(digit->GetId());
    }
    
    //Primaries
    Int_t  primMult  = 0;
    Int_t *primInts =  emcRP->GetPrimaries(primMult);
    Short_t *primList = new Short_t[primMult];
    for (Int_t ipr=0; ipr<primMult; ipr++) 
      primList[ipr] = (Short_t)(primInts[ipr]);	 
    
    // fills the ESDCaloCluster
 
    ec->SetPHOS(kTRUE);
    ec->SetPosition(xyz);                 //rec.point position in MARS
    ec->SetE(rp->Energy());         //total particle energy
    ec->SetClusterDisp(emcRP->GetDispersion()); //cluster dispersion
    ec->SetPid          (rp->GetPID()) ;        //array of particle identification
    ec->SetM02(emcRP->GetM2x()) ;               //second moment M2x
    ec->SetM20(emcRP->GetM2z()) ;               //second moment M2z
    ec->SetNExMax(emcRP->GetNExMax());          //number of local maxima
    ec->SetEmcCpvDistance(-1);                  //not yet implemented
    ec->SetClusterChi2(-1);                     //not yet implemented
    ec->SetM11(-1) ;                            //not yet implemented
 
    //Digits Lists
    TArrayS arrayAmpList(digitMult,amplList);
    TArrayS arrayTimeList(digitMult,timeList);
    TArrayS arrayIndexList(digitMult,digiList);
    ec->AddDigitAmplitude(arrayAmpList);
    ec->AddDigitTime(arrayTimeList);
    ec->AddDigitIndex(arrayIndexList);

    //Distance to the nearest bad crystal
    ec->SetDistanceToBadChannel(emcRP->GetDistanceToBadCrystal()); 
  
    //Array of MC indeces
    TArrayS arrayPrim(primMult,primList);
    ec->AddLabels(arrayPrim);

    //Array of tracks uncomment when available in future
    //TArrayS arrayTrackMatched(1);// Only one track, temporal solution.
    //arrayTrackMatched[0]= (Short_t)(matchedTrack[iClust]);
    //ec->AddTracksMatched(arrayTrackMatched);
    
    // add the track to the esd object
    esd->AddCaloCluster(ec);
    delete ec;    
  }  
}

//____________________________________________________________________________
void AliPHOSReconstructor::FillESD(AliRunLoader* runLoader,
				   AliRawReader* rawReader, AliESDEvent* esd) const
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
    Short_t *amplList  = new Short_t[digitMult];
    Short_t *digiList  = new Short_t[digitMult];

    // Convert Float_t* and Int_t* to UShort_t* to save memory
    for (Int_t iDigit=0; iDigit<digitMult; iDigit++) {
      AliPHOSDigit *digit = gime->Digit(digitsList[iDigit]);
      if(!digit) {
	AliFatal(Form("Digit not found at the expected position %d!",iDigit));
      }
      else {
	amplList[iDigit] = (Short_t)digit->GetEnergy();
	digiList[iDigit] = (Short_t)(digit->GetId());
	//timeList[iDigit] = (Short_t)(digit->GetTime());
      }
    }

    ec->SetPHOS(kTRUE);
    ec->SetPosition(xyz);                 //rec.point position in MARS
    ec->SetE(rp->Energy());         //total particle energy
    ec->SetClusterDisp(emcRP->GetDispersion()); //cluster dispersion
    ec->SetPid          (rp->GetPID()) ;        //array of particle identification
    ec->SetM02(emcRP->GetM2x()) ;               //second moment M2x
    ec->SetM20(emcRP->GetM2z()) ;               //second moment M2z
    ec->SetNExMax(emcRP->GetNExMax());          //number of local maxima

    ec->SetEmcCpvDistance(-1);                  //not yet implemented
    ec->SetClusterChi2(-1);                     //not yet implemented
    ec->SetM11(-1) ;                            //not yet implemented
 //    TArrayS arrayAmpList(digitMult,amplList);
//     TArrayS arrayTimeList(digitMult,timeList);
//     TArrayS arrayIndexList(digitMult,digiList);
//     ec->AddDigitAmplitude(arrayAmpList);
//     ec->AddDigitTime(arrayTimeList);
//     ec->AddDigitIndex(arrayIndexList);
    //Distance to the nearest bad crystal
    ec->SetDistanceToBadChannel(emcRP->GetDistanceToBadCrystal()); 
    
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

