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
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliPHOSReconstructor.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSPIDv1.h"
#include "AliPHOSTracker.h"
#include "AliRawReader.h"
#include "AliPHOSTrigger.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSRecoParam.h"
#include "AliPHOSRecoParamEmc.h"
#include "AliPHOSRecoParamCpv.h"
#include "AliPHOSDigit.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSRecParticle.h"
#include "AliPHOSRawDecoder.h"
#include "AliPHOSRawDigiProducer.h"
#include "AliPHOSPulseGenerator.h"

ClassImp(AliPHOSReconstructor)

Bool_t AliPHOSReconstructor::fgDebug = kFALSE ; 
AliPHOSRecoParam* AliPHOSReconstructor::fgkRecoParamEmc =0;  // EMC rec. parameters
AliPHOSRecoParam* AliPHOSReconstructor::fgkRecoParamCpv =0;  // CPV rec. parameters

//____________________________________________________________________________
AliPHOSReconstructor::AliPHOSReconstructor() :
  fGeom(NULL)
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

  fGeom = AliPHOSGeometry::GetInstance("IHEP","");
}

//____________________________________________________________________________
  AliPHOSReconstructor::~AliPHOSReconstructor()
{
  // dtor
  delete fGeom;
} 

//____________________________________________________________________________
void AliPHOSReconstructor::Reconstruct(TTree* digitsTree, TTree* clustersTree) const
{
  // 'single-event' local reco method called by AliReconstruction; 
  // Only the clusterization is performed,; the rest of the reconstruction is done in FillESD because the track
  // segment maker needs access to the AliESDEvent object to retrieve the tracks reconstructed by 
  // the global tracking.

  AliPHOSClusterizerv1 clu(fGeom);
  clu.SetInput(digitsTree);
  clu.SetOutput(clustersTree);
  if ( Debug() ) 
    clu.Digits2Clusters("deb all") ; 
  else 
    clu.Digits2Clusters("") ;
}

//____________________________________________________________________________
void AliPHOSReconstructor::FillESD(TTree* digitsTree, TTree* clustersTree, 
				   AliESDEvent* esd) const
{
  // This method produces PHOS rec-particles,
  // then it creates AliESDtracks out of them and
  // write tracks to the ESD

  AliPHOSTrackSegmentMaker *tsm = new AliPHOSTrackSegmentMakerv1(fGeom);
  AliPHOSPID               *pid = new AliPHOSPIDv1              (fGeom);

  // do current event; the loop over events is done by AliReconstruction::Run()
  tsm->SetESD(esd) ; 
  tsm->SetInput(clustersTree);
  if ( Debug() ) 
    tsm->Clusters2TrackSegments("deb all") ;
  else 
    tsm->Clusters2TrackSegments("") ;
  
  pid->SetInput(clustersTree, tsm->GetTrackSegments()) ; 
  pid->SetESD(esd) ; 
  if ( Debug() ) 
    pid->TrackSegments2RecParticles("deb all") ;
  else 
    pid->TrackSegments2RecParticles("") ;

	
  // This function creates AliESDtracks from AliPHOSRecParticles
  //         and
  // writes them to the ESD

  TClonesArray *recParticles  = pid->GetRecParticles();
  Int_t nOfRecParticles = recParticles->GetEntries();
  
  esd->SetNumberOfPHOSClusters(nOfRecParticles) ; 
  esd->SetFirstPHOSCluster(esd->GetNumberOfCaloClusters()) ;

  AliDebug(2,Form("%d rec. particles, option %s",nOfRecParticles,GetOption()));


  //#########Calculate trigger and set trigger info###########
 
  AliPHOSTrigger tr ;
  //   tr.SetPatchSize(1);//create 4x4 patches
  tr.Trigger();
  
  Float_t maxAmp2x2  = tr.Get2x2MaxAmplitude();
  Float_t maxAmpnxn  = tr.GetnxnMaxAmplitude();
  Float_t ampOutOfPatch2x2  = tr.Get2x2AmpOutOfPatch() ;
  Float_t ampOutOfPatchnxn  = tr.GetnxnAmpOutOfPatch() ;

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
  fGeom->RelToAbsNumbering(iRelId2x2, iAbsId2x2);
  fGeom->RelToAbsNumbering(iRelIdnxn, iAbsIdnxn);
  fGeom->RelPosInAlice(iAbsId2x2, pos2x2);
  fGeom->RelPosInAlice(iAbsIdnxn, posnxn);

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
  
  // Read digits array
  TBranch *branch = digitsTree->GetBranch("PHOS");
  if (!branch) { 
    AliError("can't get the branch with the PHOS digits !");
    return;
  }
  TClonesArray *fDigitsArr    = new TClonesArray("AliPHOSDigit",100);
  branch->SetAddress(&fDigitsArr);
  branch->GetEntry(0);

  // Get the clusters array
  TBranch *emcbranch = clustersTree->GetBranch("PHOSEmcRP");
  if (!emcbranch) { 
    AliError("can't get the branch with the PHOS EMC clusters !");
    return;
  }

  TObjArray *fEmcRecPoints = new TObjArray(100) ;
  emcbranch->SetAddress(&fEmcRecPoints);
  emcbranch->GetEntry(0);

  //Fill CaloClusters 
  const Float_t kBigShort = std::numeric_limits<short int>::max() - 1;
  const Float_t nsec100   = 1e9*100.; // units of 0.01 ns
  const Float_t gev500    = 500.;     // units of GeV/500

  for (Int_t recpart = 0 ; recpart < nOfRecParticles ; recpart++) {
    AliPHOSRecParticle * rp = dynamic_cast<AliPHOSRecParticle*>(recParticles->At(recpart));
    if (Debug()) 
      rp->Print();
    // Get track segment and EMC rec.point associated with this rec.particle
    AliPHOSTrackSegment *ts    = static_cast<AliPHOSTrackSegment *>(tsm->GetTrackSegments()->At(rp->GetPHOSTSIndex()));

    AliPHOSEmcRecPoint  *emcRP = static_cast<AliPHOSEmcRecPoint *>(fEmcRecPoints->At(ts->GetEmcIndex()));
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
      AliPHOSDigit *digit = static_cast<AliPHOSDigit *>(fDigitsArr->At(digitsList[iDigit]));
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
    ec->SetPosition(xyz);                       //rec.point position in MARS
    ec->SetE(rp->Energy());                     //total particle energy
    ec->SetClusterDisp(emcRP->GetDispersion()); //cluster dispersion
    ec->SetPid(rp->GetPID()) ;                  //array of particle identification
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
    delete [] primList;
    delete [] amplList;
    delete [] timeList;
    delete [] digiList;    
  }
  fDigitsArr   ->Delete();
  delete fDigitsArr;
  fEmcRecPoints->Delete();
  delete fEmcRecPoints;
  delete tsm;
  delete pid;
}

//____________________________________________________________________________
AliTracker* AliPHOSReconstructor::CreateTracker() const
{
  // creates the PHOS tracker
  return new AliPHOSTracker();
}

//____________________________________________________________________________
void  AliPHOSReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{
  // Converts raw data to
  // PHOS digits
  // Works on a single-event basis

  rawReader->Reset() ; 

  AliPHOSRawDecoder dc(rawReader);
  TString option = GetOption();
  if (option.Contains("OldRCUFormat"))
    dc.SetOldRCUFormat(kTRUE);
  else
    dc.SetOldRCUFormat(kFALSE);
  
  dc.SubtractPedestals(fgkRecoParamEmc->SubtractPedestals());
  
  TClonesArray *digits = new TClonesArray("AliPHOSDigit",1);
  digits->SetName("DIGITS");
  Int_t bufsize = 32000;
  digitsTree->Branch("PHOS", &digits, bufsize);

  AliPHOSRawDigiProducer pr;
  pr.MakeDigits(digits,&dc);

  //ADC counts -> GeV
  for(Int_t i=0; i<digits->GetEntries(); i++) {
    AliPHOSDigit* digit = (AliPHOSDigit*)digits->At(i);
    digit->SetEnergy(digit->GetEnergy()/AliPHOSPulseGenerator::GeV2ADC());
  }
  
  // Clean up digits below the noise threshold
  // Assuming the digit noise to be 4 MeV, we suppress digits within
  // 3-sigma of the noise.
  // This parameter should be passed via AliPHOSRecoParamEmc later

  const Double_t emcDigitThreshold = 0.012;
  for(Int_t i=0; i<digits->GetEntries(); i++) {
    AliPHOSDigit* digit = (AliPHOSDigit*)digits->At(i);
    if(digit->GetEnergy() < emcDigitThreshold)
      digits->RemoveAt(i) ;
  }
  digits->Compress() ;  

  //!!!!for debug!!!
  Int_t modMax=-111;
  Int_t colMax=-111;
  Int_t rowMax=-111;
  Float_t eMax=-333;
  //!!!for debug!!!

  Int_t relId[4];
  for(Int_t iDigit=0; iDigit<digits->GetEntries(); iDigit++) {
    AliPHOSDigit* digit = (AliPHOSDigit*)digits->At(iDigit);
    if(digit->GetEnergy()>eMax) {
      fGeom->AbsToRelNumbering(digit->GetId(),relId);
      eMax=digit->GetEnergy();
      modMax=relId[0];
      rowMax=relId[2];
      colMax=relId[3];
    }
  }

  AliDebug(1,Form("Digit with max. energy:  modMax %d colMax %d rowMax %d  eMax %f\n\n",
		  modMax,colMax,rowMax,eMax));

  digitsTree->Fill();
  digits->Delete();
  delete digits;
}
