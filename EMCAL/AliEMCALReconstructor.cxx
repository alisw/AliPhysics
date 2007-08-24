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
#include "AliEMCALReconstructor.h"

#include "AliRun.h"
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDtrack.h"
#include "AliRunLoader.h"
#include "AliEMCALLoader.h"
#include "AliEMCALRawUtils.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALPID.h"
#include "AliEMCALTrigger.h"
#include "AliRawReader.h"

ClassImp(AliEMCALReconstructor)

//____________________________________________________________________________
AliEMCALReconstructor::AliEMCALReconstructor() 
  : fDebug(kFALSE) 
{
  // ctor
} 

//____________________________________________________________________________
AliEMCALReconstructor::AliEMCALReconstructor(const AliEMCALReconstructor & rec)
  : AliReconstructor(rec),
    fDebug(rec.fDebug)
{
  //copy ctor
}

//____________________________________________________________________________
AliEMCALReconstructor::~AliEMCALReconstructor()
{
  // dtor
} 

//____________________________________________________________________________
void AliEMCALReconstructor::Reconstruct(AliRunLoader* runLoader) const 
{
  // method called by AliReconstruction; 
  // Only the clusterization is performed,; the rest of the reconstruction is done in FillESD because the track
  // segment maker needs access to the AliESD object to retrieve the tracks reconstructed by 
  // the global tracking.
 
  TString headerFile(runLoader->GetFileName()) ; 
  TString branchName(runLoader->GetEventFolder()->GetName() ) ;  
  
  AliEMCALClusterizerv1 clu(headerFile, branchName);
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(runLoader->GetDetectorLoader("EMCAL"));    

  Int_t nEvents   = runLoader->GetNumberOfEvents();
  runLoader->LoadDigits("EMCAL");
  for (Int_t ievent = 0; ievent < nEvents; ievent++) {
    runLoader->GetEvent(ievent);
    if ( Debug() ) 
      clu.ExecuteTask("deb all") ; 
    else 
      clu.ExecuteTask("pseudo") ;  
  }   
  // Unload the Digits and RecPoints
  emcalLoader->UnloadDigits() ; 
  emcalLoader->UnloadRecPoints() ; 
}

//____________________________________________________________________________
void AliEMCALReconstructor::Reconstruct(AliRunLoader* runLoader, AliRawReader* rawReader) const 

{
  // Reconstruction loop for Raw Data processing
  //
  // Only the clusterization is performed
  // Track-cluster matching is done in FillESD because the track
  // segment maker needs access to the AliESD object to retrieve 
  // the tracks reconstructed by the global tracking.
  
  TString headerFile(runLoader->GetFileName()) ; 
  TString branchName(runLoader->GetEventFolder()->GetName()) ;  

  static AliEMCALRawUtils rawUtils;
  AliEMCALClusterizerv1 clu(headerFile, branchName);
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(runLoader->GetDetectorLoader("EMCAL"));    
 
  Int_t iEvent = 0;
  rawReader->Reset() ; 
  while (rawReader->NextEvent()) {  
    runLoader->GetEvent(iEvent++);

    rawUtils.Raw2Digits(rawReader);

    if ( Debug() ) 
      clu.ExecuteTask("deb pseudo all") ; 
    else 
      clu.ExecuteTask("pseudo") ;  
    // Unload the RecPoints
    emcalLoader->UnloadRecPoints() ; 
  }
}

//____________________________________________________________________________
void AliEMCALReconstructor::FillESD(AliRunLoader* runLoader, AliESDEvent* esd) const
{
  // Called by AliReconstruct after Reconstruct() and global tracking and vertexing 
  const double timeScale = 1.e+11; // transition constant from sec to 0.01 ns 

  Int_t eventNumber = runLoader->GetEventNumber() ;

  // Creates AliESDCaloCluster from AliEMCALRecPoints 
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(runLoader->GetDetectorLoader("EMCAL"));
  runLoader->LoadRecPoints("EMCAL");
  runLoader->GetEvent(eventNumber);
  TObjArray *clusters = emcalLoader->RecPoints();
  Int_t nClusters = clusters->GetEntries(), nClustersNew=0;
  AliDebug(1,Form("Event %d: %d clusters",eventNumber,nClusters));
  //  Int_t nRP=0, nPC=0; // in input
  esd->SetFirstEMCALCluster(esd->GetNumberOfCaloClusters()); // Put after Phos clusters 

  //######################################################
  //#########Calculate trigger and set trigger info###########
  //######################################################
 
  AliEMCALTrigger tr ;
  //   tr.SetPatchSize(1);//create 4x4 patches
  tr.Trigger();
  
  Float_t maxAmp2x2  = tr.Get2x2MaxAmplitude();
  Float_t maxAmpnxn  = tr.GetnxnMaxAmplitude();
  Float_t ampOutOfPatch2x2  = tr.Get2x2AmpOutOfPatch() ;
  Float_t ampOutOfPatchnxn  = tr.GetnxnAmpOutOfPatch() ;

  AliEMCALGeometry * geom = 0;
  if (runLoader->GetAliRun() && runLoader->GetAliRun()->GetDetector("EMCAL"))
    geom = dynamic_cast<AliEMCAL*>(runLoader->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
  if (geom == 0) 
    geom = AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaulGeometryName());

  Int_t iSM2x2      = tr.Get2x2SuperModule();
  Int_t iSMnxn      = tr.GetnxnSuperModule();
  Int_t iCellPhi2x2 = tr.Get2x2CellPhi();
  Int_t iCellPhinxn = tr.GetnxnCellPhi();
  Int_t iCellEta2x2 = tr.Get2x2CellEta();
  Int_t iCellEtanxn = tr.GetnxnCellEta();

  AliDebug(2, Form("Trigger 2x2 max amp %f, out amp %f, SM %d, iphi %d ieta %d",  maxAmp2x2, ampOutOfPatch2x2, iSM2x2,iCellPhi2x2, iCellEta2x2));
  AliDebug(2, Form("Trigger 4x4 max amp %f , out amp %f, SM %d, iphi %d, ieta %d",  maxAmpnxn, ampOutOfPatchnxn, iSMnxn,iCellPhinxn, iCellEtanxn));

  TVector3    pos2x2(-1,-1,-1);
  TVector3    posnxn(-1,-1,-1);

  Int_t iAbsId2x2 = geom->GetAbsCellIdFromCellIndexes( iSM2x2, iCellPhi2x2, iCellEta2x2) ;
  Int_t iAbsIdnxn = geom->GetAbsCellIdFromCellIndexes( iSMnxn, iCellPhinxn, iCellEtanxn) ;
  geom->GetGlobal(iAbsId2x2, pos2x2);
  geom->GetGlobal(iAbsIdnxn, posnxn);
  
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

  esd->AddEMCALTriggerPosition(triggerPosition);
  esd->AddEMCALTriggerAmplitudes(triggerAmplitudes);
  
  //######################################################
  //#######################TRACK MATCHING###############
  //######################################################
  //Fill list of integers, each one is index of track to which the cluster belongs.

  // step 1 - initialize array of matched track indexes
  Int_t *matchedTrack = new Int_t[nClusters];
  for (Int_t iclus = 0; iclus < nClusters; iclus++)
    matchedTrack[iclus] = -1;  // neg. index --> no matched track
  
  // step 2, change the flag for all matched clusters found in tracks
  Int_t iemcalMatch = -1;
  Int_t endtpc = esd->GetNumberOfTracks();
  for (Int_t itrack = 0; itrack < endtpc; itrack++) {
    AliESDtrack * track = esd->GetTrack(itrack) ; // retrieve track
    iemcalMatch = track->GetEMCALcluster();
    if(iemcalMatch >= 0) matchedTrack[iemcalMatch] = itrack;
  } 
  
  //########################################
  //##############Fill CaloClusters#############
  //########################################


  for (Int_t iClust = 0 ; iClust < nClusters ; iClust++) {
    const AliEMCALRecPoint * clust = emcalLoader->RecPoint(iClust);
    //if(clust->GetClusterType()== AliESDCaloCluster::kClusterv1) nRP++; else nPC++;
    if (Debug()) clust->Print();
    // Get information from EMCAL reconstruction points
    Float_t xyz[3];
    TVector3 gpos;
    clust->GetGlobalPosition(gpos);
    for (Int_t ixyz=0; ixyz<3; ixyz++) 
      xyz[ixyz] = gpos[ixyz];
    
    Int_t digitMult = clust->GetMultiplicity();
    Short_t *amplList = new Short_t[digitMult];
    Short_t *timeList = new Short_t[digitMult];
    Short_t *digiList = new Short_t[digitMult];
    Float_t *amplFloat = clust->GetEnergiesList();
    Float_t *timeFloat = clust->GetTimeList();
    Int_t   *digitInts = clust->GetAbsId();
    Float_t elipAxis[2];
    clust->GetElipsAxis(elipAxis);
    
    // Convert Float_t* and Int_t* to Short_t* to save memory
    // Problem : we should recalculate a cluster characteristics when discard digit(s)
    Int_t newdigitMult = 0; 
    for (Int_t iDigit=0; iDigit<digitMult; iDigit++) {
      if (amplFloat[iDigit] > 0) {
	amplList[newdigitMult] = (UShort_t)(amplFloat[iDigit]*500);
        // Time in units of 0.01 ns = 10 ps
        if(timeFloat[iDigit] < 65536./timeScale) 
	  timeList[newdigitMult] = (UShort_t)(timeFloat[iDigit]*timeScale);
        else
          timeList[newdigitMult] = 65535;
	digiList[newdigitMult] = (UShort_t)(digitInts[iDigit]);
        newdigitMult++;
      }
      else if (clust->GetClusterType() != AliESDCaloCluster::kPseudoCluster)
        Warning("FillESD()","Negative or 0 digit amplitude in cluster");
    }
    
    if(newdigitMult > 0) { // accept cluster if it has some digit
      nClustersNew++;
      if(newdigitMult != digitMult) { // some digits were deleted
        Short_t *amplListNew = new Short_t[newdigitMult];
        Short_t *timeListNew = new Short_t[newdigitMult];
        Short_t *digiListNew = new Short_t[newdigitMult];
        for (Int_t iDigit=0; iDigit<newdigitMult; iDigit++) {
          amplListNew[iDigit] = amplList[iDigit];
          timeListNew[iDigit] = timeList[iDigit];
          digiListNew[iDigit] = digiList[iDigit];
        }
	
        delete [] amplList;
        delete [] timeList;
        delete [] digiList;
	
        amplList = amplListNew;
        timeList = timeListNew;
        digiList = digiListNew;
      }
      
      //Primaries
      Int_t  parentMult  = 0;
      Int_t *parentInts =  clust->GetParents(parentMult);
      Short_t *parentList = new Short_t[parentMult];
      for (Int_t ipr=0; ipr<parentMult; ipr++) 
	parentList[ipr] = (Short_t)(parentInts[ipr]);	 
      
    
      // fills the ESDCaloCluster
      AliESDCaloCluster * ec = new AliESDCaloCluster() ; 
      ec->SetEMCAL(kTRUE);
      ec->SetClusterType(clust->GetClusterType());
      ec->SetPosition(xyz);
      ec->SetE(clust->GetEnergy());
      TArrayS arrayAmpList(newdigitMult,amplList);
      TArrayS arrayTimeList(newdigitMult,timeList);
      TArrayS arrayIndexList(newdigitMult,digiList);
      ec->AddDigitAmplitude(arrayAmpList);
      ec->AddDigitTime(arrayTimeList);
      ec->AddDigitIndex(arrayIndexList);
    
      if(clust->GetClusterType()== AliESDCaloCluster::kClusterv1){

        ec->SetClusterDisp(clust->GetDispersion());
        ec->SetClusterChi2(-1); //not yet implemented
        ec->SetM02(elipAxis[0]*elipAxis[0]) ;
        ec->SetM20(elipAxis[1]*elipAxis[1]) ;
        ec->SetM11(-1) ;        //not yet implemented
	
       TArrayS arrayTrackMatched(1);// Only one track, temporal solution.
       arrayTrackMatched[0]= (Short_t)(matchedTrack[iClust]);
       ec->AddTracksMatched(arrayTrackMatched);
	
       TArrayS arrayParents(parentMult,parentList);
       ec->AddLabels(arrayParents);
      } 
      
      // add the cluster to the esd object
      esd->AddCaloCluster(ec);
      delete ec;
      delete [] parentList;
    } else { // no new ESD cluster
      
    }
    delete [] amplList;
    delete [] timeList;
    delete [] digiList;

  } // cycle on clusters

  delete [] matchedTrack;

  esd->SetNumberOfEMCALClusters(nClustersNew);
  //if(nClustersNew != nClusters) 
  //printf(" ##### nClusters %i -> new %i ##### \n", nClusters, nClustersNew );
  
  //Fill ESDCaloCluster with PID weights
  AliEMCALPID *pid = new AliEMCALPID;
  //pid->SetPrintInfo(kTRUE);
  pid->SetReconstructor(kTRUE);
  pid->RunPID(esd);
  delete pid;
}


