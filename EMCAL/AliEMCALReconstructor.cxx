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

#include "AliESD.h"
#include "AliRunLoader.h"
#include "AliEMCALLoader.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALRecPoint.h"
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
  clu.SetEventRange(0, -1) ; // do all the events
  if ( Debug() ) 
    clu.ExecuteTask("deb all") ; 
  else 
    clu.ExecuteTask("pseudo") ;  
 
}

//____________________________________________________________________________
void AliEMCALReconstructor::Reconstruct(AliRunLoader* runLoader, AliRawReader* rawreader) const 
{
  // method called by AliReconstruction; 
  // Only the clusterization is performed,; the rest of the reconstruction is done in FillESD because the track
  // segment maker needs access to the AliESD object to retrieve the tracks reconstructed by 
  // the global tracking.
  // Here we reconstruct from Raw Data
  
  rawreader->Reset() ; 
  TString headerFile(runLoader->GetFileName()) ; 
  TString branchName(runLoader->GetEventFolder()->GetName()) ;  

  AliEMCALClusterizerv1 clu(headerFile, branchName);
  clu.SetEventRange(0, -1) ; // do all the events
  if ( Debug() ) 
    clu.ExecuteTask("deb pseudo all") ; 
  else 
    clu.ExecuteTask("pseudo") ;  

}

//____________________________________________________________________________
void AliEMCALReconstructor::FillESD(AliRunLoader* runLoader, AliESD* esd) const
{
  // Called by AliReconstruct after Reconstruct() and global tracking and vertxing 

  Int_t eventNumber = runLoader->GetEventNumber() ;

  TString headerFile(runLoader->GetFileName()) ; 
  TString branchName(runLoader->GetEventFolder()->GetName()) ;  
  // Creates AliESDCaloCluster from AliEMCALRecPoints 
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));
  rl->LoadRecPoints();
  rl->LoadKinematics(); // To get the primary label
  rl->LoadDigits();     // To get the primary label
  rl->LoadHits();       // To get the primary label
  rl->GetEvent(eventNumber);
  TObjArray *clusters = emcalLoader->RecPoints();
  Int_t nClusters = clusters->GetEntries();
  esd->SetNumberOfEMCALClusters(nClusters) ; 
  esd->SetFirstEMCALCluster(esd->GetNumberOfCaloClusters()) ; 

  for (Int_t iClust = 0 ; iClust < nClusters ; iClust++) {
    const AliEMCALRecPoint * clust = emcalLoader->RecPoint(iClust);

    if (Debug()) clust->Print();
    AliESDCaloCluster * ec = new AliESDCaloCluster() ; 
    // fills the ESDCaloCluster
    Float_t xyz[3];
    TVector3 gpos;
    clust->GetGlobalPosition(gpos);
    for (Int_t ixyz=0; ixyz<3; ixyz++) 
      xyz[ixyz] = gpos[ixyz];

    Int_t digitMult = clust->GetMultiplicity();
    UShort_t *amplList = new UShort_t[digitMult];
    UShort_t *timeList = new UShort_t[digitMult];
    UShort_t *digiList = new UShort_t[digitMult];
    Float_t *amplFloat = clust->GetEnergiesList();
    Float_t *timeFloat = clust->GetTimeList();
    Int_t   *digitInts = clust->GetAbsId();
    Float_t *elipAxis  = new Float_t();
    clust->GetElipsAxis(elipAxis);

   // Convert Float_t* and Int_t* to UShort_t* to save memory
    Int_t newdigitMult = digitMult ;
    for (Int_t iDigit=0; iDigit<digitMult; iDigit++) {
      if(timeFloat[iDigit] < 65536/1e9*100){
	amplList[iDigit] = (UShort_t)(amplFloat[iDigit]*500);
	timeList[iDigit] = (UShort_t)(timeFloat[iDigit]*1e9*100); //Time in units of 100 ns = 0.1 ps
	digiList[iDigit] = (UShort_t)(digitInts[iDigit]);
      }
      else
	newdigitMult = 	newdigitMult - 1 ;
    } 
 
    ec->SetClusterType(clust->GetClusterType());
    ec->SetGlobalPosition(xyz);
    ec->SetClusterEnergy(clust->GetEnergy());
    ec->SetDigitAmplitude(amplList); //energies
    ec->SetDigitTime(timeList);      //times
    ec->SetDigitIndex(digiList);     //indices
    ec->SetNumberOfDigits(newdigitMult);
    if(clust->GetClusterType()== AliESDCaloCluster::kClusterv1){
      ec->SetClusterDisp(clust->GetDispersion());
      ec->SetClusterChi2(-1); //not yet implemented
      ec->SetM02(elipAxis[0]*elipAxis[0]) ;
      ec->SetM20(elipAxis[1]*elipAxis[1]) ;
      ec->SetM11(-1) ;        //not yet implemented
      ec->SetPrimaryIndex(clust->GetPrimaryIndex());
    }
    // add the cluster to the esd object
    esd->AddCaloCluster(ec);
    delete ec;
  }
}
