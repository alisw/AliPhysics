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

//-----------------------------------------------------------------------------
/// \class AliMUONReconstructor
///
/// Implementation of AliReconstructor for MUON subsystem.
///
/// The behavior of the MUON reconstruction can be changed, besides
/// the usual methods found in AliReconstruction (e.g. to disable tracking)
/// by using AliReconstruction::SetOption("MUON",options)
/// where options should be a space separated string.
///
/// Valid options are :
///
/// SAVEDIGITS : if you want to save in the TreeD the *calibrated* digits
///     that are used for the clustering
///
/// SIMPLEFIT : use the AliMUONClusterFinderSimpleFit clusterizer
///
/// AZ : use the AliMUONClusterFinderAZ clusterizer (default)
///
/// MLEM : another implementation of AZ, where preclustering is external
/// MLEMV3 : MLEM with preclustering=PRECLUSTERV2
/// MLEMV3 : MLEM with preclustering=PRECLUSTERV3
///
/// PRECLUSTER : use only AliMUONPreClusterFinder. Only for debug as
/// the produced clusters do not have a position, hence the tracking will not
/// work
/// PRECLUSTERV2 : another version of the preclustering
/// PRECLUSTERV3 : yet another version of the preclustering
///
/// COG : use AliMUONClusterFinderCOG clusterizer. Not really a production
/// option either, as center-of-gravity is generally not a good estimate
/// of the cluster position...
///
/// NOCLUSTERING : bypass completely the clustering stage
///
/// NOSTATUSMAP : disable the computation and usage of the pad status map. Only
/// for debug !
///
/// NOLOCALRECONSTRUCTION : for debug, to disable local reconstruction (and hence
/// "recover" old behavior)
///
/// TRIGGERDISABLE : disable the treatment of MUON trigger
///
/// DIGITSTOREV1 : use the V1 implementation of the digitstore 
/// DIGITSTOREV2R : use the V2R implementation of the digitstore 
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONReconstructor.h"

#include "AliCDBManager.h"
#include "AliLog.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONClusterFinderCOG.h"
#include "AliMUONClusterFinderMLEM.h"
#include "AliMUONClusterFinderSimpleFit.h"
#include "AliMUONClusterFinderAZ.h"
#include "AliMUONClusterReconstructor.h"
#include "AliMUONClusterStoreV1.h"
#include "AliMUONConstants.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONDigitMaker.h"
#include "AliMUONDigitStoreV1.h"
#include "AliMUONDigitStoreV2R.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONPreClusterFinder.h"
#include "AliMUONPreClusterFinderV2.h"
#include "AliMUONPreClusterFinderV3.h"
#include "AliMUONTracker.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONTriggerChamberEff.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONTriggerStoreV1.h"
#include "AliMUONVClusterFinder.h"
#include "AliMpCDB.h"
#include "AliRawReader.h"
#include "AliCodeTimer.h"
#include <Riostream.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TTree.h>
//#include "AliCodeTimer.h"
/// \cond CLASSIMP
ClassImp(AliMUONReconstructor)
/// \endcond 

//_____________________________________________________________________________
AliMUONReconstructor::AliMUONReconstructor() : 
AliReconstructor(),
fCrateManager(0x0),
fDigitMaker(0x0),
fTransformer(new AliMUONGeometryTransformer()),
fDigitStore(0x0),
fTriggerCircuit(0x0),
fCalibrationData(0x0),
fDigitCalibrator(0x0),
fClusterReconstructor(0x0),
fClusterStore(0x0),
fTriggerStore(0x0),
fTrackStore(0x0),
fTrigChamberEff(0x0)
{
  /// normal ctor

  // Load mapping
  if ( ! AliMpCDB::LoadDDLStore() ) {
    AliFatal("Could not access mapping from OCDB !");
  }
  
  // Load geometry data
  fTransformer->LoadGeometryData();
  
  
}

//_____________________________________________________________________________
AliMUONReconstructor::~AliMUONReconstructor()
{
  /// dtor
  delete fDigitMaker;
  delete fDigitStore;
  delete fTransformer;
  delete fCrateManager;
  delete fTriggerCircuit;
  delete fCalibrationData;
  delete fDigitCalibrator;
  delete fClusterReconstructor;
  delete fClusterStore;
  delete fTriggerStore;
  delete fTrackStore;
  delete fTrigChamberEff;
}

//_____________________________________________________________________________
void
AliMUONReconstructor::Calibrate(AliMUONVDigitStore& digitStore) const
{
  /// Calibrate the digitStore
  if (!fDigitCalibrator)
  {
    CreateCalibrator();
  }
  AliCodeTimerAuto(Form("%s::Calibrate(AliMUONVDigitStore*)",fDigitCalibrator->ClassName()))
  fDigitCalibrator->Calibrate(digitStore);  
}

//_____________________________________________________________________________
void
AliMUONReconstructor::Clusterize(const AliMUONVDigitStore& digitStore,
                                 AliMUONVClusterStore& clusterStore) const
{
  /// Creates clusters from digits.

  TString sopt(GetOption());
  sopt.ToUpper();
  if ( sopt.Contains("NOCLUSTERING") ) return;
  
  if  (!fClusterReconstructor)
  {
    CreateClusterReconstructor();
  }
  
  AliCodeTimerAuto(Form("%s::Digits2Clusters(const AliMUONVDigitStore&,AliMUONVClusterStore&)",
                        fClusterReconstructor->ClassName()))
  fClusterReconstructor->Digits2Clusters(digitStore,clusterStore);  
}

//_____________________________________________________________________________
AliMUONVClusterStore*
AliMUONReconstructor::ClusterStore() const
{
  /// Return (and create if necessary) the cluster container
  if (!fClusterStore) 
  {
    fClusterStore = new AliMUONClusterStoreV1;
  }
  return fClusterStore;
}

//_____________________________________________________________________________
void
AliMUONReconstructor::ConvertDigits(AliRawReader* rawReader, 
                                    AliMUONVDigitStore* digitStore,
                                    AliMUONVTriggerStore* triggerStore) const
{
  /// Convert raw data into digit and trigger stores
  CreateDigitMaker();
  
  AliCodeTimerStart(Form("%s::Raw2Digits(AliRawReader*,AliMUONVDigitStore*,AliMUONVTriggerStore*)",
                    fDigitMaker->ClassName()))
  fDigitMaker->Raw2Digits(rawReader,digitStore,triggerStore);
  AliCodeTimerStop(Form("%s::Raw2Digits(AliRawReader*,AliMUONVDigitStore*,AliMUONVTriggerStore*)",
                         fDigitMaker->ClassName()))
  Calibrate(*digitStore);
}

//_____________________________________________________________________________
void 
AliMUONReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{
   /// convert raw data into a digit tree

  Bool_t alone = ( TriggerStore() == 0 );
  
  Bool_t ok = DigitStore()->Connect(*digitsTree,alone);
  if ( TriggerStore() ) 
  {
    ok = ok && TriggerStore()->Connect(*digitsTree,kFALSE);
  }
  
  if (!ok)
  {
    AliError("Could not make branches on TreeD");
  }
  else
  {
    ConvertDigits(rawReader,DigitStore(),TriggerStore());
    digitsTree->Fill();
    DigitStore()->Clear();
  }
}

//_____________________________________________________________________________
AliMUONTriggerCrateStore*
AliMUONReconstructor::CrateManager() const
{
  /// Return (and create if necessary) the trigger crate store
  if (fCrateManager) return fCrateManager;
  fCrateManager = new AliMUONTriggerCrateStore;
  fCrateManager->ReadFromFile();
  return fCrateManager;
}

//_____________________________________________________________________________
void
AliMUONReconstructor::CreateDigitMaker() const
{
  /// Create (and create if necessary) the digit maker
  if (fDigitMaker) return;

  AliCodeTimerAuto("")

  fDigitMaker = new AliMUONDigitMaker;
}

//_____________________________________________________________________________
void 
AliMUONReconstructor::CreateTriggerCircuit() const
{
  /// Return (and create if necessary) the trigger circuit object
  if (fTriggerCircuit) return;

  AliCodeTimerAuto("")

  fTriggerCircuit = new AliMUONTriggerCircuit(fTransformer);

}

//_____________________________________________________________________________
void
AliMUONReconstructor::CreateTriggerChamberEff() const
{
  /// Create (and create if necessary) the trigger chamber efficiency class
  if (fTrigChamberEff) return;

  AliCodeTimerAuto("")

  fTrigChamberEff = new AliMUONTriggerChamberEff(fTransformer,fDigitMaker,kTRUE);
  //fTrigChamberEff->SetDebugLevel(1);
}

//_____________________________________________________________________________
AliTracker* 
AliMUONReconstructor::CreateTracker() const
{
  /// Create the MUONTracker object
  /// The MUONTracker is passed the GetOption(), i.e. our own options
  
  CreateTriggerCircuit();
  CreateDigitMaker();
  CreateTriggerChamberEff();
  
  AliMUONTracker* tracker = new AliMUONTracker(fDigitMaker,fTransformer,fTriggerCircuit,fTrigChamberEff);
  tracker->SetOption(GetOption());
  
  return tracker;
}

//_____________________________________________________________________________
void
AliMUONReconstructor::CreateClusterReconstructor() const
{
  /// Create cluster reconstructor, depending on GetOption()
  
  AliCodeTimerAuto("")

  AliDebug(1,"");
  
  AliMUONVClusterFinder* clusterFinder(0x0);
  
  TString opt(GetOption());
  opt.ToUpper();
  
  if ( strstr(opt,"PRECLUSTERV2") )
  {
    clusterFinder = new AliMUONPreClusterFinderV2;
  }    
  else if ( strstr(opt,"PRECLUSTERV3") )
  {
    clusterFinder = new AliMUONPreClusterFinderV3;
  }  
  else if ( strstr(opt,"PRECLUSTER") )
  {
    clusterFinder = new AliMUONPreClusterFinder;
  }  
  else if ( strstr(opt,"COG") )
  {
    clusterFinder = new AliMUONClusterFinderCOG(new AliMUONPreClusterFinder);
  }  
  else if ( strstr(opt,"SIMPLEFITV3") )
  {
    clusterFinder = new AliMUONClusterFinderSimpleFit(new AliMUONClusterFinderCOG(new AliMUONPreClusterFinderV3));
  }
  else if ( strstr(opt,"SIMPLEFIT") )
  {
    clusterFinder = new AliMUONClusterFinderSimpleFit(new AliMUONClusterFinderCOG(new AliMUONPreClusterFinder));
  }
  else if ( strstr(opt,"MLEM:DRAW") )
  {
    clusterFinder = new AliMUONClusterFinderMLEM(kTRUE,new AliMUONPreClusterFinder);
  }
  else if ( strstr(opt,"MLEMV3") )
  {
    clusterFinder = new AliMUONClusterFinderMLEM(kFALSE,new AliMUONPreClusterFinderV3);
  } 
  else if ( strstr(opt,"MLEMV2") )
  {
    clusterFinder = new AliMUONClusterFinderMLEM(kFALSE,new AliMUONPreClusterFinderV2);
  } 
  else if ( strstr(opt,"MLEM") )
  {
    clusterFinder = new AliMUONClusterFinderMLEM(kFALSE,new AliMUONPreClusterFinder);
  } 
  else if ( strstr(opt,"AZ") )
  {
    clusterFinder = new AliMUONClusterFinderAZ;
  }
  else
  {
    // default is currently AZ
    clusterFinder = new AliMUONClusterFinderAZ;
  }
  
  if ( clusterFinder ) 
  {
    AliInfo(Form("Will use %s for clusterizing",clusterFinder->ClassName()));
  }
  
  fClusterReconstructor = new AliMUONClusterReconstructor(clusterFinder,fTransformer);
}

//_____________________________________________________________________________
void
AliMUONReconstructor::CreateCalibrator() const
{
  /// Create the calibrator
  
  AliCodeTimerAuto("")
  
  Int_t runNumber = AliCDBManager::Instance()->GetRun();

  AliInfo("Calibration will occur.");
  
  fCalibrationData = new AliMUONCalibrationData(runNumber);
  if ( !fCalibrationData->IsValid() )
  {
    AliError("Could not retrieve calibrations !");
    delete fCalibrationData;
    fCalibrationData = 0x0;
    return;
  }    
  
  // Check that we get all the calibrations we'll need
  if ( !fCalibrationData->Pedestals() ||
       !fCalibrationData->Gains() ||
       !fCalibrationData->HV() )
  {
    AliFatal("Could not access all required calibration data");
  }
  
  TString opt(GetOption());
  opt.ToUpper();
  Bool_t statusMap(kTRUE);
  
  if ( strstr(opt,"NOSTATUSMAP") )
  {
    AliWarning("Disconnecting status map : SHOULD BE USED FOR DEBUG ONLY. NOT FOR PRODUCTION !!!");
    statusMap = kFALSE; 
  }
  fDigitCalibrator = new AliMUONDigitCalibrator(*fCalibrationData,statusMap);
}

//_____________________________________________________________________________
AliMUONVDigitStore*
AliMUONReconstructor::DigitStore() const
{
  /// Return (and create if necessary) the digit container
  if (!fDigitStore) 
  {
    TString sopt(GetOption());
    sopt.ToUpper();
    
    AliInfo(Form("Options=%s",sopt.Data()));
    
    if ( sopt.Contains("DIGITSTOREV1") )
    {
      fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV1");
    }
    else if ( sopt.Contains("DIGITSTOREV2R") ) 
    {
      fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV2R");
    }
    else if ( sopt.Contains("DIGITSTOREV2S") ) 
    {
      fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV2S");
    }
    
    if (!fDigitStore) fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV2R");
    
    AliInfo(Form("Will use %s to store digits during reconstruction",fDigitStore->ClassName()));
  }
  return fDigitStore;
}

//_____________________________________________________________________________
void
AliMUONReconstructor::FillTreeR(AliMUONVTriggerStore* triggerStore,
                                AliMUONVClusterStore* clusterStore,
                                TTree& clustersTree) const
{
  /// Write the trigger and cluster information into TreeR
  
  AliCodeTimerAuto("")

  AliDebug(1,"");
  
  Bool_t ok(kFALSE);
  if ( triggerStore ) 
  {
    Bool_t alone = ( clusterStore ? kFALSE : kTRUE );
    ok = triggerStore->Connect(clustersTree,alone);
    if (!ok)
    {
      AliError("Could not create triggerStore branches in TreeR");
    }
  }
  
  if ( clusterStore ) 
  {
    Bool_t alone = ( triggerStore ? kFALSE : kTRUE );
    ok = clusterStore->Connect(clustersTree,alone);
    if (!ok)
    {
      AliError("Could not create triggerStore branches in TreeR");
    }    
  }
  
  if (ok) // at least one type of branches created successfully
  {
    clustersTree.Fill();
  }
}

//_____________________________________________________________________________
Bool_t 
AliMUONReconstructor::HasDigitConversion() const
{
  /// We *do* have digit conversion, but we might advertise it only 
  /// if we want to save the digits.
  
  TString opt(GetOption());
  opt.ToUpper();
  if ( opt.Contains("SAVEDIGITS" ) && !opt.Contains("NOLOCALRECONSTRUCTION") )
  {
    return kTRUE;
  }
  else
  {
    return kFALSE;
  }
}

//_____________________________________________________________________________
void 
AliMUONReconstructor::Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const
{
  /// This method is called by AliReconstruction if HasLocalReconstruction()==kTRUE AND
  /// HasDigitConversion()==kFALSE
  
  if ( !clustersTree ) 
  {
    AliError("clustersTree is 0x0 !");
    return;
  }
  
  ConvertDigits(rawReader,DigitStore(),TriggerStore());
  Clusterize(*(DigitStore()),*(ClusterStore()));
    
  FillTreeR(TriggerStore(),ClusterStore(),*clustersTree);
}

//_____________________________________________________________________________
void 
AliMUONReconstructor::Reconstruct(TTree* digitsTree, TTree* clustersTree) const
{
  /// This method is called by AliReconstruction if HasLocalReconstruction()==kTRUE
  /// AND HasDigitConversion()==kTRUE
  
//  AliCodeTimerAuto("(TTree*,TTree*)")
  
  AliDebug(1,"");
  
  if (!digitsTree || !clustersTree) 
  {
    AliError(Form("Tree is null : digitsTree=%p clustersTree=%p",
                  digitsTree,clustersTree));
    return;
  }

  if (!fDigitStore)
  {
    fDigitStore = AliMUONVDigitStore::Create(*digitsTree);
    if (!fDigitStore)
    {
      AliError(Form("Could not get DigitStore from %s",digitsTree->GetName()));
    }
    else
    {
      AliInfo(Form("Created %s from %s",fDigitStore->ClassName(),digitsTree->GetName()));
    }
  }
  if (!fTriggerStore)
  {
    fTriggerStore = AliMUONVTriggerStore::Create(*digitsTree);
    if (!fTriggerStore)
    {
      AliError(Form("Could not get TriggerStore from %s",digitsTree->GetName()));
    }
    else
    {
      AliInfo(Form("Created %s from %s",fTriggerStore->ClassName(),digitsTree->GetName()));
    }
  }
  
  if (!fTriggerStore && !fDigitStore)
  {
    AliError("No store at all. Nothing to do.");
    return;
  }
  
  // insure we start with empty stores
  if ( fDigitStore ) 
  {
    fDigitStore->Clear(); 
    Bool_t alone = ( fTriggerStore ? kFALSE : kTRUE );
    Bool_t ok = fDigitStore->Connect(*digitsTree,alone);
    if (!ok)
    {
      AliError("Could not connect digitStore to digitsTree");
      return;
    }
  }
  if ( fTriggerStore ) 
  {
    fTriggerStore->Clear();
    Bool_t alone = ( fDigitStore ? kFALSE : kTRUE );
    Bool_t ok = fTriggerStore->Connect(*digitsTree,alone);
    if (!ok)
    {
      AliError("Could not connect triggerStore to digitsTree");
      return;
    }
  }
  
  digitsTree->GetEvent(0);
  
  if ( fDigitStore ) 
  {
    // Insure we got calibrated digits (if we reconstruct from pure simulated,
    // i.e. w/o going through raw data, this will be the case)
    TIter next(fDigitStore->CreateIterator());
    AliMUONVDigit* digit = static_cast<AliMUONVDigit*>(next());
    if (!digit->IsCalibrated())
    {
      Calibrate(*fDigitStore);
    }
    Clusterize(*fDigitStore,*(ClusterStore()));
  }
    
  FillTreeR(fTriggerStore,ClusterStore(),*clustersTree);
}

//_____________________________________________________________________________
AliMUONVTriggerStore*
AliMUONReconstructor::TriggerStore() const
{
  /// Return (and create if necessary and allowed) the trigger container
  TString sopt(GetOption());
  sopt.ToUpper();
  
  if (sopt.Contains("TRIGGERDISABLE"))
  {
    delete fTriggerStore;
    fTriggerStore = 0x0;
  }
  else
  {
    if (!fTriggerStore)
    {
      fTriggerStore = new AliMUONTriggerStoreV1;
    }
  }
  return fTriggerStore;
}
