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

/* $Id: AliReconstruction.cxx 63911 2013-08-19 16:46:41Z hristov $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for running the reconstruction                                      //
//                                                                           //
// Clusters and tracks are created for all detectors and all events by       //
// typing:                                                                   //
//                                                                           //
//   AliReconstruction rec;                                                  //
//   rec.Run();                                                              //
//                                                                           //
// The Run method returns kTRUE in case of successful execution.             //
//                                                                           //
// If the input to the reconstruction are not simulated digits but raw data, //
// this can be specified by an argument of the Run method or by the method   //
//                                                                           //
//   rec.SetInput("...");                                                    //
//                                                                           //
// The input formats and the corresponding argument are:                     //
// - DDL raw data files: directory name, ends with "/"                       //
// - raw data root file: root file name, extension ".root"                   //
// - raw data DATE file: DATE file name, any other non-empty string          //
// - MC root files     : empty string, default                               //
//                                                                           //
// By default all events are reconstructed. The reconstruction can be        //
// limited to a range of events by giving the index of the first and the     //
// last event as an argument to the Run method or by calling                 //
//                                                                           //
//   rec.SetEventRange(..., ...);                                            //
//                                                                           //
// The index -1 (default) can be used for the last event to indicate no      //
// upper limit of the event range.                                           //
//                                                                           //
// In case of raw-data reconstruction the user can modify the default        //
// number of events per digits/clusters/tracks file. In case the option      //
// is not used the number is set 1. In case the user provides 0, than        //
// the number of events is equal to the number of events inside the          //
// raw-data file (i.e. one digits/clusters/tracks file):                     //
//                                                                           //
//   rec.SetNumberOfEventsPerFile(...);                                      //
//                                                                           //
//                                                                           //
// The name of the galice file can be changed from the default               //
// "galice.root" by passing it as argument to the AliReconstruction          //
// constructor or by                                                         //
//                                                                           //
//   rec.SetGAliceFile("...");                                               //
//                                                                           //
// The local reconstruction can be switched on or off for individual         //
// detectors by                                                              //
//                                                                           //
//   rec.SetRunLocalReconstruction("...");                                   //
//                                                                           //
// The argument is a (case sensitive) string with the names of the           //
// detectors separated by a space. The special string "ALL" selects all      //
// available detectors. This is the default.                                 //
//                                                                           //
// The reconstruction of the primary vertex position can be switched off by  //
//                                                                           //
//   rec.SetRunVertexFinder(kFALSE);                                         //
//                                                                           //
// The tracking and the creation of ESD tracks can be switched on for        //
// selected detectors by                                                     //
//                                                                           //
//   rec.SetRunTracking("...");                                              //
//                                                                           //
// Uniform/nonuniform field tracking switches (default: uniform field)       //
//                                                                           //
//   rec.SetUniformFieldTracking(); ( rec.SetUniformFieldTracking(kFALSE); ) //
//                                                                           //
// The filling of additional ESD information can be steered by               //
//                                                                           //
//   rec.SetFillESD("...");                                                  //
//                                                                           //
// Again, for both methods the string specifies the list of detectors.       //
// The default is "ALL".                                                     //
//                                                                           //
// The call of the shortcut method                                           //
//                                                                           //
//   rec.SetRunReconstruction("...");                                        //
//                                                                           //
// is equivalent to calling SetRunLocalReconstruction, SetRunTracking and    //
// SetFillESD with the same detector selecting string as argument.           //
//                                                                           //
// The reconstruction requires digits or raw data as input. For the creation //
// of digits and raw data have a look at the class AliSimulation.            //
//                                                                           //
// The input data of a detector can be replaced by the corresponding HLT     //
// data by calling (usual detector string)                                   //
// SetUseHLTData("...");                                                     //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include <TArrayD.h>
#include <TArrayF.h>
#include <TArrayS.h>
#include <TChain.h>
#include <TFile.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TMap.h>
#include <TObjArray.h>
#include <TPRegexp.h>
#include <TParameter.h>
#include <TPluginManager.h>
#include <TProof.h>
#include <TProofOutputFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <THashTable.h>
#include <TGrid.h>
#include <TMessage.h>
#include <TUrl.h>
#include <TRandom.h>
#include <THashList.h>

#include "AliAlignObj.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCTPRawStream.h"
#include "AliCascadeVertexer.h"
#include "AliCentralTrigger.h"
#include "AliCodeTimer.h"
#include "AliDAQ.h"
#include "AliDetectorRecoParam.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDPmdTrack.h"
#include "AliESDTagCreator.h"
#include "AliESDVertex.h"
#include "AliESDcascade.h"
#include "AliESDfriend.h"
#include "AliESDkink.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrack.h"
#include "AliEventInfo.h"
#include "AliGRPObject.h"
#include "AliGRPRecoParam.h"
#include "AliGenEventHeader.h"
#include "AliGeomManager.h"
#include "AliGlobalQADataMaker.h" 
#include "AliHeader.h"
#include "AliLog.h"
#include "AliMagF.h"
#include "AliMultiplicity.h"
#include "AliPID.h"
#include "AliPlaneEff.h"
#include "AliQAv1.h"
#include "AliQADataMakerRec.h" 
#include "AliQAManager.h"
#include "AliRawVEvent.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawHLTManager.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderFile.h"
#include "AliRawReaderRoot.h"
#include "AliRecoInputHandler.h"
#include "AliReconstruction.h"
#include "AliReconstructor.h"
#include "AliRun.h"
#include "AliRunInfo.h"
#include "AliRunLoader.h"
#include "AliSysInfo.h" // memory snapshots
#include "AliTrackPointArray.h"
#include "AliTracker.h"
#include "AliTriggerClass.h"
#include "AliTriggerCluster.h"
#include "AliTriggerIR.h"
#include "AliTriggerConfiguration.h"
#include "AliV0vertexer.h"
#include "AliVertexer.h"
#include "AliTrackleter.h"
#include "AliVertexerTracks.h"
#include "AliTriggerRunScalers.h"
#include "AliTriggerScalersRecord.h"
#include "AliCTPTimeParams.h" 
#include "AliESDHLTDecision.h"
#include "AliTriggerInput.h"
#include "AliLHCData.h"
#include "ARVersion.h"
#include <RVersion.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/resource.h>
ClassImp(AliReconstruction)

using std::endl;

//_____________________________________________________________________________
const char* AliReconstruction::fgkStopEvFName = "_stopEvent_";
const char* AliReconstruction::fgkDetectorName[AliReconstruction::kNDetectors] = {"ITS", "TPC", "TRD",
"TOF", "PHOS", 
"HMPID", "EMCAL", "MUON", "FMD", "ZDC", "PMD", "T0", "VZERO", "ACORDE","AD","FIT","MFT", "HLT"};

//_____________________________________________________________________________
AliReconstruction::AliReconstruction(const char* gAliceFilename) :
  TSelector(),
  fRunVertexFinder(kTRUE),
  fRunVertexFinderTracks(kTRUE),
  fRunMuonTracking(kFALSE),
  fRunMFTTrackingMU(kFALSE),
  fRunV0Finder(kTRUE),
  fRunCascadeFinder(kTRUE),
  fRunMultFinder(kTRUE),
  fStopOnError(kTRUE),
  fStopOnMissingTriggerFile(kTRUE),
  fWriteAlignmentData(kFALSE),
  fWriteESDfriend(kFALSE),
  fWriteHLTESD(kTRUE),
  fFillTriggerESD(kTRUE),
  //
  fWriteThisFriend(kFALSE),
  fSkipFriendsForLargeZ(kFALSE),
  fMaxFriendTracks(3000),
  fFractionFriends(0.03),
  fFractionHLTESD(0.05),
  fSkipFriendsCutZ(50),
  //
  fSkipIncompleteDAQ(kTRUE),
  fCleanESD(kTRUE),
  fV0DCAmax(3.),
  fV0CsPmin(0.),
  fDmax(50.),
  fZmax(50.),

  fCosmicAlias("kCosmic"),
  fLaserAlias("kCalibLaser"),
  fRunLocalReconstruction("ALL"),
  fRunTracking("ALL"),
  fFillESD("ALL"),
  fDeleteRecPoints(""),
  fDeleteDigits(""),
  fLoadCDB(""),
  fUseTrackingErrorsForAlignment(""),
  fGAliceFileName(gAliceFilename),
  fRawInput(""),
  fESDOutput(""),
  fProofOutputFileName(""),
  fProofOutputLocation(""),
  fProofOutputDataset(kFALSE),
  fProofOutputArchive(""),
  fEquipIdMap(""),
  fFirstEvent(0),
  fLastEvent(-1),
  fNumberOfEventsPerFile((UInt_t)-1),
  fOptions(),
  fLoadAlignFromCDB(kTRUE),
  fLoadAlignData("ALL"),
  fUseHLTData(),
  fRunInfo(NULL),
  fEventInfo(),
  fRunScalers(NULL),
  fCTPTimeParams(NULL),  
  fCTPTimeAlign(NULL),  

  fRunLoader(NULL),
  fRawReader(NULL),
  fParentRawReader(NULL),

  fRecoParam(),

  fSPDTrackleter(NULL),

  fDiamondProfileSPD(NULL),
  fDiamondProfile(NULL),
  fDiamondProfileTPC(NULL),
  fListOfCosmicTriggers(NULL), //RS for BWD comp.
  fAlias2Trigger(NULL),
  
  fGRPData(NULL),

  fAlignObjArray(NULL),
  fCDBUri(),
  fQARefUri(),
  fSpecCDBUri(), 
  fCheckRecoCDBvsSimuCDB(),
  fInitCDBCalled(kFALSE),
  fCDBSnapshotMode(kFALSE),
  fSetRunNumberFromDataCalled(kFALSE),
  fQADetectors("ALL"), 
  fQATasks("ALL"), 
  fRunQA(kTRUE),  
  fRunGlobalQA(kTRUE),
  fSameQACycle(kFALSE),
  fInitQACalled(kFALSE), 
  fWriteQAExpertData(kTRUE), 
  fRunPlaneEff(kFALSE),

  fESDpid(NULL),

  fesd(NULL),
  fhltesd(NULL),
  fesdf(NULL),
  fesdfDummy(NULL),
  ffile(NULL),
  ffileF(NULL),
  ftree(NULL),
  ftreeF(NULL),
  fhlttree(NULL),
  ftVertexer(NULL),
  fIsNewRunLoader(kFALSE),
  fRunAliEVE(kFALSE),
  fChain(NULL),
  fNall(0),
  fNspecie(0),
  fSspecie(0),
  fNhighPt(0),
  fShighPt(0),
  //
  fTreeBuffSize(30000000),
  fMemCountESD(0),
  fMemCountESDF(0),
  fMemCountESDHLT(0),
  //
  fUpgradeModule(""),
  fAnalysisMacro(),
  fAnalysis(0),
  fRecoHandler(0),
  fDeclTriggerClasses(""),
  fStopped(kFALSE),
  fMaxRSS(0),
  fMaxVMEM(0),
  fNAbandonedEv(0)
{
// create reconstruction object with default parameters
  AliGeomManager::Destroy();
  
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    fReconstructor[iDet] = NULL;
    fUpgradeMask[iDet]=kFALSE;
    fLoader[iDet] = NULL;
    fTracker[iDet] = NULL;
  }
  for (Int_t iDet = 0; iDet < AliQAv1::kNDET; iDet++) {
    fQACycles[iDet] = 999999 ;
    fQAWriteExpert[iDet] = kFALSE ; 
  }
  fBeamInt[0][0]=fBeamInt[0][1]=fBeamInt[1][0]=fBeamInt[1][1] = -1;
  //
  AddCheckRecoCDBvsSimuCDB("TPC/Calib/RecoParam"); // check for similarity in the sim and rec
  //
  AliPID pid;
}

//_____________________________________________________________________________
AliReconstruction::AliReconstruction(const AliReconstruction& rec) :
  TSelector(),
  fRunVertexFinder(rec.fRunVertexFinder),
  fRunVertexFinderTracks(rec.fRunVertexFinderTracks),
  fRunMuonTracking(rec.fRunMuonTracking),
  fRunMFTTrackingMU(rec.fRunMFTTrackingMU),
  fRunV0Finder(rec.fRunV0Finder),
  fRunCascadeFinder(rec.fRunCascadeFinder),
  fRunMultFinder(rec.fRunMultFinder),
  fStopOnError(rec.fStopOnError),
  fStopOnMissingTriggerFile(rec.fStopOnMissingTriggerFile),
  fWriteAlignmentData(rec.fWriteAlignmentData),
  fWriteESDfriend(rec.fWriteESDfriend),
  fWriteHLTESD(rec.fWriteHLTESD),
  fFillTriggerESD(rec.fFillTriggerESD),
  //
  fWriteThisFriend(rec.fWriteThisFriend),
  fSkipFriendsForLargeZ(rec.fSkipFriendsForLargeZ),
  fMaxFriendTracks(rec.fMaxFriendTracks),
  fFractionFriends(rec.fFractionFriends),
  fFractionHLTESD(rec.fFractionHLTESD),
  fSkipFriendsCutZ(rec.fSkipFriendsCutZ),
  //
  fSkipIncompleteDAQ(rec.fSkipIncompleteDAQ),
  fCleanESD(rec.fCleanESD),
  fV0DCAmax(rec.fV0DCAmax),
  fV0CsPmin(rec.fV0CsPmin),
  fDmax(rec.fDmax),
  fZmax(rec.fZmax),

  fCosmicAlias(rec.fCosmicAlias),
  fLaserAlias(rec.fLaserAlias),
  fRunLocalReconstruction(rec.fRunLocalReconstruction),
  fRunTracking(rec.fRunTracking),
  fFillESD(rec.fFillESD),
  fDeleteRecPoints(""),
  fDeleteDigits(""),
  fLoadCDB(rec.fLoadCDB),
  fUseTrackingErrorsForAlignment(rec.fUseTrackingErrorsForAlignment),
  fGAliceFileName(rec.fGAliceFileName),
  fRawInput(rec.fRawInput),
  fESDOutput(rec.fESDOutput),
  fProofOutputFileName(rec.fProofOutputFileName),
  fProofOutputLocation(rec.fProofOutputLocation),
  fProofOutputDataset(rec.fProofOutputDataset),
  fProofOutputArchive(rec.fProofOutputArchive),
  fEquipIdMap(rec.fEquipIdMap),
  fFirstEvent(rec.fFirstEvent),
  fLastEvent(rec.fLastEvent),
  fNumberOfEventsPerFile(rec.fNumberOfEventsPerFile),
  fOptions(),
  fLoadAlignFromCDB(rec.fLoadAlignFromCDB),
  fLoadAlignData(rec.fLoadAlignData),
  fUseHLTData(rec.fUseHLTData),
  fRunInfo(NULL),
  fEventInfo(),
  fRunScalers(NULL),
  fCTPTimeParams(NULL),
  fCTPTimeAlign(NULL),

  fRunLoader(NULL),
  fRawReader(NULL),
  fParentRawReader(NULL),

  fRecoParam(rec.fRecoParam),

  fSPDTrackleter(NULL),

  fDiamondProfileSPD(rec.fDiamondProfileSPD),
  fDiamondProfile(rec.fDiamondProfile),
  fDiamondProfileTPC(rec.fDiamondProfileTPC),
  fListOfCosmicTriggers(NULL), //RS for BWD comp.
  fAlias2Trigger(NULL),
  
  fGRPData(NULL),

  fAlignObjArray(rec.fAlignObjArray),
  fCDBUri(rec.fCDBUri),
  fQARefUri(rec.fQARefUri),
  fSpecCDBUri(), 
  fCheckRecoCDBvsSimuCDB(),
  fInitCDBCalled(rec.fInitCDBCalled),
  fCDBSnapshotMode(rec.fCDBSnapshotMode),
  fSetRunNumberFromDataCalled(rec.fSetRunNumberFromDataCalled),
  fQADetectors(rec.fQADetectors), 
  fQATasks(rec.fQATasks), 
  fRunQA(rec.fRunQA),  
  fRunGlobalQA(rec.fRunGlobalQA),
  fSameQACycle(rec.fSameQACycle),
  fInitQACalled(rec.fInitQACalled),
  fWriteQAExpertData(rec.fWriteQAExpertData), 
  fRunPlaneEff(rec.fRunPlaneEff),

  fESDpid(NULL),

  fesd(NULL),
  fhltesd(NULL),
  fesdf(NULL),
  fesdfDummy(NULL),
  ffile(NULL),
  ffileF(NULL),
  ftree(NULL),
  ftreeF(NULL),
  fhlttree(NULL),
  ftVertexer(NULL),
  fIsNewRunLoader(rec.fIsNewRunLoader),
  fRunAliEVE(kFALSE),
  fChain(NULL),
  fNall(0),
  fNspecie(0),
  fSspecie(0),
  fNhighPt(0),
  fShighPt(0),
  //
  fTreeBuffSize(rec.fTreeBuffSize),
  fMemCountESD(0),
  fMemCountESDF(0),
  fMemCountESDHLT(0),
  //
  fUpgradeModule(""),
  fAnalysisMacro(rec.fAnalysisMacro),
  fAnalysis(0),
  fRecoHandler(0),
  fDeclTriggerClasses(rec.fDeclTriggerClasses),
  fStopped(kFALSE),
  fMaxRSS(0),
  fMaxVMEM(0),
  fNAbandonedEv(0)
{
// copy constructor

  for (Int_t i = 0; i < rec.fOptions.GetEntriesFast(); i++) {
    if (rec.fOptions[i]) fOptions.Add(rec.fOptions[i]->Clone());
  }
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    fReconstructor[iDet] = NULL;
    fUpgradeMask[iDet] = kFALSE;
    fLoader[iDet] = NULL;
    fTracker[iDet] = NULL;
  }  
  
  for (Int_t iDet = 0; iDet < AliQAv1::kNDET; iDet++) {
    fQACycles[iDet] = rec.fQACycles[iDet];
    fQAWriteExpert[iDet] = rec.fQAWriteExpert[iDet] ; 
  }

  for (Int_t i = 0; i < rec.fSpecCDBUri.GetEntriesFast(); i++) {
    if (rec.fSpecCDBUri[i]) fSpecCDBUri.Add(rec.fSpecCDBUri[i]->Clone());
  }

  for (Int_t i = 0; i < rec.fCheckRecoCDBvsSimuCDB.GetEntriesFast(); i++) {
    if (rec.fCheckRecoCDBvsSimuCDB[i]) fCheckRecoCDBvsSimuCDB.AddLast(rec.fCheckRecoCDBvsSimuCDB[i]->Clone());
  }

  for (int i=2;i--;) for (int j=2;j--;) fBeamInt[i][j] = rec.fBeamInt[i][j];

}

//_____________________________________________________________________________
AliReconstruction& AliReconstruction::operator = (const AliReconstruction& rec)
{
// assignment operator
// Used in PROOF mode
// Be very careful while modifing it!
// Simple rules to follow:
// for persistent data members - use their assignment operators
// for non-persistent ones - do nothing or take the default values from constructor
// TSelector members should not be touched
  if(&rec == this) return *this;

  fRunVertexFinder       = rec.fRunVertexFinder;
  fRunVertexFinderTracks = rec.fRunVertexFinderTracks;
  fRunMuonTracking       = rec.fRunMuonTracking;
  fRunMFTTrackingMU      = rec.fRunMFTTrackingMU;
  fRunV0Finder           = rec.fRunV0Finder;
  fRunCascadeFinder      = rec.fRunCascadeFinder;
  fRunMultFinder         = rec.fRunMultFinder;
  fStopOnError           = rec.fStopOnError;
  fStopOnMissingTriggerFile = rec.fStopOnMissingTriggerFile;
  fWriteAlignmentData    = rec.fWriteAlignmentData;
  fWriteESDfriend        = rec.fWriteESDfriend;
  fWriteHLTESD           = rec.fWriteHLTESD;
  fFillTriggerESD        = rec.fFillTriggerESD;
  //
  //
  fWriteThisFriend = rec.fWriteThisFriend;
  fSkipFriendsForLargeZ = rec.fSkipFriendsForLargeZ;
  fMaxFriendTracks = rec.fMaxFriendTracks;
  fFractionFriends = rec.fFractionFriends;
  fFractionHLTESD = rec.fFractionHLTESD;
  fSkipFriendsCutZ = rec.fSkipFriendsCutZ;
  //
  fSkipIncompleteDAQ = rec.fSkipIncompleteDAQ;
  fCleanESD  = rec.fCleanESD;
  fV0DCAmax  = rec.fV0DCAmax;
  fV0CsPmin  = rec.fV0CsPmin;
  fDmax      = rec.fDmax;
  fZmax      = rec.fZmax;

  fCosmicAlias                   = rec.fCosmicAlias;
  fLaserAlias                    = rec.fLaserAlias;

  fRunLocalReconstruction        = rec.fRunLocalReconstruction;
  fRunTracking                   = rec.fRunTracking;
  fFillESD                       = rec.fFillESD;
  fDeleteRecPoints               = rec.fDeleteRecPoints;
  fDeleteDigits                  = rec.fDeleteDigits;
  fLoadCDB                       = rec.fLoadCDB;
  fUseTrackingErrorsForAlignment = rec.fUseTrackingErrorsForAlignment;
  fGAliceFileName                = rec.fGAliceFileName;
  fRawInput                      = rec.fRawInput;
  fESDOutput                     = rec.fESDOutput;
  fProofOutputFileName           = rec.fProofOutputFileName;
  fProofOutputLocation           = rec.fProofOutputLocation;
  fProofOutputDataset            = rec.fProofOutputDataset;
  fProofOutputArchive            = rec.fProofOutputArchive;
  fEquipIdMap                    = rec.fEquipIdMap;
  fFirstEvent                    = rec.fFirstEvent;
  fLastEvent                     = rec.fLastEvent;
  fNumberOfEventsPerFile         = rec.fNumberOfEventsPerFile;

  for (Int_t i = 0; i < rec.fOptions.GetEntriesFast(); i++) {
    if (rec.fOptions[i]) fOptions.Add(rec.fOptions[i]->Clone());
  }

  fLoadAlignFromCDB              = rec.fLoadAlignFromCDB;
  fLoadAlignData                 = rec.fLoadAlignData;
  fUseHLTData                    = rec.fUseHLTData;

  delete fRunInfo; fRunInfo = NULL;
  if (rec.fRunInfo) fRunInfo = new AliRunInfo(*rec.fRunInfo);

  fEventInfo                     = rec.fEventInfo;

  delete fRunScalers; fRunScalers = NULL;
  if (rec.fRunScalers) fRunScalers = new AliTriggerRunScalers(*rec.fRunScalers); 

  delete fCTPTimeParams; fCTPTimeParams = NULL;
  if (rec.fCTPTimeParams) fCTPTimeParams = new AliCTPTimeParams(*rec.fCTPTimeParams);
  delete fCTPTimeAlign; fCTPTimeAlign = NULL;
  if (rec.fCTPTimeAlign) fCTPTimeAlign = new AliCTPTimeParams(*rec.fCTPTimeAlign);

  fRunLoader       = NULL;
  fRawReader       = NULL;
  fParentRawReader = NULL;

  fRecoParam = rec.fRecoParam;

  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    fUpgradeMask[iDet] = kFALSE;
    delete fReconstructor[iDet]; fReconstructor[iDet] = NULL;
    delete fLoader[iDet]; fLoader[iDet] = NULL;
    delete fTracker[iDet]; fTracker[iDet] = NULL;
  }
  
  for (Int_t iDet = 0; iDet < AliQAv1::kNDET; iDet++) {
    fQACycles[iDet] = rec.fQACycles[iDet];
    fQAWriteExpert[iDet] = rec.fQAWriteExpert[iDet] ;
  } 

  delete fSPDTrackleter; fSPDTrackleter = NULL;
    
  delete fDiamondProfileSPD; fDiamondProfileSPD = NULL;
  if (rec.fDiamondProfileSPD) fDiamondProfileSPD = new AliESDVertex(*rec.fDiamondProfileSPD);
  delete fDiamondProfile; fDiamondProfile = NULL;
  if (rec.fDiamondProfile) fDiamondProfile = new AliESDVertex(*rec.fDiamondProfile);
  delete fDiamondProfileTPC; fDiamondProfileTPC = NULL;
  if (rec.fDiamondProfileTPC) fDiamondProfileTPC = new AliESDVertex(*rec.fDiamondProfileTPC);

  delete fListOfCosmicTriggers; fListOfCosmicTriggers = NULL;
  if (rec.fListOfCosmicTriggers) fListOfCosmicTriggers = (THashTable*)((rec.fListOfCosmicTriggers)->Clone());
  //
  delete fAlias2Trigger; fAlias2Trigger = NULL;
  if (rec.fAlias2Trigger) fAlias2Trigger = (THashList*)((rec.fAlias2Trigger)->Clone());

  delete fGRPData; fGRPData = NULL;
  //  if (rec.fGRPData) fGRPData = (TMap*)((rec.fGRPData)->Clone());
  if (rec.fGRPData) fGRPData = (AliGRPObject*)((rec.fGRPData)->Clone());

  delete fAlignObjArray; fAlignObjArray = NULL;

  fCDBUri        = "";
  fQARefUri      = rec.fQARefUri;
  fSpecCDBUri.Delete();
  fCheckRecoCDBvsSimuCDB.Delete();
  //
  for (Int_t i = 0; i < rec.fCheckRecoCDBvsSimuCDB.GetEntriesFast(); i++) {
    if (rec.fCheckRecoCDBvsSimuCDB[i]) fCheckRecoCDBvsSimuCDB.AddLast(rec.fCheckRecoCDBvsSimuCDB[i]->Clone());
  }
  //
  fInitCDBCalled               = rec.fInitCDBCalled;
  fCDBSnapshotMode             = rec.fCDBSnapshotMode;
  fSetRunNumberFromDataCalled  = rec.fSetRunNumberFromDataCalled;
  fQADetectors                 = rec.fQADetectors;
  fQATasks                     = rec.fQATasks; 
  fRunQA                       = rec.fRunQA;  
  fRunGlobalQA                 = rec.fRunGlobalQA;
  fSameQACycle                 = rec.fSameQACycle;
  fInitQACalled                = rec.fInitQACalled;
  fWriteQAExpertData           = rec.fWriteQAExpertData;
  fRunPlaneEff                 = rec.fRunPlaneEff;
  for (int i=2;i--;) for (int j=2;j--;) fBeamInt[i][j] = rec.fBeamInt[i][j];
  fESDpid  = NULL;
  fesd     = NULL;
  fhltesd  = NULL;
  fesdf    = NULL;
  fesdfDummy    = NULL;
  ffile    = NULL;
  ffileF   = NULL;
  ftree    = NULL;
  ftreeF   = NULL;
  fhlttree = NULL;
  ftVertexer = NULL;
  fIsNewRunLoader = rec.fIsNewRunLoader;
  fRunAliEVE = kFALSE;
  fChain = NULL;
  fNall = 0;
  fNspecie = 0;
  fSspecie = 0;
  fNhighPt = 0;
  fShighPt = 0;
  //
  fTreeBuffSize = rec.fTreeBuffSize;
  fMemCountESD = 0;
  fMemCountESDF = 0;
  fMemCountESDHLT = 0;
  //
  fUpgradeModule="";
  fAnalysisMacro = rec.fAnalysisMacro;
  fAnalysis = 0;
  fRecoHandler = 0;
  fDeclTriggerClasses = rec.fDeclTriggerClasses;

  return *this;
}

//_____________________________________________________________________________
AliReconstruction::~AliReconstruction()
{
// clean up

  CleanUp();
  if (fListOfCosmicTriggers) {
    fListOfCosmicTriggers->Delete();
    delete fListOfCosmicTriggers;
  }
  if (fAlias2Trigger) {
    fAlias2Trigger->Delete();
    delete fAlias2Trigger;
  }
  delete fGRPData;
  delete fRunScalers;
  delete fCTPTimeParams;
  delete fCTPTimeAlign;
  fOptions.Delete();
  if (fAlignObjArray) {
    fAlignObjArray->Delete();
    delete fAlignObjArray;
  }
  fSpecCDBUri.Delete();
  fCheckRecoCDBvsSimuCDB.Delete();
  AliCodeTimer::Instance()->Print();
}

//_____________________________________________________________________________
void AliReconstruction::InitQA()
{
  //Initialize the QA and start of cycle 
  AliCodeTimerAuto("",0);
  
  if (fInitQACalled) return;
  fInitQACalled = kTRUE;
  
  if (fGRPData) AliQADataMaker::SetCloningRequest( fGRPData->GetQATrigClasses(), fGRPData->GetQACloningRequest());


  AliQAManager * qam = AliQAManager::QAManager(AliQAv1::kRECMODE) ; 
  qam->SetSaveData(kTRUE); 
  qam->SetCycleLength(AliQAv1::kITS, 5) ; 
  if (fWriteQAExpertData)
    qam->SetWriteExpert() ; 
 
  if (qam->IsDefaultStorageSet()) {
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliWarning("Default QA reference storage has been already set !");
    AliWarning(Form("Ignoring the default storage declared in AliReconstruction: %s",fQARefUri.Data()));
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    fQARefUri = qam->GetDefaultStorage()->GetURI();
  } else {
    if (fQARefUri.Length() > 0) {
    	AliDebug(2,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    	AliDebug(2, Form("Default QA reference storage is set to: %s", fQARefUri.Data()));
    	AliDebug(2, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
      } else {
        fQARefUri="local://$ALICE_ROOT/QAref";
        AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        AliWarning("Default QA refeference storage not yet set !!!!");
        AliWarning(Form("Setting it now to: %s", fQARefUri.Data()));
        AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    		
      }
    qam->SetDefaultStorage(fQARefUri);
  }
  
  if (fRunQA) {
  qam->SetActiveDetectors(fQADetectors) ;
  qam->SetActiveOnlineDetectors(fRunInfo->GetActiveDetectors());
    
  for (Int_t det = 0 ; det < AliQAv1::kNDET ; det++) {
    qam->SetCycleLength(AliQAv1::DETECTORINDEX_t(det), fQACycles[det]) ;  
    qam->SetWriteExpert(AliQAv1::DETECTORINDEX_t(det)) ;
  }
  if (!fRawReader && !fInput && IsInTasks(AliQAv1::kRAWS))
    fQATasks.ReplaceAll(Form("%d",AliQAv1::kRAWS), "") ;
  qam->SetTasks(fQATasks) ; 
  qam->InitQADataMaker(AliCDBManager::Instance()->GetRun()) ; 
  }
  if (fRunGlobalQA) {
    Bool_t sameCycle = kFALSE ;
    AliQADataMaker *qadm = qam->GetQADataMaker(AliQAv1::kGLOBAL);
    AliInfo(Form("Initializing the global QA data maker"));
    if (IsInTasks(AliQAv1::kRECPOINTS)) {
      qadm->StartOfCycle(AliQAv1::kRECPOINTS, AliCDBManager::Instance()->GetRun(), sameCycle) ; 
      TObjArray **arr=qadm->Init(AliQAv1::kRECPOINTS);
      AliTracker::SetResidualsArray(arr);
      sameCycle = kTRUE ; 
    }
    if (IsInTasks(AliQAv1::kESDS)) {
      qadm->StartOfCycle(AliQAv1::kESDS, AliCDBManager::Instance()->GetRun(), sameCycle) ; 
      qadm->Init(AliQAv1::kESDS);
    }
  }
    AliSysInfo::AddStamp("InitQA") ; 
}

//_____________________________________________________________________________
void AliReconstruction::MergeQA(const char *fileName)
{
  //Initialize the QA and start of cycle 
  AliCodeTimerAuto("",0) ;
  AliQAManager::QAManager()->Merge(AliCDBManager::Instance()->GetRun(),fileName) ; 
  AliSysInfo::AddStamp("MergeQA") ; 
}
  
//_____________________________________________________________________________
void AliReconstruction::InitCDB()
{
// activate a default CDB storage
// First check if we have any CDB storage set, because it is used 
// to retrieve the calibration and alignment constants
  AliCodeTimerAuto("",0);

  if (fInitCDBCalled) return;
  fInitCDBCalled = kTRUE;

  AliCDBManager* man = AliCDBManager::Instance();
  if (man->IsDefaultStorageSet())
  {
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliWarning("Default CDB storage has been already set !");
    AliWarning(Form("Ignoring the default storage declared in AliReconstruction: %s",fCDBUri.Data()));
    AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    fCDBUri = man->GetDefaultStorage()->GetURI();
  }
  else {
    if (fCDBUri.Length() > 0) 
    {
    	AliDebug(2,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    	AliDebug(2, Form("Default CDB storage is set to: %s", fCDBUri.Data()));
    	AliDebug(2, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	man->SetDefaultStorage(fCDBUri);
    } 
    else if (!man->GetRaw()){
	fCDBUri="local://$ALICE_ROOT/OCDB";
    	AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    	AliWarning("Default CDB storage not yet set !!!!");
    	AliWarning(Form("Setting it now to: %s", fCDBUri.Data()));
    	AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	man->SetDefaultStorage(fCDBUri);
    }
    else {    
    	AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	AliWarning("Default storage will be set after setting the Run Number!!!");
    	AliWarning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");    		
    }
  }

  // Now activate the detector specific CDB storage locations
  for (Int_t i = 0; i < fSpecCDBUri.GetEntriesFast(); i++) {
    TObject* obj = fSpecCDBUri[i];
    if (!obj) continue;
    AliDebug(2, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    AliDebug(2, Form("Specific CDB storage for %s is set to: %s",obj->GetName(),obj->GetTitle()));
    AliDebug(2, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    man->SetSpecificStorage(obj->GetName(), obj->GetTitle());
  }
  AliSysInfo::AddStamp("InitCDB");
}

//_____________________________________________________________________________
void AliReconstruction::SetCDBSnapshotMode(const char* snapshotFileName) {
    if (!AliCDBManager::Instance()->SetSnapshotMode(snapshotFileName))
      AliFatal("Setting CDB snapshot mode failed.");
    fCDBSnapshotMode = kTRUE;
}

//_____________________________________________________________________________
void AliReconstruction::SetDefaultStorage(const char* uri) {
// Store the desired default CDB storage location
// Activate it later within the Run() method

  fCDBUri = uri;

}

//_____________________________________________________________________________
void AliReconstruction::SetQARefDefaultStorage(const char* uri) {
  // Store the desired default CDB storage location
  // Activate it later within the Run() method
  
  fQARefUri = uri;
  AliQAv1::SetQARefStorage(fQARefUri.Data()) ;
  
}

//_____________________________________________________________________________
void AliReconstruction::SetSpecificStorage(const char* calibType, const char* uri) {
// Store a detector-specific CDB storage location
// Activate it later within the Run() method

  AliCDBPath aPath(calibType);
  if(!aPath.IsValid()){
	// if calibType is not wildcard but it is a valid detector, add "/*" to make it a valid path
	for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
		if(!strcmp(calibType, fgkDetectorName[iDet])) {
			aPath.SetPath(Form("%s/*", calibType));
			AliInfo(Form("Path for specific storage set to %s", aPath.GetPath().Data()));
			break;
		}
        }
	if(!aPath.IsValid()){
  		AliError(Form("Not a valid path or detector: %s", calibType));
  		return;
	}
  }

//  // check that calibType refers to a "valid" detector name
//  Bool_t isDetector = kFALSE;
//  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
//    TString detName = fgkDetectorName[iDet];
//    if(aPath.GetLevel0() == detName) {
//    	isDetector = kTRUE;
//	break;
//    }
//  }
//
//  if(!isDetector) {
//	AliError(Form("Not a valid detector: %s", aPath.GetLevel0().Data()));
//	return;
//  }

  TObject* obj = fSpecCDBUri.FindObject(aPath.GetPath().Data());
  if (obj) fSpecCDBUri.Remove(obj);
  fSpecCDBUri.Add(new TNamed(aPath.GetPath().Data(), uri));

}

//_____________________________________________________________________________
void AliReconstruction::AddCheckRecoCDBvsSimuCDB(const char* cdbpath,const char* comment) 
{
  // require the cdb item to be the same in the rec as in the sim
  // Activate it later within the Run() method
  TString newent = cdbpath;
  if (newent.IsNull()) return;
  TIter nextit(&fCheckRecoCDBvsSimuCDB);
  TNamed* cdbent=0;
  while ((cdbent=(TNamed*)nextit())) {
    TString str = cdbent->GetName();
    if (str==newent) {
      AliInfo(Form("%s is already in the list to check",cdbpath));
      return;
    }
  }
  fCheckRecoCDBvsSimuCDB.AddLast(new TNamed(cdbpath,comment));
  //
}

//_____________________________________________________________________________
void AliReconstruction::RemCheckRecoCDBvsSimuCDB(const char* cdbpath) 
{
  // require the cdb item to be the same in the rec as in the sim
  // Activate it later within the Run() method
  TString newent = cdbpath;
  if (newent.IsNull()) return;
  TIter nextit(&fCheckRecoCDBvsSimuCDB);
  TNamed* cdbent=0;
  while ((cdbent=(TNamed*)nextit())) {
    TString str = cdbent->GetName();
    if (str==newent) {
      AliInfo(Form("Removing %s from the list to check",cdbpath));
      delete fCheckRecoCDBvsSimuCDB.Remove(cdbent);
      fCheckRecoCDBvsSimuCDB.Compress();
      return;
    }
  }
  AliInfo(Form("%s is not in the list to check",cdbpath));
  //
}

//_____________________________________________________________________________
Bool_t AliReconstruction::SetRunNumberFromData()
{
  // The method is called in Run() in order
  // to set a correct run number.
  // In case of raw data reconstruction the
  // run number is taken from the raw data header

  if (fSetRunNumberFromDataCalled) return kTRUE;
  fSetRunNumberFromDataCalled = kTRUE;
  
  AliCDBManager* man = AliCDBManager::Instance();
  if(fRawReader) {
    if(fRawReader->NextEvent()) {
      if(man->GetRun() > 0) {
  	AliWarning("Run number is taken from raw-event header! Ignoring settings in AliCDBManager!");
      } 
      man->SetRun(fRawReader->GetRunNumber());
      GetEventInfo();
      fRawReader->RewindEvents();
    }
    else {
      if(man->GetRun() > 0) {
	AliWarning("No raw-data events are found ! Using settings in AliCDBManager !");
      }
      else {
	AliWarning("Neither raw events nor settings in AliCDBManager are found !");
	return kFALSE;
      }
    }
  }
  else {
    AliRunLoader *rl = AliRunLoader::Open(fGAliceFileName.Data());
    if (!rl) {
      AliError(Form("No run loader found in file %s", fGAliceFileName.Data()));
      return kFALSE;
    }
    else {
      rl->LoadHeader();
      // read run number from gAlice
      if(rl->GetHeader()) {
	man->SetRun(rl->GetHeader()->GetRun());
	rl->UnloadHeader();
	delete rl;
      }
      else {
	AliError("Neither run-loader header nor RawReader objects are found !");
	delete rl;
	return kFALSE;
      }
    }
  }

  man->Print();  
  
  return kTRUE;
}

//_____________________________________________________________________________
void AliReconstruction::SetCDBLock() {
  // Set CDB lock: from now on it is forbidden to reset the run number
  // or the default storage or to activate any further storage!
  
  AliCDBManager::Instance()->SetLock(1);
}

//_____________________________________________________________________________
void AliReconstruction::MatchUpgradeDetector() {
  // Translates detector name in a boolean.
  // The boolean is used in GetReconstructor to load the 
  // upgrade reconstructor instead of the standard one.
   for(Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if(fUpgradeModule.Contains(fgkDetectorName[iDet])) fUpgradeMask[iDet]=kTRUE;
   }
}
//_____________________________________________________________________________
Bool_t AliReconstruction::MisalignGeometry(const TString& detectors)
{
  // Read the alignment objects from CDB.
  // Each detector is supposed to have the
  // alignment objects in DET/Align/Data CDB path.
  // All the detector objects are then collected,
  // sorted by geometry level (starting from ALIC) and
  // then applied to the TGeo geometry.
  // Finally an overlaps check is performed.

  // Load alignment data from CDB and fill fAlignObjArray 
  if(fLoadAlignFromCDB){
  	
    TString detStr = detectors;
    TString loadAlObjsListOfDets = "";
    
    for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
      if(!IsSelected(fgkDetectorName[iDet], detStr)) continue;
      if(!strcmp(fgkDetectorName[iDet],"HLT")) continue;
      
      if(AliGeomManager::GetNalignable(fgkDetectorName[iDet]) != 0)
      {
	loadAlObjsListOfDets += fgkDetectorName[iDet];
	loadAlObjsListOfDets += " ";
      }
    } // end loop over detectors
    
    if(AliGeomManager::GetNalignable("GRP") != 0)
      loadAlObjsListOfDets.Prepend("GRP "); //add alignment objects for non-sensitive modules
    AliGeomManager::ApplyAlignObjsFromCDB(loadAlObjsListOfDets.Data());
  }else{
    // Check if the array with alignment objects was
    // provided by the user. If yes, apply the objects
    // to the present TGeo geometry
    if (fAlignObjArray) {
      if (gGeoManager && gGeoManager->IsClosed()) {
	if (AliGeomManager::ApplyAlignObjsToGeom(*fAlignObjArray) == kFALSE) {
	  AliError("The misalignment of one or more volumes failed!"
		   "Compare the list of simulated detectors and the list of detector alignment data!");
	  return kFALSE;
	}
      }
      else {
	AliError("Can't apply the misalignment! gGeoManager doesn't exist or it is still opened!");
	return kFALSE;
      }
    }
  }
  
  if (fAlignObjArray) {
    fAlignObjArray->Delete();
    delete fAlignObjArray; fAlignObjArray=NULL;
  }

  return kTRUE;
}

//_____________________________________________________________________________
void AliReconstruction::SetGAliceFile(const char* fileName)
{
// set the name of the galice file

  fGAliceFileName = fileName;
}

//_____________________________________________________________________________
void AliReconstruction::SetInput(const char* input) 
{
  // In case the input string starts with 'mem://', we run in an online mode
  // and AliRawReaderDateOnline object is created. In all other cases a raw-data
  // file is assumed. One can give as an input:
  // mem://: - events taken from DAQ monitoring libs online
  //  or
  // mem://<filename> - emulation of the above mode (via DATE monitoring libs)
  if (input) fRawInput = input;
}

//_____________________________________________________________________________
void AliReconstruction::SetOutput(const char* output) 
{
  // Set the output ESD filename
  // 'output' is a normalt ROOT url
  // The method is used in case of raw-data reco with PROOF
  if (output) fESDOutput = output;
}

//_____________________________________________________________________________
void AliReconstruction::SetOption(const char* detector, const char* option)
{
// set options for the reconstruction of a detector

  TObject* obj = fOptions.FindObject(detector);
  if (obj) fOptions.Remove(obj);
  fOptions.Add(new TNamed(detector, option));
}

//_____________________________________________________________________________
void AliReconstruction::SetRecoParam(const char* detector, AliDetectorRecoParam *par)
{
  // Set custom reconstruction parameters for a given detector
  // Single set of parameters for all the events

  // First check if the reco-params are global
  if(!strcmp(detector, "GRP")) {
    par->SetAsDefault();
    fRecoParam.AddDetRecoParam(kNDetectors,par);
    return;
  }

  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if(!strcmp(detector, fgkDetectorName[iDet])) {
      par->SetAsDefault();
      fRecoParam.AddDetRecoParam(iDet,par);
      break;
    }
  }

}

//_____________________________________________________________________________
Bool_t AliReconstruction::InitGRP() {
  //------------------------------------
  // Initialization of the GRP entry 
  //------------------------------------
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");

  if (entry) {

    TMap* m = dynamic_cast<TMap*>(entry->GetObject());  // old GRP entry

    if (m) {
       AliInfo("Found a TMap in GRP/GRP/Data, converting it into an AliGRPObject");
       m->Print();
       fGRPData = new AliGRPObject();
       fGRPData->ReadValuesFromMap(m);
    }

    else {
       AliInfo("Found an AliGRPObject in GRP/GRP/Data, reading it");
       fGRPData = dynamic_cast<AliGRPObject*>(entry->GetObject());  // new GRP entry
       entry->SetOwner(0); // RS don't manipulate with GRP
       entry->SetObject(0); // others will reread the object
       AliCDBManager::Instance()->UnloadFromCache("GRP/GRP/Data");
    }

    //    FIX ME: The unloading of GRP entry is temporarily disabled
    //    because ZDC and VZERO are using it in order to initialize
    //    their reconstructor objects. In the future one has to think
    //    of propagating AliRunInfo to the reconstructors.
    //    AliCDBManager::Instance()->UnloadFromCache("GRP/GRP/Data");
  }

  if (!fGRPData) {
     AliError("No GRP entry found in OCDB!");
     return kFALSE;
  }

  TString lhcState = fGRPData->GetLHCState();
  if (lhcState==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the LHC state ! Using UNKNOWN");
    lhcState = "UNKNOWN";
  }

  TString beamType = fGRPData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
    beamType = "UNKNOWN";
  }

  Float_t beamEnergy = fGRPData->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
    beamEnergy = 0;
  }

  TString runType = fGRPData->GetRunType();
  if (runType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the run type ! Using UNKNOWN");
    runType = "UNKNOWN";
  }

  Int_t activeDetectors = fGRPData->GetDetectorMask();
  if (activeDetectors==AliGRPObject::GetInvalidUInt()) {
    AliError("GRP/GRP/Data entry:  missing value for the detector mask ! Using 1074790399");
    activeDetectors = 1074790399;
  }
  AliDebug(1, Form("activeDetectors = %d", activeDetectors));

  fRunInfo = new AliRunInfo(lhcState, beamType, beamEnergy, runType, activeDetectors);
  fRunInfo->Dump();


  // Process the list of active detectors
  if (activeDetectors) {
    UInt_t detMask = activeDetectors;
    AliDebug(1, Form("Detector List = %s", fRunLocalReconstruction.Data()));
    fRunLocalReconstruction = MatchDetectorList(fRunLocalReconstruction,detMask);
    AliDebug(1, Form("Detector List = %s", fRunLocalReconstruction.Data()));
    fRunTracking = MatchDetectorList(fRunTracking,detMask);
    fFillESD = MatchDetectorList(fFillESD,detMask);
    if (fRunQA || fRunGlobalQA) fQADetectors = MatchDetectorList(fQADetectors,detMask);
    AliInfo(Form("fQADetectors=%s",fQADetectors.Data()));
    fDeleteRecPoints = MatchDetectorList(fDeleteRecPoints,detMask);
    fDeleteDigits    = MatchDetectorList(fDeleteDigits,detMask);
    fLoadCDB.Form("%s %s %s %s",
		  fRunLocalReconstruction.Data(),
		  fRunTracking.Data(),
		  fFillESD.Data(),
		  fQADetectors.Data());
    fLoadCDB = MatchDetectorList(fLoadCDB,detMask);
    if (!((detMask >> AliDAQ::DetectorID("ITSSPD")) & 0x1) &&
	!((detMask >> AliDAQ::DetectorID("ITSSDD")) & 0x1) &&
	!((detMask >> AliDAQ::DetectorID("ITSSSD")) & 0x1) ) {
      // switch off the vertexer
      AliInfo("SPD,SDD,SSD is not in the list of active detectors. Vertexer and Trackleter are switched off.");
      fRunVertexFinder = kFALSE;
      fRunMultFinder = kFALSE;
    }
    if (!((detMask >> AliDAQ::DetectorID("TRG")) & 0x1)) {
      // switch off the reading of CTP raw-data payload
      if (fFillTriggerESD) {
	AliInfo("CTP is not in the list of active detectors. CTP data reading switched off.");
	fFillTriggerESD = kFALSE;
      }
    }
  }

  AliInfo("===================================================================================");
  AliInfo(Form("Running local reconstruction for detectors: %s",fRunLocalReconstruction.Data()));
  AliInfo(Form("Running tracking for detectors: %s",fRunTracking.Data()));
  AliInfo(Form("Filling ESD for detectors: %s",fFillESD.Data()));
  AliInfo(Form("Quality assurance is active for detectors: %s",fQADetectors.Data()));
  AliInfo(Form("CDB and reconstruction parameters are loaded for detectors: %s",fLoadCDB.Data()));
  AliInfo("===================================================================================");

  //*** Dealing with the magnetic field map
  if ( TGeoGlobalMagField::Instance()->IsLocked() ) {
    if (TGeoGlobalMagField::Instance()->GetField()->TestBit(AliMagF::kOverrideGRP)) {
      AliInfo("ExpertMode!!! GRP information will be ignored !");
      AliInfo("ExpertMode!!! Running with the externally locked B field !");
    }
    else {
      AliInfo("Destroying existing B field instance!");
      delete TGeoGlobalMagField::Instance();
    }    
  }
  if ( !TGeoGlobalMagField::Instance()->IsLocked() ) {
    // Construct the field map out of the information retrieved from GRP.
    Bool_t ok = kTRUE;
    // L3
    Float_t l3Current = fGRPData->GetL3Current((AliGRPObject::Stats)0);
    if (l3Current == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 current !");
      ok = kFALSE;
    }
    
    Char_t l3Polarity = fGRPData->GetL3Polarity();
    if (l3Polarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 polarity !");
      ok = kFALSE;
    }

    // Dipole
    Float_t diCurrent = fGRPData->GetDipoleCurrent((AliGRPObject::Stats)0);
    if (diCurrent == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole current !");
      ok = kFALSE;
    }

    Char_t diPolarity = fGRPData->GetDipolePolarity();
    if (diPolarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole polarity !");
      ok = kFALSE;
    }

    // read special bits for the polarity convention and map type
    Int_t  polConvention = fGRPData->IsPolarityConventionLHC() ? AliMagF::kConvLHC : AliMagF::kConvDCS2008;
    Bool_t uniformB = fGRPData->IsUniformBMap();

    if (ok) { 
      AliMagF* fld = AliMagF::CreateFieldMap(TMath::Abs(l3Current) * (l3Polarity ? -1:1), 
					     TMath::Abs(diCurrent) * (diPolarity ? -1:1), 
					     polConvention,uniformB,beamEnergy, beamType.Data());
      if (fld) {
	TGeoGlobalMagField::Instance()->SetField( fld );
	TGeoGlobalMagField::Instance()->Lock();
	AliInfo("Running with the B field constructed out of GRP !");
      }
      else AliFatal("Failed to create a B field map !");
    }
    else AliFatal("B field is neither set nor constructed from GRP ! Exitig...");
  }
  
  //*** Get the diamond profiles from OCDB
  entry = AliCDBManager::Instance()->Get("GRP/Calib/MeanVertexSPD");
  if (entry) {
    fDiamondProfileSPD = dynamic_cast<AliESDVertex*> (entry->GetObject());  
  } else {
     AliError("No SPD diamond profile found in OCDB!");
  }

  entry = AliCDBManager::Instance()->Get("GRP/Calib/MeanVertex");
  if (entry) {
    fDiamondProfile = dynamic_cast<AliESDVertex*> (entry->GetObject());  
  } else {
     AliError("No diamond profile found in OCDB!");
  }

  entry = AliCDBManager::Instance()->Get("GRP/Calib/MeanVertexTPC");
  if (entry) {
    fDiamondProfileTPC = dynamic_cast<AliESDVertex*> (entry->GetObject());  
  } else {
     AliError("No TPC diamond profile found in OCDB!");
  }

  return kTRUE;
} 

//_____________________________________________________________________________
Bool_t AliReconstruction::LoadCDB()
{
  // Load CDB entries for all active detectors.
  // By default we load all the entries in <det>/Calib
  // folder.

  AliCodeTimerAuto("",0);

  AliCDBManager::Instance()->Get("GRP/CTP/Config");

  AliCDBManager::Instance()->Get("GRP/Calib/LHCClockPhase");

  TString detStr = fLoadCDB;
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
    AliCDBManager::Instance()->GetAll(Form("%s/Calib/*",fgkDetectorName[iDet]));
    AliCDBManager::Instance()->GetAll(Form("%s/Trigger/*",fgkDetectorName[iDet]));
  }

  // Temporary fix - one has to define the correct policy in order
  // to load the trigger OCDB entries only for the detectors that
  // in the trigger or that are needed in order to put correct
  // information in ESD
  AliCDBManager::Instance()->GetAll("TRIGGER/*/*");
  AliCDBManager::Instance()->GetAll("HLT/*/*");

  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/Calib/CosmicTriggers");
  if (entry) {
    fListOfCosmicTriggers = dynamic_cast<THashTable*>(entry->GetObject());
    entry->SetOwner(0);
  }
  //
  if (!fListOfCosmicTriggers) {
    AliWarning("Can not get list of cosmic triggers from OCDB! Cosmic event specie will rely on aliases if defined");
  }


  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliReconstruction::LoadTriggerScalersCDB()
{
  // Load CTP scalers from OCDB.
  // The scalers are checked for consistency.

  AliCodeTimerAuto("",0);

  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Scalers");

  if (entry) { 
   
       AliInfo("Found an AliTriggerRunScalers in GRP/CTP/Scalers, reading it");
       fRunScalers = dynamic_cast<AliTriggerRunScalers*> (entry->GetObject());
       entry->SetOwner(0);
       if (fRunScalers && (fRunScalers->CorrectScalersOverflow() == 0)) AliInfo("32bit Trigger counters corrected for overflow");

  }
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliReconstruction::LoadCTPTimeParamsCDB()
{
  // Load CTP timing information (alignment)
  // from OCDB.

  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/CTPtiming");
  if (!entry) return kFALSE;

  AliInfo("Found an AliCTPTimeParams in GRP/CTP/CTPtiming, reading it");
  fCTPTimeParams = dynamic_cast<AliCTPTimeParams*> (entry->GetObject());
  entry->SetOwner(0);

  AliCDBEntry* entry2 = AliCDBManager::Instance()->Get("GRP/CTP/TimeAlign");
  if (!entry2) return kFALSE;

  AliInfo("Found an AliCTPTimeParams in GRP/CTP/TimeAlign, reading it");
  fCTPTimeAlign = dynamic_cast<AliCTPTimeParams*> (entry2->GetObject());
  entry2->SetOwner(0);

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::ReadIntensityInfoCDB()
{
  // Load LHC DIP data
  AliCDBEntry* entry    = AliCDBManager::Instance()->Get("GRP/GRP/LHCData");
  AliCDBEntry* entryCTP = AliCDBManager::Instance()->Get("GRP/CTP/Config");
  //
  if (!entry || !entryCTP) {
    AliError(Form("Failed to extract CDB objects GRP/GRP/LHCData: %p or GRP/CTP/Config: %p",entry,entryCTP));
    return kFALSE;
  }
  // extract BC masks
  enum {kA,kB,kC,kE,kNMasks};
  AliTriggerConfiguration* conf = (AliTriggerConfiguration*)entryCTP->GetObject();
  const TObjArray& clArr = conf->GetClasses();
  TObjArray masks(kNMasks);
  TIter next(&clArr);
  AliTriggerClass* trClass = 0;
  int nFound = 0;
  masks.SetOwner(kFALSE);
  //
  while ( (trClass=(AliTriggerClass*)next()) ) {
    TString trName = trClass->GetName();
    int ind = trName.Index("-"); // prefix in front of A,B,C,E
    if (ind<1) continue;   // anomaly
    //
    trName = trName.Data() + ind;
    AliTriggerBCMask* bcMask = trClass->GetBCMask();
    if (!bcMask) continue;
    UInt_t which = 0;
    if      (trName.BeginsWith("-A-"))  which |= 0x1<<kA;
    else if (trName.BeginsWith("-B-"))  which |= 0x1<<kB;
    else if (trName.BeginsWith("-C-"))  which |= 0x1<<kC;
    else if (trName.BeginsWith("-E-"))  which |= 0x1<<kE;
    else if (trName.BeginsWith("-AC-")) which |= (0x1<<kA) | (0x1<<kC);
    else if (trName.BeginsWith("-CA-")) which |= (0x1<<kA) | (0x1<<kC);
    else { AliWarning(Form("Unknown trigger type %s\n",trClass->GetName())); continue;}
    //
    for (int ip=kNMasks;ip--;) {
      if ( !(which&(0x1<<ip)) || masks[ip] ) continue; // does not match or already done
      masks[ip] = (TObject*)bcMask;
      nFound++;
    }
    if (nFound==kNMasks) break;
  }  
  //  
  AliInfo("Reading mean bunch intensities from GRP/GRP/LHCData");
  AliLHCData* dipData = dynamic_cast<AliLHCData*> (entry->GetObject());
  //
  for (int ib=2;ib--;) {
    double intI,intNI;
    if (dipData && (dipData->GetMeanIntensity(ib,intI,intNI,&masks)>=0)) {
      fBeamInt[ib][0] = intI;
      fBeamInt[ib][1] = intNI;	
      AliInfo(Form("Mean intensity for beam %d: Interacting:%.2e Non-Interacting:%.2e",ib,intI,intNI));
    }
  }
  return kTRUE;
  //
}


//_____________________________________________________________________________
Bool_t AliReconstruction::Run(const char* input)
{
  // Run Run Run
  AliCodeTimerAuto("",0);

  InitRun(input);
  if (GetAbort() != TSelector::kContinue) return kFALSE;

  TChain *chain = NULL;
  if (fRawReader && (chain = fRawReader->GetChain())) {
    ProcessTriggerAliases();
    Long64_t nEntries = (fLastEvent < 0) ? (TChain::kBigNumber) : (fLastEvent - fFirstEvent + 1);
    // Proof mode
    if (gProof) {
      // Temporary fix for long raw-data runs (until socket timeout handling in PROOF is revised)
      gProof->Exec("gEnv->SetValue(\"Proof.SocketActivityTimeout\",-1)", kTRUE);

      if (gGrid)
	gProof->Exec("TGrid::Connect(\"alien://\")",kTRUE);

      TMessage::EnableSchemaEvolutionForAll(kTRUE);
      gProof->Exec("TMessage::EnableSchemaEvolutionForAll(kTRUE)",kTRUE);

      gProof->AddInput(this);

      if (!ParseOutput()) return kFALSE;

      gProof->SetParameter("PROOF_MaxSlavesPerNode", 9999);
      chain->SetProof();
      chain->Process("AliReconstruction","",nEntries,fFirstEvent);
    }
    else {
      chain->Process(this,"",nEntries,fFirstEvent);
    }
  }
  else {
    Begin(NULL);
    if (GetAbort() != TSelector::kContinue) return kFALSE;
    SlaveBegin(NULL);
    if (GetAbort() != TSelector::kContinue) return kFALSE;
    //******* The loop over events
    AliInfo("Starting looping over events");
    Int_t iEvent = 0;
    while ((iEvent < fRunLoader->GetNumberOfEvents()) ||
	   (fRawReader && fRawReader->NextEvent())) {
      //
      // check if process has enough resources 
      if (!HasEnoughResources(iEvent)) break;
      if (!ProcessEvent(iEvent)) {
        Abort("ProcessEvent",TSelector::kAbortFile);
        return kFALSE;
      }
      CleanProcessedEvent();
      iEvent++;
    }
    if (!iEvent) AliWarning("No events passed trigger selection");
    SlaveTerminate();
    if (GetAbort() != TSelector::kContinue) return kFALSE;
    Terminate();
    if (GetAbort() != TSelector::kContinue) return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
void AliReconstruction::InitRawReader(const char* input)
{
  // Init raw-reader and
  // set the input in case of raw data

  AliCodeTimerAuto("",0);

  if (input) fRawInput = input;
  fRawReader = AliRawReader::Create(fRawInput.Data());
  if (!fRawReader) {
    if (fRawInput.IsNull()) {
      AliInfo("Reconstruction will run over digits");
    }
    else {
      AliFatal("Can not create raw-data reader ! Exiting..."); 
    }
  }

  if (!fEquipIdMap.IsNull() && fRawReader)
    fRawReader->LoadEquipmentIdsMap(fEquipIdMap);

  if (!fUseHLTData.IsNull()) {
    // create the RawReaderHLT which performs redirection of HLT input data for
    // the specified detectors
    AliRawReader* pRawReader=AliRawHLTManager::CreateRawReaderHLT(fRawReader, fUseHLTData.Data());
    if (pRawReader) {
      fParentRawReader=fRawReader;
      fRawReader=pRawReader;
    } else {
      AliError(Form("can not create Raw Reader for HLT input %s", fUseHLTData.Data()));
    }
  }
  AliSysInfo::AddStamp("CreateRawReader");
}

//_____________________________________________________________________________
void AliReconstruction::InitRun(const char* input)
{
  // Initialization of raw-reader,
  // run number, CDB etc.
  AliCodeTimerAuto("",0);
  AliSysInfo::AddStamp("Start");

  // Initialize raw-reader if any
  InitRawReader(input);

  // Initialize the CDB storage
  InitCDB();

  // Set run number in CDBManager (if it is not already set by the user)
  if (!SetRunNumberFromData()) {
    Abort("SetRunNumberFromData", TSelector::kAbortProcess);
    return;
  }

  // Set CDB lock: from now on it is forbidden to reset the run number
  // or the default storage or to activate any further storage!
  SetCDBLock();
  //
}

//_____________________________________________________________________________
void AliReconstruction::Begin(TTree *)
{
  // Initialize AlReconstruction before
  // going into the event loop
  // Should follow the TSelector convention
  // i.e. initialize only the object on the client side
  AliCodeTimerAuto("",0);

  AliReconstruction *reco = NULL;
  if (fInput) {
    if ((reco = (AliReconstruction*)fInput->FindObject("AliReconstruction"))) {
      *this = *reco;
    }
    AliSysInfo::AddStamp("ReadInputInBegin");
  }

  // Import ideal TGeo geometry and apply misalignment
  if (!AliGeomManager::GetGeometry()) {
    TString geom(gSystem->DirName(fGAliceFileName));
    geom += "/geometry.root";
    AliGeomManager::LoadGeometry(geom.Data());
    if (!gGeoManager) {
      Abort("LoadGeometry", TSelector::kAbortProcess);
      return;
    }
    AliSysInfo::AddStamp("LoadGeom");
    TString detsToCheck=fRunLocalReconstruction;
    if(!AliGeomManager::CheckSymNamesLUT(detsToCheck.Data())) {
      Abort("CheckSymNamesLUT", TSelector::kAbortProcess);
      return;
    }
    AliSysInfo::AddStamp("CheckGeom");
  }

  Bool_t toCDBSnapshot=kFALSE;
  TString snapshotFileOut(""); // we could use fSnapshotFileName if we are not interested
  // in reading from and writing to a snapshot file at the same time
  if(TString(getenv("OCDB_SNAPSHOT_CREATE")) == TString("kTRUE")){
      toCDBSnapshot=kTRUE;
      //fFromCDBSnapshot=kFALSE;
      TString snapshotFile(getenv("OCDB_SNAPSHOT_FILENAME"));
      if(!(snapshotFile.IsNull() || snapshotFile.IsWhitespace()))
	  snapshotFileOut = snapshotFile;
      else
	  snapshotFileOut="OCDB.root";
  }

  TString detStr(fLoadAlignData);
  if (!toCDBSnapshot) {
    if (!MisalignGeometry(fLoadAlignData)) {
      Abort("MisalignGeometry", TSelector::kAbortProcess);
      return;
    }
  } else {
    // when creating the snapshot, load the CDB alignment objects without applying them to the geometry
    for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
      if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
      if (!strcmp(fgkDetectorName[iDet],"HLT")) continue;
      if (AliGeomManager::GetNalignable(fgkDetectorName[iDet]) != 0)
      {
        TString detAlignPath = fgkDetectorName[iDet];
        detAlignPath += "/Align/Data";
        AliCDBManager::Instance()->Get(detAlignPath);
      }
    } // end loop over detectors
    if(AliGeomManager::GetNalignable("GRP") != 0)
      AliCDBManager::Instance()->Get("GRP/Align/Data");
  }

  const TMap* cdbCache = AliCDBManager::Instance()->GetEntryCache();
  if(!toCDBSnapshot) {
    if(cdbCache->Contains("GRP/Geometry/Data"))
      AliCDBManager::Instance()->UnloadFromCache("GRP/Geometry/Data");
    AliCDBManager::Instance()->UnloadFromCache("*/Align/*");
  }
  AliSysInfo::AddStamp("MisalignGeom");

  if (!InitGRP()) {
    Abort("InitGRP", TSelector::kAbortProcess);
    return;
  }
  AliSysInfo::AddStamp("InitGRP");
  if(!toCDBSnapshot)
      if(cdbCache->Contains("GRP/Calib/CosmicTriggers"))
	  AliCDBManager::Instance()->UnloadFromCache("GRP/Calib/CosmicTriggers");

  if(!fCDBSnapshotMode || toCDBSnapshot){
      if (!LoadCDB()) {
	  Abort("LoadCDB", TSelector::kAbortProcess);
	  return;
      }
      AliSysInfo::AddStamp("LoadCDB"); 
  }

  if (!LoadTriggerScalersCDB()) {
    Abort("LoadTriggerScalersCDB", TSelector::kAbortProcess);
    return;
  }
  AliSysInfo::AddStamp("LoadTriggerScalersCDB");

  if (!LoadCTPTimeParamsCDB()) {
    Abort("LoadCTPTimeParamsCDB", TSelector::kAbortProcess);
    return;
  }
  AliSysInfo::AddStamp("LoadCTPTimeParamsCDB");

  if (!ReadIntensityInfoCDB()) {
    Abort("ReadIntensityInfoCDB", TSelector::kAbortProcess);
    return;
  }
  AliSysInfo::AddStamp("ReadIntensityInfoCDB");

  // Read the reconstruction parameters from OCDB
  if (!InitRecoParams()) {
    AliWarning("Not all detectors have correct RecoParam objects initialized");
  }
  AliSysInfo::AddStamp("InitRecoParams");

  if(toCDBSnapshot)
  {
      AliCDBManager::Instance()->DumpToSnapshotFile(snapshotFileOut.Data(),kFALSE);
      exit(0);
  }

  if (fInput && gProof) {
    if (reco) *reco = *this;

    gGeoManager->SetName("Geometry");
    gProof->AddInputData(gGeoManager,kTRUE);
    gGeoManager = NULL;
    gProof->AddInputData(const_cast<TMap*>(AliCDBManager::Instance()->GetEntryCache()),kTRUE);
    fInput->Add(new TParameter<Int_t>("RunNumber",AliCDBManager::Instance()->GetRun()));
    AliMagF *magFieldMap = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
    magFieldMap->SetName("MagneticFieldMap");
    gProof->AddInputData(magFieldMap,kTRUE);
    if (fAnalysis) {
      fAnalysis->SetName("Analysis");
      gProof->AddInputData(fAnalysis,kTRUE);
    }  
  }

}

//_____________________________________________________________________________
void AliReconstruction::SlaveBegin(TTree*)
{
  // Initialization related to run-loader,
  // vertexer, trackers, recontructors
  // In proof mode it is executed on the slave
  AliCodeTimerAuto("",0);
  TProofOutputFile *outProofFile = NULL;
  if (fInput) {
    if (AliDebugLevel() > 0) fInput->Print();
    if (AliDebugLevel() > 10) fInput->Dump();
    if (AliReconstruction *reco = (AliReconstruction*)fInput->FindObject("AliReconstruction")) {
      *this = *reco;
    }
    if (TGeoManager *tgeo = (TGeoManager*)fInput->FindObject("Geometry")) {
      gGeoManager = tgeo;
      AliGeomManager::SetGeometry(tgeo);
    }
    if (TMap *entryCache = (TMap*)fInput->FindObject("CDBEntryCache")) {
      Int_t runNumber = -1;
      if (TProof::GetParameter(fInput,"RunNumber",runNumber) == 0) {
	AliCDBManager *man = AliCDBManager::Instance(entryCache,runNumber);
	man->SetCacheFlag(kTRUE);
	man->SetLock(kTRUE);
	man->Print();
      }
    }
    if (AliMagF *map = (AliMagF*)fInput->FindObject("MagneticFieldMap")) {
      AliMagF *newMap = new AliMagF(*map);
      if (!newMap->LoadParameterization()) {
	Abort("AliMagF::LoadParameterization", TSelector::kAbortProcess);
	return;
      }
      TGeoGlobalMagField::Instance()->SetField(newMap);
      TGeoGlobalMagField::Instance()->Lock();
    }
    if (!fAnalysis) {
       // Attempt to get the analysis manager from the input list
       fAnalysis = (AliAnalysisManager*)fInput->FindObject("Analysis");
       if (fAnalysis) AliInfo("==== Analysis manager retrieved from input list ====");
    }   
    if (TNamed *outputFileName = (TNamed*)fInput->FindObject("PROOF_OUTPUTFILE"))
      fProofOutputFileName = outputFileName->GetTitle();
    if (TNamed *outputLocation = (TNamed*)fInput->FindObject("PROOF_OUTPUTFILE_LOCATION"))
      fProofOutputLocation = outputLocation->GetTitle();
    if (fInput->FindObject("PROOF_OUTPUTFILE_DATASET"))
      fProofOutputDataset = kTRUE;
    if (TNamed *archiveList = (TNamed*)fInput->FindObject("PROOF_OUTPUTFILE_ARCHIVE"))
      fProofOutputArchive = archiveList->GetTitle();
    if (!fProofOutputFileName.IsNull() &&
	!fProofOutputLocation.IsNull() &&
	fProofOutputArchive.IsNull()) {
      if (!fProofOutputDataset) {
	outProofFile = new TProofOutputFile(fProofOutputFileName.Data(),"M");
	outProofFile->SetOutputFileName(Form("%s%s",fProofOutputLocation.Data(),fProofOutputFileName.Data()));
      }
      else {
	outProofFile = new TProofOutputFile(fProofOutputFileName.Data(),"DROV",fProofOutputLocation.Data());
      }
      if (AliDebugLevel() > 0) outProofFile->Dump();
      fOutput->Add(outProofFile);
    }
    AliSysInfo::AddStamp("ReadInputInSlaveBegin");
  }
  // Check if analysis was requested in the reconstruction event loop
  if (!fAnalysis) {
    // Attempt to connect in-memory singleton
    fAnalysis = AliAnalysisManager::GetAnalysisManager();
    if (fAnalysis) AliInfo(Form("==== Analysis manager <%s> found in memory ====", fAnalysis->GetName()));
    // Check if an analysis macro was specified
    if (!fAnalysis && !fAnalysisMacro.IsNull()) {
      // Run specified analysis macro
      gROOT->ProcessLine(Form(".x %s",fAnalysisMacro.Data()));
      fAnalysis = AliAnalysisManager::GetAnalysisManager();
      if (!fAnalysis) AliError(Form("No analysis manager produced by analysis macro %s", fAnalysisMacro.Data()));
      else AliInfo(Form("==== Analysis manager <%s> produced by analysis macro <%s> ====", 
                        fAnalysis->GetName(), fAnalysisMacro.Data()));
    }
  }
  
  // get the run loader
  if (!InitRunLoader()) {
    Abort("InitRunLoader", TSelector::kAbortProcess);
    return;
  }
  AliSysInfo::AddStamp("LoadLoader");
 
  CheckRecoCDBvsSimuCDB();

  ftVertexer = new AliVertexerTracks(AliTracker::GetBz());

  // get trackers
  if (!fRunTracking.IsNull() && !CreateTrackers(fRunTracking)) {
    Abort("CreateTrackers", TSelector::kAbortProcess);
    return;
  }      
  AliSysInfo::AddStamp("CreateTrackers");

  // create the ESD output file and tree
  if (!outProofFile) {
    ffile = TFile::Open("AliESDs.root", "RECREATE");
    ffile->SetCompressionLevel(2);
    if (!ffile->IsOpen()) {
      Abort("OpenESDFile", TSelector::kAbortProcess);
      return;
    }
  }
  else {
    AliInfo(Form("Opening output PROOF file: %s/%s",
		 outProofFile->GetDir(), outProofFile->GetFileName()));
    if (!(ffile = outProofFile->OpenFile("RECREATE"))) {
      Abort(Form("Problems opening output PROOF file: %s/%s",
		 outProofFile->GetDir(), outProofFile->GetFileName()),
	    TSelector::kAbortProcess);
      return;
    }
  }

  ftree = new TTree("esdTree", "Tree with ESD objects");
  fesd = new AliESDEvent();
  fesd->CreateStdContent();
  // add a so far non-std object to the ESD, this will
  // become part of the std content
  fesd->AddObject(new AliESDHLTDecision);

  fesd->WriteToTree(ftree);
  if (fWriteESDfriend) {
    ffileF = TFile::Open("AliESDfriends.root", "RECREATE");
    ftreeF = new TTree("esdFriendTree", "Tree with ESD Friend objects");
    fesdf  = new AliESDfriend();
    fesdf->SetESDIndicesStored(kTRUE); // new sparse format
    fesdfDummy  = new AliESDfriend();
    fesdf->SetESDIndicesStored(kTRUE);
    fesdfDummy->SetSkipBit(kTRUE); // empty copy
    ftreeF->Branch("ESDfriend.","AliESDfriend", &fesdf);
    fesd->AddObject(fesdf);
    ffile->cd();
  }
  ftree->GetUserInfo()->Add(fesd);

  if (fWriteHLTESD) fhlttree = new TTree("HLTesdTree", "Tree with HLT ESD objects");
  fhltesd = new AliESDEvent();
  fhltesd->CreateStdContent();
  // read the ESD template from CDB
  // HLT is allowed to put non-std content to its ESD, the non-std
  // objects need to be created before invocation of WriteToTree in
  // order to create all branches. Initialization is done from an
  // ESD layout template in CDB
  AliCDBManager* man = AliCDBManager::Instance();
  AliCDBEntry* hltESDConfig = man->Get("HLT/Calib/esdLayout");
  if(!hltESDConfig){
      AliError(Form("Error getting \"HLT/Calib/esdLayout\""));
      return;
  }
  AliESDEvent* pESDLayout=dynamic_cast<AliESDEvent*>(hltESDConfig->GetObject());
  if (pESDLayout) {
    // init all internal variables from the list of objects
    pESDLayout->GetStdContent();
    
    // copy content and create non-std objects
    *fhltesd=*pESDLayout;
    fhltesd->Reset();
  } else {
    AliError(Form("error setting hltEsd layout from \"HLT/Calib/esdLayout\": invalid object type"));
  }
  if (fWriteHLTESD) {
    fhltesd->WriteToTree(fhlttree);
    fhlttree->GetUserInfo()->Add(fhltesd);
  }
  else {
    AliInfo("HLT ESD tree writing is disabled");
  }
  
  ProcInfo_t procInfo;
  gSystem->GetProcInfo(&procInfo);
  AliInfo(Form("Current memory usage %ld %ld", procInfo.fMemResident, procInfo.fMemVirtual));
  
  // PID
  fESDpid = new AliESDpid();

  //QA
  //Initialize the QA and start of cycle 
  if (fRunQA || fRunGlobalQA) 
    InitQA() ; 

  //Initialize the Plane Efficiency framework
  if (fRunPlaneEff && !InitPlaneEff()) {
    Abort("InitPlaneEff", TSelector::kAbortProcess);
    return;
  }

  if (strcmp(gProgName,"alieve") == 0)
    fRunAliEVE = InitAliEVE();
  // If we have an analysis manager, connect the AliRecoInputHandler here  
  if (fAnalysis) {
    if (!dynamic_cast<AliRecoInputHandler*>(fAnalysis->GetInputEventHandler())) {
       AliError("Analysis manager used in reconstruction should use AliRecoInputHandler - \
                 \n  ->Replacing with AliRecoInputHandler instance.");
       delete fAnalysis->GetInputEventHandler();
    }
    // Set the event and other data pointers
    fRecoHandler = new AliRecoInputHandler();
//    fRecoHandler->Init(ftree, "LOCAL");
    fRecoHandler->SetEvent(fesd);
    fRecoHandler->SetESDfriend(fesdf);
    fRecoHandler->SetHLTEvent(fhltesd);
    if (fWriteHLTESD) {
      fRecoHandler->SetHLTTree(fhlttree);
    }
    fAnalysis->SetInputEventHandler(fRecoHandler);
    // Enter external loop mode
    fAnalysis->SetExternalLoop(kTRUE);
    // Initialize analysis
    fAnalysis->SlaveBegin(ftree);
    fAnalysis->StartAnalysis("local", (TTree*)0);
    // Connect ESD tree with the input container
    fAnalysis->GetCommonInputContainer()->SetData(ftree);
  }  
  //
  ProcessTriggerAliases();
  //
  if (!fWriteHLTESD) AliInfo("Writing of HLT ESD tree is disabled");
  else AliInfoF("%.2f%% of HLT ESD events will be stored",fFractionHLTESD*100.0);
  AliInfoF("%.2f%% of ESD friends  will be stored",fFractionFriends*100.0);

  return;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::Process(Long64_t entry)
{
  // run the reconstruction over a single entry
  // from the chain with raw data
  AliCodeTimerAuto("",0);

  TTree *currTree = fChain->GetTree();
  AliRawVEvent *event = NULL;
  currTree->SetBranchAddress("rawevent",&event);
  currTree->GetEntry(entry);
  fRawReader = new AliRawReaderRoot(event,entry);
  // check if process has enough resources 
  if (!HasEnoughResources(entry)) return kFALSE;
  fStatus = ProcessEvent(fRunLoader->GetNumberOfEvents());
  delete fRawReader;
  fRawReader = NULL;
  delete event;

  if (!fStatus) Abort("ProcessEvent",TSelector::kAbortFile);  
  CleanProcessedEvent();
  return fStatus;
}

//_____________________________________________________________________________
void AliReconstruction::Init(TTree *tree)
{
  // Implementation of TSelector::Init()
  // method
  if (tree == 0) {
    AliError("The input tree is not found!");
    return;
  }
  fChain = tree;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::ProcessEvent(Int_t iEvent)
{
  // run the reconstruction over a single event
  // The event loop is steered in Run method


  static Long_t oldMres=0;
  static Long_t oldMvir=0;
  static Float_t oldCPU=0;
  static Long_t aveDMres=0;
  static Long_t aveDMvir=0;
  static Float_t aveDCPU=0;

  AliCodeTimerAuto("",0);

  if (fRawReader){
    // Store DAQ detector pattern and attributes
    fesd->SetDAQDetectorPattern(fRawReader->GetDetectorPattern()[0]);
    fesd->SetDAQAttributes(fRawReader->GetAttributes()[2]);
    if (fesd->IsIncompleteDAQ() && fSkipIncompleteDAQ) {
      AliInfo(" ");
      AliInfoF("Abandoning incomplete event with RAW_ID=%d reconstruction: DAQ attr: 0x%08x",iEvent,fesd->GetDAQAttributes());
      AliInfoF("ATTENTION: %d events were already abandoned: next stored event will have ID:%d",fNAbandonedEv,iEvent-fNAbandonedEv);
      AliInfo(" ");
      fNAbandonedEv++;
      return kTRUE;
    }
  }

  iEvent -= fNAbandonedEv;


  AliSysInfo::AddStamp(Form("StartEv_%d",iEvent), 0,0,iEvent);

  if (iEvent >= fRunLoader->GetNumberOfEvents()) {
    fRunLoader->SetEventNumber(iEvent);
    if (fRawReader)
      fRunLoader->GetHeader()->Reset(fRawReader->GetRunNumber(), 
				     iEvent, iEvent);
    fRunLoader->TreeE()->Fill();

    if (fRawReader && fRawReader->UseAutoSaveESD())
      fRunLoader->TreeE()->AutoSave("SaveSelf");
  }

  if ((iEvent < fFirstEvent) || ((fLastEvent >= 0) && (iEvent > fLastEvent))) {
    return kTRUE;
  }


  fRunLoader->GetEvent(iEvent);

  // Fill Event-info object
  GetEventInfo();
  fRecoParam.SetEventSpecie(fRunInfo,fEventInfo,fListOfCosmicTriggers);
  
  ProcInfo_t procInfo;
  if(iEvent==fFirstEvent) {
    gSystem->GetProcInfo(&procInfo);
    oldMres=procInfo.fMemResident;
    oldMvir=procInfo.fMemVirtual;
    oldCPU=procInfo.fCpuUser+procInfo.fCpuSys;
  }
  AliInfo(Form("================================= Processing event %d of type %-10s ==================================", iEvent,fRecoParam.PrintEventSpecie()));
  fEventInfo.Print();


  AliSysInfo::AddStamp(Form("StartReco_%d",iEvent), 0,0,iEvent);
  // Set the reco-params
  {
    TString detStr = fLoadCDB;
    for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
      if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
      AliReconstructor *reconstructor = GetReconstructor(iDet);
      if (reconstructor && fRecoParam.GetDetRecoParamArray(iDet)) {
        const AliDetectorRecoParam *par = fRecoParam.GetDetRecoParam(iDet);
        reconstructor->SetRecoParam(par);
	reconstructor->GetPidSettings(fESDpid);
	reconstructor->SetEventInfo(&fEventInfo);
        if (fRunQA) {
          AliQAManager::QAManager()->SetEventInfo(&fEventInfo) ;
          AliQAManager::QAManager()->SetRecoParam(iDet, par) ; 
          if (par) AliQAManager::QAManager()->SetEventSpecie(AliRecoParam::Convert(par->GetEventSpecie())) ;
        }
      }
    }
    //
    if (fRunQA || fRunGlobalQA) AliQADataMaker::SetEventTrigClasses(fEventInfo.GetTriggerClasses()); // RS: select which histo clones are to be filled
    //
    if (fRunQA) {
      const AliDetectorRecoParam *grppar = fRecoParam.GetDetRecoParam(kNDetectors);
      AliQAManager::QAManager()->SetRecoParam(AliQAv1::kGLOBAL, grppar) ; 
      AliQAManager::QAManager()->SetEventSpecie(AliRecoParam::Convert(grppar->GetEventSpecie())) ;
    }
  }

    // QA on single raw 
  if (fRunQA && IsInTasks(AliQAv1::kRAWS)) {
    AliQAManager::QAManager()->SetEventSpecie(fRecoParam.GetEventSpecie()) ;
    AliQAManager::QAManager()->RunOneEvent(fRawReader) ;  
    AliSysInfo::AddStamp(Form("RawQA_%d",iEvent), 0,0,iEvent);
  }

    // fill Event header information from the RawEventHeader
    if (fRawReader){
      FillRawEventHeaderESD(fesd);
      FillRawEventHeaderESD(fhltesd);
    }
    else { // fill MC info
      FillMCEventHeaderESD(fesd);
      FillMCEventHeaderESD(fhltesd);
    }

    fesd->SetRunNumber(fRunLoader->GetHeader()->GetRun());
    fhltesd->SetRunNumber(fRunLoader->GetHeader()->GetRun());
    
    ((AliESDRun*)fesd->GetESDRun())->SetDetectorsInDAQ(fRunInfo->GetDetectorMask());
    ((AliESDRun*)fhltesd->GetESDRun())->SetDetectorsInDAQ(fRunInfo->GetDetectorMask());
    ((AliESDRun*)fesd->GetESDRun())->SetDetectorsInReco(AliDAQ::DetectorPatternOffline(fFillESD.Data()));
    ((AliESDRun*)fhltesd->GetESDRun())->SetDetectorsInReco(AliDAQ::DetectorPatternOffline(fFillESD.Data()));

    fesd->SetEventNumberInFile(fRunLoader->GetHeader()->GetEventNrInRun());
    fhltesd->SetEventNumberInFile(fRunLoader->GetHeader()->GetEventNrInRun());

    fesd->SetEventSpecie(fRecoParam.GetEventSpecie());
    fhltesd->SetEventSpecie(fRecoParam.GetEventSpecie());
    
    // Set magnetic field from the tracker
    fesd->SetMagneticField(AliTracker::GetBz());
    fhltesd->SetMagneticField(AliTracker::GetBz());
    //
    AliESDRun *esdRun,*esdRunH;
    esdRun  = (AliESDRun*)fesd->GetESDRun();
    esdRunH = (AliESDRun*)fhltesd->GetESDRun();
    esdRun->SetBeamEnergyIsSqrtSHalfGeV();
    esdRunH->SetBeamEnergyIsSqrtSHalfGeV();
    //
    for (int ib=2;ib--;) for (int it=2;it--;) {
	esdRun->SetMeanIntensity(ib,it, fBeamInt[ib][it]); 
	esdRunH->SetMeanIntensity(ib,it, fBeamInt[ib][it]); 
      }
    //
    fesd->SetBeamEnergy(fGRPData->GetBeamEnergy());
    fesd->SetBeamType(fGRPData->GetBeamType().Data());
    fesd->SetBeamParticle(fGRPData->GetSingleBeamType(0).Atoi(),0);
    fesd->SetBeamParticle(fGRPData->GetSingleBeamType(1).Atoi(),1);
    fhltesd->SetBeamEnergy(fGRPData->GetBeamEnergy());
    fhltesd->SetBeamType(fGRPData->GetBeamType().Data());
    fhltesd->SetBeamParticle(fGRPData->GetSingleBeamType(0).Atoi(),0);
    fhltesd->SetBeamParticle(fGRPData->GetSingleBeamType(1).Atoi(),1);
    //
    AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
    if (fld) { // set info needed for field initialization
      fesd->SetCurrentL3(fld->GetCurrentSol());
      fesd->SetCurrentDip(fld->GetCurrentDip());
      fesd->SetUniformBMap(fld->IsUniform());
      fesd->SetBInfoStored();
      //
      fhltesd->SetCurrentL3(fld->GetCurrentSol());
      fhltesd->SetCurrentDip(fld->GetCurrentDip());
      fhltesd->SetUniformBMap(fld->IsUniform());
      fhltesd->SetBInfoStored();
    }

    // set CTPstart timestamp of the run
    if (fRunScalers&&fRunScalers->GetScalersRecords()->GetEntriesFast()) {
      esdRun->SetCTPStart(fRunScalers->GetScalersRecord(0)->GetTimeStamp());
      esdRunH->SetCTPStart(fRunScalers->GetScalersRecord(0)->GetTimeStamp());
    }

    for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
      if (fTracker[iDet]) { // some trackers need details about the run/time
	fTracker[iDet]->SetTimeStamp(fesd->GetTimeStamp());
	fTracker[iDet]->SetRunNumber(fesd->GetRunNumber());
      }
      //RS also some reconstructructors may need the time stamp
      if (fReconstructor[iDet]) fReconstructor[iDet]->SetTimeStamp(fesd->GetTimeStamp());
    }

    //
    // run full HLT reconstruction first
    //
    {
      TString detectors=fRunLocalReconstruction;
      if (IsSelected("HLT", detectors) &&
	  !RunLocalEventReconstruction("HLT")) {
	if (fStopOnError) {CleanUp(); return kFALSE;}
      }
      detectors=fFillESD;
      // run HLT on hltesd
      if (IsSelected("HLT", detectors) &&
	  !FillESD(fhltesd, "HLT")) {
	if (fStopOnError) {CleanUp(); return kFALSE;}
      }
    }

    // local single event reconstruction
    if (!fRunLocalReconstruction.IsNull()) {
      TString detectors=fRunLocalReconstruction;
      // the logic for selection and correct sequence of reconstruction relies on the
      // full list of detectors. Keyword 'ALL' should have been replaced at this point.
      if (detectors.Contains("ALL")) {
	AliFatal("Keyword 'ALL' needs to be replaced by the full list of detectors in "
		 "fRunLocalReconstruction. This should have been done by the framework");
      }
      detectors.ReplaceAll("HLT", "");
      if (!RunLocalEventReconstruction(detectors)) {
        if (fStopOnError) {
          CleanUp(); 
          return kFALSE;
        }
      }
    }

  
    //
    // Set most probable pt, for B=0 tracking
    // Get the global reco-params. They are atposition 16 inside the array of detectors in fRecoParam
    const AliGRPRecoParam *grpRecoParam = dynamic_cast<const AliGRPRecoParam*>(fRecoParam.GetDetRecoParam(kNDetectors));
    if (grpRecoParam) AliExternalTrackParam::SetMostProbablePt(grpRecoParam->GetMostProbablePt());
    
    // Fill raw-data error log into the ESD
    if (fRawReader) FillRawDataErrorLog(iEvent,fesd);

    AliSysInfo::AddStamp(Form("FillHeadErrs_%d",iEvent), 0,0,iEvent);

    // vertex finder
    if (fRunVertexFinder) {
      if (!RunVertexFinder(fesd)) {
	if (fStopOnError) {CleanUp(); return kFALSE;}
      }
      AliSysInfo::AddStamp(Form("VtxFinder_%d",iEvent), 0,0,iEvent);
    }

    // For Plane Efficiency: run the SPD trackleter
    if (fRunPlaneEff && fSPDTrackleter) {
      if (!RunSPDTrackleting(fesd)) {
        if (fStopOnError) {CleanUp(); return kFALSE;}
      }
      AliSysInfo::AddStamp(Form("TrackletEff_%d",iEvent), 0,0,iEvent);
    }

    // Muon tracking
    if (!fRunTracking.IsNull()) {
      if (fRunMuonTracking) {
	if (!RunMuonTracking(fesd)) {
	  if (fStopOnError) {CleanUp(); return kFALSE;}
	}
      }
      AliSysInfo::AddStamp(Form("TrackingMUON_%d",iEvent), 0,0,iEvent);      
    }

    //---------------- AU From here...

    // MFT tracking of MUON tracks
    if (!fRunTracking.IsNull()) {
      if (fRunMFTTrackingMU && fRunMuonTracking) {
	if (!RunMFTTrackingMU(fesd)) {
	  if (fStopOnError) {CleanUp(); return kFALSE;}
	}
      }
      AliSysInfo::AddStamp(Form("TrackingMFT_MUON_%d",iEvent), 0,0,iEvent);      
    }

    //---------------- ...to here

    // barrel tracking
    if (!fRunTracking.IsNull()) {
      if (!RunTracking(fesd,*fESDpid)) {
	if (fStopOnError) {CleanUp(); return kFALSE;}
      }
    }

    // fill ESD
    if (!fFillESD.IsNull()) {
      TString detectors=fFillESD;
      // the logic for selection and correct sequence of reconstruction relies on the
      // full list of detectors. Keyword 'ALL' should have been replaced at this point.
      if (detectors.Contains("ALL")) {
	AliFatal("Keyword 'ALL' needs to be replaced by the full list of detectors in "
		 "fFillESD. This should have been done by the framework");
      }
      // remove HLT as this has been executed at the beginning of the event reconstruction
      detectors.ReplaceAll("HLT", "");
      if (!FillESD(fesd, detectors)) {
	if (fStopOnError) {CleanUp(); return kFALSE;}
      }
    }

    ffile->cd();

    //
    // Propagate track to the beam pipe  (if not already done by ITS)
    //
    const Int_t ntracks = fesd->GetNumberOfTracks();
    const Double_t kRadius  = 2.8; //something less than the beam pipe radius

    TObjArray trkArray;
    UShort_t selectedIdx[ntracks];

    for (Int_t itrack=0; itrack<ntracks; itrack++){
      const Double_t kMaxStep = 5;   //max step over the material
      Bool_t ok;

      AliESDtrack *track = fesd->GetTrack(itrack);
      if (!track) continue;

      AliExternalTrackParam *tpcTrack =
           (AliExternalTrackParam *)track->GetTPCInnerParam();
      ok = kFALSE;
      if (tpcTrack)
	ok = AliTracker::
	  PropagateTrackToBxByBz(tpcTrack,kRadius,track->GetMassForTracking(),kMaxStep,kFALSE);

      if (ok) {
	Int_t n=trkArray.GetEntriesFast();
        selectedIdx[n]=track->GetID();
        trkArray.AddLast(tpcTrack);
      }

      //Tracks refitted by ITS should already be at the SPD vertex
      if (track->IsOn(AliESDtrack::kITSrefit)) continue;

      AliTracker::
         PropagateTrackToBxByBz(track,kRadius,track->GetMassForTracking(),kMaxStep,kFALSE);
      Double_t x[3]; track->GetXYZ(x);
      Double_t b[3]; AliTracker::GetBxByBz(x,b);
      track->RelateToVertexBxByBz(fesd->GetPrimaryVertexSPD(), b, kVeryBig);

    }
    AliSysInfo::AddStamp(Form("RelToSPDVtx_%d",iEvent), 0,0,iEvent);      
    //
    // Improve the reconstructed primary vertex position using the tracks
    //
    Bool_t runVertexFinderTracks = fRunVertexFinderTracks;
    if(fesd->GetPrimaryVertexSPD()) {
      TString vtitle = fesd->GetPrimaryVertexSPD()->GetTitle();
      if(vtitle.Contains("cosmics")) {
	runVertexFinderTracks=kFALSE;
      }
    }

    if (runVertexFinderTracks) {
       // TPC + ITS primary vertex
       ftVertexer->SetITSMode();
       ftVertexer->SetConstraintOff();
       // get cuts for vertexer from AliGRPRecoParam
       Bool_t constrSPD=kFALSE;
       if (grpRecoParam) {
	 Int_t nCutsVertexer = grpRecoParam->GetVertexerTracksNCuts();
	 Double_t *cutsVertexer = new Double_t[nCutsVertexer];
	 grpRecoParam->GetVertexerTracksCutsITS(cutsVertexer,nCutsVertexer);
	 ftVertexer->SetCuts(cutsVertexer,nCutsVertexer);
	 delete [] cutsVertexer; cutsVertexer = NULL; 
	 if(grpRecoParam->GetVertexerTracksConstraintITS()) { 
	   if(fDiamondProfile && fDiamondProfile->GetXRes()<kRadius){
	     ftVertexer->SetVtxStart(fDiamondProfile); // apply constraint only if sigmax is smaller than the beam pipe radius 
	   }else{
	     if(fDiamondProfileSPD && fDiamondProfileSPD->GetXRes()<kRadius){
	       ftVertexer->SetVtxStart(fDiamondProfileSPD);
	       constrSPD=kTRUE;
	     }
	   }
	 } 
       }
       Bool_t usePureSA = ftVertexer->GetITSpureSA();
       if (!fesd->GetNumberOfTPCClusters()) {
	 AliInfo("No TPC clusters: VertexerTracks will run with ITSpureSA tracks");
	 ftVertexer->SetITSpureSA(kTRUE);
       }
       AliESDVertex *pvtx=ftVertexer->FindPrimaryVertex(fesd);
       ftVertexer->SetITSpureSA(usePureSA); // restore use settings
       if (pvtx) {
	 if(constrSPD){
	   TString title=pvtx->GetTitle();
	   title.Append("SPD");
	   pvtx->SetTitle(title);
	 }
          if (pvtx->GetStatus()) {
             fesd->SetPrimaryVertexTracks(pvtx);
             for (Int_t i=0; i<ntracks; i++) {
	         AliESDtrack *t = fesd->GetTrack(i);
                 Double_t x[3]; t->GetXYZ(x);
                 Double_t b[3]; AliTracker::GetBxByBz(x,b);
                 t->RelateToVertexBxByBz(pvtx, b, kVeryBig);
             } 
          }
	  delete pvtx; pvtx=NULL;
       }
       AliSysInfo::AddStamp(Form("VtxTrk_%d",iEvent), 0,0,iEvent);      

       // TPC-only primary vertex
       ftVertexer->SetTPCMode();
       ftVertexer->SetConstraintOff();
       // get cuts for vertexer from AliGRPRecoParam
       if (grpRecoParam) {
	 Int_t nCutsVertexer = grpRecoParam->GetVertexerTracksNCuts();
	 Double_t *cutsVertexer = new Double_t[nCutsVertexer];
	 grpRecoParam->GetVertexerTracksCutsTPC(cutsVertexer,nCutsVertexer);
	 ftVertexer->SetCuts(cutsVertexer,nCutsVertexer);
	 delete [] cutsVertexer; cutsVertexer = NULL; 
	 if(fDiamondProfileTPC && grpRecoParam->GetVertexerTracksConstraintTPC()) { 
	   if(fDiamondProfileTPC->GetXRes()<kRadius) ftVertexer->SetVtxStart(fDiamondProfileTPC); // apply constraint only if sigmax is smaller than the beam pipe radius 
	 } 
       }
       pvtx=ftVertexer->FindPrimaryVertex(&trkArray,selectedIdx);
       if (pvtx) {
          if (pvtx->GetStatus()) {
             fesd->SetPrimaryVertexTPC(pvtx);
             for (Int_t i=0; i<ntracks; i++) {
	         AliESDtrack *t = fesd->GetTrack(i);
                 Double_t x[3]; t->GetXYZ(x);
                 Double_t b[3]; AliTracker::GetBxByBz(x,b);
                 t->RelateToVertexTPCBxByBz(pvtx, b, kVeryBig);
             } 
          }
	  delete pvtx; pvtx=NULL;
       }
       AliSysInfo::AddStamp(Form("VtxTPC_%d",iEvent), 0,0,iEvent);      

    }
    
    if(fDiamondProfile && fDiamondProfile->GetXRes()<kRadius) fesd->SetDiamond(fDiamondProfile);
    else fesd->SetDiamond(fDiamondProfileSPD);

    if (fRunV0Finder) {
       // V0 finding
       AliV0vertexer vtxer;
       // get cuts for V0vertexer from AliGRPRecoParam
       if (grpRecoParam) {
	 Int_t nCutsV0vertexer = grpRecoParam->GetVertexerV0NCuts();
	 Double_t cutsV0vertexer[nCutsV0vertexer];
	 grpRecoParam->GetVertexerV0Cuts(cutsV0vertexer);
	 vtxer.SetCuts(cutsV0vertexer);
       }
       vtxer.Tracks2V0vertices(fesd);
       AliSysInfo::AddStamp(Form("V0Finder_%d",iEvent), 0,0,iEvent); 

       if (fRunCascadeFinder) {
          // Cascade finding
          AliCascadeVertexer cvtxer;
	  // get cuts for CascadeVertexer from AliGRPRecoParam
	  if (grpRecoParam) {
	    Int_t nCutsCascadeVertexer = grpRecoParam->GetVertexerCascadeNCuts();
	    Double_t cutsCascadeVertexer[nCutsCascadeVertexer];
	    grpRecoParam->GetVertexerCascadeCuts(cutsCascadeVertexer);
	    cvtxer.SetCuts(cutsCascadeVertexer);
	  }
          cvtxer.V0sTracks2CascadeVertices(fesd);
	  AliSysInfo::AddStamp(Form("CascadeFinder_%d",iEvent), 0,0,iEvent); 
       }
    }

    // AdC+FN
    if (fReconstructor[3])
      GetReconstructor(3)->FillEventTimeWithTOF(fesd,fESDpid);

    // combined PID
    //    fESDpid->MakePID(fesd);

    if (fFillTriggerESD) {
      if (!FillTriggerESD(fesd)) {
	if (fStopOnError) {CleanUp(); return kFALSE;}
      }
    }
    // Always fill scalers
    if (!FillTriggerScalers(fesd)) {
       if (fStopOnError) {CleanUp(); return kFALSE;}
    }

    AliSysInfo::AddStamp(Form("FillVaria_%d",iEvent), 0,0,iEvent); 

    // write ESD
    UInt_t specie = fesd->GetEventSpecie();
    Bool_t keepAll = (specie==AliRecoParam::kCosmic || specie==AliRecoParam::kCalib);
    if (fCleanESD && (!keepAll) ) {
      CleanESD(fesd);
      AliSysInfo::AddStamp(Form("CleanESD_%d",iEvent), 0,0,iEvent); 
    }
    // 
    // RS run updated trackleter: since we want to mark the clusters used by tracks and also mark the 
    // tracks interpreted as primary, this step should be done in the very end, when full 
    // ESD info is available (particulalry, V0s)
    // vertex finder
    if (fRunMultFinder) {
      if (!RunMultFinder(fesd)) {
	if (fStopOnError) {CleanUp(); return kFALSE;}
      }
      AliSysInfo::AddStamp(Form("MultFinder_%d",iEvent), 0,0,iEvent); 
    }

  if (fRunQA && IsInTasks(AliQAv1::kESDS)) {
    AliQAManager::QAManager()->SetEventSpecie(fRecoParam.GetEventSpecie()) ;
    AliQAManager::QAManager()->RunOneEvent(fesd, fhltesd) ; 
    AliSysInfo::AddStamp(Form("RunQA_%d",iEvent), 0,0,iEvent); 
  }
  if (fRunGlobalQA) {
    AliQADataMaker *qadm = AliQAManager::QAManager()->GetQADataMaker(AliQAv1::kGLOBAL);
    if (qadm)
      qadm->SetEventSpecie(fRecoParam.GetEventSpecie()) ;
    if (qadm && IsInTasks(AliQAv1::kESDS))
      qadm->Exec(AliQAv1::kESDS, fesd);
    AliSysInfo::AddStamp(Form("RunGlobQA_%d",iEvent), 0,0,iEvent);     
  }

  // copy HLT decision from HLTesd to esd
  // the most relevant information is stored in a reduced container in the esd,
  // while the full information can be found in the HLTesd
  TObject* pHLTSrc=fhltesd->FindListObject(AliESDHLTDecision::Name());
  TObject* pHLTTgt=fesd->FindListObject(AliESDHLTDecision::Name());
  if (pHLTSrc && pHLTTgt) {
    pHLTSrc->Copy(*pHLTTgt);
  }
  //
  // Perform analysis of this event if requested
  // RS: Should be done before WriteESDfriend, since the latter may clean the esdfriend
  if (fAnalysis) {
    fRecoHandler->BeginEvent(iEvent);
    fAnalysis->ExecAnalysis();
    fRecoHandler->FinishEvent();
    AliSysInfo::AddStamp(Form("Analysis_%d",iEvent), 0,0,iEvent);     
  }  
  //
  if (fWriteThisFriend) {
    fesd->GetESDfriend(fesdf);
    AliSysInfo::AddStamp(Form("CreateFriend_%d",iEvent), 0,0,iEvent);     
  }
  //
  Long64_t nbf;
  nbf = ftree->Fill();
  if (fTreeBuffSize>0 && ftree->GetAutoFlush()<0 && (fMemCountESD += nbf)>fTreeBuffSize ) { // default limit is still not reached
    nbf = ftree->GetZipBytes();
    if (nbf>0) nbf = -nbf;
    else       nbf = ftree->GetEntries();
    ftree->SetAutoFlush(nbf);
    AliInfo(Form("Calling ftree->SetAutoFlush(%lld) | W:%lld T:%lld Z:%lld",
		 nbf,fMemCountESD,ftree->GetTotBytes(),ftree->GetZipBytes()));        
  }
  AliSysInfo::AddStamp(Form("ESDFill_%d",iEvent), 0,0,iEvent);     
  //
  if (fWriteESDfriend) {
    WriteESDfriend();
    AliSysInfo::AddStamp(Form("WriteFriend_%d",iEvent), 0,0,iEvent);     
  }
  //
  // Auto-save the ESD tree in case of prompt reco @P2
  if (fRawReader && fRawReader->UseAutoSaveESD()) {
    ftree->AutoSave("SaveSelf");
    if (fWriteESDfriend) ftreeF->AutoSave("SaveSelf");
  }
  // write HLT ESD
  if (fFractionHLTESD<1.0 && fhltesd->GetESDRun()->TestBit(AliESDRun::kESDDownscaledOnline)) {
    AliInfo("HLT ESD online prescaling flag is ON: disabling offline prescaling");
    fFractionHLTESD = 1.0;
  }
  //
  if (fWriteHLTESD) {
    if (gRandom->Rndm()>fFractionHLTESD) {
      AliInfo("HLT ESD for this event will be empty");
      fhltesd->Reset();
    } 
    nbf = fhlttree->Fill();
    if (fTreeBuffSize>0 && fhlttree->GetAutoFlush()<0 && (fMemCountESDHLT += nbf)>fTreeBuffSize ) { // default limit is still not reached
      nbf = fhlttree->GetZipBytes();
      if (nbf>0) nbf = -nbf;
      else       nbf = fhlttree->GetEntries();
      fhlttree->SetAutoFlush(nbf);
      AliInfo(Form("Calling fhlttree->SetAutoFlush(%lld) | W:%lld T:%lld Z:%lld",
		   nbf,fMemCountESDHLT,fhlttree->GetTotBytes(),fhlttree->GetZipBytes()));        
    }
  }

  gSystem->GetProcInfo(&procInfo);
  Long_t dMres=(procInfo.fMemResident-oldMres)/1024;
  Long_t dMvir=(procInfo.fMemVirtual-oldMvir)/1024;
  Float_t dCPU=procInfo.fCpuUser+procInfo.fCpuSys-oldCPU;
  aveDMres+=(dMres-aveDMres)/(iEvent-fFirstEvent+1);
  aveDMvir+=(dMvir-aveDMvir)/(iEvent-fFirstEvent+1);
  aveDCPU+=(dCPU-aveDCPU)/(iEvent-fFirstEvent+1);
  AliInfo(Form("======================= End Event %d: Res %ld(%3ld <%3ld>) Vir %ld(%3ld <%3ld>) CPU %5.2f <%5.2f> ===================",
	       iEvent, procInfo.fMemResident/1024, dMres, aveDMres, procInfo.fMemVirtual/1024, dMvir, aveDMvir, dCPU, aveDCPU));
  oldMres=procInfo.fMemResident;
  oldMvir=procInfo.fMemVirtual;
  oldCPU=procInfo.fCpuUser+procInfo.fCpuSys;
  
  if (fTracker[AliQAv1::kTPC]) {
    fTracker[AliQAv1::kTPC]->UnloadClusters();
    AliSysInfo::AddStamp(Form("TUnloadCluster%s_%d",fgkDetectorName[AliQAv1::kTPC],iEvent), AliQAv1::kTPC,4, iEvent);
  }
  if (fLoader[AliQAv1::kTPC]) {
    fLoader[AliQAv1::kTPC]->UnloadRecPoints();
    AliSysInfo::AddStamp(Form("RUnloadCluster%s_%d",fgkDetectorName[AliQAv1::kTPC],iEvent), AliQAv1::kTPC,5, iEvent);
  }
    
  return kTRUE;
}

void AliReconstruction::CleanProcessedEvent()
{
    //
    //    if (fWriteESDfriend) ResetFriends();
    ResetFriends(); // RS: Important: first reset friends to be sure objects owned by trackers are released

    fesd->Reset();
    fhltesd->Reset();

    for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
      if (fReconstructor[iDet]) fReconstructor[iDet]->FinishEvent();
    }
 
    AliInfo("======================= End Event ===================");
    

    for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
      if (fReconstructor[iDet]) {
	fReconstructor[iDet]->SetRecoParam(NULL);
	fReconstructor[iDet]->SetEventInfo(NULL);
      }
      if (fTracker[iDet]) fTracker[iDet]->SetEventInfo(NULL);
    }
	
  if (fRunQA || fRunGlobalQA) 
    AliQAManager::QAManager()->Increment() ; 

  DeleteRecPoints(fDeleteRecPoints);
  DeleteDigits(fDeleteDigits);
  //

}

//_____________________________________________________________________________
void AliReconstruction::SlaveTerminate()
{
  // Finalize the run on the slave side
  // Called after the exit
  // from the event loop
  AliCodeTimerAuto("",0);
  // If analysis was done during reconstruction, we need to call SlaveTerminate for it
  if (fAnalysis) {
     fAnalysis->PackOutput(fOutput);
     fAnalysis->SetSkipTerminate(kTRUE);
     fAnalysis->Terminate();
  }   

  if (fIsNewRunLoader) { // galice.root didn't exist
    fRunLoader->WriteHeader("OVERWRITE");
    fRunLoader->WriteTrigger("OVERWRITE");
    fRunLoader->CdGAFile();
    fRunLoader->Write(0, TObject::kOverwrite);
  }

  TMap *cdbMapCopy = new TMap();
  cdbMapCopy->SetName("cdbMap");
  TList *cdbListCopy = new TList();
  cdbListCopy->SetName("cdbList");
  // create map/list accounting for eventual snapshot
  AliCDBManager::Instance()->CreateMapListCopy(*cdbMapCopy,*cdbListCopy);
  
  ftree->GetUserInfo()->Add(cdbMapCopy);	 
  ftree->GetUserInfo()->Add(cdbListCopy);

   // Add the AliRoot version that created this file
   TString sVersion("aliroot ");
   sVersion += ALIROOT_VERSION;
   sVersion += ":";
   sVersion += ALIROOT_REVISION;
   sVersion += "; root ";
#ifdef ROOT_SVN_BRANCH
   sVersion += ROOT_SVN_BRANCH;
#elif defined(ROOT_GIT_BRANCH)
   sVersion += ROOT_GIT_BRANCH;
#else
   sVersion += "?";
#endif
   sVersion += ":";
#ifdef ROOT_SVN_REVSION
   sVersion += ROOT_SVN_REVISION;
#elif defined(ROOT_GIT_COMMIT)
   sVersion += ROOT_GIT_COMMIT;
#else 
   sVersion += "?";
#endif
   sVersion += "; metadata ";
   sVersion += getenv("PRODUCTION_METADATA");
		    

   TNamed * alirootVersion = new TNamed("alirootVersion",sVersion.Data());
   ftree->GetUserInfo()->Add(alirootVersion); // The list becomes owner of alirootVersion

  ffile->cd();

  // we want to have only one tree version number
  ftree->Write(ftree->GetName(),TObject::kOverwrite);
  if (fWriteHLTESD) fhlttree->Write(fhlttree->GetName(),TObject::kOverwrite);

  if (fWriteESDfriend) {
    ffileF->cd();
    ftreeF->Write(ftreeF->GetName(),TObject::kOverwrite);
  }

// Finish with Plane Efficiency evaluation: before of CleanUp !!!
  if (fRunPlaneEff && !FinishPlaneEff()) {
   AliWarning("Finish PlaneEff evaluation failed");
  }

  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if (fReconstructor[iDet]) fReconstructor[iDet]->Terminate();
  }
  // End of cycle for the in-loop  

  if (fRunQA || fRunGlobalQA) {
    AliQAManager::QAManager()->EndOfCycle() ;
    if (fInput &&
	!fProofOutputLocation.IsNull() &&
	fProofOutputArchive.IsNull() &&
	!fProofOutputDataset) {
      TString qaOutputFile(Form("%sMerged.%s.Data.root",
				fProofOutputLocation.Data(),
				AliQAv1::GetQADataFileName()));
      TProofOutputFile *qaProofFile = new TProofOutputFile(Form("Merged.%s.Data.root",
								AliQAv1::GetQADataFileName()));
      qaProofFile->SetOutputFileName(qaOutputFile.Data());
      if (AliDebugLevel() > 0) qaProofFile->Dump();
      fOutput->Add(qaProofFile);
      MergeQA(qaProofFile->GetFileName());
    }
    else {
      MergeQA();
    }
  }

  gROOT->cd();
  CleanUp();

  if (fInput) {
    if (!fProofOutputFileName.IsNull() &&
	!fProofOutputLocation.IsNull() &&
	fProofOutputDataset &&
	!fProofOutputArchive.IsNull()) {
      TProofOutputFile *zipProofFile = new TProofOutputFile(fProofOutputFileName.Data(),
							    "DROV",
							    fProofOutputLocation.Data());
      if (AliDebugLevel() > 0) zipProofFile->Dump();
      fOutput->Add(zipProofFile);
      TString fileList(fProofOutputArchive.Data());
      fileList.ReplaceAll(","," ");
      TString command;
#if ROOT_SVN_REVISION >= 30174
      command.Form("zip -n root %s/%s %s",zipProofFile->GetDir(kTRUE),zipProofFile->GetFileName(),fileList.Data());
#else
      command.Form("zip -n root %s/%s %s",zipProofFile->GetDir(),zipProofFile->GetFileName(),fileList.Data());
#endif
      AliInfo(Form("Executing: %s",command.Data()));
      gSystem->Exec(command.Data());
    }
  }
}
    
//_____________________________________________________________________________
void AliReconstruction::Terminate()
{
  // Create tags for the events in the ESD tree (the ESD tree is always present)
  // In case of empty events the tags will contain dummy values
  AliCodeTimerAuto("",0);

  // Do not call the ESD tag creator in case of PROOF-based reconstruction
  if (!fInput) {
    AliESDTagCreator *esdtagCreator = new AliESDTagCreator();
    esdtagCreator->CreateESDTags(fFirstEvent,fLastEvent,fGRPData, AliQAv1::Instance()->GetQA(), AliQAv1::Instance()->GetEventSpecies(), AliQAv1::kNDET, AliRecoParam::kNSpecies);
    delete esdtagCreator;
  }

  // Cleanup of CDB manager: cache and active storages!
  AliCDBManager::Instance()->ClearCache();
}

//_____________________________________________________________________________
Bool_t AliReconstruction::RunLocalEventReconstruction(const TString& detectors)
{
// run the local reconstruction

  static Int_t eventNr=0;
  AliCodeTimerAuto("",0)

  TString detStr = detectors;
  // execute HLT reconstruction first since other detector reconstruction
  // might depend on HLT data
  // key 'HLT' is removed from detStr by IsSelected
  if (IsSelected("HLT", detStr)) {
    AliReconstructor* reconstructor = GetReconstructor(kNDetectors-1);
    if (reconstructor) {
      // there is no AliLoader for HLT, see
      // https://savannah.cern.ch/bugs/?35473
      AliInfo("running reconstruction for HLT");
      if (fRawReader) {
	AliInfo("reconstructor->Reconstruct(fRawReader, NULL)");
        reconstructor->Reconstruct(fRawReader, NULL);
      } 
      else {
	AliInfo("reconstructor->Reconstruct(dummy, NULL)");
        TTree* dummy=NULL;
        reconstructor->Reconstruct(dummy, NULL);
      }
    }
    AliSysInfo::AddStamp(Form("LRecHLT_%d",eventNr), -1,1,eventNr);
  }

  AliInfo(Form("kNDetectors = %d",kNDetectors));

  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
    AliDebug(1, Form("Detector: %s", fgkDetectorName[iDet]));
    AliReconstructor* reconstructor = GetReconstructor(iDet);
    if (!reconstructor) continue;
    AliLoader* loader = fLoader[iDet];
    if (!loader) {
      AliWarning(Form("No loader is defined for %s!",fgkDetectorName[iDet]));
      continue;
    }
    // conversion of digits
    if (fRawReader && reconstructor->HasDigitConversion()) {
      AliInfo(Form("converting raw data digits into root objects for %s", 
		   fgkDetectorName[iDet]));
//      AliCodeTimerAuto(Form("converting raw data digits into root objects for %s", 
//                            fgkDetectorName[iDet]),0);
      loader->LoadDigits("update");
      loader->CleanDigits();
      loader->MakeDigitsContainer();
      TTree* digitsTree = loader->TreeD();
      reconstructor->ConvertDigits(fRawReader, digitsTree);
      loader->WriteDigits("OVERWRITE");
      loader->UnloadDigits();
    }
    // local reconstruction
    AliInfo(Form("running reconstruction for %s", fgkDetectorName[iDet]));
    //AliCodeTimerAuto(Form("running reconstruction for %s", fgkDetectorName[iDet]),0);
    AliDebug(1, "Loading Rec Points");
    loader->LoadRecPoints("update");
    AliDebug(1, "Cleaning Rec Points");
    loader->CleanRecPoints();
    AliDebug(1, "Making Rec Points Container");
    loader->MakeRecPointsContainer();
    TTree* clustersTree = loader->TreeR();
    if (fRawReader && !reconstructor->HasDigitConversion()) {
      reconstructor->Reconstruct(fRawReader, clustersTree);
    } 
    else {
      AliDebug(1, "Loading Digits");
      loader->LoadDigits("read");
      TTree* digitsTree = loader->TreeD();
      AliDebug(1, Form("Digits Tree = %p",digitsTree));
      if (!digitsTree) {
        AliError(Form("Can't get the %s digits tree", fgkDetectorName[iDet]));
        if (fStopOnError) 
          return kFALSE;
      } 
      else {
	AliDebug(1, "Digits -> Clusters");
        reconstructor->Reconstruct(digitsTree, clustersTree);
        if (fRunQA && IsInTasks(AliQAv1::kDIGITSR)) {
          AliQAManager::QAManager()->SetEventSpecie(fRecoParam.GetEventSpecie()) ;
          AliQAManager::QAManager()->RunOneEventInOneDetector(iDet, digitsTree) ; 
        }
      }
      loader->UnloadDigits();
    }
    if (fRunQA && IsInTasks(AliQAv1::kRECPOINTS)) {
      AliQAManager::QAManager()->SetEventSpecie(fRecoParam.GetEventSpecie()) ;
      AliQAManager::QAManager()->RunOneEventInOneDetector(iDet, clustersTree) ; 
    }
    loader->WriteRecPoints("OVERWRITE");
    loader->UnloadRecPoints();
    AliSysInfo::AddStamp(Form("LRec%s_%d",fgkDetectorName[iDet],eventNr), iDet,1,eventNr);
  }
  if (!IsSelected("CTP", detStr)) AliDebug(10,"No CTP");
  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    AliError(Form("the following detectors were not found: %s",
                  detStr.Data()));
    if (fStopOnError) 
      return kFALSE;
  }
  eventNr++;
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliReconstruction::RunSPDTrackleting(AliESDEvent*& esd)
{
// run the SPD trackleting (for SPD efficiency purpouses)

  AliCodeTimerAuto("",0)

  Double_t vtxPos[3] = {0, 0, 0};
  Double_t vtxErr[3] = {0.0, 0.0, 0.0};
/*
  TArrayF m
/
cVertex(3);
  // if(MC)
  if (fRunLoader->GetHeader() && fRunLoader->GetHeader()->GenEventHeader()) {
    fRunLoader->GetHeader()->GenEventHeader()->PrimaryVertex(mcVertex);
    for (Int_t i = 0; i < 3; i++) vtxPos[i] = mcVertex[i];
  }
*/
  const AliESDVertex *vertex = esd->GetVertex();
  if(!vertex){
    AliWarning("Vertex not found");
    return kFALSE;
  }
  vertex->GetXYZ(vtxPos);
  vertex->GetSigmaXYZ(vtxErr);
  if (fSPDTrackleter) {
    AliInfo("running the SPD Trackleter for Plane Efficiency Evaluation");

    // load clusters
    fLoader[0]->LoadRecPoints("read");
    TTree* tree = fLoader[0]->TreeR();
    if (!tree) {
      AliError("Can't get the ITS cluster tree");
      return kFALSE;
    }
    fSPDTrackleter->LoadClusters(tree);
    fSPDTrackleter->SetVertex(vtxPos, vtxErr);
    // run trackleting
    if (fSPDTrackleter->Clusters2Tracks(esd) != 0) {
      AliWarning("AliITSTrackleterSPDEff Clusters2Tracks failed");
     // fLoader[0]->UnloadRecPoints();
      return kFALSE;
    }
//fSPDTrackleter->UnloadRecPoints();
  } else {
    AliWarning("SPDTrackleter not available");
    return kFALSE;
  }
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::RunVertexFinder(AliESDEvent*& esd)
{
// run the barrel tracking

  AliCodeTimerAuto("",0)

  AliVertexer *vertexer = CreateVertexer();
  if (!vertexer) return kFALSE;

  AliInfo(Form("running the ITS vertex finder: %s",vertexer->ClassName()));
  AliESDVertex* vertex = NULL;
  if (fLoader[0]) {
    fLoader[0]->LoadRecPoints();
    TTree* cltree = fLoader[0]->TreeR();
    if (cltree) {
      if(fDiamondProfileSPD) vertexer->SetVtxStart(fDiamondProfileSPD);
      vertex = vertexer->FindVertexForCurrentEvent(cltree);
    }
    else {
      AliError("Can't get the ITS cluster tree");
    }
    fLoader[0]->UnloadRecPoints();
  }
  else {
    AliError("Can't get the ITS loader");
  }
  if(!vertex){
    AliWarning("Vertex not found");
    vertex = new AliESDVertex();
    vertex->SetName("default");
  }
  else {
    vertex->SetName("reconstructed");
  }

  Double_t vtxPos[3];
  Double_t vtxErr[3];
  vertex->GetXYZ(vtxPos);
  vertex->GetSigmaXYZ(vtxErr);

  esd->SetPrimaryVertexSPD(vertex);
  AliESDVertex *vpileup = NULL;
  Int_t novertices = 0;
  vpileup = vertexer->GetAllVertices(novertices);
  if(novertices>1){
    for (Int_t kk=1; kk<novertices; kk++)esd->AddPileupVertexSPD(&vpileup[kk]);
  }
  /*
  // if SPD multiplicity has been determined, it is stored in the ESD
  AliMultiplicity *mult = vertexer->GetMultiplicity();
  if(mult)esd->SetMultiplicity(mult);
  */
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if (fTracker[iDet]) fTracker[iDet]->SetVertex(vtxPos, vtxErr);
  }  
  delete vertex;

  delete vertexer;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::RunMultFinder(AliESDEvent*& esd)
{
  // run the trackleter for multiplicity study

  AliCodeTimerAuto("",0)

  AliTrackleter *trackleter = CreateMultFinder();
  if (!trackleter) return kFALSE;

  AliInfo(Form("running the ITS multiplicity finder: %s",trackleter->ClassName()));

  if (fLoader[0]) {
    fLoader[0]->LoadRecPoints();
    TTree* cltree = fLoader[0]->TreeR();
    if (cltree) {
      trackleter->Reconstruct(esd,cltree);
      AliMultiplicity *mult = trackleter->GetMultiplicity();
      if(mult) esd->SetMultiplicity(mult);
    }
    else {
      AliError("Can't get the ITS cluster tree");
    }
    fLoader[0]->UnloadRecPoints();
  }
  else {
    AliError("Can't get the ITS loader");
  }

  delete trackleter;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::RunMuonTracking(AliESDEvent*& esd)
{
// run the muon spectrometer tracking

  AliCodeTimerAuto("",0)

  if (!fRunLoader) {
    AliError("Missing runLoader!");
    return kFALSE;
  }
  Int_t iDet =  GetDetIndex("MUON"); // for MUON

  // Get a pointer to the MUON reconstructor
  AliReconstructor *reconstructor = GetReconstructor(iDet);
  if (!reconstructor) return kFALSE;

  
  TString detName = fgkDetectorName[iDet];
  AliDebug(1, Form("%s tracking", detName.Data()));
  AliTracker *tracker =  reconstructor->CreateTracker();
  if (!tracker) {
    AliWarning(Form("couldn't create a tracker for %s", detName.Data()));
    return kFALSE;
  }
     
  // read RecPoints
  fLoader[iDet]->LoadRecPoints("read");  

  tracker->LoadClusters(fLoader[iDet]->TreeR());
  
  Int_t rv = tracker->Clusters2Tracks(esd);
  
  fLoader[iDet]->UnloadRecPoints();

  tracker->UnloadClusters();
  
  if ( rv )
  {
    AliError(Form("%s Clusters2Tracks failed", fgkDetectorName[iDet]));
    return kFALSE;
  }
  
  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliReconstruction::RunMFTTrackingMU(AliESDEvent*& esd) {

  // AU

  // run the global muon tracking: matching the MUON tracks with the MFT clusters

  AliCodeTimerAuto("",0)

  if (!fRunLoader) {
    AliError("Missing runLoader!");
    return kFALSE;
  }
  Int_t iDet = GetDetIndex("MFT"); // for MFT

  // Get a pointer to the MFT reconstructor
  AliReconstructor *reconstructor = GetReconstructor(iDet);
  if (!reconstructor) return kFALSE;
  
  TString detName = fgkDetectorName[iDet];
  AliDebug(1, Form("%s tracking for muon tracks", detName.Data()));
  AliTracker *tracker = reconstructor->CreateTracker();
  if (!tracker) {
    AliWarning(Form("couldn't create a Muon tracker for %s", detName.Data()));
    return kFALSE;
  }
     
  // read RecPoints
  fLoader[iDet]->LoadRecPoints("read");  

  tracker->LoadClusters(fLoader[iDet]->TreeR());
  
  Int_t rv = tracker->Clusters2Tracks(esd);
  
  fLoader[iDet]->UnloadRecPoints();

  tracker->UnloadClusters();
  
  if (rv) {
    AliError(Form("%s Clusters2Tracks failed", fgkDetectorName[iDet]));
    return kFALSE;
  }
  
  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliReconstruction::RunTracking(AliESDEvent*& esd,AliESDpid &PID)
{
// run the barrel tracking
  static Int_t eventNr=0;
  AliCodeTimerAuto("",0)

  AliInfo("running tracking");

  // Set the event info which is used
  // by the trackers in order to obtain
  // information about read-out detectors,
  // trigger etc.
  AliDebug(1, "Setting event info");
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if (!fTracker[iDet]) continue;
    fTracker[iDet]->SetEventInfo(&fEventInfo);
  }

  //Fill the ESD with the T0 info (will be used by the TOF) 
  if (fReconstructor[11] && fLoader[11]) {
    fLoader[11]->LoadRecPoints("READ");
    TTree *treeR = fLoader[11]->TreeR();
    if (treeR) {
      GetReconstructor(11)->FillESD((TTree *)NULL,treeR,esd);
    }
  }

  // pass 1: TPC + ITS inwards
  for (Int_t iDet = 1; iDet >= 0; iDet--) {
    if (!fTracker[iDet]) continue;
    AliDebug(1, Form("%s tracking", fgkDetectorName[iDet]));

    // load clusters
    fLoader[iDet]->LoadRecPoints("read");
    AliSysInfo::AddStamp(Form("RLoadCluster%s_%d",fgkDetectorName[iDet],eventNr),iDet,1, eventNr);
    TTree* tree = fLoader[iDet]->TreeR();
    if (!tree) {
      AliError(Form("Can't get the %s cluster tree", fgkDetectorName[iDet]));
      return kFALSE;
    }
    fTracker[iDet]->LoadClusters(tree);
    AliSysInfo::AddStamp(Form("TLoadCluster%s_%d",fgkDetectorName[iDet],eventNr), iDet,2, eventNr);
    // run tracking
    if (fTracker[iDet]->Clusters2TracksHLT(esd, fhltesd) != 0) {
      AliError(Form("%s Clusters2Tracks failed", fgkDetectorName[iDet]));
      return kFALSE;
    }
    AliSysInfo::AddStamp(Form("Tracking0%s_%d",fgkDetectorName[iDet],eventNr), iDet,3,eventNr);
    // preliminary PID in TPC needed by the ITS tracker
    if (iDet == 1) {
      esd->SetNumberOfTPCClusters(fTracker[iDet]->GetNumberOfClusters());
      GetReconstructor(1)->FillESD((TTree*)NULL, (TTree*)NULL, esd);
      PID.MakePIDForTracking(esd);
      AliSysInfo::AddStamp(Form("MakePID0%s_%d",fgkDetectorName[iDet],eventNr), iDet,4,eventNr);
    } 
  }

  // pass 2: ALL backwards

  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if (!fTracker[iDet]) continue;
    AliDebug(1, Form("%s back propagation", fgkDetectorName[iDet]));

    // load clusters
    if (iDet > 1) {     // all except ITS, TPC
      TTree* tree = NULL;
      fLoader[iDet]->LoadRecPoints("read");
      AliSysInfo::AddStamp(Form("RLoadCluster0%s_%d",fgkDetectorName[iDet],eventNr), iDet,1, eventNr);
      tree = fLoader[iDet]->TreeR();
      if (!tree) {
        AliError(Form("Can't get the %s cluster tree", fgkDetectorName[iDet]));
        return kFALSE;
      }
      fTracker[iDet]->LoadClusters(tree); 
      AliSysInfo::AddStamp(Form("TLoadCluster0%s_%d",fgkDetectorName[iDet],eventNr), iDet,2, eventNr);
    }

    // run tracking
    if (iDet>1) // start filling residuals for the "outer" detectors
      if (fRunGlobalQA) {
        AliTracker::SetFillResiduals(fRecoParam.GetEventSpecie(), kTRUE);     
        TObjArray ** arr = AliTracker::GetResidualsArray() ; 
	if (arr) {
	  AliRecoParam::EventSpecie_t es=fRecoParam.GetEventSpecie();
	  TObjArray * elem = arr[AliRecoParam::AConvert(es)];
	  if ( elem && (! elem->At(0)) ) {
	    AliQADataMaker *qadm = AliQAManager::QAManager()->GetQADataMaker(AliQAv1::kGLOBAL);
	    if (qadm) qadm->InitRecPointsForTracker() ; 
	  }
	}
	//	AliSysInfo::AddStamp(Form("QAInitResid%s_%d",fgkDetectorName[iDet],eventNr), iDet,0, eventNr);
      }
    if (fTracker[iDet]->PropagateBack(esd) != 0) {
      AliError(Form("%s backward propagation failed", fgkDetectorName[iDet]));
      //      return kFALSE;
    }
    AliSysInfo::AddStamp(Form("Tracking1%s_%d",fgkDetectorName[iDet],eventNr), iDet,3, eventNr);

    // unload clusters
    if (iDet > 3) {     // all except ITS, TPC, TRD and TOF
      fTracker[iDet]->UnloadClusters();
      fLoader[iDet]->UnloadRecPoints();
    }
    // updated PID in TPC needed by the ITS tracker -MI
    if (iDet == 1) {
      //GetReconstructor(1)->FillESD((TTree*)NULL, (TTree*)NULL, esd);
      //AliESDpid::MakePID(esd);
      PID.MakePIDForTracking(esd);
      AliSysInfo::AddStamp(Form("MakePID1%s_%d",fgkDetectorName[iDet],eventNr), iDet,4,eventNr);
    }

  }
  //stop filling residuals for the "outer" detectors
  if (fRunGlobalQA) AliTracker::SetFillResiduals(fRecoParam.GetEventSpecie(), kFALSE);     

  // At this moment we decide if the friends need to be stored, since TPC stores heavy calib objects there
  fWriteThisFriend = fWriteESDfriend ? DecideFriendsStorage() : kFALSE;
  if (!fWriteThisFriend) fesd->SetNTPCFriend2Store(0); // don't store TPC seeds

  // pass 3: TRD + TPC + ITS refit inwards

  for (Int_t iDet = 2; iDet >= 0; iDet--) {
    if (!fTracker[iDet]) continue;
    AliDebug(1, Form("%s inward refit", fgkDetectorName[iDet]));

    // run tracking
    if (iDet<2) // start filling residuals for TPC and ITS
      if (fRunGlobalQA) {
        AliTracker::SetFillResiduals(fRecoParam.GetEventSpecie(), kTRUE);     
        TObjArray ** arr = AliTracker::GetResidualsArray() ; 
	if (arr) {
	  AliRecoParam::EventSpecie_t es=fRecoParam.GetEventSpecie();
	  TObjArray * elem = arr[AliRecoParam::AConvert(es)];
	  if ( elem && (! elem->At(0)) ) {
	    AliQADataMaker *qadm = AliQAManager::QAManager()->GetQADataMaker(AliQAv1::kGLOBAL);
	    if (qadm) qadm->InitRecPointsForTracker() ; 
	  }
	}
      }
    
    if (fTracker[iDet]->RefitInward(esd) != 0) {
      AliError(Form("%s inward refit failed", fgkDetectorName[iDet]));
      //      return kFALSE;
    }
    // run postprocessing
    if (fTracker[iDet]->PostProcess(esd) != 0) {
      AliError(Form("%s postprocessing failed", fgkDetectorName[iDet]));
      //      return kFALSE;
    }
    AliSysInfo::AddStamp(Form("Tracking2%s_%d",fgkDetectorName[iDet],eventNr), iDet,3, eventNr);
  }

  // write space-points to the ESD in case alignment data output
  // is switched on
  if (fWriteAlignmentData && fWriteESDfriend) {
    WriteAlignmentData(esd);
    AliSysInfo::AddStamp(Form("WrtAlignData_%d",eventNr), 0,0, eventNr);
  }
  
  for (Int_t iDet = 3; iDet >= 0; iDet--) {
    if (!fTracker[iDet]) continue;
    // 
    // RS: TPC clusters will be unloaded in the end of process event, since used in the friends calib objects
    // unload clusters
    if (iDet!=AliQAv1::kTPC) {
      fTracker[iDet]->UnloadClusters();
      AliSysInfo::AddStamp(Form("TUnloadCluster%s_%d",fgkDetectorName[iDet],eventNr), iDet,4, eventNr);
      fLoader[iDet]->UnloadRecPoints();
      AliSysInfo::AddStamp(Form("RUnloadCluster%s_%d",fgkDetectorName[iDet],eventNr), iDet,5, eventNr);
    }
  }
  // stop filling residuals for TPC and ITS
  if (fRunGlobalQA) AliTracker::SetFillResiduals(fRecoParam.GetEventSpecie(), kFALSE);     

  eventNr++;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::CleanESD(AliESDEvent *esd){
  //
  // Remove the data which are not needed for the physics analysis.
  //

  Int_t nTracks=esd->GetNumberOfTracks();
  Int_t nV0s=esd->GetNumberOfV0s();
  AliInfo
  (Form("Number of ESD tracks and V0s before cleaning: %d %d",nTracks,nV0s));
  TObjArray tracks2remove;
  Float_t cleanPars[]={fV0DCAmax,fV0CsPmin,fDmax,fZmax};
  Bool_t rc=esd->Clean(cleanPars,&tracks2remove);
  if (rc) { // physically destroy removed tracks
    for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
      if (!fTracker[iDet]) continue;
      if (fTracker[iDet]->OwnsESDObjects()) fTracker[iDet]->CleanESDTracksObjects(&tracks2remove);
    }
    tracks2remove.Delete();
  }
  nTracks=esd->GetNumberOfTracks();
  nV0s=esd->GetNumberOfV0s();
  AliInfo
  (Form("Number of ESD tracks and V0s after cleaning %d %d",nTracks,nV0s));

  return rc;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::FillESD(AliESDEvent*& esd, const TString& detectors)
{
// fill the event summary data

  AliCodeTimerAuto("",0)
    static Int_t eventNr=0; 
  TString detStr = detectors;
  
  AliSysInfo::AddStamp(Form("FillESDb%d",eventNr), -19,-19, eventNr);
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
  if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
    AliReconstructor* reconstructor = GetReconstructor(iDet);
    if (!reconstructor) continue;
    AliDebug(1, Form("filling ESD for %s", fgkDetectorName[iDet]));
    TTree* clustersTree = NULL;
    if (fLoader[iDet]) {
      fLoader[iDet]->LoadRecPoints("read");
      clustersTree = fLoader[iDet]->TreeR();
      if (!clustersTree) {
	AliError(Form("Can't get the %s clusters tree", 
		      fgkDetectorName[iDet]));
	if (fStopOnError) return kFALSE;
      }
    }
    if (fRawReader && !reconstructor->HasDigitConversion()) {
      reconstructor->FillESD(fRawReader, clustersTree, esd);
    } else {
      TTree* digitsTree = NULL;
      if (fLoader[iDet]) {
	fLoader[iDet]->LoadDigits("read");
	digitsTree = fLoader[iDet]->TreeD();
	if (!digitsTree) {
	  AliError(Form("Can't get the %s digits tree", 
			fgkDetectorName[iDet]));
	  if (fStopOnError) return kFALSE;
	}
      }
      reconstructor->FillESD(digitsTree, clustersTree, esd);
      if (fLoader[iDet]) fLoader[iDet]->UnloadDigits();
    }
    if (fLoader[iDet]) {
      fLoader[iDet]->UnloadRecPoints();
    }
  }
  
  if (!IsSelected("CTP", detStr)) AliDebug(10,"No CTP");
  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    AliError(Form("the following detectors were not found: %s", 
                  detStr.Data()));
    if (fStopOnError) return kFALSE;
  }
  AliSysInfo::AddStamp(Form("FillESDe%d",eventNr), -20,-20, eventNr);
  eventNr++;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::FillTriggerESD(AliESDEvent*& esd)
{
  // Reads the trigger decision which is
  // stored in Trigger.root file and fills
  // the corresponding esd entries

  AliCodeTimerAuto("",0)
  
  AliInfo("Filling trigger information into the ESD");

  if (fRawReader) {
    AliCTPRawStream input(fRawReader);
    if (!input.Next()) {
      AliWarning("No valid CTP (trigger) DDL raw data is found ! The trigger info is taken from the event header!");
    }
    else {
      if (esd->GetTriggerMask() != input.GetClassMask())
	AliError(Form("Invalid trigger pattern found in CTP raw-data: %llx %llx",
		      input.GetClassMask(),esd->GetTriggerMask()));
      if (esd->GetTriggerMaskNext50() != input.GetClassMaskNext50())
	AliError(Form("Invalid trigger pattern found in CTP raw-data Next50: %llx %llx",
		      input.GetClassMaskNext50(),esd->GetTriggerMaskNext50()));
      if (esd->GetOrbitNumber() != input.GetOrbitID())
	AliError(Form("Invalid orbit id found in CTP raw-data: %x %x",
		      input.GetOrbitID(),esd->GetOrbitNumber()));
      if (esd->GetBunchCrossNumber() != input.GetBCID())
	AliError(Form("Invalid bunch-crossing id found in CTP raw-data: %x %x",
		      input.GetBCID(),esd->GetBunchCrossNumber()));
      AliESDHeader* esdheader = esd->GetHeader();
      esdheader->SetL0TriggerInputs(input.GetL0Inputs());
      esdheader->SetL1TriggerInputs(input.GetL1Inputs());
      esdheader->SetL2TriggerInputs(input.GetL2Inputs());
      // IR
      //      UInt_t orbit=input.GetOrbitID();
      for(Int_t i=0 ; i<input.GetNIRs() ; i++ ) {
	esdheader->AddTriggerIR(input.GetIR(i));
      }
       AliCentralTrigger* rlCTP = fRunLoader->GetTrigger();
       if (rlCTP) {
	 rlCTP->SetL0TriggerInputs(input.GetL0Inputs());
	 rlCTP->SetL1TriggerInputs(input.GetL1Inputs());
	 rlCTP->SetL2TriggerInputs(input.GetL2Inputs());
       }
    }
    if (fIsNewRunLoader) fRunLoader->TreeCT()->Fill();
  }
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliReconstruction::FillTriggerScalers(AliESDEvent*& esd)
{
  //Scalers
  //fRunScalers->Print();
  if(fRunScalers && fRunScalers->CheckRunScalers()){
     AliTimeStamp* timestamp = new AliTimeStamp(esd->GetOrbitNumber(), esd->GetPeriodNumber(), esd->GetBunchCrossNumber());
     //AliTimeStamp* timestamp = new AliTimeStamp(10308000, 0, (ULong64_t)486238);
     AliESDHeader* esdheader = fesd->GetHeader();
     for(Int_t i=0;i<50;i++){
          if((1ull<<i) & esd->GetTriggerMask()){
          AliTriggerScalersESD* scalesd = fRunScalers->GetScalersForEventClass( timestamp, i+1);
          if(scalesd)esdheader->SetTriggerScalersRecord(scalesd);
        }
     }
     for(Int_t i=0;i<50;i++){
          if((1ull<<i) & esd->GetTriggerMaskNext50()){
          AliTriggerScalersESD* scalesd = fRunScalers->GetScalersForEventClass( timestamp, i+51);
          if(scalesd)esdheader->SetTriggerScalersRecord(scalesd);
        }
     }
     const AliTriggerScalersRecordESD* scalrecEvent = fRunScalers->GetScalersDeltaForEvent( timestamp);
     const AliTriggerScalersRecordESD* scalrecRun = fRunScalers->GetScalersDeltaForRun();
     if (scalrecEvent) esdheader->SetTriggerScalersDeltaEvent(scalrecEvent);
     if (scalrecRun) esdheader->SetTriggerScalersDeltaRun(scalrecRun);
  }
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::FillRawEventHeaderESD(AliESDEvent*& esd)
{
  // 
  // Filling information from RawReader Header
  // 

  if (!fRawReader) return kFALSE;

  AliInfo("Filling information from RawReader Header");

  esd->SetBunchCrossNumber(fRawReader->GetBCID());
  esd->SetOrbitNumber(fRawReader->GetOrbitID());
  esd->SetPeriodNumber(fRawReader->GetPeriod());

  esd->SetTimeStamp(fRawReader->GetTimestamp());  
  esd->SetEventType(fRawReader->GetType());

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::FillMCEventHeaderESD(AliESDEvent*& esd)
{
  // 
  // Filling information from MC Header
  // 
  AliHeader* headerMC = fRunLoader->GetHeader();
  if (!headerMC || !esd) return kFALSE;

  AliInfo("Filling information from MC Header");
  esd->SetTimeStamp(headerMC->GetTimeStamp());  

  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliReconstruction::IsSelected(TString detName, TString& detectors) const
{
// check whether detName is contained in detectors
// if yes, it is removed from detectors

  // check if all detectors are selected
  if ((detectors.CompareTo("ALL") == 0) ||
      detectors.BeginsWith("ALL ") ||
      detectors.EndsWith(" ALL") ||
      detectors.Contains(" ALL ")) {
    detectors = "ALL";
    return kTRUE;
  }

  // search for the given detector
  Bool_t result = kFALSE;
  if ((detectors.CompareTo(detName) == 0) ||
      detectors.BeginsWith(detName+" ") ||
      detectors.EndsWith(" "+detName) ||
      detectors.Contains(" "+detName+" ")) {
    detectors.ReplaceAll(detName, "");
    result = kTRUE;
  }

  // clean up the detectors string
  while (detectors.Contains("  ")) detectors.ReplaceAll("  ", " ");
  while (detectors.BeginsWith(" ")) detectors.Remove(0, 1);
  while (detectors.EndsWith(" ")) detectors.Remove(detectors.Length()-1, 1);

  return result;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::InitRunLoader()
{
// get or create the run loader

  if (gAlice) delete gAlice;
  gAlice = NULL;

  TFile *gafile = TFile::Open(fGAliceFileName.Data());
  //  if (!gSystem->AccessPathName(fGAliceFileName.Data())) { // galice.root exists
  if (gafile) { // galice.root exists
    gafile->Close();
    delete gafile;

    // load all base libraries to get the loader classes
    TString libs = gSystem->GetLibraries();
    for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
      TString detName = fgkDetectorName[iDet];
      if (libs.Contains("lib" + detName + "base.so")) continue;
      gSystem->Load("lib" + detName + "base");
    }
    fRunLoader = AliRunLoader::Open(fGAliceFileName.Data());
    if (!fRunLoader) {
      AliError(Form("no run loader found in file %s", fGAliceFileName.Data()));
      CleanUp();
      return kFALSE;
    }

    fRunLoader->CdGAFile();
    fRunLoader->LoadgAlice();

    //PH This is a temporary fix to give access to the kinematics
    //PH that is needed for the labels of ITS clusters
    fRunLoader->LoadHeader();
    fRunLoader->LoadKinematics();

  } else {               // galice.root does not exist
    if (!fRawReader) {
      AliError(Form("the file %s does not exist", fGAliceFileName.Data()));
    }
    fRunLoader = AliRunLoader::Open(fGAliceFileName.Data(),
				    AliConfig::GetDefaultEventFolderName(),
				    "recreate");
    if (!fRunLoader) {
      AliError(Form("could not create run loader in file %s", 
		    fGAliceFileName.Data()));
      CleanUp();
      return kFALSE;
    }
    fIsNewRunLoader = kTRUE;
    fRunLoader->MakeTree("E");
    fRunLoader->MakeTree("GG");

    if (fNumberOfEventsPerFile > 0)
      fRunLoader->SetNumberOfEventsPerFile(fNumberOfEventsPerFile);
    else
      fRunLoader->SetNumberOfEventsPerFile((UInt_t)-1);
  }

  return kTRUE;
}

//_____________________________________________________________________________
AliReconstructor* AliReconstruction::GetReconstructor(Int_t iDet)
{
// get the reconstructor object and the loader for a detector

  if (fReconstructor[iDet]) {
    if (fRecoParam.GetDetRecoParamArray(iDet) && !AliReconstructor::GetRecoParam(iDet)) {
      const AliDetectorRecoParam *par = fRecoParam.GetDetRecoParam(iDet);
      fReconstructor[iDet]->SetRecoParam(par);
      fReconstructor[iDet]->SetRunInfo(fRunInfo);
    }
    return fReconstructor[iDet];
  }

  // load the reconstructor object
  TPluginManager* pluginManager = gROOT->GetPluginManager();
  TString detName = fgkDetectorName[iDet];
  TString recName = "Ali" + detName + "Reconstructor";

  if (!fIsNewRunLoader && !fRunLoader->GetLoader(detName+"Loader") && (detName != "HLT")) return NULL;

  AliReconstructor* reconstructor = NULL;
  // first check if a plugin is defined for the reconstructor
  TPluginHandler* pluginHandler = 
    pluginManager->FindHandler("AliReconstructor", detName);
  // if not, add a plugin for it
  if (!pluginHandler) {
    AliDebug(1, Form("defining plugin for %s", recName.Data()));
    TString libs = gSystem->GetLibraries();
    if (libs.Contains("lib" + detName + "base.so") ||
	(gSystem->Load("lib" + detName + "base") >= 0)) {
      pluginManager->AddHandler("AliReconstructor", detName, 
				recName, detName + "rec", recName + "()");
    } else {
      pluginManager->AddHandler("AliReconstructor", detName, 
				recName, detName, recName + "()");
    }
    pluginHandler = pluginManager->FindHandler("AliReconstructor", detName);
  }
  if (pluginHandler && (pluginHandler->LoadPlugin() == 0)) {
    reconstructor = (AliReconstructor*) pluginHandler->ExecPlugin(0);
  }

   // check if the upgrade reconstructor should be used instead of the standard one
  if(fUpgradeMask[iDet]) {
    if(reconstructor) delete reconstructor;
    TClass *cl = new TClass(Form("Ali%sUpgradeReconstructor",fgkDetectorName[iDet]));
    reconstructor = (AliReconstructor*)(cl->New());
   }

  if (reconstructor) {
    TObject* obj = fOptions.FindObject(detName.Data());
    if (obj) reconstructor->SetOption(obj->GetTitle());
    reconstructor->SetRunInfo(fRunInfo);
    reconstructor->SetHLTESD(fhltesd);
    reconstructor->Init();
    fReconstructor[iDet] = reconstructor;
  }

  // get or create the loader
  if (detName != "HLT") {
    fLoader[iDet] = fRunLoader->GetLoader(detName + "Loader");
    if (!fLoader[iDet]) {
      AliConfig::Instance()
	->CreateDetectorFolders(fRunLoader->GetEventFolder(), 
				detName, detName);
      // first check if a plugin is defined for the loader
      pluginHandler = 
	pluginManager->FindHandler("AliLoader", detName);
      // if not, add a plugin for it
      if (!pluginHandler) {
	TString loaderName = "Ali" + detName + "Loader";
	AliDebug(1, Form("defining plugin for %s", loaderName.Data()));
	pluginManager->AddHandler("AliLoader", detName, 
				  loaderName, detName + "base", 
				  loaderName + "(const char*, TFolder*)");
	pluginHandler = pluginManager->FindHandler("AliLoader", detName);
      }
      if (pluginHandler && (pluginHandler->LoadPlugin() == 0)) {
	fLoader[iDet] = 
	  (AliLoader*) pluginHandler->ExecPlugin(2, detName.Data(), 
						 fRunLoader->GetEventFolder());
      }
      if (!fLoader[iDet]) {   // use default loader
	fLoader[iDet] = new AliLoader(detName, fRunLoader->GetEventFolder());
      }
      if (!fLoader[iDet]) {
	AliWarning(Form("couldn't get loader for %s", detName.Data()));
	if (fStopOnError) return NULL;
      } else {
	fRunLoader->AddLoader(fLoader[iDet]);
	fRunLoader->CdGAFile();
	if (gFile && !gFile->IsWritable()) gFile->ReOpen("UPDATE");
	fRunLoader->Write(0, TObject::kOverwrite);
      }
    }
  }
      
  if (fRecoParam.GetDetRecoParamArray(iDet) && !AliReconstructor::GetRecoParam(iDet)) {
    const AliDetectorRecoParam *par = fRecoParam.GetDetRecoParam(iDet);
    if (reconstructor) {
      reconstructor->SetRecoParam(par);
      reconstructor->SetRunInfo(fRunInfo);
    }
  }
  return reconstructor;
}

//_____________________________________________________________________________
AliVertexer* AliReconstruction::CreateVertexer()
{
// create the vertexer
// Please note that the caller is the owner of the
// vertexer

  AliVertexer* vertexer = NULL;
  AliReconstructor* itsReconstructor = GetReconstructor(0);
  if (itsReconstructor && ((fRunLocalReconstruction.Contains("ITS")) || 
			   fRunTracking.Contains("ITS") || fFillESD.Contains("ITS") )) {
    vertexer = itsReconstructor->CreateVertexer();
  }
  if (!vertexer) {
    AliWarning("couldn't create a vertexer for ITS");
  }

  return vertexer;
}

//_____________________________________________________________________________
AliTrackleter* AliReconstruction::CreateMultFinder()
{
// create the ITS trackleter for mult. estimation
// Please note that the caller is the owner of the
// trackleter

  AliTrackleter* trackleter = NULL;
  AliReconstructor* itsReconstructor = GetReconstructor(0);
  if (itsReconstructor && ((fRunLocalReconstruction.Contains("ITS")) || 
			   fRunTracking.Contains("ITS") || fFillESD.Contains("ITS") )) {
    trackleter = itsReconstructor->CreateMultFinder();
  }
  else {
    AliWarning("ITS is not in reconstruction, switching off RunMultFinder");
    fRunMultFinder = kFALSE;
  }

  return trackleter;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::CreateTrackers(const TString& detectors)
{
// create the trackers
	AliInfo("Creating trackers");

  TString detStr = detectors;
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
    AliReconstructor* reconstructor = GetReconstructor(iDet);
    if (!reconstructor) continue;
    TString detName = fgkDetectorName[iDet];
    if (detName == "MUON") {
      fRunMuonTracking = kTRUE;
      continue;
    }
    if (detName == "MFT") {           // AU    
      fRunMFTTrackingMU = kTRUE;      // AU
      continue;			      // AU
    }                                 // AU

    fTracker[iDet] = reconstructor->CreateTracker();
    if (!fTracker[iDet] && (iDet < 7)) {
      AliWarning(Form("couldn't create a tracker for %s", detName.Data()));
      if (fStopOnError) return kFALSE;
    }
    AliSysInfo::AddStamp(Form("LTracker%s",fgkDetectorName[iDet]), iDet,0);
  }

  return kTRUE;
}

//_____________________________________________________________________________
void AliReconstruction::CleanUp()
{
// delete trackers and the run loader and close and delete the file
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if (fReconstructor[iDet]) fReconstructor[iDet]->SetRecoParam(NULL);
    delete fReconstructor[iDet];
    fReconstructor[iDet] = NULL;
    fLoader[iDet] = NULL;
    delete fTracker[iDet];
    fTracker[iDet] = NULL;
  }

  delete fRunInfo;
  fRunInfo = NULL;

  delete fSPDTrackleter;
  fSPDTrackleter = NULL;

  delete ftVertexer;
  ftVertexer = NULL;
  
  delete fRunLoader;
  fRunLoader = NULL;
  delete fRawReader;
  fRawReader = NULL;
  delete fParentRawReader;
  fParentRawReader=NULL;

  delete fESDpid;
  fESDpid = NULL;

  if (ffile) {
    ffile->Close();
    delete ffile;
    ffile = NULL;
  }

  if (AliQAManager::QAManager())
    AliQAManager::QAManager()->ShowQA() ; 
  //  AliQAManager::Destroy() ;
  delete fAnalysis; 
  fAnalysis = NULL;
}

void AliReconstruction::WriteAlignmentData(AliESDEvent* esd)
{
  // Write space-points which are then used in the alignment procedures
  // For the moment only ITS, TPC, TRD and TOF
  if (!fWriteThisFriend) {
    AliInfo("Will not write alignment data since the friends are not stored");
    return;
  }
  Int_t ntracks = esd->GetNumberOfTracks();
  for (Int_t itrack = 0; itrack < ntracks; itrack++)
    {
      AliESDtrack *track = esd->GetTrack(itrack);
      if (track->GetFriendNotStored()) continue; // no need of alignment data for friends which are not stored
      Int_t nsp = 0;
      Int_t idx[200];
      for (Int_t i=0; i<200; ++i) idx[i] = -1; //PH avoid uninitialized values
      for (Int_t iDet = 5; iDet >= 0; iDet--) {// TOF, TRD, TPC, ITS clusters
          nsp += (iDet==GetDetIndex("TRD")) ? track->GetTRDntracklets():track->GetNcls(iDet);

          if (iDet==GetDetIndex("ITS")) { // ITS "extra" clusters
             track->GetClusters(iDet,idx);
             for (Int_t i=6; i<12; i++) if(idx[i] >= 0) nsp++;
          }  
      }

      if (nsp) {
	AliTrackPointArray *sp = new AliTrackPointArray(nsp);
	track->SetTrackPointArray(sp);
	sp->SetBit(AliTrackPointArray::kTOFBugFixed);
	Int_t isptrack = 0;
	for (Int_t iDet = 5; iDet >= 0; iDet--) {
	  AliTracker *tracker = fTracker[iDet];
	  if (!tracker) continue;
	  Int_t nspdet = (iDet==GetDetIndex("TRD")) ? track->GetTRDtracklets(idx):track->GetClusters(iDet,idx);

	  if (iDet==GetDetIndex("ITS")) // ITS "extra" clusters             
             for (Int_t i=6; i<12; i++) if(idx[i] >= 0) nspdet++;

	  if (nspdet <= 0) continue;
	  AliTrackPoint p;
	  Int_t isp = 0;
	  Int_t isp2 = 0;
	  while (isp2 < nspdet) {
	    Bool_t isvalid=kTRUE;

            Int_t index=idx[isp++];
            if (index < 0) continue;

            TString dets = fgkDetectorName[iDet];
            if ((fUseTrackingErrorsForAlignment.CompareTo(dets) == 0) ||
            fUseTrackingErrorsForAlignment.BeginsWith(dets+" ") ||
            fUseTrackingErrorsForAlignment.EndsWith(" "+dets) ||
            fUseTrackingErrorsForAlignment.Contains(" "+dets+" ")) {
              isvalid = tracker->GetTrackPointTrackingError(index,p,track);
	    } else {
	      isvalid = tracker->GetTrackPoint(index,p); 
	    } 
	    isp2++;
	    if (!isvalid) continue;
	    if (iDet==GetDetIndex("ITS") && (isp-1)>=6) p.SetExtra();
	    sp->AddPoint(isptrack,&p); isptrack++;
	  }
	}	
      }
    }
}

//_____________________________________________________________________________
void AliReconstruction::FillRawDataErrorLog(Int_t iEvent, AliESDEvent* esd)
{
  // The method reads the raw-data error log
  // accumulated within the rawReader.
  // It extracts the raw-data errors related to
  // the current event and stores them into
  // a TClonesArray inside the esd object.

  if (!fRawReader) return;

  for(Int_t i = 0; i < fRawReader->GetNumberOfErrorLogs(); i++) {

    AliRawDataErrorLog *log = fRawReader->GetErrorLog(i);
    if (!log) continue;
    if (iEvent != log->GetEventNumber()) continue;

    esd->AddRawDataErrorLog(log);
  }

}

//_____________________________________________________________________________
// void AliReconstruction::CheckQA()
// {
// check the QA of SIM for this run and remove the detectors 
// with status Fatal
  
//	TString newRunLocalReconstruction ; 
//	TString newRunTracking ;
//	TString newFillESD ;
//	 
//	for (Int_t iDet = 0; iDet < AliQAv1::kNDET; iDet++) {
//		TString detName(AliQAv1::GetDetName(iDet)) ;
//		AliQAv1 * qa = AliQAv1::Instance(AliQAv1::DETECTORINDEX_t(iDet)) ;       
//      if ( qa->IsSet(AliQAv1::DETECTORINDEX_t(iDet), AliQAv1::kSIM, specie, AliQAv1::kFATAL)) {
//        AliInfo(Form("QA status for %s %s in Hits and/or SDIGITS  and/or Digits was Fatal; No reconstruction performed", 
//                   detName.Data(), AliRecoParam::GetEventSpecieName(es))) ;
//			} else {
//			if ( fRunLocalReconstruction.Contains(AliQAv1::GetDetName(iDet)) || 
//					fRunLocalReconstruction.Contains("ALL") )  {
//				newRunLocalReconstruction += detName ; 
//				newRunLocalReconstruction += " " ; 			
//			}
//			if ( fRunTracking.Contains(AliQAv1::GetDetName(iDet)) || 
//					fRunTracking.Contains("ALL") )  {
//				newRunTracking += detName ; 
//				newRunTracking += " " ; 			
//			}
//			if ( fFillESD.Contains(AliQAv1::GetDetName(iDet)) || 
//					fFillESD.Contains("ALL") )  {
//				newFillESD += detName ; 
//				newFillESD += " " ; 			
//			}
//		}
//	}
//	fRunLocalReconstruction = newRunLocalReconstruction ; 
//	fRunTracking            = newRunTracking ; 
//	fFillESD                = newFillESD ; 
// }

//_____________________________________________________________________________
Int_t AliReconstruction::GetDetIndex(const char* detector)
{
  // return the detector index corresponding to detector
  Int_t index = -1 ; 
  for (index = 0; index < kNDetectors ; index++) {
    if ( strcmp(detector, fgkDetectorName[index]) == 0 )
	break ; 
  }	
  return index ; 
}
//_____________________________________________________________________________
Bool_t AliReconstruction::FinishPlaneEff() {
 //
 // Here execute all the necessary operationis, at the end of the tracking phase,
 // in case that evaluation of PlaneEfficiencies was required for some detector.
 // E.g., write into a DataBase file the PlaneEfficiency which have been evaluated.
 //
 // This Preliminary version works only FOR ITS !!!!!
 // other detectors (TOF,TRD, etc. have to develop their specific codes)
 //
 //  Input: none
 //  Return: kTRUE if all operations have been done properly, kFALSE otherwise
 //
 Bool_t ret=kFALSE;
 TString detStr = fLoadCDB;
 //for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
 for (Int_t iDet = 0; iDet < 1; iDet++) { // for the time being only ITS
   if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
   if(fTracker[iDet] && fTracker[iDet]->GetPlaneEff()) {
      AliPlaneEff *planeeff=fTracker[iDet]->GetPlaneEff();
      TString name=planeeff->GetName();
      name+=".root";
      TFile* pefile = TFile::Open(name, "RECREATE");
      ret=(Bool_t)planeeff->Write();
      pefile->Close();
      if(planeeff->GetCreateHistos()) {
        TString hname=planeeff->GetName();
        hname+="Histo.root";
        ret*=planeeff->WriteHistosToFile(hname,"RECREATE");
      }
   }
   if(fSPDTrackleter) {
     AliPlaneEff *planeeff=fSPDTrackleter->GetPlaneEff();
      TString name="AliITSPlaneEffSPDtracklet.root";
      TFile* pefile = TFile::Open(name, "RECREATE");
      ret=(Bool_t)planeeff->Write();
      pefile->Close();
      AliESDEvent *dummy=NULL;
      ret=(Bool_t)fSPDTrackleter->PostProcess(dummy); // take care of writing other files
   }
 }
 return ret;
}
//_____________________________________________________________________________
Bool_t AliReconstruction::InitPlaneEff() {
//
 // Here execute all the necessary operations, before of the tracking phase,
 // for the evaluation of PlaneEfficiencies, in case required for some detectors.
 // E.g., read from a DataBase file a first evaluation of the PlaneEfficiency
 // which should be updated/recalculated.
 //
 // This Preliminary version will work only FOR ITS !!!!!
 // other detectors (TOF,TRD, etc. have to develop their specific codes)
 //
 //  Input: none
 //  Return: kTRUE if all operations have been done properly, kFALSE otherwise
 //

  fSPDTrackleter = NULL;
  TString detStr = fLoadCDB;
  if (IsSelected(fgkDetectorName[0], detStr)) {
    AliReconstructor* itsReconstructor = GetReconstructor(0);
    if (itsReconstructor) {
      fSPDTrackleter = itsReconstructor->CreateTrackleter(); // this is NULL unless required in RecoParam
    }
    if (fSPDTrackleter) {
      AliInfo("Trackleter for SPD has been created");
    }
  }
 return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::InitAliEVE()
{
  // This method should be called only in case 
  // AliReconstruction is run
  // within the alieve environment.
  // It will initialize AliEVE in a way
  // so that it can visualize event processed
  // by AliReconstruction.
  // The return flag shows whenever the
  // AliEVE initialization was successful or not.

  TString macroStr(getenv("ALIEVE_ONLINE_MACRO"));

  if (macroStr.IsNull())
    macroStr.Form("%s/EVE/macros/alieve_online.C",gSystem->ExpandPathName("$ALICE_ROOT"));

  AliInfo(Form("Loading AliEVE macro: %s",macroStr.Data()));

  if (gROOT->LoadMacro(macroStr.Data()) != 0) return kFALSE;

  gROOT->ProcessLine("if (!AliEveEventManager::GetMaster()){new AliEveEventManager();AliEveEventManager::GetMaster()->AddNewEventCommand(\"alieve_online_on_new_event()\");gEve->AddEvent(AliEveEventManager::GetMaster());};");
  gROOT->ProcessLine("alieve_online_init()");

  return kTRUE;
}
  
//_____________________________________________________________________________
void AliReconstruction::RunAliEVE()
{
  // Runs AliEVE visualisation of
  // the current event.
  // Should be executed only after
  // successful initialization of AliEVE.

  AliInfo("Running AliEVE...");
  gROOT->ProcessLine(Form("AliEveEventManager::GetMaster()->SetEvent((AliRunLoader*)%p,(AliRawReader*)%p,(AliESDEvent*)%p,(AliESDfriend*)%p);",fRunLoader,fRawReader,fesd,fesdf));
  gSystem->Run();
}

//_____________________________________________________________________________
Bool_t AliReconstruction::SetRunQA(TString detAndAction) 
{
	// Allows to run QA for a selected set of detectors
	// and a selected set of tasks among RAWS, DIGITSR, RECPOINTS and ESDS
	// all selected detectors run the same selected tasks
	
	if (!detAndAction.Contains(":")) {
		AliError( Form("%s is a wrong syntax, use \"DetectorList:ActionList\" \n", detAndAction.Data()) ) ;
		fRunQA = kFALSE ;
		return kFALSE ; 		
	}
	Int_t colon = detAndAction.Index(":") ; 
	fQADetectors = detAndAction(0, colon) ; 
	fQATasks   = detAndAction(colon+1, detAndAction.Length() ) ; 
	if (fQATasks.Contains("ALL") ) {
		fQATasks = Form("%d %d %d %d", AliQAv1::kRAWS, AliQAv1::kDIGITSR, AliQAv1::kRECPOINTS, AliQAv1::kESDS) ; 
	} else {
		fQATasks.ToUpper() ; 
		TString tempo("") ; 
		if ( fQATasks.Contains("RAW") ) 
			tempo = Form("%d ", AliQAv1::kRAWS) ; 
		if ( fQATasks.Contains("DIGIT") ) 
			tempo += Form("%d ", AliQAv1::kDIGITSR) ; 
		if ( fQATasks.Contains("RECPOINT") ) 
			tempo += Form("%d ", AliQAv1::kRECPOINTS) ; 
		if ( fQATasks.Contains("ESD") ) 
			tempo += Form("%d ", AliQAv1::kESDS) ; 
		fQATasks = tempo ; 
		if (fQATasks.IsNull()) {
			AliInfo("No QA requested\n")  ;
			fRunQA = kFALSE ;
			return kTRUE ; 
		}
	}	
	TString tempo(fQATasks) ; 
	tempo.ReplaceAll(Form("%d", AliQAv1::kRAWS), AliQAv1::GetTaskName(AliQAv1::kRAWS)) 	;
	tempo.ReplaceAll(Form("%d", AliQAv1::kDIGITSR), AliQAv1::GetTaskName(AliQAv1::kDIGITSR)) ;	
	tempo.ReplaceAll(Form("%d", AliQAv1::kRECPOINTS), AliQAv1::GetTaskName(AliQAv1::kRECPOINTS)) ;	
	tempo.ReplaceAll(Form("%d", AliQAv1::kESDS), AliQAv1::GetTaskName(AliQAv1::kESDS)) ; 	
	AliInfo( Form("QA will be done on \"%s\" for \"%s\"\n", fQADetectors.Data(), tempo.Data()) ) ;  
	fRunQA = kTRUE ;
	return kTRUE; 
} 

//_____________________________________________________________________________
Bool_t AliReconstruction::InitRecoParams() 
{
  // The method accesses OCDB and retrieves all
  // the available reco-param objects from there.

  Bool_t isOK = kTRUE;

  if (fRecoParam.GetDetRecoParamArray(kNDetectors)) {
    AliInfo("Using custom GRP reconstruction parameters");
  }
  else {
    AliInfo("Loading GRP reconstruction parameter objects");

    AliCDBPath path("GRP","Calib","RecoParam");
    AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());
    if(!entry){ 
      AliWarning("Couldn't find GRP RecoParam entry in OCDB");
      isOK = kFALSE;
    }
    else {
      TObject *recoParamObj = entry->GetObject();
      if (dynamic_cast<TObjArray*>(recoParamObj)) {
	// GRP has a normal TobjArray of AliDetectorRecoParam objects
	// Registering them in AliRecoParam
	fRecoParam.AddDetRecoParamArray(kNDetectors,dynamic_cast<TObjArray*>(recoParamObj));
      }
      else if (dynamic_cast<AliDetectorRecoParam*>(recoParamObj)) {
	// GRP has only onse set of reco parameters
	// Registering it in AliRecoParam
	AliInfo("Single set of GRP reconstruction parameters found");
	dynamic_cast<AliDetectorRecoParam*>(recoParamObj)->SetAsDefault();
	fRecoParam.AddDetRecoParam(kNDetectors,dynamic_cast<AliDetectorRecoParam*>(recoParamObj));
      }
      else {
	AliError("No valid GRP RecoParam object found in the OCDB");
	isOK = kFALSE;
      }
      entry->SetOwner(0);
    }
  }

  TString detStr = fLoadCDB;
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {

    if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;

    if (fRecoParam.GetDetRecoParamArray(iDet)) {
      AliInfo(Form("Using custom reconstruction parameters for detector %s",fgkDetectorName[iDet]));
      continue;
    }

    AliInfo(Form("Loading reconstruction parameter objects for detector %s",fgkDetectorName[iDet]));
  
    AliCDBPath path(fgkDetectorName[iDet],"Calib","RecoParam");
    AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());
    if(!entry){ 
      AliWarning(Form("Couldn't find RecoParam entry in OCDB for detector %s",fgkDetectorName[iDet]));
      isOK = kFALSE;
    }
    else {
      TObject *recoParamObj = entry->GetObject();
      if (dynamic_cast<TObjArray*>(recoParamObj)) {
	// The detector has a normal TobjArray of AliDetectorRecoParam objects
	// Registering them in AliRecoParam
	fRecoParam.AddDetRecoParamArray(iDet,dynamic_cast<TObjArray*>(recoParamObj));
      }
      else if (dynamic_cast<AliDetectorRecoParam*>(recoParamObj)) {
	// The detector has only onse set of reco parameters
	// Registering it in AliRecoParam
	AliInfo(Form("Single set of reconstruction parameters found for detector %s",fgkDetectorName[iDet]));
	dynamic_cast<AliDetectorRecoParam*>(recoParamObj)->SetAsDefault();
	fRecoParam.AddDetRecoParam(iDet,dynamic_cast<AliDetectorRecoParam*>(recoParamObj));
      }
      else {
	AliError(Form("No valid RecoParam object found in the OCDB for detector %s",fgkDetectorName[iDet]));
	isOK = kFALSE;
      }
      entry->SetOwner(0);
      //      FIX ME: We have to disable the unloading of reco-param CDB
      //      entries because QA framework is using them. Has to be fix in
      //      a way that the QA takes the objects already constructed in
      //      this method.
      //      AliCDBManager::Instance()->UnloadFromCache(path.GetPath());
    }
  }

  if (AliDebugLevel() > 0) fRecoParam.Print();

  return isOK;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::GetEventInfo() 
{
  // Fill the event info object
  // ...
  AliCodeTimerAuto("",0)
  
  AliCentralTrigger *aCTP = NULL;
  fEventInfo.Reset();
  if (fRawReader) {
    fEventInfo.SetEventType(fRawReader->GetType());

    ULong64_t mask = fRawReader->GetClassMask();
    ULong64_t maskNext50 = fRawReader->GetClassMaskNext50();
    fEventInfo.SetTriggerMask(mask);
    fEventInfo.SetTriggerMaskNext50(maskNext50);
    UInt_t clmask = fRawReader->GetDetectorPattern()[0];
    fEventInfo.SetTriggerCluster(AliDAQ::ListOfTriggeredDetectors(clmask));

    aCTP = new AliCentralTrigger();
    TString configstr("");
    if (!aCTP->LoadConfiguration(configstr)) { // Load CTP config from OCDB
      AliError("No trigger configuration found in OCDB! The trigger configuration information will not be used!");
      delete aCTP;
      return kFALSE;
    }
    aCTP->SetClassMask(mask);
    aCTP->SetClusterMask(clmask);

    if (fRunLoader) {
      AliCentralTrigger* rlCTP = fRunLoader->GetTrigger();
      if (rlCTP) {
	rlCTP->SetClassMask(mask);
	rlCTP->SetClusterMask(clmask);
      }
    }
  }
  else {
    fEventInfo.SetEventType(AliRawEventHeaderBase::kPhysicsEvent);
    if (fRunLoader && (!fRunLoader->LoadTrigger())) {
      aCTP = fRunLoader->GetTrigger();
      fEventInfo.SetTriggerMask(aCTP->GetClassMask());
      fEventInfo.SetTriggerMaskNext50(aCTP->GetClassMaskNext50());
      // get inputs from actp - just get
      AliESDHeader* esdheader = fesd->GetHeader();
      esdheader->SetL0TriggerInputs(aCTP->GetL0TriggerInputs());
      esdheader->SetL1TriggerInputs(aCTP->GetL1TriggerInputs());
      esdheader->SetL2TriggerInputs(aCTP->GetL2TriggerInputs());
      fEventInfo.SetTriggerCluster(AliDAQ::ListOfTriggeredDetectors(aCTP->GetClusterMask()));
    }
    else {
      if (fStopOnMissingTriggerFile) AliFatal("No trigger can be loaded! Stopping reconstruction!");
      AliWarning("No trigger can be loaded! The trigger information will not be used!");
      return kFALSE;
    }
  }

  AliTriggerConfiguration *config = aCTP->GetConfiguration();
  if (!config) {
    AliError("No trigger configuration has been found! The trigger configuration information will not be used!");
    if (fRawReader) delete aCTP;
    return kFALSE;
  }
  UChar_t clustmask = 0;
  TString trclasses;
  ULong64_t trmask = fEventInfo.GetTriggerMask();
  ULong64_t trmaskNext50 = fEventInfo.GetTriggerMaskNext50();
  const TObjArray& classesArray = config->GetClasses();
  Int_t nclasses = classesArray.GetEntriesFast();
  for( Int_t iclass=0; iclass < nclasses; iclass++ ) {
    AliTriggerClass* trclass = (AliTriggerClass*)classesArray.At(iclass);
    Int_t trindex = trclass->GetIndex()-1;
    if (fesd) fesd->SetTriggerClass(trclass->GetName(),trindex);
    Bool_t match = trindex<50 ? (trmask & (1ull << trindex)) : (trmaskNext50 & (1ull << (trindex-50)));
    if (!match) continue;
    trclasses += " ";
    trclasses += trclass->GetName();
    trclasses += " ";
    clustmask |= trclass->GetCluster()->GetClusterMask();
    if (TriggerMatches2Alias(trclass->GetName(),fCosmicAlias)) fEventInfo.SetCosmicTrigger(kTRUE);
    else if (TriggerMatches2Alias(trclass->GetName(),fLaserAlias))  fEventInfo.SetCalibLaserTrigger(kTRUE);
    else fEventInfo.SetBeamTrigger(kTRUE);
  }
  fEventInfo.SetTriggerClasses(trclasses);

  // Write names of active trigger inputs in ESD Header
  const TObjArray& inputsArray = config->GetInputs(); 
  Int_t ninputs = inputsArray.GetEntriesFast();
  for( Int_t iinput=0; iinput < ninputs; iinput++ ) {
    AliTriggerInput* trginput = (AliTriggerInput*)inputsArray.At(iinput);
    if (trginput && trginput->GetMask()>0) {
      Int_t inputIndex = (Int_t)TMath::Nint(TMath::Log2(trginput->GetMask()));
      AliESDHeader* headeresd = 0x0;
      if (fesd) headeresd = fesd->GetHeader();
      if (headeresd) {
	Int_t trglevel = (Int_t)trginput->GetLevel();
	if (trglevel == 0) headeresd->SetActiveTriggerInputs(trginput->GetInputName(), inputIndex);
	if (trglevel == 1) headeresd->SetActiveTriggerInputs(trginput->GetInputName(), inputIndex+24);
	if (trglevel == 2) headeresd->SetActiveTriggerInputs(trginput->GetInputName(), inputIndex+48);
      }
    }
  }

  // Set the information in ESD
  if (fesd) {
    fesd->SetTriggerMask(trmask);
    fesd->SetTriggerMaskNext50(trmaskNext50);
    fesd->SetTriggerCluster(clustmask);
  }

  if (!aCTP->CheckTriggeredDetectors()) {
    if (fRawReader) delete aCTP;
    return kFALSE;
  }    

  if (fRawReader) delete aCTP;

  // We have to fill also the HLT decision here!!
  // ...
  // check if event has cosmic or laser alias
  

  return kTRUE;
}

const char *AliReconstruction::MatchDetectorList(const char *detectorList, UInt_t detectorMask)
{
  // Match the detector list found in the rec.C or the default 'ALL'
  // to the list found in the GRP (stored there by the shuttle PP which
  // gets the information from ECS)
  static TString resultList;
  TString detList = detectorList;

  resultList = "";

  for(Int_t iDet = 0; iDet < (AliDAQ::kNDetectors-1); iDet++) {
    if ((detectorMask >> iDet) & 0x1) {
      TString det = AliDAQ::OfflineModuleName(iDet);
      if ((detList.CompareTo("ALL") == 0) ||
	  ((detList.BeginsWith("ALL ") ||
	    detList.EndsWith(" ALL") ||
	    detList.Contains(" ALL ")) &&
	   !(detList.BeginsWith("-"+det+" ") ||
	     detList.EndsWith(" -"+det) ||
	     detList.Contains(" -"+det+" "))) ||
	  (detList.CompareTo(det) == 0) ||
	  detList.BeginsWith(det+" ") ||
	  detList.EndsWith(" "+det) ||
	  detList.Contains( " "+det+" " )) {
	if (!resultList.EndsWith(det + " ")) {
	  resultList += det;
	  resultList += " ";
	}
      }	       
    }
  }

  // HLT
  if ((detectorMask >> AliDAQ::kHLTId) & 0x1) {
    TString hltDet = AliDAQ::OfflineModuleName(AliDAQ::kNDetectors-1);
    if ((detList.CompareTo("ALL") == 0) ||
	((detList.BeginsWith("ALL ") ||
	  detList.EndsWith(" ALL") ||
	  detList.Contains(" ALL ")) &&
	 !(detList.BeginsWith("-"+hltDet+" ") ||
	   detList.EndsWith(" -"+hltDet) ||
	   detList.Contains(" -"+hltDet+" "))) ||
	(detList.CompareTo(hltDet) == 0) ||
	detList.BeginsWith(hltDet+" ") ||
	detList.EndsWith(" "+hltDet) ||
	detList.Contains( " "+hltDet+" " )) {
      resultList += hltDet;
    }
  }

  return resultList.Data();

}

//______________________________________________________________________________
void AliReconstruction::Abort(const char *method, EAbort what)
{
  // Abort processing. If what = kAbortProcess, the Process() loop will be
  // aborted. If what = kAbortFile, the current file in a chain will be
  // aborted and the processing will continue with the next file, if there
  // is no next file then Process() will be aborted. Abort() can also  be
  // called from Begin(), SlaveBegin(), Init() and Notify(). After abort
  // the SlaveTerminate() and Terminate() are always called. The abort flag
  // can be checked in these methods using GetAbort().
  //
  // The method is overwritten in AliReconstruction for better handling of
  // reco specific errors 

  if (!fStopOnError) return;

  CleanUp();

  TString whyMess = method;
  whyMess += " failed! Aborting...";

  AliError(whyMess.Data());

  fAbort = what;
  TString mess = "Abort";
  if (fAbort == kAbortProcess)
    mess = "AbortProcess";
  else if (fAbort == kAbortFile)
    mess = "AbortFile";

  Info(mess.Data(), "%s", whyMess.Data());
}

//______________________________________________________________________________
Bool_t AliReconstruction::ProcessEvent(void* event)
{
  // Method that is used in case the event loop
  // is steered from outside, for example by AMORE
  // 'event' is a pointer to the DATE event in the memory

  if (fRawReader) delete fRawReader;
  fRawReader = new AliRawReaderDate(event);
  fStatus = ProcessEvent(fRunLoader->GetNumberOfEvents());  
  delete fRawReader;
  fRawReader = NULL;

  return fStatus;
}

//______________________________________________________________________________
Bool_t AliReconstruction::ParseOutput()
{
  // The method parses the output file
  // location string in order to steer
  // properly the selector

  TPMERegexp re1("(\\w+\\.zip#\\w+\\.root):([,*\\w+\\.root,*]+)@dataset://(\\w++)");
  TPMERegexp re2("(\\w+\\.root)?@?dataset://(\\w++)");

  if (re1.Match(fESDOutput) == 4) {
    // root archive with output files stored and regustered
    // in proof dataset
    gProof->AddInput(new TNamed("PROOF_OUTPUTFILE",re1[1].Data()));
    gProof->AddInput(new TNamed("PROOF_OUTPUTFILE_LOCATION",re1[3].Data()));
    gProof->AddInput(new TNamed("PROOF_OUTPUTFILE_DATASET",""));
    gProof->AddInput(new TNamed("PROOF_OUTPUTFILE_ARCHIVE",re1[2].Data()));
    AliInfo(Form("%s files will be stored within %s in dataset %s",
		 re1[2].Data(),
		 re1[1].Data(),
		 re1[3].Data()));
  }
  else if (re2.Match(fESDOutput) == 3) {
    // output file stored and registered
    // in proof dataset
    gProof->AddInput(new TNamed("PROOF_OUTPUTFILE",(re2[1].IsNull()) ? "AliESDs.root" : re2[1].Data()));
    gProof->AddInput(new TNamed("PROOF_OUTPUTFILE_LOCATION",re2[2].Data()));
    gProof->AddInput(new TNamed("PROOF_OUTPUTFILE_DATASET",""));
    AliInfo(Form("%s will be stored in dataset %s",
		 (re2[1].IsNull()) ? "AliESDs.root" : re2[1].Data(),
		 re2[2].Data()));
  }
  else {
    if (fESDOutput.IsNull()) {
      // Output location not given.
      // Assuming xrootd has been already started and
      // the output file has to be sent back
      // to the client machine
      TString esdUrl(Form("root://%s/%s/",
			  TUrl(gSystem->HostName()).GetHostFQDN(),
			  gSystem->pwd()));
      gProof->AddInput(new TNamed("PROOF_OUTPUTFILE","AliESDs.root"));
      gProof->AddInput(new TNamed("PROOF_OUTPUTFILE_LOCATION",esdUrl.Data()));
      AliInfo(Form("AliESDs.root will be stored in %s",
		   esdUrl.Data()));
    }
    else {
      // User specified an output location.
      // Ones has just to parse it here
      TUrl outputUrl(fESDOutput.Data());
      TString outputFile(gSystem->BaseName(outputUrl.GetFile()));
      gProof->AddInput(new TNamed("PROOF_OUTPUTFILE",outputFile.IsNull() ? "AliESDs.root" : outputFile.Data()));
      TString outputLocation(outputUrl.GetUrl());
      outputLocation.ReplaceAll(outputFile.Data(),"");
      gProof->AddInput(new TNamed("PROOF_OUTPUTFILE_LOCATION",outputLocation.Data()));
      AliInfo(Form("%s will be stored in %s",
		   outputFile.IsNull() ? "AliESDs.root" : outputFile.Data(),
		   outputLocation.Data()));
    }
  }

  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliReconstruction::IsHighPt() const {
  // Selection of events containing "high" pT tracks
  // If at least one track is found within 1.5 and 100 GeV (pT)
  // that was reconstructed by both ITS and TPC, the event is accepted

  // Track cuts
  const Double_t pTmin = 1.5;
  const Double_t pTmax = 100;
  ULong_t mask = 0;
  //  mask |= (AliESDtrack::kITSrefit);
  //  mask |= (AliESDtrack::kTPCrefit);
  mask |= (AliESDtrack::kITSin);
  mask |= (AliESDtrack::kTPCin);
  const Double_t pTminCosmic = 5.;
  const Double_t pTmaxCosmic = 100;
  ULong_t maskCosmic = 0;
  Int_t cosmicCount=0;
  //  maskCosmic |= (AliESDtrack::kTPCrefit);
  maskCosmic |= (AliESDtrack::kTPCin);

  Bool_t isOK = kFALSE;

  if (fesd && fesd->GetEventType()==AliRawEventHeaderBase::kPhysicsEvent) {
    // Check if this ia a physics event (code 7)
    Int_t ntrk = fesd->GetNumberOfTracks();
    for (Int_t itrk=0; itrk<ntrk; ++itrk) {
	  
      AliESDtrack * trk = fesd->GetTrack(itrk);
      if (trk 
	  && trk->Pt() > pTmin 
	  && trk->Pt() < pTmax
	  && (trk->GetStatus() & mask) == mask ) {
	
	isOK = kTRUE;
	break;
      }
      if (trk 
	  && trk->GetInnerParam()
	  && trk->GetInnerParam()->Pt() > pTminCosmic 
	  && trk->GetInnerParam()->Pt() < pTmaxCosmic
	  && (trk->GetStatus() & maskCosmic) == maskCosmic ) {
	
	cosmicCount++;
	break;
      }
    }
    if (cosmicCount>1) isOK=kTRUE;
  }
  return isOK;
}

//______________________________________________________________________________
Bool_t AliReconstruction::IsCosmicOrCalibSpecie() const {
  // Select cosmic or calibration events

  Bool_t isOK = kFALSE;

  if (fesd && fesd->GetEventType()==AliRawEventHeaderBase::kPhysicsEvent) {
      // Check if this ia a physics event (code 7)
      
      UInt_t specie = fesd->GetEventSpecie();
      if (specie==AliRecoParam::kCosmic || specie==AliRecoParam::kCalib) {
	isOK = kTRUE;
      }
  }
  return isOK;
}


//______________________________________________________________________________
void AliReconstruction::ResetFriends() 
{
  // destroy friends taking care to release objects eventually owned by detectors
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if (!fTracker[iDet]) continue;
    if (fTracker[iDet]->OwnsESDObjects()) fTracker[iDet]->CleanESDFriendsObjects(fesd);
  }
  //
  if (!fesdf) return;
  // RS now friends for writing contain shallow copy of non-persistent friend tracks of AliESDtracks
  fesdf->ResetSoft(); // soft reset: Tracks are Cleared rather then Deleted
  //  fesdf->~AliESDfriend();
  //  new (fesdf) AliESDfriend(); // Reset...
}

//______________________________________________________________________________
void AliReconstruction::WriteESDfriend() 
{
  // Fill the ESD friend in the tree. 
  if (!fWriteThisFriend) ftreeF->SetBranchAddress("ESDfriend.",&fesdfDummy);
  //
  Long64_t nbf = ftreeF->Fill();
  if (fTreeBuffSize>0 && ftreeF->GetAutoFlush()<0 && (fMemCountESDF += nbf)>fTreeBuffSize ) { // default limit is still not reached
    nbf = ftreeF->GetZipBytes();
    if (nbf>0) nbf = -nbf;
    else       nbf = ftreeF->GetEntries();
    ftreeF->SetAutoFlush(nbf);
    AliInfo(Form("Calling ftreeF->SetAutoFlush(%lld) | W:%lld T:%lld Z:%lld",
		 nbf,fMemCountESDF,ftreeF->GetTotBytes(),ftreeF->GetZipBytes()));        
  }
  if (!fWriteThisFriend) ftreeF->SetBranchAddress("ESDfriend.",&fesdf); // restore real friend
}

//_________________________________________________________________
void AliReconstruction::DeleteDigits(const TString& detectors)
{
  // delete requested digit files produced at current event
  static int iEvent = 0;
  if (detectors.IsNull()) return;
  TString detStr = detectors;
  AliInfo(Form("Deleting Digits: %s",detectors.Data()));

  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if(!IsSelected(fgkDetectorName[iDet], detStr)) continue;
    unlink(Form("%s.Digits.root",fgkDetectorName[iDet]));
  }
  AliSysInfo::AddStamp(Form("DelDigits_%d",iEvent), 0,0,iEvent);
  iEvent++;
}

//_________________________________________________________________
void AliReconstruction::DeleteRecPoints(const TString& detectors)
{
  // delete requested recpoint files produced at current event
  static int iEvent = 0;
  if (detectors.IsNull()) return;
  TString detStr = detectors;
  AliInfo(Form("Deleting Recpoints: %s",detectors.Data()));
  //
  for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
    if(!IsSelected(fgkDetectorName[iDet], detStr)) continue;
    unlink(Form("%s.RecPoints.root",fgkDetectorName[iDet]));
  }
  AliSysInfo::AddStamp(Form("DelRecPoints_%d",iEvent), 0,0,iEvent);
  iEvent++;
}

//_________________________________________________________________
void AliReconstruction::SetStopOnResourcesExcess(Int_t vRSS,Int_t vVMEM)
{
  // require checking the resources left and stopping on excess
  // if 0  : no check is done
  // if >0 : stop reconstruction if exceeds this value
  // if <0 : use as margin to system limits
  //
  const int kKB2MB = 1024;
  const int kInfMem = 9999999;
  //
  struct rlimit r;
  int pgSize = getpagesize();
  //
  if (vRSS>0) {
    fMaxRSS = vRSS;
    AliInfo(Form("Setting max. RSS usage to user value %d MB",fMaxRSS));
  }
  else if (vRSS<0) {
    getrlimit(RLIMIT_RSS,&r);
    fMaxRSS = r.rlim_max==RLIM_INFINITY ? kInfMem : int(r.rlim_max*pgSize/kKB2MB/kKB2MB) + vRSS;
    AliInfo(Form("Setting max. RSS usage to system hard limit %d%s MB (%d margin)",fMaxRSS,r.rlim_max==RLIM_INFINITY ? "(inf)":"",-vRSS));
  }
  else {AliInfo("No check on RSS memory usage will be applied");}
  //
  if (vVMEM>0) {
    fMaxVMEM = vVMEM;
    AliInfo(Form("Setting max. VMEM usage to user value %d MB",fMaxVMEM));
  }
  else if (vVMEM<0) {
    getrlimit(RLIMIT_AS,&r);
    fMaxVMEM = r.rlim_max==RLIM_INFINITY ? kInfMem : int(r.rlim_max*pgSize/kKB2MB/kKB2MB) + vVMEM;
    AliInfo(Form("Setting max. VMEM usage to system hard limit %d%s MB (%d margin)",fMaxVMEM,r.rlim_max==RLIM_INFINITY ? "(inf)":"",-vVMEM));
  }
  else {AliInfo("No check on RSS memory usage will be applied");}
  //  
}

//_________________________________________________________________
Bool_t AliReconstruction::HasEnoughResources(int ev)
{
  // check if process consumed more than allowed resources
  const int kKB2MB = 1024;
  Bool_t res = kTRUE;
  if (!fMaxRSS && !fMaxVMEM) return res;
  //
  ProcInfo_t procInfo;
  gSystem->GetProcInfo(&procInfo);
  if (procInfo.fMemResident/kKB2MB > fMaxRSS)  res = kFALSE;
  if (procInfo.fMemVirtual/kKB2MB  > fMaxVMEM) res = kFALSE;  
  //
  if (!res) {
    AliInfo(Form("Job exceeded allowed limits: RSS:%d (%d) VMEM:%d (%d), will stop",
		 int(procInfo.fMemResident/kKB2MB),fMaxRSS,
		 int(procInfo.fMemVirtual/kKB2MB) ,fMaxVMEM));
    //
    unlink(Form("%s",fgkStopEvFName));
    ofstream outfile(fgkStopEvFName);
    outfile << ev << std::endl;
    outfile.close();
    fStopped = kTRUE;
  }
  return res;
}

Bool_t AliReconstruction::HasNextEventAfter(Int_t eventId)
{
	 return ( (eventId < fRunLoader->GetNumberOfEvents()) ||
	   (fRawReader && fRawReader->NextEvent()) );
}

//_________________________________________________________________
void AliReconstruction::CheckRecoCDBvsSimuCDB()
{
  // if some CDB entries must be the same in the simulation
  // and reconstruction, check here
  int nent = fCheckRecoCDBvsSimuCDB.GetEntriesFast();
  AliInfo(Form("Check %d entries for matching between sim and rec",nent));
  //
  // get simulation CDB
  fRunLoader->CdGAFile();
  TMap*  cdbMapSim  = (TMap*)gDirectory->Get("cdbMap");
  TList* cdbListSim = (TList*)gDirectory->Get("cdbList");
  if (!(cdbMapSim && cdbListSim)) {
    AliInfo(Form("No CDBMap/List found in %s, nothing to check",fGAliceFileName.Data()));
    return;
  }
  // read the requested objects to make sure they will appear in the reco list
  for (Int_t i=0;i<nent;i++) {
    TNamed* cdbent = (TNamed*) fCheckRecoCDBvsSimuCDB[i];
    if (!cdbent) continue;
    AliCDBManager::Instance()->Get(cdbent->GetName());
  }
  // get default path for simulation
  TPair* pair;
  TObjString* stro;
  pair = (TPair*)cdbMapSim->FindObject("default");
  if (!pair) {AliFatal("Did not find default storage used for simulations"); return;}
  TString defSimStore = ((TObjString*)pair->Value())->GetString();
  RectifyCDBurl(defSimStore);
  //
  // get reconstruction CDB
  TMap cdbMapRecP, *cdbMapRec = &cdbMapRecP;
  cdbMapRec->SetName("cdbMap");
  TList cdbListRecP, *cdbListRec = &cdbListRecP;
  cdbListRec->SetName("cdbList");
  // create map/list accounting for eventual snapshot
  AliCDBManager::Instance()->CreateMapListCopy(cdbMapRecP,cdbListRecP);
  //
  // get default path for reconstruction
  pair = (TPair*)cdbMapRec->FindObject("default");
  if (!pair) {AliFatal("Did not find default storage used for reconstruction"); return;}
  TString defRecStore = ((TObjString*)pair->Value())->GetString();
  RectifyCDBurl(defRecStore);
  //
  for (Int_t i=0;i<nent;i++) {
    TNamed* cdbent = (TNamed*) fCheckRecoCDBvsSimuCDB[i];
    if (!cdbent) continue;
    //
    AliInfo(Form("#%d Checking %s",i,cdbent->GetName()));
    //
    // find cdbID used for sim
    TString idSim="",storSim="";
    TIter nextSim(cdbListSim);
    while ((stro=(TObjString*)nextSim())) {
      if (stro->GetString().Contains(cdbent->GetName())) {
	idSim = stro->GetString();
	break;
      }
    }    
    // find the storage used for sim
    // check in the simuCDB special paths
    pair = (TPair*)cdbMapSim->FindObject(cdbent->GetName());
    if (pair) { // specific path is used
      storSim = ((TObjString*)pair->Value())->GetString();
      RectifyCDBurl(storSim);
    }
    else storSim = defSimStore;  // default storage list is used
    //
    if (!idSim.IsNull()) AliInfo(Form("Sim. used %s from %s",idSim.Data(), storSim.Data()));
    else                 AliInfo("Sim. did not use this object");
    //
    // find cdbID used for rec
    TString idRec="",storRec="";
    TIter nextRec(cdbListRec);
    while ((stro=(TObjString*)nextRec())) {
      if (stro->GetString().Contains(cdbent->GetName())) {
	idRec = stro->GetString();
	break;
      }
    }
    //
    // find storage used for the rec
    pair = (TPair*)cdbMapRec->FindObject(cdbent->GetName());
    if (pair) {  // specific path is used
      storRec = ((TObjString*)pair->Value())->GetString();
      RectifyCDBurl(storRec);
    }
    else storRec = defRecStore; // default storage list is used
    //
    if (!idRec.IsNull()) AliInfo(Form("Rec. used %s from %s",idRec.Data(), storRec.Data()));
    else                 AliInfo("Rec. did not use this object");
    //
    if (!idSim.IsNull() && !idRec.IsNull() && ((idSim!=idRec) || (storSim!=storRec)) ) 
      AliFatal("Different objects were used in sim and rec");
  }
  
}

//_________________________________________________________
void AliReconstruction::RectifyCDBurl(TString& url)
{
  // TBD RS
  // remove everything but the url
  TString sbs;
  if (!(sbs=url("\\?User=[^?]*")).IsNull())                url.ReplaceAll(sbs,"");
  if (!(sbs=url("\\?DBFolder=[^?]*")).IsNull())            url.ReplaceAll("?DB","");
  if (!(sbs=url("\\?SE=[^?]*")).IsNull())                  url.ReplaceAll(sbs,"");
  if (!(sbs=url("\\?CacheFolder=[^?]*")).IsNull())         url.ReplaceAll(sbs,"");
  if (!(sbs=url("\\?OperateDisconnected=[^?]*")).IsNull()) url.ReplaceAll(sbs,"");
  if (!(sbs=url("\\?CacheSize=[^?]*")).IsNull())           url.ReplaceAll(sbs,"");  
  if (!(sbs=url("\\?CleanupInterval=[^?]*")).IsNull())     url.ReplaceAll(sbs,"");  
  Bool_t slash=kFALSE,space=kFALSE;
  while ( (slash=url.EndsWith("/")) || (space=url.EndsWith(" ")) ) {
    if (slash) url = url.Strip(TString::kTrailing,'/');
    if (space) url = url.Strip(TString::kTrailing,' ');
  }
  //url.ToLower();
  //
}

//_________________________________________________________
void AliReconstruction::ProcessTriggerAliases()
{
  // load trigger aliases, attach them to the reader
  //
  fDeclTriggerClasses.Clear();
  AliCentralTrigger *aCTP = NULL;
  if (fRawReader) {
    aCTP = new AliCentralTrigger();
    TString configstr("");
    if (!aCTP->LoadConfiguration(configstr)) { // Load CTP config from OCDB
      AliError("No trigger configuration found in OCDB! The trigger configuration information will not be used!");
      delete aCTP;
      return;
    }
  }
  else if (fRunLoader && (!fRunLoader->LoadTrigger())) {
    aCTP = fRunLoader->GetTrigger();
  }
  else {
    if (fStopOnMissingTriggerFile) AliFatal("No trigger can be loaded! Stopping reconstruction!");
    AliWarning("No trigger can be loaded! The trigger information will not be used!");
    return;
  }
  //  
  AliTriggerConfiguration *config = aCTP->GetConfiguration();
  if (!config) {
    AliError("No trigger configuration has been found! The trigger configuration information will not be used!");
    if (fRawReader) delete aCTP;
    return;
  }
  //
  // here we have list of active triggers
  const TObjArray& classesArray = config->GetClasses();
  Int_t nclasses = classesArray.GetEntriesFast();
  // 
  fAlias2Trigger = new THashList();
  //
  AliCDBEntry * entry = AliCDBManager::Instance()->Get("GRP/CTP/Aliases");
  if (entry) {
    THashList * lst = dynamic_cast<THashList*>(entry->GetObject());
    if (lst) {
      lst->Sort(kSortDescending); // to avoid problems with substrungs
      if (fRawReader) fRawReader->LoadTriggerAlias(lst);
      // Now declare all the triggers present in the aliases
      TIter iter(lst);
      TNamed *nmd = 0;
      while((nmd = dynamic_cast<TNamed*>(iter.Next()))) { // account aliases of this trigger >>
	fDeclTriggerClasses += " ";
	fDeclTriggerClasses += nmd->GetName();
	//
	if (!classesArray.FindObject(nmd->GetName())) continue;
	TString aliasList(nmd->GetTitle());
	TObjArray* arrAliases = aliasList.Tokenize(',');
	Int_t nAliases = arrAliases->GetEntries();
	// Loop on aliases for the current trigger
	for(Int_t i=0; i<nAliases; i++){
	  TObjString *alias = (TObjString*) arrAliases->At(i);
	  // Find the current alias in the hash list. If it is not there, add TNamed entry
	  TNamed * inlist = (TNamed*)fAlias2Trigger->FindObject((alias->GetString()).Data());
	  if (!inlist) {
	    inlist = new TNamed((alias->GetString()).Data(),Form(" %s ",nmd->GetName()));
	    fAlias2Trigger->Add(inlist);
	  }
	  else {
	    TString tt(inlist->GetTitle());
	    tt += "||";
	    tt += Form(" %s ",nmd->GetName());
	    inlist->SetTitle(tt.Data());
	  }
	}
	delete arrAliases;
      }  // account aliases of this trigger <<
    }
    else AliError("Cannot cast the object with trigger aliases to THashList!");
  }
  else AliError("No OCDB ebtry for the trigger aliases!");
  //
  AliInfo("Aliases defined:");
  fAlias2Trigger->Print();
  //
  if (fRawReader) {
    // active classes mentioned in the alias will be converted to their masks
    for( Int_t iclass=0; iclass < nclasses; iclass++ ) {
      AliTriggerClass* trclass = (AliTriggerClass*)classesArray.At(iclass);
      if (!trclass) continue;
      int trindex = trclass->GetIndex()-1;
      fRawReader->LoadTriggerClass(trclass->GetName(),trindex);      
    }
    //
    // nullify all remaining triggers mentioned in the alias
    if (!fDeclTriggerClasses.IsNull()) {
      TObjArray *tokens = fDeclTriggerClasses.Tokenize(" ");
      Int_t ntokens = tokens->GetEntriesFast();
      for (Int_t itoken = 0; itoken < ntokens; ++itoken) {
	fRawReader->LoadTriggerClass((((TObjString*)tokens->At(itoken))->String()).Data(),-1);
      }
      delete tokens;
    }
    //
    // make sure no unparsed triggers names left in the requested triggers
    if (TPRegexp("[A-Za-z]").MatchB(fRawReader->GetParsedTriggerExpression())) {
      AliFatalF("Unknown triggers found in requested list: %s",fRawReader->GetParsedTriggerExpression().Data());
    }
  }
  //
}

//___________________________________________________
Bool_t AliReconstruction::TriggerMatches2Alias(const char* trigName, const char* alias)
{
  // check if trigger matches to alias
  TString trName = trigName;
  if (!fAlias2Trigger) return kFALSE;
  TNamed* al = (TNamed*)fAlias2Trigger->FindObject(alias);
  if (!al) return kFALSE;
  TString altrig = al->GetTitle();
  return altrig.Contains(Form(" %s ",trigName));
  //
}

//___________________________________________________
Bool_t AliReconstruction::DecideFriendsStorage()
{
  // Decide if and what should be stored in the friends
  // Decision is taken on 
  // 1) even level: unconditionally for cosmic and calib event,
  // priority for events with high pt tracks
  // 2) optionally on track level: priority according to kPriority.. flags below
  const Double_t pTminHigh = 1.5, pTmaxHigh = 100;
  const ULong_t kMaskHighPt = AliESDtrack::kITSout|AliESDtrack::kTPCin;
  const ULong_t kPriorityFlag[] =
  {    AliESDtrack::kTPCin
      ,AliESDtrack::kITSpureSA|AliESDtrack::kITSout
      ,AliESDtrack::kITSout|AliESDtrack::kTPCin
      ,AliESDtrack::kITSout|AliESDtrack::kTPCin|AliESDtrack::kTRDout
      ,AliESDtrack::kITSout|AliESDtrack::kTPCin|AliESDtrack::kTOFout
      ,AliESDtrack::kITSout|AliESDtrack::kTPCin|AliESDtrack::kTRDout|AliESDtrack::kTOFout
      ,AliESDtrack::kITSout|AliESDtrack::kTPCin|AliESDtrack::kPHOSmatch
      ,AliESDtrack::kITSout|AliESDtrack::kTPCin|AliESDtrack::kEMCALmatch
      ,AliESDtrack::kITSout|AliESDtrack::kTPCin|AliESDtrack::kHMPIDout
      };
  const int kNPrio = sizeof(kPriorityFlag)/sizeof(ULong_t);
  const ULong_t highPtMask = AliESDtrack::kITSin|AliESDtrack::kTPCin;
  float fracHighPt = 0, nHighPtCheck = 0;
  // 
  Bool_t isSelected = kFALSE;
  int ntrk = fesd->GetNumberOfTracks();
  fesd->SetNTPCFriend2Store(ntrk); // by default - all tracks
  //
  // Unconditionally store Field off, Cosmic and Calib events
  double bz = AliTrackerBase::GetBz();
  if (TMath::Abs(bz)<0.5) {
    isSelected=kTRUE;
    AliInfo("AliESDfriends event stored: non-standard field event");
    return isSelected;
  }
  if (fEventInfo.HasCosmicTrigger() || fEventInfo.HasCalibLaserTrigger()) {
    //      if (IsCosmicOrCalibSpecie()) {    
    isSelected = kTRUE; // Selection of calib or cosmic events
    AliInfo("AliESDfriends event stored: calibration or cosmic event");
    return isSelected;
  }
  //
  // analyze event on tracks level
  static TArrayS weights;
  static TArrayI indices;
  Bool_t sparsify = 
    (fMaxFriendTracks>0 && (ntrk > fMaxFriendTracks)) ||
    (fSkipFriendsForLargeZ && fSkipFriendsCutZ>0); // applied to no-its tracks
  //
  if (sparsify && ntrk>weights.GetSize()) { // arrays are used only in case of sparsification
    weights.Set(ntrk+100);
    indices.Set(ntrk+100);
  }
  Short_t *pWeights = weights.GetArray(); // for fast access
  Int_t   *pIndices = indices.GetArray();
  //
  for (Int_t itrk=0; itrk<ntrk; ++itrk) {
    //	  
    AliESDtrack * trk = fesd->GetTrack(itrk);
    ULong_t status = trk->GetStatus();
    //
    // tag high pt tracks
    Bool_t isHighPt = kFALSE;
    if ( (status&kMaskHighPt) == kMaskHighPt ) {      
      double pt = trk->Pt();
      isHighPt =  (pt>pTminHigh && pt<pTmaxHigh);
      nHighPtCheck++;
      if (isHighPt) fracHighPt++;
    }
    //
    if (sparsify) {
      int pri = 0;    // define track priority
      for (pri=kNPrio;pri--;) {
	ULong_t priFlags=kPriorityFlag[pri],priFlagsOK = priFlags&status;
	// stupidly, from the from the status flags it is not clear if there are TRD tracklets matched
	if (priFlags&AliESDtrack::kTRDout && !trk->GetTRDntracklets()) continue;
	if (priFlags&AliESDtrack::kTPCout && trk->GetTPCNcls()<16) continue; // disregard TPC tracks with low Nclus
	if (priFlagsOK==kPriorityFlag[pri]) break;
      }
      pWeights[itrk] = ++pri;      // to have it positive
      if (isHighPt) pWeights[itrk] += kNPrio;  // high pt tracks have highest priority
      else {
	if (fSkipFriendsCutZ>0 && !(status&AliESDtrack::kITSin) && (status&AliESDtrack::kTPCin)) {
	  double z;
	  if (!trk->GetInnerParam()->GetZAt(0,bz,z) || TMath::Abs(z)>fSkipFriendsCutZ) pWeights[itrk] = -kNPrio;
	}
      }
    }
    //
  }
  if (nHighPtCheck>0) fracHighPt = fracHighPt/nHighPtCheck; 
  //
  // decide if we store the friend event at all
  int nkept=0, nkeptTPC=0;
  if (gRandom->Rndm()/(1.+fracHighPt) < fFractionFriends) {
    isSelected = kTRUE;
    if (sparsify) {
      TMath::Sort(ntrk,pWeights,pIndices); // sort in order of decreasing priorities
      // flag low-priority tracks not not store friends
      //
      for (int itrk=0;itrk<ntrk;itrk++) {
	int ind = pIndices[itrk];
	AliESDtrack* trk = fesd->GetTrack(ind);
	if (nkept>=fMaxFriendTracks) { // no room anymore
	  trk->SetFriendNotStored(kTRUE); 
	  //AliInfoF("FrOff %d:%d reached %d W: %d",itrk,ind, nkept,pWeights[ind]);
	  continue; 
	}
	// did not reach the limit yet
	if (pWeights[ind]>0 || !fSkipFriendsForLargeZ) { // in principle, worth storing
	  nkept++;
	  if (trk->GetTPCNcls()>15) nkeptTPC++;
	  //AliInfoF("FrON  %d:%d saved %d/%d W: %d",itrk,ind, nkept,nkeptTPC,pWeights[ind]);
	  continue;
	}
	// negative weight and strict rejection of such track requested
	//AliInfoF("FrOff  %d:%d saved %d/%d W: %d",itrk,ind, nkept,nkeptTPC,pWeights[ind]);
	fesd->GetTrack(ind)->SetFriendNotStored(kTRUE); // unconditional rejections
	//
      }
    }
    else nkept = nkeptTPC = ntrk;
  }
  //
  if (isSelected) AliInfoF("AliESDfriends event stored: %d (with TPC: %d) out of %d friend tracks will be kept",
			   nkept,nkeptTPC,ntrk);
  else            AliInfo("AliESDfriends will not be stored for this event");
  //
  fesd->SetNTPCFriend2Store(nkeptTPC);
  //
  return isSelected;
}
