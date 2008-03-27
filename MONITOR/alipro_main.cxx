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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This the source code of AliPRO (ALICE Prompt Reconstruction Online)       //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>

#ifdef ALI_DATE
#include "event.h"
#include "monitor.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TGeoManager.h"
#include "TGrid.h"
#include "AliRun.h"
#include "AliCDBManager.h"
#include "AliLog.h"
#include "AliITSRecoParam.h"
#include "AliITSReconstructor.h"
#include "AliTPCRecoParam.h"
#include "AliTPCReconstructor.h"
#include "AliMUONRecoParam.h"
#include "AliMUONReconstructor.h"
#include "AliPHOSRecoParam.h"
#include "AliPHOSReconstructor.h"
#include "AliMagFMaps.h"
#include "AliTracker.h"
#include "AliRecoParam.h"
#include "AliReconstruction.h"
#include "AliExternalTrackParam.h"
#endif

/* Main routine
      Arguments: 
      1- monitoring data source
*/
#ifdef ALI_DATE
int main(int argc, char **argv) {

  int status;
  
  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }


  /* open result file */
  FILE *fp=NULL;
  fp=fopen("./result.txt","a");
  if (fp==NULL) {
    printf("Failed to open file\n");
    return -1;
  }
  

  /* define data source : this is argument 1 */  
  status=monitorSetDataSource( argv[1] );
  if (status!=0) {
    printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
    return -1;
  }


  /* declare monitoring program */
  status=monitorDeclareMp( __FILE__ );
  if (status!=0) {
    printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }


  /* define wait event timeout - 1s max */
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);
  

  /* log start of process */
  printf("AliPRO program started\n");  


  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

  
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  // First version of the reconstruction
  // script for the FDR'08
  /////////////////////////////////////////////////////////////////////////////////////////

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://./LocalCDB");
  // Get run number from environment / GRP-in-memory
  man->SetRun(25984);
  
  // ITS settings
  AliITSRecoParam * itsRecoParam = AliITSRecoParam::GetCosmicTestParam();
  itsRecoParam->SetClusterErrorsParam(2);
  itsRecoParam->SetFindV0s(kFALSE);
  itsRecoParam->SetAddVirtualClustersInDeadZone(kFALSE);
  itsRecoParam->SetUseAmplitudeInfo(kFALSE);
  // In case we want to switch off a layer
  //  itsRecoParam->SetLayerToSkip(<N>);
  itsRecoParam->SetLayerToSkip(4);
  itsRecoParam->SetLayerToSkip(5);
  itsRecoParam->SetLayerToSkip(2);
  itsRecoParam->SetLayerToSkip(3);
  AliITSReconstructor::SetRecoParam(itsRecoParam);

  // TPC settings
  AliLog::SetClassDebugLevel("AliTPCclustererMI",2);
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kTRUE);
  tpcRecoParam->SetTimeInterval(60,940);
  tpcRecoParam->Dump();
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(1);

  // PHOS settings
  AliPHOSRecoParam* recEmc = new AliPHOSRecoParamEmc();
  recEmc->SetSubtractPedestals(kTRUE);
  recEmc->SetMinE(0.05);
  recEmc->SetClusteringThreshold(0.10);
  AliPHOSReconstructor::SetRecoParamEmc(recEmc);

  // T0 settings
  AliLog::SetModuleDebugLevel("T0", 10);

  // MUON settings
  AliLog::SetClassDebugLevel("AliMUONRawStreamTracker",3);
  AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetLowFluxParam();
  muonRecoParam->CombineClusterTrackReco(kTRUE);
  muonRecoParam->SetCalibrationMode("NOGAIN");
  //muonRecoParam->SetClusteringMode("PEAKFIT");
  //muonRecoParam->SetClusteringMode("PEAKCOG");
  muonRecoParam->Print("FULL");
  AliRecoParam::Instance()->RegisterRecoParam(muonRecoParam);
 
  // Tracking settings
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 0., 10., 2);
  AliTracker::SetFieldMap(field,1);
  Double_t mostProbPt=0.35;
  AliExternalTrackParam::SetMostProbablePt(mostProbPt);

  // Instantiation of AliReconstruction
  // AliReconstruction settings
  AliReconstruction rec;

  rec.SetUniformFieldTracking(kFALSE);
  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData();

  rec.SetRunReconstruction("ALL");
  rec.SetUseTrackingErrorsForAlignment("ITS");
  
  // In case some detectors have to be switched off...
  //  rec.SetRunLocalReconstruction("ALL");
  //  rec.SetRunTracking("ALL");
  //  rec.SetFillESD("ALL");
  // Disable vertex finder for the moment
  rec.SetRunVertexFinder(kFALSE);
  
  // To be enabled if some equipment IDs are not set correctly by DAQ
  //  rec.SetEquipmentIdMap("EquipmentIdMap.data");
  
  // Detector options if any
  rec.SetOption("ITS","cosmics,onlyITS");
  rec.SetOption("MUON","SAVEDIGITS");
  rec.SetOption("TPC","OldRCUFormat");
  rec.SetOption("PHOS","OldRCUFormat");
  
  // To be enabled when CTP readout starts
  rec.SetFillTriggerESD(kFALSE);
  
  // all events in one single file
  rec.SetNumberOfEventsPerFile(-1);
  
  // switch off cleanESD
  rec.SetCleanESD(kFALSE);
  
  rec.SetRunQA(kFALSE);
  rec.SetRunGlobalQA(kFALSE);

  void* eventPtr = NULL;
  rec.InitRun(NULL,&eventPtr);

  /* main loop (infinite) */
  for(;;) {
    struct eventHeaderStruct *event;
    eventTypeType eventT;
    
    /* get next event (blocking call until timeout) */
    status=monitorGetEventDynamic(&eventPtr);
    event=(eventHeaderStruct*)eventPtr;
    if (status==MON_ERR_EOF) {
      printf ("End of File detected\n");
      break; /* end of monitoring file has been reached */
    }
    
    if (status!=0) {
      printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
      break;
    }
    
    /* retry if got no event */
    if (event==NULL) {
      continue;
    }
    
    
    /* use event - here, just write event id to result file */
    eventT=event->eventType;
    
    if (eventT==PHYSICS_EVENT) {
      fprintf(fp,"Run #%lu, event size: %lu, BC:%u, Orbit:%u, Period:%u\n",
	      (unsigned long)event->eventRunNb,
	      (unsigned long)event->eventSize,
	      EVENT_ID_GET_BUNCH_CROSSING(event->eventId),
	      EVENT_ID_GET_ORBIT(event->eventId),
	      EVENT_ID_GET_PERIOD(event->eventId)
	      );
      nevents_physics++;
      
      AliLog::Flush();

      rec.AddEventAndRun();
     
    }
    nevents_total++;

    /* free resources */
    free(event);
    
    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      printf("EOR event detected\n");
      break;
    }
  }

  rec.FinishRun();

  /* write report */
  fprintf(fp,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  /* close result file */
  fclose(fp);

  return status;
}
#else
int main(int /*argc*/, char** /*argv*/)
{
  printf("This program was compiled without DATE\n");

  return 1;
}
#endif
