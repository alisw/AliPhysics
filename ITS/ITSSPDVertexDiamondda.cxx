/*
  Contact: cvetan.cheshkov@cern.ch
  Link: http://alisoft.cern.ch/viewvc/trunk/ITS/ITSSPDVertexDiamondda.cxx?root=AliRoot&view=log , /afs/cern.ch/user/c/cheshkov/public/ITS/VD_da_test.date , /afs/cern.ch/user/c/cheshkov/public/08000058338016.30.root.date.gz
  Reference Run: 58338
  Run Type: PHYSICS
  DA Type: MON
  Number of events needed: 100
  Input Files: GRP/Geometry/Data , ITS/Align/Data , spd_noisy_ocdb , spd_dead_ocdb , TRIGGER/SPD/PITConditions (all the files are taken from SPD daqDetDB)
  Output Files: SPDVertexDiamondDA.root
  Trigger types used: PHYSICS, SPD-F0 
*/

#define OUTPUT_FILE "SPDVertexDiamondDA.root"

extern "C" {
#include "daqDA.h"
}

#include "event.h"
#include "monitor.h"

#ifdef ALI_AMORE
#include <AmoreDA.h>
//int amore::da::Updated(char const*) {}
#endif

#include <TPluginManager.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TPaveText.h>
#include <TSystem.h>
#include <TGeoGlobalMagField.h>

#include "AliLog.h"
#include "AliMagF.h"
#include "AliRawReaderDate.h"
#include "AliCDBManager.h"
#include "AliITSMeanVertexer.h"

int main(int argc, char **argv) {

  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()"); 

  int status;
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }
  // get mean vertex quality cuts
  status = daqDA_DB_getFile("ITSSPD_VertexQualityTuning_DA.config","./ITSSPD_VertexQualityTuning_DA.config");
  if (status) {
    printf("Failed to get config file (ITSSPD_VertexQualityTuning_DA.config) from DAQ DB, status=%d\n", status);
    return -1;
  }
  /* open the config file and retrieve running parameters */
  Float_t errX  = -1;
  Float_t r     = -1;
  UInt_t minClInner = 999 ;
  UInt_t maxClInner = 999 ;
  Int_t nEvFirstLoop = 0;
  Int_t nEvAUTOSAVE = 0; 
  char name[10][10];

  FILE *fpConfig = fopen("ITSSPD_VertexQualityTuning_DA.config","r");
  fscanf(fpConfig,"%s %f\n %s %f\n %s %d\n %s %d \n %s %d \n %s %d ",&name[0], &errX, &name[1], &r, &name[2], &minClInner,&name[3], &maxClInner, &name[4],&nEvFirstLoop,&name[5],&nEvAUTOSAVE);
  fclose(fpConfig);

  printf("\n\n Mean Vertex quality cuts : \n- errX = %f\n- r = %f\n- minSPD0 = %d maxSPD0 = %d\n- nEventsFirstLoop = %d nEventsAUTOSAVE = %d \n\n\n",errX,r,minClInner,maxClInner,nEvFirstLoop,nEvAUTOSAVE);

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
  printf("Vertex-Diamond SPD DA started\n");  

  /* init some counters */
  int nevents_with_vertex = 0;
  int nevents_physics=0;
  int nevents_total=0;

  struct eventHeaderStruct *event;
  eventTypeType eventT;

  // Get run number
  if (getenv("DATE_RUN_NUMBER")==0) {
    printf("DATE_RUN_NUMBER not properly set.\n");
    return -1;
  }
  int runNr = atoi(getenv("DATE_RUN_NUMBER"));

  // Get the necessary OCDB files from the DAQ detector DB
  if (gSystem->AccessPathName("localOCDB/GRP/Geometry/Data",kFileExists)) {
    if (gSystem->mkdir("localOCDB/GRP/Geometry/Data",kTRUE) != 0) {
      printf("Failed to create directory: localOCDB/GRP/Geometry/Data");
      return -1;
    }
  }
  status = daqDA_DB_getFile("GRP/Geometry/Data","localOCDB/GRP/Geometry/Data/Run0_999999999_v0_s0.root");
  if (status) {
    printf("Failed to get geometry file (GRP/Geometry/Data) from DAQdetDB, status=%d\n", status);
    return -1;
  }

  if (gSystem->AccessPathName("localOCDB/ITS/Align/Data",kFileExists)) {
    if (gSystem->mkdir("localOCDB/ITS/Align/Data",kTRUE) != 0) {
      printf("Failed to create directory: localOCDB/ITS/Align/Data");
      return -1;
    }
  }
  status = daqDA_DB_getFile("ITS/Align/Data","localOCDB/ITS/Align/Data/Run0_999999999_v0_s0.root");
  if (status) {
    printf("Failed to get its-alignment file (ITS/Align/Data) from DAQdetDB, status=%d\n", status);
    return -1;
  }

  if (gSystem->AccessPathName("localOCDB/ITS/Calib/SPDNoisy",kFileExists)) {
    if (gSystem->mkdir("localOCDB/ITS/Calib/SPDNoisy",kTRUE) != 0) {
      printf("Failed to create directory: localOCDB/ITS/Calib/SPDNoisy");
      return -1;
    }
  }
  status = daqDA_DB_getFile("spd_noisy_ocdb","localOCDB/ITS/Calib/SPDNoisy/Run0_999999999_v0_s0.root");
  if (status) {
    printf("Failed to get spd file (spd_noisy_ocdb) from DAQdetDB, status=%d\n", status);
    return -1;
  }

  if (gSystem->AccessPathName("localOCDB/ITS/Calib/SPDDead",kFileExists)) {
    if (gSystem->mkdir("localOCDB/ITS/Calib/SPDDead",kTRUE) != 0) {
      printf("Failed to create directory: localOCDB/ITS/Calib/SPDDead");
      return -1;
    }
  }
  status = daqDA_DB_getFile("spd_dead_ocdb","localOCDB/ITS/Calib/SPDDead/Run0_999999999_v0_s0.root");
  if (status) {
    printf("Failed to get spd file (spd_dead_ocdb) from DAQdetDB, status=%d\n", status);
    return -1;
  }

  if (gSystem->AccessPathName("localOCDB/ITS/Calib/SPDSparseDead",kFileExists)) {
    if (gSystem->mkdir("localOCDB/ITS/Calib/SPDSparseDead",kTRUE) != 0) {
      printf("Failed to create directory: localOCDB/ITS/Calib/SPDSparseDead");
      return -1;
    }
  }

  status = daqDA_DB_getFile("spd_sparsedead_ocdb","localOCDB/ITS/Calib/SPDSparseDead/Run0_999999999_v0_s0.root");
  if (status) {
    printf("Failed to get spd file (spd_sparsedead_ocdb) from DAQdetDB, status=%d\n", status);
    return -1;
  }

  if (gSystem->AccessPathName("localOCDB/TRIGGER/SPD/PITConditions",kFileExists)) {
    if (gSystem->mkdir("localOCDB/TRIGGER/SPD/PITConditions",kTRUE) != 0) {
      printf("Failed to create directory: localOCDB/TRIGGER/SPD/PITConditions");
      return -1;
    }
  }
  status = daqDA_DB_getFile("TRIGGER/SPD/PITConditions","localOCDB/TRIGGER/SPD/PITConditions/Run0_999999999_v0_s0.root");
  if (status) {
    printf("Failed to get spd trigger file (TRIGGER/SPD/PITConditions) from DAQdetDB, status=%d\n", status);
    return -1;
  }
 
  status = daqDA_DB_getFile("mfchebKGI_sym.root","localOCDB/mfchebKGI_sym.root");
  if (status) {
    printf("Failed to get spd file (mfchebKGI_sym.root) from DAQdetDB, status=%d\n", status);
    return -1;
  }

  // Global initializations

  // The B filed is required in AliITSClusterFinderV2SPD
  // for the Lorentz angle correction. B set to 0.      
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 0., 0., AliMagF::k5kGUniform,AliMagF::kBeamTypeAA,-1,2,15,"localOCDB/mfchebKGI_sym.root"));  
  AliLog::SetGlobalLogLevel(AliLog::kError);
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://localOCDB");
  man->SetRun(runNr);

  // Init mean vertexer
  AliITSMeanVertexer *mv = new AliITSMeanVertexer();
  if (!mv->Init()) {
    printf("Initialization of mean vertexer object failed ! Check the log for details");
    return -1;
  }

  mv->SetCutOnErrX(errX);
  mv->SetCutOnR(r);
  mv->SetCutOnCls(minClInner,maxClInner);

  // Initialization of AMORE sender
#ifdef ALI_AMORE
  amore::da::AmoreDA vtxAmore(amore::da::AmoreDA::kSender);
#endif
  /* main loop (infinite) */
  for(;;) {
    
    /* check shutdown condition */
    if (daqDA_checkShutdown()) {break;}
    
    /* get next event (blocking call until timeout) */
    status=monitorGetEventDynamic((void **)&event);
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

    nevents_total++;
    eventT=event->eventType;
    switch (event->eventType){
      
      /* START OF RUN */
    case START_OF_RUN:
      break;
      /* END START OF RUN */
      
      /* END OF RUN */
    case END_OF_RUN:
      break;
      
    case PHYSICS_EVENT:
      nevents_physics++;
      AliRawReader *rawReader = new AliRawReaderDate((void*)event);

      // Run mean-vertexer reco
      if (mv->Reconstruct(rawReader)) nevents_with_vertex++;
      if(nevents_physics < nEvFirstLoop) continue;
      // Auto save
      if ((nevents_physics%nEvAUTOSAVE) == 0) {
        TH2F *histo = ((TH2F*)mv->GetVertexXY());
        histo->SetStats(kFALSE);
        histo->SetTitle("");
        histo->GetListOfFunctions()->SetOwner(kTRUE);
	if(histo->GetListOfFunctions()->GetEntries()<1) histo->GetListOfFunctions()->Add(new TPaveText(-5,4.5,5.,6.2,"br"));
	for(Int_t i=0; i<histo->GetListOfFunctions()->GetEntries(); i++){
	  TString funcName = histo->GetListOfFunctions()->At(i)->ClassName();
	  if(funcName.Contains("TPaveText")){
	    TPaveText *p = (TPaveText*)histo->GetListOfFunctions()->At(i);
	    p->Clear();
	    p->AddText(Form("%f events with vertex (%i out of %i processed events)",((Double_t)nevents_with_vertex)/((Double_t)nevents_physics),nevents_with_vertex,nevents_physics));
	    p->AddText(Form("%f events with good vertex (%i out of %i events with vertex)",histo->GetEntries()/((Double_t)nevents_with_vertex),(Int_t)histo->GetEntries(),nevents_with_vertex));
	  }
	}
	mv->WriteVertices(OUTPUT_FILE);

#ifdef ALI_AMORE
	// send the histos to AMORE pool
	printf("AMORE send status: %d\n",vtxAmore.Send(mv->GetVertexXY()->GetName(),mv->GetVertexXY()));
	printf("AMORE send status: %d\n",vtxAmore.Send(mv->GetVertexZ()->GetName(),mv->GetVertexZ()));
#endif
      }

      delete rawReader;
    }
    
    /* free resources */
    free(event);
    
    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      printf("EOR event detected\n");
      break;
    }
  }

  TH2F *histo = ((TH2F*)mv->GetVertexXY());
  histo->SetStats(kFALSE);
  histo->SetTitle("");
  histo->GetListOfFunctions()->SetOwner(kTRUE);
  if(histo->GetListOfFunctions()->GetEntries()<1) histo->GetListOfFunctions()->Add(new TPaveText(-5,4.5,5.,6.2,"br"));
  for(Int_t i=0; i<histo->GetListOfFunctions()->GetEntries(); i++){
    TString funcName = histo->GetListOfFunctions()->At(i)->ClassName();
    if(funcName.Contains("TPaveText")){
      TPaveText *p = (TPaveText*)histo->GetListOfFunctions()->At(i);
      p->Clear(); 
      p->AddText(Form("%f events with vertex (%i out of %i processed events)",((Double_t)nevents_with_vertex)/((Double_t)nevents_physics),nevents_with_vertex,nevents_physics));
      p->AddText(Form("%f events with good vertex (%i out of %i events with vertex)",histo->GetEntries()/((Double_t)nevents_with_vertex),(Int_t)histo->GetEntries(),nevents_with_vertex));
    }
  }
  mv->WriteVertices(OUTPUT_FILE);

#ifdef ALI_AMORE
  // send the histos to AMORE pool
  printf("AMORE send status: %d\n",vtxAmore.Send(mv->GetVertexXY()->GetName(),mv->GetVertexXY()));
  printf("AMORE send status: %d\n",vtxAmore.Send(mv->GetVertexZ()->GetName(),mv->GetVertexZ()));
#endif

  delete mv;

  /* write report */
  printf("Run #%s, received %d events with vertex, out of %d physics and out of %d total events\n",getenv("DATE_RUN_NUMBER"),nevents_with_vertex,nevents_physics,nevents_total);

  status=0;

  /* export file to FXS */
  if (daqDA_FES_storeFile(OUTPUT_FILE, "VertexDiamond")) {
    status=-2;
  }
  
  return status;
}
