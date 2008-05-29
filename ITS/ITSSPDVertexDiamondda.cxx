/*
Contact: cvetan.cheshkov@cern.ch
Link: missing
Run Type: PHYSICS
DA Type: MON
Number of events needed: 10000
Input Files:
Output Files:
Trigger types used: PHYSICS
*/

#define X_LIMIT 2.0
#define Y_LIMIT 2.0
#define Z_LIMIT 50.0
#define X_DELTA 0.02
#define Y_DELTA 0.02
#define Z_DELTA 0.2
#define OUTPUT_FILE "SPDVertexDiamondDA.root"
#define CDB_STORAGE "local://$ALICE_ROOT"
#define N_EVENTS_AUTOSAVE 50
#define NFITPARAMS 5

extern "C" {
#include "daqDA.h"
}

#include "event.h"
#include "monitor.h"

#ifdef ALI_AMORE
#include <AmoreDA.h>
int amore::da::Updated(char const*) {}
#endif

#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TF2.h>
#include <TFile.h>
#include <TPluginManager.h>
#include <TROOT.h>
#include <TFitter.h>

#include "AliRawReaderDate.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliESDVertex.h"
#include "AliITSDetTypeRec.h"
#include "AliITSInitGeometry.h"
#include "AliITSVertexer3DTapan.h"

TF2 *fitFcn = 0x0;

AliESDVertex* FitVertexDiamond(TH2F *hXY, TH1F *hZ);

int main(int argc, char **argv) {

  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()"); 
  TFitter *minuitFit = new TFitter(NFITPARAMS);
  TVirtualFitter::SetFitter(minuitFit);

  int status;
  if (argc<2) {
    printf("Wrong number of arguments\n");
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

  // Histograms initialization
  TH2F *hXY = new TH2F("hXY","Vertex Diamond (Y vs X)",
		       2*(Int_t)(X_LIMIT/X_DELTA),-X_LIMIT,X_LIMIT,
		       2*(Int_t)(Y_LIMIT/Y_DELTA),-Y_LIMIT,Y_LIMIT);
  TH1F *hZ  = new TH1F("hZ"," Longitudinal Vertex Profile",
		       2*(Int_t)(Z_LIMIT/Z_DELTA),-Z_LIMIT,Z_LIMIT);

  // Global initializations
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(CDB_STORAGE);
  man->SetRun(runNr);
  AliGeomManager::LoadGeometry("geometry.root");
  AliGeomManager::ApplyAlignObjsFromCDB("ITS");

  // ITS initializations
  AliITSInitGeometry initgeom;
  AliITSgeom *geom = initgeom.CreateAliITSgeom();
  printf("Geometry name: %s\n",(initgeom.GetGeometryName()).Data());

  AliITSDetTypeRec *detTypeRec = new AliITSDetTypeRec();
  detTypeRec->SetITSgeom(geom);
  detTypeRec->SetDefaults();
  detTypeRec->SetDefaultClusterFindersV2(kTRUE);

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

      // Run SPD cluster finder
      TTree* clustersTree = new TTree("TreeR", "Reconstructed Points Container"); //make a tree
      detTypeRec->DigitsToRecPoints(rawReader,clustersTree,"SPD");

      // Run vertex-finder
      AliITSVertexer3DTapan *vertexer = new AliITSVertexer3DTapan(geom,1000);
      vertexer->LoadClusters(clustersTree);
      AliESDVertex *vtx = new AliESDVertex();
      vertexer->FindVertexForCurrentEvent(vtx);

      if (TMath::Abs(vtx->GetChi2()) < 0.1) {
	// Fill the vertex into the histos
	nevents_with_vertex++;
	hXY->Fill(vtx->GetXv(),vtx->GetYv());
	hZ->Fill(vtx->GetZv());

	// Auto save
	if ((nevents_with_vertex%N_EVENTS_AUTOSAVE) == 0) {
	  TFile outFile(OUTPUT_FILE, "update");
	  AliESDVertex *fitVtx = FitVertexDiamond(hXY,hZ);
	  if (fitVtx) {
	    fitVtx->Write(fitVtx->GetName(),TObject::kOverwrite);
	    TH1 *fithXY = fitFcn->CreateHistogram();
	    fithXY->Write(fithXY->GetName(),TObject::kOverwrite);
	    delete fithXY;
	  }
	  hXY->Write(hXY->GetName(),TObject::kOverwrite);
	  hZ->Write(hZ->GetName(),TObject::kOverwrite);
	  outFile.Close();
	  delete fitVtx;
	}
      }

      delete vtx;
      delete vertexer;

      delete clustersTree;
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
  
  if (detTypeRec) delete detTypeRec;

  // Store the final histograms
  TFile outFile(OUTPUT_FILE, "update");
  if (nevents_with_vertex > N_EVENTS_AUTOSAVE) { 
    // Fit XY & Z histograms
    AliESDVertex *fitVtx = FitVertexDiamond(hXY,hZ);
    if (fitVtx) {
      fitVtx->Write(fitVtx->GetName(),TObject::kOverwrite);
      TH1 *fithXY = fitFcn->CreateHistogram();
      fithXY->Write(fithXY->GetName(),TObject::kOverwrite);
      delete fithXY;
    }
    delete fitVtx;
  }
  hXY->Write(hXY->GetName(),TObject::kOverwrite);
  hZ->Write(hZ->GetName(),TObject::kOverwrite);
  outFile.Close();

#ifdef ALI_AMORE
  // send the histos to AMORE pool
  printf("AMORE send status: %d",vtxAmore.Send(hXY->GetName(),hXY));
#endif

  delete minuitFit;
  TVirtualFitter::SetFitter(0);

  /* write report */
  printf("Run #%s, received %d events with vertex, out of %d physics and out of %d total events\n",getenv("DATE_RUN_NUMBER"),nevents_with_vertex,nevents_physics,nevents_total);

  status=0;

  /* export file to FXS */
  if (daqDA_FES_storeFile(OUTPUT_FILE, "VertexDiamond")) {
    status=-2;
  }
  
  return status;
}

Double_t fitFunction(Double_t *x, Double_t *par) {

  Double_t t1 =   x[0] - par[1];
  Double_t t2 =   x[1] - par[2];

  return par[0]*TMath::Exp(-0.5*(t1*t1/(par[3]*par[3])+t2*t2/(par[4]*par[4])));
}

AliESDVertex* FitVertexDiamond(TH2F *hXY, TH1F *hZ)
{

  if (!fitFcn) {
    fitFcn = new TF2("fitFcn",fitFunction,
		     -X_LIMIT,X_LIMIT,
		     -Y_LIMIT,Y_LIMIT,NFITPARAMS);
    fitFcn->SetNpx(2*(Int_t)(X_LIMIT/X_DELTA));
    fitFcn->SetNpy(2*(Int_t)(Y_LIMIT/Y_DELTA));
    fitFcn->SetParameters(hXY->GetMaximum(),0,0,hXY->GetRMS(1),hXY->GetRMS(2));
    fitFcn->Update();
  }

  if (hXY->Fit("fitFcn","L0V+") != 0) {
    printf("XY fit failed!");
    return 0x0;
  }
  
  Double_t pos[3],poserr[3];
  pos[0] = fitFcn->GetParameter(1);
  pos[1] = fitFcn->GetParameter(2);
  poserr[0] = fitFcn->GetParameter(3);
  poserr[1] = fitFcn->GetParameter(4);

  // Could be replaced by something more robust...
  pos[2] = hZ->GetMean();
  poserr[2] = hZ->GetRMS();
 
  return new AliESDVertex(pos,poserr);
}
