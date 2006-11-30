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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This is the class for perfoming the monitoring process.                  //
//  It checks if a raw data file exists, loops over the events in the raw    //
//  data file, reconstructs TPC and ITS clusters and tracks, fills the       //
//  monitor histograms and sends the updated histograms to the clients.      //
//  Then the raw data file is deleted and it waits for a new file.           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMessage.h>
#include <TROOT.h>
#include <TServerSocket.h>
#include <TSocket.h>

#include "AliLog.h"
#include "AliESD.h"
#include "AliITSclustererV2.h"
#include "AliITSgeom.h"
#include "AliITStrackerV2.h"
#include "AliLoader.h"
#include "AliMonitorHLT.h"
#include "AliMonitorHLTHough.h"
#include "AliMonitorITS.h"
#include "AliMonitorProcess.h"
#include "AliMonitorTPC.h"
#include "AliMonitorV0s.h"
#include "AliRawReaderRoot.h"
#include "AliRun.h"
#include "AliTPC.h"
#include "AliTPCclustererMI.h"
#include "AliTPCtrackerMI.h"
#include "AliV0vertexer.h"

#include <AliHLTStandardIncludes.h>
#include <AliHLTMemHandler.h>
#include <AliHLTClusterFitter.h>
#include <AliHLTFitter.h>
#include <AliHLTHough.h>
#include <AliHLTHoughBaseTransformer.h>
#include <AliHLTStandardIncludes.h>
#include <AliHLTTrack.h>
#include <AliHLTTrackArray.h>
#include <AliHLTTransform.h>
#include <AliHLTVertex.h>
#include <AliLevel3.h>

ClassImp(AliMonitorProcess)


const Int_t AliMonitorProcess::fgkPort = 9327;


//_____________________________________________________________________________
AliMonitorProcess::AliMonitorProcess(
#if ROOT_VERSION_CODE <= 199169   // 3.10/01
				     const char* /*alienHost*/,
#else
				     const char* alienHost,
#endif
				     const char* alienDir,
				     const char* selection,
				     const char* fileNameGalice):
  fSelection(selection),
  fGrid(NULL),
  fAlienDir(alienDir),
  fRunLoader(NULL),
  fTPCParam(NULL),
  fITSgeom(NULL),
  fLogicalFileName(""),
  fFileName(""),
  fHLT(NULL),
  fHLTHough(NULL),

  fRunNumber(0),
  fSubRunNumber(0),
  fNEvents(0),
  fNEventsMin(1),
  fWriteHistoList(kFALSE),

  fTopFolder(NULL),
  fMonitors(),
  fFile(NULL),
  fTree(NULL),

  fServerSocket(NULL),
  fSockets(),
  fDisplaySocket(NULL),

  fStatus(kStopped),
  fStopping(kFALSE),

  fInterruptHandler(NULL)
{
// initialize the monitoring process and the monitor histograms

  fSelection = selection;

#if ROOT_VERSION_CODE <= 199169   // 3.10/01
  fGrid = TGrid::Connect("alien", gSystem->Getenv("USER"));
#else
  fGrid = TGrid::Connect(alienHost, gSystem->Getenv("USER"));
#endif
  if (!fGrid || fGrid->IsZombie() || !fGrid->IsConnected()) {
    delete fGrid;
    AliFatal("could not connect to alien");
  }
#if ROOT_VERSION_CODE <= 199169   // 3.10/01
  fGrid->cd(alienDir);
#endif

  fRunLoader = AliRunLoader::Open(fileNameGalice);
  if (!fRunLoader) AliFatal(Form("could not get run loader from file %s",
				 fileNameGalice));

  fRunLoader->CdGAFile();
  fTPCParam = AliTPC::LoadTPCParam(gFile);
  if (!fTPCParam) AliFatal("could not load TPC parameters");

  fRunLoader->LoadgAlice();
  gAlice = fRunLoader->GetAliRun();
  if (!gAlice) AliFatal("no gAlice object found");
  fITSgeom = (AliITSgeom*)gDirectory->Get("AliITSgeom");
  if (!fITSgeom) AliFatal("could not load ITS geometry");

  // Init TPC parameters for HLT
  Bool_t isinit=AliHLTTransform::Init(const_cast<char*>(fileNameGalice),kTRUE);
  if(!isinit){
    AliFatal("Could not create transform settings, please check log for error messages!");
  }

  fTopFolder = new TFolder("Monitor", "monitor histograms");
  fTopFolder->SetOwner(kTRUE);

  if (IsSelected("TPC")) fMonitors.Add(new AliMonitorTPC(fTPCParam));
  if (IsSelected("ITS")) fMonitors.Add(new AliMonitorITS(fITSgeom));
  if (IsSelected("V0s")) fMonitors.Add(new AliMonitorV0s);
  if (IsSelected("HLTConfMap")) fMonitors.Add(new AliMonitorHLT(fTPCParam));
  if (IsSelected("HLTHough")) fMonitors.Add(new AliMonitorHLTHough(fTPCParam));

  for (Int_t iMonitor = 0; iMonitor < fMonitors.GetEntriesFast(); iMonitor++) {
    ((AliMonitor*) fMonitors[iMonitor])->CreateHistos(fTopFolder);
  }

  TIterator* iFolder = fTopFolder->GetListOfFolders()->MakeIterator();
  while (TFolder* folder = (TFolder*) iFolder->Next()) folder->SetOwner(kTRUE);
  delete iFolder;

  fFile = TFile::Open("monitor_tree.root", "RECREATE");
  if (!fFile || !fFile->IsOpen()) {
    AliFatal("could not open file for tree");
  }
  fTree = new TTree("MonitorTree", "tree for monitoring");
  for (Int_t iMonitor = 0; iMonitor < fMonitors.GetEntriesFast(); iMonitor++) {
    ((AliMonitor*) fMonitors[iMonitor])->CreateBranches(fTree);
  }
  gROOT->cd();

  fServerSocket = new TServerSocket(fgkPort, kTRUE);
  fServerSocket->SetOption(kNoBlock, 1);
  CheckForConnections();

  fInterruptHandler = new AliMonitorInterruptHandler(this);
  gSystem->AddSignalHandler(fInterruptHandler);
}

//_____________________________________________________________________________
AliMonitorProcess::AliMonitorProcess(const AliMonitorProcess& process) :
  TObject(process),

  fSelection(""),
  fGrid(NULL),
  fAlienDir(""),
  fRunLoader(NULL),
  fTPCParam(NULL),
  fITSgeom(NULL),
  fLogicalFileName(""),
  fFileName(""),
  fHLT(NULL),
  fHLTHough(NULL),

  fRunNumber(0),
  fSubRunNumber(0),
  fNEvents(0),
  fNEventsMin(1),
  fWriteHistoList(kFALSE),

  fTopFolder(NULL),
  fMonitors(),
  fFile(NULL),
  fTree(NULL),

  fServerSocket(NULL),
  fSockets(),
  fDisplaySocket(NULL),

  fStatus(kStopped),
  fStopping(kFALSE),

  fInterruptHandler(NULL)

{
  AliFatal("copy constructor not implemented");
}

//_____________________________________________________________________________
AliMonitorProcess& AliMonitorProcess::operator = (const AliMonitorProcess&
						  /*process*/)
{
  AliFatal("assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliMonitorProcess::~AliMonitorProcess()
{
// clean up

  fMonitors.Delete();
  delete fTopFolder;
  delete fRunLoader;

  delete fServerSocket;
  fSockets.Delete();
  delete fDisplaySocket;

#if ROOT_VERSION_CODE <= 199169   // 3.10/01
  fGrid->Close();
#endif
  delete fGrid;

  fFile->Close();
  delete fFile;
  gSystem->Unlink("monitor_tree.root");

  delete fHLT;
  delete fHLTHough;

  gSystem->RemoveSignalHandler(fInterruptHandler);
  delete fInterruptHandler;
}


//_____________________________________________________________________________
const char* AliMonitorProcess::GetRevision()
{
  return "$Revision$";
}


//_____________________________________________________________________________
void AliMonitorProcess::SetStatus(EStatus status)
{
// set the current status and process system events

  fStatus = status;
  gSystem->ProcessEvents();
}


//_____________________________________________________________________________
void AliMonitorProcess::Run()
{
// run the monitor process:
//  check for a raw data file, process the raw data file and delete it

  fStopping = kFALSE;

  while (!fStopping) {
    SetStatus(kWaiting);
    while (!CheckForNewFile()) {
      CheckForConnections();
      SetStatus(kWaiting);
      if (fStopping) break;
      gSystem->Sleep(10);
    }
    if (fStopping) break;

    ProcessFile();
  }

  WriteHistos();

  fStopping = kFALSE;
  SetStatus(kStopped);
}


//_____________________________________________________________________________
void AliMonitorProcess::Stop()
{
// set the fStopping flag to terminate the monitor process after the current
// event was processed

  if (GetStatus() != kStopped) fStopping = kTRUE;
}


//_____________________________________________________________________________
void AliMonitorProcess::ProcessFile(const char* fileName)
{
// create a file with monitor histograms for a single file

  if (GetStatus() != kStopped) {
    AliError("ProcessFile can not be called"
	     " while the monitor process is running");
    return;
  }

  fFileName = fileName;
  Int_t nEventMin = fNEventsMin;
  fNEventsMin = 1;
  ProcessFile();
  WriteHistos();
  fNEventsMin = nEventMin;
  SetStatus(kStopped);
}


//_____________________________________________________________________________
Bool_t AliMonitorProcess::CheckForNewFile()
{
// check whether a new file was registered in alien

#if ROOT_VERSION_CODE < ROOT_VERSION(5,0,0)
#if ROOT_VERSION_CODE <= 199169   // 3.10/01
  TGridResult* result = fGrid->Ls();
#else
  TDatime datime;
  char dirName[256];
  sprintf(dirName, "%s/adc-%d", fAlienDir.Data(), datime.GetDate());
  char findName[256];
  sprintf(findName, "*.root");
  Grid_ResultHandle_t handle = fGrid->Find(dirName, findName);
  if (!handle) {
    AliError(Form("could not open alien directory %s",
		  dirName));
    return kFALSE;
  }
  TGridResult* result = fGrid->CreateGridResult(handle);
#endif
  Long_t maxDate = -1;
  Long_t maxTime = -1;
  TString fileName;

#if ROOT_VERSION_CODE <= 199169   // 3.10/01
  while (const char* entry = result->Next()) {
#else
  while (Grid_Result_t* resultEntry = result->Next()) {
    const char* entry = resultEntry->name.c_str();
#endif
    if (strrchr(entry, '/')) entry = strrchr(entry, '/')+1;
    // entry = host_date_time.root
    TString entryCopy(entry);
    char* p = const_cast<char*>(entryCopy.Data());
    if (!strtok(p, "_") || !p) continue;  // host name
    char* dateStr = strtok(NULL, "_");
    if (!dateStr || !p) continue;
    char* timeStr = strtok(NULL, ".");
    if (!timeStr || !p) continue;
    Long_t date = atoi(dateStr);
    Long_t time = atoi(timeStr);

    if ((date > maxDate) || ((date == maxDate) && (time > maxTime))) {
      maxDate = date;
      maxTime = time;
      fileName = entry;
    }
  }

  delete result;
  if (maxDate < 0) return kFALSE;  // no files found
  if (fLogicalFileName.CompareTo(fileName) == 0) return kFALSE;  // no new file

  fLogicalFileName = fileName;
#if ROOT_VERSION_CODE <= 199169   // 3.10/01
  result = fGrid->GetPhysicalFileNames(fLogicalFileName.Data());
  fFileName = result->Next();
#else
  fileName = dirName + ("/" + fLogicalFileName);
  handle = fGrid->GetPhysicalFileNames(fileName.Data());
  if (!handle) {
    AliError(Form("could not get physical file names for %s",
		  fileName.Data()));
    return kFALSE;
  }
  result = fGrid->CreateGridResult(handle);
  result->Reset();
  Grid_Result_t* resultEntry = result->Next();
  if (!resultEntry) {
    AliError(Form("could not get physical file names for %s",
		  fileName.Data()));
    return kFALSE;
  }
  fFileName = resultEntry->name2.c_str();
  fFileName.ReplaceAll("castor:/", "rfio:/");
#endif
  delete result;
#else
  Error("CheckForNewFile", "needs to be ported to new TGrid");
#endif
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMonitorProcess::ProcessFile()
{
// loop over all events in the raw data file, run the reconstruction
// and fill the monitor histograms

  Int_t nEvents = GetNumberOfEvents(fFileName);
  if (nEvents <= 0) return kFALSE;
  AliDebug(1, Form("found %d event(s) in file %s",
		   nEvents, fFileName.Data()));
  if (IsSelected("HLTConfMap")) CreateHLT(fFileName);
  if (IsSelected("HLTHough")) CreateHLTHough(fFileName);

  // loop over the events
  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    CheckForConnections();
    SetStatus(kReading);
    fRunLoader->SetEventNumber(0);
    AliRawReaderRoot rawReader(fFileName, iEvent);
    if (fStopping) break;
    if (rawReader.GetRunNumber() != fRunNumber) {
      WriteHistos();
      StartNewRun();
      fRunNumber = rawReader.GetRunNumber();
      fEventNumber[0] = rawReader.GetEventId()[0];
      fEventNumber[1] = rawReader.GetEventId()[1];
      fSubRunNumber = 0;
      if (fStopping) break;
    }

    // monitor only central physics events
    if (rawReader.GetType() != 7) continue;
    if ((rawReader.GetAttributes()[0] & 0x02) == 0) continue;
    AliInfo(Form("run: %d  event: %d %d\n", rawReader.GetRunNumber(),
		 rawReader.GetEventId()[0], rawReader.GetEventId()[1]));

    AliESD esd;
    if (IsSelected("TPC")) {
      CheckForConnections();
      if (!ReconstructTPC(&rawReader, &esd)) return kFALSE;
      if (fStopping) break;
    }
    if (IsSelected("ITS")) {
      CheckForConnections();
      if (!ReconstructITS(&rawReader, &esd)) return kFALSE;
      if (fStopping) break;
    }
    if (IsSelected("V0s")) {
      CheckForConnections();
      if (!ReconstructV0s(&esd)) return kFALSE;
      if (fStopping) break;
    }
    if (IsSelected("HLTConfMap")) {
      CheckForConnections();
      if (!ReconstructHLT(iEvent)) return kFALSE;
      if (fStopping) break;
    }
    if (IsSelected("HLTHough")) {
      CheckForConnections();
      if (!ReconstructHLTHough(iEvent)) return kFALSE;
      if (fStopping) break;
    }

    if (fDisplaySocket) fDisplaySocket->Send("new event");

    AliDebug(1, "filling histograms...");
    for (Int_t iMonitor = 0; iMonitor < fMonitors.GetEntriesFast(); iMonitor++) {
      CheckForConnections();
      SetStatus(kFilling);
      ((AliMonitor*) fMonitors[iMonitor])->FillHistos(fRunLoader, &rawReader,
						      &esd);
      if (fStopping) break;
    }
    if (fStopping) break;

    AliDebug(1, "updating histograms...");
    CheckForConnections();
    SetStatus(kUpdating);
    TIterator* iFolder = fTopFolder->GetListOfFolders()->MakeIterator();
    while (TFolder* folder = (TFolder*) iFolder->Next()) {
      TIterator* iHisto = folder->GetListOfFolders()->MakeIterator();
      while (AliMonitorPlot* histo = (AliMonitorPlot*) iHisto->Next()) {
	histo->Update();
      }
      delete iHisto;
    }
    delete iFolder;
    if (fStopping) break;

    AliDebug(1, "filling the tree...");
    fTree->Fill();

    AliDebug(1, "broadcasting histograms...");
    CheckForConnections();
    BroadcastHistos();

    fNEvents++;
    if (fStopping) break;
  }

  if (fHLT) {
    delete fHLT;
    fHLT = NULL;
  }
  if (fHLTHough) {
    delete fHLTHough;
    fHLTHough = NULL;
  }

  return kTRUE;
}

//_____________________________________________________________________________
void AliMonitorProcess::Reset()
{
// write the current histograms to a file and reset them

  if (fSubRunNumber == 0) fSubRunNumber++;
  WriteHistos();
  StartNewRun();
  fSubRunNumber++;
}


//_____________________________________________________________________________
UInt_t AliMonitorProcess::GetEventPeriodNumber() const
{
// get the period number from the event id

  return (fEventNumber[1] >> 4);
}

//_____________________________________________________________________________
UInt_t AliMonitorProcess::GetEventOrbitNumber() const
{
// get the orbit number from the event id

  return ((fEventNumber[1] & 0x000F) << 20) + (fEventNumber[0] >> 12);
}

//_____________________________________________________________________________
UInt_t AliMonitorProcess::GetEventBunchNumber() const
{
// get the bunch number from the event id

  return (fEventNumber[0] % 0x0FFF);
}

//_____________________________________________________________________________
Int_t AliMonitorProcess::GetNumberOfEvents(const char* fileName) const
{
// determine the number of events in the given raw data file

  Int_t nEvents = -1;

  TFile* file = TFile::Open(fileName);
  if (!file || !file->IsOpen()) {
    AliError(Form("could not open file %s", fileName));
    if (file) delete file;
    return -1;
  }

  TTree* tree = (TTree*) file->Get("RAW");
  if (!tree) {
    AliError("could not find tree with raw data");
  } else {
    nEvents = (Int_t) tree->GetEntries();
  }
  file->Close();
  delete file;

  return nEvents;
}

//_____________________________________________________________________________
Bool_t AliMonitorProcess::ReconstructTPC(AliRawReader* rawReader, AliESD* esd)
{
// find TPC clusters and tracks

  SetStatus(kRecTPC);

  AliLoader* tpcLoader = fRunLoader->GetLoader("TPCLoader");
  if (!tpcLoader) {
    AliError("no TPC loader found");
    return kFALSE;
  }
  gSystem->Unlink("TPC.RecPoints.root");

  // cluster finder
  AliDebug(1, "reconstructing clusters...");
  tpcLoader->LoadRecPoints("recreate");
  AliTPCclustererMI clusterer(fTPCParam);
  tpcLoader->MakeRecPointsContainer();
  clusterer.SetOutput(tpcLoader->TreeR());
  clusterer.Digits2Clusters(rawReader);
  tpcLoader->WriteRecPoints("OVERWRITE");

  // track finder
  AliDebug(1, "reconstructing tracks...");
  AliTPCtrackerMI tracker(fTPCParam);
  tracker.LoadClusters(tpcLoader->TreeR());
  tracker.Clusters2Tracks(esd);
  tracker.UnloadClusters();
  tpcLoader->UnloadRecPoints();

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMonitorProcess::ReconstructITS(AliRawReader* rawReader, AliESD* esd)
{
// find ITS clusters and tracks

  SetStatus(kRecITS);

  AliLoader* itsLoader = fRunLoader->GetLoader("ITSLoader");
  if (!itsLoader) {
    AliError("no ITS loader found");
    return kFALSE;
  }
  gSystem->Unlink("ITS.RecPoints.root");

  // cluster finder
  AliDebug(1, "reconstructing clusters...");
  itsLoader->LoadRecPoints("recreate");
  AliITSclustererV2 clusterer(fITSgeom);
  itsLoader->MakeRecPointsContainer();
  clusterer.Digits2Clusters(rawReader);

  // track finder
  AliDebug(1, "reconstructing tracks...");
  AliITStrackerV2 tracker(fITSgeom);
  tracker.LoadClusters(itsLoader->TreeR());
  tracker.Clusters2Tracks(esd);
  tracker.UnloadClusters();

  itsLoader->UnloadRecPoints();
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMonitorProcess::ReconstructV0s(AliESD* esd)
{
// find V0s

  SetStatus(kRecV0s);

  // V0 finder
  AliDebug(1, "reconstructing V0s...");
  AliV0vertexer vertexer;
  Double_t vtx[3];
  esd->GetVertex()->GetXYZ(vtx);
  vertexer.SetVertex(vtx);
  vertexer.Tracks2V0vertices(esd);

  return kTRUE;
}

//_____________________________________________________________________________
void AliMonitorProcess::CreateHLT(const char* fileName)
{

// create the HLT (Level3) object

  if (fHLT) delete fHLT;

  char name[256];
  strcpy(name, fileName);
  fHLT = new AliLevel3(name);
  fHLT->Init("./", AliLevel3::kRaw, 1);

  fHLT->SetClusterFinderParam(-1, -1, kTRUE);

  Int_t phiSegments = 50;
  Int_t etaSegments = 100;
  Int_t trackletlength = 3;
  Int_t tracklength = 20;//40 or 5
  Int_t rowscopetracklet = 2;
  Int_t rowscopetrack = 10;
  Double_t minPtFit = 0;
  Double_t maxangle = 0.1745;
  Double_t goodDist = 5;
  Double_t maxphi = 0.1;
  Double_t maxeta = 0.1;
  Double_t hitChi2Cut = 15;//100 or 15
  Double_t goodHitChi2 = 5;//20 or 5
  Double_t trackChi2Cut = 10;//50 or 10
  fHLT->SetTrackerParam(phiSegments, etaSegments,
			trackletlength, tracklength,
			rowscopetracklet, rowscopetrack,
			minPtFit, maxangle, goodDist, hitChi2Cut,
			goodHitChi2, trackChi2Cut, 50, maxphi, maxeta, kTRUE);

  fHLT->WriteFiles("./hlt/");
}

//_____________________________________________________________________________
void AliMonitorProcess::CreateHLTHough(const char* fileName)
{

// create the HLT Hough transform (L3Hough) object

  if (fHLTHough) delete fHLTHough;

  char name[256];
  strcpy(name, fileName);

  fHLTHough = new AliHLTHough();
  fHLTHough->SetThreshold(4);
  fHLTHough->SetTransformerParams(140,150,0.5,-1);
  fHLTHough->SetPeakThreshold(9000,-1);// or 6000
  fHLTHough->Init("./", kFALSE, 50, kFALSE,0,name);
  fHLTHough->SetAddHistograms();
  //  fHLTHough->GetMaxFinder()->SetThreshold(14000);

}

//_____________________________________________________________________________
Bool_t AliMonitorProcess::ReconstructHLT(Int_t iEvent)
{
// run the HLT cluster and track finder

  SetStatus(kRecHLT);

  gSystem->Exec("rm -rf hlt");
  gSystem->MakeDirectory("hlt");
  if (!fHLT) return kFALSE;

  fHLT->ProcessEvent(0, 35, iEvent);

  // remove the event number from the file names
  char command[256];
  sprintf(command, "rename points_%d points hlt/*.raw", iEvent);
  gSystem->Exec(command);
  sprintf(command, "rename tracks_tr_%d tracks_tr hlt/*.raw", iEvent);
  gSystem->Exec(command);
  sprintf(command, "rename tracks_gl_%d tracks_gl hlt/*.raw", iEvent);
  gSystem->Exec(command);
  sprintf(command, "rename tracks_%d tracks hlt/*.raw", iEvent);
  gSystem->Exec(command);
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMonitorProcess::ReconstructHLTHough(Int_t iEvent)
{
// run the HLT Hough transformer

  SetStatus(kRecHLT);

  gSystem->Exec("rm -rf hlt/hough");
  gSystem->MakeDirectory("hlt/hough");
  gSystem->Exec("rm -rf hlt/fitter");
  gSystem->MakeDirectory("hlt/fitter");
  if (!fHLTHough) return kFALSE;

  //  fHLTHough->Process(0, 35);
  // Loop over TPC sectors and process the data
  for(Int_t i=0; i<=35; i++)
    {
      fHLTHough->ReadData(i,iEvent);
      fHLTHough->Transform();
      //      if(fHLTHough->fAddHistograms)
      fHLTHough->AddAllHistograms();
      fHLTHough->FindTrackCandidates();
      fHLTHough->AddTracks();
    }
  fHLTHough->WriteTracks("./hlt/hough");

  // Run cluster fitter
  AliHLTClusterFitter *fitter = new AliHLTClusterFitter("./hlt");

  // Set debug flag for the cluster fitter
  //  fitter->Debug();

  // Setting fitter parameters
  fitter->SetInnerWidthFactor(1,1.5);
  fitter->SetOuterWidthFactor(1,1.5);
  fitter->SetNmaxOverlaps(5);

  //fitter->SetChiSqMax(5,kFALSE); //isolated clusters
  fitter->SetChiSqMax(5,kTRUE);  //overlapping clusters

  Int_t rowrange[2] = {0,AliHLTTransform::GetNRows()-1};

  // Takes input from global hough tracks produced by HT
  fitter->LoadSeeds(rowrange,kFALSE,iEvent);

  UInt_t ndigits;

  for(Int_t islice = 0; islice <= 35; islice++)
    {
      for(Int_t ipatch = 0; ipatch < AliHLTTransform::GetNPatches(); ipatch++)
	{
	  // Read digits
	  fHLTHough->GetMemHandler(ipatch)->Free();
	  fHLTHough->GetMemHandler(ipatch)->Init(islice,ipatch);
	  AliHLTDigitRowData *digits = (AliHLTDigitRowData *)fHLTHough->GetMemHandler(ipatch)->AliAltroDigits2Memory(ndigits,iEvent);

	  fitter->Init(islice,ipatch);
	  fitter->SetInputData(digits);
	  fitter->FindClusters();
	  fitter->WriteClusters();
 	}
    }

  // Refit of the clusters
  AliHLTVertex vertex;
  //The seeds are the input tracks from circle HT
  AliHLTTrackArray *tracks = fitter->GetSeeds();
  AliHLTFitter *ft = new AliHLTFitter(&vertex,1);

  ft->LoadClusters("./hlt/fitter/",iEvent,kFALSE);
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliHLTTrack *track = tracks->GetCheckedTrack(i);
      if(!track) continue;
      if(track->GetNHits() < 20) continue;
      ft->SortTrackClusters(track);
      ft->FitHelix(track);
      track->UpdateToFirstPoint();
    }
  delete ft;

  //Write the final tracks
  fitter->WriteTracks(20);

  delete fitter;

  // remove the event number from the file names
  char command[256];
  sprintf(command, "rename tracks_%d tracks hlt/hough/*.raw", iEvent);
  gSystem->Exec(command);
  sprintf(command, "rename tracks_%d tracks hlt/fitter/*.raw", iEvent);
  gSystem->Exec(command);
  sprintf(command, "rename points_%d points hlt/fitter/*.raw", iEvent);
  gSystem->Exec(command);
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMonitorProcess::WriteHistos()
{
// write the monitor tree and the monitor histograms to the file
// "monitor_<run number>[_<sub_run_number>].root"
// if at least fNEventsMin events were monitored

  SetStatus(kWriting);

  // rename tree file and create a new one
  fFile->cd();
  fTree->Write();
  fFile->Close();
  delete fFile;

  char fileName[256];
  sprintf(fileName, "monitor_tree_%d.root", fRunNumber);
  if (fSubRunNumber > 0) {
    sprintf(fileName, "monitor_tree_%d_%d.root", fRunNumber, fSubRunNumber);
  }
  if (fNEvents < fNEventsMin) {
    gSystem->Unlink("monitor_tree.root");
  } else {
    gSystem->Rename("monitor_tree.root", fileName);
  }

  fFile = TFile::Open("monitor_tree.root", "RECREATE");
  if (!fFile || !fFile->IsOpen()) {
    AliFatal("could not open file for tree");
  }
  fTree = new TTree("MonitorTree", "tree for monitoring");
  for (Int_t iMonitor = 0; iMonitor < fMonitors.GetEntriesFast(); iMonitor++) {
    ((AliMonitor*) fMonitors[iMonitor])->CreateBranches(fTree);
  }
  gROOT->cd();

  // write the histograms
  if (fNEvents < fNEventsMin) return kTRUE;

  if (!fWriteHistoList) {
    TIterator* iFolder = fTopFolder->GetListOfFolders()->MakeIterator();
    while (TFolder* folder = (TFolder*) iFolder->Next()) {
      TIterator* iHisto = folder->GetListOfFolders()->MakeIterator();
      while (AliMonitorPlot* histo = (AliMonitorPlot*) iHisto->Next()) {
	histo->ResetList();
      }
      delete iHisto;
    }
    delete iFolder;
  }

  Bool_t result = kTRUE;
  sprintf(fileName, "monitor_%d.root", fRunNumber);
  if (fSubRunNumber > 0) {
    sprintf(fileName, "monitor_%d_%d.root", fRunNumber, fSubRunNumber);
  }
  TFile* file = TFile::Open(fileName, "recreate");
  if (!file || !file->IsOpen()) {
    AliError(Form("could not open file %s", fileName));
    result = kFALSE;
  } else {
    fTopFolder->Write();
    file->Close();
  }
  if (file) delete file;

  return result;
}

//_____________________________________________________________________________
void AliMonitorProcess::StartNewRun()
{
// reset the histograms for a new run

  SetStatus(kResetting);
  TIterator* iFolder = fTopFolder->GetListOfFolders()->MakeIterator();
  while (TFolder* folder = (TFolder*) iFolder->Next()) {
    TIterator* iHisto = folder->GetListOfFolders()->MakeIterator();
    while (AliMonitorPlot* histo = (AliMonitorPlot*) iHisto->Next()) {
      histo->Reset();
    }
    delete iHisto;
  }
  delete iFolder;

  fNEvents = 0;
}


//_____________________________________________________________________________
void AliMonitorProcess::CheckForConnections()
{
// check if new clients want to connect and add them to the list of sockets

  TSocket* socket;
  while ((socket = fServerSocket->Accept()) != (TSocket*)-1) {
    socket->SetOption(kNoBlock, 1);
    char socketType[256];
    if (socket->Recv(socketType, 255) <= 0) {
      gSystem->Sleep(1000);
      if (socket->Recv(socketType, 255) <= 0) {
	TInetAddress adr = socket->GetInetAddress();
	AliError(Form("no socket type received - "
		      "disconnect client: %s (%s), port %d",
		      adr.GetHostName(), adr.GetHostAddress(), adr.GetPort()));
	delete socket;
	continue;
      }
    }
    if (strcmp(socketType, "client") == 0) {
      fSockets.Add(socket);
      TInetAddress adr = socket->GetInetAddress();
      AliInfo(Form("new client: %s (%s), port %d",
		   adr.GetHostName(), adr.GetHostAddress(), adr.GetPort()));
      if (fNEvents > 0) BroadcastHistos(socket);
    } else if (strcmp(socketType, "display") == 0) {
      if (fDisplaySocket) {
	fDisplaySocket->Close();
	delete fDisplaySocket;
      }
      fDisplaySocket = socket;
      TInetAddress adr = socket->GetInetAddress();
      AliInfo(Form("new display: %s (%s), port %d",
		   adr.GetHostName(), adr.GetHostAddress(), adr.GetPort()));
    } else {
      TInetAddress adr = socket->GetInetAddress();
      AliError(Form("unknown socket type %s - "
		    "disconnect client: %s (%s), port %d", socketType,
		    adr.GetHostName(), adr.GetHostAddress(), adr.GetPort()));
      delete socket;
      continue;
    }
  }

  // remove finished or invalid clients
  for (Int_t iSocket = 0; iSocket < fSockets.GetEntriesFast(); iSocket++) {
    socket = (TSocket*) fSockets[iSocket];
    if (!socket) continue;
    char controlMessage[256];
    if (socket->Recv(controlMessage, 255)) {
      if (strcmp(controlMessage, "disconnect") == 0) {
	TInetAddress adr = socket->GetInetAddress();
	AliInfo(Form("disconnect client: %s (%s), port %d",
		     adr.GetHostName(), adr.GetHostAddress(), adr.GetPort()));
	delete fSockets.RemoveAt(iSocket);
	continue;
      }
    }
    if (!socket->IsValid()) {
      // remove invalid sockets from the list
      TInetAddress adr = socket->GetInetAddress();
      AliError(Form("disconnect invalid client: %s (%s), port %d",
		    adr.GetHostName(), adr.GetHostAddress(), adr.GetPort()));
      delete fSockets.RemoveAt(iSocket);
    }
  }
  fSockets.Compress();
}

//_____________________________________________________________________________
void AliMonitorProcess::BroadcastHistos(TSocket* toSocket)
{
// send the monitor histograms to the clients

  SetStatus(kBroadcasting);
  TMessage message(kMESS_OBJECT);
  message.WriteObject(fTopFolder);

  for (Int_t iSocket = 0; iSocket < fSockets.GetEntriesFast(); iSocket++) {
    TSocket* socket = (TSocket*) fSockets[iSocket];
    if (!socket) continue;
    if (toSocket && (socket != toSocket)) continue;

    // send control message
    if (!socket->IsValid() || (socket->Send("histograms") <= 0)) {
      TInetAddress adr = socket->GetInetAddress();
      AliError(Form("connection to client failed - "
		    "disconnect client: %s (%s), port %d",
		    adr.GetHostName(), adr.GetHostAddress(), adr.GetPort()));
      delete fSockets.RemoveAt(iSocket);
    }

    // receive control message
    char controlMessage[256];
    Int_t result = socket->Recv(controlMessage, 255);
    if (result <= 0) {
      gSystem->Sleep(1000);  // wait one second and try again
      result = socket->Recv(controlMessage, 255);
    }
    if (result <= 0) {
      TInetAddress adr = socket->GetInetAddress();
      AliError(Form("no response from client - "
		    "disconnect client: %s (%s), port %d",
		    adr.GetHostName(), adr.GetHostAddress(), adr.GetPort()));
      delete fSockets.RemoveAt(iSocket);
      continue;
    }
    if (strcmp(controlMessage, "ok") != 0) {
      TInetAddress adr = socket->GetInetAddress();
      AliError(Form("no \"ok\" message from client - "
		    "disconnect client: %s (%s), port %d",
		    adr.GetHostName(), adr.GetHostAddress(), adr.GetPort()));
      delete fSockets.RemoveAt(iSocket);
      continue;
    }

    socket->SetOption(kNoBlock, 0);
    if (socket->Send(message) < 0) {
      // remove the socket from the list if there was an error
      TInetAddress adr = socket->GetInetAddress();
      AliError(Form("sending histograms failed - "
		    "disconnect client: %s (%s), port %d",
		    adr.GetHostName(), adr.GetHostAddress(), adr.GetPort()));
      delete fSockets.RemoveAt(iSocket);
    } else {
      gSystem->Sleep(100);
      socket->SetOption(kNoBlock, 1);
    }
  }
  fSockets.Compress();
}


//_____________________________________________________________________________
AliMonitorProcess::AliMonitorInterruptHandler::AliMonitorInterruptHandler
  (AliMonitorProcess* process):
  TSignalHandler(kSigUser1, kFALSE),
  fProcess(process)
{
// constructor: set process
}

//_____________________________________________________________________________
Bool_t AliMonitorProcess::AliMonitorInterruptHandler::Notify()
{
// interrupt signal -> stop process

  AliInfo("the monitoring process will be stopped.");
  fProcess->Stop();
  return kTRUE;
}
