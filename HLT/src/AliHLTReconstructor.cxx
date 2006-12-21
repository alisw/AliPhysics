// $Id$

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for HLT reconstruction                                              //
// <Cvetan.Cheshkov@cern.ch>                                                 //
// <loizides@ikf.uni-frankfurt.de>                                           //
///////////////////////////////////////////////////////////////////////////////

// very ugly but it has to work fast
#ifdef use_reconstruction

#include <Riostream.h>
#include <TSystem.h>
#include <TArrayF.h>

#include <AliRunLoader.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliESD.h>
#include <AliESDHLTtrack.h>

#include "AliHLTStandardIncludes.h"
#include "AliHLTLogging.h"
#include "AliLevel3.h"
#include "AliHLTEvaluate.h"
#include "AliHLTReconstructor.h"
#include "AliHLTTransform.h"
#include "AliHLTHough.h"
#include "AliHLTFileHandler.h"
#include "AliHLTTrack.h"
#include "AliHLTHoughTrack.h"
#include "AliHLTTrackArray.h"

#include "AliRun.h"
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliHLTITStracker.h"
#include "AliHLTTPCtracker.h"
#include "MUON/src/AliRoot/AliHLTMUONTracker.h"
#include "MUON/src/AliRoot/AliHLTMUONHitReconstructor.h"
#include "AliRawReader.h"
#if __GNUC__== 3
using namespace std;
#endif

ClassImp(AliHLTReconstructor)

AliHLTReconstructor::AliHLTReconstructor(): AliReconstructor() 
{ 
  //constructor
#ifndef use_logging
  AliHLTLog::fgLevel=AliHLTLog::kWarning;
#endif
  fDoTracker=1;
  fDoHough=0;
  fDoBench=0;
  fDoCleanUp=1;
}

AliHLTReconstructor::AliHLTReconstructor(Bool_t doTracker, Bool_t doHough): AliReconstructor() 
{ 
  //constructor
#ifndef use_logging
  AliHLTLog::fgLevel=AliHLTLog::kWarning;
#endif
  fDoTracker=doTracker;
  fDoHough=doHough;
  fDoBench=0;
  fDoCleanUp=1;
}

AliHLTReconstructor::~AliHLTReconstructor()
{ 
  //deconstructor
  if(fDoCleanUp){
    char name[256];
    gSystem->Exec("rm -rf hlt");
    sprintf(name, "rm -f confmap_*.root confmap_*.dat");
    gSystem->Exec(name);
    gSystem->Exec("rm -rf hough");
    sprintf(name, "rm -f hough_*.root hough_*.dat");
    gSystem->Exec(name);
  }
}

void AliHLTReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
  // do the standard and hough reconstruction chain
  if(!runLoader) {
    LOG(AliHLTLog::kFatal,"AliHLTReconstructor::Reconstruct","RunLoader")
      <<" Missing RunLoader! 0x0"<<ENDLOG;
    return;
  }
  gSystem->Exec("rm -rf hlt");
  gSystem->MakeDirectory("hlt");
  gSystem->Exec("rm -rf hough");
  gSystem->MakeDirectory("hough");

  Bool_t isinit=AliHLTTransform::Init(runLoader);
  if(!isinit){
    LOG(AliHLTLog::kError,"AliHLTReconstructor::Reconstruct","Transformer")
     << "Could not create transform settings, please check log for error messages!" << ENDLOG;
    return;
  }

  Int_t nEvents = runLoader->GetNumberOfEvents();

  for(Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    runLoader->GetEvent(iEvent);

    if(fDoTracker) ReconstructWithConformalMapping(runLoader,iEvent);
    if(fDoHough) ReconstructWithHoughTransform(runLoader,iEvent);
  }
}

void AliHLTReconstructor::ReconstructWithConformalMapping(AliRunLoader* runLoader,Int_t iEvent) const
{
  // reconstruct with conformal mapper
  AliLevel3 *fHLT = new AliLevel3(runLoader);
  fHLT->Init("./", AliLevel3::kRunLoader, 1);

  Int_t phiSegments = 50;
  Int_t etaSegments = 100;
  Int_t trackletlength = 3;
  Int_t tracklength = 10;
  Int_t rowscopetracklet = 2;
  Int_t rowscopetrack = 10;
  Double_t minPtFit = 0;
  Double_t maxangle = 0.1745;
  Double_t goodDist = 5;
  Double_t maxphi = 0.1;
  Double_t maxeta = 0.1;
  Double_t hitChi2Cut = 20;
  Double_t goodHitChi2 = 5;
  Double_t trackChi2Cut = 10;
  Double_t xyerror = -1;
  Double_t zerror =  -1;
  
  fHLT->SetClusterFinderParam(xyerror,zerror,kTRUE);
  fHLT->SetTrackerParam(phiSegments, etaSegments, 
			trackletlength, tracklength,
			rowscopetracklet, rowscopetrack,
			minPtFit, maxangle, goodDist, hitChi2Cut,
			goodHitChi2, trackChi2Cut, 50, maxphi, maxeta, kTRUE);
  fHLT->SetTrackerParam(phiSegments, etaSegments, 
			trackletlength, tracklength,
			rowscopetracklet, rowscopetrack,
			minPtFit, maxangle, goodDist, hitChi2Cut,
			goodHitChi2, trackChi2Cut, 50, maxphi, maxeta, kFALSE);
  fHLT->SetMergerParameters(2,3,0.003,0.1,0.05);
  fHLT->DoMc();
  fHLT->DoNonVertexTracking(); /*2 tracking passes, last without vertex contraint.*/
  fHLT->WriteFiles("./hlt/");  
  fHLT->ProcessEvent(0, 35, iEvent);
  if(fDoBench){
    char filename[256];
    sprintf(filename, "confmap_%d",iEvent);
    fHLT->DoBench(filename);
  }

  delete fHLT;
}

void AliHLTReconstructor::ReconstructWithHoughTransform(AliRunLoader* runLoader,Int_t iEvent) const
{
  //reconstruct with hough
  //not used anymore, Hough tracking is moved out of the local
  //reconstruction chain
  Float_t ptmin = 0.1*AliHLTTransform::GetSolenoidField();

  Float_t zvertex = 0;
  TArrayF mcVertex(3); 
  AliHeader * header = runLoader->GetHeader();
  if (header) {
    AliGenEventHeader * genHeader = header->GenEventHeader();
    if (genHeader) genHeader->PrimaryVertex(mcVertex);
  }
  zvertex = mcVertex[2];

  LOG(AliHLTLog::kInformational,"AliHLTReconstructor::Reconstruct","HoughTransform")
    <<" Hough Transform will run with ptmin="<<ptmin<<" and zvertex="<<zvertex<<ENDLOG;

  AliHLTHough *hough = new AliHLTHough();
    
  hough->SetThreshold(4);
  hough->CalcTransformerParams(ptmin);
  hough->SetPeakThreshold(70,-1);
  hough->SetRunLoader(runLoader);
  hough->Init("./", kFALSE, 100, kFALSE,4,0,0,zvertex);
  hough->SetAddHistograms();

  for(Int_t slice=0; slice<=35; slice++)
    {
      //cout<<"Processing slice "<<slice<<endl;
      hough->ReadData(slice,iEvent);
      hough->Transform();
      hough->AddAllHistogramsRows();
      hough->FindTrackCandidatesRow();
      //hough->WriteTracks(slice,"./hough");
      hough->AddTracks();
    }
  hough->WriteTracks("./hough");
  
  if(fDoBench){
    char filename[256];
    sprintf(filename, "hough_%d",iEvent);
    hough->DoBench(filename);
  }
  delete hough;
}

void AliHLTReconstructor::FillESD(AliRunLoader* runLoader, 
				  AliESD* esd) const
{
  //fill the esd file with found tracks
  Int_t iEvent = runLoader->GetEventNumber();

  if(fDoTracker) FillESDforConformalMapping(esd,iEvent);
  if(fDoHough) FillESDforHoughTransform(esd,iEvent);
}

void AliHLTReconstructor::FillESDforConformalMapping(AliESD* esd,Int_t iEvent) const
{
  //fill esd with tracks from conformal mapping
  Int_t slicerange[2]={0,35};
  Int_t good = (int)(0.4*AliHLTTransform::GetNRows());
  Int_t nclusters = (int)(0.4*AliHLTTransform::GetNRows());
  Int_t nminpointsontracks = (int)(0.3*AliHLTTransform::GetNRows());
  Float_t ptmin = 0.;
  Float_t ptmax = 0.;
  Float_t maxfalseratio = 0.1;
  
  AliHLTEvaluate *fHLTEval = new AliHLTEvaluate("./hlt",nclusters,good,ptmin,ptmax,slicerange);
  fHLTEval->SetMaxFalseClusters(maxfalseratio);
  fHLTEval->LoadData(iEvent,kTRUE);
  fHLTEval->AssignPIDs();
  fHLTEval->AssignIDs();
  AliHLTTrackArray *fTracks = fHLTEval->GetTracks();
  if(!fTracks){
    delete fHLTEval;
    return;
  }
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliHLTTrack *tpt = (AliHLTTrack *)fTracks->GetCheckedTrack(i);
      if(!tpt) continue; 
      if(tpt->GetNumberOfPoints() < nminpointsontracks) continue;
      
      AliESDHLTtrack *esdtrack = new AliESDHLTtrack() ; 

      esdtrack->SetRowRange(tpt->GetFirstRow(),tpt->GetLastRow());
      esdtrack->SetNHits(tpt->GetNHits());
      esdtrack->SetFirstPoint(tpt->GetFirstPointX(),tpt->GetFirstPointY(),tpt->GetFirstPointZ());
      esdtrack->SetLastPoint(tpt->GetLastPointX(),tpt->GetLastPointY(),tpt->GetLastPointZ());
      esdtrack->SetPt(tpt->GetPt());
      esdtrack->SetPsi(tpt->GetPsi());
      esdtrack->SetTgl(tpt->GetTgl());
      esdtrack->SetCharge(tpt->GetCharge());
      esdtrack->SetMCid(tpt->GetMCid());
      esdtrack->SetSector(tpt->GetSector());
      esdtrack->SetPID(tpt->GetPID());
      esdtrack->ComesFromMainVertex(tpt->ComesFromMainVertex());

      esd->AddHLTConfMapTrack(esdtrack);
      delete esdtrack;
    }
  delete fHLTEval;
}

void AliHLTReconstructor::FillESDforHoughTransform(AliESD* esd,Int_t iEvent) const
{
  //fill esd with tracks from hough
  char filename[256];
  sprintf(filename,"./hough/tracks_%d.raw",iEvent);
  
  AliHLTFileHandler *tfile = new AliHLTFileHandler();
  if(!tfile->SetBinaryInput(filename)){
    LOG(AliHLTLog::kError,"AliHLTReconstructor::FillESDforHoughTransform","Input file")
      <<" Missing file "<<filename<<ENDLOG;
    return;
  }
  
  AliHLTTrackArray *fTracks = new AliHLTTrackArray("AliHLTHoughTrack");
  tfile->Binary2TrackArray(fTracks);
  tfile->CloseBinaryInput();
  delete tfile;
  if(!fTracks) return; 
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliHLTHoughTrack *tpt = (AliHLTHoughTrack *)fTracks->GetCheckedTrack(i);
      if(!tpt) continue; 
      
      AliESDHLTtrack *esdtrack = new AliESDHLTtrack() ; 

      esdtrack->SetRowRange(tpt->GetFirstRow(),tpt->GetLastRow());
      esdtrack->SetNHits(tpt->GetNHits());
      esdtrack->SetFirstPoint(tpt->GetFirstPointX(),tpt->GetFirstPointY(),tpt->GetFirstPointZ());
      esdtrack->SetLastPoint(tpt->GetLastPointX(),tpt->GetLastPointY(),tpt->GetLastPointZ());
      esdtrack->SetPt(tpt->GetPt());
      esdtrack->SetPsi(tpt->GetPsi());
      esdtrack->SetTgl(tpt->GetTgl());
      esdtrack->SetCharge(tpt->GetCharge());
      esdtrack->SetMCid(tpt->GetMCid());
      esdtrack->SetWeight(tpt->GetWeight());
      esdtrack->SetSector(tpt->GetSector());
      esdtrack->SetBinXY(tpt->GetBinX(),tpt->GetBinY(),tpt->GetSizeX(),tpt->GetSizeY());
      esdtrack->SetPID(tpt->GetPID());
      esdtrack->ComesFromMainVertex(tpt->ComesFromMainVertex());

      esd->AddHLTHoughTrack(esdtrack);
      delete esdtrack;
    }

  delete fTracks;
}

AliTracker* AliHLTReconstructor::CreateTracker(AliRunLoader* runLoader) const
{
  //Create HLT trackers for TPC and ITS

  TString opt = GetOption();
  if(!opt.CompareTo("TPC")) {
    // Create Hough tracker for TPC
    return new AliHLTTPCtracker(runLoader);
  }
  if(!opt.CompareTo("ITS")) {
    // Create ITS tracker
    AliITSgeom* geom = GetITSgeom(runLoader);
    if (!geom) return NULL;
    return new AliHLTITStracker(geom);
  }
  if(!opt.CompareTo("MUON")) {
    return new AliHLTMUONTracker(runLoader);
  }

  return NULL;
}

//_____________________________________________________________________________
AliITSgeom* AliHLTReconstructor::GetITSgeom(AliRunLoader* runLoader) const
{
// get the ITS geometry

  if (!runLoader->GetAliRun()) runLoader->LoadgAlice();
  if (!runLoader->GetAliRun()) {
    Error("GetITSgeom", "couldn't get AliRun object");
    return NULL;
  }
  AliITS* its = (AliITS*) runLoader->GetAliRun()->GetDetector("ITS");
  if (!its) {
    Error("GetITSgeom", "couldn't get ITS detector");
    return NULL;
  }
  if (!its->GetITSgeom()) {
    Error("GetITSgeom", "no ITS geometry available");
    return NULL;
  }
  return its->GetITSgeom();
}

void AliHLTReconstructor::FillDHLTRecPoint(AliRawReader* rawReader, Int_t nofEvent, Int_t dcCut = 0) const
{
  // Hit recontruction for dimuon-HLT
  AliHLTMUONHitReconstructor kdHLTRec(rawReader);

  Int_t iEvent = 0 ;
  kdHLTRec.SetDCCut(dcCut);
  TString lookupTablePath = getenv("ALICE_ROOT");
  lookupTablePath += "/HLT/MUON/src/AliRoot/Lut";
  kdHLTRec.Init(lookupTablePath.Data(),lookupTablePath.Data());

   while(rawReader->NextEvent() && iEvent < nofEvent){
     AliInfo(Form("Event : %d",iEvent));
     kdHLTRec.WriteDHLTRecHits(iEvent);  
     iEvent++;
   }

}

#endif
