// $Id$

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for HLT reconstruction                                              //
// <Cvetan.Cheshkov@cern.ch>                                                 //
///////////////////////////////////////////////////////////////////////////////

// very ugly but it has to work fast
#ifdef use_reconstruction

#include <Riostream.h>
#include <TSystem.h>
#include <TArrayF.h>

#include "AliL3StandardIncludes.h"
#include "AliL3Logging.h"
#include "AliLevel3.h"
#include "AliL3Evaluate.h"
#include "AliHLTReconstructor.h"
#include "AliL3Transform.h"
#include "AliL3Hough.h"
#include "AliL3FileHandler.h"
#include "AliL3Track.h"
#include "AliL3HoughTrack.h"
#include "AliL3TrackArray.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliESD.h"
#include "AliESDHLTtrack.h"

#if __GNUC__== 3
using namespace std;
#endif


ClassImp(AliHLTReconstructor)

void AliHLTReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
  if(!runLoader) {
    LOG(AliL3Log::kFatal,"AliHLTReconstructor::Reconstruct","RunLoader")
      <<" Missing RunLoader! 0x0"<<ENDLOG;
    return;
  }
  gSystem->Exec("rm -rf hlt");
  gSystem->MakeDirectory("hlt");
  gSystem->Exec("rm -rf hough");
  gSystem->MakeDirectory("hough");

  Bool_t isinit=AliL3Transform::Init(runLoader);
  if(!isinit){
    LOG(AliL3Log::kError,"AliHLTReconstructor::Reconstruct","Transformer")
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

  char filename[256];
  sprintf(filename, "confmap_%d",iEvent);
  fHLT->DoBench(filename);

  delete fHLT;
}

void AliHLTReconstructor::ReconstructWithHoughTransform(AliRunLoader* runLoader,Int_t iEvent) const
{
  Float_t ptmin = 0.1*AliL3Transform::GetSolenoidField();

  Float_t zvertex = 0;
  TArrayF mcVertex(3); 
  runLoader->GetHeader()->GenEventHeader()->PrimaryVertex(mcVertex);
  zvertex = mcVertex[2];

  LOG(AliL3Log::kInformational,"AliHLTReconstructor::Reconstruct","HoughTransform")
    <<" Hough Transform will run with ptmin="<<ptmin<<" and zvertex="<<zvertex<<ENDLOG;

  AliL3Hough *hough = new AliL3Hough();
    
  hough->SetThreshold(4);
  hough->SetTransformerParams(76,140,ptmin,-1);
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
  
  char filename[256];
  sprintf(filename, "hough_%d",iEvent);
  hough->DoBench(filename);

  delete hough;
}

void AliHLTReconstructor::FillESD(AliRunLoader* runLoader, 
				  AliESD* esd) const
{
  Int_t iEvent = runLoader->GetEventNumber();

  if(fDoTracker) FillESDforConformalMapping(esd,iEvent);
  if(fDoHough) FillESDforHoughTransform(esd,iEvent);

#if 0
  char name[256];
  gSystem->Exec("rm -rf hlt");
  sprintf(name, "rm -f confmap_%d.root confmap_%d.dat",iEvent,iEvent);
  gSystem->Exec(name);
  gSystem->Exec("rm -rf hough");
  sprintf(name, "rm -f hough_%d.root hough_%d.dat",iEvent,iEvent);
  gSystem->Exec(name);
#endif
}

void AliHLTReconstructor::FillESDforConformalMapping(AliESD* esd,Int_t iEvent) const
{
  //Assign MC labels for found tracks
  Int_t slicerange[2]={0,35};
  Int_t good = (int)(0.4*AliL3Transform::GetNRows());
  Int_t nclusters = (int)(0.4*AliL3Transform::GetNRows());
  Float_t ptmin = 0.;
  Float_t ptmax = 0.;
  Float_t maxfalseratio = 0.1;
  
  AliL3Evaluate *fHLTEval = new AliL3Evaluate("./hlt",nclusters,good,ptmin,ptmax,slicerange);
  fHLTEval->SetMaxFalseClusters(maxfalseratio);

  fHLTEval->LoadData(iEvent,kTRUE);
  fHLTEval->AssignPIDs();
  fHLTEval->AssignIDs();
  AliL3TrackArray *fTracks = fHLTEval->GetTracks();
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliL3Track *tpt = (AliL3Track *)fTracks->GetCheckedTrack(i);
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
  char filename[256];
  sprintf(filename,"./hough/tracks_%d.raw",iEvent);
  
  AliL3FileHandler *tfile = new AliL3FileHandler();
  if(!tfile->SetBinaryInput(filename)){
    Error("FillESD","Inputfile ",filename," does not exist");
    return;
  }
  
  AliL3TrackArray *fTracks = new AliL3TrackArray("AliL3HoughTrack");
  tfile->Binary2TrackArray(fTracks);
  tfile->CloseBinaryInput();
  delete tfile;
  
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliL3HoughTrack *tpt = (AliL3HoughTrack *)fTracks->GetCheckedTrack(i);
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

#endif
