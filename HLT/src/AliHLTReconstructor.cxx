// $Id$

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for HLT reconstruction                                              //
// <Cvetan.Cheshkov@cern.ch>                                                 //
// <loizides@ikf.uni-frankfurt.de>                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TSystem.h>
#include <TArrayF.h>
#include <TObjString.h>

#include <AliESDEvent.h>
#include <AliESDHLTtrack.h>

#include "AliHLTStandardIncludes.h"
#include "AliHLTLogging.h"
//#include "AliLevel3.h"
//#include "AliHLTEvaluate.h"
#include "AliHLTReconstructor.h"
#include "AliHLTTransform.h"
#include "AliHLTHough.h"
#include "AliHLTFileHandler.h"
#include "AliHLTTrack.h"
#include "AliHLTHoughTrack.h"
#include "AliHLTTrackArray.h"

#include "AliLog.h"
#include "AliITS.h"
#include "AliHLTITStracker.h"
#include "AliHLTTPCtracker.h"
//#include "MUON/src/AliRoot/AliHLTMUONTracker.h"
//#include "MUON/src/AliRoot/AliHLTMUONHitReconstructor.h"
#include "AliRawReader.h"
#include "AliHLTSystem.h"

#if __GNUC__== 3
using namespace std;
#endif

/** HLT default component libraries */
const char* kHLTDefaultLibs[]= {
  "libAliHLTUtil.so", 
  "libAliHLTTPC.so", 
  //  "libAliHLTSample.so",
  "libAliHLTPHOS.so",
  NULL
};

ClassImp(AliHLTReconstructor)

AliHLTReconstructor::AliHLTReconstructor()
  : 
  AliReconstructor(),
  fDoHough(0),
  fDoTracker(1),
  fDoBench(0),
  fDoCleanUp(0),
  fpSystem(NULL)
{ 
  //constructor
#ifndef use_logging
  AliHLTLog::fgLevel=AliHLTLog::kWarning;
#endif
}

AliHLTReconstructor::AliHLTReconstructor(Bool_t doTracker, Bool_t doHough)
  : 
  AliReconstructor(),
  fDoHough(doHough),
  fDoTracker(doTracker),
  fDoBench(0),
  fDoCleanUp(0),
  fpSystem(new AliHLTSystem)
{
  //constructor
#ifndef use_logging
  AliHLTLog::fgLevel=AliHLTLog::kWarning;
#endif
}

AliHLTReconstructor::~AliHLTReconstructor()
{ 
  //destructor
  if(fDoCleanUp){
    char name[256];
    gSystem->Exec("rm -rf hlt");
    sprintf(name, "rm -f confmap_*.root confmap_*.dat");
    gSystem->Exec(name);
    gSystem->Exec("rm -rf hough");
    sprintf(name, "rm -f hough_*.root hough_*.dat");
    gSystem->Exec(name);
  }
  if (fpSystem) {
    delete fpSystem;
  }
  fpSystem=NULL;
}

void AliHLTReconstructor::Init()
{
  // init the reconstructor
  if (!fpSystem) fpSystem=new AliHLTSystem;
  if (!fpSystem) {
    AliError("can not create AliHLTSystem object");
    return;
  }
  if (fpSystem->CheckStatus(AliHLTSystem::kError)) {
    AliError("HLT system in error state");
    return;
  }

  TString libs("");
  TString option = GetOption();
  TObjArray* pTokens=option.Tokenize(" ");
  option="";
  if (pTokens) {
    int iEntries=pTokens->GetEntries();
    for (int i=0; i<iEntries; i++) {
      TString token=(((TObjString*)pTokens->At(i))->GetString());
      if (token.Contains("loglevel=")) {
	TString param=token.ReplaceAll("loglevel=", "");
	if (param.IsDigit()) {
	  fpSystem->SetGlobalLoggingLevel((AliHLTComponentLogSeverity)param.Atoi());
	} else if (param.BeginsWith("0x") &&
		   param.Replace(0,2,"",0).IsHex()) {
	  int severity=0;
	  sscanf(param.Data(),"%x", &severity);
	  fpSystem->SetGlobalLoggingLevel((AliHLTComponentLogSeverity)severity);
	} else {
	  AliWarning("wrong parameter for option \'loglevel=\', (hex) number expected");
	}
      } else if (token.Contains("alilog=off")) {
	fpSystem->SwitchAliLog(0);
      } else if (token.BeginsWith("lib") && token.EndsWith(".so")) {
	libs+=token;
	libs+=" ";
      } else {
	if (option.Length()>0) option+=" ";
	option+=token;
      }
    }
    delete pTokens;
  }
  
  Bool_t bForceLibLoad=0;
  if (bForceLibLoad=(libs.IsNull())) {
    const char** deflib=kHLTDefaultLibs;
    while (*deflib) {
      libs+=*deflib++;
      libs+=" ";
    }
  }
  if ((bForceLibLoad || !fpSystem->CheckStatus(AliHLTSystem::kLibrariesLoaded)) &&
      (fpSystem->LoadComponentLibraries(libs.Data())<0)) {
    AliError("error while loading HLT libraries");
    return;
  }

  if (!fpSystem->CheckStatus(AliHLTSystem::kReady)) {
    typedef int (*AliHLTSystemSetOptions)(AliHLTSystem* pInstance, const char* options);
    gSystem->Load("libHLTinterface.so");
    AliHLTSystemSetOptions pFunc=(AliHLTSystemSetOptions)(gSystem->DynFindSymbol("libHLTinterface.so", "AliHLTSystemSetOptions"));
    if (pFunc) {
      if ((pFunc)(fpSystem, option.Data())<0) {
      AliError("error setting options for HLT system");
      return;	
      }
    } else if (option.Length()>0) {
      AliError(Form("version of HLT system does not support the options \'%s\'", option.Data()));
      return;
    }
    if ((fpSystem->Configure())<0) {
      AliError("error during HLT system configuration");
      return;
    }
  }
}

// void AliHLTReconstructor::Reconstruct(AliRunLoader* runLoader, AliRawReader* rawReader) const 
// {
//   // reconstruction of real data if rawReader!=NULL
//   if(!runLoader) {
//     AliError("Missing RunLoader! 0x0");
//     return;
//   }

//   Int_t nEvents = runLoader->GetNumberOfEvents();
//   int iResult=0;

//   if (fpSystem) {
//     if (fpSystem->CheckStatus(AliHLTSystem::kError)) {
//       AliError("HLT system in error state");
//       return;
//     }
//     if ((iResult=fpSystem->Reconstruct(nEvents, runLoader, rawReader))>=0) {
//     }
//   }
// }

// void AliHLTReconstructor::ReconstructWithConformalMapping(AliRunLoader* runLoader,Int_t iEvent) const
// {
//   // reconstruct with conformal mapper
//   /*
//   AliLevel3 *fHLT = new AliLevel3(runLoader);
//   fHLT->Init("./", AliLevel3::kRunLoader, 1);

//   Int_t phiSegments = 50;
//   Int_t etaSegments = 100;
//   Int_t trackletlength = 3;
//   Int_t tracklength = 10;
//   Int_t rowscopetracklet = 2;
//   Int_t rowscopetrack = 10;
//   Double_t minPtFit = 0;
//   Double_t maxangle = 0.1745;
//   Double_t goodDist = 5;
//   Double_t maxphi = 0.1;
//   Double_t maxeta = 0.1;
//   Double_t hitChi2Cut = 20;
//   Double_t goodHitChi2 = 5;
//   Double_t trackChi2Cut = 10;
//   Double_t xyerror = -1;
//   Double_t zerror =  -1;
  
//   fHLT->SetClusterFinderParam(xyerror,zerror,kTRUE);
//   fHLT->SetTrackerParam(phiSegments, etaSegments, 
// 			trackletlength, tracklength,
// 			rowscopetracklet, rowscopetrack,
// 			minPtFit, maxangle, goodDist, hitChi2Cut,
// 			goodHitChi2, trackChi2Cut, 50, maxphi, maxeta, kTRUE);
//   fHLT->SetTrackerParam(phiSegments, etaSegments, 
// 			trackletlength, tracklength,
// 			rowscopetracklet, rowscopetrack,
// 			minPtFit, maxangle, goodDist, hitChi2Cut,
// 			goodHitChi2, trackChi2Cut, 50, maxphi, maxeta, kFALSE);
//   fHLT->SetMergerParameters(2,3,0.003,0.1,0.05);
//   fHLT->DoMc();
//   fHLT->DoNonVertexTracking(); // 2 tracking passes, last without vertex contraint.
//   fHLT->WriteFiles("./hlt/");  
//   fHLT->ProcessEvent(0, 35, iEvent);
//   if(fDoBench){
//     char filename[256];
//     sprintf(filename, "confmap_%d",iEvent);
//     fHLT->DoBench(filename);
//   }

//   delete fHLT;
//   */
// }

// void AliHLTReconstructor::ReconstructWithHoughTransform(AliRunLoader* runLoader,Int_t iEvent) const
// {
//   //reconstruct with hough
//   //not used anymore, Hough tracking is moved out of the local
//   //reconstruction chain
//   Float_t ptmin = 0.1*AliHLTTransform::GetSolenoidField();

//   Float_t zvertex = 0;
//   TArrayF mcVertex(3); 
//   AliHeader * header = runLoader->GetHeader();
//   if (header) {
//     AliGenEventHeader * genHeader = header->GenEventHeader();
//     if (genHeader) genHeader->PrimaryVertex(mcVertex);
//   }
//   zvertex = mcVertex[2];

//   AliInfo(Form("Hough Transform will run with ptmin=%f and zvertex=%f", ptmin, zvertex));

//   AliHLTHough *hough = new AliHLTHough();
    
//   hough->SetThreshold(4);
//   hough->CalcTransformerParams(ptmin);
//   hough->SetPeakThreshold(70,-1);
//   hough->SetRunLoader(runLoader);
//   hough->Init("./", kFALSE, 100, kFALSE,4,0,0,zvertex);
//   hough->SetAddHistograms();

//   for(Int_t slice=0; slice<=35; slice++)
//     {
//       //cout<<"Processing slice "<<slice<<endl;
//       hough->ReadData(slice,iEvent);
//       hough->Transform();
//       hough->AddAllHistogramsRows();
//       hough->FindTrackCandidatesRow();
//       //hough->WriteTracks(slice,"./hough");
//       hough->AddTracks();
//     }
//   hough->WriteTracks("./hough");
  
//   if(fDoBench){
//     char filename[256];
//     sprintf(filename, "hough_%d",iEvent);
//     hough->DoBench(filename);
//   }
//   delete hough;
// }

// void AliHLTReconstructor::FillESD(AliRunLoader* runLoader, 
// 				  AliESDEvent* esd) const
// {
//   //fill the esd file with found tracks
//   if(!runLoader) {
//     AliError("Missing RunLoader! 0x0");
//     return;
//   }
//   Int_t iEvent = runLoader->GetEventNumber();
//   if (fpSystem) {
//     if (fpSystem->CheckStatus(AliHLTSystem::kError)) {
//       AliError("HLT system in error state");
//       return;
//     }
//     if (!fpSystem->CheckStatus(AliHLTSystem::kReady)) {
//       AliError("HLT system in wrong state");
//       return;
//     }
//     fpSystem->FillESD(iEvent, runLoader, esd);
//   }
//   /*
//   if(fDoTracker) FillESDforConformalMapping(esd,iEvent);
//   if(fDoHough) FillESDforHoughTransform(esd,iEvent);
//   */
// }

// void AliHLTReconstructor::FillESDforConformalMapping(AliESDEvent* esd,Int_t iEvent) const
// {
//   //fill esd with tracks from conformal mapping
//   /*
//   Int_t slicerange[2]={0,35};
//   Int_t good = (int)(0.4*AliHLTTransform::GetNRows());
//   Int_t nclusters = (int)(0.4*AliHLTTransform::GetNRows());
//   Int_t nminpointsontracks = (int)(0.3*AliHLTTransform::GetNRows());
//   Float_t ptmin = 0.;
//   Float_t ptmax = 0.;
//   Float_t maxfalseratio = 0.1;
  
//   AliHLTEvaluate *fHLTEval = new AliHLTEvaluate("./hlt",nclusters,good,ptmin,ptmax,slicerange);
//   fHLTEval->SetMaxFalseClusters(maxfalseratio);
//   fHLTEval->LoadData(iEvent,kTRUE);
//   fHLTEval->AssignPIDs();
//   fHLTEval->AssignIDs();
//   AliHLTTrackArray *fTracks = fHLTEval->GetTracks();
//   if(!fTracks){
//     delete fHLTEval;
//     return;
//   }
//   for(Int_t i=0; i<fTracks->GetNTracks(); i++)
//     {
//       AliHLTTrack *tpt = (AliHLTTrack *)fTracks->GetCheckedTrack(i);
//       if(!tpt) continue; 
//       if(tpt->GetNumberOfPoints() < nminpointsontracks) continue;
      
//       AliESDHLTtrack *esdtrack = new AliESDHLTtrack() ; 

//       esdtrack->SetRowRange(tpt->GetFirstRow(),tpt->GetLastRow());
//       esdtrack->SetNHits(tpt->GetNHits());
//       esdtrack->SetFirstPoint(tpt->GetFirstPointX(),tpt->GetFirstPointY(),tpt->GetFirstPointZ());
//       esdtrack->SetLastPoint(tpt->GetLastPointX(),tpt->GetLastPointY(),tpt->GetLastPointZ());
//       esdtrack->SetPt(tpt->GetPt());
//       esdtrack->SetPsi(tpt->GetPsi());
//       esdtrack->SetTgl(tpt->GetTgl());
//       esdtrack->SetCharge(tpt->GetCharge());
//       esdtrack->SetMCid(tpt->GetMCid());
//       esdtrack->SetSector(tpt->GetSector());
//       esdtrack->SetPID(tpt->GetPID());
//       esdtrack->ComesFromMainVertex(tpt->ComesFromMainVertex());

//       esd->AddHLTConfMapTrack(esdtrack);
//       delete esdtrack;
//     }
//   delete fHLTEval;
//   */
// }

// void AliHLTReconstructor::FillESDforHoughTransform(AliESDEvent* esd,Int_t iEvent) const
// {
//   //fill esd with tracks from hough
//   char filename[256];
//   sprintf(filename,"./hough/tracks_%d.raw",iEvent);
  
//   AliHLTFileHandler *tfile = new AliHLTFileHandler();
//   if(!tfile->SetBinaryInput(filename)){
//     AliError(Form("Missing file %s", filename));
//     return;
//   }
  
//   AliHLTTrackArray *fTracks = new AliHLTTrackArray("AliHLTHoughTrack");
//   tfile->Binary2TrackArray(fTracks);
//   tfile->CloseBinaryInput();
//   delete tfile;
//   if(!fTracks) return; 
//   for(Int_t i=0; i<fTracks->GetNTracks(); i++)
//     {
//       AliHLTHoughTrack *tpt = (AliHLTHoughTrack *)fTracks->GetCheckedTrack(i);
//       if(!tpt) continue; 
      
//       AliESDHLTtrack *esdtrack = new AliESDHLTtrack() ; 

//       esdtrack->SetRowRange(tpt->GetFirstRow(),tpt->GetLastRow());
//       esdtrack->SetNHits(tpt->GetNHits());
//       esdtrack->SetFirstPoint(tpt->GetFirstPointX(),tpt->GetFirstPointY(),tpt->GetFirstPointZ());
//       esdtrack->SetLastPoint(tpt->GetLastPointX(),tpt->GetLastPointY(),tpt->GetLastPointZ());
//       esdtrack->SetPt(tpt->GetPt());
//       esdtrack->SetPsi(tpt->GetPsi());
//       esdtrack->SetTgl(tpt->GetTgl());
//       esdtrack->SetCharge(tpt->GetCharge());
//       esdtrack->SetMCid(tpt->GetMCid());
//       esdtrack->SetWeight(tpt->GetWeight());
//       esdtrack->SetSector(tpt->GetSector());
//       esdtrack->SetBinXY(tpt->GetBinX(),tpt->GetBinY(),tpt->GetSizeX(),tpt->GetSizeY());
//       esdtrack->SetPID(tpt->GetPID());
//       esdtrack->ComesFromMainVertex(tpt->ComesFromMainVertex());

//       esd->AddHLTHoughTrack(esdtrack);
//       delete esdtrack;
//     }

//   delete fTracks;
// }

/* The following functions are deprecated and need to be removed.
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
    return new AliHLTITStracker(0);
  }
  if(!opt.CompareTo("MUON")) {
    return new AliHLTMUONTracker(runLoader);
  }

  return NULL;
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

*/
