/* $Id: AliPhiCorrelationsQATask.cxx 47416 2011-02-15 13:55:09Z jgrosseo $ */

#include "AliPhiCorrelationsQATask.h"

#include <TH2F.h>
#include <TParticle.h>
#include <TFile.h>

#include <AliLog.h>
#include <AliESDVertex.h>
#include <AliESDEvent.h>
#include <AliMCEvent.h>
#include <AliStack.h>
#include <AliHeader.h>
#include <AliAnalysisManager.h>
#include <AliCentrality.h>

#include "AliESDtrackCuts.h"

ClassImp(AliPhiCorrelationsQATask)

AliPhiCorrelationsQATask::AliPhiCorrelationsQATask(const char* opt) :
  AliAnalysisTaskSE("AliPhiCorrelationsQATask"),
  fOutput(0),
  fOption(opt),
  fEsdTrackCuts(0),
  fEsdTrackCuts2(0),
  fCheckITS(0),
  fGlobalTracks(0),
  fCentralityCorrelation(0),
  fDCAPrimaries(0),
  fDCASecondaries(0),
  fUseUncheckedCentrality(kFALSE)
{
  //
  // Constructor. Initialization of pointers
  //

  // Define input and output slots here
  DefineOutput(1, TList::Class());
  
  AliLog::SetClassDebugLevel("AliPhiCorrelationsQATask", AliLog::kWarning);
}

AliPhiCorrelationsQATask::~AliPhiCorrelationsQATask()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
    fOutput = 0;
  }
}

void AliPhiCorrelationsQATask::UserCreateOutputObjects()
{
  // create result objects and add to output list

  Printf("AliPhiCorrelationsQATask::CreateOutputObjects");
  //AliLog::SetClassDebugLevel("AliPhysicsSelection", AliLog::kDebug);

  fOutput = new TList;
  fOutput->SetOwner();

  fCentralityCorrelation = new TH2F("fCentralityCorrelation", ";v0 centr;spd centr", 100, 0, 100.001, 100, 0, 100.001);
  fOutput->Add(fCentralityCorrelation);

  fEsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fEsdTrackCuts->SetMinNClustersTPC(70);
  // disable DCA cut
  fEsdTrackCuts->SetMaxDCAToVertexZ();
  fEsdTrackCuts->SetMaxDCAToVertexXY();
  fEsdTrackCuts->SetDCAToVertex2D(kFALSE);
  fEsdTrackCuts->SetName("cuts_quality_only");
  fEsdTrackCuts->DefineHistograms();
  fOutput->Add(fEsdTrackCuts);
  
  fEsdTrackCuts2 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fEsdTrackCuts2->SetMinNClustersTPC(70);
  fEsdTrackCuts2->SetName("cuts_quality_dca");
  fEsdTrackCuts2->DefineHistograms();
  fOutput->Add(fEsdTrackCuts2);
  
  fGlobalTracks = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
  fGlobalTracks->SetName("global_cuts");
  fGlobalTracks->DefineHistograms();
  fOutput->Add(fGlobalTracks);
        
  fCheckITS = new AliESDtrackCuts;
  fCheckITS->SetRequireITSRefit(kTRUE);
  fCheckITS->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                      AliESDtrackCuts::kAny);
  fCheckITS->SetDCAToVertex2D(kFALSE);
  fCheckITS->SetRequireSigmaToVertex(kFALSE);
  fCheckITS->SetName("check_its");
  fCheckITS->DefineHistograms();
  fOutput->Add(fCheckITS);
  
  fDCAPrimaries = new TH2F("fDCAPrimaries", ";dca_xy;dca_z", 1000, -5, 5, 1000, -5, 5);
  fDCASecondaries = (TH2F*) fDCAPrimaries->Clone("fDCASecondaries");
  
  fOutput->Add(fDCAPrimaries);
  fOutput->Add(fDCASecondaries);

  PostData(1, fOutput);
}

void AliPhiCorrelationsQATask::UserExec(Option_t*)
{
  // process the event
  
  //Printf("AliPhiCorrelationsQATask::UserExec");

  AliESDEvent* esd = (AliESDEvent*) fInputEvent;
  if (!esd)
    AliFatal("No input event");
    
  // vertex cut
  const AliESDVertex* vtxESD = esd->GetPrimaryVertex();
  const AliESDVertex *vtxSPD = esd->GetPrimaryVertexSPD();
  if (!vtxESD)
    return;
  if (!vtxSPD)
    return;
    
  Double_t vtx[3];
  vtxESD->GetXYZ(vtx);
  if (TMath::Abs(vtx[2]) > 10)
    return;

  AliCentrality *centralityObj = esd->GetCentrality();
  if (centralityObj)
  {
    Float_t v0Centrality = -1;
    Float_t spdCentrality = -1; 
      
    if (fUseUncheckedCentrality)
    {
      v0Centrality = centralityObj->GetCentralityPercentileUnchecked("V0M");
      spdCentrality = centralityObj->GetCentralityPercentileUnchecked("CL1");
    }
    else
    {
      v0Centrality = centralityObj->GetCentralityPercentile("V0M");
      spdCentrality = centralityObj->GetCentralityPercentile("CL1");
    }
    
    fCentralityCorrelation->Fill(v0Centrality, spdCentrality);
  }
      
  AliStack* stack = 0;
  if (fMCEvent)
  {
    stack = fMCEvent->Stack();
    if (!stack)
      AliFatal("Stack is 0");
  }
  
  Int_t count = 0;
    
  for (Int_t i=0; i<esd->GetNumberOfTracks(); i++)
  {
    AliESDtrack* esdTrack = esd->GetTrack(i); 
    
    if (0)
    {
      if (!fGlobalTracks->AcceptTrack(esdTrack))
        continue;
    }
    
    if (!fEsdTrackCuts->AcceptTrack(esdTrack))
      continue;
      
    count++;
    
    //if (v0Centrality < 0)
    //  continue;
      
    Float_t b[2];
    Float_t bCov[3];
    esdTrack->GetImpactParameters(b, bCov);
    
    Bool_t primary = kTRUE;
    if (stack)
      primary = stack->IsPhysicalPrimary(TMath::Abs(esdTrack->GetLabel()));
    
    if (primary)
      fDCAPrimaries->Fill(b[0], b[1]);
    else
      fDCASecondaries->Fill(b[0], b[1]);
      
    // fill histograms of second object
    if (!fEsdTrackCuts2->AcceptTrack(esdTrack))
      continue;
      
    if (fCheckITS->AcceptTrack(esdTrack))
    {
      if (!fGlobalTracks->AcceptTrack(esdTrack))
        continue;
    }
    
    /*
    // create a tpc only track
    AliESDtrack *track = AliESDtrackCuts::GetTPCOnlyTrack(esd,esdTrack->GetID());
    if (!track) 
      continue;

    if (track->Pt()>0.)
    {
      // only constrain tracks above threshold
      AliExternalTrackParam exParam;
      // take the B-feild from the ESD, no 3D fieldMap available at this point
      Bool_t relate = false;
      relate = track->RelateToVertex(vtxSPD,esd->GetMagneticField(),kVeryBig,&exParam);
      if(!relate){
        delete track;
        continue;
      }
      track->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),exParam.GetCovariance());
    }
    
    Float_t b2[2];
    track->GetImpactParameters(b2,bCov);
    
    //Printf("%.2f %.2f %.2f %.2f", b[0], b2[0], b[1], b2[1]);
    
    delete track;
    */
  }      
  
  //Printf("%.2f: %d %d tracks (out of %d)", v0Centrality, centralityObj->GetQuality(), count, esd->GetNumberOfTracks());
}

void AliPhiCorrelationsQATask::Terminate(Option_t *)
{
  // terminate
  
  Printf("In AliPhiCorrelationsQATask::Terminate...");
  
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput)
  {
    Printf("ERROR: fOutput not available");
    return;
  }

  fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("cuts_quality_only"));
  if (!fEsdTrackCuts)
  {
    Printf("ERROR: fEsdTrackCuts not available");
    return;
  }
  
  fEsdTrackCuts2 = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("cuts_quality_dca"));
  if (!fEsdTrackCuts2)
  {
    Printf("ERROR: fEsdTrackCuts2 not available");
    return;
  }
  
  fGlobalTracks = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("global_cuts"));
  if (!fGlobalTracks)
  {
    Printf("ERROR: fGlobalTracks not available");
    return;
  }
  
  TFile* file = TFile::Open("track_cuts.root", "RECREATE");
  fEsdTrackCuts->SaveHistograms();
  fEsdTrackCuts2->SaveHistograms();
  fGlobalTracks->SaveHistograms();
  file->Close();
}
