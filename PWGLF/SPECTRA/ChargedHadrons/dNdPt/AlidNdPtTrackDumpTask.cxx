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

#include "iostream"

#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include "TChain.h"
#include "TTreeStream.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TRandom3.h"

#include "AliHeader.h"  
#include "AliGenEventHeader.h"  
#include "AliInputEventHandler.h"  
#include "AliStack.h"  
#include "AliTrackReference.h"  

#include "AliPhysicsSelection.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliTracker.h"
#include "AliGeomManager.h"

#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliMultiplicity.h"

#include "AliESDtrackCuts.h"
#include "AliMCEventHandler.h"
#include "AlidNdPt.h"
#include "AlidNdPtEventCuts.h"
#include "AlidNdPtAcceptanceCuts.h"

#include "AlidNdPtTrackDumpTask.h"
#include "AliKFParticle.h"
#include "AliESDv0.h"

using namespace std;

ClassImp(AlidNdPtTrackDumpTask)

//_____________________________________________________________________________
AlidNdPtTrackDumpTask::AlidNdPtTrackDumpTask(const char *name) 
  : AliAnalysisTaskSE(name)
  , fESD(0)
  , fMC(0)
  , fESDfriend(0)
  , fOutput(0)
  , fPitList(0)
  , fUseMCInfo(kFALSE)
  , fUseESDfriends(kFALSE)
  , fdNdPtEventCuts(0)
  , fdNdPtAcceptanceCuts(0)
  , fdNdPtRecAcceptanceCuts(0)
  , fEsdTrackCuts(0)
  , fTrigger(AliTriggerAnalysis::kMB1) 
  , fAnalysisMode(AlidNdPtHelper::kTPC) 
  , fTreeSRedirector(0)
  , fCentralityEstimator(0)
  , fLowPtTrackDownscaligF(0)
  , fLowPtV0DownscaligF(0)
  , fProcessAll(kFALSE)
{
  // Constructor

  // Define input and output slots here
  DefineOutput(1, TList::Class());

}

//_____________________________________________________________________________
AlidNdPtTrackDumpTask::~AlidNdPtTrackDumpTask()
{
  if(fOutput) delete fOutput;  fOutput =0; 
  if(fTreeSRedirector) delete fTreeSRedirector;  fTreeSRedirector =0; 

  if(fdNdPtEventCuts) delete fdNdPtEventCuts; fdNdPtEventCuts=NULL; 
  if(fdNdPtAcceptanceCuts) delete fdNdPtAcceptanceCuts; fdNdPtAcceptanceCuts=NULL;
  if(fdNdPtRecAcceptanceCuts) delete fdNdPtRecAcceptanceCuts; fdNdPtRecAcceptanceCuts=NULL;  
  if(fEsdTrackCuts) delete fEsdTrackCuts; fEsdTrackCuts=NULL;
}

//____________________________________________________________________________
Bool_t AlidNdPtTrackDumpTask::Notify()
{
  static Int_t count = 0;
  count++;
  TTree *chain = (TChain*)GetInputData(0);
  if(chain)
  {
    Printf("Processing %d. file: %s", count, chain->GetCurrentFile()->GetName());
  }

return kTRUE;
}

//_____________________________________________________________________________
void AlidNdPtTrackDumpTask::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  fOutput = new TList;
  fOutput->SetOwner();

  //
  // create temporary file for output tree
  fTreeSRedirector = new TTreeSRedirector("jotwinow_Temp_Trees.root");

  PostData(1, fOutput);
}

//_____________________________________________________________________________
void AlidNdPtTrackDumpTask::UserExec(Option_t *) 
{
  //
  // Called for each event
  //

  // ESD event
  fESD = (AliESDEvent*) (InputEvent());
  if (!fESD) {
    Printf("ERROR: ESD event not available");
    return;
  }

  // MC event
  if(fUseMCInfo) {
    fMC = MCEvent();
    if (!fMC) {
      Printf("ERROR: MC event not available");
      return;
    }
  }

  if(fUseESDfriends) {
    fESDfriend = static_cast<AliESDfriend*>(fESD->FindListObject("AliESDfriend"));
      if(!fESDfriend) {
        Printf("ERROR: ESD friends not available");
    }
  }

  //
  if(fProcessAll) { 
    ProcessAll(fESD,fMC,fESDfriend); // all track stages and MC
    ProcessV0(fESD,fMC,fESDfriend);
    ProcessLaser(fESD,fMC,fESDfriend);
    ProcessdEdx(fESD,fMC,fESDfriend);
    if(IsUseMCInfo())
      ProcessMCEff(fESD,fMC,fESDfriend);
  }
  else {
    Process(fESD,fMC,fESDfriend);
    ProcessV0(fESD,fMC,fESDfriend);
    ProcessLaser(fESD,fMC,fESDfriend);
    ProcessdEdx(fESD,fMC,fESDfriend);
    if(IsUseMCInfo())
      ProcessMCEff(fESD,fMC,fESDfriend);
  }

  // Post output data.
  PostData(1, fOutput);
}

//_____________________________________________________________________________
void AlidNdPtTrackDumpTask::Process(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent, AliESDfriend *const /*esdFriend*/)
{
  //
  // Select real events with high-pT tracks 
  //
  if(!esdEvent) {
    AliDebug(AliLog::kError, "esdEvent not available");
    return;
  }

  // get selection cuts
  AlidNdPtEventCuts *evtCuts = GetEventCuts(); 
  AlidNdPtAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts  || !esdTrackCuts) {
    Printf("ERROR cuts not available");
    return;
  }

  // trigger selection
  Bool_t isEventTriggered = kTRUE;
  AliPhysicsSelection *physicsSelection = NULL;
  AliTriggerAnalysis* triggerAnalysis = NULL;

  // 
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return;
  }
   
  // get file name
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { 
    Printf("ERROR: Could not receive input chain");
    return;
  }
  TObjString fileName(chain->GetCurrentFile()->GetName());

  // trigger
  if(evtCuts->IsTriggerRequired())  
  {
    // always MB
    isEventTriggered = inputHandler->IsEventSelected() & AliVEvent::kMB;

    physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if(!physicsSelection) return;
    //SetPhysicsTriggerSelection(physicsSelection);

    if (isEventTriggered && (GetTrigger() == AliTriggerAnalysis::kV0AND)) {
      // set trigger (V0AND)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
      if(!triggerAnalysis) return;
      isEventTriggered = triggerAnalysis->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }
  }

  // centrality determination
  Float_t centralityF = -1;
  AliCentrality *esdCentrality = esdEvent->GetCentrality();
  centralityF = esdCentrality->GetCentralityPercentile(fCentralityEstimator.Data());

  // use MC information
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);

  Int_t multMCTrueTracks = 0;
  if(IsUseMCInfo())
  {
    //
    if(!mcEvent) {
      AliDebug(AliLog::kError, "mcEvent not available");
      return;
    }
    // get MC event header
    header = mcEvent->Header();
    if (!header) {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }
    // MC particle stack
    stack = mcEvent->Stack();
    if (!stack) {
      AliDebug(AliLog::kError, "Stack not available");
      return;
    }

    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);

    // multipliticy of all MC primary tracks
    // in Zv, pt and eta ranges)
    multMCTrueTracks = AlidNdPtHelper::GetMCTrueTrackMult(mcEvent,evtCuts,accCuts);

  } // end bUseMC

  // get reconstructed vertex  
  //const AliESDVertex* vtxESD = 0; 
  const AliESDVertex* vtxESD = 0; 
  if(GetAnalysisMode() == AlidNdPtHelper::kTPC) {
        vtxESD = esdEvent->GetPrimaryVertexTPC();
  }
  else if(GetAnalysisMode() == AlidNdPtHelper::kTPCITS) {
     vtxESD = esdEvent->GetPrimaryVertexTracks();
  }
  else {
    	return;
  }

  if(!vtxESD) return;

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD); 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());


  // check event cuts
  if(isEventOK && isEventTriggered)
  {

    //
    // get IR information
    //
    AliESDHeader *esdHeader = 0;
    esdHeader = esdEvent->GetHeader();
    if(!esdHeader) return;
    //Int_t ir1 = esdHeader->GetTriggerIREntries(); //all ir-s
    //Int_t ir2 = esdHeader->GetTriggerIREntries(-1,1); // int2 set, 180 ms time interval

    // Use when Peter commit the changes in the header
    Int_t ir1 = 0;
    Int_t ir2 = 0;

    //
    Double_t vert[3] = {0}; 
    vert[0] = vtxESD->GetXv();
    vert[1] = vtxESD->GetYv();
    vert[2] = vtxESD->GetZv();
    Int_t mult = vtxESD->GetNContributors();
    Double_t bz = esdEvent->GetMagneticField();
    Double_t runNumber = esdEvent->GetRunNumber();
    Double_t evtTimeStamp = esdEvent->GetTimeStamp();
    Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();

    // high pT tracks
    for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++)
    {
      AliESDtrack *track = esdEvent->GetTrack(iTrack);
      if(!track) continue;
      if(track->Charge()==0) continue;
      if(!esdTrackCuts->AcceptTrack(track)) continue;
      if(!accCuts->AcceptTrack(track)) continue;
      
      // downscale low-pT tracks
      Double_t scalempt= TMath::Min(track->Pt(),10.);
      Double_t downscaleF = gRandom->Rndm();
      downscaleF *= fLowPtTrackDownscaligF;
      if(TMath::Exp(2*scalempt)<downscaleF) continue;
      //printf("TMath::Exp(2*scalempt) %e, downscaleF %e \n",TMath::Exp(2*scalempt), downscaleF);

      AliExternalTrackParam * tpcInner = (AliExternalTrackParam *)(track->GetTPCInnerParam());
      if (!tpcInner) continue;
      // transform to the track reference frame 
      Bool_t isOK = kFALSE;
      isOK = tpcInner->Rotate(track->GetAlpha());
      isOK = tpcInner->PropagateTo(track->GetX(),esdEvent->GetMagneticField());
      if(!isOK) continue;

      // Dump to the tree 
      // vertex
      // TPC-ITS tracks
      //
      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"dNdPtTree"<<
        "fileName.="<<&fileName<<
        "runNumber="<<runNumber<<
        "evtTimeStamp="<<evtTimeStamp<<
        "evtNumberInFile="<<evtNumberInFile<<
        "Bz="<<bz<<
	"vertX="<<vert[0]<<
	"vertY="<<vert[1]<<
	"vertZ="<<vert[2]<<
	"IRtot="<<ir1<<
	"IRint2="<<ir2<<
        "mult="<<mult<<
        "esdTrack.="<<track<<
        "centralityF="<<centralityF<<
        "\n";
    }
  }
  
  PostData(1, fOutput);
}


//_____________________________________________________________________________
void AlidNdPtTrackDumpTask::ProcessLaser(AliESDEvent *const esdEvent, AliMCEvent * const /*mcEvent*/, AliESDfriend *const /*esdFriend*/)
{
  //
  // Process laser events
  //
  if(!esdEvent) {
    AliDebug(AliLog::kError, "esdEvent not available");
    return;
  }

  // get file name
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { 
    Printf("ERROR: Could not receive input chain");
    return;
  }
  TObjString fileName(chain->GetCurrentFile()->GetName());

  // laser events 
  const AliESDHeader* esdHeader = esdEvent->GetHeader();
  if(esdHeader && esdHeader->GetEventSpecie()==AliRecoParam::kCalib) 
  {
    Int_t countLaserTracks = 0;
    for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++)
    {
      AliESDtrack *track = esdEvent->GetTrack(iTrack);
      if(!track) continue;

      if(track->GetTPCInnerParam()) countLaserTracks++;
    }
       
    if(countLaserTracks > 100) 
    {      
      Double_t runNumber = esdEvent->GetRunNumber();
      Double_t evtTimeStamp = esdEvent->GetTimeStamp();
      Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();
      Double_t bz = esdEvent->GetMagneticField();

      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"Laser"<<
        "fileName.="<<&fileName<<
        "runNumber="<<runNumber<<
        "evtTimeStamp="<<evtTimeStamp<<
        "evtNumberInFile="<<evtNumberInFile<<
        "Bz="<<bz<<
        "multTPCtracks="<<countLaserTracks<<
        "\n";
    }
  }
}

//_____________________________________________________________________________
void AlidNdPtTrackDumpTask::ProcessAll(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent, AliESDfriend *const esdFriend)
{
  //
  // Select real events with high-pT tracks
  // Calculate and stor in the output tree:
  //  TPC constrained tracks
  //  InnerParams constrained tracks
  //  TPC-ITS tracks
  //  ITSout-InnerParams tracks
  //  chi2 distance between TPC constrained and TPC-ITS tracks
  //  chi2 distance between TPC refitted constrained and TPC-ITS tracks
  //  chi2 distance between ITSout and InnerParams tracks
  //  MC information: 
  //   track references at ITSin, TPCin; InnerParam at first TPC track reference, 
  //   particle ID, mother ID, production mechanism ...
  // 

  if(!esdEvent) {
    AliDebug(AliLog::kError, "esdEvent not available");
    return;
  }

  // get selection cuts
  AlidNdPtEventCuts *evtCuts = GetEventCuts(); 
  AlidNdPtAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts  || !esdTrackCuts) {
    AliDebug(AliLog::kError, "cuts not available");
    return;
  }

  // trigger selection
  Bool_t isEventTriggered = kTRUE;
  AliPhysicsSelection *physicsSelection = NULL;
  AliTriggerAnalysis* triggerAnalysis = NULL;

  // 
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return;
  }
   
  // get file name
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { 
    Printf("ERROR: Could not receive input chain");
    return;
  }
  TObjString fileName(chain->GetCurrentFile()->GetName());

  // trigger
  if(evtCuts->IsTriggerRequired())  
  {
    // always MB
    isEventTriggered = inputHandler->IsEventSelected() & AliVEvent::kMB;

    physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if(!physicsSelection) return;
    //SetPhysicsTriggerSelection(physicsSelection);

    if (isEventTriggered && (GetTrigger() == AliTriggerAnalysis::kV0AND)) {
      // set trigger (V0AND)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
      if(!triggerAnalysis) return;
      isEventTriggered = triggerAnalysis->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }
  }

  // centrality determination
  Float_t centralityF = -1;
  AliCentrality *esdCentrality = esdEvent->GetCentrality();
  centralityF = esdCentrality->GetCentralityPercentile(fCentralityEstimator.Data());

  // use MC information
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);

  Int_t multMCTrueTracks = 0;
  if(IsUseMCInfo())
  {
    //
    if(!mcEvent) {
      AliDebug(AliLog::kError, "mcEvent not available");
      return;
    }
    // get MC event header
    header = mcEvent->Header();
    if (!header) {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }
    // MC particle stack
    stack = mcEvent->Stack();
    if (!stack) {
      AliDebug(AliLog::kError, "Stack not available");
      return;
    }

    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);

    // multipliticy of all MC primary tracks
    // in Zv, pt and eta ranges)
    multMCTrueTracks = AlidNdPtHelper::GetMCTrueTrackMult(mcEvent,evtCuts,accCuts);

  } // end bUseMC

  // get reconstructed vertex  
  //const AliESDVertex* vtxESD = 0; 
  const AliESDVertex* vtxESD = 0; 
  if(GetAnalysisMode() == AlidNdPtHelper::kTPC) {
        vtxESD = esdEvent->GetPrimaryVertexTPC();
  }
  else if(GetAnalysisMode() == AlidNdPtHelper::kTPCITS) {
     vtxESD = esdEvent->GetPrimaryVertexTracks();
  }
  else {
    	return;
  }

  if(!vtxESD) return;

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD); 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());


  // check event cuts
  if(isEventOK && isEventTriggered)
  {
    //
    // get IR information
    //
    AliESDHeader *esdHeader = 0;
    esdHeader = esdEvent->GetHeader();
    if(!esdHeader) return;
    //Int_t ir1 = esdHeader->GetTriggerIREntries(); //all ir-s
    //Int_t ir2 = esdHeader->GetTriggerIREntries(-1,1); // int2 set, 180 ms time interval
    //
    Int_t ir1 = 0;
    Int_t ir2 = 0;

    //
    Double_t vert[3] = {0}; 
    vert[0] = vtxESD->GetXv();
    vert[1] = vtxESD->GetYv();
    vert[2] = vtxESD->GetZv();
    Int_t mult = vtxESD->GetNContributors();
    Double_t bz = esdEvent->GetMagneticField();
    Double_t runNumber = esdEvent->GetRunNumber();
    Double_t evtTimeStamp = esdEvent->GetTimeStamp();
    Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();

    // high pT tracks
    for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++)
    {
      AliESDtrack *track = esdEvent->GetTrack(iTrack);
      if(!track) continue;
      if(track->Charge()==0) continue;
      if(!esdTrackCuts->AcceptTrack(track)) continue;
      if(!accCuts->AcceptTrack(track)) continue;
      
      // downscale low-pT tracks
      Double_t scalempt= TMath::Min(track->Pt(),10.);
      Double_t downscaleF = gRandom->Rndm();
      downscaleF *= fLowPtTrackDownscaligF;
      if(TMath::Exp(2*scalempt)<downscaleF) continue;
      //printf("TMath::Exp(2*scalempt) %e, downscaleF %e \n",TMath::Exp(2*scalempt), downscaleF);

      // Dump to the tree 
      // vertex
      // TPC constrained tracks
      // InnerParams constrained tracks
      // TPC-ITS tracks
      // ITSout-InnerParams tracks
      // chi2 distance between TPC constrained and TPC-ITS tracks
      // chi2 distance between TPC refitted constrained and TPC-ITS tracks
      // chi2 distance between ITSout and InnerParams tracks
      // MC information
      
      Double_t x[3]; track->GetXYZ(x);
      Double_t b[3]; AliTracker::GetBxByBz(x,b);

      //
      // Transform TPC inner params to track reference frame
      //
      Bool_t isOKtpcInner = kFALSE;
      AliExternalTrackParam * tpcInner = (AliExternalTrackParam *)(track->GetTPCInnerParam());
      if (tpcInner) {
        // transform to the track reference frame 
        isOKtpcInner = tpcInner->Rotate(track->GetAlpha());
        isOKtpcInner = tpcInner->PropagateTo(track->GetX(),esdEvent->GetMagneticField());
      }

      //
      // Constrain TPC inner to vertex
      // clone TPCinner has to be deleted
      //
      Bool_t isOKtpcInnerC = kFALSE;
      AliExternalTrackParam * tpcInnerC = new AliExternalTrackParam(*(track->GetTPCInnerParam()));
      if (tpcInnerC) {
        isOKtpcInnerC = ConstrainTPCInner(tpcInnerC,vtxESD,b);
        isOKtpcInnerC = tpcInnerC->Rotate(track->GetAlpha());
        isOKtpcInnerC = tpcInnerC->PropagateTo(track->GetX(),esdEvent->GetMagneticField());
      }

      //
      // Constrain TPC refitted tracks at inner TPC wall (InnerParams) to vertex  
      // Clone track InnerParams has to be deleted
      //
      Bool_t isOKtrackInnerC = kFALSE;
      AliExternalTrackParam * trackInnerC =  new AliExternalTrackParam(*(track->GetInnerParam()));
      if (trackInnerC) {
        isOKtrackInnerC = ConstrainTrackInner(trackInnerC,vtxESD,track->GetMass(),b);
        isOKtrackInnerC = trackInnerC->Rotate(track->GetAlpha());
        isOKtrackInnerC = trackInnerC->PropagateTo(track->GetX(),esdEvent->GetMagneticField());
      } 
      
      //
      // calculate chi2 between vi and vj vectors
      // with covi and covj covariance matrices
      // chi2ij = (vi-vj)^(T)*(covi+covj)^(-1)*(vi-vj)
      //
      TMatrixD deltaT(5,1), deltaTtrackC(5,1);
      TMatrixD delta(1,5),  deltatrackC(1,5);
      TMatrixD covarM(5,5), covarMtrackC(5,5);
      TMatrixD chi2(1,1);
      TMatrixD chi2trackC(1,1);

      if(isOKtpcInnerC && isOKtrackInnerC) 
      {
        for (Int_t ipar=0; ipar<5; ipar++) {
          deltaT(ipar,0)=tpcInnerC->GetParameter()[ipar]-track->GetParameter()[ipar];
	  delta(0,ipar)=tpcInnerC->GetParameter()[ipar]-track->GetParameter()[ipar];

          deltaTtrackC(ipar,0)=trackInnerC->GetParameter()[ipar]-track->GetParameter()[ipar];
	  deltatrackC(0,ipar)=trackInnerC->GetParameter()[ipar]-track->GetParameter()[ipar];

          for (Int_t jpar=0; jpar<5; jpar++) {
	    Int_t index=track->GetIndex(ipar,jpar);
	    covarM(ipar,jpar)=track->GetCovariance()[index]+tpcInnerC->GetCovariance()[index];
	    covarMtrackC(ipar,jpar)=track->GetCovariance()[index]+trackInnerC->GetCovariance()[index];
          }
        }

        // chi2 distance TPC constrained and TPC+ITS
        TMatrixD covarMInv = covarM.Invert();
        TMatrixD mat2 = covarMInv*deltaT;
        chi2 = delta*mat2; 
        //chi2.Print();

        // chi2 distance TPC refitted constrained and TPC+ITS
        TMatrixD covarMInvtrackC = covarMtrackC.Invert();
        TMatrixD mat2trackC = covarMInvtrackC*deltaTtrackC;
        chi2trackC = deltatrackC*mat2trackC; 
        //chi2trackC.Print();
      }


      //
      // Propagate ITSout to TPC inner wall 
      // and calculate chi2 distance to track (InnerParams)
      //
      const Double_t kTPCRadius=85; 
      const Double_t kStep=3; 

      // clone track InnerParams has to be deleted
      Bool_t isOKtrackInnerC2 = kFALSE;
      AliExternalTrackParam *trackInnerC2 = new AliExternalTrackParam(*(track->GetInnerParam()));
      if (trackInnerC2) {
        isOKtrackInnerC2 = AliTracker::PropagateTrackToBxByBz(trackInnerC2,kTPCRadius,track->GetMass(),kStep,kFALSE);
      }

      Bool_t isOKouterITSc = kFALSE;
      AliExternalTrackParam *outerITSc = NULL;
      TMatrixD chi2OuterITS(1,1);

      if(esdFriend && esdFriend->TestSkipBit()==kFALSE) 
      {
        // propagate ITSout to TPC inner wall
        AliESDfriendTrack *friendTrack = esdFriend->GetTrack(iTrack);

        if(friendTrack) 
	{
          outerITSc = new AliExternalTrackParam(*friendTrack->GetITSOut());
          if(outerITSc) 
	  {
            isOKouterITSc = AliTracker::PropagateTrackToBxByBz(outerITSc,kTPCRadius,track->GetMass(),kStep,kFALSE);
            isOKouterITSc = outerITSc->Rotate(trackInnerC2->GetAlpha());
            isOKouterITSc = outerITSc->PropagateTo(trackInnerC2->GetX(),esdEvent->GetMagneticField());

	    //
            // calculate chi2 between outerITS and innerParams
	    // cov without z-coordinate at the moment
	    //
            TMatrixD deltaTouterITS(4,1);
            TMatrixD deltaouterITS(1,4);
            TMatrixD covarMouterITS(4,4);

            if(isOKtrackInnerC2 && isOKouterITSc) {
	      Int_t kipar = 0;
	      Int_t kjpar = 0;
              for (Int_t ipar=0; ipar<5; ipar++) {
		if(ipar!=1) {
                  deltaTouterITS(kipar,0)=outerITSc->GetParameter()[ipar]-trackInnerC2->GetParameter()[ipar];
	          deltaouterITS(0,kipar)=outerITSc->GetParameter()[ipar]-trackInnerC2->GetParameter()[ipar];
		}

                kjpar=0;
                for (Int_t jpar=0; jpar<5; jpar++) {
	          Int_t index=outerITSc->GetIndex(ipar,jpar);
		  if(ipar !=1 || jpar!=1) {
	            covarMouterITS(kipar,kjpar)=outerITSc->GetCovariance()[index]+trackInnerC2->GetCovariance()[index];
		  }
                  if(jpar!=1)  kjpar++;
		}
	        if(ipar!=1) kipar++;
	      }

              // chi2 distance ITSout and InnerParams
              TMatrixD covarMInvouterITS = covarMouterITS.Invert();
              TMatrixD mat2outerITS = covarMInvouterITS*deltaTouterITS;
              chi2OuterITS = deltaouterITS*mat2outerITS; 
              //chi2OuterITS.Print();
	    } 
          }
        }
      }

      //
      // MC info
      //
      TParticle *particle=NULL, *particleTPC=NULL, *particleITS=NULL;
      TParticle *particleMother=NULL, *particleMotherTPC=NULL, *particleMotherITS=NULL;
      Int_t mech=-1, mechTPC=-1, mechITS=-1;
      Bool_t isPrim=kFALSE, isPrimTPC=kFALSE, isPrimITS=kFALSE;
      Bool_t isFromStrangess=kFALSE, isFromStrangessTPC=kFALSE, isFromStrangessITS=kFALSE;
      Bool_t isFromConversion=kFALSE, isFromConversionTPC=kFALSE, isFromConversionITS=kFALSE;
      Bool_t isFromMaterial=kFALSE, isFromMaterialTPC=kFALSE, isFromMaterialITS=kFALSE;

      AliTrackReference *refTPCIn = NULL;
      AliTrackReference *refITS = NULL;

      Bool_t isOKtrackInnerC3 = kFALSE;
      AliExternalTrackParam *trackInnerC3 = new AliExternalTrackParam(*(track->GetInnerParam()));

      if(IsUseMCInfo()) 
      {
        if(!stack) return;

        //
        // global track
	//
        Int_t label = TMath::Abs(track->GetLabel()); 
        particle = stack->Particle(label);
        if(particle && particle->GetPDG() && particle->GetPDG()->Charge()!=0.)
	{
	  particleMother = GetMother(particle,stack);
          mech = particle->GetUniqueID();
          isPrim = stack->IsPhysicalPrimary(label);
	  isFromStrangess  = IsFromStrangeness(label,stack);
	  isFromConversion = IsFromConversion(label,stack);
          isFromMaterial   = IsFromMaterial(label,stack);
	}

        //
	// TPC track
	//
	Int_t labelTPC = TMath::Abs(track->GetTPCLabel()); 
        particleTPC = stack->Particle(labelTPC);
        if(particleTPC && particleTPC->GetPDG() && particleTPC->GetPDG()->Charge()!=0.)
	{
	  particleMotherTPC = GetMother(particleTPC,stack);
          mechTPC = particleTPC->GetUniqueID();
          isPrimTPC = stack->IsPhysicalPrimary(labelTPC);
	  isFromStrangessTPC  = IsFromStrangeness(labelTPC,stack);
	  isFromConversionTPC = IsFromConversion(labelTPC,stack);
          isFromMaterialTPC   = IsFromMaterial(labelTPC,stack);
	}

        //
        // store first track reference
	// for TPC track
	//
        TParticle *part=0;
        TClonesArray *trefs=0;
        Int_t status = mcEvent->GetParticleAndTR(track->GetTPCLabel(), part, trefs);

	if(status>0 && part && trefs && part->GetPDG() && part->GetPDG()->Charge()!=0.) 
	{
	  Int_t nTrackRef = trefs->GetEntries();
	  //printf("nTrackRef %d \n",nTrackRef);

          Int_t countITS = 0;
	  for (Int_t iref = 0; iref < nTrackRef; iref++) 
          {
            AliTrackReference *ref = (AliTrackReference *)trefs->At(iref);

             // Ref. in the middle ITS 
            if(ref && ref->DetectorId()==AliTrackReference::kITS)
            {
	      if(!refITS && countITS==2) {
	        refITS = ref;
	        //printf("refITS %p \n",refITS);
	      }
	      countITS++;
            }

            // TPC
            if(ref && ref->DetectorId()==AliTrackReference::kTPC)
            {
	      if(!refTPCIn) {
	        refTPCIn = ref;
	        //printf("refTPCIn %p \n",refTPCIn);
	        //break;
	      }
            }
	  }

          // transform inner params to TrackRef
	  // reference frame
          if(refTPCIn && trackInnerC3) 
	  {
	    Double_t kRefPhi = TMath::ATan2(refTPCIn->Y(),refTPCIn->X());
            isOKtrackInnerC3 = trackInnerC3->Rotate(kRefPhi);
            isOKtrackInnerC3 = AliTracker::PropagateTrackToBxByBz(trackInnerC3,refTPCIn->R(),track->GetMass(),kStep,kFALSE);
	  }
        }

        //
	// ITS track
	//
	Int_t labelITS = TMath::Abs(track->GetITSLabel()); 
        particleITS = stack->Particle(labelITS);
        if(particleITS && particleITS->GetPDG() && particleITS->GetPDG()->Charge()!=0.)
	{
	  particleMotherITS = GetMother(particleITS,stack);
          mechITS = particleITS->GetUniqueID();
          isPrimITS = stack->IsPhysicalPrimary(labelITS);
	  isFromStrangessITS  = IsFromStrangeness(labelITS,stack);
	  isFromConversionITS = IsFromConversion(labelITS,stack);
          isFromMaterialITS   = IsFromMaterial(labelITS,stack);
        }
      }

      //
      Bool_t dumpToTree=kFALSE;
      
      if(isOKtpcInnerC  && isOKtrackInnerC) dumpToTree = kTRUE;
      if(fUseESDfriends && isOKtrackInnerC2 && isOKouterITSc) dumpToTree = kTRUE;
      if(fUseMCInfo     && isOKtrackInnerC3) dumpToTree = kTRUE;

      //
      if(fTreeSRedirector && dumpToTree) 
      {
        (*fTreeSRedirector)<<"dNdPtTree"<<
          "fileName.="<<&fileName<<
          "runNumber="<<runNumber<<
          "evtTimeStamp="<<evtTimeStamp<<
          "evtNumberInFile="<<evtNumberInFile<<
          "Bz="<<bz<<
	  "vertX="<<vert[0]<<
	  "vertY="<<vert[1]<<
	  "vertZ="<<vert[2]<<
	  "IRtot="<<ir1<<
	  "IRint2="<<ir2<<
          "mult="<<mult<<
          "esdTrack.="<<track<<
          "extTPCInnerC.="<<tpcInnerC<<
          "extInnerParamC.="<<trackInnerC<<
          "extInnerParam.="<<trackInnerC2<<
          "extOuterITS.="<<outerITSc<<
          "extInnerParamRef.="<<trackInnerC3<<
          "refTPCIn.="<<refTPCIn<<
          "refITS.="<<refITS<<
          "chi2TPCInnerC="<<chi2(0,0)<<
          "chi2InnerC="<<chi2trackC(0,0)<<
          "chi2OuterITS="<<chi2OuterITS(0,0)<<
          "centralityF="<<centralityF<<
          "particle.="<<particle<<
       	  "particleMother.="<<particleMother<<
          "mech="<<mech<<
          "isPrim="<<isPrim<<
	  "isFromStrangess="<<isFromStrangess<<
	  "isFromConversion="<<isFromConversion<<
          "isFromMaterial="<<isFromMaterial<<
          "particleTPC.="<<particleTPC<<
       	  "particleMotherTPC.="<<particleMotherTPC<<
          "mechTPC="<<mechTPC<<
          "isPrimTPC="<<isPrimTPC<<
	  "isFromStrangessTPC="<<isFromStrangessTPC<<
	  "isFromConversionTPC="<<isFromConversionTPC<<
          "isFromMaterialTPC="<<isFromMaterialTPC<<
          "particleITS.="<<particleITS<<
       	  "particleMotherITS.="<<particleMotherITS<<
          "mechITS="<<mechITS<<
          "isPrimITS="<<isPrimITS<<
	  "isFromStrangessITS="<<isFromStrangessITS<<
	  "isFromConversionITS="<<isFromConversionITS<<
          "isFromMaterialITS="<<isFromMaterialITS<<
          "\n";
        }
      
	if(tpcInnerC) delete tpcInnerC;
	if(trackInnerC) delete trackInnerC;
	if(trackInnerC2) delete trackInnerC2;
	if(outerITSc) delete outerITSc;
	if(trackInnerC3) delete trackInnerC3;
    }
  }
  
  PostData(1, fOutput);
}


//_____________________________________________________________________________
void AlidNdPtTrackDumpTask::ProcessMCEff(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent, AliESDfriend *const /*esdFriend*/)
{
  //
  // Fill tree for efficiency studies MC only

  if(!esdEvent) {
    AliDebug(AliLog::kError, "esdEvent not available");
    return;
  }

   if(!mcEvent) {
    AliDebug(AliLog::kError, "mcEvent not available");
    return;
  }

  // get selection cuts
  AlidNdPtEventCuts *evtCuts = GetEventCuts(); 
  AlidNdPtAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts  || !esdTrackCuts) {
    AliDebug(AliLog::kError, "cuts not available");
    return;
  }

  // trigger selection
  Bool_t isEventTriggered = kTRUE;
  AliPhysicsSelection *physicsSelection = NULL;
  AliTriggerAnalysis* triggerAnalysis = NULL;

  // 
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return;
  }
   
  // get file name
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { 
    Printf("ERROR: Could not receive input chain");
    return;
  }
  TObjString fileName(chain->GetCurrentFile()->GetName());

  // trigger
  if(evtCuts->IsTriggerRequired())  
  {
    // always MB
    isEventTriggered = inputHandler->IsEventSelected() & AliVEvent::kMB;

    physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if(!physicsSelection) return;

    if (isEventTriggered && (GetTrigger() == AliTriggerAnalysis::kV0AND)) {
      // set trigger (V0AND)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
      if(!triggerAnalysis) return;
      isEventTriggered = triggerAnalysis->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }
  }

  // centrality determination
  Float_t centralityF = -1;
  AliCentrality *esdCentrality = esdEvent->GetCentrality();
  centralityF = esdCentrality->GetCentralityPercentile(fCentralityEstimator.Data());

  // use MC information
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);

  Int_t multMCTrueTracks = 0;
  if(IsUseMCInfo())
  {
    //
    if(!mcEvent) {
      AliDebug(AliLog::kError, "mcEvent not available");
      return;
    }
    // get MC event header
    header = mcEvent->Header();
    if (!header) {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }
    // MC particle stack
    stack = mcEvent->Stack();
    if (!stack) {
      AliDebug(AliLog::kError, "Stack not available");
      return;
    }

    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);

    // multipliticy of all MC primary tracks
    // in Zv, pt and eta ranges)
    multMCTrueTracks = AlidNdPtHelper::GetMCTrueTrackMult(mcEvent,evtCuts,accCuts);

  } // end bUseMC

  // get reconstructed vertex  
  //const AliESDVertex* vtxESD = 0; 
  const AliESDVertex* vtxESD = 0; 
  if(GetAnalysisMode() == AlidNdPtHelper::kTPC) {
        vtxESD = esdEvent->GetPrimaryVertexTPC();
  }
  else if(GetAnalysisMode() == AlidNdPtHelper::kTPCITS) {
     vtxESD = esdEvent->GetPrimaryVertexTracks();
  }
  else {
    	return;
  }

  if(!vtxESD) return;

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD); 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());

  // check event cuts
  if(isEventOK && isEventTriggered)
  {
    if(IsUseMCInfo()) 
    {
      if(!stack) return;

      //
      // MC info
      //
      TParticle *particle=NULL;
      TParticle *particleMother=NULL;
      Int_t mech=-1;

      // reco event info
      Double_t vert[3] = {0}; 
      vert[0] = vtxESD->GetXv();
      vert[1] = vtxESD->GetYv();
      vert[2] = vtxESD->GetZv();
      Int_t mult = vtxESD->GetNContributors();
      Double_t bz = esdEvent->GetMagneticField();
      Double_t runNumber = esdEvent->GetRunNumber();
      Double_t evtTimeStamp = esdEvent->GetTimeStamp();
      Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();

      // loop over MC stack
      for (Int_t iMc = 0; iMc < stack->GetNtrack(); ++iMc) 
      {
         particle = stack->Particle(iMc);
         if (!particle)
         continue;

         // only charged particles
         if(!particle->GetPDG()) continue;
         Double_t charge = particle->GetPDG()->Charge()/3.;
         if (TMath::Abs(charge) < 0.001)
         continue;

         // only primary particles
         Bool_t prim = stack->IsPhysicalPrimary(iMc);
         if(!prim) continue;

         // downscale low-pT particles
         Double_t scalempt= TMath::Min(particle->Pt(),10.);
         Double_t downscaleF = gRandom->Rndm();
         downscaleF *= fLowPtTrackDownscaligF;
         if(TMath::Exp(2*scalempt)<downscaleF) continue;

         // is particle in acceptance
         if(!accCuts->AcceptTrack(particle)) continue;
       
         // check if particle reconstructed
         Bool_t isRec = kFALSE;
         Int_t  trackIndex = -1;
         for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++)
         {
           
           AliESDtrack *track = esdEvent->GetTrack(iTrack);
           if(!track) continue;
           if(track->Charge()==0) continue;
           if(esdTrackCuts->AcceptTrack(track) && accCuts->AcceptTrack(track)) 
           {
             Int_t label =  TMath::Abs(track->GetLabel());
             if(label == iMc) {
               isRec = kTRUE;
               trackIndex = iTrack;
               break;
             }
           } 
         }

         // Store information in the output tree
         AliESDtrack *recTrack = NULL; 
         if(trackIndex>-1)  { 
           recTrack = esdEvent->GetTrack(trackIndex); 
         } else {
           recTrack = new AliESDtrack(); 
         } 

	 particleMother = GetMother(particle,stack);
         mech = particle->GetUniqueID();

         //MC particle track length
         Double_t tpcTrackLength = 0.;
         AliMCParticle *mcParticle = (AliMCParticle*) mcEvent->GetTrack(iMc);
         if(mcParticle) {
           Int_t counter;
           tpcTrackLength = mcParticle->GetTPCTrackLength(bz,0.05,counter,3.0);
         } 


         //
         if(fTreeSRedirector) {
           (*fTreeSRedirector)<<"MCEffTree"<<
           "fileName.="<<&fileName<<
           "runNumber="<<runNumber<<
           "evtTimeStamp="<<evtTimeStamp<<
           "evtNumberInFile="<<evtNumberInFile<<
           "Bz="<<bz<<
	   "vertX="<<vert[0]<<
	   "vertY="<<vert[1]<<
	   "vertZ="<<vert[2]<<
           "mult="<<mult<<
           "esdTrack.="<<recTrack<<
           "isRec="<<isRec<<
           "tpcTrackLength="<<tpcTrackLength<<
           "particle.="<<particle<<
       	   "particleMother.="<<particleMother<<
           "mech="<<mech<<
           "\n";
         }

         if(trackIndex <0 && recTrack) delete recTrack; recTrack=0;
      }
    }
  }
  
  PostData(1, fOutput);
}

//_____________________________________________________________________________
Bool_t AlidNdPtTrackDumpTask::IsHighDeDxParticle(AliESDtrack * track) {
  //
  // check if particle is Z > 1 
  //
  if (track->GetTPCNcls() < 60) return kFALSE;
  Double_t mom = track->GetInnerParam()->GetP();
  if (mom < 0.2) return kFALSE; // protection against unexpected behavior of Aleph parameterization
  Float_t dca[2], bCov[3];
  track->GetImpactParameters(dca,bCov);
  //

  Double_t triggerDeDx = 4*AliExternalTrackParam::BetheBlochAleph((mom*2)/(0.938*3),1.0288,31.9806,5.04114e-11,2.13096,2.38541);

  if (track->GetTPCsignal() > triggerDeDx && track->GetTPCsignal()<1000 && TMath::Abs(dca[0])<3.) return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________
void AlidNdPtTrackDumpTask::ProcessV0(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent, AliESDfriend *const /*esdFriend*/)
{
  //
  // Select real events with V0 (K0s and Lambda) high-pT candidates
  //
  if(!esdEvent) {
    AliDebug(AliLog::kError, "esdEvent not available");
    return;
  }

  // get selection cuts
  AlidNdPtEventCuts *evtCuts = GetEventCuts(); 
  AlidNdPtAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts  || !esdTrackCuts) {
    AliDebug(AliLog::kError, "cuts not available");
    return;
  }

  // trigger selection
  Bool_t isEventTriggered = kTRUE;
  AliPhysicsSelection *physicsSelection = NULL;
  AliTriggerAnalysis* triggerAnalysis = NULL;

  // 
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return;
  }
   
  // get file name
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { 
    Printf("ERROR: Could not receive input chain");
    return;
  }
  TObjString fileName(chain->GetCurrentFile()->GetName());

  // trigger
  if(evtCuts->IsTriggerRequired())  
  {
    // always MB
    isEventTriggered = inputHandler->IsEventSelected() & AliVEvent::kMB;

    physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if(!physicsSelection) return;
    //SetPhysicsTriggerSelection(physicsSelection);

    if (isEventTriggered && (GetTrigger() == AliTriggerAnalysis::kV0AND)) {
      // set trigger (V0AND)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
      if(!triggerAnalysis) return;
      isEventTriggered = triggerAnalysis->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }
  }

  // centrality determination
  Float_t centralityF = -1;
  AliCentrality *esdCentrality = esdEvent->GetCentrality();
  centralityF = esdCentrality->GetCentralityPercentile(fCentralityEstimator.Data());


  // get reconstructed vertex  
  //const AliESDVertex* vtxESD = 0; 
  const AliESDVertex* vtxESD = 0; 
  if(GetAnalysisMode() == AlidNdPtHelper::kTPC) {
        vtxESD = esdEvent->GetPrimaryVertexTPC();
  }
  else if(GetAnalysisMode() == AlidNdPtHelper::kTPCITS) {
     vtxESD = esdEvent->GetPrimaryVertexTracks();
  }
  else {
    	return;
  }

  if(!vtxESD) return;

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD); 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());

  // check event cuts
  if(isEventOK && isEventTriggered) {
  //
  // Dump the pt downscaled V0 into the tree
  // 
  Int_t ntracks = esdEvent->GetNumberOfTracks();
  Int_t nV0s = esdEvent->GetNumberOfV0s();
  Int_t run = esdEvent->GetRunNumber();
  Int_t time= esdEvent->GetTimeStamp();
  Int_t evNr=esdEvent->GetEventNumberInFile();
  Double_t bz = esdEvent->GetMagneticField();


  for (Int_t iv0=0; iv0<nV0s; iv0++){
    AliESDv0 * v0 = esdEvent->GetV0(iv0);
    if (!v0) continue;
    AliESDtrack * track0 = esdEvent->GetTrack(v0->GetIndex(0));
    AliESDtrack * track1 = esdEvent->GetTrack(v0->GetIndex(1));
    if (!track0) continue;
    if (!track1) continue;
    if (track0->GetSign()<0) {
      track1 = esdEvent->GetTrack(v0->GetIndex(0));
      track0 = esdEvent->GetTrack(v0->GetIndex(1));
    }
    //
    Bool_t isDownscaled = IsV0Downscaled(v0);
    if (isDownscaled) continue;
    AliKFParticle kfparticle; //
    Int_t type=GetKFParticle(v0,esdEvent,kfparticle);
    if (type==0) continue;   

    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"V0s"<<
      "isDownscaled="<<isDownscaled<<
      "Bz="<<bz<<
      "fileName.="<<&fileName<<
      "runNumber="<<run<<
      "evtTimeStamp="<<time<<
      "evtNumberInFile="<<evNr<<
      "type="<<type<<
      "ntracks="<<ntracks<<
      "v0.="<<v0<<
      "kf.="<<&kfparticle<<
      "track0.="<<track0<<
      "track1.="<<track1<<
      "centralityF="<<centralityF<<
      "\n";
  }
  }
  PostData(1, fOutput);
}

//_____________________________________________________________________________
void AlidNdPtTrackDumpTask::ProcessdEdx(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent, AliESDfriend *const /*esdFriend*/)
{
  //
  // Select real events with large TPC dEdx signal
  //
  if(!esdEvent) {
    AliDebug(AliLog::kError, "esdEvent not available");
    return;
  }

  // get selection cuts
  AlidNdPtEventCuts *evtCuts = GetEventCuts(); 
  AlidNdPtAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts  || !esdTrackCuts) {
    AliDebug(AliLog::kError, "cuts not available");
    return;
  }

  // get file name
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { 
    Printf("ERROR: Could not receive input chain");
    return;
  }
  TObjString fileName(chain->GetCurrentFile()->GetName());

  // trigger
  Bool_t isEventTriggered = kTRUE;
  AliPhysicsSelection *physicsSelection = NULL;
  AliTriggerAnalysis* triggerAnalysis = NULL;

  // 
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return;
  }

  if(evtCuts->IsTriggerRequired())  
  {
    // always MB
    isEventTriggered = inputHandler->IsEventSelected() & AliVEvent::kMB;

    physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if(!physicsSelection) return;

    if (isEventTriggered && (GetTrigger() == AliTriggerAnalysis::kV0AND)) {
      // set trigger (V0AND)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
      if(!triggerAnalysis) return;
      isEventTriggered = triggerAnalysis->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }
  }

  // get reconstructed vertex  
  const AliESDVertex* vtxESD = 0; 
  if(GetAnalysisMode() == AlidNdPtHelper::kTPC) {
        vtxESD = esdEvent->GetPrimaryVertexTPC();
  }
  else if(GetAnalysisMode() == AlidNdPtHelper::kTPCITS) {
     vtxESD = esdEvent->GetPrimaryVertexTracks();
  }
  else {
    	return;
  }
  if(!vtxESD) return;

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD); 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());


  // check event cuts
  if(isEventOK && isEventTriggered)
  {
    Double_t vert[3] = {0}; 
    vert[0] = vtxESD->GetXv();
    vert[1] = vtxESD->GetYv();
    vert[2] = vtxESD->GetZv();
    Int_t mult = vtxESD->GetNContributors();
    Double_t bz = esdEvent->GetMagneticField();
    Double_t runNumber = esdEvent->GetRunNumber();
    Double_t evtTimeStamp = esdEvent->GetTimeStamp();
    Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();

    // large dEdx
    for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++)
    {
      AliESDtrack *track = esdEvent->GetTrack(iTrack);
      if(!track) continue;
      if(track->Charge()==0) continue;
      if(!esdTrackCuts->AcceptTrack(track)) continue;
      if(!accCuts->AcceptTrack(track)) continue;

      if(!IsHighDeDxParticle(track)) continue;
      
      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"dEdx"<<
      "fileName.="<<&fileName<<
      "runNumber="<<runNumber<<
      "evtTimeStamp="<<evtTimeStamp<<
      "evtNumberInFile="<<evtNumberInFile<<
      "Bz="<<bz<<
      "vertX="<<vert[0]<<
      "vertY="<<vert[1]<<
      "vertZ="<<vert[2]<<
      "mult="<<mult<<
      "esdTrack.="<<track<<
      "\n";
    }
  }
}

//_____________________________________________________________________________
Int_t   AlidNdPtTrackDumpTask::GetKFParticle(AliESDv0 *const v0, AliESDEvent * const event, AliKFParticle & kfparticle)
{
  //
  // Create KF particle in case the V0 fullfill selection criteria
  //
  // Selection criteria
  //  0. algorithm cut
  //  1. track cut
  //  3. chi2 cut
  //  4. rough mass cut
  //  5. Normalized pointing angle cut
  //
  const Double_t cutMass=0.2;
  const Double_t kSigmaDCACut=3;
  //
  // 0.) algo cut - accept only on the fly
  //
  if (v0->GetOnFlyStatus() ==kFALSE) return 0;     
  //
  // 1.) track cut
  // 
  AliESDtrack * track0 = event->GetTrack(v0->GetIndex(0));
  AliESDtrack * track1 = event->GetTrack(v0->GetIndex(1));
  /*
    TCut cutD="abs(track0.fD/sqrt(track0.fCdd))>2&&abs(track1.fD/sqrt(track1.fCdd))>2";
    TCut cutTheta="abs(track0.fP[3])<1&&abs(track1.fP[3])<1";
    TCut cutNcl="track0.GetTPCClusterInfo(2,1)>100&&track1.GetTPCClusterInfo(2,1)>100";
  */  
  if (TMath::Abs(track0->GetTgl())>1) return 0;
  if (TMath::Abs(track1->GetTgl())>1) return 0;
  if ((track0->GetTPCClusterInfo(2,1))<100) return 0;
  if ((track1->GetTPCClusterInfo(2,1))<100) return 0;
  //if ((track0->GetITSclusters(0))<2) return 0;
  //if ((track1->GetITSclusters(0))<2) return 0; 
  Float_t pos0[2]={0}, cov0[3]={0};
  Float_t pos1[2]={0}, cov1[3]={0};
  track0->GetImpactParameters(pos0,cov0);
  track0->GetImpactParameters(pos1,cov1);
  //
  if (TMath::Abs(pos0[0])<kSigmaDCACut*TMath::Sqrt(cov0[0])) return 0;
  if (TMath::Abs(pos1[0])<kSigmaDCACut*TMath::Sqrt(cov1[0])) return 0;
  // 
  //
  // 3.) Chi2 cut
  //
  Double_t chi2KF = v0->GetKFInfo(2,2,2);
  if (chi2KF>25) return 0;
  //
  // 4.) Rough mass cut - 0.200 GeV
  //
  static Double_t masses[2]={-1};
  if (masses[0]<0){
    masses[0] = TDatabasePDG::Instance()->GetParticle("K_S0")->Mass();
    masses[1] = TDatabasePDG::Instance()->GetParticle("Lambda0")->Mass();
  }
  Double_t mass00=  v0->GetEffMass(0,0);
  Double_t mass22=  v0->GetEffMass(2,2);
  Double_t mass42=  v0->GetEffMass(4,2);
  Double_t mass24=  v0->GetEffMass(2,4);
  Bool_t massOK=kFALSE;
  Int_t type=0;
  Int_t ptype=0;
  Double_t dmass=1;
  Int_t p1=0, p2=0;
  if (TMath::Abs(mass00-0)<cutMass) {
    massOK=kTRUE; type+=1; 
    if (TMath::Abs(mass00-0)<dmass) {
      ptype=1;
      dmass=TMath::Abs(mass00-0);      
      p1=0; p2=0;
    } 
  }
  if (TMath::Abs(mass24-masses[1])<cutMass) {
    massOK=kTRUE; type+=2; 
    if (TMath::Abs(mass24-masses[1])<dmass){
      dmass = TMath::Abs(mass24-masses[1]);
      ptype=2;
      p1=2; p2=4;
    }
  }
  if (TMath::Abs(mass42-masses[1])<cutMass) {
    massOK=kTRUE; type+=4;
    if (TMath::Abs(mass42-masses[1])<dmass){
      dmass = TMath::Abs(mass42-masses[1]);
      ptype=4;
      p1=4; p2=2;
    }
  }
  if (TMath::Abs(mass22-masses[0])<cutMass) {
    massOK=kTRUE; type+=8;
    if (TMath::Abs(mass22-masses[0])<dmass){
      dmass = TMath::Abs(mass22-masses[0]);
      ptype=8;
      p1=2; p2=2;
    }
  }
  if (type==0) return 0;
  //
  const Int_t spdg[5]={kPositron,kMuonPlus,kPiPlus, kKPlus, kProton};
  const AliExternalTrackParam *paramP = v0->GetParamP();
  const AliExternalTrackParam *paramN = v0->GetParamN();
  if (paramP->GetSign()<0){
    paramP=v0->GetParamP();
    paramN=v0->GetParamN();
  }
  //Double_t *pparam1 = (Double_t*)paramP->GetParameter();
  //Double_t *pparam2 = (Double_t*)paramN->GetParameter();
  //
  AliKFParticle kfp1( *paramP, spdg[p1]  );
  AliKFParticle kfp2( *paramN, -1 *spdg[p2]  );
  AliKFParticle V0KF;
  (V0KF)+=kfp1;
  (V0KF)+=kfp2;
  kfparticle=V0KF;
  //
  // Pointing angle
  //
  Double_t  errPhi    = V0KF.GetErrPhi();
  Double_t  pointAngle= TMath::ACos(v0->GetV0CosineOfPointingAngle());
  if (pointAngle/errPhi>10) return 0;  
  //
  return ptype;  
}

//_____________________________________________________________________________
Bool_t AlidNdPtTrackDumpTask::IsV0Downscaled(AliESDv0 *const v0)
{
  //
  // Downscale randomly low pt V0
  //
  //return kFALSE;
  Double_t maxPt= TMath::Max(v0->GetParamP()->Pt(), v0->GetParamN()->Pt());
  Double_t scalempt= TMath::Min(maxPt,10.);
  Double_t downscaleF = gRandom->Rndm();
  downscaleF *= fLowPtV0DownscaligF;
  
  //printf("V0 TMath::Exp(2*scalempt) %e, downscaleF %e \n",TMath::Exp(2*scalempt), downscaleF);
  if (TMath::Exp(2*scalempt)<downscaleF) return kTRUE;
  return kFALSE;

  /*
    TH1F his1("his1","his1",100,0,10);
    TH1F his2("his2","his2",100,0,10);
    {for (Int_t i=0; i<10000; i++){
       Double_t rnd=gRandom->Exp(1);
       Bool_t isDownscaled =TMath::Exp(rnd)<100*gRandom->Rndm();
       his1->Fill(rnd); 
       if (!isDownscaled) his2->Fill(rnd); 
    }}

   */

}


//_____________________________________________________________________________
/*
Bool_t AlidNdPtTrackDumpTask::MakeTPCInnerC(AliESDtrack *const track, const AliESDVertex* vtx, Double_t b[3])
{
//
// return TPC inner constrain object (must be deleted)
//

if(!track) return NULL;
if(!vtx) return NULL;
  
  AliExternalTrackParam * tpcInnerC = NULL;
  Bool_t isOK = kFALSE;

  tpcInnerC = new AliExternalTrackParam(*(track->GetTPCInnerParam()));
  if (!tpcInnerC) return NULL;
  isOK = ConstrainTPCInner(tpcInnerC,vtx,b);
  isOK = tpcInnerC->Rotate(track->GetAlpha());
  isOK = tpcInnerC->PropagateTo(track->GetX(),b[2]);
  if(!isOK) {
    delete tpcInnerC;
    return NULL;
  }

  if(fTreeSRedirector) {
    (*fTreeSRedirector)<<"dNdPtTree"<<
    (*fTreeSRedirector)<<"dNdPtTree"<<
    "esdTrack.="<<track<<
    "extTPCInnerC.="<<tpcInnerC<<
    "chi2TPCInnerC="<<chi2TPCInnerC;
  }

return tpcInnerC;
}
*/

//_____________________________________________________________________________
Bool_t AlidNdPtTrackDumpTask::ConstrainTPCInner(AliExternalTrackParam *const tpcInnerC, const AliESDVertex* vtx, Double_t b[3])
{
 // Constrain TPC inner params constrained
 //
      if(!tpcInnerC) return kFALSE; 
      if(!vtx) return kFALSE;

      Double_t dz[2],cov[3];
      //AliESDVertex *vtx= (AliESDVertex *)esdEvent->GetPrimaryVertex();
      //if(!tpcInnerC->PropagateToDCA(vtx, esdEvent->GetMagneticField(), 3, dz, cov)) return kFALSE; 
      //if(!tpcInnerC->PropagateToDCA(vtx, Bz, 3, dz, cov)) return kFALSE; 
      if(!tpcInnerC->PropagateToDCABxByBz(vtx, b, 3, dz, cov)) return kFALSE; 


      Double_t covar[6]; vtx->GetCovMatrix(covar);
      Double_t p[2]={tpcInnerC->GetParameter()[0]-dz[0],tpcInnerC->GetParameter()[1]-dz[1]};
      Double_t c[3]={covar[2],0.,covar[5]};
      Double_t chi2C=tpcInnerC->GetPredictedChi2(p,c);
      if (chi2C>kVeryBig) return kFALSE; 

      if(!tpcInnerC->Update(p,c)) return kFALSE;

  return kTRUE;
}


//_____________________________________________________________________________
/*
AliExternalTrackParam * AlidNdPtTrackDumpTask::CalculateChi2(AliESDtrack *const track, AliExternalTrackParam *const trackParam)
{
//
// calculate chi2 between vi and vj vectors
// with covi and covj covariance matrices
// chi2ij = (vi-vj)^(T)*(covi+covj)^(-1)*(vi-vj)
//

if(!track) return 0;
if(!trackParam) return 0;

  TMatrixD deltaT(5,1);
  TMatrixD delta(1,5);
  TMatrixD covarM(5,5);

  for (Int_t ipar=0; ipar<5; ipar++) {
    deltaT(ipar,0)=trackParam->GetParameter()[ipar]-track->GetParameter()[ipar];
    delta(0,ipar)=trackParam->GetParameter()[ipar]-track->GetParameter()[ipar];
    for (Int_t jpar=0; jpar<5; jpar++) {
      Int_t index=track->GetIndex(ipar,jpar);
      covarM(ipar,jpar)=track->GetCovariance()[index]+trackParam->GetCovariance()[index];
    }
  }

  // chi2 distance 
  TMatrixD covarMInv = covarM.Invert();
  TMatrixD mat2 = covarMInv*deltaT;
  TMatrixD chi2 = delta*mat2; 
  //chi2.Print();

return ((Double_t)chi(0,0));
}
*/

//_____________________________________________________________________________
/*
AliExternalTrackParam * AlidNdPtTrackDumpTask::CalculateChi2(AliExternalTrackParam *const trackParam1, AliExternalTrackParam *const trackParam2)
{
//
// calculate chi2 between vi and vj vectors
// with covi and covj covariance matrices
// chi2ij = (vi-vj)^(T)*(covi+covj)^(-1)*(vi-vj)
//

if(!track) return 0;
if(!trackParam) return 0;

TMatrixD deltaT(5,1);
TMatrixD delta(1,5);
TMatrixD covarM(5,5);

  for (Int_t ipar=0; ipar<5; ipar++) {
    deltaT(ipar,0)=trackParam1->GetParameter()[ipar]-trackParam2->GetParameter()[ipar];
    delta(0,ipar)=trackParam1->GetParameter()[ipar]-trackParam2->GetParameter()[ipar];
    for (Int_t jpar=0; jpar<5; jpar++) {
      Int_t index=track->GetIndex(ipar,jpar);
      covarM(ipar,jpar)=trackParam1->GetCovariance()[index]+trackParam2->GetCovariance()[index];
    }
  }

  // chi2 distance 
  TMatrixD covarMInv = covarM.Invert();
  TMatrixD mat2 = covarMInv*deltaT;
  TMatrixD chi2 = delta*mat2; 
  //chi2.Print();

return ((Double_t)chi(0,0));
}
*/


//_____________________________________________________________________________
/*
AliExternalTrackParam * AlidNdPtTrackDumpTask::MakeTrackInnerC(AliESDtrack *const track, const AliESDVertex* vtx, Double_t b[3])
{
//
// Constrain TPC refitted tracks at inner TPC wall (InnerParams) 
// to vertex
//
if(!track) return NULL;
if(!vtx) return NULL;

  // clone track InnerParams has to be deleted
  AliExternalTrackParam * trackInnerC =  new AliExternalTrackParam(*(track->GetInnerParam()));
  if (!trackInnerC) return NULL;
  Bool_t isOK = ConstrainTrackInner(trackInnerC,vtx,track->GetMass(),b);
  isOK = trackInnerC->Rotate(track->GetAlpha());
  isOK = trackInnerC->PropagateTo(track->GetX(),b[2]);
  if(!isOK) {
    if(trackInnerC) delete trackInnerC;
    return NULL;
  }
      
       
  // Dump to tree
  if(fTreeSRedirector) 
  {
    (*fTreeSRedirector)<<"dNdPtTree"<<
    "extInnerParamC.="<<trackInnerC<<
    "chi2InnerC="<<chi2InnerC<<
    "\n";
  }

return trackInnerC;
}
*/

//_____________________________________________________________________________
Bool_t AlidNdPtTrackDumpTask::ConstrainTrackInner(AliExternalTrackParam *const trackInnerC, const AliESDVertex* vtx, Double_t mass, Double_t b[3])
{
 // Constrain TPC inner params constrained
 //
      if(!trackInnerC) return kFALSE; 
      if(!vtx) return kFALSE;

      const Double_t kRadius  = 2.8; 
      const Double_t kMaxStep = 1.0; 

      Double_t dz[2],cov[3];

      //AliESDVertex *vtx= (AliESDVertex *)esdEvent->GetPrimaryVertex();
      //if(!trackInnerC->PropagateToDCA(vtx, esdEvent->GetMagneticField(), 3, dz, cov)) return kFALSE; 
      //if(!trackInnerC->PropagateToDCA(vtx, Bz, 3, dz, cov)) return kFALSE; 

      if(!AliTracker::PropagateTrackToBxByBz(trackInnerC,kRadius,mass,kMaxStep,kFALSE)) return kFALSE;
      if(!trackInnerC->PropagateToDCABxByBz(vtx, b, 3, dz, cov)) return kFALSE; 

      //
      Double_t covar[6]; vtx->GetCovMatrix(covar);
      Double_t p[2]={trackInnerC->GetParameter()[0]-dz[0],trackInnerC->GetParameter()[1]-dz[1]};
      Double_t c[3]={covar[2],0.,covar[5]};
      Double_t chi2C=trackInnerC->GetPredictedChi2(p,c);
      if (chi2C>kVeryBig) return kFALSE; 

      if(!trackInnerC->Update(p,c)) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
/*
AliExternalTrackParam * AlidNdPtTrackDumpTask::PropagateInnerParam(AliESDtrack *const track, Double_t radius, Double_t step)
{
//
// Propagate InnerParams to inner wall 
//
if(!track) return NULL;

  Bool_t isOK = kFALSE;

  // clone track InnerParams has to be deleted
  AliExternalTrackParam *trackInnerC2 = new AliExternalTrackParam(*(track->GetInnerParam()));
  if (!trackInnerC2) return NULL; 
  isOK = AliTracker::PropagateTrackToBxByBz(trackInnerC2,radius,track->GetMass(),step,kFALSE);
  if(!isOK) {
    delete trackInnerC2;
    return NULL;
  }

return trackInnerC2;
}
*/


//_____________________________________________________________________________
/*
AliExternalTrackParam * AlidNdPtTrackDumpTask::PropagateITSOut(Int_t iTrack, AliExternalTrackParam * trackParam, AliESDfriend *const esdFriend, Double_t b[3], Double_t radius, Double_t step)
{
//
// Propagate ITSout to TPC inner wall 
//

if(!track) return NULL;
if(!esdFriend) return NULL;

  Bool_t isOK = kFALSE;
  AliExternalTrackParam *outerITSc = 0;

  if(esdFriend && esdFriend->TestSkipBit()==kFALSE) 
  {
    // propagate ITSout to TPC inner wall
    AliESDfriendTrack *friendTrack = esdFriend->GetTrack(iTrack);
    if(friendTrack) 
    {
      if( (outerITSc = new AliExternalTrackParam(*friendTrack->GetITSOut())) ) 
      {
        if(AliTracker::PropagateTrackToBxByBz(outerITSc,radius,track->GetMass(),step,kFALSE))
        {
          // transform ITSout to the track InnerParams reference frame 
          isOK = outerITSc->Rotate(trackInnerC2->GetAlpha());
          isOK = outerITSc->PropagateTo(trackInnerC2->GetX(),b[2]);
          if(!isOK) {
	    if(trackInnerC2) delete trackInnerC2;
            if(outerITSc) delete outerITSc;
	    return kFALSE;
          }
              
          //
          // calculate chi2 between outerITS and innerParams
	  // cov without z-coordinate at the moment
	  // 
          TMatrixD deltaTouterITS(4,1);
          TMatrixD deltaouterITS(1,4);
          TMatrixD covarMouterITS(4,4);

	  Int_t kipar = 0;
	  Int_t kjpar = 0;
          for (Int_t ipar=0; ipar<5; ipar++) {
	    if(ipar!=1) {
              deltaTouterITS(kipar,0)=outerITSc->GetParameter()[ipar]-trackInnerC2->GetParameter()[ipar];
	      deltaouterITS(0,kipar)=outerITSc->GetParameter()[ipar]-trackInnerC2->GetParameter()[ipar];
	    }

            kjpar=0;
            for (Int_t jpar=0; jpar<5; jpar++) {
	      Int_t index=outerITSc->GetIndex(ipar,jpar);
	      if(ipar !=1 || jpar!=1) {
	        covarMouterITS(kipar,kjpar)=outerITSc->GetCovariance()[index]+trackInnerC2->GetCovariance()[index];
	      }
              if(jpar!=1)  kjpar++;
	    }
	    if(ipar!=1) kipar++;
	  }

          // chi2 distance ITSout and InnerParams
          TMatrixD covarMInvouterITS = covarMouterITS.Invert();
          TMatrixD mat2outerITS = covarMInvouterITS*deltaTouterITS;
          TMatrixD chi2OuterITS = deltaouterITS*mat2outerITS; 
          //chi2OuterITS.Print();
          chi2ITSout = chi2OuterITS(0,0);

          // dump to the tree
          if(fTreeSRedirector) 
          {
            (*fTreeSRedirector)<<"dNdPtTree"<<
            "extInnerParam.="<<trackInnerC2<<
            "extOuterITS.="<<outerITSc<<
            "chi2OuterITS="<<chi2ITSout;
          }
        } 
      }
    }
  }

  if(trackInnerC2) delete trackInnerC2;
  if(outerITSc) delete outerITSc;

return kTRUE;
}
*/

//_____________________________________________________________________________
/*
Bool_t AlidNdPtTrackDumpTask::UseMCInfoAndDump(AliMCEvent *const mcEvent,  AliESDtrack *const track, AliStack *const stack)
{
//
// MC info
//
if(!mcEvent) return kFALSE;
if(!track) return kFALSE;
if(!stack) return kFALSE;


  const Double_t kStep=3; 


  TParticle *particle=NULL, *particleTPC=NULL, *particleITS=NULL;
  TParticle *particleMother=NULL, *particleMotherTPC=NULL, *particleMotherITS=NULL;
  static Int_t mech=-1, mechTPC=-1, mechITS=-1;
  static Bool_t isPrim=kFALSE, isPrimTPC=kFALSE, isPrimITS=kFALSE;
  static Bool_t isFromStrangess=kFALSE, isFromStrangessTPC=kFALSE, isFromStrangessITS=kFALSE;
  static Bool_t isFromConversion=kFALSE, isFromConversionTPC=kFALSE, isFromConversionITS=kFALSE;
  static Bool_t isFromMaterial=kFALSE, isFromMaterialTPC=kFALSE, isFromMaterialITS=kFALSE;

  AliTrackReference *refTPCIn = NULL;
  AliTrackReference *refITS = NULL;

  AliExternalTrackParam *trackInnerC3 = new AliExternalTrackParam(*(track->GetInnerParam()));
  if (!trackInnerC3) return kFALSE; 

  //
  // global track
  //
  Int_t label = TMath::Abs(track->GetLabel()); 
  particle = stack->Particle(label);
  if(particle && particle->GetPDG() && particle->GetPDG()->Charge()!=0.) 
  {
     particleMother = GetMother(particle,stack);
     mech = particle->GetUniqueID();
     isPrim = stack->IsPhysicalPrimary(label);
     isFromStrangess  = IsFromStrangeness(label,stack);
     isFromConversion = IsFromConversion(label,stack);
     isFromMaterial   = IsFromMaterial(label,stack);
  }

  //
  // TPC track
  //
  Int_t labelTPC = TMath::Abs(track->GetTPCLabel()); 
  particleTPC = stack->Particle(labelTPC);
  if(particleTPC && particleTPC->GetPDG() && particleTPC->GetPDG()->Charge()!=0.)
  {
    particleMotherTPC = GetMother(particleTPC,stack);
    mechTPC = particleTPC->GetUniqueID();
    isPrimTPC = stack->IsPhysicalPrimary(labelTPC);
    isFromStrangessTPC  = IsFromStrangeness(labelTPC,stack);
    isFromConversionTPC = IsFromConversion(labelTPC,stack);
    isFromMaterialTPC   = IsFromMaterial(labelTPC,stack);
  }

  //
  // store first track reference
  // for TPC track
  //
  TParticle *part=0;
  TClonesArray *trefs=0;
  Int_t status = mcEvent->GetParticleAndTR(track->GetTPCLabel(), part, trefs);

  if(status>0 && part && trefs && part->GetPDG() && part->GetPDG()->Charge()!=0.) 
  {
    Int_t nTrackRef = trefs->GetEntries();
    //printf("nTrackRef %d \n",nTrackRef);

    Int_t countITS = 0;
    for (Int_t iref = 0; iref < nTrackRef; iref++) 
    {
      AliTrackReference *ref = (AliTrackReference *)trefs->At(iref);
      //printf("ref %p \n",ref);
      //if(ref) printf("ref->DetectorId() %d \n",ref->DetectorId());
      //printf("AliTrackReference::kTPC  %d \n",AliTrackReference::kTPC);

      // Ref. in the middle ITS 
      if(ref && ref->DetectorId()==AliTrackReference::kITS)
      {
        if(!refITS && countITS==2) {
          refITS = ref;
         //printf("refITS %p \n",refITS);
        }
        countITS++;
      }

      // TPC
      if(ref && ref->DetectorId()==AliTrackReference::kTPC)
      {
        if(!refTPCIn) {
          refTPCIn = ref;
          //printf("refTPCIn %p \n",refTPCIn);
          //break;
        }
      }
    }

    // transform inner params to TrackRef
    // reference frame
    if(refTPCIn && trackInnerC3) {
      Double_t fRefPhi = TMath::ATan2(refTPCIn->Y(),refTPCIn->X());
      Bool_t isOK = kFALSE;
      isOK = trackInnerC3->Rotate(fRefPhi);
      isOK = AliTracker::PropagateTrackToBxByBz(trackInnerC3,refTPCIn->R(),track->GetMass(),kStep,kFALSE);
      if(!isOK){
        delete trackInnerC3;
        return kFALSE;
      }
    }
  }

  //
  // ITS track
  //
  Int_t labelITS = TMath::Abs(track->GetITSLabel()); 
  particleITS = stack->Particle(labelITS);
  if(particleITS && particleITS->GetPDG() && particleITS->GetPDG()->Charge()!=0.)
  {
    particleMotherITS = GetMother(particleITS,stack);
    mechITS = particleITS->GetUniqueID();
    isPrimITS = stack->IsPhysicalPrimary(labelITS);
    isFromStrangessITS  = IsFromStrangeness(labelITS,stack);
    isFromConversionITS = IsFromConversion(labelITS,stack);
    isFromMaterialITS   = IsFromMaterial(labelITS,stack);
  }

  // dump to tree
  if(fTreeSRedirector) 
  {
    (*fTreeSRedirector)<<"dNdPtTree"<<
    "extInnerParamRef.="<<trackInnerC3<<
    "refTPCIn.="<<refTPCIn<<
    "refITS.="<<refITS<<
    "particle.="<<particle<<
    "particleMother.="<<particleMother<<
    "mech="<<mech<<
    "isPrim="<<isPrim<<
    "isFromStrangess="<<isFromStrangess<<
    "isFromConversion="<<isFromConversion<<
    "isFromMaterial="<<isFromMaterial<<
    "particleTPC.="<<particleTPC<<
    "particleMotherTPC.="<<particleMotherTPC<<
    "mechTPC="<<mechTPC<<
    "isPrimTPC="<<isPrimTPC<<
    "isFromStrangessTPC="<<isFromStrangessTPC<<
    "isFromConversionTPC="<<isFromConversionTPC<<
    "isFromMaterialTPC="<<isFromMaterialTPC<<
    "particleITS.="<<particleITS<<
    "particleMotherITS.="<<particleMotherITS<<
    "mechITS="<<mechITS<<
    "isPrimITS="<<isPrimITS<<
    "isFromStrangessITS="<<isFromStrangessITS<<
    "isFromConversionITS="<<isFromConversionITS<<
    "isFromMaterialITS="<<isFromMaterialITS<<
    "\n";
  }

    if(trackInnerC3) delete trackInnerC3;

return kTRUE;
}
*/

//_____________________________________________________________________________
/*
Bool_t AlidNdPtTrackDumpTask::DumpEventInfo() 
{
//
// Dump run and event information to tree
// 
      if(fTreeSRedirector) 
      {
        (*fTreeSRedirector)<<"dNdPtTree"<<
        "\n";
      }

return kTRUE;
}
*/

//_____________________________________________________________________________
TParticle *AlidNdPtTrackDumpTask::GetMother(TParticle *const particle, AliStack *const stack) 
{
  if(!particle) return NULL;
  if(!stack) return NULL;

  Int_t motherLabel = TMath::Abs(particle->GetMother(0));  
  TParticle* mother = NULL; 
  mother = stack->Particle(motherLabel); 

return mother;
}

//_____________________________________________________________________________
Bool_t AlidNdPtTrackDumpTask::IsFromConversion(const Int_t label, AliStack *const stack) 
{
  Bool_t isFromConversion = kFALSE;

  if(stack) {
    TParticle* particle = stack->Particle(label);

    if(particle && particle->GetPDG() && particle->GetPDG()->Charge()!=0) 
    {
       Int_t mech = particle->GetUniqueID(); // production mechanism 
       Bool_t isPrim = stack->IsPhysicalPrimary(label);

       Int_t motherLabel = TMath::Abs(particle->GetMother(0));  
       TParticle* mother = stack->Particle(motherLabel); 
       if(mother) {
          Int_t motherPdg = mother->GetPdgCode();

          if(!isPrim && mech==5 && motherPdg==kGamma) { 
            isFromConversion=kTRUE; 
          }
       }
    } 
  } 

  return isFromConversion;
}

//_____________________________________________________________________________
Bool_t AlidNdPtTrackDumpTask::IsFromMaterial(const Int_t label, AliStack *const stack) 
{
  Bool_t isFromMaterial = kFALSE;

  if(stack) {
    TParticle* particle = stack->Particle(label);

    if(particle && particle->GetPDG() && particle->GetPDG()->Charge()!=0) 
    {
       Int_t mech = particle->GetUniqueID(); // production mechanism 
       Bool_t isPrim = stack->IsPhysicalPrimary(label);

       Int_t motherLabel = TMath::Abs(particle->GetMother(0));  
       TParticle* mother = stack->Particle(motherLabel); 
       if(mother) {
          if(!isPrim && mech==13) { 
            isFromMaterial=kTRUE; 
          }
       }
     } 
  } 

  return isFromMaterial;
}

//_____________________________________________________________________________
Bool_t AlidNdPtTrackDumpTask::IsFromStrangeness(const Int_t label, AliStack *const stack) 
{
  Bool_t isFromStrangeness = kFALSE;

  if(stack) {
    TParticle* particle = stack->Particle(label);

    if(particle && particle->GetPDG() && particle->GetPDG()->Charge()!=0) 
    {
       Int_t mech = particle->GetUniqueID(); // production mechanism 
       Bool_t isPrim = stack->IsPhysicalPrimary(label);

       Int_t motherLabel = TMath::Abs(particle->GetMother(0));  
       TParticle* mother = stack->Particle(motherLabel); 
       if(mother) {
          Int_t motherPdg = mother->GetPdgCode();

          // K+-, lambda, antilambda, K0s decays
          if(!isPrim && mech==4 && 
	      (TMath::Abs(motherPdg)==kKPlus || TMath::Abs(motherPdg)==kLambda0 || motherPdg==kK0Short))
          {
            isFromStrangeness = kTRUE;
          } 
       }
    } 
  } 

  return isFromStrangeness;
}


//_____________________________________________________________________________
void AlidNdPtTrackDumpTask::FinishTaskOutput() 
{
  //
  // Called one at the end 
  // locally on working node
  //

  // must be deleted to store trees
  if(fTreeSRedirector)  delete fTreeSRedirector; fTreeSRedirector=0;

  // open temporary file and copy trees to the ouptut container

  TChain* chain = 0;
  TTree* tree1 = 0;
  TTree* tree2 = 0;
  TTree* tree3 = 0;
  TTree* tree4 = 0;
  TTree* tree5 = 0;
  //
  chain = new TChain("dNdPtTree");
  if(chain) { 
    chain->Add("jotwinow_Temp_Trees.root");
    tree1 = chain->CopyTree("1");
    delete chain; chain=0; 
  }
  if(tree1) tree1->Print();

  //
  chain = new TChain("V0s");
  if(chain) { 
    chain->Add("jotwinow_Temp_Trees.root");
    tree2 = chain->CopyTree("1");
    delete chain; chain=0; 
  }
  if(tree2) tree2->Print();

  //
  chain = new TChain("dEdx");
  if(chain) { 
    chain->Add("jotwinow_Temp_Trees.root");
    tree3 = chain->CopyTree("1");
    delete chain; chain=0; 
  }
  if(tree3) tree3->Print();

  //
  chain = new TChain("Laser");
  if(chain) { 
    chain->Add("jotwinow_Temp_Trees.root");
    tree4 = chain->CopyTree("1");
    delete chain; chain=0; 
  }
  if(tree4) tree4->Print();

  //
  chain = new TChain("MCEffTree");
  if(chain) { 
    chain->Add("jotwinow_Temp_Trees.root");
    tree5 = chain->CopyTree("1");
    delete chain; chain=0; 
  }
  if(tree5) tree5->Print();


  OpenFile(1);

  if(tree1) fOutput->Add(tree1);
  if(tree2) fOutput->Add(tree2);
  if(tree3) fOutput->Add(tree3);
  if(tree4) fOutput->Add(tree4);
  if(tree5) fOutput->Add(tree5);
  
  // Post output data.
  PostData(1, fOutput);
}

//_____________________________________________________________________________
void AlidNdPtTrackDumpTask::Terminate(Option_t *) 
{
  // Called one at the end 
  /*
  fOutputSummary = dynamic_cast<TTree*> (GetOutputData(1));
  if(fOutputSummary) delete fOutputSummary; fOutputSummary=0;
  TChain* chain = new TChain("dNdPtTree");
  if(!chain) return;
  chain->Add("jotwinow_HighPt_TrackAndV0_Trees.root");
  TTree *tree = chain->CopyTree("1");
  if (chain) { delete chain; chain=0; }
  if(!tree) return;
  tree->Print();
  fOutputSummary = tree;

  if (!fOutputSummary) {
    Printf("ERROR: AlidNdPtTrackDumpTask::Terminate(): Output data not avaiable %p \n", GetOutputData(1));
    return;
  }
  */

  PostData(1, fOutput);

}
