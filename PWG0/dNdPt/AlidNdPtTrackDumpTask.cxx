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

#include "TChain.h"
#include "TTreeStream.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TFile.h"
#include "TMatrixD.h"

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
#include "dNdPt/AlidNdPt.h"
#include "dNdPt/AlidNdPtEventCuts.h"
#include "dNdPt/AlidNdPtAcceptanceCuts.h"

#include "dNdPt/AlidNdPtTrackDumpTask.h"

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
  , fdNdPtEventCuts(0)
  , fdNdPtAcceptanceCuts(0)
  , fdNdPtRecAcceptanceCuts(0)
  , fEsdTrackCuts(0)
  , fTrigger(AliTriggerAnalysis::kMB1) 
  , fAnalysisMode(AlidNdPtHelper::kTPC) 
  , fOutputSummary(0)
  , fTreeSRedirector(0)
  , fCentralityEstimator(0)
{
  // Constructor

  // Define input and output slots here
  DefineOutput(0, TTree::Class());
  //DefineOutput(1, TList::Class());

}

//_____________________________________________________________________________
AlidNdPtTrackDumpTask::~AlidNdPtTrackDumpTask()
{
  if(fOutput) delete fOutput;  fOutput =0; 
  //if(fOutputSummary) delete fOutputSummary;  fOutputSummary =0; 
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
  // create output tree
  //
  fTreeSRedirector = new TTreeSRedirector("dNdPtOutliersAnalysisPbPb.root");

  PostData(0, fOutputSummary);
  //PostData(1, fOutput);
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

  fESDfriend = static_cast<AliESDfriend*>(fESD->FindListObject("AliESDfriend"));
  if(!fESDfriend) {
    Printf("ERROR: ESD friends not available");
  }

  //
  Process(fESD,fMC,fESDfriend);

  // Post output data.
  PostData(0, fOutputSummary);
  //PostData(1, fOutput);
}

//_____________________________________________________________________________
void AlidNdPtTrackDumpTask::Process(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent, AliESDfriend *const esdFriend)
{
  //
  // Process real and/or simulated events
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

    for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++)
    {
      AliESDtrack *track = esdEvent->GetTrack(iTrack);
      if(!track) continue;
      if(track->Charge()==0) continue;
      if(!esdTrackCuts->AcceptTrack(track)) continue;
      if(!accCuts->AcceptTrack(track)) continue;

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
      Bool_t isOK = kFALSE;

      //
      // Constrain TPC-only tracks (TPCinner) to vertex
      //
      AliExternalTrackParam * tpcInner = (AliExternalTrackParam *)(track->GetTPCInnerParam());
      if (!tpcInner) continue;
      // transform to the track reference frame 
      isOK = tpcInner->Rotate(track->GetAlpha());
      isOK = tpcInner->PropagateTo(track->GetX(),esdEvent->GetMagneticField());
      if(!isOK) continue;

      // clone TPCinner has to be deleted
      AliExternalTrackParam * tpcInnerC = new AliExternalTrackParam(*(track->GetTPCInnerParam()));
      if (!tpcInnerC) continue;
 
      // constrain TPCinner 
      //isOK = ConstrainTPCInner(tpcInnerC,vtxESD,esdEvent->GetMagneticField());
      isOK = ConstrainTPCInner(tpcInnerC,vtxESD,b);

      // transform to the track reference frame 
      isOK = tpcInnerC->Rotate(track->GetAlpha());
      isOK = tpcInnerC->PropagateTo(track->GetX(),esdEvent->GetMagneticField());

      if(!isOK) {
        if(tpcInnerC) delete tpcInnerC;
	continue;
      }


      //
      // Constrain TPC refitted tracks at inner TPC wall (InnerParams) 
      // to vertex
      //
      // clone track InnerParams has to be deleted
      AliExternalTrackParam * trackInnerC =  new AliExternalTrackParam(*(track->GetInnerParam()));
      if (!trackInnerC) continue;
 
      // constrain track InnerParams 
      isOK = ConstrainTrackInner(trackInnerC,vtxESD,track->GetMass(),b);

      // transform to the track reference frame 
      isOK = trackInnerC->Rotate(track->GetAlpha());
      isOK = trackInnerC->PropagateTo(track->GetX(),esdEvent->GetMagneticField());

      if(!isOK) {
        if(trackInnerC) delete trackInnerC;
	continue;
      }
      
      //
      // calculate chi2 between vi and vj vectors
      // with covi and covj covariance matrices
      // chi2ij = (vi-vj)^(T)*(covi+covj)^(-1)*(vi-vj)
      //
      TMatrixD deltaT(5,1), deltaTtrackC(5,1);
      TMatrixD delta(1,5),  deltatrackC(1,5);
      TMatrixD covarM(5,5), covarMtrackC(5,5);

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
      TMatrixD chi2 = delta*mat2; 
      //chi2.Print();

      // chi2 distance TPC refitted constrained and TPC+ITS
      TMatrixD covarMInvtrackC = covarMtrackC.Invert();
      TMatrixD mat2trackC = covarMInvtrackC*deltaTtrackC;
      TMatrixD chi2trackC = deltatrackC*mat2trackC; 
      //chi2trackC.Print();


      //
      // Propagate ITSout to TPC inner wall 
      // and calculate chi2 distance to track (InnerParams)
      //
      const Double_t kTPCRadius=85; 
      const Double_t kStep=3; 

      // clone track InnerParams has to be deleted
      AliExternalTrackParam *trackInnerC2 = new AliExternalTrackParam(*(track->GetInnerParam()));
      if (!trackInnerC2) continue;
      if(!AliTracker::PropagateTrackToBxByBz(trackInnerC2,kTPCRadius,track->GetMass(),kStep,kFALSE))
      {
	  if(trackInnerC2) delete trackInnerC2;
	  continue;
      }

      AliExternalTrackParam *outerITSc = new AliExternalTrackParam();
      if(!outerITSc) continue;

      TMatrixD chi2OuterITS(1,1);
      chi2OuterITS(0,0) = 0;

      if(esdFriend && esdFriend->TestSkipBit()==kFALSE) 
      {
        // propagate ITSout to TPC inner wall
        AliESDfriendTrack *friendTrack = esdFriend->GetTrack(iTrack);

        if(friendTrack) 
	{
          if( (outerITSc = new AliExternalTrackParam(*friendTrack->GetITSOut())) ) 
	  {
	    if(AliTracker::PropagateTrackToBxByBz(outerITSc,kTPCRadius,track->GetMass(),kStep,kFALSE))
	    {
              // transform ITSout to the track InnerParams reference frame 
	      isOK = kFALSE;
              isOK = outerITSc->Rotate(trackInnerC2->GetAlpha());
              isOK = outerITSc->PropagateTo(trackInnerC2->GetX(),esdEvent->GetMagneticField());

              if(!isOK) {
                if(outerITSc) delete outerITSc;
	        if(trackInnerC2) delete trackInnerC2;
		continue;
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
      AliExternalTrackParam *trackInnerC3 = new AliExternalTrackParam(*(track->GetInnerParam()));
      if(!trackInnerC3) continue;

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
          if(refTPCIn && trackInnerC3) 
	  {
	    Double_t kRefPhi = TMath::ATan2(refTPCIn->Y(),refTPCIn->X());
	    isOK = kFALSE;
            isOK = trackInnerC3->Rotate(kRefPhi);
            isOK = AliTracker::PropagateTrackToBxByBz(trackInnerC3,refTPCIn->R(),track->GetMass(),kStep,kFALSE);

            if(!isOK){
	      if(trackInnerC3) delete trackInnerC3;
	      if(refTPCIn) delete refTPCIn;
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
      }

      //
      Double_t vert[3] = {0}; 
      vert[0] = vtxESD->GetXv();
      vert[1] = vtxESD->GetYv();
      vert[2] = vtxESD->GetZv();
      Int_t mult = vtxESD->GetNContributors();
      Double_t bz = esdEvent->GetMagneticField();
      Double_t runNumber = esdEvent->GetRunNumber();

      //
      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"dNdPtTree"<<
        "runNumber="<<runNumber<<
        "Bz="<<bz<<
        "vertX="<<vert[0]<<
        "vertY="<<vert[1]<<
        "vertZ="<<vert[2]<<
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

	if(tpcInnerC) delete tpcInnerC;
	if(trackInnerC) delete trackInnerC;
	if(trackInnerC2) delete trackInnerC2;
	if(outerITSc) delete outerITSc;

	if(trackInnerC3) delete trackInnerC3;
    }
  }

  PostData(0, fOutputSummary);
  //PostData(1, fOutput);
}


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

  // Post output data.
  PostData(1, fOutput);
  //PostData(0, fOutputSummary);
}

//_____________________________________________________________________________
void AlidNdPtTrackDumpTask::Terminate(Option_t *) 
{
  // Called one at the end 
  
  // check output data
  fOutputSummary = dynamic_cast<TTree*> (GetOutputData(0));
  if(fOutputSummary) delete fOutputSummary; fOutputSummary=0;
  if(fTreeSRedirector)  delete fTreeSRedirector; fTreeSRedirector=0;

  TChain* chain = new TChain("dNdPtTree");
  if(!chain) return;
  chain->Add("dNdPtOutliersAnalysisPbPb.root");
  TTree *tree = chain->CopyTree("1");
  if (chain) { delete chain; chain=0; }
  if(!tree) return;
  tree->Print();
  fOutputSummary = tree;

  if (!fOutputSummary) {
    Printf("ERROR: AlidNdPtTrackDumpTask::Terminate(): Output data not avaiable GetOutputData(0)==0x0 ..." );
    return;
  }
  


  PostData(0, fOutputSummary);
  //PostData(1, fOutput);
}
