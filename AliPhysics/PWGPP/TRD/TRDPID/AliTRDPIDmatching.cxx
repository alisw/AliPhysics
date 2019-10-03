/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: Yvonne Pachmayer <pachmay@physi.uni-heidelberg.de>             *
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
//
// The task:
// stores TPC TRD matching quantities in a THnSparse
// output can then be used for matching functions
//
// Author:
// Yvonne Pachmayer <pachmay@physi.uni-heidelberg.de>
//


#include "AliTRDPIDmatching.h"
#include "AliTRDPIDTree.h"
#include "TChain.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliMCEvent.h"
#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliPIDResponse.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDv0KineCuts.h"
#include "THnSparse.h"

using namespace std;

ClassImp(AliTRDPIDmatching)


//________________________________________________________________________
AliTRDPIDmatching::AliTRDPIDmatching(const char *name): AliTRDPIDTree(name),
  fESD(0), fMC(0), fTHntrdmatch(0), fOutputContainer(0), fesdTrackCuts(0), fHasMC(0)
{
    //
    // Constructor
    //

    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());

}

//________________________________________________________________________
AliTRDPIDmatching::~AliTRDPIDmatching()
{
    //
    // dtor
    //

  delete fOutputContainer;
  fOutputContainer = 0;
  
}


//________________________________________________________________________
void AliTRDPIDmatching::UserCreateOutputObjects()
{
    //
    // Output: THnSparse binning and esd track cuts
    //

    // V0 kine cuts
    fV0cuts = new AliESDv0KineCuts();
    // V0 PID Obj arrays
    fV0electrons = new TObjArray;
    fV0pions     = new TObjArray;
    fV0protons   = new TObjArray;

    fesdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "cuts");
    fesdTrackCuts->SetEtaRange(-0.9,+0.9);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetRequireTPCRefit(kTRUE);
    fesdTrackCuts->SetMaxDCAToVertexXY(3);
    fesdTrackCuts->SetMaxDCAToVertexZ(3.0);
    fesdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
    fesdTrackCuts->SetRequireITSRefit(kTRUE);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);

    const Int_t kNdim = 5;
    //                       pid,  eta,  phi,  p,  trdtracklets
    Int_t bins[kNdim] =    {  7,    180,  280,      500,    7};
    Double_t xmin[kNdim] = {  0,   -0.9,    0.,    -10.,    0};
    Double_t xmax[kNdim] = {  7,    0.9,    7.,     10.,    7};

    fTHntrdmatch= new THnSparseF("trdmatching", "TRD matching; pid;eta;phi;p;ntrdtracklets", kNdim, bins, xmin, xmax);

    fOutputContainer = new TObjArray(1);
    fOutputContainer->SetName(GetName());
    fOutputContainer->SetOwner(kTRUE);
    fOutputContainer->Add(fTHntrdmatch);

    PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliTRDPIDmatching::UserExec(Option_t *)
{
    //
    // calls the Process function
    //
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH)
    {
	Printf("ERROR: Could not get ESDInputHandler \n");
    }
    else fESD = (AliESDEvent*)esdH->GetEvent();

    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());

    if(!mcH){
//	Printf("Error: Could not get MCtruthEventHandler \n");
	fHasMC = kFALSE;
    } else {
        fMC=mcH->MCEvent();
	fHasMC = kTRUE;
    }

    FillV0PIDlist();
    ProcessData(fESD);

    PostData(1, fOutputContainer);
    // Clear the V0 PID arrays
    ClearV0PIDlist();
}



//________________________________________________________________________
void AliTRDPIDmatching::ProcessData(AliESDEvent *const esdEvent)
{
    //
    // Process function; called for each event
    //

    AliMCParticle *mctrack = NULL;
    TParticle* mcparticle = NULL;

    if (!esdEvent) {
	Printf("ERROR: esdEvent not available");
	return;
    }

    if(!((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))) Printf("ESD inputhandler not available \n");

    AliESDpid* fESDpid = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();

    if (!fESDpid) {
	Printf("ERROR: fESDpid not available");
	return;
    }

  



    for (Int_t iTracks = 0; iTracks < esdEvent->GetNumberOfTracks(); iTracks++)
    {
	AliESDtrack *trackESD = esdEvent->GetTrack(iTracks);
	if(!trackESD) continue;
	if(!fesdTrackCuts->AcceptTrack(trackESD)){continue;}

	const AliExternalTrackParam *paramIn = trackESD->GetInnerParam();
	Float_t precin=-1;
	if(paramIn) precin=paramIn->Pt();
	else continue;
        if(precin<0.3) continue; // lower momentum cut

	Double_t eta=trackESD->Eta();
	Double_t phi=trackESD->Phi();
	Int_t ntrdtracklets=trackESD->GetTRDntracklets();
	Int_t charge =trackESD->Charge();

	Int_t v0id=0;
	if(fHasMC==kFALSE)
	{
	    Int_t v0tagAllCharges = TMath::Abs(GetV0tag(iTracks));
	    if (v0tagAllCharges == -99) {
		AliError(Form("Problem with V0 tag list (requested V0 track for track %d from %d (list status %d))!", iTracks, esdEvent->GetNumberOfTracks(),
			      fV0tags != 0x0));
		v0id=0;
	    }

	    Bool_t isV0el=0;
	    Bool_t isV0pi=0;
	    Bool_t isV0pr=0;

            if(v0tagAllCharges == 11)   isV0el=1;
            if(v0tagAllCharges == 211)  isV0pi=1;
            if(v0tagAllCharges == 2212) isV0pr=1;

	    if(isV0el) { // electron
		v0id=1;
	    }
	    if(isV0pi) { // pion
		v0id=3;
	    }
	    if(isV0pr) { // proton
		v0id=5;
	    }
	}

	Int_t mcpid=0;
	if(fHasMC==kTRUE)
	{
	    if(!(mctrack = dynamic_cast<AliMCParticle *>(fMC->GetTrack(TMath::Abs(trackESD->GetLabel()))))) continue;
	    mcparticle = mctrack->Particle();
	    if(mcparticle)
	    {
		mcpid = mcparticle->GetPdgCode();
		if(TMath::Abs(mcpid)==11)   v0id=1;
		if(TMath::Abs(mcpid)==13)   v0id=2;
		if(TMath::Abs(mcpid)==211)  v0id=3;
		if(TMath::Abs(mcpid)==321)  v0id=4;
		if(TMath::Abs(mcpid)==2212) v0id=5;
	    }
	}



	Double_t contentSignal[5];
	contentSignal[0]=v0id;
	contentSignal[1]=eta;
	contentSignal[2]=phi;
	contentSignal[3]=precin*charge;
	contentSignal[4]=ntrdtracklets;

	fTHntrdmatch->Fill(contentSignal);

    } //track loop


    PostData(1, fOutputContainer);
}


//________________________________________________________________________
Int_t  AliTRDPIDmatching::CompareFloat(Float_t f1, Float_t f2) const
{
    //compares if the Float_t f1 is equal with f2 and returns 1 if true and 0 if false
    Float_t precision = 0.00001;
    if (((f1 - precision) < f2) &&
	((f1 + precision) > f2))
    {
	return 1;
    }
    else
    {
	return 0;
    }
}

//________________________________________________________________________
void AliTRDPIDmatching::Terminate(const Option_t *)
{
    //
    // Terminate function
    //

}
