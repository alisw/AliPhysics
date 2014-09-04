/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------
//         AliAnalysisTaskContMC class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliVTrack.h"
#include "AliAODMCParticle.h"
#include "AliVParticle.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskContMC.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliHelperPID.h"
#include "AliCentrality.h"
#include "TProof.h"
#include "AliPID.h"
#include "AliVEvent.h"
#include "AliPIDResponse.h"
#include "AliStack.h"
#include <TMCProcess.h>

#include <iostream>

using namespace std;

ClassImp(AliAnalysisTaskContMC) 

//________________________________________________________________________
AliAnalysisTaskContMC::AliAnalysisTaskContMC(const char *name) : AliAnalysisTaskSE(name), fAOD(0), fNSigmaPID(0), fIsMC(0), fOutput(0), fHistID(0)
{
  // Default constructor
  

  DefineInput(0, TChain::Class());
  //DefineOutput(1, AliHelperPID::Class());
  DefineOutput(1, TList::Class());
  
}
//________________________________________________________________________
//________________________________________________________________________
void AliAnalysisTaskContMC::UserCreateOutputObjects()
{
  Printf("\n\n\n\n\n\n In CreateOutput Object:");
  
  //create output objects
  Printf("NSigma: %.1f",fNSigmaPID->GetNSigmaCut());
  Printf("IsMC: %d",fNSigmaPID->GetisMC());
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("list");
  
  fHistID = new TH3F("fHistID","fHistID",6,-1.5,4.5,4,-1.5,2.5,38,0.2,4);
  fOutput->Add(fHistID);
  
  if(!fNSigmaPID)AliFatal("PID object should be set in the steering macro");
  fOutput->Add(fNSigmaPID);
  
  //PostData(1, fNSigmaPID  );
  PostData(1, fOutput  );
  
}

//________________________________________________________________________
void AliAnalysisTaskContMC::UserExec(Option_t *)
{
  const Int_t npart=4;
  const Int_t pdgcode[npart+2]={-1,211,321,2212,11,13};
  const Int_t partid[npart+2]={-1,0,1,2,3,4};
  // main event loop
  fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
  if (!fAOD) {
    AliWarning("ERROR: AliAODEvent not available \n");
    return;
  }
  
  if (strcmp(fAOD->ClassName(), "AliAODEvent"))
    {
      AliFatal("Not processing AODs");
    }
  //MC Loop
  TClonesArray *arrayMC = 0;
  Printf("fIsMC: %d",fIsMC);
  if (fIsMC)
    {
      arrayMC = (TClonesArray*) fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
      if (!arrayMC) {
	AliFatal("Error: MC particles branch not found!\n");
      }
    }
  Double_t centrality = 0;

  AliCentrality *centralityObj = ((AliVAODHeader*)fAOD->GetHeader())->GetCentralityP();
  if (centralityObj)
    {
      centrality = centralityObj->GetCentralityPercentileUnchecked("V0M");
      AliInfo(Form("Centrality is %f", centrality));
    }
  else
    {
      Printf("WARNING: Centrality object is 0");
      centrality = -1;
    }
  if (centrality < 0)
    return;
  

  //vertex selection
  Int_t nVertex = ((AliAODEvent*)fAOD)->GetNumberOfVertices();
  if( nVertex > 0 ) { 
    AliAODVertex* vertex = (AliAODVertex*)((AliAODEvent*)fAOD)->GetPrimaryVertex();
    //Int_t nTracksPrim = vertex->GetNContributors();
    Double_t zVertex = vertex->GetZ();
    //10 cm cut
    if(TMath::Abs(zVertex)>10) return;
    //AliInfo(Form(" Vertex in = %f with %d particles by  %s data ...",zVertex,nTracksPrim,vertex->GetName()));
    // Reject TPC only vertex
    TString name(vertex->GetName());
    if (name.CompareTo("PrimaryVertex") && name.CompareTo("SPDVertex"))return;
  }  

  //Int_t count=0;  
  //track loop
  for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));
    if(!track) AliFatal("Not a standard AOD");
    if(!(track->TestFilterBit(32)))continue;
    if (TMath::Abs(track->Eta()) > .8 || track->Pt() < .2) continue;
    
    Int_t pdg=-999;
    Int_t isph=-999;
    if (fIsMC && arrayMC) {
      AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(TMath::Abs(track->GetLabel()));
      if (!partMC) { 
	AliError("Cannot get MC particle");
	continue; 
      }
      pdg=TMath::Abs(partMC->GetPdgCode());
      isph=partMC->IsPhysicalPrimary();
    }
    if(!isph)continue;
    //step 1, TOF Matching
    UInt_t status;
    status=track->GetStatus();
    if((status&AliVTrack::kTOFout)==0 || (status&AliVTrack::kTIME)==0)continue;
    
    //step 2, combined PID
    Int_t IDTPCTOF=fNSigmaPID->GetParticleSpecies(track,1);
    if(IDTPCTOF==999)IDTPCTOF=-1;
    Int_t IDMC=-1;
    for(Int_t ipart=0;ipart<npart+2;ipart++)if(TMath::Abs(pdg)==pdgcode[ipart])IDMC=partid[ipart];  
    fHistID->Fill(IDMC,IDTPCTOF,track->Pt());
  } // end loop on tracks
  
  PostData(1,fOutput);
  
}

//_________________________________________________________________
void   AliAnalysisTaskContMC::Terminate(Option_t *)
{
  // Terminate analysis
  //
  fOutput = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  } 
  fHistID = dynamic_cast<TH3F*>(fOutput->FindObject("fHistID"));
  fNSigmaPID = dynamic_cast<AliHelperPID*>(fOutput->FindObject("fNSigmaPID"));
}
