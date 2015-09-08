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
 * appeuear in the supporting documentation. The authors make no claims   *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
//
//               Lc->pK0s analysis code
//
//  Input: AOD
//  Output: TTree or THnSparse (mass vs pT vs Centrality)
//
//  Cuts:
//  TTree: very loose cut
//  THnSparse: One THnSparse is created per cut. One cut is specified by
//  an array of bits, each bit corresponds to a cut in "Cut" function.
//  Use "AddCutStream" function to add a cut. 
//
//-------------------------------------------------------------------------
//
//                 Authors: Y.S Watanabe(a)
//  (a) CNS, the University of Tokyo
//  Contatcs: wyosuke@cns.s.u-tokyo.ac.jp
//-------------------------------------------------------------------------

#include <TSystem.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliExternalTrackParam.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSELc2pK0sfromAODtracks.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliTOFPIDResponse.h"
#include "AliAODPidHF.h"
#include "AliInputEventHandler.h"
#include "AliESDtrackCuts.h"
#include "AliNeutralTrackParam.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliExternalTrackParam.h"
#include "AliESDtrack.h"
#include "AliCentrality.h"
#include "AliVertexerTracks.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSELc2pK0sfromAODtracks);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSELc2pK0sfromAODtracks::AliAnalysisTaskSELc2pK0sfromAODtracks() : 
  AliAnalysisTaskSE(),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fListCuts(0),
  fCEvents(0),
  fHTrigger(0),
  fHCentrality(0),
  fAnalCuts(0),
  fIsEventSelected(kFALSE),
  fWriteVariableTree(kFALSE),
  fVariablesTree(0),
  fReconstructPrimVert(kFALSE),
  fIsMB(kFALSE),
  fIsSemi(kFALSE),
  fIsCent(kFALSE),
  fIsINT7(kFALSE),
  fIsEMC7(kFALSE),
  fCandidateVariables(),
  fVtx1(0),
  fV1(0),
  fBzkG(0),
  fCentrality(0),
  fTriggerCheck(0),
  fHistoLcK0SpMass(0),
  fHistoBachPt(0),
  fHistod0Bach(0),
  fHistod0V0(0),
  fHistod0d0(0),
  fHistoV0CosPA(0),
  fHistoProbProton(0),
  fHistoDecayLength(0),
  fHistoK0SMass(0)
{
  //
  // Default Constructor. 
  //
}

//___________________________________________________________________________
AliAnalysisTaskSELc2pK0sfromAODtracks::AliAnalysisTaskSELc2pK0sfromAODtracks(const Char_t* name,
									     AliRDHFCutsLctopK0sfromAODtracks* analCuts, 
									     Bool_t writeVariableTree) :
  AliAnalysisTaskSE(name),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fListCuts(0),
  fCEvents(0),
  fHTrigger(0),
  fHCentrality(0),
  fAnalCuts(analCuts),
  fIsEventSelected(kFALSE),
  fWriteVariableTree(writeVariableTree),
  fVariablesTree(0),
  fReconstructPrimVert(kFALSE),
  fIsMB(kFALSE),
  fIsSemi(kFALSE),
  fIsCent(kFALSE),
  fIsINT7(kFALSE),
  fIsEMC7(kFALSE),
  fCandidateVariables(),
  fVtx1(0),
  fV1(0),
  fBzkG(0),
  fCentrality(0),
  fTriggerCheck(0),
  fHistoLcK0SpMass(0),
  fHistoBachPt(0),
  fHistod0Bach(0),
  fHistod0V0(0),
  fHistod0d0(0),
  fHistoV0CosPA(0),
  fHistoProbProton(0),
  fHistoDecayLength(0),
  fHistoK0SMass(0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSELc2pK0sfromAODtracks","Calling Constructor");

  DefineOutput(1,TList::Class());  //conters
  DefineOutput(2,TList::Class());
  if(writeVariableTree){
    DefineOutput(3,TTree::Class());  //My private output
  }else{
    DefineOutput(3,TList::Class());  //conters
  }
}

//___________________________________________________________________________
AliAnalysisTaskSELc2pK0sfromAODtracks::~AliAnalysisTaskSELc2pK0sfromAODtracks() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskSELc2pK0sfromAODtracks","Calling Destructor");

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }

  if (fOutputAll) {
    delete fOutputAll;
    fOutputAll = 0;
  }

  if (fListCuts) {
    delete fListCuts;
    fListCuts = 0;
  }


  if (fAnalCuts) {
    delete fAnalCuts;
    fAnalCuts = 0;
  }

  if (fVariablesTree) {
    delete fVariablesTree;
    fVariablesTree = 0;
  }

}

//_________________________________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::Init() {
  //
  // Initialization
  //
  //

  fIsEventSelected=kFALSE;

  if (fDebug > 1) AliInfo("Init");

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("ListCuts");
  fListCuts->Add(new AliRDHFCutsLctopK0sfromAODtracks(*fAnalCuts));
  PostData(2,fListCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::UserExec(Option_t *)
{
  //
  // UserExec
  //

  if (!fInputEvent) {
    AliError("NO EVENT FOUND!");
    return;
  }
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

  fCEvents->Fill(1);

  //------------------------------------------------
  // First check if the event has proper vertex and B
  //------------------------------------------------
	
  fVtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if (!fVtx1) return;

  Double_t pos[3],cov[6];
  fVtx1->GetXYZ(pos);
  fVtx1->GetCovarianceMatrix(cov);
  fV1 = new AliESDVertex(pos,cov,100.,100,fVtx1->GetName());

  fBzkG = (Double_t)aodEvent->GetMagneticField(); 
  AliKFParticle::SetField(fBzkG);
  if (TMath::Abs(fBzkG)<0.001) {
    delete fV1;
    return;
  }
  fCEvents->Fill(2);

  //------------------------------------------------
  // Event selection 
  //------------------------------------------------
  Bool_t fIsTriggerNotOK = fAnalCuts->IsEventRejectedDueToTrigger();
  if(!fIsTriggerNotOK) fCEvents->Fill(3);
  fIsEventSelected = fAnalCuts->IsEventSelected(aodEvent); // better to initialize before CheckEventSelection call
  if(!fIsEventSelected) {
    delete fV1;
    return;
  }
  fCEvents->Fill(4);

  fIsMB=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kMB)==(AliVEvent::kMB);
  fIsSemi=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kSemiCentral)==(AliVEvent::kSemiCentral);
  fIsCent=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kCentral)==(AliVEvent::kCentral); 
  fIsINT7=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kINT7)==(AliVEvent::kINT7);  
  fIsEMC7=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kEMC7)==(AliVEvent::kEMC7);   
  fTriggerCheck = fIsMB+2*fIsSemi+4*fIsCent+8*fIsINT7+16*fIsEMC7;
  if(fIsMB) fHTrigger->Fill(1);
  if(fIsSemi) fHTrigger->Fill(2);
  if(fIsCent) fHTrigger->Fill(3);
  if(fIsINT7) fHTrigger->Fill(4);
  if(fIsEMC7) fHTrigger->Fill(5);
  if(fIsMB|fIsSemi|fIsCent) fHTrigger->Fill(7);
  if(fIsINT7|fIsEMC7) fHTrigger->Fill(8);
  if(fIsMB&fIsSemi) fHTrigger->Fill(10);
  if(fIsMB&fIsCent) fHTrigger->Fill(11);
  if(fIsINT7&fIsEMC7) fHTrigger->Fill(12);

  AliCentrality *cent = aodEvent->GetCentrality();
  fCentrality = cent->GetCentralityPercentile("V0M");
  fHCentrality->Fill(fCentrality);

  //------------------------------------------------
  // Check if the event has v0 candidate
  //------------------------------------------------
  //Int_t nv0 = aodEvent->GetNumberOfV0s();
  fCEvents->Fill(5);

  //------------------------------------------------
  // MC analysis setting
  //------------------------------------------------
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader=0;
  if (fUseMCInfo) {
    // MC array need for maching
    mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcArray) {
      AliError("Could not find Monte-Carlo in AOD");
      return;
    }
    fCEvents->Fill(6); // in case of MC events
  
    // load MC header
    mcHeader = (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader) {
      AliError("AliAnalysisTaskSELc2pK0sfromAODtracks::UserExec: MC header branch not found!\n");
      return;
    }
    fCEvents->Fill(7); // in case of MC events
  
    Double_t zMCVertex = mcHeader->GetVtxZ();
    if (TMath::Abs(zMCVertex) > fAnalCuts->GetMaxVtxZ()) {
      AliDebug(2,Form("Event rejected: abs(zVtxMC)=%f > fAnalCuts->GetMaxVtxZ()=%f",zMCVertex,fAnalCuts->GetMaxVtxZ()));
      return;
    } else {
      fCEvents->Fill(17); // in case of MC events
    }
  }

  //------------------------------------------------
  // Main analysis done in this function
  //------------------------------------------------
  MakeAnalysis(aodEvent,mcArray);


  PostData(1,fOutput);
  if(fWriteVariableTree){
    PostData(3,fVariablesTree);
  }else{
    PostData(3,fOutputAll);
  }

  fIsEventSelected=kFALSE;

  delete fV1;
  return;
}

//________________________________________ terminate ___________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::Terminate(Option_t*)
{    
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  //AliInfo("Terminate","");
  AliAnalysisTaskSE::Terminate();
  
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    AliError("fOutput not available");
    return;
  }

  if(!fWriteVariableTree){
    fOutputAll = dynamic_cast<TList*> (GetOutputData(3));
    if (!fOutputAll) {     
      AliError("fOutputAll not available");
      return;
    }
  }

  return;
}

//___________________________________________________________________________
void AliAnalysisTaskSELc2pK0sfromAODtracks::UserCreateOutputObjects() 
{ 
  //
  // UserCreateOutputObject
  //
  //AliInfo(Form("CreateOutputObjects of task %s\n", GetName()));

  //------------------------------------------------
  // output object setting
  //------------------------------------------------
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("chist0");
  DefineGeneralHistograms(); // define general histograms
  PostData(1,fOutput);

  DefineTreeVariables();
  if (fWriteVariableTree) {
    PostData(3,fVariablesTree);
  }else{
    fOutputAll = new TList();
    fOutputAll->SetOwner();
    fOutputAll->SetName("anahisto");
    DefineAnalysisHistograms(); // define general histograms
    PostData(3,fOutputAll);
  }


  return;
}

//-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2pK0sfromAODtracks::MakeAnalysis
(
 AliAODEvent *aodEvent, TClonesArray *mcArray
 )
{
  //
  // Main Analysis part
  //

  Int_t nV0s= aodEvent->GetNumberOfV0s();
  if (nV0s==0) {
    return;
  }
  Int_t nTracks= aodEvent->GetNumberOfTracks();
  if (nTracks==0) {
    return;
  }

  //------------------------------------------------
  // Arrays to store MC matching information
  //------------------------------------------------
	Int_t usedmclab[20];//Used Lc Label: Assuming there are less than 20 Lc/evt
	Int_t nusedmclab[20];//Number of times the Lc label is used: Assuming there are less than 20 Lc/evt
	for(Int_t i=0;i<20;i++) {
		usedmclab[i]=-9999;
		nusedmclab[i]=0;
	}

  //------------------------------------------------
  // V0 loop 
  //------------------------------------------------
  for (Int_t iv0 = 0; iv0<nV0s; iv0++) {
    AliAODv0 *v0 = aodEvent->GetV0(iv0);
    if(!v0) continue;
    if(!fAnalCuts->SingleV0Cuts(v0,fVtx1)) continue;

    AliAODTrack *cptrack =  (AliAODTrack*)(v0->GetDaughter(0));
    AliAODTrack *cntrack =  (AliAODTrack*)(v0->GetDaughter(1));

    //------------------------------------------------
    // track loop 
    //------------------------------------------------
    for (Int_t itrk = 0; itrk<nTracks; itrk++) {
      AliAODTrack *trk = (AliAODTrack*)aodEvent->GetTrack(itrk);
      if(trk->GetID()<0) continue;
      if(!fAnalCuts->SingleTrkCuts(trk)) continue;

      Int_t cpid = cptrack->GetID();
      Int_t cnid = cntrack->GetID();
      Int_t lpid = trk->GetID();
      if((cpid==lpid)||(cnid==lpid)) continue;

      if(!fAnalCuts->SelectWithRoughCuts(v0,trk)) continue;

      AliAODVertex *secVert = ReconstructSecondaryVertex(v0,trk,aodEvent);
      if(!secVert) continue;

      AliAODRecoCascadeHF *lcobj = MakeCascadeHF(v0,trk,aodEvent,secVert);
      if(!lcobj) {
	continue;
      }

      AliAODMCParticle *mclc = 0;
      AliAODMCParticle *mcproton = 0;
      AliAODMCParticle *mck0s = 0;
      Int_t mclablc = 0;
      Int_t nmclablc = 0;
      if(fUseMCInfo)
      {
        Int_t pdgDg[2]={2212,310};
        Int_t pdgDgv0[2]={211,211};
        mclablc = lcobj->MatchToMC(4122,pdgDg[1],pdgDg,pdgDgv0,mcArray,kTRUE);
        if(mclablc>-1){
          mclc = (AliAODMCParticle*) mcArray->At(mclablc);
          for(Int_t ia=0;ia<20;ia++){
            if(usedmclab[ia]==mclablc){
              nusedmclab[ia]++;
              nmclablc = nusedmclab[ia];
              break;
            }
            if(usedmclab[ia]==-9999){
              usedmclab[ia]=mclablc;
              nusedmclab[ia]++;
              nmclablc = nusedmclab[ia];
              break;
            }
          }
          Int_t mcprotonlabel = mclc->GetDaughter(0);
          if(mcprotonlabel>=0)
            mcproton=(AliAODMCParticle*) mcArray->At(mcprotonlabel);
          Int_t mck0slabel = mclc->GetDaughter(1);
          if(mck0slabel>=0)
            mck0s=(AliAODMCParticle*) mcArray->At(mck0slabel);
        }
      }

      FillROOTObjects(lcobj,mclc,mcproton,mck0s,nmclablc);

      lcobj->GetSecondaryVtx()->RemoveDaughters();
      lcobj->UnsetOwnPrimaryVtx();
      delete lcobj;lcobj=NULL;
      delete secVert;
    }
  }

}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2pK0sfromAODtracks::FillROOTObjects(AliAODRecoCascadeHF *lcobj, AliAODMCParticle *mcpart, AliAODMCParticle *mcproton, AliAODMCParticle *mck0s, Int_t mcnused) 
{
  //
  // Fill histograms or tree depending on fWriteVariableTree 
  //

  AliAODTrack *trk = lcobj->GetBachelor();
  AliAODv0 *v0 = lcobj->Getv0();

  fCandidateVariables[ 0] = lcobj->InvMassLctoK0sP();
  fCandidateVariables[ 1] = lcobj->Px();
  fCandidateVariables[ 2] = lcobj->Py();
  fCandidateVariables[ 3] = lcobj->Pz();
  fCandidateVariables[ 4] = v0->MassK0Short();
  fCandidateVariables[ 5] = lcobj->PxProng(0);
  fCandidateVariables[ 6] = lcobj->PyProng(0);
  fCandidateVariables[ 7] = lcobj->PzProng(0);
  fCandidateVariables[ 8] = lcobj->PxProng(1);
  fCandidateVariables[ 9] = lcobj->PyProng(1);
  fCandidateVariables[10] = lcobj->PzProng(1);
  fCandidateVariables[11] = fVtx1->GetX();
  fCandidateVariables[12] = fVtx1->GetY();
  fCandidateVariables[13] = fVtx1->GetZ();
  fCandidateVariables[14] = fCentrality;
  fCandidateVariables[15] = lcobj->DecayLengthXY();
  fCandidateVariables[16] = (Float_t) fAnalCuts->CalculateLcCosPAXY(lcobj);

  Double_t nSigmaTPCpr=-9999.;
  Double_t nSigmaTOFpr=-9999.;
  Double_t probProton=-9999.;
  if(fAnalCuts->GetIsUsePID())
    {
			nSigmaTPCpr = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kProton);
			nSigmaTOFpr = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kProton);
			if(fAnalCuts->GetPidHF()->GetUseCombined()){
				probProton = fAnalCuts->GetProtonProbabilityTPCTOF(trk);
			}
      fCandidateVariables[17] = nSigmaTPCpr;
      fCandidateVariables[18] = nSigmaTOFpr;
      fCandidateVariables[19] = probProton;
    }

	fCandidateVariables[20] = -9999;
	fCandidateVariables[21] = -9999;
	fCandidateVariables[22] = -9999;
	fCandidateVariables[23] = -9999;
	fCandidateVariables[24] = -9999;
	fCandidateVariables[25] = -9999;
	fCandidateVariables[26] = -9999;
	if(fUseMCInfo){
		if(mcpart){
			fCandidateVariables[20] = mcpart->Label();
			fCandidateVariables[21] = mcnused;
			fCandidateVariables[22] = mcpart->GetPdgCode();
			Double_t mcprimvertx = mcpart->Xv();
			Double_t mcprimverty = mcpart->Yv();
			Double_t mcsecvertx = mcproton->Xv();
			Double_t mcsecverty = mcproton->Yv();
			Double_t recosecvertx = lcobj->GetSecondaryVtx()->GetX();
			Double_t recosecverty = lcobj->GetSecondaryVtx()->GetY();
			fCandidateVariables[23] = TMath::Sqrt((mcsecvertx-mcprimvertx)*(mcsecvertx-mcprimvertx)+(mcsecverty-mcprimverty)*(mcsecverty-mcprimverty));
			fCandidateVariables[24] = TMath::Sqrt((recosecvertx-mcprimvertx)*(recosecvertx-mcprimvertx)+(recosecverty-mcprimverty)*(recosecverty-mcprimverty));
			Double_t vecx_vert = recosecvertx-mcprimvertx;
			Double_t vecy_vert = recosecverty-mcprimverty;
			Double_t vecl_vert = TMath::Sqrt(vecx_vert*vecx_vert+vecy_vert*vecy_vert);
			Double_t vecx_mom = lcobj->Px();
			Double_t vecy_mom = lcobj->Py();
			Double_t vecl_mom = lcobj->Pt();
			if(vecl_vert>0.&&vecl_mom>0.)
				fCandidateVariables[25] = (vecx_vert*vecx_mom+vecy_vert*vecy_mom)/vecl_vert/vecl_mom;
			fCandidateVariables[26] = mcpart->Pt();
		}
	}

  if(fWriteVariableTree)
    fVariablesTree->Fill();
  else{
	if(fAnalCuts->IsSelected(lcobj,AliRDHFCuts::kCandidate))
	  {
	    Double_t cont[3];
	    cont[0] = lcobj->InvMassLctoK0sP();
	    cont[1] = lcobj->Pt();
	    cont[2] = fCentrality;
	    fHistoLcK0SpMass->Fill(cont);

      fHistoBachPt->Fill(trk->Pt());
      fHistod0Bach->Fill(lcobj->Getd0Prong(0));
      fHistod0V0->Fill(lcobj->Getd0Prong(1));
      fHistod0d0->Fill(lcobj->Getd0Prong(0)*lcobj->Getd0Prong(1));
      fHistoV0CosPA->Fill(lcobj->CosV0PointingAngle());
      fHistoProbProton->Fill(probProton);
      fHistoDecayLength->Fill(lcobj->DecayLengthXY()*(fAnalCuts->CalculateLcCosPAXY(lcobj)));
      fHistoK0SMass->Fill(v0->MassK0Short());
	  }
  }
  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2pK0sfromAODtracks::DefineTreeVariables() 
{
  //
  // Define tree variables
  //

  const char* nameoutput = GetOutputSlot(3)->GetContainer()->GetName();
  fVariablesTree = new TTree(nameoutput,"Candidates variables tree");
  Int_t nVar = 27;
  fCandidateVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="InvMassLc2pK0s";
  fCandidateVariableNames[ 1]="LcPx";
  fCandidateVariableNames[ 2]="LcPy";
  fCandidateVariableNames[ 3]="LcPz";
  fCandidateVariableNames[ 4]="massK0Short";
  fCandidateVariableNames[ 5]="V0Px";
  fCandidateVariableNames[ 6]="V0Py";
  fCandidateVariableNames[ 7]="V0Pz";
  fCandidateVariableNames[ 8]="BachPx";
  fCandidateVariableNames[ 9]="BachPy";
  fCandidateVariableNames[10]="BachPz";
  fCandidateVariableNames[11]="PrimVertx";
  fCandidateVariableNames[12]="PrimVerty";
  fCandidateVariableNames[13]="PrimVertz";
  fCandidateVariableNames[14]="Centrality";
  fCandidateVariableNames[15]="DecayLengthXY";
  fCandidateVariableNames[16]="LcCosPAXY";
  fCandidateVariableNames[17]="nSigmaTPCpr";
  fCandidateVariableNames[18]="nSigmaTOFpr";
  fCandidateVariableNames[19]="probProton";
  fCandidateVariableNames[20]="mclcID";
  fCandidateVariableNames[21]="mcnused";
  fCandidateVariableNames[22]="mcpdgcode";
  fCandidateVariableNames[23]="mcdecaylength";
  fCandidateVariableNames[24]="mcdecaylength_secsmear";
  fCandidateVariableNames[25]="mclccospaxy";
  fCandidateVariableNames[26]="mclcpt";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

////__________________________________________________________________________
void  AliAnalysisTaskSELc2pK0sfromAODtracks::DefineGeneralHistograms() {
  //
  // This is to define general histograms
  //

  fCEvents = new TH1F("fCEvents","conter",18,-0.5,17.5);
  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetBinLabel(1,"X1");
  fCEvents->GetXaxis()->SetBinLabel(2,"Analyzed events");
  fCEvents->GetXaxis()->SetBinLabel(3,"AliAODVertex exists");
  fCEvents->GetXaxis()->SetBinLabel(4,"TriggerOK");
  fCEvents->GetXaxis()->SetBinLabel(5,"IsEventSelected");
  fCEvents->GetXaxis()->SetBinLabel(6,"CascadesHF exists");
  fCEvents->GetXaxis()->SetBinLabel(7,"MCarray exists");
  fCEvents->GetXaxis()->SetBinLabel(8,"MCheader exists");
  fCEvents->GetXaxis()->SetBinLabel(9,"triggerClass!=CINT1");
  fCEvents->GetXaxis()->SetBinLabel(10,"triggerMask!=kAnyINT");
  fCEvents->GetXaxis()->SetBinLabel(11,"triggerMask!=kAny");
  fCEvents->GetXaxis()->SetBinLabel(12,"vtxTitle.Contains(Z)");
  fCEvents->GetXaxis()->SetBinLabel(13,"vtxTitle.Contains(3D)");
  fCEvents->GetXaxis()->SetBinLabel(14,"vtxTitle.Doesn'tContain(Z-3D)");
  fCEvents->GetXaxis()->SetBinLabel(15,Form("zVtx<=%2.0fcm",fAnalCuts->GetMaxVtxZ()));
  fCEvents->GetXaxis()->SetBinLabel(16,"!IsEventSelected");
  fCEvents->GetXaxis()->SetBinLabel(17,"triggerMask!=kAnyINT || triggerClass!=CINT1");
  fCEvents->GetXaxis()->SetBinLabel(18,Form("zVtxMC<=%2.0fcm",fAnalCuts->GetMaxVtxZ()));
  //fCEvents->GetXaxis()->SetTitle("");
  fCEvents->GetYaxis()->SetTitle("counts");

  fHTrigger = new TH1F("fHTrigger","counter",18,-0.5,17.5);
  fHTrigger->SetStats(kTRUE);
  fHTrigger->GetXaxis()->SetBinLabel(1,"X1");
  fHTrigger->GetXaxis()->SetBinLabel(2,"kMB");
  fHTrigger->GetXaxis()->SetBinLabel(3,"kSemiCentral");
  fHTrigger->GetXaxis()->SetBinLabel(4,"kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(5,"kINT7");
  fHTrigger->GetXaxis()->SetBinLabel(6,"kEMC7");
  //fHTrigger->GetXaxis()->SetBinLabel(7,"Space");
  fHTrigger->GetXaxis()->SetBinLabel(8,"kMB|kSemiCentral|kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(9,"kINT7|kEMC7");
  fHTrigger->GetXaxis()->SetBinLabel(11,"kMB&kSemiCentral");
  fHTrigger->GetXaxis()->SetBinLabel(12,"kMB&kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(13,"kINT7&kEMC7");

  fHCentrality = new TH1F("fHCentrality","conter",100,0.,100.);

  fOutput->Add(fCEvents);
  fOutput->Add(fHTrigger);
  fOutput->Add(fHCentrality);

  return;
}
//__________________________________________________________________________
void  AliAnalysisTaskSELc2pK0sfromAODtracks::DefineAnalysisHistograms() 
{
  //
  // Define analyis histograms
  //
	
  //------------------------------------------------
  // Basic histogram
  //------------------------------------------------
  Int_t bins_base[3]=		{80				,20		,10};
  Double_t xmin_base[3]={2.286-0.2,0		,0.00};
  Double_t xmax_base[3]={2.286+0.2,20.	,100};
  fHistoLcK0SpMass = new THnSparseF("fHistoLcK0SpMass","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoLcK0SpMass);


  //------------------------------------------------
  // checking histograms
  //------------------------------------------------
  fHistoBachPt = new TH1F("fHistoBachPt","Bachelor p_{T}",100,0.0,5.0);
  fOutputAll->Add(fHistoBachPt);
  fHistod0Bach = new TH1F("fHistod0Bach","Bachelor d_{0}",100,-0.5,0.5);
  fOutputAll->Add(fHistod0Bach);
  fHistod0V0 = new TH1F("fHistod0V0","V_{0} d_{0}",100,-0.5,0.5);
  fOutputAll->Add(fHistod0V0);
  fHistod0d0 = new TH1F("fHistod0d0","d_{0} d_{0}",100,-0.5,0.5);
  fOutputAll->Add(fHistod0d0);
  fHistoV0CosPA=new TH1F("fHistoV0CosPA","V0->Second vertex cospa",100,-1.,1.0);
  fOutputAll->Add(fHistoV0CosPA);
  fHistoProbProton=new TH1F("fHistoProbProton","ProbProton",100,0.,1.0);
  fOutputAll->Add(fHistoProbProton);
  fHistoDecayLength=new TH1F("fHistoDecayLength","Decay Length",100,-0.1,0.1);
  fOutputAll->Add(fHistoDecayLength);
  fHistoK0SMass=new TH1F("fHistoK0SMass","K0S mass",100,0.497-0.05,0.497+0.05);
  fOutputAll->Add(fHistoK0SMass);

  return;
}

//________________________________________________________________________
AliAODRecoCascadeHF* AliAnalysisTaskSELc2pK0sfromAODtracks::MakeCascadeHF(AliAODv0 *v0, AliAODTrack *part, AliAODEvent * aod, AliAODVertex *secVert) 
{
  //
  // Create AliAODRecoCascadeHF object from the argument
  //

  if(!v0) return 0x0;
  if(!part) return 0x0;
  if(!aod) return 0x0;

  //------------------------------------------------
  // PrimaryVertex
  //------------------------------------------------
  AliAODVertex *primVertexAOD;
  Bool_t unsetvtx = kFALSE;
  if(fReconstructPrimVert){
    primVertexAOD = CallPrimaryVertex(v0,part,aod);
    if(!primVertexAOD){
      primVertexAOD = fVtx1;
    }else{
      unsetvtx = kTRUE;
    }
  }else{
    primVertexAOD = fVtx1;
  }
  if(!primVertexAOD) return 0x0;
  Double_t posprim[3]; primVertexAOD->GetXYZ(posprim);

  //------------------------------------------------
  // DCA between tracks
  //------------------------------------------------
  AliESDtrack *esdtrack = new AliESDtrack((AliVTrack*)part);

  AliNeutralTrackParam *trackV0=NULL;
  const AliVTrack *trackVV0 = dynamic_cast<const AliVTrack*>(v0);
  if(trackVV0)  trackV0 = new AliNeutralTrackParam(trackVV0);

  Double_t xdummy, ydummy;
  Double_t dca = esdtrack->GetDCA(trackV0,fBzkG,xdummy,ydummy);


  //------------------------------------------------
  // Propagate all tracks to the secondary vertex and calculate momentum there
  //------------------------------------------------
	
  Double_t d0z0bach[2],covd0z0bach[3];
  if(sqrt(pow(secVert->GetX(),2)+pow(secVert->GetY(),2))<1.){
    part->PropagateToDCA(secVert,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
    trackV0->PropagateToDCA(secVert,fBzkG,kVeryBig);
  }else{
    part->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
    trackV0->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig);
  }
  Double_t momv0_new[3]={-9999,-9999,-9999.};
  trackV0->GetPxPyPz(momv0_new);

  Double_t px[2],py[2],pz[2];
  px[0] = part->Px(); py[0] = part->Py(); pz[0] = part->Pz(); 
  px[1] = momv0_new[0]; py[1] = momv0_new[1]; pz[1] = momv0_new[2]; 

  //------------------------------------------------
  // d0
  //------------------------------------------------
  Double_t d0[3],d0err[3];

  part->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
  d0[0]= d0z0bach[0];
  d0err[0] = TMath::Sqrt(covd0z0bach[0]);

  Double_t d0z0v0[2],covd0z0v0[3];
  trackV0->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0v0,covd0z0v0);
  d0[1]= d0z0v0[0];
  d0err[1] = TMath::Sqrt(covd0z0v0[0]);

  //------------------------------------------------
  // Create AliAODRecoCascadeHF
  //------------------------------------------------
  Short_t charge = part->Charge();
  AliAODRecoCascadeHF *theCascade = new AliAODRecoCascadeHF(secVert,charge,px,py,pz,d0,d0err,dca);
  if(!theCascade)  
    {
      if(unsetvtx) delete primVertexAOD; primVertexAOD=NULL;
      if(esdtrack) delete esdtrack;
      if(trackV0) delete trackV0;
      return 0x0;
    }
  theCascade->SetOwnPrimaryVtx(primVertexAOD);
  UShort_t id[2]={(UShort_t)part->GetID(),(UShort_t)trackV0->GetID()};
  theCascade->SetProngIDs(2,id);

  theCascade->GetSecondaryVtx()->AddDaughter(part);
  theCascade->GetSecondaryVtx()->AddDaughter(v0);

  if(unsetvtx) delete primVertexAOD; primVertexAOD=NULL;
  if(esdtrack) delete esdtrack;
  if(trackV0) delete trackV0;

  return theCascade;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSELc2pK0sfromAODtracks::CallPrimaryVertex(AliAODv0 *v0, AliAODTrack *trk, AliAODEvent* aod)
{
  //
  // Make an array of tracks which should not be used in primary vertex calculation and 
  // Call PrimaryVertex function
  //

  TObjArray *TrackArray = new TObjArray(3);
  
  AliESDtrack *cptrk1 = new AliESDtrack((AliVTrack*)trk);
  TrackArray->AddAt(cptrk1,0);
  
  AliESDtrack *cascptrack = new AliESDtrack((AliVTrack*)v0->GetDaughter(0));
  TrackArray->AddAt(cascptrack,1);
  AliESDtrack *cascntrack = new AliESDtrack((AliVTrack*)v0->GetDaughter(1));
  TrackArray->AddAt(cascntrack,2);
  
  AliAODVertex *newvert  = PrimaryVertex(TrackArray,aod);
  
  for(Int_t i=0;i<3;i++)
    {
      AliESDtrack *tesd = (AliESDtrack*)TrackArray->UncheckedAt(i);
      delete tesd;
    }
  TrackArray->Clear();
  delete TrackArray;
  
  return newvert;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSELc2pK0sfromAODtracks::PrimaryVertex(const TObjArray *trkArray,
								   AliVEvent *event)
{
  //
  //Used only for pp
  //copied from AliAnalysisVertexingHF (except for the following 3 lines)
  //

  Bool_t fRecoPrimVtxSkippingTrks = kTRUE;
  Bool_t fRmTrksFromPrimVtx = kFALSE;

  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;
  
  //vertexESD = new AliESDVertex(*fV1);
  

  if(!fRecoPrimVtxSkippingTrks && !fRmTrksFromPrimVtx) { 
    // primary vertex from the input event
    
    vertexESD = new AliESDVertex(*fV1);
    
  } else {
    // primary vertex specific to this candidate
    
    Int_t nTrks = trkArray->GetEntriesFast();
    AliVertexerTracks *vertexer = new AliVertexerTracks(event->GetMagneticField());
    
    if(fRecoPrimVtxSkippingTrks) { 
      // recalculating the vertex
      
      if(strstr(fV1->GetTitle(),"VertexerTracksWithConstraint")) {
	Float_t diamondcovxy[3];
	event->GetDiamondCovXY(diamondcovxy);
	Double_t pos[3]={event->GetDiamondX(),event->GetDiamondY(),0.};
	Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
	AliESDVertex *diamond = new AliESDVertex(pos,cov,1.,1);
	vertexer->SetVtxStart(diamond);
	delete diamond; diamond=NULL;
	if(strstr(fV1->GetTitle(),"VertexerTracksWithConstraintOnlyFitter")) 
	  vertexer->SetOnlyFitter();
      }
      Int_t skipped[1000];
      Int_t nTrksToSkip=0,id;
      AliExternalTrackParam *t = 0;
      for(Int_t i=0; i<nTrks; i++) {
	t = (AliExternalTrackParam*)trkArray->UncheckedAt(i);
	id = (Int_t)t->GetID();
	if(id<0) continue;
	skipped[nTrksToSkip++] = id;
      }
      // TEMPORARY FIX
      // For AOD, skip also tracks without covariance matrix
      Double_t covtest[21];
      for(Int_t j=0; j<event->GetNumberOfTracks(); j++) {
	AliVTrack *vtrack = (AliVTrack*)event->GetTrack(j);
	if(!vtrack->GetCovarianceXYZPxPyPz(covtest)) {
	  id = (Int_t)vtrack->GetID();
	  if(id<0) continue;
	  skipped[nTrksToSkip++] = id;
	}
      }
      for(Int_t ijk=nTrksToSkip; ijk<1000; ijk++) skipped[ijk]=-1;
      //
      vertexer->SetSkipTracks(nTrksToSkip,skipped);
      vertexESD = (AliESDVertex*)vertexer->FindPrimaryVertex(event); 
      
    } else if(fRmTrksFromPrimVtx && nTrks>0) { 
      // removing the prongs tracks
      
      TObjArray rmArray(nTrks);
      UShort_t *rmId = new UShort_t[nTrks];
      AliESDtrack *esdTrack = 0;
      AliESDtrack *t = 0;
      for(Int_t i=0; i<nTrks; i++) {
	t = (AliESDtrack*)trkArray->UncheckedAt(i);
	esdTrack = new AliESDtrack(*t);
	rmArray.AddLast(esdTrack);
	if(esdTrack->GetID()>=0) {
	  rmId[i]=(UShort_t)esdTrack->GetID();
	} else {
	  rmId[i]=9999;
	}
      }
      Float_t diamondxy[2]={static_cast<Float_t>(event->GetDiamondX()),static_cast<Float_t>(event->GetDiamondY())};
      vertexESD = vertexer->RemoveTracksFromVertex(fV1,&rmArray,rmId,diamondxy);
      delete [] rmId; rmId=NULL;
      rmArray.Delete();
      
    }
    
    delete vertexer; vertexer=NULL;
    if(!vertexESD) return vertexAOD;
    if(vertexESD->GetNContributors()<=0) { 
      //AliDebug(2,"vertexing failed"); 
      delete vertexESD; vertexESD=NULL;
      return vertexAOD;
    }
    
    
  }
  
  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  delete vertexESD; vertexESD=NULL;
  
  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF);
  
  return vertexAOD;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSELc2pK0sfromAODtracks::ReconstructSecondaryVertex(AliAODv0 *v0, AliAODTrack *part, AliAODEvent * aod) 
{
  //
  // Reconstruct secondary vertex from trkArray (Copied from AliAnalysisVertexingHF)
  //
  //------------------------------------------------
  // PrimaryVertex
  //------------------------------------------------
  AliAODVertex *primVertexAOD;
  Bool_t unsetvtx = kFALSE;
  if(fReconstructPrimVert){
    primVertexAOD = CallPrimaryVertex(v0,part,aod);
    if(!primVertexAOD){
      primVertexAOD = fVtx1;
    }else{
      unsetvtx = kTRUE;
    }
  }else{
    primVertexAOD = fVtx1;
  }
  if(!primVertexAOD) return 0x0;

  //------------------------------------------------
  // Secondary vertex
  //------------------------------------------------

  Double_t LcPx = part->Px()+v0->Px();
  Double_t LcPy = part->Py()+v0->Py();
  Double_t LcPt = TMath::Sqrt(LcPx*LcPx+LcPy*LcPy);

  Double_t d0z0[2],covd0z0[3];
  part->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  Double_t x0 = primVertexAOD->GetX();
  Double_t y0 = primVertexAOD->GetY();
  Double_t px0 = LcPx/LcPt;
  Double_t py0 = LcPy/LcPt;
  Double_t tx[3];
  part->GetXYZ(tx);
  Double_t x1 = tx[0];
  Double_t y1 = tx[1];
  part->GetPxPyPz(tx);
  Double_t px1 = tx[0];
  Double_t py1 = tx[1];
  Double_t pt1 = sqrt(px1*px1+py1*py1);
  px1 = px1/pt1;
  py1 = py1/pt1;

  Double_t dx = x0 - x1;
  Double_t dy = y0 - y1;

  Double_t Delta = -px0*py1+py0*px1;
  Double_t a0 = -9999.;
  if(Delta!=0)
    {
      a0 = 1./Delta * (py1 * dx - px1 * dy);
    }
  Double_t neovertx = x0 + a0 * px0;
  Double_t neoverty = y0 + a0 * py0;
  Double_t z0 = primVertexAOD->GetZ();
  Double_t neovertz = z0 + TMath::Abs(a0)*part->Pz()/part->Pt();

  if(unsetvtx) delete primVertexAOD; primVertexAOD=NULL;

  Double_t pos[3],cov[6],chi2perNDF;
  pos[0]=neovertx;
  pos[1]=neoverty;
  pos[2]=neovertz;
  cov[0]=0.0;
  cov[1]=0.0;
  cov[2]=0.0;
  cov[3]=0.0;
  cov[4]=0.0;
  cov[5]=0.0;
  chi2perNDF=0.0;
  AliAODVertex *secVert = new AliAODVertex(pos,cov,chi2perNDF);
  if(!secVert){
    return 0x0;
  }
  return secVert;
}
