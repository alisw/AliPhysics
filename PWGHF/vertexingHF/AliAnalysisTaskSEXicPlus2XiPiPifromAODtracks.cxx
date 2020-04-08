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
//               XicPlus2XiPiPi analysis code
//
//  Input: AOD
//  Output: TTree or THnSparse (mass vs pT vs Centrality)
//
//  Cuts:
//  TTree: very loose cut
//  THnSparse: One THnSparse is created per cut. One cut is specified by
//  an array of bits, each bit corresponds to a cut in "Cut" function.
//
//-------------------------------------------------------------------------
//
//                 Authors: Y.S Watanabe(a), G. Luparello (b)
//                
//  (a) CNS, the University of Tokyo
//  (b) INFN, Trieste
//
//  Contatcs: wyosuke@cns.s.u-tokyo.ac.jp, gluparel@cern.ch
//-------------------------------------------------------------------------

#include <TSystem.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <THnSparse.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
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
#include "AliAODRecoCascadeHF3Prong.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliAODPidHF.h"
#include "AliInputEventHandler.h"
#include "AliESDtrackCuts.h"
#include "AliNeutralTrackParam.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliExternalTrackParam.h"
#include "AliCentrality.h"
#include "AliVertexerTracks.h"
#include "AliNormalizationCounter.h"
#include "AliVertexingHFUtils.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks() : 
  AliAnalysisTaskSE(),
  fUseMCInfo(kFALSE),
  fFillSignalOnly(kFALSE),
  fFillBkgOnly(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fListCuts(0),
  fCEvents(0),
  fHTrigger(0),
  fHCentrality(0),
  fAnalCuts(0),
  fIsEventSelected(kFALSE),
  fWriteVariableTree(kFALSE),
  fFillSparse(kFALSE),
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
  //fTriggerCheck(0),
  fHistoXicMass(0x0),
  fSparseXicMass(0x0),
  fHistoMCSpectrumAccXic(0),
  fHistoDcaPi1Pi2(0),
  fHistoDcaPi1Casc(0),
  fHistoDcaPi2Casc(0),
  fHistoLikeDecayLength(0),
  fHistoLikeDecayLengthXY(0),
  fHistoXicCosPA(0),
  fHistoXiMass(0),
  fHistoCascDcaXiDaughters(0),
  fHistoCascDcaV0Daughters(0),
  fHistoCascDcaV0ToPrimVertex(0),
  fHistoCascDcaPosToPrimVertex(0),
  fHistoCascDcaNegToPrimVertex(0),
  fHistoCascDcaBachToPrimVertex(0),
  fHistoCascCosPAXiPrim(0),
  fHistoXiPt(0),
  fHistoPiPt(0),
  fHistoPid0(0),
  fHistonSigmaTPCpi(0),
  fHistonSigmaTOFpi(0),
  fHistoProbPion(0),
  fHistoXiMassvsPtRef1(0),
  fHistoXiMassvsPtRef2(0),
  fHistoXiMassvsPtRef3(0),
  fHistoXiMassvsPtRef4(0),
  fHistoXiMassvsPtRef5(0),
  fHistoXiMassvsPtRef6(0),
  fHistoPiPtRef(0),
  fHistoPiEtaRef(0),
  fQAHistoNSelectedTracks(0),
  fQAHistoNSelectedCasc(0),
  fQAHistoDCApi1pi2(0),
  fQAHistoAODPrimVertX(0),
  fQAHistoAODPrimVertY(0),
  fQAHistoAODPrimVertZ(0),
  fQAHistoRecoPrimVertX(0),
  fQAHistoRecoPrimVertY(0),
  fQAHistoRecoPrimVertZ(0),
  fQAHistoSecondaryVertexX(0),
  fQAHistoSecondaryVertexY(0),
  fQAHistoSecondaryVertexZ(0),
  fQAHistoSecondaryVertexXY(0),
  fCounter(0)
{
  //
  // Default Constructor. 
  //
}

//___________________________________________________________________________
AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks(const Char_t* name,
											   AliRDHFCutsXicPlustoXiPiPifromAODtracks* analCuts, 
											   Bool_t writeVariableTree, Bool_t fillSparse) :
  AliAnalysisTaskSE(name),
  fUseMCInfo(kFALSE),
  fFillSignalOnly(kFALSE),
  fFillBkgOnly(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fListCuts(0),
  fCEvents(0),
  fHTrigger(0),
  fHCentrality(0),
  fAnalCuts(analCuts),
  fIsEventSelected(kFALSE),
  fWriteVariableTree(writeVariableTree),
  fFillSparse(fillSparse),
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
  //fTriggerCheck(0),
  fHistoXicMass(0x0),
  fSparseXicMass(0x0),
  fHistoMCSpectrumAccXic(0),
  fHistoDcaPi1Pi2(0),
  fHistoDcaPi1Casc(0),
  fHistoDcaPi2Casc(0),
  fHistoLikeDecayLength(0),
  fHistoLikeDecayLengthXY(0),
  fHistoXicCosPA(0),
  fHistoXiMass(0),
  fHistoCascDcaXiDaughters(0),
  fHistoCascDcaV0Daughters(0),
  fHistoCascDcaV0ToPrimVertex(0),
  fHistoCascDcaPosToPrimVertex(0),
  fHistoCascDcaNegToPrimVertex(0),
  fHistoCascDcaBachToPrimVertex(0),
  fHistoCascCosPAXiPrim(0),
  fHistoXiPt(0),
  fHistoPiPt(0),
  fHistoPid0(0),
  fHistonSigmaTPCpi(0),
  fHistonSigmaTOFpi(0),
  fHistoProbPion(0),
  fHistoXiMassvsPtRef1(0),
  fHistoXiMassvsPtRef2(0),
  fHistoXiMassvsPtRef3(0),
  fHistoXiMassvsPtRef4(0),
  fHistoXiMassvsPtRef5(0),
  fHistoXiMassvsPtRef6(0),
  fHistoPiPtRef(0),
  fHistoPiEtaRef(0),
  fQAHistoNSelectedTracks(0),
  fQAHistoNSelectedCasc(0),
  fQAHistoDCApi1pi2(0),
  fQAHistoAODPrimVertX(0),
  fQAHistoAODPrimVertY(0),
  fQAHistoAODPrimVertZ(0),
  fQAHistoRecoPrimVertX(0),
  fQAHistoRecoPrimVertY(0),
  fQAHistoRecoPrimVertZ(0),
  fQAHistoSecondaryVertexX(0),
  fQAHistoSecondaryVertexY(0),
  fQAHistoSecondaryVertexZ(0),
  fQAHistoSecondaryVertexXY(0),
  fCounter(0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks","Calling Constructor");
  
  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());
  if(writeVariableTree){
    DefineOutput(3,TTree::Class());
  }else{
    DefineOutput(3,TList::Class());
  }
  DefineOutput(4,AliNormalizationCounter::Class());
}

//___________________________________________________________________________
AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::~AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks","Calling Destructor");
  
  
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
  
  if(fCandidateVariables){
    delete fCandidateVariables;
    fCandidateVariables = 0;
  }

  if(fHistoXicMass) delete fHistoXicMass;
  if(fSparseXicMass) delete fSparseXicMass;

  if(fHistoMCSpectrumAccXic) delete fHistoMCSpectrumAccXic;
  
  if(fCounter){
    delete fCounter;
    fCounter=0;
  }
}

//_________________________________________________
void AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::Init() {
  //
  // Initialization
  //
  //
  
  //Copied from $ALICE_ROOT/PWGHF/vertexingHF/ConfigVertexingHF.C
  
  fIsEventSelected=kFALSE;
  
  if (fDebug > 1) AliInfo("Init");
  
  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("ListCuts");
  fListCuts->Add(new AliRDHFCutsXicPlustoXiPiPifromAODtracks(*fAnalCuts));
  PostData(2,fListCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::UserExec(Option_t *)
{
  //
  // UserExec code
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
  fQAHistoAODPrimVertX->Fill(pos[0]);
  fQAHistoAODPrimVertY->Fill(pos[1]);
  fQAHistoAODPrimVertZ->Fill(pos[2]);
  
  fV1 = new AliESDVertex(pos,cov,100.,100,fVtx1->GetName());
  
  fBzkG = (Double_t)aodEvent->GetMagneticField(); 
  AliKFParticle::SetField(fBzkG);
  if (TMath::Abs(fBzkG)<0.001) {
    delete fV1;
    return;
  }
  fCEvents->Fill(2);

  fCounter->StoreEvent(aodEvent, fAnalCuts, fUseMCInfo);
  
  //------------------------------------------------
  // Event selection 
  //------------------------------------------------
  Bool_t fIsTriggerNotOK = fAnalCuts->IsEventRejectedDueToTrigger();
  if(!fIsTriggerNotOK) fCEvents->Fill(3);
  
  fIsEventSelected = fAnalCuts->IsEventSelected(aodEvent); 
  if(!fIsEventSelected) {
    //cout<<"Why: "<<fAnalCuts->GetWhyRejection()<<endl;
    delete fV1;
    return;
  }
  
  //cout<<fabs(aodEvent->GetPrimaryVertex()->GetZ()-aodEvent->GetPrimaryVertexSPD()->GetZ())<<endl;
  
  fCEvents->Fill(4);
  
  fIsMB=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kMB)==(AliVEvent::kMB);
  fIsSemi=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kSemiCentral)==(AliVEvent::kSemiCentral);
  fIsCent=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kCentral)==(AliVEvent::kCentral); 
  fIsINT7=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kINT7)==(AliVEvent::kINT7);  
  fIsEMC7=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kEMC7)==(AliVEvent::kEMC7);
  //this trigger check is not used: commentig out for the moment everywhere in the task
  //  fTriggerCheck = fIsMB+2*fIsSemi+4*fIsCent+8*fIsINT7+16*fIsEMC7;
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
  // MC analysis setting
  //------------------------------------------------
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader=0;
  if (fUseMCInfo) {
    // MC array need for maching
    mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcArray) {
      AliError("Could not find Monte-Carlo in AOD");
      delete fV1;
      return;
    }
    fCEvents->Fill(6); // in case of MC events
    
    // load MC header
    mcHeader = (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader) {
      AliError("AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::UserExec: MC header branch not found!\n");
      delete fV1;
      return;
    }
    fCEvents->Fill(7); // in case of MC events
    
    Double_t zMCVertex = mcHeader->GetVtxZ();
    if (TMath::Abs(zMCVertex) > fAnalCuts->GetMaxVtxZ()) {
      AliDebug(2,Form("Event rejected: abs(zVtxMC)=%f > fAnalCuts->GetMaxVtxZ()=%f",zMCVertex,fAnalCuts->GetMaxVtxZ()));
      delete fV1;
      return;
    } else {

      //DOUBT: isEventSelected before or after LoopOverGenParticles?
      
      LoopOverGenParticles(mcArray);
      fCEvents->Fill(17); // in case of MC events
    }
  }
  
  //------------------------------------------------
  // Check if the event has cascade candidate
  //------------------------------------------------
  Int_t ncasc = aodEvent->GetNumberOfCascades();
  Int_t nselecasc = 0.;
  for(Int_t ic=0;ic<ncasc;ic++){
    AliAODcascade *casc = aodEvent->GetCascade(ic);
    if(!fAnalCuts) continue;
    if(fAnalCuts->SingleCascadeCuts(casc,pos )) nselecasc++; 
  }
  
  if(nselecasc==0){
    delete fV1;
    return;
  }
  
  fCEvents->Fill(5); //counter of the events with selected cascades
  
  //------------------------------------------------
  // Main analysis done in this function
  //------------------------------------------------
  MakeAnalysis(aodEvent, mcArray);
  
  PostData(1,fOutput);
  if(fWriteVariableTree){
    PostData(3,fVariablesTree);
  }else{
    PostData(3,fOutputAll);
  }
  PostData(4,fCounter);
  fIsEventSelected=kFALSE;
  
  delete fV1;
  return;
}

//________________________________________ terminate ___________________________
void AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::Terminate(Option_t*)
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
void AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::UserCreateOutputObjects() 
{ 
  //
  // UserCreateOutputObjects
  //

  AliInfo(Form("CreateOutputObjects of task %s\n", GetName()));
  
  //------------------------------------------------
  // output object setting
  //------------------------------------------------
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("chist0");
  DefineGeneralHistograms(); // define general histograms
  PostData(1,fOutput);
  

  if (fWriteVariableTree) {
    DefineTreeVariables();
    PostData(3,fVariablesTree);
  }else{
    fOutputAll = new TList();
    fOutputAll->SetOwner();
    fOutputAll->SetName("anahisto");
    DefineAnalysisHistograms(); // define general histograms
    PostData(3,fOutputAll);
  }

  fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(4)->GetContainer()->GetName()));
  fCounter->Init();
  PostData(4,fCounter);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::MakeAnalysis
(
 AliAODEvent *aodEvent, TClonesArray *mcArray
 )
{
  //
  // Main analysis part called from UserExec
  //

  Int_t nCascades= aodEvent->GetNumberOfCascades();
  if (nCascades==0) {
    return;
  }
  Int_t nTracks= aodEvent->GetNumberOfTracks();
  if (nTracks==0) {
    return;
  }
  
  //------------------------------------------------
  // Select good track before hand to save time
  //------------------------------------------------
  Bool_t  seleTrkFlags[nTracks];
  Int_t nSeleTrks=0;
  SelectTrack(aodEvent,nTracks,nSeleTrks,seleTrkFlags); //select candidates pions
  fQAHistoNSelectedTracks->Fill(nSeleTrks);
  
  Bool_t  seleCascFlags[nCascades];
  Int_t     nSeleCasc=0;
  SelectCascade(aodEvent,nCascades,nSeleCasc,seleCascFlags);
  fQAHistoNSelectedCasc->Fill(nSeleCasc);
 
  Int_t usedmclab[20];//Used Xic Label: Assuming there are less than 20 Xic/evt
  Int_t nusedmclab[20];//Number of times the Xic label is used: Assuming there are less than 20 Xic/evt
  for(Int_t i=0;i<20;i++) {
    usedmclab[i]=-9999;
    nusedmclab[i]=0;
  }
  
  //------------------------------------------------
  // Cascade loop 
  //------------------------------------------------
  for (Int_t icasc = 0; icasc<nCascades; icasc++) {
    if(!seleCascFlags[icasc]) continue;
    fCEvents->Fill(18);
    
    AliAODcascade *casc = aodEvent->GetCascade(icasc);
    if(!casc) continue;

    
    AliAODTrack *cptrack =  (AliAODTrack*)(casc->GetDaughter(0));
    AliAODTrack *cntrack =  (AliAODTrack*)(casc->GetDaughter(1));
    AliAODTrack *cbtrack =  (AliAODTrack*)(casc->GetDecayVertexXi()->GetDaughter(0));
    if(!cptrack || !cntrack || !cbtrack) continue;
    
    Int_t cpid = cptrack->GetID();
    Int_t cnid = cntrack->GetID();
    Int_t cbid = cbtrack->GetID();

    if(cptrack->Charge()==0) continue;
    if(cntrack->Charge()==0) continue;
    if(cbtrack->Charge()==0) continue;

    Short_t charge_casc = cptrack->Charge() + cntrack->Charge() + cbtrack->Charge();
    
    //------------------------------------------------
    // Track1 loop 
    //------------------------------------------------
    for (Int_t itrk1 = 0; itrk1<nTracks-1; itrk1++) {
      if(!seleTrkFlags[itrk1]) continue;
      AliAODTrack *trk1 = (AliAODTrack*)aodEvent->GetTrack(itrk1);
      if(!trk1||trk1->GetID()<0) continue;
      
      //------------------------------------------------
      // Track2 loop 
      //------------------------------------------------
      for (Int_t itrk2 = itrk1+1; itrk2<nTracks; itrk2++) {
	if(!seleTrkFlags[itrk2]) continue;
	AliAODTrack *trk2 = (AliAODTrack*)aodEvent->GetTrack(itrk2);
	if(!trk2||trk2->GetID()<0) continue;
	
	if(!SelectLikeSign(trk1,trk2)) continue;
	fCEvents->Fill(19);
	
	Int_t lpid1 = trk1->GetID();
	Int_t lpid2 = trk2->GetID();
	if((cpid==lpid1)||(cpid==lpid2)||(cnid==lpid1)||(cnid==lpid2)||(cbid==lpid1)||(cbid==lpid2)) continue;
	
	
	Short_t charge_like1 = trk1->Charge();
	Short_t charge_like2 = trk2->Charge();
	Bool_t ok_charge = kFALSE;
	if((charge_casc==-1)&&(charge_like1==1)&&(charge_like2==1)) ok_charge = kTRUE;
	if((charge_casc==1)&&(charge_like1==-1)&&(charge_like2==-1)) ok_charge = kTRUE;
	if(!ok_charge) continue;
	fCEvents->Fill(20);
	//------------------------------------------------
	// Roughly select good candidates to speed up
	//------------------------------------------------
	if(!fAnalCuts->SelectWithRoughCuts(casc,trk1,trk2)) continue;
	fCEvents->Fill(21);
	//------------------------------------------------
	// Secondary vertex calculation
	//------------------------------------------------
	Double_t dispersion;
	AliAODVertex *secVert = CallReconstructSecondaryVertex(trk1,trk2,dispersion);
	if(!secVert) continue;
	fCEvents->Fill(22);
	fQAHistoSecondaryVertexX->Fill(secVert->GetX());
	fQAHistoSecondaryVertexY->Fill(secVert->GetY());
	fQAHistoSecondaryVertexZ->Fill(secVert->GetZ());
	fQAHistoSecondaryVertexXY->Fill(secVert->GetX()*secVert->GetX()+secVert->GetY()*secVert->GetY());
	
	AliAODRecoCascadeHF3Prong *xicobj = MakeCascadeHF3Prong(casc,trk1,trk2,aodEvent,secVert,dispersion);
	if(!xicobj) {
	  delete secVert;
	  continue;
	}
	fCEvents->Fill(23);
	
	AliAODMCParticle *mcxic = 0;
	AliAODMCParticle *mcdaughter1 = 0;
	AliAODMCParticle *mcdaughter2 = 0;
	AliAODMCParticle *mcdaughterxi = 0;
	AliAODMCParticle *mcdaughterLambda = 0;
	AliAODMCParticle *mcdaughterPionFromLambda =0;
	AliAODMCParticle *mcdaughterProtonFromLambda=0;
	AliAODMCParticle *mcdaughterPionFromXi=0;
	Int_t nFound=0;
	
	Int_t mclabxic = 0;
	Int_t nmclabxic = 0;
	Bool_t isXic = kFALSE;
	Int_t checkOrigin=-9999;
	
	if(fUseMCInfo)
	{
	  //Note: this can be enclosed in a MatchRecoCandToMC() function
	  //Input will be the AliAODRecoCascadeHF3Prong and mcArray and return the AODMCParticle
	  
	  Int_t pdgDg[3]={211,3312,211};
	  Int_t pdgDgcasc[2]={211,3122};
	  Int_t pdgDgv0[2]={2212,211};
	  mclabxic = xicobj->MatchToMC(4232,pdgDg[1],pdgDg,pdgDgcasc,pdgDgv0,mcArray);
	  if(mclabxic>-1){
	    mcxic = (AliAODMCParticle*) mcArray->At(mclabxic);
	    if (mcxic){
	      checkOrigin=AliVertexingHFUtils::CheckOrigin(mcArray,mcxic,kTRUE);
	      Bool_t isInAcc=kTRUE;
	      if(fAnalCuts){
		if(!fAnalCuts->IsInFiducialAcceptance(mcxic->Pt(),mcxic->Y())) isInAcc=kFALSE;
	      }else{
		if(TMath::Abs(mcxic->Y())>0.8) isInAcc=kFALSE;
	      }
	      
	      for(Int_t ia=0;ia<20;ia++){   //this is to check how many Xic there are per event.. I guess always only one. I am not sure that this is really needed
		if(usedmclab[ia]==mclabxic){
		  nusedmclab[ia]++;
		  nmclabxic = nusedmclab[ia];
		  break;
		}
		if(usedmclab[ia]==-9999){
		  usedmclab[ia]=mclabxic;
		  nusedmclab[ia]++;
		  nmclabxic = nusedmclab[ia];
		  break;
		}
	      }
	      Int_t pi_counter = 0;
	      for(Int_t idau=mcxic->GetDaughterFirst();idau<mcxic->GetDaughterLast()+1;idau++)
		{
		  //cout<<idau<<endl;
		  if(idau<0) break;
		  AliAODMCParticle *mcpart = (AliAODMCParticle*) mcArray->At(idau);
		  Int_t pdgcode = TMath::Abs(mcpart->GetPdgCode());
		  if(pdgcode==211 && pi_counter==0){
		    mcdaughter1 = mcpart;
		    pi_counter=1;
		    nFound++;
		  }else if(pdgcode==211 && pi_counter==1){
		    mcdaughter2 = mcpart;
		    pi_counter=2;
		    nFound++;
		  }else if(pdgcode==3312){
		    mcdaughterxi = mcpart;
		    Int_t nXiDau= mcpart->GetNDaughters();
		    if(nXiDau!=2) break; //no correct Xi decay
		    Int_t indFirstXiDau=mcpart->GetDaughterLabel(0);
		    for(Int_t XiDau=0; XiDau<2; XiDau++){
		      Int_t indXiDau=indFirstXiDau+XiDau;
		      if(indXiDau<0) break;		      
		      AliAODMCParticle* Xidau=dynamic_cast<AliAODMCParticle*>(mcArray->At(indXiDau));
		      if(!Xidau) break;
		      Int_t pdgXidau=Xidau->GetPdgCode();
		      if(TMath::Abs(pdgXidau)==3122){
			Int_t nLambdaDau=Xidau->GetNDaughters();
			if (nLambdaDau!=2) break;
			Int_t indFirstLambdaDau=Xidau->GetDaughterLabel(0);
			for(Int_t LambdaDau=0; LambdaDau<2; LambdaDau++){
			  Int_t indLambdaDau=indFirstLambdaDau+LambdaDau;
			  AliAODMCParticle* Lambdadau=dynamic_cast<AliAODMCParticle*>(mcArray->At(indLambdaDau));
			  if(!Lambdadau) break;
			  Int_t pdgLambdadau=Lambdadau->GetPdgCode();
			  if(TMath::Abs(pdgLambdadau)==2212){
			    mcdaughterProtonFromLambda=Lambdadau;
			    nFound++;
			  } else if(TMath::Abs(pdgLambdadau)==211){
			    mcdaughterPionFromLambda=Lambdadau;
			    nFound++;
			  }
			  mcdaughterLambda=Xidau;
			}
		      }	else if(TMath::Abs(pdgXidau)==211) {
			mcdaughterPionFromXi=Xidau;
			nFound++;
		      }
		    }
		  }
		} //end of loop over xic daughters
	      if(nFound==5){
		//this is actually already done at previous steps SelectSingleTrk and SelectCascades -> this is not needed here
		//if(isInAcc){
		//		  if ((TMath::Abs(mcdaughter1->Eta())>0.8) || (TMath::Abs(mcdaughter2->Eta())>0.8) || (TMath::Abs(mcdaughterProtonFromLambda->Eta())>0.9) || (TMath::Abs(mcdaughterPionFromLambda->Eta())>0.9) || (TMath::Abs(mcdaughterPionFromXi->Eta())>0.9)){ //to check also the y acceptance of Xi and Lambda?
		  //  isInAcc=kFALSE;
		  //}
		//}
		if(isInAcc){
		  fHistoMCSpectrumAccXic->Fill(mcxic->Pt(),kReco,checkOrigin);
		  isXic=kTRUE;
		  if(fAnalCuts->IsSelected(xicobj,AliRDHFCuts::kCandidate)){
		    fHistoMCSpectrumAccXic->Fill(mcxic->Pt(),kRecoPID,checkOrigin);
		  }
		  fAnalCuts->SetUsePID(kFALSE);
		  if(fAnalCuts->IsSelected(xicobj,AliRDHFCuts::kCandidate)){
		    fHistoMCSpectrumAccXic->Fill(mcxic->Pt(),kRecoCuts,checkOrigin);
		  }
		  fAnalCuts->SetUsePID(kTRUE);
		} else {
		  fHistoMCSpectrumAccXic->Fill(mcxic->Pt(),kReco08,checkOrigin);
		  if(fAnalCuts->IsSelected(xicobj,AliRDHFCuts::kCandidate)){
		    fHistoMCSpectrumAccXic->Fill(mcxic->Pt(),kRecoPID08,checkOrigin);
		  }
		  fAnalCuts->SetUsePID(kFALSE);
		  if(fAnalCuts->IsSelected(xicobj,AliRDHFCuts::kCandidate)){
		    fHistoMCSpectrumAccXic->Fill(mcxic->Pt(),kRecoCuts08,checkOrigin);
		  }
		  fAnalCuts->SetUsePID(kTRUE);
		}
	      }
	    }
	    //cout<<"pointer: "<<mcdaughter1<<" "<<mcdaughter2<<" "<<mcdaughterxi<<endl;
	    //cout<<"1: "<<mcdaughter1->GetPdgCode()<<endl;
	    //cout<<"2: "<<mcdaughter2->GetPdgCode()<<endl;
	    //cout<<"3: "<<mcdaughterxi->GetPdgCode()<<endl; 
	  }
	} //if useMCinfo
       
	//if(!fAnalCuts->IsInFiducialAcceptance(xicobj->Pt(),xicobj->Y(4232))) continue;
	FillROOTObjects(xicobj,mcxic,mcdaughter1,mcdaughter2,mcdaughterxi,nmclabxic,isXic,checkOrigin); //AliAODRecoCascadeHF3Prong,
	
	xicobj->GetSecondaryVtx()->RemoveDaughters();
	xicobj->UnsetOwnPrimaryVtx();
	delete xicobj;xicobj=NULL;
	delete secVert;
	
      }//trk2
    }//trk1
  }//casc
}

//________________________________________________________________________
void AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::FillROOTObjects(AliAODRecoCascadeHF3Prong *xicobj, AliAODMCParticle *mcpart, AliAODMCParticle *mcdaughter1, AliAODMCParticle *mcdaughter2, AliAODMCParticle *mcdaughterxi, Int_t mcnused, Bool_t isXic, Int_t checkOrigin) 
{
  //
  // Fill histogram or Tree depending on fWriteVariableTree flag
  //

  AliAODTrack *part1 = xicobj->GetBachelor1();
  AliAODTrack *part2 = xicobj->GetBachelor2();
  AliAODcascade *casc = xicobj->GetCascade();

  Double_t dca[3];
  xicobj->GetDCAs(dca);

  Double_t nSigmaTPCpi1=-9999.;
  Double_t nSigmaTPCpi2=-9999.;
  Double_t nSigmaTOFpi1=-9999.;
  Double_t nSigmaTOFpi2=-9999.;
  Double_t probPion1=-9999.;
  Double_t probPion2=-9999.;
      
  if(fAnalCuts->GetIsUsePID()) {
    nSigmaTPCpi1 = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(part1,AliPID::kPion);    
    nSigmaTPCpi2 = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(part2,AliPID::kPion);
    nSigmaTOFpi1 = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(part1,AliPID::kPion);      
    nSigmaTOFpi2 = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(part2,AliPID::kPion);
	
    probPion1 =  fAnalCuts->GetPionProbabilityTPCTOF(part1);
    probPion2 =  fAnalCuts->GetPionProbabilityTPCTOF(part2);
    
  }

      
  if(fWriteVariableTree){
    if ((fUseMCInfo && fFillSignalOnly && !fFillBkgOnly && !isXic) || (fUseMCInfo && !fFillSignalOnly && fFillBkgOnly && isXic)) return;
    else {
	  
      //  if (!fUseMCInfo || (fUseMCInfo && fFillSignalOnly && isXiC) || (fUseMCInfo && !fFillSignalOnly)){
      fCandidateVariables[ 0] = xicobj->InvMassPiXiPi();
      fCandidateVariables[ 1] = xicobj->Pt();
      fCandidateVariables[ 2] = part1->Px();
      fCandidateVariables[ 3] = part1->Py();
      fCandidateVariables[ 4] = part1->Pz();
      fCandidateVariables[ 5] = part2->Px();
      fCandidateVariables[ 6] = part2->Py();
      fCandidateVariables[ 7] = part2->Pz();
      fCandidateVariables[ 8] = casc->MassXi();
      fCandidateVariables[ 9] = TMath::Sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)); 
      if(xicobj->GetCharge()==-1)  fCandidateVariables[10] = casc->MassAntiLambda();
      else fCandidateVariables[10] = casc->MassLambda();
      fCandidateVariables[11]= casc->Pt();
      if(casc->ChargeXi()==-1) fCandidateVariables[12]= TMath::Sqrt(pow( casc->MomPosX(),2)+pow( casc->MomPosY(),2)); ///AGGIUSta
      else fCandidateVariables[12]=TMath::Sqrt(pow( casc->MomNegX(),2)+pow( casc->MomNegY(),2));
      
      fCandidateVariables[13] = fVtx1->GetX();
      fCandidateVariables[14] = fVtx1->GetY();
      fCandidateVariables[15] = fVtx1->GetZ();
      fCandidateVariables[16] = xicobj->GetOwnPrimaryVtx()->GetX();
      fCandidateVariables[17] = xicobj->GetOwnPrimaryVtx()->GetY();
      fCandidateVariables[18] = xicobj->GetOwnPrimaryVtx()->GetZ();
      
      fCandidateVariables[19] = xicobj->CascDcaXiDaughters();
      fCandidateVariables[20] = xicobj->CascDcaV0Daughters();
      fCandidateVariables[21] = xicobj->CascDecayLength();
      fCandidateVariables[22] = xicobj->CascCosPointingAngle();
      fCandidateVariables[23] = xicobj->CascDcaV0ToPrimVertex();
      fCandidateVariables[24] = xicobj->CascDcaPosToPrimVertex();
      fCandidateVariables[25] = xicobj->CascDcaNegToPrimVertex();
      fCandidateVariables[26] = xicobj->CascDcaBachToPrimVertex();
      fCandidateVariables[27] = xicobj->CascDecayLengthV0();
      fCandidateVariables[28] = xicobj->CascCosPointingAngleV0();
      
      fCandidateVariables[29] = dca[0];
      fCandidateVariables[30] = dca[1];
      fCandidateVariables[31] = dca[2];
      fCandidateVariables[32] = xicobj->Getd0Prong(0);
      fCandidateVariables[33] = xicobj->Getd0Prong(2);
      fCandidateVariables[34] = xicobj->Getd0Prong(1);
      
      fCandidateVariables[35] = xicobj->DecayLength();
      fCandidateVariables[36] = xicobj->DecayLengthXY();
      fCandidateVariables[37] = xicobj->XicCosPointingAngle();
      fCandidateVariables[38] = xicobj->BachelorsCosPointingAngle();
      
      if(fAnalCuts->GetIsUsePID())
	{
	  fCandidateVariables[39] = nSigmaTPCpi1;
	  fCandidateVariables[40] = nSigmaTPCpi2;
	  fCandidateVariables[41] = nSigmaTOFpi1;
	  fCandidateVariables[42] = nSigmaTOFpi2;
	  fCandidateVariables[43] = probPion1;
	  fCandidateVariables[44] = probPion2;
	}
      
      fCandidateVariables[45] = -9999;
      fCandidateVariables[46] = -9999;
      fCandidateVariables[47] = -9999;
      fCandidateVariables[48] = -9999;
      fCandidateVariables[49] = -9999;
      fCandidateVariables[50] = -9999;
      fCandidateVariables[51] = -9999;
      fCandidateVariables[52] = -9999;
      fCandidateVariables[53] = -9999;
      fCandidateVariables[54] = -9999;
      fCandidateVariables[58] = -9999;
      fCandidateVariables[60] = -9999;
      fCandidateVariables[61] = -9999;
      fCandidateVariables[62] = -9999;
 
 
      if(fUseMCInfo){
	if(mcpart){
	  fCandidateVariables[45] = mcpart->Label();
	  fCandidateVariables[46] = mcnused;
	  fCandidateVariables[47] = mcpart->GetPdgCode();
	  fCandidateVariables[51] = mcpart->Pt();
	  if(mcdaughter1&&mcdaughter2&&mcdaughterxi){
	    Double_t mcprimvertx = mcpart->Xv();
	    Double_t mcprimverty = mcpart->Yv();
	    Double_t mcsecvertx = mcdaughter1->Xv();
	    Double_t mcsecverty = mcdaughter1->Yv();
	    Double_t recosecvertx = xicobj->GetSecondaryVtx()->GetX();
	    Double_t recosecverty = xicobj->GetSecondaryVtx()->GetY();
	    fCandidateVariables[48] = TMath::Sqrt((mcsecvertx-mcprimvertx)*(mcsecvertx-mcprimvertx)+(mcsecverty-mcprimverty)*(mcsecverty-mcprimverty));
	    fCandidateVariables[49] = TMath::Sqrt((recosecvertx-mcprimvertx)*(recosecvertx-mcprimvertx)+(recosecverty-mcprimverty)*(recosecverty-mcprimverty));
	    Double_t vecx_vert = recosecvertx-mcprimvertx;
	    Double_t vecy_vert = recosecverty-mcprimverty;
	    Double_t vecl_vert = TMath::Sqrt(vecx_vert*vecx_vert+vecy_vert*vecy_vert);
	    Double_t vecx_mom = xicobj->Px();
	    Double_t vecy_mom = xicobj->Py();
	    Double_t vecl_mom = xicobj->Pt();
	    if(vecl_vert>0.&&vecl_mom>0.)
	      fCandidateVariables[50] = (vecx_vert*vecx_mom+vecy_vert*vecy_mom)/vecl_vert/vecl_mom;
	    fCandidateVariables[52] = mcdaughter1->Pt();
	    fCandidateVariables[53] = mcdaughter2->Pt();
	    fCandidateVariables[54] = mcdaughterxi->Pt();
	    fCandidateVariables[58] = static_cast<Float_t>(checkOrigin);
	    fCandidateVariables[60] = mcsecvertx;
	    fCandidateVariables[61] = mcsecverty;
	    fCandidateVariables[62] = mcdaughter1->Zv();
	  }
	}
      }
      
      fCandidateVariables[55] = xicobj->GetSecondaryVtx()->GetX();
      fCandidateVariables[56] = xicobj->GetSecondaryVtx()->GetY();
      fCandidateVariables[57] = xicobj->GetSecondaryVtx()->GetZ();
      fCandidateVariables[59] = xicobj->YXicPlus();
      
    }//close if to check mc fill only signal
    fVariablesTree->Fill();
  }//fWriteTree

  //NOTE: Important: need to add a selection on the Y of the XiCPlus
  if(fFillSparse){ //if you prefer to fill the histograms
    Double_t trueXic=-9999;
    if(isXic) trueXic=1.;
    else trueXic=0.;
    Double_t point[22]={xicobj->InvMassPiXiPi(), xicobj->Pt(), dca[2],dca[0], dca[1], xicobj->DecayLength(), xicobj->DecayLengthXY(), xicobj->XicCosPointingAngle(), xicobj->CascDcaXiDaughters(), xicobj->CascDcaV0Daughters(), xicobj->CascDcaV0ToPrimVertex(), xicobj->CascDcaPosToPrimVertex(), xicobj->CascDcaNegToPrimVertex(),xicobj->CascDcaBachToPrimVertex(),xicobj->CascCosPointingAngle(),xicobj->CascMassXi(), TMath::Sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)), nSigmaTPCpi1, nSigmaTPCpi2,  nSigmaTOFpi1,  nSigmaTOFpi2, trueXic};
    
    fSparseXicMass->Fill(point);
  }
    //NOTE: Important: need to add a selection on the Y of the XiCPlus
  if(!fFillSparse && !fWriteVariableTree){
    if(fAnalCuts->IsSelected(xicobj,AliRDHFCuts::kCandidate)){ 
      Double_t cont[3];
      cont[0] = xicobj->InvMassPiXiPi();
      cont[1] = xicobj->Pt();
      cont[2] = fCentrality;
      fHistoXicMass->Fill(cont);
      
      //these histos are filled only for candidates which pass IsSelected
      fHistoDcaPi1Pi2->Fill(dca[2] );
      fHistoDcaPi1Casc->Fill(dca[0]);
      fHistoDcaPi2Casc->Fill(dca[1]);
      fHistoLikeDecayLength->Fill(xicobj->DecayLength());
      fHistoLikeDecayLengthXY->Fill(xicobj->DecayLengthXY());
      fHistoXicCosPA->Fill(xicobj->XicCosPointingAngle());
      fHistoXiMass->Fill(xicobj->CascMassXi());
      fHistoCascDcaXiDaughters->Fill(xicobj->CascDcaXiDaughters());
      fHistoCascDcaV0Daughters->Fill(xicobj->CascDcaV0Daughters());
      fHistoCascDcaV0ToPrimVertex->Fill(xicobj->CascDcaV0ToPrimVertex());
      fHistoCascDcaPosToPrimVertex->Fill(xicobj->CascDcaPosToPrimVertex());
      fHistoCascDcaNegToPrimVertex->Fill(xicobj->CascDcaNegToPrimVertex());
      fHistoCascDcaBachToPrimVertex->Fill(xicobj->CascDcaBachToPrimVertex());
      fHistoCascCosPAXiPrim->Fill(xicobj->CascCosPointingAngle());
      fHistoXiPt->Fill(xicobj->PtProng(1));
      fHistoPiPt->Fill(xicobj->PtProng(0));
      fHistoPiPt->Fill(xicobj->PtProng(2));
      fHistoPid0->Fill(xicobj->Getd0Prong(0));
      fHistoPid0->Fill(xicobj->Getd0Prong(2));
      fHistonSigmaTPCpi->Fill(nSigmaTPCpi1);
      fHistonSigmaTPCpi->Fill(nSigmaTPCpi2);
      fHistonSigmaTOFpi->Fill(nSigmaTOFpi1);
      fHistonSigmaTOFpi->Fill(nSigmaTOFpi2);
      fHistoProbPion->Fill(probPion1);
      fHistoProbPion->Fill(probPion2);
    }
  }
 return;
}

//________________________________________________________________________
void AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::DefineTreeVariables() 
{
  //
  // This is to define tree variables
  //
  const char* nameoutput = GetOutputSlot(3)->GetContainer()->GetName();
  fVariablesTree = new TTree(nameoutput,"Candidates variables tree");
  Int_t nVar = 63;
  fCandidateVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="InvMassXic";
  fCandidateVariableNames[ 1]="XicPt";
  fCandidateVariableNames[ 2]="Pi1Px";
  fCandidateVariableNames[ 3]="Pi1Py";
  fCandidateVariableNames[ 4]="Pi1Pz";
  fCandidateVariableNames[ 5]="Pi2Px";
  fCandidateVariableNames[ 6]="Pi2Py";
  fCandidateVariableNames[ 7]="Pi2Pz";
  fCandidateVariableNames[ 8]="MassXi"; //
  fCandidateVariableNames[ 9]="XiPt"; //
  fCandidateVariableNames[10]="MassLambda";
  fCandidateVariableNames[11]="LambdaPt"; 
  fCandidateVariableNames[12]="ProtonPt";
    
  fCandidateVariableNames[13]="PrimVtxX";
  fCandidateVariableNames[14]="PrimVtxY";
  fCandidateVariableNames[15]="PrimVtxZ";
  fCandidateVariableNames[16]="NewPrimVtxX";
  fCandidateVariableNames[17]="NewPrimVtxY";
  fCandidateVariableNames[18]="NewPrimVtxZ";

  fCandidateVariableNames[19]="CascDcaXiDaughters"; //
  fCandidateVariableNames[20]="CascDcaV0Daughters"; //
  fCandidateVariableNames[21]="CascDecayLengthXi";
  fCandidateVariableNames[22]="CascCosPointingAngleXi"; //
  fCandidateVariableNames[23]="CascDcaV0ToPrimVertex"; //
  fCandidateVariableNames[24]="CascDcaPosToPrimVertex"; //
  fCandidateVariableNames[25]="CascDcaNegToPrimVertex"; //
  fCandidateVariableNames[26]="CascDcaBachToPrimVertex"; //
  fCandidateVariableNames[27]="CascDecayLengthV0";
  fCandidateVariableNames[28]="CascCosPointingAngleV0";

  fCandidateVariableNames[29]="DcaPi1Casc"; //
  fCandidateVariableNames[30]="DcaPi2Casc"; //
  fCandidateVariableNames[31]="DcaPi1Pi2"; //
  
  fCandidateVariableNames[32]="Pi1d0"; //
  fCandidateVariableNames[33]="Pi2d0"; 
  fCandidateVariableNames[34]="Cascd0";

  fCandidateVariableNames[35]="DecayLength"; //for pions only?
  fCandidateVariableNames[36]="DecayLengthXY"; //for pions only??
  fCandidateVariableNames[37]="XicCosPA"; //
  fCandidateVariableNames[38]="BachelorsCosPAXY";
  
  fCandidateVariableNames[39]="nSigmaTPCpi1"; //
  fCandidateVariableNames[40]="nSigmaTPCpi2";
  fCandidateVariableNames[41]="nSigmaTOFpi1";//
  fCandidateVariableNames[42]="nSigmaTOFpi2";
  fCandidateVariableNames[43]="probPion1";
  fCandidateVariableNames[44]="probPion2";

  fCandidateVariableNames[45]="mcxicID";
  fCandidateVariableNames[46]="mcnused";
  fCandidateVariableNames[47]="mcpdgcode";
  fCandidateVariableNames[48]="mcdecaylength";
  fCandidateVariableNames[49]="mcdecaylength_secsmear";
  fCandidateVariableNames[50]="mcxiccospaxy";
  fCandidateVariableNames[51]="mcxicpt";
  fCandidateVariableNames[52]="mcpi1pt";
  fCandidateVariableNames[53]="mcpi2pt";
  fCandidateVariableNames[54]="mcxipt";

  fCandidateVariableNames[55]="XicVertX";
  fCandidateVariableNames[56]="XicVertY";
  fCandidateVariableNames[57]="XicVertZ";
  
  fCandidateVariableNames[58]="isPrompt";
  fCandidateVariableNames[59]="Xicy";

  fCandidateVariableNames[60]="mcxicsecvertx";
  fCandidateVariableNames[61]="mcxicsecverty";
  fCandidateVariableNames[62]="mcxicsecvertz";
  
  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

//__________________________________________________________________________
void  AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::DefineGeneralHistograms() {
  //
  // This is to define general histograms
  //

  fCEvents = new TH1F("fCEvents","conter",25,-0.5,24.5);
  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetBinLabel(1," ");
  fCEvents->GetXaxis()->SetBinLabel(2,"Analyzed events");
  fCEvents->GetXaxis()->SetBinLabel(3,"AliAODVertex exists");
  fCEvents->GetXaxis()->SetBinLabel(4,"Trigger not OK");
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
  fCEvents->GetXaxis()->SetBinLabel(19,"CascSelected");
  fCEvents->GetXaxis()->SetBinLabel(20,"#pi #pi pair selected");
  fCEvents->GetXaxis()->SetBinLabel(21,"correct charge combination");
  fCEvents->GetXaxis()->SetBinLabel(22,"#Xi_{c} with raugh sel");
  fCEvents->GetXaxis()->SetBinLabel(23,"Reco #pi #pi sec. vtx");
  fCEvents->GetXaxis()->SetBinLabel(24,"HF3Prongs");
	
  //fCEvents->GetXaxis()->SetTitle("");
  fCEvents->GetYaxis()->SetTitle("counts");

  fHTrigger = new TH1F("fHTrigger","counter",18,-0.5,17.5);
  fHTrigger->SetStats(kTRUE);
  fHTrigger->GetXaxis()->SetBinLabel(1," ");
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

  Double_t binx[101];
  for(Int_t ib=0;ib<101;ib++){
    binx[ib] = 1.322-0.05 + 0.1/100.*(Double_t)ib ;
  }
  Double_t biny[21]={0.0,0.60,0.80,0.90,1.00,1.1,1.2,1.3,1.4,1.5,1.7,1.9,2.2,2.6,3.1,3.9,4.9,6.0,7.2,8.5,10.};
  fHistoXiMassvsPtRef1 = new TH2F("fHistoXiMassvsPtRef","Reference #Xi spectrum (-1.5 < #eta < -1.0)",100,binx,20,biny);
  fHistoXiMassvsPtRef2 = new TH2F("fHistoXiMassvsPtRef2","Reference #Xi spectrum (-1.0 < #eta < -0.5)",100,binx,20,biny);
  fHistoXiMassvsPtRef3 = new TH2F("fHistoXiMassvsPtRef3","Reference #Xi spectrum (-0.5 < #eta < -0.)",100,binx,20,biny);
  fHistoXiMassvsPtRef4 = new TH2F("fHistoXiMassvsPtRef4","Reference #Xi spectrum (0. < #eta < 0.5)",100,binx,20,biny);
  fHistoXiMassvsPtRef5 = new TH2F("fHistoXiMassvsPtRef5","Reference #Xi spectrum (0.5 < #eta < 1.)",100,binx,20,biny);
  fHistoXiMassvsPtRef6 = new TH2F("fHistoXiMassvsPtRef6","Reference #Xi spectrum (1. < #eta < 1.5)",100,binx,20,biny);
  fHistoPiPtRef = new TH1F("fHistoPiPtRef","Reference #pi spectrum",20,0.,10.);
  fHistoPiEtaRef = new TH1F("fHistoPiEtaRef","Reference #eta distributions of #pi ",50,-1,1.);

  fQAHistoNSelectedTracks = new TH1F("fQAHistoNSelectedTracks", "Number of tracks selected as pion candidates",100, 0, 100);
  fQAHistoNSelectedCasc = new TH1F("fQAHistoNSelectedCasc", "Number of tracks selected as cascades",20, 0, 20);

  fQAHistoDCApi1pi2 = new TH1F("fQAHistoDCApi1pi2","DCA #pi - #pi", 100,0.,0.5);

  fQAHistoAODPrimVertX = new TH1F("fQAHistoAODPrimVertX", "X coord of primary vertex", 100,0., 0.1);
  fQAHistoAODPrimVertY = new TH1F("fQAHistoAODPrimVertY", "Y coord of rimary vertex", 500,0., 0.5);
  fQAHistoAODPrimVertZ = new TH1F("fQAHistoAODPrimVertZ", "Z coord of primary vertex", 140,-12., 12.);
  
  fQAHistoRecoPrimVertX = new TH1F("fQAHistoRecoPrimVertX", "X coord of primary vertex without XiC decay tracks", 100,0., 0.1);
  fQAHistoRecoPrimVertY = new TH1F("fQAHistoRecoPrimVertY", "Y coord of primary vertex without XiC decay tracks", 500,0., 0.5);
  fQAHistoRecoPrimVertZ = new TH1F("fQAHistoRecoPrimVertZ", "Z coord of primary vertex without XiC decay tracks", 140,-12, 12);

  fQAHistoSecondaryVertexX = new TH1F("fQAHistoSecondaryVertexX", "X coord of secondary vertex", 1000, -2, 2);
  fQAHistoSecondaryVertexY = new TH1F("fQAHistoSecondaryVertexY", "Y coord of secondary vertex", 1000, -2, 2);
  fQAHistoSecondaryVertexZ = new TH1F("fQAHistoSecondaryVertexZ", "Z coord of secondary vertex", 1000, -20, 20);
  fQAHistoSecondaryVertexXY = new TH1F("fQAHistoSecondaryVertexXY", "XY coord of secondary vertex", 1000, -20, 20);
  
				  
  fOutput->Add(fCEvents);
  fOutput->Add(fHTrigger);
  fOutput->Add(fHCentrality);
  fOutput->Add(fHistoXiMassvsPtRef1);
  fOutput->Add(fHistoXiMassvsPtRef2);
  fOutput->Add(fHistoXiMassvsPtRef3);
  fOutput->Add(fHistoXiMassvsPtRef4);
  fOutput->Add(fHistoXiMassvsPtRef5);
  fOutput->Add(fHistoXiMassvsPtRef6);
  fOutput->Add(fHistoPiPtRef);
  fOutput->Add(fHistoPiEtaRef);
  fOutput->Add(fQAHistoNSelectedTracks);
  fOutput->Add(fQAHistoNSelectedCasc);
  fOutput->Add(fQAHistoDCApi1pi2);
  fOutput->Add(fQAHistoAODPrimVertX);
  fOutput->Add(fQAHistoAODPrimVertY);
  fOutput->Add(fQAHistoAODPrimVertZ);
  fOutput->Add(fQAHistoRecoPrimVertX);
  fOutput->Add(fQAHistoRecoPrimVertY);
  fOutput->Add(fQAHistoRecoPrimVertZ);
  fOutput->Add(fQAHistoSecondaryVertexX);
  fOutput->Add(fQAHistoSecondaryVertexY);
  fOutput->Add(fQAHistoSecondaryVertexZ);
  fOutput->Add(fQAHistoSecondaryVertexXY);
  
  if(fUseMCInfo) {
    fHistoMCSpectrumAccXic=new TH3F("fHistoMCSpectrumAccXic","fHistoMCSpectrumAccXic",250,0,50,26,-0.5,12.5,2,3.5,5.5);
    fOutput->Add(fHistoMCSpectrumAccXic);
  }

    
  return;
}

//__________________________________________________________________________
void  AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::DefineAnalysisHistograms() 
{
  //
  // Define histograms
  //
  
  //------------------------------------------------
  // Basic histograms
  //------------------------------------------------
  Int_t bins_base[3]=		{80,20,10};
  Double_t xmin_base[3]={2.468-0.2,0,0.00};
  Double_t xmax_base[3]={2.468+0.2,20.,100};
  fHistoXicMass = new THnSparseF("fHistoXicMass","fHistoXicMass;mass:pt:centrality",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoXicMass);

  if (fFillSparse){
    Int_t bins_sparse[22]=   {      100,  20,  50,  50,  50,  50,  50,  50, 50, 50, 50, 50, 50, 50,  20,         50,   50,  50,  50,  50,  50,    2};
    Double_t xmin_sparse[22]={2.468-0.2,  0.,  0.,  0.,  0.,  0.,  0., -1., 0., 0. , 0.,  0.,  0., 0., 0.9, 1.322-0.05,   0., -5., -5., -5., -5., -0.5};
    Double_t xmax_sparse[22]={2.468+0.2, 20., 0.5,  5.,  5.,  3.,  3.,  1., 2., 1.5, 5., 0.5, 0.5, 1.,  1., 1.322+0.05,  10.,  5.,  5.,  5.,  5.,  1.5};
    fSparseXicMass = new THnSparseF("fSparseXicMass","fSparseXicMass;mass:pt:dcapi1pi2:dcapi1Xi:dcapi2Xi:decayLengthPi1Pi2:DecayLengthPi1P2XY:XicCosPA:DCAXidau:DCAV0dau:CascV0toPV:DCAPosToPrim:DCANegToPrim:DCABachToPrim:CosPAXiPrim:MassXi:XiPt:nSigmaTPCPi1:nSigmaTPCPi2:nSigmaTOFPi1:nSigmaTOFp2:isXic",22,bins_sparse,xmin_sparse, xmax_sparse);
    fOutputAll->Add(fSparseXicMass);
  }
  else {
  //------------------------------------------------
  // Checking histograms
  //------------------------------------------------
  fHistoDcaPi1Pi2 = new TH1F("fHistoDcaPi1Pi2","DCA (#pi_{1}-#pi_{2})",100,0.0,0.5);
  fOutputAll->Add(fHistoDcaPi1Pi2);
  fHistoDcaPi1Casc = new TH1F("fHistoDcaPi1Casc","DCA (#pi_{1}-#Xi)",100,0.0,1.0);
  fOutputAll->Add(fHistoDcaPi1Casc);
  fHistoDcaPi2Casc = new TH1F("fHistoDcaPi2Casc","DCA (#pi_{2}-#Xi)",100,0.0,1.0);
  fOutputAll->Add(fHistoDcaPi2Casc);
  fHistoLikeDecayLength = new TH1F("fHistoLikeDecayLength","Decay Length (#pi-#pi)",100,0.,0.2);
  fOutputAll->Add(fHistoLikeDecayLength);
  fHistoLikeDecayLengthXY = new TH1F("fHistoLikeDecayLengthXY","Decay Length (#pi-#pi)",100,0.,0.2);
  fOutputAll->Add(fHistoLikeDecayLengthXY);
  fHistoXicCosPA = new TH1F("fHistoXicCosPA","#Xi_{c} cos(pa) ",100,-1.0,1.0);
  fOutputAll->Add(fHistoXicCosPA);
  fHistoXiMass=new TH1F("fHistoXiMass","#Xi^{-} Mass",100,1.322-0.05,1.322+0.05);
  fOutputAll->Add(fHistoXiMass);
  fHistoCascDcaXiDaughters=new TH1F("fHistoCascDcaXiDaughters","DCA #Xi daughters ",100,0.0,1.0);
  fOutputAll->Add(fHistoCascDcaXiDaughters);
  fHistoCascDcaV0Daughters=new TH1F("fHistoCascDcaV0Daughters","DCA #V0 daughters ",100,0.0,1.0);
  fOutputAll->Add(fHistoCascDcaV0Daughters);
  fHistoCascDcaV0ToPrimVertex=new TH1F("fHistoCascDcaV0ToPrimVertex","DCA V0 daughters to primary vertex ",100,0.0,1.0);
  fOutputAll->Add(fHistoCascDcaV0ToPrimVertex);
  fHistoCascDcaPosToPrimVertex=new TH1F("fHistoCascDcaPosToPrimVertex","DCA Pos-Prim ",100,0.0,1.0);
  fOutputAll->Add(fHistoCascDcaPosToPrimVertex);
  fHistoCascDcaNegToPrimVertex=new TH1F("fHistoCascDcaNegToPrimVertex","DCA Neg-Prim ",100,0.0,1.0);
  fOutputAll->Add(fHistoCascDcaNegToPrimVertex);
  fHistoCascDcaBachToPrimVertex=new TH1F("fHistoCascDcaBachToPrimVertex","DCA Bach-Prim ",100,0.0,1.0);
  fOutputAll->Add(fHistoCascDcaBachToPrimVertex);
  fHistoCascCosPAXiPrim=new TH1F("fHistoCascCosPAXiPrim","#Xi CosPA (prim)",100,0.8,1.0);
  fOutputAll->Add(fHistoCascCosPAXiPrim);
  fHistoXiPt=new TH1F("fHistoXiPt","#Xi^{-} p_{T}",100,0.,10.);
  fOutputAll->Add(fHistoXiPt);
  fHistoPiPt=new TH1F("fHistoPiPt","#pi p_{T}",100,0.,10);
  fOutputAll->Add(fHistoPiPt);
  fHistoPid0=new TH1F("fHistoPid0","#pi d_{0}",100,-0.1,0.1);
  fOutputAll->Add(fHistoPid0);
  fHistonSigmaTPCpi=new TH1F("fHistonSigmaTPCpi","n#sigma (TPC, pion)",100,-10.,10.);
  fOutputAll->Add(fHistonSigmaTPCpi);
  fHistonSigmaTOFpi=new TH1F("fHistonSigmaTOFpi","n#sigma (TOF, pion)",100,-10.,10.);
  fOutputAll->Add(fHistonSigmaTOFpi);
  fHistoProbPion=new TH1F("fHistoProbPion","Bayse Prob",100,0.0,1.);
  fOutputAll->Add(fHistoProbPion);
  }
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::SelectTrack( const AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks,Bool_t *seleFlags)
{
  //
  // Select good tracks using fAnalCuts (AliRDHFCuts object) and return the array of their ids
  //
  
  //const Int_t entries = event->GetNumberOfTracks();
  if(trkEntries==0) return;
  
  nSeleTrks=0;
  for(Int_t i=0; i<trkEntries; i++) {
    seleFlags[i] = kFALSE;
    
    AliVTrack *track;
    track = (AliVTrack*)event->GetTrack(i);
    
    if(track->GetID()<0) continue;
    Double_t covtest[21];
    if(!track->GetCovarianceXYZPxPyPz(covtest)) continue;
    
    AliAODTrack *aodt = (AliAODTrack*)track;

    if(!fAnalCuts) continue;
    if(fAnalCuts->SingleTrkCuts(aodt)){
      seleFlags[i]=kTRUE;
      nSeleTrks++;
      fHistoPiPtRef->Fill(aodt->Pt());
      fHistoPiEtaRef->Fill(aodt->Eta());
    }
  } // end loop on tracks
}

//________________________________________________________________________
void AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::SelectCascade( const AliVEvent *event,Int_t nCascades,Int_t &nSeleCasc, Bool_t *seleCascFlags)
{
  //
  // Select good cascade using fAnalCuts (AliRDHFCuts object) and return the array of their ids
  //

  Double_t primVtx[3];
  fVtx1->GetXYZ(primVtx);

  nSeleCasc = 0;
  for(Int_t icasc=0;icasc<nCascades;icasc++)
    {
      seleCascFlags[icasc] = kFALSE;
      AliAODcascade *casc = ((AliAODEvent*)event)->GetCascade(icasc);
  
      if(!fAnalCuts) continue;
      if(fAnalCuts->SingleCascadeCuts(casc,primVtx)){
	seleCascFlags[icasc] = kTRUE;
	nSeleCasc++;
      }
      if(fAnalCuts->SingleCascadeCutsRef(casc,primVtx)) //the only difference wrt SingleCascadeCuts is the cut on the Xi mass
      {
        Double_t rapxi = casc->RapXi();
        if(rapxi>-1.5&&rapxi<-1.0){
  fHistoXiMassvsPtRef1->Fill(casc->MassXi(),sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)));
        }
        if(rapxi>-1.0&&rapxi<-0.5){
  fHistoXiMassvsPtRef2->Fill(casc->MassXi(),sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)));
        }
        if(rapxi>-0.5&&rapxi<0.0){
  fHistoXiMassvsPtRef3->Fill(casc->MassXi(),sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)));
        }
        if(rapxi>0.0&&rapxi<0.5){
  fHistoXiMassvsPtRef4->Fill(casc->MassXi(),sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)));
        }
        if(rapxi>0.5&&rapxi<1.0){
  fHistoXiMassvsPtRef5->Fill(casc->MassXi(),sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)));
        }
        if(rapxi>1.0&&rapxi<1.5){
  fHistoXiMassvsPtRef6->Fill(casc->MassXi(),sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)));
        }
      }
    } //loop over cascades
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::SelectLikeSign(AliAODTrack *trk1, AliAODTrack *trk2)
{
  //
  // Select LikeSign tracks
  //

  if(trk1->Charge()!=trk2->Charge()) return kFALSE;
  if(trk1->GetID()==trk2->GetID()) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::CallPrimaryVertex(AliAODcascade *casc, AliAODTrack *trk1, AliAODTrack *trk2, AliAODEvent* aod)
{
  //
  // Make an array of tracks which should not be used in primary vertex calculation and 
  // Call PrimaryVertex function
  //

  TObjArray *twoTrackArrayPlusXi = new TObjArray(5);
  
  AliESDtrack *cptrk1 = new AliESDtrack((AliVTrack*)trk1);
  twoTrackArrayPlusXi->AddAt(cptrk1,0);
  AliESDtrack *cptrk2 = new AliESDtrack((AliVTrack*)trk2);
  twoTrackArrayPlusXi->AddAt(cptrk2,1);
  
  AliESDtrack *cascptrack = new AliESDtrack((AliVTrack*)casc->GetDaughter(0));
  twoTrackArrayPlusXi->AddAt(cascptrack,2);
  AliESDtrack *cascntrack = new AliESDtrack((AliVTrack*)casc->GetDaughter(1));
  twoTrackArrayPlusXi->AddAt(cascntrack,3);
  AliESDtrack *cascbtrack = new AliESDtrack((AliVTrack*)casc->GetDecayVertexXi()->GetDaughter(0));
  twoTrackArrayPlusXi->AddAt(cascbtrack,4);
  
  AliAODVertex *newvert  = PrimaryVertex(twoTrackArrayPlusXi,aod);
  
  for(Int_t i=0;i<5;i++)
    {
      AliESDtrack *tesd = (AliESDtrack*)twoTrackArrayPlusXi->UncheckedAt(i);
      delete tesd;
    }
  twoTrackArrayPlusXi->Clear();
  delete twoTrackArrayPlusXi;
  
  return newvert;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::CallReconstructSecondaryVertex(AliAODTrack *trk1, AliAODTrack *trk2, Double_t &dispersion)
{
  //
  // Make an array of tracks which I want to use for vertex calculation and call ReconstructSecondaryVertex
  //

  TObjArray *trkArray;
  trkArray= new TObjArray(2);
  
  AliESDtrack *cptrk1 = new AliESDtrack((AliVTrack*)trk1);
  trkArray->AddAt(cptrk1,0);
  AliESDtrack *cptrk2 = new AliESDtrack((AliVTrack*)trk2);
  trkArray->AddAt(cptrk2,1);
  
  Double_t xdummy, ydummy;
  Double_t dcap1p2 = cptrk1->GetDCA(cptrk2,fBzkG,xdummy,ydummy);
  fQAHistoDCApi1pi2->Fill(dcap1p2);
  AliAODVertex *secvertex=0;
  if(dcap1p2<fAnalCuts->GetProdLikeSignDcaMax()){
    secvertex = ReconstructSecondaryVertex(trkArray,dispersion);
  }
  
  for(Int_t i=0;i<2;i++)
    {
      AliESDtrack *tesd = (AliESDtrack*)trkArray->UncheckedAt(i);
      delete tesd;
    }
  trkArray->Clear();
  delete trkArray;
  
  return secvertex;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::PrimaryVertex(const TObjArray *trkArray,
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
AliAODVertex* AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::ReconstructSecondaryVertex(TObjArray *trkArray, Double_t &dispersion,Bool_t useTRefArray)
{
  //
  // Reconstruct secondary vertex from trkArray (Copied from AliAnalysisVertexingHF)
  //
  
  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;
  
  AliVertexerTracks *fVertexerTracks = new AliVertexerTracks(fBzkG);
  
  fVertexerTracks->SetVtxStart(fV1);
  vertexESD = (AliESDVertex*)fVertexerTracks->VertexForSelectedESDTracks(trkArray);
  
  delete fVertexerTracks; fVertexerTracks=NULL;

  if(!vertexESD) return vertexAOD;

  if(vertexESD->GetNContributors()!=trkArray->GetEntriesFast()) { 
    //AliDebug(2,"vertexing failed"); 
    delete vertexESD; vertexESD=NULL;
    return vertexAOD;
  }
  
  Double_t vertRadius2=vertexESD->GetX()*vertexESD->GetX()+vertexESD->GetY()*vertexESD->GetY();
  if(vertRadius2>8.){
    // vertex outside beam pipe, reject candidate to avoid propagation through material
    delete vertexESD; vertexESD=NULL;
    return vertexAOD;
  }
  
  
  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  dispersion = vertexESD->GetDispersion();
  delete vertexESD; vertexESD=NULL;
  
  Int_t nprongs= (useTRefArray ? 0 : trkArray->GetEntriesFast());
  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,nprongs);
  
  return vertexAOD;
}

//________________________________________________________________________
AliAODRecoCascadeHF3Prong* AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::MakeCascadeHF3Prong(AliAODcascade *casc, AliAODTrack *part1, AliAODTrack *part2, AliAODEvent * aod, AliAODVertex *secVert, Double_t dispersion) 
{
  //
  // Make AliAODRecoCascadeHF3Prong object from the arguments
  //
	
  //------------------------------------------------
  // PrimaryVertex
  //------------------------------------------------
  AliAODVertex *primVertexAOD;
  Bool_t unsetvtx = kFALSE;
  if(fReconstructPrimVert){
    primVertexAOD = CallPrimaryVertex(casc,part1,part2,aod);
    if(!primVertexAOD){
      primVertexAOD = fVtx1;
    }else{
      unsetvtx = kTRUE;
    }
  }else{
    primVertexAOD = fVtx1;
  }
  if(!primVertexAOD) return 0x0;
  Double_t pos[3]; primVertexAOD->GetXYZ(pos);
  fQAHistoRecoPrimVertX->Fill(pos[0]);
  fQAHistoRecoPrimVertY->Fill(pos[1]);
  fQAHistoRecoPrimVertZ->Fill(pos[2]);
  

  //------------------------------------------------
  // DCA between tracks
  //------------------------------------------------
  AliESDtrack *esdtrack1 = new AliESDtrack((AliVTrack*)part1);
  AliESDtrack *esdtrack2 = new AliESDtrack((AliVTrack*)part2);
  Double_t xyz[3], pxpypz[3], cv[21]; Short_t sign;
  xyz[0]=casc->DecayVertexXiX();
  xyz[1]=casc->DecayVertexXiY();
  xyz[2]=casc->DecayVertexXiZ();
  pxpypz[0]=casc->MomXiX();
  pxpypz[1]=casc->MomXiY();
  pxpypz[2]=casc->MomXiZ();
  casc->GetCovarianceXYZPxPyPz(cv);
  sign=casc->ChargeXi();
  AliExternalTrackParam	*trackCasc = new AliExternalTrackParam(xyz,pxpypz,cv,sign);

  Double_t xdummy, ydummy;
  Double_t dcap1p2 = esdtrack1->GetDCA(esdtrack2,fBzkG,xdummy,ydummy);
  Double_t dcap1casc = esdtrack1->GetDCA(trackCasc,fBzkG,xdummy,ydummy);
  Double_t dcap2casc = esdtrack2->GetDCA(trackCasc,fBzkG,xdummy,ydummy);
  Double_t dca[3]={dcap1casc,dcap2casc,dcap1p2};

  //------------------------------------------------
  // Propagate all tracks to the secondary vertex and calculate momentum there
  //------------------------------------------------

  Double_t d0z0[2],covd0z0[3];
  Double_t secR = TMath::Sqrt(secVert->GetX()*secVert->GetX()+secVert->GetY()*secVert->GetY());
  if(secR<1.0){
    part1->PropagateToDCA(secVert,fBzkG,kVeryBig,d0z0,covd0z0);
    part2->PropagateToDCA(secVert,fBzkG,kVeryBig,d0z0,covd0z0);
    trackCasc->PropagateToDCA(secVert,fBzkG,kVeryBig);
  }else{
    part1->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
    part2->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
    trackCasc->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig);
  }

  Double_t momxi_new[3]={-9999,-9999,-9999.};
  trackCasc->GetPxPyPz(momxi_new);

  Double_t px[3],py[3],pz[3];
  px[0] = part1->Px(); py[0] = part1->Py(); pz[0] = part1->Pz(); 
  px[1] = momxi_new[0]; py[1] = momxi_new[1]; pz[1] = momxi_new[2]; 
  px[2] = part2->Px(); py[2] = part2->Py(); pz[2] = part2->Pz(); 

  //------------------------------------------------
  // d0
  //------------------------------------------------
  Double_t d0[3],d0err[3];
  part1->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  d0[0]= d0z0[0];
  d0err[0] = TMath::Sqrt(covd0z0[0]);

  trackCasc->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  d0[1]= d0z0[0];
  d0err[1] = TMath::Sqrt(covd0z0[0]);//Do not know what to use

  part2->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  d0[2]= d0z0[0];
  d0err[2] = TMath::Sqrt(covd0z0[0]);


  //------------------------------------------------
  // Other stuff
  //------------------------------------------------
  Double_t dist12=0.0;//Vertex pi-casc not calculated
  Double_t dist23=0.0;//Vertex pi-casc not calculated
  Short_t charge=(Short_t)(part1->Charge()+trackCasc->Charge()+part2->Charge());

  //------------------------------------------------
  // Create AliAODRecoCascadeHF3Prong
  //------------------------------------------------
  AliAODRecoCascadeHF3Prong *theCascade = new AliAODRecoCascadeHF3Prong(secVert,charge,px,py,pz,d0,d0err,dca,dispersion,dist12,dist23);
  if(!theCascade)  
    {
      if(unsetvtx) delete primVertexAOD; primVertexAOD=NULL;
      if(esdtrack1) delete esdtrack1;
      if(esdtrack2) delete esdtrack2;
      if(trackCasc) delete trackCasc;
      return 0x0;
    }
  theCascade->SetOwnPrimaryVtx(primVertexAOD);
  UShort_t id[3]={(UShort_t)part1->GetID(),(UShort_t)trackCasc->GetID(),(UShort_t)part2->GetID()};
  theCascade->SetProngIDs(3,id);

  theCascade->GetSecondaryVtx()->AddDaughter(part1);
  theCascade->GetSecondaryVtx()->AddDaughter(casc);
  theCascade->GetSecondaryVtx()->AddDaughter(part2);

  if(unsetvtx) delete primVertexAOD; primVertexAOD=NULL;
  if(esdtrack1) delete esdtrack1;
  if(esdtrack2) delete esdtrack2;
  if(trackCasc) delete trackCasc;

  return theCascade;
}

//________________________________________________________________________
void AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::LoopOverGenParticles(TClonesArray *mcArray){
  for(Int_t kmc=0; kmc<mcArray->GetEntries(); kmc++){
    AliAODMCParticle *mcpart=(AliAODMCParticle*)mcArray->At(kmc);
    if(!mcpart) continue;

    Int_t pdg=TMath::Abs(mcpart->GetPdgCode());   
    Int_t arrayDauLab[5];
    if(pdg==4232){
      if(CheckXic2XiPiPi(mcArray,mcpart,arrayDauLab)==1){ //the arrayDauLab is used to check if the single particles are in acceptance when the Xic is in the acceptance.
	Int_t checkOrigin=AliVertexingHFUtils::CheckOrigin(mcArray,mcpart,kTRUE);
	if(checkOrigin==0)continue;

	Float_t ptpart=mcpart->Pt();
	Float_t ypart=mcpart->Y();

	if(TMath::Abs(ypart)<0.5){
	  fHistoMCSpectrumAccXic->Fill(ptpart,kGenLimAcc,checkOrigin);
	}
	Bool_t isInAcc=kTRUE;
	// check GenAcc level
	if(fAnalCuts){
	  if(!fAnalCuts->IsInFiducialAcceptance(ptpart,ypart)){ 
	    isInAcc=kFALSE;
	  }
	} else if (TMath::Abs(ypart)>0.8) isInAcc=kFALSE;
	if (TMath::Abs(ypart)<0.8) {
	  fHistoMCSpectrumAccXic->Fill(ptpart,kGenAccMother08,checkOrigin);
	  for(Int_t k=0;k<5;k++){
	    AliAODMCParticle *mcpartdau=(AliAODMCParticle*)mcArray->At(arrayDauLab[k]);
	    if(TMath::Abs(mcpartdau->Eta())>0.8) fHistoMCSpectrumAccXic->Fill(ptpart,kGenAcc08,checkOrigin);;    
	  }
	}
	
	if(isInAcc){
	  fHistoMCSpectrumAccXic->Fill(ptpart,kGenAccMother,checkOrigin);
	  for(Int_t k=0;k<5;k++){
	    AliAODMCParticle *mcpartdau=(AliAODMCParticle*)mcArray->At(arrayDauLab[k]);
	    if(TMath::Abs(mcpartdau->Eta())>0.8) isInAcc=kFALSE;    
	  }
	}
	if(isInAcc) fHistoMCSpectrumAccXic->Fill(ptpart,kGenAcc,checkOrigin);
      } //else continue;//CheckXic2XiPiPi
    }//Check on PDG code
  }//loop over particles
}

//________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks::CheckXic2XiPiPi(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab){
  if(!mcPart) return -1;
  
  Int_t pdgD=mcPart->GetPdgCode();
  
  if(TMath::Abs(pdgD)!=4232) return -1;

   Int_t nDau=mcPart->GetNDaughters();
   Int_t labelFirstDau = mcPart->GetDaughterLabel(0);

   Int_t nXi=0;
   Int_t nLambda=0;
   Int_t nPions=0;
   Int_t nProtons=0;
   Double_t sumPxDau=0.;
   Double_t sumPyDau=0.;
   Double_t sumPzDau=0.;
   Int_t nFound=0;

   Int_t codeRes=-1;
   if(nDau==3){
     for(Int_t iDau=0; iDau<nDau; iDau++){
       Int_t indDau = labelFirstDau+iDau;
       if(indDau<0) return -1;
       AliAODMCParticle* dau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indDau));
       if(!dau) return -1;
       Int_t pdgdau=dau->GetPdgCode();
       if(TMath::Abs(pdgdau)==3312){ //Xi
	 nXi++;
	 sumPxDau+=dau->Px();
	 sumPyDau+=dau->Py();
	 sumPzDau+=dau->Pz();
	 Int_t nXiDau=dau->GetNDaughters();
	 if(nXiDau!=2) return -1; //no correct Xi decay
	 Int_t indFirstXiDau=dau->GetDaughterLabel(0);
	 for(Int_t XiDau=0; XiDau<2; XiDau++){
	   Int_t indXiDau=indFirstXiDau+XiDau;
	   if(indXiDau<0) return -1;
	   AliAODMCParticle* Xidau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indXiDau));
	   if(!Xidau) return -1;
	   Int_t pdgXidau=Xidau->GetPdgCode();
	   if(TMath::Abs(pdgXidau)==3122){
	     nLambda++;
	     Int_t nLambdaDau=Xidau->GetNDaughters();
	     if (nLambdaDau!=2) return -1;
	     Int_t indFirstLambdaDau=Xidau->GetDaughterLabel(0);
	     for(Int_t LambdaDau=0; LambdaDau<2; LambdaDau++){
	       Int_t indLambdaDau=indFirstLambdaDau+LambdaDau;
	       AliAODMCParticle* Lambdadau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indLambdaDau));
	       if(!Lambdadau) return -1;
	       Int_t pdgLambdadau=Lambdadau->GetPdgCode();
	       if(TMath::Abs(pdgLambdadau)==2212){
		 nProtons++;
		 arrayDauLab[nFound++]=indLambdaDau;
		 if (nFound>5) return -1;
	       } else if(TMath::Abs(pdgLambdadau)==211){
		 nPions++;
		 arrayDauLab[nFound++]=indLambdaDau;
		 if (nFound>5) return -1;
	       }
	     }
	   }else if(TMath::Abs(pdgXidau)==211) {
	     nPions++;
	     arrayDauLab[nFound++]=indXiDau;
	     if (nFound>5) return -1;
	   }
	 }
       }else if(TMath::Abs(pdgdau)==211){
	 nPions++;
	 sumPxDau+=dau->Px();
	 sumPyDau+=dau->Py();
	 sumPzDau+=dau->Pz();
	 arrayDauLab[nFound++]=indDau;
	 if (nFound>5) return -1;
       }
     }
     if(nXi!=1) return -1;
     if(nLambda!=1) return -1;
     if(nPions!=4) return -1;
     if(nProtons!=1) return -1;
     if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
     if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
     if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;

     if(nXi==1 && nLambda==1 && nPions==4 && nProtons ==1) return 1;
   } else return -1;
}
	     
