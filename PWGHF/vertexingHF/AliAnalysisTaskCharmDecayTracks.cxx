/**************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: $ */

//*************************************************************************
// Class AliAnalysisTaskCharmDecayTracks
/////////////////////////////////////////////////////////////

#include <TList.h>
#include <TH1F.h>
#include <TDatabasePDG.h>
#include <TTree.h>
#include <TChain.h>

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisVertexingHF.h"
#include "AliVertexingHFUtils.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisTaskCharmDecayTracks.h"
#include "AliAnalysisUtils.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskCharmDecayTracks);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskCharmDecayTracks::AliAnalysisTaskCharmDecayTracks():
  AliAnalysisTaskSE("CharmDecayTracks"),
  fOutput(0x0),
  fHistNEvents(0x0),
  fHistNCand(0x0),
  fHistTrLab(0x0),
  fHistCluTPCDupLab(0x0),
  fHistCluTPCDupLabCorrel(0x0),
  fHistCluITSDupLabCorrel(0x0),
  fHistMomDupLab(0x0),
  fTrackTree(0x0),
  fTreeVarInt(0x0),
  fTreeVarFloat(0x0),
  fTrPar1(),
  fTrPar2(),
  fTrPar3(),
  fTrParV0(),
  fPVertexTrk(0x0),
  fSelSpecies(421),
  fDecayMode(0),
  fFilterMask(16),
  fTrCuts(0x0),
  fReadMC(kFALSE),
  fUsePhysSel(kFALSE),
  fUsePileupCut(kFALSE),
  fTriggerMask(AliVEvent::kAnyINT),
  fGoUpToQuark(kTRUE),
  fKeepNegID(kFALSE),
  fMethod(0)
{
  /// default constructor

  DefineInput(0, TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskCharmDecayTracks::~AliAnalysisTaskCharmDecayTracks()
{
  //
  /// Destructor
  //
  if(fOutput && !fOutput->IsOwner()){
    delete fHistNEvents;
    delete fHistNCand;
    delete fHistTrLab;
    delete fHistCluTPCDupLab;
    delete fHistCluTPCDupLabCorrel;
    delete fHistCluITSDupLabCorrel;
    delete fHistMomDupLab;
    delete fTrackTree;
  }

  delete fOutput;
  delete fTrCuts;
  delete [] fTreeVarInt;
  delete [] fTreeVarFloat;

}

//________________________________________________________________________
void AliAnalysisTaskCharmDecayTracks::UserCreateOutputObjects()
{
  /// Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskCharmDecayTracks::UserCreateOutputObjects() \n");
  
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");
  
  fHistNEvents = new TH1F("hNEvents", "Number of processed events",15,-0.5,14.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"All events");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"PhysSel");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"Good vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Pass zSPD-zTrk vert sel");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"|zvert|<10");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"Pileup cut");
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fHistNCand = new TH1F("hNCand", "Number of candidates",15,-0.5,14.5);
  fHistNCand->GetXaxis()->SetBinLabel(1,"Events with method 0");
  fHistNCand->GetXaxis()->SetBinLabel(2,"Events with method 1");
  fHistNCand->GetXaxis()->SetBinLabel(3,Form("Pdg %d in kine",fSelSpecies));
  fHistNCand->GetXaxis()->SetBinLabel(4,"Good decay channel");
  fHistNCand->GetXaxis()->SetBinLabel(5,"MC truth OK");
  fHistNCand->GetXaxis()->SetBinLabel(6,"Daughter particle not tracked");
  fHistNCand->GetXaxis()->SetBinLabel(7,"Daughter track null pointer");
  fHistNCand->GetXaxis()->SetBinLabel(8,"Daughter track label mismatch");
  fHistNCand->GetXaxis()->SetBinLabel(9,"Daughter track not selected");
  fHistNCand->GetXaxis()->SetBinLabel(10,"Daughter tracks OK");
  fHistNCand->GetXaxis()->SetBinLabel(11,Form("Signal cand %d in AOD",fSelSpecies));
  fHistNCand->GetXaxis()->SetBinLabel(12,"MC truth OK");
  fHistNCand->GetXaxis()->SetBinLabel(13,"Daughter track null pointer");
  fHistNCand->GetXaxis()->SetBinLabel(14,"Daughter track not selected");
  fHistNCand->GetXaxis()->SetBinLabel(15,"Daughter tracks OK");
  fOutput->Add(fHistNCand);
  
  fHistTrLab=new TH1F("hTrLab","",12,-1.5,10.5);
  fHistCluTPCDupLab=new TH2F("hCluTPCDupLab","",10,0.5,10.5,160,-0.5,159.5);
  fHistCluTPCDupLabCorrel=new TH2F("hCluTPCDupLabCorrel","",160,-0.5,159.5,160,-0.5,159.5);
  fHistCluITSDupLabCorrel=new TH2F("hCluITSDupLabCorrel","",7,-0.5,6.5,7,-0.5,6.5);
  fHistMomDupLab=new TH2F("fHistMomDupLab","",10,0.5,10.5,100,0.,5.);
  fOutput->Add(fHistTrLab);
  fOutput->Add(fHistCluTPCDupLab);
  fOutput->Add(fHistCluTPCDupLabCorrel);
  fOutput->Add(fHistCluITSDupLabCorrel);
  fOutput->Add(fHistMomDupLab);
  
  fTrackTree = new TTree("trackTree", "Tree for analysis");
  TString intVarName[kNumOfIntVar];
  fTreeVarInt = new Int_t[kNumOfIntVar];
  intVarName[0]="pdg"; // PDG code of mother
  intVarName[1]="origin"; // charm or beauty
  intVarName[2]="pdgDau1"; // PDG of daughter
  intVarName[3]="pdgDau2"; // PDG of daughter
  intVarName[4]="pdgDau3"; // PDG of daughter
  for(Int_t ivar=0; ivar<kNumOfIntVar; ivar++){
    fTrackTree->Branch(intVarName[ivar].Data(),&fTreeVarInt[ivar],Form("%s/I",intVarName[ivar].Data()));
  }
  TString floatVarName[kNumOfFloatVar];
  fTreeVarFloat = new Float_t[kNumOfFloatVar];
  floatVarName[0]="ptgenD";
  floatVarName[1]="phigenD";
  floatVarName[2]="ygenD";
  floatVarName[3]="xcollv";
  floatVarName[4]="ycollv";
  floatVarName[5]="zcollv";
  floatVarName[6]="xprodv";
  floatVarName[7]="yprodv";
  floatVarName[8]="zprodv";
  floatVarName[9]="xdecv";
  floatVarName[10]="ydecv";
  floatVarName[11]="zdecv";
  for(Int_t ivar=0; ivar<kNumOfFloatVar; ivar++){
    fTrackTree->Branch(floatVarName[ivar].Data(),&fTreeVarFloat[ivar],Form("%s/F",floatVarName[ivar].Data()));
  }
  fTrackTree->Branch("TrPar1","AliExternalTrackParam",&fTrPar1,16000,0);
  fTrackTree->Branch("TrPar2","AliExternalTrackParam",&fTrPar2,16000,0);
  fTrackTree->Branch("TrPar3","AliExternalTrackParam",&fTrPar3,16000,0);
  fTrackTree->Branch("TrParV0","AliNeutralTrackParam",&fTrParV0,16000,0);
  fTrackTree->Branch("RecoPrimaryVtx","AliAODVertex",&fPVertexTrk,16000,0);

  PostData(1,fOutput);
  PostData(2,fTrackTree);
}

//________________________________________________________________________
void AliAnalysisTaskCharmDecayTracks::UserExec(Option_t */*option*/){
  /// 
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
  }
  if(!aod){
    printf("AliAnalysisTaskCharmDecayTracks::UserExec: AOD not found!\n");
    return;
  }
  
  // Load all the branches of the DeltaAOD - needed for SelectionBit counting
  TClonesArray *arrayD0toKpi  =0;
  TClonesArray *array3Prong   =0;
  TClonesArray *arrayV0Bachel   =0;
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*)
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {

      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();

      array3Prong  =(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
      arrayD0toKpi =(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
      arrayV0Bachel=(TClonesArray*)aodFromExt->GetList()->FindObject("CascadesHF");
    }
  } else if(aod) {
    array3Prong  =(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    arrayD0toKpi =(TClonesArray*)aod->GetList()->FindObject("D0toKpi");
    arrayV0Bachel=(TClonesArray*)aod->GetList()->FindObject("CascadesHF");
  }

  if(!aod || !array3Prong || !arrayD0toKpi) {
    printf("AliAnalysisTaskCharmDecayTracks::UserExec: AOD branch not found!\n");
    return;
  }
  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex()||TMath::Abs(aod->GetMagneticField())<0.001) return;


  // Reject events with trigger mask 0 of the LHC13d3 production
  // For these events the ITS layers are skipped in the trakcing
  // and the vertex reconstruction efficiency from tracks is biased
  if(fReadMC){
    Int_t runnumber = aod->GetRunNumber();
    if(aod->GetTriggerMask()==0 &&
       (runnumber>=195344 && runnumber<=195677)){
      return;
    }
  }
  
  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;
  if(fReadMC){
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      printf("AliAnalysisTaskCharmDecayTracks::UserExec: MC particles branch not found!\n");
      return;
    }
    
    // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskCharmDecayTracks::UserExec: MC header branch not found!\n");
      return;
    }
  }

  fHistNEvents->Fill(0); // count event
  if(fUsePhysSel){
    Bool_t isPhysSel = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
    if(!isPhysSel) return;
  }
  fHistNEvents->Fill(1);

  fPVertexTrk = aod->GetPrimaryVertex();
  const AliVVertex* vtSPD = aod->GetPrimaryVertexSPD();
  TString titTrc=fPVertexTrk->GetTitle();
  if(titTrc.IsNull() || titTrc=="vertexer: 3D" || titTrc=="vertexer: Z") return;
  if (vtSPD->GetNContributors()<1) return;
  fHistNEvents->Fill(2);

  double covTrc[6],covSPD[6];
  fPVertexTrk->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);
  double dz = fPVertexTrk->GetZ()-vtSPD->GetZ();
  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
  double errTrc = TMath::Sqrt(covTrc[5]);
  double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
  if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return; // bad vertexing
  fHistNEvents->Fill(3);

  Float_t zvert=fPVertexTrk->GetZ();
  if(TMath::Abs(zvert)>10) return;
  fHistNEvents->Fill(4);

  if(fUsePileupCut){
    AliAnalysisUtils utils;
    utils.SetMinPlpContribMV(5);
    utils.SetMaxPlpChi2MV(5.);
    utils.SetMinWDistMV(15.);
    utils.SetCheckPlpFromDifferentBCMV(kTRUE);
    Bool_t isPUMV = utils.IsPileUpMV(aod);
    if(isPUMV) return;
    fHistNEvents->Fill(5);
  }

  // Post the data already here
  PostData(1,fOutput);
  PostData(2,fTrackTree);

  if(fMethod==0){
    fHistNCand->Fill(0);
    for(Int_t i=0; i<kMaxLabel; i++) fMapTrLabel[i]=-999;
    for(Int_t i=0; i<kMaxLabel; i++) fMapV0Label[i]=-999;
    MapTrackLabels(aod);
    MapV0Labels(aod,arrayMC);

    Int_t nDauTr=2;
    if(fSelSpecies==411 || fSelSpecies==431 || (fSelSpecies==4122 && fDecayMode==0)) nDauTr=3;
    fTrPar1.Reset();
    fTrPar2.Reset();
    fTrPar3.Reset();
    fTrParV0.Reset();
    for (Int_t iPart=0; iPart<arrayMC->GetEntriesFast(); iPart++) {
      AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(iPart));
      if (!mcPart) continue;
      Int_t absPdgCode=TMath::Abs(mcPart->GetPdgCode());
      if(absPdgCode!=fSelSpecies) continue;
      fHistNCand->Fill(2);
      Int_t retCode=-1;
      Int_t arrayDauLabels[4]={-1,-1,-1,-1};
      if(fSelSpecies==411) retCode=AliVertexingHFUtils::CheckDplusDecay(arrayMC,mcPart,arrayDauLabels);
      else if(fSelSpecies==421) retCode=AliVertexingHFUtils::CheckD0Decay(arrayMC,mcPart,arrayDauLabels);
      else if(fSelSpecies==431) retCode=AliVertexingHFUtils::CheckDsDecay(arrayMC,mcPart,arrayDauLabels);
      else if(fSelSpecies==4122 && fDecayMode==0) retCode=AliVertexingHFUtils::CheckLcpKpiDecay(arrayMC,mcPart,arrayDauLabels);
      else if(fSelSpecies==4122 && fDecayMode==1){
        retCode=AliVertexingHFUtils::CheckLcV0bachelorDecay(arrayMC,mcPart,arrayDauLabels);
        if(retCode!=1) continue; // reject lambda+pion decays
        Int_t labK0s=-1;
        Int_t labbachelor=-1;
        for(Int_t jd=0; jd<3; jd++){
          Int_t labTr=arrayDauLabels[jd];
          AliAODMCParticle* mcDauPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(labTr));
          if(mcDauPart){
            Int_t labMoth=mcDauPart->GetMother();
            if(labMoth>=0){
              AliAODMCParticle* mcMothPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(labMoth));
              if(mcMothPart){
                Int_t pdgMoth=TMath::Abs(mcMothPart->GetPdgCode());
                if(pdgMoth==310){
                  labK0s=labMoth;
                }else if(pdgMoth==4122){
                  // bachelor
                  labbachelor=labTr;
                }
              }
            }
          }
        }
        if(labbachelor<0) continue;
        arrayDauLabels[0]=labbachelor;
        arrayDauLabels[1]=labK0s;
        arrayDauLabels[2]=-1;
        arrayDauLabels[3]=-1;
        nDauTr=1;
      }
      if(retCode<0 || arrayDauLabels[0]==-1) continue;
      fHistNCand->Fill(3);
      Bool_t fillTree=PrepareTreeVars(mcPart,arrayMC,mcHeader);
      if(fillTree) fHistNCand->Fill(4);
      for(Int_t jd=0; jd<nDauTr; jd++){
        Int_t labTr=arrayDauLabels[jd];
        Int_t idTr=fMapTrLabel[labTr];
        if(idTr<0 || idTr>=aod->GetNumberOfTracks()){
          if(fillTree) fHistNCand->Fill(5);
          fillTree=kFALSE;
          continue;
        }
        AliAODTrack* track = dynamic_cast<AliAODTrack*>(aod->GetTrack(idTr));
        if(!track){
          if(fillTree) fHistNCand->Fill(6);
          fillTree=kFALSE;
          continue;
        }
        if(TMath::Abs(track->GetLabel())!=labTr){
          if(fillTree) fHistNCand->Fill(7);
          fillTree=kFALSE;
          continue;
        }
        if(!IsTrackSelected(track)){
          if(fillTree) fHistNCand->Fill(8);
          fillTree=kFALSE;
        }
        AliAODMCParticle* mcDauPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(labTr));
        if(mcDauPart) fTreeVarInt[2+jd]=mcDauPart->GetPdgCode();
        if(jd==0) fTrPar1.CopyFromVTrack(track);
        else if(jd==1) fTrPar2.CopyFromVTrack(track);
        else if(jd==2) fTrPar3.CopyFromVTrack(track);
      }
      if(fSelSpecies==4122 && fDecayMode==1){
        Int_t labK0=arrayDauLabels[1];
        Int_t idK0=fMapV0Label[labK0];
        if(idK0<0 || idK0>=aod->GetNumberOfV0s()){
          if(fillTree) fHistNCand->Fill(5);
          fillTree=kFALSE;
          continue;
        }
        AliAODv0 *v0 = aod->GetV0(idK0);
        if(!v0){
          if(fillTree) fHistNCand->Fill(6);
          fillTree=kFALSE;
          continue;
        }
        Int_t labV0=v0->MatchToMC(310,arrayMC);
        if(labV0!=labK0){
          if(fillTree) fHistNCand->Fill(7);
          fillTree=kFALSE;
          continue;
        }
        if(!IsV0Selected(v0)){
          if(fillTree) fHistNCand->Fill(8);
          fillTree=kFALSE;
        }
        AliAODMCParticle* mcDauPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(labK0));
        if(mcDauPart) fTreeVarInt[3]=mcDauPart->GetPdgCode();
        Double_t xyz[3],pxpypz[3],cv[21];
        v0->GetXYZ(xyz);
        pxpypz[0]=v0->Px();
        pxpypz[1]=v0->Py();
        pxpypz[2]=v0->Pz();
        v0->GetCovarianceXYZPxPyPz(cv);
        fTrParV0=AliNeutralTrackParam(xyz,pxpypz,cv,0);
      }
      if(fillTree){
        fTrackTree->Fill();
        fHistNCand->Fill(9);
      }
    }
  }else{
    fHistNCand->Fill(1);
    // vHF object is needed to call the method that refills the missing info of the candidates
    // if they have been deleted in dAOD reconstruction phase
    // in order to reduce the size of the file
    AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
    
    TClonesArray *arrayDcand=arrayD0toKpi;
    if(fSelSpecies==411 || fSelSpecies==431 || (fSelSpecies==4122 && fDecayMode==0)) arrayDcand=array3Prong;
    if(fSelSpecies==4122 && fDecayMode==1) arrayDcand=arrayV0Bachel;
    Int_t nCand=arrayDcand->GetEntriesFast();
    Int_t pdg0[2]={321,211};
    Int_t pdgp[3]={321,211,211};
    Int_t pdgs[3]={321,211,321};
    Int_t pdgl[3]={2212,321,211};
    Int_t pdglvb[2]={2212,310}; // always 1st bachelor, 2nd V0
    Int_t pdgv0[2]={211,211};

    for (Int_t iCand = 0; iCand < nCand; iCand++) {
      AliAODRecoDecayHF *d=(AliAODRecoDecayHF*)arrayDcand->UncheckedAt(iCand);
      if(!d) continue;
      Int_t labD=-999;
      Int_t nDauTr=2;
      fTrPar1.Reset();
      fTrPar2.Reset();
      fTrPar3.Reset();
      fTrParV0.Reset();
      if(fSelSpecies==421){
        if(!vHF->FillRecoCand(aod,(AliAODRecoDecayHF2Prong*)d))continue;
        labD = d->MatchToMC(421,arrayMC,2,pdg0);
      }else if(fSelSpecies==411 || fSelSpecies==431 ||  fSelSpecies==4122){
        nDauTr=3;
        if(!vHF->FillRecoCand(aod,(AliAODRecoDecayHF3Prong*)d))continue;
        if(fSelSpecies==411){
          if(!d->HasSelectionBit(AliRDHFCuts::kDplusCuts)) continue;
          labD=d->MatchToMC(411,arrayMC,3,pdgp);
        }else if(fSelSpecies==431){
          if(!d->HasSelectionBit(AliRDHFCuts::kDsCuts)) continue;
          labD = d->MatchToMC(431,arrayMC,3,pdgs);
        }else if(fSelSpecies==4122 && fDecayMode==0){
          if(!d->HasSelectionBit(AliRDHFCuts::kLcCuts)) continue;
          labD = d->MatchToMC(4122,arrayMC,3,pdgl);
        }else if(fSelSpecies==4122 && fDecayMode==1){
          AliAODRecoCascadeHF* dcasc = dynamic_cast<AliAODRecoCascadeHF*>(arrayDcand->UncheckedAt(iCand));
          if(!dcasc) continue;
          if(!dcasc->CheckCascadeFlags()) continue;
          labD = dcasc->MatchToMC(4122,pdglvb[1],pdglvb,pdgv0,arrayMC,kTRUE);
          nDauTr=1;
        }
      }
      if(labD>=0){ 
        fHistNCand->Fill(10);
        AliAODMCParticle *partD = (AliAODMCParticle*)arrayMC->At(labD);
        if(!partD) continue;
        Bool_t fillTree=PrepareTreeVars(partD,arrayMC,mcHeader);
        if(fillTree) fHistNCand->Fill(11);
        for(Int_t jd=0; jd<nDauTr; jd++){
          AliAODTrack* track = (AliAODTrack*)d->GetDaughter(jd);
          if(!track){
            if(fillTree) fHistNCand->Fill(12);
            fillTree=kFALSE;
          }else{
            if(!IsTrackSelected(track)){
              if(fillTree) fHistNCand->Fill(13);
              fillTree=kFALSE;
            }
            Int_t labTr=TMath::Abs(track->GetLabel());
            AliAODMCParticle* mcDauPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(labTr));
            if(mcDauPart) fTreeVarInt[2+jd]=mcDauPart->GetPdgCode();
            if(jd==0) fTrPar1.CopyFromVTrack(track);
            else if(jd==1) fTrPar2.CopyFromVTrack(track);
            else if(jd==2) fTrPar3.CopyFromVTrack(track);
          }
        }
        if(fSelSpecies==4122 && fDecayMode==1){
          AliAODv0* v0=dynamic_cast<AliAODv0*>(((AliAODRecoCascadeHF*)d)->Getv0());
          if(!v0){
            if(fillTree) fHistNCand->Fill(12);
            fillTree=kFALSE;
          }else{
            if(!IsV0Selected(v0)){
              if(fillTree) fHistNCand->Fill(13);
              fillTree=kFALSE;
            }
            Int_t labV0=v0->MatchToMC(310,arrayMC,2,pdgv0);
            if(labV0>=0){
              AliAODMCParticle* mcDauPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(labV0));
              if(mcDauPart) fTreeVarInt[3]=mcDauPart->GetPdgCode();
            }else{
              fTreeVarInt[3]=-1;
              fillTree=kFALSE;
            }
            Double_t xyz[3],pxpypz[3],cv[21];
            v0->GetXYZ(xyz);
            pxpypz[0]=v0->Px();
            pxpypz[1]=v0->Py();
            pxpypz[2]=v0->Pz();
            v0->GetCovarianceXYZPxPyPz(cv);
            fTrParV0=AliNeutralTrackParam(xyz,pxpypz,cv,0);
          }
        }
        if(fillTree){
          fTrackTree->Fill();
          fHistNCand->Fill(14);
        }
      }
    } 
    delete vHF;
  }
  
  PostData(1,fOutput);
  PostData(2,fTrackTree);
    
  return;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCharmDecayTracks::IsTrackSelected(AliAODTrack* track){
  /// track selection cuts

  if(!track) return kFALSE;
  if(track->Charge()==0) return kFALSE;
  if(track->GetID()<0&&!fKeepNegID)return kFALSE;
  if(fFilterMask>0){
    if(!(track->TestFilterMask(fFilterMask))) return kFALSE;
  }
  if(fTrCuts && !fTrCuts->IsSelected(track)) return kFALSE;
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskCharmDecayTracks::IsV0Selected(AliAODv0* v0){
  /// track selection cuts

  if (!v0) return kFALSE;
  Bool_t onFlyStatus=v0->GetOnFlyStatus();
  if(onFlyStatus==kTRUE) return kFALSE;
  AliAODTrack *pTrack=(AliAODTrack *)v0->GetDaughter(0); //0->Positive Daughter
  AliAODTrack *nTrack=(AliAODTrack *)v0->GetDaughter(1); //1->Negative Daughter
  if(!pTrack || !nTrack) return kFALSE;
  if(pTrack->GetID()<0 || nTrack->GetID()<0) return kFALSE;
  if(pTrack->Charge() == nTrack->Charge()) return kFALSE;
  return kTRUE;
}



//_________________________________________________________________
void AliAnalysisTaskCharmDecayTracks::Terminate(Option_t */*option*/)
{
  /// Terminate analysis
  //
  if(fDebug > 1) printf("AliAnalysisTaskCharmDecayTracks: Terminate() \n");
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("hNEvents"));
  if(fHistNEvents){
    printf("Number of analyzed events = %d\n",(Int_t)fHistNEvents->GetBinContent(2));
  }else{
    printf("ERROR: fHistNEvents not available\n");
    return;
  }
  return;
}
//________________________________________________________________________
void AliAnalysisTaskCharmDecayTracks::MapTrackLabels(AliAODEvent* aod){
  /// Fill array of correspondence track lables <-> id
  //
  Int_t nTracks=aod->GetNumberOfTracks();

  for(Int_t it=0; it<nTracks; it++) {
    AliAODTrack *tr=dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
    if(!tr) continue;
    if(tr->GetID()<0) continue;
    if(tr->GetStatus()&AliESDtrack::kITSpureSA) continue;
    if(!(tr->GetStatus()&AliESDtrack::kITSin)) continue;
    Int_t lab=TMath::Abs(tr->GetLabel());
    if(lab<kMaxLabel){
      if(fMapTrLabel[lab]>=0) continue; // tracks with this label were already found
      Int_t countSplit=1;
      Int_t ntpclu=tr->GetTPCncls();
      Int_t nitsclu=tr->GetITSNcls();
      Int_t itBest=it;
      Double_t mom=tr->P();
      for(Int_t it2=it+1; it2<nTracks; it2++) {
        AliAODTrack *tr2=dynamic_cast<AliAODTrack*>(aod->GetTrack(it2));
        if(!tr2) continue;
        if(tr2->GetID()<0) continue;
        if(tr2->GetStatus()&AliESDtrack::kITSpureSA) continue;
        if(!(tr2->GetStatus()&AliESDtrack::kITSin)) continue;
        Int_t lab2=TMath::Abs(tr2->GetLabel());
        Int_t ntpclu2=tr2->GetTPCncls();
        Int_t nitsclu2=tr2->GetITSNcls();
        if(lab2==lab){
          if(countSplit==1){
            fHistMomDupLab->Fill(countSplit,tr->P());
            fHistCluTPCDupLab->Fill(countSplit,tr->GetTPCncls());
          }
          countSplit++;
          fHistCluTPCDupLab->Fill(countSplit,ntpclu2);
          fHistMomDupLab->Fill(countSplit,tr2->P());
          fHistCluTPCDupLabCorrel->Fill(ntpclu,ntpclu2);
          if(ntpclu2>=ntpclu) fHistCluITSDupLabCorrel->Fill(nitsclu,nitsclu2);
          else fHistCluITSDupLabCorrel->Fill(nitsclu2,nitsclu);
          // cases of two tracks with same label and similar number of TPC and ITS clusters are mainly loopers
          // we keep the leg wit higher total momentum, which should be the primary leg
          if(tr2->P()>mom){
            mom=tr2->P();
            itBest=it2;
          }
        }
      }
      fHistTrLab->Fill(countSplit);
      fMapTrLabel[lab]=itBest;
    }else{
      fHistTrLab->Fill(-1);
      printf("Label %d exceeds upper limit\n",lab);
    }
  }
  return;
}
//________________________________________________________________________
void AliAnalysisTaskCharmDecayTracks::MapV0Labels(AliAODEvent* aod, TClonesArray* arrayMC){
  /// Fill array of correspondence V0 lables <-> id
  //
  Int_t nv0s = aod->GetNumberOfV0s();
  for (Int_t iV0 = 0; iV0 < nv0s; iV0++){
    AliAODv0 *v0 = aod->GetV0(iV0);
    if (!v0) continue;
    if(!IsV0Selected(v0)) continue;
    Int_t lab=v0->MatchToMC(310,arrayMC);
    if(lab>=0 && lab<kMaxLabel){
      if(fMapV0Label[lab]>=0) continue; // V0 with this label were already found
      fMapV0Label[lab]=iV0;
    }
  }
  return;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskCharmDecayTracks::PrepareTreeVars(AliAODMCParticle* partD, TClonesArray* arrayMC, AliAODMCHeader* mcHeader){
  /// fill MC truth info in the tree
  
  for(Int_t j=0; j<kNumOfIntVar; j++) fTreeVarInt[j] = 0;
  for(Int_t j=0; j<kNumOfFloatVar; j++) fTreeVarFloat[j]=-9999.;
  if(!partD) return kFALSE;
  Int_t pdgCode=partD->GetPdgCode();
  Int_t orig=AliVertexingHFUtils::CheckOrigin(arrayMC,partD,kTRUE);
  Double_t ptgen=partD->Pt();
  Double_t phigen=partD->Phi();
  Double_t ygen=partD->Y();
  Double_t xori=partD->Xv();
  Double_t yori=partD->Yv();
  Double_t zori=partD->Zv();
  
  Int_t iDau=partD->GetDaughterFirst();
  if(iDau<0) return kFALSE;
  AliAODMCParticle *dauD=(AliAODMCParticle*)arrayMC->At(iDau);
  if(!dauD) return kFALSE;
  Double_t xdec=dauD->Xv();
  Double_t ydec=dauD->Yv();
  Double_t zdec=dauD->Zv();
  
  fTreeVarInt[0] = pdgCode;
  fTreeVarInt[1] = orig;
  fTreeVarFloat[0] = ptgen;
  fTreeVarFloat[1] = phigen;
  fTreeVarFloat[2] = ygen;
  fTreeVarFloat[3] = mcHeader->GetVtxX();
  fTreeVarFloat[4] = mcHeader->GetVtxY();
  fTreeVarFloat[5] = mcHeader->GetVtxZ();
  fTreeVarFloat[6] = xori;
  fTreeVarFloat[7] = yori;
  fTreeVarFloat[8] = zori;
  fTreeVarFloat[9] = xdec;
  fTreeVarFloat[10] = ydec;
  fTreeVarFloat[11] = zdec;

  return kTRUE;
}
