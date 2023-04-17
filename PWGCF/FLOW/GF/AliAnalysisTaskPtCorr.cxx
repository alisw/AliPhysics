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

#include "AliAnalysisTaskPtCorr.h"
#include "AliEventCuts.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "TList.h"
#include "TProfile.h"
#include "AliMultSelection.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenHepMCEventHeader.h"

using namespace std;
using namespace TMath;

ClassImp(AliAnalysisTaskPtCorr)

AliAnalysisTaskPtCorr::AliAnalysisTaskPtCorr() : AliAnalysisTaskSE(),
    fEventCuts(),
    fCentEst(0),
    fRunNo(0),
    fSystFlag(0),
    fContSubfix(0),
    fIsMC(kFALSE),
    fMCEvent(0),
    fCorrList(0),
    fQAList(0),
    fEfficiencyList(0),
    fEfficiencies(0),
    fPowerEfficiencies(0),
    fWeightSubfix(""),
    fGFWSelection(0),
    fGFWnTrackSelection(0),
    fV0MAxis(0),
    fMultiAxis(0),
    fMultiBins(0), 
    fNMultiBins(0),
    fPtAxis(0), 
    fPtBins(0), 
    fNPtBins(0),
    fEta(0.8),
    fEtaNch(0.8),
    fEtaGap(-1),
    fPUcut(1500),
    fRndm(0),
    fNbootstrap(10),
    fUseWeightsOne(false),
    fUseRecNchForMC(false),
    fPileupOff(false),
    fUseNch(false),
    fUseTPConly(0),
    fConstEff(0.8),
    fSigmaEff(0.05),
    mpar(6),
    wp(0),
    wpP(0),
    wpN(0),
    fEventWeight(PtSpace::kWperms),
    fV0MMulti(0),
    fpt(0),
    fEventCount(0),
    fNchTrueVsRec(0),
    fV0MvsMult(0),
    fPtMoms(0),
    fPtDistB(0),
    fPtDistA(0),
    fPtDistC(0),
    fPtDCA(0),
    fPtVsNTrk(0),
    fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
    fOnTheFly(false),
    fImpactParameter(0),
    EventNo(0)
{};
//_____________________________________________________________________________
AliAnalysisTaskPtCorr::AliAnalysisTaskPtCorr(const char *name, bool IsMC, bool isOnTheFly, unsigned int fl_eff, TString ContSubfix) : AliAnalysisTaskSE(name),
    fEventCuts(),
    fCentEst(0),
    fRunNo(0),
    fSystFlag(0),
    fContSubfix(0),
    fIsMC(IsMC),
    fMCEvent(0),
    fCorrList(0),
    fQAList(0),
    fEfficiencyList(0),
    fEfficiencies(0),
    fPowerEfficiencies(0),
    fWeightSubfix(""),
    fGFWSelection(0),
    fGFWnTrackSelection(0),
    fV0MAxis(0),
    fMultiAxis(0),
    fMultiBins(0), 
    fNMultiBins(0),
    fPtAxis(0), 
    fPtBins(0), 
    fNPtBins(0),
    fEta(0.8),
    fEtaNch(0.8),
    fEtaGap(-1),
    fPUcut(1500),
    fRndm(0),
    fNbootstrap(10),
    fUseWeightsOne(false),
    fUseRecNchForMC(false),
    fPileupOff(false),
    fUseNch(false),
    fUseTPConly(0),
    fConstEff(0.8),
    fSigmaEff(0.05),
    mpar(6),
    wp(0),
    wpP(0),
    wpN(0),
    fEventWeight(PtSpace::kWperms),
    fV0MMulti(0),
    fpt(0),
    fEventCount(0),
    fNchTrueVsRec(0),
    fV0MvsMult(0),
    fPtMoms(0),
    fPtDistB(0),
    fPtDistA(0),
    fPtDistC(0),
    fPtDCA(0),
    fPtVsNTrk(0),
    fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
    fOnTheFly(isOnTheFly),
    fImpactParameter(0),
    EventNo(0)
{ 
    SetContSubfix(ContSubfix);
    fCentEst = new TString("V0M");
    SetEffFlags(fl_eff);
    if(!(fIsMC&&!((eff_flags&PtCorrFlags::realeffin)==PtCorrFlags::realeffin)) &&!fOnTheFly)
    {
        DefineInput(1, TList::Class());
    }
    if(fOnTheFly)
    {
      vector<double> b = {0.0,3.72,5.23,7.31,8.88,10.20,11.38,12.47,13.50,14.51,100.0};
      vector<double> cent = {0.0,5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,100.0};
      for(size_t i(0); i<b.size(); ++i) centralitymap[b[i]]=cent[i];
    }
    DefineOutput(1, TList::Class());
    DefineOutput(2,TList::Class());
};
//_____________________________________________________________________________
AliAnalysisTaskPtCorr::~AliAnalysisTaskPtCorr()
{
    if(fCorrList) delete fCorrList;
};
void AliAnalysisTaskPtCorr::NotifyRun() {
    if(fOnTheFly) return;
    Bool_t dummy = fEventCuts.AcceptEvent(InputEvent());
    if(!fPileupOff) fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);

    //Then override PU cut if required:
    if(fGFWSelection->GetSystFlagIndex()==22 && !fPileupOff)
      fEventCuts.fESDvsTPConlyLinearCut[0] = fPUcut;
}
void AliAnalysisTaskPtCorr::UserCreateOutputObjects()
{
    if(!fGFWSelection) SetSystFlag(0);
    if(fUseTPConly) fGFWSelection->fFilterBit = 128;
    fGFWSelection->PrintSetup();
    fSystFlag = fGFWSelection->GetSystFlagIndex();
    if(fGFWSelection->GetSystFlagIndex() == 20) SetCentralityEstimator("CL0");
    else if(fGFWSelection->GetSystFlagIndex() == 21) SetCentralityEstimator("CL1");
    OpenFile(1);
    const int temp_NV0MBinsDefault = 10;
    double temp_V0MBinsDefault[11] = {0,5,10,20,30,40,50,60,70,80,90}; //Last bin to include V0M beyond anchor point
    if(!fV0MAxis) SetV0MBins(temp_NV0MBinsDefault,temp_V0MBinsDefault);
    double *l_V0MBinsDefault=GetBinsFromAxis(fV0MAxis);
    int l_NV0MBinsDefault=fV0MAxis->GetNbins();
    if(!fMultiAxis)
    {
      printf("Multiplicity axis not set. Using defaults bins\n"); 
      SetMultiplicityBins(l_NV0MBinsDefault,l_V0MBinsDefault);
    }
    fMultiBins = GetBinsFromAxis(fMultiAxis);
    fNMultiBins = fMultiAxis->GetNbins();
    const int l_NPtBinsDefault = 31;
    Double_t l_PtBinsDefault[l_NPtBinsDefault+1] = {0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0};   
    if(!fPtAxis) SetPtBins(l_NPtBinsDefault,l_PtBinsDefault);
    fPtBins = GetBinsFromAxis(fPtAxis);
    fNPtBins = fPtAxis->GetNbins();
    fV0MMulti = new TH1D("V0M_Multi","V0M_Multi",l_NV0MBinsDefault,l_V0MBinsDefault);
    fRndm = new TRandom(0);
    if(!(fIsMC&&!((eff_flags&PtCorrFlags::realeffin)==PtCorrFlags::realeffin)) &&!fOnTheFly)
    { 
        fEfficiencyList = (TList*)GetInputData(1);
        fEfficiencies = new TH1D*[l_NV0MBinsDefault];
        fPowerEfficiencies = new TH2D*[l_NV0MBinsDefault];
        for(int i=0;i<l_NV0MBinsDefault;i++) {
          if(eff_flags&PtCorrFlags::powereff)
          {
            fPowerEfficiencies[i] = (TH2D*)fEfficiencyList->FindObject(Form("Eff_Cent%i%s",i,fGFWSelection->GetSystPF()));
            if(!fPowerEfficiencies[i])
            {
              if(!i) AliFatal("Could not fetch efficiency!\n");
              printf("Could not find efficiency for V0M bin no. %i! Cloning the previous efficiency instead...\n",i);
              fPowerEfficiencies[i] = (TH2D*)fPowerEfficiencies[i-1]->Clone(Form("Eff_Cent%i%s",i,fGFWSelection->GetSystPF()));
            }
          }
          else
          {
            fEfficiencies[i] = (TH1D*)fEfficiencyList->FindObject(Form("EffRescaled_Cent%i%s",i,fGFWSelection->GetSystPF()));
            if(!fEfficiencies[i]) {
            if(!i) AliFatal("Could not fetch efficiency!\n");
            printf("Could not find efficiency for V0M bin no. %i! Cloning the previous efficiency instead...\n",i);
            fEfficiencies[i] = (TH1D*)fEfficiencies[i-1]->Clone(Form("EffRescaled_Cent%i%s",i,fGFWSelection->GetSystPF()));
            };
          }
        }
    };
    
    fCorrList = new TList(); fCorrList->SetOwner(1);
    fpt = new AliPtContainer("ptcont","ptcont",fNMultiBins,fMultiBins,mpar,fEtaGap>=0);
    fCorrList->Add(fpt);
    if(fNbootstrap) {
      fpt->InitializeSubsamples(fNbootstrap);
    }
    fpt->SetEventWeight(fEventWeight);
    fCorrList->Add(fV0MMulti);
    if(fIsMC) {
      fNchTrueVsRec = new TH2D("NchTrueVsRec",";Nch (MC-true); Nch (MC-reco)",fNMultiBins,fMultiBins,fNMultiBins,fMultiBins);
      fCorrList->Add(fNchTrueVsRec);
    }
    fV0MvsMult = new TH2D("MultVsV0M","MultVsV0M",103,0,103,fNMultiBins,fMultiBins[0],fMultiBins[fNMultiBins]);
    fCorrList->Add(fV0MvsMult);
    PostData(1,fCorrList);
    fQAList = new TList();
    fQAList->SetOwner(1);
    fEventCuts.AddQAplotsToList(fQAList,kTRUE);
    int nEventCutLabel = 6; 
    fEventCount = new TH1D("fEventCount","Event counter",nEventCutLabel,0,nEventCutLabel);
    TString eventCutLabel[6]={"Input","Centrality","Trigger","AliEventCuts","Vertex","Tracks"};
    for(int i=0;i<nEventCutLabel;++i) fEventCount->GetXaxis()->SetBinLabel(i+1,eventCutLabel[i].Data());
    fQAList->Add(fEventCount);
    double powers[7] = {0,1,2,3,4,5,6};
    int Npows = 6;
    fPtMoms = new TH3D("ptMoments","ptMoments",fNPtBins,fPtBins,Npows,powers,temp_NV0MBinsDefault,l_V0MBinsDefault);
    fQAList->Add(fPtMoms);
    fPtDistB = new TH1D("ptDistB","Before;p_{T};counts",fNPtBins,fPtBins);
    fQAList->Add(fPtDistB);
    fPtDistA = new TH1D("ptDistA","After;p_{T};counts",fNPtBins,fPtBins);
    fQAList->Add(fPtDistA);
    fPtDistC = new TH1D("ptDistC","Corrected;p_{T};counts",fNPtBins,fPtBins);
    fQAList->Add(fPtDistC);
    double binsDCA[61] = {-3.00, -2.90, -2.80, -2.70, -2.60, -2.50, -2.40, -2.30, -2.20, -2.10, -2.00, -1.90, -1.80, -1.70, -1.60, -1.50, -1.40, -1.30, -1.20, -1.10, -1.00, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00};
    int NbinsDCA = 60;
    fPtDCA = new TH3D("ptDCA","ptDCA;pt;dcaxy;dcaz",fNPtBins,fPtBins,NbinsDCA,binsDCA,NbinsDCA,binsDCA);
    fQAList->Add(fPtDCA);
    Int_t NNtrkBins=100;
    Double_t *binsNtrk = new Double_t[NNtrkBins+1];
    for(Int_t i=0;i<=NNtrkBins; i++) binsNtrk[i] = 30*i+0.5;
    fPtVsNTrk = new TH2D("NtracksPt",";NtracksPt;N_{trk};p_{T} (GeV/#it{c})",NNtrkBins,binsNtrk,fNPtBins,fPtBins);
    fQAList->Add(fPtVsNTrk);
    PostData(2,fQAList);
    fEventCuts.OverrideAutomaticTriggerSelection(fTriggerType,true);
    fGFWnTrackSelection = new AliGFWCuts();
    fGFWnTrackSelection->SetupCuts(0);
    fGFWnTrackSelection->SetEta(fEtaNch);
    printf("User output objects created!\n"); 
}
void AliAnalysisTaskPtCorr::UserExec(Option_t *)
{
  EventNo++;
  if(fOnTheFly) {
    fMCEvent = getMCEvent();
    double l_cent = getCentrality();
    if(l_cent<0) return;
    wp.resize(10,vector<double>(10));
    wpP.resize(10,vector<double>(10));
    wpN.resize(10,vector<double>(10));
    double trackXYZ[3];
    double ptMin = fPtBins[0];
    double ptMax = fPtBins[fNPtBins];
    int NTracks = fMCEvent->GetNumberOfPrimaries();
    if(NTracks < 1) { return; }
    for(Int_t iTrack(0); iTrack < NTracks; iTrack++) 
    {
        AliMCParticle* track = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(iTrack));
        if(!track) { continue; }
        double l_eta = track->Eta();          
        if (TMath::Abs(l_eta) > fEta) continue;
        double l_pt = track->Pt();
        if (l_pt<ptMin || l_pt>ptMax) continue;
        if(fEtaGap >= 0 && l_eta > fEtaGap) FillWPCounter(wpP,1,l_pt);
        if(fEtaGap >= 0 && l_eta < -fEtaGap) FillWPCounter(wpN,1,l_pt);
        FillWPCounter(wp,1,l_pt);
    }
    if(wp[1][0]==0) return;
    fV0MMulti->Fill(l_cent);
    double l_rnd = fRndm->Rndm();
    //Fill recursive pt-correlations
    fpt->FillRecursive(wp,0);
    //Test with explicit ck and Skew calculation
    fpt->FillCk(wp,l_cent,l_rnd);
    fpt->FillSkew(wp,l_cent,l_rnd);
    //Fill subevent profiles with appropriate wp arrays
    if(fEtaGap>=0) {
      fpt->FillRecursive(wpP,1); fpt->FillRecursive(wpN,2);
    }
    fpt->FillRecursiveProfiles(l_cent,l_rnd);
    wp.clear();
    wpP.clear();
    wpN.clear();
    PostData(1,fCorrList);  
  }
  else {
    AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;
    if(fIsMC) {
      fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
      if (!fMCEvent) return;
    }
    fEventCount->Fill("Input",1);
    double l_cent = getCentrality();
    if(l_cent<0) return;
    fEventCount->Fill("Centrality",1);
    if(!CheckTrigger(l_cent)) return;
    fEventCount->Fill("Trigger",1);
    double vtxXYZ[] = {0.,0.,0.};
    if(!AcceptAODEvent(fAOD, vtxXYZ)) return;
    if(!fGFWSelection->AcceptVertex(fAOD)) return;
    fEventCount->Fill("Vertex",1);
    FillPtCorr(fAOD,l_cent,vtxXYZ);
  }
}
bool AliAnalysisTaskPtCorr::CheckTrigger(Double_t lCent) {
  unsigned int fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(!(fTriggerType&fSelMask)) { return kFALSE; }; //printf("Returning from the generic check\n");
  if(fSelMask&(fTriggerType&(AliVEvent::kINT7+AliVEvent::kMB))) {return kTRUE; }; //printf("Passed by MB trigger!\n");
  if((fSelMask&fTriggerType&AliVEvent::kCentral) && lCent>10) {return kFALSE; }; //printf("Returnning from kCent case\n");
  if((fSelMask&fTriggerType&AliVEvent::kSemiCentral) && (lCent<30 || lCent>50)) {return kFALSE; }; //printf("Returning from kSC case\n");
  return kTRUE;
};
double AliAnalysisTaskPtCorr::getCentrality()
{
  if(fOnTheFly)
  {
    vector<double> b = {0.0,3.72,5.23,7.31,8.88,10.20,11.38,12.47,13.50,14.51,100.0};
    vector<double>::iterator it = upper_bound(b.begin(),b.end(),fImpactParameter);
    double l_cent = (centralitymap[b[it-b.begin()]]+centralitymap[b[it-b.begin()-1]])/2.0;
    return l_cent;
  }
  else 
  {
    AliMultSelection *l_MultSel = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
    if(!l_MultSel) { printf("MultSelection not found\n"); return -1.0; }
    double l_cent = l_MultSel->GetMultiplicityPercentile(fCentEst->Data());
    return l_cent;
  }
  return -1.0;
}
void AliAnalysisTaskPtCorr::Terminate(Option_t *)
{

}
AliMCEvent *AliAnalysisTaskPtCorr::getMCEvent() {
  AliMCEvent* ev = dynamic_cast<AliMCEvent*>(MCEvent());
  if(!ev) { AliFatal("MC event not found!"); return 0; }
  AliGenEventHeader *header = dynamic_cast<AliGenEventHeader*>(ev->GenEventHeader());
  if(!header) { AliFatal("MC event not generated!"); return 0; }
  AliCollisionGeometry* headerH;
  TString genName;
  TList *ltgen = (TList*)ev->GetCocktailList();
  if (ltgen) {
  for(auto&& listObject: *ltgen){
    genName = Form("%s",listObject->GetName());
    if (genName.Contains("Hijing")) {
      headerH = dynamic_cast<AliCollisionGeometry*>(listObject);
      break;
      }
    }
  }
  else headerH = dynamic_cast<AliCollisionGeometry*>(ev->GenEventHeader());
  if(headerH){
      fImpactParameter = headerH->ImpactParameter();
  }
  return ev;
}
bool AliAnalysisTaskPtCorr::AcceptAODTrack(AliAODTrack *tr, double *ltrackXYZ, const double &ptMin, const double &ptMax, double *vtxp)
{
  if(tr->Pt()<ptMin) return kFALSE;
  if(tr->Pt()>ptMax) return kFALSE;
  if(ltrackXYZ && vtxp) {
    tr->GetXYZ(ltrackXYZ);
    ltrackXYZ[0] = ltrackXYZ[0]-vtxp[0];
    ltrackXYZ[1] = ltrackXYZ[1]-vtxp[1];
    ltrackXYZ[2] = ltrackXYZ[2]-vtxp[2];
  } else return kFALSE; //DCA cut is a must for now
  return fGFWSelection->AcceptTrack(tr,fSystFlag==1?0:ltrackXYZ,0,kFALSE);
};
bool AliAnalysisTaskPtCorr::AcceptAODTrack(AliAODTrack *mtr, double *ltrackXYZ, const double &ptMin, const double &ptMax, double *vtxp, int &nTot) 
{
  if(mtr->Pt()<ptMin) return kFALSE;
  if(mtr->Pt()>ptMax) return kFALSE;
  if(ltrackXYZ && vtxp) {
    mtr->GetXYZ(ltrackXYZ);
    ltrackXYZ[0] = ltrackXYZ[0]-vtxp[0];
    ltrackXYZ[1] = ltrackXYZ[1]-vtxp[1];
    ltrackXYZ[2] = ltrackXYZ[2]-vtxp[2];
  } else return kFALSE; 
  if(fGFWnTrackSelection->AcceptTrack(mtr,ltrackXYZ,0,kFALSE)) nTot++;
  return fGFWSelection->AcceptTrack(mtr,fSystFlag==1?0:ltrackXYZ,0,kFALSE); 
};
bool AliAnalysisTaskPtCorr::AcceptAODEvent(AliAODEvent *fAOD, Double_t *inVtxXYZ)
{
    if(!fEventCuts.AcceptEvent(fAOD)) return 0; 
    fEventCount->Fill("AliEventCuts",1);
    const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(fAOD->GetPrimaryVertex());
    if(!vtx || vtx->GetNContributors() < 1)
    return kFALSE;
    const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(fAOD->GetPrimaryVertexSPD());
    Double_t dMaxResol = 0.25; // suggested from DPG
    Double_t cov[6] = {0};
    vtxSPD->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if ( vtxSPD->IsFromVertexerZ() && (zRes > dMaxResol)) return kFALSE;
    const Double_t aodVtxZ = vtx->GetZ();
    if(TMath::Abs(aodVtxZ) > 10)
    return kFALSE;
    vtx->GetXYZ(inVtxXYZ);
    return kTRUE;
};
void AliAnalysisTaskPtCorr::FillWPCounter(vector<vector<double>> &inarr, double w, double p)
{
  for(int i=0;i<=mpar;++i)
  {
    for(int j=0;j<=mpar;++j)
    {
      inarr[i][j] += pow(w,i)*pow(p,j);
    }
  }
  return;
}
void AliAnalysisTaskPtCorr::FillWPCounter(vector<vector<double>> &inarr, vector<double> w, double p)
{
  for(int i=0;i<=mpar;++i)
  {
    for(int j=0;j<=mpar;++j)
    {
      double ww = (i==j)?w[j]:pow(w[1],i);
      inarr[i][j] += ww*pow(p,j);
    }
  }
  return;
}
void AliAnalysisTaskPtCorr::FillPtCorr(AliAODEvent* fAOD, const double &l_cent, double *vtxXYZ)
{          
    wp.resize(10,vector<double>(10));
    wpP.resize(10,vector<double>(10));
    wpN.resize(10,vector<double>(10));
    AliAODTrack *track;
    double trackXYZ[3];
    double ptMin = fPtBins[0];
    double ptMax = fPtBins[fNPtBins];
    int iCent = fV0MMulti->FindBin(l_cent);
    if(!iCent || iCent>fV0MMulti->GetNbinsX()) return;
    iCent--;
    int nTracks=0;
    double rnd_eff = 0.0;
    if(fIsMC)
    {
      int nTracksMC=0;
      int nTracksRec=0;
      if(fUseRecNchForMC) nTracksRec = GetNTracks(fAOD,ptMin,ptMax,vtxXYZ);
      TClonesArray *tca = (TClonesArray*)fInputEvent->FindListObject("mcparticles");
      Int_t nPrim = tca->GetEntries();
      if(nPrim<1) return;
      AliAODMCParticle *part;
      vector<int> NtracksPt(fNPtBins);
      for(Int_t ipart = 0; ipart < nPrim; ipart++) {
        part = (AliAODMCParticle*)tca->At(ipart);
        if(!part->IsPhysicalPrimary()) continue;
        if(part->Charge()==0.) continue;
        double l_eta = part->Eta();
        if (TMath::Abs(l_eta) > fEta) continue;
        double l_pt = part->Pt();
        fPtDistB->Fill(l_pt);
        if (l_pt<ptMin || l_pt>ptMax) continue;
        double wNUE = 1.0;
        if(eff_flags&(PtCorrFlags::consteff|PtCorrFlags::gausseff)) {
          if(eff_flags&PtCorrFlags::realeffin){
            wNUE = fEfficiencies[iCent]->GetBinContent(fEfficiencies[iCent]->FindBin(l_pt));
            if(wNUE==0.0) continue;
          }
          else wNUE = fConstEff;
          if(eff_flags&PtCorrFlags::consteff) {
            rnd_eff = fRndm->Rndm(); 
            if(rnd_eff > wNUE) continue;
          }
          if(eff_flags&PtCorrFlags::gausseff) { 
            wNUE = fRndm->Gaus(wNUE,fSigmaEff); 
            rnd_eff = fRndm->Rndm(); 
            if(rnd_eff > wNUE) continue;
          }
          wNUE = 1./wNUE;
        }
        if(TMath::Abs(l_eta)<fEtaNch) { nTracksMC++; NtracksPt[fPtDistA->GetXaxis()->FindBin(l_pt)-1]++; }
        if(fEtaGap >= 0 && l_eta > fEtaGap) FillWPCounter(wpP,1,l_pt);
        if(fEtaGap >= 0 && l_eta < -fEtaGap) FillWPCounter(wpN,1,l_pt);
        FillWPCounter(wp,wNUE,l_pt);
        for(int i=0;i<6;++i) fPtMoms->Fill(pow(l_pt,i+1),i+0.5,l_cent);
        fPtDistA->Fill(l_pt);
        fPtDCA->Fill(l_pt,TMath::Sqrt(trackXYZ[0]*trackXYZ[0]+trackXYZ[1]*trackXYZ[1]),trackXYZ[2]);
      }
      nTracks = fUseRecNchForMC?nTracksRec:nTracksMC;
      if(fUseRecNchForMC) fNchTrueVsRec->Fill(nTracksMC,nTracksRec);
      for(int i(0);i<fNPtBins;++i) fPtVsNTrk->Fill(NtracksPt[i],fPtBins[i]+0.001); 
    }
    else
    {
      nTracks=GetNTracks(fAOD,ptMin,ptMax,vtxXYZ);
      AliAODTrack* lTrack;
      for(int iTrack(0); iTrack<fAOD->GetNumberOfTracks();iTrack++)
      {
          track = (AliAODTrack*)fAOD->GetTrack(iTrack);
          if(!track) continue;
          double l_eta = track->Eta();
          double l_pt = track->Pt();
          fPtDistB->Fill(l_pt);
          if(TMath::Abs(l_eta) > fEta) continue;
          double trackXYZ[] = {0.,0.,0.};
          if(!AcceptAODTrack(track,trackXYZ,ptMin,ptMax,vtxXYZ)) continue;
          if(eff_flags & PtCorrFlags::powereff)
          {
            vector<double> wNUE(7,1.0);
            for(int i=0;i<6;++i) {
              wNUE[i+1] = fPowerEfficiencies[iCent]->GetBinContent(fPowerEfficiencies[iCent]->GetXaxis()->FindBin(l_pt),fPowerEfficiencies[iCent]->GetYaxis()->FindBin(i+0.5)); 
              if(wNUE[i+1]==0) continue; 
              wNUE[i+1] = (eff_flags&PtCorrFlags::noeff)?1.0:1.0/wNUE[i+1];
            }
            if(fEtaGap >= 0 && l_eta > fEtaGap) FillWPCounter(wpP,wNUE,l_pt);
            if(fEtaGap >= 0 && l_eta < -fEtaGap) FillWPCounter(wpN,wNUE,l_pt);
            FillWPCounter(wp,wNUE,l_pt);
          }
          else
          {
            double wNUE = fEfficiencies[iCent]->GetBinContent(fEfficiencies[iCent]->FindBin(l_pt));
            if(wNUE==0.0) continue;
            if(eff_flags&PtCorrFlags::flateff) {
              wNUE = fEfficiencies[iCent]->GetMinimum(0.00001)/wNUE;
              rnd_eff = fRndm->Rndm(); 
              if(rnd_eff > wNUE) continue;
            }
            if(eff_flags&PtCorrFlags::noeff) wNUE = 1.0;
            wNUE = 1.0/wNUE;
            
            if(fEtaGap >= 0 && l_eta > fEtaGap) FillWPCounter(wpP,wNUE,l_pt);
            if(fEtaGap >= 0 && l_eta < -fEtaGap) FillWPCounter(wpN,wNUE,l_pt);
            FillWPCounter(wp,wNUE,l_pt);
            fPtDistA->Fill(l_pt);
            fPtDistC->Fill(l_pt,wNUE);
          }          
          for(int i=0;i<6;++i) fPtMoms->Fill(pow(l_pt,i+1),i+0.5,l_cent);
          fPtDCA->Fill(l_pt,TMath::Sqrt(trackXYZ[0]*trackXYZ[0]+trackXYZ[1]*trackXYZ[1]),trackXYZ[2]);
      }
    }
    if(wp[1][0]==0) return;
    fEventCount->Fill("Tracks",1);
    fV0MMulti->Fill(l_cent);
    fV0MvsMult->Fill(l_cent,nTracks);
    double l_mult = fUseNch?(1.0*nTracks):l_cent;
    double l_rnd = fRndm->Rndm();
    //Fill recursive pt-correlations
    fpt->FillRecursive(wp);
    //Test with explicit ck and Skew calculation
    fpt->FillCk(wp,l_mult,l_rnd);
    fpt->FillSkew(wp,l_mult,l_rnd);
    //Fill subevent profiles with appropriate wp arrays
    if(fEtaGap>=0) {
      fpt->FillRecursive(wpP,1); fpt->FillRecursive(wpN,2);
    }
    fpt->FillRecursiveProfiles(l_mult,l_rnd);
    wp.clear();
    wpP.clear();
    wpN.clear();
    PostData(1,fCorrList);
}
int AliAnalysisTaskPtCorr::GetNTracks(AliAODEvent* fAOD, const Double_t &ptmin, const Double_t &ptmax, Double_t *vtxp) 
{
  Double_t ltrackXYZ[3];
  AliAODTrack *lTrack;
  Int_t nTracks=0;
  for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
    if(!lTrack) continue;
    if(!AcceptAODTrack(lTrack,ltrackXYZ,ptmin,ptmax,vtxp,nTracks)) continue;
  };
  return nTracks;
}
double *AliAnalysisTaskPtCorr::GetBinsFromAxis(TAxis *inax) {
  Int_t lBins = inax->GetNbins();
  Double_t *retBins = new Double_t[lBins+1];
  for(Int_t i=0;i<lBins;i++)
    retBins[i] = inax->GetBinLowEdge(i+1);
  retBins[lBins] = inax->GetBinUpEdge(lBins);
  return retBins;
}
void AliAnalysisTaskPtCorr::SetPtBins(int nPtBins, double *PtBins) {
  if(fPtAxis) delete fPtAxis;
  fPtAxis = new TAxis(nPtBins, PtBins);
}
void AliAnalysisTaskPtCorr::SetMultiplicityBins(int nMultiBins, double *multibins) {
  if(fMultiAxis) delete fMultiAxis;
  fMultiAxis = new TAxis(nMultiBins, multibins);
}
void AliAnalysisTaskPtCorr::SetV0MBins(int nMultiBins, double *multibins) {
  if(fV0MAxis) delete fV0MAxis;
  fV0MAxis = new TAxis(nMultiBins, multibins);
}
