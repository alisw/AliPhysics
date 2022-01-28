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

using namespace std;
using namespace TMath;

int AliAnalysisTaskPtCorr::fFactorial[9] = {1,1,2,6,24,120,720,5040,40320};
int AliAnalysisTaskPtCorr::fSign[9] = {1,-1,1,-1,1,-1,1,-1,1};

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
    fEfficiencyList(0),
    fEfficiencies(0),
    fWeightSubfix(""),
    fGFWSelection(0),
    fV0MAxis(0),
    fMultiAxis(0),
    fMultiBins(0), 
    fNMultiBins(0),
    fPtAxis(0), 
    fPtBins(0), 
    fNPtBins(0),
    fEta(0.8),
    fPtMin(0.2),
    fPtMax(3.0),
    fRndm(0),
    fNbootstrap(10),
    fUseWeightsOne(false),
    mpar(6),
    fV0MMulti(0),
    fptcorr(0),
    pfmpt(0),
    corr(7, 0.0),
    sumw(7, 0.0),
    fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
    fOnTheFly(false),
    fImpactParameter(0)
{};
//_____________________________________________________________________________
AliAnalysisTaskPtCorr::AliAnalysisTaskPtCorr(const char *name, bool IsMC, TString ContSubfix) : AliAnalysisTaskSE(name),
    fEventCuts(),
    fCentEst(0),
    fRunNo(0),
    fSystFlag(0),
    fContSubfix(0),
    fIsMC(IsMC),
    fMCEvent(0),
    fCorrList(0),
    fEfficiencyList(0),
    fEfficiencies(0),
    fWeightSubfix(""),
    fGFWSelection(0),
    fV0MAxis(0),
    fMultiAxis(0),
    fMultiBins(0), 
    fNMultiBins(0),
    fPtAxis(0), 
    fPtBins(0), 
    fNPtBins(0),
    fEta(0.8),
    fPtMin(0.2),
    fPtMax(3.0),
    fRndm(0),
    fNbootstrap(10),
    fUseWeightsOne(false),
    mpar(6),
    fV0MMulti(0),
    fptcorr(0),
    pfmpt(0),
    corr(7, 0.0),
    sumw(7, 0.0),
    fTriggerType(AliVEvent::kMB+AliVEvent::kINT7),
    fOnTheFly(false),
    fImpactParameter(0)
{ 
    SetContSubfix(ContSubfix);
    fCentEst = new TString("V0M");
    if(!fIsMC)
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
};
//_____________________________________________________________________________
AliAnalysisTaskPtCorr::~AliAnalysisTaskPtCorr()
{
    if(fCorrList) delete fCorrList;
};
void AliAnalysisTaskPtCorr::NotifyRun() {
    Bool_t dummy = fEventCuts.AcceptEvent(InputEvent());
    fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);

    //Then override PU cut if required:
    if(fGFWSelection->GetSystFlagIndex()==22)
      fEventCuts.fESDvsTPConlyLinearCut[0] = 1500.;
}
void AliAnalysisTaskPtCorr::UserCreateOutputObjects()
{
    if(!fGFWSelection) SetSystFlag(0);
    fGFWSelection->PrintSetup();
    fSystFlag = fGFWSelection->GetSystFlagIndex();
    if(fGFWSelection->GetSystFlagIndex() == 20) SetCentralityEstimator("CL0");
    else if(fGFWSelection->GetSystFlagIndex() == 21) SetCentralityEstimator("CL1");
    OpenFile(1);
    const int temp_NV0MBinsDefault = 10;
    double temp_V0MBinsDefault[12] = {0,5,10,20,30,40,50,60,70,80,90,101}; //Last bin to include V0M beyond anchor point
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
    const int l_NPtBinsDefault = 18;
    Double_t l_PtBinsDefault[l_NPtBinsDefault+1] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,4.5,5.0};   
    if(!fPtAxis) SetPtBins(l_NPtBinsDefault,l_PtBinsDefault);
    fPtBins = GetBinsFromAxis(fPtAxis);
    fNPtBins = fPtAxis->GetNbins();
    fV0MMulti = new TH1D("V0M_Multi","V0M_Multi",l_NV0MBinsDefault,l_V0MBinsDefault);
    fRndm = new TRandom(0);
    corr.resize(mpar+1);
    sumw.resize(mpar+1);
    if(!fIsMC)
    { 
        fEfficiencyList = (TList*)GetInputData(1);
        fEfficiencies = new TH1D*[l_NV0MBinsDefault];
        for(int i=0;i<l_NV0MBinsDefault;i++) {
            fEfficiencies[i] = (TH1D*)fEfficiencyList->FindObject(Form("EffRescaled_Cent%i%s",i,fGFWSelection->GetSystPF()));
            if(!fEfficiencies[i]) {
            if(!i) AliFatal("Could not fetch efficiency!\n");
            printf("Could not find efficiency for V0M bin no. %i! Cloning the previous efficiency instead...\n",i);
            fEfficiencies[i] = (TH1D*)fEfficiencies[i-1]->Clone(Form("EffRescaled_Cent%i%s",i,fGFWSelection->GetSystPF()));
            };
        }
    };
    fCorrList = new TList(); fCorrList->SetOwner(1);
    fptcorr = new AliProfileBS*[9];
    for(int i(0);i<mpar;++i) 
    {
      fCorrList->Add(new AliProfileBS(Form("corr_%ipar",i+1),Form("corr_%ipar",i+1),fNMultiBins,fMultiBins));
      fptcorr[i] = (AliProfileBS*)fCorrList->At(i);
    }
    pfmpt = new AliProfileBS("meanpt","meanpt",fNMultiBins,fMultiBins);
    fCorrList->Add(pfmpt);
    if(fNbootstrap) {
      for(int i(0);i<mpar;++i) fptcorr[i]->InitializeSubsamples(fNbootstrap);             
      pfmpt->InitializeSubsamples(fNbootstrap);
    }
    fCorrList->Add(fV0MMulti);
    PostData(1,fCorrList);
    printf("User output objects created!\n"); 
}
void AliAnalysisTaskPtCorr::UserExec(Option_t *)
{
    AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;
    if(fIsMC) {
        fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
        if (!fMCEvent) return;
    }

    AliVEvent* fEvent;
    if(fOnTheFly)
    {
      fEvent = dynamic_cast<AliVEvent*>(MCEvent());
      if(!fEvent) { printf("Event not found!\n"); return; }
      if(!AcceptMCEvent(fEvent)) return;
    }
    
    double l_cent = getCentrality();
    if(l_cent<0) return;
    if(!fOnTheFly && !CheckTrigger(l_cent)) return;
    double vtxXYZ[] = {0.,0.,0.};
    if(!fOnTheFly && !AcceptAODEvent(fAOD, vtxXYZ)) return;
    double vtxZ = (fOnTheFly)?fEvent->GetPrimaryVertex()->GetZ():fAOD->GetPrimaryVertex()->GetZ();
    if(!fOnTheFly && !fGFWSelection->AcceptVertex(fAOD)) return;
    (fOnTheFly)?FillPtCorr(fEvent,vtxZ,l_cent,vtxXYZ):FillPtCorr(fAOD,vtxZ,l_cent,vtxXYZ);
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
bool AliAnalysisTaskPtCorr::AcceptMCEvent(AliVEvent* inev)
{
  AliMCEvent* ev = dynamic_cast<AliMCEvent*>(inev);
  if(!ev) { AliFatal("MC event not found!\n"); return kFALSE; }

  AliGenEventHeader *header = dynamic_cast<AliGenEventHeader*>(ev->GenEventHeader());
  if(!header) { AliFatal("MC event not generated!\n"); return kFALSE; }
  const AliVVertex *vertex = ev->GetPrimaryVertex();
  if(!ev) { AliError("Vertex of MC not found!\n"); }
  if(TMath::Abs(vertex->GetZ()) > 10) return kFALSE;

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
  else
    headerH = dynamic_cast<AliCollisionGeometry*>(ev->GenEventHeader());
  if(headerH){
      fImpactParameter = headerH->ImpactParameter();
  }

  return kTRUE;
}
bool AliAnalysisTaskPtCorr::AcceptAODTrack(AliAODTrack *tr, Double_t *ltrackXYZ, const Double_t &ptMin, const Double_t &ptMax, Double_t *vtxp)
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
bool AliAnalysisTaskPtCorr::AcceptAODEvent(AliAODEvent *ev, Double_t *inVtxXYZ)
{
    const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(ev->GetPrimaryVertex());
    if(!vtx || vtx->GetNContributors() < 1)
    return kFALSE;
    const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(ev->GetPrimaryVertexSPD());
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
template<typename T> void AliAnalysisTaskPtCorr::FillWPCounter(T& inarr, double w, double pt)
{
  for(int i = 0; i < mpar; ++i)
  {
    inarr[i][0] += pow(w,i+1);
    inarr[i][1] += pow(w,i+1)*pow(pt,i+1);
  }
  return;
}
void AliAnalysisTaskPtCorr::FillPtCorr(AliVEvent* ev, const double &VtxZ, const double &l_cent, double *vtxXYZ)
{
    double wp[8][2];
    AliAODTrack *track;
    double trackXYZ[3];
    double ptMin = fPtBins[0];
    double ptMax = fPtBins[fNPtBins];
    int iCent = fV0MMulti->FindBin(l_cent);
    if(!iCent || iCent>fV0MMulti->GetNbinsX()) return;
    iCent--;
    if(fOnTheFly)
    {
      AliMCEvent* mcev = dynamic_cast<AliMCEvent*>(ev);
      int NTracks = mcev->GetNumberOfPrimaries();
      if(NTracks < 1) { return; }
      for(Int_t iTrack(0); iTrack < NTracks; iTrack++) 
      {
          AliMCParticle* track = dynamic_cast<AliMCParticle*>(mcev->GetTrack(iTrack));
          if(!track) { continue; }

          if(!(mcev->IsPhysicalPrimary(iTrack))) continue;
          if(track->Charge() == 0) continue;
          double l_eta = track->Eta();          
          if (TMath::Abs(l_eta) > fEta) continue;
          double l_pt = track->Pt();
          if (l_pt<0.2 || l_pt>3.) continue;
          FillWPCounter(wp,1,l_pt);
      }    
    }
    else if(fIsMC)
    {
      TClonesArray *tca = (TClonesArray*)fInputEvent->FindListObject("mcparticles");
      Int_t nPrim = tca->GetEntries();
      AliAODMCParticle *part;
      for(Int_t ipart = 0; ipart < nPrim; ipart++) {
        part = (AliAODMCParticle*)tca->At(ipart);
        if (!part->IsPhysicalPrimary()) continue;
        if (part->Charge()==0.) continue;
        double l_eta = part->Eta();
        if (TMath::Abs(l_eta) > fEta) continue;
        double l_pt = part->Pt();
        if (l_pt<0.2 || l_pt>3.) continue;
        FillWPCounter(wp,1,l_pt);
      }
    }
    else
    {
      for(int iTrack(0); iTrack<ev->GetNumberOfTracks();iTrack++)
      {
          track = (AliAODTrack*)ev->GetTrack(iTrack);
          if(!track) continue;
          double l_eta = track->Eta();
          if(l_eta>0.8) continue;
          double trackXYZ[] = {0.,0.,0.};
          if(!AcceptAODTrack(track,trackXYZ,ptMin,ptMax,vtxXYZ)) continue;
          double l_pt = track->Pt();
          double wNUE = fEfficiencies[iCent]->GetBinContent(fEfficiencies[iCent]->FindBin(l_pt));
          if(wNUE==0.0) continue;
          double l_phi = track->Phi();
          wNUE = 1.0/wNUE;
          FillWPCounter(wp,wNUE,l_pt);
      }
    }
    if(wp[0][0]==0) return;
    fV0MMulti->Fill(l_cent);
    double l_rnd = fRndm->Rndm();
    pfmpt->FillProfile(l_cent,wp[0][1]/wp[0][0],wp[0][0],l_rnd);
    getMomentumCorrelation(wp);
    FillCorrelationProfiles(l_cent,l_rnd);
    PostData(1,fCorrList);
}
void AliAnalysisTaskPtCorr::FillCorrelationProfiles(const double &l_cent, double &rn)
{
  for(int m=1;m<=mpar;++m)
  { 
    if(sumw[m]==0.0) continue;
    if(isnan(corr[m]/sumw[m]) || isinf(corr[m]/sumw[m])) continue;
    fptcorr[m-1]->FillProfile(l_cent,corr[m]/sumw[m],sumw[m],rn);
  }
  return;
}
template<typename T> void AliAnalysisTaskPtCorr::getMomentumCorrelation(T& wp)
{
  double sumNum = 0;
  double sumDenum = 0;
  std::vector<double> valNum;
  std::vector<double> valDenum;
  corr.clear();
  sumw.clear();
  corr[0] = 1;
  sumw[0] = 1;
  for(int m(1); m<=mpar; ++m)
  {
    for(int k(1);k<=m;++k)
    {
      valNum.push_back(fSign[k-1]*corr[m-k]*(fFactorial[m-1]/fFactorial[m-k])*wp[k-1][1]);
      valDenum.push_back(fSign[k-1]*sumw[m-k]*(fFactorial[m-1]/fFactorial[m-k])*wp[k-1][0]);
    }
    sumNum = OrderedAddition(valNum, m);
    sumDenum = OrderedAddition(valDenum, m);

    valNum.clear();
    valDenum.clear();
    
    corr[m] = sumNum;
    sumw[m] = sumDenum;
  }
  return;
}
double AliAnalysisTaskPtCorr::OrderedAddition(std::vector<double> vec, int size)
{
  double sum = 0;
  std::sort(vec.begin(), vec.end());

  for(int i = 0; i < size; i++)
  {
    sum += vec[i];
  }
  return sum;
}
int AliAnalysisTaskPtCorr::GetAnalysisStage(TString instr) 
{
  if(instr.Contains("ptcorr")) return 1;
  if(instr.Contains("meanpt")) return 2;
  return 0;
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
  printf("Multi axis set\n");
}
void AliAnalysisTaskPtCorr::SetV0MBins(int nMultiBins, double *multibins) {
  if(fV0MAxis) delete fV0MAxis;
  fV0MAxis = new TAxis(nMultiBins, multibins);
}
