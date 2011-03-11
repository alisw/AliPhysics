// $Id$

#include "AliAnalysisTaskEMCALPi0PbPb.h"
#include <TAxis.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliCentrality.h"
#include "AliEMCALGeoUtils.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliLog.h"

ClassImp(AliAnalysisTaskEMCALPi0PbPb)

//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::AliAnalysisTaskEMCALPi0PbPb() 
  : AliAnalysisTaskSE(),
    fCentVar(),
    fCentFrom(0),
    fCentTo(100),
    fVtxZMin(-7),
    fVtxZMax(+7),
    fUseQualFlag(1),
    fClusName(),
    fOutput(0),
    fEsdEv(0),
    fAodEv(0),
    fRecPoints(0),
    fEsdClusters(0),
    fEsdCells(0),
    fAodClusters(0),
    fAodCells(0),
    fPtRanges(0),
    fHCuts(0x0),
    fHVertexZ(0x0),
    fHCent(0x0),
    fHCentQual(0x0),
    fHClustEtaPhi(0x0),
    fHClustEnergyPt(0x0),
    fHClustEnergySigma(0x0),
    fHClustNTowEnergyRatio(0x0),
    fHPionEtaPhi(0x0),
    fHPionMggPt(0x0),
    fHPionMggAsym(0x0)
{
  // ROOT constructor.
}

//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::AliAnalysisTaskEMCALPi0PbPb(const char *name) 
  : AliAnalysisTaskSE(name),
    fCentVar("V0M"),
    fCentFrom(0),
    fCentTo(100),
    fVtxZMin(-7),
    fVtxZMax(+7),
    fUseQualFlag(1),
    fClusName(),
    fOutput(0),
    fEsdEv(0),
    fAodEv(0),
    fRecPoints(0),
    fEsdClusters(0),
    fEsdCells(0),
    fAodClusters(0),
    fAodCells(0),
    fPtRanges(0),
    fHCuts(0x0),
    fHVertexZ(0x0),
    fHCent(0x0),
    fHCentQual(0x0),
    fHClustEtaPhi(0x0),
    fHClustEnergyPt(0x0),
    fHClustEnergySigma(0x0),
    fHClustNTowEnergyRatio(0x0),
    fHPionEtaPhi(0x0),
    fHPionMggPt(0x0),
    fHPionMggAsym(0x0)
{
  // Constructor.

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.,SPDVertex.,TPCVertex.,EMCALCells. "
               "AOD:header,vertices,emcalCells";
}

//________________________________________________________________________
AliAnalysisTaskEMCALPi0PbPb::~AliAnalysisTaskEMCALPi0PbPb()
{
  // Destructor.

  delete fOutput; fOutput = 0;
  delete fPtRanges; fPtRanges = 0;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::UserCreateOutputObjects()
{
  // Create user objects here.

  fOutput = new TList();
  fOutput->SetOwner();

  fHCuts = new TH1F("hEventCuts","",4,0.5,4.5);
  fHCuts->GetXaxis()->SetBinLabel(1,"All (PS)");
  fHCuts->GetXaxis()->SetBinLabel(2,Form("%s: %.0f-%.0f",fCentVar.Data(),fCentFrom,fCentTo));
  fHCuts->GetXaxis()->SetBinLabel(3,"QFlag");
  fHCuts->GetXaxis()->SetBinLabel(4,Form("zvtx: %.0f-%.0f",fVtxZMin,fVtxZMax));
  fOutput->Add(fHCuts);
  fHVertexZ = new TH1F("hVertexZBeforeCut","",100,-25,25);
  fHVertexZ->SetXTitle("z [cm]");
  fOutput->Add(fHVertexZ);
  fHCent = new TH1F("hCentBeforeCut","",101,-1,100);
  fHCent->SetXTitle(fCentVar.Data());
  fHCentQual = new TH1F("hCentAfterCut","",101,-1,100);
  fHCentQual->SetXTitle(fCentVar.Data());
  fOutput->Add(fHCentQual);

  fHClustEtaPhi = new TH2F("hClustEtaPhi","",100,-0.8,0.8,100,1.2,2.2);
  fHClustEtaPhi->SetXTitle("#eta");
  fHClustEtaPhi->SetYTitle("#varphi");
  fOutput->Add(fHClustEtaPhi);
  fHClustEnergyPt = new TH2F("hClustEnergyPt","",250,0,50,250,0,50);
  fHClustEnergyPt->SetXTitle("E [GeV]");
  fHClustEnergyPt->SetYTitle("p_{T}^{} [GeV/c]");
  fOutput->Add(fHClustEnergyPt);
  fHClustEnergySigma = new TH2F("hClustEnergySigma","",100,0,100,100,0,100);
  fHClustEnergySigma->SetXTitle("E*#sigma_{max}^{} [GeV*cm]");
  fHClustEnergySigma->SetYTitle("#sigma_{max}^{} [cm]");
  //fOutput->Add(fHClustEnergySigma);
  fHClustNTowEnergyRatio = new TH2F("hClustNTowEnergyRatio","",15,-0.5,14.5,101,-0.05,1.05);
  fHClustNTowEnergyRatio->SetXTitle("N towers");
  fHClustNTowEnergyRatio->SetYTitle("Energy ratio");
  //fOutput->Add(fHClustNTowEnergyRatio);

  fHPionEtaPhi = new TH2F("hPionEtaPhi","",100,-0.8,0.8,100,1.2,2.2);
  fHPionEtaPhi->SetXTitle("#eta^{#gamma#gamma}");
  fHPionEtaPhi->SetYTitle("#varphi^{#gamma#gamma}");
  fOutput->Add(fHPionEtaPhi);
  fHPionMggPt = new TH2F("hPionMggPt","",100,0,2,100,0,20.0);
  fHPionMggPt->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
  fHPionMggPt->SetYTitle("p_{T}^{#gamma#gamma} [GeV/c]");
  fOutput->Add(fHPionMggPt);
  fHPionMggAsym = new TH2F("hPionMggAsym","",100,0,2,100,0,1);
  fHPionMggAsym->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
  fHPionMggAsym->SetYTitle("Z_{#gamma#gamma} [GeV]");
  fOutput->Add(fHPionMggAsym);
  const Int_t nbins = 19; 
  Double_t xbins[nbins] = {0.5,1,1.5,2,2.5,3,3.5,4,4.5,6,7,8,9,10,12.5,15,20,25,50};
  fPtRanges = new TAxis(nbins-1,xbins);
  for (Int_t i = 0; i<=nbins; ++i) {
    fHPionInvMasses[i] = new TH1F(Form("HPionInvMass%d",i),"",100,0,2);
    fHPionInvMasses[i]->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    if (i==0)
      fHPionInvMasses[i]->SetTitle(Form("0 < p_{T}^{#gamma#gamma} <%.1f",xbins[0]));
    else if (i==nbins)
      fHPionInvMasses[i]->SetTitle(Form("p_{T}^{#gamma#gamma} > 50"));
    else 
      fHPionInvMasses[i]->SetTitle(Form("%.1f < p_{T}^{#gamma#gamma} <%.1f",xbins[i-1],xbins[i]));
    //cout << i << ": " << fHPionInvMasses[i]->GetTitle() << endl;
    fOutput->Add(fHPionInvMasses[i]);
  }

  PostData(1, fOutput); 
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::UserExec(Option_t *) 
{
  // Called for each event.

  if (!InputEvent())
    return;

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  fEsdEv = dynamic_cast<AliESDEvent*>(InputEvent());
  if (fEsdEv) {
    am->LoadBranch("AliESDRun.");
    am->LoadBranch("AliESDHeader.");
  } else {
    fAodEv = dynamic_cast<AliAODEvent*>(InputEvent());
    am->LoadBranch("header");
  }

  Int_t cut = 1;
  fHCuts->Fill(cut++);

  const AliCentrality *centP = InputEvent()->GetCentrality();
  Double_t cent = centP->GetCentralityPercentileUnchecked(fCentVar);
  fHCent->Fill(cent);
  if (cent<fCentFrom||cent>fCentTo)
    return;

  fHCuts->Fill(cut++);
  
  if (fUseQualFlag) {
    if (centP->GetQuality()>0)
      return;
  }

  fHCentQual->Fill(cent);
  fHCuts->Fill(cut++);

  if (fEsdEv) {
    am->LoadBranch("PrimaryVertex.");
    am->LoadBranch("SPDVertex.");
    am->LoadBranch("TPCVertex.");
  } else {
    fAodEv = dynamic_cast<AliAODEvent*>(InputEvent());
    am->LoadBranch("vertices");
  }

  const AliVVertex *vertex = InputEvent()->GetPrimaryVertex();
  if (!vertex)
    return;

  fHVertexZ->Fill(vertex->GetZ());

  if(vertex->GetZ()<fVtxZMin||vertex->GetZ()>fVtxZMax)
    return;

  fHCuts->Fill(cut++);

  fRecPoints   = 0; // will be set if fClusName is given and AliAnalysisTaskEMCALClusterizeFast is used
  fEsdClusters = 0; // will be set if ESD input used and if fRecPoints are not set of if clusters are attached
  fEsdCells    = 0; // will be set if ESD input used
  fAodClusters = 0; // will be set if AOD input used and if fRecPoints are not set of if clusters are attached
                    //             or if fClusName is given and AliAnalysisTaskEMCALClusterizeFast in AOD output mode
  fAodCells    = 0; // will be set if AOD input used

  // deal with special output from AliAnalysisTaskEMCALClusterizeFast first
  Bool_t clusattached = 0;
  Bool_t recalibrated = 0;
  if (1 && !fClusName.IsNull()) {
    AliAnalysisTaskEMCALClusterizeFast *cltask = 0;
    TObjArray *ts = am->GetTasks();
    cltask = dynamic_cast<AliAnalysisTaskEMCALClusterizeFast*>(ts->FindObject(fClusName));
    if (cltask && cltask->GetClusters()) {
      fRecPoints = cltask->GetClusters();
      clusattached = cltask->GetAttachClusters();
      if (cltask->GetCalibData()!=0)
        recalibrated = kTRUE;
    }
  }
  if (1 && AODEvent() && !fClusName.IsNull()) {
    TList *l = AODEvent()->GetList();
    TClonesArray *clus = 0;
    if (l) {
      clus = dynamic_cast<TClonesArray*>(l->FindObject(fClusName));
      fAodClusters = clus;
    }
  }

  if (fEsdEv) { // ESD input mode
    if (1 && (!fRecPoints||clusattached)) {
      if (!clusattached)
        am->LoadBranch("CaloClusters");
      TList *l = fEsdEv->GetList();
      TClonesArray *clus = 0;
      if (l) {
        clus = dynamic_cast<TClonesArray*>(l->FindObject("CaloClusters"));
        fEsdClusters = clus;
      }
    }
    if (1) {
      if (!recalibrated)
        am->LoadBranch("EMCALCells.");
      fEsdCells = fEsdEv->GetEMCALCells();
    }
  } else if (fAodEv) { // AOD input mode
    if (1 && (!fAodClusters || clusattached)) {
      if (!clusattached)
        am->LoadBranch("caloClusters");
      TList *l = fAodEv->GetList();
      TClonesArray *clus = 0;
      if (l) {
        clus = dynamic_cast<TClonesArray*>(l->FindObject("caloClusters"));
        fAodClusters = clus;
      }
    }
    if (1) {
      if (!recalibrated)
        am->LoadBranch("emcalCells");
      fAodCells = fAodEv->GetEMCALCells();
    }
  } else {
    AliFatal("Impossible to not have either pointer to ESD or AOD event");
  }

  if (1) {
    AliDebug(2,Form("fRecPoints   set: %p", fRecPoints));
    AliDebug(2,Form("fEsdClusters set: %p", fEsdClusters));
    AliDebug(2,Form("fEsdCells    set: %p", fEsdCells));
    AliDebug(2,Form("fAodClusters set: %p", fAodClusters));
    AliDebug(2,Form("fAodCells    set: %p", fAodCells));
  }

  FillCellHists();
  FillClusHists();
  FillPionHists();

  PostData(1, fOutput);
}      

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::Terminate(Option_t *) 
{
  // Terminate called at the end of analysis.
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::FillCellHists()
{
  // Fill histograms related to cell properties.
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::FillClusHists()
{
  // Fill histograms related to clusters.

  Double_t vertex[3] = {0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  // Fill histograms related to cluster properties.
  TLorentzVector clusterVec;
  if (fEsdClusters) {
    Int_t nclus = fEsdClusters->GetEntries();
    for(Int_t i = 0; i<nclus; ++i) {
      AliESDCaloCluster *clus = static_cast<AliESDCaloCluster*>(fEsdClusters->At(i));
      if (!clus)
        continue;
      if (!clus->IsEMCAL()) 
        continue;
      clus->GetMomentum(clusterVec,vertex);
      
      fHClustEtaPhi->Fill(clusterVec.Eta(),clusterVec.Phi());
      fHClustEnergyPt->Fill(clusterVec.E(),clusterVec.Pt());
      //fHClustEnergySigma->Fill();
      //fHClustNTowEnergyRatio->Fill();
    }
  } else if (fAodClusters) {
    Int_t nclus = fAodClusters->GetEntries();
    for (Int_t i = 0; i<nclus; ++i) {
      AliAODCaloCluster *clus = static_cast<AliAODCaloCluster*>(fAodClusters->At(i));
      if (!clus)
        continue;
      if (!clus->IsEMCAL()) 
        continue;
      clus->GetMomentum(clusterVec,vertex);
      
      fHClustEtaPhi->Fill(clusterVec.Eta(),clusterVec.Phi());
      fHClustEnergyPt->Fill(clusterVec.E(),clusterVec.Pt());
      //fHClustEnergySigma->Fill();
      //fHClustNTowEnergyRatio->Fill();
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPi0PbPb::FillPionHists()
{
  // Fill histograms related to pions.

  Double_t vertex[3] = {0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  TLorentzVector clusterVec1;
  TLorentzVector clusterVec2;
  TLorentzVector pionVec;

  if (fEsdClusters) {
    Int_t nclus = fEsdClusters->GetEntries();
    for (Int_t i = 0; i<nclus; ++i) {
      AliESDCaloCluster *clus1 = static_cast<AliESDCaloCluster*>(fEsdClusters->At(i));
      if (!clus1)
        continue;
      if (!clus1->IsEMCAL()) 
        continue;
      clus1->GetMomentum(clusterVec1,vertex);
      for (Int_t j = i+1; j<nclus; ++j) {
        AliESDCaloCluster *clus2 = static_cast<AliESDCaloCluster*>(fEsdClusters->At(j));
        if (!clus2)
          continue;
        if (!clus2->IsEMCAL()) 
          continue;
        clus2->GetMomentum(clusterVec2,vertex);
        pionVec = clusterVec1 + clusterVec2;
        Double_t pionZgg = TMath::Abs(clusterVec1.E()-clusterVec2.E())/pionVec.E();
        fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi()); 
        fHPionMggPt->Fill(pionVec.M(),pionVec.Pt()); 
        fHPionMggAsym->Fill(pionVec.M(),pionZgg); 
        Int_t bin = fPtRanges->FindBin(pionVec.Pt());
        fHPionInvMasses[bin]->Fill(pionVec.M());
      }
    }
  } else if (fAodClusters) {
    Int_t nclus = fAodClusters->GetEntries();
    for (Int_t i = 0; i<nclus; ++i) {
      AliAODCaloCluster *clus1 = static_cast<AliAODCaloCluster*>(fAodClusters->At(i));
      if (!clus1)
        continue;
      if (!clus1->IsEMCAL()) 
        continue;
      clus1->GetMomentum(clusterVec1,vertex);
      for (Int_t j = i+1; j<nclus; ++j) {
        AliAODCaloCluster *clus2 = static_cast<AliAODCaloCluster*>(fAodClusters->At(j));
        if (!clus2)
          continue;
        if (!clus2->IsEMCAL()) 
          continue;
        clus2->GetMomentum(clusterVec2,vertex);
        pionVec = clusterVec1 + clusterVec2;
        Double_t pionZgg = TMath::Abs(clusterVec1.E()-clusterVec2.E())/pionVec.E();
        fHPionEtaPhi->Fill(pionVec.Eta(),pionVec.Phi()); 
        fHPionMggPt->Fill(pionVec.M(),pionVec.Pt()); 
        fHPionMggAsym->Fill(pionVec.M(),pionZgg); 
        Int_t bin = fPtRanges->FindBin(pionVec.Pt());
        fHPionInvMasses[bin]->Fill(pionVec.M());
      }
    }
  }
}
