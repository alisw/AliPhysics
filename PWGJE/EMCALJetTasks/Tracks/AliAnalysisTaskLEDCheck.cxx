#include <THashList.h>
#include "THistManager.h"
#include <TLorentzVector.h>
#include <TMath.h>
#include <TVector.h>

#include "AliEMCALGeometry.h"
#include "AliAnalysisTaskLEDCheck.h"
#include "AliVCaloCells.h"
#include "AliVCluster.h"

ClassImp(AliAnalysisTaskLEDCheck)

AliAnalysisTaskLEDCheck::AliAnalysisTaskLEDCheck():
  AliAnalysisTaskSE(),
  fHistos(NULL),
  fGeometry(NULL)
{
}

AliAnalysisTaskLEDCheck::AliAnalysisTaskLEDCheck(const char *name):
  AliAnalysisTaskSE(name),
  fHistos(NULL),
  fGeometry(NULL)
{
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskLEDCheck::~AliAnalysisTaskLEDCheck() {
}

void AliAnalysisTaskLEDCheck::UserCreateOutputObjects(){
  fHistos = new THistManager("LEDCheck");

  fHistos->CreateTH1("fNumberOfFiringCells", "Number of non-zero cells per event", 10000, 0., 10000.);
  fHistos->CreateTH1("fNumberOfFiringClusters", " Number of clusters per event", 1000, 0., 1000.);
  fHistos->CreateTH2("fNumberOfFiringCellsSM", "Number of non-zero cells per event and supermodule", 20, -0.5, 19.5, 10000, 0., 10000.);
  fHistos->CreateTH2("fNumberOfFiringClustersSM", " Number of clusters per event", 20, -0.5, 19.5, 1000, 0., 1000.);
  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskLEDCheck::UserExec(Option_t *){
  if(!fGeometry){
    fGeometry = AliEMCALGeometry::GetInstance();
    if(!fGeometry)
      fGeometry = AliEMCALGeometry::GetInstanceFromRunNumber(InputEvent()->GetRunNumber());
  }
  // Look only at INT7-B, EMC7-B and DMC7-B events
  TString triggers(InputEvent()->GetFiredTriggerClasses());
  if(!(triggers.Contains("INT7-B") || triggers.Contains("EMC7-B") || triggers.Contains("DMC7-B"))) return;

  const int kNsupermodule = 20;
  Int_t ncell = 0, nclusters = 0;
  Int_t ncellsupermodule[kNsupermodule], nclustersupermodule[kNsupermodule];
  memset(ncellsupermodule, 0, sizeof(Int_t) * kNsupermodule);
  memset(nclustersupermodule, 0, sizeof(Int_t) * kNsupermodule);
  AliVCaloCells *emcalcells = InputEvent()->GetEMCALCells();
  for(int icell = 0; icell < emcalcells->GetNumberOfCells(); icell++){
    if(emcalcells->GetCellAmplitude(icell) > 0){
      ncell++;
      Int_t supermoduleID = 0;
      double globalpos[3];
      fGeometry->GetGlobal(icell, globalpos);
      TVector3 globalvec(globalpos[0], globalpos[1], globalpos[2]);
      Double_t phi = globalvec.Phi();
      if(phi < 0) phi += 2 * TMath::Pi();
      fGeometry->SuperModuleNumberFromEtaPhi(globalvec.Eta(), phi, supermoduleID);
      if(supermoduleID > -1)
        ncellsupermodule[supermoduleID]++;
    }
  }
  if(ncell)
    fHistos->FillTH1("fNumberOfFiringCells", ncell);
  AliVCluster *clust = NULL;
  for(Int_t icluster = 0; icluster < InputEvent()->GetNumberOfCaloClusters(); icluster++){
    clust = InputEvent()->GetCaloCluster(icluster);
    if(clust->IsEMCAL()){
      nclusters++;
      TLorentzVector position;
      Double_t vertex[3];
      InputEvent()->GetPrimaryVertexSPD()->GetXYZ(vertex);
      clust->GetMomentum(position, vertex);
      Int_t supermoduleID(-1);
      Double_t phi = position.Phi();
      if(phi < 0) phi += 2 * TMath::Pi();
      fGeometry->SuperModuleNumberFromEtaPhi(position.Eta(), phi, supermoduleID);
      if(supermoduleID > -1)
        nclustersupermodule[supermoduleID]++;
    }
  }
  if(nclusters)
    fHistos->FillTH1("fNumberOfFiringClusters", nclusters);

  for(int ism = 0; ism < 20; ism++){
    if(ncellsupermodule[ism]) fHistos->FillTH2("fNumberOfFiringCellsSM", ism, ncellsupermodule[ism]);
    if(nclustersupermodule[ism]) fHistos->FillTH2("fNumberOfFiringClustersSM", ism, nclustersupermodule[ism]);
  }
  PostData(1, fHistos->GetListOfHistograms());
}
