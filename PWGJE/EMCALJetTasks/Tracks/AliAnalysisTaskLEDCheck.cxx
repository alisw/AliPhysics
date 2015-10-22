#include <THashList.h>
#include "THistManager.h"

#include "AliAnalysisTaskLEDCheck.h"
#include "AliVCaloCells.h"
#include "AliVCluster.h"

ClassImp(AliAnalysisTaskLEDCheck)

AliAnalysisTaskLEDCheck::AliAnalysisTaskLEDCheck():
  AliAnalysisTaskSE(),
  fHistos(NULL)
{
}

AliAnalysisTaskLEDCheck::AliAnalysisTaskLEDCheck(const char *name):
  AliAnalysisTaskSE(name),
  fHistos(NULL)
{
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskLEDCheck::~AliAnalysisTaskLEDCheck() {
}

void AliAnalysisTaskLEDCheck::UserCreateOutputObjects(){
  fHistos = new THistManager("LEDCheck");

  fHistos->CreateTH1("fNumberOfFiringCells", "Number of non-zero cells per event", 10000, 0., 10000.);
  fHistos->CreateTH1("fNumberOfFiringClusters", " Number of clusters per event", 1000, 0., 1000.);
  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskLEDCheck::UserExec(Option_t *){

  // Look only at INT7-B, EMC7-B and DMC7-B events
  TString triggers(InputEvent()->GetFiredTriggerClasses());
  if(!(triggers.Contains("INT7-B") || triggers.Contains("EMC7-B") || triggers.Contains("DMC7-B"))) return;

  Int_t ncell = 0, nclusters = 0;
  AliVCaloCells *emcalcells = InputEvent()->GetEMCALCells();
  for(int icell = 0; icell < emcalcells->GetNumberOfCells(); icell++){
    if(emcalcells->GetCellAmplitude(icell) > 0) ncell++;
  }
  if(ncell)
    fHistos->FillTH1("fNumberOfFiringCells", ncell);
  AliVCluster *clust = NULL;
  for(Int_t icluster = 0; icluster < InputEvent()->GetNumberOfCaloClusters(); icluster++){
    clust = InputEvent()->GetCaloCluster(icluster);
    if(clust->IsEMCAL()) nclusters++;
  }
  if(nclusters)
    fHistos->FillTH1("fNumberOfFiringClusters", nclusters);
  PostData(1, fHistos->GetListOfHistograms());
}
