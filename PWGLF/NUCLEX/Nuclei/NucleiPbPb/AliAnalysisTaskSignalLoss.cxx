/// \author Luca Barioglio <luca.barioglio@cern.ch>, University and INFN Torino
/// \date Mar 30, 2018

#include "AliAnalysisTaskSignalLoss.h"

#include "AliVEvent.h"
#include "AliAODMCParticle.h"

#include <TClonesArray.h>
#include <TChain.h>

ClassImp(AliAnalysisTaskSignalLoss);

AliAnalysisTaskSignalLoss::AliAnalysisTaskSignalLoss(const char* task_name):
  AliAnalysisTaskSE(task_name),
  fOutputList{nullptr},
  fHistAccEvents{nullptr},
  fHistTrueINELgt0Events{nullptr},
  fHistGenPartsAcc{nullptr},
  fHistGenPartsINELgt0{nullptr},
  fPtBins{0},
  fCentBins{0},
  fPDGcodes{0},
  fNspecies{1},
  fRequireYmin{-0.5},
  fRequireYmax{0.5}
{
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

AliAnalysisTaskSignalLoss::~AliAnalysisTaskSignalLoss(){
  if(fOutputList) delete fOutputList;
}

void AliAnalysisTaskSignalLoss::UserCreateOutputObjects(){
  fOutputList = new TList();
  fOutputList->SetOwner(true);

  const int nPtBins = fPtBins.GetSize() - 1;
  const int nCentBins = fCentBins.GetSize() - 1;

  const float *ptBins = fPtBins.GetArray();
  const float *centBins = fCentBins.GetArray();

  TArrayF arrSpeciesBins;
  arrSpeciesBins.Set(fNspecies+1);
  for(int iSp=0; iSp<=fNspecies; iSp++){
    arrSpeciesBins[iSp]=iSp-0.5;
  }
  const float *speciesBins = arrSpeciesBins.GetArray();
  const int nSpeciesBins = arrSpeciesBins.GetSize()-1;

  fHistAccEvents = new TH1F("fHistAccEvents",";Centrality (%); Accepted events", nCentBins, centBins);
  fHistTrueINELgt0Events = new TH1F("fHistTrueINELgt0Events",";Centrality; True INEL>0 events",nCentBins,centBins);
  fOutputList->Add(fHistAccEvents);
  fOutputList->Add(fHistTrueINELgt0Events);

  char matter_letter[2] = {'M','A'};

  for(int iS=0; iS<2; iS++){
    fHistGenPartsAcc[iS] = new TH3F(Form("fHistGenPartsAcc%c",matter_letter[iS]),";Centrality (%); #it{p}_{T} (GeV/c); counts",nCentBins,centBins,nPtBins,ptBins,nSpeciesBins,speciesBins);
    fHistGenPartsINELgt0[iS] = new TH3F(Form("fHistGenPartsINELgt0%c",matter_letter[iS]),";Centrality (%); #it{p}_{T} (GeV/c); species",nCentBins,centBins,nPtBins,ptBins,nSpeciesBins,speciesBins);
    fOutputList->Add(fHistGenPartsAcc[iS]);
    fOutputList->Add(fHistGenPartsINELgt0[iS]);
  }
  PostData(1,fOutputList);
}

void AliAnalysisTaskSignalLoss::UserExec(Option_t*){

  AliVEvent *ev = InputEvent();

  bool bIsEventTrueINELgt0 = fEventCuts.IsTrueINELgtZero(ev,true);
  bool bIsEventAccepted = fEventCuts.AcceptEvent(ev);
  float centrality = fEventCuts.GetCentrality();

  if(bIsEventTrueINELgt0) fHistTrueINELgt0Events->Fill(centrality);
  if(bIsEventAccepted) fHistAccEvents->Fill(centrality);

  TClonesArray* stack = nullptr;
  stack = (TClonesArray*)ev->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(!stack)
    ::Fatal("AliAnalysisTaskSignalLoss::UserExec","No MC particle array provided");

  for(int iMC=0; iMC < stack->GetEntriesFast(); iMC++){
    AliAODMCParticle* part = (AliAODMCParticle*)stack->UncheckedAt(iMC);
    const int pdgCode = std::abs(part->GetPdgCode());
    const int iCharge = part->Charge() > 0 ? 0 : 1;

    int iSpecies = -1;
    for(int iPDG=0; iPDG<fNspecies; iPDG++){
      if(pdgCode == fPDGcodes[iPDG]) iSpecies = iPDG;
    }
    if(iSpecies<0) continue;
    if (part->Y() > fRequireYmax || part->Y() < fRequireYmin) continue;
    if (!part->IsPhysicalPrimary()) continue;
    if(bIsEventTrueINELgt0) fHistGenPartsINELgt0[iCharge]->Fill(centrality,part->Pt(),iSpecies);
    if(bIsEventAccepted) fHistGenPartsAcc[iCharge]->Fill(centrality,part->Pt(),iSpecies);
  }

  PostData(1,fOutputList);

}

void AliAnalysisTaskSignalLoss::Terminate(Option_t*){
  return;
}

void AliAnalysisTaskSignalLoss::SetPtBins(int nbins, float *bins) {
  fPtBins.Set(nbins + 1, bins);
}

void AliAnalysisTaskSignalLoss::SetCentBins(int nbins, float *bins) {
  fCentBins.Set(nbins + 1, bins);
}

void AliAnalysisTaskSignalLoss::SetPDGcodes(int nbins, int *bins){
  fNspecies = nbins;
  fPDGcodes.Set(nbins,bins);
}
