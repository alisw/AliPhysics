#include <sstream>
#include "AliAnalysisTaskEmcalOutliersGen.h"
#include "AliAnalysisManager.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliMCEvent.h"
#include "AliParticleContainer.h"
#include "AliVParticle.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalOutliersGen)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalOutliersGen::AliAnalysisTaskEmcalOutliersGen():
    AliAnalysisTaskEmcalJet(),
    fHistJetPt(nullptr)
{

}

AliAnalysisTaskEmcalOutliersGen::AliAnalysisTaskEmcalOutliersGen(const char *name):
    AliAnalysisTaskEmcalJet(name, kTRUE),
    fHistJetPt(nullptr)
{
    SetMakeGeneralHistograms(true);
    SetGetPtHardBinFromPath(kFALSE);
    SetIsPythia(true);
}

void AliAnalysisTaskEmcalOutliersGen::UserCreateOutputObjects() {
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

    fHistJetPt = new TH1F("fHistJetPt", "Jet pt", 1000, 0., 1000.);
    fOutput->Add(fHistJetPt);

    PostData(1, fOutput);
}

bool AliAnalysisTaskEmcalOutliersGen::Run() {
    auto particles = GetParticleContainer("mcparticlesSelected");
    //auto particles = GetParticleContainer("mcparticles");
    auto jets = GetJetContainer("partjets");

    double partptmin[21] = {0., 20., 25., 30., 40., 50., 70., 80., 100., 120., 140., 180., 200., 250., 270., 300., 350., 380., 420., 450., 600.};
    //std::cout << "Using pt-hard cut " << partptmin[fPtHardBin] << " for pt-hard bin "  << fPtHardBin << std::endl;

    for(auto j : jets->accepted()){
        fHistJetPt->Fill(j->Pt());
        if(j->Pt() > partptmin[fPtHardBin]) {
            std::cout << "Outlier jet found, pt > " << j->Pt() << " GeV/c, NEF = " << j->NEF() << ", " 
                                                    << j->GetNumberOfTracks() << " / " << j->GetNumberOfClusters() << " constituents" << std::endl;
            for(int i = 0; i < j->GetNumberOfTracks(); i++){
                auto part = j->TrackAt(i, particles->GetArray());
                auto z = j->GetZ(part->Px(), part->Py(), part->Pz());
                auto mother = part->GetMother() > -1 ? MCEvent()->GetTrack(part->GetMother()) : nullptr;
                std::cout << "Particle " << i << ": pt " << part->Pt() << " GeV/c, z " << z << ", pdg " << part->PdgCode();
                if(mother) {
                    std::cout << ", mother pdg " << mother->PdgCode() << ", mother pt " << mother->Pt();
                    auto grandmother =  mother->GetMother() > -1 ? MCEvent()->GetTrack(mother->GetMother()) : nullptr;
                    if(grandmother) {
                        std::cout << ", grandmother pdg " << grandmother->PdgCode() << ",  grandmother pt " << grandmother->Pt();

                    } else {
                        std::cout << ", no grand mother";
                    }

                } else {
                    std::cout << ", no mother";
                }
                std::cout << std::endl;
            }
        }
    }

    return true;
}

AliAnalysisTaskEmcalOutliersGen *AliAnalysisTaskEmcalOutliersGen::AddTaskEmcalOutliersGen(const char *name) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScale", "No analysis manager available");
    return nullptr;
  } 

  auto task = new AliAnalysisTaskEmcalOutliersGen(name);
  mgr->AddTask(task);

  auto partcont = task->AddMCParticleContainer("mcparticlesSelected");
  //auto partcont = task->AddMCParticleContainer("mcparticles");
  partcont->SetMinPt(0.);
  auto contjet = task->AddJetContainer(AliJetContainer::kFullJet, AliJetContainer::antikt_algorithm, AliJetContainer::E_scheme, 0.2,
                                                     AliJetContainer::kTPCfid, partcont, nullptr);
  contjet->SetName("partjets");
  contjet->SetMaxTrackPt(1000.);

  std::stringstream outnamebuilder;
  outnamebuilder << mgr->GetCommonFileName() << ":Outliers";

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("OutlierHists", TList::Class(), AliAnalysisManager::kOutputContainer, outnamebuilder.str().data()));

  return task;
}