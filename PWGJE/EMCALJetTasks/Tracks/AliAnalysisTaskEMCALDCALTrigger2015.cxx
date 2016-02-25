#include <map>
#include <vector>

#include <TArrayD.h>
#include <TClonesArray.h>
#include <THashList.h>
#include <THistManager.h>
#include <TLorentzVector.h>

#include "AliAnalysisUtils.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVTrack.h"

#include "AliAnalysisTaskEMCALDCALTrigger2015.h"

#if __cplusplus < 201103L
/*
 * Old C++
 */
#ifndef nullptr
#define nullptr NULL
#endif
#endif

ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEMCALDCALTrigger2015)

namespace EMCalTriggerPtAnalysis {

const TString AliAnalysisTaskEMCALDCALTrigger2015::fgkTriggerClasses[11] = {"INT7", "EMC7", "DMC7", "EJ1", "EJ2", "EG1", "EG2", "DJ1", "DJ2", "DG1", "DG2"};
const TString AliAnalysisTaskEMCALDCALTrigger2015::fgkBeamDirs[4] = {"A", "B", "C", "E"};

AliAnalysisTaskEMCALDCALTrigger2015::AliAnalysisTaskEMCALDCALTrigger2015() :
  AliAnalysisTaskSE(),
  fHistos(nullptr),
  fGeometry(nullptr),
  fClusterContainer(nullptr),
  fPatchContainer(nullptr)
{
}

AliAnalysisTaskEMCALDCALTrigger2015::AliAnalysisTaskEMCALDCALTrigger2015(const char *name) :
  AliAnalysisTaskSE(name),
  fHistos(nullptr),
  fGeometry(nullptr),
  fClusterContainer(nullptr),
  fPatchContainer(nullptr)
{
  DefineOutput(1, TList::Class());
}

void AliAnalysisTaskEMCALDCALTrigger2015::UserCreateOutputObjects(){
  TArrayD energybinning, etabinning, supermodulebinning;
  CreateEnergyBinning(energybinning);
  CreateLinearBinning(etabinning, 100, -1, 1);
  CreateLinearBinning(supermodulebinning, 20, -0.5, 19.5);

  fHistos = new THistManager("histos");
  for(const TString *trgit = fgkTriggerClasses; trgit < fgkTriggerClasses + sizeof(fgkTriggerClasses)/sizeof(const TString); ++trgit){
    for(const TString *beamit = fgkBeamDirs; beamit < fgkBeamDirs + sizeof(fgkBeamDirs)/sizeof(const TString); ++beamit){
      fHistos->CreateTH1(Form("hEventCount%s%s", trgit->Data(), beamit->Data()), Form("Event counter for trigger class %s-%s", trgit->Data(), beamit->Data()), 1, 0.5, 1.5);
      fHistos->CreateTH1(Form("hVertexDistBefore%s%s", trgit->Data(), beamit->Data()), Form("Vertex distribution in trigger class %s-%s before event selection", trgit->Data(), beamit->Data()), 100, -40., 40.);
      fHistos->CreateTH1(Form("hVertexDistAfter%s%s", trgit->Data(), beamit->Data()), Form("Vertex distribution in trigger class %s-%s after event selection", trgit->Data(), beamit->Data()), 100, -40., 40.);
      fHistos->CreateTH1(Form("hClusterEnergyUncalib%s%s", trgit->Data(), beamit->Data()), Form("(Uncalibrated) cluster energy distribution in trigger class %s-%s", trgit->Data(), beamit->Data()), energybinning);
      fHistos->CreateTH1(Form("hClusterEnergyCalib%s%s", trgit->Data(), beamit->Data()), Form("(Calibrated) cluster energy distribution in trigger class %s-%s", trgit->Data(), beamit->Data()), energybinning);
      fHistos->CreateTH2(Form("hClusterEnergyEtaUncalib%s%s", trgit->Data(), beamit->Data()), Form("(Uncalibrated) cluster energy vs eta distribution in trigger class %s-%s", trgit->Data(), beamit->Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hClusterEnergyEtaCalib%s%s", trgit->Data(), beamit->Data()), Form("(Calibrated) cluster energy vs eta distribution in trigger class %s-%s", trgit->Data(), beamit->Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hClusterEnergySupermoduleUncalib%s%s", trgit->Data(), beamit->Data()), Form("(Uncalibrated) cluster energy vs sm distribution in trigger class %s-%s", trgit->Data(), beamit->Data()), supermodulebinning, energybinning);
      fHistos->CreateTH2(Form("hClusterEnergySupermoduleCalib%s%s", trgit->Data(), beamit->Data()), Form("(Calibrated) cluster energy vs sm distribution in trigger class %s-%s", trgit->Data(), beamit->Data()), supermodulebinning, energybinning);
      /*
      fHistos->CreateTH1(Form("hPatchEnergyOnline%s%s", trgit->Data(), beamit->Data()), Form("(Online) patch energy distribution in trigger class %s-%s", trgit->Data(), beamit->Data()), energybinning);
      fHistos->CreateTH1(Form("hPatchEnergyOffline%s%s", trgit->Data(), beamit->Data()), Form("(Offline) patch energy distribution in trigger class %s-%s", trgit->Data(), beamit->Data()), energybinning);
      fHistos->CreateTH2(Form("hPatchEnergyEtaOnline%s%s", trgit->Data(), beamit->Data()), Form("(Online) patch energy vs eta distribution in trigger class %s-%s", trgit->Data(), beamit->Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hPatchEnergyEtaOffline%s%s", trgit->Data(), beamit->Data()), Form("(Offline) patch energy vs eta distribution in trigger class %s-%s", trgit->Data(), beamit->Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hPatchEnergySupermoduleOnline%s%s", trgit->Data(), beamit->Data()), Form("(Online) patch energy vs sm distribution in trigger class %s-%s", trgit->Data(), beamit->Data()), supermodulebinning, energybinning);
      fHistos->CreateTH2(Form("hPatchEnergySupermoduleOffline%s%s", trgit->Data(), beamit->Data()), Form("(Offline) patch energy vs sm distribution in trigger class %s-%s", trgit->Data(), beamit->Data()), supermodulebinning, energybinning);
      */
      for(int ism = 0; ism < 20; ism++){
        fHistos->CreateTH2(Form("hClusterEnergyEtaSupermodule%dCalib%s%s", ism, trgit->Data(), beamit->Data()), Form("(Calibrated) cluster energy vs eta for sm %d in trigger %s-%s", ism, trgit->Data(), beamit->Data()), etabinning, energybinning);
        fHistos->CreateTH2(Form("hClusterEnergyEtaSupermodule%dUncalib%s%s", ism, trgit->Data(), beamit->Data()), Form("(Uncalibrated) cluster energy vs eta for sm %d in trigger %s-%s", ism, trgit->Data(), beamit->Data()), etabinning, energybinning);
        /*
        fHistos->CreateTH2(Form("hPatchEnergyEtaSupermodule%dOnline%s%s", ism, trgit->Data(), beamit->Data()), Form("(Online) patch energy vs eta for sm %d in trigger %s-%s", ism, trgit->Data(), beamit->Data()), etabinning, energybinning);
        fHistos->CreateTH2(Form("hPatchEnergyEtaSupermodule%dOffline%s%s", ism, trgit->Data(), beamit->Data()), Form("(Offline) patch energy vs eta for sm %d in trigger %s-%s", ism, trgit->Data(), beamit->Data()), etabinning, energybinning);
        */
      }
    }
  }

  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEMCALDCALTrigger2015::UserExec(Option_t * /*option*/){
  if(!fGeometry){
    fGeometry = AliEMCALGeometry::GetInstance();
    if(!fGeometry)
      fGeometry = AliEMCALGeometry::GetInstanceFromRunNumber(InputEvent()->GetRunNumber());
  }

  fClusterContainer = static_cast<TClonesArray *>(fInputEvent->FindListObject(fClusterContainerName.Data()));
  fPatchContainer = static_cast<TClonesArray *>(fInputEvent->FindListObject("EmcalTriggers"));

  std::vector<TString> triggerclassesSelected;

  TString classes = fInputEvent->GetFiredTriggerClasses();

  for(const TString *trgit = fgkTriggerClasses; trgit < fgkTriggerClasses + sizeof(fgkTriggerClasses)/sizeof(TString); ++ trgit){
    for(const TString *beamit = fgkBeamDirs; beamit < fgkBeamDirs + sizeof(fgkBeamDirs)/sizeof(TString); ++beamit){
      if(classes.Contains(Form("%s-%s", trgit->Data(), beamit->Data())))
        triggerclassesSelected.push_back(Form("%s%s", trgit->Data(), beamit->Data()));

    }
  }

  if(!triggerclassesSelected.size()) return;

  const AliVVertex *spdvertex = fInputEvent->GetPrimaryVertexSPD();
  if(!spdvertex) return;
  for(std::vector<TString>::iterator trgit = triggerclassesSelected.begin(); trgit < triggerclassesSelected.end(); ++trgit){
    fHistos->FillTH1(Form("hVertexDistBefore%s", trgit->Data()), spdvertex->GetZ());
  }
  /*
  if(spdvertex->GetNContributors() < 1) return;
  if(InputEvent()->IsPileupFromSPD(5, 0.8, 3., 2., 5.)) return;
  if(TMath::Abs(spdvertex->GetZ()) > 10) return;
  int ntrack(0);
  for(int itrk = 0; itrk < InputEvent()->GetNumberOfTracks(); itrk++){
    AliVTrack *trk = static_cast<AliVTrack *>(InputEvent()->GetTrack(itrk));
    if(!(trk->GetStatus() & AliVTrack::kITSpureSA)) continue;
    if(TMath::Abs(trk->Pt()) > 1.) ntrack++;
  }
  if(ntrack < 3) return;
  */

  for(std::vector<TString>::iterator trgit = triggerclassesSelected.begin(); trgit < triggerclassesSelected.end(); ++trgit){
    fHistos->FillTH1(Form("hEventCount%s", trgit->Data()), 1.);
    fHistos->FillTH1(Form("hVertexDistAfter%s", trgit->Data()), spdvertex->GetZ());
  }

  // Loop over uncalibrated clusters
  AliVCluster *clust = nullptr;
  for(int icls = 0; icls < fInputEvent->GetNumberOfCaloClusters(); icls++){
    clust = fInputEvent->GetCaloCluster(icls);
    if(!clust->IsEMCAL()) continue;
    if(clust->GetIsExotic()) continue;
    for(std::vector<TString>::iterator trgit = triggerclassesSelected.begin(); trgit != triggerclassesSelected.end(); ++trgit){
      ProcessCluster(*trgit, clust, false);
    }
  }

  // Loop over calibrated clusters
  if(fClusterContainer){
    for(TIter clustit = TIter(fClusterContainer).Begin(); clustit != TIter::End(); ++clustit){
      clust = dynamic_cast<AliVCluster *>(*clustit);
      if(!clust) continue;
      if(!clust->IsEMCAL()) continue;
      if(clust->GetIsExotic()) continue;
      for(std::vector<TString>::iterator trgit = triggerclassesSelected.begin(); trgit != triggerclassesSelected.end(); ++trgit){
        ProcessCluster(*trgit, clust, true);
      }
    }
  }

  // Loop over trigger patches
  /*
  if(fPatchContainer){
    for(TIter patchiter = TIter(fPatchContainer).Begin(); patchiter != TIter::End(); ++patchiter){
      AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(*patchiter);
      if(!patch) continue;
      if(patch->IsOfflineSimple() && patch->IsGammaLowSimple()){
        for(std::vector<TString>::iterator trgit = triggerclassesSelected.begin(); trgit != triggerclassesSelected.end(); ++trgit){
          ProcessPatch(*trgit, patch, false);
        }
      } else if(!patch->IsOfflineSimple() && patch->IsLevel0()){
        for(std::vector<TString>::iterator trgit = triggerclassesSelected.begin(); trgit != triggerclassesSelected.end(); ++trgit){
          ProcessPatch(*trgit, patch, true);
        }
      }
    }
  }
  */
  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEMCALDCALTrigger2015::ProcessCluster(const TString &triggerclass, const AliVCluster * const clust, bool isCalib){
  if(!clust->IsEMCAL()) return;
  double energy = clust->GetNonLinCorrEnergy();
  fHistos->FillTH1(Form("hClusterEnergy%s%s", isCalib ? "Calib" : "Uncalib", triggerclass.Data()), energy);
  Double_t vertex[3];
  fInputEvent->GetPrimaryVertexSPD()->GetXYZ(vertex);
  TLorentzVector clustervec;
  clust->GetMomentum(clustervec, vertex);
  fHistos->FillTH2(Form("hClusterEnergyEta%s%s", isCalib ? "Calib" : "Uncalib", triggerclass.Data()), clustervec.Eta(), energy);
  Int_t supermoduleID(-1);
  if(fGeometry->SuperModuleNumberFromEtaPhi(clustervec.Eta(), clustervec.Phi(), supermoduleID)){
    fHistos->FillTH1(Form("hClusterEnergySupermodule%s%s", isCalib ? "Calib" : "Uncalib", triggerclass.Data()), supermoduleID, energy);
    fHistos->FillTH2(Form("hClusterEnergyEtaSupermodule%d%s%s", supermoduleID,  isCalib ? "Calib" : "Uncalib", triggerclass.Data()), clustervec.Eta(), energy);
  }
}


void AliAnalysisTaskEMCALDCALTrigger2015::ProcessPatch(const TString &triggerclass, const AliEMCALTriggerPatchInfo * const patch, bool isOnline){
  fHistos->FillTH1(Form("hPatchEnergy%s%s", isOnline ? "Online" : "Offline", triggerclass.Data()), patch->GetPatchE());
  fHistos->FillTH2(Form("hPatchEnergyEta%s%s", isOnline ? "Online" : "Offline", triggerclass.Data()), patch->GetEtaCM(), patch->GetPatchE());
  Int_t supermoduleID(-1);
  if(fGeometry->SuperModuleNumberFromEtaPhi(patch->GetEtaCM(), patch->GetPhiCM(), supermoduleID)){
    fHistos->FillTH1(Form("hPatchEnergySupermodule%s%s", isOnline ? "Online" : "Offline", triggerclass.Data()), supermoduleID, patch->GetPatchE());
    fHistos->FillTH2(Form("hPatchEnergyEtaSupermodule%d%s%s", supermoduleID,  isOnline ? "Online" : "Offline", triggerclass.Data()), patch->GetEtaCM(), patch->GetPatchE());
  }
}

void AliAnalysisTaskEMCALDCALTrigger2015::CreateEnergyBinning(TArrayD& binning) const {
  std::vector<double> mybinning;
  std::map<double,double> definitions;
  definitions.insert(std::pair<double, double>(1, 0.05));
  definitions.insert(std::pair<double, double>(2, 0.1));
  definitions.insert(std::pair<double, double>(4, 0.2));
  definitions.insert(std::pair<double, double>(7, 0.5));
  definitions.insert(std::pair<double, double>(16, 1));
  definitions.insert(std::pair<double, double>(32, 2));
  definitions.insert(std::pair<double, double>(40, 4));
  definitions.insert(std::pair<double, double>(50, 5));
  definitions.insert(std::pair<double, double>(100, 10));
  definitions.insert(std::pair<double, double>(200, 20));
  double currentval = 0.;
  mybinning.push_back(currentval);
  for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
    double limit = id->first, binwidth = id->second;
    while(currentval < limit){
      currentval += binwidth;
      mybinning.push_back(currentval);
    }
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

void AliAnalysisTaskEMCALDCALTrigger2015::CreateLinearBinning(TArrayD& binning, int nbins, double min, double max) const {
  double binwidth = (max-min)/static_cast<double>(nbins);
  binning.Set(nbins+1);
  binning[0] = min;
  double currentlimit = min + binwidth;
  for(int ibin = 0; ibin < nbins; ibin++){
    binning[ibin+1] = currentlimit;
    currentlimit += binwidth;
  }
}


} /* namespace EMCalTriggerPtAnalysis */
