#include <TClonesArray.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TNtuple.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TRandom.h>

#include "AliVCluster.h"
//#include "AliAODCaloCluster.h"
//#include "AliESDCaloCluster.h"
#include "AliVTrack.h"
#include "AliAODMCParticle.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
//#include "AliClusterContainer.h"
//#include "AliPicoTrack.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSDKL.h"
#include "AliAnalysisTaskSDKLResponse.h"

//ClassImp(AliAnalysisTaskSDKLResponse)

//________________________________________________________________________
AliAnalysisTaskSDKLResponse::AliAnalysisTaskSDKLResponse() :
  AliAnalysisTaskSDKL("AliAnalysisTaskSDKLResponse"),
  fhResponse(),
  fhPtDeltaPt(),
  fhPtDeltaZg(),
  fhPtDeltaRg(),
  fhPtDeltaMg(),
  fhResponseBackSub(),
  fhPtDeltaPtBackSub(0),
  fhPtDeltaZgBackSub(),
  fhPtDeltaRgBackSub(),
  fhPtDeltaMgBackSub(),
  fhResponseDet(),
  fhPtDeltaPtDet(0),
  fhPtDeltaZgDet(),
  fhPtDeltaRgDet(),
  fhPtDeltaMgDet(),
//  fhRho(0),
//  fhRhoSparse(0),
  fhPtDeltaPtAreaBackSub(0),
  fhPtDeltaPtAreaBackSubSparse(0),
  fTreeDL(0),
  fTreeDLUEBS(0),
  fJetsCont1(0),
  fJetsCont2(0),
  fTracksCont1(0),
  fTracksCont2(0),
  fFractionEventsDumpedToTree(0),
  fRandom(0)
{
  // Default constructor.
  // SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskSDKLResponse::AliAnalysisTaskSDKLResponse(const char *name, Int_t const backgroption, Double_t const fractioneventsfortree) :
  AliAnalysisTaskSDKL(name, backgroption),
  fhResponse(),
  fhPtDeltaPt(0),
  fhPtDeltaZg(),
  fhPtDeltaRg(),
  fhPtDeltaMg(),
  fhResponseBackSub(),
  fhPtDeltaPtBackSub(0),
  fhPtDeltaZgBackSub(),
  fhPtDeltaRgBackSub(),
  fhPtDeltaMgBackSub(),
  fhResponseDet(),
  fhPtDeltaPtDet(0),
  fhPtDeltaZgDet(),
  fhPtDeltaRgDet(),
  fhPtDeltaMgDet(),
//  fhRho(0),
//  fhRhoSparse(0),
  fhPtDeltaPtAreaBackSub(0),
  fhPtDeltaPtAreaBackSubSparse(0),
  fTreeDL(0),
  fTreeDLUEBS(0),
  fJetsCont1(0),
  fJetsCont2(0),
  fTracksCont1(0),
  fTracksCont2(0),
  fFractionEventsDumpedToTree(fractioneventsfortree),
  fRandom(0)
{
  // Standard constructor.
  // SetMakeGeneralHistograms(kTRUE);
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskSDKLResponse::~AliAnalysisTaskSDKLResponse()
{
  if (fTreeDL)      delete fTreeDL;
  if (fTreeDLUEBS)  delete fTreeDLUEBS;
  if (fRandom)      delete fRandom;
  // Destructor.
}

AliAnalysisTaskSDKLResponse* AliAnalysisTaskSDKLResponse::AddTaskSoftDropResponse(
  const char *ntracks,
  const char *njets1,
  const char *njets2,
  const char *nrho,
  Int_t       nCentBins,
  Double_t    jetradius,
  Double_t    jetptcut,
  Double_t    jetareacut,
  const char *type,
  Int_t       backgroption,
  Int_t       leadhadtype,
  Double_t    fractioneventsfortree,
  const char *taskname
)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetSample", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetSample", "This task requires an input event handler");
    return NULL;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString trackName(ntracks);

  TString name(taskname);
  if (strcmp(njets1,"")) {
    name += "_";
    name += njets1;
  }
  if (strcmp(njets2,"")) {
    name += "_";
    name += njets2;
  }
  if (strcmp(nrho,"")) {
    name += "_";
    name += nrho;
  }
  if (strcmp(type,"")) {
    name += "_";
    name += type;
  }

  Printf("name: %s",name.Data());

  AliAnalysisTaskSDKLResponse* jetTask = new AliAnalysisTaskSDKLResponse(name, backgroption, fractioneventsfortree);
  //jetTask->SetCentRange(0.,100.);
  //jetTask->SetNCentBins(nCentBins);
  //AliMultSelection *multSelection =static_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
  //if(multSelection) centrality = multSelection->GetMultiplicityPercentile("V0M");

//  AliParticleContainer* partCont = 0;
//  if (trackName == "mcparticles") {
//    partCont = jetTask->AddMCParticleContainer(ntracks);
//  }
//  else if (trackName == "usedefault" || trackName == "tracks" || trackName == "Tracks") {
//    partCont = jetTask->AddTrackContainer(ntracks);
//  }
//  else if (!trackName.IsNull()) {
//    partCont = new AliParticleContainer(trackName);
//  }

//  AliParticleContainer* partCont_pl = jetTask->AddParticleContainer("mcparticles");
  AliParticleContainer* partCont = jetTask->AddTrackContainer("tracks");
  AliParticleContainer* partCont_emb = jetTask->AddTrackContainer("tracks");

  TString strType(type);
  AliJetContainer *jetCont1 = jetTask->AddJetContainer(njets1,strType,jetradius);
  if(jetCont1) {
    jetCont1->SetRhoName(nrho);
//    jetCont1->ConnectParticleContainer(partCont);
    //jetCont->SetZLeadingCut(0.98,0.98);
    jetCont1->SetPercAreaCut(jetareacut);
    jetCont1->SetJetPtCut(jetptcut);
    jetCont1->SetLeadingHadronType(leadhadtype);
  }

  AliJetContainer *jetCont2 = jetTask->AddJetContainer(njets2,strType,jetradius);
  if(jetCont2) {
    jetCont2->SetRhoName(nrho);
    //jetCont2->ConnectParticleContainer(partCont);
    jetCont2->SetPercAreaCut(jetareacut);
    jetCont2->SetJetPtCut(jetptcut);
    jetCont2->SetLeadingHadronType(leadhadtype);
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(jetTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname1(name);
  contname1 += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname1.Data(),
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));

  TString contname2(name);
  contname2 += "_treeDL";
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contname2.Data(),
                                                            TTree::Class(),AliAnalysisManager::kOutputContainer,
                                                            Form("%s", AliAnalysisManager::GetCommonFileName()));

  TString contname3(name);
  contname3 += "_treeDLUEBS";
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(contname3.Data(),
                                                            TTree::Class(),AliAnalysisManager::kOutputContainer,
                                                            Form("%s", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  mgr->ConnectOutput (jetTask, 1, coutput1 );
  mgr->ConnectOutput (jetTask, 2, coutput2 );
  mgr->ConnectOutput (jetTask, 3, coutput3 );

  return jetTask;
}

//________________________________________________________________________
void AliAnalysisTaskSDKLResponse::UserCreateOutputObjects() {
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  OpenFile(1);

  const Int_t nbins = 11;

  //dist eshare iscl pt1   zg1  rg1  mg1/pt1 nsd1                pt2  zg2  rg2 mg2/pt2
  //                        0    1  2    3     4    5    6        7    8    9   10
  Int_t bins[nbins]    = { 10,  10, 2,  200,  20,  20,  20,     200,  20,  20,  20};
  Double_t xmin[nbins] = {0.0, 0.0, 0,  0.0, 0.0, 0.0, 0.0,      0., 0.0, 0.0, 0.0};
  Double_t xmax[nbins] = {0.8, 1.0, 2, 200., 0.5, 0.4, 0.3,    200., 0.5, 0.4, 0.3};

  fhResponse[0] = new THnSparseD("hResponse_1", "hResponse_1", nbins, bins, xmin, xmax); //by splits
  fhResponse[1] = new THnSparseD("hResponse_2", "hResponse_2", nbins, bins, xmin, xmax);

  fhResponseBackSub[0] = new THnSparseD("hResponseBackSub_1", "hResponseBackSub_1", nbins, bins, xmin, xmax);
  fhResponseBackSub[1] = new THnSparseD("hResponseBackSub_2", "hResponseBackSub_2", nbins, bins, xmin, xmax);

  fhResponseDet[0] = new THnSparseD("hResponseDet_1", "hResponseDet_1", nbins, bins, xmin, xmax);
  fhResponseDet[1] = new THnSparseD("hResponseDet_2", "hResponseDet_2", nbins, bins, xmin, xmax);

  xmax[5] = 0.3; //rg1
  xmax[9] = 0.3; //rg2
  xmax[6] = 0.15; //mg1/pt
  xmax[10] = 0.15; //mg2/pt
  fhResponse[2] = new THnSparseD("hResponse_3", "hResponse_3", nbins, bins, xmin, xmax);
  fhResponseBackSub[2] = new THnSparseD("hResponseBackSub_3", "hResponseBackSub_3", nbins, bins, xmin, xmax);
  fhResponseDet[2] = new THnSparseD("hResponseDet_3", "hResponseDet_3", nbins, bins, xmin, xmax);

  xmax[6] = 0.1; //mg1/pt
  xmax[10] = 0.1; //mg2/pt
  fhResponse[3] = new THnSparseD("hResponse_4", "hResponse_4", nbins, bins, xmin, xmax);
  fhResponseBackSub[3] = new THnSparseD("hResponseBackSub_4", "hResponseBackSub_4", nbins, bins, xmin, xmax);
  fhResponseDet[3] = new THnSparseD("hResponseDet_4", "hResponseDet_4", nbins, bins, xmin, xmax);

  for (int s = 0; s < 4; s++) {
    fOutput->Add( fhResponse[s] );
    fOutput->Add( fhResponseBackSub[s] );
    fOutput->Add( fhResponseDet[s] );
  }

  fhPtDeltaPt = new TH2F("fhPtDeltaPt","fhPtDeltaPt",15,0.,150.,200,-20.,20.);
  fOutput->Add(fhPtDeltaPt);

  fhPtDeltaPtBackSub = new TH2F("fhPtDeltaPtBackSub","fhPtDeltaPtBackSub",15,0.,150.,200,-20.,20.);
  fOutput->Add(fhPtDeltaPtBackSub);

  fhPtDeltaPtDet = new TH2F("fhPtDeltaPtDet","fhPtDeltaPtDet",15,0.,150.,200,-20.,20.);
  fOutput->Add(fhPtDeltaPtDet);

  auto r_l = -0.05;
  auto r_r =  0.05;
  for (int s = 0; s < 4; s++) {

    TString dzg_name = "fhPtDeltaZg";
    dzg_name += (s+1);
    fhPtDeltaZg[s] = new TH2F(dzg_name,dzg_name,15,0.,150.,200,-0.2,0.2);

    TString drg_name = "fhPtDeltaRg";
    drg_name += (s+1);
    fhPtDeltaRg[s] = new TH2F(drg_name,drg_name,15,0.,150.,200,r_l,r_r);

    TString dmg_name = "fhPtDeltaMg";
    dmg_name += (s+1);
    fhPtDeltaMg[s] = new TH2F(dmg_name,dmg_name,15,0.,150.,400,-0.3,0.3);

    fOutput->Add( fhPtDeltaZg[s] );
    fOutput->Add( fhPtDeltaRg[s] );
    fOutput->Add( fhPtDeltaMg[s] );

    dzg_name += "BackSub";
    fhPtDeltaZgBackSub[s] = new TH2F(dzg_name,dzg_name,15,0.,150.,200,-0.2,0.2);
    drg_name += "BackSub";
    fhPtDeltaRgBackSub[s] = new TH2F(drg_name,drg_name,15,0.,150.,200,r_l,r_r);
    dmg_name += "BackSub";
    fhPtDeltaMgBackSub[s] = new TH2F(dmg_name,dmg_name,15,0.,150.,400,-0.3,0.3);

    fOutput->Add( fhPtDeltaZgBackSub[s] );
    fOutput->Add( fhPtDeltaRgBackSub[s] );
    fOutput->Add( fhPtDeltaMgBackSub[s] );

    dzg_name = "fhPtDeltaZgDet";
    fhPtDeltaZgDet[s] = new TH2F(dzg_name,dzg_name,15,0.,150.,200,-0.2,0.2);
    drg_name = "fhPtDeltaRgDet";
    fhPtDeltaRgDet[s] = new TH2F(drg_name,drg_name,15,0.,150.,200,r_l,r_r);
    dmg_name = "fhPtDeltaMgDet";
    fhPtDeltaMgDet[s] = new TH2F(dmg_name,dmg_name,15,0.,150.,400,-0.3,0.3);

    fOutput->Add( fhPtDeltaZgDet[s] );
    fOutput->Add( fhPtDeltaRgDet[s] );
    fOutput->Add( fhPtDeltaMgDet[s] );

  }

  fhRho = new TH1F("fhRho","fhRho",1000,0,100);
  fOutput->Add(fhRho);

  fhRhoSparse = new TH1F("fhRhoSparse","fhRhoSparse",1000,0,100);
  fOutput->Add(fhRhoSparse);

  fhPtDeltaPtAreaBackSub = new TH2F("fhPtDeltaPtAreaBackSub","fhPtDeltaPtAreaBackSub",20,0.,200.,200,-20.,20.);
  fOutput->Add(fhPtDeltaPtAreaBackSub);

  fhPtDeltaPtAreaBackSubSparse = new TH2F("fhPtDeltaPtAreaBackSubSparse","fhPtDeltaPtAreaBackSubSparse",20,0.,200.,200,-20.,20.);
  fOutput->Add(fhPtDeltaPtAreaBackSubSparse);

  fTreeDL = new TNtuple("JetTrackTreeDL", "jet-track tree dl", "pt:eta:phi:jetm");
  PostData(2, fTreeDL);

  fTreeDLUEBS = new TNtuple("JetTrackTreeDLUEBS", "jet-track tree dl+ue backgr sub", "pt:eta:phi:jetm");
  PostData(3, fTreeDLUEBS);

  fRandom = new TRandom;

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSDKLResponse::FillHistograms() {

  auto const area_nominal = TMath::Pi() * 0.4 * 0.4;
  bool isDumpEventToTree = fRandom->Uniform() < fFractionEventsDumpedToTree;

  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4, fastjet::E_scheme);
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(0.9, 1));
  fastjet::Selector sel_jets = fastjet::SelectorAbsEtaMax(0.9 - 0.4); //max_eta_jet

  // Fill histograms.

  std::vector<mjet> mjet_cont_pl; //probe PL
  std::vector<mjet> mjet_cont_dlue; //response DL+UE
  FillMjetContainer(fJetsCont1, mjet_cont_pl);
  FillMjetContainer(fJetsCont2, mjet_cont_dlue);

//  if (fJetsCont1) {
//    for (auto jet : fJetsCont1->accepted()) {
//    }
//  }

  std::vector<fastjet::PseudoJet> event_dl;
  AddTracksToEvent(fTracksCont2, event_dl); //only pythia det level tracks
  fastjet::ClusterSequence cs_dl(event_dl, jet_def);
  std::vector<fastjet::PseudoJet> jets_dl = sorted_by_pt(sel_jets(cs_dl.inclusive_jets()));

  std::vector<fastjet::PseudoJet> jets_dl_filtered;
  for (auto jet : jets_dl) {
    if ( jet.pt() < 1. ) continue;
    if ( jet.has_area() ) {
      auto jarea = jet.area_4vector().perp();
      if (jarea < (0.6 * area_nominal)) continue;
    }
    jets_dl_filtered.push_back(jet);
  }

  std::vector<mjet> mjet_cont_dl;
  FillMjetContainer(jets_dl_filtered, mjet_cont_dl);

  //PL-DL, no UE
  std::vector<fastjet::PseudoJet> jets_dl_matched;
  for (int i = 0; i < mjet_cont_pl.size(); i++) {
    for (int j = 0; j < mjet_cont_dl.size(); j++) {
      auto mjet1 = mjet_cont_pl[i];
      auto mjet2 = mjet_cont_dl[j];
      auto dist = CalcDist(mjet1, mjet2);
      if (dist > 0.8) continue; //loose cut
      auto eshare = CalcEnergyShare(mjet1, mjet2);
      if (eshare < 0.1) continue; //loose cut

      //delayed calculation of the splits
      //deep declustering is performed
      //only for the (roughly) matched jets
      if (mjet1.pointerAJet) mjet1.splits = ReclusterFindHardSplits( mjet1.pointerAJet );
      if (mjet2.pointerPJet) mjet2.splits = ReclusterFindHardSplits( *(mjet2.pointerPJet) );

      auto iscl = IsClosestPair(i,j,mjet_cont_pl,mjet_cont_dl);
      FillResponseFromMjets(mjet1, mjet2, fhResponseDet, dist, eshare, iscl);
      if ( (dist < 0.4) && (eshare > 0.5) && iscl ) { //strict cuts
        FillDeltasFromMjets(mjet1, mjet2, fhPtDeltaPtDet, fhPtDeltaZgDet, fhPtDeltaRgDet, fhPtDeltaMgDet);
        if (mjet2.pointerPJet) jets_dl_matched.push_back( *(mjet2.pointerPJet) );
      }
    }
  }

  //dump only matched jets
  if (isDumpEventToTree) {
    FillTree(jets_dl_matched, fTreeDL);
    PostData(2, fTreeDL);
  }

  //FULL EVENT
  std::vector<fastjet::PseudoJet> event_full;
  AddTracksToEvent(fTracksCont1, event_full);
  AddTracksToEvent(fTracksCont2, event_full);

  Double_t rho;
  Double_t rho_sparse;
  std::vector<fastjet::PseudoJet> event_backsub = GetBackSubEvent(event_full, rho, rho_sparse, fbcoption);
  fhRho->Fill(rho);
  fhRhoSparse->Fill(rho_sparse);

  fastjet::ClusterSequenceArea clust_seq_backsub(event_backsub, jet_def, area_def);
  std::vector<fastjet::PseudoJet> jets_backsub = sorted_by_pt( sel_jets(clust_seq_backsub.inclusive_jets()) );

  std::vector<fastjet::PseudoJet> jets_backsub_filtered;
  for (auto jet : jets_backsub) {
    if ( jet.pt() < 1. ) continue;
    if ( jet.has_area() ) {
      auto jarea = jet.area_4vector().perp();
      if (jarea < (0.6 * area_nominal)) continue;
    }
    jets_backsub_filtered.push_back(jet);
  }

  std::vector<mjet> mjet_container_backsub;
  FillMjetContainer(jets_backsub_filtered, mjet_container_backsub);

  //background-subtracted jets
  std::vector<fastjet::PseudoJet> jets_backsub_filtered_matched;
  for (int i = 0; i < mjet_cont_pl.size(); i++) {
    for (int j = 0; j < mjet_container_backsub.size(); j++) {
      auto mjet1 = mjet_cont_pl[i];
      auto mjet2 = mjet_container_backsub[j];
      auto dist = CalcDist(mjet1, mjet2);
      if (dist > 0.8) continue;   //loose cut
      auto eshare = CalcEnergyShare(mjet1, mjet2);
      if (eshare < 0.1) continue; //loose cut

      //delayed calculation
      if (mjet1.pointerAJet) mjet1.splits = ReclusterFindHardSplits( mjet1.pointerAJet );
      if (mjet2.pointerPJet) mjet2.splits = ReclusterFindHardSplits( *(mjet2.pointerPJet) );

      auto iscl = IsClosestPair(i,j,mjet_cont_pl,mjet_container_backsub);
      FillResponseFromMjets(mjet1, mjet2, fhResponseBackSub, dist, eshare, iscl);
      if ( (dist < 0.4) && (eshare > 0.5) && iscl ) { //strict cuts
        FillDeltasFromMjets(mjet1, mjet2, fhPtDeltaPtBackSub, fhPtDeltaZgBackSub, fhPtDeltaRgBackSub, fhPtDeltaMgBackSub);
        if (mjet2.pointerPJet) jets_backsub_filtered_matched.push_back( *(mjet2.pointerPJet) );
      }
    }
  }

  //and now dump only matched jets
  if (isDumpEventToTree) {
    FillTree(jets_backsub_filtered_matched, fTreeDLUEBS);
    PostData(3, fTreeDLUEBS);
  }

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskSDKLResponse::ExecOnce() {


  AliAnalysisTaskEmcalJet::ExecOnce();

  fJetsCont1           = GetJetContainer(0);
  fJetsCont2           = GetJetContainer(1);
//  fTracksCont          = GetParticleContainer(0); //old

//  fTracksCont1 = GetParticleContainer(0);
//  fTracksCont1->SetClassName("AliAODMCParticle");

//  fTracksCont0 = GetParticleContainer(0);
//  fTracksCont0->SetClassName("AliAODMCParticle");

  fTracksCont1 = GetParticleContainer(0);
  fTracksCont1->SetClassName("AliVTrack");

  fTracksCont2 = GetParticleContainer(1);
  if (fTracksCont2) fTracksCont2->SetClassName("AliVTrack");

  if (fJetsCont1 && fJetsCont1->GetArray() == 0) fJetsCont1 = 0;
  if (fJetsCont2 && fJetsCont2->GetArray() == 0) fJetsCont2 = 0;
//  if (fTracksCont0 && fTracksCont0->GetArray() == 0) fTracksCont0 = 0;
  if (fTracksCont1 && fTracksCont1->GetArray() == 0) fTracksCont1 = 0;
  if (fTracksCont2 && fTracksCont2->GetArray() == 0) fTracksCont2 = 0;

}

int AliAnalysisTaskSDKLResponse::FillResponseFromMjets(mjet const & mjet1, mjet const & mjet2, THnSparse* hr[4], Float_t const dist, Float_t const eshare, Bool_t const iscl) {

  std::vector<split> const & splits_1 = mjet1.splits;
  std::vector<split> const & splits_2 = mjet2.splits;
  auto pt1 = mjet1.pt;
  auto pt2 = mjet2.pt;

  auto nsd_1 = splits_1.size();
  float zg_1[4] = {0.};
  float rg_1[4] = {0.};
  float mg_1[4] = {-1.};
  for (int i = 0; (i < 4) && (i < nsd_1); i++) {
    zg_1[i] = splits_1[i].z;
    rg_1[i] = splits_1[i].r;
    mg_1[i] = splits_1[i].m;
  }

  auto nsd_2 = splits_2.size();
  float zg_2[4] = {0.};
  float rg_2[4] = {0.};
  float mg_2[4] = {-1.};
  for (int i = 0; (i < 4) && (i < nsd_2); i++) {
    zg_2[i] = splits_2[i].z;
    rg_2[i] = splits_2[i].r;
    mg_2[i] = splits_2[i].m;
  }

  //dist eshare iscl pt1   zg1  rg1  mg1/pt1 nsd1                pt2  zg2  rg2 mg2/pt2
  //                          0    1  2    3     4    5    6        7    8    9   10
//  Int_t bins[nbins]    = { 10,  10, 2,  200,  20,  20,  20,     200,  20,  20,  20};
//  Double_t xmin[nbins] = {0.0, 0.0, 0,  0.0, 0.0, 0.0, 0.0,      0., 0.0, 0.0, 0.0};
//  Double_t xmax[nbins] = {0.8, 1.0, 2, 200., 0.5, 0.4, 0.3,    200., 0.5, 0.4, 0.3};

  Double_t xvalue[4][11]; //4 hard splits
  for (int s = 0; s < 4; s++) {
    xvalue[s][0] = dist;
    xvalue[s][1] = eshare;
    xvalue[s][2] = 0.1 + iscl;

    xvalue[s][3] = pt1;
    xvalue[s][4] = zg_1[s];
    xvalue[s][5] = rg_1[s];
    xvalue[s][6] = mg_1[s] / pt1;

    xvalue[s][7]  = pt2;
    xvalue[s][8]  = zg_2[s];
    xvalue[s][9]  = rg_2[s];
    xvalue[s][10] = mg_2[s] / pt2;

    hr[s]->Fill( xvalue[s] );
  }

  return 0;

}


int AliAnalysisTaskSDKLResponse::FillDeltasFromMjets(mjet const & mjet1, mjet const & mjet2, TH2F *hptdpt, TH2F *h1[4], TH2F *h2[4], TH2F *h3[4]) {

  auto dist = CalcDist(mjet1, mjet2);
  if ( dist > (0.4*0.8) ) return 1; //tight cut - only "matched" jets

  auto eshare = CalcEnergyShare(mjet1, mjet2);
  if (eshare < 0.5) return 1; //eshare cut

  std::vector<split> const & splits_1 = mjet1.splits;
  std::vector<split> const & splits_2 = mjet2.splits;
  auto pt1 = mjet1.pt;
  auto pt2 = mjet2.pt;

  int nsd_1 = splits_1.size();
  float zg_1[4] = {0.};
  float rg_1[4] = {0.};
  float mg_1[4] = {-1.};
  for (int i = 0; (i < 4) && (i < nsd_1); i++) {
    zg_1[i] = splits_1[i].z;
    rg_1[i] = splits_1[i].r;
    mg_1[i] = splits_1[i].m;
  }

  int nsd_2 = splits_2.size();
  float zg_2[4] = {0.};
  float rg_2[4] = {0.};
  float mg_2[4] = {-1.};
  for (int i = 0; (i < 4) && (i < nsd_2); i++) {
    zg_2[i] = splits_2[i].z;
    rg_2[i] = splits_2[i].r;
    mg_2[i] = splits_2[i].m;
  }

  hptdpt->Fill(pt1, pt2 - pt1);
  for (int s = 0; s < 4; s++) {
    if ( (zg_1[s] > 0.05) && (zg_2[s] > 0.05) ) {
      h1[s]->Fill(pt1, zg_2[s] - zg_1[s]);
      h2[s]->Fill(pt1, rg_2[s] - rg_1[s]);
      h3[s]->Fill(pt1, mg_2[s] - mg_1[s]);
    }//only tagged jets
  }//splits iter

  return 0;

}

void AliAnalysisTaskSDKLResponse::FillMjetContainer(AliJetContainer *jet_container, std::vector <mjet> & mjet_container) {
  if (jet_container) {
    for ( auto jet : jet_container->accepted() ) {

      double pt_scalar = 0.0;
      std::vector<double> track_pts;
      std::vector<int> track_labels;

      UShort_t ntracks = jet->GetNumberOfTracks();
      for (int j = 0; j < ntracks; j++) {
        auto jtrack = jet->Track(j);
        auto track_pt = jtrack->Pt();
        pt_scalar += track_pt;
        auto mclabel = jtrack->GetLabel();
        if (mclabel > -1) { //keep only MC labels
          track_pts.push_back( track_pt );
          track_labels.push_back( mclabel );
        }
      }
      std::vector<split> empty;
      mjet_container.push_back( mjet{ jet->Pt(), pt_scalar, jet->Eta(), jet->Phi(), jet->Area(), track_pts, track_labels, empty, jet, nullptr } );
    }
  }
}

void AliAnalysisTaskSDKLResponse::FillMjetContainer(std::vector <fastjet::PseudoJet> & jet_container, std::vector <mjet> & mjet_container) {

  for ( int k = 0; k < jet_container.size(); k++ ) {
    auto jet = jet_container[k];
    if ( jet.perp() < 1.e-6 ) continue; //skip ghost-only jets
    double area = 0.0;
    double pt_scalar = 0.0;
    std::vector<double> track_pts;
    std::vector<int> track_labels;

    for (auto c : jet.constituents() ) {
      pt_scalar += c.pt();
      auto mclabel = c.user_index();
      if (mclabel > -1) { //keep only MC labels
        track_pts.push_back( c.pt() );
        track_labels.push_back( mclabel );
      }
    }

    if ( jet.has_area() ) area = jet.area_4vector().perp();
    std::vector<split> empty;
    mjet_container.push_back( mjet{ jet.perp(), pt_scalar, jet.eta(), jet.phi(), area, track_pts, track_labels, empty, nullptr, &(jet_container[k]) } );
  }

}

//________________________________________________________________________
Bool_t AliAnalysisTaskSDKLResponse::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
void AliAnalysisTaskSDKLResponse::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

Float_t AliAnalysisTaskSDKLResponse::CalcDist(mjet const & mjet1, mjet const & mjet2) {

  auto phi1 = mjet1.phi;
  auto eta1 = mjet1.eta;
  auto phi2 = mjet2.phi;
  auto eta2 = mjet2.eta;

  auto phi_diff = fabs(phi1 - phi2);
  if ( phi_diff > TMath::Pi() ) phi_diff = 2.0 * TMath::Pi() - phi_diff;
  auto eta_diff = fabs(eta1 - eta2);

  auto dist = sqrt(phi_diff*phi_diff + eta_diff*eta_diff);

  return dist;

}

Bool_t AliAnalysisTaskSDKLResponse::IsClosestPair(int const idx1, int const idx2, std::vector <mjet> const & mjet_container1, std::vector <mjet> const & mjet_container2) {

  double min_dist12 = 1000.0;
  int min_dist12_index2 = -7;

  auto size2 = mjet_container2.size();
  for (int k = 0; k < size2; k++) {
    auto cur_dist = CalcDist( mjet_container1[idx1], mjet_container2[k] );
    if (cur_dist < min_dist12) {
      min_dist12 = cur_dist;
      min_dist12_index2 = k;
    }
  }

  if ( idx2 != min_dist12_index2 ) return kFALSE;

  double min_dist21 = 1000.0;
  int min_dist21_index1 = -7;

  auto size1 = mjet_container1.size();
  for (int l = 0; l < size1; l++) {
    auto cur_dist = CalcDist( mjet_container2[idx2], mjet_container1[l] );
    if (cur_dist < min_dist21) {
      min_dist21 = cur_dist;
      min_dist21_index1 = l;
    }
  }

  if ( idx1 != min_dist21_index1 ) return kFALSE;

  return kTRUE;

}

Float_t AliAnalysisTaskSDKLResponse::CalcEnergyShare(mjet const & mjet1, mjet const & mjet2) {

  auto track_pts2 = mjet2.track_pts;
  auto track_labels1 = mjet1.track_labels;
  auto track_labels2 = mjet2.track_labels;

  double pt_tot_from_jet1 = 0.0;
  for (int j = 0; j < track_labels2.size(); j++) {
    if ( std::find( track_labels1.begin(), track_labels1.end(), track_labels2[j] ) != track_labels1.end() ) {
      pt_tot_from_jet1 += track_pts2[j];
    }
  }

  return pt_tot_from_jet1/mjet2.pt_scalar;

}