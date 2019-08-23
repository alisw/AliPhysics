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

  fhDist = new TH1F("fhDist", "fhDist", 100, 0., 2.0);
  fOutput->Add(fhDist);

  const Int_t nbins = 15;

  //dist    pt1   zg1  rg1  mg1  mg1/pt1 nsd1 sdstep1 pt2 zg2 rg2 mg2 mg2/pt2 nsd2 sdstep2
  //                      0    1     2    3    4    5     6    7    8     9    10   11   12    13   14
  Int_t bins[nbins]    = {10,  30,   20,  20,  20,  20,   10,  10,  30,   20,  20,  20,  20,   10,  10};
  Double_t xmin[nbins] = {0.,  0.,   0.,  0.,  0.,  0.00, 0.,  0.,  0.,   0.,  0.,  0.,  0.00, 0.,  0.};
  Double_t xmax[nbins] = {0.8, 150., 0.5, 0.4, 20., 0.30, 10., 10., 150., 0.5, 0.4, 20., 0.30, 10., 10.};

  fhResponse[0] = new THnSparseD("hResponse_1", "hResponse_1", nbins, bins, xmin, xmax); //by splits
  fhResponse[1] = new THnSparseD("hResponse_2", "hResponse_2", nbins, bins, xmin, xmax);

  fhResponseBackSub[0] = new THnSparseD("hResponseBackSub_1", "hResponseBackSub_1", nbins, bins, xmin, xmax);
  fhResponseBackSub[1] = new THnSparseD("hResponseBackSub_2", "hResponseBackSub_2", nbins, bins, xmin, xmax);

  fhResponseDet[0] = new THnSparseD("hResponseDet_1", "hResponseDet_1", nbins, bins, xmin, xmax);
  fhResponseDet[1] = new THnSparseD("hResponseDet_2", "hResponseDet_2", nbins, bins, xmin, xmax);

  xmax[3] = 0.3; //rg1
  xmax[10] = 0.3; //rg2
  xmax[5] = 0.15; //mg1/pt
  xmax[12] = 0.15; //mg2/pt
  fhResponse[2] = new THnSparseD("hResponse_3", "hResponse_3", nbins, bins, xmin, xmax);
  fhResponseBackSub[2] = new THnSparseD("hResponseBackSub_3", "hResponseBackSub_3", nbins, bins, xmin, xmax);
  fhResponseDet[2] = new THnSparseD("hResponseDet_3", "hResponseDet_3", nbins, bins, xmin, xmax);

  xmax[5] = 0.1; //mg1/pt
  xmax[12] = 0.1; //mg2/pt
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

  fhPtDeltaPtAreaBackSub = new TH2F("fhPtDeltaPtAreaBackSub","fhPtDeltaPtAreaBackSub",15,0.,150.,200,-20.,20.);
  fOutput->Add(fhPtDeltaPtAreaBackSub);

  fhPtDeltaPtAreaBackSubSparse = new TH2F("fhPtDeltaPtAreaBackSubSparse","fhPtDeltaPtAreaBackSubSparse",15,0.,150.,200,-20.,20.);
  fOutput->Add(fhPtDeltaPtAreaBackSubSparse);

  fTreeDL = new TNtuple("JetTrackTreeDL", "jet-track tree dl", "pt:eta:phi:jetm");
  PostData(2, fTreeDL);

  fTreeDLUEBS = new TNtuple("JetTrackTreeDLUEBS", "jet-track tree dl+ue backgr sub", "pt:eta:phi:jetm");
  PostData(3, fTreeDLUEBS);

  fRandom = new TRandom;

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSDKLResponse::FillHistograms()
{
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4, fastjet::pt_scheme);
  fastjet::AreaDefinition area_def( fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(0.9,1) );
  fastjet::Selector sel_jets = fastjet::SelectorAbsEtaMax(0.9 - 0.4); //max_eta_jet

  // Fill histograms.

  std::vector <mjet> mjet_cont_1; //probe PL
  std::vector <mjet> mjet_cont_2; //response DL+UE
  FillMjetContainer(fJetsCont1, mjet_cont_1);
  FillMjetContainer(fJetsCont2, mjet_cont_2);

  std::vector <fastjet::PseudoJet> event_dl;
  FillAllTracks(fTracksCont2, nullptr, event_dl); //only pythia det level tracks
  fastjet::ClusterSequence cs_dl(event_dl, jet_def);
  std::vector <fastjet::PseudoJet> jets_dl = sorted_by_pt( sel_jets( cs_dl.inclusive_jets() ) );
  std::vector<mjet> mjet_cont_dl;
  FillMjetContainer(jets_dl, mjet_cont_dl);

  bool isDumpEventToTree = fRandom->Uniform() < fFractionEventsDumpedToTree;

  if (isDumpEventToTree) {
    FillTree(jets_dl, fTreeDL);
    PostData(2, fTreeDL);
  }

  //PL-DL, no UE
  for ( auto mjet1 : mjet_cont_1 ) {
    for ( auto mjet2 : mjet_cont_dl ) {
      FillResponseFromMjets(mjet1, mjet2, fhResponseDet);
      FillDeltasFromMjets(mjet1, mjet2, fhPtDeltaPtDet, fhPtDeltaZgDet, fhPtDeltaRgDet, fhPtDeltaMgDet);
    }
  }

  std::vector <fastjet::PseudoJet> event_full;
  FillAllTracks(fTracksCont1, fTracksCont2, event_full);

  fastjet::ClusterSequenceArea cs_full(event_full, jet_def, area_def);
  std::vector<fastjet::PseudoJet> jets_full = sorted_by_pt( sel_jets( cs_full.inclusive_jets() ) );
//  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4, fastjet::pt_scheme);
//  fastjet::ClusterSequence cs(full_event, jet_def);
//  std::vector <fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

  std::vector<mjet> mjet_cont_full;
  FillMjetContainer(jets_full, mjet_cont_full);

  Double_t rho;
  Double_t rho_sparse;
  std::vector <fastjet::PseudoJet> event_backsub = GetBackSubEvent(event_full, rho, rho_sparse, fbcoption);
  fhRho->Fill(rho);
  fhRhoSparse->Fill(rho_sparse);

  //area subtraction
  //background-subtracted jets
  for ( auto mjet1 : mjet_cont_1 ) {
    for ( auto mjet2 : mjet_cont_full ) {
      auto dist = CalcDist(mjet1, mjet2);
      if ( dist > (0.4*0.8) ) continue;
      auto pt1 = mjet1.pt;
      auto pt2 = mjet2.pt;
      auto area2 = mjet2.area;
      fhPtDeltaPtAreaBackSub->Fill( pt1, (pt2 - rho * area2) - pt1);
      fhPtDeltaPtAreaBackSubSparse->Fill( pt1, (pt2 - rho_sparse * area2) - pt1);
    }
  }

  //0.9 -> max_eta
//  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4, fastjet::pt_scheme);
//  fastjet::AreaDefinition area_def( fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(0.9,1) );
  fastjet::ClusterSequenceArea clust_seq_backsub(event_backsub, jet_def, area_def);
//  fastjet::Selector sel_jets = fastjet::SelectorAbsEtaMax(0.9 - 0.4); //max_eta_jet
  std::vector<fastjet::PseudoJet> jets_backsub = sorted_by_pt( sel_jets( clust_seq_backsub.inclusive_jets() ) );

  if (isDumpEventToTree) {
    FillTree(jets_backsub, fTreeDLUEBS);
    PostData(3, fTreeDLUEBS);
  }

  std::vector<mjet> mjet_container_backsub;

  FillMjetContainer(jets_backsub, mjet_container_backsub);

//  std::cout<<"back sub DL"<<std::endl;
//  for ( auto j : mjet_containerBackSub ) {
//  }

  for ( auto mjet1 : mjet_cont_1 ) {
    for ( auto mjet2 : mjet_cont_2 ) {
      FillResponseFromMjets(mjet1, mjet2, fhResponse);
      FillDeltasFromMjets(mjet1, mjet2, fhPtDeltaPt, fhPtDeltaZg, fhPtDeltaRg, fhPtDeltaMg);
    }
  }

  //background-subtracted jets
  for ( auto mjet1 : mjet_cont_1 ) {
    for ( auto mjet2 : mjet_container_backsub ) {
      FillResponseFromMjets(mjet1, mjet2, fhResponseBackSub);
      FillDeltasFromMjets(mjet1, mjet2, fhPtDeltaPtBackSub, fhPtDeltaZgBackSub, fhPtDeltaRgBackSub, fhPtDeltaMgBackSub);
    }
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

  fTracksCont1 = GetParticleContainer(0);
  fTracksCont1->SetClassName("AliVTrack");

  fTracksCont2 = GetParticleContainer(1);
  if (fTracksCont2) fTracksCont2->SetClassName("AliVTrack");

  if (fJetsCont1 && fJetsCont1->GetArray() == 0) fJetsCont1 = 0;
  if (fJetsCont2 && fJetsCont2->GetArray() == 0) fJetsCont2 = 0;
  if (fTracksCont1 && fTracksCont1->GetArray() == 0) fTracksCont1 = 0;
  if (fTracksCont2 && fTracksCont2->GetArray() == 0) fTracksCont2 = 0;

}

int AliAnalysisTaskSDKLResponse::FillResponseFromMjets(mjet const & mjet1, mjet const & mjet2, THnSparse* hr[4]) {

  auto dist = CalcDist(mjet1, mjet2);
  if (dist > 0.7) return 1; //wide distance cut-off for matching

  std::vector<split> const & splits_1 = mjet1.splits;
  std::vector<split> const & splits_2 = mjet2.splits;
  auto pt1 = mjet1.pt;
  auto pt2 = mjet2.pt;

  auto nsd_1 = splits_1.size();
  float zg_1[4] = {0.};
  float rg_1[4] = {0.};
  float mg_1[4] = {-1.};
  int sd_step_1[4] = {0};
  for (int i = 0; (i < 4) && (i < nsd_1); i++) {
    zg_1[i] = splits_1[i].z;
    rg_1[i] = splits_1[i].r;
    mg_1[i] = splits_1[i].m;
    sd_step_1[i] = splits_1[i].sd_step;
  }

  auto nsd_2 = splits_2.size();
  float zg_2[4] = {0.};
  float rg_2[4] = {0.};
  float mg_2[4] = {-1.};
  int sd_step_2[4] = {0};
  for (int i = 0; (i < 4) && (i < nsd_2); i++) {
    zg_2[i] = splits_2[i].z;
    rg_2[i] = splits_2[i].r;
    mg_2[i] = splits_2[i].m;
    sd_step_2[i] = splits_2[i].sd_step;
  }

//  fhDist->Fill(dist);
  Double_t xvalue[4][15]; //4 hard splits
  for (int s = 0; s < 4; s++) {
    xvalue[s][0] = dist;

    xvalue[s][1] = pt1;
    xvalue[s][2] = zg_1[s];
    xvalue[s][3] = rg_1[s];
    xvalue[s][4] = mg_1[s];
    xvalue[s][5] = mg_1[s] / pt1; //todo, check

    xvalue[s][6] = nsd_1;
    xvalue[s][7] = sd_step_1[s];

    xvalue[s][8] = pt2;
    xvalue[s][9] = zg_2[s];
    xvalue[s][10] = rg_2[s];
    xvalue[s][11] = mg_2[s];
    xvalue[s][12] = mg_2[s] / pt2; //todo, check

    xvalue[s][13] = nsd_2;
    xvalue[s][14] = sd_step_2[s]; //testing, one value only

    hr[s]->Fill(xvalue[s]); //response - split - ordering
  }

  return 0;

}


int AliAnalysisTaskSDKLResponse::FillDeltasFromMjets(mjet const & mjet1, mjet const & mjet2, TH2F *hptdpt, TH2F *h1[4], TH2F *h2[4], TH2F *h3[4]) {

  auto dist = CalcDist(mjet1, mjet2);
  if ( dist > (0.4*0.8) ) return 1; //tight cut - only "matched" jets

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
      mjet_container.push_back( mjet{jet->Pt(), jet->Eta(), jet->Phi(), jet->Area(), AliAnalysisTaskSDKL::ReclusterFindHardSplits(jet)} );
    }
  }
}

void AliAnalysisTaskSDKLResponse::FillMjetContainer(std::vector <fastjet::PseudoJet> const & jet_container, std::vector <mjet> & mjet_container) {
  //don't fill ghosts
  for ( auto jet : jet_container ) {
    if ( jet.perp() < 1.e-6 ) continue; //skip ghost-only jets
    double area = 0.0;
    if ( jet.has_area() ) area = jet.area_4vector().perp();
    mjet_container.push_back( mjet{jet.perp(), jet.eta(), jet.phi(), area, AliAnalysisTaskSDKL::ReclusterFindHardSplits(jet)} );
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