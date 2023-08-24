#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
//#include <TH3F.h>
#include <THnSparse.h>
//#include <TTree.h>
#include <TNtuple.h>
#include <TList.h>
#include <TRandom.h>
//#include <TLorentzVector.h>

#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliAODMCParticle.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSDKL.h"

ClassImp(AliAnalysisTaskSDKL)

//________________________________________________________________________
AliAnalysisTaskSDKL::AliAnalysisTaskSDKL(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fhAll(0),
  fhAllBackSub(0),
  fhRho(0),
  fhRhoSparse(0),
  fTree(0),
  fTreeBackSub(0),
  fJetsCont(0),
  fTracksCont(0),
  fbcoption(0),
  fCSOption(0),
  fCSAlpha(0.0),
  fCSDeltaRmax(0.25),
  fCSubtractor(0),
  fCSubtractorCS(0)
{
  // Default constructor.
  // SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskSDKL::AliAnalysisTaskSDKL(const char *name, Int_t const backgroption) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fhAll(0),
  fhAllBackSub(0),
  fhRho(0),
  fhRhoSparse(0),
  fTree(0),
  fTreeBackSub(0),
  fJetsCont(0),
  fTracksCont(0),
  fbcoption(backgroption),
  fCSOption(0),
  fCSAlpha(0.0),
  fCSDeltaRmax(0.25),
  fCSubtractor(0),
  fCSubtractorCS(0)
{
  // Standard constructor.
  // SetMakeGeneralHistograms(kTRUE);
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskSDKL::~AliAnalysisTaskSDKL()
{
  // Destructor.
  if (fTree) delete fTree;
  if (fTreeBackSub) delete fTreeBackSub;
}

AliAnalysisTaskSDKL* AliAnalysisTaskSDKL::AddTaskSoftDrop(
  const char *ntracks,
  const char *njets,
  const char *nrho,
  Int_t       nCentBins,
  Double_t    jetradius,
  Double_t    jetptcut,
  Double_t    jetareacut,
  const char *type,
  Int_t       backgroption,
  Int_t       leadhadtype,
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
  if (strcmp(njets,"")) {
    name += "_";
    name += njets;
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

  AliAnalysisTaskSDKL* jetTask = new AliAnalysisTaskSDKL(name, backgroption);
  //jetTask->SetCentRange(0.,100.);
  //jetTask->SetNCentBins(nCentBins);
  //AliMultSelection *multSelection = static_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
  //if(multSelection) centrality = multSelection->GetMultiplicityPercentile("V0M");

  AliParticleContainer* partCont = 0;
  if (trackName == "mcparticles") {
    partCont = jetTask->AddMCParticleContainer("mcparticles");
  }
  else if (trackName == "usedefault" || trackName == "tracks" || trackName == "Tracks") {
    partCont = jetTask->AddTrackContainer("tracks");
  }
//  else if (!trackName.IsNull()) {
//    partCont = new AliParticleContainer(trackName);
//  }

  TString strType(type);
  AliJetContainer *jetCont = jetTask->AddJetContainer(njets,strType,jetradius);
  if(jetCont) {
    jetCont->SetRhoName(nrho);
//    jetCont->ConnectParticleContainer(partCont);
//    jetCont->SetZLeadingCut(0.98,0.98);
    jetCont->SetPercAreaCut(jetareacut);
    jetCont->SetJetPtCut(jetptcut);
    jetCont->SetLeadingHadronType(leadhadtype);
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
  contname2 += "_tree";
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contname2.Data(),
                                                            TTree::Class(),AliAnalysisManager::kOutputContainer,
                                                            Form("%s", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  mgr->ConnectOutput (jetTask, 1, coutput1 );
  mgr->ConnectOutput (jetTask, 2, coutput2 );

  return jetTask;
}

//________________________________________________________________________
void AliAnalysisTaskSDKL::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  const Int_t nzgbins = 10;
  Double_t zg_min = 0.0;
  Double_t zg_max = 0.5;

  const Int_t nrgbins = 10;
  Double_t rg_min = 0.0;
  Double_t rg_max = 0.4;

  const Int_t nmgbins = 10;
  Double_t mg_min = 0.0;
  Double_t mg_max = 0.3;

                      //pt    n_sd nsteps_1  zg1       rg1       mg1      zg2       rg2       mg2      zg3       rg3       mg3      zg4       rg4       mg4
  Int_t     bins[15] = {15,   10,  10,       nzgbins,  nrgbins,  nmgbins, nzgbins,  nrgbins,  nmgbins, nzgbins,  nrgbins,  nmgbins, nzgbins,  nrgbins,  nmgbins};
  Double_t  xmin[15] = {0.,   0,   0,        zg_min,   rg_min,   mg_min,  zg_min,   rg_min,   mg_min,  zg_min,   rg_min,   mg_min,  zg_min,   rg_min,   mg_min};
  Double_t  xmax[15] = {150., 10,  10,       zg_max,   rg_max,   mg_max,  zg_max,   rg_max,   0.2,     zg_max,   0.3,      0.1,     zg_max,   0.3,      0.1};
  fhAll = new THnSparseD("hAll", "hAll", 15, bins, xmin, xmax);
  fOutput->Add(fhAll);

  fhAllBackSub = new THnSparseD("hAllBackSub", "hAllBackSub", 15, bins, xmin, xmax);
  fOutput->Add(fhAllBackSub);

  fhRho = new TH1F("fhRho","fhRho",1000,0,1000);
  fOutput->Add(fhRho);

  fhRhoSparse = new TH1F("fhRhoSparse","fhRhoSparse",1000,0,1000);
  fOutput->Add(fhRhoSparse);

  fTree = new TNtuple("JetTrackTree", "jet track tree", "pt:eta:phi:jetm");
//  fOutput->Add(fTree);

  fTreeBackSub = new TNtuple("JetTrackTreeBackSub", "jet track tree back sub", "pt:eta:phi:jetm");
  PostData(2, fTreeBackSub);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.

}

//________________________________________________________________________
Bool_t AliAnalysisTaskSDKL::FillHistograms() {
  // Fill histograms.
  float max_jet_pt = 0.0;

  if (fJetsCont) {
    for (auto jet : fJetsCont->accepted()) {
      auto jet_pt = jet->Pt();
      if (jet_pt > max_jet_pt) max_jet_pt = jet_pt;
    }
  }

  //hardest jet is too soft
  if (max_jet_pt < 5.) return kTRUE;
  //tighter cut for AA
  if ( (1==fbcoption) && (max_jet_pt < 30.) ) return kTRUE;

  //fill full event
  std::vector <fastjet::PseudoJet> event_full;
  AddTracksToEvent(fTracksCont, event_full);

  if (fbcoption >= 0) {
    //get backgr-subtracted jets
    Double_t rho;
    Double_t rho_sparse;
    InitializeSubtractor(event_full, rho, rho_sparse, fbcoption);
    fhRho->Fill(rho);
    fhRhoSparse->Fill(rho_sparse);
  }

  std::vector<fastjet::PseudoJet> jets_backsub = GetBackSubJets(event_full);

  std::vector<fastjet::PseudoJet> jets_backsub_filtered;
  FilterJets(jets_backsub, jets_backsub_filtered, 5.0);

//  for (auto & jet : jets_backsub_filtered) {
//    std::vector<split> splits = ReclusterFindHardSplits(jet);
//    FillSparseFromSplits( fhAllBackSub, splits, jet.pt() );
//  }

  FillTree(jets_backsub_filtered, fTreeBackSub);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.

  if (fCSubtractor)   delete fCSubtractor;
  if (fCSubtractorCS) delete fCSubtractorCS;

  fCSubtractor = 0;
  fCSubtractorCS = 0;

  return kTRUE;
}

std::vector<fastjet::PseudoJet> AliAnalysisTaskSDKL::GetBackSubJets(std::vector <fastjet::PseudoJet> const & event_full) {

  fastjet::Selector sel_jets = fastjet::SelectorAbsEtaMax(0.9 - 0.4); //max_eta_jet
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4, fastjet::E_scheme);
  fastjet::AreaDefinition area_def( fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(0.9,1) ); //0.9 -> max_eta

  std::vector<fastjet::PseudoJet> jets_backsub;

  if (fbcoption < 0) {

    fCSubtractorCS = new fastjet::ClusterSequenceArea(event_full, jet_def, area_def);
    jets_backsub = sel_jets( fCSubtractorCS->inclusive_jets() );

    return jets_backsub;
  }

  if (fCSubtractor) {

    if (0 == fCSOption) { //full event: subtract then cluster
      std::vector<fastjet::PseudoJet> event_backsub = fCSubtractor->subtract_event(event_full, 0.9);
      fCSubtractorCS = new fastjet::ClusterSequenceArea(event_backsub, jet_def, area_def);
      jets_backsub = sel_jets( fCSubtractorCS->inclusive_jets() );
    }
    else if (1 == fCSOption) { //jet-by-jet: cluster then subtract
      fCSubtractorCS = new fastjet::ClusterSequenceArea(event_full, jet_def, area_def);
      std::vector<fastjet::PseudoJet> jets_full = sel_jets( fCSubtractorCS->inclusive_jets() );
      for (auto & j : jets_full) {
        jets_backsub.push_back( fCSubtractor->result(j) );
      }
    }

  }

  return jets_backsub;

}

std::vector<split> AliAnalysisTaskSDKL::ReclusterFindHardSplits(AliEmcalJet *jet) {

  std::vector <fastjet::PseudoJet> particles;
  UShort_t ntracks = jet->GetNumberOfTracks();
  for (int j = 0; j < ntracks; j++) {
    particles.push_back(fastjet::PseudoJet(jet->Track(j)->Px(), jet->Track(j)->Py(), jet->Track(j)->Pz(), jet->Track(j)->E()));
  }

  return ReclusterFindHardSplits(particles);

}

std::vector<split> AliAnalysisTaskSDKL::ReclusterFindHardSplits(fastjet::PseudoJet const & jet) {

  //remove ghost particles
  std::vector <fastjet::PseudoJet> particles;
  for (auto constituent : jet.constituents() ) {
    if ( constituent.perp() > 1.e-5 ) particles.push_back( constituent );
  }

  return ReclusterFindHardSplits(particles);

}

std::vector<split> AliAnalysisTaskSDKL::ReclusterFindHardSplits(std::vector <fastjet::PseudoJet> const & particles) {

  fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 0.4, fastjet::E_scheme); //scheme
  fastjet::ClusterSequence cs(particles, jet_def);
  std::vector <fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

  std::vector<split> splits;
  if (jets.size() > 0) {
    splits = FindHardSplits( jets[0] );
  }

  return splits;

}

std::vector<split> AliAnalysisTaskSDKL::FindHardSplits(fastjet::PseudoJet const & jet) {

  auto hardest_subjet = jet;

  fastjet::PseudoJet jet1;
  fastjet::PseudoJet jet2;

  std::vector<split> splits;

  int sd_step = 0;

  while ( hardest_subjet.has_parents(jet1, jet2) ) {

    sd_step++;
    hardest_subjet = jet1;

    Float_t pt1 = jet1.pt();
    Float_t pt2 = jet2.pt();

    Float_t z = pt2/(pt1 + pt2);
    Float_t r = TMath::Sqrt( jet1.plain_distance(jet2) );
    Float_t m = hardest_subjet.m();

    if (z > 0.1) {
      splits.push_back( split{z,r,m,sd_step} );
    }

  }

  return splits;

}


void AliAnalysisTaskSDKL::FillSparseFromSplits(THnSparse *histo, std::vector<split> const & splits, const double jet_pt) {

  int nsd = splits.size();

  float zg[4] = {0.0};
  float rg[4] = {0.0};
  float mg[4] = {-1.0};
  int sd_step[4] = {0};
  for (int i = 0; (i < 4) && (i < nsd); i++) {
    zg[i] = splits[i].z;
    rg[i] = splits[i].r;
    mg[i] = splits[i].m;
    sd_step[i] = splits[i].sd_step;
  }

  //pt n_sd nsteps_1 zg1 rg1 mg1 zg2 rg2 mg2 zg3 rg3 mg3 zg4 rg4 mg4
  Double_t xvalue[15];
  xvalue[0] = jet_pt;
  xvalue[1] = splits.size();

  xvalue[2] = sd_step[0]; //for the moment nsd_steps only for the first split

  xvalue[3] = zg[0];
  xvalue[4] = rg[0];
  xvalue[5] = mg[0]/jet_pt;

  xvalue[6] = zg[1];
  xvalue[7] = rg[1];
  xvalue[8] = mg[1]/jet_pt;

  xvalue[9]  = zg[2];
  xvalue[10] = rg[2];
  xvalue[11] = mg[2]/jet_pt;

  xvalue[12] = zg[3];
  xvalue[13] = rg[3];
  xvalue[14] = mg[3]/jet_pt;

  histo->Fill( xvalue );

}

int AliAnalysisTaskSDKL::InitializeSubtractor(std::vector <fastjet::PseudoJet> const & event_full, Double_t & rho, Double_t & rho_sparse, Int_t opt) {

  double max_eta = 0.9;
  double max_eta_jet = 0.5;

  fastjet::Selector sel_jets = fastjet::SelectorAbsEtaMax(max_eta_jet);

  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(max_eta, 1));

  fastjet::JetDefinition jet_def_kt(fastjet::kt_algorithm, 0.4);
  fastjet::ClusterSequenceArea cs_kt(event_full, jet_def_kt, area_def);
  std::vector <fastjet::PseudoJet> jets_kt = sorted_by_pt(sel_jets(cs_kt.inclusive_jets()));

  std::vector <fastjet::PseudoJet> jets_akt;
  //option 0 (pPb)
  //tag matched kt-antikt jets
  if (0 == opt) {

    fastjet::JetDefinition jet_def_akt(fastjet::antikt_algorithm, 0.4);
    fastjet::ClusterSequenceArea cs_akt(event_full, jet_def_akt, area_def);
    jets_akt = sorted_by_pt(sel_jets(cs_akt.inclusive_jets()));

    for (auto & jet_akt : jets_akt) {
      auto jet_akt_pt = jet_akt.perp();
      if (jet_akt_pt < 5.) continue;
      for (auto & jet_kt : jets_kt) {
        auto dist = jet_akt.delta_R(jet_kt);
        if (dist < 0.4 * 0.6) {
          jet_kt.set_user_index(7); //tagging kt jets matched to anti-kt
        }
      }
    }

    jets_akt.clear();
  }

  //option 1 (PbPb)
  //exclude two hardest kt jets
  if (1 == opt) {
    if (jets_kt.size() < 3) return 1; //vector of kt-jets has only 2 elements
    else {
      jets_kt[0].set_user_index(7);
      jets_kt[1].set_user_index(7);
    }
  }

  std::vector <Double_t> vec_pt_area;
  Double_t area_sum = 0.0;
  for (auto & jet : jets_kt) {

    if ( 7==jet.user_index() ) continue;

    //exclude ghost-only jets
    if ( jet.has_area() && (jet.pt() > 0.001) )  {
      auto jarea = jet.area_4vector().perp();
      if (jarea > 0) {
        vec_pt_area.push_back( jet.perp() / jarea ); //test
        area_sum += jarea;
      }
    }

  }

  sort( vec_pt_area.begin(), vec_pt_area.end() );
  auto size_vec_pt_area = vec_pt_area.size();
  if (size_vec_pt_area > 0) {
    if (size_vec_pt_area % 2 == 0) {
      rho = (vec_pt_area[size_vec_pt_area / 2 - 1] + vec_pt_area[size_vec_pt_area / 2]) / 2.0;
    } else {
      rho = vec_pt_area[size_vec_pt_area / 2];
    }
  }
  vec_pt_area.clear();

  Double_t const full_acceptance = 2.0 * 0.9 * 2.0 * TMath::Pi();
  auto acc_correction = area_sum/full_acceptance;

  rho_sparse = rho * acc_correction;

  if (rho > 0.001 && rho_sparse > 0.001) {
    //provide external rho and rhom
    Double_t rho_C;
    if (0==opt) rho_C = rho_sparse; //pA
    else        rho_C = rho;        //AA
    Double_t const rho_m_C = 0.0; //massless tracks
    fCSubtractor = new fastjet::contrib::ConstituentSubtractor(rho_C, rho_m_C, 0.0, 0.25);
    fCSubtractor->set_alpha(fCSAlpha);
    fCSubtractor->set_max_standardDeltaR(fCSDeltaRmax);
    fCSubtractor->set_ghost_area(0.01);
    return 0;
  }
  else return 1;

}

void AliAnalysisTaskSDKL::AddTracksToEvent(AliParticleContainer* cont, std::vector <fastjet::PseudoJet> & event, Double_t const efficiency, TRandom* fRandom) {

  if (cont) {
    cont->ResetCurrentID();
    for ( auto track : cont->accepted() ) {

      if ( (efficiency < 1.0) && fRandom ) {
        if (fRandom->Uniform() > efficiency) continue;
      }

      Double_t track_p[3];
      track->PxPyPz(track_p);
      Double_t ptot = sqrt( track_p[0]*track_p[0] + track_p[1]*track_p[1] + track_p[2]*track_p[2] );
      fastjet::PseudoJet pj(track_p[0], track_p[1], track_p[2], ptot); //massless particle
      pj.set_user_index( track->GetLabel() );
      event.push_back( pj );
    }
  }

}

void AliAnalysisTaskSDKL::FillTree(std::vector<fastjet::PseudoJet> const & jets, TNtuple* tree) {

  for (auto & jet : jets) {

    if (fbcoption < 0) {
      if (jet.pt() < 5.) continue; //hard-coded jet-pt cut
    }
    else {
      if (jet.pt() < 30.) continue; //hard-coded jet-pt cut
    }

    int nconst = 0;
    for (auto c : jet.constituents() ) {
      if ( c.pt() > 1.e-5 ) nconst++;
    }
    tree->Fill(jet.pt(), jet.eta(), jet.phi(), nconst);
    for (auto c : jet.constituents() ) {
      if ( c.pt() > 1.e-5 ) {
        tree->Fill(c.pt(), c.eta(), c.phi(), -7);
      }
    }

  }

}

void AliAnalysisTaskSDKL::FillTree(AliJetContainer *jets, TNtuple* tree) {

  if (jets) {
    for (auto jet : jets->accepted()) {

      if (fbcoption < 0) {
        if (jet->Pt() < 5.) continue; //hard-coded jet-pt cut
      }
      else {
        if (jet->Pt() < 30.) continue; //hard-coded jet-pt cut
      }

      UShort_t ntracks = jet->GetNumberOfTracks();
      tree->Fill(jet->Pt(), jet->Eta(), jet->Phi(), ntracks);
      for (int j = 0; j < ntracks; j++) {
        auto jtrack = jet->Track(j);
        if (jtrack->Pt() > 1.e-5) {
          tree->Fill(jtrack->Pt(), jtrack->Eta(), jtrack->Phi(), -7);
        }
      }

    }
  }

}

//________________________________________________________________________
void AliAnalysisTaskSDKL::ExecOnce() {

  AliAnalysisTaskEmcalJet::ExecOnce();

  fJetsCont   = GetJetContainer(0);
  fTracksCont = GetParticleContainer(0);

//  TString trackName(ntracks);
//  if (trackName == "mcparticles") {
//    fTracksCont->SetClassName("AliAODMCParticle");
//  }
//  else if (trackName == "usedefault" || trackName == "tracks" || trackName == "Tracks") {
  fTracksCont->SetClassName("AliVTrack");
//  }

  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskSDKL::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
void AliAnalysisTaskSDKL::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

void AliAnalysisTaskSDKL::FilterJets(std::vector<fastjet::PseudoJet> const & jets, std::vector<fastjet::PseudoJet> & jets_filtered, Float_t pt_cut) {
  for (auto & jet : jets) {
    if ( jet.pt() < pt_cut ) continue;
    if ( fabs(jet.eta()) > 0.5 ) continue;
    if ( jet.has_area() ) {
      auto jarea = jet.area_4vector().perp();
      if (jarea < 0.3016) continue;
    }
    jets_filtered.push_back(jet);
  }
}
