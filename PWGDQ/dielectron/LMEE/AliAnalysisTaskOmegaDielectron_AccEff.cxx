#include "AliAnalysisTaskOmegaDielectron_AccEff.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskOmegaDielectron_AccEff);
/// \endcond

/***************************************************************************/ /**
                                                                               * ROOT I/O Constructor.
                                                                               ******************************************************************************/
AliAnalysisTaskOmegaDielectron_AccEff::AliAnalysisTaskOmegaDielectron_AccEff()
    : AliAnalysisTaskSE(),
      // General member variables
      fOutputList(nullptr), fEvent(nullptr), fMCEvent(nullptr), fEventCuts(),
      fESDtrackCuts(nullptr), fCutMode(100),
      // Toggles
      fIsESD(kTRUE), fIsMC(kFALSE), fUseCent(kFALSE),
      // Cut Parameters
      fTriggerMask(AliVEvent::kMB | AliVEvent::kINT7), fMinEta(-0.8),
      fMaxEta(0.8), fMinPt(0.2), fMaxPt(10.0), fMaxZv(10.0), fMinCent(0.0),
      fMaxCent(100.0),
      // Arrays for Binning
      fBinsMult(nullptr), fBinsCent(nullptr), fBinsPt(nullptr), fBinsy(nullptr),
      fBinsEta(nullptr), fBinsZv(nullptr), fBinsPtReso(nullptr),
      //pdg codes:
      felectron_pdg(11),fpositron_pdg(-11),
      fmother_pdg(223),
      // fmother_pdg(333),
      //storage Vectors
      v_elec_true_omega(0x0),
      v_posi_true_omega(0x0),
      v_elec_motherID_true_omega(0),
      v_posi_motherID_true_omega(0),
      // Event-Histograms
      fHistEventSelection(nullptr),
      fHistMC_Omegas_Rapidity(nullptr),
      fHistMC_ele1_posi2_OmegaDielDeacay(nullptr),
      fHistMC_Omegas_gen(nullptr),
      fHistMC_Omegas_gen_DaughtersinAcc(nullptr),
      fHist_rec_true_Ele_Omegas_Mothers(nullptr),
      fHist_rec_true_Pos_Omegas_Mothers(nullptr),
      fHist_rec_true_Dielec(nullptr)
{
  // ROOT IO constructor, don't allocate memory here!
}

/***************************************************************************/ /**
                                                                               * Constructor.
                                                                               ******************************************************************************/
AliAnalysisTaskOmegaDielectron_AccEff::AliAnalysisTaskOmegaDielectron_AccEff(const char *name)
    : AliAnalysisTaskSE(name),
      // General member variables
      fOutputList(nullptr), fEvent(nullptr), fMCEvent(nullptr), fEventCuts(),
      fESDtrackCuts(nullptr), fCutMode(100),
      // Toggles
      fIsESD(kTRUE),
      // fIsMC(kTRUE),
      fIsMC(kFALSE), fUseCent(kFALSE),
      // Cut Parameters
      fTriggerMask(AliVEvent::kMB | AliVEvent::kINT7), fMinEta(-0.8),
      fMaxEta(0.8), fMinPt(0.2), fMaxPt(10.0), fMaxZv(10.0), fMinCent(0.0),
      fMaxCent(100.0),
      // Arrays for Binning
      fBinsMult(nullptr), fBinsCent(nullptr), fBinsPt(nullptr), fBinsy(nullptr),
      fBinsEta(nullptr), fBinsZv(nullptr), fBinsPtReso(nullptr),
      //pdg codes:
      felectron_pdg(11),fpositron_pdg(-11),
      fmother_pdg(223),
      // fmother_pdg(333),
      //storage Vectors
      v_elec_true_omega(0x0),
      v_posi_true_omega(0x0),
      v_elec_motherID_true_omega(0),
      v_posi_motherID_true_omega(0),
      // Event-Histograms
      fHistEventSelection(nullptr),
      fHistMC_Omegas_Rapidity(nullptr),
      fHistMC_ele1_posi2_OmegaDielDeacay(nullptr),
      fHistMC_Omegas_gen(nullptr),
      fHistMC_Omegas_gen_DaughtersinAcc(nullptr),
      fHist_rec_true_Ele_Omegas_Mothers(nullptr),
      fHist_rec_true_Pos_Omegas_Mothers(nullptr),
      fHist_rec_true_Dielec(nullptr)
{
  // Set default binning
  Double_t binsMultDefault[2] = {0., 10000.};
  Double_t binsCentDefault[2] = {0., 100.};
  Double_t binsPtDefault[53] = {
      0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45, 0.5, 0.55, 0.6,
      0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95, 1.0,  1.1, 1.2,  1.3,
      1.4,  1.5,  1.6,  1.7,  1.8,  1.9,  2.0,  2.2,  2.4, 2.6,  2.8,
      3.0,  3.2,  3.4,  3.6,  3.8,  4.0,  4.5,  5.0,  5.5, 6.0,  6.5,
      7.0,  8.0,  9.0,  10.0, 20.0, 30.0, 40.0, 50.0, 60.0};
  Double_t binsEtaDefault[19] = {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3,
                                 -0.2, -0.1, 0.,   0.1,  0.2,  0.3,  0.4,
                                 0.5,  0.6,  0.7,  0.8,  0.9};
  Double_t binsyDefault[19] = {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3,
                                -0.2, -0.1, 0.,   0.1,  0.2,  0.3,  0.4,
                                0.5,  0.6,  0.7,  0.8,  0.9};

  Double_t binsZvDefault[13] = {-30., -25., -20., -15., -10., -5., 0.,
                                5.,   10.,  15.,  20.,  25.,  30.};


  // binning for relative pT resolution
  const Int_t nBinsPtReso = 300;
  Double_t binsPtReso[nBinsPtReso + 1];
  SetFixedBinEdges(binsPtReso, 0., 0.3, nBinsPtReso);
  SetBinsPtReso(nBinsPtReso, binsPtReso);

  SetBinsMult(1, binsMultDefault);
  SetBinsCent(1, binsCentDefault);
  SetBinsPt(52, binsPtDefault);
  SetBinsEta(18, binsEtaDefault);
  SetBinsy(18, binsyDefault);
  SetBinsZv(12, binsZvDefault);

  DefineOutput(1, TList::Class());
}

/***************************************************************************/ /**
                                                                               * Function executed once before the event loop. Create histograms here.
                                                                               ******************************************************************************/
void AliAnalysisTaskOmegaDielectron_AccEff::UserCreateOutputObjects() {

  OpenFile(1, "recreate");
  fOutputList = new TList();
  fOutputList->SetOwner();

  // Control histogram to check the effect of event cuts
  fHistEventSelection = new TH1F("fHistEventSelection", "fHistEventSelection [all : selected]", 2, 0.5, 2.5);
  fHistEventSelection->GetYaxis()->SetTitle("#it{N}_{events}");
  fHistEventSelection->GetXaxis()->SetBinLabel(1, "all");
  fHistEventSelection->GetXaxis()->SetBinLabel(2, "selected");
  fOutputList->Add(fHistEventSelection);

  fHistMC_Omegas_Rapidity = new TH1F("fHistMC_Omegas_Rapidity", "Omega - Rapidity", 200, -10, 10);
  fHistMC_Omegas_Rapidity->GetYaxis()->SetTitle("#it{N}_{events}");
  fHistMC_Omegas_Rapidity->GetXaxis()->SetTitle("y");
  fOutputList->Add(fHistMC_Omegas_Rapidity);

  // if(fIsMC)
  // {
  fHistMC_ele1_posi2_OmegaDielDeacay = new TH1F("fHistMC_ele1_posi2_OmegaDielDeacy", "Omega - Dielectron Deacay", 2, 0.5, 2.5);
  fHistMC_ele1_posi2_OmegaDielDeacay->GetYaxis()->SetTitle("N");
  fHistMC_ele1_posi2_OmegaDielDeacay->GetXaxis()->SetBinLabel(1, "electrons");
  fHistMC_ele1_posi2_OmegaDielDeacay->GetXaxis()->SetBinLabel(2, "positrons");
  fOutputList->Add(fHistMC_ele1_posi2_OmegaDielDeacay);
  // }

  fHistMC_Omegas_gen = new TH2F("fHistMC_Omegas_gen", "fHistMC_Omegas_gen;#it{p}_{T} (GeV/#it{c});y", 50, 0, 10, 50, -10, 10);
  fOutputList->Add(fHistMC_Omegas_gen);
  fHistMC_Omegas_gen_DaughtersinAcc = new TH2F("fHistMC_Omegas_gen_DaughtersinAcc", "fHistMC_Omegas_gen_DaughtersinAcc;#it{p}_{T} (GeV/#it{c});y", 50, 0, 10, 50, -10, 10);
  fOutputList->Add(fHistMC_Omegas_gen_DaughtersinAcc);

  fHist_rec_true_Ele_Omegas_Mothers = new TH1F("Reconstructed_True_Electrons_Omegas_Mothers_Pt", "Reconstructed_True_Electrons_Omegas_Mothers_Pt;#it{p}_{T} (GeV/#it{c})", 50, 0, 10);
  fHistMC_Omegas_Rapidity->GetYaxis()->SetTitle("#it{N}_{events}");
  fOutputList->Add(fHist_rec_true_Ele_Omegas_Mothers);
  fHist_rec_true_Pos_Omegas_Mothers = new TH1F("Reconstructed_True_Positrons_Omegas_Mothers_Pt", "Reconstructed_True_Positrons_Omegas_Mothers_Pt;#it{p}_{T} (GeV/#it{c})", 50, 0, 10);
  fHistMC_Omegas_Rapidity->GetYaxis()->SetTitle("#it{N}_{events}");
  fOutputList->Add(fHist_rec_true_Pos_Omegas_Mothers);
  fHist_rec_true_Dielec = new TH2D("Reconstructed_True_Omegas_DielectronD_InvMass_Pt", "Reconstructed_True_Omegas_DielectronD_InvMass_Pt;Inv_Mass;#it{p}_{T,ee} (GeV/#it{c})", 50, 0, 10, 50, 0, 10);
  fHistMC_Omegas_Rapidity->GetYaxis()->SetTitle("#it{N}_{events}");
  fOutputList->Add(fHist_rec_true_Dielec);



  // override event automatic event selection settings


  ///??????????????///, 60, -0.5, 6 - 0.5);
  fEventCuts.SetMaxVertexZposition(fMaxZv);
  if (fUseCent)
    fEventCuts.SetCentralityRange(fMinCent, fMaxCent);
  fEventCuts.OverrideAutomaticTriggerSelection(fTriggerMask);

    ///??????????????///
  if (fIsESD)
    InitESDTrackCuts();

  PostData(1, fOutputList);
}

/// Destructor
AliAnalysisTaskOmegaDielectron_AccEff::~AliAnalysisTaskOmegaDielectron_AccEff() {
  if (fESDtrackCuts) {
    delete fESDtrackCuts;
    fESDtrackCuts = nullptr;
  }
}

/***************************************************************************/ /**
                                                                               * Function executed for each event.
                                                                               ******************************************************************************/
void AliAnalysisTaskOmegaDielectron_AccEff::UserExec(Option_t *) {

  fEvent = InputEvent();
  if (!fEvent) {
    AliError("fEvent not available\n");
    return;
  }

  // if(fIsMC){
  // AliMCEvent* fMCEvent = MCEvent();
  fMCEvent = MCEvent();
  if (!fMCEvent) {
    AliError("fMCEvent not available\n");
    return;
  }
  // }


  ///-------------------- Fill Events in Histogramm :: All;Selected;???? --------------------///

  fHistEventSelection->Fill(1.0); // all events
  if (!fEventCuts.AcceptEvent(fEvent)) // https://github.com/alisw/AliPhysics/blob/master/OADB/AliEventCuts.cxx
    return;
  fHistEventSelection->Fill(2.0); // selected events


  ///-------------------- Loop over Generated MC Particles ---OMEGAS --------------------///
  // if (fIsMC){
  double parity = 0.8;

  for(Int_t iGenPart = 0; iGenPart < fMCEvent->GetNumberOfTracks();iGenPart++) {
    AliMCParticle* mcGenParticle  = (AliMCParticle*)fMCEvent->GetTrack(iGenPart);
    if(!mcGenParticle){
      AliError("mcGenParticle  not available\n");
      continue;
    }

    Int_t pdgcode_gen = mcGenParticle->PdgCode();
    Int_t mother_id_gen = ( mcGenParticle->GetLabel() );

    if(pdgcode_gen == fmother_pdg){
      fHistMC_Omegas_Rapidity->Fill(mcGenParticle->Y());
    }


    if ( TMath::Abs( mcGenParticle->Y() ) >= parity ) continue;

    if(pdgcode_gen == fmother_pdg){
      if( CheckDielectronDecay(mcGenParticle) ){
        fHistMC_Omegas_gen->Fill(mcGenParticle->Pt(), mcGenParticle->Y());
        // std::cout << "mother_id Gen: " << mother_id_gen << '\n';
      }

      if( CheckDielectronDecay_DaughterinAcc(mcGenParticle) ){
        fHistMC_Omegas_gen_DaughtersinAcc->Fill(mcGenParticle->Pt(), mcGenParticle->Y());
        // std::cout << "mother_id Gen with daughters within /eta/<0.8 : " << mother_id_gen << '\n';
      }

      else{
        continue;
      }
    }
  }




  ///-------------------- Loop over AOD Tracks --------------------///

  AliAODTrack *track = nullptr;
  // Create a vector containing the electron Track
  v_elec_true_omega.clear();
  v_elec_motherID_true_omega.clear();
  v_posi_true_omega.clear();
  v_posi_motherID_true_omega.clear();

  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); iTrack++) {
    track = static_cast<AliAODTrack*>(fEvent->GetTrack(iTrack));
    if (!track || !track->TestFilterBit(128)){  // filter bit 128 denotes TPC-only tracks, use only them for the analysis
      // AliErrorF("Could not receive track %d\n", iTrack);
      continue;
    }
    track = static_cast<AliAODTrack*>(track);
    if (!AcceptKinematics(track)){
      continue;
    }



    ///--- Find original particle in MC-Stack ---///
    // if(fIsMC){
    Int_t mcLabel =  TMath::Abs(track->GetLabel()); // negative label means bad quality track
    AliMCParticle *mcParticle = (AliMCParticle *)fMCEvent->GetTrack(mcLabel);
    if (!mcParticle) {
      AliError("mcParticle not available\n");
      continue;
    }
    if ( !mcParticle->IsPhysicalPrimary() ){
      // AliError("mcParticle not physical primary\n");
      continue;
    }

    Int_t pdgcode = mcParticle->PdgCode();
    AliMCParticle* mcParticle_mother  = (AliMCParticle*)fMCEvent->GetTrack( mcParticle->GetMother() );
    Int_t mother_id = ( mcParticle_mother->GetLabel() );
    Int_t pdgCode_mother = ( mcParticle_mother->PdgCode() );

////***************ELECTRONS*****************////
    if(pdgCode_mother == fmother_pdg){
      // std::cout << "Mother is Omega" << '\n';

      if(pdgcode == felectron_pdg){
        if(CheckDielectronDecay(mcParticle_mother)){
          fHistMC_ele1_posi2_OmegaDielDeacay->Fill(1.0);
          v_elec_true_omega.push_back(track);
          v_elec_motherID_true_omega.push_back(mother_id);
          std::cout << "--------::  electron track: " << iTrack << "  electron label: " << mcLabel << "  electron mother: " << mother_id << '\n';
        }
      }
////***************POSITRONS*****************////
      if(pdgcode == fpositron_pdg){
        if(CheckDielectronDecay(mcParticle_mother)){
          fHistMC_ele1_posi2_OmegaDielDeacay->Fill(2.0);
          v_posi_true_omega.push_back(track);
          v_posi_motherID_true_omega.push_back(mother_id);
          std::cout << "--------::  positron track: " << iTrack << "  positron label: " << mcLabel << "  positron mother: " << mother_id << '\n';
        }
      }
      else{
        continue;
      }
    }

  }

  // //Check vectors
  // if (v_elec_true_omega.size() >0){
  //   for( auto const & item : v_elec_true_omega)
  //       std::cout<<item<<", BLAAAAAAAAAAAAAAAAA"<< '\n';
  // }
  // if (v_posi_true_omega.size() >0){
  //   for( auto const & item_pos : v_posi_true_omega)
  //       std::cout<<item_pos<<", BLAAAAAAAAAAAAAAAAAMMMMMMMMMMM"<< '\n';
  // }


  if (v_posi_true_omega.size() == 0) return;
  if (v_elec_true_omega.size() == 0) return;
  else{
    for(int i=0; i< (Int_t)v_elec_true_omega.size(); i++){
    // for(auto const & i : v_elec_true_omega){
      AliAODTrack * ele_from_same_mother_track= (AliAODTrack *)v_elec_true_omega.at(i);
      int mother_electron= (Int_t)v_elec_motherID_true_omega.at(i);

      for(int j=0; j< (Int_t)v_posi_true_omega.size(); j++){
      // for(auto const & j : v_posi_true_omega){
        AliAODTrack * pos_from_same_mother_track= (AliAODTrack *)v_posi_true_omega.at(j);
        int mother_positron= (Int_t)v_posi_motherID_true_omega.at(j);

        if(mother_electron==mother_positron){
          fHist_rec_true_Ele_Omegas_Mothers->Fill(ele_from_same_mother_track->Pt());
          fHist_rec_true_Pos_Omegas_Mothers->Fill(pos_from_same_mother_track->Pt());

          AliMCParticle *mcParticle_ele_trail = (AliMCParticle *)fMCEvent->GetTrack(TMath::Abs(ele_from_same_mother_track->GetLabel()));
          std::cout << " Electron pdg: " << mcParticle_ele_trail->PdgCode() << "  mclabel: " << TMath::Abs(ele_from_same_mother_track->GetLabel()) << '\n';
          AliMCParticle *mcParticle_pos_trail = (AliMCParticle *)fMCEvent->GetTrack(TMath::Abs(pos_from_same_mother_track->GetLabel()));
          std::cout << " Positron pdg: " << mcParticle_pos_trail->PdgCode() << "  mclabel: " << TMath::Abs(pos_from_same_mother_track->GetLabel()) << '\n';

          // Double_t electron_mass=0.00051099895;
          Double_t electron_mass = (AliPID::ParticleMass(AliPID::kElectron) );
          // std::cout << "ELECTRON  pt: " << Pt_elec << "  eta: " << Eta_elec << "  phi: " << Phi_elec << " m: " << electron_mass << '\n';
          // std::cout << "POSITRON  pt: " << Pt_posi << "  eta: " << Eta_posi << "  phi: " << Phi_posi << " m: " << electron_mass << '\n';
          TLorentzVector Lvec1;
          Lvec1.SetPtEtaPhiM(ele_from_same_mother_track->Pt(), ele_from_same_mother_track->Eta(), ele_from_same_mother_track->Phi(), electron_mass);
          TLorentzVector Lvec2;
          Lvec2.SetPtEtaPhiM(pos_from_same_mother_track->Pt(), pos_from_same_mother_track->Eta(), pos_from_same_mother_track->Phi(), electron_mass);

          TLorentzVector LvecM = Lvec1 + Lvec2;
          Double_t mass = LvecM.M();
          Double_t pairpt = LvecM.Pt();
          Double_t weight = 1;
          std::cout << "  mass: " << mass << "  pairpt: " << pairpt << '\n';

          fHist_rec_true_Dielec->Fill(mass, pairpt);
        }
        else{
          continue;
        }
      }
    }
  }






  PostData(1, fOutputList);
}

/***************************************************************************/
/*** Function executed after all events were processed.*********************
*********************************************************/
void AliAnalysisTaskOmegaDielectron_AccEff::Terminate(Option_t *) {}

/***************************************************************************/ /**
* Function to select primary charged particles.
******************************************************************************/
Bool_t AliAnalysisTaskOmegaDielectron_AccEff::IsChargedPrimary(Int_t mcLabel) {
  AliMCParticle *mcParticle = (AliMCParticle *)fMCEvent->GetTrack(mcLabel);
  if (!mcParticle->IsPhysicalPrimary()){
    return kFALSE;
  }
  if (!mcParticle) {
    AliError("mcGenParticle  not available\n");
    return kFALSE;
  }
  if (!(TMath::Abs(mcParticle->Charge()) > 0.01))
    return kFALSE;
  return kTRUE;
}

Bool_t AliAnalysisTaskOmegaDielectron_AccEff::CheckDielectronDecay(AliMCParticle *particle) {
  // Get the daughters:
  Int_t number_daughters_bool = particle->GetNDaughters();
  Int_t LabelFirstDaughter_bool = particle->GetDaughterFirst();
  Int_t LabelLastDaughter_bool = particle->GetDaughterLast();

  bool ele_from_same_mother_bool = false;
  bool pos_from_same_mother_bool = false;

  if(number_daughters_bool != 2){
    return false;
  }
  else{
    for (int daughter_i_bool = LabelFirstDaughter_bool; daughter_i_bool <= LabelLastDaughter_bool; daughter_i_bool++) {
      int pdgCode_daughter_bool = fMCEvent->GetTrack(daughter_i_bool)->PdgCode();
      if (pdgCode_daughter_bool == felectron_pdg) {
        ele_from_same_mother_bool = true;
      }
      if (pdgCode_daughter_bool == fpositron_pdg) {
        pos_from_same_mother_bool = true;
      }
    }
    // if (ele_from_same_mother_bool == true && pos_from_same_mother_bool == true && number_daughters_bool==2){
    if (ele_from_same_mother_bool && pos_from_same_mother_bool ){
      return true; // true -> if daughters ARE Dieceltrons!
      // (number_daughters==2 Wichtig! )
    }
    else
      return false; // false -> if daughters are not e+ + e-
  }
}

Bool_t AliAnalysisTaskOmegaDielectron_AccEff::CheckDielectronDecay_DaughterinAcc(AliMCParticle *particle) {
  // Get the daughters:
  Int_t number_daughters_bool = particle->GetNDaughters();
  Int_t LabelFirstDaughter_bool = particle->GetDaughterFirst();
  Int_t LabelLastDaughter_bool = particle->GetDaughterLast();

  AliMCParticle *daughter_elec = nullptr;
  AliMCParticle *daughter_pos = nullptr;

  bool ele_from_same_mother_bool = false;
  bool pos_from_same_mother_bool = false;


  if(number_daughters_bool != 2){
    return false;
  }
  else{
    for (int daughter_i_bool = LabelFirstDaughter_bool; daughter_i_bool <= LabelLastDaughter_bool; daughter_i_bool++) {
      int pdgCode_daughter_bool = fMCEvent->GetTrack(daughter_i_bool)->PdgCode();
      if (pdgCode_daughter_bool == felectron_pdg) {
        ele_from_same_mother_bool = true;
        daughter_elec = (AliMCParticle*)fMCEvent->GetTrack(daughter_i_bool);
      }
      if (pdgCode_daughter_bool == fpositron_pdg) {
        pos_from_same_mother_bool = true;
        daughter_pos = (AliMCParticle*)fMCEvent->GetTrack(daughter_i_bool);
      }
    }
    if (ele_from_same_mother_bool == true && pos_from_same_mother_bool == true && AcceptKinematics(daughter_elec) && AcceptKinematics(daughter_pos) ){ // eta und pt cut
    // if (ele_from_same_mother_bool == true && pos_from_same_mother_bool == true && EtaCut(daughter_elec) && EtaCut(daughter_pos) ){ // eta cut
      return true; // true -> if daughters ARE Dieceltrons! within Acceptance
      // (number_daughters==2 Wichtig! )
    }
    else
      return false; // false -> if daughters are not e+ + e-
  }
}

/***************************************************************************/ /**
                                                                               * Function to select if the track or mc particle is in defined kinematic range.
                                                                               ******************************************************************************/
Bool_t AliAnalysisTaskOmegaDielectron_AccEff::AcceptKinematics(AliVParticle *particle) {
  if (!particle)
    return kFALSE;

  Double_t eta = particle->Eta();
  Double_t pt = particle->Pt();

  if (eta <= fMinEta + PRECISION)
    return kFALSE;
  if (eta >= fMaxEta - PRECISION)
    return kFALSE;
  if (pt <= fMinPt + PRECISION)
    return kFALSE;
  if (pt >= fMaxPt - PRECISION)
    return kFALSE;
  return kTRUE;
}

Bool_t AliAnalysisTaskOmegaDielectron_AccEff::EtaCut(AliVParticle *particle) {
  if (!particle)
    return kFALSE;

  Double_t eta = particle->Eta();

  if (eta <= fMinEta + PRECISION)
    return kFALSE;
  if (eta >= fMaxEta - PRECISION)
    return kFALSE;
  return kTRUE;
}






/***************************************************************************/ /**
                                                                               * Function to select tracks with required quality.
                                                                               ******************************************************************************/
Bool_t AliAnalysisTaskOmegaDielectron_AccEff::AcceptTrackQuality(AliVTrack *track) {
  if (fIsESD) {
    AliESDtrack *esdTrack = dynamic_cast<AliESDtrack *>(track);
    if (!fESDtrackCuts->AcceptTrack(esdTrack))
      return kFALSE;
  }
  return kTRUE;
}

/***************************************************************************/ /**
                                                                               * Function to obtain V0M centrality.
                                                                               ******************************************************************************/
Double_t AliAnalysisTaskOmegaDielectron_AccEff::GetCentrality(AliVEvent *event) {
  Double_t centrality = -1;
  AliMultSelection *multSelection =
      (AliMultSelection *)fEvent->FindListObject("MultSelection");
  if (!multSelection) {
    AliError("No MultSelection found!");
    return 999;
  }
  centrality = multSelection->GetMultiplicityPercentile("V0M");
  if (centrality > 100) {
    AliError("Centrality determination does not work proprely!");
    return 999;
  }
  return centrality;
}

/***************************************************************************/ /**
                                                                               * Function to initialize the ESD track cuts object.
                                                                               ******************************************************************************/
void AliAnalysisTaskOmegaDielectron_AccEff::InitESDTrackCuts() {

  if (fESDtrackCuts)
    delete fESDtrackCuts;
  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
  if (!fESDtrackCuts) {
    AliError("fESDtrackCuts not available");
    return;
  }

  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(4.0);
  fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.4);
  fESDtrackCuts->SetRequireITSRefit(kTRUE);
  fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                          AliESDtrackCuts::kAny);
  fESDtrackCuts->SetMaxChi2PerClusterITS(36.);
  fESDtrackCuts->SetDCAToVertex2D(kFALSE);
  fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);
  fESDtrackCuts->SetMaxDCAToVertexZ(2.0);
  fESDtrackCuts->SetMaxDCAToVertexXYPtDep(
      "0.0182+0.0350/pt^1.01"); // 7*(0.0026+0.0050/pt^1.01)
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36.); // tpcc cut
  fESDtrackCuts->SetCutGeoNcrNcl(3, 130, 1.5, 0.85,
                                 0.7); // Geometrical-Length Cut

  // cut-variations:
  if (fCutMode == 101) {
    fESDtrackCuts->SetMaxChi2PerClusterITS(25.);
  }
  if (fCutMode == 102) {
    fESDtrackCuts->SetMaxChi2PerClusterITS(49.);
  }

  if (fCutMode == 103) {
    fESDtrackCuts->SetMaxChi2PerClusterTPC(3.0);
  }
  if (fCutMode == 104) {
    fESDtrackCuts->SetMaxChi2PerClusterTPC(5.0);
  }

  if (fCutMode == 105) {
    fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
  }
  if (fCutMode == 106) {
    fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
  }

  if (fCutMode == 107) {
    fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.2);
  }
  if (fCutMode == 108) {
    fESDtrackCuts->SetMaxFractionSharedTPCClusters(1.0);
  }

  if (fCutMode == 109) {
    fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(25.);
  }
  if (fCutMode == 110) {
    fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(49.);
  }

  if (fCutMode == 111) {
    fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0104+0.0200/pt^1.01");
  }
  if (fCutMode == 112) {
    fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0260+0.0500/pt^1.01");
  }

  if (fCutMode == 113) {
    fESDtrackCuts->SetMaxDCAToVertexZ(1.0);
  }
  if (fCutMode == 114) {
    fESDtrackCuts->SetMaxDCAToVertexZ(5.0);
  }

  if (fCutMode == 115) {
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                            AliESDtrackCuts::kOff);
  }

  if (fCutMode == 116) {
    fESDtrackCuts->SetCutGeoNcrNcl(3, 120, 1.5, 0.85, 0.7);
  }
  if (fCutMode == 117) {
    fESDtrackCuts->SetCutGeoNcrNcl(3, 140, 1.5, 0.85, 0.7);
  }

  if (fCutMode == 118) {
    fESDtrackCuts->SetCutGeoNcrNcl(4, 130, 1.5, 0.85, 0.7);
  }
  if (fCutMode == 119) {
    fESDtrackCuts->SetCutGeoNcrNcl(2, 130, 1.5, 0.85, 0.7);
  }
}

/***************************************************************************/ /**
                                                                               * Function to get array of equidistant bin edges between lower and upper edge.
                                                                               ******************************************************************************/
void AliAnalysisTaskOmegaDielectron_AccEff::SetFixedBinEdges(Double_t *array,
                                             Double_t lowerEdge,
                                             Double_t upperEdge, Int_t nBins) {
  for (Int_t i = 0; i <= nBins; i++) {
    array[i] = lowerEdge + i * (upperEdge - lowerEdge) / nBins;
  }
}

/***************************************************************************/ /**
                                                                               * Function to create THnSparseF histogram with the specified axes.
                                                                               ******************************************************************************/
THnSparseF *AliAnalysisTaskOmegaDielectron_AccEff::CreateHistogram(string name,
                                                   vector<string> axes) {
  Int_t nAxes = axes.size();
  if (nAxes > MAX_HISTO_DIM)
    return nullptr;

  Int_t nBins[MAX_HISTO_DIM] = {0};
  Double_t lowerBounds[MAX_HISTO_DIM] = {0.0};
  Double_t upperBounds[MAX_HISTO_DIM] = {0.0};

  string title = name + " [";
  // first figure out number of bins and dimensions
  for (Int_t i = 0; i < nAxes; i++) {
    TArrayD *binEdges = GetBinEdges(axes[i]);
    nBins[i] = binEdges->GetSize() - 1;
    lowerBounds[i] = binEdges->GetAt(0);
    upperBounds[i] = binEdges->GetAt(binEdges->GetSize() - 1);
    title += axes[i];
    if (i < nAxes - 1)
      title += " : ";
    else
      title += "]";
  }
  // create histogram
  THnSparseF *histogram = new THnSparseF(name.c_str(), title.c_str(), nAxes,
                                         nBins, lowerBounds, upperBounds);

  // set histogram axes
  for (Int_t i = 0; i < nAxes; i++) {
    TArrayD *binEdges = GetBinEdges(axes[i]);
    histogram->SetBinEdges(i, binEdges->GetArray());
    histogram->GetAxis(i)->SetTitle(GetAxisTitle(axes[i]).c_str());
  }
  histogram->Sumw2();
  return histogram;
}

/***************************************************************************/ /**
                                                                               * Function to obtain the correct binning for the respective axis.
                                                                               ******************************************************************************/
TArrayD *AliAnalysisTaskOmegaDielectron_AccEff::GetBinEdges(string &axisName) {
  if (axisName.find("sigmapt") != string::npos)
    return fBinsPtReso;
  else if (axisName.find("deltapt") != string::npos)
    return fBinsPtReso;
  else if (axisName.find("pt") != string::npos)
    return fBinsPt;
  else if (axisName.find("eta") != string::npos)
    return fBinsEta;
  else if (axisName.find("y") != string::npos)
    return fBinsy;
  else if (axisName.find("mult") != string::npos)
    return fBinsMult;
  else if (axisName.find("cent") != string::npos)
    return fBinsCent;
  else if (axisName.find("zv") != string::npos)
    return fBinsZv;
  else
    return nullptr;
}

/***************************************************************************/ /**
                                                                               * Function to get the correct title for each histogram axis.
                                                                               ******************************************************************************/
string AliAnalysisTaskOmegaDielectron_AccEff::GetAxisTitle(string &axisName) {
  if (axisName == "pt")
    return "#it{p}_{T} (GeV/#it{c})";
  else if (axisName == "deltapt")
    return "#Delta(#it{p}_{T}) / #it{p}^{ meas}_{T}";
  else if (axisName == "mult")
    return "Multiplicity";
  else if (axisName == "cent")
    return "Centrality (%)";
  else if (axisName == "eta_meas")
    return "#eta^{ meas}";
  else if (axisName == "eta_true")
    return "#eta^{ true}";
  else if (axisName == "pt_meas")
    return "#it{p}^{ meas}_{T} (GeV/#it{c})";
  else if (axisName == "pt_true")
    return "#it{p}^{ true}_{T} (GeV/#it{c})";
  else if (axisName == "mult_meas")
    return "#it{N}^{ meas}_{ch}";
  else if (axisName == "mult_true")
    return "#it{N}^{ true}_{ch}";
  else if (axisName == "sigmapt")
    return "#sigma(#it{p}^{ meas}_{T}) / #it{p}^{ meas}_{T}";
  else
    return "dummyTitle";
}

/***************************************************************************/ /**
                                                                               * Function to fill a histogram.
                                                                               ******************************************************************************/
void AliAnalysisTaskOmegaDielectron_AccEff::FillHisto(THnSparseF *histo,
                                      array<Double_t, MAX_HISTO_DIM> values) {
  histo->Fill(values.data());
}

/***************************************************************************/ /**
                                                                               * Function to set variable binning for multiplicity.
                                                                               ******************************************************************************/
void AliAnalysisTaskOmegaDielectron_AccEff::SetBinsMult(vector<Int_t> multSteps,
                                        vector<Int_t> multBinWidth) {
  if (multSteps.size() != multBinWidth.size()) {
    AliError("SetBinsMult:: Vectors need to have same size!");
    return;
  }

  Int_t nMultSteps = multSteps.size();
  Int_t nBinsMult = 1; // for mult=0 bin
  for (Int_t multBins : multSteps)
    nBinsMult += multBins;
  Double_t *multBinEdges = new Double_t[nBinsMult + 1]; // edges need one more

  multBinEdges[0] = -0.5;
  multBinEdges[1] = 0.5;
  Int_t startBin = 1;
  Int_t endBin = 1;
  for (Int_t multStep = 0; multStep < nMultSteps; multStep++) {
    endBin += multSteps[multStep];
    for (Int_t multBin = startBin; multBin < endBin; multBin++) {
      multBinEdges[multBin + 1] =
          multBinEdges[multBin] + multBinWidth[multStep];
    }
    startBin = endBin;
  }
  SetBinsMult(nBinsMult, multBinEdges);
}

/***************************************************************************/ /**
                                                                               * Function to set maxMult single multiplicity steps.
                                                                               ******************************************************************************/
void AliAnalysisTaskOmegaDielectron_AccEff::SetBinsMult(Int_t maxMult) {
  return SetBinsMult({maxMult}, {1});
}

/***************************************************************************/ /**
                                                                               * Function to add this task to a train.
                                                                               ******************************************************************************/
AliAnalysisTaskOmegaDielectron_AccEff *
AliAnalysisTaskOmegaDielectron_AccEff::AddTaskMultDepSpec(TString controlstring,
                                          Int_t cutModeLow, Int_t cutModeHigh) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMultDepSpec", "No analysis manager found.");
    return nullptr;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMultDepSpec", "No input event handler found.");
    return nullptr;
  }
  TString type =
      mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Bool_t isMC =
      (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() !=
       nullptr);

  // Default cut settings:
  UInt_t triggerMask = AliVEvent::kMB | AliVEvent::kINT7;

  Double_t cutVertexZ = 10.0;

  Double_t cutCentLow = 0.0;
  Double_t cutCentHigh = 100.0;

  Double_t cutPtLow = 0.15;
  Double_t cutPtHigh = 50.0;

  Double_t cutEtaLow = -0.8;
  Double_t cutEtaHigh = 0.8;

  Int_t maxMult = 100;
  string colsys = "pp";

  Bool_t useCent = kFALSE;
  if (controlstring.Contains("useCent"))
    useCent = kTRUE;
  Double_t centBinEdges[9] = {0., 5., 10., 20., 40., 60., 80., 90., 100.};

  // colison system specific settings
  if (controlstring.Contains("pp")) {
    colsys = "pp";
    maxMult = 100;
  } else if (controlstring.Contains("pPb")) {
    colsys = "pPb";
    maxMult = 200;
  } else if (controlstring.Contains("XeXe")) {
    colsys = "XeXe";
    maxMult = 3500;
  } else if (controlstring.Contains("PbPb")) {
    colsys = "PbPb";
    maxMult = 4500;
  }

  AliAnalysisTaskOmegaDielectron_AccEff *returnTask = nullptr;

  char taskName[100] = "";

  for (Int_t cutMode = cutModeLow; cutMode <= cutModeHigh; cutMode++) {
    sprintf(taskName, "multDepSpec_%s_cutMode_%d", colsys.c_str(), cutMode);

    AliAnalysisTaskOmegaDielectron_AccEff *task = new AliAnalysisTaskOmegaDielectron_AccEff(taskName);
    if (cutMode == cutModeLow)
      returnTask = task; // return one of the tasks

    task->SetCutMode(cutMode);
    task->SetTriggerMask(triggerMask);

    task->SetIsMC(isMC);
    if (type.Contains("ESD"))
      task->SetUseESD();
    else
      task->SetUseAOD();

    task->SetBinsMult(maxMult);

    if (useCent) {
      task->SetMinCent(cutCentLow);
      task->SetMaxCent(cutCentHigh);
      task->SetBinsCent(8, centBinEdges);
    }

    // kinematic cuts:
    task->SetMinEta(cutEtaLow);
    task->SetMaxEta(cutEtaHigh);
    task->SetMinPt(cutPtLow);
    task->SetMaxPt(cutPtHigh);
    task->SetMaxZv(cutVertexZ);

    // hang task in train
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(
        task, 1,
        mgr->CreateContainer(taskName, TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             "AnalysisResults.root"));
  }
  return returnTask;
}
