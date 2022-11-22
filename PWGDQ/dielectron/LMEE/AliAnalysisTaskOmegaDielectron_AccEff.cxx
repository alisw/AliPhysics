#include "AliAnalysisTaskOmegaDielectron_AccEff.h"
#include <iostream>

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskOmegaDielectron_AccEff);
/// \endcond

/***************************************************************************/ /**
* ROOT I/O Constructor.
******************************************************************************/
AliAnalysisTaskOmegaDielectron_AccEff::AliAnalysisTaskOmegaDielectron_AccEff()
: AliAnalysisTaskSE(),
// General member variables
fOutputList(nullptr), fEvent(nullptr), fMCEvent(nullptr), fEventCuts(nullptr),
fFilter_TrackCuts(), fFilter_PID(),
fPIDResponse(nullptr),
fMinEta(-0.8), fMaxEta(0.8), fMinPt(0.2), fMaxPt(10.0),
//pdg codes:
felectron_pdg(11),fpositron_pdg(-11),
fmother_pdg(223),
// fmother_pdg(333),
//storage Vectors
v_elec_true_omega(0x0),
v_posi_true_omega(0x0),
v_elec_true_omega_MCPart(0x0),
v_posi_true_omega_MCPart(0x0),
v_elec_motherID_true_omega(0),
v_posi_motherID_true_omega(0),
// Event-Histograms
fHistEventSelection(nullptr),
fHist_MC_Omegas_Rapidity(nullptr),
fHist_MC_Omegas_gen(nullptr),
fHist_MC_Omegas_gen_DaughtersinAcc(nullptr),
fHist_MC_elec_gen(nullptr),
fHist_MC_posi_gen(nullptr),
fHist_MC_elec_gen_inAcc(nullptr),
fHist_MC_posi_gen_inAcc(nullptr),
fHist_elec_rec_inAcc(nullptr),
fHist_elec_rec_inAcc_Track(nullptr),
fHist_elec_rec_inAcc_Track_PID(nullptr),
fHist_posi_rec_inAcc(nullptr),
fHist_posi_rec_inAcc_Track(nullptr),
fHist_posi_rec_inAcc_Track_PID(nullptr),
fHist_MC_Omegas_withoutCuts(nullptr),
fHist_MC_Omegas_TrackCuts(nullptr),
fHist_MC_Omegas_TrackPID(nullptr),
fHist_Rec_Omegas_withoutCuts(nullptr),
fHist_Rec_Omegas_TrackCuts(nullptr),
fHist_Rec_Omegas_TrackPID(nullptr)
{
  // ROOT IO constructor, don't allocate memory here!
}

/***************************************************************************/ /**
* Constructor.
******************************************************************************/
AliAnalysisTaskOmegaDielectron_AccEff::AliAnalysisTaskOmegaDielectron_AccEff(const char *name)
: AliAnalysisTaskSE(name),
// General member variables
fOutputList(nullptr), fEvent(nullptr), fMCEvent(nullptr), fEventCuts(nullptr),
fFilter_TrackCuts(), fFilter_PID(),
fPIDResponse(nullptr),
fMinEta(-0.8), fMaxEta(0.8), fMinPt(0.2), fMaxPt(10.0),
//pdg codes:
felectron_pdg(11),fpositron_pdg(-11),
fmother_pdg(223),
// fmother_pdg(333),
//storage Vectors
v_elec_true_omega(0x0),
v_posi_true_omega(0x0),
v_elec_true_omega_MCPart(0x0),
v_posi_true_omega_MCPart(0x0),
v_elec_motherID_true_omega(0),
v_posi_motherID_true_omega(0),
// Event-Histograms
fHistEventSelection(nullptr),
fHist_MC_Omegas_Rapidity(nullptr),
fHist_MC_Omegas_gen(nullptr),
fHist_MC_Omegas_gen_DaughtersinAcc(nullptr),
fHist_MC_elec_gen(nullptr),
fHist_MC_posi_gen(nullptr),
fHist_MC_elec_gen_inAcc(nullptr),
fHist_MC_posi_gen_inAcc(nullptr),
fHist_elec_rec_inAcc(nullptr),
fHist_elec_rec_inAcc_Track(nullptr),
fHist_elec_rec_inAcc_Track_PID(nullptr),
fHist_posi_rec_inAcc(nullptr),
fHist_posi_rec_inAcc_Track(nullptr),
fHist_posi_rec_inAcc_Track_PID(nullptr),
fHist_MC_Omegas_withoutCuts(nullptr),
fHist_MC_Omegas_TrackCuts(nullptr),
fHist_MC_Omegas_TrackPID(nullptr),
fHist_Rec_Omegas_withoutCuts(nullptr),
fHist_Rec_Omegas_TrackCuts(nullptr),
fHist_Rec_Omegas_TrackPID(nullptr)
{
  DefineOutput(1, TList::Class());
}

/***************************************************************************/
/*** Function executed once before the event loop. Create histograms here.**
****************************************************************************/
void AliAnalysisTaskOmegaDielectron_AccEff::UserCreateOutputObjects() {

  OpenFile(1, "recreate");
  fOutputList = new TList();
  fOutputList->SetOwner();

  // Control histogram to check the effect of event cuts
  fHistEventSelection = new TH1F("fHistEventSelection", "fHistEventSelection [all : selected]", 2, 0., 3.);
  fHistEventSelection->GetYaxis()->SetTitle("#it{N}_{events}");
  fHistEventSelection->GetXaxis()->SetBinLabel(1, "in");
  fHistEventSelection->GetXaxis()->SetBinLabel(2, "out");
  fOutputList->Add(fHistEventSelection);

  // Histogram to check Rapidity distribution
  fHist_MC_Omegas_Rapidity = new TH1F("fHist_MC_Omegas_Rapidity", "Omega - Rapidity;y;#it{N}_{events}", 200, -10, 10);
  fOutputList->Add(fHist_MC_Omegas_Rapidity);


  /// --- --- --- MC histograms - generated/accepted Omegas --- --- --- ///
  fHist_MC_Omegas_gen = new TH3D("fHist_MC_Omegas_gen", "fHist_MC_Omegas_gen;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});y",  100, 0, 10, 100, 0, 10, 20, -1, 1);
  fOutputList->Add(fHist_MC_Omegas_gen);
  fHist_MC_Omegas_gen_DaughtersinAcc = new TH3D("fHist_MC_Omegas_gen_DaughtersinAcc", "fHist_MC_Omegas_gen_DaughtersinAcc;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});y", 100, 0, 10, 100, 0, 10, 20, -1, 1);
  fOutputList->Add(fHist_MC_Omegas_gen_DaughtersinAcc);

  // --- single electron histos: pt - eta - phi
  fHist_MC_elec_gen = new TH3D("fHist_MC_elec_gen", "fHist_MC_elec_gen;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", 100, 0, 10, 20, -1, 1, 64, 0, 6.4); // 6.4 == 2*PI()
  fOutputList->Add(fHist_MC_elec_gen);
  fHist_MC_elec_gen_inAcc = new TH3D("fHist_MC_elec_gen_inAcc", "fHist_MC_elec_gen_inAcc;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", 100, 0, 10, 20, -1, 1, 64, 0, 6.4); // 6.4 == 2*PI()
  fOutputList->Add(fHist_MC_elec_gen_inAcc);
  fHist_elec_rec_inAcc = new TH3D("fHist_elec_rec_inAcc", "fHist_elec_rec_inAcc;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", 100, 0, 10, 20, -1, 1, 64, 0, 6.4); // 6.4 == 2*PI()
  fOutputList->Add(fHist_elec_rec_inAcc);
  fHist_elec_rec_inAcc_Track = new TH3D("fHist_elec_rec_inAcc_Track", "fHist_elec_rec_inAcc_Track;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", 100, 0, 10, 20, -1, 1, 64, 0, 6.4); // 6.4 == 2*PI()
  fOutputList->Add(fHist_elec_rec_inAcc_Track);
  fHist_elec_rec_inAcc_Track_PID = new TH3D("fHist_elec_rec_inAcc_Track_PID", "fHist_elec_rec_inAcc_Track_PID;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", 100, 0, 10, 20, -1, 1, 64, 0, 6.4); // 6.4 == 2*PI()
  fOutputList->Add(fHist_elec_rec_inAcc_Track_PID);

  // --- single positron histos: pt - eta - phi
  fHist_MC_posi_gen = new TH3D("fHist_MC_posi_gen", "fHist_MC_posi_gen;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", 100, 0, 10, 20, -1, 1, 64, 0, 6.4); // 6.4 == 2*PI()
  fOutputList->Add(fHist_MC_posi_gen);
  fHist_MC_posi_gen_inAcc = new TH3D("fHist_MC_posi_gen_inAcc", "fHist_MC_posi_gen_inAcc;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", 100, 0, 10, 20, -1, 1, 64, 0, 6.4); // 6.4 == 2*PI()
  fOutputList->Add(fHist_MC_posi_gen_inAcc);
  fHist_posi_rec_inAcc = new TH3D("fHist_posi_rec_inAcc", "fHist_posi_rec_inAcc;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", 100, 0, 10, 20, -1, 1, 64, 0, 6.4); // 6.4 == 2*PI()
  fOutputList->Add(fHist_posi_rec_inAcc);
  fHist_posi_rec_inAcc_Track = new TH3D("fHist_posi_rec_inAcc_Track", "fHist_posi_rec_inAcc_Track;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", 100, 0, 10, 20, -1, 1, 64, 0, 6.4); // 6.4 == 2*PI()
  fOutputList->Add(fHist_posi_rec_inAcc_Track);
  fHist_posi_rec_inAcc_Track_PID = new TH3D("fHist_posi_rec_inAcc_Track_PID", "fHist_posi_rec_inAcc_Track_PID;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", 100, 0, 10, 20, -1, 1, 64, 0, 6.4); // 6.4 == 2*PI()
  fOutputList->Add(fHist_posi_rec_inAcc_Track_PID);

  /// --- --- --- histograms -> reconstructed Omega --- --- --- ///
  // m - pt - y
  fHist_MC_Omegas_withoutCuts = new TH3D("fHist_MC_Omegas_withoutCuts", "fHist_MC_Omegas_withoutCuts;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T,ee} (GeV/#it{c});y", 100, 0, 10, 100, 0, 10, 20, -1, 1);
  fOutputList->Add(fHist_MC_Omegas_withoutCuts);
  fHist_MC_Omegas_TrackCuts = new TH3D("fHist_MC_Omegas_TrackCuts", "fHist_MC_Omegas_TrackCuts;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T,ee} (GeV/#it{c});y", 100, 0, 10, 100, 0, 10, 20, -1, 1);
  fOutputList->Add(fHist_MC_Omegas_TrackCuts);
  fHist_MC_Omegas_TrackPID = new TH3D("fHist_MC_Omegas_TrackPID", "fHist_MC_Omegas_TrackPID;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T,ee} (GeV/#it{c});y",  100, 0, 10, 100, 0, 10, 20, -1, 1);
  fOutputList->Add(fHist_MC_Omegas_TrackPID);
  fHist_Rec_Omegas_withoutCuts = new TH3D("fHist_Rec_Omegas_withoutCuts", "fHist_Rec_Omegas_withoutCuts;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T,ee} (GeV/#it{c});y", 100, 0, 10, 100, 0, 10, 20, -1, 1);
  fOutputList->Add(fHist_Rec_Omegas_withoutCuts);
  fHist_Rec_Omegas_TrackCuts = new TH3D("fHist_Rec_Omegas_TrackCuts", "fHist_Rec_Omegas_TrackCuts;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T,ee} (GeV/#it{c});y", 100, 0, 10, 100, 0, 10, 20, -1, 1);
  fOutputList->Add(fHist_Rec_Omegas_TrackCuts);
  fHist_Rec_Omegas_TrackPID = new TH3D("fHist_Rec_Omegas_TrackPID", "fHist_Rec_Omegas_TrackPID;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T,ee} (GeV/#it{c});y", 100, 0, 10, 100, 0, 10, 20, -1, 1);
  fOutputList->Add(fHist_Rec_Omegas_TrackPID);

  /// --- --- --- EventCuts (for AOD/ ESD) --- --- --- ///
  fEventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  fEventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD);
  fEventCuts->SetRequireVertex();
  fEventCuts->SetVertexZ(-10.,10.);
  fEventCuts->SetMinVtxContributors(1);

  InitCuts();

  PostData(1, fOutputList);
}

/// Destructor
AliAnalysisTaskOmegaDielectron_AccEff::~AliAnalysisTaskOmegaDielectron_AccEff() {

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

  AliInputEventHandler *eventHandler = nullptr;
  eventHandler = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  if (!fPIDResponse) SetPIDResponse( eventHandler->GetPIDResponse() );
  AliDielectronVarManager::SetPIDResponse(fPIDResponse);

  fMCEvent = MCEvent();
  AliDielectronVarManager::SetEvent(fMCEvent);
  if (!fMCEvent) {
    AliError("fMCEvent not available\n");
    return;
  }


  ///-------------------- Fill Events in Histogramm :: All;Selected--------------------///
  fHistEventSelection->Fill(1.0); // all events
  if (!fEventCuts->IsSelected(fEvent)){ // https://github.com/alisw/AliPhysics/blob/master/OADB/AliEventCuts.cxx
    return;
  }
  fHistEventSelection->Fill(2.0); // selected events


  ///-------------------- Loop over Generated MC Particles ---OMEGAS --------------------///
  for(Int_t iGenPart = 0; iGenPart < fMCEvent->GetNumberOfTracks();iGenPart++) {
    AliMCParticle* mcGenParticle  = (AliMCParticle*)fMCEvent->GetTrack(iGenPart);
    if(!mcGenParticle){
      AliError("mcGenParticle not available\n");
      continue;
    }

    Int_t pdgcode_gen = mcGenParticle->PdgCode();

    if(pdgcode_gen == fmother_pdg){  // Check if MC Particle is an omega
      fHist_MC_Omegas_Rapidity->Fill(mcGenParticle->Y());

      if( CheckDielectronDecay(mcGenParticle, kFALSE) ){ // Check if :: w->e+e- and omega within y<+-0.8
        fHist_MC_Omegas_gen->Fill(mcGenParticle->M(), mcGenParticle->Pt(), mcGenParticle->Y());
      }

      if( CheckDielectronDecay(mcGenParticle, kTRUE) ){ // Check if :: w->e+e-and omega within y<+-0.8 --> with /n_e/<0.8 & 0.2<p_T,e < 10GeV/c
        fHist_MC_Omegas_gen_DaughtersinAcc->Fill(mcGenParticle->M(), mcGenParticle->Pt(), mcGenParticle->Y());
      }
    }

    // for the single electron efficiency:
    if (mcGenParticle->IsPhysicalPrimary() ){
      if(pdgcode_gen == felectron_pdg){  // Check if MC Particle is an electron
        fHist_MC_elec_gen->Fill(mcGenParticle->Pt(), mcGenParticle->Eta(), mcGenParticle->Phi());
        if (AcceptKinematics(mcGenParticle)){
          fHist_MC_elec_gen_inAcc->Fill(mcGenParticle->Pt(), mcGenParticle->Eta(), mcGenParticle->Phi());
        }
      }
      if(pdgcode_gen == fpositron_pdg){  // Check if MC Particle is an positron
        fHist_MC_posi_gen->Fill(mcGenParticle->Pt(), mcGenParticle->Eta(), mcGenParticle->Phi());
        if (AcceptKinematics(mcGenParticle)){
          fHist_MC_posi_gen_inAcc->Fill(mcGenParticle->Pt(), mcGenParticle->Eta(), mcGenParticle->Phi());
        }
      }
    }
    // else{continue;}
  }


  // SelectionMask for the track and PID Cuts
  UInt_t selectedMask_Track =( 1 << fFilter_TrackCuts->GetCuts()->GetEntries())-1; // cutting logic taken from AliDielectron::FillTrackArrays()   // apply track cuts
  UInt_t selectedMask_PID =( 1 << fFilter_PID->GetCuts()->GetEntries())-1;

  ///-------------------- Loop over AOD Tracks - electrons// positrons --------------------///
  AliAODTrack *track = nullptr;
  // Create a vector containing the electron/positron tracks
  v_elec_true_omega.clear();
  v_elec_true_omega_MCPart.clear();
  v_elec_motherID_true_omega.clear();
  v_posi_true_omega.clear();
  v_posi_true_omega_MCPart.clear();
  v_posi_motherID_true_omega.clear();

  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); iTrack++) {
    track = static_cast<AliAODTrack*>(fEvent->GetTrack(iTrack));

    if (!track || !track->TestFilterBit(16)){ // SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); // 1<<4 // -> FilterBit 16
      // filter bit 128 denotes TPC-only tracks, use only them for the analysis
      continue;
    }

    if (!AcceptKinematics(track)){ // Check if electron/positron  /n_e/<0.8 & 0.2<p_T,e < 10GeV/c
      continue;
    }

    ///--- Find original particle in MC-Stack ---///
    Int_t mcLabel =  TMath::Abs(track->GetLabel()); // negative label means bad quality track
    AliMCParticle *mcParticle = (AliMCParticle *)fMCEvent->GetTrack(mcLabel);
    if (!mcParticle) {
      AliError("mcParticle not available\n");
      continue;
    }
    if ( !mcParticle->IsPhysicalPrimary() ){
      continue;
    }

    Int_t pdgcode = mcParticle->PdgCode();
    AliMCParticle* mcParticle_mother  = (AliMCParticle*)fMCEvent->GetTrack( mcParticle->GetMother() );
    Int_t mother_id = ( mcParticle_mother->GetLabel() );
    Int_t pdgCode_mother = ( mcParticle_mother->PdgCode() );

    UInt_t cutMask = fFilter_TrackCuts->IsSelected(track);
    UInt_t cutMask_PID = fFilter_PID->IsSelected(track);
    if(pdgcode == felectron_pdg){
      fHist_elec_rec_inAcc->Fill(track->Pt(), track->Eta(), track->Phi());
      if (cutMask == selectedMask_Track) {
        fHist_elec_rec_inAcc_Track->Fill(track->Pt(), track->Eta(), track->Phi());
        if (cutMask_PID == selectedMask_PID) {
          fHist_elec_rec_inAcc_Track_PID->Fill(track->Pt(), track->Eta(), track->Phi());
        }
      }
    }
    if(pdgcode == fpositron_pdg){
      fHist_posi_rec_inAcc->Fill(track->Pt(), track->Eta(), track->Phi());
      if (cutMask == selectedMask_Track) {
        fHist_posi_rec_inAcc_Track->Fill(track->Pt(), track->Eta(), track->Phi());
        if (cutMask_PID == selectedMask_PID) {
          fHist_posi_rec_inAcc_Track_PID->Fill(track->Pt(), track->Eta(), track->Phi());
        }
      }
    }

    ////***************ELECTRONS*****************////
    if(pdgCode_mother == fmother_pdg){                        // Check if mother of MC Particle (Track) is an omega
      if(pdgcode == felectron_pdg){                           // Check if MC Particle (Track) is an electron
        if(CheckDielectronDecay(mcParticle_mother, kFALSE)){          // Check if omega has a 2-body dielectron decay
          v_elec_true_omega.push_back(track);
          v_elec_true_omega_MCPart.push_back(mcParticle_mother);
          v_elec_motherID_true_omega.push_back(mother_id);
        }
      }
      ////***************POSITRONS*****************////
      if(pdgcode == fpositron_pdg){
        if(CheckDielectronDecay(mcParticle_mother, kFALSE)){
          v_posi_true_omega.push_back(track);
          v_posi_true_omega_MCPart.push_back(mcParticle_mother);
          v_posi_motherID_true_omega.push_back(mother_id);
        }
      }

      else{
        continue;
      }
    }

  }



  ///-------------------- Pair electron + positron of the same Omega-mother  --------------------///
  ///--------------------  Go over the selected electron/positron Tracks    --------------------///

  if (v_posi_true_omega.size() == 0) return;
  if (v_elec_true_omega.size() == 0) return;
  else{
    /// ~~~~~~~ electron vectors ~~~~~~~ ///
    for(int i=0; i< (Int_t)v_elec_true_omega.size(); i++){
      AliAODTrack * ele_from_same_mother_track= (AliAODTrack *)v_elec_true_omega.at(i);
      AliMCParticle *mcParticle_same_mother = (AliMCParticle *)v_elec_true_omega_MCPart.at(i); // Omega mother!
      int mother_electron= (Int_t)v_elec_motherID_true_omega.at(i);

      /// ~~~~~~~ positron vectors ~~~~~~~ ///
      for(int j=0; j< (Int_t)v_posi_true_omega.size(); j++){
        AliAODTrack * pos_from_same_mother_track= (AliAODTrack *)v_posi_true_omega.at(j);
        // AliMCParticle *mcParticle_p_same_mother = (AliMCParticle *)v_posi_true_omega_MCPart.at(i); // Omega mother!
        int mother_positron= (Int_t)v_posi_motherID_true_omega.at(j);

        if(mother_electron==mother_positron){
          Double_t electron_mass = (AliPID::ParticleMass(AliPID::kElectron) );
          //Tracks:
          TLorentzVector Lvec1;
          Lvec1.SetPtEtaPhiM(ele_from_same_mother_track->Pt(), ele_from_same_mother_track->Eta(), ele_from_same_mother_track->Phi(), electron_mass);
          TLorentzVector Lvec2;
          Lvec2.SetPtEtaPhiM(pos_from_same_mother_track->Pt(), pos_from_same_mother_track->Eta(), pos_from_same_mother_track->Phi(), electron_mass);

          TLorentzVector LvecM = Lvec1 + Lvec2;
          Double_t mass = LvecM.M();
          Double_t pairpt = LvecM.Pt();
          Double_t pairy = LvecM.Rapidity();
          Double_t weight = 1;
          // std::cout << "  mass: " << mass << "  pairpt: " << pairpt << '\n';
          fHist_Rec_Omegas_withoutCuts->Fill(mass, pairpt, pairy);
          fHist_MC_Omegas_withoutCuts->Fill(mcParticle_same_mother->M(), mcParticle_same_mother->Pt(), mcParticle_same_mother->Y());

          /// ~~~~~~ ***** pass Track cuts *****  ~~~~~~ ///
          UInt_t cutMask_electron = fFilter_TrackCuts->IsSelected(ele_from_same_mother_track);
          UInt_t cutMask_positron = fFilter_TrackCuts->IsSelected(pos_from_same_mother_track);
          if (cutMask_electron == selectedMask_Track && cutMask_positron == selectedMask_Track) {
            fHist_Rec_Omegas_TrackCuts->Fill(mass, pairpt, pairy);
            // omega histo
            // cout << mcParticle_same_mother->M()<< '\n';
            fHist_MC_Omegas_TrackCuts->Fill(mcParticle_same_mother->M(), mcParticle_same_mother->Pt(), mcParticle_same_mother->Y());


            /// ~~~~~~ ***** pass PID cuts *****  ~~~~~~ ///
            UInt_t cutMask_PID_electron = fFilter_PID->IsSelected(ele_from_same_mother_track);
            UInt_t cutMask_PID_positron = fFilter_PID->IsSelected(pos_from_same_mother_track);
            if (cutMask_PID_electron == selectedMask_PID && cutMask_PID_positron == selectedMask_PID) {
              fHist_Rec_Omegas_TrackPID->Fill(mass, pairpt, pairy);
              // omega histo
              fHist_MC_Omegas_TrackPID->Fill(mcParticle_same_mother->M(), mcParticle_same_mother->Pt(), mcParticle_same_mother->Y());
            }
          }

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

// /***************************************************************************/ /**
// * Function to Check if the (mother)particle decays via dielectrons
// * with Bool to check if the daughters are within the kinematic acceptance
// ******************************************************************************/
Bool_t AliAnalysisTaskOmegaDielectron_AccEff::CheckDielectronDecay(AliMCParticle *motherparticle, Bool_t checkacc) {
  // Check if the decay of "particle" has only 2 daughters::
  Int_t number_daughters = motherparticle->GetNDaughters();
  if(number_daughters != 2){
    return false;
  }
  // Check the Rapidity range of motherparticle:
  double rapidi = 0.8;
  if ( TMath::Abs( motherparticle->Y() ) >= rapidi ){
    return false;
  } // omega :: /y/<0.8

  // Get the daughters:
  Int_t LabelFirstDaughter = motherparticle->GetDaughterFirst();
  Int_t LabelLastDaughter = motherparticle->GetDaughterLast();
  AliMCParticle *daughter1 = (AliMCParticle*)fMCEvent->GetTrack(LabelFirstDaughter);
  int pdgCode_daughter1 = daughter1->PdgCode();
  AliMCParticle *daughter2 = (AliMCParticle*)fMCEvent->GetTrack(LabelLastDaughter);
  int pdgCode_daughter2 = daughter2->PdgCode();

  //checkacc==true-> check if daughters in Acceptence
  if (checkacc){
    if (!(AcceptKinematics(daughter1) && AcceptKinematics(daughter2))){
      return false;
    }
  }

  /// ***** OPTION 1 ***** ///
  if ( TMath::Abs(pdgCode_daughter1) == TMath::Abs(pdgCode_daughter2) ){    // check if daughters have the same absolute pdgcode
    if (pdgCode_daughter1 * pdgCode_daughter2 == (-121)){                   // the multipication electron * positron pdg = 11*(-11) = -121 --> 121 divides only by 121, 11, 1 & there is no particle with PDGcode 121...
      return true;
    }
    else
    return false;
  }

  else {
    return false;
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

// /***************************************************************************/ /**
// * Function to initialize the cuts.
// ******************************************************************************/
void AliAnalysisTaskOmegaDielectron_AccEff::InitCuts() {
  std::cout << "InitCuts()" <<std::endl;

  // TrackCuts (in one)
  if (fFilter_TrackCuts)
  delete fFilter_TrackCuts;
  fFilter_TrackCuts = new AliAnalysisFilter("Analysis Filter _ fFilter_TrackCuts");

  if (!fFilter_TrackCuts) {
    AliError("fFilter_TrackCuts not available");
    return;
  }
  fFilter_TrackCuts->AddCuts(SetupTrackCuts());


  // PIDCuts
  if (fFilter_PID)
  delete fFilter_PID;
  fFilter_PID = new AliAnalysisFilter("Analysis Filter _ fFilter_PID");

  if (!fFilter_PID) {
    AliError("fFilter_PID not available");
    return;
  }
  fFilter_PID->AddCuts(SetPIDcuts());
}

// /***************************************************************************/ /**
// * Track cuts.
// ******************************************************************************/
AliAnalysisCuts *AliAnalysisTaskOmegaDielectron_AccEff::SetupTrackCuts(){


  std::cout << "SetupTrackCuts()" <<std::endl;
  // AliAnalysisCuts* trackCuts=0x0;

  AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
  trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); // 1<<4
  trackCutsDiel->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);
  trackCutsDiel->SetRequireITSRefit(kTRUE);
  trackCutsDiel->SetRequireTPCRefit(kTRUE);


  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
  trackCutsAOD->AddCut(AliDielectronVarManager::kEta, -0.8,   0.8);
  trackCutsAOD->AddCut(AliDielectronVarManager::kPt, 0.2,   10.);
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0); // needs in UserExec() -->  AliDielectronVarManager::SetEvent(fMCEvent);
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0); // needs in UserExec() -->  AliDielectronVarManager::SetEvent(fMCEvent);

  // Number of ITS clusters
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 999.0);

  // Chi2 per ITS cluster
  // trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    -999.,   4.5);
  trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    -999.,   5.5);

  // Number of shared ITS clusters
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS, 1.0, 6.0, kTRUE);

  // Min TPCcls
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 999.0);

  //Min TPC cross rows
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 999.0);

  // Min ratio TPC cross row over findable
  // trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 999.);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 999.);

  // Shared TPC clusters
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSFracTPC,     -999., 0.4);

  // chi2 per TPC cls
  trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    -999.,   4.0);


  AliDielectronCutGroup* trackCuts = new AliDielectronCutGroup("Trackcuts","Trackcuts",AliDielectronCutGroup::kCompAND);
  trackCuts->AddCut(trackCutsAOD);
  trackCuts->AddCut(trackCutsDiel);


  return trackCuts;

}


// /***************************************************************************/ /**
// * PID cuts.
// ******************************************************************************/
AliAnalysisCuts *AliAnalysisTaskOmegaDielectron_AccEff::SetPIDcuts(){

  // AliAnalysisCuts *fancyCut = 0x0;

  AliDielectronPID *mastermind_TPC = new AliDielectronPID("mastermind_TPC","mastermind_TPC");
  AliDielectronPID *mastermind_TOF = new AliDielectronPID("mastermind_TOF","mastermind_TOF");

  //TPC electrons: includes electrons and exclude all possible other contributions using the TPC
  mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3.5,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  // mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kMuon,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);

  // //TOF electrons: includes all electrons, exlcludes Pions using the TPC
  mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3., 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 3.5 , 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);

  // Combine
  AliDielectronCutGroup* mastermind_cg = new AliDielectronCutGroup("mastermind_cg","mastermind_cg",AliDielectronCutGroup::kCompOR);
  mastermind_cg->AddCut(mastermind_TPC);
  mastermind_cg->AddCut(mastermind_TOF);

  return mastermind_cg;

}
