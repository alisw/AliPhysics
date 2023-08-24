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
fOutputList(nullptr), fCutVariList(0x0), fEvent(nullptr), fMCEvent(nullptr), fEventCuts(nullptr), fEventCuts_VertexZ(nullptr),
fFilter_TrackCuts(), fFilter_PID(),
fPIDResponse(nullptr),
fMinEta(-0.8), fMaxEta(0.8), fMinPt(0.2), fMaxPt(10.0),
fsysUnc(kFALSE),
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
fHist_MC_Omegas_gen(),
fHist_MC_Omegas_gen_DaughtersinAcc(),
fHist_Rec_Omegas_TrackPID(),

fHistVertex(),
fHist_MC_Omegas_Rapidity(),
fHist_MC_elec_gen(),
fHist_MC_posi_gen(),
fHist_MC_elec_gen_inAcc(),
fHist_MC_posi_gen_inAcc(),
fHist_elec_rec_inAcc(),
fHist_elec_rec_inAcc_Track(),
fHist_elec_rec_inAcc_Track_PID(),
fHist_posi_rec_inAcc(),
fHist_posi_rec_inAcc_Track(),
fHist_posi_rec_inAcc_Track_PID(),
fHist_MC_Omegas_withoutCuts(),
fHist_MC_Omegas_TrackCuts(),
fHist_MC_Omegas_TrackPID(),
fHist_Rec_Omegas_withoutCuts(),
fHist_Rec_Omegas_TrackCuts()
{
  // ROOT IO constructor, don't allocate memory here!
}

/***************************************************************************/ /**
* Constructor.
******************************************************************************/
AliAnalysisTaskOmegaDielectron_AccEff::AliAnalysisTaskOmegaDielectron_AccEff(const char *name)
: AliAnalysisTaskSE(name),
// General member variables
fOutputList(nullptr), fCutVariList(0x0), fEvent(nullptr), fMCEvent(nullptr), fEventCuts(nullptr), fEventCuts_VertexZ(nullptr),
fFilter_TrackCuts(), fFilter_PID(),
fPIDResponse(nullptr),
fMinEta(-0.8), fMaxEta(0.8), fMinPt(0.2), fMaxPt(10.0),
fsysUnc(kFALSE),
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
fHist_MC_Omegas_gen(),
fHist_MC_Omegas_gen_DaughtersinAcc(),
fHist_Rec_Omegas_TrackPID(),

fHistVertex(),
fHist_MC_Omegas_Rapidity(),
fHist_MC_elec_gen(),
fHist_MC_posi_gen(),
fHist_MC_elec_gen_inAcc(),
fHist_MC_posi_gen_inAcc(),
fHist_elec_rec_inAcc(),
fHist_elec_rec_inAcc_Track(),
fHist_elec_rec_inAcc_Track_PID(),
fHist_posi_rec_inAcc(),
fHist_posi_rec_inAcc_Track(),
fHist_posi_rec_inAcc_Track_PID(),
fHist_MC_Omegas_withoutCuts(),
fHist_MC_Omegas_TrackCuts(),
fHist_MC_Omegas_TrackPID(),
fHist_Rec_Omegas_withoutCuts(),
fHist_Rec_Omegas_TrackCuts()
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


  std::vector<std::string> vector_EventCut_Names{"withoutEventCuts","EventCuts_without_VertexZ_pm10","withEventCuts"};
  int bins_pTe=800, bins_pTomega=160, bins_mee=200, bins_eta=20, bins_y=20, bins_phi=64;
  double edges_pTe[2]={0.0,8.0}, edges_pTomega[2]={0.0,8.0}, edges_mee[2]={0.0,2.0}, edges_eta[2]={-1.0,1.0}, edges_y[2]={-1.0,1.0}, edges_phi[2]={0.0,6.4};


  for (unsigned int i_Ev = 0; i_Ev < 3 ; ++i_Ev){
    TH3D* th3_tmp_fHist_MC_Omegas_gen = new TH3D(Form("Hist_MC_Omegas_gen_%s",vector_EventCut_Names.at(i_Ev).c_str()), Form("Hist_MC_Omegas_gen_%s;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});y",vector_EventCut_Names.at(i_Ev).c_str()), bins_mee,edges_mee[0],edges_mee[1], bins_pTomega,edges_pTomega[0],edges_pTomega[1], bins_y,edges_y[0],edges_y[1]);
    th3_tmp_fHist_MC_Omegas_gen->Sumw2();
    fHist_MC_Omegas_gen.push_back(th3_tmp_fHist_MC_Omegas_gen);
    fOutputList->Add(th3_tmp_fHist_MC_Omegas_gen);
  }

  for (unsigned int i_cv = 0; i_cv < fFilter_TrackCuts.size(); ++i_cv){

    TH3D* th3_tmp_fHist_Rec_Omegas_TrackPID = new TH3D(Form("Hist_Rec_Omegas_TrackPID_Cut%d", i_cv+1), Form("Hist_Rec_Omegas_TrackPID_Cut%d;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});y", i_cv+1),  bins_mee,edges_mee[0],edges_mee[1], bins_pTomega,edges_pTomega[0],edges_pTomega[1], bins_y,edges_y[0],edges_y[1]);
    th3_tmp_fHist_Rec_Omegas_TrackPID->Sumw2();
    fHist_Rec_Omegas_TrackPID.push_back(th3_tmp_fHist_Rec_Omegas_TrackPID);
    fOutputList->Add(th3_tmp_fHist_Rec_Omegas_TrackPID);


    if(!fsysUnc){
      // --- single electron histos: pt - eta - phi
      //th3_tmp_fHist_elec_rec_inAcc -- control Histo -- for every CutVariation the same
      TH3D* th3_tmp_fHist_elec_rec_inAcc = new TH3D(Form("Hist_elec_rec_inAcc_Cut%d", i_cv+1), Form("Hist_elec_rec_inAcc_Cut%d;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", i_cv+1),bins_pTe,edges_pTe[0],edges_pTe[1], bins_eta,edges_eta[0],edges_eta[1], bins_phi,edges_phi[0],edges_phi[1]); // 6.4 == 2*PI()
      th3_tmp_fHist_elec_rec_inAcc->Sumw2();
      fHist_elec_rec_inAcc.push_back(th3_tmp_fHist_elec_rec_inAcc);
      fOutputList->Add(th3_tmp_fHist_elec_rec_inAcc);
      TH3D* th3_tmp_fHist_elec_rec_inAcc_Track = new TH3D(Form("Hist_elec_rec_inAcc_Track_Cut%d", i_cv+1), Form("Hist_elec_rec_inAcc_Track_Cut%d;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", i_cv+1),bins_pTe,edges_pTe[0],edges_pTe[1], bins_eta,edges_eta[0],edges_eta[1], bins_phi,edges_phi[0],edges_phi[1]); // 6.4 == 2*PI()
      th3_tmp_fHist_elec_rec_inAcc_Track->Sumw2();
      fHist_elec_rec_inAcc_Track.push_back(th3_tmp_fHist_elec_rec_inAcc_Track);
      fOutputList->Add(th3_tmp_fHist_elec_rec_inAcc_Track);
      TH3D* th3_tmp_fHist_elec_rec_inAcc_Track_PID = new TH3D(Form("Hist_elec_rec_inAcc_Track_PID_Cut%d", i_cv+1), Form("Hist_elec_rec_inAcc_Track_PID_Cut%d;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", i_cv+1),bins_pTe,edges_pTe[0],edges_pTe[1], bins_eta,edges_eta[0],edges_eta[1], bins_phi,edges_phi[0],edges_phi[1]); // 6.4 == 2*PI()
      th3_tmp_fHist_elec_rec_inAcc_Track_PID->Sumw2();
      fHist_elec_rec_inAcc_Track_PID.push_back(th3_tmp_fHist_elec_rec_inAcc_Track_PID);
      fOutputList->Add(th3_tmp_fHist_elec_rec_inAcc_Track_PID);

      // --- single positron histos: pt - eta - phi
      //th3_tmp_fHist_posi_rec_inAcc -- control Histo -- for every CutVariation the same
      TH3D* th3_tmp_fHist_posi_rec_inAcc = new TH3D(Form("Hist_posi_rec_inAcc_Cut%d", i_cv+1), Form("Hist_posi_rec_inAcc_Cut%d;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", i_cv+1),bins_pTe,edges_pTe[0],edges_pTe[1], bins_eta,edges_eta[0],edges_eta[1], bins_phi,edges_phi[0],edges_phi[1]); // 6.4 == 2*PI()
      th3_tmp_fHist_posi_rec_inAcc->Sumw2();
      fHist_posi_rec_inAcc.push_back(th3_tmp_fHist_posi_rec_inAcc);
      fOutputList->Add(th3_tmp_fHist_posi_rec_inAcc);
      TH3D* th3_tmp_fHist_posi_rec_inAcc_Track = new TH3D(Form("Hist_posi_rec_inAcc_Track_Cut%d", i_cv+1), Form("Hist_posi_rec_inAcc_Track_Cut%d;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", i_cv+1),bins_pTe,edges_pTe[0],edges_pTe[1], bins_eta,edges_eta[0],edges_eta[1], bins_phi,edges_phi[0],edges_phi[1]); // 6.4 == 2*PI()
      th3_tmp_fHist_posi_rec_inAcc_Track->Sumw2();
      fHist_posi_rec_inAcc_Track.push_back(th3_tmp_fHist_posi_rec_inAcc_Track);
      fOutputList->Add(th3_tmp_fHist_posi_rec_inAcc_Track);
      TH3D* th3_tmp_fHist_posi_rec_inAcc_Track_PID = new TH3D(Form("Hist_posi_rec_inAcc_Track_PID_Cut%d", i_cv+1), Form("Hist_posi_rec_inAcc_Track_PID_Cut%d;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", i_cv+1),bins_pTe,edges_pTe[0],edges_pTe[1], bins_eta,edges_eta[0],edges_eta[1], bins_phi,edges_phi[0],edges_phi[1]); // 6.4 == 2*PI()
      th3_tmp_fHist_posi_rec_inAcc_Track_PID->Sumw2();
      fHist_posi_rec_inAcc_Track_PID.push_back(th3_tmp_fHist_posi_rec_inAcc_Track_PID);
      fOutputList->Add(th3_tmp_fHist_posi_rec_inAcc_Track_PID);

      /// --- --- --- histograms -> reconstructed Omega --- --- --- ///
      // m - pt - y
      TH3D* th3_tmp_fHist_MC_Omegas_TrackCuts = new TH3D(Form("Hist_MC_Omegas_TrackCuts_Cut%d", i_cv+1), Form("Hist_MC_Omegas_TrackCuts_Cut%d;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});y", i_cv+1),  bins_mee,edges_mee[0],edges_mee[1], bins_pTomega,edges_pTomega[0],edges_pTomega[1], bins_y,edges_y[0],edges_y[1]);
      th3_tmp_fHist_MC_Omegas_TrackCuts->Sumw2();
      fHist_MC_Omegas_TrackCuts.push_back(th3_tmp_fHist_MC_Omegas_TrackCuts);
      fOutputList->Add(th3_tmp_fHist_MC_Omegas_TrackCuts);
      TH3D* th3_tmp_fHist_MC_Omegas_TrackPID = new TH3D(Form("Hist_MC_Omegas_TrackPID_Cut%d", i_cv+1), Form("Hist_MC_Omegas_TrackPID_Cut%d;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});y", i_cv+1),  bins_mee,edges_mee[0],edges_mee[1], bins_pTomega,edges_pTomega[0],edges_pTomega[1], bins_y,edges_y[0],edges_y[1]);
      th3_tmp_fHist_MC_Omegas_TrackPID->Sumw2();
      fHist_MC_Omegas_TrackPID.push_back(th3_tmp_fHist_MC_Omegas_TrackPID);
      fOutputList->Add(th3_tmp_fHist_MC_Omegas_TrackPID);

      TH3D* th3_tmp_fHist_Rec_Omegas_TrackCuts = new TH3D(Form("Hist_Rec_Omegas_TrackCuts_Cut%d", i_cv+1), Form("Hist_Rec_Omegas_TrackCuts_Cut%d;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});y", i_cv+1),  bins_mee,edges_mee[0],edges_mee[1], bins_pTomega,edges_pTomega[0],edges_pTomega[1], bins_y,edges_y[0],edges_y[1]);
      th3_tmp_fHist_Rec_Omegas_TrackCuts->Sumw2();
      fHist_Rec_Omegas_TrackCuts.push_back(th3_tmp_fHist_Rec_Omegas_TrackCuts);
      fOutputList->Add(th3_tmp_fHist_Rec_Omegas_TrackCuts);
    }

  }

  fHist_MC_Omegas_gen_DaughtersinAcc = new TH3D("fHist_MC_Omegas_gen_DaughtersinAcc", "fHist_MC_Omegas_gen_DaughtersinAcc;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});y",  bins_mee,edges_mee[0],edges_mee[1], bins_pTomega,edges_pTomega[0],edges_pTomega[1], bins_y,edges_y[0],edges_y[1]);
  fOutputList->Add(fHist_MC_Omegas_gen_DaughtersinAcc);

  if(!fsysUnc){
    fHistVertex             = new TH2D("zVertex", "zVertex; Event selection ; z(cm)", 3, 0.5, 3.5, 300, -15.0, 15.0);
    fHistVertex->GetXaxis()->SetBinLabel(1, vector_EventCut_Names.at(0).c_str());
    fHistVertex->GetXaxis()->SetBinLabel(2, vector_EventCut_Names.at(1).c_str());
    fHistVertex->GetXaxis()->SetBinLabel(3, vector_EventCut_Names.at(2).c_str());
    fOutputList->Add(fHistVertex);

    // Histogram to check Rapidity distribution
    fHist_MC_Omegas_Rapidity = new TH1F("fHist_MC_Omegas_Rapidity", "Omega - Rapidity;y;#it{N}_{events}", 200, -10, 10);
    fOutputList->Add(fHist_MC_Omegas_Rapidity);

    /// --- --- --- MC histograms - generated/accepted Omegas --- --- --- ///

    fHist_MC_elec_gen = new TH3D("fHist_MC_elec_gen", "fHist_MC_elec_gen;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", bins_pTe,edges_pTe[0],edges_pTe[1], bins_eta,edges_eta[0],edges_eta[1], bins_phi,edges_phi[0],edges_phi[1]); // 6.4 == 2*PI()
    fOutputList->Add(fHist_MC_elec_gen);
    fHist_MC_elec_gen_inAcc = new TH3D("fHist_MC_elec_gen_inAcc", "fHist_MC_elec_gen_inAcc;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", bins_pTe,edges_pTe[0],edges_pTe[1], bins_eta,edges_eta[0],edges_eta[1], bins_phi,edges_phi[0],edges_phi[1]); // 6.4 == 2*PI()
    fOutputList->Add(fHist_MC_elec_gen_inAcc);

    fHist_MC_posi_gen = new TH3D("fHist_MC_posi_gen", "fHist_MC_posi_gen;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", bins_pTe,edges_pTe[0],edges_pTe[1], bins_eta,edges_eta[0],edges_eta[1], bins_phi,edges_phi[0],edges_phi[1]); // 6.4 == 2*PI()
    fOutputList->Add(fHist_MC_posi_gen);
    fHist_MC_posi_gen_inAcc = new TH3D("fHist_MC_posi_gen_inAcc", "fHist_MC_posi_gen_inAcc;#it{p}_{T,e} (GeV/#it{c});#eta;#phi", bins_pTe,edges_pTe[0],edges_pTe[1], bins_eta,edges_eta[0],edges_eta[1], bins_phi,edges_phi[0],edges_phi[1]); // 6.4 == 2*PI()
    fOutputList->Add(fHist_MC_posi_gen_inAcc);


    fHist_MC_Omegas_withoutCuts = new TH3D("fHist_MC_Omegas_withoutCuts", "fHist_MC_Omegas_withoutCuts;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T,ee} (GeV/#it{c});y",  bins_mee,edges_mee[0],edges_mee[1], bins_pTomega,edges_pTomega[0],edges_pTomega[1], bins_y,edges_y[0],edges_y[1]);
    fOutputList->Add(fHist_MC_Omegas_withoutCuts);

    fHist_Rec_Omegas_withoutCuts = new TH3D("fHist_Rec_Omegas_withoutCuts", "fHist_Rec_Omegas_withoutCuts;#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T,ee} (GeV/#it{c});y", bins_mee,edges_mee[0],edges_mee[1], bins_pTomega,edges_pTomega[0],edges_pTomega[1], bins_y,edges_y[0],edges_y[1]);
    fOutputList->Add(fHist_Rec_Omegas_withoutCuts);
  }


  /// --- --- --- EventCuts (for AOD/ ESD) --- --- --- ///
  fEventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && ncontrib>0");
  fEventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD);
  fEventCuts->SetRequireVertex();
  fEventCuts->SetMinVtxContributors(1);

  //fEventCuts_VertexZ=new AliDielectronEventCuts("eventCuts_VertexZ","|vtxZ|<10");
  fEventCuts_VertexZ=new AliDielectronEventCuts("eventCuts_VertexZ","Vertex Track && |vtxZ|<10 && ncontrib>0");
  fEventCuts_VertexZ->SetVertexType(AliDielectronEventCuts::kVtxSPD);
  fEventCuts_VertexZ->SetRequireVertex();
  fEventCuts_VertexZ->SetVertexZ(-10.,10.);
  fEventCuts_VertexZ->SetMinVtxContributors(1);

  //TO CHECK if cuts are set via config-file-> no InitCuts();
  // InitCuts();

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

  // for signal loss correction due to trigger -- >
  for (unsigned int i_EventC = 0; i_EventC < 3 ; ++i_EventC){
    const AliVVertex* vtx = fEvent->GetPrimaryVertex();
    Double_t vtxZGlobal = -99.;
    if (vtx) {
      vtxZGlobal = vtx->GetZ();
    }
    if (i_EventC == 1){
      if (!fEventCuts->IsSelected(fEvent)){ return; }
    }
    if (i_EventC == 2){
      if (!fEventCuts_VertexZ->IsSelected(fEvent)){ return; }
    }

    if(!fsysUnc){
      ///-------------------- Fill Events in Histogramm :: All;Selected--------------------///
      fHistVertex->Fill(i_EventC+1.0, vtxZGlobal); // all events Vertex z
    }



    ///-------------------- Loop over Generated MC Particles ---OMEGAS --------------------///
    for(Int_t iGenPart = 0; iGenPart < fMCEvent->GetNumberOfTracks();iGenPart++) {
      AliMCParticle* mcGenParticle  = (AliMCParticle*)fMCEvent->GetTrack(iGenPart);
      if(!mcGenParticle){
        AliError("mcGenParticle not available\n");
        continue;
      }

      Int_t pdgcode_gen = mcGenParticle->PdgCode();
      
      if(pdgcode_gen == fmother_pdg){  // Check if MC Particle is an omega

        if( CheckDielectronDecay(mcGenParticle, kFALSE) ){ // Check if :: w->e+e- and omega within y<+-0.8
          fHist_MC_Omegas_gen.at(i_EventC)->Fill(mcGenParticle->M(), mcGenParticle->Pt(), mcGenParticle->Y());
        }

        if(i_EventC == 2){ // only fill these Histos with the full event selection
          if( CheckDielectronDecay(mcGenParticle, kTRUE) ){ // Check if :: w->e+e-and omega within y<+-0.8 --> with /n_e/<0.8 & 0.2<p_T,e < 10GeV/c
            fHist_MC_Omegas_gen_DaughtersinAcc->Fill(mcGenParticle->M(), mcGenParticle->Pt(), mcGenParticle->Y());
          }
          if(!fsysUnc){
            fHist_MC_Omegas_Rapidity->Fill(mcGenParticle->Y());
          }
        }
      }
      if(!fsysUnc && i_EventC == 2){// only fill these Histos with the full event selection && ...
        // for the single electron efficiency:
        if (mcGenParticle->IsPhysicalPrimary() ){
          if(pdgcode_gen == felectron_pdg){  // Check if MC Particle is an electron
            fHist_MC_elec_gen->Fill(mcGenParticle->Pt(), mcGenParticle->Eta(), mcGenParticle->Phi());
            if (AcceptKinematics(mcGenParticle)){
              fHist_MC_elec_gen_inAcc->Fill(mcGenParticle->Pt(), mcGenParticle->Eta(), mcGenParticle->Phi());
            }
          }
          if(pdgcode_gen == fpositron_pdg){  // Check if MC Particle is a positron
            fHist_MC_posi_gen->Fill(mcGenParticle->Pt(), mcGenParticle->Eta(), mcGenParticle->Phi());
            if (AcceptKinematics(mcGenParticle)){
              fHist_MC_posi_gen_inAcc->Fill(mcGenParticle->Pt(), mcGenParticle->Eta(), mcGenParticle->Phi());
            }
          }
        }
      }
    }
  }
  
  if (!fEventCuts_VertexZ->IsSelected(fEvent)){ return; }


  ///-------------------- Loop over AOD Tracks - electrons// positrons --------------------///
  const int size_Cuts = fFilter_TrackCuts.size();
  // SelectionMask for the track and PID Cuts
  UInt_t selectedMask_Track[size_Cuts], selectedMask_PID[size_Cuts];
  UInt_t cutMask[size_Cuts], cutMask_PID[size_Cuts];
  UInt_t cutMask_electron[size_Cuts], cutMask_positron[size_Cuts];
  UInt_t cutMask_PID_electron[size_Cuts], cutMask_PID_positron[size_Cuts];

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
      continue; // filter bit 128 denotes TPC-only tracks, use only them for the analysis
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
    //loop for Cut Variations:

    for (unsigned int i_cutV_1 = 0; i_cutV_1 < fFilter_TrackCuts.size(); ++i_cutV_1){
      selectedMask_Track[i_cutV_1] =( 1 << fFilter_TrackCuts.at(i_cutV_1)->GetCuts()->GetEntries())-1; // cutting logic taken from AliDielectron::FillTrackArrays()   // apply track cuts
      selectedMask_PID[i_cutV_1] =( 1 << fFilter_PID.at(i_cutV_1)->GetCuts()->GetEntries())-1;
      cutMask[i_cutV_1] = fFilter_TrackCuts.at(i_cutV_1)->IsSelected(track);
      cutMask_PID[i_cutV_1] = fFilter_PID.at(i_cutV_1)->IsSelected(track);
      if(!fsysUnc){
        if(pdgcode == felectron_pdg){
          fHist_elec_rec_inAcc.at(i_cutV_1)->Fill(track->Pt(), track->Eta(), track->Phi());
          if (cutMask[i_cutV_1] == selectedMask_Track[i_cutV_1]) {
            fHist_elec_rec_inAcc_Track.at(i_cutV_1)->Fill(track->Pt(), track->Eta(), track->Phi());
            if (cutMask_PID[i_cutV_1] == selectedMask_PID[i_cutV_1]) {
              fHist_elec_rec_inAcc_Track_PID.at(i_cutV_1)->Fill(track->Pt(), track->Eta(), track->Phi());
            }
          }
        }
        if(pdgcode == fpositron_pdg){
          fHist_posi_rec_inAcc.at(i_cutV_1)->Fill(track->Pt(), track->Eta(), track->Phi());
          if (cutMask[i_cutV_1] == selectedMask_Track[i_cutV_1]) {
            fHist_posi_rec_inAcc_Track.at(i_cutV_1)->Fill(track->Pt(), track->Eta(), track->Phi());
            if (cutMask_PID[i_cutV_1] == selectedMask_PID[i_cutV_1]) {
              fHist_posi_rec_inAcc_Track_PID.at(i_cutV_1)->Fill(track->Pt(), track->Eta(), track->Phi());
            }
          }
        }
      }
    }


    if(pdgCode_mother == fmother_pdg){                        // Check if mother of MC Particle (Track) is an omega
      ////***************ELECTRONS*****************////
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
          if(!fsysUnc){
            fHist_Rec_Omegas_withoutCuts->Fill(mass, pairpt, pairy);
            fHist_MC_Omegas_withoutCuts->Fill(mcParticle_same_mother->M(), mcParticle_same_mother->Pt(), mcParticle_same_mother->Y());
          }

          /// ~~~~~~ ***** pass Track cuts *****  ~~~~~~ ///
          //loop for Cut Variations:
          for (unsigned int i_cutV_2 = 0; i_cutV_2 < fFilter_TrackCuts.size(); ++i_cutV_2){
            cutMask_electron[i_cutV_2] = fFilter_TrackCuts.at(i_cutV_2)->IsSelected(ele_from_same_mother_track);
            cutMask_positron[i_cutV_2] = fFilter_TrackCuts.at(i_cutV_2)->IsSelected(pos_from_same_mother_track);
            if (cutMask_electron[i_cutV_2] == selectedMask_Track[i_cutV_2] && cutMask_positron[i_cutV_2] == selectedMask_Track[i_cutV_2]) {
              if(!fsysUnc){
                fHist_Rec_Omegas_TrackCuts.at(i_cutV_2)->Fill(mass, pairpt, pairy);
                fHist_MC_Omegas_TrackCuts.at(i_cutV_2)->Fill(mcParticle_same_mother->M(), mcParticle_same_mother->Pt(), mcParticle_same_mother->Y());
              }


              /// ~~~~~~ ***** pass PID cuts *****  ~~~~~~ ///
              cutMask_PID_electron[i_cutV_2] = fFilter_PID.at(i_cutV_2)->IsSelected(ele_from_same_mother_track);
              cutMask_PID_positron[i_cutV_2] = fFilter_PID.at(i_cutV_2)->IsSelected(pos_from_same_mother_track);
              if (cutMask_PID_electron[i_cutV_2] == selectedMask_PID[i_cutV_2] && cutMask_PID_positron[i_cutV_2] == selectedMask_PID[i_cutV_2]) {
                fHist_Rec_Omegas_TrackPID.at(i_cutV_2)->Fill(mass, pairpt, pairy);
                // omega histo
                if(!fsysUnc){
                  fHist_MC_Omegas_TrackPID.at(i_cutV_2)->Fill(mcParticle_same_mother->M(), mcParticle_same_mother->Pt(), mcParticle_same_mother->Y());
                }
              }
            }

            else{continue;}
          }
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
