///
/// \file AliFemtoPionLambdaCutMonitor.cxx
///


#include "AliFemtoPionLambdaCutMonitor.h"
#include "AliFemtoModelHiddenInfo.h"

#include "AliFemtoEvent.h"


static const double PionMass = 0.13956995;

#include <TH1F.h>
#include <TH2F.h>
#include <TH2I.h>
#include <TObjArray.h>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoPionLambdaCutMonitor);
  /// \endcond
#endif

AliFemtoPionLambdaCutMonitor::Event::Event(const bool passing,
                                           const bool is_mc_analysis,
                                           const bool suffix_output):
  AliFemtoCutMonitor()
  , _centrality(NULL)
  , _multiplicity(NULL)
  , _vertex_z(NULL)
  , _vertex_xy(NULL)
  , _collection_size_pass(NULL)
  , _collection_size_fail(NULL)
  , _prev_ev(NULL)
  , _prev_pion_coll_size(0)
  , _prev_lam_coll_size(0)
{

  const char *title_suffix = (passing ? " (PASS)" : " (FAIL)");

  const TString pf(suffix_output ? passing ? "_P" : "_F" : "");

  _centrality = new TH1F(
    "centrality" + pf,
    TString::Format("Event Centrality%s", title_suffix),
    100, 0, 100.0
  );
  _centrality->Sumw2();

  _multiplicity = new TH1F(
    "multiplicity" + pf,
    "Event Multiplicity",
    100, 0, 10000.0
  );
  _multiplicity->Sumw2();

  _vertex_z = new TH1F(
    "VertexZ" + pf,
    TString::Format("Vertex Z Distribution%s;z (cm);dN/dz", title_suffix),
    128, -15.0f, 15.0f
  );
  _vertex_z->Sumw2();

  _vertex_xy = new TH2F(
    "VertexXY" + pf,
    TString::Format("Vertex XY Distribution%s;x (cm);y (cm); dN/(dx $\\cdot$ dy)", title_suffix),
    48, -0.3f, 0.3f,
    48, -0.3f, 0.3f
  );
  _vertex_xy->Sumw2();

  // only create _collection_size histograms if this is the passing event cut monitor
  if (passing) {
    _collection_size_pass = new TH2I("collection_size_p",
                                     "Size of Particle Collection in Passing Events;"
                                     "# pions;"
                                     "# lambdas;",
                                     100, -0.5, 2000.5,
                                     11, -0.5, 10.5);
    _collection_size_pass->Sumw2();
    _collection_size_fail = (TH2I*)_collection_size_pass->Clone("collection_size_f");
  } else {
     _collection_size_pass = NULL;
     _collection_size_fail = NULL;
  }
}

TList*
AliFemtoPionLambdaCutMonitor::Event::GetOutputList()
{
  TList *olist = new TList();
  TCollection *output = olist;

  output->Add(_centrality);
  output->Add(_multiplicity);
  output->Add(_vertex_z);
  output->Add(_vertex_xy);

  if (_collection_size_pass) {
     output->Add(_collection_size_pass);
     output->Add(_collection_size_fail);
  }

  return olist;
}

void
AliFemtoPionLambdaCutMonitor::Event::Fill(const AliFemtoEvent* ev)
{
  const Float_t centrality = ev->CentralityV0();
  const Int_t multiplicty = ev->NumberOfTracks();
  const AliFemtoThreeVector vertex = ev->PrimVertPos();

  _centrality->Fill(centrality);
  _multiplicity->Fill(multiplicty);
  _vertex_z->Fill(vertex.z());
  _vertex_xy->Fill(vertex.x(), vertex.y());

  if (_collection_size_pass && ev == _prev_ev) {
    _collection_size_pass->Fill(_prev_pion_coll_size, _prev_lam_coll_size);
    _prev_ev = NULL;
  }
}

void AliFemtoPionLambdaCutMonitor::Event::EventBegin(const AliFemtoEvent* ev)
{
  if (_collection_size_pass == NULL) {
    return;
  }
  _prev_ev = ev;
}

void AliFemtoPionLambdaCutMonitor::Event::EventEnd(const AliFemtoEvent* ev)
{
  // this cut monitor does not monitor collection size
  if (_collection_size_pass == NULL) {
    return;
  }

  // We were not called with the previous event - must have gone to failed
  if (_prev_ev != NULL) {
    _collection_size_fail->Fill(_prev_pion_coll_size, _prev_lam_coll_size);
    _prev_ev = NULL;
  }
}


void
AliFemtoPionLambdaCutMonitor::Event::Fill(const AliFemtoParticleCollection *coll_1,
                                          const AliFemtoParticleCollection *coll_2)
{
  _prev_lam_coll_size = coll_2->size();
  _prev_pion_coll_size = coll_1->size();
}

AliFemtoPionLambdaCutMonitor::Pion::Pion(const bool passing,
                                         const TString& typestr,
                                         const bool is_mc_analysis,
                                         const bool suffix_output):
  AliFemtoCutMonitor()
  , fYPt(NULL)
  , fPtPhi(NULL)
  , fEtaPhi(NULL)
  , fdEdX(NULL)
  , fMinv(NULL)
{
  // Build 'standard' format for histogram titles
  //  <ParticleType> <Title> <Pass/Fail>; <AxisInfo>
  const TString title_format = TString::Format("%s %%s %s; %%s",
                                               typestr.Data(),
                                               (passing ? "(PASS)" : "(FAIL)"));
  const TString pf(suffix_output ? passing ? "_P" : "_F" : "");

  fYPt = new TH2F(
    "eta_Pt" + pf,
    TString::Format(title_format,
                    "$\\eta  vs  p_{T}$",
             /*X*/  "$\\eta$;"
             /*Y*/  "p_{T} (GeV);"
             /*Z*/  "dN/(p_{T} $\\cdot$ \\eta)"),
     140, -1.4, 1.4,
     100, 0, 3.0);
  fYPt->Sumw2();

  fPtPhi = new TH2F(
    "PtPhi" + pf,
    TString::Format(title_format,
                    "Pt vs Phi",
                    "\\phi (rads);"
                    "p_{T} (GeV);"),
    144, -TMath::Pi(), TMath::Pi(),
    144,  0.0, 3.0);
  fPtPhi->Sumw2();

  fEtaPhi = new TH2F(
    "EtaPhi" + pf,
    TString::Format(title_format,
                    "\\eta vs Phi",
                    "\\phi (rads);"
                    "\\eta;"),
    144, -TMath::Pi(), TMath::Pi(),
    144, -1.4, 1.4);
  fEtaPhi->Sumw2();

  fdEdX = new TH2F(
    "dEdX" + pf,
    Form(title_format,
         "dE/dx vs p",
         "p (GeV);"
         "dE/dx;"
         "dN/(p_{T} $\\cdot$ dE/dx)"),
     128, 0, 6.0,
     128, 0, 500.0);
  fdEdX->Sumw2();

  if (is_mc_analysis) {
    fMinv = new TH1F(
      "mc_Mass" + pf,
      TString::Format(title_format, "M_{inv}",
                                    "M_{inv} (GeV);"
                                    "dN/dM"),
      144, 0.0120, 1.158);
    fMinv->Sumw2();
  }

}


TList*
AliFemtoPionLambdaCutMonitor::Pion::GetOutputList()
{
  TList *olist = new TList();
  TCollection *output = olist;

  output->Add(fYPt);
  output->Add(fPtPhi);
  output->Add(fEtaPhi);
  output->Add(fdEdX);
  if (fMinv) {
    output->Add(fMinv);
  }

  return olist;
};

void AliFemtoPionLambdaCutMonitor::Pion::Fill(const AliFemtoTrack* track)
{
  const float pz = track->P().z(),
              pt = track->Pt(),
               p = track->P().Mag(),
             phi = track->P().Phi();

  const double energy = ::sqrt(track->P().Mag2() + PionMass * PionMass),
                  eta = 0.5 * ::log((energy + pz) / (energy - pz));


  if (fMinv) {
    fMinv->Fill(track->GetMass());
  }

  fYPt->Fill(eta, pt);
  fPtPhi->Fill(phi, pt);
  fEtaPhi->Fill(phi, eta);
  fdEdX->Fill(p, track->TPCsignal());
}

AliFemtoPionLambdaCutMonitor::Lambda::Lambda(const bool passing,
                                             const TString& typestr,
                                             const AliFemtoAnalysisPionLambda::LambdaType ltype,
                                             const bool is_mc_analysis,
                                             const bool suffix_output):
  AliFemtoCutMonitor()
  , fLambdaType(ltype)
  , _minv(NULL)
  , _ypt(NULL)
  , _dedx_p_pro(NULL)
  , _dedx_p_pi(NULL)
  , fCosPointingAngle(NULL)
  , fMCTrue_minv(NULL)
  , fMCTrue_ypt(NULL)
{
  // Build 'standard' format for histogram titles
  //  <ParticleType> <Title> <Pass/Fail>; <AxisInfo>
  const TString title_format = TString::Format("%s %%s %s; %%s",
                                         typestr.Data(),
                                         (passing ? "(PASS)" : "(FAIL)"));
  const TString pf(suffix_output ? passing ? "_P" : "_F" : "");

  _minv = new TH1F(
    "V0_Minv" + pf,
    TString::Format(title_format,
                    "M_{inv}",
                    "M_{inv} (GeV); dN/dM"),
    576, 1.070, 1.140);
  _minv->Sumw2();

  _ypt = new TH2F(
    "V0_YPt" + pf,
    TString::Format(title_format,
                    "\\eta vs p_{T}",
                    "\\eta;"
                    "p_{T} (GeV);"
                    "dN/(p_{T} $\\cdot$ \\eta)"),
    140, -1.4, 1.4,
    100, 0.0, 3.0);
  _ypt->Sumw2();

  _dedx_p_pro = new TH2F(
    "dEdX_Pro" + pf,
    TString::Format(title_format,
                    "dE/dx vs p (Proton Daughter)",
                    "p (GeV);"
                    "dE/dx;"
                    "dN/(p_{T} $\\cdot$ dE/dx)"),
     128, 0, 6.0,
     128, 0, 500);
  _dedx_p_pro->Sumw2();

  _dedx_p_pi = new TH2F(
    "dEdX_Pi" + pf,
    TString::Format(title_format,
                    "dE/dx vs p (Pion Daughter)",
                    "p (GeV);"
                    "dE/dx;"
                    "dN/(p_{T} $\\cdot$ dE/dx)"),
     128, 0, 6.0,
     128, 0, 500.0);
  _dedx_p_pi->Sumw2();

  fCosPointingAngle = new TH1F(
    "CosPointingAngle",
    TString::Format(title_format,
                    "Cosine Pointing Angle",
                    "Cos(\\Theta);"),
     60, 0.98, 1.002);
  fCosPointingAngle->Sumw2();

  if (is_mc_analysis) {
    fMCTrue_minv = _minv = new TH1F(
      "mc_V0_Minv" + pf,
      TString::Format(title_format,
                      "(MC) M_{inv}",
                      "M_{inv} (GeV); dN/dM"),
      576, 1.070, 1.140);
    fMCTrue_minv->Sumw2();

    fMCTrue_ypt = new TH2F(
      "mc_V0_YPt" + pf,
      TString::Format(title_format,
                      "(MC) \\eta vs p_{T}",
                      "\\eta;"
                      "p (GeV);"
                      "dN/(p_{T} $\\cdot$ \\eta)"),
      140, -1.4, 1.4,
      100, 0.0, 3.0);
    fMCTrue_minv->Sumw2();
  }
}


void
AliFemtoPionLambdaCutMonitor::Lambda::Fill(const AliFemtoV0* track)
{
  const bool type_is_lambda = (fLambdaType == AliFemtoAnalysisPionLambda::kLambda);

  const Float_t eta = track->EtaV0(),
                 pt = track->PtV0(),

               minv = (type_is_lambda) ? track->MassLambda() : track->MassAntiLambda(),

             p_pion = (type_is_lambda) ? track->MomNeg().Mag() : track->MomPos().Mag(),
           p_proton = (type_is_lambda) ? track->MomPos().Mag() : track->MomNeg().Mag(),

//         eta_pion = (type_is_lambda) ? track->EtaNeg() : track->EtaPos(),
//       eta_proton = (type_is_lambda) ? track->EtaPos() : track->EtaNeg(),

          dedx_pion = (type_is_lambda) ? track->DedxNeg() : track->DedxPos(),
        dedx_proton = (type_is_lambda) ? track->DedxPos() : track->DedxNeg();
/*
       nsig_tpc_pro = (type_is_lambda) ? track->PosNSigmaTPCP() : track->NegNSigmaTPCP(),
       nsig_tof_pro = (type_is_lambda) ? track->PosNSigmaTOFP() : track->NegNSigmaTOFP();

      nsig_tpc_pion = (type_is_lambda) ? track->NegNSigmaTPCPi() : track->PosNSigmaTPCPi(),
      nsig_tof_pion = (type_is_lambda) ? track->NegNSigmaTOFPi() : track->PosNSigmaTOFPi();
*/

  _minv->Fill(minv);
  _ypt->Fill(eta, pt);
  _dedx_p_pro->Fill(p_proton, dedx_proton);
  _dedx_p_pi->Fill(p_pion, dedx_pion);
  fCosPointingAngle->Fill(track->CosPointingAngle());

  if (fMCTrue_minv) {
    const AliFemtoModelHiddenInfo *mc_data = dynamic_cast<AliFemtoModelHiddenInfo*>(track->GetHiddenInfo());
    if (mc_data) {
      fMCTrue_minv->Fill(mc_data->GetMass());
      fMCTrue_ypt->Fill(mc_data->GetTrueMomentum()->PseudoRapidity(),
                        mc_data->GetTrueMomentum()->Perp());
    }
  }


}


TList*
AliFemtoPionLambdaCutMonitor::Lambda::GetOutputList()
{
  TList *olist = new TList();
  TCollection *output = olist;

  output->Add(_minv);
  output->Add(_ypt);
  output->Add(_dedx_p_pro);
  output->Add(_dedx_p_pi);
  output->Add(fCosPointingAngle);

  if (fMCTrue_minv) {
    output->Add(fMCTrue_minv);
    output->Add(fMCTrue_ypt);
  }

  return olist;
}


AliFemtoPionLambdaCutMonitor::Pair::Pair(const bool passing,
                                         const TString& typestr,
                                         const bool is_mc_analysis,
                                         const bool suffix_output
                                         ):
  AliFemtoCutMonitor()
  , _minv(NULL)
  , fKt(NULL)
  , fAvgSep_pion(NULL)
  , fAvgSep_proton(NULL)
  , fMCTrue_minv(NULL)
  , fMCTrue_kstar(NULL)
{
  const TString title_format = TString::Format("%s %%s %s; %%s",
                                               typestr.Data(),
                                               (passing ? "(PASS)" : "(FAIL)"));
  const TString pf(suffix_output ? passing ? "_P" : "_F" : "");

  _minv = new TH1F(
    "Pair_Minv" + pf,
    TString::Format(title_format, "M_{inv}", "M_{inv} (GeV)"),
    288, 0.0, 8.0);
  _minv->Sumw2();

  fKt = new TH1F(
    "kt" + pf,
    TString::Format(title_format,
                    "k_{T} Distribution",
                    "k_{T} (GeV); dN/k_{T}"),
   144, 0.0, 4.0);
  fKt->Sumw2();

  fAvgSep_pion = new TH1F(
    "AvgSep_pi" + pf,
    TString::Format(title_format,
      "AvgSep Pion Daughter", "Average Separation (cm)"),
    144, 0.0, 20.0);
  fAvgSep_pion->Sumw2();

  fAvgSep_proton = new TH1F(
    "AvgSep_pro" + pf,
    TString::Format(title_format,
      "AvgSep Proton Daughter", "Average Separation (cm)"),
    144, 0.0, 20.0);
  fAvgSep_proton->Sumw2();

  if (is_mc_analysis) {
    fMCTrue_minv = new TH2F(
      "mc_Minv" + pf,
      TString::Format(title_format,
        "Minv True vs Reconstructed",
        "M_{inv}^{r} (GeV);"
        "M_{inv}^{t} (Gev);"),
      144, 1.0, 6.0,
      144, 1.0, 6.0);
    fMCTrue_minv->Sumw2();

    fMCTrue_kstar = new TH2F(
      "mc_Kstar" + pf,
      TString::Format(title_format,
        "K* True vs Reconstructed",
        "K*^{r} (GeV);"
        "K*^{t} (Gev);"),
      144, 0.0, 4.0,
      144, 0.0, 4.0);
    fMCTrue_kstar->Sumw2();
  }
}

void
AliFemtoPionLambdaCutMonitor::Pair::Fill(const AliFemtoPair *pair)
{
  const float minv = pair->MInv(),
             kstar = pair->KStar();

  _minv->Fill(minv);
  fKt->Fill(pair->KT());

  if (fMCTrue_minv) {
    const AliFemtoModelHiddenInfo *mc_1 = dynamic_cast<const AliFemtoModelHiddenInfo*>(pair->Track1()->HiddenInfo()),
                                  *mc_2 = dynamic_cast<const AliFemtoModelHiddenInfo*>(pair->Track2()->HiddenInfo());
    if (mc_1 && mc_2) {

      // const AliFemtoV0 *lambda = pair->Track1()->V0();
      // const AliFemtoTrack *pion = pair->Track2()->Track();

      // get true momentums
      const AliFemtoThreeVector &momentum_1 = *mc_1->GetTrueMomentum(),
                                &momentum_2 = *mc_2->GetTrueMomentum();

      // get true mass and calculate true energy
      const float m1 = mc_1->GetMass(),
                  e1 = TMath::Sqrt(m1 * m1 + momentum_1.Mag2()),

                  m2 = mc_2->GetMass(),
                  e2 = TMath::Sqrt(m2 * m2 + momentum_2.Mag2());

      // build momentum four-vectors
      const AliFemtoLorentzVector p1(e1, momentum_1),
                                  p2(e2, momentum_2),
                                  p_sum = p1 + p2,
                                  p_diff = p1 - p2;

      // Minv is just the magnitude of the momentum sum
      fMCTrue_minv->Fill(minv, p_sum.m());

      const float tQ = pow(m1 * m1 - m2 * m2, 2) / p_sum.m2(),
                  q2 = tQ - p_diff.m2(),

                  mc_kstar = q2 > 0
                           ? TMath::Sqrt(q2) / 2.0
                           : 0.0;

      // kstar calculation
      fMCTrue_kstar->Fill(kstar, mc_kstar);

    }
  }

}


TList*
AliFemtoPionLambdaCutMonitor::Pair::GetOutputList()
{
  TList *olist = new TList();
  TCollection *output = olist;

  output->Add(_minv);
  output->Add(fKt);

  if (fMCTrue_kstar) {
    output->Add(fMCTrue_kstar);
  }

  if (fMCTrue_minv) {
    output->Add(fMCTrue_minv);
  }

  return olist;
};

/*
void AliFemtoPionLambdaCutMonitor::_init()
{

  //
  // Event Histograms
  //
  _event_centrality_pass = new TH1F(
    "centrality_pass",
    "Event Centrality",
    100, 0, 100.0);
  _event_centrality_fail = (TH1F*)_event_centrality_pass->Clone("centrality_fail");

  _event_multiplicity_pass = new TH1F(
    "MultPass",
    "Event Multiplicty",
    100, 0, 10000.0);
  _event_multiplicity_fail = (TH1F*)_event_multiplicity_pass->Clone("MultFail");

  _event_vertex_z_pass = new TH1F(
    "VertexZPass",
    "Vertex Z Distribution;z (cm);dN/dz",
    128, -15.0f, 15.0f);
  _event_vertex_z_fail = (TH1F*)_event_vertex_z_pass->Clone("VertexZFail");

  _event_vertex_xy_pass = new TH2F(
    "VertexXYPass",
    "Vertex XY Distribution (Pass);x (cm);y (cm); dN/(dx $\\cdot$ dy)",
    48, -1.0f, 1.0f,
    48, -1.0f, 1.0f);
  _event_vertex_xy_fail = (TH2F*)_event_vertex_xy_pass->Clone("VertexXYFail");


  //
  // Pion Histograms
  //

  _pion_minv_pass = new TH1F(
    "MinvPass",
    "Pion M_{inv}; M_{inv} (GeV); dN/dM",
    144, 0.120, .158);
  _pion_minv_fail = (TH1F*)_pion_minv_pass->Clone("MinvFail");

  _pion_dedx_p_pass = new TH2F(
    "dEdX_Pt_Pass",
    "dE/dx vs p_{T}; p_{T} (GeV); dE/dx; dN/(p_{T} $\\cdot$ dE/dx)",
     128, 0, 500,
     128, 0, 1.5);
  _pion_dedx_pt_fail = (TH2F*)_pion_dedx_p_pass->Clone("dEdX_Pt_Fail");

  _pion_ypt_pass = new TH2F(
    "YPt_Pass",
    "Rapidity vs Pt; Y; p_{T} (GeV); dN/(p_{T} $\\cdot$ Y)",
     140, -1.4, 1.4,
     100, 0, 3.0);
  _pion_ypt_fail = (TH2F*)_pion_ypt_pass->Clone("YPt_Fail");


  _lambda_minv_pass = new TH1F(
    "V0_MinvPass",
    "Lambda M_{inv}",
    576, 1.07, 1.118);
  _lambda_minv_fail = (TH1F*)_lambda_minv_pass->Clone("V0_MinvFail");

  _lambda_ypt_pass = new TH2F(
    "V0_YPt_Pass",
    "Rapidity vs p_{T}; Y; p_{T} (GeV); dN/(p_{T} $\\cdot$ Y)",
    140, -1.4, 1.4,
    100, 0.0, 3.0);
  _lambda_ypt_fail = (TH2F*)_lambda_ypt_pass->Clone("V0_YPt_Fail");

  _lambda_dedx_pt_pro_pass = new TH2F(
    "dEdX_Pro_Pass",
    "dE/dx vs p_{T} (Proton Daughter); p_{T} (GeV); dE/dx",
     128, 0, 2.2,
     128, 0, 500);
  _lambda_dedx_p_pro_fail = (TH2F*)_lambda_dedx_p_pro_pass->Clone("dEdX_Pro_Fail");

  _lambda_dedx_pt_pi_pass = new TH2F(
    "dEdX_Pi_Pass",
    "dE/dx vs p_{T} (Pion Daughter); p_{T} (GeV); dE/dx",
     128, 0, 2.2,
     128, 0, 500.0);
  _lambda_dedx_pt_pi_fail = (TH2F*)_lambda_dedx_pt_pi_pass->Clone("dEdX_Pi_Fail");


  //
  // Pair Histograms
  //
  _pair_minv_pass = new TH1F(
    "Pair_Minv_Pass",
    "M_{inv} (Pass);M_{inv} (GeV)",
    288, 0.0, 8.0);
  _pair_minv_fail = (TH1F*)_pair_minv_pass->Clone("Pair_Minv_Fail");
}

AliFemtoPionLambdaCutMonitor::AliFemtoPionLambdaCutMonitor()
{
  // _init();
}

TList* AliFemtoPionLambdaCutMonitor::GetOutputList()
{
  TList *olist = new TList();
  TObjArray *event_objs = new TObjArray();
  event_objs->SetName("Event");
  event_objs->Add(_event_centrality_pass);
  event_objs->Add(_event_centrality_fail);
  event_objs->Add(_event_multiplicity_pass);
  event_objs->Add(_event_multiplicity_fail);
  event_objs->Add(_event_vertex_z_pass);
  event_objs->Add(_event_vertex_z_fail);
  event_objs->Add(_event_vertex_xy_pass);
  event_objs->Add(_event_vertex_xy_fail);
  olist->Add(event_objs);

  TObjArray *pion_objs = new TObjArray();
  pion_objs->SetName("Pion");
  pion_objs->Add(_pion_minv_pass);
  pion_objs->Add(_pion_minv_fail);
  pion_objs->Add(_pion_dedx_p_pass);
  pion_objs->Add(_pion_dedx_pt_fail);
  pion_objs->Add(_pion_ypt_pass);
  pion_objs->Add(_pion_ypt_fail);
  olist->Add(pion_objs);

  TObjArray *lambda_objs = new TObjArray();
  lambda_objs->SetName("Lambda");
  lambda_objs->Add(_lambda_minv_pass);
  lambda_objs->Add(_lambda_minv_fail);
  lambda_objs->Add(_lambda_ypt_pass);
  lambda_objs->Add(_lambda_ypt_fail);
  lambda_objs->Add(_lambda_dedx_pt_pro_pass);
  lambda_objs->Add(_lambda_dedx_pt_pro_fail);
  lambda_objs->Add(_lambda_dedx_pt_pi_pass);
  lambda_objs->Add(_lambda_dedx_pt_pi_fail);
  olist->Add(lambda_objs);

  TObjArray *pair_objs = new TObjArray();
  pair_objs->SetName("Pair");
  pair_objs->Add(_pair_minv_pass);
  pair_objs->Add(_pair_minv_fail);
  olist->Add(pair_objs);
  return olist;
}

void AliFemtoPionLambdaCutMonitor::Fill(const AliFemtoEvent* ev)
{
  const Float_t centrality = ev->CentralityV0();
  const Int_t multiplicty = ev->NumberOfTracks();
  const AliFemtoThreeVector vertex = ev->PrimVertPos();

  _event_centrality->Fill(centrality);
  _event_multiplicity->Fill(multiplicty);
  _event_vertex_z_pass->Fill(vertex.z());
  _event_vertex_xy_pass->Fill(vertex.x(), vertex.y());

}

void AliFemtoPionLambdaCutMonitor::Fill(const AliFemtoTrack* aTrack)
{}

void AliFemtoPionLambdaCutMonitor::Fill(const AliFemtoV0* aV0)
{}

void AliFemtoPionLambdaCutMonitor::Fill(const AliFemtoPair* aPair)
{}
*/
