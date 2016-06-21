///
/// \file AliFemtoCutMonitorPionPion.cxx
///


#include "AliFemtoCutMonitorPionPion.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoAvgSepCalculator.h"

#include "AliFemtoPairCutDetaDphi.h"

#include "AliFemtoEvent.h"

static const double PionMass = 0.13956995;

#include <TAxis.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TH2I.h>
#include <TObjArray.h>

#include <map>
#include <string>


AliFemtoCutMonitorPionPion::Event::Event(const bool passing,
                                         const bool is_identical_analysis,
                                         const bool is_mc_analysis,
                                         const bool suffix_output):
  AliFemtoCutMonitor()
  , fCentMult(nullptr)
  , fVertexZ(nullptr)
  , fVertexXY(nullptr)
  , _collection_size_pass(nullptr)
  , _collection_size_fail(nullptr)
  , _identical_collection_size_pass(nullptr)
  , _identical_collection_size_fail(nullptr)
  , _prev_ev(nullptr)
  , _prev_pion_coll_1_size(0)
  , _prev_pion_coll_2_size(0)
{
  const char *title_suffix = (passing ? " (PASS)" : " (FAIL)");

  const TString pf(suffix_output ? passing ? "_P" : "_F" : "");

  fCentMult = new TH2F(
    "cent_mult" + pf,
    TString::Format("Event Centrality vs Multiplicity%s;"
                    "centrality; multiplicity (N_{tracks}); N_{ev}", title_suffix),
    100, 0, 100.0,
    100, 0, 10000.0
  );
  fCentMult->Sumw2();

  fVertexZ = new TH1F(
    "VertexZ" + pf,
    TString::Format("Vertex Z Distribution%s;z (cm);N_{ev}", title_suffix),
    128, -15.0f, 15.0f
  );
  fVertexZ->Sumw2();

  fVertexXY = new TH2F(
    "VertexXY" + pf,
    TString::Format("Vertex XY Distribution%s;x (cm);y (cm); dN/(dx $\\cdot$ dy)", title_suffix),
    48, 0.05f, 0.08f,
    48, 0.31f, 0.345f
    // 48, 0.0f, 0.12f,
    // 48, 0.22f, 0.32f
  );
  fVertexXY->Sumw2();

  // only create _collection_size histograms if this is the passing event cut monitor
  if (passing) {
    // if identical - then skip this
    if (is_identical_analysis) {
      _identical_collection_size_pass = new TH1I("collection_size_p",
                                                 "Size of particle collection;"
                                                 "# pions; N_{ev}",
                                                 100, -0.5, 800.5);
      _identical_collection_size_pass->Sumw2();
      _identical_collection_size_fail = (TH1I*)_identical_collection_size_pass->Clone("collection_size_f");
    } else {
      _collection_size_pass = new TH2I("collection_size_p",
                                       "Size of Particle Collection in Passing Events;"
                                       "# pions (1);"
                                       "# pions (2);"
                                       "N_{ev}",
                                       100, -0.5, 800.5,
                                       100, -0.5, 800.5);
      _collection_size_fail = (TH2I*)_collection_size_pass->Clone("collection_size_f");
    }
  }
}


TList*
AliFemtoCutMonitorPionPion::Event::GetOutputList()
{
  TList *olist = new TList();
  TCollection *output = olist;

  output->Add(fCentMult);
  output->Add(fVertexZ);
  output->Add(fVertexXY);

  if (_collection_size_pass) {
     output->Add(_collection_size_pass);
     output->Add(_collection_size_fail);
  } else if (_identical_collection_size_pass) {
    output->Add(_identical_collection_size_pass);
    output->Add(_identical_collection_size_fail);
  }

  return olist;
}


void
AliFemtoCutMonitorPionPion::Event::Fill(const AliFemtoEvent* ev)
{
  const Float_t centrality = ev->CentralityV0();
  const Int_t multiplicty = ev->NumberOfTracks();
  const AliFemtoThreeVector vertex = ev->PrimVertPos();

  fCentMult->Fill(centrality, multiplicty);
  fVertexZ->Fill(vertex.z());
  fVertexXY->Fill(vertex.x(), vertex.y());

  if (ev == _prev_ev) {
    if (_collection_size_pass) {
      _collection_size_pass->Fill(_prev_pion_coll_1_size, _prev_pion_coll_2_size);
      _prev_ev = nullptr;
    } else if (_identical_collection_size_pass) {
      _identical_collection_size_pass->Fill(_prev_pion_coll_1_size);
    }
  }

}

void
AliFemtoCutMonitorPionPion::Event::EventBegin(const AliFemtoEvent* ev)
{
  if (_collection_size_pass == nullptr && _identical_collection_size_pass == nullptr) {
    return;
  }
  _prev_ev = ev;
}


void
AliFemtoCutMonitorPionPion::Event::EventEnd(const AliFemtoEvent* ev)
{

  // We were not called with the previous event - must have gone to failed
  if (_prev_ev != nullptr) {

    // this cut monitor does not monitor collection size
    if (_collection_size_pass != nullptr) {
      _collection_size_fail->Fill(_prev_pion_coll_1_size, _prev_pion_coll_2_size);
    } else if (_identical_collection_size_pass) {
      _identical_collection_size_fail->Fill(_prev_pion_coll_1_size);
    }

    _prev_ev = nullptr;
  }
}


void
AliFemtoCutMonitorPionPion::Event::Fill(const AliFemtoParticleCollection *coll_1,
                                        const AliFemtoParticleCollection *coll_2)
{
  _prev_pion_coll_2_size = coll_2->size();
  _prev_pion_coll_1_size = coll_1->size();
}


//
// Pion Cut Monitor
//

#define UNSPECIFIED_PARTICLE -999999999

const std::map<Int_t, std::string> code_to_label = {
  {UNSPECIFIED_PARTICLE, "other"},
  {0, "none"},
  {11, "e^{-}"},
  {-11, "e^{+}"},
  {13, "#mu^{-}"},
  {-13, "#mu^{+}"},
  {15, "#tau^{-}"},
  {-15, "#tau^{+}"},
  {-211, "#pi^{-}"},
  {211, "#pi^{+}"},
  {311, "K^{0}"},
  {321, "K^{+}"},
  {-321, "K^{-}"},
  {113, "#rho^{0}"},
  {213, "#rho^{+}"},
  {-213, "#rho^{-}"},
  {221, "#eta"},
  {223, "#omega"},
  {2212, "p^{+}"},
  {-2212, /* "\\bar{p}" */ "p^{-}"},
  {3122, "#Lambda"},
  {-3122, "$\\bar{\\Lambda}$"},
  {3222, "#Sigma^{+}"},
  {3212, "#Sigma^{0}"},
  {3112, "#Sigma^{-}"},
};

const std::vector<Int_t> codes = {
  UNSPECIFIED_PARTICLE,
  11, -11, 13, -13, -211, 211, -321, 321, -2212, 2212,
};

const std::vector<Int_t> parent_codes = {
  UNSPECIFIED_PARTICLE,  // unknown
  0,  // no parent
   -211, 211, 221, 223, 113, 213, -213, 311, 321, -321,
  3122, -3122, 3222, 3212, 3112
};



AliFemtoCutMonitorPionPion::Pion::Pion(const bool passing,
                                       const TString& typestr,
                                       const bool is_mc_analysis,
                                       const bool suffix_output):
  AliFemtoCutMonitor()
  , fYPt(nullptr)
  , fPtPhi(nullptr)
  , fEtaPhi(nullptr)
  , fChiTpcIts(nullptr)
  , fdEdX(nullptr)
  , fMC_mass(nullptr)
  , fMC_pt(nullptr)
  , fMC_type(nullptr)
  , fMC_parent(nullptr)
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
                    "#eta  vs  p_{T}",
             /*X*/  "#eta;"
             /*Y*/  "p_{T} (GeV);"
             /*Z*/  "dN/(p_{T} $\\cdot$ \\eta)"),
    140, -1.4, 1.4,
    100, 0, 3.0);
  fYPt->Sumw2();

  fPtPhi = new TH2F(
    "PtPhi" + pf,
    TString::Format(title_format,
                    "Pt vs Phi",
                    "#phi (rads);"
                    "p_{T} (GeV);"),
    144, -TMath::Pi(), TMath::Pi(),
    144,  0.0, 3.0);
  fPtPhi->Sumw2();

  fEtaPhi = new TH2F(
    "EtaPhi" + pf,
    TString::Format(title_format,
                    "#eta vs Phi",
                    "#phi (rads);"
                    "#eta;"),
    144, -TMath::Pi(), TMath::Pi(),
    144, -1.4, 1.4);
  fEtaPhi->Sumw2();

  fChiTpcIts = new TH2F(
    "ChiTpcIts" + pf,
    TString::Format(title_format,
                    "#chi^{2} / N_{cls} TPC vs ITS",
                    "TPC; ITS;"),
    144, 0.0, 0.1,
    144, 0.0, 0.1);
  fChiTpcIts->Sumw2();

  fdEdX = new TH2F(
    "dEdX" + pf,
    Form(title_format,
         "dE/dx vs p",
         "p (GeV);"
         "dE/dx;"
         "N_{tracks}"),
    128, 0, 6.0,
    128, 0, 500.0);
  fdEdX->Sumw2();

  fImpact = new TH2F(
    "impact" + pf,
    Form(title_format,
         "Track impact parameter components",
         "z (cm?); r (cm?); N_{#pi}"
        ),
    256, -0.25, 0.25,
    256, -0.01, 0.1);
  fImpact->Sumw2();

  if (is_mc_analysis) {
    fMC_mass = new TH1F(
      "mc_Mass" + pf,
      TString::Format(title_format, "M_{inv}",
                                    "M_{inv} (GeV);"
                                    "N_{#pi}"),
      144, 0.0, 1.5
      // 144, 0.0120, 0.5,
    );
    fMC_mass->Sumw2();

    fMC_pt = new TH2F(
      "mc_Pt" + pf,
      TString::Format(title_format, "p_{T}",
                                    "p_{T}^{reconstrcted};"
                                    "p_{T}^{true}"),
      144,  0.0, 3.0,
      144,  0.0, 3.0
    );
    fMC_pt->Sumw2();

    fMC_type = new TH1I(
      "mc_pdg",
      TString::Format(title_format,
        "PDG Code; Code; N_{code};"),
      codes.size(), -0.5, codes.size() - 0.5
    );

    for (Int_t bin = 0; bin < codes.size(); bin++) {
      Int_t code = codes[bin];
      auto label = code_to_label.find(code);
      fMC_type->GetXaxis()->SetBinLabel(bin + 1, label->second.c_str());
    }
    fMC_type->GetXaxis()->CenterLabels();

    fMC_parent = new TH2I(
      "mc_parent_pdg",
      TString::Format(title_format,
        "Parent PDG Code; Daughter PID; Parent;"),
        codes.size(), -0.5, codes.size() - 0.5,
        parent_codes.size(), -0.5, parent_codes.size() - 0.5
    );

    for (Int_t bin = 0; bin < codes.size(); bin++) {
      Int_t code = codes[bin];
      auto label = code_to_label.find(code);
      fMC_parent->GetXaxis()->SetBinLabel(bin + 1, label->second.c_str());
    }
    fMC_parent->GetXaxis()->CenterLabels();


    for (Int_t bin = 0; bin < parent_codes.size(); bin++) {
      Int_t code = parent_codes[bin];
      auto label = code_to_label.find(code);
      fMC_parent->GetYaxis()->SetBinLabel(bin + 1, label->second.c_str());
    }
    fMC_parent->GetYaxis()->CenterLabels();
  }

}


TList*
AliFemtoCutMonitorPionPion::Pion::GetOutputList()
{
  TList *olist = new TList();
  TCollection *output = olist;

  output->Add(fYPt);
  output->Add(fPtPhi);
  output->Add(fEtaPhi);
  output->Add(fChiTpcIts);
  output->Add(fdEdX);
  output->Add(fImpact);
  if (fMC_type) {
    output->Add(fMC_mass);
    output->Add(fMC_pt);
    output->Add(fMC_type);
    output->Add(fMC_parent);
  }

  return olist;
};


void AliFemtoCutMonitorPionPion::Pion::Fill(const AliFemtoTrack* track)
{
  const float pz = track->P().z(),
              pt = track->Pt(),
               p = track->P().Mag(),
             phi = track->P().Phi();

  const double energy = ::sqrt(track->P().Mag2() + PionMass * PionMass),
                  eta = 0.5 * ::log((energy + pz) / (energy - pz));


  const Int_t ITS_ncls = track->ITSncls(),
              TPC_ncls = track->TPCncls();


  if (fMC_mass) {
    const AliFemtoModelHiddenInfo *mc = dynamic_cast<const AliFemtoModelHiddenInfo*>(track->GetHiddenInfo());

    fMC_mass->Fill(mc->GetMass());
    fMC_pt->Fill(pt, mc->GetTrueMomentum()->Perp());

    const auto type_location = std::find(codes.begin(), codes.end(), mc->GetPDGPid());
    const Int_t type_bin = type_location == codes.end()
                         ? 0
                         : std::distance(codes.begin(), type_location);
    fMC_type->Fill(type_bin);


    auto parent_location = std::find(parent_codes.begin(), parent_codes.end(), mc->GetMotherPdgCode());
    const Int_t parent_type_bin = parent_location == parent_codes.end()
                                ? 0
                                : std::distance(parent_codes.begin(), parent_location);
    // if (parent_location == parent_codes.end() && mc->GetPDGPid() == 211) {
    // std::cout << " " << mc->GetPDGPid()
    //           << " <- " << mc->GetMotherPdgCode()
    //           << "\n";
    // }
    fMC_parent->Fill(type_bin, parent_type_bin);
  }

  fYPt->Fill(eta, pt);
  fPtPhi->Fill(phi, pt);
  fEtaPhi->Fill(phi, eta);
  fdEdX->Fill(p, track->TPCsignal());

  fChiTpcIts->Fill( (TPC_ncls > 0) ? track->TPCchi2() / TPC_ncls : 0.0,
                    (ITS_ncls > 0) ? track->ITSchi2() / ITS_ncls : 0.0);


  fImpact->Fill(track->ImpactZ(), track->ImpactD());
}


AliFemtoCutMonitorPionPion::Pair::Pair(const bool passing,
                                       const TString& typestr,
                                       const bool is_mc_analysis,
                                       const bool suffix_output):
  AliFemtoCutMonitor()
  , fMinv(nullptr)
  , fKt(nullptr)
  , fDetaDphi(nullptr)
  , fQinvDeta(nullptr)
  , fQinvDphiStar(nullptr)
  , fMCTrue_minv(nullptr)
  , fMCTrue_qinv(nullptr)
{
  const TString title_format = TString::Format("%s %%s %s; %%s",
                                               typestr.Data(),
                                               (passing ? "(PASS)" : "(FAIL)"));
  const TString pf(suffix_output ? passing ? "_P" : "_F" : "");

  fMinv = new TH1F(
    "Pair_Minv" + pf,
    TString::Format(title_format, "M_{inv}", "M_{inv} (GeV); N_{pairs}"),
    288, 0.0, 8.0);
  fMinv->Sumw2();

  fKt = new TH1F(
    "kt" + pf,
    TString::Format(title_format,
                    "k_{T} Distribution",
                    "k_{T} (GeV); N_{pairs}"),
    144, 0.0, 4.0);
  fKt->Sumw2();

  fDetaDphi = new TH2F(
    "DetaDphi" + pf,
    TString::Format(title_format,
                    "#Delta #eta* vs #Delta #phi*",
                    "#Delta #eta*; #Delta #phi*"),
    145, -0.2, 0.2,
    145, -0.2, 0.2
  );
  fDetaDphi->Sumw2();

  fQinvDeta = new TH2F(
    "QinvDeta" + pf,
    TString::Format(title_format,
                    "Q_{inv} vs #Delta #eta",
                    "Q_{inv} (GeV); #Delta #eta"),
    100, 0.0, 1.2,
    75, -0.1, 0.1
  );
  fDetaDphi->Sumw2();

  fQinvDphiStar = new TH2F(
    "QinvDphiStar" + pf,
    TString::Format(title_format,
                    "Q_{inv} vs #Delta #phi*",
                    "Q_{inv} (GeV); #Delta #phi*"),
    100, 0.0, 1.2,
    75, -0.1, 0.1
  );
  fQinvDphiStar->Sumw2();

  if (is_mc_analysis) {
    fMCTrue_minv = new TH2F(
      "mc_Minv" + pf,
      TString::Format(title_format,
        "Minv True vs Reconstructed",
        "M_{inv}^{r} (GeV);"
        "M_{inv}^{t} (Gev);"),
      144, 0.0, 4.5,
      144, 0.0, 4.5);
    fMCTrue_minv->Sumw2();

    fMCTrue_qinv = new TH2F(
      "mc_Qinv" + pf,
      TString::Format(title_format,
        "q_{inv} True vs Reconstructed",
        "q_{inv}^{r} (GeV);"
        "q_{inv}^{t} (Gev);"),
      400, 0.0, 1.0,
      400, 0.0, 1.0);
    fMCTrue_qinv->Sumw2();
  }
}

void
AliFemtoCutMonitorPionPion::Pair::Fill(const AliFemtoPair *pair)
{
  const float minv = pair->MInv(),
             kstar = pair->KStar(),
              qinv = fabs(pair->QInv());

  fMinv->Fill(minv);
  fKt->Fill(pair->KT());

  const AliFemtoTrack *track1 = pair->Track1()->Track(),
                      *track2 = pair->Track2()->Track();

  const AliFemtoThreeVector p1 = track1->P(),
                            p2 = track2->P();

  const short chg_1 = track1->Charge(),
              chg_2 = track2->Charge();

  const float delta_eta = AliFemtoPairCutDetaDphi::CalculateDEta(p1, p2),
         delta_phi_star = AliFemtoPairCutDetaDphi::CalculateDPhiStar(p1, chg_1, p2, chg_2, 1.6, 0.5006670488586);

  fDetaDphi->Fill(delta_eta, delta_phi_star);
  fQinvDeta->Fill(qinv, delta_eta);
  fQinvDphiStar->Fill(qinv, delta_phi_star);

  if (fMCTrue_minv) {
    const AliFemtoModelHiddenInfo *mc_1 = dynamic_cast<const AliFemtoModelHiddenInfo*>(pair->Track1()->HiddenInfo()),
                                  *mc_2 = dynamic_cast<const AliFemtoModelHiddenInfo*>(pair->Track2()->HiddenInfo());
    if (mc_1 == nullptr || mc_2 == nullptr) {
      return;
    }

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
                                p2(e2, momentum_2);

    const double reconstructed_qinv = qinv,
                          true_qinv = -1.0 * (p1 - p2).m(),
                          true_minv = (p1 + p2).m();

    // skip mis-identified protons
    if (mc_1->GetPDGPid() == 2212 || mc_2->GetPDGPid() == 2212) {
      return;
    }
    if (mc_1->GetPDGPid() != 211 || mc_2->GetPDGPid() != 211) {
      return;
    }

    fMCTrue_qinv->Fill(qinv, true_qinv);
    fMCTrue_minv->Fill(minv, true_minv);

//     if (0.2 < (qinv - true_qinv)) {
//         printf(" => %6d %6d\n", mc_1->GetPDGPid(), mc_2->GetPDGPid());
// //        cout << " => " << mc_1->GetPDGPid() << "  " << mc_2->GetPDGPid() << "\n";
//     }
  }
}


TList* AliFemtoCutMonitorPionPion::Pair::GetOutputList()
{
  TList *olist = new TList();
  TCollection *output = olist;

  output->Add(fMinv);
  output->Add(fKt);
  output->Add(fDetaDphi);
  output->Add(fQinvDeta);
  output->Add(fQinvDphiStar);

  if (fMCTrue_qinv) {
    output->Add(fMCTrue_qinv);
  }

  if (fMCTrue_minv) {
    output->Add(fMCTrue_minv);
  }

  return olist;
};
