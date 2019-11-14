///
/// \file AliFemtoCutMonitorPionPion.cxx
///


#include "AliFemtoCutMonitorPionPion.h"
#include "AliFemtoModelHiddenInfo.h"
// #include "AliFemtoAvgSepCalculator.h"

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
#include <vector>

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

  fVertexZ = new TH1F(
    "VertexZ" + pf,
    TString::Format("Vertex Z Distribution%s;z (cm);N_{ev}", title_suffix),
    128, -15.0f, 15.0f
  );

  fVertexXY = new TH2F(
    "VertexXY" + pf,
    TString::Format("Vertex XY Distribution%s;x (cm);y (cm); dN/(dx $\\cdot$ dy)", title_suffix),
    48, 0.05f, 0.095f,
    48, 0.31f, 0.355f
    // 48, 0.0f, 0.12f,
    // 48, 0.22f, 0.32f
  );

  // only create _collection_size histograms if this is the passing event cut monitor
  if (passing) {
    // if identical - then skip this
    if (is_identical_analysis) {
      _identical_collection_size_pass = new TH1I("collection_size_p",
                                                 "Size of particle collection;"
                                                 "# pions; N_{ev}",
                                                 100, -0.5, 1000.5);
      _identical_collection_size_fail = (TH1I*)_identical_collection_size_pass->Clone("collection_size_f");
    } else {
      _collection_size_pass = new TH2I("collection_size_p",
                                       "Size of Particle Collection in Passing Events;"
                                       "# pions (1);"
                                       "# pions (2);"
                                       "N_{ev}",
                                       100, -0.5, 1000.5,
                                       100, -0.5, 1000.5);
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

const static int UNSPECIFIED_PARTICLE = -999999999;

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
                                       const bool suffix_output,
                                       const bool wide_impact_range):
  AliFemtoCutMonitor()
  , fAllowCharge(0)
  , fYPt(nullptr)
  , fPtPhi(nullptr)
  , fEtaPhi(nullptr)
  , fChi2Tpc(nullptr)
  , fChi2Its(nullptr)
  , fChiTpcIts(nullptr)
  , fClsTpcIts(nullptr)

  , fPidProbPion(nullptr)
  , fPidProbKaon(nullptr)
  , fPidProbProton(nullptr)
  , fPidProbElectron(nullptr)

  , fdEdX(nullptr)
  , fTofVsP(nullptr)
  , fTofPionVsP(nullptr)
  , fTofKaonVsP(nullptr)
  , fTofProtonVsP(nullptr)
  , fTpcTofPionSigma(nullptr)
  , fTpcTofKaonSigma(nullptr)

  , fTofMass(nullptr)

  , fNsigPionTpc(nullptr)
  , fNsigKaonTpc(nullptr)
  , fNsigProtonTpc(nullptr)

  , fNsigPionTof(nullptr)
  , fNsigKaonTof(nullptr)
  , fNsigProtonTof(nullptr)

  , fImpact(nullptr)
  , fEtaY(nullptr)
  , fMC_mass(nullptr)
  , fMC_pt(nullptr)
  , fMC_rap(nullptr)
  , fMC_type(nullptr)
  , fMC_parent(nullptr)
{
  // Build 'standard' format for histogram titles
  //  <ParticleType> <Title> <Pass/Fail>; <AxisInfo>
  const TString title_format = TString::Format("%s %%s %s; %%s",
                                               typestr.Data(),
                                               (passing ? "(PASS)" : "(FAIL)"));
  const TString pf(suffix_output ? passing ? "_P" : "_F" : "");

  const auto hist_name = [&] (const TString &name)
    {
      return name + pf;
    };

  const auto hist_title = [&] (const char *title, const char *axes)
    {
      return TString::Format(title_format, title, axes);
    };

  fYPt = new TH2F(
    hist_name("EtaPt"),
    hist_title("#eta  vs  p_{T}",
                /*X*/  "#eta;"
                /*Y*/  "p_{T} (GeV);"
                /*Z*/  "dN/(p_{T} $\\cdot$ \\eta)"),
    140, -1.4, 1.4,
    100, 0, 3.0);

  fPtPhi = new TH2F(
    hist_name("PtPhi"),
    hist_title("Pt vs Phi",
               "#phi (rads);"
               "p_{T} (GeV);"),
    144, -TMath::Pi(), TMath::Pi(),
    144,  0.0, 3.0);

  fEtaPhi = new TH2F(
    hist_name("EtaPhi"),
    hist_title("#eta vs Phi",
              "#phi (rads);"
              "#eta;"),
    144, -TMath::Pi(), TMath::Pi(),
    144, -1.4, 1.4);

  fChi2Tpc = new TH1F(
    hist_name("Chi2Tpc"),
    hist_title("#chi^{2} / N_{cls} TPC", "TPC"),
    144, 0.0, 5.0
  );

  fChi2Its = new TH1F(
    hist_name("Chi2Its"),
    hist_title("#chi^{2} / N_{cls} ITS", "ITS"),
    144, 0.0, 5.0
  );

  fChiTpcIts = new TH2F(
    hist_name("ChiTpcIts"),
    hist_title("#chi^{2} / N_{DoF} TPC vs ITS", "TPC; ITS;"),
    144, 0.0, 6.1,
    144, 0.0, 7.1);

  fClsTpcIts = new TH2F(
    hist_name("ClsTpcIts"),
    hist_title("N-Clusters ITS vs TPC", "TPC; ITS;"),
    161, -0.5, 160.5,
    11, -0.5, 10.5);

#if PIDPROB_ENABLED
  fPidProbPion = new TH2F(
    hist_name("PidProbPion"),
    hist_title("Pion PID Probability Vs P", "p (GeV); Prob;"),
    128, 0.0, 4.0,
    140, 0.0, 1.0);

  fPidProbKaon = new TH2F(
    hist_name("PidProbKaon"),
    hist_title("Kaon PID Probability Vs P", "p (GeV); Prob;"),
    128, 0.0, 4.0,
    140, 0.0, 1.0);

  fPidProbProton = new TH2F(
    hist_name("PidProbProton"),
    hist_title("Proton PID Probability Vs P", "p (GeV); Prob;"),
    128, 0.0, 4.0,
    140, 0.0, 1.0);

  fPidProbElectron = new TH2F(
    hist_name("PidProbElectron"),
    hist_title("Electron PID Probability Vs P", "p (GeV); Prob;"),
    128, 0.0, 4.0,
    140, 0.0, 1.0);
#endif

  fdEdX = new TH2F(
    hist_name("dEdX"),
    hist_title("dE/dx vs p",
               "p (GeV);"
               "dE/dx;"
               "N_{tracks}"),
    128, 0, 4.0,
    128, 0, 500.0
  );

  fTofVsP = new TH2F(
    hist_name("TofVsP"),
    hist_title("TOF Time vs p",
               "p (GeV);"
               "TOF Time;"
               "N_{tracks}"),
    256, 0, 4.0,
    3*256, 0.0, 101000.0
  );

  fTofPionVsP = new TH2F(
    hist_name("TofPiVsP"),
    hist_title("TOF Time - T_{pi} vs p",
               "p (GeV);"
               "TOF Time;"
               "N_{tracks}"),
    256, 0, 4.0,
    3*256, -30000.0, 101000.0
  );

  fTofKaonVsP = new TH2F(
    hist_name("TofKVsP"),
    hist_title("(TOF Time - T_{kaon} vs p",
               "p (GeV);"
               "TOF Time;"
               "N_{tracks}"),
    256, 0, 4.0,
    3*256, -30000.0, 101000.0
  );

  fTofProtonVsP = new TH2F(
    hist_name("TofPVsP"),
    hist_title("(TOF Time - T_{proton} vs p",
               "p (GeV);"
               "TOF Time;"
               "N_{tracks}"),
    256, 0, 4.0,
    3*256, -30000.0, 101000.0
  );

  const int sig_nbins_per_unit = 4;
  const double sig_target_range = 50,
               sig_binsize = 1.0 / sig_nbins_per_unit,
               sig_nbins = std::ceil(sig_target_range * 2 * sig_nbins_per_unit),
               sig_max = sig_binsize * sig_nbins / 2;

  fTpcTofPionSigma = new TH2F(
    hist_name("TpcTofSigmaPion"),
    hist_title("(TOF #sigma_{pion} vs TPC #sigma_{pion} ",
               "TPC #sigma;"
               "TOF #sigma;"
               "N_{tracks}"),
    sig_nbins, -sig_max, sig_max,
    sig_nbins, -sig_max, sig_max);

  fTpcTofKaonSigma = new TH2F(
    hist_name("TpcTofSigmaKaon"),
    hist_title("(TOF #sigma_{kaon} vs TPC #sigma_{kaon} ",
               "TPC #sigma;"
               "TOF #sigma;"
               "N_{tracks}"),
    sig_nbins, -sig_max, sig_max,
    sig_nbins, -sig_max, sig_max);

  fTpcTofProtonSigma = new TH2F(
    hist_name("TpcTofSigmaProton"),
    hist_title("(TOF #sigma_{proton} vs TPC #sigma_{proton} ",
               "TPC #sigma;"
               "TOF #sigma;"
               "N_{tracks}"),
    sig_nbins, -sig_max, sig_max,
    sig_nbins, -sig_max, sig_max);

  fNsigPionTof = new TH2F(
    hist_name("TofSigmaPion"),
    hist_title("TOF #sigma_{pion} vs p",
               "p (GeV);"
               "TOF #sigma;"
               "N_{tracks}"),
    128, 0, 4.0,
    sig_nbins, -sig_max, sig_max);

  fNsigKaonTof = new TH2F(
    hist_name("TofSigmaKaon"),
    hist_title("TOF #sigma_{kaon} vs p",
               "p (GeV);"
               "TOF #sigma;"
               "N_{tracks}"),
    128, 0, 4.0,
    sig_nbins, -sig_max, sig_max);

  fNsigProtonTof = new TH2F(
    hist_name("TofSigmaProton"),
    hist_title("TOF #sigma_{proton} vs p",
               "p (GeV);"
               "TOF #sigma;"
               "N_{tracks}"),
    128, 0, 4.0,
    sig_nbins, -sig_max, sig_max);

  fNsigPionTpc = new TH2F(
    hist_name("TpcSigmaPion"),
    hist_title("TPC #sigma_{pion} vs p",
               "p (GeV);"
               "TPC #sigma;"
               "N_{tracks}"),
    128, 0, 4.0,
    sig_nbins, -sig_max / 2.0, sig_max / 2.0);

  fNsigKaonTpc = new TH2F(
    hist_name("TpcSigmaKaon"),
    hist_title("TPC #sigma_{kaon} vs p",
               "p (GeV);"
               "TPC #sigma;"
               "N_{tracks}"),
    128, 0, 4.0,
    sig_nbins, -sig_max, sig_max);

  fNsigProtonTpc = new TH2F(
    hist_name("TpcSigProton"),
    hist_title("TOF #sigma_{proton} vs p",
               "p (GeV);"
               "TOF #sigma_{proton};"
               "N_{tracks}"),
    128, 0, 4.0,
    sig_nbins, -sig_max, sig_max);


  const double impact_range = wide_impact_range ? 3.25 : 0.25;

  fImpact = new TH2F(
    hist_name("impact"),
    hist_title("Track impact parameter components",
                "z (cm); "
                "r (cm); "
                "N_{#pi} "),
    257, -impact_range, impact_range,
    129, 0, impact_range);

  fEtaY = new TH2F(
    hist_name("eta_y"),
    hist_title("Rapidity vs PseudoRapidity",
               "pseudorapidity, #eta; rapidity, y"),
    400, -2.1, 2.1,
    400, -2.1, 2.1);

  fTofMass = new TH1F(
    hist_name("MassTOF"),
    hist_title("TOF Mass", "Mass (GeV);"),
    400, 0.0, 1.1);

  if (is_mc_analysis) {
    fMC_mass = new TH1F(
      hist_name("mc_Mass"),
      hist_title("M_{inv}",
                "M_{inv} (GeV);"
                "N_{#pi}"),
      144, 0.0, 1.5);

    fMC_pt = new TH2F(
      hist_name("mc_Pt"),
      hist_title("p_{T}",
                 "p_{T}^{reconstrcted};"
                 "p_{T}^{true}"),
      144,  0.0, 3.0,
      144,  0.0, 3.0);

    fMC_type = new TH1I(
      hist_name("mc_pdg"),
      hist_title("PDG Code",
                 "Code;"
                 "N_{code};"),
      codes.size(), -0.5, codes.size() - 0.5);

    fMC_rap = new TH2F(
      hist_name("mc_rapidity"),
      hist_title("Ideal Rapidity vs Rapidity",
                 "rapidity (assumed pion mass); ideal rapidity (true mass)"),
      400, -2.1, 2.1,
      400, -2.1, 2.1);

    for (UInt_t bin = 0; bin < codes.size(); bin++) {
      Int_t code = codes[bin];
      auto label = code_to_label.find(code);
      fMC_type->GetXaxis()->SetBinLabel(bin + 1, label->second.c_str());
    }
    fMC_type->GetXaxis()->CenterLabels();

#ifndef MC_PARENT_IS_THSPARSE
    fMC_parent = new TH2I(
      hist_name("mc_parent_pdg"),
      hist_title("Parent PDG Code",
                 "Daughter PID;"
                 "Parent;"),
      codes.size(), -0.5, codes.size() - 0.5,
      parent_codes.size(), -0.5, parent_codes.size() - 0.5
    );

    for (UInt_t bin = 0; bin < codes.size(); bin++) {
      Int_t code = codes[bin];
      auto label = code_to_label.find(code);
      fMC_parent->GetXaxis()->SetBinLabel(bin + 1, label->second.c_str());
    }
    fMC_parent->GetXaxis()->CenterLabels();

    for (UInt_t bin = 0; bin < parent_codes.size(); bin++) {
      Int_t code = parent_codes[bin];
      auto label = code_to_label.find(code);
      fMC_parent->GetYaxis()->SetBinLabel(bin + 1, label->second.c_str());
    }
    fMC_parent->GetYaxis()->CenterLabels();

#else
    const int max_code = 100555;
    const static int ndim = 2;
    const Int_t nbins[ndim] = {2 * max_code + 1, 2 * max_code + 1};
    const Double_t low[ndim] = {-max_code - 0.5, -max_code - 0.5},
                  high[ndim] = {max_code + 0.5, max_code + 0.5};

    fMC_parent = new THnSparseI(
      hist_name("mc_parent_pdg"),
      hist_title("parent PDG code",
                 "daughter Code;"
                 "parent Code;"),
      ndim,
      nbins, low, high);
#endif
  }

}

AliFemtoCutMonitorPionPion::Pion::Pion(const Pion &orig):
  AliFemtoCutMonitor()
  , fAllowCharge(orig.fAllowCharge)
  , fYPt(static_cast<TH2F*>(orig.fYPt->Clone()))
  , fPtPhi(static_cast<TH2F*>(orig.fPtPhi->Clone()))
  , fEtaPhi(static_cast<TH2F*>(orig.fEtaPhi->Clone()))
  , fChi2Tpc(static_cast<TH1F*>(orig.fChi2Tpc->Clone()))
  , fChi2Its(static_cast<TH1F*>(orig.fChi2Its->Clone()))
  , fChiTpcIts(static_cast<TH2F*>(orig.fChiTpcIts->Clone()))
  , fClsTpcIts(static_cast<TH2F*>(orig.fClsTpcIts->Clone()))
  , fdEdX(static_cast<TH2F*>(orig.fdEdX->Clone()))

#if PIDPROB_ENABLED
  , fPidProbPion(static_cast<TH2F*>(orig.fPidProbPion->Clone()))
  , fPidProbKaon(static_cast<TH2F*>(orig.fPidProbKaon->Clone()))
  , fPidProbProton(static_cast<TH2F*>(orig.fPidProbProton->Clone()))
  , fPidProbElectron(static_cast<TH2F*>(orig.fPidProbElectron->Clone()))
#endif

  , fTofVsP(static_cast<TH2F*>(orig.fTofVsP->Clone()))
  , fTofPionVsP(static_cast<TH2F*>(orig.fTofPionVsP->Clone()))
  , fTofKaonVsP(static_cast<TH2F*>(orig.fTofKaonVsP->Clone()))
  , fTofProtonVsP(static_cast<TH2F*>(orig.fTofProtonVsP->Clone()))

  , fTpcTofPionSigma(static_cast<TH2F*>(orig.fTpcTofPionSigma->Clone()))
  , fTpcTofKaonSigma(static_cast<TH2F*>(orig.fTpcTofKaonSigma->Clone()))
  , fTpcTofProtonSigma(static_cast<TH2F*>(orig.fTpcTofProtonSigma->Clone()))
  , fTofMass(static_cast<TH1F*>(orig.fTofMass->Clone()))

  , fNsigPionTpc(static_cast<TH2F*>(orig.fNsigPionTpc->Clone()))
  , fNsigKaonTpc(static_cast<TH2F*>(orig.fNsigKaonTpc->Clone()))
  , fNsigProtonTpc(static_cast<TH2F*>(orig.fNsigProtonTpc->Clone()))

  , fNsigPionTof(static_cast<TH2F*>(orig.fNsigPionTof->Clone()))
  , fNsigKaonTof(static_cast<TH2F*>(orig.fNsigKaonTof->Clone()))
  , fNsigProtonTof(static_cast<TH2F*>(orig.fNsigProtonTof->Clone()))

  , fImpact(static_cast<TH2F*>(orig.fImpact->Clone()))
  , fEtaY(static_cast<TH2F*>(orig.fEtaY->Clone()))
  , fMC_mass(static_cast<TH1F*>(orig.fMC_mass ? orig.fMC_mass->Clone(): nullptr))
  , fMC_pt(static_cast<TH2F*>(orig.fMC_pt ? orig.fMC_pt->Clone(): nullptr))
  , fMC_rap(static_cast<TH2F*>(orig.fMC_rap ? orig.fMC_rap->Clone() : nullptr))
  , fMC_type(static_cast<TH1I*>(orig.fMC_type ? orig.fMC_type->Clone(): nullptr))
  , fMC_parent(static_cast<THnSparseI*>(orig.fMC_parent ? orig.fMC_parent->Clone(): nullptr))
{
}

TList*
AliFemtoCutMonitorPionPion::Pion::GetOutputList()
{
  TList *olist = new TList();
  TCollection *output = olist;

  output->Add(fYPt);
  output->Add(fPtPhi);
  output->Add(fEtaPhi);
  output->Add(fChi2Tpc);
  output->Add(fChi2Its);
  output->Add(fChiTpcIts);
  output->Add(fClsTpcIts);

#if PIDPROB_ENABLED
  output->Add(fPidProbPion);
  output->Add(fPidProbKaon);
  output->Add(fPidProbProton);
  output->Add(fPidProbElectron);
#endif

  output->Add(fdEdX);
  output->Add(fNsigPionTpc);
  output->Add(fNsigKaonTpc);
  output->Add(fNsigProtonTpc);

  output->Add(fTofMass);
  output->Add(fTofVsP);
  output->Add(fTofPionVsP);
  output->Add(fTofKaonVsP);
  output->Add(fTofProtonVsP);

  output->Add(fNsigPionTof);
  output->Add(fNsigKaonTof);
  output->Add(fNsigProtonTof);

  output->Add(fTpcTofPionSigma);
  output->Add(fTpcTofKaonSigma);
  output->Add(fTpcTofProtonSigma);

  output->Add(fImpact);
  output->Add(fEtaY);
  if (fMC_type) {
    output->Add(fMC_mass);
    output->Add(fMC_pt);
    output->Add(fMC_rap);
    output->Add(fMC_type);
    output->Add(fMC_parent);
  }

  return olist;
};


void AliFemtoCutMonitorPionPion::Pion::Fill(const AliFemtoTrack* track)
{
  if (fAllowCharge && fAllowCharge != track->Charge()) {
    return;
  }

  const float pz = track->P().z(),
              pt = track->Pt(),
               p = track->P().Mag(),
             eta = track->P().PseudoRapidity(),
             phi = track->P().Phi();

  const double energy = ::sqrt(p * p + PionMass * PionMass),
             rapidity = 0.5 * ::log((energy + pz) / (energy - pz));

  const Int_t TPC_ncls = track->TPCncls();
  const Int_t ITS_ncls = track->ITSncls();

  if (fMC_mass) {
    const auto &mc = static_cast<const AliFemtoModelHiddenInfo&>(*track->GetHiddenInfo());

    fMC_mass->Fill(mc.GetMass());
    fMC_pt->Fill(pt, mc.GetTrueMomentum()->Perp());

    auto pdg_code = mc.GetPDGPid(),
         pdg_code_parent = mc.GetMotherPdgCode();

    const auto type_location = std::find(codes.begin(), codes.end(), pdg_code);
    const Int_t type_bin = type_location == codes.end()
                         ? 0
                         : std::distance(codes.begin(), type_location);
    fMC_type->Fill(type_bin);

    #ifndef MC_PARENT_IS_THSPARSE
    auto parent_location = std::find(parent_codes.begin(), parent_codes.end(), pdg_code_parent);
    const Int_t parent_type_bin = parent_location == parent_codes.end()
                                ? 0
                                : std::distance(parent_codes.begin(), parent_location);
    fMC_parent->Fill(type_bin, parent_type_bin);
    #else
    Double_t value[] = {static_cast<double>(pdg_code), static_cast<double>(pdg_code_parent)};
    fMC_parent->Fill(value);
    #endif

    const AliFemtoThreeVector &ideal_p = *mc.GetTrueMomentum();
    const double
      ipz = ideal_p.z(),
      iE = ideal_p.MassHypothesis(mc.GetMass());

    double ideal_rapidity = iE <= ipz ? -2.1 : 0.5 * std::log((iE + ipz) / (iE - ipz));
    fMC_rap->Fill(rapidity, ideal_rapidity);
  }

  fYPt->Fill(rapidity, pt);
  fPtPhi->Fill(phi, pt);
  fEtaPhi->Fill(phi, eta);

#if PIDPROB_ENABLED
  fPidProbPion->Fill(p, track->PidProbPion());
  fPidProbKaon->Fill(p, track->PidProbKaon());
  fPidProbProton->Fill(p, track->PidProbProton());
  fPidProbElectron->Fill(p, track->PidProbElectron());
#endif

  // fTofVsP->Fill(p, track->TOFdeuteronTime());
  fTofVsP->Fill(p, track->TOFsignal());
  fTofPionVsP->Fill(p, track->TOFpionTime());
  fTofKaonVsP->Fill(p, track->TOFkaonTime());
  fTofProtonVsP->Fill(p, track->TOFprotonTime());

  fTofMass->Fill(track->MassTOF());

  fdEdX->Fill(p, track->TPCsignal());
  fTpcTofPionSigma->Fill(track->NSigmaTPCPi(), track->NSigmaTOFPi());
  fTpcTofKaonSigma->Fill(track->NSigmaTPCK(), track->NSigmaTOFK());

  fNsigPionTpc->Fill(p, track->NSigmaTPCPi());
  fNsigKaonTpc->Fill(p, track->NSigmaTPCK());
  fNsigProtonTpc->Fill(p, track->NSigmaTPCP());

  fNsigPionTof->Fill(p, track->NSigmaTOFPi());
  fNsigKaonTof->Fill(p, track->NSigmaTOFK());
  fNsigProtonTof->Fill(p, track->NSigmaTOFP());

  fChi2Tpc->Fill(TPC_ncls > 0 ? track->TPCchi2() / TPC_ncls : -1.0);
  fChi2Its->Fill(ITS_ncls > 0 ? track->ITSchi2() / ITS_ncls : -1.0);

  fChiTpcIts->Fill(track->TPCchi2perNDF(), track->ITSchi2perNDF());

  fClsTpcIts->Fill(track->TPCncls(), track->ITSncls());

  fImpact->Fill(track->ImpactZ(), track->ImpactD());
  fEtaY->Fill(eta, rapidity);
}


AliFemtoCutMonitorPionPion::Pair::Pair(const bool passing,
                                       const TString& typestr,
                                       const bool is_mc_analysis,
                                       const bool suffix_output):
  AliFemtoCutMonitor()
  , fCurrentMagneticField(0.500668)
  , fRadius(1.2)
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

  auto hist_name = [&] (const TString &name) { return name + pf; };
  auto hist_title = [&] (const char *title, const char *axes) {
    return TString::Format(title_format, title, axes);
  };

  fMinv = new TH1F(
    hist_name("Pair_Minv"),
    hist_title("M_{inv}",
               "M_{inv} (GeV); N_{pairs}"),
    288, 0.0, 8.0);

  fKt = new TH1F(
    hist_name("kt"),
    hist_title("k_{T} Distribution",
               "k_{T} (GeV); N_{pairs}"),
    144, 0.0, 4.0);

  fDetaDphi = new TH2F(
    hist_name("DetaDphi"),
    hist_title("#Delta #eta vs #Delta #phi*",
               "#Delta #eta; #Delta #phi*"),
    145, -0.2, 0.2,
    145, -0.2, 0.2);

  fQinvDeta = new TH2F(
    hist_name("QinvDeta"),
    hist_title("Q_{inv} vs #Delta #eta",
               "Q_{inv} (GeV); #Delta #eta"),
    100, 0.0, 1.2,
    75, -0.1, 0.1);

  fQinvDphiStar = new TH2F(
    hist_name("QinvDphiStar"),
    hist_title("Q_{inv} vs #Delta #phi*",
               "Q_{inv} (GeV); #Delta #phi*"),
    100, 0.0, 1.2,
    75, -0.1, 0.1);

  if (is_mc_analysis) {
    fMCTrue_minv = new TH2F(
      hist_name("mc_Minv"),
      hist_title("Minv True vs Reconstructed",
                 "M_{inv}^{gen} (Gev);"
                 "M_{inv}^{rec} (GeV);"
                 ),
      144, 0.0, 4.5,
      144, 0.0, 4.5);

    fMCTrue_qinv = new TH2F(
      hist_name("mc_Qinv"),
      hist_title("q_{inv} True vs Reconstructed",
                 "q_{inv}^{gen} (Gev);"
                 "q_{inv}^{rec} (GeV);"
                 ),
      400, 0.0, 1.0,
      400, 0.0, 1.0);
  }
}

void
AliFemtoCutMonitorPionPion::Pair::Fill(const AliFemtoPair *pair)
{
  const float minv = pair->MInv(),
              // kstar = pair->KStar(),
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
         delta_phi_star = AliFemtoPairCutDetaDphi::CalculateDPhiStar(p1, chg_1, p2, chg_2, fRadius, fCurrentMagneticField);

  // if (fabs(delta_eta) <= 0.07 && fabs(delta_phi_star) <= 0.07 && passes)
  //     printf(">> %f % 6f % 6f %p \n", pair, delta_eta, delta_phi_star, this);
    // std::cout << ">> " << pair << " " << delta_eta << ", " << delta_phi_star << "\n"; // -> " << passes << "\n";

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

    const double true_qinv = -1.0 * (p1 - p2).m(),
                 true_minv = (p1 + p2).m();

    // skip mis-identified protons
    if (mc_1->GetPDGPid() == 2212 || mc_2->GetPDGPid() == 2212) {
      return;
    }
    if (mc_1->GetPDGPid() != 211 || mc_2->GetPDGPid() != 211) {
      return;
    }

    fMCTrue_qinv->Fill(true_qinv, qinv);
    fMCTrue_minv->Fill(true_minv, minv);

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
