#include <TFile.h>
#include <TTree.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#include <TLorentzVector.h>
#pragma GCC diagnostic pop

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TProfile.h>

#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequenceArea.hh>

#ifdef FASTJET_MEDIAN
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/Subtractor.hh>
#endif // FASTJET_MEDIAN

#include <jet.h>

#define NTRACK_MAX (1U << 17)

int main(int argc, char *argv[])
{
    if (argc < 2) {
        exit(EXIT_FAILURE);
    }

    int dummyc = 1;
    char **dummyv = new char *[1];

    dummyv[0] = strdup("main");

    TApplication application("", &dummyc, dummyv);
    TCanvas canvas("canvas", "", 960 + 4, 720 + 28);
    TH1D histogram0("histogram0", "", 120, -2, 2);
    TH1D histogram1("histogram1", "", 120, -2, 2);

    canvas.SetLogy();
    histogram0.Sumw2();
    histogram0.SetMarkerStyle(20);
    histogram0.SetMarkerColor(EColor::kBlack);
    histogram0.SetLineColor(EColor::kBlack);
    histogram0.SetXTitle("\\Delta p_T\\:(\\mathrm{GeV}/c)");
    histogram1.Sumw2();
    histogram1.SetMarkerStyle(25);
    histogram1.SetMarkerColor(EColor::kBlack);
    histogram1.SetLineColor(EColor::kBlack);
    histogram1.SetXTitle("\\Delta p_T\\:(\\mathrm{GeV}/c)");

    for (int iarg = 1; iarg < argc; iarg++) {
        TFile *file = TFile::Open(argv[iarg]);

        if (file == NULL) {
            exit(EXIT_FAILURE);
        }

        TTree *_tree_event = dynamic_cast<TTree *>
            (dynamic_cast<TDirectoryFile *>
             (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));

        UInt_t ntrack;
        Float_t track_e[NTRACK_MAX];
        Float_t track_pt[NTRACK_MAX];
        Float_t track_eta[NTRACK_MAX];
        Float_t track_phi[NTRACK_MAX];
        UChar_t track_quality[NTRACK_MAX];

        _tree_event->SetBranchAddress("ntrack", &ntrack);
        _tree_event->SetBranchAddress("track_e", track_e);
        _tree_event->SetBranchAddress("track_pt", track_pt);
        _tree_event->SetBranchAddress("track_eta", track_eta);
        _tree_event->SetBranchAddress("track_phi", track_phi);
        _tree_event->SetBranchAddress("track_quality", track_quality);

#ifdef FASTJET_MEDIAN
        fastjet::JetMedianBackgroundEstimator fastjet_ue_estimator(
            fastjet::SelectorAbsRapMax(0.9),
            fastjet::JetDefinition(fastjet::kt_algorithm, 0.3),
            fastjet::AreaDefinition(
                fastjet::active_area,
                fastjet::GhostedAreaSpec(0.9)));
#endif // FASTJET_MEDIAN

        for (Long64_t i = 0; i < _tree_event->GetEntries(); i++) {
            _tree_event->GetEntry(i);

            std::vector<fastjet::PseudoJet> particle_reco;
            std::vector<point_2d_t> particle_reco_area_estimation;
            std::map<size_t, size_t> track_reco_index;
            double sum_pt = 0;

            for (UInt_t j = 0; j < ntrack; j++) {
                if ((track_quality[j] & 3) != 0 &&
                    std::isfinite(track_eta[j]) &&
                    std::isfinite(track_phi[j])) {
                    track_reco_index[j] = particle_reco.size();

                    TLorentzVector v;

                    v.SetPtEtaPhiE(track_pt[j], track_eta[j],
                                   track_phi[j], track_e[j]);
                    particle_reco.push_back(fastjet::PseudoJet(
                        v.Px(), v.Py(), v.Pz(), v.P()));
                    particle_reco_area_estimation.push_back(
                        point_2d_t(track_eta[j], track_phi[j]));
                    sum_pt += track_pt[j];
                }
            }

            for (size_t j = 0; j < particle_reco.size(); j++) {
                particle_reco[j].set_user_index(static_cast<int>(j));
            }

            std::vector<double> particle_reco_area;
            std::vector<std::set<size_t> > particle_reco_incident;

            voronoi_area_incident(particle_reco_area,
                                  particle_reco_incident,
                                  particle_reco_area_estimation);

            double sum_area = 0;

            for (size_t j = 0; j < particle_reco_area.size(); j++) {
                sum_area += particle_reco_area[j];
            }

            static const double jet_kt_d_ue_estimation = 0.3;
            const fastjet::ClusterSequenceArea
                cluster_sequence_ue_estimation(
                    particle_reco,
                    fastjet::JetDefinition(fastjet::JetDefinition(
                        fastjet::kt_algorithm, jet_kt_d_ue_estimation)),
                    fastjet::VoronoiAreaSpec());

            std::pair<std::vector<double>, std::vector<double> >
                ue_estimate = ue_estimation_truncated_mean(
                    cluster_sequence_ue_estimation, particle_reco_area);

            const fastjet::ClusterSequenceArea
                cluster_sequence(
                    particle_reco,
                    fastjet::JetDefinition(fastjet::JetDefinition(
                        fastjet::antikt_algorithm, 0.4)),
                    fastjet::AreaDefinition(
                        fastjet::active_area,
                        fastjet::GhostedAreaSpec(0.9)));
            const std::vector<fastjet::PseudoJet> jet_with_ue =
                cluster_sequence.inclusive_jets(0);

#ifdef FASTJET_MEDIAN
            fastjet_ue_estimator.set_particles(particle_reco);
#endif // FASTJET_MEDIAN

            for (std::vector<fastjet::PseudoJet>::const_iterator
                     iterator_jet = jet_with_ue.begin();
                 iterator_jet != jet_with_ue.end(); iterator_jet++) {
                const std::vector<fastjet::PseudoJet> constituent =
                    cluster_sequence.constituents(*iterator_jet);
                double pt_raw_ue = 0;

                for (std::vector<fastjet::PseudoJet>::const_iterator
                         iterator_constituent = constituent.begin();
                     iterator_constituent != constituent.end();
                     iterator_constituent++) {
                    const int index = iterator_constituent->user_index();

                    if (index >= 0 && static_cast<size_t>(index) <
                        particle_reco_area.size()) {
                        pt_raw_ue += evaluate_ue(
                            ue_estimate,
                            iterator_constituent->pseudorapidity(),
                            iterator_constituent->phi_std()) *
                            particle_reco_area[index];
                    }
                }
                if (fabs(iterator_jet->pseudorapidity()) < 0.4) {
                    if (std::isfinite(pt_raw_ue)) {
                        histogram0.Fill(iterator_jet->perp() - pt_raw_ue);
                    }
#ifdef FASTJET_MEDIAN
                    histogram1.Fill(iterator_jet->perp() -
                                    fastjet_ue_estimator.rho(*iterator_jet) *
                                    iterator_jet->area());
#endif // FASTJET_MEDIAN
                }
            }

            if (i % 100 == 0) {
                histogram0.Draw("e1x0");
                histogram1.Draw("e1x0same");
                canvas.Update();
                canvas.SaveAs("sub.eps");
            }
        }
    }

    histogram0.Draw("e1x0");
    histogram1.Draw("e1x0same");
    canvas.Update();

    application.Run();

    return EXIT_SUCCESS;
}
