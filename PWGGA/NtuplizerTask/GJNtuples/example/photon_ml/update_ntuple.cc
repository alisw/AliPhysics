#include <numeric>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TProfile.h>

#include <H5Cpp.h>

#include <emcal.h>

#define RANK 2
#define NCLUSTER_MAX        (1U << 17)

int main(int argc, char *argv[])
{
    if (argc < 3) {
        exit(EXIT_FAILURE);
    }

    KerasModel model;

    model.LoadModel("../../photon_discr.model");

    TFile *root_file_in = TFile::Open(argv[1]);

    if (root_file_in == NULL) {
        return EXIT_FAILURE;
    }

    TDirectoryFile *df_in = dynamic_cast<TDirectoryFile *>
        (root_file_in->Get("AliAnalysisTaskNTGJ"));

    if (df_in == NULL) {
        return EXIT_FAILURE;
    }

    TTree *tree_event_in = dynamic_cast<TTree *>
        (df_in->Get("_tree_event"));

    if (tree_event_in == NULL) {
        return EXIT_FAILURE;
    }

    Float_t multiplicity_v0[64];
    UInt_t ncluster;
    Float_t cluster_e[NCLUSTER_MAX];
    Float_t cluster_eta[NCLUSTER_MAX];
    Float_t cluster_phi[NCLUSTER_MAX];
    Float_t cluster_s_nphoton[NCLUSTER_MAX][4];
    UShort_t cluster_cell_id_max[NCLUSTER_MAX];
    Float_t cell_e[17664];

    tree_event_in->SetBranchAddress("multiplicity_v0", multiplicity_v0);
    tree_event_in->SetBranchAddress("ncluster", &ncluster);
    tree_event_in->SetBranchAddress("cluster_e", cluster_e);
    tree_event_in->SetBranchAddress("cluster_eta", cluster_eta);
    tree_event_in->SetBranchAddress("cluster_phi", cluster_phi);
    tree_event_in->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton);
    tree_event_in->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
    tree_event_in->SetBranchAddress("cell_e", cell_e);

    TFile *root_file_out = TFile::Open(argv[2], "recreate");

    root_file_out->mkdir("AliAnalysisTaskNTGJ");
    root_file_out->cd("AliAnalysisTaskNTGJ");

    TTree *tree_event_out = tree_event_in->CloneTree(0);

    for (Long64_t i = 0; i < tree_event_in->GetEntries(); i++) {
        tree_event_in->GetEntry(i);
        for (UInt_t j = 0; j < ncluster; j++) {
            const double pt = cluster_e[j] / cosh(cluster_eta[j]);
            AliVCluster cluster;
            AliVCaloCells cell;
            AliVVZERO v0;
            const double vertex[3] = { 0 };

            cluster._momentum.SetPtEtaPhiM(pt, cluster_eta[j],
                                           cluster_phi[j], 0);

            unsigned int sm_max;
            unsigned int nphi_max;

            to_sm_nphi(sm_max, nphi_max, cluster_cell_id_max[j]);

            unsigned int cell_id_5_5[25];

            cell_5_5(cell_id_5_5, cluster_cell_id_max[j]);

            for (Int_t k = 0; k < 25; k++) {
                unsigned int sm;
                unsigned int nphi;

                to_sm_nphi(sm, nphi, cell_id_5_5[k]);
                cell._amplitude[cluster.GetCellsAbsId()[k]] =
                    sm == sm_max ?
                    cell_e[cell_id_5_5[k]] : 0;
            }
            v0._multiplicity_sum = 0;
            for (size_t k = 0; k < 64; k++) {
                v0._multiplicity_sum += multiplicity_v0[k];
            }

            std::vector<float> activation =
                cluster_cell_keras_inference(&cluster, &cell, vertex,
                                             &v0, model);

#if 0
            if (cluster_e[j] >= 8) {
                fprintf(stderr, "%s:%d: %f %f %f %f %f\n", __FILE__,
                        __LINE__, cluster_e[j], activation[0],
                        activation[1], cluster_s_nphoton[j][1],
                        cluster_s_nphoton[j][2]);
            }
#endif
            std::copy(activation.begin(), activation.end(),
                      cluster_s_nphoton[j] + 1);
        }
        tree_event_out->Fill();
    }
    tree_event_out->AutoSave();

    delete root_file_in;
    delete root_file_out;

    return EXIT_SUCCESS;
}
