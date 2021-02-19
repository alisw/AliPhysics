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

int main(int argc, char *argv[])
{
    if (argc < 2) {
        exit(EXIT_FAILURE);
    }

    KerasModel model;

    model.LoadModel("../../photon_discr.model");

    H5::H5File hdf5_file(argv[1], H5F_ACC_RDONLY);
    H5::DataSet data_set_X = hdf5_file.openDataSet("X");
    H5::DataSet data_set_y = hdf5_file.openDataSet("y");
	H5::DataSpace data_space_X = data_set_X.getSpace();
	H5::DataSpace data_space_y = data_set_y.getSpace();
    hsize_t dim_max_X[RANK];
    hsize_t dim_max_y[RANK];

    data_space_X.getSimpleExtentDims(dim_max_X, NULL);
    data_space_y.getSimpleExtentDims(dim_max_y, NULL);

    hsize_t dim_extend_X[RANK] = { 1, dim_max_X[1] };
    hsize_t dim_extend_y[RANK] = { 1, dim_max_y[1] };
    hsize_t offset[RANK] = {0, 0};

    fprintf(stderr, "%s:%d: %llu %llu\n", __FILE__, __LINE__,
            dim_max_X[0], dim_max_y[0]);

    std::vector<float> X(dim_max_X[1], NAN);
    std::vector<float> y(dim_max_y[1], NAN);
    H5::DataSpace mem_space_X(RANK, dim_extend_X);
    H5::DataSpace mem_space_y(RANK, dim_extend_y);

    // The convention in photon_dnn.py is that odd clusters are
    // test/validation ones
    for (hsize_t i = 1; i < std::min(dim_max_X[0], dim_max_y[0]);
         i += 2) {
        offset[0] = i;
        data_space_X.selectHyperslab(H5S_SELECT_SET, dim_extend_X,
                                     offset);
        data_space_y.selectHyperslab(H5S_SELECT_SET, dim_extend_y,
                                     offset);
        data_set_X.read(&X[0], H5::PredType::NATIVE_FLOAT,
                        mem_space_X, data_space_X);
        data_set_y.read(&y[0], H5::PredType::NATIVE_FLOAT,
                        mem_space_y, data_space_y);

        const double pt = 1 / (std::pow(X[25], 2) * cosh(X[26]));
        const double eta = X[26];
        const double phi = angular_range_reduce(
            (X[27] >= 0 ?
             (X[26] >= 0 ? 3.2079436405129194 : 3.2083967655129197) :
             (X[26] >= 0 ? 1.5707958333333356 : 1.570796701388894))
            + X[27]);
        AliVCluster cluster;
        AliVCaloCells cell;
        AliVVZERO v0;
        const double vertex[3] = { 0 };

        cluster._momentum.SetPtEtaPhiM(pt, eta, phi, 0);
        for (Int_t i = 0; i < cluster.GetNCells(); i++) {
            cell._amplitude[cluster.GetCellsAbsId()[i]] =
                X[i] / std::pow(X[25], 2);
        }
        v0._multiplicity_sum = X[28];

#if 0
        fprintf(stderr, "%s:%d: %g %g %g\n", __FILE__, __LINE__,
                X[25], X[27], X[28]);
#endif

        std::vector<float> activation =
            cluster_cell_keras_inference(&cluster, &cell, vertex,
                                         &v0, model);

        fprintf(stderr, "%s:%d: %g %g %g %g\n", __FILE__, __LINE__,
                1 / std::pow(X[25], 2), y[0],
                activation[0], activation[1]);
    }

    hdf5_file.close();

    return EXIT_SUCCESS;
}
