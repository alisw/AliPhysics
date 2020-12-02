// -*- mode: c++; -*-

#ifndef EMCAL_H_
#define EMCAL_H_

#include <vector>
#include <set>
#include <AliVCluster.h>
#include <AliVCaloCells.h>
#include <AliVVZERO.h>
#include "special_function.h"
#include "keras_model.h"
#include "emcal_cell.h"

namespace {

    void cell_max_cross(Int_t &cell_id_max,
                        Double_t &cell_energy_max,
                        Double_t &energy_cross,
                        AliVCluster *c, AliVCaloCells *emcal_cell)
    {
        std::set<Int_t> cluster_cell_id;

        cell_id_max = -1;
        cell_energy_max = -INFINITY;

        for (Int_t j = 0; j < c->GetNCells(); j++) {
            const Int_t cell_id = c->GetCellsAbsId()[j];
            const Double_t cell_energy =
                emcal_cell->GetCellAmplitude(cell_id);

            cluster_cell_id.insert(cell_id);
            if (cell_energy > cell_energy_max) {
                cell_energy_max = cell_energy;
                cell_id_max = cell_id;
            }
        }

        energy_cross = NAN;

        if (cell_id_max != -1) {
            unsigned int cell_id_max_cross[4];

            cell_cross(cell_id_max_cross, cell_id_max);
            energy_cross = 0;
            for (size_t j = 0; j < 4; j++) {
                if (cluster_cell_id.find(cell_id_max_cross[j]) !=
                    cluster_cell_id.end()) {
                    energy_cross += emcal_cell->
                        GetCellAmplitude(cell_id_max_cross[j]);
                }
            }
        }
    }

    bool cell_masked(AliVCluster *c, std::vector<bool> emcal_mask)
    {
        for (Int_t j = 0; j < c->GetNCells(); j++) {
            const Int_t cell_id = c->GetCellsAbsId()[j];

            if (cell_id >= 0 &&
                cell_id < static_cast<Int_t>(emcal_mask.size()) &&
                !emcal_mask[cell_id]) {
                return true;
            }
        }

        return false;
    }

    std::vector<float>
    cluster_cell_keras_inference(AliVCluster *cluster,
                                 AliVCaloCells *cell,
                                 const double *vertex,
                                 const AliVVZERO *v0,
                                 KerasModel &model)
    {
        Int_t cell_id_max = -1;
        Double_t cell_energy_max = -INFINITY;
        Double_t cell_cross = NAN;

        cell_max_cross(cell_id_max, cell_energy_max, cell_cross,
                       cluster, cell);

        unsigned int sm_max;
        unsigned int nphi_max;

        to_sm_nphi(sm_max, nphi_max, cell_id_max);

#if 0
        if (!inside_edge(cell_id_max, 2)) {
            return std::vector<float>();
        }
#endif

        TLorentzVector p;

        cluster->GetMomentum(p, vertex);

        unsigned int cell_id_5_5[25];

        cell_5_5(cell_id_5_5, cell_id_max);

        Tensor feature(29);

        for (size_t i = 0; i < 25; i++) {
            unsigned int sm;
            unsigned int nphi;

            to_sm_nphi(sm, nphi, cell_id_5_5[i]);
            feature(i) = sm == sm_max ?
                cell->GetCellAmplitude(cell_id_5_5[i]) / p.E() :
                0;
            if (std::isnan(feature(i))) {
                feature(i) = 0;
            }
        }
        feature(25) = 1 / sqrt(p.E());
        feature(26) = p.Eta();

        static const double sm_center[20][2] = {
            {  1,  1.5707958333333356 },
            { -1,  1.570796701388894  },
            {  1,  1.9198628472222254 },
            { -1,  1.9198605902777761 },
            {  1,  2.2689267361111165 },
            { -1,  2.2689289930555647 },
            {  1,  2.6179921875000085 },
            { -1,  2.617995138888896  },
            {  1,  2.9670605902777836 },
            { -1,  2.9670579861111235 },
            {  1,  3.2079436405129194 },
            { -1,  3.2083967655129197 },
            {  1, -1.5707963541666699 },
            { -1, -1.5707968750000028 },
            {  1, -1.2217302083333352 },
            { -1, -1.2217322916666704 },
            {  1, -0.8726658854166689 },
            { -1, -0.8726617187500013 },
            {  1, -0.6317869791666666 },
            { -1, -0.6313296875       }
        };
        double min_dazimuth = M_PI;

        for (size_t i = 0; i < 20; i++) {
            if (sm_center[i][0] * p.Eta() >= 0) {
                const double dazimuth = angular_range_reduce(
                    p.Phi() - sm_center[i][1]);
                if (fabs(dazimuth) < min_dazimuth) {
                    min_dazimuth = dazimuth;
                }
            }
        }
        feature(27) = min_dazimuth;

        double sum_v0 = 0;

        for (size_t i = 0; i < 64; i++) {
            sum_v0 += v0->GetMultiplicity(i);
        }
        feature(28) = sum_v0;

#if 0
        fprintf(stderr, "%s:%d: %g %g %g\n", __FILE__, __LINE__,
                feature(25), feature(27), feature(28));
#endif

static const float mean[29] = {2.12970655e-03, 3.88862821e-03, 4.92862100e-03, 3.90173425e-03, 2.08488060e-03, 4.72605834e-03, 1.24210687e-02, 4.01937366e-02, 1.25719616e-02, 4.69736801e-03, 9.08698235e-03, 5.94846457e-02, 6.88571036e-01, 6.02749288e-02, 9.03223641e-03, 4.62681567e-03, 1.22468071e-02, 3.87016498e-02, 1.23087661e-02, 4.65176860e-03, 2.12111045e-03, 3.95621359e-03, 4.97721089e-03, 3.94566311e-03, 2.14471412e-03, 3.06629092e-01, -8.71083885e-03, 1.00243546e-03, 3.31128784e+02};
static const float std[29] = {1.60925817e-02, 2.37639789e-02, 2.78392453e-02, 2.37607714e-02, 1.58410259e-02, 2.56051980e-02, 3.41129899e-02, 6.61690012e-02, 3.42592224e-02, 2.55473647e-02, 3.53715308e-02, 8.71860087e-02, 1.73532262e-01, 8.74797031e-02, 3.50968055e-02, 2.52052099e-02, 3.41960452e-02, 6.45336062e-02, 3.40170525e-02, 2.53446139e-02, 1.57360733e-02, 2.41008513e-02, 2.80672945e-02, 2.39805095e-02, 1.62536539e-02, 3.54581177e-02, 3.99770141e-01, 7.72248954e-02, 1.65971893e+02};

        for (size_t i = 0; i < 29; i++) {
            feature(i) = (feature(i) - mean[i]) / std[i];
        }

        Tensor output_tensor;

        model.Apply(&feature, &output_tensor);

        return output_tensor.data_;
    }

}

#endif // EMCAL_H_
