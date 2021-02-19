#ifndef ALIVCLUSTER_H_
#define ALIVCLUSTER_H_

#include <TLorentzVector.h>

#define NCELL 5

namespace {

    void cell_neighbor(Int_t n_m_m[], const unsigned int n,
                       const unsigned int m = 5, unsigned int ld = 0)
    {
        if (ld == 0) {
            ld = m;
        }

        const unsigned int sm = n < 11520 ? n / 1152 :
            n < 12288 ? 10 + (n - 11520) / 384 :
            n < 16896 ? 12 + (n - 12288) / 768 :
            18 + (n - 16896) / 384;
        const unsigned int nphi =
            sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;

        for (unsigned int i = 0; i < m; i++) {
            const unsigned int i_centered = i - (m - 1) / 2;
            const unsigned int offset_i = (i_centered & 1) == 0 ?
                i_centered * nphi :
                (n & 1) == 0 ?
                ((i_centered + 2) & ~1) * nphi + 1 :
                (i_centered & ~1) * nphi - 1;
            for (unsigned int j = 0; j < m; j++) {
                const unsigned int j_times_2_centered =
                    2 * j - (m - 1);

                n_m_m[i * ld + j] =
                    n + offset_i + j_times_2_centered;
            }
        }
    }

    class AliVCluster {
    public:
        bool _initialized;
        Int_t _cell_id[NCELL * NCELL];
        TLorentzVector _momentum;
        AliVCluster(void)
        {
            cell_neighbor(_cell_id, 0, NCELL);
            cell_neighbor(_cell_id, _cell_id[NCELL * NCELL - 1], NCELL);
        }
        const Int_t *GetCellsAbsId(void) const
        {
            return _cell_id;
        }
        Int_t GetNCells(void) const
        {
            return NCELL * NCELL;
        }
        void GetMomentum(TLorentzVector &p, const double *) const
        {
            p = _momentum;
        }
    };

}

#endif // ALIVCLUSTER_H_
