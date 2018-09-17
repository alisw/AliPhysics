// -*- mode: c++; -*-

#ifndef ISOLATION_H_
#define ISOLATION_H_

#include <algorithm>

namespace {

    // ROOT ACLiC weirdness is preventing delta_vs_iso from being
    // declared const

    double frixione_iso_max_x_e_eps(std::vector<std::pair<
                                    double, double> >
                                    delta_vs_iso,
                                    double delta_0, double n)
    {
        if (delta_vs_iso.empty()) {
            return 0;
        }

        std::sort(delta_vs_iso.begin(), delta_vs_iso.end());

        double sum_iso = 0;
        double max_x_e_eps = -INFINITY;
        const double one_cos_delta_0 = 1 - cos(delta_0);

        for (std::vector<std::pair<double, double> >::
                 const_iterator iterator = delta_vs_iso.begin();
             iterator != delta_vs_iso.end(); iterator++) {
            const double delta = iterator->first;

            sum_iso += iterator->second;

            const double x_e_eps = sum_iso *
                pow(one_cos_delta_0 / (1 - cos(delta)), n);

            max_x_e_eps = std::max(max_x_e_eps, x_e_eps);
        }

        return max_x_e_eps;
    }

    double anti_frixione_iso_max_x_e_eps(std::vector<std::pair<
                                         double, double> >
                                         delta_vs_iso,
                                         double delta_0, double n)
    {
        if (delta_vs_iso.empty()) {
            return 0;
        }

        std::sort(delta_vs_iso.begin(), delta_vs_iso.end());

        double sum_iso = 0;
        double min_x_e_eps = INFINITY;
        const double one_cos_delta_0 = 1 - cos(delta_0);

        for (std::vector<std::pair<double, double> >::
                 const_iterator iterator = delta_vs_iso.begin();
             iterator != delta_vs_iso.end(); iterator++) {
            const double delta = iterator->first;

            sum_iso += iterator->second;

            const double x_e_eps = sum_iso *
                pow(one_cos_delta_0 / (1 - cos(delta)), n);

            min_x_e_eps = std::min(min_x_e_eps, x_e_eps);
        }

        return min_x_e_eps;
    }

}

#endif // ISOLATION_H_
