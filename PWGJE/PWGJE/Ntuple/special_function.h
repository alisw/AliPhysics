// -*- mode: c++; -*-

#ifndef SPECIAL_FUNCTION_H_
#define SPECIAL_FUNCTION_H_

namespace {

    void kahan_sum(double &s, double &s2, const double x)
    {
        s2 += x;

        const double t = s + s2;

        s2 += s - t;
        s = t;
    }

    double angular_range_reduce(const double x)
    {
        if (!std::isfinite(x)) {
            return x;
        }

        static const double cody_waite_x_max = 1608.4954386379741381;
        static const double two_pi_0 = 6.2831853071795649157;
        static const double two_pi_1 = 2.1561211432631314669e-14;
        static const double two_pi_2 = 1.1615423895917441336e-27;
        double ret;

        if(x >= -cody_waite_x_max && x <= cody_waite_x_max) {
            static const double inverse_two_pi =
                0.15915494309189534197;
            const double k = rint(x * inverse_two_pi);
            ret = ((x - (k * two_pi_0)) - k * two_pi_1) -
                k * two_pi_2;
        }
        else {
            long double sin_x;
            long double cos_x;

            sincosl(x, &sin_x, &cos_x);
            ret = (double)atan2l(sin_x, cos_x);
        }
        if(ret == -M_PI) {
            ret = M_PI;
        }

        return ret;
    }

}

#endif // SPECIAL_FUNCTION_H_
