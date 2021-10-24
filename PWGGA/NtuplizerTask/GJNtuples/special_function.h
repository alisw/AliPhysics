// -*- mode: c++; -*-

#ifndef SPECIAL_FUNCTION_H_
#define SPECIAL_FUNCTION_H_

namespace {

    void kahan_sum(double &s, double &s2, const double x)
    {
        // The Kahan summation algorithm. See e.g. N. J. Higham, "The
        // accuracy of floating point summation", SIAM J. Sci. Comput.
        // 14(4), 783--799 (1993), https://doi.org/10.1137/0914050 .
        s2 += x;

        const double t = s + s2;

        s2 += s - t;
        s = t;
    }

    double fast2sum(double &z, const double x, const double y)
    {
        // Fast2Sum, see T. J. Dekker, "A floating-point technique for
        // extending the available precision", Numer. Math. 18(3),
        // 224--242 (1971), https://doi.org/10.1007/BF01397083 , p.
        // 227, eq. (4.1) and p. 228, (4.3)
        z = x + y;

        return y - (z - x);
    }

    double angular_range_reduce(const double x)
    {
        if (!std::isfinite(x)) {
            return x;
        }

        // Cody and Waite range reduction with three constants. See
        // e.g. J.-M. Muller, Elementary Functions: Algorithms and
        // Implementation, 3rd Ed. (Springer, New York, 2016),
        // https://doi.org/10.1007/978-1-4899-7983-4 , section 11.2.1,
        // p. 203.

        // The worst case for small IEEE 754 double precision
        // arguments (that is until 6381956970095103 2^799, compare
        // table 11.2, p. 211 in Muller 2016) is 3205513981387887
        // 2^(-44) = 182.2..., which requires 59 additional digits for
        // accurate reduction. And 53 - ceil(59 / 2) - 1 = 22,
        // limiting the possible magnitude to 2^23 pi = 1.317... 10^7
        static const double cody_waite_x_max_3 = ldexp(M_PI, 23);
        static const double inverse_two_pi = 1 / (2 * M_PI);
        double ret;

        if (x >= -cody_waite_x_max_3 && x <= cody_waite_x_max_3) {
            // 1686629713 2^(-28)
            static const double two_pi_1 =  6.28318530693650245667;
            // 560513589 2^(-61)
            static const double two_pi_2 =  2.43084020360578856312e-10;
            // 2 pi - 14488038916154245685 2^(-61)
            static const double two_pi_3 = -1.00331152253366640471e-19;
            const double k = rint(x * inverse_two_pi);

            ret = x - k * two_pi_1;

            const double t = fast2sum(ret, ret, -k * two_pi_2);

            ret += t - k * two_pi_3;
        }
        else {
            // Fallback using the system trigonometric function

            long double sin_x;
            long double cos_x;
            
            sin_x = sinl(x);
            cos_x = cosl(x);
            ret = (double)atan2l(sin_x, cos_x);
        }

        // In all cases, avoid rounding to a duplicate -pi

        if(ret == -M_PI) {
            ret = M_PI;
        }

        return ret;
    }

}

#endif // SPECIAL_FUNCTION_H_
