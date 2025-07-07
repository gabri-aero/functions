/**
 * @file Nlm.hpp
 *
 * @brief Header file for computing normalization constants for spherical
 * harmonics
 *
 * @author Gabriel Valles
 * @date 2025-02-17
 */
#ifndef _NLM_HPP_
#define _NLM_HPP_

#include <cmath>

/**
 * @file Nlm.hpp
 * @brief Normalization constants.
 */

/**
 * @class Nlm
 *
 * @brief Computes recursively and stores normalization constants for
 * fully-normalized spherical harmonics
 *
 * This class computes the normalization constants for fully-normalized
 * spherical harmonics as defined in Heiskanen and Moritz (1967, eq. 1-91). As a
 * result, the orthogonality integrals values are:
 * \f[
 * \frac{1}{4\pi} \iint \bar{Y}_{lm}(\phi,\theta) \bar{Y}_{l'm'}(\phi,\theta) \,
 * d\sigma =
 * \delta_{ll'}\delta_{mm'}
 * \f]
 * The normalization constants take the following form:
 * \f[
 * N_{lm} = \sqrt{\frac{(2-\delta_{0m})(2l+1)(l+m)!}{(l-m)!}}
 * \f]
 * This class leverages recursive relations to compute all the normalization
 * constants up to a maximum input degree minimizing the overflow problem.
 */
class Nlm {
    double *_Nlm = nullptr; // Private attribute storing Nlm coefficients
    int l_max;

    /**
     * Function that computes global index for internal data structure.
     * @param l Degree
     * @param m Order
     * @return Global index
     */
    int lm_idx(int l, int m) { return (l * (l + 1)) / 2 + m; };

  public:
    /**
     * Default constructor
     */
    Nlm() {};
    /**
     * Class constructor
     * @param l_max Maximum degree and order to which the normalization
     * constants are computed
     */
    Nlm(int l_max) : l_max(l_max) {
        this->_Nlm = new double[((l_max + 1) * (l_max + 2)) / 2];
        for (int l = 0; l <= l_max; l++) {
            // Compute for m = 0
            _Nlm[lm_idx(l, 0)] = sqrt(2 * l + 1);
        }
        for (int m = 1; m <= l_max; m++) {
            for (int l = m; l <= l_max; l++) {
                // Recursion for m > 0
                _Nlm[lm_idx(l, m)] = _Nlm[lm_idx(l, m - 1)] *
                                     sqrt(1.0 / ((l - m + 1) * (l + m)));
            }
        }
        // Adjust k for m>0
        for (int m = 1; m <= l_max; m++) {
            for (int l = m; l <= l_max; l++) {
                _Nlm[lm_idx(l, m)] *= sqrt(2);
            }
        }
    };

    // Copy assignment operator constructor
    Nlm &operator=(const Nlm &other) {
        if (this != &other) {
            l_max = other.l_max;
            // Allocate and assign Nlm
            if (other._Nlm) {
                int Nlm_size = (l_max + 1) * (l_max + 2) / 2;
                _Nlm = new double[Nlm_size];
                std::copy(other._Nlm, other._Nlm + Nlm_size, _Nlm);
            }
        }
        return *this;
    };

    // Copy constructor
    Nlm(const Nlm &other) : l_max(other.l_max) {
        // Allocate and assign Nlm
        int Nlm_size = (l_max + 1) * (l_max + 2) / 2;
        _Nlm = new double[Nlm_size];
        std::copy(other._Nlm, other._Nlm + Nlm_size, _Nlm);
    };

    // Destructor
    ~Nlm() {
        if (_Nlm)
            delete[] _Nlm;
    };

    /**
     * @brief Getter for normalization constant
     * @param l degree
     * @param m order
     */
    double get_Nlm(int l, int m) { return _Nlm[lm_idx(l, m)]; };
};

#endif //_NLM_HPP_
