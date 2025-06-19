/**
 * @file Nlm.hpp
 *
 * @brief Header file for computing normalization constants for spherical harmonics
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
 * @brief Computes recursively and stores normalization constants for fully-normalized spherical
 * harmonics
 *
 * This class computes the normalization constants for fully-normalized spherical harmonics as
 * defined in Heiskanen and Moritz (1967, eq. 1-91). As a result, the orthogonality integrals values
 * are:
 * \f[
 * \frac{1}{4\pi} \iint \bar{Y}_{lm}(\phi,\theta) \bar{Y}_{l'm'}(\phi,\theta) \, d\sigma =
 * \delta_{ll'}\delta_{mm'}
 * \f]
 * The normalization constants take the following form:
 * \f[
 * N_{lm} = \sqrt{\frac{(2-\delta_{0m})(2l+1)(l+m)!}{(l-m)!}}
 * \f]
 * This class leverages recursive relations to compute all the normalization constants up to a
 * maximum input degree minimizing the overflow problem.
 */
class Nlm
{
    double *_Nlm = nullptr; // Private attribute storing Nlm coefficients
    int l_max;

    /**
     * Function that computes global index for internal data structure.
     * @param l Degree
     * @param m Order
     * @return Global index
     */
    int lm_idx(int l, int m);

public:
    /**
     * Default constructor
     */
    Nlm();
    /**
     * Class constructor
     * @param l_max Maximum degree and order to which the normalization constants are computed
     */
    Nlm(int l_max);

    // Copy assignment operator constructor
    Nlm &operator=(const Nlm &other);

    // Copy constructor
    Nlm(const Nlm &other);

    // Destructor
    ~Nlm();

    /**
     * @brief Getter for normalization constant
     * @param l degree
     * @param m order
     */
    double get_Nlm(int l, int m);
};

#endif //_NLM_HPP_
