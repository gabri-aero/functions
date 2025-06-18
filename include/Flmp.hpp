#ifndef _FLMP_HPP_
#define _FLMP_HPP_

/**
 * @file Flmp.hpp
 * @brief Inclination functions.
 */

#include <cmath>

#include <Eigen/Dense>
#include <chrono>
#include <complex>
#include <iostream>
#include "Fourier.hpp"
#include "Plm.hpp"

/**
 * @class Flmp
 *
 * @brief Class that computes and stores the normalized inclination functions and its
 * derivatives at a given inclination.
 *
 * This class computes the fully-normalized inclination functions and its derivatives applying a
 * real FFdouble to a disturbing potential along a great circle at the associated inclination without
 * approximation (Wagner, 1983). The same procedure can be followed for the derivative of the
 * disturbing potential w.r.t. the inclination in order to compute the derivatives of the
 * inclination function.
 *
 * Further details on the normalization can also be found in Nlm.hpp
 *
 * The class enables two different formulations found in literature, both \f$\bar{F}_{lmp}\f$ (e.g.
 * Kaula, 1966) and \f$\bar{F}_{lmk}\f$ with \f$k=l-2p\f$. The latter is more useful for gravity
 * field spectral analysis.
 *
 */
class Flmp
{
    int l_idx(int l) { return (l * (l + 1) * (2 * l + 1)) / 6; }
    double *_Flmp;
    double *_dFlmp;
    double I;
    int l_max;
    /**
     * Function that retrieves global index for a given l,m,p set
     *
     * @param l degree
     * @param m order
     * @param p p-index
     *
     * @return Global storing index associated to \f$\bar{F}_{lmp}\f$
     */
    int lmp_idx(int l, int m, int p) const;

    /**
     * Function that retrieves global index for a given l,m,k set
     *
     * @param l degree
     * @param m order
     * @param k k-index
     *
     * @return Global storing index associated to \f$\bar{F}_{lmp}\f$
     */
    int lmk_idx(int l, int m, int k) const;

public:
    /**
     * Class default constructor
     */
    Flmp();

    /**
     * Class constructor
     * @param l_max Maximum degree to which the inclination functions (and its derivatives) will be
     * computed
     * @param I Inclination at which the inclination functions (and its derivatives) are evaluated
     * @param compute_derivatives Flag to determine whether inclination functions derivatives shall
     * be computed or not
     */
    Flmp(int l_max, double I, bool compute_derivatives = false);

    /**
     * Copy assignment operator
     */
    Flmp &operator=(const Flmp &other);

    /**
     * Copy constructor
     */
    Flmp(const Flmp &other);

    /**
     * Destructor
     */
    ~Flmp();

    /**
     * Getter for maximum degree computed
     */
    int get_l_max() const;
    /**
     * Inclination function getter for l,m,p set
     * @param l Degree
     * @param m Order
     * @param p p-index
     * @return \f$\bar{F}_{lmp}\f$
     */
    double get_Flmp(int l, int m, int p) const;
    /**
     * Inclination function getter for l,m,k set
     * @param l Degree
     * @param m Order
     * @param k k-index
     * @return \f$\bar{F}_{lmk}\f$
     */
    double get_Flmk(int l, int m, int k) const;
    /**
     * Inclination function derivative getter for l,m,p set
     * @param l Degree
     * @param m Order
     * @param p k-index
     * @return \f$\bar{F}_{lmp}\f$
     */
    double get_dFlmp(int l, int m, int p) const;
    /**
     * Inclination function derivative getter for l,m,k set
     * @param l Degree
     * @param m Order
     * @param k k-index
     * @return \f$\bar{F}_{lmk}\f$
     */
    double get_dFlmk(int l, int m, int k) const;
    /**
     * Cross-track inclination function derivative getter for l,m,k set
     * @param l Degree
     * @param m Order
     * @param k k-index
     * @return \f$\bar{F}_{lmk}\f$
     */
    double get_Flmk_star(int l, int m, int k) const;
};

#endif //_FLMP_HPP_
