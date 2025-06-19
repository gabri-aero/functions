/**
 * @file Plm.hpp
 *
 * @brief Header file to define Associated Legendre Functions (ALFs) class
 *
 * @author Gabriel Valles
 * @date 2025-02-17
 */
#ifndef _PLM_HPP_
#define _PLM_HPP_

#include <math.h>

#include "Nlm.hpp"

/**
 * @class Plm
 *
 * @brief Class that computes and stores the Associated Legendre Functions (ALFs) and its
 * derivatives at a given co-latitude.
 *
 * This class computes the Associated Legendre Functions (ALFs) up to a certain degree and order
 * employing a recursive standard forward column method, so-called Fixed-Order-Increase-Degree
 * (FOID), as described by Holmes and Featherstone, 2002 (sec. 2.1).
 *
 * First the following constants are defined:
 * \f[
 * a_{lm} = \sqrt{\frac{(2l-1)(2l+1)}{(l-m)(l+m)}} \quad\quad\quad b_{lm} =
 * \sqrt{\frac{(2l+1)(l+m-1)(l-m+1)}{(l-m)(l+m)(2l-3)}}
 * \f]
 * We define that \f$ t=\cos{\theta}, u=\sin{\theta}\f$ and the algorithm is initialised with \f$
 * \bar{P}_{0,0}(\theta)=1, \bar{P}_{1,1}(\theta)=\sqrt{3}u\f$. Then, sectorial recursive
 * relationships are applied:
 * \f[
 *  \bar{P}_{mm} = u \sqrt{\frac{2m+1}{2m}} \bar{P}_{m-1,m-1}(\theta)
 * \f]
 * Next, we sequentially fix the order \f$ m\f$ and increase the degree \f$ l\f$ applying the
 * following recursive relationships:
 * \f[
 *  \bar{P}_{lm}(\theta) = a_{lm} t \bar{P}_{lm}(\theta) - b_{lm} \bar{P}_{l-2,m}
 * \f]
 *
 * Derivatives with respect to the co-latitude are of interest for multiple applications and can be
 * computed from the ALFs with the following expressions:
 * \f[
 *  f_{lm} = \sqrt{\frac{(l^2-m^2)(2l+1)}{2l-1}}
 * \f]
 * \f[
 *  \frac{d \bar{P}_{lm}(\theta)}{d\theta} = \frac{1}{u} (lt \bar{P}_{lm}(\theta)-f_{lm}
 * \bar{P}_{l-1,m}(\theta))
 * \f]
 *
 *
 */
class Plm
{
    double *_Plm = nullptr;   // Fully-normalized ALFs
    double *_dPlm = nullptr;  // Fully-normalized ALFs co-latitude derivatives
    double *_ddPlm = nullptr; // Fully-normalized ALFs co-latitude 2nd order derivatives
    double theta;             // Co-latitude
    Nlm _Nlm;                 // Normalization constants
    int l_max;                // Maximum degree of ALFs

    /**
     * Function that computes global index for internal data structure.
     * @param l degree
     * @param m order
     * @return Global index
     */
    int lm_idx(int l, int m);

public:
    /**
     * Default constructor
     */
    Plm();
    /**
     * Class constructor
     * @tparam double Data type of ALFs and derivatives
     * @param l_max Maximum degree to which the ALFs (or its derivatives) are computed. Note that to
     * compute the derivatives up to a degree L, it is necessary to compute the ALFs up to degree
     * L+2.
     * @param theta Co-latitude at which the ALFs (and its derivatives) are evaluated
     * @param derivatives Flag to indicate whether derivatives are computed or not
     * @param second_derivatives Flag to indicate whether 2nd order derivatives are computed or not
     */
    Plm(int l_max, double theta, bool derivatives = false, bool second_derivatives = false);

    // Copy constructor
    Plm(const Plm &other);

    // Copy assignment operator
    Plm &operator=(const Plm &other);

    // Destructor
    ~Plm();

    /**
     * @brief Getter for fully-normalized ALF
     * @param l degree
     * @param m order
     */
    double get_Plm_bar(int l, int m);

    /**
     * @brief Getter for unnormalized ALF
     * @param l degree
     * @param m order
     */
    double get_Plm(int l, int m);

    /**
     * @brief Getter for fully-normalized ALF derivative
     * @param l degree
     * @param m order
     */
    double get_dPlm_bar(int l, int m);

    /**
     * @brief Getter for unnormalized ALF derivative
     * @param l degree
     * @param m order
     */
    double get_dPlm(int l, int m);

    /**
     * @brief Getter for fully-normalized ALF 2nd order derivative
     * @param l degree
     * @param m order
     */
    double get_ddPlm_bar(int l, int m);

    /**
     * @brief Getter for unnormalized ALF derivative
     * @param l degree
     * @param m order
     */
    double get_ddPlm(int l, int m);

    /**
     * @brief Getter for associated colatitude
     */
    double get_theta() const;
};

#endif // _PLM_HPP_
