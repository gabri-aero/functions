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

#include "Nlm.hpp"

#include <cmath>

/**
 * @class Plm
 *
 * @brief Class that computes and stores the Associated Legendre Functions
 * (ALFs) and its derivatives at a given co-latitude.
 *
 * This class computes the Associated Legendre Functions (ALFs) up to a certain
 * degree and order employing a recursive standard forward column method,
 * so-called Fixed-Order-Increase-Degree (FOID), as described by Holmes and
 * Featherstone, 2002 (sec. 2.1).
 *
 * First the following constants are defined:
 * \f[
 * a_{lm} = \sqrt{\frac{(2l-1)(2l+1)}{(l-m)(l+m)}} \quad\quad\quad b_{lm} =
 * \sqrt{\frac{(2l+1)(l+m-1)(l-m+1)}{(l-m)(l+m)(2l-3)}}
 * \f]
 * We define that \f$ t=\cos{\theta}, u=\sin{\theta}\f$ and the algorithm is
 * initialised with \f$
 * \bar{P}_{0,0}(\theta)=1, \bar{P}_{1,1}(\theta)=\sqrt{3}u\f$. Then, sectorial
 * recursive relationships are applied:
 * \f[
 *  \bar{P}_{mm} = u \sqrt{\frac{2m+1}{2m}} \bar{P}_{m-1,m-1}(\theta)
 * \f]
 * Next, we sequentially fix the order \f$ m\f$ and increase the degree \f$ l\f$
 * applying the following recursive relationships:
 * \f[
 *  \bar{P}_{lm}(\theta) = a_{lm} t \bar{P}_{lm}(\theta) - b_{lm}
 * \bar{P}_{l-2,m}
 * \f]
 *
 * Derivatives with respect to the co-latitude are of interest for multiple
 * applications and can be computed from the ALFs with the following
 * expressions:
 * \f[
 *  f_{lm} = \sqrt{\frac{(l^2-m^2)(2l+1)}{2l-1}}
 * \f]
 * \f[
 *  \frac{d \bar{P}_{lm}(\theta)}{d\theta} = \frac{1}{u} (lt
 * \bar{P}_{lm}(\theta)-f_{lm}
 * \bar{P}_{l-1,m}(\theta))
 * \f]
 *
 *
 */
class Plm {
    int l_max;    // Maximum degree of ALFs
    Nlm _Nlm;     // Normalization constants
    double theta; // Co-latitude

    double *_Plm = nullptr;  // Fully-normalized ALFs
    double *_dPlm = nullptr; // Fully-normalized ALFs co-latitude derivatives
    double *_ddPlm =
        nullptr; // Fully-normalized ALFs co-latitude 2nd order derivatives

    /**
     * Function that computes global index for internal data structure.
     * @param l degree
     * @param m order
     * @return Global index
     */
    int lm_idx(int l, int m) { return (l * (l + 1)) / 2 + m; };

  public:
    /**
     * Default constructor
     */
    Plm() {};
    /**
     * Class constructor
     * @param l_max Maximum degree to which the ALFs (or its derivatives) are
     * computed. Note that to compute the derivatives up to a degree L, it is
     * necessary to compute the ALFs up to degree L+2. This is automatically
     * handled based on derivative flags.
     * @param theta Co-latitude at which the ALFs (and its derivatives) are
     * evaluated
     * @param derivatives Flag to indicate whether derivatives are computed or
     * not
     * @param second_derivatives Flag to indicate whether 2nd order derivatives
     * are computed or not
     */
    Plm(int l_max, double theta, bool derivatives = false,
        bool second_derivatives = false)
        : l_max(l_max), _Nlm(Nlm(l_max)), theta(theta) {
        // Allocate ALFs
        this->_Plm = new double[((l_max + 1) * (l_max + 2)) / 2];
        // Define constants for FOID recursion
        double *a = new double[((l_max + 1) * (l_max + 2)) / 2];
        double *b = new double[((l_max + 1) * (l_max + 2)) / 2];
        for (int l = 0; l <= l_max; l++) {
            for (int m = 0; m < l; m++) {
                a[lm_idx(l, m)] =
                    sqrt((2 * l - 1.0) * (2 * l + 1) / ((l - m) * (l + m)));
                b[lm_idx(l, m)] =
                    l - m != 1
                        ? sqrt(((2 * l + 1.0) * (l + m - 1) * (l - m - 1)) /
                               ((l - m) * (l + m) * (2 * l - 3)))
                        : 0;
            }
        }
        // Define cosine, sine
        double t = cos(theta);
        double u = sin(theta);
        // Define P00
        _Plm[0] = 1;
        // Define P11
        if (l_max > 0)
            _Plm[lm_idx(1, 1)] = sqrt(3) * u;
        // Recursion for sectorial polynomials
        for (int l = 2; l <= l_max; l++) {
            _Plm[lm_idx(l, l)] =
                sqrt((2 * l + 1.0) / (2 * l)) * u * _Plm[lm_idx(l - 1, l - 1)];
        }
        // Recursion for terms below diagonal
        for (int m = 0; m < l_max; m++) { // Fix order
            // Now increase degree
            int l = m + 1;
            // Terms right below the diagonal
            _Plm[lm_idx(l, m)] = a[lm_idx(l, m)] * t * _Plm[lm_idx(l - 1, m)];
            // Other terms
            for (int l = m + 2; l <= l_max; l++) {
                _Plm[lm_idx(l, m)] =
                    a[lm_idx(l, m)] * t * _Plm[lm_idx(l - 1, m)] -
                    b[lm_idx(l, m)] * _Plm[lm_idx(l - 2, m)];
            }
        }
        // Delete allocated data
        delete[] a;
        delete[] b;
        // Compute derivatives
        if (derivatives) {
            // Allocate derivatives
            _dPlm = new double[((l_max + 1) * (l_max + 2)) / 2];
            // Constant for derivatives recursion
            double *f = new double[((l_max + 1) * (l_max + 2)) / 2];
            for (int l = 0; l <= l_max; l++) {
                for (int m = 0; m <= l; m++) {
                    f[lm_idx(l, m)] =
                        sqrt((l * l - m * m) * (2 * l + 1) / (2 * l - 1.0));
                }
            }
            // Sectorial terms
            for (int m = 0; m <= l_max; m++) {
                _dPlm[lm_idx(m, m)] = m * t / u * _Plm[lm_idx(m, m)];
            }
            // Terms below diagonal
            for (int l = 1; l <= l_max; l++) {
                for (int m = 0; m < l; m++) {
                    _dPlm[lm_idx(l, m)] =
                        1.0 / u *
                        (l * t * _Plm[lm_idx(l, m)] -
                         f[lm_idx(l, m)] * _Plm[lm_idx(l - 1, m)]);
                }
            }

            // Compute 2nd order derivatives
            if (second_derivatives) {
                // Allocate 2nd order derivatives
                _ddPlm = new double[((l_max + 1) * (l_max + 2)) / 2];
                // Sectorial terms
                for (int m = 0; m <= l_max; m++) {
                    _ddPlm[lm_idx(m, m)] =
                        (m - 1) * t / u * _dPlm[lm_idx(m, m)] -
                        m * _Plm[lm_idx(m, m)];
                }
                // Terms below diagonal
                for (int l = 1; l <= l_max; l++) {
                    for (int m = 0; m < l; m++) {
                        _ddPlm[lm_idx(l, m)] =
                            1.0 / u *
                                ((l - 1) * t * _dPlm[lm_idx(l, m)] -
                                 f[lm_idx(l, m)] * _dPlm[lm_idx(l - 1, m)]) -
                            l * _Plm[lm_idx(l, m)];
                    }
                }
            } else {
                _ddPlm = nullptr;
            }
            // Delete allocated data
            delete[] f;
        } else {
            _dPlm = nullptr;
        }
    };

    // Copy constructor
    Plm(const Plm &other)
        : l_max(other.l_max), _Nlm(other._Nlm), theta(other.theta) {
        // Allocate and assign Plm
        int Plm_size = (l_max + 1) * (l_max + 2) / 2;
        _Plm = new double[Plm_size];
        std::copy(other._Plm, other._Plm + Plm_size, _Plm);
        // Allocate and assign derivatives
        if (other._dPlm) {
            _dPlm = new double[Plm_size];
            std::copy(other._dPlm, other._dPlm + Plm_size, _dPlm);
        } else {
            _dPlm = nullptr;
        }
        // Allocate and assign 2nd order derivatives
        if (other._ddPlm) {
            _ddPlm = new double[Plm_size];
            std::copy(other._ddPlm, other._ddPlm + Plm_size, _ddPlm);
        } else {
            _ddPlm = nullptr;
        }
    };

    // Copy assignment operator
    Plm &operator=(const Plm &other) {
        if (this != &other) {
            if (_Plm)
                delete[] _Plm;
            if (_dPlm)
                delete[] _dPlm;
            if (_ddPlm)
                delete[] _ddPlm;
            // Copy data
            theta = other.theta;
            _Nlm = other._Nlm;
            l_max = other.l_max;
            // Compute total size
            int Plm_size = (l_max + 1) * (l_max + 2) / 2;
            // Allocate and assign Plm
            if (other._Plm) {
                _Plm = new double[Plm_size];
                std::copy(other._Plm, other._Plm + Plm_size, _Plm);
            } else {
                _Plm = nullptr;
            }
            // Allocate and assign derivatives
            if (other._dPlm) {
                _dPlm = new double[Plm_size];
                std::copy(other._dPlm, other._dPlm + Plm_size, _dPlm);
            } else {
                _dPlm = nullptr;
            }
            // Allocate and assign 2nd order derivatives
            if (other._ddPlm) {
                _ddPlm = new double[Plm_size];
                std::copy(other._ddPlm, other._ddPlm + Plm_size, _ddPlm);
            } else {
                _ddPlm = nullptr;
            }
        }
        return *this;
    };

    // Destructor
    ~Plm() {
        if (_Plm)
            delete[] _Plm;
        if (_dPlm)
            delete[] _dPlm;
        if (_ddPlm)
            delete[] _ddPlm;
    };

    /**
     * @brief Getter for fully-normalized ALF
     * @param l degree
     * @param m order
     */
    double get_Plm_bar(int l, int m) { return _Plm[lm_idx(l, m)]; };

    /**
     * @brief Getter for unnormalized ALF
     * @param l degree
     * @param m order
     */
    double get_Plm(int l, int m) {
        return _Plm[lm_idx(l, m)] / _Nlm.get_Nlm(l, m);
    };

    /**
     * @brief Getter for fully-normalized ALF derivative
     * @param l degree
     * @param m order
     */
    double get_dPlm_bar(int l, int m) { return _dPlm[lm_idx(l, m)]; };

    /**
     * @brief Getter for unnormalized ALF derivative
     * @param l degree
     * @param m order
     */
    double get_dPlm(int l, int m) {
        return _dPlm[lm_idx(l, m)] / _Nlm.get_Nlm(l, m);
    };

    /**
     * @brief Getter for fully-normalized ALF 2nd order derivative
     * @param l degree
     * @param m order
     */
    double get_ddPlm_bar(int l, int m) { return _ddPlm[lm_idx(l, m)]; };

    /**
     * @brief Getter for unnormalized ALF derivative
     * @param l degree
     * @param m order
     */
    double get_ddPlm(int l, int m) {
        return _ddPlm[lm_idx(l, m)] / _Nlm.get_Nlm(l, m);
    };

    /**
     * @brief Getter for associated colatitude
     */
    double get_theta() const { return theta; };
};

#endif // _PLM_HPP_
