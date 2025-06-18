/**
 * @file Plm.hpp
 *
 * @brief Source file for defining Associated Legendre Functions (ALFs)
 *
 * @author Gabriel Valles
 * @date 2025-02-17
 */

#include "Plm.hpp"

int Plm::lm_idx(int l, int m) { return (l * (l + 1)) / 2 + m; }

Plm::Plm() {};

Plm::Plm(int l_max, double theta, bool derivatives, bool second_derivatives)
    : l_max(l_max), _Nlm(Nlm(l_max)), theta(theta)
{
    // Allocate ALFs
    this->_Plm = new double[((l_max + 1) * (l_max + 2)) / 2];
    // Define constants for FOID recursion
    double *a = new double[((l_max + 1) * (l_max + 2)) / 2];
    double *b = new double[((l_max + 1) * (l_max + 2)) / 2];
    for (int l = 0; l <= l_max; l++)
    {
        for (int m = 0; m < l; m++)
        {
            a[lm_idx(l, m)] = sqrt((2 * l - 1.0) * (2 * l + 1) / ((l - m) * (l + m)));
            b[lm_idx(l, m)] = l - m != 1 ? sqrt(((2 * l + 1.0) * (l + m - 1) * (l - m - 1)) /
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
    for (int l = 2; l <= l_max; l++)
    {
        _Plm[lm_idx(l, l)] = sqrt((2 * l + 1.0) / (2 * l)) * u * _Plm[lm_idx(l - 1, l - 1)];
    }
    // Recursion for terms below diagonal
    for (int m = 0; m < l_max; m++)
    { // Fix order
        // Now increase degree
        int l = m + 1;
        // Terms right below the diagonal
        _Plm[lm_idx(l, m)] = a[lm_idx(l, m)] * t * _Plm[lm_idx(l - 1, m)];
        // Other terms
        for (int l = m + 2; l <= l_max; l++)
        {
            _Plm[lm_idx(l, m)] = a[lm_idx(l, m)] * t * _Plm[lm_idx(l - 1, m)] -
                                 b[lm_idx(l, m)] * _Plm[lm_idx(l - 2, m)];
        }
    }
    // Delete allocated data
    delete[] a;
    delete[] b;
    // Compute derivatives
    if (derivatives)
    {
        // Allocate derivatives
        _dPlm = new double[((l_max + 1) * (l_max + 2)) / 2];
        // Constant for derivatives recursion
        double *f = new double[((l_max + 1) * (l_max + 2)) / 2];
        for (int l = 0; l <= l_max; l++)
        {
            for (int m = 0; m <= l; m++)
            {
                f[lm_idx(l, m)] = sqrt((l * l - m * m) * (2 * l + 1) / (2 * l - 1.0));
            }
        }
        // Sectorial terms
        for (int m = 0; m <= l_max; m++)
        {
            _dPlm[lm_idx(m, m)] = m * t / u * _Plm[lm_idx(m, m)];
        }
        // Terms below diagonal
        for (int l = 1; l <= l_max; l++)
        {
            for (int m = 0; m < l; m++)
            {
                _dPlm[lm_idx(l, m)] =
                    1.0 / u *
                    (l * t * _Plm[lm_idx(l, m)] - f[lm_idx(l, m)] * _Plm[lm_idx(l - 1, m)]);
            }
        }

        // Compute 2nd order derivatives
        if (second_derivatives)
        {
            // Allocate 2nd order derivatives
            _ddPlm = new double[((l_max + 1) * (l_max + 2)) / 2];
            // Sectorial terms
            for (int m = 0; m <= l_max; m++)
            {
                _ddPlm[lm_idx(m, m)] =
                    (m - 1) * t / u * _dPlm[lm_idx(m, m)] - m * _Plm[lm_idx(m, m)];
            }
            // Terms below diagonal
            for (int l = 1; l <= l_max; l++)
            {
                for (int m = 0; m < l; m++)
                {
                    _ddPlm[lm_idx(l, m)] = 1.0 / u *
                                               ((l - 1) * t * _dPlm[lm_idx(l, m)] -
                                                f[lm_idx(l, m)] * _dPlm[lm_idx(l - 1, m)]) -
                                           l * _Plm[lm_idx(l, m)];
                }
            }
        }
        else
        {
            _ddPlm = nullptr;
        }
        // Delete allocated data
        delete[] f;
    }
    else
    {
        _dPlm = nullptr;
    }
}

Plm::Plm(const Plm &other) : theta(other.theta), _Nlm(other._Nlm), l_max(l_max)
{
    // Allocate and assign Plm
    int Plm_size = (l_max + 1) * (l_max + 2) / 2;
    _Plm = new double[Plm_size];
    std::copy(other._Plm, other._Plm + Plm_size, _Plm);
    // Allocate and assign derivatives
    if (other._dPlm)
    {
        _dPlm = new double[Plm_size];
        std::copy(other._dPlm, other._dPlm + Plm_size, _dPlm);
    }
    else
    {
        _dPlm = nullptr;
    }
    // Allocate and assign 2nd order derivatives
    if (other._ddPlm)
    {
        _ddPlm = new double[Plm_size];
        std::copy(other._ddPlm, other._ddPlm + Plm_size, _ddPlm);
    }
    else
    {
        _ddPlm = nullptr;
    }
}

Plm &Plm::operator=(const Plm &other)
{
    if (this != &other)
    {
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
        if (other._Plm)
        {
            _Plm = new double[Plm_size];
            std::copy(other._Plm, other._Plm + Plm_size, _Plm);
        }
        else
        {
            _Plm = nullptr;
        }
        // Allocate and assign derivatives
        if (other._dPlm)
        {
            _dPlm = new double[Plm_size];
            std::copy(other._dPlm, other._dPlm + Plm_size, _dPlm);
        }
        else
        {
            _dPlm = nullptr;
        }
        // Allocate and assign 2nd order derivatives
        if (other._ddPlm)
        {
            _ddPlm = new double[Plm_size];
            std::copy(other._ddPlm, other._ddPlm + Plm_size, _ddPlm);
        }
        else
        {
            _ddPlm = nullptr;
        }
    }
    return *this;
}

Plm::~Plm()
{
    if (_Plm)
        delete[] _Plm;
    if (_dPlm)
        delete[] _dPlm;
    if (_ddPlm)
        delete[] _ddPlm;
}

double Plm::get_Plm_bar(int l, int m) { return _Plm[lm_idx(l, m)]; }

double Plm::get_Plm(int l, int m) { return _Plm[lm_idx(l, m)] / _Nlm.get_Nlm(l, m); }

double Plm::get_dPlm_bar(int l, int m) { return _dPlm[lm_idx(l, m)]; }

double Plm::get_dPlm(int l, int m) { return _dPlm[lm_idx(l, m)] / _Nlm.get_Nlm(l, m); }

double Plm::get_ddPlm_bar(int l, int m) { return _ddPlm[lm_idx(l, m)]; }

double Plm::get_ddPlm(int l, int m) { return _ddPlm[lm_idx(l, m)] / _Nlm.get_Nlm(l, m); }

double Plm::get_theta() const { return theta; }
