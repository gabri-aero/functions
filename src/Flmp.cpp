/**
 * @file Flmp.cpp
 * @brief Source file to define inclination function class
 * @date 2025-02-17
 */

#include "Flmp.hpp"

int Flmp::lmp_idx(int l, int m, int p) const
{
    return (l + 1) * (l * (2 * l + 1)) / 6 + m * (l + 1) + p;
}

int Flmp::lmk_idx(int l, int m, int k) const { return lmp_idx(l, m, (l - k) / 2); }

Flmp::Flmp() : _Flmp(nullptr), _dFlmp(nullptr), l_max(0) {};

Flmp::Flmp(int l_max, double I, bool compute_derivatives) : l_max(l_max), I(I)
{
    // Allocate inclination functions
    _Flmp = new double[l_idx(l_max + 1)];
    // Compute inclination functions
    // Determine great circle sampling
    const int N = pow(2, ceil(log2(2 * l_max + 1))); // number of samples
    double du = 2 * M_PI / N;                        // step
    std::vector<double> u(N), lam(N), theta(N);
    Eigen::VectorX<double> Tlm(N);
    std::vector<Plm> plm;
    plm.reserve(N);
    double cos_I = cos(I);
    double sin_I = sin(I);
    std::vector<double> sin_u(N), cos_u(N);
    for (int i = 0; i < N; i++)
    {
        u[i] = du * i;
        sin_u[i] = sin(u[i]);
        cos_u[i] = cos(u[i]);
        lam[i] = atan2(cos_I * sin_u[i], cos_u[i]);
        theta[i] = acos(sin_I * sin_u[i]);
        plm.emplace_back(l_max, theta[i], compute_derivatives);
    }
    // Iterate over degree and order
    std::vector<double> C(l_max + 1), S(l_max + 1);
    Eigen::VectorX<std::complex<double>> y(N);
    int l, m, i, lm = 0, h = N;
    for (l = 0; l <= l_max; l++)
    {
        for (m = 0; m <= l; m++)
        {
            // Compute unit disturbing potential along great circle
            for (i = 0; i < N; i++)
            {
                Tlm[i] = plm[i].get_Plm_bar(l, m) * (cos(m * lam[i]) + sin(m * lam[i]));
            }
            // Analyse perturbing potential with FFT
            y = rfft(Tlm);
            for (i = 0; i <= l; i++)
            {
                C[i] = 2 * y[i].real() / N;
                S[i] = -2 * y[i].imag() / N;
            }
            // Map coefficients Ci, Si to Flmp
            if (l % 2 == 0)
            {
                _Flmp[lm + l / 2] = m % 2 == 0 ? C[0] : -C[0];
            }
            if (l % 2 == m % 2)
            {
                for (i = 0; i <= l; i++)
                {
                    if (i % 2 == l % 2)
                    {
                        _Flmp[lm + (l - i) / 2] = (C[i] + S[i]) / 2;
                        _Flmp[lm + (l + i) / 2] = (C[i] - S[i]) / 2;
                    }
                }
            }
            else
            {
                for (i = 0; i <= l; i++)
                {
                    if (i % 2 == l % 2)
                    {
                        _Flmp[lm + (l + i) / 2] = -(C[i] + S[i]) / 2;
                        _Flmp[lm + (l - i) / 2] = -(C[i] - S[i]) / 2;
                    }
                }
            }
            lm += (l + 1);
        }
    }
    // Compute derivatives
    if (compute_derivatives)
    {
        // Allocate derivatives
        _dFlmp = new double[l_idx(l_max + 1)];
        // Define additional variables
        Eigen::VectorX<double> dTlm(N);
        std::vector<double> dtheta_dI(N), dlam_dI(N);
        double tan_u;
        for (int i = 0; i < N; i++)
        {
            tan_u = sin_u[i] / cos_u[i];
            dtheta_dI[i] = -sin_u[i] * cos_I / sqrt(1 - sin_I * sin_I * sin_u[i] * sin_u[i]);
            dlam_dI[i] = -sin_I * tan_u / (1 + cos_I * cos_I * tan_u * tan_u);
        }
        // Iterate over degree and order
        lm = 0;
        for (l = 0; l <= l_max; l++)
        {
            for (m = 0; m <= l; m++)
            {
                // Compute unit disturbing potential derivative along great circle
                for (i = 0; i < N; i++)
                {
                    dTlm[i] = plm[i].get_dPlm_bar(l, m) * dtheta_dI[i] *
                                  (cos(m * lam[i]) + sin(m * lam[i])) +
                              +plm[i].get_Plm_bar(l, m) *
                                  (-m * sin(m * lam[i]) + m * cos(m * lam[i])) * dlam_dI[i];
                }
                // Analyse perturbing potential with FFT
                y = rfft(dTlm);
                for (i = 0; i <= l; i++)
                {
                    C[i] = 2 * y[i].real() / N;
                    S[i] = -2 * y[i].imag() / N;
                }
                // Map coefficients Ci, Si to Flmp
                if (l % 2 == 0)
                {
                    _dFlmp[lm + l / 2] = m % 2 == 0 ? C[0] : -C[0];
                }
                if (l % 2 == m % 2)
                {
                    for (i = 0; i <= l; i++)
                    {
                        if (i % 2 == l % 2)
                        {
                            _dFlmp[lm + (l - i) / 2] = (C[i] + S[i]) / 2;
                            _dFlmp[lm + (l + i) / 2] = (C[i] - S[i]) / 2;
                        }
                    }
                }
                else
                {
                    for (i = 0; i <= l; i++)
                    {
                        if (i % 2 == l % 2)
                        {
                            _dFlmp[lm + (l + i) / 2] = -(C[i] + S[i]) / 2;
                            _dFlmp[lm + (l - i) / 2] = -(C[i] - S[i]) / 2;
                        }
                    }
                }
                lm += (l + 1);
            }
        }
    }
    else
    {
        _dFlmp = nullptr;
    }
}

Flmp &Flmp::operator=(const Flmp &other)
{
    if (this != &other)
    {
        // Copy basic attributes
        I = other.I;
        l_max = other.l_max;
        // Assign inclination functions
        if (other._Flmp != nullptr && other.l_max > 0)
        {
            int size = l_idx(other.l_max + 1);
            _Flmp = new double[size];

            for (int i = 0; i < size; ++i)
            {
                _Flmp[i] = other._Flmp[i];
            }
        }
        else
        {
            _Flmp == nullptr;
        }
        // Assign inclination functions derivatives
        if (other._dFlmp != nullptr && other.l_max > 0)
        {
            int size = l_idx(other.l_max + 1);
            _dFlmp = new double[size];

            for (int i = 0; i < size; ++i)
            {
                _dFlmp[i] = other._dFlmp[i];
            }
        }
        else
        {
            _dFlmp == nullptr;
        }
    }
    return *this;
}

Flmp::Flmp(const Flmp &other) : _Flmp(nullptr), _dFlmp(nullptr), I(other.I), l_max(other.l_max)
{
    // Assign inclination functions
    if (other._Flmp != nullptr && other.l_max > 0)
    {
        int size = l_idx(other.l_max + 1);
        _Flmp = new double[size];

        for (int i = 0; i < size; ++i)
        {
            _Flmp[i] = other._Flmp[i];
        }
    }
    // Assign inclination functions derivatives
    if (other._dFlmp != nullptr && other.l_max > 0)
    {
        int size = l_idx(other.l_max + 1);
        _dFlmp = new double[size];

        for (int i = 0; i < size; ++i)
        {
            _dFlmp[i] = other._dFlmp[i];
        }
    }
}

Flmp::~Flmp()
{
    delete[] _Flmp;
    if (_dFlmp)
        delete[] _dFlmp;
}

int Flmp::get_l_max() const { return l_max; }

double Flmp::get_Flmp(int l, int m, int p) const { return _Flmp[lmp_idx(l, m, p)]; }

double Flmp::get_Flmk(int l, int m, int k) const { return abs(k) > l ? 0 : _Flmp[lmk_idx(l, m, k)]; }

double Flmp::get_dFlmp(int l, int m, int p) const { return _dFlmp[lmp_idx(l, m, p)]; }

double Flmp::get_dFlmk(int l, int m, int k) const { return abs(k) > l ? 0 : _dFlmp[lmk_idx(l, m, k)]; }

double Flmp::get_Flmk_star(int l, int m, int k) const
{
    return 0.5 * (((k - 1) * cos(I) - m) / sin(I) * this->get_Flmk(l, m, k - 1) +
                  ((k + 1) * cos(I) - m) / sin(I) * this->get_Flmk(l, m, k + 1) +
                  -this->get_dFlmk(l, m, k - 1) + this->get_dFlmk(l, m, k + 1));
}
