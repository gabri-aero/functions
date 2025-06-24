/**
 * @file Nlm.cpp
 *
 * @brief Source file for computing normalization constants for spherical harmonics
 *
 * @author Gabriel Valles
 * @date 2025-02-17
 */

#include "Nlm.hpp"

int Nlm::lm_idx(int l, int m) { return (l * (l + 1)) / 2 + m; }

Nlm::Nlm() {};

Nlm::Nlm(int l_max) : l_max(l_max)
{
    this->_Nlm = new double[((l_max + 1) * (l_max + 2)) / 2];
    for (int l = 0; l <= l_max; l++)
    {
        // Compute for m = 0
        _Nlm[lm_idx(l, 0)] = sqrt(2 * l + 1);
    }
    for (int m = 1; m <= l_max; m++)
    {
        for (int l = m; l <= l_max; l++)
        {
            // Recursion for m > 0
            _Nlm[lm_idx(l, m)] = _Nlm[lm_idx(l, m - 1)] * sqrt(1.0 / ((l - m + 1) * (l + m)));
        }
    }
    // Adjust k for m>0
    for (int m = 1; m <= l_max; m++)
    {
        for (int l = m; l <= l_max; l++)
        {
            _Nlm[lm_idx(l, m)] *= sqrt(2);
        }
    }
}

Nlm::Nlm(const Nlm &other) : l_max(other.l_max)
{
    // Allocate and assign Nlm
    int Nlm_size = (l_max + 1) * (l_max + 2) / 2;
    _Nlm = new double[Nlm_size];
    std::copy(other._Nlm, other._Nlm + Nlm_size, _Nlm);
}

Nlm &Nlm::operator=(const Nlm &other)
{
    if (this != &other)
    {
        l_max = other.l_max;
        // Allocate and assign Nlm
        if (other._Nlm)
        {
            int Nlm_size = (l_max + 1) * (l_max + 2) / 2;
            _Nlm = new double[Nlm_size];
            std::copy(other._Nlm, other._Nlm + Nlm_size, _Nlm);
        }
    }
    return *this;
}

Nlm::~Nlm()
{
    if (_Nlm)
        delete[] _Nlm;
}

double Nlm::get_Nlm(int l, int m) { return _Nlm[lm_idx(l, m)]; }
