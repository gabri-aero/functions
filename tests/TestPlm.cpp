#include <functions/Plm.hpp>

#include "gtest/gtest.h"

TEST(Plm, Idx)
{
    int l_max = 100;
    Plm<> plm(l_max, 1);
    int i = 0;
    for (int l = 0; l <= l_max; l++)
    {
        for (int m = 0; m <= l; m++)
        {
            ASSERT_EQ(plm.lm_idx(l, m), i++);
        }
    }
}

TEST(Plm, Value)
{
    Plm<> plm(100, 65 * M_PI / 180);
    // Values retrieved from Matlab
    ASSERT_NEAR(plm.get_Plm(14, 4), -9.251507461437021e+03, 1e-10);
    ASSERT_NEAR(plm.get_Plm(97, 26), 1.765752185461010e+49, 1e36);
}

TEST(Plm, Derivatives)
{
    auto theta = 65 * M_PI / 180;
    double dtheta = 5e-5 * M_PI / 180;
    Plm<> pa(100, theta + dtheta);
    Plm<> pb(100, theta - dtheta);
    Plm<> plm(100, theta, true);
    int l = 13;
    int m = 5;
    double dPlm_num = (pa.get_Plm_bar(l, m) - pb.get_Plm_bar(l, m)) / (2 * dtheta);
    ASSERT_NEAR((plm.get_dPlm_bar(l, m) - dPlm_num) / dPlm_num, 0, 1e-7);
}

TEST(Plm, SecondDerivatives)
{
    auto theta = 65 * M_PI / 180;
    double dtheta = 5e-5 * M_PI / 180;
    Plm<> pa(100, theta + dtheta, true);
    Plm<> pb(100, theta - dtheta, true);
    Plm<> plm(100, theta, true, true);
    int l = 13;
    int m = 5;
    double ddPlm_num = (pa.get_dPlm_bar(l, m) - pb.get_dPlm_bar(l, m)) / (2 * dtheta);
    ASSERT_NEAR((plm.get_ddPlm_bar(l, m) - ddPlm_num) / ddPlm_num, 0, 1e-7);
}