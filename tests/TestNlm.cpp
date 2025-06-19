#include "gtest/gtest.h"
#include <functions/Plm.hpp>

double factorial(int n)
{
    return n == 1 ? 1 : n * factorial(n - 1);
};

TEST(Nlm, Value)
{
    Nlm<> nlm(10);
    double d0m;
    for (int l = 0; l < 5; l++)
    {
        for (int m = 0; m < l; m++)
        {
            d0m = m == 0 ? 1 : 0;
            EXPECT_NEAR(nlm.get_Nlm(l, m), sqrt((2 - d0m) * (2 * l + 1) * factorial(l - m) / factorial(l + m)), 1e-15);
        }
    }
}