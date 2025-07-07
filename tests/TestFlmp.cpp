#include <chrono>
#include <complex>

#include <functions>
#include <gtest/gtest.h>

double get_complex_Flmp(const Flmp &flmp, int l, int m, int p)
{
    auto j = std::complex<double>(0, 1.0);
    return (pow(j, l - m) * flmp.get_Flmp(l, m, p)).real();
};

TEST(Flmp, Value100)
{
    double I = 109.9 * M_PI / 180;
    Flmp flmp(100, I);

    ASSERT_NEAR(get_complex_Flmp(flmp, 15, 15, 7), 0.163727788669698, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 17, 15, 8), 0.487417791777481, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 19, 15, 9), 0.039444885080361, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 21, 15, 10), -0.334234993689438, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 23, 15, 11), 0.238101170358486, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 25, 15, 12), 0.035197122324998, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 27, 15, 13), -0.238961053270882, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 29, 15, 14), 0.250820102027528, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 31, 15, 15), -0.098284229213865, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 33, 15, 16), -0.099812590952652, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 35, 15, 17), 0.220401483107786, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 37, 15, 18), -0.203459255803049, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 39, 15, 19), 0.072853902584608, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 41, 15, 20), 0.089117362850045, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 43, 15, 21), -0.192487848426302, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 45, 15, 22), 0.186917527873700, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 47, 15, 23), -0.083106948025162, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 49, 15, 24), -0.058636531371390, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 51, 15, 25), 0.163214940273027, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 53, 15, 26), -0.179533365185972, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 55, 15, 27), 0.104101730627469, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 57, 15, 28), 0.020582796611666, 1e-10);
    ASSERT_NEAR(get_complex_Flmp(flmp, 59, 15, 29), -0.129982540091162, 1e-10);
}

TEST(Flmp, Derivative100)
{
    const int l_max = 100;
    double I = 25 * M_PI / 180;
    Flmp flmp(l_max, I, true);

    ASSERT_NEAR(std::abs(flmp.get_dFlmp(15, 15, 7)), 0.000193588834461, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(17, 15, 8)), 0.002962643282053, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(19, 15, 9)), 0.019210738800719, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(21, 15, 10)), 0.080996204022307, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(23, 15, 11)), 0.254529309868877, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(25, 15, 12)), 0.635791817300206, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(27, 15, 13)), 1.304718954007593, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(29, 15, 14)), 2.229338572015512, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(31, 15, 15)), 3.154511340102659, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(33, 15, 16)), 3.561310705132132, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(35, 15, 17)), 2.797301141098675, 1e-10);

    ASSERT_NEAR(std::abs(flmp.get_dFlmp(59, 15, 29)), 7.135563481217891, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(61, 15, 30)), 13.533758345144610, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(63, 15, 31)), 12.842455780020720, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(65, 15, 32)), 4.896633828451622, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(67, 15, 33)), 6.247154772263426, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(69, 15, 34)), 14.285109814165770, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(71, 15, 35)), 14.262965486747120, 1e-10);
    ASSERT_NEAR(std::abs(flmp.get_dFlmp(73, 15, 36)), 5.729761501008049, 1e-10);
}