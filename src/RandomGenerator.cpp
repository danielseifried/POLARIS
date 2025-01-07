#include "RandomGenerator.hpp"


void CRandomGenerator::init(ullong seed)
{
    // KISS will work just fine even w/o the init process
    KISS_state[0] = seed;
    // The seed for CONG is given by the OMP thread number
    // call CONG several times to get a nice random seed for KISS
    for(int i = 0; i < 2000; i++)
        KISS_state[0] = CONG(KISS_state[0]);

    // fill the state array for KISS
    KISS_state[1] = CONG(KISS_state[0]);
    KISS_state[2] = CONG(KISS_state[1]);
    KISS_state[3] = CONG(KISS_state[2]);
}

ullong CRandomGenerator::CONG(ullong current_state)
{
    // CONG - very simple linear congruential generator
    return 6906969069ULL * current_state + 1234567;
}

double CRandomGenerator::getRND()
{
    // KISS (Keep it Simple Stupid) is a family of pseudorandom number generators
    // introduced by George Marsaglia.
    // Source: https://de.wikipedia.org/wiki/KISS_(Zufallszahlengenerator)

    // linear congruential generator
    KISS_state[2] = CONG(KISS_state[2]);

    // Xorshift
    KISS_state[1] ^= KISS_state[1] << 13;
    KISS_state[1] ^= KISS_state[1] >> 17;
    KISS_state[1] ^= KISS_state[1] << 43;

    // Multiply-with-carry
    ullong tmp = (KISS_state[0] << 58) + KISS_state[3];
    KISS_state[3] = (KISS_state[0] >> 6);
    KISS_state[0] += tmp;
    KISS_state[3] += (KISS_state[0] < tmp);

    // Return double between 0 and 1
    return double(KISS_state[0] + KISS_state[1] + KISS_state[2]) / 18446744073709551615ULL;
}

double CRandomGenerator::getRNDnormal(double mu, double sigma)
{
    double U1, U2, W, mult;
    double X1, X2;

    do
    {
        U1 = -1 + getRND() * 2;
        U2 = -1 + getRND() * 2;
        W = pow(U1, 2) + pow(U2, 2);
    } while(W >= 1 || W == 0);

    mult = sqrt((-2 * log(W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;

    double res = mu + sigma * X1;

    if(res < 0)
        return getRNDnormal(mu, sigma);

    return res;
}
