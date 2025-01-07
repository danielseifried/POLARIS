#include "MathProbList.hpp"


uint prob_list::size()
{
    return N + 1;
}

double prob_list::X(uint i)
{
    return x[i];
}

void prob_list::resize(uint size)
{
    if(x != 0)
        delete[] x;

    N = size - 1;
    x = new double[size];

    for(uint i = 0; i < size; i++)
        x[i] = 0;
}

void prob_list::setValue(uint pos, double _x)
{
#ifdef DEBUG
    if(x == 0)
    {
        cout << "\nERROR: Spline was not initiated!" << endl;
        return;
    }
#endif
    x[pos] = _x;
}

uint prob_list::getIndex(double v)
{
    uint min = 0;

    if(v < x[0] || v > x[N] || N == 1)
        return 0;

    min = upper_bound(x, x+N+1, v) - x - 1;

    return min;
}

void prob_list::normalize(double integValue)
{
    if(integValue > 0)
        for(uint i = 0; i <= N; i++)
            x[i] /= integValue;
    else
        resize(N + 1);
}

prob_list operator-(prob_list & list1, prob_list & list2)
{
    if(list1.size() != list2.size())
        cout << "\nERROR: Probability lists have different lengths!";
    prob_list diff_list(list1.size());

    for(uint i = 0; i < list1.size(); i++)
        diff_list.setValue(i, list1.X(i) - list2.X(i));

    return diff_list;
}
