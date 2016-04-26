#include <cmath>
#include <iostream>
#include <polymul.h>

#ifdef _WIN32
#include <windows.h>
#define FREQ LARGE_INTEGER
#define TICK LARGE_INTEGER
#else
#include "sys/time.h"
#define TICK timeval
#endif

using namespace std;
using namespace polymul;

#ifdef _WIN32
static FREQ frequency;
#endif
static TICK starttick;

void start()
{
#ifdef _WIN32
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&starttick);
#else
    gettimeofday(&starttick, NULL);
#endif
}

double stop()
{
    TICK stoptick;
#ifdef _WIN32
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&stoptick);
    return (stoptick.QuadPart - starttick.QuadPart) / (double) frequency.QuadPart;
#else
    gettimeofday(&stoptick, NULL);
    return stoptick.tv_sec - starttick.tv_sec + (stoptick.tv_usec - starttick.tv_usec) * 1e-6;
#endif
}

template <class numtype, int Nvar, int Ndeg>
void benchmark(void)
{
    double dt;

    polynomial<numtype, Nvar, Ndeg> p1, p2, pt;
    polynomial<numtype, Nvar, 2 * Ndeg> p3, sum;
    for (int i = 0; i < p1.size; i++)
    {
        p1[i] = cos(i);
        p2[i] = cos(3 * i);
    }
    sum = 0;

    cout << "Benchmarking " << Nvar << ", " << Ndeg << ": " << p1.size << " variables." << endl;

    // Multiplication
    int nrep = 2e7 / p1.size;
    start();
    for (int i = 0; i < nrep; i++)
    {
        p1[0] = cos(i);
        mul(p3, p1, p2);
        //for (int j=0;j<sum.size;j++)
        //    sum[j] += p3[j];
    }
    dt = stop();
    cout << "polymul " << nrep / dt << " muls/s " << sum[0] << " " << dt << endl;

    // Taylor Multiplication
    nrep = 20e7 / p1.size;

    start();
    for (int i = 0; i < nrep; i++)
    {
        p1[0] = cos(i);
        taylormul(pt, p1, p2);
    }
    dt = stop();
    cout << "taylormul " << nrep / dt << " muls/s " << p1[0] << " " << dt << endl;

    // Taylor in place Multiplication
    nrep = 20e6 / p1.size;

    start();
    for (int i = 0; i < nrep; i++)
    {
        p1[0] = cos(i);
        taylormul(p1, p2);
    }
    dt = stop();
    cout << "taylormul inplace " << nrep / dt << " muls/s " << p1[0] << " " << dt << endl;

    // Linear transformation
    nrep = 1e7 / p1.size;
    numtype T[Nvar * Nvar];
    for (int i = 0; i < Nvar * Nvar; i++)
        T[i] = sin(i);

    start();
    for (int i = 0; i < nrep; i++)
    {
        p1[0] = cos(i);
        T[0] = sin(i);
        trans(p1, p2, T);
    }
    dt = stop();
    cout << "polytrans " << nrep / dt << " trans/s " << p1[0] << " " << dt << endl;
}

int main(void)
{
    benchmark<double, 4, 4>();
}

