#include<iostream>
#include"Movers/AutoRegressiveMover.h"
#include"EnsembleSampler.h"
using namespace MarkovChainMonteCarlo;

int main()
{
    float offsets[5] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    float phis[5] = {0.8f, 0.904761904762f, 0.9354838709677f, 0.9672131147541f, 0.990050200903734685f};
    float vars[5] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    Mover::AutoRegressiveMove<float> mover(5, offsets, phis, vars);
    EnsembleSampler<float, Mover::AutoRegressiveMove<float> > sampler(0, 100, 5, 2147483648ULL);
}
