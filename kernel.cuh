#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <map>
enum class function { Rastrigin, Michalewicz, Dejong, Schwefel };
enum class improvment { Firstimprov, Bestimprov, Worstimprov, Annealing };
enum class dimensions { d5, d10, d30, d100 };

struct Config {

    double a;
    double b;
    double mutationRate;
    int p;
    int d;
    int it;
    int bits;
    int bitsPerDim;
    int temp;
    int threads;
    int blocks;
    dimensions dimens;
    function func;
    improvment strat;
};

std::vector<double> launch(const Config& configs);
