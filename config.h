#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <map>
using bitstring = std::vector<bool>;
enum class function { Rastrigin, Michalewicz, Dejong, Schwefel };
enum class improvment { Firstimprov, Bestimprov, Worstimprov, Annealing };
enum class dimensions { d5, d10, d30, d100 };

struct Config{

    double a;
    double b;
    int p;
    int d;
    int it;
    int bits;
    int bitsPerDim;
    int temp;
    int threads;
    int blocks;
    int pop;
    dimensions dimens;
    function func;
    improvment strat;
};

struct Result {
    bitstring individual;
    double eval;
    double fitness;
};


std::vector<double> launch(const Config& configs);
