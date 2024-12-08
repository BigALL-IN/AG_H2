#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <map>
using bitstring = std::vector<bool>;
enum class function { Rastrigin, Michalewicz, Dejong, Schwefel };


struct Result {
    bitstring individual;
    double eval;
    double fitness;

};

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
    int pop;
    function func;
};




