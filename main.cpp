#include "kernel.cuh"
#include <cmath>
#include <iostream>
#include <numeric>
#include <iomanip>
#include <chrono>
#include <string.h>
int main()
{
    Config config;
    std::string theFunction;
    std::string theDim;
    std::string theStrat;
    config.func = function::Rastrigin;
    config.strat = improvment::Firstimprov;
    config.dimens = dimensions::d5;
    theStrat = "Firstimprov";
    switch (config.func)
    {
    case function::Rastrigin:
    {
        config.a = -5.12;
        config.b = 5.12;
        theFunction = "Rastrigin";
        break;
    }
    case function::Michalewicz:
    {
        config.a = 0;
        config.b = M_PI;
        theFunction = "Michalewicz";
        break;
    }
    case function::Dejong:
    {
        config.a = -5.12;
        config.b = 5.12;
        theFunction = "DeJong";
        break;
    }
    case function::Schwefel:
    {
        config.a = -500;
        config.b = 500;
        theFunction = "Schwefel";
        break;
    }

    default:
        break;
    }

    switch (config.dimens)
    {
    case dimensions::d5:
    {
        config.d = 5;
        theDim = "D5";
        break;
    }
    case dimensions::d10:
    {
        config.d = 10;
        theDim = "D10";
        break;
    }
    case dimensions::d30:
    {
        config.d = 30;
        theDim = "D30";
        break;
    }
    case dimensions::d100:
    {
        config.d = 100;
        theDim = "D100";
        break;
    }

    default:
        break;
    }

    config.it = 100;
    config.p = 5;
    int samples = 30;

    int segments = (config.b - config.a) * pow(10, config.p);
    config.bitsPerDim = static_cast<int>(std::ceil(log2(segments)));
    config.bits = config.bitsPerDim * config.d;
    config.threads = 64;
    config.blocks = (config.it + config.threads - 1) / config.threads;
    std::vector<double> result(samples);
    std::vector<double> runtime(samples);

    for (int i = 0; i < samples; i++) {
        //auto start = std::chrono::high_resolution_clock::now();
        std::vector<double> iniresult = launch(config);
       // auto end = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> duration = end - start;
        for (int i = 0; i < config.it; i++) {
            std::cout << iniresult[i] << " ";
        }
        std::cout << "\n";
        std::cout << "\n";
        std::cout << "\n";
        std::cout << "\n";
      
        //double finres = *std::min_element(iniresult.begin(), iniresult.end());
       // runtime[i] = duration.count();
        //std::cout << std::fixed << std::setprecision(config.p);
       // std::cout << i + 1 << "/" << samples << " " << "best: " << finres << "\n";
       // std::cout << "Runtime: " << runtime[i] << "\n";
        //result[i] = finres;
    }




    std::cout << "\n===================" << theFunction << " " << theDim << " " << theStrat << " Final Results : ================\n" << std::flush;
    std::cout << std::fixed << std::setprecision(config.p);
    std::cout << "Best result: " << *std::min_element(result.begin(), result.end()) << '\n' << std::flush;
    std::cout << "Best Runtime: " << *std::min_element(runtime.begin(), runtime.end()) << "\n\n" << std::flush;
    std::cout << "Average result: " << std::accumulate(result.begin(), result.end(), 0.0) / (double)result.size() << '\n' << std::flush;
    std::cout << "Average Runtime: " << std::accumulate(runtime.begin(), runtime.end(), 0.0) / (double)runtime.size() << "\n\n" << std::flush;
    std::cout << "Worst result: " << *std::max_element(result.begin(), result.end()) << '\n' << std::flush;
    std::cout << "Worst Runtime: " << *std::max_element(runtime.begin(), runtime.end()) << '\n' << std::flush;






}