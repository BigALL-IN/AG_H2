#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <thread>
#include <future>
#include <chrono>
#include <fstream>
#include <numeric>

#include "settings.h"
#include "functions.h"


std::random_device rd;
std::mt19937 gen(rd());
std::bernoulli_distribution distribution(0.5);
std::uniform_real_distribution<> realDist(0.0, 1);

int bestVersion = 0;
Config settings;

int Nbit() {
    int N = (settings.b - settings.a) * pow(10, settings.p);
    int bits = static_cast<int>(std::ceil(log2(N)));

    return bits;
}

int Cbit() {
    int bit = Nbit();
    return bit * settings.d;
}

bitstring Gen_num() {
    int bitcount = Cbit();
    bitstring vc(bitcount);
    for (int i = 0; i < vc.size(); i++) {
        vc[i] = distribution(gen);
    }
    return vc;
}
double Convert(bitstring vc)
{
    long long decval = 0;
    for (bool bit : vc) {
        decval = (decval << 1) | bit;
    }
    double finval = settings.a + ((static_cast<double>(decval) * (settings.b - settings.a)) / (pow(2, Nbit()) - 1));
    return finval;
}

double Eval(bitstring vc) {

    int vs = vc.size(); // vector size
    int cs = vs / settings.d; //chunk size
    std::vector<double> results(settings.d);
    bitstring aux;
    for (int i = 0; i < settings.d; ++i) {
        aux.clear();
        for (int j = i * cs; j < cs + cs * i; ++j) {
            aux.push_back(vc[j]);
        }
        results[i] = Convert(aux);
    }
    return Michalewicz(results);
}

//to be revised
double EvalFitness(double eval) {
    double results;

    switch (settings.func)
    {
    case function::Rastrigin:

        return (eval != 0) ? 1.0 / eval : 1e6;
        break;
    case function::Michalewicz:
        return exp(-0.025*eval);
        break;

    case function::Schwefel:
        return eval + 501;
        break;

    case function::Dejong:
        return (eval != 0) ? 1.0 / eval : 1e6;
        break;

    default:
        return -69;
        break;

    }

}
 


std::vector<double> CummulativeFitness(std::vector<Result> eval) {
    std::vector<double> helper, result;
    double sum = 0;
    for (int i = 0; i < eval.size(); i++) {
        sum += eval[i].fitness;
    }
    for (int i = 0; i < eval.size(); i++) {
        result.push_back(eval[i].fitness / sum);
    }
    for (int i = 1; i < eval.size(); i++) {
        result[i] = result[i-1] + result[i];

    }

    return result;
}

Result InitSelection(std::vector<double> eval, std::vector<Result> population) { 
    Result result;
    for (int i = 1; i < eval.size(); i++) {
        double randomDouble = realDist(gen);
        if (eval[0] >= randomDouble > 0) {
            return population[0];
        }
        if ((eval[eval.size()-1] < randomDouble)) {
            return population[population.size()-1];

        }
        if (eval[i - 1] < randomDouble <= eval[i]) {
            return population[i];
           
        }

    }

}
bitstring Swap_Ranges(bitstring C1, bitstring C2, int j) {
    for (int i = j; i < C1.size(); i++) {
        C1[i] = C2[i];
 }
    return C1;
}
std::vector<bitstring> crossOver(std::vector<Result>population, std::vector<double> q) {

    double p = realDist(gen);
    double prob1;
    int size = q.size();
    bitstring C1, C2, C1c, C2c;
    std::vector<bitstring> results;
    int j1, j2;
    
        double prob = realDist(gen); 
       
            bool swapped = false, diffrent = false;
            prob1 = realDist(gen);
            while (diffrent == false) {

                double prob2 = realDist(gen);

                for (int j = 0; j < size-1; ++j) {
                    if (q[0] >= prob1) {
                        j1 = 0;
                    }
                    if (q[0] >= prob2) {
                        j2 = 0;
                    }

                    if (q[q.size()-1] < prob1) {
                        j1 = q.size() - 1;
                    }

                    if (q[q.size() - 1] < prob1) {
                        j2 = q.size() - 1;
                    }

                    if (q[j] < prob1 && prob1 <= q[j + 1])
                        j1 = j;

                    if (q[j] < prob2 && prob2 <= q[j + 1])
                        j2 = j;

                }



                if (j1 != j2)diffrent = true;

            }

       
            C1 = population[j1].individual;
            C2 = population[j2].individual;
        


            if (prob < 0.8) {
                while (!swapped) {

                    for (int j = 0; j < C1.size(); ++j) {
                        double crossProb = realDist(gen);
                        double crossProbThreshold = 5.0 / C1.size();
                        if (crossProb <crossProbThreshold) {
                            C1c = C1;
                            C1 = Swap_Ranges(C1, C2, j);
                            C2 = Swap_Ranges(C2, C1c, j);
                            swapped = true;
                        }


                    }

                }

            }


    
            results.push_back(C1);
            results.push_back(C2);
            C1.clear();
            C2.clear();
            return results;


}

bitstring mutateInstance(bitstring candidate) {
    bool mutated = false;

    while (!mutated) {
        for (int i = 0; i < candidate.size(); ++i) {
            double p = realDist(gen);
            double mutateProb = 1.0 / candidate.size();
            if (p < mutateProb) {
                candidate[i] = !candidate[i];
                mutated = true;

            }
        }
    }

   return candidate;

}
Result initpop() {
    Result result;
    result.individual = Gen_num();
    result.eval=(Eval(result.individual));
    result.fitness = EvalFitness(result.eval);
  
    return result;
}

std::vector<Result> individual(std::vector<Result> population, std::vector<double>cummulativeFitness) {
    std::vector<Result> result(2);
    std::vector<bitstring> instances(2);
    std::vector<double> candidates(2), fitnesses(2);
    instances = crossOver(population, cummulativeFitness);
    result[0].individual = mutateInstance(instances[0]);
    result[1].individual = mutateInstance(instances[1]);
    result[0].eval = Eval(instances[0]);
    result[1].eval = Eval(instances[0]);
    result[0].fitness = EvalFitness(result[0].eval);
    result[1].fitness = EvalFitness(result[1].eval);
    return result;
}

double FindMin(std::vector<Result> population) {
    double min = population[0].eval;
    for (int i = 1; i < population.size(); i++) {
        if (min > population[i].eval) {
            min = population[i].eval;
        }
    }
    return min;
}
int main()
{

    settings.func = function::Michalewicz;
    settings.d = 5;
   
    switch (settings.func)
    {
    case function::Rastrigin:
    {
        settings.a = -5.12;
        settings.b = 5.12;
        break;
    }
    case function::Michalewicz:
    {
        settings.a = 0;
        settings.b = M_PI;
        break;
    }
    case function::Dejong:
    {
        settings.a = -5.12;
        settings.b = 5.12;  
        break;
    }
    case function::Schwefel:
    {
        settings.a = -500;
        settings.b = 500;
        break;
    }

    default:
        break;
    }
   
    settings.it = 100;
    settings.p = 5;
    settings.pop = 200;
    int samples = 1;
    int counter = 0;
     int nthreads = std::thread::hardware_concurrency();
    constexpr double inf{ std::numeric_limits<double>::infinity() };
    std::vector<double> sampleCandidates(samples);
    std::vector<double> sampleRuntimeDurations(samples);

  

    std::fill(sampleCandidates.begin(), sampleCandidates.end(), inf);

    while (counter < samples) {
        double minCandidate;
        std::vector<std::future<Result>> futures1(nthreads);
        std::vector<Result> gen0p; 
        std::vector<Result> gen0;
        std::vector<double> cummulativeFitnesses;
        int iterationCount = 0;
        int popcount = 0;

        //populate first generation with twice as many individuals in order to select from them
        while (popcount < settings.pop * 2) { 
            for (int i = 0; i < nthreads; i++) { 
                futures1[i] = std::async(initpop);
            }

            for (int i = 0; i < nthreads; i++) {
                gen0p.push_back(futures1[i].get());
            }
            popcount += nthreads; 
        }

        cummulativeFitnesses = CummulativeFitness(gen0p);
        popcount = 0; 
       
        //perform select on initial population
        while (popcount < settings.pop) {  
            for (int i = 0; i < nthreads; i++)  {
                futures1[i] = std::async(InitSelection, cummulativeFitnesses, gen0p);
            }

            for (int i = 0; i < nthreads; i++) {
                gen0.push_back(futures1[i].get());
            }
            popcount += nthreads;
        }
        double bestcandidate = FindMin(gen0);
        std::cout << "BEST CANDIDATE OF GENERATION 0: " << bestcandidate << "\n";

        cummulativeFitnesses.clear();
        cummulativeFitnesses = CummulativeFitness(gen0);
        auto start = std::chrono::high_resolution_clock::now();


        //start the loop
        while (iterationCount < settings.it) {
            popcount = 0;
            cummulativeFitnesses.clear();
            cummulativeFitnesses = CummulativeFitness(gen0);
            std::vector<std::future<std::vector<Result>>> futures(nthreads);
            std::vector<std::vector<Result>> candidates(nthreads);
           
            // generate new population
            while (popcount < settings.pop) {
                int i = 0;
                for (i = 0; i < nthreads; i++) {
                    futures[i] = std::async(individual, gen0, cummulativeFitnesses);
                }


     
                for (i = 0; i < nthreads; i++) {
                    candidates[i] = futures[i].get();
                    //place new population together with old population
                    for (int j = 0; j < 2; j++) {
                        if (settings.pop + popcount + i < settings.pop * 2 - 1) {
                            gen0.push_back(candidates[i][j]);
                        }
                    }
                    //gen0.insert(gen0.begin()+settings.pop, candidates[i].begin(), candidates[i].end());
                }
                popcount += nthreads;
            }

            //free up old initial population, move the total population in it, then free up current population to use for next run
            gen0p.clear();
            gen0p = gen0; 
            gen0.clear();
            popcount = 0;
            cummulativeFitnesses.clear();
            cummulativeFitnesses = CummulativeFitness(gen0p);
            //perform selection
            while (popcount < settings.pop) {
                for (int i = 0; i < nthreads; i++) {
                    futures1[i] = std::async(InitSelection, cummulativeFitnesses, gen0p);
                }

                for (int i = 0; i < nthreads; i++) {
                    gen0.push_back(futures1[i].get());
                }
                popcount += nthreads;
            }

minCandidate = FindMin(gen0);

            if (minCandidate < bestcandidate) {
                bestcandidate = minCandidate;
            }   

            iterationCount += 1;
            std::cout << "BEST CANDIDATE OF GENERATION " << iterationCount << ": " << minCandidate << "/BEST OVERALL CANDIDATE:" << bestcandidate << "\n";
           ;
            double lowest = *std::min_element(sampleCandidates.begin(), sampleCandidates.end());

        }



        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;


        sampleCandidates[counter] = bestcandidate;
        sampleRuntimeDurations[counter] = duration.count();

        counter++;


    }

    std::cout << "\n\n===================Final Results================\n\n\n" << std::flush;
    std::cout << "\n\n" << std::flush;;

    std::cout << "Best result: " << *std::min_element(sampleCandidates.begin(), sampleCandidates.end()) << '\n' << std::flush;
    std::cout << "Best Runtime: " << *std::min_element(sampleRuntimeDurations.begin(), sampleRuntimeDurations.end()) << "\n\n" << std::flush;
    std::cout << "Average result: " << std::accumulate(sampleCandidates.begin(), sampleCandidates.end(), 0.0) / (double)sampleCandidates.size() << '\n' << std::flush;
    std::cout << "Average Runtime: " << std::accumulate(sampleRuntimeDurations.begin(), sampleRuntimeDurations.end(), 0.0) / (double)sampleRuntimeDurations.size() << "\n\n" << std::flush;
    std::cout << "Worst result: " << *std::max_element(sampleCandidates.begin(), sampleCandidates.end()) << '\n' << std::flush;
    std::cout << "Worst Runtime: " << *std::max_element(sampleRuntimeDurations.begin(), sampleRuntimeDurations.end()) << '\n' << std::flush;

    return 0;
}
