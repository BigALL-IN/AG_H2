#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <thread>
#include <future>
#include <chrono>
#include <fstream>
#include <numeric>

#include "config.h"
#include "functions.h"




std::random_device rd;
std::mt19937 gen(rd());
std::bernoulli_distribution distribution(0.5);
std::uniform_real_distribution<> realDist(0.0, 0.1);

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
        return -eval;
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
    std::vector<double> result;
    double sum;
    for (int i = 0; i < eval.size(); i++) {
        sum += eval[i].eval;
    }
    for (int i = 0; i < eval.size(); i++) {
        result[i] = eval[i].eval / sum;
    }
    for (int i = 1; i < eval.size(); i++) {
        result[i] = eval[i].eval + eval[i - 1].eval;

    }

    return result;
}

Result InitSelection(std::vector<double> eval, std::vector<Result> population) {
    double randomDouble;
    int j = 0;
    Result result;
    for (int i = 1; i < eval.size(); i++) {
        randomDouble = realDist(gen);
        if (eval[0] >= randomDouble) {
            result = population[0];
            j++;
            break;
        }
        if (eval[eval.size()-1] < randomDouble) {
            result = population[population.size()];
            j++;
            break;
        }
        if (eval[i - 1] < randomDouble <= eval[i]) {
            result = population[i];
            break;
        }
    }

}


std::vector<bitstring> crossOver(std::vector<Result>population, std::vector<double> q) {

    double p = realDist(gen);
    double prob1, prob2;
    int size = population.size();
    bitstring C1, C2;
    std::vector<bitstring> results;
    int j1, j2;
    
        double prob = realDist(gen); 
       
            bool swapped = false, diffrent = false;
            prob1 = realDist(gen);
            while (diffrent == false) {

                prob2 = realDist(gen);

                for (int j = 0; j < size; ++j) {
                    if (q[0] >= prob1) {
                        j1 = 0;
                    }
                    if (q[0] >= prob2) {
                        j1 = 0;
                    }

                    if (q[q.size()-1] < prob1) {
                        j1 = q.size() - 1;
                    }

                    if (q[q.size() - 1] < prob1) {
                        j1 = q.size() - 1;
                    }

                    if (q[j] < prob1 && prob1 <= q[j + 1])
                        j1 = j;

                    if (q[j] < prob2 && prob2 <= q[j + 1])
                        j2 = j;

                }

                if (j1 != j2)diffrent = true;

            }

       

            std::copy(population[j1].individual[0], population[j1].individual[population[j1].individual.size() - 1], C1);
            std::copy(population[j2].individual[0], population[j2].individual[population[j2].individual.size() - 1], C1);

            if (prob < 0.8) {
                while (!swapped) {

                    for (int j = 0; j < C1.size(); ++j) {
                        double crossProb = realDist(gen);
                        if (crossProb < 1 / C1.size()) {

                            std::swap_ranges(C1[j], C1[C1.size() - 1], C2[j]);
                            swapped = true;
                        }


                    }

                }

                results.push_back(C1);
                results.push_back(C2);
                C1.clear();
                C2.clear();
                return results;
            

            }


    



}

bitstring mutateInstance(bitstring candidate) {
    bool mutated = false;

    while (!mutated) {
        for (int i = 0; i < candidate.size(); ++i) {
            double p = realDist(gen);
            if (p < (1 / candidate.size())) {
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
    std::vector<Result> result;
    std::vector<bitstring> instances;
    std::vector<double> candidates, fitnesses;
    double candidate0, candidate1;
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
}
int main()
{

    settings.func = function::Rastrigin;
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
        std::vector<Result> gen0p(settings.pop * 2); 
        std::vector<Result> gen0(settings.pop);
        std::vector<double> cummulativeFitnesses;
        int iterationCount = 0;
        int popcount = 0;

        //populate first generation with twice as many individuals in order to select from them
        while (popcount < settings.pop * 2) { 
            for (int i = 0; i < nthreads; i++) { 
                futures1[i] = std::async(initpop);
            }

            for (int i = 0; i <i; i++) { 
                gen0p[popcount+i] = futures1[i].get();
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
                gen0[i] = futures1[i].get();
            }
            popcount += nthreads;
        }
        double bestcandidate = FindMin(gen0);
        cummulativeFitnesses.clear();
        cummulativeFitnesses = CummulativeFitness(gen0);
        auto start = std::chrono::high_resolution_clock::now();


        //start the loop
        while (iterationCount < settings.it) {
            popcount = 0;
            std::vector<std::future<std::vector<Result>>> futures(nthreads);
            std::vector<std::vector<Result>> candidates(nthreads);
           
            // generate new population
            while (popcount < settings.pop) {
                for (int i = 0; i < nthreads; i++) {
                    futures[i] = std::async(individual, gen0, cummulativeFitnesses);
                }


                
                for (int i = 0; i < nthreads; i++) {
                    candidates[i] = futures[i].get();
                    //place new population together with old population
                    gen0.insert(gen0.end(), candidates[i].begin(), candidates[i].end());
                }
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
                    gen0[popcount+i] = futures1[i].get();
                }
                popcount += nthreads;
            }
            minCandidate = FindMin(gen0);


           


           // double bestcandidate = *std::min_element(evals.begin(), evals.end());

            if (minCandidate < bestcandidate) {
                bestcandidate = minCandidate;
            }

            iterationCount += 1;
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
