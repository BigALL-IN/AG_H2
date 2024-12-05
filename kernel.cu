#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <curand_kernel.h>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <numeric>
#include <iomanip>

#include "kernel.cuh"


__constant__ Config globalConfig;


__global__ void Init_states(curandState* states, long long seed) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= globalConfig.it) return;
    curand_init(seed, idx, 1000, &states[idx]);
}

__global__ void InitBitstring(curandState* states, bool* b) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= globalConfig.it) return;
    int startBit = idx * globalConfig.bits;
    for (int i = startBit; i < startBit + globalConfig.bits; i++)
    {
        b[i] = curand_uniform(&states[idx]) > 0.5f;
    }
}

__device__ void Convert(bool* bits, double* values)
{
    for (int j = 0; j < globalConfig.d; j++) {
        unsigned long long dec = 0;
        for (int i = 0; i < globalConfig.bitsPerDim; i++)
        {
            dec = (dec << 1) | bits[j * globalConfig.bitsPerDim + i];

        }
        values[j] = globalConfig.a + dec * (globalConfig.b - globalConfig.a) / ((1ull << globalConfig.bitsPerDim) - 1);
    }
}
__global__ void GenRealValues(bool* bits, double* values) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= globalConfig.it) return;
    Convert(bits + idx * globalConfig.bits, values + idx * globalConfig.d);
}

__device__ double Rastrigin(double* v, int dimensions) {

    double res = 10 * dimensions;
    for (int i = 0; i < dimensions; i++) {
        res += v[i] * v[i] - 10 * cos(2 * M_PI * v[i]);
    }
    return res;
}

__device__ double Michalewicz(double* v, int dimensions) {
    double res = 0;
    for (int i = 0; i < dimensions; i++) {
        res += sin(v[i]) * pow(sin(((i + 1) * v[i] * v[i]) / M_PI), 20);
    }
    return -res;
}

//reminder to check if this si actually dejong
__device__ double Dejong(double* v, int dimensions) {
    double res = 0;
    for (int i = 0; i < dimensions; i++) {
        res += v[i] * v[i];
    }
    return res;
}

__device__ double Schwefel(double* v, int dimensions) {
    double res = 0;
    for (int i = 0; i < dimensions; i++) {
        res += -v[i] * sin(sqrt(abs(v[i])));
    }
    return res;
}

__device__ double Eval(double* values)
{
    switch (globalConfig.func)
    {
    case function::Rastrigin:
        return Rastrigin(values, globalConfig.d);
        break;
    case function::Michalewicz:
        return Michalewicz(values, globalConfig.d);
        break;

    case function::Schwefel:
        return Schwefel(values, globalConfig.d);
        break;

    case function::Dejong:
        return Dejong(values, globalConfig.d);
        break;
    }
}

__global__ void EvalValue(double* values, double* candidates)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= globalConfig.it) return;
    candidates[idx] = Eval(idx * globalConfig.d + values);
}



__global__ void CummulativeFitness(double* fitnessscores) {
    double result;
    for (int i = 0; i < globalConfig.it; i++) {
        result += fitnessscores[i];
    }
    for (int i = 0; i < globalConfig.it; i++) {
        fitnessscores[i] = fitnessscores[i] / result;
    }
    for (int i = 1; i < globalConfig.it; i++) {
        fitnessscores[i] = fitnessscores[i] + fitnessscores[i-1];
      
    } 
  
   
}

__global__ void mutate(curandState* states, bool* binPopulation, double* eval) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= globalConfig.it) return;
    for (int i = 0; i < globalConfig.it; ++i) {
        double probM = curand_uniform_double(&states[idx]);
        if (probM <= globalConfig.mutationRate) {
           mutateInstance(states, binPopulation, idx);
         
        }
    }
}


__device__ void mutateInstance(curandState* states, bool* candidate, int idx) {
    bool mutated = false;
    int startBit = idx * globalConfig.bits;
    while (!mutated) {
        for (int i = startBit; i < startBit + globalConfig.bits; ++i) {
            double p = curand_uniform_double(&states[idx]);
            if (p < 1 / globalConfig.bits) {
                candidate[i] = !candidate[i];
                mutated = true;
                break;

            }
        }
    }
}
__device__ void EvalFitnessdevice(double* candidates, double* fitnessscores, int idx) {
    //to be changed to e^ or even have switch case
    fitnessscores[idx] = -candidates[idx];
}
__global__ void EvalFitness(double* candidates, double* fitnessscores)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= globalConfig.it) return;
    EvalFitnessdevice(candidates, fitnessscores, idx);

}



__global__  void Algorithm(bool* bitstr, double* values, double* candidates, curandState* states) {

}

std::vector<double> launch(const Config& config) {

    bool* bitstr;
    double* candidates;
    double* realValues;
    double* fitnessScores;
    curandState* states;
    std::vector<double> result(config.it);

    // Allocate device memory
    cudaMalloc(&bitstr, sizeof(bool) * config.bits * config.it);
    cudaMalloc(&candidates, sizeof(double) * config.it);
    cudaMalloc(&fitnessScores, sizeof(double) * config.it);
    cudaMalloc(&states, sizeof(curandState) * config.it);
    cudaMalloc(&realValues, sizeof(double) * config.it * config.d);
    cudaMemcpyToSymbol(globalConfig, &config, sizeof(Config));


    // Launch kernel
    Init_states << < config.blocks, config.threads >> > (states, std::random_device{}());
    InitBitstring << < config.blocks, config.threads >> > (states, bitstr);
    GenRealValues << < config.blocks, config.threads >> > (bitstr, realValues);
    EvalValue << < config.blocks, config.threads >> > (realValues, candidates);
    EvalFitness << < config.blocks, config.threads >> > (candidates, fitnessScores);
    CummulativeFitness << < 1, 1 >> > (fitnessScores);

    //Algorithm << < config.blocks, config.threads >> > (bitstr, realValues, candidates, states);


    // Copy result back to host
    cudaMemcpy(result.data(), fitnessScores, sizeof(double) * config.it, cudaMemcpyDeviceToHost);


    // Clean up device memory
    cudaFree(bitstr);
    cudaFree(candidates);
    cudaFree(states);
    cudaFree(realValues);

    return result;
}