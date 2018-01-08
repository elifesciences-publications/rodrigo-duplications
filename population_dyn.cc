#ifndef _population_dyn_
#define _population_dyn_

//For compiling
//g++ -c -g -Wall -O4 population_dyn.cc
//g++ -o population_dyn population_dyn.o
//For executing
//./population_dyn evo_time ini_freq mean_signal

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#define RANDOM ( (double) (rand()/(RAND_MAX+1.0)) ) //uniform random number 0-1
#define NORMA_RANDOM ( (double) sqrt(-2*log(RANDOM))*cos(2*3.1416*RANDOM) ) //normal random number, mean 0, std 1 (Box–Muller transform)
#define DILUTION 100
#define MAX_N_POP 100000
#define IN_NOISE 0.5
#define EX_NOISE 0
#define VAR_SIGNAL 1

using namespace std;

/////////////////////////////////////

class cell {
public:
	cell();
	void set(bool copies, double env_molecule, double t_current);
    void reset(double env_molecule, double t_current);
    void reset(double t_advance);
    double compute_volume(double t_current);
	bool has_duplication() const;
    double get_diff_exp() const;
	double get_expression() const;
    double get_fitness() const;
	double get_t_born() const;
	double get_t_doubling() const;
private:
    void compute_fitness();
	bool duplication;
    double lactose;
    double diff_exp;
    double diff_exp_dupl;
    double expression;
	double fitness;
	double t_born;
	double t_doubling;
};

cell :: cell(){
}

void cell :: set(bool copies, double env_molecule, double t_current){
	duplication = copies;
    lactose = env_molecule;
    diff_exp = 0.0;
    diff_exp_dupl = 0.0;
    t_born = t_current;
    compute_fitness();
}

void cell :: reset(double env_molecule, double t_current){
    lactose = env_molecule;
    t_born = t_current;
    compute_fitness();
}

void cell :: reset(double t_advance){
    t_born += t_advance;
    t_doubling = 1.0/fitness + t_born;
}

double cell :: compute_volume(double t_current){
    double dt = t_current-t_born;
    double vol = pow(2,fitness*dt);
    return vol;
}

void cell :: compute_fitness(){
    double f0 = 1;
    double a = 0.17;
    double k = 0.40;
    double b = 0.036;
    double h = 1.8;
    double n = 4;
    double x0 = 0.13;
    double m = 1.0;
    
    double noiseX = exp(EX_NOISE*NORMA_RANDOM);
    double noise1 = exp(IN_NOISE*NORMA_RANDOM);
    double noise2 = exp(IN_NOISE*NORMA_RANDOM);
    if (!duplication)
        expression = (1.0+diff_exp)/(1 + pow(lactose*noise1*noiseX/x0, -n));
    else
        expression = (0.5+diff_exp)/(1 + pow(lactose*noise1*noiseX/x0, -n)) + (0.5+diff_exp_dupl)/(1 + pow(lactose*noise2*noiseX/x0, -n));
    
    fitness = f0*(1.0 + a*expression*lactose/(k + lactose) - b*pow(expression, m)/(h - expression));
    if (expression>=h)
        fitness = 0.000001; //it is indeed 0
    t_doubling = 1.0/fitness + t_born;
}

bool cell :: has_duplication() const{
	return duplication;
}

double cell :: get_diff_exp() const{
    return diff_exp;
}

double cell :: get_expression() const{
	return expression;
}

double cell :: get_fitness() const{
    return fitness;
}

double cell :: get_t_born() const{
	return t_born;
}

double cell :: get_t_doubling() const{
	return t_doubling;
}

/////////////////////////////////////

vector<cell> find_mutants (vector<cell> population, int option) {
    vector<cell> mutants;
    int n_pop = population.size();
    for (int i=0; i<n_pop; i++){
        if (option==2){ //find duplicates
            if (population[i].has_duplication())
                mutants.push_back(population[i]);
        }
        else{ //find singletons
            if (!population[i].has_duplication())
                mutants.push_back(population[i]);
        }
    }
    return mutants;
}

void reset_pop (vector<cell> &population, double t) {
    int n_pop = population.size();
    for (int i=0; i<n_pop; i++)
        population[i].reset(t);
}

double volume_pop (vector<cell> population, double t) {
    double vol = 0;
    int n_pop = population.size();
    for (int i=0; i<n_pop; i++)
        vol += population[i].compute_volume(t);
    return vol;
}

double frequency_dupl (vector<cell> population, double t){
    vector<cell> subpopulation;
    subpopulation = find_mutants(population, 2);
    double vol = volume_pop(population, t);
    double vol2 = volume_pop(subpopulation, t);
    double freq = vol2/vol;
    return freq;
}

vector<cell> bottleneck (vector<cell> population, double t) {
    vector<cell> singletons;
    vector<cell> duplicates;
    singletons = find_mutants(population, 1);
    duplicates = find_mutants(population, 2);
    double vol_sing = volume_pop(singletons, t);
    double vol_dupl = volume_pop(duplicates, t);
    double vol_sing_after = (vol_sing/(vol_sing+vol_dupl))*MAX_N_POP/DILUTION;
    double vol_dupl_after = (vol_dupl/(vol_sing+vol_dupl))*MAX_N_POP/DILUTION;
    
    vector<cell> subpopulation;
    int n_pop = singletons.size();
    double vol = 0;
    for (int i=0; i<n_pop; i++){
        int pos = (int)floor(RANDOM*n_pop);
        subpopulation.push_back(singletons[pos]);
        vol += singletons[pos].compute_volume(t);
        if(vol > vol_sing_after)
            break;
    }
    n_pop = duplicates.size();
    vol = 0;
    for (int i=0; i<n_pop; i++){
        int pos = (int)floor(RANDOM*n_pop);
        subpopulation.push_back(duplicates[pos]);
        vol += duplicates[pos].compute_volume(t);
        if(vol > vol_dupl_after)
            break;
    }
    return subpopulation;
}

/////////////////////////////////////

int main (int argc, char *argv[]) {

	double evo_time = atof( argv[1] );
    double ini_freq  = atof( argv[2] );
    double mean_in_signal  = atof( argv[3] );
	srand( time(NULL) );

    double env_molecule = mean_in_signal*exp(VAR_SIGNAL*NORMA_RANDOM);
    int dilution_period = 24;
    int t_bottleneck = dilution_period;
    double generations = log(DILUTION)/log(2);
    int t_passages = 0;
    int plot_period = 75; //in passages
    double t_plot = plot_period;

	//initialize
	vector <cell> population;
	int n_pop2 = (int)floor(ini_freq*MAX_N_POP/DILUTION);
    int n_pop1 = MAX_N_POP/DILUTION - n_pop2;
	for (int i=0; i<n_pop1; i++){
		cell c;
		c.set(false, env_molecule, 0);
		population.push_back(c);
    }
    for (int i=0; i<n_pop2; i++){
        cell c;
		c.set(true, env_molecule, 0);
		population.push_back(c);
	}
    double vol = MAX_N_POP/DILUTION;
    double t_generations = 0;
    double t_advance;
    double freq = frequency_dupl(population, 0);
    cout<< t_generations <<" "<< freq <<endl;
    
	//growing
	for (int t=1; t<=1000*evo_time; t++){ //t is 1000*time
        int n_pop = population.size();
		for (int i=0; i<n_pop; i++){
            if (vol < MAX_N_POP){
                vol += population[i].compute_volume((double)t/1000) - population[i].compute_volume((double)(t-1)/1000);
                if (t >= 1000*population[i].get_t_doubling()){
                    population[i].reset(env_molecule, (double)t/1000);
                    cell c = population[i];
                    population.push_back(c);
                }
			}
		}
        if (vol >= MAX_N_POP){
            t_advance = t_bottleneck - (double)t/1000;
            t = 1000*t_bottleneck;
            reset_pop(population, t_advance);
        }
        if (t == 1000*t_bottleneck){
            t_passages++;
            t_generations += generations;
            t_bottleneck += dilution_period;
            env_molecule = mean_in_signal*exp(VAR_SIGNAL*NORMA_RANDOM);
            freq = frequency_dupl(population, (double)t/1000);
            cout<< t_generations <<" "<< freq <<endl;
            population = bottleneck(population, (double)t/1000);
            vol = volume_pop(population, (double)t/1000);
        }
        
        if (t_passages == t_plot){
            freq = frequency_dupl(population, (double)t/1000);
            cout<< t_plot*generations <<" "<< freq <<endl;
            t_plot += plot_period;
        }
	}

	return 0;

}

#endif
