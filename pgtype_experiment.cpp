/*

        Author: Matheus S. S. Macedo

    This is the C++ implementation of the mathematical model of the relationship
    between phenotype and genotype of populations across their evolution over
    generations.

    The implemented model was presented by:

    Laarits, T., P. Bordalo, and B. Lemos. "Genes under weaker stabilizing selection 
    increase network evolvability and rapid regulatory adaptation to an environmental
    shift." Journal of Evolutionary Biology 29.8 (2016): 1602-1616.

    available at: https://doi.org/10.1111/jeb.12897

*/


//code
#include <iostream>
#include <time.h>
#include <string>
#include <chrono>
#include <stdlib.h>

//math
#include <vector>
#include <math.h>
#include <random>
#include <Eigen/Dense>

//saving
#include <fstream>

using namespace std;
using namespace Eigen;

class MPopulation
{
private:
    //random number generator
    //seed is restarted whenever a distribution or random number is cast
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937_64 generator;

    //parameters
    int m_size;
    int pop_size;
    int generations;
    double epsilon;
    double l;
    int tau;
    double alpha;
    double nonzero_prob;
    double mu;
    double mut_prob;
    bool stable;
    bool get_condition = false;

    //matrix elements distribution
    string dist_config = "normal";

    //normal distribution parameters
    //double mean = 0.0;
    //double std_dev = 1.0;
    double mean;
    double std_dev;

    Eigen::MatrixXd initial_matrix;
    vector<Eigen::MatrixXd> population;
    vector<vector<int>> parents;
    Eigen::MatrixXd phenotypes;
    vector<Eigen::MatrixXd> old_phen;
    vector<double> fitness;
    vector<double> selection_prob;
    vector<double> condition_n;

    vector<double> target;
    vector<double> gene_expression;

    vector<vector<double>> mean_phenotype;
    vector<double> mean_fitness;

    //functions
    double phi(Eigen::MatrixXd phen, vector<Eigen::MatrixXd> prev_phen, int time, int v_position);
    double rho(double sum);
    double fit(Eigen::MatrixXd phen);
    Eigen::MatrixXd calc_phenotype(Eigen::MatrixXd matrix, Eigen::MatrixXd prev_phen);

    //methods
    void draw_initial_population();
    void update_pop_fitness();
    void calc_reprod_prob(int time);
    void select_parents();
    void draw_new_population();
    void update_pop_phenotypes();

    //saving
    string file_suffix = "";
    int save_freq;
    bool save_bool;
    void save(int time);
    void store_mean_fitness();
    void store_mean_phenotype();
    void store_pop_condition();

public:
    //constructor
    MPopulation(int matrix_dimensions,int generations, vector<double> target, vector<double> gene_expression);

    //methods
    void reproduction();
    void save_means();

    //setters and getters
    void setFileSuffix(string suffix);
    void setDistribution(string dist_type);
    void setDistribution(string normal_dist, double mean, double std_dev);
    void setCalculateCondition(bool y_or_n);
    //void setSaveBool(bool y_or_n);
    void setSavingFrequency(bool y_or_n, int frequency);

    //functions
    double conditionNumber(Eigen::MatrixXd matrix);

    //destructor
    ~MPopulation();
};

#include "MPopulation_constructors.H"
#include "MPopulation_savers.H"
#include "MPopulation_functions.H"
#include "MPopulation_methods.H"

void help_output()
{
    cout << "\nAvailable options are: " << endl;
    cout << "-dim: " << "Matrix dimensions, standard = 19." << endl;
    cout << "-dist: " << "Distribution to draw nonzero matrix elements and mutations from" << endl;
    cout << "       normal (normal distribution)" << endl << "       enter distribution mean and std. dev. through options -mean and -stddev" << endl;
    cout << "       standard is normal distribution with mean 0.0 and std. dev 1.0." << endl;
    cout << "-nsim: " << "Number of simulations to be executed, standard = 1." << endl;
    cout << "-gens: " << "Number of generations per simulation, standard = 12000." << endl;
    cout << "-cond: " << "Wether or not to calculate population`s condition numbers. 0 (false) or 1 (true), standard is false!" << endl;
    cout << "       This option greatly increases runtime." << endl;
    cout << "-savepop: " << "Wether or not to save the genotypes and phenotypes during runtime (standard saving frequency is of 1000 generations)" << endl;
    cout << "       0 (false) or 1 (true), standard is TRUE. Mean phenotypes and fit values will ALWAYS be saved." << endl;
    cout << "-savefreq: " << "Frequency of generations at which populations will be saved. Must be an integer." << endl;
    cout << "       Standard value is 1000. Note that if -savepop is set to 0, this will be ignored." << endl;
    cout << endl;
}

#include "declare_standard_parameters.H"

int main(int argc, char const *argv[])
{
    clock_t start, end;
    start = clock();

    #include "execution_options.H"
    for (int r = 0; r < runs; r++)
    {
        if (runs > 1)
        {
            cout << "\nSIMULATION " << r+1 << " OF " << runs << endl;
        }
        
        //initializing experiment object
        MPopulation experiment(dim, gens, target, gene_expr);
        //setting parameters
        if (r == 0)
        {
            experiment.setDistribution(dist, mean, std_dev);
            experiment.setCalculateCondition(cond);
            //experiment.setSaveBool(save_bool);
            experiment.setSavingFrequency(save_bool, save_freq);
        }
        string suffix = to_string(r);
        experiment.setFileSuffix(suffix);    
        experiment.reproduction();
        experiment.save_means();
    }

    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    cout << "\nExecution time = " << time_taken; 
    cout << " s. " << endl; 
    //system("pause");   

    return 0;
}

