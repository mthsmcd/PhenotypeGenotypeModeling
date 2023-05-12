MPopulation::MPopulation(int matrix_dimensions, int generations, vector<double> target, vector<double> gene_expression)
{
    //configuring matrix dimensions
    this->m_size = matrix_dimensions;
    this->initial_matrix.resize(m_size, m_size);

    //configuring population dimensions
    this->pop_size = 100;
    this->population.resize(pop_size, Eigen::MatrixXd(m_size, m_size));

    //target phenotype
    this->target = target;
    //initial population's phenotype and fitness
    this->phenotypes.resize(m_size, pop_size);
    this->fitness.resize(pop_size);
    //number of generations
    this->generations = generations;
    this->gene_expression = gene_expression;

    //other parameters
    this->epsilon = 1e-3;
    this->tau = 10;
    this->l = 1e-2;
    this->alpha = 0.4;
    this->mu = 0.1;
    this->nonzero_prob = 0.6;
    //stability of population
    this->stable = false;

    cout<< "Initial population = " << population.size() << endl;
    cout << "Matrix dimensions = " << initial_matrix.rows() << "x" << initial_matrix.cols() << endl;

    //mu comes from the article
    this->mut_prob = mu/(nonzero_prob*m_size*m_size);
    cout << "The probability of mutations in non-zero elements of matrices is of " << mut_prob << endl;
}

void MPopulation::draw_initial_population()
{
    this->seed = chrono::system_clock::now().time_since_epoch().count();
    this->generator.seed(this->seed);

    //if (!dist_config.compare("normal"))
    //{
    #include "normal_population.H"
    //}
}

void MPopulation::setDistribution(string dist_type)
{
    this->dist_config = dist_type;
}

void MPopulation::setDistribution(string normal_dist, double m, double std)
{
    this->dist_config = normal_dist;
    this->mean = m;
    this->std_dev = std;
    cout << "Normal distribution with mean = " << mean << " and standard deviation = " << std_dev << endl;
    cout << "The probability of non-zero elements in the matrices is of " << nonzero_prob <<endl;
}

void MPopulation::setCalculateCondition(bool y_or_n)
{
    this->get_condition = y_or_n;
    if (get_condition)
    {
        cout << "\nCalculating and saving matrices condition numbers!" << endl;
    }     
}

/*void MPopulation::setSaveBool(bool y_or_n)
{
    this->save_bool = y_or_n;
    if (save_bool)
    {
        cout << "\nSaving genotypes and phenotypes every " << save_freq << " generations!" << endl;
    }
    else
    {
        cout << "\nOnly mean phenotypes and population fitness being saved at the end of each experiment!" << endl;
    }
}*/