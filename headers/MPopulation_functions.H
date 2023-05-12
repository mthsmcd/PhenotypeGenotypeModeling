double MPopulation::phi(Eigen::MatrixXd phen, vector<Eigen::MatrixXd> prev_phen, int time, int v_position)
{
    double value = 0.0;   
    for (int t = (time - tau); t < time; t++)
    {
        Eigen::MatrixXd diff = prev_phen[t].col(v_position) - phen.col(v_position);
        diff = diff.cwiseAbs();
        value += diff.mean();
    }
    value /= tau;
    return value;
}

double MPopulation::rho(double x)
{
    double value = 1/(1 + exp(-alpha*x));
    return value;
}

double MPopulation::fit(Eigen::MatrixXd phen)
{
    double sum = 0.0;
    for (int r = 0; r < m_size; r++)
    {
        double diff = fabs(phen(r) - target[r]);
        if (diff <= l)
        {
            diff = 0.0;
        }
        sum += diff/gene_expression[r];
    }

    double value = exp(-sum);
    return value;
}

Eigen::MatrixXd MPopulation::calc_phenotype(Eigen::MatrixXd matrix, Eigen::MatrixXd prev_phen)
{
    Eigen::MatrixXd vetor(m_size,1);
    Eigen::MatrixXd aux = matrix * prev_phen;

    for (int j = 0; j < m_size; j++)
    {
        vetor(j) = rho(aux(j));
    }

    return vetor;
} 

void MPopulation::update_pop_fitness()
{
    for (int i = 0; i < pop_size; i++)
    {
        this->fitness[i] = fit(phenotypes.col(i));
    }
}

void MPopulation::update_pop_phenotypes()
{

    Eigen::MatrixXd aux_phenotypes = phenotypes;
 
    for (int i = 0; i < pop_size; i++)
    {
        phenotypes.col(i) = calc_phenotype(population[i], aux_phenotypes.col(i));
    }
}

void MPopulation::store_mean_fitness()
{
    double sum = 0.0;
    for (int i = 0; i < pop_size; i++)
    {
        sum += fitness[i];
    }

    sum /= pop_size;
    mean_fitness.push_back(sum);
}

void MPopulation::store_mean_phenotype()
{
    vector<double> mf(m_size);

    for (int i = 0; i < m_size; i++)
    {
        mf[i] = phenotypes.row(i).mean();
    }

    mean_phenotype.push_back(mf);
}

double MPopulation::conditionNumber(Eigen::MatrixXd matrix)
{
    Eigen::MatrixXd inv_matrix = matrix.inverse();
    double norm = matrix.lpNorm<Infinity>();
    double inv_norm = inv_matrix.lpNorm<Infinity>();
    double k = norm*inv_norm;
    return k;
}

void MPopulation::store_pop_condition()
{
    Eigen::VectorXd conditions(pop_size);
    for (int i = 0; i < pop_size; i++)
    {
        conditions(i) = conditionNumber(population[i]);
    }

    condition_n.push_back(conditions.mean());
}