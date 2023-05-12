//uniform integer distribution to get zeroed entries
uniform_int_distribution<int> uniform_dist(0,9);
int zero_thre = int (nonzero_prob*10);

//uniform real distribution for initial phenotypes
uniform_real_distribution<double> uniform_real_dist(0.0, 1.0);

normal_distribution<double> n_dist(this->mean,this->std_dev);

for (int i = 0; i < m_size; i++)
{
    for (int j = 0; j < m_size; j++)
    {
        int zero_prob = uniform_dist(generator);
        if (zero_prob <= zero_thre)
        {
            this->initial_matrix(i,j) =  n_dist(generator);
        }
        else
        {
            this->initial_matrix(i,j) = 0.0;
        } 
                
    }
}

for (int i = 0; i < m_size; i++)
{
    for (int j = 0; j < pop_size; j++)
    {
        this->phenotypes(i,j) = uniform_real_dist(generator);
        this->population[j] = initial_matrix;
    }
}