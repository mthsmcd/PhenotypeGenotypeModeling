int multiplier = 1e9;
double corrector = 1e-10;
int lower_lim = int (multiplier*(mut_prob + corrector));
int upper_lim = int (multiplier*((1.0 - mut_prob) + corrector));
vector<int> limites = {lower_lim, upper_lim};
discrete_distribution<int> mut_dist(limites.begin(), limites.end());

normal_distribution<double> n_dist(this->mean,this->std_dev);

for (int p = 0; p < pop_size; p++)
{
    for (int i = 0; i < m_size; i++)
    {
        for (int j = 0; j < m_size; j++)
        {
            if (temp_population[p](i,j) != 0.0)
            {
                int draw = mut_dist(generator);
                if (!draw)
                {
                    temp_population[p](i,j) = n_dist(generator);
                }
            }
        }
    }
}