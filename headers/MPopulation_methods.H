void MPopulation::reproduction()
{

    cout << "\nInitializing matrices" << endl;
    draw_initial_population();
    if (save_bool)
    {
        save(0);
    }
    if (get_condition)
    {
        store_pop_condition();
    }
    
    cout << "\nBeginning iterative process" << endl;

    for (int k = 0; k < generations; k++)
    {

        if (k == 0)
        {
            update_pop_phenotypes();
            this->old_phen.clear();
            update_pop_fitness();
            store_mean_fitness();
            store_mean_phenotype();
        }

        calc_reprod_prob(k);

        if (k == 99 && !stable)
        {
            cout << "\nSimulation failed to meet the phenotype stability criteria in 100 steps." << endl;
            cout << "Restarting simulation..." << endl;

            draw_initial_population();
            if (save_bool)
            {
                save(0);
            }
            this->mean_fitness.clear();
            this->mean_phenotype.clear();
            if (get_condition)
            {
                this->condition_n.clear();
            }      
            k = -1;
        }
        else
        { 
            select_parents();
            draw_new_population();
            update_pop_phenotypes();
            update_pop_fitness();
            store_mean_phenotype();
            store_mean_fitness();
            if (get_condition)
            {
                store_pop_condition();
            }
            
            if ((k+1) % save_freq == 0)
            {
                cout << "\nSimulation is at generation T = " << k+1 << endl;
                
                if (save_bool == true)
                {
                    save(k+1);
                }
                
                
            }   
        }
    }
}


void MPopulation::calc_reprod_prob(int time)
{

    if (time < 10)
    {
        this->selection_prob = this->fitness;
        old_phen.push_back(phenotypes);
    }
    else if (time < 100 && !stable)
    {
        old_phen.push_back(phenotypes);
        
        int cont = 0;
        for (int i = 0; i < pop_size; i++)
        {
            if (phi(phenotypes, old_phen, time, i) >= epsilon)
                {this->selection_prob[i] = 0;}
            else
            {
                this->selection_prob[i] = this->fitness[i];
                cont++;
            }
        }

        if (cont >= 2)
        {
            cout << "Population became stable in time T = " << time << endl;
            this->stable = true;
        }
        else
        {
            this->selection_prob = this->fitness;
        }
          
    }
    else
    {
        this->selection_prob = this->fitness;  
    }
}

void MPopulation::select_parents()
{
    //update generator seed
    //this->seed = chrono::system_clock::now().time_since_epoch().count();
    //this->generator.seed(seed); 
    
    //converting probabilities to integers
    vector<int> int_selection_prob(pop_size);
    int multiplier = 1000;
    double corrector = 1e-5;

    for (int i = 0; i < pop_size; i++)
    {
        int_selection_prob[i] = int(multiplier*(selection_prob[i] + corrector));
    }
    
    discrete_distribution<int> disc_dist(int_selection_prob.begin(), int_selection_prob.end());

    parents.resize(pop_size/2, vector<int>(2));

    for (int i = 0; i < parents.size(); i++)
    {
        int couple_0 = disc_dist(generator);
        int couple_1 = disc_dist(generator);

        while (couple_0 == couple_1)
        {
            couple_0 = disc_dist(generator);
            couple_1 = disc_dist(generator);
        }
        
        parents[i][0] = couple_0;
        parents[i][1] = couple_1;
    }
}

void MPopulation::draw_new_population()
{
    //update generator seed
    this->seed = chrono::system_clock::now().time_since_epoch().count();
    this->generator.seed(seed); 

    vector<int> line_selection = {1,1};
    vector<Eigen::MatrixXd> temp_population(population);
    Eigen::MatrixXd temp_matrix(m_size, m_size);

    discrete_distribution<int> disc_dist(line_selection.begin(), line_selection.end());

    int cont = 0;
    for (int i = 0; i < parents.size(); i++)
    {
        for (int k = 0; k < 2; k++)
        {
            for (int j = 0; j < temp_matrix.rows(); j++)
            {
                int chosen_ind = disc_dist(generator);
                int chosen = parents[i][chosen_ind];
                temp_matrix.row(j) = population[chosen].row(j);
            }

            temp_population[cont] = temp_matrix;
            cont++;
        }
    }
    
    //mutations
    //if (!dist_config.compare("normal"))
    //{
    #include "normal_mutation.H"
    //}
    
    //temporary mutated population`s phenotype
    Eigen::MatrixXd temp_phen(m_size,pop_size);
    Eigen::MatrixXd aux(m_size,1);
    vector<int> temp_prob(pop_size);
    multiplier = 1000;
    corrector = 1e-5;

    for (int i = 0; i < pop_size; i++)
    {
        aux = calc_phenotype(temp_population[i], phenotypes.col(i));
        temp_phen.col(i) = aux;
        temp_prob[i] = int (multiplier*(fit(temp_phen.col(i)) + corrector));
    }

    vector<int> selected(pop_size);
    discrete_distribution<int> another_dist(temp_prob.begin(), temp_prob.end());
    for (int i = 0; i < pop_size; i++)
    {
        selected[i] = another_dist(generator);
        population[i] = temp_population[selected[i]];
    }
}

MPopulation::~MPopulation()
{
}