void MPopulation::setFileSuffix(string suffix)
{
    this->file_suffix = suffix;
}

void MPopulation::setSavingFrequency(bool y_or_n, int frequency)
{
    this->save_freq = frequency;
    this->save_bool = y_or_n;
    if (save_bool)
    {
        cout << "\nSaving genotypes and phenotypes every " << save_freq << " generations!" << endl;
    }
    else
    {
        cout << "\nOnly mean phenotypes and population fitness being saved at the end of each experiment!" << endl;
    }
}

void MPopulation::save(int time)
{
    int precision = 8;
    Eigen::IOFormat output(precision, 0, " ", "\n","","","","\n");

    string file_name = "outputs/";
    file_name += "gen";
    file_name += to_string(time);
    if (!file_suffix.empty())
    {
        file_name += "_";
        file_name += this->file_suffix;
    }
    file_name += ".dat";
    ofstream myFile(file_name);

    cout << "Saving matrices into `" << file_name << "` file" << endl;
    for (int i = 0; i < pop_size; i++)
    {
        myFile << population[i].format(output) << endl;
    }

    file_name = "outputs/";
    file_name += "phen";
    file_name += to_string(time);
    if (!file_suffix.empty())
    {
        file_name += "_";
        file_name += this->file_suffix;
    }
    file_name += ".dat";
    ofstream myFile2(file_name);

    cout << "Saving phenotypes into `" << file_name << "` file" << endl;
    for (int i = 0; i < phenotypes.cols(); i++)
    {
        myFile2 << phenotypes.col(i).format(output) << endl;
    }    
}


void MPopulation::save_means()
{
    int precision = 8;
    Eigen::IOFormat output(precision, 0, " ", "\n","","","","\n");

    string file_name = "outputs/";
    file_name += "meanf";
    if (!file_suffix.empty())
    {
        file_name += "_";
        file_name += this->file_suffix;
    }
    file_name += ".dat";
    ofstream myFile(file_name);
    Eigen::VectorXd fitnesses(mean_fitness.size());

    for (int i = 0; i < fitnesses.size(); i++)
    {
        fitnesses(i) = mean_fitness[i];
    }
    
    cout << "\nSaving mean fitnesses into `" << file_name << "` file" << endl;
    myFile << fitnesses.format(output) << endl;

    file_name = "outputs/";
    file_name += "meanph";
    if (!file_suffix.empty())
    {
        file_name += "_";
        file_name += this->file_suffix;
    }
    file_name += ".dat";
    ofstream myFile2(file_name);

    Eigen::MatrixXd meanps(mean_phenotype.size(), mean_phenotype[0].size());

    for (int i = 0; i < meanps.rows(); i++)
    {
        for (int j = 0; j < meanps.cols(); j++)
        {
            meanps(i,j) = mean_phenotype[i][j];
        }
    }
    
    cout << "Saving mean phenotypes into `" << file_name << "` file" << endl;
    myFile2 << meanps.format(output) << endl;

    if (get_condition)
    {

        file_name = "outputs/";
        file_name += "meancond";
        if (!file_suffix.empty())
        {
            file_name += "_";
            file_name += this->file_suffix;
        }
        file_name += ".dat";
        ofstream myFile3(file_name);

        Eigen::VectorXd conditions(condition_n.size());
        for (int i = 0; i < condition_n.size(); i++)
        {
            conditions(i) = condition_n[i];
        }

        cout << "Saving mean phenotypes into `" << file_name << "` file" << endl;
        myFile3 << conditions.format(output) << endl;
    } 
}