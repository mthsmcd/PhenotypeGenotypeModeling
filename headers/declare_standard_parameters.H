/*

    These are the standard parameters, some of them might be changed
    by terminal options.

*/

vector<double> target = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.8, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                        0.2, 0.2, 0.2, 0.2, 0.2, 0.2};
                    
vector<double> gene_expr = {1, 1, 1, 1, 1, 1, 1, 1000, 1000, 1000, 1000, 1000, 1000,
                            1, 1, 1, 1, 1, 1};

string dist = "normal";
double mean = 0.0;
double std_dev = 1.0;
int dim = 19;
int runs = 1;
int gens = 12000;
bool cond = false;
bool save_bool = true;
int save_freq = 1000;