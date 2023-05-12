
vector<string> options = {"-help","-dist","-dim","-nsim","-gens","-cond","-savepop","-savefreq"};

if (argc > 1)
{
    for (int i = 1; i < argc; i += 2)
    {
        string aux = argv[i];
        if (!aux.compare("-help"))
        {
            help_output();
            exit(EXIT_FAILURE);
        }
        else if (!aux.compare("-dist"))
        {
            dist = argv[i+1];
        }
        else if (!aux.compare("-dim"))
        {
            dim = atoi(argv[i+1]);
        }
        else if (!aux.compare("-nsim"))
        {
            runs = atoi(argv[i+1]);
        }
        else if (!aux.compare("-gens"))
        {
            gens = atoi(argv[i+1]);
        }
        else if (!aux.compare("-mean"))
        {
            mean = atof(argv[i+1]);
        }
        else if (!aux.compare("-stddev"))
        {
            std_dev = atof(argv[i+1]);
        } 
        else if (!aux.compare("-cond"))
        {
            cond = atoi(argv[i+1]);
        }
        else if (!aux.compare("-savepop"))
        {
            save_bool = atoi(argv[i+1]);
        }
        else if (!aux.compare("-savefreq"))
        {
            save_freq = atof(argv[i+1]);
        }
        else
        {
            cout << "Invalid option " << aux << ", valid options are: ";
            for (int i = 0; i < options.size(); i++)
            {
                cout << options[i] << " ";
            }
            cout << endl << "Aborting execution!" << endl;
            exit(EXIT_FAILURE);
        }
    }
}