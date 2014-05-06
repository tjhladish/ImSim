#include "MSMSC.h"
#include <time.h>
#include <stdlib.h>
#include <cmath>

int j_max; // number of different networks or (repetitions on same network) to iterate through
bool reuse_net; // whether to reuse networks or generate a new one for each j 
int s_max;   // number of times to try each network
int n;     // size of network to use
int patient_zero_total;// number of infections to introduce, divided (and rounded) among subtypes as needed
double R_zero; //r0 value for epidemic
 
double CHI;
string RunID;

void connect_network (Network* net) {
    vector<double> dist;
    double deg_array[] = {0, 3, 45, 160, 393, 715, 883, 989, 897, 752,
                          697, 755, 856, 1085, 1224, 1452, 1578, 1622, 1711, 1584,
                          1514, 1355, 1209, 990, 805, 597, 477, 353, 232, 177,
                          126, 90, 69, 54, 36, 29, 21, 16, 10, 5,
                          8, 13, 7, 9, 3, 1, 3, 3, 2, 5,
                          0, 1, 0, 0, 0, 1};
    dist.assign(deg_array,deg_array+56);

    dist = normalize_dist(dist, sum(dist));
    net->rand_connect_user(dist);
}

map< int, vector<double> > import_historical_subtypes(const string &filename) {
    // This function assumes 3 possible reported subtypes (h1n1, h3n2, b)
    ifstream iss(filename.c_str());
    if (!iss) {
        cerr << "ERROR: " << filename << " not found." << endl;
    }

    map< int, vector<double> > subtype_by_season;

    while (iss) {
        char buffer[500];
        iss.getline(buffer,500);
        istringstream line(buffer);
        int season;
        float s0, s1, s2;
        if (line >> season >> s0 >> s1 >> s2 ) {
            vector<double> values(NumStrains);
            values[0] = s0; values[1] = s1; values[2] = s2;
            subtype_by_season[season] = values;
        } else {
            cerr << "Unparseable content: " << buffer << endl;
        }
        
    }
    iss.close();
    return subtype_by_season;
}


vector<int> sample_multinomial(vector<double> weights, int trials, MTRand* mtrand) {       
    vector<int> sample(trials, 0);
    for (int i = 0; i < trials; ++i) {
        double r = mtrand->rand();
        for (unsigned int j = 0; j < weights.size(); ++j) {
            if (weights[j] > r) {
                ++sample[j];
            } else {
                r -= weights[j];
            }
        }
    } 
    return sample;
}


void processCmdlineParameter(int argc, char* argv[] ) {
    cerr << "Arguments provided: " << argc - 1 << endl;
    cerr << "Argument order: reps seasons net_size R0 CHI p0_ct RunID\n";
    if ( argc == 8 ) {  
        j_max     = (int) atoi(argv[1]);
        s_max     = (int) atoi(argv[2]); 
        n         = (int) atoi(argv[3]); 
        R_zero    = (double) atof(argv[4]);
        CHI       = (double) atof(argv[5]);
        patient_zero_total = (int) atoi(argv[6]);
        RunID     = argv[7];
 
        cout << "# "; for(int i=0; i < argc; i++ )  cout << argv[i] <<  " ";  cout << endl;
        //need to check bounds
        if (CHI < 0)        { cerr << "Expecting positive immune decay parameter\n"; exit(-1); }

    } else {
        cerr << "Expecting 7 parameters.\n";
        exit(-1);
    }
}


int main(int argc, char* argv[]) {
    processCmdlineParameter(argc,argv);

    // Header line
    cout << "# RunID Network Season Epi_size P0_size R0\n";
    if (reuse_net == true) {
        cerr << "Same network topology will be used for all repetitions.\n";
        cerr << "Generating network toplogy . . . ";
    }
    
    const int first_season = 1992; // should agree with strains.txt file
    map< int, vector<double> > subtype_by_season = import_historical_subtypes("./raw_data/grog/strains.txt");
    
    //for (int s = 1992; s<2013; ++s) {
    //    cout << s << " "; for ( int i = 0; i<3; ++i ) cout << subtype_by_season[s][i] << " ";
    //    cout << endl;
    //}
    //return -5;

    Network* net = new Network("EpiNet", Network::Undirected);
    net->populate(n);
    connect_network(net);
    MSMSC_Sim* sim = new MSMSC_Sim(net,CHI);
    vector<double> r_zeros(3, R_zero);
    sim->calc_naive_transmissibility(r_zeros); //calculates correct T for this network topology
    vector<Node*> patients_zero;
    double new_R_zero = R_zero;
    for ( int j = 0; j < j_max; j++) {
        cerr << "Repetition: " << j << endl;
        //double Tc_actual = sim->calc_critical_transmissibility();

        for ( int season = 0; season < s_max; season++) {
            const int cluster = season + 1;
            for (int s = 0; s < (int) NumStrains; ++s) {
                sim->set_cluster( (StrainType) s, cluster);
            }
            vector<int> patient_zero_cts = sample_multinomial(subtype_by_season[first_season + season], patient_zero_total, net->get_rng());
            patients_zero = sim->rand_infect(patient_zero_cts);
            sim->run_simulation();
            cout << RunID << " " << j << " " << first_season + season << " | ";
            vector<int> final_sizes = sim->final_size_by_strain();
            cout_vector(final_sizes);
            cout << " | " << sim->epidemic_size() ; 
//            cout_vector(patient_zero_cts);
            //cout << " " << new_R_zero << " ";
            
            // now calculate what R_zero will be at the start of the next season
            //vector<double> average_tk;
            //double average_t = 0;
            //sim->calculate_average_transmissibility(average_tk, average_t);
            //new_R_zero = average_t / Tc_actual;
            cout << endl;
            //run_mixture_model(mixture, net);
            //sim.summary();
        }
        //sim->write_event_history("/tmp/eventHistory");
        sim->reset();
        //new_R_zero = R_zero;
    }
}
