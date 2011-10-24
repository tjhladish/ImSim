#include "/home/tjhladish/work/ImSim/MultiSeason_Network/MultiSeason_Sim.h"
#include <time.h>
#include <stdlib.h>

int NET_SIZE;     // size of network to use
int patient_zero_ct; // number of infections to introduce (usually 1)
double R_zero; //r0 value for epidemic

enum DistType { POI, EXP, POW, URB, CON}; // specifies which degree distribution to use: poisson, exponential, powerlaw, urban, or constant
DistType dist = URB;

// Distribution parameters.  param2 is a dummy for poisson and exponential.
double param1;
double param2;
double Ih; // Immunity halflife
int burnin;
double h; // hospitalization/doctor visit rate
string observed_data_filename;

void generate_network(Network* net, MultiSeason_Sim* sim);
void connect_network (Network* net);
void read_observed_data(string filename, char sep, map<string, vector<float> > &data);
double KS(vector<float>s1, vector<float>s2);
vector<double> report_metrics(map<string, vector<float> > s1, map<string, vector<float> > s2);
//double autocorrelation_score( vector<float>s1, vector<float>s2 );
vector<double> autocorrelation_matrix( map<string, vector<float> >& data );

//float obs_R0 = 3;   // estimated from Vynnycky et al 2007. Estimates of the reproduction numbers of Spanish influenza using morbidity data
//float obs_R  = 1.3; // from Chowell et al 2007. Seasonal Influenza in the United States, France, and Australia: Transmission and prospects for control

float controlB[] = {1.26, 0.34, 0.70, 1.75, 50.57, 1.55, 0.08, 0.42, 0.50, 3.20, 0.15, 0.49, 0.95, 0.24, 1.37, 0.17, 6.98, 0.10, 0.94, 0.38 };
float treatmentB[] = {2.37, 2.16, 14.82, 1.73, 41.04, 0.23, 1.32, 2.91, 39.41, 0.11, 27.44, 4.51, 0.51, 4.50, 0.18, 14.68, 4.66, 1.30, 2.06, 1.19};

vector<float> flatten_map(map<string, vector<float> > data) {
    vector<float> flat;
    map<string, vector<float> >::iterator it;

    for ( it = data.begin(); it != data.end(); it++ ) {
        string loc = it->first;
        for (unsigned int i = 0; i < data[loc].size(); i++) {
            flat.push_back(data[loc][i]);
        }
    }
    return flat;
}

void connect_network (Network* net) {
    if (dist == POI) {
        net->rand_connect_poisson(param1);
    } else if (dist == EXP) {
        net->rand_connect_exponential(param1);
    } else if (dist == POW) {
        net->rand_connect_powerlaw(param1, param2);
    } else if (dist == URB) {
        vector<double> dist;
        double deg_array[] = {0, 0, 1, 12, 45, 50, 73, 106, 93, 74, 68, 78, 91, 102, 127, 137, 170, 165, 181, 181, 150, 166, 154, 101, 67, 69, 58, 44, 26, 24, 17, 6, 11, 4, 0, 6, 5, 3, 1, 1, 3, 1, 1, 0, 1, 0, 2};
        dist.assign(deg_array,deg_array+47);
        dist = normalize_dist(dist, sum(dist));
        net->rand_connect_user(dist);
    } else if (dist == CON) {
        vector<double> dist(param1+1, 0);
        dist[param1] = 1;
        net->rand_connect_user(dist);
    }
}

double calculate_mean_effective_R (vector<float> epi_values, vector<float> Re_values, float threshold) {
    int n = 0;
    double running_sum = 0.0;
    for ( unsigned int i = 0; i < epi_values.size(); i++) {
        if (Re_values[i] > 1) {
        //if (epi_values[i] > threshold) {
            running_sum += Re_values[i];
            n++;
        }
    }
    if ( n == 0 ) {
        return 0.0;
    } else {
        return running_sum / n;
    }
}

void processCmdlineParameter(int argc, char* argv[] ) {
    //cerr << "Arguments provided: " << argc - 1 << endl;
    //cerr << "Argument order: burnin_seasons net_size R0 Ih p0_ct hosp_rate real_data_filename\n";
    if ( argc == 8 ) {  

        burnin     = (int) atoi(argv[1]); 
        NET_SIZE   = (int) atoi(argv[2]); 
        R_zero     = (double) atof(argv[3]);
        dist       = (DistType) URB;
        Ih         = (double) atof(argv[4]);
        patient_zero_ct =  (int) atoi(argv[5]);
        h          = (double) atof(argv[6]);

        observed_data_filename = argv[7];
 
        //cerr << "# "; for(int i=0; i < argc; i++ )  cerr << argv[i] <<  " ";  cerr << endl;
        //need to check bounds
        if (Ih < 0)        { cerr << "Expecting non-negative immunity halflife parameter\n"; exit(-1); }
        if (h <= 0 || h > 1)        { cerr << "Expecting h to be between 0 and 1 \n"; exit(-1); }

    } else {
        cerr << "Expecting 7 parameters.\n";
        cerr << "Argument order: burnin_seasons net_size R0 Ih p0_ct hosp_rate real_data_filename\n";
        exit(-1);
    }
}


int main(int argc, char* argv[]) {
    processCmdlineParameter(argc,argv);

    map<string, vector<float> > sim_data;
    map<string, vector<float> > R0_vals;
    map<string, vector<float> > obs_data; //string = location, float = incidence on [0,1]
    read_observed_data(observed_data_filename, ',', obs_data);

    // Header line
    //cerr << "# Network Season Epi_size P0_size R0\n";

    Network* net = new Network("EpiNet", Network::Undirected);
    MultiSeason_Sim* sim = new MultiSeason_Sim(net,Ih);
    generate_network(net, sim);
    double new_R_zero = R_zero;

    map<string, vector<float> >::iterator it;

    for ( it = obs_data.begin(); it != obs_data.end(); it++ ) {
        double Tc_actual = sim->calc_critical_transmissibility();
        string loc = (*it).first;
        const int obs_N = obs_data[loc].size();

        sim_data[loc].resize(obs_N);
        R0_vals[loc].resize(obs_N);

        for ( int season = 0; season < burnin + obs_N; season++) {
            sim->rand_infect(patient_zero_ct);
            sim->run_simulation();
                                    
            if (season >= burnin) {
                // epi size in percent, reduced by hospitalization factor
                const double transmitted_size = double(sim->epidemic_size() - patient_zero_ct)/(NET_SIZE - patient_zero_ct);
                sim_data[loc][season - burnin] =  h * transmitted_size;
                R0_vals[loc][season - burnin]  = new_R_zero;
                //cerr << "* ";
            }
            //cerr << burnin << " " << season << " " << new_R_zero << endl;

            // now calculate what R_zero will be at the start of the next season
            vector<double> average_tk;
            double average_t = 0;
            sim->calculate_average_transmissibility(average_tk, average_t);
            new_R_zero = average_t / Tc_actual;
        }
        sim->reset();
        new_R_zero = R_zero;
    }
    // Report parameters
    cout << "R0 " << R_zero;
    cout << " Ih " << Ih;
    cout << " h " << h;
    cout << " P0 " << patient_zero_ct;
    report_metrics(sim_data, R0_vals);
    //vector<float> r0_flat = flatten_map(R0_vals); 
    //cerr_vector(r0_flat ); cerr << endl;
    return 0;
}

void generate_network(Network* net, MultiSeason_Sim* sim) {
    net->populate(NET_SIZE);
    connect_network(net); // connect network using the parameters above
    sim->calc_naive_transmissibility(R_zero); //calculates correct T for this network topology
}

void read_observed_data(string filename, char sep, map<string, vector<float> > & obs_data) {
    //cerr << "Reading observed data\n";
    ifstream myfile(filename.c_str());

    if (myfile.is_open()) {
        string line;

        while ( getline(myfile,line) ) {
            //split string based on "," and store results into vector
            vector<string> fields;
            split(line,sep, fields);
            const char whitespace[] = " \n\t\r";

            //format check
            if (fields.size() > 3 ) {
                cerr << "Skipping line: too many fields: " << line << endl;
                continue;
            } else { 
                string loc   = strip(fields[0],whitespace);
                string year  = strip(fields[1],whitespace);
                float  count = to_float(strip(fields[2],whitespace));
                obs_data[loc].push_back(count); 
            }
        }
    }
}

double KS(vector<float>s1, vector<float>s2) {
    double D = 0;
    int s1size = s1.size();
    int s2size = s2.size();

    sort(s1.begin(), s1.end());
    sort(s2.begin(), s2.end());

    s1.push_back(numeric_limits<float>::infinity());
    s2.push_back(numeric_limits<float>::infinity());

    int i=0; int j=0;
    while( i < s1size-1 || j < s2size-1) {
        if (s1[i] < s2[j]) i++;
        else if (s1[i] > s2[j]) j++;
        else { i++; j++; }
        D = max(D, fabs(float(i)/s1size - float(j)/s2size));
     }
    return D;
}

vector<double> autocorrelation_matrix(map<string, vector<float> > &data) {
    map<string, vector<float> >::iterator it;
    double N = 0.0;
    
    // Determine total number of comparisons available
    for ( it = data.begin(); it != data.end(); it++ ) {
        string loc = it->first;
        if (data[loc].size() > 1) {
            N += data[loc].size() - 1; // n-1 per time series
        }
    }
    const double Inv_N = 1.0 / (double) N; // percent contributed by each comparison

    // Autocorrelation scores
    // indexing: 0 == small, small; 1 == large, large; 2 == small, large; 3 == large, small
    vector<double> ac_matrix(4,0);
    
    for ( it = data.begin(); it != data.end(); it++ ) {
        string loc = it->first;
        if (data[loc].size() <= 1) continue;

        const double avg = mean(data[loc]);
        for (unsigned int i = 0; i < data[loc].size() - 1; i++) {
            const float this_year = data[loc][i];
            const float next_year = data[loc][i+1];
            if ( this_year < avg and next_year < avg ) { // small, small
                ac_matrix[0] += Inv_N;
            } else if ( this_year >= avg and next_year >= avg ) { // large, large 
                ac_matrix[1] += Inv_N;
            } else if ( this_year < avg and next_year >= avg ) { // small, large 
                ac_matrix[2] += Inv_N;
            } else { // large, small
                ac_matrix[3] += Inv_N;
            }
        }
    }
    return ac_matrix;
}

vector<double> report_metrics(map<string, vector<float> > sim, map<string, vector<float> > Re_values) {
    vector<float> sim_flat = flatten_map( sim );
    vector<float> Re_flat  = flatten_map( Re_values );

    double mean_sim = mean(sim_flat);
    double median_sim = median(sim_flat);
    double skew_sim = mean_sim - median_sim;
    double sd_sim = stdev(sim_flat);
    vector<double> acm_sim = autocorrelation_matrix(sim);

    float epi_threshold = 0.05;
    double mean_Re = calculate_mean_effective_R (sim_flat, Re_flat, epi_threshold);

    cout << " mean " << mean_sim;
    cout << " median " << median_sim;
    cout << " sd " << sd_sim;
    cout << " skew " << skew_sim;
    cout << " ss " << acm_sim[0];
    cout << " ll " << acm_sim[1];
    cout << " sl " << acm_sim[2];
    cout << " ls " << acm_sim[3];
    cout << " Re " << mean_Re;
/*
    //cerr << "R0, Ih, h, P0 : mean, median, sd, skew, ss, ll, sl, ls : delta ";
    cerr << R_zero << "," << Ih << "," << h << "," << patient_zero_ct << ",  ";
    cerr << mean_sim << ",";
    cerr << median_sim << ",";
    cerr << sd_sim << ",";
    cerr << skew_sim << ",";
    cerr << acm_sim[0] << ",";
    cerr << acm_sim[1] << ",";
    cerr << acm_sim[2] << ",";
    cerr << acm_sim[3];
    cerr << endl;*/

    vector<double> variables(12);
/*    variables[0] = R_zero;
    variables[1] = Ih;
    variables[2] = h;
    variables[3] = patient_zero_ct;
    variables[4] = mean_sim;
    variables[5] = median_sim;
    variables[6] = sd_sim;
    variables[7] = skew_sim;
    variables[8] = acm_sim[0];
    variables[9] = acm_sim[1];
    variables[10] = acm_sim[2];
    variables[11] = acm_sim[3];
*/ 
    return variables;
}


