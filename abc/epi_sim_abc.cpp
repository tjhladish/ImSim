#include "../MultiSeason_Network/MultiSeason_Sim.h"
#include <time.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <map>


struct Parameters {
    string observed_data_file;
    int network_size, burnin, sample_size, threads;
    double epsilon, obs_mean, obs_median, obs_max, obs_range, obs_sd, obs_skew, obs_ss, obs_ll, obs_sl, obs_ls;
};

struct Particle { 
    double R0, Ih, h; 
    int P0;
};

struct Doubled_variances {
    double R0, Ih, P0, h;
};

/*
struct ParticleSet { 
    vector<Particle> particles;
    Particle back() { return particles.back(); }
    void resize(const int n) { particles.resize(n); }
    int size() { return (int) particles.size(); }
    Particle& operator[] ( const int nIndex) { return particles[nIndex]; }
};
*/

//vector<ParticleSet> theta; 

// Define ABC parameters
const double R0_min = 1.0;
const double R0_max = 8.0;
const double Ih_min = 1.0/12.0;
const double Ih_max = 100.0;
const int P0_min = 1;
const int P0_max = 256;
const double h_min  = 0.0;
const double h_max  = 1.0;

// Don't consider values outside of these bounds
// (because they would be nonsensical or would cause simulator to crash)
const double R0_hard_min = 0.0;
const double R0_hard_max = 1000.0;
const double Ih_hard_min = 0.00001;
const double Ih_hard_max = 10000.0;
const int P0_hard_min = 1;
const int P0_hard_max = 1000;
const double h_hard_min  = 0.0;
const double h_hard_max  = 1.0;

int NET_SIZE;     // size of network to use
int patient_zero_ct; // number of infections to introduce (usually 1)
double R_zero; //r0 value for epidemic

enum DistType { POI, EXP, POW, URB, CON}; // specifies which degree distribution to use: poisson, exponential, powerlaw, urban, or constant
DistType dist = URB;

double Ih; // Immunity halflife
int burnin;
double h; // hospitalization/doctor visit rate
string observed_data_filename;

void generate_network(Network* net, MultiSeason_Sim* sim);
void connect_network (Network* net);
void load_observed_data(string filename, char sep, map<string, vector<float> > &data);
double KS(vector<float>s1, vector<float>s2);
void report_metrics(map<string, vector<float> > s1, map<string, vector<float> > s2);
//double autocorrelation_score( vector<float>s1, vector<float>s2 );
vector<double> autocorrelation_matrix( map<string, vector<float> >& data );
Doubled_variances calculate_doubled_variances( vector<Particle> particle_set );
vector<double> calculate_predictive_prior_weights(vector<Particle> predictive_prior, int set_num); 

//float obs_R0 = 3;   // estimated from Vynnycky et al 2007. Estimates of the reproduction numbers of Spanish influenza using morbidity data
//float obs_R  = 1.3; // from Chowell et al 2007. Seasonal Influenza in the United States, France, and Australia: Transmission and prospects for control

float controlB[] = {1.26, 0.34, 0.70, 1.75, 50.57, 1.55, 0.08, 0.42, 0.50, 3.20, 0.15, 0.49, 0.95, 0.24, 1.37, 0.17, 6.98, 0.10, 0.94, 0.38 };
float treatmentB[] = {2.37, 2.16, 14.82, 1.73, 41.04, 0.23, 1.32, 2.91, 39.41, 0.11, 27.44, 4.51, 0.51, 4.50, 0.18, 14.68, 4.66, 1.30, 2.06, 1.19};
double uniform_pdf(double a, double b) { return 1.0 / fabs(b-a); }

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
    vector<double> dist;
    double deg_array[] = {0, 0, 1, 12, 45, 50, 73, 106, 93, 74, 68, 78, 91, 102, 127, 137, 170, 165, 181, 181, 150, 166, 154, 101, 67, 69, 58, 44, 26, 24, 17, 6, 11, 4, 0, 6, 5, 3, 1, 1, 3, 1, 1, 0, 1, 0, 2};
    dist.assign(deg_array,deg_array+47);
    dist = normalize_dist(dist, sum(dist));
    net->rand_connect_user(dist);
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

void process_config_file (string config_filename, Parameters &p) {

    std::map<string,string> tmp_par;
    ifstream myfile(config_filename.c_str());
    std::stringstream ss;

    if (myfile.is_open()) {
        string line;
        char sep = '\t';

        while ( getline(myfile, line) ) {
            vector<string> fields;
            split(line, sep, fields);

            const char whitespace[] = " \n\t\r";
            string key =   strip( fields[0], whitespace ); 
            string value = strip( fields[1], whitespace ); 

            tmp_par[ key ] = value;
        }
    } else {
        cerr << "Failed to open parameter file\n";
        exit(102);
    }

    p.observed_data_file = tmp_par["observed_data_file"];
    p.network_size = to_int( tmp_par["network_size"] );
    p.burnin       = to_int( tmp_par["burnin"] );
    p.sample_size  = to_int( tmp_par["sample_size"] );  
    p.epsilon      = to_double( tmp_par["epsilon"] );
    p.threads      = to_int( tmp_par["threads"] );
    p.obs_mean     = to_double( tmp_par["obs_mean"] );
    p.obs_median   = to_double( tmp_par["obs_median"] );
    p.obs_max      = to_double( tmp_par["obs_max"] );
    p.obs_range    = to_double( tmp_par["obs_range"] );
    p.obs_sd       = to_double( tmp_par["obs_sd"] );
    p.obs_skew     = to_double( tmp_par["obs_skew"] );
    p.obs_ss       = to_double( tmp_par["obs_ss"] );
    p.obs_ll       = to_double( tmp_par["obs_ll"] );
    p.obs_sl       = to_double( tmp_par["obs_sl"] );
    p.obs_ls       = to_double( tmp_par["obs_ls"] );

    return;
}


bool fileExists(string strFilename) { 
    struct stat stFileInfo; 
    bool blnReturn; 
    int intStat; 
    // Attempt to get the file attributes 
    intStat = stat(strFilename.c_str(),&stFileInfo); 
    if(intStat == 0) { 
        blnReturn = true; 
    } else { 
        blnReturn = false;
    }
    return(blnReturn);
}

vector<Particle> read_predictive_prior_file(int set_num) {
    string filename = "predictive_prior." + to_string(set_num);
// 'R0', 'Ih', 'h', 'P0', 'mean', 'median', 'max', 'range', 'sd', 'skew', 'ss', 'll', 'sl', 'ls', 'Re'   
    ifstream myfile(filename.c_str());
    std::stringstream ss;
    vector<Particle> predictive_prior;

    if (myfile.is_open()) {
        string line;
        char sep = ' ';

        while ( getline(myfile, line) ) {
            vector<string> fields;
            split(line, sep, fields);
            Particle p;

            const char whitespace[] = " \n\t\r";
            string R0_str = strip( fields[0], whitespace ); 
            string Ih_str = strip( fields[1], whitespace ); 
            string h_str  = strip( fields[2], whitespace ); 
            string P0_str = strip( fields[3], whitespace ); 

            p.R0 = to_double( R0_str );
            p.Ih = to_double( Ih_str );
            p.h  = to_double( h_str );
            p.P0 = to_int( P0_str );

            predictive_prior.push_back(p);
        }
    } else {
        cerr << "Failed to open parameter file\n";
        exit(103);
    }
    return predictive_prior;
}

vector<double> read_weights_file(int set_num) {
    string filename = "weights." + to_string(set_num);

    vector<double> weights;
    ifstream myfile(filename.c_str());
    std::stringstream ss;

    if (myfile.is_open()) {
        string line;

        while ( getline(myfile,line) ) {
            //split string based on "," and store results into vector
            const char whitespace[] = " \n\t\r";
            string val_str = strip( line, whitespace );
        
            weights.push_back( atof( val_str.c_str() ) );
        }
    }
    return weights;
}


void write_weights_file(vector<double> weights, int set_num) {
    string filename = "weights." + to_string(set_num);

    ofstream pipe(filename.c_str(), ios::out);
    for (unsigned int i = 0; i < weights.size(); i++) {
        pipe << weights[i] << endl;
    }
    pipe.close();
    return;
}


int determine_set_number() {
    bool found_file = true;
    int set_num = 0;
    while (found_file) {
        stringstream ss;
        ss << "predictive_prior." << set_num;
        if ( fileExists(ss.str()) ) {
            set_num++;
        } else {
            found_file = false;
        }
    }
    return set_num;
}


void sample_prior(vector<Particle> &particles, int sample_size, MTRand* mtrand) {
    particles.resize(sample_size);
    for (int i = 0; i < sample_size; i++) {
        particles[i].R0 = rand_uniform( R0_min, R0_max, mtrand );
        particles[i].Ih = rand_uniform( Ih_min, Ih_max, mtrand );
        particles[i].h  = rand_uniform( h_min,  h_max,  mtrand  );
        particles[i].P0 = rand_uniform_int( P0_min, P0_max, mtrand );
    }
}

double rand_trunc_normal(double mu, double sigma_squared, double min, double max, MTRand* mtrand) {
    double sigma = sqrt(sigma_squared);
    // Don't like this, but it should work
    while (1) {
        double dev = mtrand->randNorm(mu, sigma);
        if (dev >= min and dev <= max) {
            return dev;
        }
    }
}
                                

void sample_predictive_prior(int set_num, vector<Particle> &particles, int sample_size, MTRand* mtrand) {
    // Read in predictive prior from the last set.
    // We need to calculate the proper weights for the predictive prior so that we know how to sample from it.
    vector<Particle> predictive_prior = read_predictive_prior_file( set_num - 1 );
    Doubled_variances sampling_variance = calculate_doubled_variances( predictive_prior );
    
    vector<double> weights;
    if (set_num == 1) {
        // uniform weights for set 0 predictive prior
        weights.resize(predictive_prior.size(), 1.0/(double) predictive_prior.size());
    } else if ( set_num > 1 ) {
        // weights from set - 2 are needed to calculate weights for set - 1
        weights = calculate_predictive_prior_weights( predictive_prior, set_num );
    }
    write_weights_file( weights, set_num - 1 );

    particles.resize(sample_size);
    for (int i = 0; i < sample_size; i++) {
        // Select a particle index j to use from the predictive prior
        int j = rand_nonuniform_int( weights, mtrand );
        particles[i].R0 = rand_trunc_normal( predictive_prior[ j ].R0, sampling_variance.R0, R0_hard_min, R0_hard_max, mtrand ); 
        particles[i].Ih = rand_trunc_normal( predictive_prior[ j ].Ih, sampling_variance.Ih, Ih_hard_min, Ih_hard_max, mtrand ); 
        particles[i].h  = rand_trunc_normal( predictive_prior[ j ].h,  sampling_variance.h,  h_hard_min,  h_hard_max,  mtrand ); 
        particles[i].P0 = (int) (rand_trunc_normal( predictive_prior[ j ].P0, sampling_variance.P0, P0_hard_min, P0_hard_max, mtrand ) + 0.5);
    }
    return;
}

/*
vector<double> determine_weights(int particles_in_prior) {
    const double N = (double) particles_in_prior;
    vector<double> weights;
    int set_num = determine_set_number();
    switch (set_num) {
        case 0:
            break;
        case 1:
            weights.clear();
            weights.resize((int) N, 1.0/N);
            break;
        default:
            cerr << "determine weights is not fully implemented\n";
            break;
    }
        
    return weights;
}*/


Doubled_variances calculate_doubled_variances( vector<Particle> particle_set ) {
    int set_size = particle_set.size();
    vector<double> R0_vec( set_size );
    vector<double> Ih_vec( set_size );
    vector<double> h_vec( set_size );
    vector<double> P0_vec( set_size );

    for (int i = 0; i < set_size; i++) {
        R0_vec[i] = particle_set[i].R0;
        Ih_vec[i] = particle_set[i].Ih;
        h_vec[i]  = particle_set[i].h;
        P0_vec[i] = particle_set[i].P0;
    }

    Doubled_variances doubled_var;
    doubled_var.R0 = 2 * variance( R0_vec );
    doubled_var.Ih = 2 * variance( Ih_vec );
    doubled_var.P0 = 2 * variance( P0_vec );
    doubled_var.h  = 2 * variance( h_vec );

    return doubled_var;
}


vector<double> calculate_predictive_prior_weights(vector<Particle> predictive_prior, int set_num) {

    const int prior_size = predictive_prior.size();
    vector<double> predictive_prior_weights( prior_size );
    const double numerator = uniform_pdf(R0_min, R0_max)
                           * uniform_pdf(Ih_min, Ih_max)
                           * uniform_pdf(P0_min, P0_max)
                           * uniform_pdf(h_min, h_max);

    vector <double>   old_predictive_prior_weights = read_weights_file( set_num - 2 ); 
    vector<Particle>  old_predictive_prior         = read_predictive_prior_file( set_num - 2 );
    const int old_prior_size                       = old_predictive_prior.size();
    Doubled_variances old_par_2var                 = calculate_doubled_variances( old_predictive_prior );


    for (int i = 0; i < prior_size; i++) {
        double denominator = 0.0;
        for (int j = 0; j < old_prior_size; j++) {
            denominator += old_predictive_prior_weights[j]
                           //double normal_pdf(double x, double mu, double var) {
                           * normal_pdf(predictive_prior[i].R0, old_predictive_prior[j].R0, old_par_2var.R0)
                           * normal_pdf(predictive_prior[i].Ih, old_predictive_prior[j].Ih, old_par_2var.Ih)
                           * normal_pdf(predictive_prior[i].P0, old_predictive_prior[j].P0, old_par_2var.P0)
                           * normal_pdf(predictive_prior[i].h,  old_predictive_prior[j].h,  old_par_2var.h);
/*
                 denominator = 1; // DEBUGGING, DELETE LATER
                 cerr << predictive_prior[i].P0 
                 << " " << old_predictive_prior[j].P0 
                 << " " << old_par_2var.P0 
                 << " " << normal_pdf(predictive_prior[i].P0, old_predictive_prior[j].P0, old_par_2var.P0)
                 << endl;
            cerr << old_predictive_prior_weights[j]
                 << " " << normal_pdf(predictive_prior[i].R0, old_predictive_prior[j].R0, old_par_2var.R0)
                 << " " << normal_pdf(predictive_prior[i].Ih, old_predictive_prior[j].Ih, old_par_2var.Ih)
                 << " " << normal_pdf(predictive_prior[i].P0, old_predictive_prior[j].P0, old_par_2var.P0)
                 << " " << normal_pdf(predictive_prior[i].h,  old_predictive_prior[j].h,  old_par_2var.h) << endl;
                 */
        }
        predictive_prior_weights[i] = numerator / denominator;
    }

    return normalize_dist( predictive_prior_weights );
}


/*
Determine run number: what is the last predictive_prior.X file? (or does one not exist->first run)

If first time being run, no input files exist.  Sample from uniform priors, output particles.

If second time, read in predictive prior from first run, use uniform weights. Output uniform
weights as weights for first run, output second run particles.

If third+ time, read in predictive prior from last run, weights from second to last run,
calculate and output weights for last run; output particles.

*/


int main(int argc, char* argv[]) {
    if (argc < 2) { cerr << "Provide configuration file name as a parameter\n"; }
    string config_file_name( argv[1] );

    // Process configuration file
    Parameters par;
    process_config_file(config_file_name, par);

    // Determine which ABC set we're supposed to run
    const int set_num = determine_set_number();

    // Determine what parameter combinations we will try
    static MTRand mtrand;
    vector<Particle> particles;
    if (set_num == 0) {
        // sample naive priors
        sample_prior(particles, par.sample_size, &mtrand);
    } else {
        // sample predictive priors based on last set
        sample_predictive_prior(set_num, particles, par.sample_size, &mtrand);
    }

    for (unsigned int i = 0; i < particles.size(); i++) {
        cerr << particles[i].R0 << " " << particles[i].Ih << " " << particles[i].h << " " << particles[i].P0 << endl;
    }

exit(0);

    map<string, vector<float> > sim_data;
    map<string, vector<float> > R0_vals;
    map<string, vector<float> > obs_data; //string = location, float = incidence on [0,1]

    load_observed_data(par.observed_data_file, ',', obs_data);

    Network* net = new Network("EpiNet", Network::Undirected);
    MultiSeason_Sim* sim = new MultiSeason_Sim(net, Ih);
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

void load_observed_data(string filename, char sep, map<string, vector<float> > & obs_data) {
    //cerr << "Reading observed data\n";
    ifstream myfile(filename.c_str());

    if (myfile.is_open()) {
        string line;

        while ( getline(myfile,line) ) {
            //split string based on "," and store results into vector
            vector<string> fields;
            split(line, sep, fields);
            const char whitespace[] = " \n\t\r";

            //format check
            if (fields.size() > 3 ) {
                cerr << "Skipping line: too many fields: " << line << endl;
                continue;
            } else if (fields.size() < 3 ) {
                cerr << "Skipping line: too few fields: " << line << endl;
                continue;
            } else { 
                string loc   = strip(fields[0],whitespace);
                string year  = strip(fields[1],whitespace);
                float  count = to_float(strip(fields[2],whitespace));
                obs_data[loc].push_back(count); 
            }
        }
    }
    return;
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

void report_metrics(map<string, vector<float> > sim, map<string, vector<float> > Re_values) {
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

    return;
}


