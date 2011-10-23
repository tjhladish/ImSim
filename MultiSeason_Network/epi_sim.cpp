#include "MultiSeason_Sim.h"
#include <time.h>
#include <stdlib.h>

int j_max; // number of different networks or (repetitions on same network) to iterate through
bool reuse_net; // whether to reuse networks or generate a new one for each j 
int s_max;   // number of times to try each network
int n;     // size of network to use
int patient_zero_ct; // number of infections to introduce (usually 1)
double R_zero; //r0 value for epidemic
double m;   // fraction of individuals that migrate in natural mixture model
 
// specifies which degree distribution to use: poisson, exponential, powerlaw, urban, or constant
// EDGELIST means that an edgelist file will be used instead of randomly generating a network
enum DistType { POI, EXP, POW, URB, CON, EDGELIST};
enum MixType  { NO, HALF, FULL, NATURAL, HALF_SHUFFLE, FULL_SHUFFLE }; // specifies how much to shuffle the network between seasons
enum P0Type   { STATIC, DYNAMIC, CLUSTER }; // use a static number of randomly drawn nodes, use a dynamic number that depends 
                                            // on immune states of nodes, or a static number of clustered nodes

MixType mixture = NO;
DistType dist = POI;
P0Type P0_model = STATIC;

// Distribution parameters.  param2 is a dummy for poisson and exponential.
double param1;
double param2;
double Ih;
string RunID;
string EdgelistFilename;

void run_mixture_model(MixType mixture, Network* net);
void generate_network(Network* net, MultiSeason_Sim* sim);
void connect_network (Network* net);

void connect_network (Network* net) {
    if (dist == POI) {
        net->rand_connect_poisson(param1);
    } else if (dist == EXP) {
        net->rand_connect_exponential(param1);
    } else if (dist == POW) {
        net->rand_connect_powerlaw(param1, param2);
    } else if (dist == URB) {
        vector<double> dist;
        //double deg_array[] = {0, 0, 1, 12, 45, 50, 73, 106, 93, 74, 68, 78, 91, 102, 127, 137, 170, 165, 181, 181, 150, 166, 154, 101, 67, 69, 58, 44, 26, 24, 17, 6, 11, 4, 0, 6, 5, 3, 1, 1, 3, 1, 1, 0, 1, 0, 2};
        //dist.assign(deg_array,deg_array+47);

        double deg_array[] = {0, 3, 45, 160, 393, 715, 883, 989, 897, 752,
                              697, 755, 856, 1085, 1224, 1452, 1578, 1622, 1711, 1584,
                              1514, 1355, 1209, 990, 805, 597, 477, 353, 232, 177,
                              126, 90, 69, 54, 36, 29, 21, 16, 10, 5,
                              8, 13, 7, 9, 3, 1, 3, 3, 2, 5,
                              0, 1, 0, 0, 0, 1};
        dist.assign(deg_array,deg_array+56);

        dist = normalize_dist(dist, sum(dist));
        net->rand_connect_user(dist);
    } else if (dist == CON) {
        vector<double> dist(param1+1, 0);
        dist[param1] = 1;
        net->rand_connect_user(dist);
    }
}

void processCmdlineParameter(int argc, char* argv[] ) {
    cerr << "Arguments provided: " << argc - 1 << endl;
    cerr << "Argument order: reps re-use_net seasons net_size R0 dist par1 par2 Ih mixture mix_frac p0_ct p0_model RunID [EdgelistFilename]\n";
    cerr << "Mixture Model: 0 = none, 1 = half, 2 = full, 3 = natural, 4 = half_shufle, 5 = full_shuffle\n"
         << "Dist: 0 = poisson, 1 = exponential, 2 = powerlaw, 3 = urban, 4 = constant, 5 = edgelist_file\n";
    if ( argc == 15 || argc == 16 ) {  

        j_max     = (int) atoi(argv[1]);
        reuse_net = (bool) atoi(argv[2]);
        s_max     = (int) atoi(argv[3]); 
        n         = (int) atoi(argv[4]); 
        R_zero    = (double) atof(argv[5]);
        dist      = (DistType) atoi(argv[6]);
        param1    = (double) atof(argv[7]);
        param2    = (double) atof(argv[8]);
        Ih        = (double) atof(argv[9]);
        mixture   = (MixType) atoi(argv[10]);  
        m         = (double) atof(argv[11]);
        patient_zero_ct =  (int) atoi(argv[12]);
        P0_model  = (P0Type) atoi(argv[13]);
        RunID     = argv[14];
        
        if ( dist == EDGELIST and argc == 15) {
            cerr << "An edgelist filename argument must be specified if the distribution type is \"5\" (== edgelist_file)\n"; 
            exit(-1);
        }
        if ( argc == 16 ) { 
            if ( dist != EDGELIST ) {
                cerr << "An edgelist filename argument should not be specified if the distribution type is not \"5\" (== edgelist_file)\n"; 
                exit(-1);
            }
            cerr << "An edgelist filename has been specified.\n"
                 << "Network will be reused for all simulations.\n"
                 << "Network size and parameters for random network generation will be ignored.\n";
            reuse_net = true;
            EdgelistFilename = argv[15];
        }
 
        cout << "# "; for(int i=0; i < argc; i++ )  cout << argv[i] <<  " ";  cout << endl;
        //need to check bounds
        if (mixture > 5 ) { cerr << "Invalid mixture model value\n"; exit(-1); }
        if (dist > 5 )    { cerr << "Invalid distribution value\n"; exit(-1); }
        if (Ih < 0)        { cerr << "Expecting positive immune decay parameter\n"; exit(-1); }
        if (P0_model > 2) { cerr << "Invalid patient-zero model value\n"; exit(-1); }

    } else {
        cerr << "Expecting 14 or 15 parameters.  Set them in config.pm and use sim.pl to as a wrapper.\n";
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

    Network* net = new Network("EpiNet", Network::Undirected);
    MultiSeason_Sim* sim = new MultiSeason_Sim(net,Ih);
    generate_network(net, sim);
    vector<Node*> patients_zero;
    double new_R_zero = R_zero;
    for ( int j = 0; j < j_max; j++) {
        cerr << "Repetition: " << j << endl;
        if (reuse_net == false) {
            delete(net);
            delete(sim);
            net = new Network("EpiNet",Network::Undirected);
            sim = new MultiSeason_Sim(net,Ih);
            generate_network(net,sim);
        }
        double Tc_actual = sim->calc_critical_transmissibility();

        for ( int season = 0; season < s_max; season++) {
            
            if (P0_model == STATIC) {
                patients_zero = sim->rand_infect(patient_zero_ct);
            } else if (P0_model == DYNAMIC) {
                patients_zero = sim->rand_infect_using_immunity(patient_zero_ct);
            } else if (P0_model == CLUSTER) {
                patients_zero = sim->rand_infect_cluster(patient_zero_ct);
            }
            
            sim->run_simulation();
            cout << RunID << " " << j << " " << season << " " << sim->epidemic_size() << " " << patients_zero.size() << " ";
            cout << " " << new_R_zero << " ";
            
            // now calculate what R_zero will be at the start of the next season
            vector<double> average_tk;
            double average_t = 0;
            sim->calculate_average_transmissibility(average_tk, average_t);
            new_R_zero = average_t / Tc_actual;
            cout << endl;
            run_mixture_model(mixture, net);
            //sim.summary();
        }
        //sim->write_event_history("/tmp/eventHistory");
        sim->reset();
        new_R_zero = R_zero;
    }
}

void generate_network(Network* net, MultiSeason_Sim* sim) {
    if ( dist == EDGELIST ) {
        net->read_edgelist(EdgelistFilename, ',');
    } else {
        net->populate(n);
        connect_network(net); // connect network using the parameters above
    }
    sim->calc_naive_transmissibility(R_zero); //calculates correct T for this network topology
}

void run_mixture_model(MixType mixture, Network* net) { 
    if (mixture == NO) {
       return;
    } else if (mixture == FULL) {
        net->clear_edges();
        connect_network(net);
    } else if (mixture == FULL_SHUFFLE) {
        vector<int> states = net->get_states();
        shuffle(states, net->get_rng());
        vector<Node*> nodes = net->get_nodes();
        for (int i = 0; i <(signed)  nodes.size(); i++) {
            nodes[i]->set_state(states[i]);
        }
    } else if (mixture == HALF) {
        net->disconnect_edges();
        net->rand_connect_stubs(net->get_edges());
    } else if (mixture == HALF_SHUFFLE) { // topology stays same, states are shuffled within degree classes
        vector< vector<int> > states_by_degree = net->get_states_by_degree();
        for ( int d = 0; d<(signed) states_by_degree.size(); d++ ) shuffle(states_by_degree[d], net->get_rng());

        vector<int> deg_dist = net->get_deg_dist();
        vector<int> counter(deg_dist.size(), 0);
        vector<Node*> nodes = net->get_nodes();
        for ( int i = 0; i < (signed) nodes.size(); i++ ) {
            Node* node = nodes[i];
            int deg    = node->deg();
            node->set_state( states_by_degree[deg][counter[deg]] );
            counter[deg]++;
        }

    } else if (mixture == NATURAL) {
        // m is the fraction of nodes to be 'mixed',
        // so k is the total number of nodes to be mixed 
        int k = (int) (m * net->size());
        // create empty vectors that will store the nodes to be
        // mixed and the new degree series that will be assigned
        // to them
        vector<int> sample(k);
        vector<int> new_degrees(k);
        // generate the new degrees (their sum will be even,
        // ensuring that no stubs will be left over)
        net->gen_deg_series(new_degrees);

        MTRand* mtrand = net->get_rng();
        vector<double> dist = net->get_gen_deg_dist();
        rand_nchoosek( (int) net->size() - 1, sample, mtrand);

        vector<int> deltas(k);
        int delta_sum = 0;
        int node_id, deg, new_deg;
        Node* node;
        for (int i = 0; i < (signed) sample.size(); i++) {
            node_id = sample[i];
            node = net->get_node(node_id);
            deg = node->deg();
            new_deg = rand_nonuniform_int(dist, mtrand);
            deltas[i] = new_deg - deg;
            delta_sum += deltas[i];
        }
        while ( delta_sum % 2 != 0 ) {
            int i = mtrand->randInt(k - 1);
            node_id = sample[i];
            node = net->get_node(node_id);
            deg = node->deg();
            new_deg = rand_nonuniform_int(dist, mtrand);
            delta_sum -= deltas[i];
            deltas[i] = new_deg - deg;
            delta_sum += deltas[i];
        }

        vector<Edge*> stubs;

        for (int i = 0; i < (signed) sample.size(); i++) {
            node_id = sample[i];
            node = net->get_node(node_id);
            deg = node->deg();
            new_deg = deltas[i] + deg;
            while (new_deg < node->deg()) {
                Edge* edge = node->get_rand_edge();
                if (edge->get_end() != NULL) {
                    Edge* comp;// = edge->get_complement();
                    vector<Edge*> edges_out = edge->get_end()->get_edges_out(); // get the edges leaving the endpoint;
                    for (int i = 0; i < (signed) edges_out.size(); i++) {
                        if (edges_out[i] == edge) continue;
                        if (edges_out[i]->get_end() == edge->get_start()) {
                            comp = edges_out[i];
                            comp->break_end();
                            stubs.push_back( comp );
                            break;
                        }
                    }
                }
                if (edge->get_end() == NULL) { delete_element(stubs, edge); }
                edge->delete_edge();

            }
            while (new_deg > node->deg()) {
                Edge* stub = node->add_stub_out();
                stubs.push_back( stub );
            }
        }

        net->rand_connect_stubs(stubs);
    }
} 

