#ifndef MULTISEASON_SIM_H
#define MULTISEASON_SIM_H

#include "Percolation_Sim.h"
//#include <assert.h>
//#include <queue>

class MultiSeason_Sim: public Percolation_Sim
{
    //double D;                 // immune decay parameter (T_node = T * (1 - exp(-D*G))
    double Ih;                 // immunity halflife parameter (T_node = T * (1 - 2^(-G/Ih))

    
    //int current_epi_size;       // for r_zero_test; can be deleted later
    //double mean_sec_infections; // for r_zero_test
    //int infections_to_watch;    // for r_zero_test
    
    //Ruben's version: inline double immune_decay(double D, int G) {if (G == -1) return 1; return exp(-D/G);}
    //inline double immune_decay(double D, int G) {if (G == -1) return 1; return 1 - exp(-D*G);} 
    inline double immune_decay(double Ih, int G) {if (G == -1 or Ih == 0) return 1; return 1.0 - pow(2.0, -G/Ih);} 
    inline int time_since_infection(Node* node) {return node->get_state() - 1;}

    public:
        //typedef enum { S=0, I=1, R=2 } stateType; //Whatever is equal to zero is the default state
        //int epi_size;
        //vector<Node*> recovered; // we don't need this for the simulation per se, but we do need it
                                 // in order to quickly reset the network to completely susceptible

        MultiSeason_Sim():Percolation_Sim() {}
        MultiSeason_Sim(Network* net, double Ih):Percolation_Sim(net) { this->Ih = Ih;}
        
        void set_immune_decay_par( double Ih ) { this->Ih = Ih; }
        double calc_naive_transmissibility(double R_zero) { this->T = R_zero * calc_critical_transmissibility(); return T; }
        double get_naive_transmissibility() { return T; }

        double get_transmissibility(Node* end) { // uses the exponential decay model of immunity loss
            int G = time_since_infection(end);
            return (double) (T * immune_decay(Ih,G)); // Ih is the immunity halflife parameter, G is seasons since last infection 
        }

        vector<Node*> rand_infect_using_immunity (int n) { // attempt to infect n nodes, succeeding with probability 2^(-G/Ih)
            vector<Node*> nodes = net->get_nodes();
            vector<Node*> p_zeros;
            double p_naive = (double) n / nodes.size(); // expected fraction infected in naive network
            for (int i = 0; i<(signed) nodes.size(); i++) {
                Node* node = nodes[i];
                double p = p_naive * immune_decay(Ih, time_since_infection(node));
                if ( mtrand->rand() < p ){
                    node->set_state(1); 
                    infected.push_back( node );
                    p_zeros.push_back( node );
                }
            }
            return p_zeros;
        }
/*
        vector<Node*> rand_set_nodes_to_state(int n, stateType state) {
            vector<Node*> patients_zero = Percolation_Sim::rand_set_nodes_to_state(n, state);
            vector<Node*>::iterator itr;
            itr = infected.end();
            infected.insert(itr, patients_zero.begin(), patients_zero.end());
            return patients_zero;
        }
*/
        // this is to test the effects of clumped vs. scattered patients-zero
        vector<Node*> rand_infect_cluster(int n) {
            int attempts = 0;
            int max_attempts = 100;
            while (attempts < max_attempts) {
                Node* p_zero = net->get_rand_node();
                p_zero->set_state(1);
                vector<Node*> patients_zero;
                patients_zero.push_back( p_zero );
                vector<Node*> tmp_infected(1, p_zero);
                vector<Node*> new_infected;
                while ((signed) patients_zero.size() < n && tmp_infected.size() > 0) {
                    for (int i = 0; i < (signed) tmp_infected.size(); i++) {
                        Node* inode = tmp_infected[i];
                        vector<Node*> neighbors = inode->get_neighbors();
                        for (int j = 0; j < (signed) neighbors.size(); j++) {
                            Node* test = neighbors[j];
                            test->set_state( 1 );
                            new_infected.push_back( test );
                            patients_zero.push_back( test );
                            if ((signed) patients_zero.size() >= n ) { goto DONE; }
                        }
                    }
                    tmp_infected = new_infected;
                    new_infected.clear();
                }
                DONE:
                vector<Node*>::iterator itr;
                itr = infected.end();
                infected.insert(itr, patients_zero.begin(), patients_zero.end());
                if ( (signed) patients_zero.size() == n ) { 
                    return patients_zero;
                } else {
                    attempts++;
                }
            }
            cerr << "NetSim::rand_infect_cluster() failed to infect a cluster of size " << n << " after " << max_attempts << " attempts." << endl;
            exit(1);
        }

        void step_simulation () {
            time++;
            assert(infected.size() > 0);
            vector<Node*> new_infected;
            for (int i = 0; i < (signed) infected.size(); i++) {
                //current_epi_size++; // for r_zero_test
                //int sec_infection_ct = 0; // for r_zero_test
                Node* inode = infected[i];
                vector<Node*> neighbors = inode->get_neighbors();
                for (int j = 0; j < (signed) neighbors.size(); j++) {
                    Node* test = neighbors[j];
                    if ( mtrand->rand() < get_transmissibility(test) ) {
                        //sec_infection_ct++; // for r_zero_test
                        test->set_state( 1 );
                        new_infected.push_back( test );
                    }
                }
                
                //if (current_epi_size > 1 && current_epi_size <= infections_to_watch + 1) { // for r_zero_test; ignore P_zero
                    //mean_sec_infections += (double) sec_infection_ct / infections_to_watch; // for r_zero_test
                //} // for r_zero_test
                recovered.push_back(inode);
            }
            infected = new_infected;
        }


        void calculate_average_transmissibility(vector<double> &average_tk, double &average_t) {
            average_t = 0;
            vector<int> deg_dist = net->get_deg_dist();
            vector<double> tk( deg_dist.size() );    // total transmissibility by degree -- used to calc mean
		    vector<Node*> nodes = net->get_nodes();

            for (int i = 0; i < (signed) nodes.size(); i++) {
                int deg = nodes[i]->deg();
                tk[deg] += get_transmissibility(nodes[i]);
            }
            
            average_tk.resize( tk.size() ); // average transmissibility of degree k nodes
            int deg_sum = sum(net->get_deg_series()); // sum of all degrees in the network, i.e. the total number of stubs
            for (int deg = 0; deg < (signed) tk.size(); deg++) {
                average_tk[deg] = tk[deg]/ ((double) deg_dist[deg]);
                average_t += tk[deg] * deg / ((double) deg_sum);
            }
        }
		
        void run_simulation() {
            recovered.clear();
            //epi_size = 0;
            //current_epi_size = 0;      // for r_zero_test
            //mean_sec_infections = 0;   // for r_zero_test
            //infections_to_watch = 25;  // for r_zero_test
            vector<Node*> nodes = net->get_nodes();
            while (infected.size() > 0) { //let the epidemic spread until it runs out
                step_simulation();
                ///////////////////////////////
            }
        
            for ( int i = 0; i < (signed)  nodes.size(); i++ ) {
                int curr_state = nodes[i]->get_state();
                //if ( curr_state == 1 ) epi_size++;

                if ( curr_state >= 1 ) nodes[i]->set_state(curr_state + 1);
            }
            infected.clear();
//            if (current_epi_size > 1) { // for r_zero_test
//                if (current_epi_size < infections_to_watch + 1) { // for r_zero_test; ignore P_zero
//                    mean_sec_infections *= (double) infections_to_watch / current_epi_size; // for r_zero_test
//                } // for r_zero_test
//                cout << epi_size << " " << mean_sec_infections << endl; // for r_zero_test
//            } else { // for r_zero_test
//                cout << epi_size << " " << "NAN" << endl; // for r_zero_test
//            } // for r_zero_test
        }

        
        // For each degree, count the number of nodes that have each state
        vector< vector<int> > tabulate_states() {
            vector< vector<int> > states_by_degree = net->get_states_by_degree();
            vector< vector<int> > tabulated(states_by_degree.size());
            for (int d = 0; d<(signed) states_by_degree.size(); d++) tabulated[d] = tabulate_vector(states_by_degree[d]);
            return tabulated;
        }

        void reset() {
            reset_time();
            
            set_these_nodes_to_state(net->get_nodes(), 0);
            
            infected.clear();
            recovered.clear();
        }
        
};

/*

Network net("test", false );
net.populate(10000);

NetSim sim(net);
*/
#endif
