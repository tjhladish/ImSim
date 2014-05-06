#ifndef MSMS_CLUSTER_SIM_H
#define MSMS_CLUSTER_SIM_H
#include "Simulator.h"

//#include <assert.h>
//#include <queue>

typedef enum { H1N1,
               H3N2,
               B,
               NumStrains } StrainType; // must be last


class MSMSC_Sim: public Simulator {
    double CHI;                 // cross immunity between clusters
    vector<int> clusters;

    vector<Node*> infected;   // list of infected nodes; used so we don't have to check every node's state when transmitting
    vector<double> T;                 // naive transmission probability
    vector< vector<int> > node_history; // immune state of all nodes (1st dim) for all strains (2nd dim)
                                       // value is 0 for never infected, cluster number for previously infected
    vector< vector<bool> > recent_infections; // node id is 1st dim, whether infected this season is 2nd dim
    //vector<Node*> recovered; // we don't need this for the simulation per se, but we do need it
                             // in order to quickly reset the network to completely susceptible
    

    inline int last_infecting_cluster(Node* node, StrainType strain) { return node_history[node->get_id()][strain]; }

    inline double susceptibility(Node* node, StrainType challenging_strain, int challenging_cluster ) { 
        int last_cluster = last_infecting_cluster( node, challenging_strain );
        if (last_cluster == 0) {
            return 1.0;
        } else {
            return 1.0 - pow(CHI, challenging_cluster - last_cluster);
        }
    } 

   
    // this should not be called without first checking that the node has, in fact, been previously infected
    inline int clusters_since_infection(Node* node, StrainType strain) {return clusters[strain] - node_history[node->get_id()][strain];}

    public:
        MSMSC_Sim():Simulator() {}
        MSMSC_Sim(Network* net, double CHI):Simulator(net) {
            this->CHI = CHI;
            set_these_nodes_to_state(net->get_nodes(), -1);
            node_history.resize( net->size(), vector<int>(NumStrains, 0) );
            T.resize( NumStrains, 0 );
            clusters.resize( NumStrains, 1 );
            recent_infections.resize( net->size(), vector<bool>(NumStrains, false) );
        }
        
        void set_network( Network* net ) { 
            Simulator::set_network( net );
            set_these_nodes_to_state(net->get_nodes(), -1);
            node_history.clear();
            node_history.resize( net->size(), vector<int>(NumStrains, 0) ); // -1 is never infected
            T.resize( NumStrains, 0 ); // calc_naive_transmissibility() must be called to set T values
            clusters.resize( NumStrains, 1 ); // all clusters start at 1
            recent_infections.resize( net->size(), vector<bool>(NumStrains, false) );
        };

        void set_cluster( StrainType strain, int cluster ) { clusters[strain] = cluster; }
        void set_cluster_crossimmunity( double CHI ) { this->CHI = CHI; }

        vector<double> calc_naive_transmissibility(vector<double> R_zeros) {
            for (int s = 0; s < (int) NumStrains; s++) {
                this->T[s] = R_zeros[s] * calc_critical_transmissibility();
            } 
            return T; 
        }

        vector<double> get_naive_transmissibility() { return T; }

        double get_transmissibility(Node* node, StrainType strain) {
            const int id = node->get_id();
            const int cluster = clusters[strain];
            double t; // Transmission probability given immunity
            // are we challenging with influenza A when we've been infected this year with influenza A?
            if ( (strain == H1N1 or strain == H3N2) and (recent_infections[id][H1N1] == true or recent_infections[id][H3N2] == true) ) {
                t = 0.0; // Short-term heterosubtypic immunity
            } else {
                t = T[strain] * susceptibility(node, strain, cluster); 
            }
            return t;
        }

        void infect(Node* node, StrainType strain) {
            node->set_state(strain); // most recently infecting strain this year
            const int id = node->get_id();
            node_history[id][strain] = clusters[strain]; // update node's infection history
            recent_infections[id][strain] = true;
        }

        vector<Node*> rand_infect (vector<int> initial_infections_by_strain) {
            vector<StrainType> infection_types; 
            int total_p_zero = 0;
            for ( unsigned int s = 0; s < initial_infections_by_strain.size(); ++s ) {
                total_p_zero += initial_infections_by_strain[s];
                // build up a vector of strains, one for each infection
                infection_types.resize(total_p_zero, (StrainType) s);
            }
            assert(total_p_zero > 0 and total_p_zero <= net->size() );
            shuffle(infection_types, mtrand); // probably not necessary to shuffle
            vector<Node*> patients_zero = rand_choose_nodes(total_p_zero);
            for ( unsigned int i = 0; i < patients_zero.size(); ++i) {
                const StrainType strain = infection_types[i];
                infect( patients_zero[i], strain );
                infected.push_back(patients_zero[i]); // keep track to facilitate future transmission
            };
            return patients_zero;
        }


        void step_simulation () {
            time++;
            //assert(infected.size() > 0);
            vector<Node*> new_infected;
            for (int i = 0; i < (signed) infected.size(); i++) {
                Node* inode = infected[i];
                // nodes transmit with whatever the most recently infected strain is
                StrainType strain = (StrainType) inode->get_state();
                vector<Node*> neighbors = inode->get_neighbors();
                for (int j = 0; j < (signed) neighbors.size(); j++) {
                    Node* test = neighbors[j];
                    if ( mtrand->rand() < get_transmissibility(test, strain) ) {
                        infect( test, strain );
                        new_infected.push_back( test );
                    }
                }
                //recovered.push_back(inode);
            }
            infected = new_infected;
        }


        vector<int> final_size_by_strain() {
            vector<int> tally(NumStrains, 0);
            for ( unsigned int i = 0; i < recent_infections.size(); ++i ) {
                for ( int s = 0; s < NumStrains; ++s ) {
                    if ( recent_infections[i][s] == true ) ++tally[s];
                }
            } 
            return tally;
        }

        // can be greater than pop size, since there is no heterotypic immunity
        int total_infections() { return sum( final_size_by_strain() ); } 
        int epidemic_size() { return total_infections(); } 

        void calculate_average_transmissibility(vector<double> &average_tk, double &average_t, StrainType strain) {
            average_t = 0;
            vector<int> deg_dist = net->get_deg_dist();
            vector<double> tk( deg_dist.size() );    // total transmissibility by degree -- used to calc mean
            vector<Node*> nodes = net->get_nodes();

            for (int i = 0; i < (signed) nodes.size(); i++) {
                int deg = nodes[i]->deg();
                tk[deg] += get_transmissibility(nodes[i], strain);
            }
            
            average_tk.resize( tk.size() ); // average transmissibility of degree k nodes
            int deg_sum = sum(net->get_deg_series()); // sum of all degrees in the network, i.e. the total number of stubs
            for (int deg = 0; deg < (signed) tk.size(); deg++) {
                average_tk[deg] = tk[deg]/ ((double) deg_dist[deg]);
                average_t += tk[deg] * deg / ((double) deg_sum);
            }
        }

        void run_simulation() {
            vector<Node*> nodes = net->get_nodes();
            recent_infections.clear();
            recent_infections.resize( nodes.size(), vector<bool>(3,false) );

            //let the epidemic spread until it runs out
            while (infected.size() > 0) { 
                step_simulation();
            }
        
            infected.clear();
        }

        /* 
        // For each degree, count the number of nodes that have each state
        vector< vector<int> > tabulate_states() {
            vector< vector<int> > states_by_degree = net->get_states_by_degree();
            vector< vector<int> > tabulated(states_by_degree.size());
            for (unsigned int d = 0; d < states_by_degree.size(); ++d) {
                tabulated[d] = tabulate_vector(states_by_degree[d]);
            }
            return tabulated;
        }
        */

        void reset() {
            reset_time();
            
            set_these_nodes_to_state(net->get_nodes(), -1);
            node_history.clear();
            node_history.resize( net->size(), vector<int>(NumStrains, 0) );

            infected.clear();
            //recovered.clear();
        }
        
};

#endif
