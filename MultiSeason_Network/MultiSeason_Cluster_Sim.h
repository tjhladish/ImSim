#ifndef MULTISEASONCLUSTER_SIM_H
#define MULTISEASONCLUSTER_SIM_H

#include "Percolation_Sim.h"

class MultiSeason_Cluster_Sim: public Percolation_Sim
{
    double CI;                 // cross immunity between clusters
    int Cluster;

    inline double cluster_susceptibility(double CI, int last_cluster ) { 
        if (last_cluster == 0) {
            return 1.0;
        } else {
            return pow(CI, Cluster - last_cluster);
        }
    } 

    inline int last_infecting_cluster(Node* node) {return node->get_state();}

    public:
        MultiSeason_Cluster_Sim():Percolation_Sim() {}
        MultiSeason_Cluster_Sim(Network* net, double CI):Percolation_Sim(net) { this->CI = CI; }
        
        void set_crossimmunity(int ci) {CI = ci;}
        int get_crossimmunity() {return CI;}
        
        void set_cluster(int cluster) {Cluster = cluster;}
        int get_cluster() {return Cluster;}
        
        double calc_naive_transmissibility(double R_zero) { this->T = R_zero * calc_critical_transmissibility(); return T; }
        double get_naive_transmissibility() { return T; }

        double get_transmissibility(Node* end) { // uses the exponential decay model of immunity loss
            int last_cluster = last_infecting_cluster(end);
            return (double) (T * cluster_susceptibility(CI, last_cluster));
        }

        void step_simulation () {
            time++;
            vector<Node*> new_infected;
            for (int i = 0; i < (signed) infected.size(); i++) {
                Node* inode = infected[i];
                vector<Node*> neighbors = inode->get_neighbors();
                for (int j = 0; j < (signed) neighbors.size(); j++) {
                    Node* test = neighbors[j];
                    if ( last_infecting_cluster(test) != Cluster and mtrand->rand() < get_transmissibility(test) ) {
                        test->set_state( Cluster );
                        new_infected.push_back( test );
                    }
                }
                
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
            vector<Node*> nodes = net->get_nodes();
            while (infected.size() > 0) { //let the epidemic spread until it runs out
                step_simulation();
            }
        
            infected.clear();
        }

        void reset() {
            reset_time();
            
            set_these_nodes_to_state(net->get_nodes(), 0);
            
            infected.clear();
            recovered.clear();
        }
        
};

#endif
