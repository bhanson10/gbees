#ifndef GBEES_H
#define GBEES_H

#include "BST.h"
/*==============================================================================
CLASS DEFINITIONS
==============================================================================*/
class GBEES{       // GBEES Class
  public:
    BST P;
    TreeNode* dead; 
    int a_count = 0; 
    int tot_count = 1; 
    uint64_t max_key = 0; 

    void Initialize_D(Grid G, Traj lyap);
    void Initialize_vuw(Grid G, Traj lyap, TreeNode* r);
    void Initialize_ik_nodes(Grid G, TreeNode* r); 
    void Modify_pointset(Grid G, Traj lyap);
    TreeNode* create_neighbors(Grid G, TreeNode* r, double& prob_sum);
    TreeNode* delete_neighbors(Grid G, TreeNode* r);
    void normalize_prob(TreeNode* r, double prob_sum);
    void RHS_P(Grid G, Traj lyap, int rk);
    void original_prob(TreeNode* r); 
    void initial_f(Grid G, TreeNode* r); 
    void total_f(Grid G, TreeNode* r); 
    void calculating_K(Grid G, TreeNode* r, int rk);
    void update_prob(Grid G, TreeNode* r, int rk); 
    void get_size(Grid G, TreeNode* r, int& a_count, int& tot_count, uint64_t& max_key); 
    void Record_Data(std::string file_name, Grid G, TreeNode* r);
    void measurement_update(Measurement m, Grid G, TreeNode* r, double& C);
};
/*==============================================================================
Non-member Function DEFINITIONS
==============================================================================*/
uint64_t CantorPair(int state[], int n){
    uint64_t key; 
    if(n>2){
        int last = state[n-1]; int new_state[n-1];
        for (int i = 0; i < n - 1; i++){
            new_state[i] = state[i];
        }
        int x = CantorPair(new_state, n-1); int y = last;
        key = static_cast<uint64_t>((0.5)*(x+y)*(x+y+1)+y);
    }else{
        int x = state[0]; int y = state[1];
        key = static_cast<uint64_t>((0.5)*(x+y)*(x+y+1)+y);
    }
    return key; 
};

uint64_t RosenbergPair(int state[], int d, int m){
    uint64_t key; 
    if(d == 1){
        return key = static_cast<uint64_t>(state[0]);
    }else{
        int new_state[d-1];
        for (int i = 0; i < d-1; i++){
            new_state[i] = state[i]; 
        }
        int new_m = *std::max_element(new_state, new_state + d-1);
        return key = RosenbergPair(new_state, d-1, new_m) + pow(m,d) + (m - state[d-1])*(pow(m+1, d-1) - pow(m,d-1)); 
    }
};

uint64_t state_conversion(std::array<int,DIM> state, Grid G){
    int shift_state[DIM]; int m; uint64_t key; 
    for (int i = 0; i < DIM; i++){
        if(state[i]<0){
            shift_state[i] = -2*state[i]-1;
        }else{
            shift_state[i] = 2*state[i];
        }
    }
    if(DIM == 1){
        key = static_cast<uint64_t>(shift_state[0]);
    }else{
        if(G.pair==1){
            key = CantorPair(shift_state, DIM);
        }
        else if(G.pair==2){
            m = *std::max_element(shift_state, shift_state + DIM);
            key = RosenbergPair(shift_state, DIM, m);
        }
    }
    return key; 
};

double MC(double th){
    double phi;
    phi = std::max(0.0,std::min({(1+th)/2,2.0,2*th}));
    return phi;
};
/*==============================================================================
Member Function DEFINITIONS
==============================================================================*/
void GBEES::Initialize_D(Grid G, Traj lyap){
    Cell blank_c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .f = {0}, .status1 = -1, .status2 = -1, .K = {0}, .x_n = 0};
    TreeNode* dead_node = new TreeNode(-1, blank_c); dead_node->cell.i_nodes = {dead_node, dead_node, dead_node, dead_node, dead_node, dead_node}; dead_node->cell.k_nodes = {dead_node, dead_node, dead_node, dead_node, dead_node, dead_node};
    dead = dead_node; 

    std::array<int,DIM> current_state; uint64_t key; 
    for (int i = round((-3*G.std[0])/G.del[0]); i <= round((3*G.std[0])/G.del[0]); i++){
        for (int j = round((-3*G.std[1])/G.del[1]); j <= round((3*G.std[1])/G.del[1]); j++){
            for (int k = round((-3*G.std[2])/G.del[2]); k <= round((3*G.std[2])/G.del[2]); k++){
                for (int l = round((-3*G.std[3])/G.del[3]); l <= round((3*G.std[3])/G.del[3]); l++){
                    for (int m = round((-3*G.std[4])/G.del[4]); m <= round((3*G.std[4])/G.del[4]); m++){
                        for (int n = round((-3*G.std[5])/G.del[5]); n <= round((3*G.std[5])/G.del[5]); n++){
                            current_state = {i,j,k,l,m,n}; key = state_conversion(current_state, G);
                            /*
                            std::cout << i << " " << j << " " << k << " " << l << std::endl; 
                            std::cout << key << std::endl;
                            std::cin.get(); 
                            */

                            double x = 0; 
                            for(int q = 0; q < DIM; q++){
                                x += pow((current_state[q]*G.del[q]),2)/pow(G.std[q],2); 
                            }
                            
                            Cell c = {.prob = exp(-4*x/2), .v = {0}, .u = {0}, .w = {0}, .f = {0}, .state = current_state, .status1 = 0, .status2 = 0, .K = {0}, .x_n = 0};
                            TreeNode* new_node = new TreeNode(key, c); P.root = P.insertRecursive(P.root, new_node);
                        }
                    }
                }
            }
        }
    }

    Initialize_vuw(G,lyap,P.root);
    Initialize_ik_nodes(G,P.root); 
};

void GBEES::Initialize_vuw(Grid G, Traj lyap, TreeNode* r){
    if(r==NULL){ //Iterator over entire BST
        return;
    }
    Initialize_vuw(G,lyap,r->left);
    Initialize_vuw(G,lyap,r->right);  
    
    if(r->cell.status1==0){
        std::array<double,DIM> x; 
        for(int i = 0; i < DIM; i++){
            x[i] = G.del[i]*r->cell.state[i]+G.epoch[i];
        }
        double r1 = pow(pow(x[0]+lyap.mu,2)+pow(x[1],2)+pow(x[2],2),1.5); 
        double r2 = pow(pow(x[0]-1+lyap.mu,2)+pow(x[1],2)+pow(x[2],2),1.5); 

        double v1 = x[3];
        double v2 = x[4];
        double v3 = x[5];
        double v4 = 2*x[4]+x[0]-((1-lyap.mu)*(x[0]+lyap.mu)/r1)-((x[0]-1+lyap.mu)*(lyap.mu)/r2);
        double v5 = -2*x[3]+x[1]-((1-lyap.mu)*x[1]/r1)-((lyap.mu)*x[1]/r2);
        double v6 = -((1-lyap.mu)*x[2]/r1)-((lyap.mu)*x[2]/r2);
        r->cell.v = {v1,v2,v3,v4,v5,v6};
        r->cell.u = {std::min(v1,0.0),std::min(v2,0.0),std::min(v3,0.0),std::min(v4,0.0),std::min(v5,0.0),std::min(v6,0.0)};
        r->cell.w = {std::max(v1,0.0),std::max(v2,0.0),std::max(v3,0.0),std::max(v4,0.0),std::max(v5,0.0),std::max(v6,0.0)}; 
        r->cell.status1 = 1;
    }
};

void GBEES::Initialize_ik_nodes(Grid G,TreeNode* r){
    if (r==NULL){
        return; 
    }
    Initialize_ik_nodes(G,r->left);
    Initialize_ik_nodes(G,r->right);

    if(r->cell.status2 == 0){
        std::array<int, DIM> l_state = r->cell.state;
        for(int q = 0; q < DIM; q++){
            //Initializing i, k nodes
            std::array<int,DIM> i_state = l_state; i_state[q] = i_state[q]-1; uint64_t i_key = state_conversion(i_state,G);
            std::array<int,DIM> k_state = l_state; k_state[q] = k_state[q]+1; uint64_t k_key = state_conversion(k_state,G);
            TreeNode* i_node = P.recursiveSearch(P.root, i_key); TreeNode* k_node = P.recursiveSearch(P.root, k_key); 
            
            if(i_node == NULL){
                i_node = dead; r->cell.i_nodes[q] = i_node; 
            }else{
                r->cell.i_nodes[q] = i_node; i_node->cell.k_nodes[q] = r; 
            }

            if(k_node == NULL){
                k_node = dead; r->cell.k_nodes[q] = k_node; 
            }else{
                r->cell.k_nodes[q] = k_node; k_node->cell.i_nodes[q] = r; 
            }      
        }
        r->cell.status2 = 1; 
    }
};

void GBEES::Modify_pointset(Grid G, Traj lyap){
    double prob_sum = 0; 
    P.root = create_neighbors(G, P.root, prob_sum);
    Initialize_ik_nodes(G,P.root); 
    Initialize_vuw(G, lyap, P.root); 
    normalize_prob(P.root, prob_sum);
    if(!P.isBalanced(P.root)){
        P.root = P.balance(P.root);
    }
};

TreeNode* GBEES::create_neighbors(Grid G, TreeNode* r, double& prob_sum){
    if (r == NULL){
        return r;
    }
    create_neighbors(G, r->left, prob_sum);
    create_neighbors(G, r->right, prob_sum);

    if (r->cell.prob >= G.thresh){
        prob_sum += r->cell.prob; 
        bool outsider = false; 
        for(int i = 0; i < DIM; i++){
            if ((r->cell.i_nodes[i]==dead)||(r->cell.k_nodes[i]==dead)){
                outsider = true; 
            }
        };

        if(outsider == true){
            std::array<double,DIM> current_v = r->cell.v; std::array<int,DIM> current_state = r->cell.state;    

            std::array<std::array<int,2>,DIM> n_states;  

            for(int c = 0; c < DIM; c++){
                if(current_v[c] > 0){
                    n_states[c][0] = 0; n_states[c][1] = 1;
                }else if(current_v[c]==0){
                    n_states[c][0] = -1; n_states[c][1] = 1;
                }else{
                    n_states[c][0] = -1; n_states[c][1] = 0;
                }
            }

            for (int i = n_states[0][0]; i <= n_states[0][1]; i++){
                for (int j = n_states[1][0]; j <= n_states[1][1]; j++){
                    for (int k = n_states[2][0]; k <= n_states[2][1]; k++){
                        for (int l = n_states[3][0]; l <= n_states[3][1]; l++){
                            for (int m = n_states[4][0]; m <= n_states[4][1]; m++){
                                for (int n = n_states[5][0]; n <= n_states[5][1]; n++){
                                    std::array<int,DIM> state_update = {i, j, k, l, m, n}; 
                                    
                                    TreeNode* neighbor_node = r;
                                    for (int p = 0; p < DIM; p++){
                                        if(state_update[p] == 1){
                                            neighbor_node = neighbor_node->cell.k_nodes[p];
                                        }else if (state_update[p] == -1){
                                            neighbor_node = neighbor_node->cell.i_nodes[p];
                                        }

                                        if(neighbor_node == dead){
                                            std::array<int,DIM> new_state = current_state; 
                                            for(int q = 0; q < DIM; q++){
                                                new_state[q] += state_update[q]; 
                                            }
                                            uint64_t new_key = state_conversion(new_state,G);
                                            TreeNode* test_node = P.recursiveSearch(P.root, new_key); 
                                            if(test_node == NULL){
                                                Cell c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .f = {0}, .state = new_state, .status1 = 0, .status2 = 0, .K = {0}, .x_n = 0};
                                                TreeNode* new_node = new TreeNode(new_key, c);
                                                P.root = P.insertRecursive(P.root, new_node); 
                                            } 
                                            break;      
                                        }
                                    }  
                                }
                            }            
                        }
                    }
                }
            }     
        }  
    }
    return P.root; 
};

TreeNode* GBEES::delete_neighbors(Grid G, TreeNode* r){
    if (r==NULL){
        return r;
    }
    delete_neighbors(G, r->left);
    delete_neighbors(G, r->right);

    r->cell.status2 = 0; 
    if ((r->cell.prob < G.thresh)&&(r->cell.status1 == 1)){
        bool neighbors = true; std::array<double,DIM> current_v = r->cell.v; std::array<int,DIM> current_state = r->cell.state;    

        std::array<int,DIM> num_n; std::array<std::array<int,DIM>,DIM> which_n;  

        for(int c = 0; c < DIM; c++){
            if(current_v[c] < 0){
                num_n[c] = 2; 
                which_n[c][0] = current_state[c]; which_n[c][1] = current_state[c]+1;
            }else if(current_v[c]==0){
                num_n[c] = 3; 
                which_n[c][0] = current_state[c]-1; which_n[c][1] = current_state[c]; which_n[c][2] = current_state[c]+1;
            }else{
                num_n[c] = 2;
                which_n[c][0] = current_state[c]-1; which_n[c][1] = current_state[c];
            }
        }

        for (int i = 0; i < num_n[0]; i++){
            for (int j = 0; j < num_n[1]; j++){
                for (int k = 0; k < num_n[2]; k++){
                    for (int l = 0; l < num_n[3]; l++){
                        for (int m = 0; m < num_n[4]; m++){
                            for (int n = 0; n < num_n[5]; n++){
                                std::array<int,DIM> new_state = {which_n[0][i], which_n[1][j], which_n[2][k], which_n[3][l], which_n[4][m], which_n[5][n]}; uint64_t new_key = state_conversion(new_state,G);
                                TreeNode* new_node = P.recursiveSearch(P.root, new_key); 
                                if(new_node != NULL){
                                    if(new_node->cell.prob >= G.thresh){
                                        neighbors = false;
                                        break; 
                                    }
                                }    
                            }
                        }        
                    }      
                }
            }
        }   
        if(neighbors){
            P.root = P.deleteNode(P.root, r->key);
        }
    }
    
    return P.root; 
};

void GBEES::normalize_prob(TreeNode* r, double prob_sum){
    if (r==NULL){
        return;
    }
    normalize_prob(r->left, prob_sum);
    normalize_prob(r->right, prob_sum);

    double current_p = r->cell.prob;

    r->cell.prob = current_p/prob_sum; 
};

void GBEES::RHS_P(Grid G, Traj lyap, int rk){
    if (rk == 1){
        original_prob(P.root); 
    }
    initial_f(G, P.root); 
    total_f(G, P.root); 
    calculating_K(G, P.root, rk);
};

void GBEES::original_prob(TreeNode* r){
    if (r==NULL){
        return; 
    }
    original_prob(r->left);
    original_prob(r->right);

    if (r->key != -1){
        r->cell.x_n = r->cell.prob;
    }
};

void GBEES::initial_f(Grid G, TreeNode* r){
    if (r==NULL){
        return; 
    }
    initial_f(G,r->left);
    initial_f(G,r->right);

    int l_key = r->key;
    if(l_key != -1){
        r->cell.f = {0.0}; 
        for(int q = 0; q < DIM; q++){
            r->cell.f[q] = r->cell.w[q] * r->cell.prob + r->cell.u[q] * r->cell.k_nodes[q]->cell.prob;
        }
    }
};

void GBEES::total_f(Grid G, TreeNode* r){
    if (r==NULL){
        return;
    }
    total_f(G,r->left);
    total_f(G,r->right);

    int l_key = r->key; 
    if(l_key != -1){ 
        for(int q = 0; q < DIM; q++){
            TreeNode* i_node = r->cell.i_nodes[q]; 
            if ((r->cell.prob >= G.thresh)||(i_node->cell.prob >= G.thresh)){ 
                double F = G.dt*(r->cell.prob-i_node->cell.prob)/(2*G.del[q]); 
                for(int e = 0; e < DIM; e++){
                    if (e!=q){
                        TreeNode* j_node = r->cell.i_nodes[e]; 
                        TreeNode* p_node = i_node->cell.i_nodes[e];
                        
                        r->cell.f[e] -= r->cell.w[e] * i_node->cell.w[q] * F; 
                        i_node->cell.f[e] -= i_node->cell.w[e] * i_node->cell.u[q] * F; 
                        j_node->cell.f[e] -= j_node->cell.u[e] * i_node->cell.w[q] * F; 
                        p_node->cell.f[e] -= p_node->cell.u[e] * i_node->cell.u[q] * F; 
                    }                       
                }
                
                double th,t; 
                if (i_node->cell.v[q]>0){
                    TreeNode* i_i_node = i_node->cell.i_nodes[q];
                    th = (i_node->cell.prob-i_i_node->cell.prob)/(r->cell.prob-i_node->cell.prob);
                    t = i_node->cell.v[q];
                }else{
                    th = (r->cell.k_nodes[q]->cell.prob-r->cell.prob)/(r->cell.prob-i_node->cell.prob);
                    t = -i_node->cell.v[q];
                }
                
                i_node->cell.f[q] += t*(G.del[q]/G.dt - t)*F*MC(th); 
            }
        }
    }
};

void GBEES::calculating_K(Grid G, TreeNode* r, int rk){
    if (r == NULL){
        return;
    }
    calculating_K(G,r->left, rk);
    calculating_K(G,r->right, rk);

    uint64_t l_key = r->key; r->cell.K[rk-1] = 0.0; 
    if(l_key != -1){
        for(int q = 0; q < DIM; q++){
            r->cell.K[rk-1] -= (r->cell.f[q]-r->cell.i_nodes[q]->cell.f[q])/G.del[q];  
        }
    }
};

void GBEES::update_prob(Grid G, TreeNode* r, int rk){
    if (r == NULL){
        return; 
    }

    update_prob(G, r->left, rk);
    update_prob(G, r->right, rk);
    
    if(r->key != -1){
        if(rk == 1){
            r->cell.prob = r->cell.x_n + G.dt*r->cell.K[0]; 
        }else if (rk == 2){
            if(G.rk == 2){
                r->cell.prob = r->cell.x_n + (G.dt/2)*(r->cell.K[0]+r->cell.K[1]); 
            }else if(G.rk == 3){
                r->cell.prob = r->cell.x_n + (G.dt/4)*(r->cell.K[0]+r->cell.K[1]); 
            }
        }else if (rk == 3){
            r->cell.prob = r->cell.x_n + (G.dt/6)*(r->cell.K[0]+r->cell.K[1]+4*r->cell.K[2]); 
        }
    }
};

void GBEES::get_size(Grid G, TreeNode* r, int& a_count, int& tot_count, uint64_t& max_key){
    if (r == NULL){
        return; 
    }

    get_size(G, r->left, a_count, tot_count, max_key);
    get_size(G, r->right, a_count, tot_count, max_key);

    if(r->key != -1){
        tot_count += 1; 
        if(r->cell.prob >= G.thresh){
            a_count += 1; 
        }
        if (r->key > max_key){
            max_key = r->key; 
        }
    }
};

void GBEES::Record_Data(std::string file_name, Grid G, TreeNode* r){
	std::ofstream myfile; myfile.open(file_name);
    P.writeFile(myfile, G, P.root);
    myfile.close(); 
};

void GBEES::measurement_update(Measurement m, Grid G, TreeNode* r, double& C){

    if (r == NULL){
        return; 
    }

    measurement_update(m, G, r->left, C);
    measurement_update(m, G, r->right, C);

    double x = 0; 
    for(int q = 0; q < DIM; q++){
        if (m.unc[q] != 0){
            x += pow(((r->cell.state[q]*G.del[q]+G.epoch[q])-m.mean[q]),2)/pow(m.unc[q],2); 
        }
    }

    double prob = exp(-4*x/2); r->cell.prob *= prob; C += r->cell.prob; 
};
#endif // GBEES_H