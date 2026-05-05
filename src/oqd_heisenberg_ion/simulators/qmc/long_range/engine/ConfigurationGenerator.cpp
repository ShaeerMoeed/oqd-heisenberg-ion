#include "ConfigurationGenerator.h"

ConfigurationGenerator::ConfigurationGenerator(const SimulationParameters &sim_params,
                                               const ProbabilityTables &prob_tables,
                                               std::shared_ptr<spdlog::logger> &logger_ptr){

    init_config_index = sim_params.init_config_index;
    hamiltonian_type = sim_params.hamiltonian_type;
    a_parameter = sim_params.new_M_multiplier;

    logger=logger_ptr;

    // Updated and used throughout the simulation if necessary by the required methods
    cumulative_loop_size = -1;
    num_free_spins = -1;
    num_total_legs = -1;

    off_diagonal_update_seed = sim_params.off_diagonal_update_seed;

    // Will be set during equilibration
    if (hamiltonian_type == 2 or hamiltonian_type == 3) {
        average_cumulative_loop_size = 1;
        max_loop_size = 10;
        N_l = 1;
        count_non_skipped_loop_updates = 0;
        skip_loop_update = true;

        metropolis_bond_generator_seed = sim_params.metropolis_bond_generator_seed;
        metropolis_bond_generator.seed(metropolis_bond_generator_seed);

        loop_start_position_generator.seed(off_diagonal_update_seed);

        logger->info("SSE loop parameter initialization:");
        logger->info("hamiltonian_type = {}, average_cumulative_loop_size = {}, max_loop_size = {}, N_l = {}",
            hamiltonian_type, average_cumulative_loop_size, max_loop_size, N_l);
        logger->info("count_non_skipped_loop_updates = {}, skip_loop_update = {}",
            count_non_skipped_loop_updates, skip_loop_update);

    }
    else {
        off_diagonal_spin_flip_generator.seed(off_diagonal_update_seed);
    }

    num_bonds = sim_params.num_bonds;

    beta = sim_params.beta;

    N = sim_params.N;
    if (N % 2 == 0){
        N_over_2 = N/2;
    }
    else {
        N_over_2 = 0;
    }

    equilibration_steps = sim_params.equilibration_steps;
    mc_steps = sim_params.simulation_steps;

    diagonal_update_seed = sim_params.diagonal_update_seed;
    diagonal_update_generator.seed(diagonal_update_seed);

    disconnected_spin_flip_seed = sim_params.disconnected_spin_flip_seed;
    disconnected_spin_flip_generator.seed(disconnected_spin_flip_seed);

    metropolis_insert_seed = sim_params.metropolis_insert_seed;
    metropolis_insert_generator.seed(metropolis_insert_seed);

    metropolis_remove_seed = sim_params.metropolis_remove_seed;
    metropolis_remove_generator.seed(metropolis_remove_seed);
    //twist_seed = sim_params.twist_seed;

    if (hamiltonian_type == 0 or hamiltonian_type == 2 or hamiltonian_type == 3){
        exit_leg_seed = sim_params.exit_leg_seed;
        exit_leg_generator.seed(exit_leg_seed);
    }

    spin_labels = {-1,1};

    ConfigurationGenerator::setInitialConfiguration(sim_params);
    if (init_config_index == -2) {
        operator_locations = sim_params.init_operator_locations;
        n = sim_params.init_n;
        init_M = sim_params.init_M;
        M = init_M;
        num_winding = sim_params.winding;

        logger->info("SSE sampling initialization:");
        logger->info("init_config_index = {}, n = {}, init_M = {}, num_winding = {}",
            init_config_index,n, init_M, num_winding);
    }
    else {
        //init_M = std::max((int)(2 * sim_params.beta * prob_tables.max_diagonal_norm), sim_params.init_M);
        init_M = sim_params.init_M;
        M = init_M;
        n=0;
        ConfigurationGenerator::populateOperatorLocations(M);
        num_winding = 0.0;

        logger->info("SSE sampling initialization:");
        logger->info("init_config_index = {}, n = {}, init_M = {}, num_winding = {}",
            init_config_index, n, init_M, num_winding);
    }

    if (init_config_index == 0) {
        initial_config_seed = sim_params.initial_config_seed;
        initial_config_generator.seed(initial_config_seed);
    }
}

void ConfigurationGenerator::populateOperatorLocations(const int &num_fill_zeros) {

    std::vector<int> operator_locations_new;
    for (int i=0;i<n;i++){
        operator_locations_new.push_back(operator_locations.at(p_list.at(i)));
        p_list[i] = i;
    }
    operator_locations = operator_locations_new;
    for (int i=0; i<num_fill_zeros; i++){
        operator_locations.push_back(0);
    }
}

void ConfigurationGenerator::getNewOperatorLocations(const int &num_fill_zeros) {

    for (int i=0; i<num_fill_zeros; i++){
        operator_locations.push_back(0);
    }
}

void ConfigurationGenerator::setInitialConfiguration(const SimulationParameters &sim_params) {

    if (init_config_index == 0){
        std::uniform_int_distribution<> initial_config_dist(0, 1);
        for (int i=0; i<N; i++){
            spin_configuration.push_back(spin_labels.at(initial_config_dist(initial_config_generator)));
        }
    }
    else if (init_config_index == 1){
        for (int i=0; i<N; i++){
            spin_configuration.push_back(1);
        }
    }
    else if (init_config_index == -1){
        for (int i=0; i<N; i++){
            spin_configuration.push_back(-1);
        }
    }
    else if (init_config_index > 1){
        for (int i=0; i<N; i++){
            if (i % init_config_index == 0){
                spin_configuration.push_back(1);
            }
            else {
                spin_configuration.push_back(-1);
            }
        }
    }
    else if (init_config_index == -2) {
        spin_configuration = sim_params.init_spin_config;
    }
    else {
        logger->error("Invalid initial configuration index provided");
        logger->flush();
        throw std::runtime_error("Invalid initial configuration index provided\n");
    }
}

void ConfigurationGenerator::diagonalUpdatesXXZh(const ProbabilityTables &prob_tables,
                                                 const VertexTypes &vertex_types)
{

    p_list.clear();
    b_list.clear();
    std::discrete_distribution<int> diag_update_distribution(prob_tables.max_norm_probabilities.begin(),
                                                             prob_tables.max_norm_probabilities.end());

    for (int t=0; t<M; t++){
        double M_minus_n = M - n;
        if (operator_locations.at(t) == 0) {
            double u_1 = metropolis_acceptance_distribution(metropolis_insert_generator);
            double acceptance_prob = beta * prob_tables.max_diagonal_norm/M_minus_n;
            if (u_1 < acceptance_prob) {
                int b_prime = diag_update_distribution(diagonal_update_generator);
                int i_b = prob_tables.lattice_sites.at(b_prime).at(0);
                int j_b = prob_tables.lattice_sites.at(b_prime).at(1);
                std::vector<int> diag_config = {spin_configuration.at(i_b), spin_configuration.at(j_b),
                                                spin_configuration.at(i_b), spin_configuration.at(j_b)};
                int vertex_type = vertex_types.getVertexTypeIndex(diag_config);
                double u_2 = metropolis_acceptance_distribution(metropolis_bond_generator);
                double acceptance_prob_2 = prob_tables.diagonal_probabilities.at(vertex_type).at(b_prime);
                if (u_2 < acceptance_prob_2) {
                    int bond_num = b_prime + 1;
                    operator_locations.at(t) = 2 * bond_num;
                    n++;
                    p_list.push_back(t);
                    b_list.push_back(bond_num);
                }
            }
        }
        else if (operator_locations.at(t) % 2 == 1){
            p_list.push_back(t);
            int b = (operator_locations.at(t) - 1)/2;
            b_list.push_back(b);
            int bond_index = b - 1;
            int i_b = prob_tables.lattice_sites.at(bond_index).at(0);
            int j_b = prob_tables.lattice_sites.at(bond_index).at(1);
            spin_configuration.at(i_b) = -spin_configuration.at(i_b);
            spin_configuration.at(j_b) = -spin_configuration.at(j_b);
        }
        else {
            double acceptance_prob = (M_minus_n + 1.0)/(beta * prob_tables.max_diagonal_norm);
            double u_1 = metropolis_acceptance_distribution(metropolis_remove_generator);
            if (u_1 < acceptance_prob) {
                operator_locations.at(t) = 0;
                n--;
            }
            else {
                p_list.push_back(t);
                int b = operator_locations.at(t)/2;
                b_list.push_back(b);
            }
        }
    }
}

void ConfigurationGenerator::diagonalUpdatesXY(const ProbabilityTables &prob_tables)
{

    p_list.clear();
    b_list.clear();

    std::discrete_distribution<int> diag_update_distribution(prob_tables.max_norm_probabilities.begin(),
                                                             prob_tables.max_norm_probabilities.end());

    for (int t=0; t<M; t++){
        double M_minus_n = (double)M - (double)n;
        if (operator_locations.at(t) == 0) {
            double u_1 = metropolis_acceptance_distribution(metropolis_insert_generator);
            double acceptance_prob = beta * prob_tables.max_diagonal_norm/M_minus_n;
            if (u_1 < acceptance_prob) {
                int b_prime = diag_update_distribution(diagonal_update_generator);
                int bond_num = b_prime + 1;
                operator_locations.at(t) = 2*bond_num;
                n++;
                p_list.push_back(t);
                b_list.push_back(bond_num);
            }
        }
        else if (operator_locations.at(t) % 2 == 1){
            p_list.push_back(t);
            int b = (operator_locations.at(t) - 1)/2;
            b_list.push_back(b);
            int bond_index = b - 1;
            int i_b = prob_tables.lattice_sites.at(bond_index).at(0);
            int j_b = prob_tables.lattice_sites.at(bond_index).at(1);
            spin_configuration.at(i_b) = -spin_configuration.at(i_b);
            spin_configuration.at(j_b) = -spin_configuration.at(j_b);
        }
        else {
            double acceptance_prob = (M_minus_n + 1.0)/(beta * prob_tables.max_diagonal_norm);
            double u_1 = metropolis_acceptance_distribution(metropolis_remove_generator);
            if (u_1 < acceptance_prob) {
                operator_locations.at(t) = 0;
                n--;
            }
            else {
                p_list.push_back(t);
                int b = operator_locations.at(t)/2;
                b_list.push_back(b);
            }
        }
    }
}

void ConfigurationGenerator::diagonalUpdatesIsotropic(const ProbabilityTables &prob_tables)
{

    p_list.clear();
    b_list.clear();

    std::discrete_distribution<int> diag_update_distribution(prob_tables.max_norm_probabilities.begin(),
                                                             prob_tables.max_norm_probabilities.end());

    for (int t=0; t<M; t++){
        double M_minus_n = M - n;
        if (operator_locations.at(t) == 0) {
            double u_1 = metropolis_acceptance_distribution(metropolis_insert_generator);
            double acceptance_prob = beta * prob_tables.max_diagonal_norm/M_minus_n;
            if (u_1 < acceptance_prob) {
                int b_prime = diag_update_distribution(diagonal_update_generator);
                int i_b = prob_tables.lattice_sites.at(b_prime).at(0);
                int j_b = prob_tables.lattice_sites.at(b_prime).at(1);
                if (spin_configuration.at(i_b) == hamiltonian_type * spin_configuration.at(j_b)) {
                    int bond_num = b_prime + 1;
                    operator_locations.at(t) = 2*bond_num;
                    n++;
                    p_list.push_back(t);
                    b_list.push_back(bond_num);
                }
            }
        }
        else if (operator_locations.at(t) % 2 == 1){
            p_list.push_back(t);
            int b = (operator_locations.at(t) - 1)/2;
            b_list.push_back(b);
            int bond_index = b - 1;
            int i_b = prob_tables.lattice_sites.at(bond_index).at(0);
            int j_b = prob_tables.lattice_sites.at(bond_index).at(1);
            spin_configuration.at(i_b) = -spin_configuration.at(i_b);
            spin_configuration.at(j_b) = -spin_configuration.at(j_b);
        }
        else {
            double acceptance_prob = (M_minus_n + 1.0)/(beta * prob_tables.max_diagonal_norm);
            double u_1 = metropolis_acceptance_distribution(metropolis_remove_generator);
            if (u_1 < acceptance_prob) {
                operator_locations.at(t) = 0;
                n--;
            }
            else {
                p_list.push_back(t);
                int b = operator_locations.at(t)/2;
                b_list.push_back(b);
            }
        }
    }
}

void ConfigurationGenerator::initializeOffDiagonalUpdates() {

    first_vertex_leg.clear();
    last_vertex_leg.clear();
    linked_list.clear();
    vertex_configuration.clear();

    num_total_legs = num_legs_per_vertex * n;

    for (int i=0; i<N; i++){
        first_vertex_leg.push_back(-1);
        last_vertex_leg.push_back(-1);
    }

    for (int i=0; i<n; i++){
        vertex_configuration.push_back(0);
        for (int j=0; j<num_legs_per_vertex; j++){
            linked_list.push_back(0);
        }
    }
}

void ConfigurationGenerator::populateLinkedList(const ProbabilityTables &prob_tables,
                                                const VertexTypes &vertex_types) {

    ConfigurationGenerator::initializeOffDiagonalUpdates();
    for (int p=0; p<n; p++){
        int t = p_list.at(p);
        int bond_num = b_list.at(p);
        int bond_index = bond_num - 1;
        int i_b = prob_tables.lattice_sites.at(bond_index).at(0);
        int j_b = prob_tables.lattice_sites.at(bond_index).at(1);

        if (last_vertex_leg.at(i_b) == -1) {
            first_vertex_leg.at(i_b) = num_legs_per_vertex*p;
            last_vertex_leg.at(i_b) = num_legs_per_vertex*p + 2;
        }
        else{
            linked_list.at(num_legs_per_vertex*p) = last_vertex_leg.at(i_b);
            linked_list.at(last_vertex_leg.at(i_b)) = num_legs_per_vertex*p;
            last_vertex_leg.at(i_b) = num_legs_per_vertex*p + 2;
        }

        if (last_vertex_leg.at(j_b) == -1) {
            first_vertex_leg.at(j_b) = num_legs_per_vertex*p + 1;
            last_vertex_leg.at(j_b) = num_legs_per_vertex*p + 3;
        }
        else {
            linked_list.at(num_legs_per_vertex*p + 1) = last_vertex_leg.at(j_b);
            linked_list.at(last_vertex_leg.at(j_b)) = num_legs_per_vertex*p + 1;
            last_vertex_leg.at(j_b) = num_legs_per_vertex*p + 3;
        }

        if (operator_locations.at(t) % 2 == 1) {
            std::vector<int> current_config = {spin_configuration.at(i_b), spin_configuration.at(j_b),
                                               -spin_configuration.at(i_b), -spin_configuration.at(j_b)};
            vertex_configuration.at(p) = vertex_types.getVertexTypeIndex(current_config);
            spin_configuration.at(i_b) = -spin_configuration.at(i_b);
            spin_configuration.at(j_b) = -spin_configuration.at(j_b);
        }
        else {
            std::vector<int> current_config = {spin_configuration.at(i_b), spin_configuration.at(j_b),
                                               spin_configuration.at(i_b), spin_configuration.at(j_b)};
            vertex_configuration.at(p) = vertex_types.getVertexTypeIndex(current_config);
        }
    }

    for (int i=0; i<N; i++){
        if (last_vertex_leg.at(i) != -1) {
            linked_list.at(last_vertex_leg.at(i)) = first_vertex_leg.at(i);
            linked_list.at(first_vertex_leg.at(i)) = last_vertex_leg.at(i);
        }
    }
}

void ConfigurationGenerator::offDiagonalUpdatesXXZh(const ProbabilityTables &prob_tables,
                                                    const VertexTypes &vertex_types) {

    if (n == 0) {
        num_free_spins = 0;
        for (int i = 0; i < N; i++) {
            flipFreeSpin(i);
        }
    }
    else {
        populateLinkedList(prob_tables, vertex_types);
        loop_start_pos_dist = std::uniform_int_distribution<>(0, num_total_legs - 1);
        int loop_size = 0;
        cumulative_loop_size = 0;

        for (int loop_num=0; loop_num<N_l; loop_num++){

            skip_loop_update = true;
            loop_size = 0;

            int j_0 = loop_start_pos_dist(loop_start_position_generator);
            int j_current = j_0;

            for (int iter=0;iter<max_loop_size;iter++){

                int l_e = j_current % num_legs_per_vertex;
                int p = (j_current - l_e)/num_legs_per_vertex;
                int vertex_type = vertex_configuration.at(p);

                int composite_leg_index = num_legs_per_vertex*l_e;
                int row_index = num_composite_leg_indices*vertex_type + composite_leg_index;
                int loc_index = num_legs_per_vertex*(num_legs_per_vertex-1)*vertex_type
                        + (num_legs_per_vertex-1)*l_e;

                out_leg_probs.clear();
                for (int k=0; k<num_legs_per_vertex-1; k++){
                    int l_k = vertex_types.allowed_exit_legs.at(loc_index + k);
                    double prob_l_k = prob_tables.loop_update_probabilities.at(
                            row_index + l_k).at(b_list.at(p)-1);
                    out_leg_probs.push_back(prob_l_k);
                }

                std::discrete_distribution<int> exit_leg_dist(out_leg_probs.begin(),out_leg_probs.end());
                int l_x_index = exit_leg_dist(exit_leg_generator);
                int l_x = vertex_types.allowed_exit_legs.at(loc_index + l_x_index);
                vertex_configuration.at(p) = vertex_types.getFlippedSpinsVertexIndex(l_e, l_x,
                                                                                     vertex_type);

                j_current = num_legs_per_vertex*p + l_x;
                if (l_e != l_x) {
                    loop_size++;
                }
                if (j_current == j_0){
                    skip_loop_update = false;
                    cumulative_loop_size += loop_size;
                    break;
                }
                else {
                    j_current = linked_list.at(j_current);
                    if (j_current == j_0) {
                        skip_loop_update = false;
                        cumulative_loop_size += loop_size;
                        break;
                    }
                }
            }

            if (skip_loop_update) {
                cumulative_loop_size = 0;
                break;
            }
        }

        num_free_spins = 0;
        if (!skip_loop_update) {
            mapVerticesToOperatorLocations(vertex_types, prob_tables);
        }
        else {
            for (int i=0; i<N; i++) {
                if (first_vertex_leg.at(i) == -1) {
                    flipFreeSpin(i);
                }
            }
        }
    }
}

void ConfigurationGenerator::offDiagonalUpdatesXY(const ProbabilityTables &prob_tables,
                                                  const VertexTypes &vertex_types) {

    if (n == 0) {
        num_free_spins = 0;
        for (int i = 0; i < N; i++) {
            flipFreeSpin(i);
        }
    }
    else {
        populateLinkedList(prob_tables, vertex_types);
        bool loop_closed = false;
        for (int j_0 = 0; j_0 < num_total_legs; j_0++) {

            if (!connected_legs.contains(j_0)) {

                int j = j_0;
                int flip_spins = spin_labels.at(binary_dist(off_diagonal_spin_flip_generator));

                while (!loop_closed) {

                    connected_legs.insert(j);

                    int p = j / num_legs_per_vertex;
                    int l_e = j % num_legs_per_vertex;

                    int vertex_type = vertex_configuration.at(p);
                    int exit_leg_index = binary_dist(exit_leg_generator);
                    int loc_index = num_legs_per_vertex*(num_legs_per_vertex-2)*vertex_type
                            + (num_legs_per_vertex-2)*l_e;
                    int l_x = vertex_types.allowed_exit_legs.at(loc_index + exit_leg_index);

                    if (flip_spins == 1) {
                        vertex_configuration.at(p) = vertex_types.getFlippedSpinsVertexIndex(l_e, l_x,
                                                                                 vertex_type);
                    }

                    j = num_legs_per_vertex*p + l_x;
                    connected_legs.insert(j);
                    j = linked_list.at(j);

                    if (j == j_0) {
                        loop_closed = true;
                    }
                }
            }
        }
        connected_legs.clear();
        mapVerticesToOperatorLocations(vertex_types, prob_tables);
    }
}

void ConfigurationGenerator::offDiagonalUpdatesIsotropic(const ProbabilityTables &prob_tables,
                                                         const VertexTypes &vertex_types) {

    if (n == 0) {
        num_free_spins = 0;
        for (int i = 0; i < N; i++) {
            flipFreeSpin(i);
        }
    }
    else {
        populateLinkedList(prob_tables, vertex_types);
        bool loop_closed = false;
        for (int j_0 = 0; j_0 < num_total_legs; j_0++) {

            if (!connected_legs.contains(j_0)) {

                int j = j_0;
                int flip_spins = spin_labels.at(binary_dist(off_diagonal_spin_flip_generator));

                while (!loop_closed) {

                    connected_legs.insert(j);

                    int p = j / num_legs_per_vertex;
                    int l_e = j % num_legs_per_vertex;

                    int vertex_type = vertex_configuration.at(p);
                    int l_x = vertex_types.allowed_exit_legs.at(l_e);

                    if (flip_spins == 1) {
                        vertex_configuration.at(p) = vertex_types.getFlippedSpinsVertexIndex(l_e, l_x,
                                                                                 vertex_type);
                    }

                    j = num_legs_per_vertex*p + l_x;
                    connected_legs.insert(j);
                    j = linked_list.at(j);

                    if (j == j_0) {
                        loop_closed = true;
                    }
                }
            }
        }
        connected_legs.clear();
        mapVerticesToOperatorLocations(vertex_types, prob_tables);
    }

}

void ConfigurationGenerator::mapVerticesToOperatorLocations(const VertexTypes &vertex_types,
                                                            const ProbabilityTables &prob_tables){

    count_non_skipped_loop_updates++;
    num_winding = 0.0;
    for (int p = 0; p < n; p++) {
        int t = p_list.at(p);
        int b = b_list.at(p);
        operator_locations.at(t) = 2 * b + vertex_types.is_off_diag.at(vertex_configuration.at(p));
        int signed_distance = prob_tables.lattice_sites.at(b-1).at(2);
        num_winding += vertex_types.twist_mapping.at(vertex_configuration.at(p)) * signed_distance;
    }
    for (int i = 0; i < N; i++) {
        if (first_vertex_leg.at(i) == -1) {
            flipFreeSpin(i);
        } else {
            int l = first_vertex_leg.at(i) % num_legs_per_vertex;
            int p = (first_vertex_leg.at(i) - l) / num_legs_per_vertex;
            int vertex_type = vertex_configuration.at(p);
            std::vector<int> config = vertex_types.getVertexConfig(vertex_type);
            spin_configuration.at(i) = config.at(l);
        }
    }
}

void ConfigurationGenerator::flipFreeSpin(const int &i) {

    int multiplier = spin_labels.at(binary_dist(disconnected_spin_flip_generator));
    spin_configuration.at(i) = multiplier * spin_configuration.at(i);
    num_free_spins++;
}

void ConfigurationGenerator::randomSpinFlipsXXZ() {

    int flip_spins_bool = binary_dist(disconnected_spin_flip_generator);
    if (flip_spins_bool == 1){
        for (int j=0; j<N; j++){
            spin_configuration.at(j) = -1 * spin_configuration.at(j);
        }
    }
}

void ConfigurationGenerator::flipAllSpins() {
    for (int j=0; j<N; j++){
        spin_configuration.at(j) = -1 * spin_configuration.at(j);
    }
}

void ConfigurationGenerator::simulateProbabilisticLoopsXXZh(const ProbabilityTables &prob_tables,
                                                            Estimators &estimators,
                                                            const SimulationParameters &sim_params,
                                                            const VertexTypes &vertex_types) {

    int max_n = n;
    int avg_n = 0;
    double avg_num_samples = (double)std::min(1000, equilibration_steps);
    int n_sum = 0;
    int cumulative_loop_size_sum = 1;
    int M_new;

    double avg_n_d = 0.0;
    double avg_cumul_loop_size_d;
    double prob_non_skip_update;
    double N_l_d;

    logger->info("Starting burn in");

    for (int step=0; step<equilibration_steps; step++) {
        ConfigurationGenerator::diagonalUpdatesXXZh(prob_tables, vertex_types);
        if (n > max_n) {
            max_n = n;
        }
        n_list.push_back(n);
        n_sum += n;

        avg_n_d = n_sum / ((double)step+1.0);
        avg_n = (int)std::round(avg_n_d);

        //M_new = std::max((int)std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d),n+1);
        //ConfigurationGenerator::populateOperatorLocations(int(M_new - n));
        M_new = std::max((int)std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d),M);
        ConfigurationGenerator::getNewOperatorLocations(M_new-M);
        M = M_new;

        max_loop_size = 10000 * std::max(avg_n, 1);
        ConfigurationGenerator::offDiagonalUpdatesXXZh(prob_tables, vertex_types);

        cumulative_loop_size_sum += cumulative_loop_size;
        cumulative_loop_size_list.push_back(cumulative_loop_size);
        avg_cumul_loop_size_d = ((double)cumulative_loop_size_sum)/((double)step+1.0);
        prob_non_skip_update = (double)count_non_skipped_loop_updates/((double)step+1.0);
        N_l_d = 2.0*avg_n_d/(avg_cumul_loop_size_d);
        N_l = std::max((int)std::round(N_l_d), 1);

        if ((double)cumulative_loop_size_list.size() == avg_num_samples){
            break;
        }

        if ((step+1) % 1000 == 0) {
            logger->info("Burn-in step = {}, n = {}, M = {}, N_l = {}",
                step, n, M, N_l);
            logger->flush();
        }
    }

    logger->info("Burn in finished");

    avg_cumul_loop_size_d = ConfigurationGenerator::computeAverage(cumulative_loop_size_list,
                                                                   avg_num_samples);
    avg_n_d = ConfigurationGenerator::computeAverage(n_list, avg_num_samples);
    prob_non_skip_update = (double)count_non_skipped_loop_updates/(avg_num_samples);

    avg_n = (int)std::round(avg_n_d);
    average_cumulative_loop_size = (int)std::round(avg_cumul_loop_size_d);

    //M_new = std::max((int)std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d),n+1);
    //ConfigurationGenerator::populateOperatorLocations(int(M_new - n));
    M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), M);
    ConfigurationGenerator::getNewOperatorLocations(M_new-M);

    M = M_new;
    count_non_skipped_loop_updates = 0;

    N_l_d = 2.0*avg_n_d/(avg_cumul_loop_size_d * prob_non_skip_update);
    N_l = std::max((int)std::round(N_l_d), 1);

    logger->info("Starting equilibration");

    for (int step=0; step<equilibration_steps; step++){
        ConfigurationGenerator::diagonalUpdatesXXZh(prob_tables, vertex_types);

        if (n > max_n) {
            max_n = n;
        }

        avg_n_d = ((avg_n_d * avg_num_samples) - n_list.at(n_list.size() - (int)avg_num_samples)
                + n)/avg_num_samples;
        avg_n = (int)std::round(avg_n_d);
        n_list.push_back(n);

        //M_new = std::max((int)std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d),n+1);
        //ConfigurationGenerator::populateOperatorLocations(int(M_new - n));
        M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), M);
        ConfigurationGenerator::getNewOperatorLocations(M_new-M);
        M = M_new;

        max_loop_size = std::max(10000*avg_n,1);

        ConfigurationGenerator::offDiagonalUpdatesXXZh(prob_tables, vertex_types);

        avg_cumul_loop_size_d = ((avg_cumul_loop_size_d * avg_num_samples)
                                 - cumulative_loop_size_list.at(cumulative_loop_size_list.size()
                                 - (int)avg_num_samples)
                                 + cumulative_loop_size)/avg_num_samples;

        average_cumulative_loop_size = (int)std::round(avg_cumul_loop_size_d);
        cumulative_loop_size_list.push_back(cumulative_loop_size);
        prob_non_skip_update = (double)count_non_skipped_loop_updates/((double)step+1.0);
        N_l_d = 2.0*avg_n_d/(avg_cumul_loop_size_d * prob_non_skip_update);
        N_l = std::max((int)std::round(N_l_d), 1);

        if ((step+1) % 1000 == 0) {
            logger->info("Equilibration step = {}, n = {}, M = {}, N_l = {}",
                step, n, M, N_l);
            logger->flush();
        }

    }

    logger->info("Equilibration finished");
    logger->flush();

    max_loop_size = std::max(10000*avg_n,1);

    avg_cumul_loop_size_d = ((avg_cumul_loop_size_d * avg_num_samples)
                             - cumulative_loop_size_list.at(cumulative_loop_size_list.size()
                             - (int)avg_num_samples)
                             + cumulative_loop_size)/avg_num_samples;
    average_cumulative_loop_size = (int)std::round(avg_cumul_loop_size_d);
    cumulative_loop_size_list.push_back(cumulative_loop_size);
    prob_non_skip_update = (double)count_non_skipped_loop_updates/((double)equilibration_steps);
    N_l_d = 2.0*avg_n_d/(avg_cumul_loop_size_d * prob_non_skip_update);
    N_l = std::max((int)std::round(N_l_d), 1);

    logger->info("Starting estimation run");
    logger->flush();

    for (int step=0; step<mc_steps; step++){

        ConfigurationGenerator::diagonalUpdatesXXZh(prob_tables, vertex_types);

        ConfigurationGenerator::offDiagonalUpdatesXXZh(prob_tables, vertex_types);

        estimators.updateAllPropertiesProbabilistic(n, spin_configuration, sim_params,
                                       prob_tables.spectrum_offset, num_winding,
                                       skip_loop_update, prob_tables.lattice_sites);

        if (sim_params.track_spin_configs) {
            estimators.trackSpinConfigs(spin_configuration, sim_params);
        }

        if ((step+1) % 1000 == 0) {
            logger->info("Simulation step = {}, n = {}", step, n);
            logger->flush();
        }

        if ((step+1) % estimators.chunk_size == 0) {
            estimators.populateShotDataFile(sim_params);
        }
    }

    logger->info("Estimation run finished");
}

void ConfigurationGenerator::simulateProbabilisticLoopsXXZ(const ProbabilityTables &prob_tables,
                                                           Estimators &estimators,
                                                           const SimulationParameters &sim_params,
                                                           const VertexTypes &vertex_types) {

    int max_n = n;
    int avg_n = 0;
    double avg_num_samples = (double)std::min(1000, equilibration_steps);
    int n_sum = 0;
    int cumulative_loop_size_sum = 1;
    int M_new;

    double avg_n_d = 0.0;
    double avg_cumul_loop_size_d;
    double prob_non_skip_update;
    double N_l_d;

    logger->info("Starting burn in");

    for (int step=0; step<equilibration_steps; step++) {
        ConfigurationGenerator::diagonalUpdatesXXZh(prob_tables, vertex_types);
        if (n > max_n) {
            max_n = n;
        }
        n_list.push_back(n);
        n_sum += n;

        avg_n_d = n_sum / ((double)step+1.0);
        avg_n = (int)std::round(avg_n_d);

        //M_new = std::max((int)std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d),n+1);
        //ConfigurationGenerator::populateOperatorLocations(int(M_new - n));
        M_new = std::max((int)std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d),M);
        ConfigurationGenerator::getNewOperatorLocations(M_new-M);
        M = M_new;

        max_loop_size = 10000 * std::max(avg_n, 1);

        ConfigurationGenerator::offDiagonalUpdatesXXZh(prob_tables, vertex_types);
        ConfigurationGenerator::randomSpinFlipsXXZ();

        cumulative_loop_size_sum += cumulative_loop_size;
        cumulative_loop_size_list.push_back(cumulative_loop_size);
        avg_cumul_loop_size_d = ((double)cumulative_loop_size_sum)/((double)step+1.0);
        prob_non_skip_update = (double)count_non_skipped_loop_updates/((double)step+1.0);
        N_l_d = 2.0*avg_n_d/(avg_cumul_loop_size_d);
        N_l = std::max((int)std::round(N_l_d), 1);

        if ((double)cumulative_loop_size_list.size() == avg_num_samples){
            break;
        }

        if ((step+1) % 1000 == 0) {
            logger->info("Burn-in step = {}, n = {}, M = {}, N_l = {}",
                step, n, M, N_l);
            logger->flush();
        }
    }

    logger->info("Burn in finished");

    avg_cumul_loop_size_d = ConfigurationGenerator::computeAverage(cumulative_loop_size_list,
                                                                   avg_num_samples);
    avg_n_d = ConfigurationGenerator::computeAverage(n_list, avg_num_samples);
    prob_non_skip_update = (double)count_non_skipped_loop_updates/(avg_num_samples);

    avg_n = (int)std::round(avg_n_d);
    average_cumulative_loop_size = (int)std::round(avg_cumul_loop_size_d);

    //M_new = std::max((int)std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d),n+1);
    //ConfigurationGenerator::populateOperatorLocations(int(M_new - n));
    M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), M);
    ConfigurationGenerator::getNewOperatorLocations(M_new-M);

    M = M_new;
    count_non_skipped_loop_updates = 0;

    N_l_d = 2.0*avg_n_d/(avg_cumul_loop_size_d * prob_non_skip_update);
    N_l = std::max((int)std::round(N_l_d), 1);

    logger->info("Starting equilibration");

    for (int step=0; step<equilibration_steps; step++){
        ConfigurationGenerator::diagonalUpdatesXXZh(prob_tables, vertex_types);

        if (n > max_n) {
            max_n = n;
        }

        avg_n_d = ((avg_n_d * avg_num_samples) - n_list.at(n_list.size() - (int)avg_num_samples)
                + n)/avg_num_samples;
        avg_n = (int)std::round(avg_n_d);
        n_list.push_back(n);

        //M_new = std::max((int)std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d),n+1);
        //ConfigurationGenerator::populateOperatorLocations(int(M_new - n));
        M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), M);
        ConfigurationGenerator::getNewOperatorLocations(M_new-M);
        M = M_new;

        max_loop_size = std::max(10000*avg_n,1);

        ConfigurationGenerator::offDiagonalUpdatesXXZh(prob_tables, vertex_types);
        ConfigurationGenerator::randomSpinFlipsXXZ();

        avg_cumul_loop_size_d = ((avg_cumul_loop_size_d * avg_num_samples)
                                 - cumulative_loop_size_list.at(cumulative_loop_size_list.size()
                                 - (int)avg_num_samples)
                                 + cumulative_loop_size)/avg_num_samples;

        average_cumulative_loop_size = (int)std::round(avg_cumul_loop_size_d);
        cumulative_loop_size_list.push_back(cumulative_loop_size);
        prob_non_skip_update = (double)count_non_skipped_loop_updates/((double)step+1.0);
        N_l_d = 2.0*avg_n_d/(avg_cumul_loop_size_d * prob_non_skip_update);
        N_l = std::max((int)std::round(N_l_d), 1);

        if ((step+1) % 1000 == 0) {
            logger->info("Equilibration step = {}, n = {}, M = {}, N_l = {}",
                step, n, M, N_l);
            logger->flush();
        }

    }

    logger->info("Equilibration finished");
    logger->flush();

    max_loop_size = std::max(10000*avg_n,1);

    avg_cumul_loop_size_d = ((avg_cumul_loop_size_d * avg_num_samples)
                             - cumulative_loop_size_list.at(cumulative_loop_size_list.size()
                             - (int)avg_num_samples)
                             + cumulative_loop_size)/avg_num_samples;
    average_cumulative_loop_size = (int)std::round(avg_cumul_loop_size_d);
    cumulative_loop_size_list.push_back(cumulative_loop_size);
    prob_non_skip_update = (double)count_non_skipped_loop_updates/((double)equilibration_steps);
    N_l_d = 2.0*avg_n_d/(avg_cumul_loop_size_d * prob_non_skip_update);
    N_l = std::max((int)std::round(N_l_d), 1);

    logger->info("Starting estimation run");
    logger->flush();

    for (int step=0; step<mc_steps; step++){

        ConfigurationGenerator::diagonalUpdatesXXZh(prob_tables, vertex_types);

        ConfigurationGenerator::offDiagonalUpdatesXXZh(prob_tables, vertex_types);

        estimators.updateAllPropertiesProbabilistic(n, spin_configuration, sim_params,
                                       prob_tables.spectrum_offset, num_winding,
                                       skip_loop_update, prob_tables.lattice_sites);
        ConfigurationGenerator::randomSpinFlipsXXZ();

        if (sim_params.track_spin_configs) {
            estimators.trackSpinConfigs(spin_configuration, sim_params);
        }

        if ((step+1) % 1000 == 0) {
            logger->info("Simulation step = {}, n = {}", step, n);
            logger->flush();
        }

        if ((step+1) % estimators.chunk_size == 0) {
            estimators.populateShotDataFile(sim_params);
        }
    }

    logger->info("Estimation run finished");
}

void ConfigurationGenerator::simulateDeterministicIsotropic(const ProbabilityTables &prob_tables,
                                                           Estimators &estimators,
                                                           const SimulationParameters &sim_params,
                                                           const VertexTypes &vertex_types) {

    if (hamiltonian_type != 1 and hamiltonian_type != -1) {
        logger->error("hamiltonian type is not 1 or -1");
        logger->flush();
        throw std::runtime_error("hamiltonian type is not 1 or -1\n");
    }

    int max_n = n;
    int n_sum = 0;
    double avg_n_d = 0.0;
    int avg_num_samples = 1000;
    int M_new;

    logger->info("Starting burn in");

    for (int step=0; step<avg_num_samples; step++) {
        ConfigurationGenerator::diagonalUpdatesIsotropic(prob_tables);
        if (n > max_n) {
            max_n = n;
        }
        n_list.push_back(n);
        n_sum += n;

        avg_n_d = n_sum / ((double) step + 1.0);

        //M_new = std::max((int)std::round(1.25 * max_n), M);
        //M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), n + 1);
        //ConfigurationGenerator::populateOperatorLocations(int(M_new - n));
        M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), M);
        ConfigurationGenerator::getNewOperatorLocations(M_new-M);
        M = M_new;

        ConfigurationGenerator::offDiagonalUpdatesIsotropic(prob_tables, vertex_types);
        ConfigurationGenerator::flipAllSpins();

        if ((step+1) % 1000 == 0) {
            logger->info("Equilibration step = {}, n = {}, M = {}", step, n, M);
            logger->flush();
        }
    }

    logger->info("Burn in finished");

    avg_n_d = ConfigurationGenerator::computeAverage(n_list, avg_num_samples);

    //M_new = std::max((int)std::round(1.25 * max_n), M);
    //M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), n + 1);
    ConfigurationGenerator::populateOperatorLocations(int(M_new - n));
    M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), M);
    ConfigurationGenerator::getNewOperatorLocations(M_new-M);
    M = M_new;

    logger->info("Starting equilibration");

    for (int step=0; step<equilibration_steps; step++) {
        ConfigurationGenerator::diagonalUpdatesIsotropic(prob_tables);

        if (n > max_n) {
            max_n = n;
        }

        avg_n_d = ((avg_n_d * avg_num_samples) - n_list.at(n_list.size() - (int) avg_num_samples) + n)
                / avg_num_samples;
        n_list.push_back(n);

        //M_new = std::max((int)std::round(1.25 * max_n), M);
        //M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), n + 1);
        //ConfigurationGenerator::populateOperatorLocations(int(M_new - n));
        M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), M);
        ConfigurationGenerator::getNewOperatorLocations(M_new-M);
        M = M_new;

        ConfigurationGenerator::offDiagonalUpdatesXXZh(prob_tables, vertex_types);

        ConfigurationGenerator::flipAllSpins();

        if ((step+1) % 1000 == 0) {
            logger->info("Equilibration step = {}, n = {}, M = {}", step, n, M);
            logger->flush();
        }
    }

    logger->info("Equilibration finished");
    logger->info("Starting estimation run");
    logger->flush();

    for (int step=0; step<mc_steps; step++){

        ConfigurationGenerator::diagonalUpdatesIsotropic(prob_tables);

        ConfigurationGenerator::offDiagonalUpdatesIsotropic(prob_tables, vertex_types);

        estimators.updateAllPropertiesDeterministic(n, spin_configuration, sim_params,
                                                    prob_tables.spectrum_offset, num_winding,
                                                    prob_tables.lattice_sites);

        ConfigurationGenerator::flipAllSpins();

        estimators.updateAllPropertiesDeterministic(n, spin_configuration, sim_params,
                                                    prob_tables.spectrum_offset, num_winding,
                                                    prob_tables.lattice_sites);

        if ((step+1) % 1000 == 0) {
            logger->info("Simulation step = {}, n = {}", step, n);
            logger->flush();
        }

        if (sim_params.track_spin_configs) {
            estimators.trackSpinConfigs(spin_configuration, sim_params);
        }

        if ((step+1) % estimators.chunk_size == 0) {
            estimators.populateShotDataFile(sim_params);
        }
    }
    logger->info("Estimation run finished");
}

void ConfigurationGenerator::simulateDeterministicXY(const ProbabilityTables &prob_tables,
                                                            Estimators &estimators,
                                                            const SimulationParameters &sim_params,
                                                            const VertexTypes &vertex_types) {

    if (hamiltonian_type != 0) {
        logger->error("hamiltonian type is not 0");
        logger->flush();
        throw std::runtime_error("hamiltonian type is not 0\n");
    }

    int max_n = n;
    int n_sum = 0;
    double avg_n_d = 0.0;
    int avg_num_samples = 100;
    int M_new;

    logger->info("Starting burn in");

    for (int step=0; step<5000; step++) {
        ConfigurationGenerator::diagonalUpdatesXY(prob_tables);
        if (n > max_n) {
            max_n = n;
        }
        n_list.push_back(n);
        n_sum += n;

        avg_n_d = n_sum / ((double) step + 1.0);

        //M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), M);
        //ConfigurationGenerator::getNewOperatorLocations(M_new-M);
        //M = M_new;

        M_new = std::max((int)std::round(a_parameter * avg_n_d), M);
        //M_new = std::max((int)std::round(a_parameter * max_n), M);
        //M_new = std::max((int)std::round(a_parameter*(avg_n_d) + beta*prob_tables.max_diagonal_norm), M);
        if (M_new > M) {
            ConfigurationGenerator::getNewOperatorLocations(M_new-M);
        }
        M = M_new;

        //M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), n + 1);
        //ConfigurationGenerator::populateOperatorLocations(int(M_new - n));

        ConfigurationGenerator::offDiagonalUpdatesXY(prob_tables, vertex_types);
        ConfigurationGenerator::flipAllSpins();

        if (step % avg_num_samples == 0) {
            logger->info("Equilibration step = {}, n = {}, M = {}", step, n, M);
            logger->flush();
        }
    }

    logger->info("Burn in finished");

    avg_n_d = ConfigurationGenerator::computeAverage(n_list, avg_num_samples);

    M_new = std::max((int)std::round(a_parameter * avg_n_d), M);
    //M_new = std::max((int)std::round(a_parameter*(avg_n_d) + beta*prob_tables.max_diagonal_norm), M);
    if (M_new > M) {
        ConfigurationGenerator::getNewOperatorLocations(M_new-M);
    }
    M = M_new;

    //M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), M);
    //ConfigurationGenerator::getNewOperatorLocations(M_new-M);
    //M = M_new;

    //M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), n + 1);
    //ConfigurationGenerator::populateOperatorLocations(int(M_new - n));

    logger->info("Starting equilibration");

    for (int step=0; step<equilibration_steps; step++) {
        ConfigurationGenerator::diagonalUpdatesXY(prob_tables);

        if (n > max_n) {
            max_n = n;
        }

        avg_n_d = ((avg_n_d * avg_num_samples) - (double)n_list.at(n_list.size()
                - (int)avg_num_samples) + (double)n) / avg_num_samples;
        n_list.push_back(n);

        //M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), M);
        //ConfigurationGenerator::populateOperatorLocations(int(M_new - n));
        //ConfigurationGenerator::getNewOperatorLocations(M_new-M);
        //M = M_new;

        M_new = std::max((int)std::round(a_parameter * avg_n_d), M);
        //M_new = std::max((int)std::round(a_parameter*(avg_n_d) + beta*prob_tables.max_diagonal_norm), M);
        if (M_new > M) {
            ConfigurationGenerator::getNewOperatorLocations(M_new-M);
        }
        M = M_new;

        //M_new = std::max((int) std::round(sim_params.beta * prob_tables.max_diagonal_norm + avg_n_d), M);
        //ConfigurationGenerator::getNewOperatorLocations(M_new-M);
        //M = M_new;

        ConfigurationGenerator::offDiagonalUpdatesXY(prob_tables, vertex_types);

        ConfigurationGenerator::flipAllSpins();

        if ((step+1) % 1000 == 0) {
            logger->info("Equilibration step = {}, n = {}, M = {}", step, n, M);
            logger->flush();
        }
    }

    logger->info("Equilibration finished");
    logger->info("Starting estimation run");
    logger->flush();

    for (int step=0; step<mc_steps; step++){

        ConfigurationGenerator::diagonalUpdatesXY(prob_tables);

        ConfigurationGenerator::offDiagonalUpdatesXY(prob_tables, vertex_types);

        estimators.updateAllPropertiesDeterministic(n, spin_configuration, sim_params,
                                                    prob_tables.spectrum_offset, num_winding,
                                                    prob_tables.lattice_sites);

        ConfigurationGenerator::flipAllSpins();

        estimators.updateAllPropertiesDeterministic(n, spin_configuration, sim_params,
                                                    prob_tables.spectrum_offset, num_winding,
                                                    prob_tables.lattice_sites);

        if (sim_params.track_spin_configs) {
            estimators.trackSpinConfigs(spin_configuration, sim_params);
        }

        if ((step+1) % 1000 == 0) {
            logger->info("Simulation step = {}, n = {}", step, n);
            logger->flush();
        }

        if ((step+1) % estimators.chunk_size == 0) {
            estimators.populateShotDataFile(sim_params);
        }
    }
    logger->info("Estimation run finished");
}

double ConfigurationGenerator::computeAverage(std::vector<int> &vector_entries,
                                              const double &num_samples_in) {

    double sum_entries = 0.0;
    int num_samples = (int)num_samples_in;
    for (int i=0; i<num_samples; i++){
        sum_entries += vector_entries.at(vector_entries.size()-i-1);
    }
    sum_entries = sum_entries/num_samples_in;
    double avg_val = sum_entries;

    return avg_val;
}

void ConfigurationGenerator::generateConfigurations(const ProbabilityTables &prob_tables, Estimators &estimators,
                                                    const SimulationParameters &sim_params,
                                                    const VertexTypes &vertex_types) {

    if (sim_params.hamiltonian_type == 0) {
        simulateDeterministicXY(prob_tables, estimators, sim_params, vertex_types);
    }
    else if (sim_params.hamiltonian_type == 1 || sim_params.hamiltonian_type == -1) {
        simulateDeterministicIsotropic(prob_tables, estimators, sim_params, vertex_types);
    }
    else if (sim_params.hamiltonian_type == 2) {
        simulateProbabilisticLoopsXXZ(prob_tables, estimators, sim_params, vertex_types);
    }
    else if (sim_params.hamiltonian_type == 3) {
        simulateProbabilisticLoopsXXZh(prob_tables, estimators, sim_params, vertex_types);
    }
    else {
        logger->error("Hamiltonian type can not be resolved");
        logger->flush();
        throw std::runtime_error("Hamiltonian type can not be resolved.\n");
    }
}

void ConfigurationGenerator::writeFinalConfigurations(const SimulationParameters &sim_params) {

    std::string file_path = sim_params.simulation_subfolder + "/" + "Final SSE Configuration.txt";
    logger->info("Writing final configuration to: {}", file_path);
    std::ofstream ofs(file_path);

    ofs << "N" << "\t" << N;
    ofs << "M" << "\t" << M;
    ofs << "n" << "\t" << n;
    ofs << "W" << "\t" << num_winding;

    ofs << "Operator Locations List" << "\t" << operator_locations.at(0);
    for (int i=1; i < M; i++) {
        ofs << "," << operator_locations.at(i);
    }
    ofs << "\n";

    ofs << "Spin Configurations List" << "\t" << spin_configuration.at(0);
    for (int i=1; i < N; i++) {
        ofs << "," << spin_configuration.at(i);
    }
    ofs << "\n";

    ofs.close();
}