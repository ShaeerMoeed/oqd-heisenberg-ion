#ifndef cpp_qmc_SIMULATIONPARAMETERS_H
#define cpp_qmc_SIMULATIONPARAMETERS_H
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <map>
#include <spdlog/spdlog.h>

class SimulationParameters {

public:
    int N;
    double T;
    double J;
    double beta;
    double gamma;
    double xy_alpha;
    double zz_alpha;
    double ksi;
    double Delta;
    double h;
    double h_B;
    int num_bonds;
    std::shared_ptr<spdlog::logger> logger;

    int hamiltonian_type;
    std::string xy_interaction_type;
    std::string zz_interaction_type;

    std::string simulation_subfolder;
    std::string root_folder;

    int simulation_steps;
    int equilibration_steps;
    double new_M_multiplier;

    std::string loop_type;
    bool distance_dep_offset;

    //int init_M;
    int init_config_index;
    int initial_init_config_index;

    int init_M;
    int init_n;

    std::string boundary_conditions;

    bool track_spin_configs;
    bool write_final_SSE_configs;

    std::string uuid;

    std::string init_config_file_path;

    std::vector<int> init_spin_config;
    std::vector<int> init_operator_locations;
    double winding;

    uint64_t initial_config_seed;
    uint64_t diagonal_update_seed;
    uint64_t exit_leg_seed;
    uint64_t disconnected_spin_flip_seed;
    uint64_t off_diagonal_update_seed;
    uint64_t metropolis_insert_seed;
    uint64_t metropolis_bond_generator_seed;
    uint64_t metropolis_remove_seed;

    explicit SimulationParameters(std::map<std::string, std::string> &input_key_vals,
        const std::shared_ptr<spdlog::logger> &logger_ptr);

    void extractIntegerEntry(const std::string &key_str, const std::string &val_str, int &member_var,
                             const bool &enforce_minimum, const int &min_val=0) const;

    void extractIntegerEntry(const std::string &key_str, const std::string &val_str, uint64_t &member_var,
                             const bool &enforce_minimum, const uint64_t &min_val=0) const;

    void setOptionalIntegerEntry(const std::string &key_str, const std::string &val_str,
                                int &member_var, const bool &enforce_minimum,
                                const int &min_val, const int &default_val) const;

    void setOptionalIntegerEntry(const std::string &key_str, const std::string &val_str,
                                uint64_t &member_var, const bool &enforce_minimum,
                                const uint64_t &min_val, const uint64_t &default_val) const;

     void extractDoubleEntry(const std::string &key_str, const std::string &val_str, double &member_var,
                            const bool &enforce_minimum, const double &min_val=0.0);

     void extractBoolEntry(const std::string &key_str, const std::string &val_str, bool &member_var);

     void extractStringEntry(const std::string &key_str, const std::string &val_str, std::string &member_var);

     void extractListInts(const std::string &key_str, const std::string &val_str, std::vector<int> &member_var,
                         const int &list_size, const bool &enforce_minimum, const int &min_val=0);

    void extractHamiltonianType(const std::string &key_str, const std::string &val_str);

    void extractBoundaryConditions(const std::string &key_str, const std::string &val_str);

    void extractLoopType(const std::string &key_str, const std::string &val_str);

    void extractInteractionType(const std::string &key_str, const std::string &val_str, std::string &member_var);

    void extractInitialConditionsFromFile(std::string &file_path);

    void writeNumericEntry(const std::string &key_str, const int &val, std::ofstream &file_stream) const;

    void writeNumericEntry(const std::string &key_str, const uint64_t &val, std::ofstream &file_stream) const;

    void writeNumericEntry(const std::string &key_str, const double &val, std::ofstream &file_stream) const;

    void writeStringEntry(const std::string &key_str, const std::string &val, std::ofstream &file_stream) const;

    void writeBoolEntry(const std::string &key_str, const bool &val, std::ofstream &file_stream) const;

};


#endif //cpp_qmc_SIMULATIONPARAMETERS_H
