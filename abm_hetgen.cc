/*
 * Highly heterogeneous simulations of infectious diseases.
 *
 * Running from command line is quite flexible. You can specify
 * scenarios and sensitivity tests. The output is a CSV file written
 * to standard output.
 *
 * To compile:
 *
 * g++ -Wall -O3  -std=c++17 hetgen.cc -o hetgen -lpthread
 *
 * All the parameters are in the options variable in the main() function.
 *
 * Examples:
 *
 * You can run it without any options, in which case it's just the
 * default parameters, and it will run one simulation.
 *
 * ./abm_hetgen
 *
 * Here we set some parameters so that no isolation or tracing of the
 * symptomatic takes place:
 *
 * ./abm_hetgen --trace_effective=0.0 --min_isolation=0.0 --max_isolation=0.0
 *
 * Here we run it 100 times, and make it only report the final result for
 * each run. There are no deaths in this run because everyone recovers before
 * hospital.
 *
 * ./abm_hetgen --trace_effective=0.0 --min_isolation=0.0 --max_isolation=0.0 --threads=1  --runs=100 --iterations=500 --report_frequency=600 --recover_before_hospital=1.0
 *
 * To do sensitivity testing we have to "jiggle" our sensitive parameters.
 * An of the options that's of struct Jiggle can be jiggled. To specify that
 * a parameter musst be jiggled you separate the lower and upper bounds with a
 * colon. E.g. (note the min_test, isolation_period and exposed_risk parameters):
 *
 * ./abm_hetgen --trace_effective=0.0 --min_isolation=0.0 --max_isolation=0.0 --threads=1  --jiggles=1000 --runs=6 --iterations=500 --report_frequency=600 --recover_before_hospital=1.0 --prob_test_infectious_s=0.01:0.99 --mean_test=1:8 --min_test=0:2 --isolation_period=6:14 --exposed_risk=0.1:0.9
 *
 * The above runs 6,000 simulations (1,000 jiggles and 6 runs per jiggle).
 *
 * We can also specify multiple scenarios. Use the + sign at the end of a
 * set of options to specify a new scenario. The new scenario inherits all the
 * values of the previous one.
 *
 * Here we have two scenarios, no intervention and maximum intervention (with
 * full tracing and isolation):
 *
 * ./abm_hetgen --trace_effective=0.0 --min_isolation=0.0 --max_isolation=0.0 --threads=1  --jiggles=1000 --runs=6 --iterations=500 --report_frequency=600 --recover_before_hospital=1.0 --prob_test_infectious_s=0.01:0.99 --mean_test=1:8 --min_test=0:2 --isolation_period=6:14 --exposed_risk=0.1:0.9 + --trace_effective=1.0 --min_isolation=1.0 --max_isolation=1.0
 *
 * The above runs 12,000 simulations (1,000 jiggles, 6 runs of 2 scenarios per
 * jiggle).
 *
 * Here we run 3000 jiggles of 10 runs over 4 scenarios:
 *
 * ./abm_hetgen --trace_effective=0.0 --min_isolation=0.0 --max_isolation=0.0 --threads=0  --jiggles=3000 --runs=10 --iterations=500 --report_frequency=600 --recover_before_hospital=1.0 --prob_test_infectious_s=0.01:0.99 --mean_test=1:8 --min_test=0:2 --isolation_period=6:14 --exposed_risk=0.1:0.9 --asymptomatic=0.05:0.95 --infectious_a_risk=0.1:0.9 --infectious_s_risk=0.1:0.9 --k_assort=32:44  --seed=1 + --trace_effective=1.0 --min_isolation=1.0 --max_isolation=1.0 + --trace_effective=0.3 --min_isolation=0.02 --max_isolation=1.0 + --trace_effective=0.0 --min_isolation=0.4 --max_isolation=1.0 >out.csv
 *
 *
 */

#include <atomic>
#include <algorithm>
#include <any>
#include <cassert>
#include <cmath>
#include <cstdbool>
#include <functional>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include <vector>
#include <utility>

#include <boost/asio.hpp>

#define ISOLATED_BEFORE -11
#define INFECTED_BEFORE -1
#define NONE -2

#define NEGATIVE 0
#define POSITIVE 1

#define SUSCEPTIBLE 0
#define EXPOSED 1
#define INFECTIOUS_A 2
#define INFECTIOUS_S 3
#define INFECTIOUS_H 4
#define INFECTIOUS_I 5
#define RECOVERED 6
#define DEAD 7

#define POS_SUSCEPTIBLE 0.001
#define POS_EXPOSED 0.5
#define POS_INFECTIOUS_A 0.9
#define POS_INFECTIOUS_S 0.999
#define POS_INFECTIOUS_H 0.999
#define POS_INFECTIOUS_I 0.999
#define POS_RECOVERED 0.3
#define POS_DEAD 0.99

#define TEST_SUSCEPTIBLE 0.0
#define TEST_EXPOSED 0.0
#define TEST_INFECTIOUS_A 0.0
#define TEST_INFECTIOUS_S 0.3
#define TEST_INFECTIOUS_H 0.0
#define TEST_INFECTIOUS_I 0.0
#define TEST_RECOVERED 0.0
#define TEST_DEAD 0.0

std::atomic<int> global_seed(NONE);
thread_local std::random_device rd;

thread_local std::mt19937 rng((global_seed == NONE) ? rd() : global_seed++);

inline double rand_0_1() {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}

inline int rand_range(int a, int b)
{
    std::uniform_int_distribution<int> dist(a, b);
    return dist(rng);
}

inline double rand_range(double a, double b)
{
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}

int calc_threads(int threads)
{
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) threads = 1;
    }
    return threads;
}

template<typename T>
struct Jiggle {
    T l; // lower
    T u; // upper
    inline T operator()() { return l; }
    void set() {
        if (u != NONE) {
            l = rand_range(l, u);
        }
    }
    void get(Jiggle<T>& from) {
        if (from.u != NONE) {
            l = from.l;
        }
    }
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const Jiggle<T>& jiggle)
{
    if (jiggle.u == NONE) {
        os << jiggle.l;
    } else {
        os << '[' << jiggle.l << ':' << jiggle.u << ']';
    }
    return os;
}

struct Option {
    std::string name;
    std::string description;
    std::any value;
};

std::ostream& operator<<(std::ostream& os, const Option& option)
{
    os << option.name << ": " << option.description;
    os << " (default: ";
    if (option.value.type() == typeid(int *)) {
        os << *std::any_cast<int *>(option.value);
    } else if (option.value.type() == typeid(double *)) {
        os << *std::any_cast<double *>(option.value);
    } else if (option.value.type() == typeid(Jiggle<int> *)) {
        os << *std::any_cast<Jiggle<int> *>(option.value);
    } else if (option.value.type() == typeid(Jiggle<double> *)) {
        os << *std::any_cast<Jiggle<double> *>(option.value);
    } else {
        os << *std::any_cast<std::string *>(option.value);
    }
    os << ')';
    return os;
}


struct Parameters {
    int first_id = 0;
    int scenario = 0;
    int jiggle = 0;
    int run = 0;
    int id = 0;
    int num_jiggles = 1;
    int num_runs_per_jiggle = 1;
    int seed = NONE;
    int threads = 0;
    int num_iterations = 365;
    int num_agents = 10000;
    int stats_frequency = 1;
    int report_frequency = 20;
    int initial_infections = 0;
    bool verbose = true;
    std::vector<double> health = {0.999, 1.0};
    double risk_infection = 0.1;
    std::unordered_map<int, double> risk_infecting = {
        std::pair<int, double>(INFECTIOUS_A, 0.05),
        std::pair<int, double>(INFECTIOUS_S, 0.1),
        std::pair<int, double>(INFECTIOUS_H, 0.1),
        std::pair<int, double>(INFECTIOUS_I, 0.1)
    };

    std::vector<double> risk_positive = {
        POS_SUSCEPTIBLE,
        POS_EXPOSED,
        POS_INFECTIOUS_A,
        POS_INFECTIOUS_S,
        POS_INFECTIOUS_H,
        POS_INFECTIOUS_I,
        POS_RECOVERED,
        POS_DEAD
    };

    std::vector<double> prob_test = {
        TEST_SUSCEPTIBLE,
        TEST_EXPOSED,
        TEST_INFECTIOUS_A,
        TEST_INFECTIOUS_S,
        TEST_INFECTIOUS_H,
        TEST_INFECTIOUS_I,
        TEST_RECOVERED,
        TEST_DEAD
    };

    double min_isolation = 0.0;
    double max_isolation = 1.0;
    int min_before_isolate = 0;

    double trace_effective = 0.9;
    int min_before_trace = 0;
    double recover_before_hospital = 0.8;
    double recover_before_icu = 0.65;
    double infectious_h_risk = 0.15;
    double recover_before_death = 0.3;
    double infectious_i_risk = 0.15;

    Jiggle<double> exposed_risk = {0.25, NONE};
    Jiggle<double> asymptomatic = {0.75, NONE};
    Jiggle<double> infectious_a_risk = {0.5, NONE};
    Jiggle<int> k_assort = {38, NONE};
    Jiggle<int> k_unassort = {0, NONE};
    Jiggle<int> k_tracing = {38, NONE};

    Jiggle<double> pos_test_susceptible = {POS_SUSCEPTIBLE, NONE};
    Jiggle<double> pos_test_exposed = {POS_EXPOSED, NONE};
    Jiggle<double> pos_test_infectious_a = {POS_INFECTIOUS_A, NONE};
    Jiggle<double> pos_test_infectious_s = {POS_INFECTIOUS_S, NONE};
    Jiggle<double> pos_test_infectious_h = {POS_INFECTIOUS_H, NONE};
    Jiggle<double> pos_test_infectious_i = {POS_INFECTIOUS_I, NONE};
    Jiggle<double> pos_test_recovered = {POS_RECOVERED, NONE};
    Jiggle<double> pos_test_dead = {POS_DEAD, NONE};

    Jiggle<double> prob_test_susceptible = {TEST_SUSCEPTIBLE, NONE};
    Jiggle<double> prob_test_exposed = {TEST_EXPOSED, NONE};
    Jiggle<double> prob_test_infectious_a = {TEST_INFECTIOUS_A, NONE};
    Jiggle<double> prob_test_infectious_s = {TEST_INFECTIOUS_S, NONE};
    Jiggle<double> prob_test_infectious_h = {TEST_INFECTIOUS_H, NONE};
    Jiggle<double> prob_test_infectious_i = {TEST_INFECTIOUS_I, NONE};
    Jiggle<double> prob_test_recovered = {TEST_RECOVERED, NONE};
    Jiggle<double> prob_test_dead = {TEST_DEAD, NONE};

    Jiggle<int> mean_test = {2, NONE};
    Jiggle<int> min_test = {2, NONE};
    Jiggle<int> isolation_period = {10, NONE};
    Jiggle<double> infectious_s_risk = {0.2, NONE};
};

void print_help(const std::vector<Option>& options, const char *prog_name,
                const char *description)
{
    std::cout << prog_name;
    if (description)
        std::cout << ": " << description;
    std::cout << std::endl;

    std::cout << "Syntax:" << std::endl;
    std::cout << prog_name << " ";
    for (auto & option: options)
        std::cout << "[-" << option.name << "=<value>] ";
    std::cout << std::endl;

    std::cout << "\tOptions:" << std::endl;
    for (auto &option: options) {

        std::cout << "\t-" << option << std::endl;
    }
}

void process_file_options(const std::string filename,
                          std::vector<std::string>& arguments)
{
    std::ifstream ifs (filename, std::ifstream::in);

    if (ifs.fail()) {
        std::cerr << "File " << filename << " can't be opened." << std::endl;
        exit(1);
    }

    std::string s;
    while (std::getline(ifs, s)) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
        s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), s.end());
        if (s.size() > 0 && s[0] != '#') {
            arguments.push_back(s);
        }
    }
    ifs.close();
}

void process_command_options(int argc, char *argv[],
                             std::vector<std::string>& arguments)
{
    for (int i = 1; i < argc; i++) {
        arguments.push_back(std::string(argv[i]));
    }
}

void process_options(int argc, char *argv[], std::vector<Option>& options,
                     Parameters &p, std::vector<Parameters> &parameter_set,
                     const char *prog_desc = NULL)
{
    std::vector<std::string> arguments;

    process_command_options(argc, argv, arguments);

    size_t i = 0;
    while (i < arguments.size()) {
        std::string arg = arguments[i];
        if (arg == "-h" || arg == "--help" || arg == "help") {
            print_help(options, argv[0], prog_desc);
            exit(0);
        }

        if (arg == "+") {
            parameter_set.push_back(p);
            ++i;
            continue;
        }

        size_t p = arg.find('=');
        if (p == std::string::npos || p == arg.size() - 1) {
            print_help(options, argv[0], prog_desc);
            exit(1);
        }
        int start = 0;
        if (arg[0] == '-') {
            if (arg[1] == '-')
                start = 2;
            else
                start = 1;
        }
        std::string name = arg.substr(start, p-start);
        std::string value_s = arg.substr(p + 1);

        int found = 0;
        for (auto & option: options) {
            if (option.name == name) {
                found = 1;
                if (name == "input") {
                    process_file_options(value_s, arguments);
                    break;
                }
                if (option.value.type() == typeid(double*)) {
                    *(std::any_cast<double*>(option.value)) = std::stod(value_s);
                } else if (option.value.type() == typeid(int *)) {
                    *(std::any_cast<int*>(option.value))  = std::stoi(value_s);
                } else if (option.value.type() == typeid(Jiggle<double> *)) {
                    size_t p = value_s.find(':');
                    if (p == std::string::npos) {
                        *(std::any_cast<Jiggle<double>*>(option.value))  =
                            Jiggle<double>{std::stod(value_s), NONE};
                    } else {
                        std::string l = value_s.substr(0, p);
                        std::string u = value_s.substr(p + 1);
                        *(std::any_cast<Jiggle<double>*>(option.value)) =
                            Jiggle<double>{std::stod(l), std::stod(u)};
                    }
                } else if (option.value.type() == typeid(Jiggle<int> *)) {
                    size_t p = value_s.find(':');
                    if (p == std::string::npos) {
                        *(std::any_cast<Jiggle<int>*>(option.value)) =
                            Jiggle<int>{std::stoi(value_s), NONE};
                    } else {
                        std::string l = value_s.substr(0, p);
                        std::string u = value_s.substr(p + 1);
                        *(std::any_cast<Jiggle<int>*>(option.value)) =
                            Jiggle<int>{std::stoi(l), std::stoi(u)};
                    }
                }
                break;
            }
        }
        if (!found) {
            std::cerr << "Unknown option:" << arg << std::endl;
            std::cerr << "Try: " << std::endl;
            std::cerr << '\t' << argv[0] << " --help" << std::endl;
            exit(1);
        }
        ++i;
    }
    parameter_set.push_back(p);
}

struct Agent {
    int id_;
    double risk_infection_;
    std::unordered_map<int, double> risk_infecting_;
    std::pair<int, int> infector_ = std::pair<int, int>(NONE, NONE);
    std::vector< std::pair<int, int> > infected_by_me_;
    std::unordered_map< int, std::vector<int> > health_change_iters_;
    int tested_ = 0;
    int test_result_ = NEGATIVE;
    int test_res_iter_ = NONE;
    int isolation_iter_ = NONE;
    double isolated_ = 0.0;
    int health_;
    bool asymptomatic_ = true;
    bool recover_before_hospital_ = true;
    bool recover_before_icu_ = true;
    bool recover_before_death_ = true;

    Agent(Parameters & p, int id) {
        id_ = id;
        {
            std::exponential_distribution<double> dist (1.0 / p.risk_infection);
            risk_infection_ = dist(rng);
        }
        for (auto &m: p.risk_infecting) {
            std::exponential_distribution<double> dist(1.0 / std::get<1>(m));
            risk_infecting_[std::get<0>(m)] = dist(rng);
        }
        if (p.initial_infections == 0) {
            int stage = SUSCEPTIBLE;
            for (auto d: p.health) {
                if (rand_0_1() < d) {
                    health_ = stage;
                    break;
                }
                ++stage;
            }
        } else {
            health_ = SUSCEPTIBLE;
        }
        asymptomatic_ = (rand_0_1() < p.asymptomatic()) ? true : false;
        if (asymptomatic_ == false) {
            recover_before_hospital_ = rand_0_1() < p.recover_before_hospital;
            if (recover_before_hospital_ == false) {
                recover_before_icu_ = rand_0_1() < p.recover_before_icu;
                if (recover_before_icu_ == false) {
                    recover_before_death_ = rand_0_1() < p.recover_before_death;
                }
            }
        }
    }
};


struct Simulation {
    Parameters parameters_;
    std::vector<Agent *> agents_;
    int iteration_;
    int num_agents_isolated_ = 0;
    int num_isolated_ = 0;
    int num_deisolated_ = 0;
    int num_traced_ = 0;
    int total_infected_ = 0;
    int num_agents_tested_ = 0;
    int num_tests_ = 0;
    int num_positives_ = 0;
    int peak_ = 0;
    int peak_total_infections_ = 0;
    int peak_iter_ = 0;
    std::vector< std::unordered_map<std::string, double> > results_;

    void infect(Agent *from, Agent *to) {
        int from_id = INFECTED_BEFORE;
        if (from) {
            from->infected_by_me_.push_back(std::pair<int, int>
                                           (to->id_, iteration_));
            from_id = from->id_;
        }
        to->infector_ = std::pair<int, int> (from_id, iteration_);
        ++total_infected_;
    }

    void isolate(Agent *a) {
        if (a->isolation_iter_ == NONE) ++num_agents_isolated_;
        a->isolation_iter_ = iteration_ + parameters_.isolation_period();
        a->isolated_ = rand_0_1() *
            (parameters_.max_isolation - parameters_.min_isolation) +
            parameters_.min_isolation;
        ++num_isolated_;
    }

    void deisolate(Agent *a) {
        a->isolation_iter_ = - (a->isolation_iter_ + ISOLATED_BEFORE);
        a->isolated_ = 0.0;
        ++num_deisolated_;
    }

    Simulation(const Parameters &p) : parameters_(p) {};

    void init_agents(int n) {
        for (int i = 0; i < n; i++) {
            Agent *a = new Agent(parameters_, i);
            if (a->health_ > SUSCEPTIBLE && a->health_ < DEAD)
                infect(NULL, a);
            agents_.push_back(a);
        }
        if (parameters_.initial_infections) {
            std::vector<Agent *> indices = agents_;
            std::shuffle(indices.begin(), indices.end(), rng);
            for (auto a = indices.begin();
                 a < indices.begin() + parameters_.initial_infections; a++) {
                (*a)->health_ = EXPOSED;
                infect(NULL, *a);
            }
        }
    }

    void stats(int forced = false) {
        if (forced || iteration_ % parameters_.stats_frequency == 0) {
            int infections = 0;
            for (auto & a: agents_) {
                infections += (a->health_ > SUSCEPTIBLE &&
                               a->health_ < RECOVERED);
            }
            if (peak_ < infections) {
                peak_ = infections;
                peak_total_infections_ = total_infected_;
                peak_iter_ = iteration_;
            }
        }
    }

    void report(int forced = false) {
        static std::mutex io_mutex;

        if (forced || iteration_ % parameters_.report_frequency == 0) {
            int susceptible = 0;
            int exposed = 0, infectious_a = 0, infectious_s = 0;
            int infectious_h = 0, infectious_i = 0;
            int recovered = 0, dead = 0, active = 0;
            for (auto & a: agents_) {
                switch(a->health_) {
                case SUSCEPTIBLE: ++susceptible; break;
                case EXPOSED: ++exposed; break;
                case INFECTIOUS_A: ++infectious_a; break;
                case INFECTIOUS_S: ++infectious_s; break;
                case INFECTIOUS_H: ++infectious_h; break;
                case INFECTIOUS_I: ++infectious_i; break;
                case RECOVERED: ++recovered; break;
                case DEAD: ++dead; break;
                }
            }
            active = exposed + infectious_a + infectious_s + infectious_h +
                infectious_i;
            {
                std::lock_guard<std::mutex> lk(io_mutex);
                std::cout << parameters_.id << ','
                          << parameters_.scenario << ','
                          << parameters_.jiggle << ','
                          << parameters_.run << ','
                          << iteration_ << ','
                          << susceptible << ','
                          << exposed << ','
                          << infectious_a << ','
                          << infectious_s << ','
                          << infectious_h << ','
                          << infectious_i << ','
                          << recovered << ','
                          << dead << ','
                          << active << ','
                          << total_infected_ << ','
                          << num_agents_isolated_ << ','
                          << num_isolated_ << ','
                          << num_deisolated_ << ','
                          << num_traced_ << ','
                          << num_agents_tested_ << ','
                          << num_tests_ << ','
                          << num_positives_ << ','
                          << parameters_.k_assort() << ','
                          << parameters_.prob_test_infectious_s() << ','
                          << parameters_.mean_test() << ','
                          << parameters_.min_test() << ','
                          << parameters_.isolation_period() << ','
                          << parameters_.exposed_risk() << ','
                          << parameters_.asymptomatic() << ','
                          << parameters_.infectious_a_risk() << ','
                          << parameters_.infectious_s_risk() << ','
                          << parameters_.min_isolation << ','
                          << parameters_.max_isolation << ','
                          << parameters_.trace_effective
                          << std::endl;
            }
        }
    }

    /*
     * This infection event is the preferred one. Agents can only infect up to k
     * neighbours.
     */

    void event_infect_assort() {
        int neighbors = std::round( (double) parameters_.k_assort() / 2.0);
        std::vector<int> infected(agents_.size(), NONE);
        std::vector<int> indices(agents_.size());
        for (size_t i = 0; i < agents_.size(); i++) indices[i] = i;
        std::shuffle(indices.begin(), indices.end(), rng);

        for (int i = 0; (size_t) i < indices.size(); i++) {
            Agent *a = agents_[i];
            if (a->health_ > EXPOSED && a->health_ < RECOVERED) {
                int from = std::max(0, i - neighbors);
                int to = std::min(i + 1 + neighbors, (int) agents_.size());
                for (int j = from; j < to; j++) {
                    Agent *b = agents_[j];
                    if (b->health_ > SUSCEPTIBLE) continue;
                    double risk;
                    risk = std::min(1.0 - a->isolated_, 1.0 - b->isolated_) *
                        ((a->risk_infecting_[a->health_] + b->risk_infection_) /
                         2.0);
                    if (rand_0_1() < risk) {
                        infected[j] = i;
                    }
                }
            }
        }
        for (int i = 0; (size_t) i < infected.size(); i++) {
            if (infected[i] > NONE) {
                assert(agents_[infected[i]]->health_ > EXPOSED &&
                       agents_[infected[i]]->health_ < RECOVERED);
                assert(agents_[i]->health_ == SUSCEPTIBLE);
                assert(std::abs(i - infected[i]) <=
                       std::round((double)parameters_.k_assort() / 2));
                agents_[i]->health_ = EXPOSED;
                infect(agents_[infected[i]], agents_[i]);
            }
        }
    }

    /*
     * This event randomly infects agents.
     */
    void event_infect_unassort() {
        int k = std::min(parameters_.k_unassort(), (int) agents_.size());
        std::vector<int> infected(agents_.size(), NONE);
        std::vector<Agent*> agents = agents_;
        std::shuffle(agents.begin(), agents.end(), rng);
        for (size_t i = 0; i < agents.size(); i++ ) {
            auto a = agents[i];
            if (a->health_ > EXPOSED && a->health_ < RECOVERED) {
                std::vector<size_t> partners(k);
                for (auto &i: partners) i = rand_range(0, agents_.size() - 1);
                for (auto j: partners) {
                    auto b = agents_[j];
                    if (b->health_ == SUSCEPTIBLE && infected[j] == NONE) {
                        double risk;
                        risk = std::min(1.0 - a->isolated_, 1.0 - b->isolated_) *
                            ((a->risk_infecting_[a->health_] + b->risk_infection_) /
                             2.0);
                        if (rand_0_1() < risk) infected[j] = i;
                    }
                }
            }
        }
        for (size_t i = 0; i < infected.size(); i++) {
            if (infected[i] > NONE) {
                agents_[i]->health_ = EXPOSED;
                infect(agents_[infected[i]], agents_[i]);
            }
        }
    }

    void event_test() {
        for (auto a: agents_) {
            if (a->test_res_iter_ == NONE &&
                a->test_result_ == NEGATIVE &&
                rand_0_1() < parameters_.prob_test[a->health_]) {
                std::poisson_distribution<int> dist(parameters_.mean_test());
                int i = dist(rng);
                a->test_res_iter_ = iteration_ +
                    std::max(i, parameters_.min_test());
                if (rand_0_1() < parameters_.risk_positive[a->health_]) {
                    a->test_result_ = POSITIVE;
                    ++num_positives_;
                } else {
                    a->test_result_ = NEGATIVE;
                }
                ++num_tests_;
                if (a->tested_ == 0) ++num_agents_tested_;
                ++a->tested_;
            }
        }
    }

    void event_isolate() {
        if (parameters_.min_before_isolate > 0) {
            int infections = 0;
            for (auto a: agents_) infections += (a->health_ > INFECTIOUS_A);
            if (infections < parameters_.min_before_isolate)
                return;
            else
                parameters_.min_before_isolate = 0;
        }

        for (auto a: agents_) {
            if (a->test_res_iter_ == iteration_ &&
                a->test_result_ == POSITIVE &&
                a->isolated_ == 0.0) {
                isolate(a);
            }
        }
    }

    void event_deisolate() {
        for (auto a: agents_) {
            if (a->isolation_iter_ == iteration_) deisolate(a);
        }
    }

    void event_trace() {
        if (parameters_.min_before_trace > 0) {
            int infections = 0;
            for (auto a: agents_) infections += (a->health_ > INFECTIOUS_A);
            if (infections < parameters_.min_before_trace) return;
        }
        parameters_.min_before_trace = 0;
        for (auto a: agents_) {
            if (a->test_res_iter_ == iteration_ &&
                a->test_result_ == POSITIVE) {
                int neighbors = std::round((double) parameters_.k_assort() / 2.0);
                int from = std::max(0, a->id_ - neighbors);
                int to = std::min(a->id_ + 1 + neighbors, (int) agents_.size());
                for (int i = from; i < to; i++) {
                    if (i != a->id_ &&
                        agents_[i]->isolated_ == 0.0 &&
                        agents_[i]->health_ < RECOVERED) {
                        if (rand_0_1() < parameters_.trace_effective) {
                            ++num_traced_;
                            isolate(agents_[i]);
                        }
                    }
                }
            }
        }
    }

    // Get test result, i.e. reset test info. This must be done
    // after all other events dependent on test result
    void event_result() {
        for (auto a: agents_) {
            if (a->test_res_iter_ == iteration_) {
                a->test_res_iter_ = NONE;
            }
        }
    }

    void advance_infection(Agent *a, int stage_from, int stage_to, double prob,
                           bool recover) {
        if (a->health_ == stage_from) {
            if (rand_0_1() < prob) {
                if (recover) {
                    a->health_ = RECOVERED;
                    a->health_change_iters_[RECOVERED].push_back(iteration_);
                } else {
                    a->health_ = stage_to;
                    a->health_change_iters_[stage_to].push_back(iteration_);
                }
            }
        }
    }

    void event_exposed() {
        for (auto a: agents_)
            advance_infection(a, EXPOSED, INFECTIOUS_A,
                              parameters_.exposed_risk(), false);
    }

    void event_infectious_a() {
        for (auto a: agents_)
            advance_infection(a, INFECTIOUS_A, INFECTIOUS_S,
                              parameters_.infectious_a_risk(), a->asymptomatic_);
    }

    void event_infectious_s() {
        for (auto a: agents_)
            advance_infection(a, INFECTIOUS_S, INFECTIOUS_H,
                              parameters_.infectious_s_risk(),
                              a->recover_before_hospital_);
    }

    void event_infectious_h() {
        for (auto a: agents_)
            advance_infection(a, INFECTIOUS_H, INFECTIOUS_I,
                              parameters_.infectious_h_risk,
                              a->recover_before_icu_);
    }

    void event_infectious_i() {
        for (auto a: agents_)
            advance_infection(a, INFECTIOUS_I, DEAD,
                              parameters_.infectious_i_risk,
                              a->recover_before_death_);
    }

    void iterate() {
        int n = iteration_ + parameters_.num_iterations;
        for (; iteration_ < n; iteration_++){
            stats();
            report();
            if (parameters_.k_assort() > 0) event_infect_assort();
            if (parameters_.k_unassort() > 0) event_infect_unassort();
            event_test();
            if (parameters_.max_isolation > 0.0) event_isolate();
            if (parameters_.max_isolation > 0.0) event_deisolate();
            if (parameters_.trace_effective > 0.0) event_trace();
            event_result();
            // Stage advances
            event_exposed();
            event_infectious_a();
            event_infectious_s();
            event_infectious_h();
            event_infectious_i();
        }
        stats(true);
        report(true);
    }

    void simulate() {
        init_agents(parameters_.num_agents);
        iteration_ = 0;
        iterate();
    }

    ~Simulation() {
        for (size_t i = 0; i < agents_.size(); i++) delete agents_[i];
    }
};

void write_csv_header() {
    std::cout << "id,scenario,jiggle,run,"
              << "iteration,susceptible,exposed,asymptomatic,symptomatic,"
              << "hospital,icu,recover,dead,active,total_infected,"
              << "agents_isolated,isolated,deisolated,traced,agents_tested,tested,"
              << "positives,k,test_infectious_s,mean_test,"
              << "min_test,isolation_period,exposed_risk,asymp_prob,inf_a_risk,"
              << "inf_s_risk,min_isolation,max_isolation,trace_effective"
              << std::endl;
}

void run_one_sim(Simulation *s)
{
    if (s->parameters_.seed == NONE) {
        rng.seed(time(NULL));
    } else {
        rng.seed(s->parameters_.seed + s->parameters_.id);
    }
    s->simulate();
}

void run_thread(Parameters *p)
{
    static int count = 0;
    static int max_threads = calc_threads(p->threads);
    static std::vector<Simulation *> simulations;
    static std::vector<std::thread> threads;

    if (max_threads == 1) {
        if (p) {
            Simulation s(*p);
            run_one_sim(&s);
        }
    } else {
        if (p == NULL) {
            for (int i = 0; i < count; i++) {
                threads[i].join();
                delete simulations[i];
            }
        } else {
            Simulation *s = new Simulation(*p);
            simulations.push_back(s);
            threads.push_back(std::thread(run_one_sim, s));
            ++count;
            if (count == max_threads) {
                for (int i = 0; i < max_threads; i++) {
                    threads[i].join();
                    delete simulations[i];
                }
                count = 0;
                threads.clear();
                simulations.clear();
            }
        }
    }
}

void set_jiggles(Parameters &p)
{
    p.exposed_risk.set();
    p.asymptomatic.set();
    p.infectious_a_risk.set();
    p.k_assort.set();
    p.prob_test_susceptible.set();
    p.prob_test_exposed.set();
    p.prob_test_infectious_a.set();
    p.prob_test_infectious_s.set();
    p.prob_test_infectious_h.set();
    p.prob_test_infectious_i.set();
    p.prob_test_recovered.set();
    p.prob_test_dead.set();
    p.mean_test.set();
    p.min_test.set();
    p.isolation_period.set();
    p.infectious_s_risk.set();
}

void get_jiggles(Parameters &from, Parameters &to)
{
    to.exposed_risk.get(from.exposed_risk);
    to.asymptomatic.get(from.asymptomatic);
    to.infectious_a_risk.get(from.infectious_a_risk);
    to.k_assort.get(from.k_assort);
    to.prob_test_susceptible.get(from.prob_test_susceptible);
    to.prob_test[SUSCEPTIBLE] = to.prob_test_susceptible();
    to.prob_test_exposed.get(from.prob_test_exposed);
    to.prob_test[EXPOSED] = to.prob_test_exposed();
    to.prob_test_infectious_a.get(from.prob_test_infectious_a);
    to.prob_test[INFECTIOUS_A] = to.prob_test_infectious_a();
    to.prob_test_infectious_s.get(from.prob_test_infectious_s);
    to.prob_test[INFECTIOUS_S] = to.prob_test_infectious_s();
    to.prob_test_infectious_s.get(from.prob_test_infectious_s);
    to.prob_test[INFECTIOUS_H] = to.prob_test_infectious_h();
    to.prob_test_infectious_i.get(from.prob_test_infectious_i);
    to.prob_test[INFECTIOUS_I] = to.prob_test_infectious_i();
    to.prob_test_recovered.get(from.prob_test_recovered);
    to.prob_test[RECOVERED] = to.prob_test_recovered();
    to.prob_test_dead.get(from.prob_test_dead);
    to.prob_test[DEAD] = to.prob_test_dead();
    to.mean_test.get(from.mean_test);
    to.min_test.get(from.min_test);
    to.isolation_period.get(from.isolation_period);
    to.infectious_s_risk.get(from.infectious_s_risk);
}


void run_simulations(std::vector<Parameters> &parms)
{
    assert(parms.size());
    boost::asio::thread_pool pool(calc_threads(parms[0].threads));
    int id = parms[0].first_id;
    size_t num_jiggles = parms[0].num_jiggles;
    size_t num_runs_per_jiggle = parms[0].num_runs_per_jiggle;
    Parameters template_parms = parms[0];

    if (template_parms.seed != NONE) global_seed = template_parms.seed;

    for (size_t i = 0; i < num_jiggles; i++) {
        auto p = template_parms;
        set_jiggles(p);
        for (size_t j = 0; j < num_runs_per_jiggle; j++) {
            for (size_t k = 0; k < parms.size(); k++) {
                get_jiggles(p, parms[k]);
                parms[k].jiggle = i;
                parms[k].run = j;
                parms[k].scenario = k;
                parms[k].id = id++;
                boost::asio::post(pool, [p = parms[k]]() {
                    Simulation s(p);
                    s.simulate();
                });
            }
        }
    }
    pool.join();
}

int main(int argc, char *argv[])
{
    Parameters p;
    std::string s;
    std::vector<Parameters> parameter_set;
    std::vector<Option> options = {
        {
            "input",
            "Filename for input parameters",
            &s
        },
        {
            "first_id",
            "Number to give to first id",
            &p.first_id
        },
        {
            "first_id",
            "Number to give to first id",
            &p.first_id
        },
        {
            "scenario",
            "Number of scenario",
            &p.scenario
        },
        {
            "threads",
            "Number of threads to use in multithreaded runs",
            &p.threads
        },
        {
            "seed",
            "Random seed",
            &p.seed
        },
        {
            "jiggles",
            "Number of sensitivity tests",
            &p.num_jiggles
        },
        {
            "runs",
            "Number of runs per sensitivity test",
            &p.num_runs_per_jiggle
        },
        {
            "agents",
            "Number of agents in the simulation",
            &p.num_agents
        },
        {
            "iterations",
            "Number of iterations per simulation",
            &p.num_iterations
        },
        {
            "report_frequency",
            "How often (in iterations) to output the results",
            &p.report_frequency
        },
        {
            "stats_frequency",
            "How often (in iterations) to record stats for calculating R0",
            &p.stats_frequency
        },
        {
            "initial_infections",
            "If this (n) > 0, overrides health array and sets n agents to exposed",
            &p.initial_infections
        },
        {
            "k_assort",
            "Number of neighbours in assorted infection & trace matching",
            &p.k_assort
        },
        {
            "k_unassort",
            "Number of neighbours in unassorted (random) infection",
            &p.k_unassort
        },
        {
            "trace_effective",
            "Effectiveness of tracing",
            &p.trace_effective
        },
        {
            "min_before_trace",
            "Minimum number of symptomatic infections before tracing begins",
            &p.min_before_trace
        },
        {
            "min_before_isolation",
            "Minimum number of symptomatic infections before isolation begins",
            &p.min_before_isolate
        },
        {
            "min_isolation",
            "Minimum of isolation range",
            &p.min_isolation
        },
        {
            "max_isolation",
            "Maximum of isolation range",
            &p.max_isolation
        },
        {
            "isolation_period",
            "Length of isolation period",
            &p.isolation_period
        },
        {
            "recover_before_hospital",
            "Likelihood of recovering before hospital",
            &p.recover_before_hospital
        },
        {
            "mean_test",
            "Mean number of iterations for test result",
            &p.mean_test
        },
        {
            "min_test",
            "Minimum number of iterations for test result",
            &p.min_test
        },
        {
            "pos_test_susceptible",
            "Likelihood of susceptible agent getting tested",
            &p.pos_test_susceptible
        },
        {
            "pos_test_exposed",
            "Likelihood of exposed agent testing positive",
            &p.pos_test_exposed
        },
        {
            "pos_test_infectious_a",
            "Likelihood of infectious asymptomatic testing positive",
            &p.pos_test_infectious_a
        },
        {
            "pos_test_infectious_s",
            "Likelihood of infectious symptomatic testing positive",
            &p.pos_test_infectious_s
        },
        {
            "pos_test_infectious_h",
            "Likelihood of infectious hospitalised testing positive",
            &p.pos_test_infectious_h
        },
        {
            "pos_test_infectious_i",
            "Likelihood of infectious ICU testing positive",
            &p.pos_test_infectious_i
        },
        {
            "pos_test_recovered",
            "Likelihood of recovered testing positive",
            &p.pos_test_recovered
        },
        {
            "pos_test_dead",
            "Likelihood of dead testing positive",
            &p.pos_test_dead
        },
        {
            "prob_test_susceptible",
            "Likelihood of susceptible agent getting tested",
            &p.prob_test_susceptible
        },
        {
            "prob_test_exposed",
            "Likelihood of exposed agent getting tested",
            &p.prob_test_exposed
        },
        {
            "prob_test_infectious_a",
            "Likelihood of infectious asymptomatic getting tested",
            &p.prob_test_infectious_a
        },
        {
            "prob_test_infectious_s",
            "Likelihood of infectious symptomatic getting tested",
            &p.prob_test_infectious_s
        },
        {
            "prob_test_infectious_h",
            "Likelihood of infectious hospitalised getting tested",
            &p.prob_test_infectious_h
        },
        {
            "prob_test_infectious_i",
            "Likelihood of infectious ICU getting tested",
            &p.prob_test_infectious_i
        },
        {
            "prob_test_recovered",
            "Likelihood of recovered getting tested",
            &p.prob_test_recovered
        },
        {
            "prob_test_dead",
            "Likelihood of dead getting tested",
            &p.prob_test_dead
            },
        {
            "exposed_risk",
            "Probability per day of exposed changing to infectious",
            &p.exposed_risk
        },
        {
            "asymptomatic",
            "Probability of asymptomatic going straight to recovered",
            &p.asymptomatic
            },
        {
            "infectious_a_risk",
            "Probability per day of infectious_a changing health",
            &p.infectious_a_risk
            },
        {
            "infectious_s_risk",
            "Probability per day of infectious_s changing health",
            &p.infectious_s_risk
        }
    };



    process_options(argc, argv, options, p, parameter_set,
                    "Simulate infectious diseases");
    write_csv_header();
    run_simulations(parameter_set);
}
