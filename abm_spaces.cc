/* To compile optimized:
 *   g++ -O3 -DNDEBUG abm_spaces.cc -lpthread
 *   Or if you're using clang:
 *   clang++ -O3 -DNDEBUG c19p5.cc -lpthread
 *
 * To compile for debugging:
 *   g++ -g c19p5.cc -lpthread
 *
 * Examples:
 * To run the default version:
 * ./a.out
 * To run it with tracing on 100 times:
 * ./a.out  tracing=1 num_runs=100
 * To run it with tracing on 10 times and test turnaround set to 8 days:
 * ./a.out  tracing=1 num_runs=10 tat=8
 * Command line arguments:
 * seed: Random seed (uses time otherwise)
 * cores: Number of CPU cores to use (defaults to number on machine)
 * num_runs: Number of simulations to run (default 1)
 * num_agents: Number of agents (default 10000)
 * num_time: Number of days (default 270)
 * tracing: 0 = no tracing, 1 = tracing
 * test_contacts: 0 = no testing of contacts, 1 = testing of contacts
 * tat: Test turnaround time (default 2)
 * initial_infection_rate: Proportion of initially exposed agents (default 0.002)
 * report_frequency: How frequently to report (default 270)
 */

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <ctime>
#include <future>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#define NOT_INFECTED(a) ((a.exposed + a.infectious + a.recovered + a.dead) == 0)

struct Option {
    std::string name;
    std::string description;
    double* value_d;
    int* value_i;
};

struct Parameters {
    int cores = 0;
    int num_runs = 1;
    int num_agents = 10000;
    int num_time = 270;
    int seed = 0;
    int report_frequency = num_time;
    int tracing = 0;
    int test_contacts = 1;
    int tat = 2;
    int output_parameters = 0;
    int sensitivity_tests =  0;
    int uniform_age = 0;
    double asymptomatic_rate = 0.75;
    double icu_rate = 0.1;
    double icu_capacity = 0.001;
    double death_rate = 0.006;
    double stay_home_rate = 0.25;
    double tracing_efficacy = 0.8;
    double isolation_rate = 0.8;
    double house_risk = 0.008;
    double block_risk = 0.0003;
    double class_risk = 0.02;
    double taxi_risk = 0.02;
    double work_risk = 0.036;
    double initial_infection_rate = 0.002;
    double sensitivity_range = 0.0;
};


static thread_local std::mt19937 rng;

inline int rand_range_int(int a, int b)
{
    std::uniform_int_distribution<int> dist(a, b);
    return dist(rng);
}

inline double rand_0_1()
{
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}

inline double rand_range(double a, double b)
{
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}

int age_giver(int uniform_age) {
    if (uniform_age) {
        return rand_range_int(0,90);
    } else {
        std::discrete_distribution<int> distribution  {1027, 1016, 910, 820, 869,
            951, 926, 759, 597, 501, 420,
            357, 288, 217, 148, 96, 98};
        int x = distribution(rng);
        int age = (x * 5) + rand_range_int(0, 4);
        return age;
    }
}

struct Agent {
    int id_;
    int gender;
    int exposed = 0;
    int infectious = 0;
    int isolated = 0;
    int tested_positive = 0;
    int times_tested = 0;
    int days_infectious = 0;
    int days_exposed = 0;
    int time_for_infectiousness = -1;
    int time_for_symptoms = -1;
    int time_for_test = -20;
    int time_for_stay_home = -1;
    int time_for_icu = -1;
    int time_for_cure = -1;
    int time_for_end_isolation = -1;
    int time_for_end_stay_home = -1;
    int time_for_death = -1;
    int will_need_icu = 0;
    int will_die = 0;
    int will_stay_home = 0;
    int will_be_symptomatic = 0;
    int stay_home = 0;
    int age;
    int household = 0;
    int block = 0;
    int class_ = -1;
    int workplace = 0;
    int recovered = 0;
    int taxi = 0;
    int symptomatic = 0;
    int icu = 0;
    int dead = 0;
    int dead_icu_full = 0;
    int times_quarantined = 0;
    int where_infected = 0; // 1=home, 2=block, 3=class, 4=work, 5=taxi
    int birth_week;

    void init_agent(int id, double infection_rate, int uniform_age) {
        id_ = id;
        gender = rand_range_int(0,1);
        age = age_giver(uniform_age);
        birth_week = rand_range_int(1, 52);
        if (rand_0_1() < infection_rate)
            exposed = 1;
    }
};

struct Simulation {

    std::vector<Agent> agents;
    int current_run;
    int household_size = 4;
    int households_per_block = 40;
    int num_schools = 1;
    unsigned class_size = 40;
    double employment_rate = 0.65;
    double regular_taxi_takers = 0.8;
    unsigned num_workplaces;
    unsigned taxi_capacity = 12;
    int icu_capacity;
    int tat;

    int day_first_infectious = 4;
    int v_day_first_infectious = 2;
    int days_infectious_asymptomatic = 2;
    int v_days_infectious_asymptomatic = 1;
    int symptoms_to_end_infectious = 5;
    int v_symptoms_to_end_infectious = 1;
    int symptoms_to_test = 4;
    int v_symptoms_to_test = 3;
    int symptoms_to_icu = 5;
    int v_symptoms_to_icu = 1;
    int time_step = 0;
    int num_in_icu = 0;

    double asymptomatic_rate;
    double icu_rate;
    double death_rate;
    double stay_home_rate;
    double tracing_efficacy;
    double isolation_rate;
    double house_risk;
    double block_risk;
    double class_risk;
    double taxi_risk;
    double work_risk;
    int tracing;
    int test_contacts;

    int peak = 0;
    int peak_time = 0;
    int total_infections_at_peak;

    std::vector<double> death_risk_array;
    std::vector< std::vector<int> > household_indices;
    std::vector< std::vector<int> > block_indices;
    std::vector< std::vector<int> > class_indices;
    std::vector< std::vector<int> > work_indices;
    std::vector< int > taxi_indices;

    std::ostringstream out;

    void death_adjust() {
        double dr0s = 0;
        double dr10s = 0.00695 / 100;
        double dr20s = 0.0309 / 100;
        double dr30s = 0.0844 / 100;
        double dr40s = 0.161 / 100;
        double dr50s = 0.595 / 100;
        double dr60s = 1.93 / 100;
        double dr70s = 4.28 / 100;
        double dr80s = 7.8 / 100;

        double total_risk = (dr0s + dr10s + dr20s + dr30s + dr40s +
                             dr50s + dr60s + dr70s + dr80s) / 9;
        dr0s = dr0s * (death_rate / total_risk);
        dr10s = dr10s * (death_rate / total_risk);
        if (dr10s > 1) dr10s = 1;
        dr20s = dr20s * (death_rate / total_risk);
        if (dr20s > 1) dr20s = 1;
        dr30s = dr30s * (death_rate / total_risk);
        if (dr30s > 1) dr30s = 1;
        dr40s = dr40s * (death_rate / total_risk);
        if (dr40s > 1) dr40s = 1;
        dr50s = dr50s * (death_rate / total_risk);
        if (dr50s > 1) dr50s = 1;
        dr60s = dr60s * (death_rate / total_risk);
        if (dr60s > 1) dr60s = 1;
        dr70s = dr70s * (death_rate / total_risk);
        if (dr70s > 1) dr70s = 1;
        dr80s = dr80s * (death_rate / total_risk);
        if (dr80s > 1) dr80s = 1;
        death_risk_array = {dr10s, dr20s, dr30s, dr40s,
                            dr50s, dr60s, dr70s, dr80s};
    }


    int count_exposures() {
        int total = 0;
        for (auto & agent: agents)
            total += agent.exposed;
        return total;
    }

    void make_agents(int num_agents, double infection_rate, int uniform_age)  {
        for (int i = 0; i < num_agents; i++) {
            Agent a;
            a.init_agent(i, infection_rate, uniform_age);
            agents.push_back(a);
        }
        assert(count_exposures() > 0);
    }

    void make_households()  {
        int x = 0;
        int house_num = 0;
        int len = agents.size();
        while (x < len) {
            std::vector<int> indices;
            int size_next_house = rand_range_int(2, household_size * 2 + 2);
            if ( len - x <= (household_size * 2) + 2) {
                size_next_house = agents.size() - x;
            }
            for (int j = x;  j < x + size_next_house; j++) {
                agents[j].household = house_num;
                indices.push_back(j);
            }
            assert(indices.size() > 0);
            household_indices.push_back(indices);
            ++house_num;
            x += size_next_house;
        }
    }

    int check_blocks() {
        int i = 0;
        for (auto & block : block_indices) {
            if (block.size() < 2)
                return 0;
            ++i;
        }
        return 1;
    }

    void make_blocks()  {
        int num_blocks = std::round(agents.size() / households_per_block);
        for (int i = 0; i < num_blocks; i++) {
            std::vector<int> indices;
            block_indices.push_back(indices);
        }
        for (auto & agent: agents) {
            std::vector<int> indices;
            agent.block = agent.household % num_blocks;
            block_indices[agent.block].push_back(agent.id_);
        }
        if (block_indices.back().size() == 0)
            block_indices.pop_back();
        assert(check_blocks());
    }

    int check_classes() {
        for (auto & room: class_indices) {
            if(room.size() < 2) {
                std::cerr << "Too small class" << std::endl;
                return 0;
            }
        }
        for (auto & room : class_indices) {
            int r = agents[room[0]].class_;
            int a = agents[room[0]].age;
            for (auto i: room) {
                if (agents[i].class_ != r) {
                    std::cerr << "Classroom mismatch: "
                              << i << " "
                              << agents[i].class_ << " "
                              << r << std::endl;
                    return 0;
                }
                if (agents[i].age != a) {
                    std::cerr << "Age mismatch: "
                              << i << " "
                              << agents[i].age << " "
                              << a << std::endl;
                    return 0;
                }
            }
        }
        return 1;
    }

    void make_classes()  {
        std::vector<int> indices[12];

        for (auto& agent: agents) {
            if (agent.age >=6 && agent.age < 18) {
                indices[agent.age - 6].push_back(agent.id_);
            }
        }

        for (auto& indice: indices)
            std::shuffle(indice.begin(), indice.end(), rng);

        int room = -1;
        for (unsigned i = 0; i < 12; i++) {
            for (size_t j = 0; j < indices[i].size(); j++) {
                if (j % class_size == 0) {
                    if (indices[i].size() - j > class_size / 2) {
                        std::vector<int> v;
                        class_indices.push_back(v);
                        room++;
                    }
                }
                assert(agents[indices[i][j]].age >= 6 &&
                       agents[indices[i][j]].age < 18);
                agents[indices[i][j]].class_ = room;
                class_indices.back().push_back(indices[i][j]);
            }
        }

        assert(check_classes());
    }

    int check_workplaces() {
        int w = 1;
        for (auto & work: work_indices) {
            if (work.size() < 1) {
                std::cerr << "Wrong workplace size: " << w << " "
                          << work.size() << std::endl;
                return 0;
            }
            for (auto i : work) {
                if (agents[i].workplace != w) {
                    std::cerr << "Wrong workplace: " << agents[i].workplace
                              << " " << w << " " << work_indices.size()
                              << std::endl;
                    return 0;
                }
            }
            ++w;
        }
        return 1;
    }

    void make_workplaces() {
        std::vector<int> indices;

        for (auto& agent: agents) {
            if (agent.age >= 18 && agent.age < 61) {
                if (rand_0_1() < employment_rate) {
                    indices.push_back(agent.id_);
                }
            }
        }
        std::shuffle(indices.begin(), indices.end(), rng);

        int avg_workplace_size = std::round((double) indices.size() /
                                            num_workplaces);
        int workplace_size = -1;
        int workplace = 0;
        int j = 0;

        for (auto i: indices) {
            if (j > workplace_size) {
                ++workplace;
                std::vector<int> v;
                work_indices.push_back(v);
                j = 0;
                workplace_size = rand_range_int(2, 2 * avg_workplace_size + 2);
            }
            agents[i].workplace = workplace;
            work_indices.back().push_back(i);
            j++;
        }
        assert(check_workplaces());
    }

    void make_taxis() {
        for (auto& agent: agents) {
            if (agent.workplace > 0) {
                if (rand_0_1() < regular_taxi_takers) {
                    agent.taxi = 1;
                    taxi_indices.push_back(agent.id_);
                }
            }
        }
    }

    void house_transmit() {
        int num_infected;
        for (auto& indices: household_indices) {
            std::vector<int> house;
            num_infected = 0;
            for (auto i: indices) {
                if (agents[i].isolated == 0 && agents[i].icu == 0) {
                    house.push_back(i);
                    if (agents[i].infectious)
                        ++num_infected;
                }
            }
            if (num_infected) {
                double this_house_risk = house_risk * num_infected;
                for (auto i : house) {
                    if (NOT_INFECTED(agents[i])) {
                        if (rand_range(0, agents.size()) < this_house_risk) {
                            agents[i].exposed = 1;
                            agents[i].where_infected = 1;
                        }
                    }
                }
            }
        }
    }


    void block_transmit() {
        for (auto& indices: block_indices) {
            std::vector<int> current_block;
            int infections = 0;
            for (auto i: indices) {
                if (agents[i].isolated +
                    agents[i].stay_home +
                    agents[i].icu == 0) {
                    current_block.push_back(i);
                    infections += agents[i].infectious;
                }
            }
            if (current_block.size()) {
                double current_block_risk = block_risk * std::sqrt(infections);
                for (auto i: current_block) {
                    if (NOT_INFECTED(agents[i])) {
                        if (rand_range(0, agents.size()) < current_block_risk) {
                            agents[i].exposed = 1;
                            agents[i].where_infected = 2;
                        }
                    }
                }
            }
        }
    }

    void class_transmit() {
        for (auto& room_indices: class_indices) {
            int infected = 0;
            for (auto i: room_indices)
                infected += agents[i].infectious;
            double this_class_risk = class_risk * infected;
            for (auto i: room_indices) {
                if (NOT_INFECTED(agents[i]) &&
                    agents[i].isolated + agents[i].icu  == 0) {
                    if (rand_range(0, agents.size()) < this_class_risk) {
                        agents[i].exposed = 1;
                        agents[i].where_infected = 3;
                    }
                }
            }
        }
    }

    void work_transmit() {
        for (auto &indices : work_indices) {
            int infected = 0;
            for (auto i: indices) {
                if (agents[i].stay_home + agents[i].isolated +
                    agents[i].dead + agents[i].icu == 0)
                    infected += agents[i].infectious;
            }
            double this_workplace_risk = work_risk * infected;
            for (auto i: indices) {
                if (NOT_INFECTED(agents[i]) &&
                    agents[i].isolated + agents[i].icu == 0) {
                    if (rand_range(0, agents.size()) < this_workplace_risk) {
                        agents[i].exposed = 1;
                        agents[i].where_infected = 4;
                    }
                }
            }
        }
    }

    void taxi_transmit() {
        std::vector<int> indices;
        for (auto i: taxi_indices) {
            if (agents[i].isolated + agents[i].stay_home + agents[i].dead +
                agents[i].icu == 0)
                indices.push_back(i);
        }
        shuffle(indices.begin(), indices.end(), rng);
        unsigned current = 0;
        int num_taxis = (double) taxi_indices.size() / taxi_capacity;
        for (int i = 0; i < num_taxis; i++) {
            int infections = 0;
            size_t j = current;
            std::vector<int> taxi;
            for (;j < current + taxi_capacity && j < indices.size(); j++) {
                taxi.push_back(indices[j]);
                infections += agents[indices[j]].infectious;
            }
            if (infections) {
                double current_taxi_risk = taxi_risk * infections;
                for (auto i: taxi) {
                    if (NOT_INFECTED(agents[i])) {
                        if (rand_range(0, agents.size()) < current_taxi_risk) {
                            agents[i].exposed = 1;
                            agents[i].where_infected = 5;
                        }
                    }
                }
            }
            current += taxi_capacity;
        }
    }

    int get_tested(int days_infected) {
        int test_result = 0;
        int k = rand_range_int(1, 100);
        // Below figures are based on this article in AIM https://www.acpjournals.org/doi/10.7326/M20-1495
        if (days_infected == 1) {
            if (k < 1) test_result = 1;
        } else if (days_infected == 2) {
            if (k < 1) test_result = 1;
        } else if (days_infected == 3) {
            if (k < 6) test_result = 1;
        } else if (days_infected == 4) {
            if (k < 33) test_result = 1;
        } else if (days_infected == 5) {
            if (k < 63) test_result = 1;
        } else if (days_infected == 6) {
            if (k < 76) test_result = 1;
        } else if (days_infected == 7) {
            if (k < 81) test_result = 1;
        } else if (days_infected == 8) {
            if (k < 81) test_result = 1;
        } else if (days_infected == 9) {
            if (k < 80) test_result = 1;
        } else if (days_infected == 10) {
            if (k < 78) test_result = 1;
        } else if (days_infected == 11) {
            if (k < 75) test_result = 1;
        } else if (days_infected == 12) {
            if (k < 71) test_result = 1;
        } else if (days_infected == 13) {
            if (k < 67) test_result = 1;
        } else if (days_infected == 14) {
            if (k < 63) test_result = 1;
        } else if (days_infected == 15) {
            if (k < 59) test_result = 1;
        } else if (days_infected == 16) {
            if (k < 55) test_result = 1;
        } else if (days_infected == 17) {
            if (k < 51) test_result = 1;
        } else if (days_infected == 18) {
            if (k < 47) test_result = 1;
        } else if (days_infected == 19) {
            if (k < 43) test_result = 1;
        } else if (days_infected == 20) {
            if (k < 39) test_result = 1;
        } else if (days_infected == 21) {
            if (k < 35) test_result = 1;
        } else if (days_infected == 22) {
            if (k < 31) test_result = 1;
        } else if (days_infected == 23) {
            if (k < 27) test_result = 1;
        } else if (days_infected == 24) {
            if (k < 19) test_result = 1;
        } else if (days_infected == 25) {
            if (k < 15) test_result = 1;
        } else if (days_infected == 26) {
            if (k < 11) test_result = 1;
        } else if (days_infected == 27) {
            if (k < 7) test_result = 1;
        } else if (days_infected == 28) {
            if (k < 3) test_result = 1;
        }
        return test_result;
    }

    void household_tracing(struct Agent & agent) {
        std::vector<int> household;
        for (auto i: household_indices[agent.household]) {
            if (agents[i].tested_positive == 0 &&
                agents[i].time_for_test < time_step - 14 &&
                agents[i].isolated == 0) {
                assert(agents[i].id_ != agent.id_);
                household.push_back(i);
            }
        }
        shuffle(household.begin(), household.end(), rng);
        int n = std::min(household.size(), (size_t)
                         std::round(tracing_efficacy * household.size()));
        for (int i = 0; i < n; i++) {
            int j = household[i];
            assert(agents[j].household == agent.household);
            if (test_contacts == 1) {agents[j].time_for_test = time_step + tat +
                rand_range_int(1, 2);}
            if (rand_0_1() < isolation_rate) {
                agents[j].isolated = 1;
                agents[j].time_for_end_isolation = time_step + 14;
                ++agents[j].times_quarantined;
            } else {
                agents[j].stay_home = 1;
                agents[j].time_for_end_stay_home = time_step + 14;
            }
        }
    }

    void workplace_tracing(struct Agent & agent) {
        if (agent.workplace > 0) {
            std::vector<int> workplace;
            for (auto i: work_indices[agent.workplace - 1]) {
                if (agents[i].tested_positive == 0 &&
                    agents[i].time_for_test < (time_step - 14) &&
                    agents[i].isolated == 0) {
                    assert(agents[i].id_ != agent.id_);
                    assert(agents[i].workplace == agent.workplace);
                    workplace.push_back(i);
                }
            }
            shuffle(workplace.begin(), workplace.end(), rng);
            int n = std::min(workplace.size(), (size_t)
                             std::round(tracing_efficacy * workplace.size()));
            for (int i = 0; i < n; i++) {
                int j = workplace[i];
                assert(agents[j].workplace == agent.workplace);
                if (test_contacts == 1) {agents[j].time_for_test = time_step + tat +
                    rand_range_int(1, 2);}
                if (rand_0_1() < isolation_rate) {
                    agents[j].isolated = 1;
                    agents[j].time_for_end_isolation = time_step + 14;
                    ++agents[j].times_quarantined;
                } else {
                    agents[j].stay_home = 1;
                    agents[j].time_for_end_stay_home = time_step + 14;
                }
            }
        }
    }

    void class_tracing(struct Agent & agent) {
        if (agent.class_ > -1) {
            std::vector<int> room;
            for (auto i: class_indices[agent.class_]) {
                assert(agents[i].class_ == agent.class_);
                assert(agents[i].age == agent.age);
                if (agents[i].tested_positive + agents[i].isolated == 0 &&
                    agents[i].time_for_test < time_step - 14) {
                    assert(agents[i].id_ != agent.id_);
                    room.push_back(i);
                }
            }
            shuffle(room.begin(), room.end(), rng);
            int n = std::min(room.size(), (size_t)
                             std::round(tracing_efficacy * room.size()));
            assert(n <= (int) room.size());
            for (int i = 0; i < n; i++) {
                int j = room[i];
                assert(agents[j].class_ == agent.class_);
                if (test_contacts == 1) {agents[j].time_for_test = time_step + tat + rand_range_int(1, 2);}
                if (rand_0_1() < isolation_rate) {
                    agents[j].isolated = 1;
                    agents[j].time_for_end_isolation = time_step + 14;
                    ++agents[j].times_quarantined;
                } else {
                    agents[j].stay_home = 1;
                    agents[j].time_for_end_stay_home = time_step + 14;
                }
            }
        }
    }

    void disease_progress() {
        for (auto& agent: agents) {
            if (agent.dead)
                continue;
            if (agent.infectious) {
                ++agent.days_infectious;
            }
            if (agent.time_for_end_stay_home == time_step) {
                agent.stay_home = 0;
            }
            if (agent.time_for_end_isolation == time_step) {
                agent.isolated = 0;
            }
            if (agent.time_for_test == time_step &&
                agent.tested_positive + agent.dead == 0) {
                int days_infected = agent.days_exposed + agent.days_infectious -
                    tat;
                if (agent.days_exposed) {
                    agent.tested_positive = get_tested(days_infected);
                }
                ++agent.times_tested;
                if (tracing && agent.tested_positive) {
                    agent.isolated = 1;
                    agent.time_for_end_isolation = time_step + 21;
                    ++agent.times_quarantined;
                    household_tracing(agent);
                    workplace_tracing(agent);
                    class_tracing(agent);
                }
            }
            if (agent.time_for_stay_home == time_step) {
                agent.stay_home = 1;
            }

            if (agent.time_for_icu == time_step) {
                if (num_in_icu < icu_capacity) {
                    ++num_in_icu;
                    agent.icu = 1;
                } else {
                    agent.infectious = 0;
                    agent.dead = 1;
                    agent.dead_icu_full = 1;
                    continue;
                }
            }
            if (agent.time_for_cure == time_step) {
                agent.infectious = 0;
                agent.recovered = 1;
                if (agent.icu) {
                    --num_in_icu;
                    agent.icu = 0;
                }
                agent.isolated = 0;
                agent.stay_home = 0;
                continue;
            }
            if (agent.time_for_death == time_step) {
                agent.infectious = 0;
                agent.dead = 1;
                agent.icu = 0;
                --num_in_icu;
                agent.isolated = 0;
            }
        }
    }

    int death_test(int age) {
        int index = -1;
        int len = death_risk_array.size();
        // Apply the age-based death risk to the agent
        if (age < 10)
            return 0;
        else {
            index = age / 10 - 1;
            if (index >= len)
                index = len - 1;
            return (rand_0_1() < death_risk_array[index]) ? 1 : 0;
        }
    }

    int check_agent(struct Agent & agent) {
        if (agent.exposed + agent.infectious == 2) {
            std::cout << agent.id_ << " "
                      << agent.exposed << " " << agent.infectious << " "
                      << agent.where_infected << std::endl;
            return 0;
        }
        return 1;
    }

    void exposed_progress() {
        int from_i = day_first_infectious - v_day_first_infectious;
        int to_i = day_first_infectious + v_day_first_infectious;
        int from_s = days_infectious_asymptomatic -
            v_days_infectious_asymptomatic;
        int to_s = days_infectious_asymptomatic + v_days_infectious_asymptomatic;
        int from_t = symptoms_to_test - v_symptoms_to_test;
        int to_t = symptoms_to_test + v_symptoms_to_test;

        for (auto & agent: agents) {
            if (agent.dead)
                continue;

            if (agent.exposed)
                ++agent.days_exposed;

            if (agent.infectious == 0 && agent.days_exposed == 1) {
                agent.time_for_infectiousness = time_step +
                    rand_range_int(from_i, to_i);
                agent.will_be_symptomatic =
                    (rand_0_1() < asymptomatic_rate) ? 0 : 1;
                if (agent.will_be_symptomatic) {
                    agent.will_stay_home =
                        (rand_0_1() < stay_home_rate) ? 1 : 0;
                    agent.will_need_icu =
                        (rand_0_1() < icu_rate) ? 1 : 0;
                    agent.time_for_symptoms = agent.time_for_infectiousness +
                        rand_range_int(from_s, to_s);
                    agent.time_for_test = agent.time_for_symptoms +
                        rand_range_int(from_t, to_t) + tat;
                }
                if (agent.will_stay_home == 1) {
                    agent.time_for_stay_home = agent.time_for_symptoms +
                        rand_range_int(1, 5);
                }
                if (agent.will_need_icu) {
                    agent.time_for_icu = agent.time_for_symptoms +
                        rand_range_int(3, 10);
                    agent.will_die = death_test(agent.age);
                }
                if (agent.will_die == 0) {
                    agent.time_for_cure = time_step + rand_range_int(8, 17);
                    if (agent.time_for_icu > agent.time_for_cure) {
                        agent.time_for_icu = agent.time_for_cure;
                    }
                } else {
                    agent.time_for_death = agent.time_for_icu +
                        rand_range_int(1, 10);
                }
            }
            if (agent.time_for_infectiousness == time_step) {
                agent.exposed = 0;
                agent.infectious = 1;
            }
            assert(check_agent(agent));
        }
    }

    void report() {
        int total_tests = 0;
        int total_positive_tests = 0;
        int exposed = 0;
        int infectious = 0;
        int recovered = 0;
        int dead = 0;
        int dead_icu_full = 0;
        int total_infections = 0;
        int infected_initial = 0;
        int infected_home = 0;
        int infected_school = 0;
        int infected_block = 0;
        int infected_work = 0;
        int infected_taxi = 0;
        int quarantines = 0;
        double R0 = 1.0 / ( (double) (agents.size() - total_infections_at_peak)
                            / agents.size());
        double percent_detected = 0.0;
        for (auto & agent: agents) {
            total_tests += agent.times_tested;
            total_positive_tests += agent.tested_positive;
            dead += agent.dead;
            recovered += agent.recovered;
            infectious += agent.infectious;
            exposed += agent.exposed;
            dead_icu_full += agent.dead_icu_full;
            if (agent.recovered || agent.dead ||
                agent.exposed || agent.infectious)
                ++total_infections;
            infected_initial += (!NOT_INFECTED(agent) &&
                                agent.where_infected == 0);
            infected_home += (agent.where_infected == 1);
            infected_block += (agent.where_infected == 2);
            infected_school += (agent.where_infected == 3);
            infected_work += (agent.where_infected == 4);
            infected_taxi += (agent.where_infected == 5);
            quarantines += agent.times_quarantined;
        }
        if (total_infections) {
            percent_detected =  ( (double) total_positive_tests /
                                  total_infections) * 100.0;
        }
        assert(infected_initial + infected_home + infected_school +
               infected_block + infected_work +
               infected_taxi == total_infections);


        out << current_run
            << "," << time_step
            << "," << total_infections
            << "," << exposed
            << "," << infectious
            << "," << recovered
            << "," << dead
            << "," << dead_icu_full
            << "," << peak
            << "," << peak_time
            << "," << total_infections_at_peak
            << "," << infected_initial
            << "," << infected_home
            << "," << infected_school
            << "," << infected_block
            << "," << infected_work
            << "," << infected_taxi
            << "," << total_tests
            << "," << total_positive_tests
            << "," << percent_detected
            << "," << quarantines
            << "," << R0
            << "," << tracing
            << "," << tat
            << std::endl;
    }

    void run_model(int run, Parameters& parameters) {
        current_run = run;
        tracing = parameters.tracing;
        test_contacts = parameters.test_contacts;
        tat = parameters.tat;

        asymptomatic_rate = parameters.asymptomatic_rate;
        icu_rate = parameters.icu_rate;
        death_rate = parameters.death_rate;
        stay_home_rate = parameters.stay_home_rate;
        tracing_efficacy = parameters.tracing_efficacy;
        isolation_rate = parameters.isolation_rate;

        house_risk = parameters.house_risk * parameters.num_agents;
        block_risk = parameters.block_risk * parameters.num_agents;
        class_risk = parameters.class_risk * parameters.num_agents;
        taxi_risk = parameters.taxi_risk * parameters.num_agents;
        work_risk = parameters.work_risk * parameters.num_agents;

        num_workplaces =
            (double) parameters.num_agents / 28.57142857142857142857;
        icu_capacity = parameters.icu_capacity * parameters.num_agents;

        death_rate = parameters.death_rate *
            (1 / (icu_rate * (1 - asymptomatic_rate)));
        death_adjust();

        make_agents(parameters.num_agents, parameters.initial_infection_rate,
                    parameters.uniform_age);
        make_households();

        make_blocks();
        make_classes();
        make_workplaces();
        make_taxis();

        for (time_step = 0; time_step < parameters.num_time; time_step++) {
            house_transmit();
            block_transmit();
            class_transmit();
            work_transmit();
            taxi_transmit();

            disease_progress();
            exposed_progress();

            int infected = 0, recovered = 0;

            for (auto & agent: agents) {
                if (agent.infectious || agent.exposed)
                    ++infected;
                else if (agent.recovered)
                    ++recovered;
            }
            if (peak < infected) {
                peak = infected;
                peak_time = time_step;
                total_infections_at_peak = infected + recovered;
            }
            if ( (time_step + 1) % parameters.report_frequency == 0)
                report();
        }
    }
};

double peturb(double mean, double range)
{
    double from = mean * (1.0 - range);
    double to = mean * (1.0 + range);

    return rand_range(from, to);
}

void peturb_parameters(Parameters &p)
{
    double range = p.sensitivity_range;
    p.asymptomatic_rate = peturb(p.asymptomatic_rate, range);
    p.icu_rate = peturb(p.icu_rate, range);
    p.icu_capacity = peturb(p.icu_capacity, range);
    p.death_rate = peturb(p.death_rate, range);
    p.stay_home_rate = peturb(p.stay_home_rate, range);
    p.tracing_efficacy = peturb(p.tracing_efficacy, range);
    p.isolation_rate = peturb(p.isolation_rate, range);
    p.house_risk = peturb(p.house_risk, range);
    p.block_risk = peturb(p.block_risk, range);
    p.class_risk = peturb(p.class_risk, range);
    p.taxi_risk = peturb(p.taxi_risk, range);
    p.work_risk = peturb(p.work_risk, range);
    p.initial_infection_rate = peturb(p.initial_infection_rate, range);
}

std::string run_one_simulation(int current_run, Parameters parameters)
{
    Simulation s;
    if (parameters.seed > 0) {
        rng.seed(parameters.seed + 11 * current_run);
    } else {
        rng.seed(time(NULL) + 111 * current_run);
    }
    s.run_model(current_run, parameters);
    return s.out.str();
}

void calc_cores(int & cores)
{
    if (cores == 0) {
        cores = std::thread::hardware_concurrency();
        if (cores == 0)
            cores = 1;
    }
}

void print_csv_header()
{
    std::cout << "Run,"
              << "Time step,"
              << "Infections,"
              << "Exposed,"
              << "Infectious,"
              << "Recovered,"
              << "Dead,"
              << "ICU full,"
              << "Peak,"
              << "Peak time,"
              << "Peak total,"
              << "Initial,"
              << "Home,"
              << "School,"
              << "Block,"
              << "Work,"
              << "Taxi,"
              << "Tests,"
              << "Positive,"
              << "% Detected,"
              << "Quarantines,"
              << "R0,"
              << "Tracing,"
              << "TaT"
              << std::endl;
}

void run_simulations(Parameters & parameters)
{
    calc_cores(parameters.cores);
    int current_run = 0;
    print_csv_header();

    for (int i = 0; i < parameters.num_runs; i += parameters.cores) {
        std::vector<std::future<std::string> > output;
        int j = 0;
        for (; j < parameters.cores && current_run < parameters.num_runs;
             current_run++, j++) {
            Parameters p = parameters;
            if (parameters.sensitivity_range > 0.0)
                peturb_parameters(p);
            output.push_back(std::async(std::launch::async, run_one_simulation,
                                        current_run, p));
        }
        for (auto & s: output) {
            std::cout << s.get();
        }
    }
}

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
    for (auto &option: options)
        std::cout << "\t-" << option.name << ": "
                  << option.description
                  << " (default: " << ( (option.value_d == NULL) ?
                                        *option.value_i : *option.value_d)
                  << ")" << std::endl;
}

void process_options(int argc, char *argv[], std::vector<Option>& options,
                     const char *prog_desc = NULL)
{
    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg == "-h" || arg == "--help" || arg == "help") {
            print_help(options, argv[0], prog_desc);
            exit(0);
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
                if (option.value_d) {
                    *option.value_d = std::stod(value_s);
                } else {
                    *option.value_i = std::stoi(value_s);
                }
            }
        }
        if (!found) {
            std::cerr << "Unknown option:" << arg << std::endl;
            std::cerr << "Try: " << std::endl;
            std::cerr << '\t' << argv[0] << " --help" << std::endl;
            exit(1);
        }
    }
}

int main(int argc, char *argv[])
{
    Parameters parameters;
    std::vector<Option> options = {
        {
            "cores",
            "Number of cores to use in multithreaded runs",
            NULL,
            &parameters.cores
        },
        {
            "seed",
            "Random seed",
            NULL,
            &parameters.seed
        },
        {
            "num_runs",
            "Number of simulations",
            NULL,
            &parameters.num_runs
        },
        {
            "num_agents",
            "Number of agents in the simulation",
            NULL,
            &parameters.num_agents
        },
        {
            "num_time",
            "Number of days in each simulation",
            NULL,
            &parameters.num_time
        },
        {
            "tracing",
            "Whether or not to do do tracing (0=off, 1=on)",
            NULL,
            &parameters.tracing
        },
        {
            "test_contacts",
            "Whether or not to test contacts (0=off, 1=on)",
            NULL,
            &parameters.test_contacts
        },
        {
            "tat",
            "Turnaround time for a test",
            NULL,
            &parameters.tat
        },
        {
            "initial_infection_rate",
            "Proportion of agents exposed at beginning of simulation",
            &parameters.initial_infection_rate,
            NULL
        },
        {
            "house_risk",
            "Risk of house exposure",
            &parameters.house_risk,
            NULL
        },
        {
            "block_risk",
            "Risk of exposure in neighborhood",
            &parameters.block_risk,
            NULL
        },
        {
            "class_risk",
            "Risk of classroom exposure",
            &parameters.class_risk,
            NULL
        },
        {
            "taxi_risk",
            "Risk of exposure in taxi",
            &parameters.taxi_risk,
            NULL
        },
        {
            "work_risk",
            "Risk of exposure at work",
            &parameters.work_risk,
            NULL
        },
        {
            "asymptomatic_rate",
            "Asymptotic infection rate",
            &parameters.asymptomatic_rate,
            NULL
        },
        {
            "icu_rate",
            "Admitted to ICU rate",
            &parameters.icu_rate,
            NULL
        },
        {
            "icu_capacity",
            "Proportion of population that could be admitted to ICU",
            &parameters.icu_capacity,
            NULL
        },
        {
            "death_rate",
            "Infection death rate",
            &parameters.death_rate,
            NULL
        },
        {
            "uniform_age",
            "Set ages using uniform random",
            NULL,
            &parameters.uniform_age
        },
        {
            "stay_home_rate",
            "Rate at which possibly exposed who stay home",
            &parameters.stay_home_rate,
            NULL
        },
        {
            "tracing_efficacy",
            "Proportion of contacts traced",
            &parameters.tracing_efficacy,
            NULL
        },
        {
            "isolation_rate",
            "Proportion of infected who isolate",
            &parameters.isolation_rate,
            NULL
        },
        {
            "sensitivity_range",
            "Modify parameters on each iteration proportionate to this",
            &parameters.sensitivity_range,
            NULL
        },
        {
            "report_frequency",
            "How often (in simulation days) to output the results",
            NULL,
            &parameters.report_frequency
        },
        {
            "output_parameters",
            "Whether to output the parameters",
            NULL,
            &parameters.output_parameters
        }
    };

    process_options(argc, argv, options, "Simulate Covid-19");

    if (parameters.output_parameters) {
        for (auto & option: options) {
            std::cout << option.name << ": "
                      << ( (option.value_d == NULL) ?
                           *option.value_i : *option.value_d)
                      << std::endl;
        }
    }
    run_simulations(parameters);
}
