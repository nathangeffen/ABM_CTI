#!/bin/sh
"""
Example of how to use this:

from abm import Simulation
s = Simulation()
s.simulate()

You can also change parameters, e.g. this would turn off tracing and isolation
and run for 700 iterations (or days).

from abm import Simulation
s = Simulation(iterations=700, trace_effective=0.0, min_isolation=0.0,
               max_isolation=0.0)
s.simulate()

Parameters:

All user modifiable parameters are set in the simulation.__init__ method.

iterations: 365 # Number of iterations (1 iteration by default equal 1 day)
num_agents: 10000 # Number of agents
stats_frequency: 1 # Calc stats ever X iterations (lower = more accurate R0)
report_frequency: 20 # Report results ever X iterations

# There are 8 health stages from Health.SUSCEPTIBLE to Health.DEAD and
# this list is used initialise the proportion of agents in each stage.
# By default 99.9% are initiated to Health.SUSCEPTIBLE and the remainder
# to Health.EXPOSED to get the epidemic going. But you can change the
# proportions and add more stages.
health: [0.999, 1.0]

# This next parameter has a big effect on R. It is the number of adjacent
# agents next to an infectious agent (half to the left, half to the right)
# who come into contact with that infectious agent. Agents close to the
# edges of the array have fewer contacts else there'd be an array out-of-bounds
# error.
k_match: 38
# Daily risk of infection for agent exposed to infectious agent. This is
# an average risk, and is the parameter passed to the exponential random
# distribution when each agent is initialized.
risk_infection: 0.1
# This next parameter is the infectiousness of each infectious health stage.
# As with risk_infection it is the average risk that is drawn from an
# exponential distribution for each agent when it's initialized.
# To calculate the infectiousness of a contact, add this to the
# risk_infection and divide by 2 (and multiple by level of isolation).
risk_infecting = {
   Health.INFECTIOUS_A: 0.05, # Asymptomatic
   Health.INFECTIOUS_S: 0.1,  # Symptomatic
   Health.INFECTIOUS_H: 0.1,  # Hospitalised
   Health.INFECTIOUS_I: 0.1,  # ICU
}

test_likelihood: 0.3 # Probability of test per iteration of infected agent
mean_test: 2 # Mean number of test days to get test result back
min_test: 2 # Minimum number of days for test result to come back
isolation_period: 10 # Number iterations in isolation
min_isolation: 0.8 # Minimum effectiveness of isolation
max_isolation: 1.0 # Maximum effectiveness of isolation
trace_effective: 0.9 # Effectiveness of tracing

# Health stage transition parameters
exposed_risk: 0.25 # Probability of staying Health.EXPOSED per day
infectious_a_risk: 0.5 # Prob of remaining Health.ASYMPTOMATIC per day
asymptomatic: 0.75 # Prob of going direct from Health.ASYMPTOMATIC to recovered
infectious_s_risk: 0.2 # Prob of remaining Health.SYMPTOMATIC per day
recover_before_hospital: 0.8 # Prob of recovering before hospitalisation
infectious_h_risk: 0.15 # Prob of remaining in hospital per day
recover_before_icu: 0.5 # Prob of only recovering before needing icu
infectious_i_risk: 0.15 # Prob of remaining in ICU per day
recover_before_death: 0.6 # Prob of recovering from ICU rather than dying

"""

import argparse
import datetime
from multiprocessing import Pool
import random
import numpy.random as npr
from enum import IntEnum
import copy
import pickle
import csv


class Health(IntEnum):
    SUSCEPTIBLE = 0
    EXPOSED = 1
    INFECTIOUS_A = 2 # Asymptomatic
    INFECTIOUS_S = 3 # Symptomatic
    INFECTIOUS_H = 4 # Hospital
    INFECTIOUS_I = 5 # ICU
    RECOVERED = 6
    DEAD = 7

def _getarg(dct, key, missing):
    if key in dct:
        result = dct[key]
        del(dct[key])
    else:
        result = missing

    return result

class Agent:

    def __init__(self, sim, id):
        self.id = id
        self.risk_infection = random.expovariate(1.0 / sim.risk_infection)
        self.risk_infecting = {}
        for key, value in sim.risk_infecting.items():
            self.risk_infecting[key] = random.expovariate(1.0 / value)
        r = random.random()
        probs = sim.health
        for stage, prob in zip(range(len(probs)), probs):
            if r < prob:
                self.health = stage
                break
        self.infector = None
        self.infected_by_me = []
        self.test_res_iter = None # Iteration that test result comes back
        self.isolation_iter = -1
        self.isolated = 0.0
        self.asymptomatic = True if random.random() < sim.asymptomatic else False
        self.recover_before_hospital = True \
                    if (self.asymptomatic is True or
                        random.random() < sim.recover_before_hospital) else False
        self.recover_before_icu = True \
                    if (self.recover_before_hospital is True or
                        random.random() < sim.recover_before_icu) else False
        self.recover_before_death = True \
                if (self.recover_before_icu is True or
                    random.random() < sim.recover_before_death) else False
        self.inf_end_iter = None # Iteration on which no longer infectious

    def isolate(self, sim):
        self.isolation_iter = sim.iteration + sim.isolation_period
        self.isolated = random.random() * \
                     (sim.max_isolation - sim.min_isolation) + \
                     sim.min_isolation

    def deisolate(self, sim):
        self.isolation_iter = -1
        self.isolated = 0.0


class Simulation:

    def __init__(self, **kwargs):

        dct = copy.deepcopy(kwargs)
        # All user modifiable parameters are set in this method
        # simulation engine
        self.iterations = _getarg(dct, 'iterations', 365)
        self.num_agents = _getarg(dct, 'num_agents', 10000)
        self.stats_frequency = _getarg(dct,'stats_frequency', 1)
        self.report_frequency = _getarg(dct, 'report_frequency', 20)
        self.verbose = _getarg(dct,'verbose', True)
        self.name = _getarg(dct, "name", "")

        # Agent initiation

        ## This initializes 999/1000 agents to susceptible and 1/1000 to
        ## exposed.
        self.health = _getarg(dct, 'health', [0.999, 1.0])

        # Infection event
        ## Max number of neighbours to consider for infections
        ## The number is divided in half and agents to left and right are
        ## considered.
        ## Agents on the edge of the array have fewer contacts.
        self.k_match = _getarg(dct, 'k_match', 38)

        ## Daily risk of getting infected (exponential dist)
        self.risk_infection = _getarg(dct, 'risk_infection', 0.1)

        ## Daily risk of infecting (exponential dist)
        self.risk_infecting = {
            Health.INFECTIOUS_A: 0.05,
            Health.INFECTIOUS_S: 0.1,
            Health.INFECTIOUS_H: 0.1,
            Health.INFECTIOUS_I: 0.1,
        }

        # test event
        ## Likelihood per iteration of getting tested
        self.test_likelihood = _getarg(dct, 'test_likelihood', 0.3)

        ## Mean duration it takes to test (Poisson distribution)
        self.mean_test = _getarg(dct, 'mean_test', 2)

        ## Minimum number of test days
        self.min_test = _getarg(dct, 'min_test', 2)

        # Isolate event

        ## Length of isolation period in days
        self.isolation_period = _getarg(dct, 'isolation_period', 10)

        ## Minimum effectiveness of isolation
        self.min_isolation = _getarg(dct, 'min_isolation', 0.0)

        ## Maximum effectiveness of isolation
        self.max_isolation = _getarg(dct, 'max_isolation', 1.0)

        ## Minimum symptomatic infections before isolating
        self.min_before_isolate = _getarg(dct, 'min_before_isolate', 0)


        # Trace event

        ## Effectiveness of tracing (0=no tracing 1=100% tracing)
        self.trace_effective = _getarg(dct, 'trace_effective', 0.9)
        ## Minimum symptomatic infections before tracing
        self.min_before_trace = _getarg(dct, 'min_before_trace', 0)

        # Exposed event
        ## Prob of remaining exposed per day
        self.exposed_risk = _getarg(dct, 'exposed_risk', 0.25)

        # Infectious asymptomatic event

        ## Prob of only being asymptomatic
        self.asymptomatic = _getarg(dct, 'asymptomatic', 0.75)

        ## Prob of remaining infectious asymptomatic per day
        self.infectious_a_risk = _getarg(dct, 'infectious_a_risk', 0.5)

        # Infectious symptomatic event

        ## Prob of infectious_a only recovering before hospitalisation
        self.recover_before_hospital = _getarg(dct, 'recover_before_hospital',
                                               0.8)

        ## Prob of remaining infectious symptomatic per day
        self.infectious_s_risk = _getarg(dct, 'infectious_s_risk', 0.2)

        # Infectious hospitalisation event

        ## Prob of infectious_h only recovering before needing icu
        self.recover_before_icu = _getarg(dct, 'recover_before_icu', 0.65)

        ## Prob of remaining infectious hospitalised per day
        self.infectious_h_risk = _getarg(dct, 'infectious_h_risk', 0.15)

        # Infectious icu event

        ## Prob of recovering if in (or needing) icu
        self.recover_before_death = _getarg(dct, 'recover_before_death', 0.3)

        ## Prob of remaining infectious icu per day
        self.infectious_i_risk = _getarg(dct, 'infectious_i_risk', 0.15)

        if dct:
            raise TypeError('Unknown argument ' + str(dct))

        # Trackers
        self.num_isolated = 0
        self.num_deisolated = 0
        self.num_traced = 0
        self.peak = 0
        self.peak_total_infections = 0
        self.peak_iter = 0
        self.results = []

    def init_agents(self, n):
        self.agents = []
        for i in range(n):
            self.agents.append(Agent(self, i))

    def make_infected(self, infected, infector):
        self.agents[infected].health = Health.EXPOSED
        self.agents[infected].infector = (self.iteration, infector)
        self.agents[infector].infected_by_me.append((self.iteration, infected))

    def event_infect(self):
        neighbors = round(self.k_match / 2)
        agents = self.agents
        indices = list(range(len(agents)))
        npr.shuffle(indices)
        infected = [None] * len(agents)
        for i in indices:
            if agents[i].health > Health.EXPOSED and \
               agents[i].health < Health.RECOVERED:
                _from = max(0, i - neighbors)
                _to = min(i + 1 + neighbors, len(agents))
                for j in range(_from, _to):
                    if agents[j].health > Health.SUSCEPTIBLE:
                        continue
                    risk = min(1 - agents[i].isolated, 1 - agents[j].isolated) \
                           * ((agents[i].risk_infecting[agents[i].health] +
                              agents[j].risk_infection) / 2.0)
                    if random.random() < risk:
                        infected[j] = i
        for i in range(len(infected)):
            if infected[i]:
                self.make_infected(i, infected[i])

    def event_test(self):
        agents = [a for a in self.agents if
                  a.test_res_iter is None and
                  a.health > Health.INFECTIOUS_A and
                  a.health < Health.RECOVERED]
        for a in agents:
            if random.random() < self.test_likelihood:
                a.test_res_iter = self.iteration + \
                                max(npr.poisson(self.mean_test), self.min_test)


    def event_isolate(self):
        if self.min_before_isolate > 0:
           infectious = len([a for a in self.agents
                             if a.health > Health.INFECTIOUS_A])
           if infectious < self.min_before_isolate:
               return

        self.min_before_isolate = 0
        agents = [a for a in self.agents if
                  a.test_res_iter == self.iteration and
                  a.health < Health.RECOVERED]
        for a in agents:
            self.num_isolated += 1
            a.isolate(self)

    def event_deisolate(self):
        agents = [a for a in self.agents if
                  a.isolation_iter == self.iteration]
        for a in agents:
            self.num_deisolated += 1
            a.deisolate(self)

    def event_trace(self):
        if self.min_before_trace > 0:
           infectious = len([a for a in self.agents
                             if a.health > Health.INFECTIOUS_A])
           if infectious < self.min_before_trace:
               return
        self.min_before_trace = 0
        agents = [a for a in self.agents if
                  a.test_res_iter and a.test_res_iter <= self.iteration and
                  a.health > Health.INFECTIOUS_A and
                  a.health < Health.RECOVERED]
        neighbors = round(self.k_match / 2)
        for a in agents:
            _from = max(0, a.id - neighbors)
            _to = min(a.id + 1 + neighbors, len(self.agents))
            for i in range(_from, _to):
                if i != a.id:
                    if random.random() < self.trace_effective:
                        self.num_traced += 1
                        self.num_isolated += 1
                        self.agents[i].isolate(self)

    def event_exposed(self):
        agents = [a for a in self.agents if a.health == Health.EXPOSED]
        for a in agents:
            if random.random() < self.exposed_risk:
                a.health = Health.INFECTIOUS_A

    def event_infectious_a(self):
        agents = [a for a in self.agents if a.health == Health.INFECTIOUS_A]
        for a in agents:
            if random.random() < self.infectious_a_risk:
                if a.asymptomatic is True:
                    a.health = Health.RECOVERED
                    a.inf_end_iter = self.iteration
                else:
                    a.health = Health.INFECTIOUS_S

    def event_infectious_s(self):
        agents = [a for a in self.agents if a.health == Health.INFECTIOUS_S]
        for a in agents:
            if random.random() < self.infectious_s_risk:
                if a.recover_before_hospital is True:
                    a.health = Health.RECOVERED
                    a.inf_end_iter = self.iteration
                else:
                    a.health = Health.INFECTIOUS_H

    def event_infectious_h(self):
        agents = [a for a in self.agents if a.health == Health.INFECTIOUS_H]
        for a in agents:
            if random.random() < self.infectious_h_risk:
                if a.recover_before_icu is True:
                    a.health = Health.RECOVERED
                    a.inf_end_iter = self.iteration
                else:
                    a.health = Health.INFECTIOUS_I

    def event_infectious_i(self):
        agents = [a for a in self.agents if a.health == Health.INFECTIOUS_I]
        for a in agents:
            if random.random() < self.infectious_i_risk:
                a.inf_end_iter = self.iteration
                if a.recover_before_death is True:
                    a.health = Health.RECOVERED
                else:
                    a.health = Health.DEAD

    def stats(self, forced=False):
        if forced or self.iteration % self.stats_frequency == 0:
            infections = len([a for a in self.agents
                              if a.health > Health.SUSCEPTIBLE and
                              a.health < Health.RECOVERED])
            if self.peak < infections:
                self.peak = infections
                self.peak_total_infections = len([a for a in self.agents
                                            if a.health > Health.SUSCEPTIBLE])
                self.peak_iter = self.iteration

    def report(self, forced=False):
        if forced or (self.iteration % self.report_frequency) == 0:
            dead = len([a for a in self.agents if a.health == Health.DEAD])
            recovered = len([a for a in self.agents
                             if a.health == Health.RECOVERED])
            infected = len([a for a in self.agents
                            if a.health > Health.SUSCEPTIBLE and
                            a.health < Health.RECOVERED])
            if self.verbose:
                print(f"{self.name}Iteration {self.iteration}. "
                      f"Total infected {infected + recovered + dead}.")
            self.results.append({
                'iteration': self.iteration,
                'infected': infected,
                'recovered': recovered,
                'dead': dead}
            )

    def iterate(self):
        from_ = self.iteration
        to_ = self.iteration + self.iterations
        for self.iteration in range(from_, to_):
            self.stats()
            self.report()
            self.event_infect()
            self.event_test()
            self.event_isolate()
            self.event_deisolate()
            self.event_trace()
            self.event_exposed()
            self.event_infectious_a()
            self.event_infectious_s()
            self.event_infectious_h()
            self.event_infectious_i()
        self.stats(forced=True)
        self.report(forced=True)
        return self.results

    def simulate(self):
        self.init_agents(self.num_agents)
        self.iteration = 0
        return self.iterate()

    def R0(self, period=10):
        infections = [len(a.infected_by_me) for a in self.agents
                if a.inf_end_iter and a.inf_end_iter <= period]
        method_1 = sum(infections) / len(infections)
        n = len(self.agents)
        method_2 = 1.0 / ((n - self.peak_total_infections) / n)
        return (method_1, method_2)


######################################
# Scenarios for paper
######################################


def run_scenarios(num=10, scenarios=['a','b','c', 'd', 'e', 'f', 'g', 'h', 'i'],
                  filename_prefix="scenarios",
                  parms={}):

    results = {}
    simulations = {}
    unique_file = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    filename = filename_prefix + unique_file
    if 'iterations' not in parms:
        parms['iterations'] = 500

    if 'recover_before_hospital' not in parms:
        parms['recover_before_hospital'] = 1.0

    if 'verbose' not in parms:
        parms['verbose'] = True

    for scenario in scenarios:
        simulations[scenario] = []
        if parms['verbose'] is True:
            print("Scenario:", scenario)
        for i in range(num):
            simulation_parms = copy.deepcopy(parms)
            if scenario == 'a': # No interventions
                simulation_parms['trace_effective'] = 0.0
                simulation_parms['min_isolation'] = 0.0
                simulation_parms['max_isolation'] = 0.0
            elif scenario == 'b': # Max intervention
                simulation_parms['trace_effective'] = 1.0
                simulation_parms['min_isolation'] = 1.0
                simulation_parms['max_isolation'] = 1.0
            elif scenario == 'c': # Realistic intervention lower bound
                simulation_parms['trace_effective'] = 0.1
                simulation_parms['min_isolation'] = 0.2
                simulation_parms['max_isolation'] = 0.9
                simulation_parms['mean_test'] = 7
            elif scenario == 'd': # Realistic intervention upper bound
                simulation_parms['trace_effective'] = 0.6
                simulation_parms['min_isolation'] = 0.7
                simulation_parms['max_isolation'] = 1.0
            elif scenario == 'e': # Realistic intervention middle-of-the-road
                simulation_parms['trace_effective'] = 0.3
                simulation_parms['min_isolation'] = 0.4
                simulation_parms['max_isolation'] = 1.0
                simulation_parms['mean_test'] = 3
            elif scenario == 'f': # Realistic from middle of epidemic
                simulation_parms['trace_effective'] = 0.3
                simulation_parms['min_isolation'] = 0.4
                simulation_parms['max_isolation'] = 1.0
                simulation_parms['mean_test'] = 3
                simulation_parms['min_before_trace'] = 1000
                simulation_parms['min_before_isolate'] = 1000
            elif scenario == 'g': # Realistic from middle of epidemic no trace
                simulation_parms['trace_effective'] = 0.0
                simulation_parms['min_isolation'] = 0.4
                simulation_parms['max_isolation'] = 1.0
                simulation_parms['mean_test'] = 3
                simulation_parms['min_before_trace'] = 1000
                simulation_parms['min_before_isolate'] = 1000
            elif scenario == 'h': # Scenario a sensitivity test of R0
                simulation_parms['trace_effective'] = 0.0
                simulation_parms['min_isolation'] = 0.0
                simulation_parms['max_isolation'] = 0.0
                simulation_parms['k_match'] = random.randint(30, 46)
            elif scenario == 'i': # Scenario e sensitivity test of R0
                simulation_parms['trace_effective'] = 0.3
                simulation_parms['min_isolation'] = 0.4
                simulation_parms['max_isolation'] = 1.0
                simulation_parms['mean_test'] = 3
                simulation_parms['k_match'] = random.randint(30, 46)
            simulations[scenario].append(Simulation(**simulation_parms))
        results[scenario] = [s.simulate() for s in simulations[scenario]]
    output = {
        "simulations": simulations,
        "results": results
    }
    with open(filename + ".bin", "wb") as f:
        pickle.dump(output, f)
    with open(filename + ".csv", 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['scenario', 'active', 'recovered', 'dead', 'total'])
        for scenario, result in results.items():
            for r in result:
                csvwriter.writerow([scenario,
                                    r[-1]['infected'],
                                    r[-1]['recovered'],
                                    r[-1]['dead'],
                                    sum([r[-1]['infected'],
                                         r[-1]['recovered'],
                                         r[-1]['dead']])])
    return output



def _run_sim(s):
    r = s.copy()
    r.update(r['simulation'].simulate()[-1])
    del(r['simulation'])
    r.update(s['parms'])
    del(r['parms'])
    return r

def run_sensitivity(num_jiggles=10, num_sims_per_jiggle=5,
                    scenarios=['j','k','l',],
                    filename_prefix="scenarios",
                    parms_user={},
                    jiggle_parms_user={}):

    jiggle_parms = {
        'test_likelihood': [0.01, 0.99, random.uniform],
        'mean_test': [1, 8, random.randint],
        'min_test': [0, 2, random.randint],
        'isolation_period': [6, 14, random.randint],
        'exposed_risk': [0.1, 0.9, random.uniform],
        'asymptomatic': [0.05, 0.95, random.uniform],
        'infectious_a_risk': [0.1, 0.9, random.uniform],
        'infectious_s_risk': [0.05, 0.7, random.uniform],
    }

    l = len(jiggle_parms)
    jiggle_parms.update(jiggle_parms_user)
    if len(jiggle_parms) > l:
        raise TypeError('Unknown arguments ' + str(jiggle_parms_user))

    parms = {
        'iterations': 500,
        'recover_before_hospital': 1.0,
        'verbose': True,
    }
    parms.update(parms_user)

    results = {}
    unique_file = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    filename = filename_prefix + unique_file

    results = []
    for i in range(num_jiggles):
        if parms['verbose']:
            print("Jiggle:", i)
        for key, value in jiggle_parms.items():
            if len(value) == 3:
                rng = value[2]
            else:
                rng = random.uniform
            parms[key] = rng(value[0], value[1])

        simulations = []
        for j in range(num_sims_per_jiggle):
            for scenario in scenarios:
                if scenario == 'j': # No interventions
                    parms['trace_effective'] = 0.0
                    parms['min_isolation'] = 0.0
                    parms['max_isolation'] = 0.0
                elif scenario == 'k': # Max intervention
                    parms['trace_effective'] = 1.0
                    parms['min_isolation'] = 1.0
                    parms['max_isolation'] = 1.0
                elif scenario == 'l': # Realistic middle-of-the-road
                    parms['trace_effective'] = 0.3
                    parms['min_isolation'] = 0.4
                    parms['max_isolation'] = 1.0
                else:
                    continue
                parms['name'] =  "_".join([str(i), str(j), scenario, " "])
                s = Simulation(**parms)
                simulations.append({
                    'jiggle': i,
                    'run': j,
                    'scenario': scenario,
                    'parms': copy.deepcopy(parms),
                    'simulation': s
                })
        with Pool() as p:
            results = results + p.map(_run_sim, simulations)

    with open(filename + ".csv", 'w', newline='') as csvfile:
        header = list(results[0].keys())
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(header)
        for r in results:
            csvwriter.writerow(r.values())
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='abm', description='Run simulations')
    parser.add_argument("--sensitivity",
                        dest="run",
                        action='store_const',
                        const=run_sensitivity,
                        default=None,
                        help="Run sensitivity tests")
    parser.add_argument("--num_runs",
                        dest="num_runs",
                        type=int,
                        default=10,
                        help='Number of simulations per scenario or jiggle')
    parser.add_argument("--num_jiggles",
                        dest="num_jiggles",
                        type=int,
                        default=10,
                        help='Number of sensitivity jiggles')
    parser.add_argument("--scenarios",
                        dest="scenarios",
                        type=str,
                        default="abcdefghijkl",
                        help='Scenarios to run')
    parser.add_argument("--filename",
                        dest="filename",
                        type=str,
                        default="scenarios",
                        help='Output csv prefix name')
    parser.add_argument("--parameters",
                        dest="parameters",
                        type=str,
                        nargs="*",
                        default=None,
                        help='Customise parameters')
    args = parser.parse_args()

    parms = {}
    if args.parameters:
        for p in args.parameters:
            key, value = p.split("=")
            if value == "True":
                parms[key] = True
            elif value == "False":
                parms[key] = False
            elif "." in value:
                parms[key] = float(value)
            else:
                parms[key] = int(value)

    if args.run:
        run_sensitivity(args.num_jiggles, args.num_runs, list(args.scenarios),
                        args.filename, parms)
    else:
        run_scenarios(args.num_runs, args.scenarios, filename, parms)
