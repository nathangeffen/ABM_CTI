import argparse
import subprocess
import pandas as pd
import statistics
import datetime
import time

def read_file(filename='out.csv'):
    return pd.read_csv(filename)

def get_scenarios(df, num_scenarios=5):
    scenarios = []
    for i in range(num_scenarios):
        scenarios.append(df.loc[(df['scenario'] == i) &
                                (df['iteration'] == 500)])
    return scenarios


def get_jiggle_stats(scenarios):
    jiggle_stats = {'mean': [], 'median': [], 'min': [], 'max': [], 'std': []}
    for s in scenarios:
        jiggle_stats['mean'].append(s.groupby(['jiggle']).mean())
        jiggle_stats['median'].append(s.groupby(['jiggle']).median())
        jiggle_stats['min'].append(s.groupby(['jiggle']).min())
        jiggle_stats['max'].append(s.groupby(['jiggle']).max())
        jiggle_stats['std'].append(s.groupby(['jiggle']).std())
    return jiggle_stats

def cmp_scenarios(jiggle_stats, stat='mean', s_1=0, s_2=1,
                  col='total_infected'):
    wins = 0
    for i in range(len(jiggle_stats[stat][s_1])):
        wins += jiggle_stats[stat][s_1][col][i] < jiggle_stats[stat][s_2][col][i]
    return wins

def count_bad(jiggle_stats, stat='mean', s_1=0, s_2=1, threshold=0.9,
              min_size=0, col='total_infected'):
    bad = 0
    for i in range(len(jiggle_stats[stat][s_1])):
        if (stat == 'mean' and \
           # jiggle_stats[stat][s_1]['meant_test'][i] < 6 and \
           jiggle_stats[stat][s_1][col][i] < threshold * \
           jiggle_stats[stat][s_2][col][i] and \
           jiggle_stats[stat][s_2][col][i] > min_size):
                bad += 1
    return bad


def stats(filename='out.csv'):
    df = read_file(filename)
    scenarios = get_scenarios(df)
    for i in range(len(scenarios)):
        print("Mean scenario", i, ": ", scenarios[i]['total_infected'].mean())
        print("Median scenario", i, ": ",scenarios[i]['total_infected'].median())
        print("Min scenario", i, ": ", scenarios[i]['total_infected'].min())
        print("Max scenario", i, ": ", scenarios[i]['total_infected'].max())
        print("Std scenario", i, ": ", scenarios[i]['total_infected'].std())
    return scenarios

def good_vs_bad(scenarios, threshold=0.9, min_size=0):
    jiggle_stats = get_jiggle_stats(scenarios)
    n = len(jiggle_stats['mean'][0])

    print("Mean: Min v OnlyI", cmp_scenarios(jiggle_stats, s_1=0, s_2=1))
    print("Mean: Min v PoorCTI", cmp_scenarios(jiggle_stats, s_1=0, s_2=2))
    print("Mean: Min v OKCTI", cmp_scenarios(jiggle_stats, s_1=0, s_2=3))
    print("Mean: Min v MaxCTI", cmp_scenarios(jiggle_stats, s_1=0, s_2=4))

    print("Mean: Onlyiso v PoorCTI", cmp_scenarios(jiggle_stats, s_1=1, s_2=2))
    print("Mean: OnlyI v OKCTI", cmp_scenarios(jiggle_stats, s_1=1, s_2=3))
    print("Mean: OnlyI v MaxCTI", cmp_scenarios(jiggle_stats, s_1=1, s_2=4))

    print("Mean: PoorCTI v OKCTI", cmp_scenarios(jiggle_stats, s_1=2, s_2=3))
    print("Mean: PoorCTI v MaxCTI", cmp_scenarios(jiggle_stats, s_1=2, s_2=4))

    print("Mean: OKCTI v MaxCTI", cmp_scenarios(jiggle_stats, s_1=3, s_2=4))


    print("Median: Min v OnlyI", cmp_scenarios(jiggle_stats, s_1=0, s_2=1,
                                             stat='median'))
    print("Median: Min v PoorCTI", cmp_scenarios(jiggle_stats, s_1=0, s_2=2,
                                             stat='median'))
    print("Median: Min v OKCTI", cmp_scenarios(jiggle_stats, s_1=0, s_2=3,
                                             stat='median'))
    print("Median: Min v MaxCTI", cmp_scenarios(jiggle_stats, s_1=0, s_2=4,
                                             stat='median'))

    print("Median: Onlyiso v PoorCTI", cmp_scenarios(jiggle_stats, s_1=1, s_2=2,
                                             stat='median'))
    print("Median: OnlyI v OKCTI", cmp_scenarios(jiggle_stats, s_1=1, s_2=3,
                                             stat='median'))
    print("Median: OnlyI v MaxCTI", cmp_scenarios(jiggle_stats, s_1=1, s_2=4,
                                             stat='median'))

    print("Median: PoorCTI v OKCTI", cmp_scenarios(jiggle_stats, s_1=2, s_2=3,
                                             stat='median'))
    print("Median: PoorCTI v MaxCTI", cmp_scenarios(jiggle_stats, s_1=2, s_2=4,
                                             stat='median'))

    print("Median: OKCTI v MaxCTI", cmp_scenarios(jiggle_stats, s_1=3, s_2=4,
                                             stat='median'))

    print("Std: Min v OnlyI", cmp_scenarios(jiggle_stats, s_1=0, s_2=1,
                                             stat='std'))
    print("Std: Min v PoorCTI", cmp_scenarios(jiggle_stats, s_1=0, s_2=2,
                                             stat='std'))
    print("Std: Min v OKCTI", cmp_scenarios(jiggle_stats, s_1=0, s_2=3,
                                             stat='std'))
    print("Std: Min v MaxCTI", cmp_scenarios(jiggle_stats, s_1=0, s_2=4,
                                             stat='std'))

    print("Std: Onlyiso v PoorCTI", cmp_scenarios(jiggle_stats, s_1=1, s_2=2,
                                             stat='std'))
    print("Std: OnlyI v OKCTI", cmp_scenarios(jiggle_stats, s_1=1, s_2=3,
                                             stat='std'))
    print("Std: OnlyI v MaxCTI", cmp_scenarios(jiggle_stats, s_1=1, s_2=4,
                                             stat='std'))

    print("Std: PoorCTI v OKCTI", cmp_scenarios(jiggle_stats, s_1=2, s_2=3,
                                             stat='std'))
    print("Std: PoorCTI v MaxCTI", cmp_scenarios(jiggle_stats, s_1=2, s_2=4,
                                             stat='std'))

    print("Std: OKCTI v MaxCTI", cmp_scenarios(jiggle_stats, s_1=3, s_2=4,
                                             stat='std'))

    print("Bad: Min v OnlyI", count_bad(jiggle_stats, s_1=0, s_2=1,
                                      threshold=threshold,
                                      min_size=min_size))
    print("Bad: Min v PoorCTI", count_bad(jiggle_stats, s_1=0, s_2=2,
                                      threshold=threshold,
                                      min_size=min_size))
    print("Bad: Min v OKCTI", count_bad(jiggle_stats, s_1=0, s_2=3,
                                      threshold=threshold,
                                      min_size=min_size))
    print("Bad: Min v MaxCTI", count_bad(jiggle_stats, s_1=0, s_2=4,
                                      threshold=threshold,
                                      min_size=min_size))

    print("Bad: Onlyiso v PoorCTI", count_bad(jiggle_stats, s_1=1, s_2=2,
                                      threshold=threshold,
                                      min_size=min_size))
    print("Bad: OnlyI v OKCTI", count_bad(jiggle_stats, s_1=1, s_2=3,
                                      threshold=threshold,
                                      min_size=min_size))
    print("Bad: OnlyI v MaxCTI", count_bad(jiggle_stats, s_1=1, s_2=4,
                                      threshold=threshold,
                                      min_size=min_size))

    print("Bad: PoorCTI v OKCTI", count_bad(jiggle_stats, s_1=2, s_2=3,
                                      threshold=threshold,
                                      min_size=min_size))
    print("Bad: PoorCTI v MaxCTI", count_bad(jiggle_stats, s_1=2, s_2=4,
                                             threshold=threshold,
                                             min_size=min_size))

    print("Bad: OKCTI v MaxCTI", count_bad(jiggle_stats, s_1=3, s_2=4,
                                             threshold=threshold,
                                             min_size=min_size))

    print("Good: Min v OnlyI", count_bad(jiggle_stats, s_2=0, s_1=1,
                                      threshold=threshold,
                                      min_size=min_size))
    print("Good: Min v PoorCTI", count_bad(jiggle_stats, s_2=0, s_1=2,
                                      threshold=threshold,
                                      min_size=min_size))
    print("Good: Min v OKCTI", count_bad(jiggle_stats, s_2=0, s_1=3,
                                      threshold=threshold,
                                      min_size=min_size))
    print("Good: Min v MaxCTI", count_bad(jiggle_stats, s_2=0, s_1=4,
                                      threshold=threshold,
                                      min_size=min_size))

    print("Good: Onlyiso v PoorCTI", count_bad(jiggle_stats, s_2=1, s_1=2,
                                      threshold=threshold,
                                      min_size=min_size))
    print("Good: OnlyI v OKCTI", count_bad(jiggle_stats, s_2=1, s_1=3,
                                      threshold=threshold,
                                      min_size=min_size))
    print("Good: OnlyI v MaxCTI", count_bad(jiggle_stats, s_2=1, s_1=4,
                                      threshold=threshold,
                                      min_size=min_size))

    print("Good: PoorCTI v OKCTI", count_bad(jiggle_stats, s_2=2, s_1=3,
                                      threshold=threshold,
                                      min_size=min_size))
    print("Good: PoorCTI v MaxCTI", count_bad(jiggle_stats, s_2=2, s_1=4,
                                             threshold=threshold,
                                             min_size=min_size))

    print("Good: OKCTI v MaxCTI", count_bad(jiggle_stats, s_2=3, s_1=4,
                                             threshold=threshold,
                                             min_size=min_size))



def main():
    n = datetime.datetime.now()
    file_id = '-'.join([str(n.year), str(n.month), str(n.day),
                        str(n.hour), str(n.minute), str(n.second)])
    parser = argparse.ArgumentParser(description='Analyse hetgen output')
    parser.add_argument("--run_hetgen",
                        dest="run_hetgen",
                        type=bool,
                        default=False,
                        help="Whether to run hetgen")
    parser.add_argument("--compile_hetgen",
                        dest="compile_hetgen",
                        type=bool,
                        default=True,
                        help="Whether to compile hetgen")
    parser.add_argument("--jiggles",
                        dest="jiggles",
                        type=int,
                        default=10,
                        help="Number of jiggles")
    parser.add_argument("--runs",
                        dest="runs",
                        type=int,
                        default=1,
                        help="Number of runs per jiggle")
    parser.add_argument("--threads",
                        dest="threads",
                        type=int,
                        default=0,
                        help="Number of threads (0 = system determined)")
    parser.add_argument("--seed",
                        dest="seed",
                        type=int,
                        default=-2,
                        help="Random seed for hetgen (-2 uses random device)")
    parser.add_argument("--agents",
                        dest="agents",
                        type=int,
                        default=10000,
                        help="Number of agents")
    parser.add_argument("--initial_infections",
                        dest="initial_infections",
                        type=int,
                        default=10,
                        help="Initial infections")
    parser.add_argument("--k_assort",
                        dest="k_assort",
                        type=str,
                        default="32:44",
                        help="Neighbouring agents exposed to infection")
    parser.add_argument("--k_unassort",
                        dest="k_unassort",
                        type=str,
                        default="0",
                        help="Agents exposed to random infection")
    parser.add_argument("--filename",
                        dest="filename",
                        default="out-" + file_id + ".csv",
                        help="File to analyse")
    parser.add_argument("--threshold",
                        dest="threshold",
                        type=float,
                        default=0.9,
                        help="Threshold prop for classifying good v bad")
    parser.add_argument("--min_size",
                        dest="min_size",
                        type=int,
                        default=0,
                        help="Minimum infections for classifying good v bad")
    args = parser.parse_args()


    if (args.run_hetgen):
        if args.compile_hetgen:
            subprocess.run(['make', 'abm_hetgen'])
        output = open(args.filename, "w")
        t0 = time.time();
        subprocess.run([
            './abm_hetgen',
            f'--threads={args.threads}',
            f'--jiggles={args.jiggles}',
            f'--seed={args.seed}',
            f'--runs={args.runs}',
            f'--agents={args.agents}',
            f'--initial_infections={args.initial_infections}',
            '--iterations=500',
            '--report_frequency=600',
            '--recover_before_hospital=1.0',
            '--prob_test_exposed=0.1:0.3',
            '--prob_test_infectious_a=0.1:0.3',
            '--prob_test_infectious_s=0.3:0.9',
            '--mean_test=1:6',
            '--min_test=0:2',
            '--isolation_period=6:14',
            '--exposed_risk=0.1:0.9',
            '--asymptomatic=0.05:0.95',
            '--infectious_a_risk=0.1:0.9',
            '--infectious_s_risk=0.1:0.9',
            f'--k_assort={args.k_assort}',
            f'--k_unassort={args.k_unassort}',

            # No intervention
            '--trace_effective=0.0',
            '--min_isolation=0.0',
            '--max_isolation=0.0',

            # Isolation only
            '+',
            '--trace_effective=0.0',
            '--min_isolation=0.7',
            '--max_isolation=1.0',

            # Minimal CTI
            '+',
            '--trace_effective=0.1',
            '--min_isolation=0.7',
            '--max_isolation=1.0',

            # Moderate CTI
            '+',
            '--trace_effective=0.3',
            '--min_isolation=0.7',
            '--max_isolation=1.0',

            # Maximum CTI
            '+',
            '--trace_effective=1.0',
            '--min_isolation=1.0',
            '--max_isolation=1.0'

        ], text=True, stdout=output)
        t1 = time.time()
        output.close()
        print("Time take: ", t1 - t0, "seconds")
    scenarios = stats(args.filename)
    good_vs_bad(scenarios, args.threshold, args.min_size)

if __name__ == "__main__":
    main()
