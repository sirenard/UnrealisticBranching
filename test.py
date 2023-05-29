import os
import subprocess
import argparse
import math

OUT = None

def dprint(text="", end="\n"):
    text = text.replace(".",",")
    print(text, end=end)
    if OUT is not None:
        with open(OUT, 'a') as file:
            file.write(text+end)

def get_n_nodes(instance,combination):
    if "strat" in combination:
        command = CLASSIC_SOLVE_COMMAND.format(combination["strat"], instance)

    elif "out" in combination:
        if combination["out"] is None:
            command = UB_SOLVE_COMMAND.format(combination["rd"], instance)
        else:
            out = os.path.join(combination["out"], instance.split("/")[-1]+".csv")
            command = GENERATE_DS_COMMAND.format(combination["rd"], instance, out)
    elif "model" in combination:
        command = MODEL_SOLVE.format(instance, combination["model"],combination["secondstrat"], combination["threshold"], combination["md"])

    try:
        output = subprocess.check_output(command, shell=True)
    except:
        return -1, -1, -1

    lines = str(output).split("\\n")
    nnode = -1
    time = -1
    gap = -1
    for line in lines:
        if "Solving Nodes" in line:
            nnode = int(line.split()[-1])
        if "Solving Time (sec)" in line:
            time = float(line.split()[-1])
        if time != -1 and "Gap" in line:
            els = line.split()
            for i in range(len(els)-1):
                if els[i] == ":":
                    if els[i+1] != "infinite":
                        gap = float(els[i+1])

    return nnode, time, gap


def evaluate(instances, recursiondepths, outputpath=None, maxdepths=[-1], modelpath=[], thresholds=[0.5], strategies=[]):
    assert not(len(modelpath) > 0 and outputpath is not None), "outputpath and modelpath cannot be specified at the same time"
    assert not(len(modelpath)>1 and len(maxdepths)>1), "cannot test several depth and several models at the same time"
    assert len(maxdepths)==1 or len(thresholds)==1, "cannot vary the maxdepth and the thresholds"
    assert not (len(modelpath)==1 and len(strategies)>1) or maxdepths==[-1], "cannot test a model with several second strategies"

    combinations = []
    for rd in recursiondepths:
        if rd == 0:
            for strat in strategies:
                combinations.append({
                    "strat": strat,
                })
        else:
            combinations.append({
                "out": outputpath,
                "rd": rd
            })

    for model in modelpath:
        for md in maxdepths:
            for threshold in thresholds:
                combinations.append({
                    "model": model,
                    "md": md,
                    "threshold": threshold,
                    "secondstrat": strategies[0]
                })

    dprint("{:<30}".format("instance"), end="")
    for combination in combinations:
        if "strat" in combination:
            title = combination["strat"]
        elif "rd" in combination:
            title = "rd={}".format(combination["rd"])
        else:
            if len(maxdepths) > 1:
                title = "md={}".format(combination["md"])
            elif len(thresholds) > 1:
                title = "th={}".format(combination["threshold"])
            else:
                title = combination["model"].split("/")[-1]

        dprint("{:<15} {:<8} {:<8}".format(title, 't(s)', 'gap(%)'), end="")

    dprint()
    for instance in instances:
        for i, combination in enumerate(combinations):
            nNode, time, gap = get_n_nodes(instance, combination)
            if i == 0:
                dprint("{:<30}".format(instance.split("/")[-1]), end="")

            dprint("{:<15} {:<8} {:<8}".format(nNode, time, gap), end="")

        dprint()



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    argParser = argparse.ArgumentParser()
    argParser.add_argument("-e", "--executable", help="path to scip learning executable", required=True)
    argParser.add_argument("-d", "--dir", help="path to lp instances", required=True)
    argParser.add_argument("-m", "--models", help="path to models to use", nargs="+", default=[])
    argParser.add_argument("-o", "--output", help="path to output directory", default=None)
    argParser.add_argument("-f", "--froms", type=int, help="index from", default=0)
    argParser.add_argument("-t", "--to", type=int, help="index to", default=None)
    argParser.add_argument('-r','--recursiondepth', type=int, nargs='+', help='list of recursion (in addition to 0)', default=None)
    argParser.add_argument('-a','--abc', type=int, nargs=3, help='list a b c', required=True)
    argParser.add_argument('-md','--maxdepth', type=int, nargs="+", help='list of max depths tested', default=[-1])
    argParser.add_argument('-th','--thresholds', type=float, nargs="+", help='threshold values', default=[0.5])
    argParser.add_argument('-bs','--branchingstrategies', type=str, nargs="+", help='Branching strategies used', default=["vanillafullstrong"])
    argParser.add_argument("-eps", "--epsilon", type=float, help="Epsilon for exploration", default=0.8)
    argParser.add_argument("-exp", "--exploremethod", help="Name of the branching scheme used for exploration", type=str, default="pscost")
    argParser.add_argument("-tl", "--timelimit", type=float, help="if solve take more than that time, don't do it", default=10)
    argParser.add_argument("-ltl", "--leaftimelimit", type=float, help="leaftime limit (coefficient of the running time)", default=-1)
    argParser.add_argument("-nl", "--nodelimit", type=int, help="if solve take more than in node number that time, don't do it", default=-1)
    argParser.add_argument("-of", "--outputfile", type=str, help="path for prints", default=None)


    args = argParser.parse_args()

    EXECUTABLE = args.executable
    CLASSIC_SOLVE_COMMAND = EXECUTABLE + " -c 'set branching {{}} priority 99999999' -c 'read {{}}' -c 'set limits time {}' -c 'set limits nodes {}' -c 'optimize' -c 'quit'".format(args.timelimit, args.nodelimit)
    UB_SOLVE_COMMAND = EXECUTABLE + "  -c 'set branching unrealistic leaftimelimit {}'  -c 'set branching unrealistic maxdepth {}' -c 'set limits nodes {}' -c 'set branching unrealistic recursion {{}}' -c 'read {{}}' -c 'optimize' -c 'quit'".format(args.leaftimelimit, args.maxdepth[0], args.nodelimit)
    GENERATE_DS_COMMAND = EXECUTABLE + " -c 'set limits time {}'  -c 'set branching unrealisticTrained maxdepth {}' -c 'set limits nodes {}' -c 'set data scoremethod n' -c 'set branching unrealistic leaftimelimit {}' -c 'set branching {} priority 70000' -c 'set data epsilon {}' -c 'set branching unrealistic recursion {{}}' -c 'generateDataset {{}} {{}} {} {} {}' -c 'quit'".format(args.timelimit, args.maxdepth[0],args.nodelimit, args.leaftimelimit, args.exploremethod, args.epsilon, *args.abc)
    MODEL_SOLVE = EXECUTABLE + " -c 'set limits time {}' -c 'set limits nodes {}' -c 'read {{}}' -c 'loadmodel {{}} {} {} {}' -c 'set branching unrealistic priority 0' -c 'set branching {{}} priority 9999998' -c 'set branching unrealisticTrained priority 9999999' -c 'set branching unrealisticTrained maxbounddist 1'  -c 'set branching unrealisticTrained threshold {{}}' -c 'set branching unrealisticTrained maxdepth {{}}' -c 'optimize' -c 'quit'".format(args.timelimit, args.nodelimit, *args.abc)

    dir = args.dir
    instances = os.listdir(dir)
    instances.sort()
    for i in range(len(instances)):
        instances[i] = os.path.join(dir, instances[i])
    to = args.to
    if to is None:
        instances = instances[args.froms:]
    else:
        instances = instances[args.froms:to]


    recursions = []
    if args.recursiondepth is not None:
        recursions = recursions + args.recursiondepth

    models = args.models

    OUT = args.outputfile
    if OUT is not None:
        with open(OUT, "w") as file:
            file.write("init\n")



    # instances = ["/home/simon/Téléchargements/bppgen/ds1/BPP_60.lp"]
    evaluate(instances, recursions, outputpath=args.output, modelpath=models,maxdepths=args.maxdepth, thresholds=args.thresholds, strategies=args.branchingstrategies)

    #generate_ds(instances, output, 0, 0, 1)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
