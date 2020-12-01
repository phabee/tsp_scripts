import pants
import math
import random
import pandas as pd

# Ant colony optimization for TSP
# documentation see here: https://pypi.org/project/ACO-Pants/

def euclidean(a, b):
    return math.sqrt(pow(a[1] - b[1], 2) + pow(a[0] - b[0], 2))

def main():
    dat = pd.read_csv("./TSP_Instances/pr1002.tsp", sep='\t', skiprows=2, names=['nodeId', 'lat', 'lng'])
    nodes = [tuple(x) for x in dat.to_numpy()]
    # nodes.append((dat.lat.array, dat.lng.array, dat.nodeId.array))
    world = pants.World(nodes, euclidean)
    solver = pants.Solver()
    solution = solver.solve(world)
    print(solution.distance)
    print(solution.tour)    # Nodes visited in order
    print(solution.path)    # Edges taken in order

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
