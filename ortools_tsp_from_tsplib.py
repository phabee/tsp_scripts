# https://developers.google.com/optimization/routing/tsp

"""Simple travelling salesman problem between cities."""
# read instance from tsplib format

from __future__ import print_function
from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp

import pandas as pd
import scipy.spatial as sp
import numpy as np


def create_data_model(dima):
    """Stores the data for the problem."""
    data = {}
    data['distance_matrix'] = dima
    data['num_vehicles'] = 1
    data['depot'] = 0
    return data


def print_solution(manager, routing, solution, dima):
    """Prints solution on console."""
    print('Objective: {} miles'.format(solution.ObjectiveValue()))
    index = routing.Start(0)
    plan_output = 'Route for vehicle 0:\n'
    route_distance = 0
    while not routing.IsEnd(index):
        plan_output += ' {} ->'.format(manager.IndexToNode(index))
        previous_index = index
        index = solution.Value(routing.NextVar(index))
        route_distance += routing.GetArcCostForVehicle(previous_index, index, 0)
    plan_output += ' {}\n'.format(manager.IndexToNode(index))
    plan_output += 'Route distance: {} miles\n'.format(route_distance)
    print(plan_output)


def main():
    # load tsp instance
    tsp = pd.read_csv("./TSP_Instances/bier127.tsp", sep='\t', skiprows=2, names=['nodeId', 'lat', 'lng'])
    tsp = tsp.sort_values(by='nodeId', inplace=False)

    A = tsp[['lat', 'lng']].to_numpy()*100
    dima = sp.distance_matrix(A, A)

    """Entry point of the program."""
    # Instantiate the data problem.
    data = create_data_model(dima)

    # Create the routing index manager.
    manager = pywrapcp.RoutingIndexManager(len(data['distance_matrix']),
                                           data['num_vehicles'], data['depot'])

    # Create Routing Model.
    routing = pywrapcp.RoutingModel(manager)

    def distance_callback(from_index, to_index):
        """Returns the distance between the two nodes."""
        # Convert from routing variable Index to distance matrix NodeIndex.
        from_node = manager.IndexToNode(from_index)
        to_node = manager.IndexToNode(to_index)
        return data['distance_matrix'][from_node][to_node]

    transit_callback_index = routing.RegisterTransitCallback(distance_callback)

    # Define cost of each arc.
    routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

    # Setting first solution heuristic.
    search_parameters = pywrapcp.DefaultRoutingSearchParameters()
    # search_parameters.first_solution_strategy = (
        # routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC)

    search_parameters.first_solution_strategy = (
        routing_enums_pb2.FirstSolutionStrategy.AUTOMATIC)

    search_parameters.local_search_metaheuristic = (
        routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH)

    search_parameters.time_limit.seconds = 30

    # Solve the problem.
    solution = routing.SolveWithParameters(search_parameters)

    # Print solution on console.
    if solution:
        print_solution(manager, routing, solution, dima)

if __name__ == '__main__':
    main()
