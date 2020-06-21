import os
import sys
import getopt
import time
import subprocess
import argparse

import numpy as np

from toric_model import Action
from toric_model import Toric_code


def usage():
    print('usage:  sdafklasfjlksdaklfjsa')

def distance(coord1, coord2):
    x_distance = np.abs(coord1[0] - coord2[0])
    y_distance = np.abs(coord1[1] - coord2[1])
    return x_distance + y_distance

def get_non_periodic_distance(coord1, coord2):
    distance = np.abs(coord1 - coord2)
    return distance

def get_periodic_distance_for_one_axis(coord1, coord2, system_size):
    distance_to_border_coord1 = system_size-1-coord2
    distance_to_border_coord2 = system_size-1-coord1
    distance1 = coord1 + 1 + distance_to_border_coord1
    distance2 = coord2 + 1 + distance_to_border_coord2
    return min([distance1, distance2])

def get_shortest_distance(coord1, coord2, system_size):

    x_distance = min([get_periodic_distance_for_one_axis(coord1[0], coord2[0], system_size), get_non_periodic_distance(coord1[0], coord2[0])])
    y_distance = min([get_periodic_distance_for_one_axis(coord1[1], coord2[1], system_size), get_non_periodic_distance(coord1[1], coord2[1])])
    distance = x_distance + y_distance
    return distance

def generate_node_indices(edges, nbr_of_edges, nbr_of_nodes):
    start_nodes = [np.repeat(i, nbr_of_nodes-i-1) for i in range(nbr_of_nodes)]
    edges[:,0] = np.concatenate(start_nodes)

    end_nodes = [np.arange(i+1, nbr_of_nodes) for i in range(nbr_of_nodes)]
    edges[:,1] = np.concatenate(end_nodes)

    return edges

def get_distances(MWPM_edges, edges_no_periodic, edges, nbr_of_nodes):
    if (nbr_of_nodes == 2):
        MWPM_edge_indices = 0
        periodic_distances = edges[0][2]
        non_periodic_distances = edges_no_periodic[0][2]
    else:
        MWPM_edge_indices = np.where((edges_no_periodic[:, 0:2] == MWPM_edges[:, None]).all(-1))[1]
        periodic_distances = edges[MWPM_edge_indices, 2]
        non_periodic_distances = edges_no_periodic[MWPM_edge_indices, 2]

    return non_periodic_distances, periodic_distances

def generate_syndrome(toric_code, p):
    toric_code.generate_random_error(p)
    return toric_code

def generate_syndrome_smart(toric_code, n):
    toric_code.generate_n_random_errors(n)
    return toric_code

def generate_MWPM(matrix, system_size):
    defect_coords = np.array(np.nonzero(matrix)).T

    nbr_of_nodes = len(defect_coords)

    nbr_of_edges= int(nbr_of_nodes*(nbr_of_nodes-1)/2)

    edges = np.zeros((nbr_of_edges, 3))

    edges_no_periodic = np.zeros((nbr_of_edges, 3))

    edges = generate_node_indices(edges, nbr_of_edges, nbr_of_nodes)
    edges_no_periodic = generate_node_indices(edges_no_periodic, nbr_of_edges, nbr_of_nodes)

    shortest_distances = [get_shortest_distance(coord1, coord2, system_size)
            for i, coord1 in enumerate(defect_coords[:-1]) for coord2 in defect_coords[i+1:, :]]

    shortest_non_periodic_distances= [distance(coord1, coord2)
            for i, coord1 in enumerate(defect_coords[:-1]) for coord2 in defect_coords[i+1:, :]]

    edges[:, 2] = shortest_distances
    edges_no_periodic[:, 2] = shortest_non_periodic_distances


    processId = os.getpid()
    PATH = str(processId) + 'edges.TXT'
    OUTPUT_PATH = str(processId) +'output.TXT'

    header_str = "{} {}".format(nbr_of_nodes, nbr_of_edges)
    np.savetxt(PATH, edges, fmt='%i', header=header_str, comments='')

    subprocess.call(["./blossom5-v2.05.src/blossom5", "-e", PATH, "-w", OUTPUT_PATH, "-V"], stdout=open(os.devnull, 'wb'))

    MWPM_edges = np.loadtxt(OUTPUT_PATH, skiprows=1, dtype=int)
    MWPM_edges = MWPM_edges.reshape((int(nbr_of_nodes/2), 2))

    os.remove(PATH)
    os.remove(OUTPUT_PATH)

    return MWPM_edges, edges, edges_no_periodic, defect_coords

def eliminate_defect_pair(toric_code, start_coord, end_coord, diff_coord, matrix_index, system_size):
    coord = start_coord
    size_diff = (np.array([system_size, system_size])-np.abs(diff_coord))
    size_diff[size_diff == 5] = 0

    nbr_of_vertical_steps = min(size_diff[0], np.abs(diff_coord[0]))
    nbr_of_horizontal_steps = min(size_diff[1], np.abs(diff_coord[1]))

    action_index = (not matrix_index)*2 + 1

    for i in range(nbr_of_vertical_steps):
        direction = size_diff[0]>np.abs(diff_coord[0])
        coord[0] = (coord[0] + direction-(not matrix_index))%system_size
        pos = [matrix_index, coord[0], coord[1]]
        coord[0] = (coord[0] - (not direction)+(not matrix_index))%system_size

        action = Action(pos, action_index)
        toric_code.step(action)
        toric_code.syndrom('state')

    for i in range(nbr_of_horizontal_steps):
        if (diff_coord[1] < 0):
            direction = size_diff[1]<np.abs(diff_coord[1])
            coord[1] = (coord[1] + (direction)-(not matrix_index))%system_size
            pos = [int(not matrix_index), coord[0], coord[1]]
            coord[1] = (coord[1] - (not direction))+(not matrix_index)%system_size
        else:
            direction = size_diff[1]>np.abs(diff_coord[1])
            coord[1] = (coord[1] + (direction)-(not matrix_index))%system_size
            pos = [int(not matrix_index), coord[0], coord[1]]
            coord[1] = (coord[1] - (not direction)+(not matrix_index))%system_size

        action = Action(pos, action_index)
        toric_code.step(action)
        toric_code.syndrom('state')
    return toric_code

def generate_solution(MWPM_edges, defect_coords, toric_code, matrix_index, system_size):
    start_coords = defect_coords[MWPM_edges[:, 0].T, :]
    end_coords = defect_coords[MWPM_edges[:, 1].T, :]
    diff_coords = end_coords-start_coords

    if (np.sum(np.sum(toric_code.state[matrix_index])) > 0):
        for start_coord, end_coord, diff_coord in zip(start_coords, end_coords, diff_coords):
            toric_code = eliminate_defect_pair(toric_code, start_coord, end_coord, diff_coord, matrix_index, system_size)
    return toric_code

def main(args):
    #TODO: add support for arguments
    # parser = argparse.ArgumentParser()
    # parser.add_argument('test2', help='test adsfadd')
    # parser.parse_args()
    # print(args.test2)
    # try:
    #     opts, args = getopt.getopt(argv, "hn:d:p:", ["help"])
    # except getopt.GetoptError as err:
    #     print(err)
    #     usage()
    #     sys.exit(2)
    # print(opts)
    # print(args)

    # for o, a in opts:
    #     print(o)
    #     if o in ("-h", "--help"):
    #         usage()
    #     elif o == "-d":
    #         system_size = int(a)
    #         print("d = {}".format(system_size))
    #     elif o == "-n":
    #         nbr_of_iterations = int(float(a))
    #         print(nbr_of_iterations)


    p_errors = [0.05, 0.06, 0.07]
    system_size = int(5)
    nbr_of_iterations = int(float(1e2))

    print(p_errors)
    ground_state_kept_list = []


    for p in p_errors:
        print(p)
        ground_states = 0
        for _ in range(nbr_of_iterations):
            toric_code = Toric_code(system_size)
            nbr_of_vertex_nodes = 0
            nbr_of_plaquette_nodes = 0

            toric_code = generate_syndrome(toric_code, p)
            n = 3

            MWPM_edges_vertex = []
            edges_vertex = []
            edges_no_periodic_vertex = []
            defect_coords_vertex = []

            MWPM_edges_plaquette = []
            edges_plaquette = []
            edges_no_periodic_plaquette = []
            defect_coords_plaquette = []

            if np.sum(np.sum(toric_code.state[0])) > 0:
                MWPM_edges_vertex, edges_vertex, edges_no_periodic_vertex, defect_coords_vertex = generate_MWPM(toric_code.state[0], system_size)
            if np.sum(np.sum(toric_code.state[1])) > 0:
                MWPM_edges_plaquette, edges_plaquette, edges_no_periodic_plaquette, defect_coords_plaquette = generate_MWPM(toric_code.state[1], system_size)

            if len(MWPM_edges_vertex) > 0:
                toric_code = generate_solution(MWPM_edges_vertex, defect_coords_vertex, toric_code, 0, system_size)
            if len(MWPM_edges_plaquette) > 0:
                toric_code = generate_solution(MWPM_edges_plaquette, defect_coords_plaquette, toric_code, 1, system_size)


            toric_code.eval_ground_state()
            ground_states += toric_code.ground_state

        ground_state_kept_list.append(ground_states/nbr_of_iterations)

    timestamp = time.ctime()
    try:
        os.mkdir('results')
    except FileExistsError:
        pass

    PATH_ground2 = 'results/p_succes_MWPM_d={}_time={}.TXT'.format(system_size, timestamp)

    np.savetxt(PATH_ground2, ground_state_kept_list, fmt='%e', comments='')
    print(ground_state_kept_list)

if __name__ == "__main__":
    main(sys.argv[1:])
