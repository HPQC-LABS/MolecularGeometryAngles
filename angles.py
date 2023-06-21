import argparse
import numpy as np
import itertools
from collections import deque

def cyclic_permutations(obj):
    obj = list(obj)
    cyclic_perms = [obj]

    obj_dequed = deque(obj)

    while True:
        obj_dequed.rotate(1)
        if list(obj_dequed) != obj:
            cyclic_perms.append(list(obj_dequed))
        else:
            break

    return cyclic_perms

def compute_angle(coords):
    vec1 = np.array(coords[0]) - np.array(coords[1])
    vec2 = np.array(coords[2]) - np.array(coords[1])

    _, angle = compute_angle_from_vectors(vec1, vec2)

    return angle

def compute_angle_from_vectors(vec1, vec2):
    cos_angle = (
        np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    )

    angle_rads = np.arccos(cos_angle)
    angle_deg = np.rad2deg(angle_rads)

    return angle_rads, angle_deg

def compute_normal_to_plane(a, b, c):
    vec1 = a - b
    vec2 = c - b

    n = np.cross(vec1, vec2)
    return n

def compute_angle_between_planes(plane1, plane2):
    a1, b1, c1 = plane1
    a2, b2, c2 = plane2

    n1 = compute_normal_to_plane(a1, b1, c1)
    n2 = compute_normal_to_plane(a2, b2, c2)

    _, angle = compute_angle_from_vectors(n1, n2)

    return angle

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="ComputeAngles",
        description="Compute planar and dihedral angles from XYZ file",
        epilog="",
    )
    parser.add_argument("filename")
    args = parser.parse_args()

    filename = args.filename

    elements = []
    coords = []

    with open(filename) as f:
        contents = f.readlines()

        for line in contents:
            line_split = [x for x in line.rstrip('\n').split(" ") if x != '']

            if len(line_split) != 4: continue

            if not line_split[0].isalpha(): continue
            if len(line_split[0]) > 3: continue

            for i in range(1, 4):
                try:
                    float(line_split[i])
                except:
                    continue

            elements.append(line_split[0])
            coords.append([float(line_split[i]) for i in range(1, 4)])

    ids = list(np.arange(0, len(elements)))

    elements_enum = [f"{elem}{(id+1)}" for elem, id in zip(elements, ids)]

    three_combs = list(itertools.combinations(ids, 3))
    four_combs = list(itertools.combinations(ids, 4))

    planar_angles_elems = []
    planar_angles_degs = []

    for triangle in three_combs:
        angles = cyclic_permutations(triangle)

        for angle_ids in angles:
            angle_elements = [elements_enum[i] for i in angle_ids]
            angle_coords = [coords[i] for i in angle_ids]

            planar_angles_elems.append("-".join(angle_elements))

            angle = compute_angle(angle_coords)

            planar_angles_degs.append(np.round(angle, 5))

    non_repeating_planar_angles_degs = list(set(planar_angles_degs))
    non_repeating_planar_angles_degs.sort()

    planar_angles_elems = np.array(planar_angles_elems)
    planar_angles_degs = np.array(planar_angles_degs)

    planar_angles_dict = {}
    
    for angle in non_repeating_planar_angles_degs:
        ids = np.ravel(np.argwhere(planar_angles_degs == angle))

        planar_angles_dict[angle] = list(planar_angles_elems[ids.astype(int)])

    for pyramid in four_combs:
        print(pyramid)

        planes = list(itertools.combinations(pyramid, 3))
        pairs_of_planes = list(itertools.combinations(planes, 2))
        
        for pair in pairs_of_planes:
            plane1_ids, plane2_ids = pair
            
            plane1_coords = np.array([coords[i] for i in plane1_ids])
            plane2_coords = np.array([coords[i] for i in plane2_ids])

            angle = compute_angle_between_planes(plane1_coords, plane2_coords)
            
    
