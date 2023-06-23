"""
\file IMP/saxs/SolventAccessibleSurface.h \brief

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

import math
import numpy as np
from scipy import spatial

class SolventAccessibleSurface:
    def compute_accessibility(self, point, probe_radius, density):
        sphere_dots = self.create_sphere_dots(probe_radius, density)
        intersecting_count = 0
        for sphere_center, radius in sphere_dots:
            if self.is_intersecting(point, sphere_center, point[3], radius):
                intersecting_count += 1
        accessibility = 1.0 - (intersecting_count / len(sphere_dots))
        return accessibility

    def get_solvent_accessible_surface(self, points, probe_radius=1.8, density=5.0):
        return self.get_solvent_accessibility(points, probe_radius, density)

    def get_solvent_accessibility(self, ps, probe_radius=1.8, density=5.0):
        res = []
        # coordinates = [p.get_coordinates() for p in ps]
        # radii = [p.get_radius() for p in ps]
        coordinates = [p.coordinates for p in ps]
        radii = [p.radius for p in ps]

        # generate sphere dots for radii present in the ps set
        self.create_sphere_dots(ps, density)

        # init grid
        grid = self.init_grid(coordinates)

        # compute surface
        max_radius = 3.0
        for i in range(len(ps)):
            atom_radius = radii[i]
            radius = atom_radius + 2 * probe_radius + max_radius
            # query
            # bb = self.get_bounding_box(coordinates[i], radius)
            inside_rad = grid.query_ball_point(coordinates[i], radius)
            # lb, ub = grid.query(bb[0]), grid.query(bb[1])
            neighbours1, neighbours2 = [], []
            for mol_index in inside_rad:
            # for it in grid.indexes_begin(lb, ub):
            #     for vIndex in range(len(grid[it])):
                # mol_index = grid[it][vIndex]
                radius_sum1 = atom_radius + radii[mol_index]
                radius_sum2 = radius_sum1 + 2 * probe_radius
                dist2 = self.get_squared_distance(coordinates[i], coordinates[mol_index])
                if dist2 < radius_sum1 * radius_sum1:
                    neighbours1.append(mol_index)
                elif dist2 < radius_sum2 * radius_sum2:
                    neighbours2.append(mol_index)

            ratio = (atom_radius + probe_radius) / atom_radius
            spoints = np.array(self.get_sphere_dots(atom_radius))

            dotNum = 0
            for s_index in range(len(spoints)):
                probe_center = coordinates[i] + ratio * spoints[s_index]
                # check for intersection with neighbors1
                collides = False
                for n_index in range(len(neighbours1)):
                    if self.is_intersecting(probe_center, coordinates[neighbours1[n_index]],
                                            probe_radius, radii[neighbours1[n_index]]):
                        collides = True
                        break
                if not collides:  # check for intersection with neighbors2
                    for n_index in range(len(neighbours2)):
                        if self.is_intersecting(probe_center, coordinates[neighbours2[n_index]],
                                                probe_radius, radii[neighbours2[n_index]]):
                            collides = True
                            break
                if not collides:
                    dotNum += 1
            res.append(float(dotNum) / len(spoints))
        return res

    def init_grid(self, coordinates):
        grid = spatial.KDTree(coordinates)
        # # Define grid dimensions
        # size = int(max(max(coord) for coord in coordinates) - min(min(coord) for coord in coordinates) + 10)
        # grid_size = (size, size, size)  # Size of the grid in each dimension
        # # Create an empty grid
        # grid = np.zeros(grid_size)

        # # grid = defaultdict(list)
        # for i in range(len(coordinates)):
        #     # grid_index = self.get_nearest_index(grid, coordinates[i])
        #     grid[coordinates[i]] = i
        return grid


    def get_nearest_index(self, grid_coordinates, target_coordinates):
        # Convert the grid coordinates and target coordinates to NumPy arrays
        grid_coordinates = np.array(grid_coordinates)
        target_coordinates = np.array(target_coordinates)

        # Calculate the Euclidean distances between the target coordinates and all grid coordinates
        distances = np.linalg.norm(grid_coordinates - target_coordinates, axis=1)

        # Find the index corresponding to the minimum distance
        nearest_index = np.argmin(distances)

        return nearest_index


    def get_bounding_box(self, coordinates, radius=0.0):
        coordinates = [coordinates]
        min_point = [min(coord[:][i] for coord in coordinates) - radius for i in range(3)]
        max_point = [max(coord[:][i] for coord in coordinates) + radius for i in range(3)]
        return [min_point, max_point]

    def is_intersecting(self, center1, center2, radius1, radius2):
        squared_radius_sum = (radius1 + radius2) * (radius1 + radius2)
        squared_dist = self.get_squared_distance(center1, center2)
        if abs(squared_radius_sum - squared_dist) < 0.0001:
            return False
        if squared_radius_sum > squared_dist:
            return True
        return False

    def get_squared_distance(self, v1, v2):
        squared_dist = sum((x1 - x2) ** 2 for x1, x2 in zip(v1, v2))
        return squared_dist

    def create_sphere_dots(self, ps, density):
        radii2type = {}
        sphere_dots = []

        for p in ps:
            r = p.radius # [3] # radius in XYZR
            if r not in radii2type:
                type_ = len(radii2type)
                radii2type[r] = type_
                dots = self.create_sphere_dots_for_radius(r, density)
                sphere_dots.append(dots)

        self.radii2type = radii2type
        self.sphere_dots = sphere_dots

    def get_sphere_dots(self, r):
        if r not in self.radii2type:
            raise ValueError(f"SolventAccessibleSurface: can't find sphere dots for radius {r}")
        index = self.radii2type[r]
        return self.sphere_dots[index]

    def create_sphere_dots_for_radius(self, radius, density):
        res = []
        num_equat = 2 * math.pi * radius * math.sqrt(density)
        vert_count = 0.5 * num_equat

        for i in range(int(vert_count)):
            phi = (math.pi * i) / vert_count
            z = math.cos(phi)
            xy = math.sin(phi)
            horz_count = xy * num_equat
            for j in range(int(horz_count - 1)):
                teta = (2 * math.pi * j) / horz_count
                x = xy * math.cos(teta)
                y = xy * math.sin(teta)
                res.append((radius * x, radius * y, radius * z))
        return res

# def get_nearest_index(self, pt):
#     if self.get_is_dense():
#         ei = self.get_nearest_extended_index(pt)
#         return self.get_index(ei)
#     else:
#         raise ValueError("get_nearest_index only works on dense grids.")
# 
# def get_nearest_extended_index(self, pt):
#     if self.get_is_bounded():
#         ei = self.get_extended_index(pt)
#         for i in range(pt.get_dimension()):
#             ei[i] = max(0, ei[i])
#             ei[i] = min(len(self) - 1, ei[i])
#         return ei
#     else:
#         raise ValueError("get_nearest_index only works on bounded grids.")
