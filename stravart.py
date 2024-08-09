import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Polygon
from shapely.ops import unary_union, transform
from shapely.affinity import rotate, scale

def read_polygons_from_file(file_path):
    polygons = []
    current_polygon = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("Polygon"):
                if current_polygon:
                    polygons.append(Polygon(current_polygon))
                    current_polygon = []
            else:
                coords = line.strip().split(',')
                current_polygon.append((float(coords[0]), float(coords[1])))
        if current_polygon:  # Add the last polygon
            polygons.append(Polygon(current_polygon))
    return polygons

def compute_overlap(polygons, shape):
    overlaps = []
    for i, polygon in enumerate(polygons):
        intersection = polygon.intersection(shape)
        intersection_area = intersection.area
        polygon_area = polygon.area
        percentage_overlap = (intersection_area / polygon_area) * 100
        if percentage_overlap > 0.0:
            overlaps.append((i, percentage_overlap, polygon))
    return overlaps

def read_shape_coords_from_file(file_path):
    data = pd.read_csv(file_path)
    coords = data[['X', 'Y']].values.tolist()
    return coords

def export_polygon_to_file(polygon, file_path):
    with open(file_path, 'w') as file:
        for x, y in polygon.exterior.coords:
            file.write(f"{x},{y}\n")

def rotate_shape(shape, angle):
    return rotate(shape, angle, origin='centroid', use_radians=False)

def scale_shape(shape, scale_factor):
    return scale(shape, xfact=scale_factor, yfact=scale_factor, origin='centroid')

def plot_results(grid_polygons, scaled_shape, overlapping_polygons, merged_polygon):
    plt.figure(figsize=(8, 8))
    for polygon in grid_polygons:
        x, y = polygon.exterior.xy
        plt.plot(x, y, 'bo-')

    plt.plot(*scaled_shape.exterior.xy, 'ro--')

    for (index, overlap, polygon) in overlapping_polygons:
        x, y = polygon.exterior.xy
        plt.fill(x, y, color='g', alpha=0.3)

    if merged_polygon and merged_polygon.is_valid:
        x, y = merged_polygon.exterior.xy
        plt.plot(x, y, color='k', linewidth=2)

    plt.gca().set_aspect('equal', adjustable='box')
    plt.title(f"Grid Polygons with Overlaps (Rotated and Scaled)")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.grid(True)
    plt.show()

def process_shape_and_grid(shape_file, grid_file, angle_of_rotation, scaling_factor, overlapper, plot_results_flag=True):
    shape_coords = read_shape_coords_from_file(shape_file)
    shape = Polygon(shape_coords)

    rotated_shape = rotate_shape(shape, angle_of_rotation)
    scaled_shape = scale_shape(rotated_shape, scaling_factor)

    grid_polygons = read_polygons_from_file(grid_file)
    overlapping_polygons = compute_overlap(grid_polygons, scaled_shape)

    #print("Overlapping polygons and their percentage overlap:")
    #for (index, overlap, polygon) in overlapping_polygons:
        #print(f"Polygon {index + 1}: {overlap:.2f}% overlap")
    
    overlap_values = [overlap for _, overlap, _ in overlapping_polygons]
    average_overlap = sum(overlap_values) / len(overlap_values) if overlap_values else 0
    print("Average overlap:", average_overlap/100)
    print("# Overlapped polygons:", len(overlap_values)/len(grid_polygons))

    polygons_to_merge = [polygon for (index, overlap, polygon) in overlapping_polygons if overlap > overlapper]
    merged_polygon = unary_union(polygons_to_merge) if polygons_to_merge else None

    if merged_polygon and merged_polygon.is_valid:
        export_polygon_to_file(merged_polygon, 'merged_polygon.dat')

    if plot_results_flag:
        plot_results(grid_polygons, scaled_shape, overlapping_polygons, merged_polygon)

    # Computing weighted function for optimser:
    weighted_average = ((average_overlap/100)*1.0 + (len(overlap_values)/len(grid_polygons))*1.0)*0.5
    print("Average average:", weighted_average)


# Run:
process_shape_and_grid(
    shape_file='drawing-coords.dat',
    grid_file='road-coords.dat',
    angle_of_rotation=-20,  # Angle in degrees
    scaling_factor=0.8,     # Scaling factor
    overlapper=50.0,        # % overlap to include block
    plot_results_flag=0     # Set to False to skip plotting
)

