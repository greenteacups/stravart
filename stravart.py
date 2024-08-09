import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Polygon
from shapely.ops import unary_union, transform
from shapely.affinity import rotate, scale
from scipy.optimize import basinhopping

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

    overlap_values = [overlap for _, overlap, _ in overlapping_polygons]
    average_overlap = sum(overlap_values) / len(overlap_values) if overlap_values else 0
    #print("Average overlap:", average_overlap/100)
    #print("# Overlapped polygons:", len(overlap_values)/len(grid_polygons))

    polygons_to_merge = [polygon for (index, overlap, polygon) in overlapping_polygons if overlap > overlapper]
    merged_polygon = unary_union(polygons_to_merge) if polygons_to_merge else None

    # Handle export with error handling
    try:
        if merged_polygon and merged_polygon.is_valid:
            export_polygon_to_file(merged_polygon, 'merged_polygon.dat')
    except Exception as e:
        print(f"Error exporting merged polygon: {e}")
        # Set a flag or adjust the result as needed
        merged_polygon = None

    if plot_results_flag:
        plot_results(grid_polygons, scaled_shape, overlapping_polygons, merged_polygon)

    # Computing weighted function for optimizer
    #weighted_average = ((average_overlap/100)*1.0 + (len(overlap_values)/len(grid_polygons))*1.0)*0.5
    weighted_average = merged_polygon.area/shape.area
    #print("Average average:", shape.area, "test", merged_polygon.area)
    
    return weighted_average if merged_polygon else 0  # Return 0 or another value if there's an issue with the export


# Run:
""
process_shape_and_grid(
    shape_file='drawing-coords.dat',
    grid_file='road-coords.dat',
    angle_of_rotation=12, #32.56,  # Angle in degrees
    scaling_factor=0.8,     # Scaling factor
    overlapper=50.0,        # % overlap to include block
    plot_results_flag=1     # Set to False to skip plotting
)

"""
def objective_function(angle_of_rotation, *args):
    shape_file, grid_file, scaling_factor, overlapper = args
    weighted_average = process_shape_and_grid(
        shape_file=shape_file,
        grid_file=grid_file,
        angle_of_rotation=angle_of_rotation,
        scaling_factor=scaling_factor,
        overlapper=overlapper,
        plot_results_flag=False  # Disable plotting for optimization
    )
    return -weighted_average  # Minimize the negative of the weighted average

# Define parameters for optimization
shape_file = 'drawing-coords.dat'
grid_file = 'road-coords.dat'
scaling_factor = 0.8
overlapper = 50.0

# Set the initial guess and bounds for the angle_of_rotation
initial_guess = 0
bounds = [(-180, 180)]  # Angle of rotation can vary from -180 to 180 degrees

# Run basinhopping
result = basinhopping(
    objective_function,
    initial_guess,
    niter=50,  # Number of iterations
    T=1.0,  # Temperature for the random jumps
    stepsize=10.0,  # Size of the random jumps
    minimizer_kwargs={
        'method': 'L-BFGS-B',
        'bounds': bounds,
        'args': (shape_file, grid_file, scaling_factor, overlapper)
    }
)

print(f"Optimal angle of rotation: {result.x[0]:.2f} degrees")
print(f"Maximum weighted average: {-result.fun:.2f}")
"""

