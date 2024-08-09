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
    
def translate_shape(shape, x_translation, y_translation):
    return transform(lambda x, y: (x + x_translation, y + y_translation), shape)

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

def process_shape_and_grid(shape_file, grid_file, angle_of_rotation, scaling_factor, overlapper, x_translation, y_translation, plot_results_flag=True):
    shape_coords = read_shape_coords_from_file(shape_file)
    shape = Polygon(shape_coords)

    # Apply rotation, scaling, and translation
    rotated_shape = rotate_shape(shape, angle_of_rotation)
    scaled_shape = scale_shape(rotated_shape, scaling_factor)
    translated_shape = translate_shape(scaled_shape, x_translation, y_translation)

    grid_polygons = read_polygons_from_file(grid_file)
    overlapping_polygons = compute_overlap(grid_polygons, translated_shape)

    overlap_values = [overlap for _, overlap, _ in overlapping_polygons]
    average_overlap = sum(overlap_values) / len(overlap_values) if overlap_values else 0

    polygons_to_merge = [polygon for (index, overlap, polygon) in overlapping_polygons if overlap > overlapper]
    merged_polygon = unary_union(polygons_to_merge) if polygons_to_merge else None

    # Handle export with error handling
    try:
        if merged_polygon and merged_polygon.is_valid:
            export_polygon_to_file(merged_polygon, 'merged_polygon.dat')
    except Exception as e:
        print(f"Error exporting merged polygon: {e}")
        merged_polygon = None

    if plot_results_flag:
        plot_results(grid_polygons, translated_shape, overlapping_polygons, merged_polygon)

    weighted_average = merged_polygon.area / (shape.area*scaling_factor) if merged_polygon else 0
    
    return weighted_average if merged_polygon else 0  # Return 0 or another value if there's an issue with the export


# Run:
""
process_shape_and_grid(
    shape_file = 'drawing-coords.dat',
    grid_file = 'road-coords.dat',
    angle_of_rotation = 32.54,	            # Angle in degrees
    scaling_factor = 0.5,                   # Scaling factor
    overlapper = 80.0,                      # % overlap to include block
    x_translation = 0.1,
    y_translation = -0.1,
    plot_results_flag = 1                   # Set to False to skip plotting
)

"""
def objective_function(params, *args):
    shape_file, grid_file, overlapper = args
    angle_of_rotation, x_translation, y_translation, scaling_factor = params
    weighted_average = process_shape_and_grid(
        shape_file=shape_file,
        grid_file=grid_file,
        angle_of_rotation=angle_of_rotation,
        scaling_factor=scaling_factor,
        overlapper=overlapper,
        x_translation=x_translation,
        y_translation=y_translation,
        plot_results_flag=False
    )
    return -weighted_average  # Minimize the negative of the weighted average

# Define parameters for optimization
shape_file = 'drawing-coords.dat'
grid_file = 'road-coords.dat'
overlapper = 80.0

# Set the initial guess and bounds for the angle_of_rotation, translations, and scaling factor
initial_guess = [0, 0, 0, 0.8]  # [angle_of_rotation, x_translation, y_translation, scaling_factor]
bounds = [(-180, 180), (-0.1, 0.1), (-0.1, 0.1), (0.5, 1.0)]  # Angle of rotation, translations, and scaling factor

# Run basinhopping
result = basinhopping(
    objective_function,
    initial_guess,
    niter=500,
    T=1.0,
    stepsize=10.0,
    minimizer_kwargs={
        'method': 'L-BFGS-B',
        'bounds': bounds,
        'args': (shape_file, grid_file, overlapper)
    }
)

print(f"Optimal angle of rotation: {result.x[0]:.2f} degrees")
print(f"Optimal x translation: {result.x[1]:.2f}")
print(f"Optimal y translation: {result.x[2]:.2f}")
print(f"Optimal scaling factor: {result.x[3]:.2f}")
print(f"Maximum weighted average: {-result.fun:.2f}")
"""
