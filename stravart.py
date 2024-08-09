import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Polygon
from shapely.ops import unary_union

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
        if percentage_overlap > 0:
            overlaps.append((i, percentage_overlap, polygon))
    return overlaps

def read_shape_coords_from_file(file_path):
    data = pd.read_csv(file_path)
    coords = data[['X', 'Y']].values.tolist()
    return coords

# Define the shape polygon
file_path_shape = 'drawing-coords.dat'
shape_coords = read_shape_coords_from_file(file_path_shape)
shape = Polygon(shape_coords)

# Read grid polygons from file
file_path_grid = 'road-coords.dat'
grid_polygons = read_polygons_from_file(file_path_grid)

# Compute overlap
overlapping_polygons = compute_overlap(grid_polygons, shape)

# Output the list of overlapping polygons, their percentage overlaps, and intersection points
print("Overlapping polygons and their percentage overlap:")
for (index, overlap, polygon) in overlapping_polygons:
    print(f"Polygon {index + 1}: {overlap:.2f}% overlap")

# Merge polygons with more than X% overlap into a single polygon
polygons_to_merge = [polygon for (index, overlap, polygon) in overlapping_polygons if overlap > 50]
merged_polygon = unary_union(polygons_to_merge) if polygons_to_merge else None

# Plot the grid polygons and shape
plt.figure(figsize=(8, 8))
for polygon in grid_polygons:
    x, y = polygon.exterior.xy
    plt.plot(x, y, 'bo-')

plt.plot(*shape.exterior.xy, 'ro--')

for (index, overlap, polygon) in overlapping_polygons:
    x, y = polygon.exterior.xy
    plt.fill(x, y, color='g', alpha=0.3)

# Highlight the merged polygon
if merged_polygon and merged_polygon.is_valid:
    x, y = merged_polygon.exterior.xy
    plt.plot(x, y, color='k', linewidth=2)

plt.gca().set_aspect('equal', adjustable='box')
plt.title("Grid Polygons with Shaded Overlapping Areas")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)
plt.show()

