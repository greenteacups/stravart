import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Polygon
from shapely.ops import unary_union

def generate_grid(n_rows, n_cols, spacing):
    x_coords = np.arange(0, n_cols * spacing, spacing)
    y_coords = np.arange(0, n_rows * spacing, spacing)
    grid_points = np.array(np.meshgrid(x_coords, y_coords)).T.reshape(-1, 2)
    return grid_points

def generate_lines(grid_points, n_rows, n_cols):
    lines = []
    for i in range(n_rows):
        for j in range(n_cols - 1):
            start = tuple(grid_points[i * n_cols + j])
            end = tuple(grid_points[i * n_cols + j + 1])
            lines.append((start, end))
    for j in range(n_cols):
        for i in range(n_rows - 1):
            start = tuple(grid_points[i * n_cols + j])
            end = tuple(grid_points[(i + 1) * n_cols + j])
            lines.append((start, end))
    return lines

def compute_overlap(grid_points, shape, n_rows, n_cols, spacing):
    overlaps = []
    for i in range(n_rows - 1):
        for j in range(n_cols - 1):
            square_coords = [
                (j * spacing, i * spacing),
                ((j + 1) * spacing, i * spacing),
                ((j + 1) * spacing, (i + 1) * spacing),
                (j * spacing, (i + 1) * spacing)
            ]
            grid_square = Polygon(square_coords)
            square_area = grid_square.area

            intersection = grid_square.intersection(shape)
            intersection_area = intersection.area
            percentage_overlap = (intersection_area / square_area) * 100
            if percentage_overlap > 0:
                overlaps.append(((i, j), percentage_overlap, grid_square))
    
    return overlaps

def read_shape_coords_from_file(file_path):
    data = pd.read_csv(file_path)
    coords = data[['X', 'Y']].values.tolist()
    return coords

# Parameters
n_rows = 20
n_cols = 20
spacing = 0.05

# Generate grid and lines
grid_points = generate_grid(n_rows, n_cols, spacing)
grid_lines = generate_lines(grid_points, n_rows, n_cols)

# Define the shape polygon
file_path = 'drawing-coords.dat'
shape_coords = read_shape_coords_from_file(file_path)
shape = Polygon(shape_coords)

# Compute overlap
overlapping_squares = compute_overlap(grid_points, shape, n_rows, n_cols, spacing)

# Output the list of overlapping squares, their percentage overlaps, and intersection points
print("Overlapping squares and their percentage overlap:")
for (square, overlap, polygon) in overlapping_squares:
    print(f"Square {square}: {overlap:.2f}% overlap")

# Merge squares with more than X% overlap into a single polygon
polygons_to_merge = [polygon for (square, overlap, polygon) in overlapping_squares if overlap > 50]
merged_polygon = unary_union(polygons_to_merge) if polygons_to_merge else None

# Plot the grid and shape
plt.figure(figsize=(8, 8))
for grid_line in grid_lines:
    start, end = grid_line
    #plt.plot([start[0], end[0]], [start[1], end[1]], 'bo-')
    
plt.plot(*shape.exterior.xy, 'ro--')

for (square, overlap, polygon) in overlapping_squares:
    i, j = square
    plt.fill_between([j * spacing, (j + 1) * spacing], [i * spacing, i * spacing], [(i + 1) * spacing, (i + 1) * spacing], color='g', alpha=0.3)

# Highlight the merged polygon
if merged_polygon and merged_polygon.is_valid:
    x, y = merged_polygon.exterior.xy
    plt.plot(x, y, color='k', linewidth=2)

plt.gca().set_aspect('equal', adjustable='box')
plt.title("Grid of Points Connected by Lines with Shaded Overlapping Areas")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)
plt.show()

