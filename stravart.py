import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import euclidean

def generate_grid(n_rows, n_cols, spacing):
    """
    Generates a grid of points connected by lines.
    
    Parameters:
    - n_rows: Number of rows in the grid.
    - n_cols: Number of columns in the grid.
    - spacing: Distance between adjacent points.
    
    Returns:
    - grid_points: A 2D array of points in the grid.
    """
    x_coords = np.arange(0, n_cols * spacing, spacing)
    y_coords = np.arange(0, n_rows * spacing, spacing)
    grid_points = np.array(np.meshgrid(x_coords, y_coords)).T.reshape(-1, 2)
    
    return grid_points

def generate_lines(grid_points, n_rows, n_cols):
    """
    Generates a list of lines connecting the grid points.
    
    Parameters:
    - grid_points: A 2D array of points in the grid.
    - n_rows: Number of rows in the grid.
    - n_cols: Number of columns in the grid.
    
    Returns:
    - lines: A list of tuples representing lines (start_point, end_point).
    """
    lines = []
    
    # Generate horizontal lines
    for i in range(n_rows):
        for j in range(n_cols - 1):
            start = tuple(grid_points[i * n_cols + j])
            end = tuple(grid_points[i * n_cols + j + 1])
            lines.append((start, end))
    
    # Generate vertical lines
    for j in range(n_cols):
        for i in range(n_rows - 1):
            start = tuple(grid_points[i * n_cols + j])
            end = tuple(grid_points[(i + 1) * n_cols + j])
            lines.append((start, end))
    
    return lines

# Function to calculate the distance between two lines
def line_distance(line1, line2):
    (x1, y1), (x2, y2) = line1
    (a1, b1), (a2, b2) = line2
    return min(euclidean((x1, y1), (a1, b1)), euclidean((x1, y1), (a2, b2)), euclidean((x2, y2), (a1, b1)), euclidean((x2, y2), (a2, b2)))


# Compute:
n_rows = 5
n_cols = 5
spacing = 1

grid_points = generate_grid(n_rows, n_cols, spacing)
grid_lines = generate_lines(grid_points, n_rows, n_cols)

#shape_lines = [((1, 1), (2, 1)),((1, 2), (2, 2)), ((1, 1), (1, 2)),((2, 1), (2, 2))]
shape_lines = [((1.1, 1.2), (2, 0.9)),((0.9, 2.2), (2.1, 2.1)), ((1.1, 1.2), (0.9, 2.2)),((2, 0.9), (2.1, 2.1))]

# Output the list of grid lines
for grid_line in grid_lines:
    print(grid_line)
    

# Calculate distances and find the closest lines
closest_lines = []
for shape_line in shape_lines:
    distances = [line_distance(shape_line, possible_line) for possible_line in grid_lines]
    min_distance_index = np.argmin(distances)
    closest_lines.append(grid_lines[min_distance_index])

# Output the list of grid lines
print("Closest:")
for closest_line in closest_lines:
    print(closest_line)


    
# Plot the lines
plt.figure(figsize=(8, 8))
for grid_line in grid_lines:
    start, end = grid_line
    plt.plot([start[0], end[0]], [start[1], end[1]], 'bo-')
    
for line in shape_lines:
    start, end = line
    plt.plot([start[0], end[0]], [start[1], end[1]], 'ro--') 
    
for line in closest_lines:
    (x1, y1), (x2, y2) = line
    plt.plot([x1, x2], [y1, y2], 'g', label='Closest Possible Lines')

plt.gca().set_aspect('equal', adjustable='box')
plt.title("Grid of Points Connected by Lines")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)
plt.show()




