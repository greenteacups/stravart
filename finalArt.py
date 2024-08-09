import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

def read_polygon_from_file(file_path):
    coords = []
    with open(file_path, 'r') as file:
        for line in file:
            x, y = map(float, line.strip().split(','))
            coords.append((x, y))
    return coords

# Load the image
image_path = 'manhattan.png'
image = Image.open(image_path)

# Read the polygon coordinates from the file
file_path_polygon = 'merged_polygon.dat'
polygon_coords = read_polygon_from_file(file_path_polygon)

# Convert the polygon coordinates to separate lists of x and y values
x_coords, y_coords = zip(*polygon_coords)

# Plot the image and the polygon
plt.figure(figsize=(8, 8))
plt.imshow(image, extent=[0, 1, 0, 1])  # Assuming the image should be scaled to fit the [0, 1] range
plt.plot(x_coords, y_coords, 'r-', linewidth=3)  # Plot the polygon in red

# Optionally, close the polygon by connecting the last point to the first
plt.plot([x_coords[-1], x_coords[0]], [y_coords[-1], y_coords[0]], 'r-', linewidth=3)

plt.title("Polygon Overlay on Manhattan Image")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)
plt.show()

