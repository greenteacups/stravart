import matplotlib.pyplot as plt
from PIL import Image
import csv
import math

# Function to handle mouse clicks
def onclick(event):
    if event.xdata and event.ydata:
        x, y = event.xdata / img_width, event.ydata / img_height
        use_existing = False

        for (px, py) in all_coords:
            if math.sqrt((x - px)**2 + (y - py)**2) < 0.02:
                x, y = px, py
                use_existing = True
                print(f'Using existing point: ({x:.4f}, {y:.4f})')
                break

        coords.append((x, y))
        if not use_existing:
            print(f'Clicked at: ({x:.4f}, {y:.4f})')

        # Plot the point
        ax.plot(event.xdata, event.ydata, 'ro')

        # Plot lines between consecutive points within the same polygon
        if len(coords) % 4 != 0 and len(coords) > 1:
            prev_x, prev_y = coords[-2]
            prev_x *= img_width
            prev_y *= img_height
            ax.plot([prev_x, event.xdata], [prev_y, event.ydata], 'r-', linewidth=2)

        # Connect every 4 points and fill the polygon
        if len(coords) % 4 == 0:
            polygon_points = coords[-4:]
            polygon_points = [(px * img_width, py * img_height) for px, py in polygon_points]
            polygon = plt.Polygon(polygon_points, closed=True, fill=True, edgecolor='r', facecolor='r', alpha=0.3, linewidth=2)
            ax.add_patch(polygon)
            all_coords.extend(coords)  # Add the current polygon points to all_coords

            # Write the polygon coordinates to the CSV file
            with open('road-coords.dat', 'a', newline='') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow([f'Polygon {len(all_coords)//4}'])
                for px, py in polygon_points:
                    csvwriter.writerow([px / img_width, 1 - (py / img_height)])

            coords.clear()  # Reset coordinates to start a new polygon

        fig.canvas.draw()

# Load the image
image_path = 'manhattan.png'  # Replace with your image file path
image = Image.open(image_path)
img_width, img_height = image.size

# Initialize the coordinates lists
coords = []
all_coords = []

# Create or clear the file for exporting coordinates
with open('road-coords.dat', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['Polygon', 'X', 'Y'])

# Display the image
fig, ax = plt.subplots()
ax.imshow(image)

# Connect the click event handler
fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()

print('Coordinates have been saved to road-coords.dat')

