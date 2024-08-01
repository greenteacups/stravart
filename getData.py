import matplotlib.pyplot as plt
from PIL import Image
import csv

# Function to handle mouse clicks
def onclick(event):
    if event.xdata and event.ydata:
        x, y = event.xdata / img_width, event.ydata / img_height
        coords.append((x, y))
        print(f'Clicked at: ({x:.4f}, {y:.4f})')

        # Write the coordinates to the CSV file
        with open('drawing-coords.dat', 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow([x, 1-y])

        # Plot the point
        ax.plot(event.xdata, event.ydata, 'ro')

        # Plot lines between consecutive points
        if len(coords) > 1:
            prev_x, prev_y = coords[-2]
            prev_x *= img_width
            prev_y *= img_height
            ax.plot([prev_x, event.xdata], [prev_y, event.ydata], 'r-')

        fig.canvas.draw()

# Load the image
image_path = 'bluey.png'  # Replace with your image file path
image = Image.open(image_path)
img_width, img_height = image.size

# Initialize the coordinates list
coords = []

# Create or clear the file for exporting coordinates
with open('drawing-coords.dat', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['X', 'Y'])

# Display the image
fig, ax = plt.subplots()
ax.imshow(image)

# Connect the click event handler
fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()

print('Coordinates have been saved to drawing-coords.dat')

