import matplotlib.pyplot as plt
import numpy as np

# Example matrix
matrix = np.array([
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9]
])

# Use the imshow function to create a heatmap
plt.imshow(matrix, cmap='hot', interpolation='nearest')

# Adding a color bar to show the scale
plt.colorbar()

# Optional: Adding titles and labels
plt.title("Matrix Visualization")
plt.xlabel("Column")
plt.ylabel("Row")

# Display the plot
plt.show()
