import matplotlib.pyplot as plt
import numpy as np
import meshio

# # Your data
data = np.array([
    [0, 0, 0],
    [1, 0, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [1, 1, 1],
    [0, 1, 1],
    [0.5, 0, 0],
    [1, 0.5, 0],
    [0.5, 1, 0],
    [0, 0.5, 0],
    [0, 0, 0.5],
    [1, 0, 0.5],
    [0.5, 0, 1],
    [1, 1, 0.5],
    [1, 0.5, 1],
    [0, 1, 0.5],
    [0.5, 1, 1],
    [0, 0.5, 1],
    [0.5, 0.5, 0],
    [0.5, 0, 0.5],
    [1, 0.5, 0.5],
    [0.5, 1, 0.5],
    [0, 0.5, 0.5],
    [0.5, 0.5, 1],
    [0.5, 0.5, 0.5]
])

# # Extract x, y, z coordinates
# x = data[:, 0].flatten()
# y = data[:, 1].flatten()
# z = data[:, 2].flatten()

# # Create the 3D scatter plot
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(x, y, z, s=50)

# # Set labels for the axes
# ax.set_xlabel('X-axis')
# ax.set_ylabel('Y-axis')
# ax.set_zlabel('Z-axis')

# for i, txt in enumerate(range(0, len(x))):
#     ax.text(x[i], y[i], z[i], str(txt), color='red')

# # Show the plot
# plt.show()

hexa27_elem = [
            1, 13, 5, 9, 22, 15, 2, 14, 6,  
            12, 15, 20, 21, 27, 26, 10, 23, 17,
            4, 18, 8, 11, 24, 19, 3, 16, 7
]

hexa27_elem = [x - 1 for x in hexa27_elem]

hexa27_nodes = np.array(data[hexa27_elem])
shape_title = "hexahedron27"

hexa27_mesh = meshio.Mesh(hexa27_nodes, {shape_title: np.array([hexa27_elem])})

print(hexa27_mesh.points)
print(hexa27_mesh.cells)

meshio.write("meshiotest.vtk", hexa27_mesh)


