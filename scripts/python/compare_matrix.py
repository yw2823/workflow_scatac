import numpy as np
import matplotlib.pyplot as plt

# Load the matrices from .npy files
matrix1 = np.load("./output/cluster8/prediction/npy/chr12_74500000.npy")
matrix2 = np.load("./output/cluster0/prediction/npy/chr12_74500000.npy")

# Make sure the matrices have the same shape before subtracting
if matrix1.shape != matrix2.shape:
    raise ValueError("The matrices must have the same dimensions to be subtracted.")

#add psuedocount
pseudocount = 1e-6

# Compute the relative difference
result_relative = (matrix1 - matrix2) / (matrix2 + pseudocount)

# Subtract matrix2 from matrix1
result = matrix1 - matrix2

# Plot the resulting matrix
#fig, ax = plt.subplots(figsize=(8, 6))
plt.figure(figsize=(8, 6))
plt.imshow(result_relative, cmap='RdBu_r')  # you can choose any colormap you prefer
#im= ax.imshow(result_relative, cmap='RdBu_r')  # you can choose any colormap you prefer
plt.colorbar(label='Difference')
plt.title("Difference between prediction")
plt.xlabel("Genomic position")
plt.ylabel("Genomic position")
plt.savefig("matrix_difference_log.png")