import numpy as np
import matplotlib.pyplot as plt

#input matrix is a matrix of 256 bins
def chr_score(matrix, res = 10000, radius = 500000, pseudocount_coeff = 30):
    pseudocount = matrix.mean() * pseudocount_coeff
    pixel_radius = int(radius / res)
    scores = []
    for loc_i, loc in enumerate(range(len(matrix))):
        scores.append(point_score(loc, pixel_radius, matrix, pseudocount))
    return scores

# def point_score(locus, radius, matrix, pseudocount):
#     l_edge = max(locus - radius, 0)
#     r_edge = min(locus + radius, len(matrix))
#     l_mask = matrix[l_edge : locus, l_edge : locus] 
#     r_mask = matrix[locus : r_edge, locus : r_edge]
#     center_mask = matrix[l_edge : locus, locus : r_edge]
#     score = (max(l_mask.mean(), r_mask.mean()) +  pseudocount) /(center_mask.mean() + pseudocount)
#     return score

def point_score(locus, radius, matrix, pseudocount):
    l_edge = max(locus - radius, 0)
    r_edge = min(locus + radius, len(matrix))
    
    # Compute slices
    l_mask = matrix[l_edge:locus, l_edge:locus] if locus - l_edge > 0 else None
    r_mask = matrix[locus:r_edge, locus:r_edge] if r_edge - locus > 0 else None
    center_mask = matrix[l_edge:locus, locus:r_edge] if (locus - l_edge > 0 and r_edge - locus > 0) else None

    # Calculate means safely by checking if the slice exists and is non-empty
    l_mean = l_mask.mean() if l_mask is not None and l_mask.size > 0 else 0
    r_mean = r_mask.mean() if r_mask is not None and r_mask.size > 0 else 0
    center_mean = center_mask.mean() if center_mask is not None and center_mask.size > 0 else 0

    # Compute the score, ensuring we don't divide by zero
    score = (max(l_mean, r_mean) + pseudocount) / (center_mean + pseudocount)
    return score

matrix= np.load("./output/cluster0/prediction/npy/chr12_74500000.npy") 
score = chr_score(matrix)

plt.figure(figsize=(8, 6))
plt.plot(score,marker='o',color='black',markersize=2)
plt.title('Insulation score')
plt.savefig("insulation_score.png")