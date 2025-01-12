import numpy as np 

N_total = 8
N_pheno = 5
y = np.array([4.5, 2.9, 3.9, 3.5, 5.0])
sex = np.array([1, 0, 0, 1, 1])
X = np.zeros((N_pheno, 2), dtype = np.int8)
X[:,0] = sex
X[:,1] = 1 - sex
Z = np.zeros((N_pheno, N_total), dtype = np.int8)
Z[:, (N_total-N_pheno):] = np.eye(N_pheno)
A = np.array([[1.00,	0.00,	0.00,	0.500,	0.000,	0.50,	0.25,	0.250],
                [0.00,	1.00,	0.00,	0.000,	0.500,	0.50,	0.25,	0.250],
                [0.00,	0.00,	1.00,	0.000,	0.500,	0.00,	0.25,	0.500],
                [0.50,	0.00,	0.00,	1.000,	0.000,	0.25,	0.50,	0.125],
                [0.00,	0.50,	0.50,	0.000,	1.000,	0.25,	0.50,	0.375],
                [0.50,	0.50,	0.00,	0.250,	0.250,	1.00,	0.25,	0.500],
                [0.25,	0.25,	0.25,	0.500,	0.500,	0.25,	1.00,	0.250],
                [0.25,	0.25,	0.50,	0.125,	0.375,	0.50,	0.25,	1.000]])

variance_components = np.array([0.25, 0.25]) # sigma_a, sigma_e

np.savetxt('y_vec.txt', y)
np.savetxt('X_mat.txt', X)
np.savetxt('Z_mat.txt', Z)
np.savetxt('A_mat.txt', A)
np.savetxt('variance_components.txt', variance_components)

print('done.')

