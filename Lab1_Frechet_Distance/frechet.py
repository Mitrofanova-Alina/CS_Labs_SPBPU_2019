import numpy as np
import matplotlib.pyplot as plt


class Frechet:
    def __init__(self, P, Q):
        self.P = np.array(P)
        self.Q = np.array(Q)

        self.p = len(P)
        self.q = len(Q)

        self.ca = np.full((self.p, self.q), -1.0)


    def frechet_distance(self):
        dist, i, j = self.c(self.p - 1, self.q - 1)
        self.plot(i, j, dist)
        return dist, i, j

    def c(self, i, j):
        n_i = i
        n_j = j
        d = np.linalg.norm(self.P[i] - self.Q[j])
        if self.ca[i][j] > -1:
            return self.ca[i][j], n_i, n_j
        elif i == 0 and j == 0:
            self.ca[i][j] = d
        elif i > 0 and j == 0:
            self.ca[i][j], n_i, n_j = max(
                self.c(i - 1, 0), (d, i, 0)
            )
        elif i == 0 and j > 0:
            self.ca[i][j], n_i, n_j = max(
                self.c(0, j - 1), (d, 0, j)
            )
        elif i > 0 and j > 0:
            self.ca[i][j], n_i, n_j = max(min(
                self.c(i - 1, j), self.c(i - 1, j - 1), self.c(i, j - 1)),
                (d, i, j)
            )
        else:
            self.ca[i][j] = float('inf')

        return self.ca[i][j], n_i, n_j


    def plot(self, i, j, d):
        plt.figure()
        plt.plot(self.P[:, 0], self.P[:, 1], color='blue')
        plt.plot(self.Q[:, 0], self.Q[:, 1], color='orange')
        plt.plot([self.P[i][0], self.Q[j][0]], [self.P[i][1], self.Q[j][1]], color='red')
        plt.legend(['P', 'Q', 'Frechet distance = %.3f' % d])
        plt.title('Frechet distance')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.grid(True)
        plt.savefig("Frechet_dist%d.png"%i, dpi=500, format='png')
        plt.show()