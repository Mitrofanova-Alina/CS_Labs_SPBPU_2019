import numpy as np
import matplotlib.pyplot as plt


def lemniscate(c):
    n = 360
    fi = 0
    step = 2 * np.pi / n
    left = []
    right = []
    for i in range(0, n // 4):
        t = np.tan(fi)
        x = c * (t + t ** 3) / (1 + t ** 4)
        y = c * (t - t ** 3) / (1 + t ** 4)
        right.append([x, y])
        fi += step

    for i in range(n // 4, n // 2):
        t = np.tan(fi)
        x = c * (t + t ** 3) / (1 + t ** 4)
        y = c * (t - t ** 3) / (1 + t ** 4)
        left.append([x, y])
        fi += step

    # names = ["left", "right"]
    #
    # fig, ax = plt.subplots(nrows=1, ncols=1, sharey=True)
    # data = [left, right]
    #
    # i = 0
    # plt.axis('equal')
    # for points in data:
    #     x = []
    #     y = []
    #     for elem in points:
    #         x.append(elem[0])
    #         y.append(elem[1])
    #     ax.plot(x, y, label=names[i])
    #     i += 1
    # fig.legend()
    # fig.show()

    return left, right
