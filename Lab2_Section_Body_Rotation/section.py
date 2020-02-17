from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab
import matplotlib.pyplot as plt
import frechet


def lemniscate(c):
    n = 360
    fi = 0
    step = 2 * np.pi / n
    left, right = [], []
    for i in range(0, n // 4):
        p = np.tan(fi)
        x = c * (p + p ** 3) / (1 + p ** 4)
        y = c * (p - p ** 3) / (1 + p ** 4)
        right.append([x, y])
        fi += step

    for i in range(n // 4, n // 2):
        p = np.tan(fi)
        x = c * (p + p ** 3) / (1 + p ** 4)
        y = c * (p - p ** 3) / (1 + p ** 4)
        left.append([x, y])
        fi += step

    return left, right


class TorPoint:
    def __init__(self, z, r1, r2):
        self.z_ = z
        self.r1_ = r1
        self.r2_ = r2


class Point:
    def __init__(self, x, y, z):
        self.x_ = x
        self.y_ = y
        self.z_ = z


class Tor:
    def __init__(self, z, r_rotate, r_circle):
        self.z_ = z
        # distance from the center of the forming circle to the axis of rotation
        self.R_ = r_rotate
        # radius of the forming circle
        self.r_ = r_circle
        self.points_ = []

    def generate_points(self, num_points):
        self.points_ = []
        step = np.pi / num_points
        fi = - np.pi / 2

        for i in range(0, num_points):
            dr = self.r_ * np.cos(fi)
            point = TorPoint(self.r_ * np.sin(fi), self.R_ - dr, self.R_ + dr)
            self.points_.append(point)
            fi += step

    def get_tor_points(self):
        return self.points_


def points_on_plane(tor_point, x_plane):
    z = tor_point.z_
    r1 = tor_point.r1_
    r2 = tor_point.r2_

    point1 = Point(x_plane, np.sqrt(r1 ** 2 - x_plane ** 2), z)
    point2 = Point(x_plane, np.sqrt(r2 ** 2 - x_plane ** 2), z)
    point3 = Point(x_plane, - np.sqrt(r1 ** 2 - x_plane ** 2), z)
    point4 = Point(x_plane, - np.sqrt(r2 ** 2 - x_plane ** 2), z)

    return point1, point2, point3, point4


def find_intersection_tor_plane(tor, plane):
    tor_points = tor.get_tor_points()

    left1, left2, left3, left4 = [], [], [], []
    right1, right2, right3, right4 = [], [], [], []

    left_jump_1_4_flag, left_jump_2_3_flag = False, False
    right_jump_1_4_flag, right_jump_2_3_flag = False, False

    for point in tor_points:
        tmp = points_on_plane(point, plane)

        if not np.isnan(tmp[0].y_):
            if not right_jump_1_4_flag:
                right1.append(tmp[0])
            else:
                right4.append(tmp[0])
        else:
            right_jump_1_4_flag = True

        if not np.isnan(tmp[1].y_):
            if not right_jump_2_3_flag:
                right2.append(tmp[1])
            else:
                right3.append(tmp[1])
        else:
            right_jump_2_3_flag = True

        if not np.isnan(tmp[2].y_):
            if not left_jump_1_4_flag:
                left1.append(tmp[2])
            else:
                left4.append(tmp[2])
        else:
            left_jump_1_4_flag = True

        if not np.isnan(tmp[3].y_):
            if not left_jump_2_3_flag:
                left2.append(tmp[3])
            else:
                left3.append(tmp[3])
        else:
            left_jump_2_3_flag = True

    right1.reverse()
    right4.reverse()
    left1.reverse()
    left4.reverse()

    right = right2 + right3 + right4 + right1
    right.reverse()
    left = left2 + left3 + left4 + left1

    return left, right


def plot_3d(plane_cut, filename):
    fig = pylab.figure()
    axes = Axes3D(fig)
    number = len(plane_cut)
    cmap = plt.get_cmap('gnuplot')
    colors = [cmap(i) for i in np.linspace(0, 1, number)]

    i = 0
    for points_arr in plane_cut:
        for points in points_arr:
            x, y, z = [], [], []
            for elem in points:
                x.append(elem.x_)
                y.append(elem.y_)
                z.append(elem.z_)

            axes.plot(x, y, z, ".", color=colors[i])

        i += 1

    fig.savefig(filename, dpi=300, format='png', bbox_inches='tight')
    fig.show()
    plt.close(fig)


def draw_cut(points_arr, name, filename):
    plt.axis('scaled')
    fig, ax = plt.subplots(nrows=1, ncols=1, sharey='all')
    ax.set_xlabel("y")
    ax.set_ylabel("z")
    ax.set_title(name)

    for points in points_arr:
        x, y, z = [], [], []
        for elem in points:
            x.append(elem.x_)
            y.append(elem.y_)
            z.append(elem.z_)
        plt.axis('equal')
        plt.grid(True)
        ax.plot(y, z, ".", color="purple")

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])
    fig.savefig(filename, dpi=300, format='png', bbox_inches='tight')
    fig.show()
    plt.close(fig)


def plane_research(tor, planes):
    data = []
    for plane in planes:
        left, right = find_intersection_tor_plane(tor, plane)
        draw_cut([left, right], "H = %i" % plane, "tor_section_H=%i.png" % plane)
        data.append([left, right])

    plot_3d(data, "tor_all_sections.png")


def draw_line_on_plot(line, ax, color):
    data = [line]
    for i in range(0, len(data)):
        x, y = [], []
        for elem in data[i]:
            x.append(elem[0])
            y.append(elem[1])
        plt.axis('equal')
        plt.grid(True)
        ax.plot(x, y, ":", label="lemniscate", color=color)


def process_lemniscate(tor, plane, c):
    left_lemn, right_lemn = lemniscate(c)
    left_sec, right_sec = find_intersection_tor_plane(tor, plane)
    points_arr = [left_sec, right_sec]

    fig, ax = plt.subplots(nrows=1, ncols=1, sharey='all')
    ax.set_xlabel("y")
    ax.set_ylabel("z")
    ax.set_title("Сomparison of leminiscate and tor sections")
    plt.axis('equal')
    plt.grid(True)

    data_arr = []
    draw_flag = True
    for points in points_arr:
        x, y, z = [], [], []
        tmp = []
        for elem in points:
            x.append(elem.x_)
            y.append(elem.y_)
            z.append(elem.z_)
            tmp.append([elem.y_, elem.z_])
        data_arr.append(tmp)
        if draw_flag:
            ax.plot(y, z, "-", color="red", label="dist")
            ax.plot(y, z, ".", color="orange", label="section")
            draw_flag = False
        else:
            ax.plot(y, z, ".", color="orange")

    draw_line_on_plot(left_lemn + right_lemn, ax, "green")

    solver1 = frechet.Frechet(right_lemn, data_arr[1])
    solver2 = frechet.Frechet(left_lemn, data_arr[0])

    dist1, i1, j1 = solver1.frechet_distance()
    dist2, i2, j2 = solver2.frechet_distance()
    print("dist1 = ", dist1, "i1 = ", i1, "j1 = ", j1)
    print("dist2 = ", dist2, "i2 = ", i2, "j2 = ", j2)

    fig.legend(loc='upper left')
    fig.savefig("comparison.png", dpi=300, format='png', bbox_inches='tight')

    fig.show()
    plt.close(fig)
