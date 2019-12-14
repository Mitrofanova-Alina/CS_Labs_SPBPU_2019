from frechet import Frechet


def test_frechet_open():
    P = [(0, 0), (4, 2), (6, 5), (12, 6), (15, 7), (15, 10), (18, 13)]
    Q = [(1, 1), (2, 5), (7, 7), (8, 12), (13, 14), (15, 16)]

    solver = Frechet(P, Q)
    dist, i, j = solver.frechet_distance()

    print("Frechet distance from P to Q is %f (%d, %d)" % (dist, i, j))


def test_frechet_close():
    P = [(2, 2), (3, 4), (2, 7), (5, 6), (9, 8), (8, 5), (10, 1), (6, 3), (2, 2)]
    Q = [(12, 1), (10, 3), (6, 6), (9, 7), (10, 9), (12, 6), (15, 5), (13, 3), (12, 1)]

    solver = Frechet(P, Q)
    dist, i, j = solver.frechet_distance()

    print("Frechet distance from P to Q is %f (%d, %d)" % (dist, i, j))


if __name__ == "__main__":
    print("Hello, World!")
    test_frechet_open()
    test_frechet_close()
