import section

if __name__ == "__main__":
    R = 20
    r = 10
    z = 0
    c = 28
    n = 100

    tor = section.Tor(z, R, r)
    tor.generate_points(n)

    H = [0, 5, 10, 15, 20, 25, 30]
    section.plane_research(tor, H)
    section.process_lemniscate(tor, R - r, c)
