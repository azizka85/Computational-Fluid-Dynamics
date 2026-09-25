def solve_uniform(
    l, l1, c, c0, c1, r, r0,
    d, u
):
    n = len(u)

    if n > 0:
        u[0] = r0 / c0
        d[0] = d[0] / c0

        if n > 1:
            for i in range(1, n-1):
                ct = c - l * u[i-1]

                u[i] = r / ct
                d[i] = (d[i] - l * d[i-1]) / ct

            ct = c1 - l1 * u[n-2]
            d[n-1] = (d[n-1] - l1 * d[n-2]) / ct

        u[n-1] = d[n-1]

        for i in range(n-2, -1, -1):
            u[i] = d[i] - u[i]*u[i+1]
