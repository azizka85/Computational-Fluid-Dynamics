def max_abs_diff(u, u1):
    n = len(u)
    diff = 0

    for i in range(n):
        diff = max(diff, abs(u1[i] - u[i]))

    return diff
