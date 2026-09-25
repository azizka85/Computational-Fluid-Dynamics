def adjust_time_step(b, t, dt, ots):
    if t >= ots:
        return ots

    if t < ots and t + dt >= ots:
        return ots - t

    return b*dt
