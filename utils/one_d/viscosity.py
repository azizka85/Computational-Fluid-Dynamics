def linear_profile(z, ht, nus, nub):
    if z > ht:
        return nub

    return (nub*z + nus*(ht - z)) / ht
