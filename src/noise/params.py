"""Measured wall parameters used in Fig. 4 and Eq. (2). Regenerate with: python src/noise/measure_parameters.py"""
P_ESC = (0.1738 + 0.3496 + 0.1729) / 3     # probability that a flipped cell creates a pair of walls, averaged over the 3 steps
VR, VL = 0.2188, 0.2360                    # speeds of right- and left-moving walls (blocks per step)
