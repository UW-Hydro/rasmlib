"""
rasmlib analysis stats
"""
import numpy as np


def rmse(s, o):
    """
    Root Mean Squared Error
    input:
        s: simulated
        o: observed
    output:
        rmses: root mean squared error
    """
    return np.sqrt(np.mean((s - o) ** 2))


def mae(s, o):
    """
    Mean Absolute Error
    input:
        s: simulated
        o: observed
    output:
        maes: mean absolute error
    """
    return np.mean(abs(s - o))


def bias(s, o):
    """
    Bias
    input:
        s: simulated
        o: observed
    output:
        bias: bias
    """
    return np.mean(s - o)


def nash_sutcliffe(s, o):
    """
    Nash Sutcliffe efficiency coefficient
    input:
        s: simulated
        o: observed
    output:
        ns: Nash Sutcliffe efficient coefficient
    """
    return 1 - sum((s - o) ** 2) / sum((o - np.mean(o)) ** 2)


def percent_error(approx, exact):
    return abs(approx - exact) / abs(exact)
