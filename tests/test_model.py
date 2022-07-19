import math
import pytest


def logit_p_left(u_left, u_right):
    p = math.exp(u_left) / (math.exp(u_left) + math.exp(u_right))
    return p / (1 - p)


def simplified_logit_p_left(u_left, u_right):
    return 1 / math.exp(u_right - u_left)


@pytest.mark.parametrize("u_left,u_right", [
    (-10, -20),
    (-20, -10),
    (0, 0),
    (5, 5),
    (10, -5),
    (-10, 5)
])
def test_simplified_logit_equals_normal_logit(u_left, u_right):
    logit = logit_p_left(u_left, u_right)
    simplified_logit = simplified_logit_p_left(u_left, u_right)
    assert abs(logit - simplified_logit) == pytest.approx(0, abs=1e-3)
