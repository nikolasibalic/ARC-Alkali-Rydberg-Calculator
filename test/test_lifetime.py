from arc import Rubidium85


def test_Rb_D1_and_D2_lifetime():
    atom = Rubidium85()
    assert abs(27.679e-9 - atom.getStateLifetime(5, 1, 0.5)) <= 3e-11
    assert abs(26.2348e-9 - atom.getStateLifetime(5, 1, 1.5)) <= 8e-12
