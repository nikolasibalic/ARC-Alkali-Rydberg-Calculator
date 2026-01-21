from arc import PairStateInteractions, Potassium


def testPotassium():
    # 66⁢𝑆1/2 states used in Ref. S. Helmrich, A. Arias, S. Whitlock, arXiv:1605.08609 by calling

    calculation = PairStateInteractions(
        Potassium(), 66, 0, 0.5, 66, 0, 0.5, 0.5, 0.5
    )
    assert abs(calculation.getC6perturbatively(0, 0, 5, 30e9) - (-265)) <= 1
