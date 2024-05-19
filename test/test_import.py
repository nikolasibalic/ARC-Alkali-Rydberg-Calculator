from arc import (
    Cesium,
    Rubidium85,
    Rubidium87,
    Potassium39,
    Potassium40,
    Potassium41,
    Lithium6,
    Lithium7,
    Strontium88,
    Ytterbium174,
    Calcium40,
    Sodium,
)


def test_import_cs():
    _ = Cesium()


def test_import_rb():
    _ = Rubidium85()
    _ = Rubidium87()


def test_import_li():
    _ = Lithium6()
    _ = Lithium7


def test_import_K():
    _ = Potassium39()
    _ = Potassium40()
    _ = Potassium41()


def test_import_Na():
    _ = Sodium()


def test_import_Sr():
    _ = Strontium88()


def test_import_Ca():
    _ = Calcium40()


def test_import_Yt():
    _ = Ytterbium174()
