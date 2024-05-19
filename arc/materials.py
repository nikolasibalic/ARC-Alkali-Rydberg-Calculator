import numpy as np
import os
from .alkali_atom_functions import DPATH


class OpticalMaterial(object):
    """
    Abstract class implementing calculation of basic properties for optical
    materials.
    """

    #: Human-friendly name of material
    name = ""
    #: List of .csv files listing refractive index measurements
    #: first column in these files is wavelength (in mu m), the second
    #: refractive index
    sources = []
    # This array is loaded automatically based on sources list
    sourcesN = []
    #: Any notes about measured values
    sourcesComment = []
    #: Array of max and minimal wavelegth pairs [lambdaMin, lambdaMax]
    #: for each of the sources. Automatically loaded from sources list
    sourcesRange = []

    def __init__(self):
        for s in self.sources:
            self.sourcesN.append(
                np.loadtxt(
                    os.path.join(DPATH, "refractive_index_data", s),
                    skiprows=1,
                    delimiter=",",
                    unpack=True,
                )
            )
            self.sourcesRange.append(
                [self.sourcesN[-1][0].min(), self.sourcesN[-1][0].max()]
            )

    def getN(self, *args, **kwargs):
        """
        Refractive index of material
        """
        return "To-do: refractive index"

    def getRho(self):
        return "To-do: density"

    def getElectricConductance(self):
        return "To-do: electric condctance"

    def getThermalConductance(self):
        return "To-do: thermal conductance"


class Air(OpticalMaterial):
    """
    Air as an optical material at normal conditions
    """

    name = "Air (dry, normal conditions)"
    sources = [
        "Mathar-1.3.csv",
        "Mathar-2.8.csv",
        "Mathar-4.35.csv",
        "Mathar-7.5.csv",
    ]
    sourcesComment = ["vacuum", "vacuum", "vacuum", "vacuum"]

    def getN(self, vacuumWavelength=None, *args, **kwargs):
        """

        Assumes temperature: 15 °C, pressure: 101325 Pa
        """
        if vacuumWavelength is not None:
            x = vacuumWavelength
        else:
            raise ValueError("wavelength not specified for refractive index")

        if (x > 0.23) and (x < 1.690):
            return (
                1
                + 0.05792105 / (238.0185 - x ** (-2))
                + 0.00167917 / (57.362 - x ** (-2))
            )
        else:
            for i, rangeN in enumerate(self.sourcesRange):
                if (x > rangeN[0]) and (x < rangeN[1]):
                    return np.interp(
                        x, self.sourcesN[i][0], self.sourcesN[i][1]
                    )
            raise ValueError(
                "No refrative index data available for requested"
                " wavelength %.3f mum" % x
            )


class Sapphire(OpticalMaterial):
    """
    Sapphire as optical material.
    """

    name = "Sapphire"
    # data from: https://refractiveindex.info
    sources = ["Querry-o.csv", "Querry-e.csv"]
    sourcesN = []
    sourcesComment = ["o", "e"]

    def getN(
        self,
        vacuumWavelength=None,
        airWavelength=None,
        axis="ordinary",
        *args,
        **kwargs,
    ):
        """ """

        if vacuumWavelength is not None:
            air = Air()
            x = vacuumWavelength / air.getN(vacuumWavelength=vacuumWavelength)
        elif airWavelength is not None:
            x = airWavelength
        else:
            raise ValueError("wavelength not specified for refractive index")

        if (axis == "ordinary") or (axis == "o"):
            # electric field polarisation perpendicular to cristal axis
            if (x > 0.2) and (x < 5.0):
                return (
                    1
                    + 1.4313493 / (1 - (0.0726631 / x) ** 2)
                    + 0.65054713 / (1 - (0.1193242 / x) ** 2)
                    + 5.3414021 / (1 - (18.028251 / x) ** 2)
                ) ** 0.5
            else:
                for i, rangeN in enumerate(self.sourcesRange):
                    if (
                        (x > rangeN[0])
                        and (x < rangeN[1])
                        and (self.sourcesComment[i] == "o")
                    ):
                        return np.interp(
                            x, self.sourcesN[i][0], self.sourcesN[i][1]
                        )
                raise ValueError(
                    "No refrative index data available for "
                    "requested wavelength %.3f mum" % x
                )

        elif (axis == "extraordinary") or (axis == "e"):
            # electric field polarisation along cristal axis
            if (x > 0.2) or (x < 5.0):
                return (
                    1
                    + 1.5039759 / (1 - (0.0740288 / x) ** 2)
                    + 0.55069141 / (1 - (0.1216529 / x) ** 2)
                    + 6.5927379 / (1 - (20.072248 / x) ** 2)
                ) ** 0.5
            else:
                for i, rangeN in enumerate(self.sourcesRange):
                    if (
                        (x > rangeN[0])
                        and (x < rangeN[1])
                        and (self.sourcesComment[i] == "e")
                    ):
                        return np.interp(
                            x, self.sourcesN[i][0], self.sourcesN[i][1]
                        )
                raise ValueError(
                    "No refrative index data available for "
                    "requested wavelength %.3f mum" % x
                )
        else:
            raise ValueError("Uknown axis")
