from __future__ import print_function, absolute_import

import sys

if sys.version_info > (2,):
    xrange = range

import numpy as np
from .alkali_atom_functions import printStateString, C_e, C_h, pi


def htmlLiteratureOutput(v, ref):
    print(
        "<div class='lit'><p>Literature values<p>Radial part of dipole matrix element: %.3f</p>"
        % v
    )
    typeOfSource = "experimental value"
    if ref[0] == 1:
        typeOfSource = "theoretical value"
    print(
        "<p>Source: <a class='link' target='_blank' href='http://dx.doi.org/%s'>%s</a>, %s (%s) </p>"
        % (ref[4], ref[3], typeOfSource, ref[2])
    )
    print("</div>")


def rabiFrequencyWidget(atom, n1, l1, j1, n2, l2, j2, laserPower, laserWaist):
    sol = []
    inputMj = '<p>Rabi frequency $=$ <span id="rabival">0</span><p><form id="polarization" onchange="myFunction()">'
    inputMj += '<p>for driving from <select id="mj" onchange="myFunction()">'
    index = 0
    for mj1 in np.linspace(-j1, j1, int(round(2 * j1 + 1))):
        inputMj += '<option value="%d">m_j = %d/2 ' % (
            index,
            int(round(2.0 * mj1)),
        )
        arr = []
        for q in [-1, 0, 1]:
            if abs(mj1 + q) - 0.1 < j2:
                rabiFreq = atom.getRabiFrequency(
                    n1, l1, j1, mj1, n2, l2, j2, q, laserPower, laserWaist
                ) / (2 * pi)
                arr.append(
                    "$2 \\pi \\times$"
                    + printValueString(rabiFreq, "Hz", decimalPlaces=2)
                )
            else:
                arr.append("not coupled")
        sol.append(arr)
        index += 1
    inputMj += r'</select>\
    <input type="radio" name="colors" id="sigma-" value="0" >$\sigma^-$ | \
<input type="radio" name="colors" id="pi" value="1" checked>$\pi$ |\
<input type="radio" name="colors" id="sigma+" value="2" >$\sigma^+$\
 transition</p></form>'

    script = "<script id='returnscript' type='text/javascript'>"
    script = script + "var rabiFreq =" + str(sol) + "; "
    script += 'function myFunction() {\
    var mj = document.getElementById("mj").value;\
    var p = 0;\
    if (document.getElementById("sigma-").checked){\
        p=0;\
    }\
    if (document.getElementById("pi").checked){\
        p=1;    \
    }\
    if (document.getElementById("sigma+").checked){\
        p=2;    \
    }\
    document.getElementById("rabival").innerHTML = rabiFreq[mj][p] ;\
    MathJax.Hub.Queue(["Typeset",MathJax.Hub,"rabival"]);\
}\
document.getElementById("polarization").addEventListener("click", myFunction);\
myFunction();\
</script>'
    return inputMj + script


def printValueString(value, unit, decimalPlaces=3):
    prefix = ["f", "p", "n", "$\\mu$", "m", "", "k", "M", "G", "T"]
    i = 5
    sg = 1.0
    if value < 0:
        sg = -1.0
    value = abs(value)
    formatString = "%%.%df %%s%%s" % decimalPlaces

    if value > 1000:
        while (value > 1000) and (i < 9):
            value = value * 1.0e-3
            i += 1
        return formatString % (sg * value, prefix[i], unit)
    elif value < 1:
        while (value < 1) and (i > 0):
            value = value * 1.0e3
            i -= 1
        return formatString % (sg * value, prefix[i], unit)
    else:
        return formatString % (sg * value, "", unit)


def plotStarkMap(calc, units=1, xlim=[], ylim=[], filename=""):
    originalState = calc.basisStates[calc.indexOfCoupledState]
    n = originalState[0]
    l = originalState[1]
    j = originalState[2]

    ax = webPlot()

    x = []
    y = []
    yState = []

    ax.xlabel = "E field (V/cm)"

    coeff = 1.0
    ax.ylabel = "Energy/h (GHz)"

    if units == 1:
        # in cm^{-1}
        coeff = 0.03336  # conversion factor from GHz to cm^{-1}
        ax.ylabel = "Energy/(h c) (cm^{-1})"
    if ylim == []:
        ylim = [
            calc.atom.getEnergy(n, l, j) * C_e / C_h * 1e-9 * coeff - 10,
            calc.atom.getEnergy(n, l, j) * C_e / C_h * 1e-9 * coeff + 10,
        ]

    for br in xrange(len(calc.y)):
        for i in xrange(len(calc.y[br])):
            yt = calc.y[br][i] * coeff
            if yt < ylim[1] and ylim[0] < yt:
                x.append(calc.eFieldList[i])
                y.append(yt)
                yState.append(calc.highlight[br][i])

    yState = np.array(yState)
    sortOrder = yState.argsort(kind="heapsort")
    x = np.array(x)
    y = np.array(y)

    x = x[sortOrder]
    y = y[sortOrder]
    yState = yState[sortOrder]

    ct = "|< %s | \\mu > |^2" % printStateString(n, l, j)

    ax.scatter(x / 100.0, y, c=yState, cmin=0, cmax=1, ctitle=ct)

    if xlim == []:
        xlim = [min(x) / 100.0, max(x) / 100.0]

    ax.printPlot(
        xlim=xlim, ylim=ylim, filename=filename, name="starkdiv1", height=600
    )

    return 0


def plotInteractionLevels(calc, xlim=[], ylim=[], filename=""):
    ax = webPlot()
    ax.xlabel = r"R (\mu m)"
    ax.ylabel = r"\Delta E (GHz)"

    if calc.drivingFromState[0] == 0:
        # colouring is based on the contribution of the original pair state here
        ct = r"|< %s %.1f , %s %.1f | \mu > |^2$" % (
            printStateString(calc.n, calc.l, calc.j),
            calc.m1,
            printStateString(calc.nn, calc.ll, calc.jj),
            calc.m1,
        )
    else:
        # colouring is based on the coupling to different states
        ct = r"\Omega_\mu/\Omega"

    x = []
    y = []
    yState = []
    for br in xrange(len(calc.y)):
        for i in xrange(len(calc.y[br])):
            x.append(calc.r[i])
            y.append(calc.y[br][i])
            yState.append(calc.highlight[br][i])

    yState = np.array(yState)
    sortOrder = yState.argsort(kind="heapsort")
    x = np.array(x)
    y = np.array(y)

    x = x[sortOrder]
    y = y[sortOrder]
    yState = yState[sortOrder]

    ax.scatter(x, y, c=yState, cmin=0, cmax=1, ctitle=ct)

    ax.printPlot(xlim=xlim, ylim=ylim, filename=filename, name="levelintdiv")
    return


class webPlot:
    def __init__(self):
        self.traces = []
        self.layout = []
        self.traceNo = 0
        self.xlabel = ""
        self.ylabel = ""
        self.layoutx = ""
        self.layouty = ""
        self.title = ""

    def plot(self, x, y, type, name=""):
        np.set_printoptions(threshold=1e10)
        self.traceNo += 1
        temp = "{ x:" + np.array2string(x, separator=",") + ",\n"
        temp = temp + "y: " + np.array2string(y, separator=",") + ",\n"
        if type == ".":
            temp += "mode: 'markers',\n marker: {size:5},\n"
        elif type == "-":
            temp += "mode: 'lines',\n"
        temp += "name: '%s'" % name
        temp += "}"
        self.traces.append(temp)

    def semilogx(self, x, y, type, name=""):
        self.layoutx = "type:'log' ,\n\
        tickformat :'.1e',\n        "
        self.plot(x, y, type, name)

    def semilogy(self, x, y, type, name=""):
        self.layouty = "type:'log' ,\n\
        tickformat :'.1e',\n        "
        self.plot(x, y, type, name)

    def scatter(self, x, y, c=[], cmin=0, cmax=1, ctitle="", name=""):
        np.set_printoptions(threshold=1e10)
        self.traceNo += 1
        temp = (
            "{ x:"
            + np.array2string(
                x,
                separator=",",
            )
            + ",\n"
        )
        temp = temp + "y: " + np.array2string(y, separator=",") + ",\n"
        temp += "name: '%s',\n" % name
        if c != []:
            temp = temp + " text: " + np.array2string(c, separator=",") + ",\n"
        temp += "mode: 'markers',\n"
        if c != []:
            temp = (
                temp
                + "marker:{\n\
                color:"
                + np.array2string(c, separator=",")
                + ",\n\
                cmin:%f,\n\
                cmax:%f,\n\
                showscale: true,\n\
                colorbar:{\n\
                    title:'"
                % (cmin, cmax)
                + str(ctitle)
                + "',\n\
                },\n\
                size:5\n\
                 },\n"
            )
        else:
            temp = (
                temp
                + "marker:{\n\
                size:5\n\
                 },\n"
            )
        temp += "}"
        self.traces.append(temp)

    def printPlot(
        self,
        name="",
        width=600,
        height=363,
        xlim=[],
        ylim=[],
        filename="",
        scriptName="returnscript",
    ):
        d = ""
        i = 0
        while i < self.traceNo:
            if i != 0:
                d += ","
            d += self.traces[i]
            i += 1
        d = "data=[" + d + "];\n"

        xLimData = ""
        if not xlim == []:
            xLimData = "range: [%.2E,%.2E],\n" % (xlim[0], xlim[1])
        yLimData = ""
        if not ylim == []:
            yLimData = "range: [%.2E,%.2E],\n" % (ylim[0], ylim[1])

        # now layout

        l = (
            "layout = {\n\
    hovermode: 'closest',\n\
    xaxis:{\n\
        zeroline:false,\n\
        "
            + self.layoutx
            + "\
        "
            + xLimData
            + "\
title: '"
            + self.xlabel
            + "',\n\
        ticks: 'inside',\n\
        showline: true\n\
        },\n\
    yaxis:{\n\
        zeroline:false,\n\
        "
            + self.layouty
            + "\
        "
            + yLimData
            + "\
title: '"
            + self.ylabel
            + "',\n\
        ticks: 'inside' ,\n\
        showline: true  \n\
        }\n\
    };\n"
        )

        if filename == "":
            if name == "":
                name = "plotdiv"
            if self.title != "":
                print("<p>" + self.title + "</p>")
            print(
                "<div id='"
                + name
                + "' style='width:%dpx;height:%dpx;'></div>\n" % (width, height)
            )
            print("<script id='" + scriptName + "' type='text/javascript'>\n")
            print("plotarea = document.getElementById('" + name + "');\n")
            print(d)
            print(l)
            print("Plotly.plot(plotarea, data, layout);\n")
            print("</script>\n")
        else:
            f = open(filename, "w")
            if name == "":
                name = "plotdiv"
            if self.title != "":
                f.write("<p>" + self.title + "</p>")
            f.write(
                "<div id='"
                + name
                + "' style='width:%dpx;height:%dpx;'></div>\n" % (width, height)
            )
            f.write("<script id='" + scriptName + "' type='text/javascript'>\n")
            f.write("plotarea = document.getElementById('" + name + "')\n")
            f.write(d)
            f.write(l)
            f.write("Plotly.plot(plotarea, data, layout);\n")
            f.write("</script>\n")
            f.close()
