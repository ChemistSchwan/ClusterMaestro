import matplotlib.pyplot as plt
import numpy as np


def MO_Diagram(eigerFile="eiger.out", plotName="MO-diagram.pdf"):
    """
    Plots the MO-Diagram from the eiger-program, as implemented in turbomole.
    Returns the HOMO-LUMO-Gap

    Parameters
    ----------

    eigerFile: str
        Output-file from the eiger-program.
    plotName: str
        Name of the resulting plot.
    """
    getEnergy = False
    orbitalEng = []
    HOMO, LUMO = -9999, 9999
    with open(eigerFile, "r") as source:
        for line in source:
            if len(line.split()) == 0:
                getEnergy = False
            if getEnergy:
                splitLine = line.split()
                orbitalEng.append(float(splitLine[-2]))
                if splitLine[-6] == "2.000":
                    if float(splitLine[-2]) > HOMO:
                        HOMO = float(splitLine[-2])
                else:
                    if float(splitLine[-2]) < LUMO:
                        LUMO = float(splitLine[-2])
            if "Occupation" in line:
                getEnergy = True
    orbitalEng.sort()
    orbitalEng = np.array(orbitalEng)
    plt.clf()
    fig = plt.figure(figsize=(4, 3))
    fig.tight_layout()
    for i, val in enumerate(orbitalEng):
        if i < len(orbitalEng):
            try:
                if val == orbitalEng[i + 1]:
                    plt.plot([-0.35, -0.15], [val, val], c="black")
            except:
                None
            if val == orbitalEng[i - 1]:
                plt.plot([0.35, 0.15], [val, val], c="black")
            elif val == orbitalEng[i - 1] and val == orbitalEng[i + 1]:
                plt.plot([-0.1, 0.1], [val, val], c="black")
            else:
                plt.plot([-0.1, 0.1], [val, val], c="black")
        else:
            plt.plot([-0.1, 0.1], [val, val], c="black")

    plt.ylabel("Orbital Energy / eV")
    plt.xticks()
    plt.xlim(-0.5, 0.5)
    plt.ylim(-8, 2)
    plt.xticks([], [])
    plt.tight_layout()
    plt.savefig(plotName)
    plt.clf()
    return LUMO - HOMO


def vibSpectrum(
    vibspectrum="vibspectrum",
    plotName="Vibspectrum.pdf",
    waveRange=[0, 3500],
    yRange=None,
):
    """
    Plots a simple idealized vibration spectrum from the aoforce-output file of turbomole.

    Parameters
    ----------

    vibspectrum: str
        Vibspectrum-Output-file from the aoforce-program.
    plotName: str
        Name of the resulting plot.
    waveRange: 1x2 list
        xRange of the resulting plot.
    yRange: 1x2 list or None
        YRange of the resulting plot.
    """
    freqList, intList = [], []
    with open(vibspectrum, "r") as source:
        for line in source:
            if "#" not in line and "spectrum" not in line and "end" not in line:
                splitLine = line.split()
                if len(splitLine) == 6:
                    freqList.append(float(splitLine[2]))
                    intList.append(float(splitLine[3]))
    fig, ax = plt.subplots(figsize=(6, 3))
    markerline, stemlines, baseline = plt.stem(
        freqList, intList, markerfmt=" ", bottom=-0
    )
    markerline.set_markerfacecolor("none")
    stemlines.set_linewidths(1)
    baseline.set_linewidth(1)
    baseline.set_color("black")
    plt.xlabel("frequency / cm^-1")
    plt.ylabel("Intensity / km/mol")
    plt.tight_layout()
    plt.xlim(waveRange[0], waveRange[1])
    if yRange == None:
        plt.ylim(0)  # , 105)
    else:
        plt.ylim(yRange[0], yRange[1])
    plt.savefig(plotName)
    plt.clf()


def exSpectrum(exFile="exspectrum", plotName="exSpectrum.pdf", waveRange=None):
    """
    Plots a simple idealized excitation spectrum from the exspectrum-output file of turbomole.

    Parameters
    ----------

    exFile: str
        Output-file from the escf-program.
    plotName: str
        Name of the resulting plot.
    waveRange: 1x2 list or None

    """
    array = np.genfromtxt(exFile)
    fig, ax = plt.subplots(figsize=(6, 3))
    markerline, stemlines, baseline = plt.stem(
        array[:, 5], array[:, 7], markerfmt=" ", bottom=-0
    )
    markerline.set_markerfacecolor("none")
    stemlines.set_linewidths(1)
    baseline.set_linewidth(1)
    baseline.set_color("black")
    plt.xlabel("frequency / nm")
    plt.ylabel("Intensity / arb. units")
    plt.tight_layout()
    #    plt.xlim(waveRange[0], waveRange[1])
    plt.ylim(0)
    plt.savefig(plotName)
    plt.clf()


def histPlot(data, plotName="hist.pdf", bins=20, xlabel="", ylabel="", mean=False):
    """
    Plots a 1D array as a Histogramm.

    Parameters
    ----------

    data: 1D np.array
        Dataarray with the values.
    bins: int
        How many bars there will be in the plot.
    xlabel: str
        The text for the x-axis.
    ylabel: str
        The text for the y-axis.
    mean: bool
        Does some statistics stuff.
    """
    from scipy.stats import norm

    def mean(numbers):
        return float(sum(numbers)) / max(len(numbers), 1)

    fig, axs = plt.subplots(figsize=(4, 3), tight_layout=True)
    n, nbins, patches = axs.hist(data, bins=bins)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if mean:
        mu, std = norm.fit(data)
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x, mu, std)
        factor = max(n) / max(p)
        p *= factor
        plt.plot(x, p, "k", linewidth=2)
        plt.title("mean: {:.3f}, std: {:.5f}".format(mu, std))
    plt.savefig(plotName)
    plt.clf()
