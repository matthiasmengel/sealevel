# This file is part of SEALEVEL - a tool to estimates future sea-level rise
# constrained by past obervations and long-term sea-level commitment
# Copyright (C) 2016 Matthias Mengel working at PIK Potsdam

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# LICENSE.txt for more details.

import matplotlib.pylab as plt
import numpy as np
import dimarray as da
from mpl_toolkits.axes_grid1 import make_axes_locatable
from . import get_ipcc_data as ipcc
import importlib

importlib.reload(ipcc)

# plot settings
plt.rcParams["xtick.major.pad"] = 10
plt.rcParams["font.size"] = 12
plt.rcParams["lines.markeredgewidth"] = 2
plt.rcParams["legend.fontsize"] = 12
plt.rcParams["figure.figsize"] = 10, 8
plt.rcParams["figure.facecolor"] = "white"


contrib_ids = ["thermexp", "gic", "gis_sid", "gis_smb", "ant_sid", "ant_smb"]
rcpcoldict = {"RCP3PD": "#2256A6", "RCP45": "#73B2E1", "RCP85": "#EE322D"}
rcpnamedict = {"RCP3PD": "RCP26", "RCP45": "RCP45", "RCP85": "RCP85"}


def fig3(projection_data):

    labels = {
        "gic": "Mountain glaciers",
        "thermexp": "Thermal expansion",
        "gis_sid": "Greenland SID",
        "gis_smb": "Greenland SMB",
        "ant_sid": "Antarctica SID",
        "ant_smb": "Antarctica SMB",
    }

    plt.subplots_adjust(
        left=0.1, bottom=0.06, right=1.0, top=0.97, wspace=0.1, hspace=None
    )

    # print "## sl contribtioon, scenario, 5th percentiles in 2100, median in
    # 2100, 95th percentile in 2100"

    plot_period = np.arange(2000, 2101, 1)
    axs = []
    plotlb = "abcdef"
    for i, name in enumerate(contrib_ids[0:6]):
        ax = plt.subplot(3, 2, i + 1)
        # ax = plt.subplot(111)
        divider = make_axes_locatable(ax)
        axy = divider.append_axes("right", size=0.6, pad=0.0, sharey=ax, axisbg="grey")
        axy.axis("off")
        axy.axvspan(0, 10, facecolor="0.5", alpha=0.2, lw=0)
        # print name
        xloc = 10  # ipcc
        oloc = 4  # our estimates
        # single_contribs[name] = {}
        for k, scen in enumerate(["RCP85", "RCP45", "RCP3PD"]):

            contrib = projection_data[scen][name]
            # anomaly to 1986-2005
            contrib = contrib - contrib[1986:2005, :].mean(axis="time")
            # pandas handles NaNs within percentile calculus
            contrib = contrib.T.to_pandas()
            percentiles = da.from_pandas(contrib.quantile([0.05, 0.5, 0.95]))[
                :, plot_period
            ]
            ax.fill_between(
                plot_period,
                percentiles[0.05] * 1e3,
                percentiles[0.95] * 1e3,
                color=rcpcoldict[scen],
                alpha=0.4,
                lw=0.5,
            )
            ax.plot(
                plot_period,
                percentiles[0.5] * 1e3,
                lw=3,
                color=rcpcoldict[scen],
                alpha=1.0,
                label=rcpnamedict[scen],
            )

            oloc -= 1.0
            xloc -= 1.0
            med = percentiles[0.5][np.arange(2081, 2100, 1)].mean() * 1e3
            low = percentiles[0.05][np.arange(2081, 2100, 1)].mean() * 1e3
            upp = percentiles[0.95][np.arange(2081, 2100, 1)].mean() * 1e3
            axy.plot([oloc - 1, oloc + 1], [med, med], lw=3, color=rcpcoldict[scen])

            axy.fill_between(
                [oloc - 1, oloc + 1],
                [low, low],
                [upp, upp],
                color=rcpcoldict[scen],
                alpha=0.4,
                lw=0.0,
            )

            low, med, upp = ipcc.ipcc_contrib_estimates[name][scen]
            axy.fill_between(
                [xloc - 1, xloc + 1],
                [low, low],
                [upp, upp],
                color=rcpcoldict[scen],
                alpha=0.4,
                lw=0.0,
            )
            axy.plot(
                [xloc - 1, xloc + 1],
                [med, med],
                color=rcpcoldict[scen],
                lw=3,
                alpha=1.0,
            )
            axy.set_xlim(-1, 11)

        ax.text(
            0.05,
            0.9,
            plotlb[i],
            transform=ax.transAxes,
            fontdict={"family": "sans-serif", "weight": "bold", "size": 16},
        )
        ax.text(
            0.05,
            0.8,
            labels[name],
            transform=ax.transAxes,
            fontdict={"family": "sans-serif", "weight": "bold"},
        )
        if i == 0:
            axy.text(
                1.5,
                440,
                "M16",
                rotation="vertical",
                horizontalalignment="center",
                verticalalignment="center",
            )
            axy.text(
                7.5,
                440,
                "IPCC",
                rotation="vertical",
                horizontalalignment="center",
                verticalalignment="center",
            )

        l1 = ax.legend(ncol=1, loc="center left")
        l1.draw_frame(0)
        for l in l1.get_lines():
            l.set_alpha(1)

        axs.append(ax)
        ax.set_xlim(plot_period[0], plot_period[-1])

    for ax in axs[-2:]:
        ax.set_xlabel("Time in years")
    for ax in [axs[0], axs[2], axs[4]]:
        ax.set_ylabel("Sea level in mm")

    axs[3].set_ylim(0, 500)


def fig4(projection_data):

    plot_period = np.arange(2000, 2101, 1)

    plt.subplots_adjust(
        left=0.1, bottom=0.1, right=1.0, top=0.97, wspace=0.1, hspace=None
    )

    ax6 = plt.subplot(111)

    # ipcc ranges
    divider = make_axes_locatable(ax6)
    axy = divider.append_axes("right", size=0.8, pad=0.0, sharey=ax6)
    axy.axis("off")
    axy.axvspan(0, 10, facecolor="0.5", alpha=0.2, lw=0)
    total_contribs = {}
    xloc = 10
    oloc = 4

    for k, scen in enumerate(["RCP85", "RCP45", "RCP3PD"]):
        total_slr = da.zeros_like(projection_data[scen]["thermexp"])

        for i, name in enumerate(contrib_ids):
            # sum up all contributions
            single_contrib = projection_data[scen][name]
            total_slr += single_contrib

        contrib = total_slr - total_slr[1986:2005, :].mean(axis=0)
        upper_perc = da.DimArray(
            np.percentile(contrib, 95, axis=1), axes=contrib.time, dims="time"
        )
        lower_perc = da.DimArray(
            np.percentile(contrib, 5, axis=1), axes=contrib.time, dims="time"
        )
        median = da.DimArray(
            np.percentile(contrib, 50, axis=1), axes=contrib.time, dims="time"
        )

        h = ax6.fill_between(
            plot_period,
            lower_perc[plot_period] * 1e3,
            upper_perc[plot_period] * 1e3,
            color=rcpcoldict[scen],
            alpha=0.4,
            lw=0.5,
        )
        ax6.plot(
            plot_period,
            median[plot_period] * 1e3,
            lw=3,
            color=rcpcoldict[scen],
            alpha=1.0,
            label=rcpnamedict[scen],
        )

        total_contribs[scen] = np.array(
            [median[2100], lower_perc[2100], upper_perc[2100]]
        )

        oloc -= 1.0
        xloc -= 1.0
        low = lower_perc[2081:2100].mean() * 1.0e3
        med = median[2081:2100].mean() * 1.0e3
        upp = upper_perc[2081:2100].mean() * 1.0e3
        axy.plot([oloc - 1, oloc + 1], [med, med], lw=3, color=rcpcoldict[scen])
        axy.fill_between(
            [oloc - 1, oloc + 1],
            [low, low],
            [upp, upp],
            color=rcpcoldict[scen],
            alpha=0.4,
            lw=0.0,
        )
        low, med, upp = ipcc.get_ipcc_range(scen, "mean_slr_2081_2100")
        axy.fill_between(
            [xloc - 1, xloc + 1],
            [low, low],
            [upp, upp],
            color=rcpcoldict[scen],
            alpha=0.4,
            lw=0.0,
        )
        axy.plot(
            [xloc - 1, xloc + 1], [med, med], color=rcpcoldict[scen], lw=3, alpha=1.0
        )
        axy.set_xlim(-1, 11)

    axy.text(
        1.5,
        1200,
        "M16",
        rotation="vertical",
        horizontalalignment="center",
        verticalalignment="center",
    )
    axy.text(
        7.5,
        1200,
        "IPCC",
        rotation="vertical",
        horizontalalignment="center",
        verticalalignment="center",
    )

    ax6.set_xlim(plot_period[0], plot_period[-1])
    ax6.set_xlabel("Time in years")
    ax6.set_ylabel("Sea level in mm")

    l1 = ax6.legend(ncol=1, loc="center left")
    l1.draw_frame(0)
    for l in l1.get_lines():
        l.set_alpha(1)
