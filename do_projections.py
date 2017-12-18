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

""" Minimal way to do sea level projectoons using the sealevel.projection.project
    function.
"""


import settings
reload(settings)
import sealevel.projection

if settings.probablistic_climate:
    import sealevel.get_magicc_gmt_data as mag
    reload(mag)
else:
    import sealevel.get_ipcc_data

if __name__ == "__main__":

    for scen in settings.scenarios:

        print "scenario", scen

        if settings.probablistic_climate:
            # 600 member ensemble
            gmt = mag.magicc_gmt[scen]
        else:
            # single timeseries from IPCC AR5, for illustration and testing
            gmt = sealevel.get_ipcc_data.tas_data[scen]

        sealevel.projection.project_slr(scen, gmt, settings)