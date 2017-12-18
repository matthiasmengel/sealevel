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
import sealevel.projection as pr
reload(pr)
import sealevel.get_magicc_gmt_data as mag
reload(mag)

if __name__ == "__main__":

    projection_data = {}

    for scen in settings.scenarios:

        print "scenario", scen

        gmt = mag.magicc_gmt[scen]
        pr.project_slr(scen, gmt, settings)