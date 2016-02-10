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

import cProfile
import create_projections as cp
import pstats

cProfile.run('cp.do_projection("RCP3PD")', "restats")

p = pstats.Stats('restats')
# p.strip_dirs().sort_stats(-1).print_stats()

# p.sort_stats('cumulative').print_stats(10)
# p.sort_stats('time').print_stats(10)
