import cProfile
import create_projections as cp
import pstats

cProfile.run('cp.do_projection("RCP3PD")',"restats")

p = pstats.Stats('restats')
# p.strip_dirs().sort_stats(-1).print_stats()

# p.sort_stats('cumulative').print_stats(10)
# p.sort_stats('time').print_stats(10)