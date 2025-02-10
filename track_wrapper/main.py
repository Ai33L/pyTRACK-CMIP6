from .track_wrapper import *
from .track_stats import *
from .track_wrapper_generic import *
from .composite import *

__all__ = ['Case']

class Case(object):
    def __init__(self, indir, outdir):
        self.indir = indir
        self.outdir = outdir

    def accelerate(self):
        self.speed = self.speed + 2
        print (self.speed)