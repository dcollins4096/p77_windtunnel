
from dtools.starter1 import *

import dtools.taxi.taxi as taxi
reload(taxi)


car = taxi.load('p77_a1')
car.plot()
