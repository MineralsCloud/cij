'''The plotter that helps you to plot physical properties inside ``qha-cij``'s
calulator, without the worries of dealing with boundaries and units.
'''

from cij.core.calculator import Calculator


class Plotter:
    '''The plotter for QHA's calculator
    '''
    def __init__(self, calculator: Calculator):
        self.calculator = calculator
