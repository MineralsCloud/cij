def main():

    import sys
    from glob import glob
    from cij.plot.quick import plot_table

    for arg in sys.argv[1:]:
        for fname in glob(arg):
            plot_table(fname)