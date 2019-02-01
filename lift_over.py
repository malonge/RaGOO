

if __name__ == "__main__":
    import os
    import argparse

    parser = argparse.ArgumentParser(description='Lift-over gff coordinates to from contigs to RaGOO pseudomolecules')
    parser.add_argument("gff", metavar="<genes.gff>", type=str, help="Gff file to be lifted-over")
    parser.add_argument("-g", metavar="100", type=int, default=100, help="Gap size for padding in pseudomolecules (must match what was used for 'ragoo.py'.")


    # Get the command line arguments
    args = parser.parse_args()
    gff_file = args.gff


    # Needs to account for chimeric broken contigs.
    # Really, this should just be integrated into the ragoo.py command
    #    -gff should trigger chimeric breaking to avoid intervals, and should lift-over as it goes

    # Looks like ragoo.py will already update gff files after breaking chimeras.
    #   Next question is, do we then do the final lift-over automatically
    #   Or do it as a stand-alone post processesing script.