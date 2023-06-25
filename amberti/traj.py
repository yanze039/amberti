from amberti.amber import cpptraj


def extract_conf(conf, prmtop):

    script = [
        f"trajin {conf}",
        # remove the two ligands and keep the rest
        'strip ":1,2"',
        'outtraj ${s}_solvated.pdb onlyframes 1'

    # extract the first ligand
    unstrip
    strip ":2-999999"
    outtraj ${s}_M7G.pdb onlyframes 1

    # extract the second ligand
    unstrip
    strip ":1,3-999999"
    outtraj ${s}_LNA.pdb onlyframes 1
    ]