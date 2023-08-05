import rdkit
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolAlign
from rdkit.Chem.rdForceFieldHelpers import UFFGetMoleculeForceField
# BUG: only when this import exists, UFFGetMoleculeForceField can be registered.
from rdkit.Chem import ChemicalForceFields
import numpy as np

print("RDKit Version: ", rdkit.__version__)


def calc_dist(pos1, pos2, offset=0):
    return np.linalg.norm(np.array(pos1)+offset- np.array(pos2))


def get_chiral_tag(mol, from_structure=False):
    if from_structure:
        Chem.AssignAtomChiralTagsFromStructure(mol)
    chiral = Chem.FindMolChiralCenters(mol)
    chiral_dict = {x[0]:x[1] for x in chiral}
    return chiral_dict


def find_common_core(target,
                     query):
    # we expect the target doesn't have Hs, since H's can slow down the MCS search
    target = Chem.RemoveHs(target)
    mcs = rdFMCS.FindMCS(
                 [target, query], 
                 ringMatchesRingOnly=True
                    )
    m1 = Chem.MolFromSmarts(mcs.smartsString)
    query_Match = query.GetSubstructMatch(m1)

    # now we remove mapping that are ring atoms in query but not in MCS.
    is_ring = np.ones(m1.GetNumAtoms(), dtype=bool)
    not_ring_atoms= []
    for ii in range(m1.GetNumAtoms()):
        atom_i = m1.GetAtomWithIdx(ii)
        atom_q = query.GetAtomWithIdx(query_Match[ii])
        atom_i_is_in_ring = atom_i.IsInRing()
        atom_q_is_in_ring = atom_q.IsInRing()
        if atom_i_is_in_ring != atom_q_is_in_ring:
            is_ring[ii] = 0
            not_ring_atoms.append(ii)
    
    # remove the not ring atoms from mcs mol
    rw_m1 = Chem.RWMol(m1)
    not_ring_atoms.sort(reverse=True)
    for atom in not_ring_atoms:
        rw_m1.RemoveAtom(atom)
    m1 = rw_m1.GetMol()
    
    # after removing the not ring atoms, mcs are separated into fragments.
    fragments = Chem.GetMolFrags(m1, asMols=True)
    
    fragments_with_atom_n = []
    for frag in fragments:
        fragments_with_atom_n.append((frag, frag.GetNumAtoms()))
    fragments_with_atom_n = sorted(fragments_with_atom_n, key=lambda x: x[1], reverse=True)

    # pick the largest fragmentas the MCS
    target_Match = target.GetSubstructMatch(fragments_with_atom_n[0][0])
    query_Match = query.GetSubstructMatch(fragments_with_atom_n[0][0])

    mcs_atom_mapping = list(zip(target_Match, query_Match))
    
    chiral_target = get_chiral_tag(target, from_structure=False)
    chiral_query = get_chiral_tag(query, from_structure=True)

    for atom_pair in mcs_atom_mapping:
        if not atom_pair[0] in chiral_target or not atom_pair[1] in chiral_query:
            continue
        chiral_t = chiral_target[atom_pair[0]]
        chiral_q = chiral_query[atom_pair[1]]
        assert chiral_t == chiral_q, f"Atom target@{atom_pair[0]}:{chiral_t} and "\
                       f"query@{atom_pair[1]}:{chiral_q} have different chiral tags."
        
    return mcs_atom_mapping

    
def ConstrainedEmbed(mol, core, atom_mapping, force_constant=100., randomseed=-1, getForceField=UFFGetMoleculeForceField):
    
    """ generates an embedding of a molecule where part of the molecule
        is constrained to have particular coordinates

        Arguments
          - mol: the molecule to embed
          - core: the molecule to use as a source of constraints
          - randomSeed: (optional) seed for the random number generator
    """
    coordMap = {}
    coreConf = core.GetConformer()

    for i_t, i_q in atom_mapping:
        coordMap[i_t] = coreConf.GetAtomPosition(i_q)
        
    ci = rdDistGeom.EmbedMolecule(mol, 
                                  coordMap=coordMap, 
                                  clearConfs=True,
                                  enforceChirality=True,
                                  randomSeed=randomseed
                                 )
    if ci < 0:
        raise ValueError('Could not embed molecule.')

    # rotate the embedded conformation onto the core:
    rms = rdMolAlign.AlignMol(mol, core, atomMap=atom_mapping)
    ff = getForceField(mol, confId=0)
    for i_t, i_q in atom_mapping:
        p = coreConf.GetAtomPosition(i_q)
        pIdx = ff.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
        ff.AddDistanceConstraint(pIdx, i_t, 0, 0, force_constant)
        
    ff.Initialize()
    n = 6
    more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
    while more and n:
        more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
        n -= 1
    # realign
    rms = rdMolAlign.AlignMol(mol, core, atomMap=atom_mapping)
    mol.SetProp('EmbedRMS', str(rms))
    return mol


def PositionRestrainedEnergyMinimization(mol, atom_list, force_constant=100., max_attempts=5, getForceField=UFFGetMoleculeForceField):
    ff = getForceField(mol, confId=0)
    conf = mol.GetConformer()
    for i in atom_list:
        p = conf.GetAtomPosition(i)
        pIdx = ff.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
        ff.AddDistanceConstraint(pIdx, i, 0, 0, force_constant)

    ff.Initialize()
    n = max_attempts
    more = ff.Minimize(energyTol=1e-5, forceTol=1e-4)
    
    # more restrict minimization
    while more and n:
        more = ff.Minimize(energyTol=1e-5, forceTol=1e-4)
        n -= 1
    return mol


def alignHs(target, query, atom_mapping, dist_tol = 1.0, debug = False):
    conf_target = target.GetConformer()
    conf_query = query.GetConformer()
    assert query.GetNumConformers(), "No conformers are found in query molecule."
    assert target.GetNumConformers(), "No conformers are found in target molecule."
    hydrogen_mapping = []
    for hvy_atom_pair_idx in atom_mapping:
        if debug: print("processing: ", hvy_atom_pair_idx)
        hyv_atom_target_idx, hyv_atom_query_idx = hvy_atom_pair_idx
        hyv_atom_target_pos = conf_target.GetAtomPosition(hyv_atom_target_idx)
        hyv_atom_query_pos = conf_query.GetAtomPosition(hyv_atom_query_idx)
        dist =  calc_dist(hyv_atom_target_pos, hyv_atom_query_pos)
        
        if debug: print(f"distance of pairs: {dist:.2f}")
        if dist > dist_tol:
            writer = Chem.SDWriter('embed.target.H_align.fail.sdf')
            for cid in range(target.GetNumConformers()):
                writer.write(target, confId=cid)
            raise RuntimeError(f"alignment falied. distance {dist:.2f} exceeds tolerance {dist_tol:.2f}"\
                            f" for atom pair target:{hyv_atom_target_idx} and query:{hyv_atom_query_idx}")
            
            
        # otherwise, the atom pairs are within tolerance.
        conf_target.SetAtomPosition(hyv_atom_target_idx,
                                    hyv_atom_query_pos)
        displacement = np.array(hyv_atom_query_pos) - np.array(hyv_atom_target_pos)

        # now handle the hydrogen atoms.
        hyv_atom_target = target.GetAtomWithIdx(hyv_atom_target_idx)
        hyv_atom_query = query.GetAtomWithIdx(hyv_atom_query_idx)
        
        targetNumHs = hyv_atom_target.GetTotalNumHs(includeNeighbors=True)
        queryNumHs = hyv_atom_query.GetTotalNumHs(includeNeighbors=True)
        
        if debug: print(f"Target@{hyv_atom_target_idx} has {targetNumHs} Hs")
        if debug: print(f"Query@{hyv_atom_query_idx} has {queryNumHs} Hs")
        
        if queryNumHs == 0:  
            # Ground truth. no H's in query, we skip this atom.
            # H in targets can be randomly placed.
            continue
        
        neighbor_H_query = []
        for neighbor_query in hyv_atom_query.GetNeighbors():
            if neighbor_query.GetAtomicNum() == 1:
                neighbor_H_query.append(neighbor_query.GetIdx())
        assert len(neighbor_H_query) == queryNumHs
        
        neighbor_H_target = []
        for neighbor_target in hyv_atom_target.GetNeighbors():
            if neighbor_target.GetAtomicNum() == 1:
                neighbor_H_target.append(neighbor_target.GetIdx())
        assert len(neighbor_H_target) == targetNumHs
        
        if targetNumHs == queryNumHs:
            # the most likely cases, let's directly assign the pos from query to target.
            for ii in range(targetNumHs):
                conf_target.SetAtomPosition(neighbor_H_target[ii],
                                            conf_query.GetAtomPosition(neighbor_H_query[ii]))
                # an artifacial mapping, cause we force the order of H's of target as query.
                hydrogen_mapping.append((neighbor_H_target[ii], neighbor_H_query[ii]))
        
        else:
            if targetNumHs < queryNumHs:
                # this indicates that one Hydrogen atom in query mol is changed to other group.
                assert queryNumHs - targetNumHs <= 1, \
                        "multiple hydrogen atoms are alchemically changed, but this is not encouraged."
                # we assign the pos from the closest atoms.
                for ii in range(targetNumHs):
                    distance = []
                    for jj in range(len(neighbor_H_query)):
                        dist_i_j = calc_dist(
                                            conf_target.GetAtomPosition(neighbor_H_target[ii]),
                                            conf_query.GetAtomPosition(neighbor_H_query[jj]),
                                            offset=displacement
                                            )
                        distance.apppend((jj, dist_i_j))
                    distance = sorted(distance, key=lambda x: x[1])
                    conf_target.SetAtomPosition(neighbor_H_target[ii],
                                                conf_query.GetAtomPosition(distance[0][0]))
                    neighbor_H_query = [x[0] for x in distance[1:]]
                    hydrogen_mapping.append((neighbor_H_target[ii], distance[0][0]))
                
            else:
                # target > query
                # this indicates that one group in query is changed to Hydrogen. we need to find the closest 
                # H's in target for the H's in querys. and assign those remaining H's in query.
                assert targetNumHs - queryNumHs <= 1, \
                        "multiple hydrogen atoms are alchemically changed, but this is not encouraged."
                # we assign the pos from the closest atoms.
                for ii in range(queryNumHs):
                    distance = []
                    for jj in range(len(neighbor_H_target)):
                        dist_i_j = calc_dist(
                                            conf_target.GetAtomPosition(neighbor_H_target[jj]),
                                            conf_query.GetAtomPosition(neighbor_H_query[ii]),
                                            offset=displacement
                                            )
                        distance.apppend((jj, dist_i_j))
                    distance = sorted(distance, key=lambda x: x[1])
                    conf_target.SetAtomPosition(distance[0][0],
                                                conf_query.GetAtomPosition(neighbor_H_query[ii]))
                    neighbor_H_query = [x[0] for x in distance[1:]]
                    hydrogen_mapping.append((distance[0][0], neighbor_H_query[ii]))

    return target, hydrogen_mapping


def check_formal_charge(target, query, atom_mapping):

    # Check formal charges for molecules.
    formal_charge_query = Chem.GetFormalCharge(query)
    formal_charge_target_embed = Chem.GetFormalCharge(target)
    assert formal_charge_query == formal_charge_target_embed, \
            f"Formal charges of query:{formal_charge_query} and target:{formal_charge_target_embed} are different."
    
    # check the formal chages. Somethims difference in H's are due to the hydrogen adding algorithm.
    # Ideally, FEP expects the same charges in target atoms and query atoms.
    for hvy_atom_pair_idx in atom_mapping:
        hyv_atom_target_idx, hyv_atom_query_idx = hvy_atom_pair_idx
        atarget = target.GetAtomWithIdx(hyv_atom_target_idx)
        aquery = query.GetAtomWithIdx(hyv_atom_query_idx)
        chg_target = atarget.GetFormalCharge()
        chg_query = aquery.GetFormalCharge()
        if chg_target != chg_query:
            Warning(f"Formal charges of target@{hyv_atom_target_idx}:{chg_target} and "\
                            f"query@{hyv_atom_query_idx}:{chg_query} are different.")
            # raise RuntimeError("Formal charges of target and query are different.")
            atarget.SetFormalCharge(chg_query)
            atarget.SetNumExplicitHs(aquery.GetTotalNumHs())
    
    # Check formal charges for each unique atoms
    target_common_core = [x[0] for x in atom_mapping]
    target_unique_atom = [x for x in range(target.GetNumAtoms()) if x not in target_common_core]

    for atom_idx in target_unique_atom:
        charge_i = target.GetAtomWithIdx(atom_idx).GetFormalCharge()
        assert charge_i == 0, (f"Target Unique Atom@{atom_idx} has formal charge {charge_i}, you are changing the charges!!")


def ConstrainedConfGeneration(target, query, randomseed=-1):
    mcs_atom_mapping = find_common_core(target, query)
    
    check_formal_charge(target, query, mcs_atom_mapping)
    target = Chem.AddHs(target)
    target_embed = ConstrainedEmbed(target, query, mcs_atom_mapping, randomseed=randomseed)
    # target_embed_without_H = Chem.RemoveHs(target_embed)
    target_embed_H, Hmap = alignHs(target_embed, query, mcs_atom_mapping, dist_tol = 1.0, debug = False)
    mcs_atom_mapping_with_H = mcs_atom_mapping + Hmap
    target_embed_H = PositionRestrainedEnergyMinimization(target_embed_H, 
                                                          [x[0] for x in mcs_atom_mapping_with_H], 
                                                          force_constant=1000.)
    return target_embed_H


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--ftarget", type=str, default="CC1=CC=C(C=C1)C2=CC=CC=C2")
    parser.add_argument("--fquery", type=str, default="embed.query.sdf")
    parser.add_argument("--output", type=str, default="embed.target.sdf")
    parser.add_argument("--randomseed", type=int, default=-1)
    args = parser.parse_args()
    with open(args.ftarget, "r") as f:
        target_SMILES = f.read()
    target = Chem.MolFromSmiles(target_SMILES)
    query = Chem.SDMolSupplier(args.fquery, removeHs=False)[0]
    target_embed_H = ConstrainedConfGeneration(target, query, args.randomseed)

    writer = Chem.SDWriter(args.output)
    writer.write(target_embed_H, confId=0)