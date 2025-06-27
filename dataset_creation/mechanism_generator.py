import random
import copy
import json
import jsonlines
from rdchiral import template_extractor
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdmolops


class changed_atoms:
    def __init__(self, atom):
        self.atom = atom
        self.neighbor = []
        if atom.HasProp('molAtomMapNumber'):
            idx = int(atom.GetProp('molAtomMapNumber'))
            self.atom_idx = idx
        else:
            self.atom_idx = -1

    def add_neighbor(self, neighbor):
        info = [neighbor]
        if neighbor.HasProp('molAtomMapNumber'):
            idx = int(neighbor.GetProp('molAtomMapNumber'))
            info.append(idx)
        else:
            info.append(-1)
        self.neighbor.append(info)

    def get_atom_idx(self):
        return self.atom_idx

    def get_neighbor(self):
        return self.neighbor

    def get_self(self):
        return self.atom


class different_molecule_overlap:
    def __init__(self, molecule):
        self.molecule = molecule
        self.list_of_idx = []

    def combine_idx_list(self, list_of_idx):
        self.list_of_idx = self.list_of_idx + list_of_idx

    def update_idx_list(self, list_of_idx):
        self.list_of_idx.append(list_of_idx)



def atom_mapper(molecule, changed_atom, idx_list, react_radius_list, reactant_mol_list):
    bad = 0
    checkpoint_num = 0
    try:
        all_mol_smi = Chem.MolToSmiles(molecule)
        list_all_pdt = all_mol_smi.split('.')
        new_list = []
        
        if any(molecule.HasSubstructMatch(mol) and mol.HasSubstructMatch(molecule) for mol in reactant_mol_list):
            return None
        # Check whether the product is unreacted
        
        # Refers to the number of atoms that is in idx_list and in product
        reacting_atom = 0
        for smi in list_all_pdt:
            mol = Chem.MolFromSmiles(smi, sanitize=False)
            for atom in mol.GetAtoms():
                if atom.HasProp("molAtomMapNumber"):
                    if atom.GetProp('molAtomMapNumber') in idx_list:
                        reacting_atom += 1
            new_list.append(mol)

        if reacting_atom == len(idx_list) or (not new_list):
            return None
        # Determine whether we have a molecule with atom mapping

        max_num_atom_mapping = 0
        for mol in new_list:
            all_idx = []
            for atom in mol.GetAtoms():
                if atom.HasProp("molAtomMapNumber"):
                    all_idx.append(atom.GetProp("molAtomMapNumber"))
            if len(all_idx) > max_num_atom_mapping:
                max_num_atom_mapping = len(all_idx)
                molecule = mol

        remaining_atom_mapping = []
        for atom in molecule.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                remaining_atom_mapping.append(str(atom.GetProp('molAtomMapNumber')))
        # Get remaining atom index and store in remaining_atom_mapping

        checkpoint_num = 1

        if all(idx in remaining_atom_mapping for idx in react_radius_list):
            return None

        idx_to_be_mapped = list(set(react_radius_list) - set(remaining_atom_mapping))
        properly_mapped = False

        while idx_to_be_mapped and not properly_mapped:

            map_list = []
            checkpoint_num = 2

            for atom in molecule.GetAtoms():
                if atom.HasProp("molAtomMapNumber"):
                    if atom.GetProp('molAtomMapNumber') in idx_to_be_mapped:
                        idx_to_be_mapped.remove(atom.GetProp('molAtomMapNumber'))

            checkpoint_num = 3

            changed_atom_copy = []
            for map_candidate in changed_atom:
                if str(map_candidate.get_atom_idx()) in idx_to_be_mapped:
                    changed_atom_copy.append(map_candidate)
            changed_atom = changed_atom_copy

            checkpoint_num = 4

            for atom in molecule.GetAtoms():
                if atom.HasProp("molAtomMapNumber"):
                    continue
                else:
                    neighbors = atom.GetNeighbors()
                    if neighbors:
                        most_probable = []
                        prev_count = 0
                        count = 0
                        for change in neighbors:
                            for map_candidate in changed_atom:
                                if map_candidate.get_self().GetAtomicNum() != atom.GetAtomicNum():
                                    continue
                                else:
                                    count += 1
                                for comparator in map_candidate.get_neighbor():
                                    if bool(change.HasProp("molAtomMapNumber")):
                                        if int(comparator[1]) == int(change.GetProp('molAtomMapNumber')):
                                            count += 10
                                    if not atoms_are_different(comparator[0], change):
                                        count += 5
                                    if comparator[0].GetAtomicNum() == change.GetAtomicNum():
                                        count += 1
                                if count > prev_count:
                                    most_probable = []
                                    most_probable.append(atom)
                                    most_probable.append(map_candidate.get_atom_idx())
                                    most_probable.append(count)
                                    prev_count = count
                                    count = 0
                                elif count == prev_count:
                                    if count != 0:
                                        most_probable = []
                                        most_probable.append(atom)
                                        most_probable.append(map_candidate.get_atom_idx())
                                        most_probable.append(count)
                                        count = 0
                            if most_probable:
                                if most_probable not in map_list:
                                    map_list.append(most_probable)
                    else:
                        most_probable = []
                        prev_count = 0
                        count = 0
                        for map_candidate in changed_atom:
                            if atom.GetAtomicNum() == map_candidate.get_self().GetAtomicNum():
                                count += 1
                            if count > prev_count:
                                    most_probable = []
                                    most_probable.append(atom)
                                    most_probable.append(map_candidate.get_atom_idx())
                                    most_probable.append(count)
                                    prev_count = count
                                    count = 0
                            elif count == prev_count:
                                if count != 0:
                                    most_probable = []
                                    most_probable.append(atom)
                                    most_probable.append(map_candidate.get_atom_idx())
                                    most_probable.append(count)
                                    count = 0
                        if most_probable:
                            if most_probable not in map_list:
                                map_list.append(most_probable)

            checkpoint_num = 5

            if idx_to_be_mapped and (not map_list):

                if all(a.HasProp('molAtomMapNumber') for a in molecule.GetAtoms()):
                    checkpoint_num = 99
                    properly_mapped = True
                    break
                elif not changed_atom:
                    checkpoint_num = 9999
                    break
                else:
                    checkpoint_num = 999999
                    bad += 1
                    if bad < 2:
                        for atom in molecule.GetAtoms():
                            if not atom.HasProp('molAtomMapNumber'):
                                list_of_neighbor_idx = []
                                for neighbor in atom.GetNeighbors():
                                    if neighbor.HasProp('molAtomMapNumber'):
                                        list_of_neighbor_idx.append(int(neighbor.GetProp('molAtomMapNumber')))
                                for atom2 in changed_atom:
                                    if atom.GetAtomicNum() == atom2.get_self().GetAtomicNum():
                                        if any((i - 1) <= atom2.get_atom_idx() <= (i + 1) for i in list_of_neighbor_idx):
                                            atom.SetAtomMapNum(atom2.get_atom_idx())
                                            break
                    else:
                        break

            ## In case of Diels_Alder reaction

            checkpoint_num = 999

            map_list_sorted = sorted(map_list, key=lambda x: x[2], reverse=True)
            mapped = []
            mapped_idx = []

            for i in map_list_sorted:
                if mapped:
                    if (i[0] not in mapped and i[1] not in mapped_idx):
                        i[0].SetAtomMapNum(i[1])
                        mapped.append(i[0])
                        mapped_idx.append(i[1])
                elif not mapped:
                    i[0].SetAtomMapNum(i[1])
                    mapped.append(i[0])
                    mapped_idx.append(i[1])

            checkpoint_num = 6

            remaining_atom_mapping = []
            for a in molecule.GetAtoms():
                if a.HasProp("molAtomMapNumber"):
                    remaining_atom_mapping.append(a.GetProp('molAtomMapNumber'))

            checkpoint_num = 7

        return molecule
    except Exception as e:
        return None


def atoms_are_different(atom1, atom2):
    if atom1.GetAtomicNum() != atom2.GetAtomicNum(): return True  # must be true for atom mapping
    if atom1.GetFormalCharge() != atom2.GetFormalCharge(): return True
    if atom1.GetDegree() != atom2.GetDegree(): return True
    # if atom1.IsInRing() != atom2.IsInRing(): return True # do not want to check this!
    # e.g., in macrocycle formation, don't want the template to include the entire ring structure
    if atom1.GetNumRadicalElectrons() != atom2.GetNumRadicalElectrons(): return True
    if atom1.GetIsAromatic() != atom2.GetIsAromatic(): return True

    # Check bonds and nearest neighbor identity
    bonds1 = sorted([template_extractor.bond_to_label(bond) for bond in atom1.GetBonds()])
    bonds2 = sorted([template_extractor.bond_to_label(bond) for bond in atom2.GetBonds()])
    if bonds1 != bonds2: return True

    return False


def combine_molecules(overlap, reacting_mol):
    finalized_mol_list = []
    for mol1 in reacting_mol:
        sum_mol = different_molecule_overlap(mol1)
        for mol2 in overlap:
            if mol2.molecule == sum_mol.molecule:
                sum_mol.combine_idx_list(mol2.list_of_idx)
        finalized_mol_list.append(sum_mol)
    return finalized_mol_list


def database_construct(database_path, db_num, start_step, debug):
    if start_step == 1:

        sampler = 0
        with open(database_path, errors="ignore") as database:
            rxn_idx = 1
            reactions = database.readlines()
            
            if debug:
                sampler = random.randint(1, 900)
                temp_reactions = []
                for reaction in reactions:
                    if sampler == 3:
                        mol_list = []
                        if '>>' in reaction:
                            mol_list = reaction.split('>>')
                            mol_list = mol_list[0].split('.') + mol_list[1].split('.')
                        else:
                            mol_list = reaction.split('>')
                            mol_list = mol_list[0].split('.') + mol_list[1].split('.') + mol_list[2].split('.')
                        try:
                            for mol in mol_list:
                                mol = Chem.MolFromSmiles(mol)
                                if mol is None:
                                    raise ValueError
                            temp_reactions.append(reaction)
                        except:
                            pass
                    sampler = random.randint(1, 900)
                reactions = temp_reactions

            ct = len(reactions)
            
            if not debug:
                chunk_size = ct // db_num
                reactions_list = [reactions[i : i + chunk_size] for i in range(0, len(reactions), chunk_size)]
                if len(reactions_list != db_num):
                    reactions_list[-1] = reactions_list[-1] + reactions_list[-2]
                    reactions_list.pop(-2)

            for i in range(db_num):
                if debug:
                    num_rxn = 0
                    if i == db_num - 1:
                        num_rxn = ct // db_num + ct % db_num
                    else:
                        num_rxn = ct // db_num
                else:
                    reactions = reactions_list[i]
                    num_rxn = len(reactions)

                with open(f"database/database_{i}.json", 'w') as write:
                    for j in range(num_rxn):
                        reaction = reactions.pop(0).split(' ')[0].split('\t')[0].strip()
                        json.dump({"Reaction index": rxn_idx, "Reaction": reaction}, write)
                        write.write('\n')  # Add newline after each reaction dictionary
                        rxn_idx += 1
            print(rxn_idx-1)

    else:
        ct = 0
        with jsonlines.open(f"Data/Potential_Match.json", 'r') as database:
            for dict in database:
                ct += 1

        with jsonlines.open(f"Data/Potential_Match.json", 'r') as database:
            num_rxn = ct // db_num
            for i in range(db_num):
                if i == db_num - 1:
                    num_rxn = ct // db_num + ct % db_num
                with open(f"database/database_{i}.json", 'w') as write:
                    for j in range(num_rxn):
                        reaction = database.read()
                        json.dump(reaction, write)
                        write.write('\n')  # Add newline after each reaction dictionary


def duplicate_template(num_parallel):
    with open(r'database_json') as db:
        temp = db.read()
        for num in range(num_parallel):
            with open(r'mech_db/mechanism_database_%d.txt' % num, 'a') as w:
                w.write(temp)


def extract_reaction(all_rxn):
    reaction = all_rxn.readline()
    if 'index' in reaction:
        return reaction
    reaction = reaction.split(' ')[0]
    reaction = reaction.split('\t')[0]
    if reaction == '':
        return None
    return reaction


def extract_template(all_template):
    mechanism = all_template.readline()
    is_mechanism = mechanism.find('>')
    is_name = mechanism.find('_')
    if is_name != -1:
        return mechanism
    while is_mechanism == -1:
        if not bool(mechanism):
            return None
        mechanism = all_template.readline()
        is_mechanism = mechanism.find('>')
    return mechanism


def get_reacting_atom(reaction, radius):
    reactant = get_reactant_smiles(reaction)
    pdt_holder = get_product(reaction)
    reactant_match = Chem.MolFromSmiles(reactant)
    rct_ls = []
    pdt_ls = []
    rct_ls.append(reactant_match)
    pdt_mol = Chem.MolFromSmiles(pdt_holder[0])
    if pdt_mol is None:
        return None
    pdt_ls.append(pdt_mol)
    idx_list = list(map(int, template_extractor.get_changed_atoms(rct_ls, pdt_ls)[-2]))
    all_atom_idx = []
    for atom in reactant_match.GetAtoms():
        if atom.HasProp('molAtomMapNumber'):
            all_atom_idx.append(int(atom.GetProp('molAtomMapNumber')))

    react_radius_temp = []

    for i in idx_list:
        if i == 1:
            if i + radius > max(all_atom_idx):
                for j in range(1, max(all_atom_idx) + 1):
                    if j not in set(react_radius_temp):
                        react_radius_temp.append(j)
            else:
                for j in range(1, radius + 1):
                    if j not in set(react_radius_temp):
                        react_radius_temp.append(j)
        elif i == max(all_atom_idx):
            if max(all_atom_idx) - radius < 1:
                for j in range(1, i + 1):
                    if j not in set(react_radius_temp):
                        react_radius_temp.append(j)
            else:
                for j in range(max(all_atom_idx) - (radius), i + 1):
                    if j not in set(react_radius_temp):
                        react_radius_temp.append(j)
        else:
            if i + radius > max(all_atom_idx):
                if i - radius < 1:
                    for j in range(1, max(all_atom_idx) + 1):
                        if j not in set(react_radius_temp):
                            react_radius_temp.append(j)
                else:
                    for j in range(i - radius, max(all_atom_idx) + 1):
                        if j not in set(react_radius_temp):
                            react_radius_temp.append(j)

            elif i - radius < 1:
                if i + radius > max(all_atom_idx):
                    for j in range(1, max(all_atom_idx) + 1):
                        if j not in set(react_radius_temp):
                            react_radius_temp.append(j)
                else:
                    for j in range(1, i + radius):
                        if j not in set(react_radius_temp):
                            react_radius_temp.append(j)

            else:
                for j in range(i - radius, i + radius + 1):
                    if j not in set(react_radius_temp):
                        react_radius_temp.append(j)
    react_radius_temp = list(map(str, react_radius_temp))
    return react_radius_temp


def get_product(reaction):
    reaction = reaction.split('\n')[0]
    is_rxn = reaction.find('>')
    if is_rxn == -1:
        return None
    if 'smiles' in reaction:
        reaction = reaction.split(' ')[2]
    else:
        reaction = reaction.split()[0]
    locate = reaction.find('>>')
    if locate != -1:
        reaction = reaction.split('>>')[1]
        reaction = reaction.split('.')
    else:
        reaction = reaction.split('>')[2]
        reaction = reaction.split('.')
    return reaction


def get_reactant_mol_list(reaction):
    if reaction is None:
        return None

    reactant_mol_list = []
    if reaction.find('>>') == -1:
        hold = reaction.split('>')
        all_reactant_smi = ''.join((hold[0], '.', hold[1]))
        reactant_smi_list = all_reactant_smi.split('.')
        for smi in reactant_smi_list:
            reactant_mol_list.append(Chem.MolFromSmiles(smi))
        return reactant_mol_list

    if reaction.find('>>') != -1:
        all_reactant_smi = reaction.split('>')[0]
        reactant_smi_list = all_reactant_smi.split('.')
        for smi in reactant_smi_list:
            reactant_mol_list.append(Chem.MolFromSmiles(smi))
        return reactant_mol_list


def get_reactant_smiles(reaction):
    if reaction is None:
        return None
    if reaction.find('>>') == -1:
        hold = reaction.split('>')
        reactant = ''.join((hold[0], '.', hold[1]))
        return reactant
    if reaction.find('>>') != -1:
        reactant = reaction.split('>')[0]
        return reactant


def inter_molecular_rxn(reactant, template):
    try:
        template = rdChemReactions.ReactionFromSmarts(template)
        product = template.RunReactants(reactant)
        return product
    except ValueError as ex:
        if str(ex) == ("ChemicalParserException: Number of reactants provided does not match number of reactant "
                       "templates."):
            return None


def intra_molecular_rxn(reactant, template):
    template = template.split('>>')
    template[0] = ''.join(('(', template[0], ')'))
    template[1] = ''.join(('(', template[1], ')'))
    template = ''.join((template[0], '>>', template[1]))
    template = rdChemReactions.ReactionFromSmarts(template)
    product = template.RunReactants(reactant)
    if product is not None:
        return product
    else:
        return None


def multi_step_reactant_update(reactant, pdt_hold):
    rct_hold = reactant.split('.')
    updated_rct_hold = []
    rct_mol_mapping_list = []
    pdt_mol_mapping_dict = {'smi': pdt_hold[0], 'Atom mapping': []}
    for smi in rct_hold:
        rct_mol_mapping_dict = {'smi': smi, 'Atom mapping': []}
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            return None
        for atom in mol.GetAtoms():
            if atom.HasProp('molAtomMapNumber'):
                rct_mol_mapping_dict['Atom mapping'].append(atom.GetProp('molAtomMapNumber'))
        rct_mol_mapping_list.append(rct_mol_mapping_dict)
    pdt_mol = Chem.MolFromSmarts(pdt_mol_mapping_dict['smi'])
    for atom in pdt_mol.GetAtoms():
        if atom.HasProp('molAtomMapNumber'):
            pdt_mol_mapping_dict['Atom mapping'].append(atom.GetProp('molAtomMapNumber'))
    for rct_dict in rct_mol_mapping_list:
        if not set(pdt_mol_mapping_dict['Atom mapping']).intersection(set(rct_dict['Atom mapping'])):
            updated_rct_hold.append(rct_dict['smi'])
    updated_rct_hold.append(pdt_mol_mapping_dict['smi'])
    reactant = '.'.join(updated_rct_hold)
    return reactant


def overlap_solver(list_of_tuples):
    count_dict = {}  # create an empty dictionary to hold the count of each integer
    for tuple1 in list_of_tuples:
        for elem in tuple1:
            for num in elem:
                if num in count_dict:
                    count_dict[num] += 1
                else:
                    count_dict[num] = 1
    return [num for num in count_dict if count_dict[num] > 1]



def try_to_apply_mechanistic_template(reaction, mechanism, react_radius_list, multi_step):
    reactant = get_reactant_smiles(reaction)
    pdt_holder = get_product(reaction)
    reactant_match = Chem.MolFromSmiles(reactant)
    if reactant_match is None:
        return None
    if multi_step:
        reactant = multi_step_reactant_update(reactant, pdt_holder)
    if not reactant:
        return None
    temp_reactant = Chem.MolFromSmarts(mechanism.split('>')[0])

    # Important: Get the changed atoms

    changed_atom = []

    list_of_reactant = reactant.split('.')
    temp_reactant = mechanism.split('>')[0]
    temp_reactant_list = temp_reactant.split('.')
    reacting_mol = []
    overlap = []
    for element_mechanism in temp_reactant_list:
        for element_reactant in list_of_reactant:
            element_reactant = Chem.MolFromSmiles(element_reactant)
            if element_reactant is None:
                return None
            temp_mol = Chem.MolFromSmarts(element_mechanism)
            test = element_reactant.HasSubstructMatch(temp_mol)
            if test:
                molecule = different_molecule_overlap(element_reactant)
                molecule.update_idx_list(element_reactant.GetSubstructMatches(temp_mol))
                for atom in element_reactant.GetAtoms():
                    if not atom.HasProp('molAtomMapNumber'):
                        continue
                    if (atom.GetProp('molAtomMapNumber')) in react_radius_list:
                        atom_info = changed_atoms(atom)
                        neighbors = atom.GetNeighbors()
                        if neighbors:
                            for change in neighbors:
                                atom_info.add_neighbor(change)
                        changed_atom.append(atom_info)
                    elif (atom.GetProp('molAtomMapNumber')) not in react_radius_list:
                        atom.SetProp('_protected', '1')
                same_mol = False
                if reacting_mol:
                    if any(i.HasSubstructMatch(element_reactant) and element_reactant.HasSubstructMatch(i) for i in reacting_mol):
                        same_mol = True
                if not same_mol:
                    reacting_mol.append(element_reactant)
                    overlap.append(molecule)
    combined_mol_list = combine_molecules(overlap, reacting_mol)
    for mol in combined_mol_list:
        mol.list_of_idx = overlap_solver(mol.list_of_idx)
    for mol in reacting_mol:
        rdmolops.FastFindRings(mol)
        for mol2 in combined_mol_list:
            if mol == mol2.molecule:
                overlap_atom = mol2.list_of_idx
                for idx in overlap_atom:
                    atom = mol.GetAtomWithIdx(idx)
                    atom.SetProp('_protected', '1')
                    for a in atom.GetNeighbors():
                        a.SetProp('_protected', '1')
                break
    product_temp = []
    if len(reacting_mol) == 1:
        product_temp_holder = []
        reactant_mol = Chem.RWMol()
        for mol in reacting_mol:
            reactant_mol = Chem.CombineMols(reactant_mol, mol)
        rdmolops.FastFindRings(reactant_mol)
        rxn_input = [reactant_mol]
        product_temp = intra_molecular_rxn(rxn_input, mechanism)
        for product_tuple in product_temp:
            for pdt_mol in product_tuple:
                product_temp_holder.append(pdt_mol)
        product_temp = product_temp_holder
    elif len(reacting_mol) > 1:
        product_temp_holder = []
        for mol in reacting_mol:
            product_temp = intra_molecular_rxn([mol], mechanism)
            for product_tuple in product_temp:
                for pdt_mol in product_tuple:
                    product_temp_holder.append(pdt_mol)
        product_temp = inter_molecular_rxn(reacting_mol, mechanism)
        if product_temp:
            for product_tuple in product_temp:
                for pdt_mol in product_tuple:
                    product_temp_holder.append(pdt_mol)
        product_temp = product_temp_holder
    else:
        return None
    product = []
    for fil in product_temp:
        mapped_product = atom_mapper(fil, changed_atom, react_radius_list, react_radius_list, get_reactant_mol_list(reaction))
        if mapped_product:
            product.append(mapped_product)
    product_list = []
    for i in product:
        match = False
        if not product_list:
            product_list.append(i)
        for j in product_list:
            if bool(i.HasSubstructMatch(j)) and bool(j.HasSubstructMatch(i)):
                match = True
        if not match:
            if Chem.MolToSmiles(i):
                product_list.append(i)
    if len(product_list) > 50:
        return None
    return product_list
