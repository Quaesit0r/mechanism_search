import json
import mechanism_generator
import time
from rdkit import Chem


class reaction:
    def __init__(self, reaction_dict, radius):
        self.id = reaction_dict.get("Reaction index")
        self.reaction_smi = reaction_dict.get("Reaction")
        self.ground_truth_product = mechanism_generator.get_product(self.reaction_smi)[0]
        self.radius = radius
        self.react_idx_list = mechanism_generator.get_reacting_atom(self.reaction_smi, radius)
        self.dead_node = 0
        self.alive = reaction_dict.get("alive", True)
        self.subsequent_mechanisms = reaction_dict.get("subsequent mechanisms", [])

    def add_mechanism(self, mechanism_obj):
        self.subsequent_mechanisms.append(mechanism_obj)

    def check_alive(self):
        if not self.subsequent_mechanisms or self.dead_node == len(self.subsequent_mechanisms):
            self.alive = False
        return self.alive

    def reaction_filter(self):
        if self.reaction_smi.count('C') + self.reaction_smi.count('c') > 200 or not self.reaction_smi or len(mechanism_generator.get_product(self.reaction_smi)) != 1:
            return False
        else:
            return True


class mechanism_node:
    def __init__(self, mechanism_dict, parent_node):
        self.parent = parent_node
        self.step = mechanism_dict.get("Step")
        self.id = mechanism_dict.get("Mechanism index", 0)
        self.mechanism_name = mechanism_dict.get("Mechanism name", None)
        self.template = mechanism_dict.get("Template", None)
        self.mechanism_smi = mechanism_dict.get("Mechanism smi", None)
        self.radius = mechanism_dict.get("Radius", None)
        self.alive = mechanism_dict.get("alive", True)
        self.ground_truth_product = self.parent.ground_truth_product
        self.react_idx_list = self.parent.react_idx_list
        self.dead_node = 0
        self.subsequent_mechanisms = mechanism_dict.get("subsequent mechanisms", [])

    def add_mechanism(self, mechanism_obj):
        self.subsequent_mechanisms.append(mechanism_obj)

    def check_alive(self):
        if not self.subsequent_mechanisms or self.dead_node == len(self.subsequent_mechanisms):
            self.alive = False
        return self.alive

    def get_mechanism_name(self):
        return self.mechanism_name


class template_list:
    def __init__(self, mechanism_database):
        self.mechanism_holder = None
        with open(mechanism_database) as mech_db:
            self.template = json.load(mech_db)

    def get_mechanism(self):
        return self.template



def apply_all_mechanism(obj, template_obj, search_step, start_time=0):
    # Helper function to apply a mechanism and determine a match
    
    def apply_mechanism(obj, mechanism_dict, search_step, start_time):
        # Helper function to determine a match in the product molecule list
        if time.time() - start_time > 30:
            return [], False
        
        def determining_match(product_mol_list):
            finalized_obj_list = []
            for product_mol in product_mol_list:
                # Check if the product molecule matches the ground truth product and vice versa
                if product_mol.HasSubstructMatch(Chem.MolFromSmiles(obj.ground_truth_product)) and Chem.MolFromSmiles(obj.ground_truth_product).HasSubstructMatch(product_mol):
                    # Update the mechanism SMILES and create a mechanism object
                    if isinstance(obj, reaction):
                        mechanism_smi = update_mechanism_smi(obj.reaction_smi, product_mol)
                    else:
                        mechanism_smi = update_mechanism_smi(obj.mechanism_smi, product_mol, True)
                    data_input_dict = {"Mechanism name": mechanism_name, "Template": template_smi,
                                       "Mechanism smi": mechanism_smi, "Step": search_step,
                                       "Radius": obj.radius}
                    mechanism_obj = mechanism_node(data_input_dict, obj)
                    return mechanism_obj, True
                else:
                    # Update the mechanism SMILES and create a mechanism object
                    if isinstance(obj, reaction):
                        mechanism_smi = update_mechanism_smi(obj.reaction_smi, product_mol)
                    else:
                        mechanism_smi = update_mechanism_smi(obj.mechanism_smi, product_mol, True)
                    data_input_dict = {"Mechanism name": mechanism_name, "Template": template_smi,
                                       "Mechanism smi": mechanism_smi, "Step": search_step,
                                       "Radius": obj.radius}
                    mechanism_obj = mechanism_node(data_input_dict, obj)
                    finalized_obj_list.append(mechanism_obj)
            return finalized_obj_list, False

        if isinstance(obj, reaction):
            obj_smi = obj.reaction_smi
        else:
            obj_smi = obj.mechanism_smi
        mechanism_name = mechanism_dict.get("name")
        template_smi_list = mechanism_dict.get("template")
        all_product_obj_list = []

        # Apply the mechanism templates to the object's reaction SMILES
        for template_smi in template_smi_list:
            if isinstance(obj, reaction):
                # Apply the mechanism template to the reaction SMILES without considering the root node
                product_list = mechanism_generator.try_to_apply_mechanistic_template(obj_smi, template_smi, obj.react_idx_list, False)
                if product_list:
                    if len(product_list) > 100:
                        return [], False
                    product_nodes, good_match = determining_match(product_list)
                    if good_match:
                        return product_nodes, good_match
                    else:
                        all_product_obj_list = all_product_obj_list + product_nodes
                        if len(all_product_obj_list) > 100:
                            return [], False
            else:
                # Apply the mechanism template to the reaction SMILES considering the root node
                product_list = mechanism_generator.try_to_apply_mechanistic_template(obj_smi, template_smi, obj.react_idx_list, True)
                if product_list:
                    if len(product_list) > 100:
                        return [], False
                    product_nodes, good_match = determining_match(product_list)
                    if good_match:
                        return product_nodes, good_match
                    else:
                        all_product_obj_list = all_product_obj_list + product_nodes
                        if len(all_product_obj_list) > 100:
                            return [], False

        # Update the parent of each product node
        for single_product_node in all_product_obj_list:
            single_product_node.parent = obj

        return all_product_obj_list, False

    start_time = time.time()
    finalized_list = []
    bottom_layer = False

    # Check if the object is alive
    if not obj.alive:
        return None

    if isinstance(obj, reaction):
        if len(obj.subsequent_mechanisms) == 0 and obj.alive:
            bottom_layer = True
            pass
        else:
            for child_obj in obj.subsequent_mechanisms:
                # Recursively apply mechanisms to the child object
                product_node = apply_all_mechanism(child_obj, template_obj, search_step, start_time)
                if product_node:
                    match_node = child_obj
                    obj.subsequent_mechanisms = []
                    obj.subsequent_mechanisms.append(match_node)
                    return product_node
                else:
                    # Check if the object is not alive
                    obj.check_alive()
                    if not obj.alive:
                        return None

    # Recursive call to get to the desired search step
    elif isinstance(obj, mechanism_node):
        if search_step-1 != obj.step:
            for child_obj in obj.subsequent_mechanisms:
                # Recursively apply mechanisms to the child object
                product_node = apply_all_mechanism(child_obj, template_obj, search_step, start_time)
                if product_node:
                    match_node = child_obj
                    obj.subsequent_mechanisms = []
                    obj.subsequent_mechanisms.append(match_node)
                    return product_node
                else:
                    # Check if the object is not alive and it is not the root node
                    alive = obj.check_alive()
                    if not alive and not isinstance(obj, reaction):
                        obj.parent.dead_node += 1
                        return None
        else:
            bottom_layer = True

    if bottom_layer:
        # Apply all mechanisms in the list
        for mechanism_dict in template_obj:
            if time.time() - start_time > 150:
                finalized_list = []
                break
            product_node, match = apply_mechanism(obj, mechanism_dict, search_step, start_time)
            if match:
                # Add the product node to the subsequent mechanisms of the object
                obj.subsequent_mechanisms.append(product_node)
                return product_node
            else:
                finalized_list = finalized_list + product_node

        #TODO: Add functions to remove the repeating items
        finalized_node_list = []
        for node in finalized_list:
            if time.time() - start_time > 150:
                finalized_list = []
                break
            match = False
            if not finalized_node_list:
                finalized_node_list.append(node)
            else:
                for finalized_node in finalized_node_list:
                    try:
                        product = Chem.MolFromSmiles(node.mechanism_smi.split('>>')[1])
                        finalized_product = Chem.MolFromSmiles(finalized_node.mechanism_smi.split('>>')[1])
                        if product.HasSubstructMatch(finalized_product) and finalized_product.HasSubstructMatch(product):
                            match = True
                    except Exception as e:
                        continue
            if not match:
                finalized_node_list.append(node)
        finalized_list = finalized_node_list
        if len(finalized_list) > 200:
            finalized_list = []

        # Check if finalized list is empty
        if not finalized_list:
            alive = obj.check_alive()
            if not isinstance(obj, reaction) and not alive:
                # Increment the dead_node counter of the parent object
                obj.parent.dead_node += 1
        else:
            # Assign ids to the finalized mechanisms and update the subsequent mechanisms of the object
            for index, mechanism_obj in enumerate(finalized_list):
                mechanism_obj.id = index
            obj.subsequent_mechanisms = finalized_list

    return None


def good_match_processor(mechanism_obj):
    if not isinstance(mechanism_obj, reaction):
        parent = mechanism_obj.parent
        index = mechanism_obj.id
        if not isinstance(parent, reaction):
            good_match_processor(parent)
        parent.subsequent_mechanisms = parent.subsequent_mechanisms[index]
        return None
    else:
        return None


def node_decoder(dict, radius, parent_node=None):
    if "subsequent mechanisms" not in dict:
        # If "subsequent mechanism" not in dict aka not initialized, create a reaction object
        reaction_obj = reaction(dict, radius)
        return reaction_obj
    elif dict.get("Reaction"):
        # If "is root" key is present in the dictionary, create a reaction object and process subsequent mechanisms
            reaction_obj = reaction(dict, radius)
            subsequent_mechanisms_obj_list = []
            if dict.get("subsequent mechanisms"):
                for subsequent_mechanisms_dict in dict.get("subsequent mechanisms"):
                    # Recursively call node_decoder to process each subsequent mechanism
                    subsequent_mechanisms_node = node_decoder(subsequent_mechanisms_dict, radius, reaction_obj)
                    subsequent_mechanisms_obj_list.append(subsequent_mechanisms_node)
                # Set the list of subsequent mechanisms as an attribute of the reaction object
            reaction_obj.subsequent_mechanisms = subsequent_mechanisms_obj_list
            return reaction_obj
    else:
        # If neither "step" nor "is root" key is present, create a subsequent mechanism object
        subsequent_mechanisms_obj = mechanism_node(dict, parent_node)
        subsequent_mechanisms_obj_list = []
        if dict.get("subsequent mechanisms"):
            for subsequent_mechanisms_dict in dict.get("subsequent mechanisms"):
                # Recursively call node_decoder to process each subsequent mechanism
                subsequent_mechanisms_obj_child = node_decoder(subsequent_mechanisms_dict, radius, subsequent_mechanisms_obj)
                subsequent_mechanisms_obj_list.append(subsequent_mechanisms_obj_child)
        # Set the list of subsequent mechanisms as an attribute of the subsequent mechanism object
        subsequent_mechanisms_obj.subsequent_mechanisms = subsequent_mechanisms_obj_list
        return subsequent_mechanisms_obj



def node_encoder(obj):
    if isinstance(obj, reaction):
        # If the object is a reaction (root), encode its attributes
        rxn_dict = {
            "Reaction index": obj.id,
            "Reaction": obj.reaction_smi,
            "Ground truth product": obj.ground_truth_product,
            "alive": obj.alive,
            "subsequent mechanisms": []
        }
        subsequent_mechanisms = obj.subsequent_mechanisms
        # Recursively encode the subsequent mechanisms
        if isinstance(subsequent_mechanisms, list):
            for mechanism in subsequent_mechanisms:
                encoded_subsequent_mechanism = node_encoder(mechanism)
                rxn_dict["subsequent mechanisms"].append(encoded_subsequent_mechanism)
        else:
            encoded_subsequent_mechanism = node_encoder(subsequent_mechanisms)
            rxn_dict["subsequent mechanisms"].append(encoded_subsequent_mechanism)
        return rxn_dict
    else:
        # If the object is a subsequent mechanism, encode its attributes
        mech_dict = {
            "Step": obj.step,
            "ID": obj.id,
            "Mechanism name": obj.mechanism_name,
            "Template": obj.template,
            "Mechanism smi": obj.mechanism_smi,
            "Radius": obj.radius,
            "alive": obj.alive,
            "subsequent mechanisms": []
        }
        subsequent_mechanisms = obj.subsequent_mechanisms
        # Recursively encode the subsequent mechanisms
        if isinstance(subsequent_mechanisms, list):
            for mechanism in subsequent_mechanisms:
                encoded_subsequent_mechanism = node_encoder(mechanism)
                mech_dict["subsequent mechanisms"].append(encoded_subsequent_mechanism)
        else:
            encoded_subsequent_mechanisms = node_encoder(subsequent_mechanisms)
            mech_dict["subsequent mechanisms"].append(encoded_subsequent_mechanisms)
        return mech_dict


def update_mechanism_smi(original_smi, product_mol, multistep=False):
    try:
        product_smi = Chem.MolToSmiles(product_mol)
    except Exception as e:
        return None
    smi_to_be_modified = original_smi
    if '>>' in smi_to_be_modified:
        if multistep:
            smi_to_be_modified = mechanism_generator.multi_step_reactant_update(smi_to_be_modified.split('>>')[0], [smi_to_be_modified.split('>>')[1]])
        else:
            smi_to_be_modified = smi_to_be_modified.split('>>')[0]
        if not smi_to_be_modified:
            print('smi problem')
            print(original_smi)
            print(product_smi)
        smi_to_be_modified = ''.join((smi_to_be_modified, '>>', product_smi))
    else:
        smi_to_be_modified = original_smi.split('>')
        smi_to_be_modified = ''.join(('.'.join((smi_to_be_modified[0], smi_to_be_modified[1])), '>>', product_smi))
        
        # For further handling of reagent individually
        """
        else:
        if multistep:
            reactant_smi_hold, reagent_smi_hold = original_smi.split('>')[0], original_smi.split('>')[1]
            reactant_smi_hold = mechanism_generator.multi_step_reactant_update(reactant_smi_hold, [original_smi.split('>')[2]])
            reagent_smi_hold = mechanism_generator.multi_step_reactant_update(reagent_smi_hold, [original_smi.split('>')[2]])
            reactant_smi_hold = reactant_smi_hold.split('.')
            reagent_smi_hold = reagent_smi_hold.split('.')
            repetitive_obj = []
            if set(reactant_smi_hold).intersection(reagent_smi_hold):
                repetitive_obj = list(set(reactant_smi_hold).intersection(reagent_smi_hold))
            for smi in repetitive_obj:
                reagent_smi_hold.remove(smi)
            reactant_smi_hold = '.'.join(reactant_smi_hold)
            reagent_smi_hold = '.'.join(reagent_smi_hold)
            smi_to_be_modified = [reactant_smi_hold, reagent_smi_hold]
        else:
            smi_to_be_modified = original_smi.split('>')
        smi_to_be_modified = ''.join((smi_to_be_modified[0], '>', smi_to_be_modified[1], '>', product_smi))
        """
    return smi_to_be_modified

