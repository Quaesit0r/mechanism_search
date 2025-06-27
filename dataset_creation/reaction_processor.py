import traceback
import time
import reaction
import json
import jsonlines
import os
import statistic_module


class path_holder:
    def __init__(self, num):
        self.all_rxn_name = r'database/database_%d.json' % num
        self.destination_stat = r'Data/Stat_%d.txt' % num
        self.mech_database = r'mech_db/mechanism_database_%d.json' % num
        self.good_match_file = r'Data/Good_Match_%d.json' % num
        self.potential_match_file = r'Data/Potential_Match_%d.json' % num
        self.no_match_file = r'Data/No_Match_%d.json' % num

    def return_all_rxn_name(self):
        return self.all_rxn_name

    def return_destination_stat(self):
        return self.destination_stat

    def return_mech_database(self):
        return self.mech_database

    def return_good_match_file(self):
        return self.good_match_file

    def return_potential_match_file(self):
        return self.potential_match_file

    def return_no_match_file(self):
        return self.no_match_file


class mech_writer:
    def __init__(self, num, step):
        self.good_match_holder = r'Data/Good_Match_%d.json' % num
        self.potential_match_holder = r'Data/Potential_Match_%d_step_%d.json' % (num, step)
        self.no_match_holder = r'Data/No_Match_%d.json' % num
        self.good_match_file = None
        self.potential_match_file = None
        self.no_match_file = None

    def initialize(self):
        self.good_match_file = open(self.good_match_holder, 'a')
        self.potential_match_file = open(self.potential_match_holder, 'a')
        self.no_match_file = open(self.no_match_holder, 'a')

    def close(self):
        self.good_match_file.close()
        self.potential_match_file.close()
        self.no_match_file.close()

    def write_to_good_match_holder(self, rxn_node):
        rxn_dict = reaction.node_encoder(rxn_node)
        json.dump(rxn_dict, self.good_match_file)
        self.good_match_file.write('\n')

    def write_to_potential_match_holder(self, rxn_node):
        rxn_dict = reaction.node_encoder(rxn_node)
        json.dump(rxn_dict, self.potential_match_file)
        self.potential_match_file.write('\n')

    def write_to_no_match_holder(self, rxn_node):
        rxn_dict = reaction.node_encoder(rxn_node)
        json.dump(rxn_dict, self.no_match_file)
        self.no_match_file.write('\n')



def process_reaction(num, start_step, end_step, radius):
    path = path_holder(num)
    statistic = statistic_module.statistics()
    template_obj = reaction.template_list(path.return_mech_database())
    time_st = time.time()

    try:
        for search_step in range(start_step, end_step+1):
            writer = mech_writer(num, search_step)
            writer.initialize()
            # Initialization of reaction nodes
            if not os.path.exists(path.return_all_rxn_name()):
                break
            with jsonlines.open(path.return_all_rxn_name()) as all_rxn_dict:
                for rxn_dict in all_rxn_dict:
                    rxn_node = reaction.node_decoder(rxn_dict, radius)
                    if not rxn_node.reaction_filter():
                        writer.write_to_no_match_holder(rxn_node)
                        statistic.add_no_match()
                        continue
                    match = reaction.apply_all_mechanism(rxn_node, template_obj.get_mechanism(), search_step)
                    if match:
                        statistic.determining_match(match)
                        reaction.good_match_processor(match)
                        writer.write_to_good_match_holder(rxn_node)
                    elif rxn_node.check_alive():
                        writer.write_to_potential_match_holder(rxn_node)
                        statistic.add_potential_match()
                    else:
                        writer.write_to_no_match_holder(rxn_node)
                        statistic.add_no_match()

            writer.close()
            statistic.write_statistics(path.return_destination_stat(), time_st)
            statistic.reset()
            if os.path.exists(path.all_rxn_name):
                os.remove(path.all_rxn_name)
            path.all_rxn_name = writer.potential_match_holder


    except Exception as e:
        print({e})
        traceback.print_exc()
        statistic.write_statistics(path.return_destination_stat(), time_st)

    if os.path.exists(r'Data/Potential_Match_%d_step_%d.json' % (num, end_step)):
        os.rename(r'Data/Potential_Match_%d_step_%d.json' % (num, end_step), r'Data/Potential_Match_%d.json' % num)

