import time
import os
import re
import json
import jsonlines

class statistics:
    def __init__(self):
        self.max_time = 0
        self.tot_num = 0
        self.num_match = 0
        self.num_no_match = 0
        self.num_potential_match = 0
        self.num_DA = 0
        self.num_SN2 = 0
        self.num_carbo_rearr = 0
        self.num_3_3_sigma = 0
        self.num_1_5_sigma = 0
        self.num_1_7_sigma = 0
        self.num_Prot = 0
        self.num_Deprot = 0
        self.num_Nuc_add = 0
        self.num_Ele_add = 0
        self.num_Het_cleav = 0
        self.num_E2 = 0
        self.num_Elimination = 0
        self.num_Dissociation = 0
        self.num_SnAr = 0
        self.num_Isomerization = 0

    def determining_match(self, good_match_node):
        template_name = good_match_node.mechanism_name
        if template_name == ('Diels_Alder_Reaction'):
            self.num_DA += 1
        elif template_name == ('SN2_Reaction'):
            self.num_SN2 += 1
        elif template_name == ('Carbocation_rearrangement'):
            self.num_carbo_rearr += 1
        elif template_name == ('3,3_sigmatropic_rearrangement'):
            self.num_3_3_sigma += 1
        elif template_name == ('1,5_sigmatropic_rearrangement'):
            self.num_1_5_sigma += 1
        elif template_name == ('1,7_sigmatropic_rearrangement'):
            self.num_1_7_sigma += 1
        elif template_name == ('Protonation'):
            self.num_Prot += 1
        elif template_name == ('Deprotonation'):
            self.num_Deprot += 1
        elif template_name == ('Nucleophilic_addition'):
            self.num_Nuc_add += 1
        elif template_name == ('Electrophilic_addition'):
            self.num_Ele_add += 1
        elif template_name == ('Heterolytic_cleavage'):
            self.num_Het_cleav += 1
        elif template_name == ('E2_Reaction'):
            self.num_E2 += 1
        elif template_name == ('Elimination_reaction'):
            self.num_Elimination += 1
        elif template_name == ('Dissociation'):
            self.num_Dissociation += 1
        elif template_name == ('SnAr_reaction'):
            self.num_SnAr += 1
        elif template_name == ('Isomerization'):
            self.num_Isomerization += 1
        else:
            print('Unexpected good match')
        self.num_match += 1

    def write_statistics(self, Stat, time_st):
        with open(Stat, "a") as stat:
            time_ck = time.time()
            tot_time = time_ck - time_st
            total_num = self.num_match + self.num_no_match + self.num_potential_match
            stat.write("The total number of reactions examined is %d" % total_num+'\n')
            stat.write("The number of good match is %d" % self.num_match+"\n")
            stat.write("The number of no match is %d" % self.num_no_match+"\n")
            stat.write("The number of potential match is %d" % self.num_potential_match+"\n"+"\n")
            stat.write("Mechanism hits:"+"\n")
            stat.write("The number of Diels-Alder reaction hit is %d" % self.num_DA+"\n")
            stat.write("The number of SN2 reaction hit is %d" % self.num_SN2+"\n")
            stat.write("The number of carbocation rearrangement reaction hit is %d" % self.num_carbo_rearr + "\n")
            stat.write("The number of 3,3-sigmatropic rearrangement reaction hit is %d" % self.num_3_3_sigma+"\n")
            stat.write("The number of 1,5-sigmatropic rearrangement reaction hit is %d" % self.num_1_5_sigma+"\n")
            stat.write("The number of 1,7-sigmatropic rearrangement reaction hit is %d" % self.num_1_7_sigma+"\n")
            stat.write("The number of protonation reaction hit is %d" % self.num_Prot+"\n")
            stat.write("The number of deprotonation reaction hit is %d" % self.num_Deprot+"\n")
            stat.write("The number of nucleophilic reaction hit is %d" % self.num_Nuc_add+"\n")
            stat.write("The number of electrophilic reaction hit is %d" % self.num_Ele_add+"\n")
            stat.write("The number of heterolytic cleavage reaction hit is %d" % self.num_Het_cleav+"\n")
            stat.write("The number of E2 reaction hit is %d" % self.num_E2+"\n")
            stat.write("The number of elimination reaction hit is %d" % self.num_Elimination+"\n")
            stat.write("The number of dissociation reaction hit is %d" % self.num_Dissociation+"\n")
            stat.write("The number of SnAr reaction hit is %d" % self.num_SnAr+"\n")
            stat.write("The number of isomerization reaction hit is %d" % self.num_Isomerization+"\n"+"\n")
            stat.write("Runtime is %f" % tot_time+"\n"+"\n"+"\n")

    def add_potential_match(self):
        self.num_potential_match += 1

    def add_no_match(self):
        self.num_no_match += 1

    def reset(self):
        self.max_time = 0
        self.tot_num = 0
        self.num_match = 0
        self.num_no_match = 0
        self.num_potential_match = 0
        self.num_DA = 0
        self.num_SN2 = 0
        self.num_carbo_rearr = 0
        self.num_3_3_sigma = 0
        self.num_1_5_sigma = 0
        self.num_1_7_sigma = 0
        self.num_Prot = 0
        self.num_Deprot = 0
        self.num_Nuc_add = 0
        self.num_Ele_add = 0
        self.num_Het_cleav = 0
        self.num_E2 = 0
        self.num_Elimination = 0
        self.num_Dissociation = 0
        self.num_SnAr = 0
        self.num_Isomerization = 0


def calculate_total_statistics(total_num_database):
    finalized_stat = statistics()
    for database_number in range(total_num_database):
        try:
            temp_potential_match = 0
            with open(r'Data/Stat_%d.txt' % database_number) as db:
                writer = db.readline()
                while writer:
                    if 'number of good match' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[0])
                        finalized_stat.num_match += temp_num
                    elif 'number of no match' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[0])
                        finalized_stat.num_no_match += temp_num
                    elif 'number of potential match' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[0])
                        temp_potential_match = temp_num
                    elif 'Diels-Alder reaction' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[0])
                        finalized_stat.num_DA += temp_num
                    elif 'SN2 reaction' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[1])
                        finalized_stat.num_SN2 += temp_num
                    elif 'carbocation rearrangement' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[0])
                        finalized_stat.num_carbo_rearr += temp_num
                    elif '3,3-sigmatropic rearrangement' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[2])
                        finalized_stat.num_3_3_sigma += temp_num
                    elif '1,5_sigmatropic rearrangement' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[2])
                        finalized_stat.num_1_5_sigma += temp_num
                    elif '1,7_sigmatropic rearrangement' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[2])
                        finalized_stat.num_1_7_sigma += temp_num
                    elif 'deprotonation' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[0])
                        finalized_stat.num_Deprot += temp_num
                    elif 'protonation' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[0])
                        finalized_stat.num_Prot += temp_num
                    elif 'nucleophilic reaction' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[0])
                        finalized_stat.num_Nuc_add += temp_num
                    elif 'electrophilic reaction' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[0])
                        finalized_stat.num_Ele_add += temp_num
                    elif 'heterolytic cleavage reaction' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[0])
                        finalized_stat.num_Het_cleav += temp_num
                    elif 'E2 reaction' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[1])
                        finalized_stat.num_E2 += temp_num
                    elif 'elimination reaction' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[0])
                        finalized_stat.num_Elimination += temp_num
                    elif 'dissociation reaction' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[0])
                        finalized_stat.num_Dissociation += temp_num
                    elif 'SnAr reaction' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[0])
                        finalized_stat.num_SnAr += temp_num
                    elif 'isomerization' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[0])
                        finalized_stat.num_Isomerization += temp_num
                    elif 'Runtime' in writer:
                        temp_num = re.findall(r'\d+', writer)
                        temp_num = int(temp_num[0])
                        if temp_num > finalized_stat.max_time:
                            finalized_stat.max_time = temp_num
                    writer = db.readline()

            finalized_stat.num_potential_match += temp_potential_match
            os.remove(r'Data/Stat_%d.txt' % database_number)

        except FileNotFoundError:
            print('Stat file not found, skip')

    finalized_stat.tot_num = finalized_stat.num_match + finalized_stat.num_potential_match + finalized_stat.num_no_match
    return finalized_stat

# TODO: Make this function more general
def finalize_statistics(num_parallel):

    for num in range(num_parallel):
        try:
            with jsonlines.open(r'Data/Good_Match_%d.json' % num) as db:
                with open(r'Data/Good_Match.json', 'a') as match:
                    for good_match_dict in db:
                        json.dump(good_match_dict, match)
                        match.write('\n')
            os.remove(r'Data/Good_Match_%d.json' % num)
        except FileNotFoundError:
            print('Good match file not found, skip')

        try:
            with jsonlines.open(r'Data/Potential_Match_%d.json' % num) as db:
                with open(r'Data/Potential_Match.json', 'a') as match:
                    for potential_match_dict in db:
                        json.dump(potential_match_dict, match)
                        match.write('\n')
            os.remove(r'Data/Potential_Match_%d.json' % num)
        except FileNotFoundError:
            print('Potential match file not found, skip')

        try:
            with jsonlines.open(r'Data/No_Match_%d.json' % num) as db:
                with open(r'Data/No_Match.json', 'a') as match:
                    for no_match_dict in db:
                        json.dump(no_match_dict, match)
                        match.write('\n')
            os.remove(r'Data/No_Match_%d.json' % num)
        except FileNotFoundError:
            print('No match file not found, skip')

        if os.path.exists(r'database/database_%d.json' % num):
            os.remove(r'database/database_%d.json' % num)
        os.remove(r'mech_db/mechanism_database_%d.json' % num)


def write_finalized_statistics(finalized_stat, Stat):
    with open(Stat, "a") as stat:
        tot_time = finalized_stat.max_time
        stat.write("The total number of reactions examined is %d" % finalized_stat.tot_num+'\n')
        stat.write("The number of good match is %d" % finalized_stat.num_match+"\n")
        stat.write("The number of no match is %d" % finalized_stat.num_no_match+"\n")
        stat.write("The number of potential match is %d" % finalized_stat.num_potential_match+"\n"+"\n")
        stat.write("Mechanism hits:"+"\n")
        stat.write("The number of Diels-Alder reaction hit is %d" % finalized_stat.num_DA+"\n")
        stat.write("The number of SN2 reaction hit is %d" % finalized_stat.num_SN2+"\n")
        stat.write("The number of carbocation rearrangement reaction hit is %d" % finalized_stat.num_carbo_rearr + "\n")
        stat.write("The number of 3,3-sigmatropic rearrangement reaction hit is %d" % finalized_stat.num_3_3_sigma+"\n")
        stat.write("The number of 1,5-sigmatropic rearrangement reaction hit is %d" % finalized_stat.num_1_5_sigma+"\n")
        stat.write("The number of 1,7-sigmatropic rearrangement reaction hit is %d" % finalized_stat.num_1_7_sigma+"\n")
        stat.write("The number of protonation reaction hit is %d" % finalized_stat.num_Prot+"\n")
        stat.write("The number of deprotonation reaction hit is %d" % finalized_stat.num_Deprot+"\n")
        stat.write("The number of nucleophilic reaction hit is %d" % finalized_stat.num_Nuc_add+"\n")
        stat.write("The number of electrophilic reaction hit is %d" % finalized_stat.num_Ele_add+"\n")
        stat.write("The number of heterolytic cleavage reaction hit is %d" % finalized_stat.num_Het_cleav+"\n")
        stat.write("The number of E2 reaction hit is %d" % finalized_stat.num_E2+"\n")
        stat.write("The number of elimination reaction hit is %d" % finalized_stat.num_Elimination+"\n")
        stat.write("The number of dissociation reaction hit is %d" % finalized_stat.num_Dissociation+"\n")
        stat.write("The number of SnAr reaction hit is %d" % finalized_stat.num_SnAr+"\n")
        stat.write("The number of isomerization reaction hit is %d" % finalized_stat.num_Isomerization+"\n"+"\n")
        stat.write("Runtime is %f" % tot_time+"\n"+"\n"+"\n")