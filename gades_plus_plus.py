#!/usr/bin/env python3

import requests, json, os, rdflib, datetime
from sem_sim_list import *
from statistics import mean


ALPHA = 0.6

MATRIX_FILE = "patient_similarity.txt"
MATRIX_FILE_TO_PLOT = "to_plot.txt"
ENTITIES_FILE = "entities.txt"
LISTS_FILE = "list.txt"

SIM_OUTPUT = 'tmp/sim_output'
SIM_INPUT = 'tmp/sim_input'

DICTIONARY = "input_dictionary.txt"
DICTIONARY_QUERY = "query_dictionary.txt"
DICTIONARY_CONSTRUCT = "construct_dictionary.txt"

class Entity:    
    ids = None
    monthsInTreatment = 0

    ecog = []  
    biopsy = [] 
    stage = [] 
    drugs = [] 
    comorbidities = [] 
    chemotherapyList = [] 
    tkiList = [] 
    immunotherapyList = []
    familial_antecedent_list = []
    non_oncologycal_treatme = []

    immunotherapy = False 
    chemotherapy = False 
    tki = False 
    surgery = False 
    antiangiogenic = False 
    radiationtherapy = False
    familial_antecedent = False
    
    systemicProgression = False 
    localProgression = False 
    brainMetastasis = False 

    def __init__(self, id):
        self.ids = id
        self.ecog = []
        self.biopsy = []
        self.stage = []
        self.drugs = []
        self.comorbidities = []
        self.chemotherapyList = []
        self.tkiList = []
        self.immunotherapyList = []
        self.familial_antecedent_list = []
        self.non_oncologycal_treatme = []
        
    def __str__(self):
        return "++++++++++++\nPatient: "+self.ids+"\nComorbidities: "+str(self.comorbidities)+"\nDrugs: "+str(self.drugs)+"\nChemotherapy: "+str(self.chemotherapy)+"\nChemotherapy: "+str(self.chemotherapyList)+"\ntki: "+str(self.tki)+"\ntkiList: "+str(self.tkiList)+"\nImmunotherapy: "+str(self.immunotherapy)+"\nImmunotherapyList: "+str(self.immunotherapyList)+"\nAntiangiogenic: "+str(self.antiangiogenic)+"\nradiationtherapy : "+str(self.radiationtherapy)+"\nSurgery :"+str(self.surgery)+"\nSystemicProgression: "+str(self.systemicProgression)+"\nLocalProgression: "+str(self.localProgression)+"\nBrainMetastasis: "+str(self.brainMetastasis)+"\necog: "+str(self.ecog)+"\nStage: "+str(self.stage)+"\nBiopsy: "+str(self.biopsy)+"\nMonthsInTreatment: "+str(self.monthsInTreatment)

def sim_dates(d1, d2):
    d1 = abs(d1)
    d2 = abs(d2)
    return 1 - abs(d1 - d2)/max(d1, d2)

def get_entity_cancer_drugs(e):
    drugs = set()
    for c in e.chemotherapyList:
        drugs.add(c)
    for c in e.immunotherapyList:
        drugs.add(c)
    for c in e.tkiList:
        drugs.add(c)
    for c in e.non_oncologycal_treatme:
        drugs.add(c)
    return drugs

def get_treatments_set(e):
    l = []
    if e.immunotherapy:
        l.append("immunotherapy")
    if e.chemotherapy:
        l.append("chemotherapy")
    if e.tki:
        l.append("tki")
    if e.surgery:
        l.append("surgery")
    if e.antiangiogenic:
        l.append("antiangiogenic")
    if e.radiationtherapy:
        l.append("radiationtherapy")
    return l

def get_progression_set(e):
    l = []
    if e.systemicProgression:
        l.append("systemicProgression")
    if e.localProgression:
        l.append("localProgression")
    if e.brainMetastasis:
        l.append("brainMetastasis")
    return l

def get_antecedent_set(e):
    antecedent = set()
    for a in e.familial_antecedent_list:
        antecedent.add(a)
    return antecedent

def similarity_with_set(e1, e2):
    if e1.ids == e2.ids:
        return 1.0

    if len(e1.ecog) != 0 and len(e2.ecog) != 0:
        if len(e1.ecog) != 1 and len(e2.ecog) != 1:
            ecog = sim_subseq(e1.ecog, e2.ecog)
        else:
            if len(e1.ecog) == 1 and len(e2.ecog) == 1:
                #print("Warning using soft tf-idf", (e1.ecog, e2.ecog))
                ecog = jaro_winkler(e1.ecog, e2.ecog)
            else:
                ecog = 0.0
    else:
        ecog = 0.0

    if len(e1.stage) != 0 and len(e2.stage) != 0:
        if len(e1.stage) != 1 and len(e2.stage) != 1:
            stage = sim_subseq(e1.stage, e2.stage)
        else:
            if len(e1.stage) == 1 and len(e2.stage) == 1:
                #print("Warning using jaro winkler", (e1.stage, e2.stage))
                stage = jaro_winkler(e1.stage, e2.stage)
            else:
                stage = 0.0
    else:
        stage = 0.0
        
    if len(e1.biopsy) != 0 and len(e2.biopsy) != 0:
        if len(e1.biopsy) != 1 and len(e2.biopsy) != 1:
            biopsy = sim_subseq(e1.biopsy, e2.biopsy)
        else:
            if len(e1.biopsy) == 1 and len(e2.biopsy) == 1:
                #print("Warning using jaro winkler", (e1.biopsy, e2.biopsy))
                #print("biopsy")
                biopsy = sim_jaccard(e1.biopsy, e2.biopsy)
            else:
                biopsy = 0.0
    else:
        biopsy = 0.0

    if e1.monthsInTreatment > 0.0 or e2.monthsInTreatment > 0.0:
        months = sim_dates(e1.monthsInTreatment, e2.monthsInTreatment)
    else:
        months = 0

    if (len(e1.comorbidities) > 0) and  (len(e2.comorbidities) > 0):
        #print("comorbidities")
        sim_pair_comorbidities = sim_jaccard(e1.comorbidities, e2.comorbidities)
    else:
        sim_pair_comorbidities = 0.0

    if (len(e1.drugs) > 0) and  (len(e2.drugs) > 0):
        #print("drugs")
        sim_pair_drugs = sim_jaccard(e1.drugs, e2.drugs)
    else:
        sim_pair_drugs = 0.0

    cancer_drugs_e1 = get_entity_cancer_drugs(e1)
    cancer_drugs_e2 = get_entity_cancer_drugs(e2)

    if (len(cancer_drugs_e1) > 0) and (len(cancer_drugs_e2) > 0):
        #print("cancer_drugs")
        sim_pair_drugs_c = sim_jaccard(cancer_drugs_e1, cancer_drugs_e2)
    else:
        sim_pair_drugs_c = 0.0

    treatments_e1 = get_treatments_set(e1)
    treatments_e2 = get_treatments_set(e2)
    #print("treatments")
    if (len(treatments_e1) > 0) and (len(treatments_e2) > 0):
        treatments_pair_sim = sim_jaccard(treatments_e1, treatments_e2)
    else:
        treatments_pair_sim = 0.0

    #print(e1)
    progression_e1 = get_progression_set(e1)
    progression_e2 = get_progression_set(e2)
    if (len(progression_e1) > 0) and (len(progression_e2) > 0):
        progression_pair_sim = sim_jaccard(progression_e1, progression_e2)
    else:
        progression_pair_sim = 0.0

    antecedent_e1 = get_antecedent_set(e1)
    antecedent_e2 = get_antecedent_set(e2)
    if (len(antecedent_e1) > 0) and (len(antecedent_e2) > 0):
        sim_pair_antecedent_f = sim_jaccard(antecedent_e1, antecedent_e2)
    else:
        sim_pair_antecedent_f = 0.0

    #evolution = [ecog, stage, biopsy, sim_pair_comorbidities, sim_pair_drugs, sim_pair_drugs_c, treatments_pair_sim, progression_pair_sim, sim_pair_antecedent_f]
    evolution = [ecog, stage, biopsy, sim_pair_comorbidities, sim_pair_drugs_c, treatments_pair_sim,
                 progression_pair_sim, sim_pair_antecedent_f]
    sim = (((1.0-ALPHA)*mean(evolution)) + months*ALPHA)/2.0
    if sim < 0:
        print(evolution)
        print(months)
        print(e1.monthsInTreatment)
        print(e2.monthsInTreatment)
        return -1
    return sim


############################################################################


def compute_all_similarities_with_set(d_entities):
    #patient_id_list = list(d_entities.keys())

    with open(MATRIX_FILE, "w") as fm:
        with open(MATRIX_FILE_TO_PLOT, "w") as fp:
            fm.write("patient_id1" + "," + "patient_id2" + "," + "similarityValue" + "\n")
            fp.write("patient_id1" + "," + "patient_id2" + "," + "similarityValue" + "\n")
            results = []
            cont = 0
            n = len(d_entities)
            for i in range(n):
                e1 = d_entities[i]
                for j in range(i, n):
                    if cont > 0 and cont % 10 == 0:
                        print("\nNumber of pairs so far: " + str(cont))
                    e2 = d_entities[j]
                    """""
                    print("P1")
                    print(e1.ids)
                    print(e1.familial_antecedent_list)
                    print("P2")
                    print(e2.ids)
                    print(e2.familial_antecedent_list)
                    """
                    sim = similarity_with_set(e1, e2)
                    if sim < 0:
                        print(e1.ids)
                        print(e2.ids)
                    results.append((e1.ids, e2.ids, sim))
                    cont += 1
                    fm.write(str(e1.ids) + "," + str(e2.ids) + "," + str(sim) + "\n")
                    if e1.ids != e2.ids:
                        fp.write(str(e1.ids) + "," + str(e2.ids) + "," + str(sim) + "\n")
        return results


def load_patients_from_file(filename):
    with open(filename) as fd:
        header = fd.readline().rstrip().split(",")
        #print(header)
        #print("Number of concepts: "+str(len(header)))
        d_entities = []
        for line in fd:
            tok = line.rstrip().split(",")
            #assert( len(header) == len(tok) )
            if len(header) != len(tok):
                #print(len(header), len(tok))
                sys.exit(1)
            entity = Entity(tok[0])
            entity.chemotherapy = True if tok[12] == "True" else False
            entity.tki = True if tok[13] == "True" else False
            entity.immunotherapy = True if tok[14] == "True" else False
            entity.antiangiogenic = True if tok[15] == "True" else False
            entity.radiationtherapy = True if tok[16] == "True" else False
            entity.surgery = True if tok[17] == "True" else False
            entity.familial_antecedent = True if tok[18] == "True" else False
            entity.systemicProgression = True if tok[19] == "True" else False
            entity.localProgression = True if tok[20] == "True" else False
            entity.brainMetastasis = True if tok[21] == "True" else False
            entity.monthsInTreatment = int(tok[22])

            d_entities.append(entity)
        return d_entities

def get_chemotherapy_List(pat_chemotherapy, chemotherapy, d_entities):
    p_ch = open(pat_chemotherapy)
    chemoth_list = open(chemotherapy)
    length_p_ch = len(p_ch.readlines())
    length_chemoth_list = len(chemoth_list.readlines())
    for i in range(len(d_entities)):
        p_ch.seek(1)
        p_ch.readline()
        if (d_entities[i].chemotherapy == False):
            continue
        for j in range(length_p_ch):
            tok = p_ch.readline().rstrip().split(",")
            if d_entities[i].ids == tok[0]:
                chemoth_list.seek(1)
                chemoth_list.readline()
                for k in range(length_chemoth_list):
                    tok_k = chemoth_list.readline().rstrip().split(",")
                    if tok[1] == tok_k[1]:
                        d_entities[i].chemotherapyList.append(tok[1])
                        break
    p_ch.close()
    chemoth_list.close()

def get_familial_antecedents_List(pat_familial_antecedents, familial_antecedents, d_entities):
    p_familial_antecedents = open(pat_familial_antecedents)
    familial_antecedents_list = open(familial_antecedents)
    length_p_familial_antecedents = len(p_familial_antecedents.readlines())
    length_familial_antecedents_list = len(familial_antecedents_list.readlines())
    #print("length_p_familial_antecedents")
    #print(length_p_familial_antecedents)
    #print("length_familial_antecedents_list")
    #print(length_familial_antecedents_list)
    for i in range(len(d_entities)):
        p_familial_antecedents.seek(1)
        p_familial_antecedents.readline()
        if(d_entities[i].familial_antecedent == False):
            continue
        for j in range(length_p_familial_antecedents):
            tok = p_familial_antecedents.readline().rstrip().split(",")
            #print("p_familial_antecedents")
            #print(tok)
            if d_entities[i].ids == tok[0]:
                #print("d_entities[i].ids == tok[0]")
                #print(d_entities[i].ids)
                familial_antecedents_list.seek(1)
                familial_antecedents_list.readline()
                for k in range(length_familial_antecedents_list):
                    tok_k = familial_antecedents_list.readline().rstrip().split(",")
                    #print("familial_antecedents_list")
                    #print(tok_k)
                    if tok[1] == tok_k[1]:
                        d_entities[i].familial_antecedent_list.append(tok[1])
                        #print("antecedent ADD")
                        #print(tok[1])
                        break
        #print(d_entities[i].ids)
        #print(d_entities[i].familial_antecedent_list)
    p_familial_antecedents.close()
    familial_antecedents_list.close()

def get_immunotherapy_List(pat_immunotherapy, immunotherapy, d_entities):
    p_immun = open(pat_immunotherapy)
    immunoth_list = open(immunotherapy)
    length_p_immun = len(p_immun.readlines())
    length_immunoth_list = len(immunoth_list.readlines())
    for i in range(len(d_entities)):
        p_immun.seek(1)
        p_immun.readline()
        if (d_entities[i].immunotherapy == False):
            continue
        for j in range(length_p_immun):
            tok = p_immun.readline().rstrip().split(",")
            if d_entities[i].ids == tok[0]:
                immunoth_list.seek(1)
                immunoth_list.readline()
                for k in range(length_immunoth_list):
                    tok_k = immunoth_list.readline().rstrip().split(",")
                    if tok[1] == tok_k[1]:
                        d_entities[i].immunotherapyList.append(tok[1])
                        break
    p_immun.close()
    immunoth_list.close()

def get_tki_List(pat_tki, tki, d_entities):
    p_tki = open(pat_tki)
    tki_list = open(tki)
    length_p_tki = len(p_tki.readlines())
    length_tki_list = len(tki_list.readlines())
    for i in range(len(d_entities)):
        p_tki.seek(1)
        p_tki.readline()
        if (d_entities[i].tki == False):
            continue
        for j in range(length_p_tki):
            tok = p_tki.readline().rstrip().split(",")
            if d_entities[i].ids == tok[0]:
                tki_list.seek(1)
                tki_list.readline()
                for k in range(length_tki_list):
                    tok_k = tki_list.readline().rstrip().split(",")
                    if tok[1] == tok_k[1]:
                        d_entities[i].tkiList.append(tok[1])
                        break
    p_tki.close()
    tki_list.close()

def get_comorbidity_list(pat_comorbidities, comorbidities, d_entities):
    p_comorbidities = open(pat_comorbidities)
    comorbidities_list = open(comorbidities)
    length_p_comorbidities = len(p_comorbidities.readlines())
    length_comorbidities_list = len(comorbidities_list.readlines())
    for i in range(len(d_entities)):
        p_comorbidities.seek(1)
        p_comorbidities.readline()
        for j in range(length_p_comorbidities):
            tok = p_comorbidities.readline().rstrip().split(",")
            if d_entities[i].ids == tok[0]:
                comorbidities_list.seek(1)
                comorbidities_list.readline()
                for k in range(length_comorbidities_list):
                    tok_k = comorbidities_list.readline().rstrip().split(",")
                    if tok[1] == tok_k[1]:
                        d_entities[i].comorbidities.append(tok[1])
                        break
    p_comorbidities.close()
    comorbidities_list.close()

def get_drugs_list(pat_drug_groups, d_entities):
    p_drug_groups = open(pat_drug_groups)
    length_p_drug_groups = len(p_drug_groups.readlines())
    for i in range(len(d_entities)):
        p_drug_groups.seek(1)
        p_drug_groups.readline()
        for j in range(length_p_drug_groups):
            tok = p_drug_groups.readline().rstrip().split(",")
            if d_entities[i].ids == tok[0]:
                d_entities[i].drugs.append(tok[1])
                #print("drug:")
                #print(tok[1])
    p_drug_groups.close()

def get_non_oncologycal_treatment(pat_oncologycal_treatment, oncologycal_treatment, d_entities):
    p_oncologycal_treatment = open(pat_oncologycal_treatment)
    oncologycal_treatment_list = open(oncologycal_treatment)
    length_p_oncologycal_treatment = len(p_oncologycal_treatment.readlines())
    length_oncologycal_treatment_list = len(oncologycal_treatment_list.readlines())
    for i in range(len(d_entities)):
        p_oncologycal_treatment.seek(1)
        p_oncologycal_treatment.readline()
        for j in range(length_p_oncologycal_treatment):
            tok = p_oncologycal_treatment.readline().rstrip().split(",")
            if d_entities[i].ids == tok[0]:
                oncologycal_treatment_list.seek(1)
                oncologycal_treatment_list.readline()
                for k in range(length_oncologycal_treatment_list):
                    tok_k = oncologycal_treatment_list.readline().rstrip().split(",")
                    if tok[1] == tok_k[1]:
                        d_entities[i].non_oncologycal_treatme.append(tok[1])
                        break
    p_oncologycal_treatment.close()
    oncologycal_treatment_list.close()

def get_biopsy_list(pat_biopsia, d_entities):
    p_biopsia = open(pat_biopsia)
    length_p_biopsia = len(p_biopsia.readlines())
    d_list = {}
    for i in range(len(d_entities)):
        d_list[i] = []
        p_biopsia.seek(1)
        p_biopsia.readline()
        for j in range(length_p_biopsia):
            tok = p_biopsia.readline().rstrip().split(",")
            if d_entities[i].ids == tok[0]:
                d_list[i].append({'date': tok[2], 'value': tok[4]})

    for p in d_list:
        d_list[p].sort(key=lambda x: datetime.datetime.strptime(x['date'], '%Y-%m-%d'))
        for x in d_list[p]:
            d_entities[p].biopsy.append(x['value'])
    p_biopsia.close()

def get_performance_status_list(pat_performance_status, d_entities):
    p_performance_status = open(pat_performance_status)
    length_p_performance_status = len(p_performance_status.readlines())
    d_ecog_list = {}
    for i in range(len(d_entities)):
        d_ecog_list[i] = []
        p_performance_status.seek(1)
        p_performance_status.readline()
        for j in range(length_p_performance_status):
            tok = p_performance_status.readline().rstrip().split(",")
            if d_entities[i].ids == tok[0]:
                d_ecog_list[i].append({'date': tok[3], 'value': tok[2]})

    for p in d_ecog_list:
        d_ecog_list[p].sort(key=lambda x: datetime.datetime.strptime(x['date'], '%Y-%m-%d'))
        for x in d_ecog_list[p]:
            d_entities[p].ecog.append(x['value'])
    p_performance_status.close()

def get_stage_list(pat_stage, d_entities):
    p_stage = open(pat_stage)
    length_p_stage = len(p_stage.readlines())
    d_stage_list = {}
    for i in range(len(d_entities)):
        d_stage_list[i] = []
        p_stage.seek(1)
        p_stage.readline()
        for j in range(length_p_stage):
            tok = p_stage.readline().rstrip().split(",")
            if d_entities[i].ids == tok[0]:
                d_stage_list[i].append({'date': tok[3], 'value': tok[2]})

    for p in d_stage_list:
        d_stage_list[p].sort(key=lambda x: datetime.datetime.strptime(x['date'], '%Y-%m-%d'))
        for x in d_stage_list[p]:
            d_entities[p].stage.append(x['value'])
    p_stage.close()

def save_list(d_entities):
    with open(LISTS_FILE, "w") as fm:
        n = len(d_entities)
        fm.write("ids" + "\t" + "chemotherapyList" + "\t" + "immunotherapyList" + "\t" + "tkiList" + "\t"
                 + "comorbidities" + "\t" + "drugs" + "\t" + "biopsy" + "\t" + "ecog" + "\t" +"stage" + "\t" +
                 "familial_antecedent_list" + "\t" + "non_oncologycal_treatme" + "\n")
        for i in range(n):
            e1 = d_entities[i]
            fm.write(str(e1.ids) + "\t" + str(e1.chemotherapyList) + "\t" + str(e1.immunotherapyList) + "\t" + str(
                e1.tkiList)+ "\t" + str(e1.comorbidities) + "\t" + str(e1.drugs) + "\t" + str(e1.biopsy) + "\t" +
                     str(e1.ecog) + "\t" +str(e1.stage) + "\t" + str(e1.familial_antecedent_list) + "\t" +
                     str(e1.non_oncologycal_treatme) + "\n")

def create_filter_dictionary ():
    with open(DICTIONARY, "w") as fd:
        fd.write("gender" + "," + "<http://project-iasis.eu/Gender" +"\n")
        fd.write("age" + "," + "http://www.w3.org/2001/XMLSchema#" + "\n")
        fd.write("tumorStage" + "," + "<http://project-iasis.eu/LCPatient/TumorStage" + "\n")
        fd.write("immunotherapy" + "," + "<http://project-iasis.eu/vocab/immunotherapy" + "\n")
        fd.write("chemotherapy" + "," + "<http://project-iasis.eu/vocab/chemotherapy" + "\n")
        fd.write("tki" + "," + "<http://project-iasis.eu/vocab/tki" + "\n")
        fd.write("surgery" + "," + "<http://project-iasis.eu/vocab/surgery" + "\n")
        fd.write("antiangiogentic" + "," + "<http://project-iasis.eu/vocab/antiangiogentic" + "\n")
        fd.write("radiationtherapy" + "," + "<http://project-iasis.eu/vocab/radiationtherapy" + "\n")
        fd.write("familiarAntecedents" + "," + "<http://project-iasis.eu/vocab/familiarAntecedents" + "\n")
        fd.write("systemicProgression" + "," + "<http://project-iasis.eu/vocab/systemicProgression" + "\n")
        fd.write("localProgression" + "," + "<http://project-iasis.eu/vocab/localProgression" + "\n")
        fd.write("brainMetastasis" + "," + "<http://project-iasis.eu/vocab/brainMetastasis" + "\n")
        fd.write("survivalMonths" + "," + "<http://project-iasis.eu/vocab/survivalMonths")

def create_query_dictionary ():
    with open(DICTIONARY_QUERY, "w") as fd:
        fd.write("gender" + "," + "?p1 <http://project-iasis.eu/vocab/gender> ?gender."+"\n")
        fd.write("age" + "," + "?p1 <http://project-iasis.eu/vocab/hasDiagnosis> ?d1.***?d1 <http://project-iasis.eu/vocab/hasDiagnosisAge> ?age." + "\n")
        fd.write("tumorStage" + "," + "?p1 <http://project-iasis.eu/vocab/hasDiagnosis> ?d1.***?d1 <http://project-iasis.eu/vocab/hasDiagnosisTumorStage> ?tumorStage." + "\n")
        fd.write("performanceStatus" + "," + "?p1 <http://project-iasis.eu/vocab/hasDiagnosis> ?d1.***?d1 <http://project-iasis.eu/vocab/hasDiagnosisPerformanceStatus> ?PerformanceStatus." + "\n")
        fd.write("date" + "," + "?p1 <http://project-iasis.eu/vocab/hasDiagnosis> ?d1.***?d1 <http://project-iasis.eu/vocab/hasDiagnosisDate> ?date." + "\n")
        fd.write("immunotherapy" + "," + "?p1 <http://project-iasis.eu/vocab/immunotherapy> ?immunotherapy." + "\n")
        fd.write("chemotherapy" + "," + "?p1 <http://project-iasis.eu/vocab/chemotherapy> ?chemotherapy." + "\n")
        fd.write("tki" + "," + "?p1 <http://project-iasis.eu/vocab/tki> ?tki." + "\n")
        fd.write("surgery" + "," + "?p1 <http://project-iasis.eu/vocab/surgery> ?surgery." + "\n")
        fd.write("antiangiogentic" + "," + "?p1 <http://project-iasis.eu/vocab/antiangiogentic> ?antiangiogentic." + "\n")
        fd.write("radiationtherapy" + "," + "?p1 <http://project-iasis.eu/vocab/radiationtherapy> ?radiationtherapy." + "\n")
        fd.write("familiarAntecedents" + "," + "?p1 <http://project-iasis.eu/vocab/familiarAntecedents> ?familiarAntecedents." + "\n")
        fd.write("systemicProgression" + "," + "?p1 <http://project-iasis.eu/vocab/systemicProgression> ?systemicProgression." + "\n")
        fd.write("localProgression" + "," + "?p1 <http://project-iasis.eu/vocab/localProgression> ?localProgression." + "\n")
        fd.write("brainMetastasis" + "," + "?p1 <http://project-iasis.eu/vocab/brainMetastasis> ?brainMetastasis." + "\n")
        fd.write("survivalMonths" + "," + "?p1 <http://project-iasis.eu/vocab/survivalMonths> ?survivalMonths.")

def create_construct_dictionary ():
    with open(DICTIONARY_CONSTRUCT, "w") as fd:
        fd.write("gender" + "," + "?p1 <http://project-iasis.eu/vocab/gender> ?gender."+"\n")
        fd.write("age" + "," + "?d1 <http://project-iasis.eu/vocab/hasDiagnosisAge> ?age." + "\n")
        fd.write("tumorStage" + "," + "?d1 <http://project-iasis.eu/vocab/hasDiagnosisTumorStage> ?tumorStage." + "\n")
        fd.write("performanceStatus" + "," + "?d1 <http://project-iasis.eu/vocab/hasDiagnosisPerformanceStatus> ?PerformanceStatus." + "\n")
        fd.write("date" + "," + "?d1 <http://project-iasis.eu/vocab/hasDiagnosisDate> ?date." + "\n")
        fd.write("immunotherapy" + "," + "?p1 <http://project-iasis.eu/vocab/immunotherapy> ?immunotherapy." + "\n")
        fd.write("chemotherapy" + "," + "?p1 <http://project-iasis.eu/vocab/chemotherapy> ?chemotherapy." + "\n")
        fd.write("tki" + "," + "?p1 <http://project-iasis.eu/vocab/tki> ?tki." + "\n")
        fd.write("surgery" + "," + "?p1 <http://project-iasis.eu/vocab/surgery> ?surgery." + "\n")
        fd.write("antiangiogentic" + "," + "?p1 <http://project-iasis.eu/vocab/antiangiogentic> ?antiangiogentic." + "\n")
        fd.write("radiationtherapy" + "," + "?p1 <http://project-iasis.eu/vocab/radiationtherapy> ?radiationtherapy." + "\n")
        fd.write("familiarAntecedents" + "," + "?p1 <http://project-iasis.eu/vocab/familiarAntecedents> ?familiarAntecedents." + "\n")
        fd.write("systemicProgression" + "," + "?p1 <http://project-iasis.eu/vocab/systemicProgression> ?systemicProgression." + "\n")
        fd.write("localProgression" + "," + "?p1 <http://project-iasis.eu/vocab/localProgression> ?localProgression." + "\n")
        fd.write("brainMetastasis" + "," + "?p1 <http://project-iasis.eu/vocab/brainMetastasis> ?brainMetastasis." + "\n")
        fd.write("survivalMonths" + "," + "?p1 <http://project-iasis.eu/vocab/survivalMonths> ?survivalMonths.")

def main(*args):

    d_entities = load_patients_from_file(args[0])
    get_chemotherapy_List(args[1], args[2], d_entities)
    get_immunotherapy_List(args[3], args[4], d_entities)
    get_tki_List(args[5], args[6], d_entities)
    get_comorbidity_list(args[7], args[8], d_entities)
    get_drugs_list(args[9], d_entities)
    get_biopsy_list(args[10], d_entities)
    get_performance_status_list(args[11], d_entities)
    get_stage_list(args[12], d_entities)
    get_familial_antecedents_List(args[13], args[14], d_entities)
    get_non_oncologycal_treatment(args[15], args[16], d_entities)

    save_list(d_entities)

    similarities = compute_all_similarities_with_set(d_entities)
    print(similarities)
    print("Building the semEP-Node input files")
    #create_semEP_node_input(similarities)
    create_filter_dictionary()
    create_query_dictionary()
    create_construct_dictionary()
    print("Done Gades Plus Plus")

if __name__ == "__main__":
    main(*sys.argv[1:])
