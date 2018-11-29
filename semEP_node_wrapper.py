#!/usr/bin/env python3
#
# Description: Wrapper for SemEP-Node
#
# Copying:  MIT License
import glob
import sys, os, numpy



from PIL import Image
from SPARQLWrapper import SPARQLWrapper, JSON
from time import time
import json
from os import listdir
from os.path import isfile, join
from distributions import compare_distributions

import rpy2.robjects as ro

import rpy2.robjects.numpy2ri

rpy2.robjects.numpy2ri.activate()

############################
#
# Solver constants
#
############################

MATRIX_FILE = "matrix.tsv"
ENTITIES_FILE = "entities.txt"
ENTITIES = "entities"
FORMAT = ".json"
PERCENTILE = "80"
SOLVER_NAME = "semEP-nodeV1"
TMP_FILE = ".tmp_nt_results_query"
EMPTY_JSON = "[]"
INPUT_DICTIONARY = "input/input_dictionary.txt"
DICTIONARY_QUERY = "input/query_dictionary.txt"
DICTIONARY_CONSTRUCT = "input/construct_dictionary.txt"
boolean_value = ["biomarkers", "bioEgfr", "bioAlk", "bioRos1", "chemotherapy", "tki", "immunotherapy",
                 "antiangiogenic", "radiationtherapy", "surgery", "familialAntecedents", "systemicProgression",
                 "localProgression", "brainMetastasis"]
TO_PLOT = "output/file_to_plot.csv"

############################
#
# Query constants
#
############################


BASIC_QUERY_PATIENTS1 = """PREFIX iASIS_vocab: <http://project-iasis.eu/vocab/>
PREFIX iASIS: <http://project-iasis.eu/>
PREFIX rdf:	<http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

 Select distinct ?p1 where {
 
"""
BASIC_QUERY_CONSTRUCT = """PREFIX iASIS_vocab: <http://project-iasis.eu/vocab/>
PREFIX iASIS: <http://project-iasis.eu/>
PREFIX rdf:	<http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

 construct { 
 
"""

BASIC_QUERY_SIMILARITY = """PREFIX iASIS_vocab: <http://project-iasis.eu/vocab/>
PREFIX iASIS: <http://project-iasis.eu/>
PREFIX rdf:	<http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

 Select distinct ?p1 ?p2 ?s where {

"""
TO_SIMILARITY = """?r <http://project-iasis.eu/vocab/isBetweenPatients> ?p1.
         ?r <http://project-iasis.eu/vocab/isBetweenPatients> ?p2.
         ?r <http://project-iasis.eu/vocab/similarityValue> ?s.

"""
CLOSE = "}"

FILTER = 'FILTER'

AND = "&&"
SPC = " "
OPENP = "("
CLOSEP = ")"
IN = " IN "
DOT = '.'

FROM_AGE_END = '"^^xsd:integer)'

NEWLINE = "\n"

F_INI = "?"

F_AGE_INI = '(xsd:integer(?'


def processing_similarity_data(graph):
    index_sim = {}
    patients = set()
    sim_value = []
    with open(ENTITIES_FILE, "w") as fd:
        for (e1, e2, sim) in graph:
            key = e1 + "**" + e2
            index_sim[key] = sim
            key = e2 + "**" + e1
            index_sim[key] = sim
            if e1 not in patients:
                key = e1 + "**" + e1
                index_sim[key] = "1.0"
            patients.add(e1)
            if e2 not in patients:
                key = e2 + "**" + e2
                index_sim[key] = "1.0"
            patients.add(e2)
        elements = list(patients)

        with open(MATRIX_FILE, "w") as fm:
            s = ""
            j = 0
            size_population = len(elements)
            fm.write(str(size_population) + "\n")
            fd.write(str(size_population) + "\n")
            for e1 in elements:
                j += 1
                fd.write(str(e1) + "\n")
                for e2 in range(size_population):
                    key = e1 + "**" + elements[e2]
                    sim = index_sim[key]
                    s += str(sim) + " "
                    if e2 >= j:
                        sim_value.append(float(sim))
                s = s[:-1] + "\n"
            fm.write(s)

    return sim_value, size_population


###################################
#
# Computation of the distribution
#
###################################

def compute_distributions(clusters, q_results):
    cluster_dist = {}
    population_dist = {}
    sub = "http://project-iasis.eu/LCPatient/"
    predicate_list = set()
    cluster_list = set()
    for c in clusters:
        cluster_list.add(c)
        cluster_dist[c] = {}
        for id_patient in clusters[c]:
            id_p = id_patient.lstrip(sub)
            (cluster_dist[c], predicate_list, population_dist) = get_construct_data(q_results, id_p, cluster_dist[c],
                                                                                    predicate_list, population_dist)
    return cluster_dist, predicate_list, cluster_list, population_dist


############################
#
# Filter generation
#
############################


def get_query_parameter(filter1, sentence_where, sentence_construct, set_selection, set_parameter, construct_dicc):
    set_intersection = set_selection & set_parameter
    # print("set_intersection: "+str(list(set_intersection)))
    set_diference = set_parameter - set_intersection
    # print("set_diference: " + str(list(set_diference)))
    parameter = list(set_diference)
    remaining_sentences, set_diference = create_sentence_construct(parameter, construct_dicc)

    query = BASIC_QUERY_CONSTRUCT
    query += sentence_construct + "}" + "\n"
    query += "where {" + "\n" + sentence_where
    query += remaining_sentences
    query += filter1
    query += CLOSE
    return query


def get_query(filter, sentence_where):
    query = BASIC_QUERY_PATIENTS1
    query += sentence_where
    query += filter
    # query += selection
    query += CLOSE
    return query


def get_query_similarity(filter1, sentence_where, to_similarity_where, to_similarity_filter, ):
    query_similarity = BASIC_QUERY_SIMILARITY
    query_similarity += sentence_where
    query_similarity += to_similarity_where
    query_similarity += TO_SIMILARITY
    query_similarity += filter1
    query_similarity += to_similarity_filter
    query_similarity += "FILTER (?p1 !=?p2)"
    query_similarity += CLOSE
    return query_similarity


def get_number_filter(age_range, k):
    fstr = SPC + SPC + SPC + SPC + FILTER + OPENP + SPC
    if len(age_range) != 2:
        print("Error in the age filter")
        print(age_range)
        sys.exit(1)
    if age_range[0] > age_range[1]:
        print("Error in the age range")
        print(age_range)
        sys.exit(1)
    to_similarity_filter = fstr

    fstr += F_AGE_INI + k + ') > "' + str(age_range[0]) + FROM_AGE_END
    fstr += SPC + AND + SPC
    fstr += F_AGE_INI + k + ') < "' + str(age_range[1]) + FROM_AGE_END + CLOSEP + DOT + NEWLINE
    to_similarity_filter += F_AGE_INI + k + '2) > "' + str(age_range[0]) + FROM_AGE_END
    to_similarity_filter += SPC + AND + SPC
    to_similarity_filter += F_AGE_INI + k + '2) < "' + str(age_range[1]) + FROM_AGE_END + CLOSEP + DOT + NEWLINE
    # print("Filter "+str(k))
    # print(fstr)

    return fstr, to_similarity_filter


def get_string_filters(k, value_list, input_dictionary):
    fstr = SPC + SPC + SPC + SPC + FILTER + OPENP + SPC
    to_similarity_filter = fstr
    fstr += F_INI + k + IN + OPENP
    to_similarity_filter += F_INI + k + "2" + IN + OPENP
    if k in boolean_value:
        for i in value_list:
            fstr += '"' + str(i) + '", '
            to_similarity_filter += '"' + str(i) + '", '
        fstr = fstr.strip()
        fstr = fstr.rstrip(",")
        fstr += CLOSEP + SPC + CLOSEP + DOT + NEWLINE
        to_similarity_filter = to_similarity_filter.strip()
        to_similarity_filter = to_similarity_filter.rstrip(",")
        to_similarity_filter += CLOSEP + SPC + CLOSEP + DOT + NEWLINE
    else:
        for i in value_list:
            fstr += str(input_dictionary[str(k)]) + "/" + str(i) + ">, "
            to_similarity_filter += str(input_dictionary[str(k)]) + "/" + str(i) + ">, "
        fstr = fstr.strip()
        fstr = fstr.rstrip(",")
        fstr += CLOSEP + SPC + CLOSEP + DOT + NEWLINE
        to_similarity_filter = to_similarity_filter.strip()
        to_similarity_filter = to_similarity_filter.rstrip(",")
        to_similarity_filter += CLOSEP + SPC + CLOSEP + DOT + NEWLINE
    # print("Filter "+str(k))
    # print(fstr)
    return fstr, to_similarity_filter


def create_sentence_construct(parameter, construct_dicc):
    pstr = ""
    set_parameter = set()
    for par1 in parameter:
        pstr += str(construct_dicc[str(par1)]) + "\n"
        set_parameter.add(par1)
    return pstr, set_parameter


def load_input(request, input_dictionary, query_dictionary, construct_dicc):
    filter = ""
    sentence = set()
    sentence_where = ""
    to_similarity_where = ""
    to_similarity_filter = ""
    sentence_construct = ""
    set_selection = set()
    for i in request:
        if i == "top_clusters":
            top_cluster = int(request[i][0])
            continue
        if i == "parameter":
            (sentence_construct, set_parameter) = create_sentence_construct(request[i], construct_dicc)
            break
        for k in request[i]:
            if (k == "age"):
                set_selection.add(k)
                age_range = [request["selection"]["age"]["from"], request["selection"]["age"]["to"]]
                # filter += get_number_filter(age_range, k)
                (a, b) = get_number_filter(age_range, k)
                filter += a
                to_similarity_filter += b
                s = str(query_dictionary[str(k)])
                (sentence_where, sentence, to_similarity_where) = aux1(s, sentence_where, sentence, to_similarity_where)
            else:
                set_selection.add(k)
                list = request[i][k]
                # filter += get_string_filters(k, list, input_dictionary)
                (a, b) = get_string_filters(k, list, input_dictionary)
                filter += a
                to_similarity_filter += b
                s = str(query_dictionary[str(k)])
                (sentence_where, sentence, to_similarity_where) = aux1(s, sentence_where, sentence, to_similarity_where)
    # print(sentence_where)
    # print(to_similarity_where)
    # print(to_similarity_filter)
    # print(sentence_construct)
    return filter, sentence_where, to_similarity_where, to_similarity_filter, sentence_construct, set_selection, set_parameter, top_cluster


def aux1(s, sentence_where, sentence, to_similarity_where):
    value = s.split("***")
    for l in value:
        if l not in sentence:
            sentence_where += l + "\n"
            to_similarity_where += aux2(l) + "\n"
        sentence.add(l)
    return sentence_where, sentence, to_similarity_where


def aux2(l):
    temp = l.split(" ")
    temp[0] = temp[0].rstrip("1")
    temp[0] = temp[0] + "2"
    temp[2] = temp[2].rstrip(".")
    temp[2] = temp[2].rstrip("1")
    temp[2] += "2."
    to_similarity_where = temp[0] + " " + temp[1] + " " + temp[2]
    return to_similarity_where


############################
#
# SemEP data
#
############################


def compute_percentile(sim, percentile):

    return numpy.percentile(sim, percentile)


############################
#
# Query proccesing
#
############################


def sendSPARQL(sparql_ins, origin_query, num_paging, len_max):
    cur_offset = 0
    num_paged_result = sys.maxsize
    list_result = []
    while (num_paged_result >= num_paging and cur_offset < len_max):
        cur_query = origin_query + "LIMIT" + " " + str(num_paging) + " " + "OFFSET" + " " + str(cur_offset)
        sparql_ins.setQuery(cur_query)
        sparql_ins.setReturnFormat(JSON)
        results = sparql_ins.query().convert()
        cur_list_result = results["results"]["bindings"]
        list_result = list_result + cur_list_result
        num_paged_result = len(cur_list_result)
        cur_offset = cur_offset + num_paging
    return list_result


def get_triples_data(qresults):
    entities = set()
    for cur_result in qresults:
        subj = cur_result["p1"]["value"]
        pred = cur_result["p2"]["value"]
        obj = cur_result["s"]["value"]
        # print((subj, pred, obj))
        entities.add((subj, pred, obj))
    return entities


############################
#
# Building SemEP-Node data
#
############################

def get_dicc_clusters(onlyfiles):
    cont = 1
    dicc_clusters = {}
    for filename in onlyfiles:
        key = "Cluster-" + str(cont)
        dicc_clusters[key] = []
        with open(filename) as fd:
            for line in fd:
                e = line.rstrip()
                dicc_clusters[key].append(e)
        cont = cont + 1
    return dicc_clusters


def call_semEP(threshold):
    EXE_SEM_EP = "sudo docker run -it --rm -v "
    DIR_SEM_EP = ":/data kemele/semepnode:23-05-2018 semEP-node"
    current_path = os.path.dirname(os.path.realpath(__file__))
    th = "{:.4f}".format(float(threshold))

    commd = EXE_SEM_EP + "" + current_path + "" + DIR_SEM_EP + " " + ENTITIES_FILE + " " + MATRIX_FILE + " " + str(th)
    # print("Command to execute " + commd)
    os.system(commd)
    pattern = ENTITIES_FILE.replace(".txt", "")

    results_folder = glob.glob(current_path + "/" + pattern + "-*")

    print("sudo chown -R rivas:rivas " + results_folder[0] + "*")

    os.system("sudo chown -R rivas:rivas " + results_folder[0] + "*")

    onlyfiles = [os.path.join(results_folder[0], f) for f in listdir(results_folder[0]) if
                 isfile(join(results_folder[0], f))]
    dicc_clusters = get_dicc_clusters(onlyfiles)

    for r, d, f in os.walk(results_folder[0]):
        for files in f:
            os.remove(os.path.join(r, files))
        os.removedirs(r)
    return dicc_clusters


def generate_results(top_clusters, dict_top_cluster, population):
    results = {}
    results["TopDifferentClusters"] = top_clusters
    results["DistributionsInTopClusters"] = dict_top_cluster
    results["Population"] = population
    result_json = json.dumps(results, indent=4)
    # print("Resutls JSON:")
    # print(result_json)
    return result_json


def get_patient_graph(filter, sentence_where, end_point):
    query = get_query(filter, sentence_where)
    print("QUERY patient:", query)
    sparql_ins = SPARQLWrapper(end_point)
    print("Waiting for the SPARQL Endpoint")
    qresults = sendSPARQL(sparql_ins=sparql_ins, origin_query=query, num_paging=50000000, len_max=sys.maxsize)
    # print(qresults)

    if len(qresults) == 0:
        # return EMPTY_JSON
        return None


def get_parameter_projection(filter1, sentence_where, sentence_construct, set_selection, set_parameter, construct_dicc,
                             end_point):
    query = get_query_parameter(filter1, sentence_where, sentence_construct, set_selection, set_parameter,
                                construct_dicc)
    # print("QUERY_CONSTRUCT:", query)
    sparql_ins = SPARQLWrapper(end_point)
    print("Waiting for the SPARQL Endpoint")
    c_results = sendSPARQL(sparql_ins=sparql_ins, origin_query=query, num_paging=50000000, len_max=sys.maxsize)
    # print(qresults)

    if len(c_results) == 0:
        return EMPTY_JSON
    # graph = get_construct_data(c_results)
    return c_results


def get_matrix(filter1, sentence_where, to_similarity_where, to_similarity_filter, end_point):
    query_similarity = get_query_similarity(filter1, sentence_where, to_similarity_where, to_similarity_filter)
    # print("Query_matrix:", query_similarity)
    sparql_ins = SPARQLWrapper(end_point)
    print("Waiting for the SPARQL Endpoint")
    qresults = sendSPARQL(sparql_ins=sparql_ins, origin_query=query_similarity, num_paging=100000,
                          len_max=sys.maxsize)
    # print(qresults)

    if len(qresults) == 0:
        # return EMPTY_JSON
        return None
    # print("print SIMILARITY-PATIENT ")
    # print(qresults)

    graph = get_triples_data(qresults)
    return graph


def get_construct_data(c_results, id, dict_predicate, predicate_list, population_dist):
    for cur_result in c_results:
        subj = cur_result["s"]["value"]
        if subj.__contains__(id):
            pred = cur_result["p"]["value"]
            if pred not in dict_predicate:
                predicate_list.add(pred)
                dict_predicate[pred] = {}
                obj = cur_result["o"]["value"]
                dict_predicate[pred][obj] = 0
                dict_predicate[pred][obj] += 1
            else:
                obj = cur_result["o"]["value"]
                if obj not in dict_predicate[pred]:
                    dict_predicate[pred][obj] = 0
                    dict_predicate[pred][obj] += 1
                else:
                    dict_predicate[pred][obj] += 1
            if pred not in population_dist:
                population_dist[pred] = {}
                population_dist[pred][obj] = 0
                population_dist[pred][obj] += 1
            else:
                if obj not in population_dist[pred]:
                    population_dist[pred][obj] = 0
                    population_dist[pred][obj] += 1
                else:
                    population_dist[pred][obj] += 1
    return dict_predicate, predicate_list, population_dist


def get_all_parameter(population):
    list_parameter = []
    for parameter in population:
        key1 = parameter
        for parameter_2 in population[parameter]:
            key2 = parameter_2
            key = key1 + "***" + key2
            list_parameter.append(key)
    return list_parameter


def case_true_false(key1, key2):
    column = ""
    if key2 in ["True", "False"]:
        array = key1.split("/")
        column = array[len(array) - 1]
    column += key2
    return column


def modify_column(column):
    array = column.split("/")
    if len(array) > 2:
        a = array[len(array) - 2]
        b = array[len(array) - 1]
        column = a + b
    return column


def get_file_plot(population, dict_top_cluster):
    file_csv = {}
    name = "All Patients"
    list_parameter = get_all_parameter(population)
    for sample in dict_top_cluster:
        file_csv[sample] = {}
    file_csv[name] = {}
    for l_p in list_parameter:
        string = l_p.split("***")
        key1 = string[0]
        key2 = string[1]
        for cluster in dict_top_cluster:
            if key2 in dict_top_cluster[cluster][key1]:
                column = case_true_false(key1, key2)
                column = modify_column(column)
                file_csv[cluster][column] = dict_top_cluster[cluster][key1][key2]
                file_csv[name][column] = population[key1][key2]
            else:
                column = case_true_false(key1, key2)
                column = modify_column(column)
                file_csv[cluster][column] = 0
                file_csv[name][column] = population[key1][key2]
    return file_csv


def get_csv_to_plot(file_csv):
    with open(TO_PLOT, "w") as fp:
        # Write CSV Header
        aux = "Comparison,"
        for header in file_csv:
            for p in file_csv[header]:
                aux += str(p) + ","
            aux = aux.rstrip(",")
            fp.write(aux + "\n")
            break

        for line in file_csv:
            string = line + ","
            for i in file_csv[line]:
                string += str(file_csv[line][i]) + ","
            string = string.rstrip(",")
            fp.write(string + "\n")


def plot_heatmap(file):
    codigo_r = """
    heatmap_plot <- function(dataset) {
      library(pheatmap)
      setEPS()
      postscript("output/Plot.eps")
      profiling <- read.csv(dataset)
      m <- as.matrix(profiling[, -1])
      rownames(m) <- profiling$Comparison
      n.ind <- dim(profiling)[1]
      n_cluster<-dim(profiling$Comparison)[1]
      pheatmap(m, scale = "none", clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean", clustering_method = "average",
               cutree_rows = n.ind, fontsize = 11)
      dev.off()
    }
    """
    #pdf(file="output/Plot.pdf")
    ro.r(codigo_r)
    code_py = ro.globalenv['heatmap_plot']
    read = code_py(file)

    imagen = Image.open("output/Plot.eps")
    imagen.show()
    #pdf_file = open("output/Plot.pdf")
    #read_pdf = PyPDF2.PdfFileReader(pdf_file)

    return read


##################################
#
# Running the SemEP-Node Wrapper
#
##################################

def run_wrapper(end_point, request):
    (query_dictionary, input_dictionary, construct_dicc) = load_files()
    (filter1, sentence_where, to_similarity_where, to_similarity_filter, sentence_construct, set_selection,
     set_parameter, top_cluster) = load_input(request, input_dictionary, query_dictionary, construct_dicc)
    start_time = time()
    # graph = get_patient_graph(filter1, sentence_where, end_point)
    graph2 = get_matrix(filter1, sentence_where, to_similarity_where, to_similarity_filter, end_point)
    elapsed_time = time() - start_time
    print("Elapsed time (get_matrix(): " + str(elapsed_time))
    sim, size_population = processing_similarity_data(graph2)
    threshold = compute_percentile(sim, PERCENTILE)

    clusters = call_semEP(threshold)
    c_results = get_parameter_projection(filter1, sentence_where, sentence_construct, set_selection, set_parameter,
                                         construct_dicc, end_point)
    # print("num_Clusters")
    # print(len(clusters))
    distributions, predicate_list, cluster_list, population_dist = compute_distributions(clusters, c_results)
    top_clusters, dict_top_cluster, population = compare_distributions(distributions, clusters, top_cluster,
                                                                       set_parameter, population_dist, size_population)

    file_csv = get_file_plot(population, dict_top_cluster)

    get_csv_to_plot(file_csv)
    aa = plot_heatmap(TO_PLOT)

    resutls_json = generate_results(top_clusters, dict_top_cluster, population)
    return resutls_json


def load_files():
    q_dicc = open(DICTIONARY_QUERY)
    length_q_dicc = len(q_dicc.readlines())
    query_dicc = {}
    q_dicc.seek(0)

    i_dicc = open(INPUT_DICTIONARY)
    length_i_dicc = len(i_dicc.readlines())
    input_dicc = {}
    i_dicc.seek(0)

    c_dicc = open(DICTIONARY_CONSTRUCT)
    length_c_dicc = len(c_dicc.readlines())
    construct_dicc = {}
    c_dicc.seek(0)

    length = max(length_q_dicc, length_i_dicc, length_c_dicc)
    # print(length)
    for i in range(length):
        tok = q_dicc.readline().rstrip().split(",")
        if len(tok) > 1:
            query_dicc[tok[0]] = tok[1]
        tok = i_dicc.readline().rstrip().split(",")
        if len(tok) > 1:
            input_dicc[tok[0]] = tok[1]
        tok = c_dicc.readline().rstrip().split(",")
        if len(tok) > 1:
            construct_dicc[tok[0]] = tok[1]
    return query_dicc, input_dicc, construct_dicc
