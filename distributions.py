#!/usr/bin/env python3

import sys
import copy
from scipy import stats

MIN_NUM_ELEM = 3


def get_list_value(sample):
    list_value = []
    for cluster in sample:
        temp = []
        for key1 in cluster:
            for key2 in cluster[key1]:
                temp.append(cluster[key1][key2])
        list_value.append(temp)
    return list_value


def get_global_list_value(population):
    list_value = []
    for key1 in population:
        for key2 in population[key1]:
            list_value.append(population[key1][key2])
    return list_value


def compute_ks_2(population, sample):
    ks_cluster_pvalues = []
    list_value = get_list_value(sample)
    global_list = get_global_list_value(population)
    i = 1
    for c in list_value:
        values = c
        pvalue = stats.ks_2samp(global_list, values)[1]
        ks_cluster_pvalues.append(("Cluster-" + str(i), pvalue))
        i += 1
    ks_cluster_pvalues.sort(key=lambda tup: tup[1])
    return ks_cluster_pvalues


def clusters_top(ks_clusters_pvalues, clusters, max_clusters, sample):
    top_clusters = []
    dict_top_cluster = {}
    for (c, pvalue) in ks_clusters_pvalues:
        if len(clusters[c]) >= MIN_NUM_ELEM:
            if len(top_clusters) < max_clusters:
                top_clusters.append(c)
                pos = int(c.split("-")[1]) - 1
                dict_top_cluster[c] = sample[pos]
            else:
                break
    # print(dict_top_cluster)
    return top_clusters, dict_top_cluster


def get_age_range(population_dist, key):
    dict_age = {}
    for age in population_dist[key]:
        if int(age) > 75:
            if "Age >75" not in dict_age:
                dict_age["Age >75"] = {}
                dict_age["Age >75"] = 0
                dict_age["Age >75"] += population_dist[key][age]
            else:
                dict_age["Age >75"] += population_dist[key][age]
        elif int(age) > 65 & int(age) < 76:
            if "Age 66-75" not in dict_age:
                dict_age["Age 66-75"] = {}
                dict_age["Age 66-75"] = 0
                dict_age["Age 66-75"] += population_dist[key][age]
            else:
                dict_age["Age 66-75"] += population_dist[key][age]
        elif int(age) > 55 & int(age) < 66:
            if "Age 56-65" not in dict_age:
                dict_age["Age 56-65"] = {}
                dict_age["Age 56-65"] = 0
                dict_age["Age 56-65"] += population_dist[key][age]
            else:
                dict_age["Age 56-65"] += population_dist[key][age]
        elif int(age) > 45 & int(age) < 56:
            if "Age 46-55" not in dict_age:
                dict_age["Age 46-55"] = {}
                dict_age["Age 46-55"] = 0
                dict_age["Age 46-55"] += population_dist[key][age]
            else:
                dict_age["Age 46-55"] += population_dist[key][age]
        elif int(age) > 35 & int(age) < 46:
            if "Age 36-45" not in dict_age:
                dict_age["Age 36-45"] = {}
                dict_age["Age 36-45"] = 0
                dict_age["Age 36-45"] += population_dist[key][age]
            else:
                dict_age["Age 36-45"] += population_dist[key][age]
        else:
            if "Age <35" not in dict_age:
                dict_age["Age <35"] = {}
                dict_age["Age <35"] = 0
                dict_age["Age <35"] += population_dist[key][age]
            else:
                dict_age["Age <35"] += population_dist[key][age]
    return dict_age


def graphic_with_age(set_parameter, uri, population_dist):
    population = {}
    if "age" in set_parameter:
        age_range = get_age_range(population_dist, uri)
        for parameter in population_dist:
            if parameter == uri:
                population[parameter] = {}
                for age in age_range:
                    population[parameter][age] = age_range[age]
                continue
            population[parameter] = population_dist[parameter]
        return population
    return population_dist


def graph_cluster(distributions, uri):
    top_cluster = []
    for cluster in distributions:
        sample = {}
        for parameter in distributions[cluster]:
            if parameter == uri:
                c_age = get_age_range(distributions[cluster], uri)
                sample[parameter] = {}
                for age in c_age:
                    sample[parameter][age] = c_age[age]
                continue
            sample[parameter] = distributions[cluster][parameter]
        top_cluster.append(sample)
    return top_cluster


def normarlization_all_parametres(population_dist, size):
    copy_population = copy.deepcopy(population_dist)
    for predicate in copy_population:
        for obj in copy_population[predicate]:
            copy_population[predicate][obj] = copy_population[predicate][obj] / size
    return copy_population


def get_size(population):
    total_patients = 0
    for predicate in population:
        for obj in population[predicate]:
            total_patients += population[predicate][obj]
        return total_patients


def normarlization_clusters(diff_cluster):
    clusters = []
    for cl in diff_cluster:
        size = get_size(cl)
        clusters.append(normarlization_all_parametres(cl, size))
    return clusters


def compare_distributions(distributions, clusters, top_cluster, set_parameter, population_dist, size_population):

    population_dist = graphic_with_age(set_parameter, "http://project-iasis.eu/vocab/hasDiagnosisAge", population_dist)
    diff_cluster = graph_cluster(distributions, "http://project-iasis.eu/vocab/hasDiagnosisAge")

    population = normarlization_all_parametres(population_dist, size_population)
    sample = normarlization_clusters(diff_cluster)

    ks_2 = compute_ks_2(population, sample)

    top_clusters, dict_top_cluster = clusters_top(ks_2, clusters, top_cluster, sample)


    return top_clusters, dict_top_cluster, population, diff_cluster, population_dist


def main(*args):
    print("Done")
    # print(stats.ks_2samp(x, y))


if __name__ == "__main__":
    main(*sys.argv[1:])
