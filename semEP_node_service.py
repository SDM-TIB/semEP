#!/usr/bin/env python3
#
# Description: POST service for SemEP-Node
#
# Author: Guillermo Palma (palma@l3s.de)
# Copying:  MIT License

import sys, os
from flask import Flask, jsonify, abort, request, send_from_directory, send_file
from semEP_node_wrapper import run_wrapper
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

KG = "http://node2.research.tib.eu:19191/sparql"
#KG = "http://node2.research.tib.eu:9191/sparql"
#KG = "http://localhost:18890/sparql"
#KG = os.environ["IASISKG_ENDPOINT"]
M = "M"
F = "F"

STAGE_1 = "C0027651"
STAGE_2 = "C0152013"
STAGE_3 = "C0024121"
STAGE_4 = "C4304521"
STAGE_5 = "C0242379"
STAGE_6 = "C0001418"

SMOKER      = "Smoker"
NONSMOKER   = "NonSmoker"

EGFR        = "EGFR"
NONMUTATION = "NonMutation"

def check_gender(lgender):
    if len(lgender) == 0:
        abort(400)
    for e in lgender:
        if e != 'C0043210' and e != 'C0025266':
            logger.error("Error in the gender request")
            abort(400)
            
def check_age_range(age):
    from_age, to_age = age[0], age[1]
    if to_age < from_age:
        logger.error("Error in the age range request")
        abort(400)
""""
def check_stage(tumor_stage):
    if len(tumor_stage) == 0:
        logger.error("Error, Tumor stage empty")
        abort(400)
    for e in tumor_stage:
        if STAGE_1 != e and  STAGE_2 != e and STAGE_3 != e and STAGE_4 != e:
            logger.error("Error in the tumor stage request")
            abort(400)
"""
def check_toxic_habits(toxic_habits):
    if len(toxic_habits) == 0:
        logger.error("Error, toxic habits empty")
        abort(400)
    for e in toxic_habits:
        if SMOKER != e and  NONSMOKER != e:
            logger.error("Error in the toxic habits request")
            abort(400)

def check_mutation_type(mutation_type):
    if len(mutation_type) == 0:
        logger.error("Error, mutation type empty")
        abort(400)
    for e in mutation_type:
        if EGFR != e and NONMUTATION != e:
            logger.error("Error in the mutation type request")
            abort(400)

def get_selection_parameters(input_request):
    gender = input_request["selection"]["gender"]
    logger.info("gender" + str(gender))
    check_gender(gender)

    age_range = [input_request["selection"]["age"]["from"], input_request["selection"]["age"]["to"]] 
    logger.info("age_range" + str(age_range))
    check_age_range(age_range)
    
    tumor_stage = input_request["selection"]["tumorStage"]
    logger.info("tumor_stage" + str(tumor_stage))
    #check_stage(tumor_stage)

    return (gender, age_range, tumor_stage)

app = Flask(__name__)

@app.route('/semepnode', methods=['POST'])
def run_semep_node_service():
    if ((not request.json) or (not "parameter" in request.json) or (not "selection" in request.json)):
        logger.error("Error in the format of the input")
        abort(400)
    logger.info("Hello, this is the SemEP-Node service")
    endpoint = KG
    logger.info("Endpoint: " + str(endpoint))
    #selection = get_selection_parameters(request.json)
    #logger.info("Selection: " + str(selection))
    parameters = request.json["parameter"]
    if len(parameters) == 0:
        logger.error("List of parameters is empty")
        abort(400)
    #logger.info("Parameters " + str(parameters))
    

    result = run_wrapper(endpoint, request.json)
    #logger.info("Sending the results in JSON format")

    response = make_response(result, 200)
    #response.mimetype = "application/javascript"
    return response

def run_service():
    app.run(debug=True)

def main(*args):
    if len(args) == 1:
        myhost = args[0]
    else:
        myhost = "0.0.0.0"
    app.run(debug=False, host=myhost)
    
if __name__ == '__main__':
     main(*sys.argv[1:])
