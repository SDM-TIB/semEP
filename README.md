#  SemEP-node service

## 1.  About

Service for detection of communities of Lung Cancer Patients in the iASiS KG

## 2. Usage

Setup the variable with the iASiS endpoint.

`export IASISKG_ENDPOINT=<Endpoint address>`

For example:

`export IASISKG_ENDPOINT=http://node2.research.tib.eu:19191/sparql`

Execute

`./run_semEP-Node_service.sh`

### Example of calls to service

curl -H "Content-Type: application/json" -X POST -d '{"top_clusters": [4], "selection": {"gender": ["C0025266", "C0043210"], "age": {"from":20, "to": 90}, "tumorStage": ["C0027651","C0152013", "C0024121", "C4304521", "C0242379", "C0001418"], "surgery": ["True", "False"]}, "parameter": ["tki", "gender" ,"age", "tumorStage", "familiarAntecedents", "brainMetastasis"]}' http://10.115.83.140:5000/semepnode --output output/output.json

# Note
If you use Postman tool::

Request: 
http://10.115.83.140:5000/semepnode

HTTP (POST)

Body:
{"top_clusters": [4], "selection": {"gender": ["C0025266", "C0043210"], "age": {"from":20, "to": 90}, "tumorStage": ["C0027651","C0152013", "C0024121", "C4304521", "C0242379", "C0001418"], "surgery": ["True", "False"]}, "parameter": ["tki", "gender" ,"age", "tumorStage", "familiarAntecedents", "brainMetastasis"]}

Apply by: Send and Download

## Documentation:

The documentation of the API is available on:

https://docs.google.com/document/d/1elKzAjmcOPcRw1fmV3GULIkqVrVGAgxn1mkI4pEvwhc/edit#heading=h.gjdgxs

