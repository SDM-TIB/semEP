#  SemEP-node service

## 1.  About

Service for detection of communities of Lung Cancer Patients in the iASiS KG

## 2. Usage

Setup the variable with the iASiS endpoint.

`export IASISKG_ENDPOINT=<Endpoint address>`

For example:

`export IASISKG_ENDPOINT=http://10.114.113.14:8181/sparql`

Execute

`./run_semEP-Node_service.sh`

### Example of calls to service

curl -H "Content-Type: application/json" -X POST -d '{"selection": {"gender": ["M", "F"], "age": {"from": 20, "to": 90}, "tumorStage": ["I","III", "IV", "II"]}, "parameter": ["gender" ,"age", "tumorStage", "toxicHabits", "mutationType"]}' http://0.0.0.0:5001/semepnode --output example.json

## Documentation:

The documentation of the API is available on:

https://docs.google.com/document/d/1elKzAjmcOPcRw1fmV3GULIkqVrVGAgxn1mkI4pEvwhc/edit#heading=h.gjdgxs

