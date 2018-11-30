#!/bin/bash

echo "***** Running SemEP-Node service on port 18890 ***** "
gunicorn3 -w 4 -b 10.115.83.140:5000 --timeout 9000 semEP_node_service:app
