FROM ubuntu:17.10

WORKDIR /SemepNodeService
ADD . /SemepNodeService

RUN apt-get --assume-yes update
RUN apt-get --assume-yes upgrade
RUN apt-get --assume-yes install python3 python3-numpy python3-flask python3-sparqlwrapper gunicorn3 python3-requests python3-pip python3-scipy
RUN pip3 install py_stringmatching

# Make port 5001 available to the world outside this container
EXPOSE 5001

# Define environment variable
ENV NAME SemEP

# Run app.py when the container launches
CMD ./run_semEP-Node_service.sh

