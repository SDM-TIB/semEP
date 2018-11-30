FROM ubuntu:18.04

WORKDIR /SemepNodeService
ADD . /SemepNodeService

RUN apt-get --assume-yes update
RUN apt-get --assume-yes upgrade
RUN apt-get install -y tzdata
RUN apt-get --assume-yes install python3 python3-pip
RUN apt-get install r-base --assume-yes
RUN pip3 install -r requirements.txt
RUN R --vanilla -e 'install.packages("pheatmap",repos="http://cran.us.r-project.org")'
# Make port 5001 available to the world outside this container
EXPOSE 5000

# Define environment variable
ENV NAME SemEP

# Run app.py when the container launches
ENTRYPOINT python3 /SemepNodeService/semEP_node_service.py
