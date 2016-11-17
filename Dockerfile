# run as:
#   docker build -f Dockerfile -t csi .

# base everything on a recent Ubuntu
FROM debian:latest

# get system packages up to date then install a basic scientific python
RUN apt-get update && apt-get -y upgrade && \
    apt-get -y install python3 python3-dev \
         python3-numpy python3-scipy python3-matplotlib \
	 python3-seaborn python3-setuptools wget libhdf5-dev

#pip in the dependencies that need to be all fancy and new, or just don't exist in apt-get
#note h5py needing the HDF5_DIR thing as libhdf5-dev installs somewhere it doesn't see normally
RUN easy_install3 pip
RUN pip install pandas --upgrade
RUN pip install cython
RUN HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial/ pip install h5py

# add and configure our code
RUN mkdir /scripts
COPY csi/ /scripts/csi/
COPY html/ /scripts/html/

# grab the expression filter utility
RUN apt-get -y install git
RUN git clone https://github.com/cyversewarwick/expression_filter
RUN cp /expression_filter/scripts/expression_filter.py /scripts/csi

#set up analysis crash text file
RUN apt-get -y install git
RUN git clone https://github.com/cyversewarwick/analysis_crash.git

MAINTAINER Sam Mason <sam@samason.uk>
ENTRYPOINT ["bash", "/scripts/csi/csi_tarwrapper.sh"]
CMD ["--help"]
