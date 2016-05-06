# run as:
#   docker build -f Dockerfile -t csi .

# base everything on a recent Ubuntu
FROM debian:latest

# get system packages up to date then install a basic scientific python
RUN apt-get update && apt-get -y upgrade && \
    apt-get -y install python3 python3-dev \
         python3-numpy python3-scipy python3-matplotlib \
	 python3-seaborn python3-setuptools wget libhdf5-dev

#install new pandas and h5py by hand. also needs cython

#pandas
RUN wget https://pypi.python.org/packages/source/p/pandas/pandas-0.17.1.tar.gz#md5=1e18b9a5496ec92752b3cb6674bbe987
RUN tar -zxvf pandas-0.17.1.tar.gz && rm pandas-0.17.1.tar.gz
RUN cd pandas-0.17.1 && python3 setup.py install && cd ..

#cython (for h5py)
RUN wget http://cython.org/release/Cython-0.23.4.tar.gz
RUN tar -zxvf Cython-0.23.4.tar.gz && rm Cython-0.23.4.tar.gz
RUN cd Cython-0.23.4 && python3 setup.py install && cd ..

#h5py (note the HDF5_DIR - libhdf5-dev installs somewhere where the setup doesn't see it by default)
RUN wget https://pypi.python.org/packages/source/h/h5py/h5py-2.5.0.tar.gz#md5=6e4301b5ad5da0d51b0a1e5ac19e3b74
RUN tar -zxvf h5py-2.5.0.tar.gz && rm h5py-2.5.0.tar.gz
RUN cd h5py-2.5.0 && HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial/ python3 setup.py install && cd ..

# add and configure our code
RUN mkdir /scripts
COPY csi/ /scripts/csi/
COPY html/ /scripts/html/
MAINTAINER Sam Mason <sam@samason.uk>
ENTRYPOINT ["bash", "/scripts/csi/csi_wrapper.sh"]
CMD ["--help"]
