FROM ubuntu
RUN apt-get -y update; \
  DEBIAN_FRONTEND=noninteractive apt-get --no-install-recommends install -y \
  libopenblas-dev gfortran build-essential python-software-properties python-dev \
  python-numpy python-scipy python-matplotlib python-pip python-h5py wget unzip \
  openmpi-bin libopenmpi-dev openmpi-doc zlib1g-dev python-tk; \
  wget https://bootstrap.pypa.io/ez_setup.py -O - | python; \
  pip install nugridpy; \
  echo 'download, unpack and configure openmpi'; \
  wget http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-1.8.3/src/hdf5-1.8.3.tar.gz > /dev/null; \
  tar -xzvf hdf5-1.8.3.tar.gz > /dev/null; mkdir hdf5; cd hdf5-1.8.3; \
  echo 'start configure'; \
  ./configure --prefix=/hdf5 F77=gfortran > /dev/null 2>&1; echo 'start make'; \
  make > /dev/null 2>&1; echo 'make done'; \
  make check > /dev/null 2>&1; echo 'make check'; make install> /dev/null 2>&1; \
  echo 'check hdf5'; ls /hdf5; cd ../; \
  wget https://github.com/NuGrid/SE-library/archive/master.zip > /dev/null; \
  unzip master.zip; mkdir SE; cd NuSE-master/SE; \
  ./configure --prefix=/SE F77=gfortran --with-hdf5=/hdf5 > /dev/null 2>&1; \
  echo 'start make SE'; make > /dev/null 2>&1; echo 'start make check SE'; \
  make check > /dev/null 2>&1; \
  echo 'start make install SE'; make install; echo 'check SE'; ls /SE; cd ../
ADD . /home/NuPPN/
