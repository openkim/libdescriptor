# base image
FROM ghcr.io/openkim/developer-platform

# set root user
USER root

# Install clang-12 (taken from enzyme dev-container Dockerfile)
RUN apt-get -q update \
    && apt-get install -y --no-install-recommends ca-certificates software-properties-common curl gnupg2 git\
    && curl -fsSL https://apt.llvm.org/llvm-snapshot.gpg.key|apt-key add - \
    && apt-add-repository "deb http://apt.llvm.org/`lsb_release -cs`/ llvm-toolchain-`lsb_release -cs`-12 main" || true \
    && apt-get -q update \
    && apt-get install -y --no-install-recommends zlib1g-dev lldb ninja-build llvm-12-dev clang-format clang-12 libclang-12-dev libomp-12-dev lld-12\
    && python -m pip install --upgrade pip setuptools \
    && python -m pip install lit==12.0.1 pathlib2 \
    && touch /usr/lib/llvm-12/bin/yaml-bench

# Install Enzyme
WORKDIR /opt/enzyme
RUN git clone https://github.com/enzymead/enzyme \
    && cd enzyme/enzyme \
    && mkdir build \
    && cd build \
    && CC=clang-12 CXX=clang++-12 cmake .. -DLLVM_DIR=/usr/lib/llvm-12/cmake -DLLVM_EXTERNAL_LIT=/usr/local/bin/lit \
    && make \
    && make install

# Switch back to openkim env
WORKDIR /home/openkim
USER openkim

# Download and compile libdescriptor
# providing sudo password from commandline, as anyway it is out in open!
# using build-type release, as enzyme fails with debug
RUN git clone https://github.com/ipcamit/colabfit-descriptor-library \
    && cd colabfit-descriptor-library \
    && mkdir build && cd build \
    && CC=clang-12 CXX=clang++-12 cmake .. -DENZYME_LIB=/usr/local/lib -DCMAKE_BUILD_TYPE=Release \
    && make \
    && echo openkim | sudo -S make install

# Install 

