FROM quay.io/fenicsproject/stable:latest

LABEL maintainer="johannes.hoermann@imtek.uni-freiburg.de"

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    LANGUAGE=C.UTF-8 \
    PYTHONUNBUFFERED=1

RUN apt-get --yes update && \
    apt-get -y install \
        python3-pip \
        gmsh && \
    python3 -m pip install meshio matscipy[cli,fenics] && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

USER fenics
WORKDIR $FENICS_HOME
ENTRYPOINT ["/sbin/my_init","--quiet","--","/sbin/setuser","fenics","/bin/bash","-l","-c"]
CMD ["/bin/bash","-i"]
