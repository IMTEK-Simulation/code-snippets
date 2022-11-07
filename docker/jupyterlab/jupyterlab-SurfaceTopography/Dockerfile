# Build from within repository root 
#
# Modify conda enviornment running jupyter notebooks
FROM jupyter/scipy-notebook:latest

COPY --chown=${NB_UID}:${NB_GID} jupyterlab-SurfaceTopography/conda-requirements.in "/home/${NB_USER}/"
RUN pip install pip-tools
RUN pip-compile "/home/${NB_USER}/conda-requirements.in" > "/home/${NB_USER}/conda-requirements.txt"
RUN pip install --quiet --no-cache-dir --requirement "/home/${NB_USER}/conda-requirements.txt" && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

# FROM imkteksim/dtool-jupyter:latest
# Install from requirements.txt file
USER root
RUN apt-get --yes update && \
    apt-get --yes install \
        clang cmake curl g++ gdb git m4 wget \
        libboost-test-dev \ 
        libcurl4-openssl-dev \
        libeigen3-dev \
        libfftw3-dev \
        libgmp-dev \
        libnetcdf-dev \
        libopenblas-base \
        libopenblas-dev \
        python3-pip \
        python3-dev \
        python3-netcdf4 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
USER ${NB_UID}

# Workaround for building packages and introducing an ipykernel outside of the jupyterlab conda environment
RUN PATH=$(echo "$PATH" | sed -e 's|/opt/conda[^:]*:||g') && ( \
        pip install --no-binary numpy SurfaceTopography && \
        pip install ipykernel && \
        python3 -m ipykernel install --user --name SurfaceTopography --display-name "Python (SurfaceTopography)" && \
        echo "c.MultiKernelManager.default_kernel_name = 'SurfaceTopography'" >> ${HOME}/.jupyter/jupyter_notebook_config.py) && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

# Workaround for modifying system python outside of the jupyterlab conda environment
COPY --chown=${NB_UID}:${NB_GID} jupyterlab-SurfaceTopography/requirements.txt "/home/${NB_USER}/"
RUN PATH="/home/${NB_USER}/.local/bin:$(echo "$PATH" | sed -e 's|/opt/conda[^:]*:||g')" && ( \
        pip install --quiet --no-cache-dir --requirement "/home/${NB_USER}/requirements.txt" )

WORKDIR /home/jovyan/work

# Custom fonts config
COPY jupyterlab-SurfaceTopography/start-wrapper.sh /usr/local/bin/
RUN mkdir -p /home/${NB_USER}/.config/fontconfig
COPY jupyterlab-SurfaceTopography/etc/fonts.conf /home/${NB_USER}/.config/fontconfig/fonts.conf
RUN mkdir -p /home/${NB_USER}/.fonts
RUN rm -rf /home/${NB_USER}/.cache/matplotlib
CMD ["start-wrapper.sh", "start-notebook.sh"]
