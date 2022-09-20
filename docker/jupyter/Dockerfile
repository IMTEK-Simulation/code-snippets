# Start from a core stack version
FROM jupyter/scipy-notebook:latest
# Install from requirements.txt file
COPY --chown=${NB_UID}:${NB_GID} requirements.txt /tmp/
RUN pip install --quiet --no-cache-dir --requirement /tmp/requirements.txt && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

COPY start-wrapper.sh /usr/local/bin/

# Custom fonts config
RUN mkdir -p /home/${NB_USER}/.config/fontconfig
COPY etc/fonts.conf /home/${NB_USER}/.config/fontconfig/fonts.conf
COPY etc/fonts /home/${NB_USER}/.fonts
RUN rm -rf /home/${NB_USER}/.cache/matplotlib
RUN fc-cache --force
RUN python3 -c 'import matplotlib.font_manager; matplotlib.font_manager._load_fontmanager(try_read_cache=False)'
CMD ["start-wrapper.sh", "start-notebook.sh"]
