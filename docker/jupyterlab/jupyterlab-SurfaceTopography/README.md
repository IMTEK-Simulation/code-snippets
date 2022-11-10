# Sample setup jupyterlab-SurfaceTopography

## Changes to the container

This docker image bases on the [Jupyter Docker Stacks](https://jupyter-docker-stacks.readthedocs.io/en/latest/).
These base images use condato run a JupyterLab instance. Extensions to the Jupyter ecosystem (i.e. `jupytext`) 
must hence go into the `conda-requirements.in` file. 

`SurfaceTopography` does not work well within a conda environment. Instead, this image introduces a separate 
ipython kernel `SurfaceTopography` that just uses the container-internal system Python and a few system libraries,
i.e. `numpy`. Anything needed for production in connection with `SurfaceTopography` goes in `requirements.in` or
into `requirements.txt` directly. When modifying `requirements.in`, regenerate `requirments.txt`with pinned versions
as described below.

If in need to quickly modify the `SurfaceTopography` kernel environment at runtime, pay attention to
breaking out of the default conda enviornment in the shell

* `conda deactivate`, and
* removing any conda-related artifacts from the `$PATH`, i.e. with `export PATH=$(echo "$PATH" | sed -e 's|/opt/conda[^:]*:||g')`

The latter is necessary to avoid conda compiler interfering with custom build steps when compiling against system libraries.

## Python requirements

To (re-)generate a `requirements.txt` with fixed package versions, use

```console
python3 -m venv venv
source venv/bin/activate

pip install --upgrade pip
pip install wheel
pip install pip-tools

pip-compile requirements.in > requirements.txt
```

## Build container

Build container from parent directory (`docker/jupyterlab`) with

    docker build -t imteksim/jupyterlab-surfacetopography -f jupyterlab-SurfaceTopography/Dockerfile .

or write directly to tar for transfer with

    docker build --output type=tar,dest=imteksim-jupyterlab-SurfaceTopography.tar -t imteksim/jupyterlab-surfacetopography -f jupyterlab-SurfaceTopography/Dockerfile .

and load at target destination with

    docker load --input imteksim-jupyterlab-SurfaceTopography.tar

## Certificates with acme.sh

Identify the public hostname used for issuing SSL certificates, i.e.

    export HOSTNAME="SOME-UUID.fr.bw-cloud-instance.org"

for a public bw-cloud instance and a seecure directory to hold certificates, i.e.

    export CERTDIR=${HOME}/acme.sh

and use

    mkdir -p acme.sh

    # zerossl
    docker run -v "${CERTDIR}:/acme.sh" --rm neilpang/acme.sh:latest \
      acme.sh --register-account -m ${MY_EMAIL_ADDRESS}
    docker run -v "${CERTDIR}:/acme.sh" -p 80:80 --rm neilpang/acme.sh:latest \
      acme.sh --issue -d "${HOSTNAME}" --standalone --force

    # or, letsencrypt
    docker run -v "${CERTDIR}:/acme.sh" --rm neilpang/acme.sh:latest \
      acme.sh --server letsencrypt --register-account -m ${MY_EMAIL_ADDRESS}
    docker run -v "${CERTDIR}:/acme.sh" -p 80:80 --rm neilpang/acme.sh:latest \
      acme.sh --server letsencrypt --issue -d "${HOSTNAME}" --standalone --force

    sudo chown -R $USER:$USER $HOME/acme.sh

to (re-)issue certiicates.

## Fonts

Microsoft fonts can be installed with

    sudo apt-get install -y fontconfig
    sudo apt-get install -y ttf-mscorefonts-installer
    sudo fc-cache -f

Check with

    $ fc-match Arial
    LiberationSans-Regular.ttf: "Liberation Sans" "Regular"

Find location of installed fonts with

    $ dpkg-query -L ttf-mscorefonts-installer
    /.
    /usr
    /usr/lib
    /usr/lib/msttcorefonts
    /usr/lib/msttcorefonts/update-ms-fonts
    /usr/share
    /usr/share/doc
    /usr/share/doc/ttf-mscorefonts-installer
    /usr/share/doc/ttf-mscorefonts-installer/README.Debian
    /usr/share/doc/ttf-mscorefonts-installer/changelog.gz
    /usr/share/doc/ttf-mscorefonts-installer/copyright
    /usr/share/fonts
    /usr/share/fonts/truetype
    /usr/share/lintian
    /usr/share/lintian/overrides
    /usr/share/lintian/overrides/ttf-mscorefonts-installer
    /usr/share/package-data-downloads
    /usr/share/package-data-downloads/ttf-mscorefonts-installer
    /var
    /var/lib
    /var/lib/msttcorefonts

## Run Jupyter Lab

Identify a persistent workspace directory, i.e.
    
    export WORKDIR=${HOME}/work

and, optionally, a directory of additional fonts to make available to matplotlib, i.e. as above

    export FONTDIR=/usr/share/fonts/truetype/msttcorefonts

Run public Jupyter Lab on default https port with fixed password

    docker run -d --rm -p 443:443 \
        -v ${WORKDIR}:/home/jovyan/work \
        -v ${FONTDIR}:/fonts \
        -v ${CERTDIR}/${HOSTNAME}:/etc/ssl/notebook \
        imteksim/jupyterlab-surfacetopography \
        start-notebook.sh \
        --ServerApp.keyfile=/etc/ssl/notebook/${HOSTNAME}.key \
        --ServerApp.certfile=/etc/ssl/notebook/${HOSTNAME}.cer \
        --ServerApp.password='argon2:$argon2id$v=19$m=10240,t=10,p=8$Jg1gFoiE0ybmBCVG+1s6uA$Jxz4/1591Z7so7JK2M1lRA' \
        --ServerApp.ip='*' \
        --ServerApp.port=443

Note that the hash of the password is passed.

    argon2:$argon2id$v=19$m=10240,t=10,p=8$Jg1gFoiE0ybmBCVG+1s6uA$Jxz4/1591Z7so7JK2M1lRA

corresponds to password `imtek-simulation` and has been generated with

    from notebook.auth import passwd
    passwd()

Add

    --user root -e GRANT_SUDO=yes -e NB_GID=100 -e NB_USER=jovyan

to the docker command line arguments if you need passwordless `sudo` within the container.

## Access Jupyter Lab web interface

Make sure `https` port 443 is reachable, navigate browser to `https://${HOSTNAME}` and log in with password specified on Jupyter Lab launch.
