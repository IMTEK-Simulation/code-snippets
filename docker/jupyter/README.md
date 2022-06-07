# JupyterLab with jupytext and other custom extensions

## Build & run


    docker build -t imkteksim/dtool-jupyter .
    docker run -it -p 8888:8888 -v ${HOME}/.config/dtool:/home/jovyan/.config/dtool -v ${HOME}:/home/jovyan/work imkteksim/dtool-jupyter

## Creation protocol

Describes how to extract imports from python scripts and dump as pip-processible `requirements.txt`.

```console
python3 -m venv ~/venv/celan-venv
source ~/venv/clean-venv/bin/activate

pip install --upgrade pip
pip install wheel
pip install pip-tools

pip freeze > initial.txt
cat initial.txt | sed 's/==.*//g' > initial.in

# extract requirements from some scripts with pipreqs
cat pipreqs_requirements.txt | sed 's/==.*//g' > pipreqs_requirements.in

cp pipreqs_requirements.in edited_requirements.in  # and edit manually

cat initial.in edited_requirements.in

cat initial.in edited_requirements.in > requirements.in

# remove pkg_resources from requirements.in,
# manually add packages, i.e. jupytext, openpyxl, ...

pip-compile requirements.in
```
