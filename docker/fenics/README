# Build and run with podman

To build with podman, evoke

```console
$ podman build -f Dockerfile
STEP 1: ...
...
...
STEP 9: COMMIT
--> 6eab0eda254
6eab0eda254bce72ef7676ac6278dd383321d6e51efaa56a84db2c9ff13a73e2
```

and tag the resulting image appropriately, i.e. as

```console
$ podman tag 6eab0eda254bce fenics_and_gmsh:latest
```

to easily identify it later on.

For a jupyter instance at a random free port of your choice, in this sample at 8765, with

```console
$ podman run -dt --rm --publish 8765:8765 -v $(pwd):/home/fenics/shared localhost/fenics_and_gmsh:latest 'jupyter notebook --port=8765'
655544dc3351b65a49b70b8906cbe3a1dca9728711c797711b6c5b722fae62c9
```

to make the content of your current working directory (`pwd`) on the host available within the container at `/home/fenics/shared`.

Show log output using the container's id returned by the run command with

```console
$ podman logs 655544dc335
[I 20:58:34.305 NotebookApp] Writing notebook server cookie secret to /home/fenics/.local/share/jupyter/runtime/notebook_cookie_secret
[I 20:58:34.502 NotebookApp] JupyterLab extension loaded from /usr/local/lib/python3.6/dist-packages/jupyterlab
[I 20:58:34.502 NotebookApp] JupyterLab application directory is /usr/local/share/jupyter/lab
[I 20:58:34.504 NotebookApp] Serving notebooks from local directory: /home/fenics
[I 20:58:34.504 NotebookApp] The Jupyter Notebook is running at:
[I 20:58:34.504 NotebookApp] http://localhost:8765/?token=cbd5069612c7556f3909a4cff02cb8abe0be5a88f4ec5649
[I 20:58:34.504 NotebookApp]  or http://127.0.0.1:8765/?token=cbd5069612c7556f3909a4cff02cb8abe0be5a88f4ec5649
[I 20:58:34.504 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[W 20:58:34.508 NotebookApp] No web browser found: could not locate runnable browser.
[C 20:58:34.508 NotebookApp] 
    
    To access the notebook, open this file in a browser:
        file:///home/fenics/.local/share/jupyter/runtime/nbserver-15-open.html
    Or copy and paste one of these URLs:
        http://localhost:8765/?token=cbd5069612c7556f3909a4cff02cb8abe0be5a88f4ec5649
     or http://127.0.0.1:8765/?token=cbd5069612c7556f3909a4cff02cb8abe0be5a88f4ec5649
```

and use token link to connect with notebook server as usual. 

Building and running with docker should work similarly.
