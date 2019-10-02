# Environments

Each process has its own environment, either managed by `nextflow` via `conda` or `docker`.

To create and run the environments for yourself, follow these commands:

## conda

```sh
conda env create -f assemblyWithMiniasm/environment.yaml -n miniasm
conda activate miniasm
```

## docker

```sh
cd assemblyWithMiniasm
docker build -t willrowe/miniasm .
docker push willrowe/miniasm
docker run -it willrowe/miniasm
```

> replace willrowe with your docker hub username
