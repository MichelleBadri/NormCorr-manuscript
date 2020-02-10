## Docker container recipe for reproducing results in this repo
FROM continuumio/miniconda3:4.7.12
LABEL maintainer="Zachary D Kurtz <zdkurtz@gmail.com>"

COPY environment.yml ./
RUN conda env update -n base --file environment.yml --prune

ENTRYPOINT [ "/opt/conda/bin/tini", "--" ]
CMD [ "/bin/bash" ]
