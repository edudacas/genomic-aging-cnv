# Candidate build recipe. Release requires successful build, native tests,
# resolved dependency manifests and recorded image ID. Never use an unpinned base.
ARG BASE_IMAGE
FROM ${BASE_IMAGE}
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates python3 python3-venv python3-dev build-essential gfortran \
    samtools libhts-dev r-base r-base-dev r-cran-dplyr r-cran-data.table \
    r-cran-ggplot2 r-cran-argparse r-cran-stringr r-cran-patchwork r-cran-glue \
    r-cran-r.utils libcurl4-openssl-dev libssl-dev libxml2-dev
WORKDIR /opt/cnv
COPY . /opt/cnv
RUN python3 runtime/bootstrap_local.py /opt/cnv-runtime
ENV PATH=/opt/cnv-runtime/python/bin:/opt/cnv-runtime/bin:$PATH
ENV R_LIBS_USER=/opt/cnv-runtime/R-library R_LIBS_SITE=/opt/cnv-runtime/R-library
RUN python3 runtime/validate_native.py --root /opt/cnv-runtime --output /opt/cnv-validation
RUN python3 runtime/check_event_signals.py /opt/cnv-validation
RUN python3 runtime/freeze_runtime.py --root /opt/cnv-runtime --output /opt/cnv-runtime/FILES.json
RUN dpkg-query -W > /opt/cnv-runtime/os-packages.tsv
ENTRYPOINT ["/opt/cnv-runtime/python/bin/python", "/opt/cnv/cnv_cli.py"]
CMD ["--help"]
