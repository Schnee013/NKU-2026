FROM python:3.12.14-slim@sha256:f77ac9e44ae96ef2c90b8053ea08c31f8be030f824196b0ae4db6d462c84e51f
ARG REQUIREMENTS=pathway-requirements.txt
RUN --mount=type=secret,id=proxy_ca \
    if [ "$REQUIREMENTS" = "docking-requirements.txt" ]; then \
      if [ -f /run/secrets/proxy_ca ]; then \
        apt-get -o Acquire::https::CaInfo=/run/secrets/proxy_ca update; \
      else apt-get update; fi && \
      apt-get install -y --no-install-recommends libxrender1 libxext6 libx11-6 libexpat1 && \
      rm -rf /var/lib/apt/lists/*; \
    fi
COPY ${REQUIREMENTS} /requirements.txt
RUN --mount=type=secret,id=proxy_ca \
    if [ -f /run/secrets/proxy_ca ]; then \
      PIP_CERT=/run/secrets/proxy_ca python -m pip install --disable-pip-version-check --no-cache-dir --root-user-action=ignore -r /requirements.txt; \
    else \
      python -m pip install --disable-pip-version-check --no-cache-dir --root-user-action=ignore -r /requirements.txt; \
    fi
WORKDIR /app
ENV PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/matplotlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
