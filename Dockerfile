# Dockerfile -- the single image every Nextflow process runs in.
#
# Two separate environments, each on its own interpreter:
#   - Pacybara (rothlab/pacybara at PACYBARA_COMMIT) and its toolchain in a
#     conda env built from Pacybara's own pacybara_env.yml (Python 3.10, R,
#     bwa, bowtie2, samtools, emboss, muscle). Its bin/ is first on PATH
#     (the env is never activated), so pacybara_simplex.sh's
#     `python3`/`Rscript` and its scripts' `#!/usr/bin/env` shebangs all
#     resolve there.
#   - This repo's package in a uv venv (Python 3.13). Only its
#     `pacybara-make-params` console script is put on PATH; that script's
#     shebang points at the venv's own interpreter, so the two Pythons never
#     meet.
#
# linux/amd64 only: not every bioconda package in pacybara_env.yml has an
# arm64 build. On Apple silicon, build with --platform linux/amd64.
FROM mambaorg/micromamba:2.3.0

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates git procps \
    && rm -rf /var/lib/apt/lists/*

# ── Pacybara ────────────────────────────────────────────────────────────────
# Bump deliberately: rebuild and rerun a known library when changing it.
ARG PACYBARA_COMMIT=b5e39bd0e2731283eee2ee88c0e0852ddf83fe33
# Pacybara's GitHub-only R dependencies (install.sh installs whatever their
# default branch is; pinned here instead).
ARG YOGITOOLS_COMMIT=35265272206d2e3daed1b8b0cc11da4c17166d49
ARG YOGISEQ_COMMIT=2363474193fc6e6d474a45579fe13f9729ad74f8
ARG HGVSPARSER_COMMIT=7a3a5c81ac0037a4514225924aa108d663b1fb90

RUN git clone https://github.com/rothlab/pacybara.git /opt/pacybara \
    && git -C /opt/pacybara checkout --quiet "${PACYBARA_COMMIT}"

# install.sh is interactive, so its steps are reproduced here: the conda env
# from pacybara_env.yml, then its R packages. The CRAN ones come from
# conda-forge as prebuilt binaries (no compiler needed in the image); only
# the three GitHub packages are installed from source, without touching
# their already-installed dependencies.
RUN micromamba create -y -n pacybara -f /opt/pacybara/pacybara_env.yml \
    && micromamba install -y -n pacybara -c conda-forge \
        r-argparser r-hash r-bitops r-pbmcapply r-remotes \
    && micromamba clean -a -y

ENV PATH="/opt/conda/envs/pacybara/bin:${PATH}"

RUN Rscript -e " \
    options(warn = 2); \
    for (spec in c('jweile/yogitools@${YOGITOOLS_COMMIT}', \
                   'jweile/yogiseq@${YOGISEQ_COMMIT}', \
                   'VariantEffect/hgvsParseR@${HGVSPARSER_COMMIT}')) \
      remotes::install_github(spec, dependencies = FALSE, upgrade = 'never'); \
    for (p in c('argparser', 'hash', 'bitops', 'pbmcapply', \
                'yogitools', 'yogiseq', 'hgvsParseR')) \
      library(p, character.only = TRUE)"

# install.sh's last step: copy the scripts onto PATH.
RUN cp /opt/pacybara/src/*.sh /opt/pacybara/src/*.R /opt/pacybara/src/*.py /usr/local/bin/ \
    && chmod 755 /usr/local/bin/*.sh /usr/local/bin/*.R /usr/local/bin/*.py

# ── This repo's package ─────────────────────────────────────────────────────
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/local/bin/
# Keep uv's managed Python out of /root, so tasks run as the invoking user
# (nextflow.config's `-u $(id -u):$(id -g)`) can still execute it.
ENV UV_PYTHON_INSTALL_DIR=/opt/uv-python \
    UV_LINK_MODE=copy \
    UV_COMPILE_BYTECODE=1

WORKDIR /opt/pacybara-workflow
# .python-version must be present before the first `uv sync`, so uv
# provisions exactly 3.13. Dependencies first, cached independently of
# source changes (there are none today, but the layer split costs nothing).
COPY pyproject.toml uv.lock .python-version ./
RUN uv sync --frozen --no-install-project --no-dev
COPY README.md ./
COPY src/ src/
RUN uv sync --frozen --no-dev \
    && ln -s /opt/pacybara-workflow/.venv/bin/pacybara-make-params /usr/local/bin/

# pacybara_simplex.sh runs under `set -u` and reads $CONDA_DEFAULT_ENV,
# which is unset when the env is on PATH rather than activated. "base" takes
# its "No conda environment detected" branch -- harmless, since every tool
# already resolves through PATH.
ENV CONDA_DEFAULT_ENV=base

WORKDIR /
USER mambauser
ENTRYPOINT []
CMD []
