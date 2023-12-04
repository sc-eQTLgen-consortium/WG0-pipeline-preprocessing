
################## BASE IMAGE ######################

FROM ubuntu:22.04

################## METADATA ######################

LABEL base_image="ubuntu:22.04"
LABEL version="1.0.0"
LABEL software="WG0 Pipeline"
LABEL about.summary="WG0 sceQTLGen Consortium Prepreocessing Pipeline"
LABEL about.documentation="https://github.com/sc-eQTLgen-consortium/WG0-pipeline-preprocessing"
LABEL about.tags="Genomics"

# Build syntax: docker build ./ -t wg0-pipeline-preprocessing:2023.11.28.0 --progress=plain > build.log 2>&1
# Total build takes 6 minutes and has a size of 6.66 GB.
# Use dive wg0-pipeline-preprocessing:2023.11.28.0 to investigate memory usage.

################## MAINTAINER ######################

MAINTAINER Martijn Vochteloo <m.vochteloo@umcg.nl>

################## INSTALLATION ######################

ADD . /tmp/repo
WORKDIR /tmp/repo

ENV PATH=/opt:/usr/games:/opt/conda/envs/py37/bin/:/opt/conda/bin:/opt/cellranger-7.0.1:$PATH
ENV SHELL=/bin/bash
ENV LC_ALL=C
ENV LANG=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN echo 'alias python=python3' >> ~/.bashrc

# Needed to prevent asking for geographic location when installing things.
RUN export TZ=Europe/Amsterdam \
    && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

RUN apt-get update -y \
    && apt-get upgrade -y \
    # ca-certificates openssl
    && apt-get install -y --no-install-recommends ca-certificates \
    # libpsl5 wget
    && apt-get install -y --no-install-recommends wget \
    # git git-man libbrotli1 libcurl3-gnutls liberror-perl libexpat1
    # libgdbm-compat4 libgdbm6 libldap-2.5-0 libnghttp2-14 libperl5.34 librtmp1
    # libsasl2-2 libsasl2-modules-db libssh-4 perl perl-modules-5.34
    && apt-get install -y --no-install-recommends git \
    # unzip
    && apt-get install -y --no-install-recommends unzip

##################################
############# CELLRANGER #########
##################################

RUN cd /opt \
    && wget http://regmedsrv1.wustl.edu/Public_SPACE/litd/Public_html/pkg/cellranger-7.0.1.tar.gz \
    && tar -xzvf cellranger-7.0.1.tar.gz \
    && rm cellranger-7.0.1.tar.gz

##################################
############# CELLBENDER #########
##################################

# Reduce conda size by preventing Python from recreating a corresponding bytecode cache file (*.pyc) at runtime.
ENV PYTHONDONTWRITEBYTECODE=true

# Install Python. Also required for bedtools2.
# libexpat1 libmpdec3 libpython3-stdlib libpython3.10-minimal
# libpython3.10-stdlib libreadline8 libsqlite3-0 media-types python3
# python3-minimal python3.10 python3.10-minimal readline-common
RUN apt-get install -y --no-install-recommends python3

# Install miniconda for the virtual environment.
# https://github.com/ContinuumIO/docker-images/blob/main/miniconda3/debian/Dockerfile
RUN cd /opt \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-py37_23.1.0-1-Linux-x86_64.sh -O miniconda.sh -q \
    && mkdir -p /opt \
    && bash miniconda.sh -b -p /opt/conda \
    && rm miniconda.sh \
    && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
    && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && /opt/conda/bin/conda clean -afy

# Create and activate virtual environment
RUN eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)" \
    && conda create -n py37 python=3.7.16 \
    # None
    && /opt/conda/envs/py37/bin/pip install anaconda==0.0.1.1 \
    # numexpr-2.8.6 numpy-1.21.6 packaging-23.2
    && /opt/conda/envs/py37/bin/pip install tables==3.7.0 \
    # nvidia-cublas-cu11-11.10.3.66 nvidia-cuda-nvrtc-cu11-11.7.99  \
    # nvidia-cuda-runtime-cu11-11.7.99 nvidia-cudnn-cu11-8.5.0.96  \
    # typing-extensions-4.7.1
    && /opt/conda/envs/py37/bin/pip install torch==1.13.1 \
    # charset-normalizer-3.3.2 idna-3.4 pillow-9.5.0 requests-2.31.0  \
    # urllib3-2.0.7
    && /opt/conda/envs/py37/bin/pip install torchvision==0.14.1 \
    # Using Commits on Aug 8, 2023.
    # MarkupSafe-2.1.3 Send2Trash-1.8.2 anndata-0.8.0 anyio-3.7.1 argon2-cffi-23.1.0 argon2-cffi-bindings-21.2.0
    # attrs-23.1.0 backcall-0.2.0 beautifulsoup4-4.12.2 bleach-6.0.0 cellbender-0.3.0 cffi-1.15.1 click-8.1.7
    # comm-0.1.4 cycler-0.11.0 debugpy-1.7.0 decorator-5.1.1 defusedxml-0.7.1 entrypoints-0.4 exceptiongroup-1.2.0
    # fastjsonschema-2.19.0 fonttools-4.38.0 h5py-3.8.0 importlib-resources-5.12.0 importlib_metadata-6.7.0
    # ipykernel-6.16.2 ipython-7.34.0 ipython-genutils-0.2.0 ipywidgets-8.1.1 jedi-0.19.1 jinja2-3.1.2 jsonschema-4.17.3
    # jupyter-1.0.0 jupyter-client-7.4.9 jupyter-console-6.6.3 jupyter-core-4.12.0 jupyter-server-1.24.0
    # jupyter_contrib_core-0.4.2 jupyter_contrib_nbextensions-0.7.0 jupyter_highlight_selected_word-0.2.0
    # jupyter_nbextensions_configurator-0.6.3 jupyterlab-pygments-0.2.2 jupyterlab-widgets-3.0.9 kiwisolver-1.4.5
    # llvmlite-0.39.1 loompy-3.0.7 lxml-4.9.3 matplotlib-3.5.3 matplotlib-inline-0.1.6 mistune-0.8.4 natsort-8.4.0
    # nbclassic-1.0.0 nbclient-0.7.4 nbconvert-6.5.4 nbformat-5.8.0 nest-asyncio-1.5.8 notebook-6.5.6 notebook-shim-0.2.3
    # numba-0.56.4 numpy-groupies-0.9.22 opt-einsum-3.3.0 pandas-1.3.5 pandocfilters-1.5.0 parso-0.8.3 pexpect-4.8.0
    # pickleshare-0.7.5 pkgutil-resolve-name-1.3.10 prometheus-client-0.17.1 prompt-toolkit-3.0.41 psutil-5.9.6
    # ptyprocess-0.7.0 pycparser-2.21 pygments-2.17.2 pyparsing-3.1.1 pyro-api-0.1.2 pyro-ppl-1.8.6 pyrsistent-0.19.3
    # python-dateutil-2.8.2 pytz-2023.3.post1 pyyaml-6.0.1 pyzmq-24.0.1 qtconsole-5.4.4 qtpy-2.4.1 scipy-1.7.3 six-1.16.0
    # sniffio-1.3.0 soupsieve-2.4.1 terminado-0.17.1 tinycss2-1.2.1 tornado-6.2 tqdm-4.66.1 traitlets-5.9.0 wcwidth-0.2.12
    # webencodings-0.5.1 websocket-client-1.6.1 widgetsnbextension-4.0.9 zipp-3.15.0
    && /opt/conda/envs/py37/bin/pip install --no-cache-dir -U git+https://github.com/broadinstitute/CellBender.git@4990df713f296256577c92cab3314daeeca0f3d7 \
    && /opt/conda/envs/py37/bin/pip cache purge

RUN conda clean -y --all

#######################################
################ WG1-CODE #############
#######################################

# Always get our own newest software. IMPORTANT: make sure you use the correct branch here.
# Delete all the non python / R files since we won't need them anyway.
RUN cd /opt \
    && GITHUB_BRANCH=scMetaBrain \
    && wget https://github.com/sc-eQTLgen-consortium/WG0-pipeline-preprocessing/archive/refs/heads/${GITHUB_BRANCH}.zip \
    && unzip -q ${GITHUB_BRANCH}.zip \
    && rm ${GITHUB_BRANCH}.zip \
    && mv WG0-pipeline-preprocessing-${GITHUB_BRANCH} WG0-pipeline-preprocessing \
    && cd WG0-pipeline-preprocessing \
    && find . -type f ! \( -iname \*.py -o -iname \*.R \) -delete \
    && find . -type d -empty -delete \
    && find . -type f \( -iname \*.py -o -iname \*.R \) -exec chmod 777 {} \;

####################################
################ CLEAN #############
####################################

RUN apt-get clean \
    && apt-get autoremove -y
