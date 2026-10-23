# NGS Pipeline Dockerfile
# Multi-stage build for minimal production image

# ============================================================================
# BASE STAGE - Build dependencies
# ============================================================================
FROM ubuntu:22.04 AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    unzip \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python3 \
    python3-pip \
    openjdk-11-jdk-headless \
    && rm -rf /var/lib/apt/lists/*

# Set Java home
ENV JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
ENV PATH=$JAVA_HOME/bin:$PATH

# ============================================================================
# TOOLS STAGE - Install bioinformatics tools
# ============================================================================
FROM builder AS tools

# Install conda
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda clean -afy

ENV PATH=/opt/conda/bin:$PATH

# Configure conda channels
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# Install bioinformatics tools via conda
RUN conda install -y \
    fastqc=0.11.9 \
    multiqc=1.14 \
    cutadapt=4.4 \
    fastp=0.23.2 \
    trimmomatic=0.39 \
    bwa=0.7.17 \
    bowtie2=2.5.1 \
    samtools=1.17 \
    gatk4=4.5.0.0 \
    freebayes=1.3.6 \
    bcftools=1.17 \
    snpeff=5.1d \
    vep=108.0 \
    bedtools=2.31.0 \
    htslib=1.17 \
    python=3.11 \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    && conda clean -afy

# Install additional Python packages
RUN pip install --no-cache-dir \
    multiqc[pdf]==1.14 \
    pandas==2.1.0 \
    numpy==1.24.3 \
    scipy==1.11.0

# ============================================================================
# VEP CACHE SETUP (Optional - download during build)
# ============================================================================
# VEP cache is large (~20GB), typically mounted at runtime
# Uncomment below to pre-download (increases image size significantly)
# RUN vep_install -a cf -s homo_sapiens -y GRCh38 -c /opt/vep_cache --CONVERT

# ============================================================================
# SNPEFF DATABASE SETUP
# ============================================================================
RUN mkdir -p /opt/snpEff/data && \
    echo "database.dir=/opt/snpEff/data" > /opt/snpEff/snpEff.config

# ============================================================================
# FINAL STAGE - Minimal runtime image
# ============================================================================
FROM ubuntu:22.04 AS runtime

# Install runtime dependencies only
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-11-jre-headless \
    python3 \
    python3-pip \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libcurl4 \
    libssl3 \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Copy tools from tools stage
COPY --from=tools /opt/conda /opt/conda
COPY --from=tools /opt/snpEff /opt/snpEff

ENV PATH=/opt/conda/bin:$PATH
ENV JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
ENV PATH=$JAVA_HOME/bin:$PATH

# Create non-root user
RUN groupadd -r ngsuser -g 1000 && \
    useradd -r -g ngsuser -u 1000 -m -d /home/ngsuser -s /bin/bash ngsuser

# Set up working directories
RUN mkdir -p /data /results /work /resources && \
    chown -R ngsuser:ngsuser /data /results /work /resources /opt/snpEff

USER ngsuser
WORKDIR /work

# Default entrypoint
ENTRYPOINT ["nextflow"]
CMD ["run", "/pipeline/main.nf", "--help"]

# Labels
LABEL maintainer="NGS Pipeline Team" \
      version="1.0.0" \
      description="NGS Analysis Pipeline Docker Image" \
      org.opencontainers.image.source="https://github.com/ngs-pipeline/ngs-pipeline" \
      org.opencontainers.image.documentation="https://github.com/ngs-pipeline/ngs-pipeline/blob/main/README.md"