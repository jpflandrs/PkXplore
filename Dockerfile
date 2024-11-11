#FromPlatformFlagConstDisallowed: FROM --platform flag should not use constant value "linux/amd64"
FROM julia:latest 

# Create user and set up directories
RUN useradd --create-home --shell /bin/bash genie
RUN mkdir /home/genie/app
COPY . /home/genie/app
WORKDIR /home/genie/app
COPY Project.toml /home/genie/app/

# Install system dependencies and clean up in one RUN to reduce layers
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    mafft \
    r-base \
    r-base-dev \
    wget \
    unzip \
    g++ \
    make 
RUN rm -rf /var/lib/apt/lists/*

# Install R packages RUN apt-get update && apt-get install -y r-base r-base-dev
RUN R -e "install.packages(\"ape\")" 
RUN R -e "install.packages(\"phangorn\")"
   
RUN apt-get update && apt-get install -y wget \
    unzip \
    g++ \
    make 
# Download, compile and install rapidnj
 
RUN    wget -O rapidnj.zip https://github.com/somme89/rapidNJ/archive/refs/heads/master.zip
RUN    unzip rapidnj.zip 
WORKDIR rapidNJ-master 
RUN    make 
RUN    cp bin/* /usr/local/bin/ 
WORKDIR /home/genie/app  
RUN    rm -rf rapid*

# Set ownership
RUN chown -R genie:genie /home/

# Switch to genie user
USER genie

# Configure ports
EXPOSE 8000
EXPOSE 80

# Set environment variables LegacyKeyValueFormat: "ENV key=value" should be used instead of legacy "ENV key value" format
ENV JULIA_DEPOT_PATH="/home/genie/.julia"
ENV JULIA_REVISE="off"
ENV GENIE_ENV="prod"
ENV GENIE_HOST="0.0.0.0"
ENV PORT="8000"
ENV WSPORT="8000"
ENV EARLYBIND="true"

# Install Julia packages
RUN julia -e "using Pkg; Pkg.activate(\".\"); Pkg.instantiate(); Pkg.precompile();"

ENTRYPOINT ["julia", "--project", "-e", "using GenieFramework; Genie.loadapp(); up(async=false);"]
