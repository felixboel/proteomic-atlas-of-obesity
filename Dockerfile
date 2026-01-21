FROM rocker/r-ver:4.3.2

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    ca-certificates \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /srv/shiny-app

COPY renv.lock renv.lock
COPY renv ./renv
COPY .Rprofile .Rprofile

RUN R -e "install.packages('renv', repos='https://cloud.r-project.org'); renv::restore(prompt = FALSE)"

COPY app ./app
COPY data/processed ./data/processed

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/srv/shiny-app/app', host='0.0.0.0', port=3838)"]
