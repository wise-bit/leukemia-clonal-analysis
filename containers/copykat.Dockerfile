# Latest official R runtime for parent image
FROM r-base:latest

# Installing required R packages
RUN R -e "install.packages('Seurat')"
RUN R -e "install.packages('devtools')
RUN R -e "devtools::install_github('navinlabcode/copykat')"

# Copy R script and other necessary files into container
COPY scripts-r/copykat-workflow.R /usr/src/T-ALL/copykat-workflow.R

# Set working directory
WORKDIR /usr/src/T-ALL/

# Run R script when the container launches
CMD ["Rscript", "copykat-workflow.R"]
