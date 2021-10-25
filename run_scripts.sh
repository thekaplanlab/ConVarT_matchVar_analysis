#!/bin/bash

# Make scripts executable
chmod +x downloads.sh
chmod +x topmed.sh

# Run local scripts
Rscript helper.R
./downloads.sh
./topmed.sh
Rscript variant_tables.R
Rscript matchvar_tables.R
Rscript write_sql.R
Rscript matchVar_analysis.R
