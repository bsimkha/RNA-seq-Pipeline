Before you begin make sure these dependenceies are poperly installed
  - NextFlow (Refuqires java version 17-25)
  - FastP

Step 1: Edit the config.yaml files to proper file paths

Step 2 (Optional): Run the Preflight.R file to make sure all the pair ended reads have proper mates
  - Run using Rscript Preflight.R command
  - Once you have established all the reads have pairs, proceed to Step 3

Step 3: Run main.nf with nextflow run main.nf script

