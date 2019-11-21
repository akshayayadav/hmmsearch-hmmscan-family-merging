### Container for HMMSEARCH-HMMSCAN based family merging workflow

This container contains the application for merging gene families using two-way HMMSEARCH-HMMSCAN search procedure. Each given family is searched against a combined database of all family sequences and the desired outgroup sequences to find the candidate missing sequences for the family. The candidate missing sequences reverse searched against the combined database of family HMMs and the sequences that find the original family as the best match are predicted missing sequences for the family. Finally, the predicted missing sequences from all the famililies are used to merge families.

 1. #### Downloading the container
  ```
  docker pull akshayayadav/hmmsearch-hmmscan-family-merging
  ```

 2. #### Preparing the data
 Create a data  directory *<my_directory>* with a user defined name. This directory **MUST** contain 3 items: *family_fasta* directory containing all the fasta files for the user-defined gene families, a *family-outgroup-sequences.fa* fasta file containing sequences from all the families and sequences from all the desired outgroups and a *outgroups.list* file listing the names of all the outgroup species (one per line) present in the *family-outgroup-sequences.fa* file. All the outgroup sequence headers in the *family-outgroup-sequences.fa* file **MUST** be prefixed with the respective outgroup species names present in *outgroups.list* file.
 

 3. #### Running the analysis
 Once the data is ready execute the workflow using the following command
 ```
docker run -v <absolute_path_to_the_data_directory>:/data akshayayadav/hmmsearch-hmmscan-family-merging run_analysis.sh
 ```
