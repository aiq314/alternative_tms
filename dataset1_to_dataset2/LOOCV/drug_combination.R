# Load the gtools package
install.packages("gtools")
install.packages("comprehenr")
library(gtools)
library(comprehenr)


# Define dict of drug_risks
elements <- c("mexiletine","loratadine","metoprolol","ranolazine","diltiazem","tamoxifen","nitrendipine",
              "droperidol","cisapride","risperidone","astemizole","clozapine","domperidone",
              "sotalol","dofetilide","disopyramide","bepridil","vandetanib")

combine_drug <- function(num_of_drugs) {
  # Number of elements to choose (excluding L1 and L2)
  # rand_elements <- num_of_drugs-2 # 8 drugs = 2 fixed and 6 random
  rand_elements <- num_of_drugs # All random
  
  # Generate combinations with repetitions for the remaining elements
  combinations <- combinations(n = length(elements), r = rand_elements, v = elements, repeats.allowed = FALSE)
  
  # # Filter combinations to ensure no more than 3 L's in each combination
  # #filtered_combinations <- subset(combinations, rowSums(combinations == "L1" | combinations == "L2") <= 3)
  check_counts <- function(combination) {
  	 # 8 drugs = L<2; 2<I<3; 2<H<3 
    max_drugs <- num_of_drugs / 3
    # sum(combination %in% c("tamoxifen", "loratadine")) >= floor(max_drugs)-2 &&
    # sum(combination %in% c("tamoxifen", "loratadine")) < ceiling(max_drugs)-1 &&
    sum(combination %in% c("mexiletine","loratadine","metoprolol","ranolazine","diltiazem","tamoxifen","nitrendipine")) >= floor(max_drugs) &&
    sum(combination %in% c("mexiletine","loratadine","metoprolol","ranolazine","diltiazem","tamoxifen","nitrendipine")) <= ceiling(max_drugs) &&
    sum(combination %in% c("droperidol","cisapride","risperidone","astemizole","clozapine","domperidone")) >= floor(max_drugs) &&
    sum(combination %in% c("droperidol","cisapride","risperidone","astemizole","clozapine","domperidone")) <= ceiling(max_drugs) &&
    sum(combination %in% c("sotalol","dofetilide","disopyramide","bepridil","vandetanib")) >= floor(max_drugs) &&
    sum(combination %in% c("sotalol","dofetilide","disopyramide","bepridil","vandetanib")) <= ceiling(max_drugs)
  }
  # Filter combinations to ensure no more than 3 H in each combination
  filtered_combinations <- subset(combinations, apply(combinations, 1, check_counts))
  # # Add L1 and L2 to the filtered combinations
  # final_combinations <- cbind("verapamil","ranolazine",filtered_combinations)
  final_combinations <- cbind(filtered_combinations)
  return(final_combinations)
}