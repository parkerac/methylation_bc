#install.packages(c("tidyverse", R.utils"))
#BiocManager::install("minfi")
#BiocManager::install("IlluminaHumanMethylation450kmanifest")
#BiocManager::install("GEOquery")

library(tidyverse)
library(R.utils)
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(GEOquery)

dir.create("BC_Files")

# DOWNLOAD RAW DATA --------------------------------------------------------------------------------
# Define function for unpacking files
prepare_files = function(accession) {
  dir.create(accession)
  setwd(accession)
  untar(paste0("../", accession, ".tar"))
  idat_files = list.files(pattern = "idat.gz$", full = TRUE)
  sapply(idat_files, gunzip, overwrite = TRUE)
}

# GSE84207
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84207&format=file", "GSE84207_RAW.tar")
prepare_files("GSE84207_RAW")

# GSE72308
setwd("../")
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE72308&format=file", "GSE72308_RAW.tar")
prepare_files("GSE72308_RAW")

# GSE141338
setwd("../")
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141338&format=file", "GSE141338_RAW.tar")
prepare_files("GSE141338_RAW")

# GSE89278
setwd("../")
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE89278&format=file", "GSE89278_RAW.tar")
prepare_files("GSE89278_RAW")


# SEPARATE CONTROLS --------------------------------------------------------------------------------
setwd("../")

# Define function for downloading series matrix files
download_series_matrix = function(accession, link) {
  file_name = str_c("BC_Files/", accession, "_series_matrix.txt.gz")
  download.file(link, file_name)
  gunzip(file_name)
}

download_series_matrix("GSE89278", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89278/matrix/GSE89278_series_matrix.txt.gz")

# Define function for gathering metadata for controls
get_metadata = function(accession, skip, n_max) {
  sm_file = str_c("BC_Files/", accession, "_series_matrix.txt")
  sm_data = read_tsv(sm_file, col_names = TRUE, skip = skip, n_max = n_max)
  return(sm_data)
}

# Skip header lines that only have 2 columns
# Get data up to the row that identifies controls
GSE89278_metadata = get_metadata("GSE89278", 39, 1)
GSE89278_metadata[2,] = as.list(names(GSE89278_metadata))

# Define function for moving control files to a separate folder
move_controls = function(accession, metadata, identifier) {
  # Identify control files
  names(metadata) = metadata[1,]
  metadata = metadata[nrow(metadata),]
  keep = names(metadata)[str_detect(metadata, identifier)]
  print(keep)
  base_dir = str_c(accession, "_RAW")
  idat_files = list.files(base_dir)
  files_to_keep = idat_files[str_detect(idat_files, paste(keep, collapse = "|"))]
  print(files_to_keep)
  
  # Make directory for control files
  control_dir = str_c(accession, "_Control")
  dir.create(control_dir)
  
  # Move control files to that directory
  new_locations = str_c(control_dir, "/", files_to_keep)
  current_locations = str_c(base_dir, "/", files_to_keep)
  file.rename(current_locations, new_locations)
}

move_controls("GSE89278", GSE89278_metadata, "CONTROL")


# RUN MINFI --------------------------------------------------------------------------------
# Define function for running minfi
run_minfi = function(accession) {
  dir = str_c(accession)
  rg_set = read.metharray.exp(dir, verbose = TRUE)
  m_set = preprocessRaw(rg_set)
  ratio_set = ratioConvert(m_set, what = "both", keepCN = TRUE)
  beta = getBeta(ratio_set)
  beta_frame = as.data.frame(beta, row.names = rownames(beta))
  output_name = str_c("BC_Files/", accession, "_beta.csv")
  write.csv(beta_frame, output_name)
}

run_minfi("GSE141338_RAW")
run_minfi("GSE84207_RAW")
run_minfi("GSE72308_RAW")
run_minfi("GSE89278_RAW_Control")


# SUMMARIZE AT GENE LEVEL --------------------------------------------------------------------------------
setwd("BC_Files/")

# Download this file to your local computer: "https://www.dropbox.com/s/fffnki8j9e3c958/GPL16304-47833.txt.gz?dl=0"
# Load file for mapping probes to genes
metadata450K = read_tsv("GPL16304-47833.txt", comment="#") %>%
  dplyr::select(ID, Distance_closest_TSS, Closest_TSS_gene_name) %>%
  dplyr::rename(Probe=ID, TSS_Distance=Distance_closest_TSS, Gene=Closest_TSS_gene_name) %>% 
  arrange(Gene, TSS_Distance)

# Define function for summarizing beta values
perform_summary = function(accession) {
  # Load data
  beta_file = str_c(accession, "_beta_all.csv")
  rawData450K = read_csv(beta_file) %>% 
    dplyr::rename(Probe=X1)
  
  # Set variables for separating data into 10 patient sections
  num_patients = ncol(rawData450K)
  sections = list()
  start_index = 2
  stop_index = 11
  probes = rawData450K[,1]
  i = 1
  
  # Separate the data to avoid overloading memory
  while (stop_index < num_patients) {
    section = rawData450K[,start_index:stop_index]
    complete_section = cbind(probes, section)
    sections[[i]] = as.data.frame(complete_section)
    start_index = start_index + 10
    stop_index = stop_index + 10
    i = i + 1
  }
  
  # Add final section of data
  section = rawData450K[,start_index:num_patients]
  complete_section = cbind(probes, section)
  sections[[i]] = as.data.frame(complete_section)
  
  # Set variables for collecting data
  output = list()
  j = 1
  
  # Process each subsection of data
  for (section in sections) {
    data450K = inner_join(metadata450K, section, by = "Probe")
    tidyData450K = gather(data450K, Patient_ID, beta, -Probe, -TSS_Distance, -Gene, na.rm = TRUE)
    
    # Filter data to include only probes closest to TSS
    maxTSSDistance = 300
    filteredData450K = dplyr::filter(tidyData450K, TSS_Distance > -maxTSSDistance) %>%
      dplyr::filter(TSS_Distance < maxTSSDistance)
    
    # Summarize to mean value per gene
    filteredData450K = filteredData450K %>% dplyr::select(Patient_ID, Probe, Gene, beta) %>% group_by(Patient_ID, Gene) %>% summarise(Value=mean(beta)) %>% ungroup()
    output[[j]] = filteredData450K
    j = j + 1
  }
  
  # Combine data and save
  final_output = reduce(output, rbind)
  output_name = str_c(accession, "_Methylation_Analysis_Data.txt")
  write_tsv(final_output, output_name)
}

perform_summary("GSE141338_RAW")
perform_summary("GSE84207_RAW")
perform_summary("GSE72308_RAW")
perform_summary("GSE89278_RAW_Control")


# IDENTIFY SUBTYPE CATEGORIES --------------------------------------------------------------------------------
download_series_matrix("GSE141338", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141338/matrix/GSE141338_series_matrix.txt.gz")
gse = getGEO(filename = "GSE141338_series_matrix.txt")

metadata = data.frame(Category = gse$title, Accession = gse$geo_accession)
write_tsv(metadata, "classifications.tsv")

