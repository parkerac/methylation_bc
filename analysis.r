source('utils.r')


"TEST 1"

"Prepare Data"
raw_normal = read_tsv('GSE89278_RAW_Control_Methylation_Analysis_Data.txt')
raw_general = read_tsv('GSE72308_RAW_Methylation_Analysis_Data.txt')

normal_test_data = mutate(raw_normal, Category = 'Normal')
general_test_data = mutate(raw_general, Category = 'General')


processed_test_normal = process_data(normal_test_data, 'Normal')
processed_general = process_data(general_test_data, 'General')


"Compare First Dataset to Normal"
general_vs_normal = compareTwo(processed_test_normal, processed_general)
general_vs_normal = rename(general_vs_normal, general_diff = diff)
significant_genes_test_1 = compareTwoCategories(processed_test_normal, normal_test_data, processed_general, general_test_data, .01, "general_diff")

#Confirm Results
"Prepare Data"
raw_confirm = read_tsv('GSE84207_RAW_Methylation_Analysis_Data.txt')
confirm_test_data = mutate(raw_confirm, Category = 'Confirm')
processed_confirm = process_data(confirm_test_data, 'Confirm')

"Compare Second Dataset to Normal"
confirm_vs_normal = compareTwo(processed_test_normal, processed_confirm)
confirm_vs_normal = rename(confirm_vs_normal, confirm_diff = diff)
significant_genes_test_2 =  compareTwoCategories(processed_test_normal, normal_test_data, processed_confirm, confirm_test_data, .01, "confirm_diff")

"Combine Results"
significant_genes_test_1 = rename(significant_genes_test_1, p_value_test_1 = p_value)
significant_genes_test_2 = rename(significant_genes_test_2, p_value_test_2 = p_value)

joined_final_results = inner_join(significant_genes_test_1,significant_genes_test_2)
joined = inner_join(joined_final_results, general_vs_normal)
joined_final_results_confirmation_test = inner_join(joined, confirm_vs_normal)
print(joined_final_results_confirmation_test)
significant_genes_final = pull(joined_final_results_confirmation_test, Gene)

"TEST 2 - sub-types"

classes = read_tsv('classifications.tsv')
raw_multiple_groups_data = read_tsv('GSE141338_RAW_Methylation_Analysis_Data.txt') %>% separate(Patient_ID, c("Accession", "Patient_ID"), "_", extra = "merge")
pre_name_joined_data = inner_join(classes, raw_multiple_groups_data)
joined_data = rename_categories(pre_name_joined_data)  %>% mutate(Category = as.factor(Category))

"Prepare Data"
HER2 = filter(joined_data, Category == 'HER2')
LuminalA = filter(joined_data, Category == 'Luminal A')
LuminalB = filter(joined_data, Category == 'Luminal B')
luminalHER2 = filter(joined_data, Category =='Luminal-HER2')
normalBreast = filter(joined_data, Category == 'Normal')
TN = filter(joined_data, Category =='TN')
inVitro = filter(joined_data, Category == 'In Vitro')

"Process Data"
processedHER2 = process_data(HER2, 'HER2')
processedLuminalA = process_data(LuminalA, 'LuminalA')
processedLuminalB = process_data(LuminalB, 'LuminalB')
processedLuminalHER2 = process_data(luminalHER2, 'LuminalHER2')
processedNormal = process_data(normalBreast, 'Normal')
processedTN = process_data(TN, 'TN')
processedInVitro = process_data(inVitro, 'inVitro')

processedAllSubTypeDataSet = list(processedHER2, processedLuminalA, 
                                processedLuminalB, processedLuminalHER2,
                                processedNormal, processedTN, processedInVitro) %>%
                                reduce(rbind)

"Compare all versus normal"
HER2_NORMAL = compareTwo(processedHER2, processedNormal)
TN_NORMAL = compareTwo(processedTN, processedNormal)
LuminalA_NORMAL = compareTwo(processedLuminalA, processedNormal)
LuminalB_NORMAL = compareTwo(processedLuminalB, processedNormal)
LuminalHER2_NORMAL = compareTwo(processedLuminalHER2, processedNormal)
InVitro_NORMAL = compareTwo(processedInVitro, processedNormal)

HER2_NORMAL = rename(HER2_NORMAL, HER_diff = diff)
TN_NORMAL = rename(TN_NORMAL, TN_diff = diff)
LuminalA_NORMAL = rename(LuminalA_NORMAL, LuminalA_diff = diff)
LuminalB_NORMAL = rename(LuminalB_NORMAL, LuminalB_diff = diff)
LuminalHER2_NORMAL = rename(LuminalHER2_NORMAL, LuminalHER2_diff = diff)

"Discover Interesting Genes"
all_joined = list(HER2_NORMAL, TN_NORMAL, LuminalA_NORMAL, LuminalB_NORMAL, LuminalHER2_NORMAL) %>%
  reduce(inner_join)
view(all_joined)
genes = all_joined %>% pull(Gene)
common_genes_all = filter(joined_data, Gene %in% genes) %>% filter(Category != "In Vitro")

"STATISTICAL ANALYSIS"
print(paste0("ALPHA ", .05 / length(genes)))
significant_genes_all_types = performStatisticalAnalysis(common_genes_all, genes, .05 / length(genes))
print(significant_genes_all_types)

"TEST GENES NOT IN ORIGNIAL SIGNIFICANT SET"

genes_not_in_original = NULL
for (gene in genes){
  if (!(gene %in% significant_genes_final)){
    print(paste0(c(gene, " not found significant in original test")))
    genes_not_in_original = c(genes_not_in_original, gene)
  }
  
}
print(genes_not_in_original)
significant_genes_not_significant_in_test_1 = genes_not_in_original
"Perform Statistical Analysis against original data"

"Test Data Set"

filtered_normal_processed = filter(processed_test_normal, Gene %in% genes)
filtered_general_processed = filter(processed_general, Gene %in% genes)
significant_genes_test_3_general =  compareTwoCategories(processed_test_normal, 
                                            normal_test_data, 
                                            filtered_general_processed,
                                            general_test_data, .01, "General_diff", .0)

"Test Confirm Data Set"
filtered_confirm_processed = filter(processed_confirm, Gene %in% genes)
significant_genes_test_3_confirm =  compareTwoCategories(processed_test_normal, 
                                            normal_test_data, 
                                            filtered_confirm_processed,
                                            confirm_test_data, .01, "Confirm_diff", .0)
compared_to_large = inner_join(significant_genes_test_3_general, significant_genes_test_3_confirm)
print(compared_to_large)
write_tsv(compared_to_large, 'out/compated_to_large.tsv')
# general_vs_normal = compareTwo(processed_test_normal, processed_general, threshold=0)
# print(filter(general_vs_normal, Gene %in% genes_not_in_original))
# confirm_vs_normal = compareTwo(processed_test_normal, processed_confirm, threshold=0)
# print(filter(confirm_vs_normal, Gene %in% genes_not_in_original))


"COMMON GENES WITHOUT HER2"
all_joined_not_her2 = list(TN_NORMAL, LuminalA_NORMAL, LuminalB_NORMAL, LuminalHER2_NORMAL) %>%
  reduce(inner_join)
view(all_joined_not_her2)
genes = all_joined_not_her2 %>% pull(Gene)
common_genes_not_her2 = filter(joined_data, Gene %in% genes) %>% filter(Category != "HER2") %>% filter(Category != "In Vitro")

"STATISTICAL ANALYSIS"
print(paste0("ALPHA ", .05 / length(genes)))
significant_genes_all_types_not_her2 = performStatisticalAnalysis(common_genes_not_her2, genes, .05 / length(genes))
print(significant_genes_all_types_not_her2)

significant_genes_not_significant_in_test_1_not_her2 = genes_not_in_original


filtered_normal_processed = filter(processed_test_normal, Gene %in% genes)
filtered_general_processed = filter(processed_general, Gene %in% genes)
significant_genes_test_3_not_her2_general =  compareTwoCategories(processed_test_normal, 
                                            normal_test_data, 
                                            filtered_general_processed,
                                            general_test_data, .01, "General_diff", .0) %>% rename(p_value_general = p_value)

"Test Confirm Data Set"
filtered_confirm_processed = filter(processed_confirm, Gene %in% genes)
significant_genes_test_3_not_her2_confirm =  compareTwoCategories(processed_test_normal, 
                                            normal_test_data, 
                                            filtered_confirm_processed,
                                            confirm_test_data, .01, "Confirm_diff", .0) %>% rename(p_value_confirm = p_value)
compared_to_large = inner_join(significant_genes_test_3_not_her2_general, significant_genes_test_3_not_her2_confirm)
print(compared_to_large)
write_tsv(compared_to_large, 'out/compated_to_large_not_her2.tsv')
"Each SUBTYPE VS LARGE NORMAL DATASET"

"TN"
significant_genes_TN_large_dataset =  compareTwoCategories(processed_test_normal, 
                                            normal_test_data, 
                                            processedTN,
                                            TN, .05, "TN_diff")
"HER2"
significant_genes_HER2_large_dataset =  compareTwoCategories(processed_test_normal, 
                                            normal_test_data, 
                                            processedHER2,
                                            HER2, .05, "HER2_diff")

"LuminalA"
significant_genes_LuminalA_large_dataset =  compareTwoCategories(processed_test_normal, 
                                            normal_test_data, 
                                            processedLuminalA,
                                            LuminalA, .05, "LuminalA_diff")

"LuminalB"
significant_genes_LuminalB_large_dataset =  compareTwoCategories(processed_test_normal, 
                                            normal_test_data, 
                                            processedLuminalB,
                                            LuminalB, .05, "LuminalB_diff")

"luminalHER2"
significant_genes_luminalHER2_large_dataset =  compareTwoCategories(processed_test_normal, 
                                            normal_test_data, 
                                            processedLuminalHER2,
                                            luminalHER2, .05, "luminalHER2_diff")                                            

"EACH SUBTYPE VS SUBTYPE DATASET"

"TN"
significant_genes_TN =  compareTwoCategories(processedNormal, 
                                            normalBreast, 
                                            processedTN,
                                            TN, .05, "TN_diff")
"HER2"
significant_genes_HER2 =  compareTwoCategories(processedNormal, 
                                            normalBreast, 
                                            processedHER2,
                                            HER2, .05, "HER2_diff")

"LuminalA"
significant_genes_LuminalA =  compareTwoCategories(processedNormal, 
                                            normalBreast, 
                                            processedLuminalA,
                                            LuminalA, .05, "LuminalA_diff")

"LuminalB"
significant_genes_LuminalB =  compareTwoCategories(processedNormal, 
                                            normalBreast, 
                                            processedLuminalB,
                                            LuminalB, .05, "LuminalB_diff")

"luminalHER2"
significant_genes_luminalHER2 =  compareTwoCategories(processedNormal, 
                                            normalBreast, 
                                            processedLuminalHER2,
                                            luminalHER2, .05, "luminalHER2_diff")   


"TEST 2 - GWAS"
gwas_genes = read_tsv('genes.tsv')
significant_gwas_genes_summary = inner_join(joined_final_results_confirmation_test, gwas_genes)
significant_gwas_genes = pull(significant_gwas_genes_summary, Gene)
print(significant_gwas_genes)

general_gwas_genes = filter(general_test_data, Gene %in% significant_gwas_genes)
normal_gwas_genes = filter(normal_test_data, Gene %in% significant_gwas_genes)
confirm_gwas_genes = filter(confirm_test_data, Gene %in% significant_gwas_genes)

"GWAS VS SUBTYPES"
"Prepare Data"
subtypes_gwas_genes = filter(joined_data, Gene %in% significant_gwas_genes)
normal_test_data = filter(normal_test_data, Gene %in% significant_gwas_genes)
HER2 = filter(subtypes_gwas_genes, Category == 'HER2')
LuminalA = filter(subtypes_gwas_genes, Category == 'Luminal A')
LuminalB = filter(subtypes_gwas_genes, Category == 'Luminal B')
luminalHER2 = filter(subtypes_gwas_genes, Category =='Luminal-HER2')
normalBreast = filter(subtypes_gwas_genes, Category == 'Normal')
TN = filter(subtypes_gwas_genes, Category =='TN')
inVitro = filter(subtypes_gwas_genes, Category == 'In Vitro')

"Process Data"
processedHER2GWAS = process_data(HER2, "HER2")
processedLuminalAGWAS = process_data(LuminalA, 'LUMINALA')
processedLuminalBGWAS = process_data(LuminalB, 'LUMINALB')
processedLuminalHER2GWAS = process_data(luminalHER2, 'LUMINALHER2')
processedNormalGWAS = process_data(normalBreast, 'NORMAL')
processedTNGWAS = process_data(TN, 'TN')
processedInVitroGWAS = process_data(inVitro, 'INVITRO')

processedAllSubTypeDataSetGWAS = list(processedHER2GWAS, processedLuminalAGWAS, 
                                processedLuminalBGWAS, processedLuminalHER2GWAS,
                                processedNormalGWAS, processedTNGWAS, processedInVitroGWAS) %>%
                                reduce(rbind)

"Compare all versus large normal"
HER2_NORMAL_GWAS = compareTwo(processedHER2, processedNormal)
TN_NORMAL_GWAS = compareTwo(processedTN, processedNormal)
LuminalA_NORMAL_GWAS = compareTwo(processedLuminalA, processedNormal)
LuminalB_NORMAL_GWAS = compareTwo(processedLuminalB, processedNormal)
LuminalHER2_NORMAL_GWAS = compareTwo(processedLuminalHER2, processedNormal)
InVitro_NORMAL_GWAS = compareTwo(processedInVitro, processedNormal)

#TODO OUT THESE TO A FILE

print(unique(pull(subtypes_gwas_genes, Category)))
print(processedAllSubTypeDataSetGWAS)
test_data = rbind(subtypes_gwas_genes, normal_test_data %>% mutate(Accession = "") %>% mutate(Category = "LargeNormal"))
print(test_data)
"STATISTICAL ANALYSIS"
print(paste0("ALPHA ", .05 / length(significant_gwas_genes)))
significant_genes_all_types_gwas = performStatisticalAnalysis(test_data, significant_gwas_genes, .05 / length(significant_gwas_genes))
print(significant_genes_all_types_gwas)
write_tsv(significant_genes_all_types_gwas, 'out/significant_genes_all_types_gwas.tsv')

"Each SUBTYPE VS LARGE NORMAL DATASET"
print(unique(pull(TN, Category)))
"TN"
significant_genes_TN_large_dataset =  compareTwoCategories(processed_test_normal, 
                                            normal_test_data, 
                                            processedTN,
                                            TN, .05, "TN_diff")
"HER2"
significant_genes_HER2_large_dataset =  compareTwoCategories(processed_test_normal, 
                                            normal_test_data, 
                                            processedHER2,
                                            HER2, .05, "HER2_diff")

"LuminalA"
significant_genes_LuminalA_large_dataset =  compareTwoCategories(processed_test_normal, 
                                            normal_test_data, 
                                            processedLuminalA,
                                            LuminalA, .05, "LuminalA_diff")

"LuminalB"
significant_genes_LuminalB_large_dataset =  compareTwoCategories(processed_test_normal, 
                                            normal_test_data, 
                                            processedLuminalB,
                                            LuminalB, .05, "LuminalB_diff")

"luminalHER2"
significant_genes_luminalHER2_large_dataset =  compareTwoCategories(processed_test_normal, 
                                            normal_test_data, 
                                            processedLuminalHER2,
                                            luminalHER2, .05, "luminalHER2_diff")                                            

"EACH SUBTYPE VS SUBTYPE DATASET"

"TN"
significant_genes_TN =  compareTwoCategories(processedNormal, 
                                            normalBreast, 
                                            processedTN,
                                            TN, .05, "TN_diff")
"HER2"
significant_genes_HER2 =  compareTwoCategories(processedNormal, 
                                            normalBreast, 
                                            processedHER2,
                                            HER2, .05, "HER2_diff")

"LuminalA"
significant_genes_LuminalA =  compareTwoCategories(processedNormal, 
                                            normalBreast, 
                                            processedLuminalA,
                                            LuminalA, .05, "LuminalA_diff")

"LuminalB"
significant_genes_LuminalB =  compareTwoCategories(processedNormal, 
                                            normalBreast, 
                                            processedLuminalB,
                                            LuminalB, .05, "LuminalB_diff")

"luminalHER2"
significant_genes_luminalHER2 =  compareTwoCategories(processedNormal, 
                                            normalBreast, 
                                            processedLuminalHER2,
                                            luminalHER2, .05, "luminalHER2_diff")   


"CREATE FILES FOR LATER ANALYSIS"

generate_one_v_all_diff_file = function(data1, filename, exclude) {
        data1_HER2 = compareTwo(data1, processedHER2)
        data1_LUM_A = compareTwo(data1, processedLuminalA)
        data1_LUM_B = compareTwo(data1, processedLuminalB)
        data1_TN = compareTwo(data1, processedTN)
        data1_LUM_HER2 = compareTwo(data1, processedLuminalHER2)
        data1_InVitro = compareTwo(data1, processedInVitro)
        data1_NORMAL = compareTwo(data1, processedNormal)

        data1_NORMAL = rename(data1_NORMAL, NORM_diff = diff)
        data1_TN = rename(data1_TN, TN_diff = diff)
        data1_LUM_A = rename(data1_LUM_A, LuminalA_diff = diff)
        data1_LUM_B = rename(data1_LUM_B, LuminalB_diff = diff)
        data1_LUM_HER2 = rename(data1_LUM_HER2, Luminal_HER2_diff = diff)
        data1_HER2 = rename(data1_HER2, HER2_diff = diff)
        data1_InVitro = rename(data1_InVitro, InVitro_diff = diff)
        data_joined = list(data1_NORMAL, data1_TN, data1_LUM_A, data1_LUM_B, data1_LUM_HER2, data1_HER2,data1_InVitro ) %>% reduce(full_join)

        write_tsv(
        data_joined,
        filename)
}

#Results from confirmation test
write_tsv(
  joined_final_results_confirmation_test,
  'out/results_confirmation_test.tsv'
)

write_tsv(
  significant_genes_all_types,
  'out/significant_genes_all_types.tsv'
)

#Save Processed Data
write_tsv(
  processedAllSubTypeDataSet,
  'out/processedAllSubTypeDataSet.tsv'
)

write_tsv(
  processedAllSubTypeDataSetGWAS,
  'out/processedGWASSubTypeDataSetVsLargeNormal.tsv'
)

#Save all common genes compared to normal
write_tsv(
  common_genes_all,
  'out/common_genes_all_vs_normal.tsv')


write_tsv(
  common_genes_not_her2, 
  'out/common_genes_not_her2.tsv'
)

write_tsv(
  significant_gwas_genes_summary,
  'out/significant_gwas_genes_summary.tsv'
)

write_tsv(
  significant_genes_all_types_not_her2, 
  'out/significant_genes_all_types_not_her2.tsv'
)
write_tsv(
  general_gwas_genes,
  'out/general_gwas_genes.tsv'
)
write_tsv(
  normal_gwas_genes,
  'out/normal_gwas_genes.tsv'
)
write_tsv(
  confirm_gwas_genes,
  'out/confirm_gwas_genes.tsv'
)
write_tsv(
  subtypes_gwas_genes,
  'out/subtypes_gwas_genes.tsv'
)
#HER2 vs all

generate_one_v_all_diff_file(processedHER2, "out/her2_vs_all.tsv", "HER2")

#LuminalA vs all

generate_one_v_all_diff_file(processedLuminalA, "out/lum_a_vs_all.tsv", "LUM_A")

#LuminalB vs all

generate_one_v_all_diff_file(processedLuminalB, "out/lum_b_vs_all.tsv", "LUM_B")

#Luminal HER2 vs all

generate_one_v_all_diff_file(processedLuminalHER2, "out/lum_her2_vs_all.tsv", "LUM_HER2")

#Normal vs all

generate_one_v_all_diff_file(processedNormal, "out/normal_vs_all.tsv", "NORMAL")

#TN vs all

generate_one_v_all_diff_file(processedTN, "out/TN_vs_all.tsv", "TN")

#Invitro vs all

generate_one_v_all_diff_file(processedInVitro, "out/invitro_vs_all.tsv", "INVITRO")


"EACH SUBTYPE SIGNIFICANT GENES"

write_tsv(
  significant_genes_TN,
  'out/significant_genes_TN.tsv'
)
write_tsv(
  significant_genes_HER2,
  'out/significant_genes_HER2.tsv'
)
write_tsv(
  significant_genes_LuminalA,
  'out/significant_genes_LuminalA.tsv'
)
write_tsv(
  significant_genes_LuminalB,
  'out/significant_genes_LuminalB.tsv'
)
write_tsv(
  significant_genes_luminalHER2,
  'out/significant_genes_luminalHER2.tsv'
)
                                  
"EACH SUBTYPE SIGNIFICANT GENES"

write_tsv(
  significant_genes_TN_large_dataset,
  'out/significant_genes_TN_large_dataset.tsv'
)
write_tsv(
  significant_genes_HER2_large_dataset,
  'out/significant_genes_HER2_large_dataset.tsv'
)
write_tsv(
  significant_genes_LuminalA_large_dataset,
  'out/significant_genes_LuminalA_large_dataset.tsv'
)
write_tsv(
  significant_genes_LuminalB_large_dataset,
  'out/significant_genes_LuminalB_large_dataset.tsv'
)
write_tsv(
  significant_genes_luminalHER2_large_dataset,
  'out/significant_genes_luminalHER2_large_dataset.tsv'
)

missing_genes = tibble(Gene=significant_genes_not_significant_in_test_1)
missing_genes_not_her_2 = tibble(Gene=significant_genes_not_significant_in_test_1_not_her2)

write_tsv(
  missing_genes,
  'out/common_genes_all_missing_from_test1.tsv'
)

write_tsv(
  missing_genes_not_her_2,
  'out/common_genes_all_minus_her2_missing_from_test1.tsv'
)


"All Gene no threshold"
normal_test_data = mutate(raw_normal, Category = 'Normal')
general_test_data = mutate(raw_general, Category = 'General')
confirm_test_data = mutate(raw_confirm, Category = 'Confirm')

"Compare First Dataset to Normal"
general_vs_normal_all = compareTwo(processed_test_normal, processed_general, 0)
general_vs_normal_all = rename(general_vs_normal_all, general_diff = diff)
significant_genes_test_1_all = compareTwoCategories(processed_test_normal, normal_test_data, processed_general, general_test_data, -1, "general_diff", 0)

#Confirm Results

"Compare Second Dataset to Normal"
confirm_vs_normal_all = compareTwo(processed_test_normal, processed_confirm, 0)
confirm_vs_normal_all = rename(confirm_vs_normal_all, confirm_diff = diff)
significant_genes_test_2_all =  compareTwoCategories(processed_test_normal, normal_test_data, processed_confirm, confirm_test_data, -1, "confirm_diff", 0)

"Combine Results"
significant_genes_test_1_all = rename(significant_genes_test_1_all, p_value_test_1 = p_value)
significant_genes_test_2_all = rename(significant_genes_test_2_all, p_value_test_2 = p_value)

joined_final_results_all = inner_join(significant_genes_test_1_all,significant_genes_test_2_all)
joined_all = inner_join(joined_final_results_all, general_vs_normal_all)
joined_final_results_confirmation_test_all = inner_join(joined_all, confirm_vs_normal_all)

write_tsv(joined_final_results_confirmation_test_all, 'out/joined_final_results_confirmation_test_all.tsv')

print(paste0(c("Number of Genes", length(unique(pull(joined_final_results_confirmation_test_all, Gene))))))
corrected_alpha = .05/length(unique(pull(joined_final_results_confirmation_test_all, Gene)))
print(paste0(c("Corrected alpha", corrected_alpha)))
