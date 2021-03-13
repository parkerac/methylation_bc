source('utils.r')


"TEST 1"

"Prepare Data"
raw_normal = read_tsv('GSE89278_RAW_Control_Methylation_Analysis_Data.txt')
raw_general = read_tsv('GSE72308_RAW_Methylation_Analysis_Data.txt')

normal_test_data = mutate(raw_normal, Category = 'Normal')
general_test_data = mutate(raw_general, Category = 'General')


processed_test_normal = process_data(raw_normal)
processed_general = process_data(raw_general)


"Compare First Dataset to Normal"
general_vs_normal = compareTwo(processed_test_normal, processed_general)
general_vs_normal = rename(general_vs_normal, general_diff = diff)
significant_genes_test_1 = compareTwoCategories(processed_test_normal, normal_test_data, processed_general, general_test_data, .01, "general_diff")

#Confirm Results
"Prepare Data"
raw_confirm = read_tsv('GSE84207_RAW_Methylation_Analysis_Data.txt')
confirm_test_data = mutate(raw_confirm, Category = 'Confirm')
processed_confirm = process_data(raw_confirm)

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

"TEST 2 - GWAS"
gwas_genes = read_tsv('genes.tsv')
significant_gwas_genes_summary = inner_join(joined_final_results_confirmation_test, gwas_genes)
significant_gwas_genes = pull(significant_gwas_genes_summary, Gene)
print(significant_gwas_genes)

general_gwas_genes = filter(general_test_data, Gene %in% significant_gwas_genes)
normal_gwas_genes = filter(normal_test_data, Gene %in% significant_gwas_genes)
confirm_gwas_genes = filter(confirm_test_data, Gene %in% significant_gwas_genes)

"TEST 3 - sub-types"

"Prepare Data"
classes = read_tsv('classifications.tsv')
raw_multiple_groups_data = read_tsv('GSE141338_RAW_Methylation_Analysis_Data.txt') %>% separate(Patient_ID, c("Accession", "Patient_ID"), "_", extra = "merge")
pre_name_joined_data = inner_join(classes, raw_multiple_groups_data)
joined_data = rename_categories(pre_name_joined_data)  %>% mutate(Category = as.factor(Category))

subtypes_gwas_genes = filter(joined_data, Gene %in% significant_gwas_genes)

HER2 = filter(joined_data, Category == 'HER2')
LuminalA = filter(joined_data, Category == 'Luminal A')
LuminalB = filter(joined_data, Category == 'Luminal B')
luminalHER2 = filter(joined_data, Category =='Luminal-HER2')
normalBreast = filter(joined_data, Category == 'Normal')
TN = filter(joined_data, Category =='TN')
inVitro = filter(joined_data, Category == 'In Vitro')

"Process Data"
processedHER2 = process_data(HER2)
processedLuminalA = process_data(LuminalA)
processedLuminalB = process_data(LuminalB)
processedLuminalHER2 = process_data(luminalHER2)
processedNormal = process_data(normalBreast)
processedTN = process_data(TN)
processedInVitro = process_data(inVitro)

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
filtered_normal_3 = filter(normal_test_data, Gene %in% genes) 
filtered_general_3 = filter(general_test_data, Gene %in% genes) 
joined_test_3 = bind_rows(filtered_normal_3, filtered_general_3) %>% mutate(Category = as.factor(Category))

corrected_alpha = .01/length(genes)
print(paste0("ALPHA ", corrected_alpha))
significant_genes_test_3 = performTwoWayStatisticalAnalysis(joined_test_3, genes, corrected_alpha)
print(significant_genes_test_3)

"Test Confirm Data Set"
filtered_confirm_3 = filter(confirm_test_data, Gene %in% genes) 
joined_test_4 = bind_rows(filtered_normal_3, filtered_confirm_3) %>% mutate(Category = as.factor(Category))

corrected_alpha = .01/length(genes)
print(paste0("ALPHA ", corrected_alpha))
significant_genes_test_4 = performTwoWayStatisticalAnalysis(joined_test_4, genes, corrected_alpha)
print(significant_genes_test_4)

general_vs_normal = compareTwo(processed_test_normal, processed_general, threshold=0)
print(filter(general_vs_normal, Gene %in% genes_not_in_original))
confirm_vs_normal = compareTwo(processed_test_normal, processed_confirm, threshold=0)
print(filter(confirm_vs_normal, Gene %in% genes_not_in_original))


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

genes_not_in_original = NULL
for (gene in genes){
  if (!(gene %in% significant_genes_final)){
    print(paste0(c(gene, " not found significant in original test")))
    genes_not_in_original = c(genes_not_in_original, gene)
  }
  
}

significant_genes_not_significant_in_test_1_not_her2 = genes_not_in_original

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
  'out/test_1_results/results_confirmation_test.tsv'
)

write_tsv(
  significant_genes_all_types,
  'out/test_3_results/significant_genes_all_types.tsv'
)

#Save Processed Data
write_tsv(
  processedHER2,
  'out/test_3_results/sumarized_data/HER2.tsv')

write_tsv(
  processedLuminalA,
  'out/test_3_results/sumarized_data/LuminalA.tsv')

write_tsv(
  processedLuminalB,
  'out/test_3_results/sumarized_data/LuminalB.tsv')

write_tsv(
  processedLuminalHER2,
  'out/test_3_results/sumarized_data/LuminalHER2.tsv')

write_tsv(
  processedNormal,
  'out/test_3_results/sumarized_data/Normal.tsv')

write_tsv(
  processedTN,
  'out/test_3_results/sumarized_data/TN.tsv')

write_tsv(
  processedInVitro,
  'out/test_3_results/sumarized_data/InVitro.tsv')

#Save all common genes compared to normal
write_tsv(
  common_genes_all,
  'out/test_3_results/common_genes_all_vs_normal.tsv')


write_tsv(
  common_genes_not_her2, 
  'out/test_3_results/common_genes_not_her2.tsv'
)

write_tsv(
  significant_gwas_genes_summary,
  'out/test_2_results/significant_gwas_genes_summary.tsv'
)

write_tsv(
  significant_genes_all_types_not_her2, 
  'out/test_3_results/significant_genes_all_types_not_her2.tsv'
)
write_tsv(
  general_gwas_genes,
  'out/test_2_results/general_gwas_genes.tsv'
)
write_tsv(
  normal_gwas_genes,
  'out/test_2_results/normal_gwas_genes.tsv'
)
write_tsv(
  confirm_gwas_genes,
  'out/test_2_results/confirm_gwas_genes.tsv'
)
write_tsv(
  subtypes_gwas_genes,
  'out/test_2_results/subtypes_gwas_genes.tsv'
)
#HER2 vs all

generate_one_v_all_diff_file(processedHER2, "out/test_3_results/other_out_files/her2_vs_all.tsv", "HER2")

#LuminalA vs all

generate_one_v_all_diff_file(processedLuminalA, "out/test_3_results/other_out_files/lum_a_vs_all.tsv", "LUM_A")

#LuminalB vs all

generate_one_v_all_diff_file(processedLuminalB, "out/test_3_results/other_out_files/lum_b_vs_all.tsv", "LUM_B")

#Luminal HER2 vs all

generate_one_v_all_diff_file(processedLuminalHER2, "out/test_3_results/other_out_files/lum_her2_vs_all.tsv", "LUM_HER2")

#Normal vs all

generate_one_v_all_diff_file(processedNormal, "out/test_3_results/other_out_files/normal_vs_all.tsv", "NORMAL")

#TN vs all

generate_one_v_all_diff_file(processedTN, "out/test_3_results/other_out_files/TN_vs_all.tsv", "TN")

#Invitro vs all

generate_one_v_all_diff_file(processedInVitro, "out/test_3_results/other_out_files/invitro_vs_all.tsv", "INVITRO")


"EACH SUBTYPE SIGNIFICANT GENES"

write_tsv(
  significant_genes_TN,
  'out/test_3_results/significant_genes_TN.tsv'
)
write_tsv(
  significant_genes_HER2,
  'out/test_3_results/significant_genes_HER2.tsv'
)
write_tsv(
  significant_genes_LuminalA,
  'out/test_3_results/significant_genes_LuminalA.tsv'
)
write_tsv(
  significant_genes_LuminalB,
  'out/test_3_results/significant_genes_LuminalB.tsv'
)
write_tsv(
  significant_genes_luminalHER2,
  'out/test_3_results/significant_genes_luminalHER2.tsv'
)
                                  
"EACH SUBTYPE SIGNIFICANT GENES"

write_tsv(
  significant_genes_TN_large_dataset,
  'out/test_3_results/significant_genes_TN_large_dataset.tsv'
)
write_tsv(
  significant_genes_HER2_large_dataset,
  'out/test_3_results/significant_genes_HER2_large_dataset.tsv'
)
write_tsv(
  significant_genes_LuminalA_large_dataset,
  'out/test_3_results/significant_genes_LuminalA_large_dataset.tsv'
)
write_tsv(
  significant_genes_LuminalB_large_dataset,
  'out/test_3_results/significant_genes_LuminalB_large_dataset.tsv'
)
write_tsv(
  significant_genes_luminalHER2_large_dataset,
  'out/test_3_results/significant_genes_luminalHER2_large_dataset.tsv'
)

missing_genes = tibble(Gene=significant_genes_not_significant_in_test_1)
missing_genes_not_her_2 = tibble(Gene=significant_genes_not_significant_in_test_1_not_her2)

write_tsv(
  missing_genes,
  'out/test_3_results/common_genes_all_missing_from_test1.tsv'
)

write_tsv(
  missing_genes_not_her_2,
  'out/test_3_results/common_genes_all_minus_her2_missing_from_test1.tsv'
)
