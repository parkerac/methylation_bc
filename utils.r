# install.packages('PMCMRplus')
library(PMCMR)
library('tidyverse')
#First four functions assume data in format Category, Acession, PatientID, Value, Gene
#Last one assume data in format from processData
rename_categories_general_cancer = function(data, normalStr = ' normal breast '){
        return(mutate(data, Category = 
        ifelse(str_detect(Category, normalStr), "Normal", "Cancer")))
}

rename_categories = function(data){
        return(mutate(data, Category = 
        ifelse(str_detect(Category, ' HER2 '), "HER2", 
        ifelse(str_detect(Category, ' luminal A '), "Luminal A", 
        ifelse(str_detect(Category, ' luminal B '), "Luminal B", 
        ifelse(str_detect(Category, ' luminal-HER2 '), "Luminal-HER2", 
        ifelse(str_detect(Category, ' normal breast '), "Normal", 
        ifelse(str_detect(Category, ' TN '), "TN", 
        ifelse(str_detect(Category, 'in vitro '), "In Vitro", "NA")))))))))
}


performTwoWayStatisticalAnalysis = function(data, genes, corrected_alpha=.01, process_data, colname){
        results = tibble(Gene = character(), p_value = numeric(), diff = numeric())
        for (gene in genes){
                stat_data = filter(data, Gene == gene) %>% mutate(Category = as.factor(Category))
                # view(stat_data)
                if (dim(stat_data)[1] == 0){
                }
                else {
                        p_value = wilcox.test(Value ~ Category, data=stat_data)$p.value
                        if (p_value < corrected_alpha){
                                diff = filter(process_data, Gene == gene)$diff
                                results = results %>% add_row(Gene = gene, p_value = p_value, diff = diff)
                        }
                }
        }
        
        results = rename(results, !!colname := diff)

        return(results)
}

compareTwoCategories = function(data1, rawData1, data2, rawData2, alpha, colname, threshold = .3){
        data1_vs_data2 = compareTwo(data1, data2, threshold)
        genes_test = pull(data1_vs_data2, Gene)

        filtered_data1 = filter(rawData1, Gene %in% genes_test) 
        filtered_data2 = filter(rawData2, Gene %in% genes_test) 
        joined_test = bind_rows(filtered_data1, filtered_data2) %>% mutate(Category = as.factor(Category))

        "Perform Statistical Analysis"
        if (alpha < 0){
                corrected_alpha = 100
        }
        else {
                corrected_alpha = alpha/length(genes_test)
        }
        print(paste0("ALPHA ", corrected_alpha))
        significant_genes_test = performTwoWayStatisticalAnalysis(joined_test, genes_test, corrected_alpha,data1_vs_data2, colname)
        print(significant_genes_test)
        return(significant_genes_test)
}

addToTable = function(tbl, data, gene, alpha = .01) {
        table = tbl
        for (name1 in rownames(data$p.value)) {
                for (name2 in colnames(data$p.value)){
                        p_value = data$p.value[name1, name2]
                        if(is.numeric(p_value) && !is.na(p_value) && p_value < alpha){
                                table = table %>% add_row(Gene = gene, 
                                                Category1 = name1, 
                                                Category2 = name2, 
                                                corrected_p_value = p_value)
                        }
                }
        }     
        return(table)
}

performStatisticalAnalysis = function(data, genes, alpha = .05){
        table = tibble(Gene = character(), 
              Category1 = character(), 
              Category2 = character(), 
              corrected_p_value = numeric())
        print(paste0("ALPHA ", alpha))
        for (gene in genes){
                print(paste("Statistical test for ", gene))
                result = statisticalAnalysis(data, gene)
                if (result == "NOT STATISTICALY DIFFERENT"){
                        print(result)
                        
                }
                else {
                        table = addToTable(table, result, gene, alpha)
                        print(result)
                }
        }
        return(table)
}

# Performs an ANOVA test if data is normal. If data is not the same, it them performs a TukeyHSD analysis
statisticalAnalysis = function(data, gene, alpha = .05){
        stat_data = filter(data, Gene == gene)

        if (kruskal.test(Value ~ Category, data=stat_data)$p.value < alpha) {
                results = posthoc.kruskal.dunn.test(Value ~ Category, data=stat_data, dist="Tukey", p.adjust.method="fdr")
                return(results)
        }
        else {
                return("NOT STATISTICALY DIFFERENT")
        }
}



process_data = function(data, category) {
        data %>% group_by(Gene) %>% 
                summarise(average=mean(Value), median=median(Value), 
                Category = category,
                min=min(Value), max=max(Value),
                Q1 = quantile(Value, .25, na.rm = TRUE), 
                Q3 = quantile(Value, .75, na.rm = TRUE),
                var = var(Value),
                IQR = IQR(Value, na.rm = TRUE)) %>% 
                arrange(desc(median)) -> processed_data
        return(processed_data)
}

compareTwo= function(data1, data2, threshold = .3){
        d1 = data1 %>% select(Gene, median) %>% rename(d1_med = median)
        d2 = data2 %>% select(Gene, median) %>% rename(d2_med = median)
        result = inner_join(d1, d2) %>% mutate(diff = d1_med - d2_med)
        return(arrange(filter(select(result, Gene, diff), abs(diff) > threshold), desc(diff)))
}


