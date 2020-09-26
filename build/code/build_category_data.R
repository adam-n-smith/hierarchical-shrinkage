library(readxl)
library(tidyverse)
library(lubridate)
library(fastDummies)
library(janitor)

hd_path = "/Volumes/LaCie/IRI Academic Data/"
categories = c("beer","carbbev","coffee","coldcer","fzdinent",
               "fzpizza","hotdog","margbutr","mayo","milk",
               "mustketc","peanbutr","saltsnck","spagsauc",
               "soup","yogurt")

years = 2011

# weeks data
weeks_table = read_xls(paste0(hd_path,"Year11/beer/IRI week translation_2008_2017.xls"),
                 col_types=c("numeric","skip","date","skip","skip","skip")) %>%
  rename("WEEK"=`IRI Week`,"WEEK_END"=`Calendar week ending on`) %>%
  mutate(WEEK_END=as.Date(WEEK_END),
         MONTH = month(WEEK_END),
         QUARTER = quarter(WEEK_END),
         YEAR = year(WEEK_END))

# ----------------------------------------------- #
# load UPC tables
# ----------------------------------------------- #

upc_table = NULL
for(category in categories){
  
  for(year in years){
    
    # short year name
    yr = as.numeric(substr(year,3,4))
    
    # check to see if there is UPC file for given year
    base_files = list.files(hd_path)
    folder = base_files[str_detect(base_files,"parsed stub files") & str_detect(base_files,as.character(yr))]
    
    if(length(folder)>0){
      
      # upc table
      path = paste0(hd_path,folder)
      files = list.files(path)
      wchfile = files[str_detect(files,category)]
      cat_upc_table = read_xlsx(paste0(path,"/",wchfile), col_types="guess") %>%
        rename("LARGE_CATEGORY"=L1,"SMALL_CATEGORY"=L2,"PARENT_COMPANY"=L3,
               "VENDOR"=L4,"BRAND"=L5,"UPC_DESCRIPTION"=L9) %>%
        select("SY","GE","VEND","ITEM","LARGE_CATEGORY","SMALL_CATEGORY","PARENT_COMPANY",
               "VENDOR","BRAND","UPC","UPC_DESCRIPTION","VOL_EQ")
      
      # bind together
      upc_table = bind_rows(upc_table,cat_upc_table)
      rm(cat_upc_table)
      
    }
    
  }
  
}

upc_table = upc_table %>% 
  # remove duplicates
  group_by(SY,GE,VEND,ITEM) %>%
  mutate(n=1:n()) %>%
  filter(n==n()) %>%
  select(-n) %>%
  # rename large category column
  mutate(LARGE_CATEGORY = str_remove(LARGE_CATEGORY,"CATEGORY - "))

# ----------------------------------------------- #
# load store tables
# ----------------------------------------------- #

store_table = read_table(paste0(hd_path,"Year",yr,"/",categories[1],"/Delivery_Stores"))

markets = "PITTSFIELD"

chain = store_table %>%
  filter(Market_Name %in% markets) %>%
  top_n(1,EST_ACV) %>%
  pull(MskdName)

stores = store_table %>%
  filter(Market_Name %in% markets, MskdName %in% chain) %>%
  pull(IRI_KEY)

# ----------------------------------------------- #
# load transaction data
# ----------------------------------------------- #

for(category in categories){

  cat_data_store = NULL
  for(year in years){
    
    # short year name
    yr = as.numeric(substr(year,3,4))

    # file path for transaction data
    txn_path = paste0(hd_path,"Year",yr,"/",category)
    
    # transaction data
    filename = list.files(txn_path)[str_detect(list.files(txn_path),"groc")]
    cat_data_store_tmp = read_table2(paste0(txn_path,"/",filename), col_names=TRUE, 
                                 cols(SY=col_character(),
                                      GE=col_character(),
                                      VEND=col_character(),
                                      ITEM=col_character())) %>%
      filter(IRI_KEY %in% stores) %>%
      left_join(upc_table,by=c("SY","GE","VEND","ITEM")) %>%
      left_join(store_table,by="IRI_KEY")
      
    cat_data_store = bind_rows(cat_data_store,cat_data_store_tmp)
    
    rm(cat_data_store_tmp)
    
  }
  
  # save category file
  assign(paste0("data_store_",category),cat_data_store)
  rm(cat_data_store)
  print(category)
}

# bind together all data
data_store = NULL
for(category in categories){
  data_store = rbind(data_store,get(paste0("data_store_",category)))
}

# save
save(data_store,upc_table,years,store_table,markets,stores,weeks_table,
     file="build/output/IRI_category_data.RData")
