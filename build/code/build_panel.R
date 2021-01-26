library(tidyverse)
library(fastDummies)

load("build/output/IRI_category_data.RData")

# ----------------------------------------------- #
# manipulate
# ----------------------------------------------- #

# weeks
minweek = as.numeric(min(data_store$WEEK))
maxweek = as.numeric(max(data_store$WEEK))

data_store %>%
  group_by(LARGE_CATEGORY) %>%
  summarise(SUBCAT = n_distinct(SMALL_CATEGORY),
            BRAND = n_distinct(BRAND),
            SALES = sum(DOLLARS),.groups="drop") %>%
  arrange(-SALES) %>%
  mutate(SHARE = SALES/sum(SALES),
         CUMSHARE = cumsum(SHARE))

# ----------------------------------------------- #
# 
# ----------------------------------------------- #

prod_def = c("LARGE_CATEGORY","SMALL_CATEGORY","BRAND")

# weeks in the training data
weeks = weeks_table %>%
  filter(YEAR==2011) %>%
  pull(WEEK)

# compute pack size by category
size_table = data_store %>%
  mutate(NUMBER = as.numeric(str_extract(word(UPC_DESCRIPTION,-1),"\\d\\.\\d+|\\d+")),
         UNITS = str_sub(UPC_DESCRIPTION,-2,-1),
         CAT_SIZE = NUMBER/VOL_EQ) %>%
  group_by(LARGE_CATEGORY,UNITS) %>%
  count(CAT_SIZE) %>%
  filter(n==max(n)) %>%
  mutate(CAT_SIZE = ifelse(UNITS %in% c("PT","LB"),CAT_SIZE*16,CAT_SIZE)) %>%
  ungroup() %>%
  select(-c(n,UNITS)) %>%
  distinct()

# compute category/subcategory shares and remove small subcategories
cat_weight_table = data_store %>%
  mutate(SALES = sum(DOLLARS)) %>%
  group_by(LARGE_CATEGORY) %>%
  mutate(CAT_SALES = sum(DOLLARS),
         CAT_SHARE = CAT_SALES/SALES,
         SUBCAT_COUNT = n_distinct(SMALL_CATEGORY)) %>%
  group_by(LARGE_CATEGORY,SMALL_CATEGORY) %>%
  summarise(CAT_SALES = CAT_SALES[1],
            CAT_SHARE = CAT_SHARE[1],
            SUBCAT_COUNT = SUBCAT_COUNT[1],
            CAT_N = max(round(100*CAT_SHARE),SUBCAT_COUNT),
            SUBCAT_SALES = sum(DOLLARS),
            SUBCAT_SHARE = SUBCAT_SALES / CAT_SALES,
            N = max(2,round(CAT_N*SUBCAT_SHARE)), .groups="drop") %>%
  filter(SUBCAT_SHARE>0.01) %>%
  select(LARGE_CATEGORY,SMALL_CATEGORY,N)
small_categories = cat_weight_table$SMALL_CATEGORY

# ----------------------------------------------- #
# select products
# ----------------------------------------------- #

# product summary
product_table = data_store %>%
  filter(SMALL_CATEGORY %in% small_categories) %>%
  filter(UNITS>0) %>%
  group_by_at(c("IRI_KEY",prod_def)) %>%
  mutate(PRICE = DOLLARS/UNITS,
         LAG_PRICE = lag(PRICE,1)) %>%
  group_by_at(prod_def) %>%
  summarise(revenue = sum(DOLLARS,na.rm=TRUE),
            n_stores = n_distinct(IRI_KEY),
            n_weeks = n_distinct(WEEK),
            price_levels = n_distinct(round(PRICE,2)),
            price_change = mean(abs(PRICE-LAG_PRICE)>0.01,na.rm=TRUE),
            UPCs = list(UPC), .groups="drop") %>%
  # remove missing
  filter(n_weeks == length(weeks)) %>%
  filter(price_levels > 1) %>%
  # filter(price_change > 0.05) %>%
  filter(n_stores==length(stores)) %>%
  # ranking
  left_join(cat_weight_table,by=c("LARGE_CATEGORY","SMALL_CATEGORY")) %>%
  group_by(LARGE_CATEGORY,SMALL_CATEGORY) %>%
  mutate(RANK=rank(-revenue)) %>%
  filter(RANK<=3*N) %>%
  # create product IDs
  group_by_at(prod_def) %>%
  mutate(PRODUCT = cur_group_id()) %>%
  ungroup() %>%
  select(PRODUCT,LARGE_CATEGORY,SMALL_CATEGORY,BRAND,
         revenue,n_weeks,price_change,price_levels,UPCs)

# share of revenue
upcs = unique(unlist(product_table$UPCs))
coverage_table = data_store %>%
  filter(WEEK %in% weeks) %>%
  filter(IRI_KEY %in% stores) %>%
  mutate(TOTAL=sum(DOLLARS)) %>%
  filter(UPC %in% upcs) %>%
  summarise(SHARE=sum(DOLLARS)/TOTAL[1])

# trees
tree_names = product_table %>%
  select(LARGE_CATEGORY,SMALL_CATEGORY,BRAND) %>%
  distinct()
tree = product_table %>%
  select(LARGE_CATEGORY,SMALL_CATEGORY,BRAND) %>%
  mutate(LARGE_CATEGORY = as.numeric(factor(LARGE_CATEGORY,levels=unique(LARGE_CATEGORY))),
         SMALL_CATEGORY = as.numeric(factor(SMALL_CATEGORY,levels=unique(SMALL_CATEGORY))),
         BRAND = as.numeric(factor(paste(SMALL_CATEGORY,BRAND),levels=unique(paste(SMALL_CATEGORY,BRAND))))) %>%
  as.matrix()

# ----------------------------------------------- #
# construct panel
# ----------------------------------------------- #

# compute UPC weight table
weight_table = data_store %>%
  filter(UPC %in% upcs) %>%
  group_by(LARGE_CATEGORY,SMALL_CATEGORY,BRAND) %>%
  mutate(TOTAL_SALES=sum(DOLLARS)) %>%
  group_by(LARGE_CATEGORY,SMALL_CATEGORY,BRAND,UPC) %>%
  summarise(WT=sum(DOLLARS)/TOTAL_SALES[1], .groups="drop")

# compute panel in long format
panel = data_store %>%
  filter(IRI_KEY %in% stores, UPC %in% upcs) %>%
  filter(UNITS>0) %>%
  left_join(size_table,by="LARGE_CATEGORY") %>%
  left_join(weight_table,by=c("LARGE_CATEGORY","SMALL_CATEGORY","BRAND","UPC")) %>%
  mutate(PRICE = DOLLARS/UNITS,
         FEATURE = ifelse(F=="NONE",0,1),
         DISPLAY = ifelse(D==0,0,1)) %>%
  # create rows for missing UPC prices
  select(IRI_KEY,WEEK,UPC,UNITS,PRICE,FEATURE,DISPLAY,UPC_DESCRIPTION,BRAND,VENDOR,SMALL_CATEGORY,
         LARGE_CATEGORY,VOL_EQ,CAT_SIZE,WT) %>%
  group_by(across(-c(WEEK,UNITS,PRICE,FEATURE,DISPLAY))) %>%
  complete(WEEK=seq(minweek,maxweek)) %>%
  # impute missing values
  fill(PRICE,.direction="downup") %>%
  mutate(across(everything(),~replace_na(.,0))) %>%
  group_by_at(c("IRI_KEY","WEEK",prod_def)) %>%
  # normalize weights
  mutate(WT = WT/sum(WT)) %>%
  # compute weighted average of marketing variables
  summarise(PRICE = sum(WT*PRICE/VOL_EQ),
            UNITS = sum(UNITS*VOL_EQ),
            FEATURE = sum(WT*FEATURE),
            DISPLAY = sum(WT*DISPLAY),
            UPCs = list(UPC), .groups="drop")
  # summarise(PRICE = sum(WT*PRICE/(VOL_EQ*CAT_SIZE)),
  #           UNITS = sum(UNITS*VOL_EQ*CAT_SIZE),
  #           FEATURE = sum(WT*FEATURE),
  #           DISPLAY = sum(WT*DISPLAY),
  #           UPCs = list(UPC), .groups="drop")

# ----------------------------------------------- #
# construct price, demand, and promotion tables
# ----------------------------------------------- #

# prices, quantities, promotions
prices = panel %>%
  inner_join(product_table[,c("PRODUCT",prod_def)], by = c("LARGE_CATEGORY", "SMALL_CATEGORY", "BRAND")) %>%
  mutate(PRICE=log(PRICE)) %>%
  select(IRI_KEY,PRODUCT,WEEK,PRICE) %>%
  pivot_wider(names_from=PRODUCT,values_from=PRICE,names_prefix="PRICE_")
# check UNITS ###################################
demand = panel %>%
  inner_join(product_table[,c("PRODUCT",prod_def)], by = c("LARGE_CATEGORY", "SMALL_CATEGORY", "BRAND")) %>%
  mutate(UNITS = ifelse(UNITS==0,log(0.1),log(UNITS))) %>%
  select(IRI_KEY,PRODUCT,WEEK,UNITS) %>%
  arrange(IRI_KEY,PRODUCT,WEEK) %>%
  pivot_wider(names_from=PRODUCT,values_from=UNITS,names_prefix="UNITS_")
feature = panel %>%
  inner_join(product_table[,c("PRODUCT",prod_def)], by = c("LARGE_CATEGORY", "SMALL_CATEGORY", "BRAND")) %>%
  mutate(IRI_KEY,FEATURE = ifelse(is.nan(FEATURE),0,FEATURE)) %>%
  group_by(PRODUCT) %>%
  select(IRI_KEY,PRODUCT,WEEK,FEATURE) %>%
  pivot_wider(names_from=PRODUCT,values_from=FEATURE,names_prefix="FEATURE_")
display = panel %>%
  inner_join(product_table[,c("PRODUCT",prod_def)], by = c("LARGE_CATEGORY", "SMALL_CATEGORY", "BRAND")) %>%
  mutate(DISPLAY = ifelse(is.nan(DISPLAY),0,DISPLAY)) %>%
  group_by(PRODUCT) %>%
  select(IRI_KEY,PRODUCT,WEEK,DISPLAY) %>%
  pivot_wider(names_from=PRODUCT,values_from=DISPLAY,names_prefix="DISPLAY_") 
seasonal = panel %>%
  distinct(WEEK) %>%
  left_join(weeks_table,by="WEEK") %>%
  dummy_cols(select_columns="QUARTER",
             remove_first_dummy = TRUE) %>%
  left_join(demand,by="WEEK") %>%
  select(starts_with("QUARTER_")) 

# save 
save(panel,demand,prices,feature,display,seasonal,
     product_table,tree_names,tree,coverage_table,markets,
     file="build/output/store_panel.RData")
