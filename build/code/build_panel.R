library(tidyverse)
library(fastDummies)

load("build/output/IRI_category_data.RData")

# ----------------------------------------------- #
# preliminaries
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

prod_def = c("LARGE_CATEGORY","SMALL_CATEGORY","BRAND")

# weeks in the training data
weeks = weeks_table %>%
  filter(YEAR %in% years) %>%
  pull(WEEK)

# compute pack size by category
size_table = data_store %>%
  mutate(NUMBER = as.numeric(str_extract(word(UPC_DESCRIPTION,-1),"\\d\\.\\d+|\\d+")),
         UOM = str_sub(UPC_DESCRIPTION,-2,-1),
         CAT_SIZE = NUMBER/VOL_EQ) %>%
  drop_na() %>%
  group_by(LARGE_CATEGORY,UOM) %>%
  count(CAT_SIZE) %>%
  filter(n==max(n)) %>%
  mutate(CAT_SIZE = ifelse(UOM %in% c("PT","LB"),CAT_SIZE*16,CAT_SIZE)) %>%
  ungroup() %>%
  select(-c(n,UOM)) %>%
  distinct()

# ----------------------------------------------- #
# select products
# ----------------------------------------------- #

# small categories
small_categories = data_store %>%
  left_join(weeks_table) %>%
  group_by(LARGE_CATEGORY,SMALL_CATEGORY) %>%
  summarise(REVENUE = sum(DOLLARS, na.rm=TRUE),
            nYEAR = n_distinct(YEAR),
            nBRAND = n_distinct(BRAND)) %>%
  group_by(LARGE_CATEGORY) %>%
  mutate(SHARE = REVENUE/sum(REVENUE)) %>%
  filter(SHARE>0.02, nYEAR==length(years)) %>%
  pull(SMALL_CATEGORY)

# product
product_table = data_store %>%
  filter(SMALL_CATEGORY %in% small_categories) %>%
  filter(UNITS>0) %>%
  group_by_at(c("IRI_KEY",prod_def)) %>%
  mutate(PRICE = DOLLARS/UNITS) %>%
  group_by_at(prod_def) %>%
  summarise(revenue = sum(DOLLARS,na.rm=TRUE),
            n_stores = n_distinct(IRI_KEY),
            min_week = min(WEEK),
            max_week = max(WEEK),
            n_weeks = n_distinct(IRI_KEY,WEEK),
            price_levels = n_distinct(round(PRICE,2)),
            UPCs = list(UPC),
            nUPCs = n_distinct(UPC), .groups="drop") %>%
  group_by_at(c("LARGE_CATEGORY","SMALL_CATEGORY")) %>%
  mutate(total_revenue = sum(revenue),
         share = revenue/total_revenue) %>%
  arrange(LARGE_CATEGORY,SMALL_CATEGORY,-share) %>%
  mutate(cum_share = cumsum(share),
         rank = 1:n()) %>%
  # apply filtering criteria
  filter(cum_share<0.85 | rank<=2) %>%
  filter(n_weeks>0.95*length(stores)*length(weeks)) %>%
  # create product IDs
  group_by_at(prod_def) %>%
  mutate(PRODUCT = cur_group_id()) %>%
  ungroup() %>%
  select(PRODUCT,LARGE_CATEGORY,SMALL_CATEGORY,BRAND,nUPCs,UPCs,
         revenue,min_week,max_week,n_weeks,price_levels,share,cum_share)

# share of revenue
upcs = unique(unlist(product_table$UPCs))
coverage = data_store %>%
  filter(WEEK %in% weeks) %>%
  filter(IRI_KEY %in% stores) %>%
  mutate(TOTAL=sum(DOLLARS)) %>%
  filter(UPC %in% upcs) %>%
  summarise(SHARE=sum(DOLLARS)/TOTAL[1]) %>%
  pull(SHARE)

# trees
tree_names = product_table %>%
  select(LARGE_CATEGORY,SMALL_CATEGORY,BRAND) %>%
  distinct() %>%
  mutate(PRODUCT = 1:n())
tree_names_clean = tree_names %>%
  mutate(SMALL_CATEGORY = to_any_case(SMALL_CATEGORY,case="title",parsing_option=0),
         SMALL_CATEGORY = str_remove(SMALL_CATEGORY,"\\([^(]*$"),
         BRAND = to_any_case(BRAND,case="title",parsing_option=0))
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
  left_join(weeks_table) %>%
  filter(UPC %in% upcs) %>%
  group_by(LARGE_CATEGORY,SMALL_CATEGORY,BRAND,YEAR) %>%
  mutate(TOTAL_SALES=sum(DOLLARS)) %>%
  group_by(LARGE_CATEGORY,SMALL_CATEGORY,BRAND,UPC,YEAR) %>%
  summarise(WT=sum(DOLLARS)/TOTAL_SALES[1], .groups="drop")

# compute panel in long format
panel = data_store %>%
  filter(IRI_KEY %in% stores, UPC %in% upcs) %>%
  filter(UNITS>0) %>%
  left_join(size_table,by="LARGE_CATEGORY") %>%
  mutate(PRICE = DOLLARS/UNITS,
         FEATURE = ifelse(F=="NONE",0,1),
         DISPLAY = ifelse(D==0,0,1),
         SIZE = VOL_EQ*CAT_SIZE) %>%
  # bring in and normalize weights
  left_join(weeks_table,by="WEEK") %>%
  left_join(weight_table,by=c("LARGE_CATEGORY","SMALL_CATEGORY","BRAND","UPC","YEAR")) %>%
  mutate(WT = replace_na(WT,0)) %>%
  group_by_at(c("IRI_KEY","WEEK",prod_def)) %>%
  mutate(WT = WT/sum(WT)) %>%
  # compute weighted average of marketing variables
  # summarise(PRICE = sum(WT*PRICE/SIZE),
  #           UNITS = sum(UNITS*SIZE),
  #           FEATURE = sum(WT*FEATURE),
  #           DISPLAY = sum(WT*DISPLAY)) %>%
  summarise(PRICE = sum(WT*PRICE/VOL_EQ),
            UNITS = sum(UNITS*VOL_EQ),
            FEATURE = sum(WT*FEATURE),
            DISPLAY = sum(WT*DISPLAY)) %>%
  # imputation for missing obs
  group_by_at(c("IRI_KEY",prod_def)) %>%
  complete(WEEK=seq(minweek,maxweek)) %>%
  group_by_at(c("IRI_KEY",prod_def)) %>%
  fill(PRICE,.direction="downup") %>%
  mutate(UNITS = replace_na(UNITS,min(UNITS[!is.na(UNITS)])/2),
         FEATURE = replace_na(FEATURE,0),
         DISPLAY = replace_na(DISPLAY,0))%>%
  ungroup()

# ----------------------------------------------- #
# construct price, demand, and promotion tables
# ----------------------------------------------- #

# prices, quantities, promotions
prices = panel %>%
  inner_join(product_table[,c("PRODUCT",prod_def)], by = prod_def) %>%
  mutate(PRICE = log(PRICE)) %>%
  select(IRI_KEY,PRODUCT,WEEK,PRICE) %>%
  pivot_wider(names_from=PRODUCT,values_from=PRICE,names_prefix="PRICE_")
demand = panel %>%
  inner_join(product_table[,c("PRODUCT",prod_def)], by = prod_def) %>%
  mutate(UNITS = log(UNITS)) %>%
  select(IRI_KEY,PRODUCT,WEEK,UNITS) %>%
  arrange(IRI_KEY,PRODUCT,WEEK) %>%
  pivot_wider(names_from=PRODUCT,values_from=UNITS,names_prefix="UNITS_")
feature = panel %>%
  inner_join(product_table[,c("PRODUCT",prod_def)], by = prod_def) %>%
  group_by(PRODUCT) %>%
  select(IRI_KEY,PRODUCT,WEEK,FEATURE) %>%
  pivot_wider(names_from=PRODUCT,values_from=FEATURE,names_prefix="FEATURE_")
display = panel %>%
  inner_join(product_table[,c("PRODUCT",prod_def)], by = prod_def) %>%
  group_by(PRODUCT) %>%
  select(IRI_KEY,PRODUCT,WEEK,DISPLAY) %>%
  pivot_wider(names_from=PRODUCT,values_from=DISPLAY,names_prefix="DISPLAY_") 
seasonal = panel %>%
  distinct(WEEK) %>%
  left_join(weeks_table,by="WEEK") %>%
  dummy_cols(select_columns="QUARTER",
             remove_first_dummy = TRUE) %>%
  left_join(demand,by="WEEK") %>%
  select(starts_with(c("MONTH","QUARTER_","YEAR","HOLIDAY"))) %>%
  mutate(SUMMER = ifelse(MONTH>=6 & MONTH<=9,1,0))

# save 
save(panel,demand,prices,feature,display,seasonal,
     product_table,tree,tree_names,tree_names_clean,coverage,markets,
     file="build/output/store_panel.RData")
