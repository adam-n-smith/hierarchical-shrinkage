library(lubridate)
library(scales)
library(broom)

# prices
panel %>%
  left_join(weeks_table, by="WEEK") %>%
  ggplot(aes(x=WEEK_END,y=PRICE,color=interaction(SMALL_CATEGORY,BRAND))) + 
  geom_line(show.legend=FALSE) + 
  labs(x="") +
  scale_x_date(date_breaks = "4 months", labels = date_format("%m-%Y")) +
  facet_wrap(LARGE_CATEGORY~., scales="free") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90))

# quantities
panel %>%
  left_join(weeks_table, by="WEEK") %>%
  ggplot(aes(x=WEEK_END,y=log(UNITS),color=interaction(SMALL_CATEGORY,BRAND))) + 
  geom_line(show.legend=FALSE) + 
  labs(x="") +
  scale_x_date(date_breaks = "4 months", labels = date_format("%m-%Y")) +
  facet_wrap(LARGE_CATEGORY~., scales="free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90))

# price vs. demand
panel %>%
  ggplot(aes(x=log(PRICE), y=log(UNITS),color=SMALL_CATEGORY)) + 
  geom_point(show.legend=FALSE) + 
  facet_wrap(LARGE_CATEGORY~., scales="free") +
  theme_minimal() 
panel %>%
  ggplot(aes(x=log(PRICE), y=log(UNITS),color=interaction(SMALL_CATEGORY,BRAND))) + 
  geom_point(show.legend=FALSE) + 
  facet_wrap(LARGE_CATEGORY~., scales="free") + 
  theme_minimal() + 
  theme(strip.text.x = element_text(face="bold"))

# regressions
reg = panel %>%
  nest_by(LARGE_CATEGORY,SMALL_CATEGORY,BRAND) %>%
  mutate(fit = list(lm(log(UNITS)~log(PRICE)+FEATURE+DISPLAY,data=data))) %>%
  summarise(tidy(fit)) 
reg %>%
  filter(str_detect(term,"PRICE")) %>%
  ggplot(aes(x=estimate,fill=p.value<0.05)) + 
  geom_histogram() +
  xlim(-10,10)
