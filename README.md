---
title: "Liver Transplant Survival Analysis"
author: "Attoumassa Samaké, Vixra Keo, Mounir M'barki, Antoine Bedouch"
output:
  pdf_document: default
  html_notebook: default
---


```{r, echo=FALSE, out.width="75%", include=FALSE}
# Importing the packages if not installed yet
# install.packages("survival", "survminer", "dplyr", "lubridate", "ggsurvfit", "gtsummary", "tidycmprsk", "finalfit")
library(survival)
library(survminer)
library(dplyr)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(finalfit)
```


# Introduction : 

For our survival analysis work, we chose a dataset available in the survival R package call "transplant". This dataset provides information about subjects from the Mayo Clinic Rochester liver transplant waiting list from 1990 and 1999. It contains 815 patients, from which 636 underwent a transplantation, 66 died while waiting, 37 withdrew from the waiting list and 76 were censored. This data initially comes from a scientific study (Kim et al., Deaths on the liver transplant waiting list: An analysis of competing risks, 2006), though the final dataset used by the study is more complete (it includes data on the patients disease)

The dateset contains 6 variables :

  * age : age at the addition of the waiting list
  * sex : male or female
  * abo : which is the blood type (A, B, AB and 0)
  * year : the year the patient entered in the waiting list
  * futime : the time (in days) between the registry on the waiting list and the final outcome
  * event : 
    + ltx = The patient was transplanted during the study
    + death = The patient died while on the waiting list
    + withdraw = The patient withdrew from the waiting list
    + censored = The outcome is censored (for most cases, no outcome at the end of the study in 1999)

This dataset reports on the evolution of liver transplants. It aims to answer crucial questions about wait times, the demographic composition of people on the waiting list, and how the patient medical information are impacting the death and transplantation outcome.

An essential factor influencing transplant dynamics is blood type. The dataset highlights the impact of blood type compatibility on the process.Donor livers from type O individuals can be utilized by recipients with A, B, AB, or O blood types, while an AB liver is exclusive to an AB recipient.

Our goal will be to study the survival of two outcomes:
  * the transplant outcome `ltx`
  * the death while waiting `death`

We will consider the patients who withdrew from the study as right-censored data. 


# Data preparation 


First let's see what our dataset looks like: 

```{r, echo=FALSE, out.width="75%"}
data = transplant
summary(data)
```

We note that 18 values for the age are missing. We could just ignore these patients in our analysis. Instead we have decided to impute these missing values with the average age of 50 years old.

For the age non parametric analysis, we split the patients in two age groups:

  * old (age >= 50)
  * young (age < 50)


Finally, we created in our dataframe specific boolean and integer columns containing outcome data for future analysis

```{r, echo=FALSE, out.width="75%"}
age_moyen <- mean(data$age, na.rm = TRUE)
data$age[is.na(data$age)] <- as.integer(age_moyen)
print(paste("average age =", as.integer(age_moyen)))

data <- data %>% mutate(age_group = ifelse(age >=50, "old", "young"))
data$age_group <- factor(data$age_group)

data$death <- (data$event == "death")
data$ltx <- (data$event == "ltx")
data <- data %>% mutate(
  status = as.factor(recode(event, 'death'=2, 'ltx' = 1, 'censored' = 0, 'withdraw' = 0)),
  status_ltx = as.factor(recode(event, 'death'=2, 'ltx' = 1, 'censored' = 0, 'withdraw' = 0)),
  status_death = as.factor(recode(event, 'death'=1, 'ltx' = 2, 'censored' = 0, 'withdraw' = 0))
)
```


# Non parametric modeling

Before doing any modeling, we can look at the empirical cumulative incidence function to see the relative incidence of each outcome: transplant or death. This curves tells us what is the probability to have a given outcome after a given number of days.

```{r, echo=FALSE, out.width="75%"}
cuminc_obj <- tidycmprsk::cuminc(Surv(futime, status) ~ 1, data=data)

ggcuminc(cuminc_obj, outcome = c(1, 2)) +
  labs(
    x = "Days",
    y = "Cumulative incidence"
  ) +
  add_confidence_interval()

print(paste("Median number of days before the patient is transplanted:",median(data[data$ltx,]$futime)))
print(paste("Median number of days before the patient dies:",median(data[data$death,]$futime)))
```
Legend:
The outcome 1 is the transplant incidence
The outcome 2 is the death incidence

We observe that in our dataset the majority outcome is to be transplanted (~78%) while the probability to die while on the waiting list is much lower (~8%). We also observe that most outcomes happen quite quickly after the registry on the list (108.5 days if the outcome is transplantation, 65.5 days if the outcome is death).



## Kaplan-Meier Method and Log Rank Test

The simplest survival modeling is the Kaplan-Meier model. It is  non-parametric method used to estimate the survival function from observed survival times in a sample. It is a good starting point to get an idea of the survival rates of our dataset. 


Note that the Kaplan-Meier analysis only considers one outcome at a time. Therefore when using Kaplan-Meier on a given outcome (say death), it considers that no alternate outcome existed. In other words, it considers that other outcomes are censored (say transplant) and therefore will overestimate the hazard of the studied outcome. If we study the death outcome, then it will consider that the patients that received a transplant have merely been censored from the study, and therefore have an equal chance of dying as other patients who haven't been transplanted. We will talk more about this in the last section.

We will plot the survival probability of the transplant and death case considering the different covariates : sex, age_group and blood type. 
The survival probability represents the likelihood that an event of interest has not occurred up to a certain point in time.

### Kaplan-Meier modeling of transplantation hasard:

```{r, echo=FALSE, out.width="75%"}
#Examine the predictive value of the patient's sex
data$fustat <- ifelse(data$event == "ltx", 1, 0)
#Fit survival data using the Kaplan-Meier method :
surv_object <- Surv(time = data$futime, event = data$fustat)
fit1 <- survfit(surv_object ~ sex, data )
ggsurvplot(fit1, data = data, pval = TRUE, title = "Transplant probability by gender ")
```


This plot compares the 2 groups male and female considering the outcome transplantation. 
The p value of 0.6 which is greater than 0.05 indicates that the difference between the 2 groups is not statistically significant. 
We can then conclude that the gender has no impact for the transplantation outcome. 


```{r, echo=FALSE, out.width="75%"}
#Examine the predictive value of the patient's age_group
data$fustat <- ifelse(data$event == "ltx", 1, 0)
surv_object <- Surv(time = transplant$futime, event = data$fustat)
fit2 <- survfit(surv_object ~ age_group,data )
ggsurvplot(fit2, data = data, pval = TRUE, title="Transplant probability by age")
```


We now consider the transplantation outcome by age. As described previously, the sample has been divided into 2 groups : patients younger than 50 years old and older than 50.
The p-value of 0.38 means that there is no significant difference between the 2 groups for transplanation.

```{r, echo=FALSE, out.width="75%"}
#Examine the predictive value of the patient's blood type
data$fustat <- ifelse(data$event == "ltx", 1, 0)
surv_object <- Surv(time = data$futime, event = data$fustat)
fit3 <- survfit(surv_object ~ abo,data)
ggsurvplot(fit3, data = data, pval = TRUE, title="Transplant probability by blood type")
```


The p-value which is 0.001, indicates a significant result if considering that p < 0.05 denotes statistical significance. In this study, blood groups have a significant impact on the patient's transplant waiting time. Overall, the curves indicate that patients with blood groups AB and A have a shorter waiting time compared to others, and patients with type O wait longer. 
This is explained by the fact that livers from donors with blood group O can be used by patients with blood groups A, B, AB, or O. This is what we observe when looking at the curves: the O survival is always above the other blood types. This means that patients with O blood type wait the longest for a transplant.


However, one may wonder if the waiting time has an impact on the patient's survival: let's analyse the death outcome.


```{r, echo=FALSE, out.width="75%"}
#Examine the predictive value of the patient's sex
data$fustat <- ifelse(data$event == "death", 1, 0)
surv_object <- Surv(time = data$futime, event = data$fustat)
fit10 <- survfit(surv_object ~ sex,data )
ggsurvplot(fit10, data = data, pval = TRUE, title="Survival probability considering the gender")
```


The p-value of 0.17 is greater than 0.05 indicates that there is no significant difference between the 2 genders.

For age and blood type, we observe the similar results: neither the age nor the blood type have an significant impact in the survival probability in the death outcome (p=0.34 and p=0.94 respectively). 



```{r, echo=FALSE, out.width="75%"}
#Examine the predictive value of the patient's age_group
data$fustat <- ifelse(data$event == "death", 1, 0)
surv_object <- Surv(time = data$futime, event = data$fustat)
fit20 <- survfit(surv_object ~ age_group,data )
ggsurvplot(fit20, data = data, pval = TRUE, title="Survival probability considering the age")
```





```{r, echo=FALSE, out.width="75%"}
#Examine the predictive value of the patient's blood type
data$fustat <- ifelse(data$event == "death", 1, 0)
surv_object <- Surv(time = data$futime, event = data$fustat)
fit30 <- survfit(surv_object ~ abo,data )
ggsurvplot(fit30, data = data, pval = TRUE, title="Survival probability considering the blood type")
```


The analysis using Kaplan Meier indicates that:

  * For the transplant survival rate, only the blood type has a significant impact
  * For the death survival rate, none of the studied covariate have a significant impact

# Advanced analysis

We have seen a basic analysis with the Kaplan-Meier model, it however presents some limitations. Most notably, it can only compare two discrete groups at a time. Secondly, it fails to take into account the multiple outcome as we've seen previously. In section, we will see more advanced models.

First the Cox Proportional Hazards (CoxPH) modeling. The main advantage is that it is capable to perform regressions, and in particular to completely model the age covariate. Like Kaplan-Meier however, it can only analyse one outcome.

Finally, we will see Competing Risks Regression (CRR) which takes into considerations all outcomes.

# Cox Model

The Cox Proportional Hazard model works by computing a baseline hazard $h_0$ and then fitting the proportional hazards $h_i$ for each studied covariate. If the proportional hazard is larger than 1, then it means the covariate presents a higher risk of having the studied outcome and vice versa compared to the baseline.

First, let's see the survival rates  for the transplant outcome:

We use the CoxPH model to model the impact of our explanatory covariates : sex, abo, age.
The model computes the Hazard Ratios (HR) and the associated Confidence interval (CI) and p-value. We consider the covariate to be statistically significant if the p-value $p > 0.05$.

```{r, echo=FALSE, out.width="75%"}
dependent_transplant <- "Surv(futime, ltx)"
explanatory   <- c("sex", "abo", "age")


cox_fit_transplant <- coxphmulti(data, dependent_transplant, explanatory)
transplcox_fit_transplant_df <- fit2df(cox_fit_transplant)
names(transplcox_fit_transplant_df) <- c("explanatory", "transplant HR (LCI, UCI, p)")
transplcox_fit_transplant_df
```

Let's analyse how each covariate impacts the transplant hazard:

  * Sex: the HR is close to 1 and the p-value is quite large. This means that there is no significant difference in the transplant chances between male and female patients.
  * Blood type: Our model compares the Hazard Ratio with respect to the A blood type. The B and O blood types are significantly different from the A blood type (p=0.009 and p<0.001 respectively). We note that for both of them, the HR < 1, meaning that the transplantation hazard is higher than compared to blood type A. This makes sense for group O because there are less donors that are compatible than for group A. For blood type B, the lower prevalence of B type people in the US population could account for the lower hazard rate.
  * There seems to be absolutely no impact of age on transplant probability (HR = 0.99, p-value=0.189). We find this result surprising.


When focusing on death survival:

```{r, echo=FALSE, out.width="75%"}
dependent_death <- "Surv(futime, death)"
explanatory   <- c("sex", "abo", "age")


cox_fit_death <- coxphmulti(data, dependent_death, explanatory)
deathcox_fit_transplant_df <- fit2df(cox_fit_death)
names(deathcox_fit_transplant_df) <- c("explanatory", "death HR (LCI, UCI, p)")

deathcox_fit_transplant_df
```

Let's analyse how each covariate impacts the death hazard:

  * Sex: the HR is lower than 1 (0.68) and the p-value greater than 0.05. This means that in our study, women are on average less likely to die while on the waiting list but this result is not statistically significative. We cannot conclude there are significant differences in the death chances between male and female patients.
  * Blood type: Our model compares the Hazard Ratio with respect to the A blood type. The B, AB and O blood types are significantly equivalent from the A blood type (p > 0.7). This means that there is no significant difference in the death chances between each blood type patients.
  * There seems to be absolutely no impact of age on death probability (HR = 1.02, p-value=0.155).

Overall, the death CoxPH analysis concludes that none of the studied covariates significantly impact the death hazard.

# Competing Risk Regression

As mentioned previously, we have multiple competing outcomes. It is important to consider that the death and transplant outcomes are competing. If we don't, we then consider that the outcome we are not studying is considered as censored. This is an incorrect reasoning as the alternative outcome should not be considered as a censored one. 

Say we study the death while waiting for a transplant. If we consider that transplants are censored observations, then the model will overestimate the death risk. Therefore the above analysis isn't entirely rigorous.

A solution to this is to consider Competing Risks Regression (CRR) first proposed by Fine and Gray in 1999. This model takes into account that we are in a competing risk setting.


We will model our data for the covariates `age`, `sex` and `abo` using the CRR model.

First focusing the the transplant hazard:

```{r, echo=FALSE, out.width="75%"}
dependent_transplant <- "Surv(futime, status_ltx)"

crr_fit_transplant <- crrmulti(data, dependent_transplant, explanatory)

transplant_result_df <- cbind(fit2df(cox_fit_transplant), fit2df(crr_fit_transplant)$HR)
names(transplant_result_df) <- c("explanatory", "transplant CoxPH HR (LCI, UCI, p)", "transplant CRR HR (LCI, UCI, p)")
transplant_result_df
```

For each covariate, we are comparing the computed HR and associated p-value of the previously computed the CoxPH model with the newly computed CRR model.

Our first observation is that the newly HR are very similar to the ones computed by CoxPH with variations of at most a few percents. We note that the computed p-values change slightly, for age and sex, they become closer to 0. This means that the CRR model finds that the effect of sex and age to be slightly more impactful on the transplant hazard ratio. However, it doesn't make them statistically significant either.

Overall, using the CRR model gives slightly different HR and p values. However, it doesn't change the previous conclusions: the statistically significative covariates remain significative. 



Secondly, focusing the the transplant hazard:

```{r, echo=FALSE, out.width="75%"}
dependent_death <- "Surv(futime, status_death)"

crr_fit_death <- crrmulti(data, dependent_death, explanatory)

death_result_df <- cbind(fit2df(cox_fit_death), fit2df(crr_fit_death)$HR)
names(death_result_df) <- c("explanatory", "death CoxPH HR (LCI, UCI, p)", "death CRR HR (LCI, UCI, p)")
death_result_df
```

Similarly to the transplant hazard, we observe some slight changes in the HR and p values when compared to the previously computed CoxPH model.

The most notable changes are for the B and O blood types where in both cases the HR increases by a factor of about 1.5. However, like in the Cox model, none of the covariate impact significantly the death rate. The lowest observed p-value is of 0.094 for the age.

Once again, using CRR instead of CoxPH yields some slight changes in the computed HR and p-values, but doesn't change the conclusions we've drawn previously.

The use of CRR doesn't seem to be all that relevant in our analysis, however, we are convinced that it is methodologically more rigourous to use it over CoxPH and Kaplan-Meier models. Especially when we see no drawback to use CRR over other models.

# Conclusion

In this work, we have studied the survival analysis of the patients waiting for a liver transplant. Our focus was to study the impact of various medical factors (sex, age and blood type) on two possible outcomes:

   * The patient receives a liver transplant
   * The patient dies while on the waiting list

We have tested different models of varying complexity. The Kaplan-Meier model, the Cox Proportional Hazards model and the Competing Risks Regression model.
All of the previous models arrived to the same conclusion:

  * For the transplant outcome:
    + Age and Sex doesn't significantly impact the transplant hazard
    + Blood types AB doesn't significantly impact the transplant hazard compared to blood type A (baseline)
    + Blood types O and B have a significant impact on the transplant hazard: their transplantation hazard is lower compared to blood type A
  * For the death outcome:
    + No covariates have a signficant impact on the death hazard

We believe that the CRR model is the most rigorous and should be used to infer predicted survival data.
