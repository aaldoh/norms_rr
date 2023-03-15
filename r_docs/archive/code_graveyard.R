#where code comes to die

#--------------------------
contrast_interact = list(interaction     = c(1, -1, -1,  1,  0), #interaction
                         dytext_vs_dyvis = c(1, -1,  0,  0,  0)) #text vs visual for dynamic only

model_lm = lm(attitude_T1 ~ condition, data = no_out)

lm_emmeans = emmeans(model_lm, "condition", weights = "cells")

contrast(lm_emmeans, contrast_interact, adjust="none")

ggplot2::autoplot(goggles_lm, which = c(1, 3, 2, 4), colour = "#5c97bf", smooth.colour = "#ef4836", alpha = 0.5, size = 1) + 
  theme_minimal() ## assumptions


h1_emmeans <- list()
for (i in c("attitude_T1", "interest_T1", "intentionuni_T1")) {
  model <- lm(as.formula(paste(outcome_var, "~ norm_info")), no_out)
  h1_emmeans[[i]] <- emmeans(model, "norm_info", contr = contrasts_norm, infer = T)
} 


h1.models <- map(outcomes_T1, ~ lm(substitute(i ~ norm_info, list(i = as.name(.))), data = no_out) %>%
                   emmeans(., "norm_info", contr = contrasts_norm, infer = T) %>%
                   .$contrasts %>%
                   summary()) %>%
  bind_rows(.id = "outcome")

map_df(set_names(names(ds)[1:3]), ~ 
         lm(formula(paste0(.x, "~", 
                           paste(names(ds)[4:5], collapse=" + "))), data = ds) %>%
         tidy, .id = "Dep_Variable")


#--------------------------
init <- mice(no_out, maxit = 0)
meth <- init$meth
pred <- init$pred
pred[, names(no_out)[c(1:5,32, 54:64)]] <- 0
pred[names(no_out)[c(1:5,32, 54:64)],] <- 0
meth[names(meth) %in% names(no_out)[c(1:5,32, 54:64)]] <- ""

###### split by group
impute1 <- no_out %>%
  split(.$norm_info) %>%
  map(~ mice(.x, m = 1, maxit = 1, meth = meth, pred = pred))

impute1_mapped <- impute1 %>%
  map(~ mice::complete(.x, "long", include = T))
impute1_mapped_long <- bind_rows(impute1_mapped)
impute1_mapped_df <- impute1 %>%
  map_df(~ mice::complete(.x, "long", include = T))

###### native miceadd function
impute2 <- mice(no_out, m = 1, maxit = 1, meth = meth, pred = pred, group = group)
impute2_long <- mice::complete(impute2, "long", include = T)
impute_test <- complete(impute2, 1)

###### overall imputation
impute3 <- mice(no_out, m = 1, maxit = 1, meth = meth, pred = pred)
impute3_long <- mice::complete(impute3, "long", include = T)

###### split by group
impute4 <- no_out %>%
  split(.$norm_info) %>%
  map(~ mice(.x, m = 1, maxit = 1, meth = meth, pred = pred))

impute4_rbind1 <- rbind(impute4$dynamic, impute4$static)
impute4_rbindfull <- rbind(impute4_rbind1, impute4$none)
impute4_long <- mice::complete(impute4_rbindfull, "long", include = T)

init <- mice(no_out, maxit = 0)
meth <- init$meth
pred <- init$pred
pred[, names(no_out)[c(1:5,32, 54:64)]] <- 0
pred[names(no_out)[c(1:5,32, 54:64)],] <- 0
meth[names(meth) %in% names(no_out)[c(1:5,32, 54:64)]] <- ""

imputationFunction <- as.list(meth[meth != ""])
meth[meth != ""] <- "bygroup"
group <- imputationFunction
group[] <- "norm_info"

pred[, "norm_info"] <- 0


imp <- mice::mice(no_out, meth = meth, pred = pred, m = 1, maxit = 1,
                  group = group, imputationFunction = imputationFunction)

imp_long <-  mice::complete(imp, "long", include = F)



summary(impute1_mapped_long$interest_T2) # mapped mice then mapped long
summary(impute1_mapped_df$interest_T2) # mapped mice as df
summary(impute2_long$interest_T2) # using miceadd group function
summary(impute3_long$interest_T2) # baseline - general without split
summary(imp_long$interest_T2) # the convoluted internet response

as.mids(impute1_mapped)
as.mids(impute1_mapped_df)
as.mids(impute2_long)
as.mids(impute3_long)

library(mice)
data( data.ma01, package="miceadds")
dat <- data.ma01

# use sub-dataset
dat <- dat[ dat$idschool <=1006, ]
V <- ncol(dat)
# create initial predictor matrix and imputation methods
predictorMatrix <- matrix( 1, nrow=V, ncol=V)
diag(predictorMatrix) <- 0
rownames(predictorMatrix) <- colnames(predictorMatrix) <- colnames(dat)
predictorMatrix[, c("idstud", "studwgt","urban" ) ] <- 0
method <- rep("norm", V)
names(method) <- colnames(dat)

#** groupwise imputation of variable books
method["books"] <- "bygroup"
# specify name of the grouping variable ('idschool') and imputation method ('norm')
group <- list( "books"="idschool" )
imputationFunction <- list("books"="norm" )

#** conduct multiple imputation in mice
imp <- mice::mice( dat, method=method, predictorMatrix=predictorMatrix,
                   m=1, maxit=1, group=group, imputationFunction=imputationFunction )