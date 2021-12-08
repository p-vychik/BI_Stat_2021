library(data.table)
# library(vegan)
library(dplyr)
library(Metrics)
library(kernlab)

# provide path to folder with .csv data
input_folder <- "/home/aither/institute-for-bioinformatics/r-studio/BI_Stat_2021/pca/superconduct/"

# import function
data_import <- function(input_folder, extension){
  files <- list.files(paste(input_folder), pattern = paste0(".+",extension,"$"))
  print(files)
  bind_cols(lapply(paste0(input_folder,"/",files), fread, drop = "material"))
}
source_data <- data_import(input_folder, "csv")

# remove duplicated critical_temp variable and rename column
source_data <- subset(source_data, select = -ncol(source_data))
names(source_data)
names(source_data)[82] <- "critical_temp"

# check structure
str(source_data)

# split data to training and evaluating parts
test_sep <- rbinom(nrow(source_data), 1, 0.5)
training_set <- source_data[test_sep == 0, ]
testing_set <- source_data[test_sep == 1, ]

# standardize testing data using training set

for (i in colnames(testing_set)){
  data_sd <- sd(training_set[[i]])
  if (i != "critical_temp" && data_sd != 0){
    testing_set[[i]] <- (testing_set[[i]] - mean(training_set[[i]])) / data_sd
  }
}

# build lm for critical_temp with all other variables as predictors
lm_model1 <- lm(critical_temp~., data = training_set)
# got here R-adjusted ~ 0.77 
summary(lm_model1)
mae(training_set$critical_temp, predict(lm_model1, training_set))
# mean absolute error increased on testing_set
mae(testing_set$critical_temp, predict(lm_model1, testing_set))


# PCA
pc <- prcomp(training_set[,-82])
summary(pc)

# take first eight components for scores calculations
pc_scores <- scores(pc, display = "species", choices=c(1:8), scaling = 0)
testing_set_tf <- as.data.frame(as.matrix(testing_set[,-82]) %*% pc_scores)
testing_set_tf['critical_temp'] <- testing_set$critical_temp
# build model
pca_model <- lm(critical_temp~., data = testing_set_tf)
summary(pca_model)
# evaluate
mae(testing_set$critical_temp, predict(pca_model, testing_set_tf))

# KPCA
kpca <- kpca(~., data = training_set[,-82], kernel = 'rbfdot', features = 8)
training_set_kpca <- as.data.frame(predict(kpca, training_set))
training_set_kpca['critical_temp'] <- training_set$critical_temp
kpca_model <- lm(critical_temp~., data = training_set_kpca)
summary(kpca_model)
# Adjusted R-squared:  0.02813  - significant decrease, too few components added
