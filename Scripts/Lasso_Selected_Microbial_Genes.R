
############################## A. prepare data set for lasso  ####################################


library(stringr)
library(readr)
library(dplyr)
library(factoextra) 
library(readxl)
library(caret)


## 1. load annotated genes  
genes<- genes<- read.delim("genes_filtered_by_annot.txt")
dim(genes) #151583
## exclude the individual not including in the rest analysis
genes<- genes[,which(!colnames(genes) %in% "X986")]


## 2. load metadata 
metadata<-read.delim("metadata.txt", sep="")
metadata$Id<-paste0("X",metadata$Id)
metadata<-metadata[which(metadata$Id %in% colnames(genes)),]
## check that the two files agreed
identical(metadata$Id,colnames(genes))



##3. convert TTo to a binary factor
stress.genes<- metadata[which(metadata$Tto=="stress"),]   
control.genes<- metadata[which(metadata$Tto=="control"),] 


## standardize the data using the scale() function in R
genes<-as.data.frame(t(genes))
str(genes)
genes<-  scale(genes)
genes<-as.data.frame(genes)

## convert data to binary
metadata$Tto2<-ifelse(metadata$Tto=="stress",1,0)
metadata$Tto2<- as.factor(metadata$Tto2)
metadata$Tto<- as.factor(metadata$Tto)


##4. correct trait by Sex.  Use residuals for the model as new response!
phenoad= metadata$Tto
model <- glm(phenoad ~ metadata$Sex, family=binomial(link='logit'))
residuals1<- residuals(model)
residuals1<-scale(residuals1)
residuals1<-as.numeric(residuals1)


## 4. Create partitions to ensure no bias selection 
## Use same seed across the analysis to ensure reproducibility 
set.seed(42)

inTraining <- createDataPartition(residuals1, p = .80, list = FALSE)
genes<-cbind.data.frame(y1=residuals1,genes)
training <- genes[ inTraining,]
dim(training)
testing  <- genes[-inTraining,]
dim(testing)
rownames(testing)
write.table(rownames(testing),"testing.txt",row.names=TRUE)





################################## B. apply lasso ###############################
library(glmnet)
library(coefplot)
 

# Extract the response variable (assumed to be in the first column)
response<-training[,1]
# Extract the predictor variables (all columns except the first)
df<-training[,-1]
# Perform LASSO regression with 15-fold cross-validation using glmnet
# alpha = 1 specifies LASSO (L1 penalty)
cv_lasso_fit <- cv.glmnet(x = as.matrix(df), y = response, alpha = 1, nfolds = 15)
# Plot the mean cross-validated error curve as a function of lambda
plot(cv_lasso_fit)
# Identify the value of lambda that gives the minimum mean cross-validated error
cv_lasso_fit$lambda.min

# Display all coefficients from the cross-validated LASSO model
coef(cv_lasso_fit)
# Extract the names of predictors with non-zero coefficients at lambda.min
all.coef1<-rownames(coef(cv_lasso_fit, s = 'lambda.min'))[coef(cv_lasso_fit, s = 'lambda.min')[,1]!= 0] 
# Extract all coefficients 
all.coef<-extract.coef(cv_lasso_fit)



# 2. Fit LASSO model on the full training data (without cross-validation)
lasso_fit <- glmnet(x = as.matrix(df), y = response, alpha = 1 )
# Predict on the test set using the lambda.min from cross-validation
glm.probs <- predict(lasso_fit, s = cv_lasso_fit$lambda.min,as.matrix(testing[-1])) 
# Calculate the correlation between the observed response and predicted values
cor(as.numeric(testing[,1]),glm.probs) # 0.8

# Extract coefficients from the fitted LASSO model using the optimal lambda
lasso_coef = predict(lasso_fit, type = "coefficients", s = cv_lasso_fit$lambda.min) 

# Select and display only the non-zero coefficients (i.e., variables selected by LASSO)
coef2<-lasso_coef[lasso_coef != 0] 


write.table(all.coef1,"lasso.coef.after.adjusted.txt",row.names=TRUE)



