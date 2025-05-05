#This file contains the companion functions for the blood_trajectory_analysis.Rmd

# Define a filtering function to be applied to each chunk for loading and filtering mimic 3
mimic3_filter <- function(x, pos) {
  subset(x, HADM_ID %in% mimic_patients_all$hadm_id)
}

# Define a filter function to be applied to each chunk
mimic4_filter <- function(x, pos) {
  subset(x, hadm_id %in% mimic_patients_all$hadm_id)
}

# Custom function for detecting and deleting extreme values
extreme_deletion<-function(x, prob = c(0.2,0.8), factor = 1.5){
  down <- quantile(x, prob[1])- factor*IQR(x)
  up <- quantile(x, prob[2])+ factor*IQR(x)
  
  x <- ifelse(x>up, NA,x)
  x <- ifelse(x<down, NA,x)
  
  return(x)
}

## Custom function for Exact fisher analysis with p.value simulation 
my_fisher <- function(data,variable, by, conf.level, ...){
  data <- data[c(variable, by)] %>% dplyr::filter(complete.cases(.))
  fisher.test(data[[variable]], factor(data[[by]]), simulate.p.value = T) %>%
    broom::tidy()
}

## Helper function to run the LCGA search
my_LCGA<-function(dataset, G=2:5){
  
  dataset<-as.data.frame(dataset)
  var_name<-dataset$label[1]
  
  ### initial 1 class model
  cat("Model: m1p1", "\n")
  
  file<-paste0("models/expLCGA_",var_name,"_linear_m1p1.Rds")
  if(!file%in%model_files){
    m1 <- lcmm(value_clean ~ time_days,random =~ 1, 
               subject = 'subject_id_num', data = dataset,maxiter = 200)
    
    saveRDS(m1, file = file)
  } else{
    m1<-readRDS(file)
  }
  
  file<-paste0("models/expLCGA_",var_name,"_beta_m1p1.Rds")
  if(!file%in%model_files){
    m1_beta <- lcmm(value_clean ~ time_days,random =~ 1, 
                    subject = 'subject_id_num',
                    data = dataset, link="beta",maxiter = 200)
    
    saveRDS(m1_beta, file = file)
  }else{
    m1_beta<-readRDS(file)
  }
  
  file<-paste0("models/expLCGA_",var_name,"_spline_m1p1.Rds")
  if(!file%in%model_files){
    m1_splines <- lcmm(value_clean ~ time_days,random =~ 1, 
                       subject = 'subject_id_num',
                       data = dataset, link="3-quant-splines",maxiter = 200)
    
    saveRDS(m1_splines, file = file)
  }else{
    m1_splines<-readRDS(file)
  }
  
  cat("Model: m1p2", "\n")
  
  file<-paste0("models/expLCGA_",var_name,"_linear_m1p2.Rds")
  if(!file%in%model_files){
    m1p2 <- lcmm(value_clean ~ poly(time_days, 2),random =~ 1, 
                 subject = 'subject_id_num', data = dataset,maxiter = 200)
    
    saveRDS(m1p2, file = file)
  }else{
    m1p2<-readRDS(file)
  }
  
  file<-paste0("models/expLCGA_",var_name,"_beta_m1p2.Rds")
  if(!file%in%model_files){
    m1p2_beta <- lcmm(value_clean ~ poly(time_days, 2),random =~ 1, 
                      subject = 'subject_id_num', data = dataset, 
                      link = "beta",maxiter = 200)
    
    saveRDS(m1p2_beta, file = file)
  }else{
    m1p2_beta<-readRDS(file)
  }
  
  file<-paste0("models/expLCGA_",var_name,"_splines_m1p2.Rds")
  if(!file%in%model_files){
    m1p2_splines <- lcmm(value_clean ~ poly(time_days, 2),random =~ 1, 
                         subject = 'subject_id_num', data = dataset, 
                         link = "3-quant-splines",maxiter = 200)
    
    saveRDS(m1p2_splines, file = file)
  }else{
    m1p2_splines<-readRDS(file)
  }
  
  cat("Model: m1p3", "\n")
  
  file<-paste0("models/expLCGA_",var_name,"_linear_m1p3.Rds")
  if(!file%in%model_files){
    m1p3 <- lcmm(value_clean ~ poly(time_days, 3),random =~ 1, 
                 subject = 'subject_id_num', data = dataset,maxiter = 200)
    
    saveRDS(m1p3, file = file)
  }else{
    m1p3<-readRDS(file)
  }
  
  file<-paste0("models/expLCGA_",var_name,"_beta_m1p3.Rds")
  if(!file%in%model_files){
    m1p3_beta <- lcmm(value_clean ~ poly(time_days, 3),random =~ 1, 
                      subject = 'subject_id_num', data = dataset, 
                      link = "beta",maxiter = 200)
    
    saveRDS(m1p3_beta, file = file)
  }else{
    m1p3_beta<-readRDS(file)
  }
  
  file<-paste0("models/expLCGA_",var_name,"_splines_m1p3.Rds")
  if(!file%in%model_files){
    m1p3_splines <- lcmm(value_clean ~ poly(time_days, 3),random =~ 1, 
                         subject = 'subject_id_num', data = dataset, 
                         link = "3-quant-splines",maxiter = 200)
    
    saveRDS(m1p3_splines, file = file)
  }else{
    m1p3_splines<-readRDS(file)
  }
  
  
  for (k in G){
    cat("Model: k",k,"p1", "\n")
    
    file<-paste0("models/expLCGA_",var_name,"_linear_k",k,"p1.Rds")
    if(!file%in%model_files){
      temp_k <- lcmm(value_clean ~ time_days,random =~ 1, 
                     subject = 'subject_id_num', 
                     data = dataset, mixture = ~time_days, 
                     ng = k, B = m1,maxiter = 200)
      
      saveRDS(temp_k, file = file)
    }
    
    file<-paste0("models/expLCGA_",var_name,"_beta_k",k,"p1.Rds")
    if(!file%in%model_files){
      temp_k_beta <- lcmm(value_clean ~ time_days,random =~ 1, 
                          subject = 'subject_id_num',
                          data = dataset, mixture = ~time_days, 
                          ng = k, B = m1_beta, link = "beta",maxiter = 200)
      
      saveRDS(temp_k_beta, file = file)
    }
    
    file<-paste0("models/expLCGA_",var_name,"_splines_k",k,"p1.Rds")
    if(!file%in%model_files){
      temp_k_splines <- lcmm(value_clean ~ time_days,random =~ 1, 
                             subject = 'subject_id_num', 
                             link = "3-quant-splines",
                             data = dataset, mixture = ~time_days, ng = k, 
                             B = m1_splines,maxiter = 200)
      
      saveRDS(temp_k_splines, file = file)
    }
    
    cat("Model: k",k,"p2", "\n")
    
    file<-paste0("models/expLCGA_",var_name,"_linear_k",k,"p2.Rds")
    if(!file%in%model_files){
      temp_kp2 <- lcmm(value_clean ~ poly(time_days, 2),random =~ 1, 
                       subject = 'subject_id_num', 
                       data = dataset, mixture = ~poly(time_days, 2), 
                       ng = k, B = m1p2,maxiter = 200)
      
      saveRDS(temp_kp2, file = file)
    }
    
    file<-paste0("models/expLCGA_",var_name,"_beta_k",k,"p2.Rds")
    if(!file%in%model_files){
      temp_kp2_beta <- lcmm(value_clean ~ poly(time_days, 2),random =~ 1, 
                            subject = 'subject_id_num', 
                            data = dataset, mixture = ~poly(time_days, 2), 
                            ng = k, B = m1p2_beta, link = "beta",maxiter = 200)
      
      saveRDS(temp_kp2_beta, file = file)
    }
    
    file <- paste0("models/expLCGA_",var_name,"_splines_k",k,"p2.Rds")
    if(!file%in%model_files){
      temp_kp2_splines <- lcmm(value_clean ~ poly(time_days, 2),random =~ 1, 
                               subject = 'subject_id_num', 
                               data = dataset, mixture = ~poly(time_days, 2), 
                               ng = k, B = m1p2_splines, link = "3-quant-splines",
                               maxiter = 200)
      
      saveRDS(temp_kp2_splines, file = file)
    }
    
    cat("Model: k",k,"p3", "\n")
    
    file<-paste0("models/expLCGA_",var_name,"_linear_k",k,"p3.Rds")
    if(!file%in%model_files){
      temp_kp3 <- lcmm(value_clean ~ poly(time_days, 3),
                       random =~ 1, subject = 'subject_id_num', 
                       data = dataset,
                       mixture = ~poly(time_days, 3),ng = k, 
                       B = m1p3,maxiter = 200)
      
      saveRDS(temp_kp3, file = file)
    }
    
    file<-paste0("models/expLCGA_",var_name,"_beta_k",k,"p3.Rds")
    if(!file%in%model_files){
      temp_kp3_beta <- lcmm(value_clean ~ poly(time_days, 3),
                            random =~ 1, subject = 'subject_id_num', 
                            data = dataset, link = "beta",
                            mixture = ~poly(time_days, 3),ng = k, 
                            B = m1p3_beta,maxiter = 200)
      
      saveRDS(temp_kp3_beta, file = file)
    }
    
    file<-paste0("models/expLCGA_",var_name,"_splines_k",k,"p3.Rds")
    if(!file%in%model_files){
      temp_kp3_splines <- lcmm(value_clean ~ poly(time_days, 3),
                               random =~ 1, subject = 'subject_id_num', 
                               data = dataset, link = "3-quant-splines",
                               mixture = ~poly(time_days, 3),ng = k, 
                               B = m1p3_splines,maxiter = 200)
      
      saveRDS(temp_kp3_splines, file = file)
    }
  }
}

summary_metrics<-function(x, i, .modeltype = 6, .analyte = 4, .linktype = 5){
  
  ## x: a model object fitted with lcmm
  ## i: the name of the model
  temp_x<-as.data.frame(summarytable(x, which = c("G", "conv", "npm", "AIC", 
                                                  "BIC", "entropy","ICL", "%class")))
  
  for (c in c("%class2", "%class3", "%class4", "%class5")){
    if (!c%in%colnames(temp_x)){
      temp_x[,c]<-NA
    }
  }
  
  if( temp_x[1,]$G != 1){
    pprob<-x$pprob
    pprob$PPA<-apply(pprob[,3:ncol(pprob)], 1, max)
    temp_x$APPA = mean(pprob$PPA)
  }else{
    temp_x$APPA = 1
  }
  print(i)
  modeltype<-str_split(i, "_", simplify = T)[,.modeltype]
  
  temp_x$analyte<-str_split(i, "_", simplify = T)[,.analyte]
  temp_x$linktype<-str_split(i, "_", simplify = T)[,.linktype]
  temp_x$modeltype<-str_sub(modeltype, 1, 4)
  temp_x$poly = str_sub(modeltype, 3,4)
  
  return(temp_x)
}

#Extracts the model selection metrics for each marker
model_metrics<-function(lab_name, models_list){
  
  lab_list<-models_list[[lab_name]]
  
  table_m_selection_l<-mapply(
    summary_metrics(x, i),
    lab_list, 
    names(lab_list), SIMPLIFY = F)
  
  table_m_selection<-bind_rows(table_m_selection_l)
  
  table_m_selection<-table_m_selection%>%
    mutate(linktype = ifelse(linktype == "splines", "spline", linktype),
           linktype = factor(linktype, levels = c("linear", "beta", "spline")))
  
  table_m_selection_plot<-table_m_selection%>%
    filter(conv %in% c(1,2))%>%
    select(G, BIC, ICL, APPA, poly, linktype)%>%
    pivot_longer(-c(G, poly, linktype))%>%
    mutate(name = factor(name, levels = c("BIC", "ICL", "APPA")))
  
  figure<-ggplot(table_m_selection_plot, aes(G, value, color = poly,
                                             shape = linktype, linetype = linktype))+
    geom_point()+
    geom_line()+
    facet_wrap(~name, scales = "free")+
    xlab("Number of classes")+
    ylab(NULL)+
    theme_minimal()+
    ggtitle(unique(table_m_selection$analyte))
  
  return(list("lab_name" = lab_name, "table"= table_m_selection, 
              "fig" = figure))
}

# Custom function for the fitting of the GMM models based on the provided parameters
my_GMM<-function(analyte, link, k, d, maxiter, save = TRUE, .return = FALSE, 
                 .nproc =1, seed = 500){
  
  ## Custom function for the fitting of the GMM models based on the provided parameters
  ## @ analyte: string with the name of the analyte
  ## @ link: string with the name of the link function
  ## @ k: number of classes to fit
  ## @ d: degree of the polynomial for the time trend
  
  # analyte<-"Glucose"
  set.seed(500)
  dataset<-as.data.frame(modeling_set_list[[analyte]])
  
  file<-paste0("models_GMM/GMM_",analyte,"_", link, "_m",k,"p",d,".Rds")
  
  if (d == 1){
    model_form<-"value_clean ~ time_days"
    random<-"~ time_days"
  }else{
    model_form<-paste0("value_clean ~ splines::ns(time_days,",d,")")
    random<-paste0("~ splines::ns(time_days,",d,")")
  }
  
  if (k == 1){
    if(!file%in%model_files){
      
      sink("GMM fitting log.txt",append = T)
      cat(Sys.time(),analyte, link, k, d,"\n")
      sink()
      
      m1 <- lcmm(as.formula(model_form),random =as.formula(random), 
                 subject = 'subject_id_num', data = dataset,link = link,
                 maxiter = maxiter, nproc = .nproc)
      if (save) saveRDS(m1, file = file)
      if (.return) return(m1)
    }
  }else{
    
    
    error<-try(
      m1<-readRDS(paste0("models_GMM/GMM_",analyte,"_", 
                         link, "_m1p",d,".Rds"))
    )
    
    if (class(error) == "try_error"){
      sink("GMM fitting log.txt",append = T)
      cat(Sys.time(), analyte, link, k, d,"\n")
      sink()
    }
    
    if(!file%in%model_files){
      
      sink("test.txt",append = T)
      cat(Sys.time(), analyte, link, k, d,"\n")
      sink()
      
      m2 <- lcmm(as.formula(model_form),random =as.formula(random), 
                 subject = 'subject_id_num', data = dataset,link = link,
                 ng =k, mixture = as.formula(random), B = m1, 
                 maxiter = maxiter)
      
      if (save) saveRDS(m2, file = file)
      if (.return) return(m2)
    }
  }
}

## Helper function to plot the GMM trajectories
plot_GMM<-function(gmm_model, analyte, k, ndraws = 100, .margins = 5){
  
  class_size<-table(gmm_model$pprob$class)
  
  gmm_model$call$data<-new_dataset
  pred<-predictY(gmm_model, new_dataset, draws = T, ndraws = ndraws,
                 var.time = "time_days")
  
  if(gmm_model$conv ==1){
    if(k == 1){
      pred_df<-do.call(cbind, pred)%>%
        pivot_longer(-time_days)%>%
        mutate(class = "class1",
               quantile = str_split (name, "_", simplify = T)[,2])
    }else{
      pred_df<-do.call(cbind, pred)%>%
        pivot_longer(-time_days)%>%
        mutate(class = str_split(name, "_", simplify = T)[,3],
               quantile = str_split (name, "_", simplify = T)[,2])
    }
    
    p<-pred_df%>%
      filter(quantile == 50)%>%
      ggplot(aes(time_days, value, color = class, group = class))+
      geom_path(size = 1)+
      geom_path(data = pred_df%>%filter(quantile == 2.5), linetype = "dashed")+
      geom_path(data = pred_df%>%filter(quantile == 97.5), linetype = "dashed")+
      scale_color_manual(values = c("firebrick1","steelblue1", "chartreuse3"),
                         labels = paste0("class ", 1:as.numeric(k), " (n=", class_size,")"))+
      xlab("Days after hospital arrival")+
      ylab("Marker value")+
      theme_minimal(base_size = 16)+
      ggtitle(analyte)+
      theme(plot.margin = unit(rep(.margins, 4),units = "mm"), 
            legend.position = "bottom", legend.title = element_blank(),
            legend.text = element_text(size = 16))+
      guides(color=guide_legend(ncol=2))
    
  }else{
    pred_df<-do.call(cbind, pred)%>%
      pivot_longer(-time_days)%>%
      mutate(class = str_split(name, "_", simplify = T)[,2])
    
    p<-pred_df%>%
      ggplot(aes(time_days, value, color = class, group = class))+
      geom_path(size = 1)+
      xlab("Days after hospital arrival")+
      ylab("Marker value")+
      scale_color_manual(values = c("firebrick1","steelblue1", "chartreuse3"), 
                         labels = paste0("class ", 1:as.numeric(k), " (n=", class_size,")"))+
      theme_minimal(base_size = 16)+
      ggtitle(analyte)+
      theme(plot.margin = unit(rep(.margins, 4), units = "mm"), 
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 16))+
      guides(color=guide_legend(ncol=2))
  }
  
  return(p)
}
#for extracting prediction performance
pred_performance<-function(fit, X, Y, positive = "yes", negative = "no", roc = TRUE){
  
  # posterior prob of learner
  prediction<-predict(fit, X, type = "prob")
  prediction$truth<-relevel(Y, ref = positive)
  
  # posterior prob of majority class
  prevalence_positive<-mean(Y==positive)
  
  # Precision - Recall 
  pred_prROC<-pr_curve(prediction,truth,!!sym(positive))
  pred_prAUC<-pr_auc(prediction, truth, !!sym(positive))
  
  
  # Sensitivity - specificity
  pred_ROC<-roc_curve(prediction,truth,!!sym(positive))
  pred_AUC<-roc_auc(prediction, truth, !!sym(positive))
  
  
  return(list("pred_prROC" = pred_prROC, "pred_prAUC" = pred_prAUC,
              "prevalence_positive" = prevalence_positive,
              "pred_ROC" = pred_ROC, "pred_AUC" = pred_AUC))
}

# to run the ML experiments
run_experiment<-function(prediction_df, target.name, predictor_start = 3,
                         seed = 555, reps = 25, yes_name = "yes",no_name = "no",
                         n.yes.train, n.no.train, resampling = T, 
                         method = "glmnet"){
  
  # prediction_df = predict_df 
  # target.name = "is_SCI"
  # n.yes.train = n_yes*0.8
  # n.no.train = n_no*0.8
  # resampling = FALSE
  # reps = 2
  # method = "glmnet"
  
  set.seed(seed)
  
  ## Training with resampling and repeated model building
  exp_fit_list<-list()
  exp_test_fit<-list()
  exp_dummy_fit<-list()
  exp_varimp<-list()
  
  for (i in 1:reps) {
    print(i)
    
    yes_class<-prediction_df%>%
      filter(!!sym(target.name) == yes_name)
    yes_index<-sample(1:nrow(yes_class), n.yes.train, replace = F)
    yes_train_df<-yes_class[yes_index, ]
    yes_test_df<-yes_class[-yes_index, ]
    
    no_class<-prediction_df%>%
      filter(!!sym(target.name) == no_name)
    no_index<-sample(1:nrow(no_class), n.no.train, replace = F)
    no_train_df<-no_class[no_index, ]
    no_test_df<-no_class[-no_index, ]
    
    ## Testing data preparation
    test_df<-rbind(yes_test_df, no_test_df)
    dmy <- dummyVars(" ~ .", data = test_df[,predictor_start:ncol(test_df)])
    
    X_test <- data.frame(predict(dmy, newdata = test_df))
    Y_test<-as.factor(unlist(test_df[,target.name]))
    
    
    if(resampling){
      new_yes_train_df<-yes_train_df[sample(1:nrow(yes_train_df), nrow(no_train_df),
                                            replace = T), ]
      
      temp_train_df<-rbind(new_yes_train_df, no_train_df)
      train_df<-temp_train_df[,predictor_start:ncol(temp_train_df)]
      
      # dmy <- dummyVars(" ~ .", data = train_df)
      X <- data.frame(predict(dmy, newdata = train_df))
      X <- X[colnames(X)%in%colnames(X_test)]
      
      Y <- as.factor(unlist(temp_train_df[, target.name]))
      
    }else{
      temp_train_df<-rbind(yes_train_df, no_train_df)
      train_df<-temp_train_df[,predictor_start:ncol(temp_train_df)]
      
      dmy <- dummyVars(" ~ .", data = train_df)
      X <- data.frame(predict(dmy, newdata = train_df))
      X <- X[colnames(X)%in%colnames(X_test)]
      
      Y <- as.factor(unlist(temp_train_df[, target.name]))
    }
    
    folds <- 5
    cvIndex <- createFolds(factor(Y), folds, returnTrain = T)
    
    if (method == "lr"){
      
      tc <- trainControl(
        index = cvIndex,
        method = 'cv',
        # repeats = 25,
        number = folds
      )
      
      fit <- train( x = X, y = Y, method = "glm", family = "binomial",
                    maxit= 100, trControl = tc)
      
      exp_fit_list[[i]] <- pred_performance(
        fit, X, Y, positive = yes_name,  negative = no_name
      )
      
      exp_test_fit[[i]] <- pred_performance(
        fit, X_test, Y_test, positive = yes_name, negative = no_name
      )
      
    }else{
      tc <- trainControl(index = cvIndex,
                         method = 'cv', 
                         number = folds)
      
      tuneGrid <- expand.grid(alpha = seq(0, 1, length = 10), 
                              lambda = seq(0.0001, 1, length = 10))
      
      fit<-try({
        train( x = X, y = Y, method = method,trControl = tc, tuneGrid = tuneGrid,metric = "Kappa")
      })
      
      if (class(fit) == "try-error"){
        exp_fit_list[[i]]<-NA
        next()
      }
      
      exp_fit_list[[i]]<-pred_performance(fit, X, Y, 
                                          positive = yes_name, negative = no_name)
      
      exp_test_fit[[i]]<-pred_performance(fit, X_test, Y_test,
                                          positive = yes_name, negative = no_name)
      
      exp_varimp[[i]] <- varImp(fit)$importance 
    }
    
  }
  
  return(list("train_list" = exp_fit_list, "test_list" = exp_test_fit, "var_imp" = exp_varimp))
}

experiment_performance<-function(exp_list, rep = 25){
  res_list<-mapply(function(i, name_i) {
    
    # Loop over 25 repetitions for each cutoff
    temp_df <- lapply(1:rep, function(j) {
      df <- data.frame(
        roc_auc = unlist(c(
          i$train_list[[j]]$pred_AUC[, 3],  # Train ROC AUC
          i$test_list[[j]]$pred_AUC[, 3]    # Test ROC AUC
        )),
        pr_auc = unlist(c(
          i$train_list[[j]]$pred_prAUC[, 3],  # Train PR AUC
          i$test_list[[j]]$pred_prAUC[, 3]    # Test PR AUC
        )),
        pr_BL = unlist(c(
          i$train_list[[j]]$prevalence_positive,  # Train class balance
          i$test_list[[j]]$prevalence_positive    # Test class balance
        ))
      )
      
      df$per_type <- c("in-train", "out-train")  
      df$cutoff <- name_i                        
      
      return(df)
    })
    
    # Bind rows from all reps into a single data frame
    temp_df <- do.call(rbind, temp_df)
    temp_df
  }, exp_list, names(exp_list), SIMPLIFY = FALSE)
  
  # Combine all performance data across cutoff values into a single long-format dataframe
  res<-do.call(rbind, res_list)
  
  return(res)
}



get_varimp_plot <- function(day, dat, mycol, ex, va) {
  
  # Extract list of variable importance vectors for a given day
  varimp_list <- lapply(dat[[as.character(day)]][["var_imp"]], function(x) {
    v <- as.numeric(x$Overall)
    if (length(v) == 134) v <- v[-132]
    if (length(v) == 137) v <- v[-135]
    return(v)
  })
  
  # Combine into matrix
  varimp_mat <- do.call(cbind, varimp_list)
  
  # Set biomarker names (row names)
  biomarker_names <- rownames(dat[[as.character(day)]][["var_imp"]][[1]])
  if (length(biomarker_names) == 134) biomarker_names <- biomarker_names[-132]
  if (length(biomarker_names) == 137) biomarker_names <- biomarker_names[-135]
  biomarker_names <- sub("^X\\.", "", biomarker_names)
  
  rownames(varimp_mat) <- biomarker_names
  
  # Transpose: rows = samples, columns = biomarkers
  varimp_df <- as.data.frame(t(varimp_mat))
  
  # Pivot to long format
  long_df <- pivot_longer(varimp_df, cols = everything(),
                          names_to = "Biomarker", values_to = "Importance")
  
  # Summary statistics
  stats_df <- long_df %>%
    group_by(Biomarker) %>%
    summarise(Mean = mean(Importance, na.rm = T), SD = sd(Importance,na.rm = T),
              N = length(!is.na(Importance)),
              .groups = "drop") %>%
    mutate(Color = case_when(
      Mean < 20 ~ "< 20",
      Mean > 20 & Mean < 40 ~ "M",
      Mean > 40 ~ ">40"
    ), SE = SD/sqrt(N))
  
  # Bar plot for biomarkers with Mean > 40
  p <- stats_df %>%
    filter(Color == ">40") %>%
    ggplot(aes(x = reorder(Biomarker, Mean), y = Mean)) +
    geom_bar(stat = "identity", fill = mycol) +
    geom_errorbar(aes(ymin = Mean - 1.96 * SE, ymax = Mean + 1.96 * SE),
                  width = 0.1, colour = "black") +
    labs(x = NULL, y = "Mean Importance",
         title = paste0("Experiment ", ex, ": Day ", day, " ", va)) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, 
                                 face = "bold", color = "black"),
      axis.text.y = element_text(size = 15, face = "bold", color = "black")
    ) +
    scale_y_continuous(limits = c(0, 110))+
    coord_flip()

  
  return(p)
}



varimp_function <- function(dat = NULL, ex = NULL, type = NULL, va = NULL, mycol = NULL) {
  
  days <- c(1, 3, 7, 14, 21)
  plot_list <- lapply(days, get_varimp_plot, dat = dat, mycol = mycol, ex = ex, va = va)
  
  varimp_figure <- ggarrange(plotlist = plot_list,
                             ncol = 1,
                             align = "v",
                             common.legend = TRUE,
                             legend = "none")
  
  ggsave(plot = varimp_figure,
         filename = paste0("figures/Variable Importance plots/varimp_figure_", ex, "_", type, ".png"),
         width = 4, height = 8)
  
  return(varimp_figure)
}
