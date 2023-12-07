library("BOIN")
setwd("C:/Users/katie.wang/OneDrive - Karyopharm Therapeutics/Documents/BOIN")


################Get Boundary table format###############
boin_decision_table<-function(target, ncohort, cohortsize,
                                     n.earlystop, cutoff.eli,
                                     extrasafe, offset ){
      bound=get.boundary(target = target,ncohort = ncohort,
                         cohortsize=cohortsize,
                         n.earlystop=n.earlystop, 
                         p.saf=0.6*target,
                         p.tox=1.4*target,
                         cutoff.eli=cutoff.eli, 
                         extrasafe=extrasafe,
                         offset=offset)
      decision_table = c()
      
      
      for (i in 1:n.earlystop){
            column = rep(0, n.earlystop+1)
            for(j in 0: (n.earlystop)){
                 
                  if (j <= (bound$full_boundary_tab[2,i]) ){
                        column[j+1]="E"
                        
                  }else if(j > (bound$full_boundary_tab[2,i])& 
                           j < (bound$full_boundary_tab[3,i])){
                        column[j+1]="S"
                  }
                  else if(j >= (bound$full_boundary_tab[3,i])){
                        column[j+1]= "D"
                        
                        if(i >=3 &!is.na(bound$full_boundary_tab[4,i])){
                              if(j>= (bound$full_boundary_tab[4,i])){
                                    column[j+1]= "DE"
                              }
                        }else{
                              if(i ==1 & j >=2){
                                    column[j+1]= ""   
                              }else if(i == 2 & j>=3){
                                    column[j+1]= "" 
                              }
                        
                               
                        }

                  }
            }
        
            decision_table = cbind(decision_table, column)
           
      }
      decision_table = as.data.frame(decision_table)
      colnames(decision_table)= 1:n.earlystop
      rownames(decision_table) = 0:n.earlystop
      return(decision_table)
}

############ OUTPUT DECISION BOUNDARY #####################

get_boin_decision<-function(file_dir, target, ncohort, cohortsize,
                                     n.earlystop, cutoff.eli,
                                     extrasafe, offset ){
      bound=get.boundary(target = target,
                         ncohort = ncohort,
                         cohortsize=cohortsize,
                         n.earlystop=n.earlystop, 
                         p.saf=0.6*target,
                         p.tox=1.4*target,
                         cutoff.eli=cutoff.eli, 
                         extrasafe=extrasafe,
                         offset=offset)
      
      decision_table = boin_decision_table(target = target,
                                               ncohort = ncohort,
                                               cohortsize=cohortsize,
                                               n.earlystop=n.earlystop,
                                               cutoff.eli=cutoff.eli,
                                               extrasafe=extrasafe,
                                             offset=offset)
      
      
      parameters = c("target", "ncohort", "cohortsize","n.earlystop",
                     "cutoff.eli","extrasafe","offset")
      
      
      meaning = c("the target toxicity rate",
                  "the total number of times to add patients",
                  "number of patients in each cohort",
                  "the max number of patients assigned to each dose",
                  "The cutoff to eliminate the overly toxic dose for safety",
                  "Set extrasafe=TRUE to impose a stricter stopping rule",
                  "A small positive number(between0and0.5) to control how strict the stopping rule is")
      
      values = c(target, ncohort, cohortsize,n.earlystop,
                 cutoff.eli, extrasafe, offset)
      df = data.frame(parameters, meaning, values)
      
      
      
      #write out the decision table and decision boundary
      write.table(df,file_dir , na = "NA", append = TRUE,
                  col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
      
      write.table('\n', file_dir, append = TRUE,
                  col.names = FALSE, row.names = FALSE, sep = ",", quote = FALSE)
     
       write.table(bound$full_boundary_tab,file_dir , na = "NA", append = TRUE,
                  col.names = FALSE, row.names = TRUE, sep = ",", quote = FALSE)

      write.table('\n', file_dir, append = TRUE,
                  col.names = FALSE, row.names = FALSE, sep = ",", quote = FALSE)
      
      
      write.table(decision_table, file_dir, na = "NA", append = TRUE, 
                  col.names = NA, row.names = TRUE, sep = ",", quote = FALSE)
      
      write.table('\n', file_dir, append = TRUE, 
                  col.names = FALSE, row.names = FALSE, sep = ",", quote = FALSE)
}


################# Summary Simulation Results ############################
summary.boin_csv<- function (file_dir, object, ...)
{
      if (!is.null(object$boundary_tab)) {
            if (!is.na(object$lambda_e))
                  cat("Escalate dose if the observed DLT rate at the current dose <= ",
                      object$lambda_e, "\n")
            if (!is.na(object$lambda_d))
                  cat("Deescalate dose if the observed DLT rate at the current dose > ",
                      object$lambda_d, "\n\n")
            if (!is.null(object$boundary_tab)) {
                  cat("This is equivalent to the following decision boundaries\n")
                  print(object$boundary_tab)
            }
            if (!is.null(object$full_boundary_tab)) {
                  cat("\n")
                  cat("A more completed version of the decision boundaries is given by\n")
                  print(object$full_boundary_tab)
            }
            if (!is.null(object$stop_boundary)) {
                  cat("\n")
                  cat("In addition to the default stopping rule (i.e., stop the trial if the lowest dose is eliminated), \n")
                  cat("the following more strict stopping safety rule will be used for extra safety: \n")
                  cat(" stop the trial if (1) the number of patients treated at the lowest dose >= 3 AND",
                      "\n", "(2) Pr(the DLT rate of the lowest dose >",
                      object$target, "| data) > ", object$cutoff, ",\n",
                      "which corresponds to the following stopping boundaries:\n")
                  print(object$stop_boundary)
            }
            else {
                  cat("\n")
                  cat("Default stopping rule: stop the trial if the lowest dose is eliminated.\n")
            }
      }
      if (!is.null(object$next_subtrial)) {
            if (is.na(object$next_subtrial) == TRUE) {
                  cat("No additional next subtrials are needed!!\n")
            }
            else {
                  cat("Next subtrial includes doses: ", "\n")
                  cat("\t\t", object$next_subtrial, "\n\n")
                  cat("The starting dose for this subtrial is:\n",
                      "\t\t", paste("(", object$starting_dose[1], ", ",
                                    object$starting_dose[2], ")", sep = ""), "\n")
            }
      }
      if (!is.null(object$next_dc)) {
            if (is.na(object$next_dc[1]) == TRUE) {
                  cat("The trial experienced an early stopping.")
            }
            else {
                  cat("The recommended dose combination for the next cohort of patients is (",
                      object$next_dc[1], ", ", object$next_dc[2], ").",
                      "\n")
            }
      }
      if (!is.null(object$MTD)) {
            if (length(object$MTD) == 1) {
                  if (object$MTD == 99) {
                        cat("All tested doses are overly toxic. No MTD should be selected! \n\n")
                  }
                  else {
                        cat("The MTD is dose level ", object$MTD, "\n\n")
                  }
                  cat("Dose    Posterior DLT             95%                  \n",
                      sep = "")
                  cat("Level     Estimate         Credible Interval   Pr(toxicity>",
                      object$target, "|data)\n", sep = "")
                  for (i in 1:nrow(object$p_est)) {
                        cat(" ", i, "        ", as.character(object$p_est[i,
                                                                          2]), "         ", as.character(object$p_est[i,
                                                                                                                      3]), "         ", as.character(object$p_overdose[i]),
                            "\n")
                  }
                  cat("NOTE: no estimate is provided for the doses at which no patient was treated.\n")
            }
            if (length(object$MTD) >= 2) {
                  if (length(object$MTD) == 2) {
                        if (object$MTD[1, 1] == 99 && object$MTD[1, 2] ==
                            99) {
                              cat("All tested doses are overly toxic. No MTD is selected! \n")
                        }
                        else cat("The MTD is dose combination (", object$MTD[1,
                                                                             1], ", ", object$MTD[1, 2], ") \n\n")
                  }
                  else {
                        if (length(object$MTD) == 0) {
                              cat("All tested doses are overly toxic. No MTD is selected! \n")
                        }
                        else {
                              cat("The MTD contour includes dose combinations ",
                                  paste("(", object$MTD[, 1], ", ", object$MTD[,
                                                                               2], ")", sep = ""), "\n\n")
                        }
                  }
                  cat("Isotonic estimates of toxicity probabilities and 95% confidence intervals for combinations are \n")
                  # for (i in 1:dim(object$p_est_CI)[1]) {
                  #   cat(formatC(object$p_est_CI[i, ], digits = 2, format = "f",
                  #               width = 5), sep = "  ", "\n")
                  # }
                  print(noquote(object$p_est_CI))
                  cat("\n")
                  cat("NOTE: no estimate is provided for the doses at which no patient was treated.\n\n")
            }
      }
      if (!is.null(object$percentstop)) {
            if (!is.null(object$overdose60)) {
                  cat("True DLT rate(%):\n")
                  cat(formatC(object$p.true, digits = 2, format = "f"),
                      sep = "  ", "\n")
                  cat("selection percentage at each dose level (%):\n")
                  cat(formatC(object$selpercent, digits = 1, format = "f"),
                      sep = "  ", "\n")

                  cat("average number of patients treated at each dose level:\n")
                  cat(formatC(object$npatients, digits = 1, format = "f"),
                      sep = "  ", "\n")
                  cat("average number of toxicity observed at each dose level:\n")
                  cat(formatC(object$ntox, digits = 1, format = "f"),
                      sep = "  ", "\n")
                  cat("average DLT percentage per dose:", formatC(object$perct_DLT_per_dose, digits=2, 
                                                                  format="f"), "\n")
                  cat("average percentage of patients per dose:", formatC(apply(object$perct_pts_per_dose,2,mean), digits=4, 
                                                                          format="f"), "\n")
                  cat("average number of toxicities:", formatC(object$totaltox,
                                                               digits = 1, format = "f"), "\n")
                  cat("average number of patients:", formatC(object$totaln,
                                                             digits = 1, format = "f"), "\n")
                  cat("average trial_duration(month):", formatC(object$average_trial_duration, digits=1, format="f"), "\n")
                  cat("25% quantile trial_duration(month):", formatC(object$twentyFive_Perct_trial_duration, digits=2, format="f"), "\n")
                  cat("75% quantile trial_duration(month):", formatC(object$senventyFive_Perct_trial_duration, digits=2, format="f"), "\n")
                  
                  
                  cat("risk of poor allocation:", formatC(object$risk_poor_allocation, digits=2, format="f"), "\n")
                  cat("risk of high toxicity:", formatC(object$risk_high_toxicity, digits=2, format="f"), "\n")
                  
                  cat("percentage of early stopping due to toxicity:",
                      formatC(object$percentstop, digits = 1, format = "f"),
                      "% \n")
                  cat("risk of overdosing (>60% of patients treated above the MTD):",
                      formatC(object$overdose60, digits = 1, format = "f"),
                      "% \n")
                  cat("risk of overdosing (>80% of patients treated above the MTD):",
                      formatC(object$overdose80, digits = 1, format = "f"),
                      "% \n")
                  cat("\n")
            }
            else {
                  cat("True DLT rate(%):\n")
                  cat(formatC(object$p.true, digits = 2, format = "f"),
                      sep = "  ", "\n")
                  cat("selection percentage at each dose level (%):\n")
                  cat(formatC(object$selpercent, digits = 1, format = "f"),
                      sep = "  ", "\n")

                  cat("average number of patients treated at each dose level:\n")
                  cat(formatC(object$npatients, digits = 1, format = "f"),
                      sep = "  ", "\n")
                  cat("average number of toxicity observed at each dose level:\n")
                  cat(formatC(object$ntox, digits = 1, format = "f"),
                      sep = "  ", "\n")
                  cat("average DLT percentage per dose:", formatC(object$perct_DLT_per_dose, digits=2, 
                                                                  format="f"), "\n")
                  cat("average percentage of patients per dose:", formatC(apply(object$perct_pts_per_dose,2,mean), digits=4, 
                                                                          format="f"), "\n")
                  cat("average number of toxicities:", formatC(object$totaltox,
                                                               digits = 1, format = "f"), "\n")
                  cat("average number of patients:", formatC(object$totaln,
                                                             digits = 1, format = "f"), "\n")
                  cat("average trial_duration(month):", formatC(object$average_trial_duration, digits=1, format="f"), "\n")
                  cat("25% quantile trial_duration(month):", formatC(object$twentyFive_Perct_trial_duration, digits=2, format="f"), "\n")
                  cat("75% quantile trial_duration(month):", formatC(object$senventyFive_Perct_trial_duration, digits=2, format="f"), "\n")
                  cat("risk of poor allocation:", formatC(object$risk_poor_allocation, digits=2, format="f"), "\n")
                  cat("risk of high toxicity:", formatC(object$risk_high_toxicity, digits=2, format="f"), "\n")
                  
                  cat("percentage of early stopping due to toxicity:",
                      formatC(object$percentstop, digits = 1, format = "f"),
                      "% \n")
                  cat("\n")
            }
      }
      if (!is.null(object$npercent) | !is.null(object$npercent.contour)) {
            if (!is.null(object$npercent.contour)) {
                  cat("true DLT rate of dose combinations:\n")
                  for (i in 1:dim(object$p.true)[1]) cat(formatC(object$p.true[i,
                  ], digits = 2, format = "f", width = 5), sep = "  ",
                  "\n")
                  cat("\n")
                  cat("selection percentage at each dose combination (%):\n")
                  for (i in 1:dim(object$p.true)[1]) cat(formatC(object$selpercent[i,
                  ], digits = 2, format = "f", width = 5), sep = "  ",
                  "\n")

                  cat("\n")
                  cat("average number of patients treated at each dose combination:\n")
                  for (i in 1:dim(object$p.true)[1]) cat(formatC(object$npatients[i,
                  ], digits = 2, format = "f", width = 5), sep = "  ",
                  "\n")
                  cat("\n")
                  cat("average number of toxicity observed at each dose combination:\n")
                  for (i in 1:dim(object$p.true)[1]) cat(formatC(object$ntox[i,
                  ], digits = 2, format = "f", width = 5), sep = "  ",
                  "\n")
                  cat("\n")
                  cat("average DLT percentage per dose:", formatC(object$perct_DLT_per_dose, digits=2, 
                                                                  format="f"), "\n")
                  cat("\n")
                  cat("average percentage of patients per dose:", formatC(apply(object$perct_pts_per_dose,2,mean), digits=4, 
                                                                          format="f"), "\n")
                  cat("\n")
                  cat("average number of toxicities:", formatC(object$totaltox,
                                                               digits = 1, format = "f"), "\n")
                  cat("\n")
                  cat("average number of patients:", formatC(object$totaln,
                                                             digits = 1, format = "f"), "\n")
                  cat("\n")
                  cat("average trial_duration(month):", formatC(object$average_trial_duration, digits=1, format="f"), "\n")
                  cat("\n")
                  
                  cat("25% quantile trial_duration(month):", formatC(object$twentyFive_Perct_trial_duration, digits=2, format="f"), "\n")
                  cat("\n")
                  cat("75% quantile trial_duration(month):", formatC(object$senventyFive_Perct_trial_duration, digits=2, format="f"), "\n")
                  cat("\n")
                  
                  cat("risk of poor allocation:", formatC(object$risk_poor_allocation, digits=2, format="f"), "\n")
                  cat("\n")
                  cat("risk of high toxicity:", formatC(object$risk_high_toxicity, digits=2, format="f"), "\n")
                  car("\n")
                  cat("percentage of patients treated at MTD contour:",
                      object$npercent.contour, "\n")
                  cat("percentage of patients treated above MTD contour:",
                      formatC(object$npercent.above.contour, digits = 1,
                              format = "f"), "\n")
                  cat("percentage of patients treated below MTD contour:",
                      formatC(object$npercent.below.contour, digits = 1,
                              format = "f"), "\n")
                  cat("percentage of correct selection of the MTD contour:",
                      formatC(object$pcs.contour, digits = 1, format = "f"),
                      "\n")
            }
            else {
                  cat("true DLT rate of dose combinations:\n")
                  for (i in 1:dim(object$p.true)[1]) {
                        cat(formatC(object$p.true[i, ], digits = 2, format = "f",
                                    width = 5), sep = "  ", "\n")
                  }
                  cat("\n")
                  cat("selection percentage at each dose combination (%):\n")
                  for (i in 1:dim(object$p.true)[1]) {
                        cat(formatC(object$selpercent[i, ], digits = 2,
                                    format = "f", width = 5), sep = "  ", "\n")
                  }
                  cat("\n")
                  

                  cat("average number of patients treated at each dose combination:\n")
                  for (i in 1:dim(object$p.true)[1]) {
                        cat(formatC(object$npatients[i, ], digits = 2,
                                    format = "f", width = 5), sep = "  ", "\n")
                  }
                  cat("\n")
                  cat("average number of toxicity observed at each dose combination:\n")
                  for (i in 1:dim(object$p.true)[1]) {
                        cat(formatC(object$ntox[i, ], digits = 2, format = "f",
                                    width = 5), sep = "  ", "\n")
                  }
                  cat("average DLT percentage per dose:", formatC(object$perct_DLT_per_dose, digits=2, 
                                                                  format="f"), "\n")
                  cat("\n")
                  cat("average percentage of patients per dose:", formatC(apply(object$perct_pts_per_dose,3,mean), digits=4, 
                                                                          format="f"), "\n")
                  cat("\n")
                  cat("average trial_duration(month):", formatC(object$average_trial_duration, digits=1, format="f"), "\n")
                  cat("\n")
                  cat("25% quantile trial_duration(month):", formatC(object$twentyFive_Perct_trial_duration, digits=2, format="f"), "\n")
                  cat("\n")
                  cat("75% quantile trial_duration(month):", formatC(object$senventyFive_Perct_trial_duration, digits=2, format="f"), "\n")
                  cat("\n")
                  
                  cat("risk of poor allocation:", formatC(object$risk_poor_allocation, digits=2, format="f"), "\n")
                  cat("\n")
                  cat("risk of high toxicity:", formatC(object$risk_high_toxicity, digits=2, format="f"), "\n")
                  car("\n")
                  
                  cat("average number of toxicities:", formatC(object$totaltox,
                                                               digits = 1, format = "f"), "\n")
                  cat("average number of patients:", formatC(object$totaln,
                                                             digits = 1, format = "f"), "\n")
                  # cat("selection percentage of MTD:", formatC(object$pcs,
                  #                                             digits = 1, format = "f"), "\n")
                  # cat("percentage of patients treated at MTD:", formatC(object$npercent,
                  #                                                       digits = 1, format = "f"), "\n")
                  cat("percentage of early stopping due to toxicity:",
                      formatC(object$percentstop, digits = 1, format = "f"),
                      "% \n")
            }
      }
      
      #create a data frame for csv output
      df_title = data.frame(information_name = c("true DLT rate: ",
                                                 "selection percentage at each dose level (%): ", 
                                                 "average number of patients treated at each dose level: ",
                                                 "average DLT percentage per dose: ",
                                                 "average percentage of patients per dose: ",
                                                 "average number of DLT: ",
                                                 "average number of patients: ",
                                                 "average trial_duration(month): ",
                                                 "25% quantile trial_duration(month): ",
                                                 "75% quantile trial_duration(month): ",
                                                 "risk of poor allocation: ",
                                                 "risk of high toxicity: "))
      df_content = data.frame(content =
                                    c(paste(formatC(object$p.true, digits=2, format="f"),collapse = ","),
                                      paste(formatC(object$selpercent, digits=1, format="f"),collapse = ","),
                                      paste(formatC(object$npatients, digits=1, format="f"),collapse = ","),
                                      
                                      paste(formatC(object$perct_DLT_per_dose, digits=2, format="f"),collapse = ","),
                                      paste(formatC(apply(object$perct_pts_per_dose,2,mean), digits=2, format="f"),collapse = ","),
                                      paste(formatC(object$totaltox, digits=1, format="f"),collapse = ","),
                                      paste(formatC(object$totaln, digits=1, format="f"),collapse = ","),
                                      paste(formatC(object$average_trial_duration, digits=1, format="f"),collapse = ","),
                                      paste(formatC(object$twentyFive_Perct_trial_duration, digits=2, format="f"),collapse = ","),
                                      paste(formatC(object$senventyFive_Perct_trial_duration, digits=2, format="f"),collapse = ","),
                                      paste(formatC(object$risk_poor_allocation, digits=2, format="f"),collapse = ","),
                                      paste(formatC(object$risk_high_toxicity, digits=2, format="f"),collapse = ",")))
      
      result = cbind(df_title, df_content)
      
      write.table(result, file_dir, na = "NA", append = TRUE, 
                  col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
      write.table('\n', file_dir, append = TRUE, 
                  col.names = FALSE, row.names = FALSE, sep = ",", quote = FALSE)    
}

################################## Calculate trial duration #####################################

calculate_duration <-function(newcohort, newcohort_dlt, monthly_pts_enroll,follow_up_time ){
      
      
      start_date = c()
      #calculate enrollment duration per cohort
      #if number of patients enrolled per month smaller than cohort size
      len_cohort = length(newcohort)
      if(monthly_pts_enroll < len_cohort){
            num_follow_up_times_needed = ceiling(length(newcohort)/monthly_pts_enroll)
            time_period_lower_boundry = 1
            time_period_higher_boundry = follow_up_time
            
            #allocate all cohort patients to different month
            for(i in 1: num_follow_up_times_needed){
                  if(i == num_follow_up_times_needed){
                        num_pts_left = length(newcohort) - (i-1)*monthly_pts_enroll
                        #generate start date for pts in last follow_up_time
                        startdate = floor(runif(num_pts_left, 
                                                min=time_period_lower_boundry, max=time_period_higher_boundry))
                  }
                  else{
                        #generate start date for pts in each follow_up_time
                        startdate = floor(runif(monthly_pts_enroll, 
                                                min=time_period_lower_boundry, max=time_period_higher_boundry))
                  }
                  start_date = c(start_date, startdate)
                  time_period_lower_boundry = time_period_higher_boundry+1
                  time_period_higher_boundry = time_period_higher_boundry + follow_up_time
                  
            }
            #if number of patients enrolled per month larger or equal cohort size
      }else{
            # we can only add cohort_size pts for current dose
            start_date = floor(runif(len_cohort, min=1, max=follow_up_time))
      }
      
      
      #generate end_date and calculate current cohort duration
      end_date = rep(0, length(newcohort))
      new_cohort_duration = 0
      
      #if no DLT in current cohort
      if(newcohort_dlt==0){
            end_date = start_date + follow_up_time
            
      }
      #if >=1 DLT in current cohort
      else{
            for (i in 1: cohortsize){
                  if(newcohort[i] ==1){
                        #generate patients follow up time within follow_up_time
                        dlt_follow_up_time = floor(runif(1, min=1, max=follow_up_time))
                        end_date[i] = start_date[i] + dlt_follow_up_time
                  }else{
                        end_date[i] = start_date[i] + follow_up_time}
            }
            
      }
      
      #calculate the trail duration for current cohort
      new_cohort_duration = max(end_date) - start_date[which.min(start_date)]
      return(new_cohort_duration)
      
}








################################### Simulation #################################################
#no gap between each cohort
#number of pats enrolled each month >=1
boin_simulation <-
      function (target,
                p.true,
                ncohort,
                cohortsize,
                n.earlystop = 100,
                startdose = 1,
                titration = FALSE,
                p.saf = 0.6 * target,
                p.tox = 1.4 *
                      target,
                cutoff.eli = 0.95,
                extrasafe = FALSE,
                offset = 0.05,
                boundMTD = FALSE,
                ntrial = 1000,
                follow_up_time = 30, 
                monthly_pts_enroll = 3,
                seed = 6)
      {
            if (target < 0.05) {
                  stop("the target is too low!")
                  
            }
            if (target > 0.6) {
                  stop("the target is too high!")
                  
            }
            if ((target - p.saf) < (0.1 * target)) {
                  stop("the probability deemed safe cannot be higher than or too close to the target!")
            }
            if ((p.tox - target) < (0.1 * target)) {
                  stop("the probability deemed toxic cannot be lower than or too close to the target!")
            }
            if (offset >= 0.5) {
                  stop("the offset is too large!")
            }
            if (n.earlystop <= 6) {
                  warning(
                        "the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18."
                  )
            }
            set.seed(seed)
            if (cohortsize == 1)
                  titration = FALSE
            lambda_e = log((1 - p.saf) / (1 - target)) / log(target * (1 - p.saf) /
                                                            (p.saf * (1 - target)))
            lambda_d = log((1 - target) / (1 - p.tox)) / log(p.tox * (1 - target) /
                                                            (target * (1 - p.tox)))
            ndose = length(p.true)
            npts = ncohort * cohortsize
            Y = matrix(rep(0, ndose * ntrial), ncol = ndose)
            N = matrix(rep(0, ndose * ntrial), ncol = ndose)
            dselect = rep(0, ntrial)
            
            
            if (cohortsize > 1) {
                  temp = get.boundary(
                        target,
                        ncohort,
                        cohortsize,
                        n.earlystop = ncohort * cohortsize,
                        p.saf,
                        p.tox,
                        cutoff.eli,
                        extrasafe
                  )$full_boundary_tab
            }
            #cohortsize =1, the $full_boundary_tab is NULL so we need $boundary_tab
            else {
                  temp = get.boundary(
                        target,
                        ncohort,
                        cohortsize,
                        n.earlystop = ncohort * cohortsize,
                        p.saf,
                        p.tox,
                        cutoff.eli,
                        extrasafe
                  )$boundary_tab
            }
            #escalate dose level
            b.e = temp[2,]
            #de-escalate dose level
            b.d = temp[3,]
            #eliminate dose level
            b.elim = temp[4,]
            all_trials_duration = rep(0, ntrial)
            
            check_toxicity = 0
            pat_in_mtd = rep(0, ntrial)
            
            #run simulation
            for (trial in 1:ntrial) {
                  y <- rep(0, ndose)
                  n <- rep(0, ndose)
                  earlystop = 0
                  d = startdose
                  elimi = rep(0, ndose)
                  ft = TRUE #flag used to determine whether or not to add cohortsize-1 patients to a dose for the first time when titration is triggered.
                  trial_duration = 0
                 
                  #use cohort size =1 for initial dose escalation
                  if (titration) {
                        #random generate number of ndose(6) numbers from 0~1
                        #if number < p.true, means that dose has DLT
                        z <- (runif(ndose) < p.true)
                        #if all random generated num > p.true, 
                        #no DLT,assign each dose with 1 patient
                        #d = last dose
                        if (sum(z) == 0) {
                              d = ndose
                              n[1:ndose] = 1
                        }
                        #else: find the first num < p.true, 
                        #assign each dose with 1 pats before num
                        #d = first dose with dlt
                        else {
                              d = which(z == 1)[1]
                              n[1:d] = 1
                              y[d] = 1
                        }
                        #calculate trial duration time
                        startDate = floor(runif(1, 
                                                min=1, max=follow_up_time))
                        duration = startDate + follow_up_time*(d-1)
                        trial_duration = trial_duration + duration
                        #dlt pat need dlt_follow_up time
                        if(d != ndose){
                              dlt_follow_up_time = floor(runif(1, min=1, max=follow_up_time))
                              trial_duration = trial_duration + dlt_follow_up_time
                        }
                  }
                  #for each cohort
                  for (i in 1:ncohort) {
                        #if use titration and does d has 1 pats
                        #and < cohort size
                        if (titration && n[d] < cohortsize && ft) {
                              ft = FALSE
                              #add toxicity pats in left cohortsize pts to y[d]
                              newcohort = runif(cohortsize - 1) < p.true[d]
                              newcohort_dlt = sum(newcohort)
                              
                              y[d] = y[d] + newcohort_dlt
                              n[d] = n[d] + cohortsize - 1
                        }
                        else {
                              #generate newcohort number of pts with/without dlt
                              newcohort = runif(cohortsize) < p.true[d]
                              
                              #if total pats > npts, only add part of cohort pts
                              if ((sum(n) + cohortsize) >= npts) {
                                    nremain = npts - sum(n)
                                    
                                    newcohort = newcohort[1:nremain]
                                    newcohort_dlt = sum(newcohort)
                                    
                                    y[d] = y[d] + newcohort_dlt
                                    
                                    n[d] = n[d] + nremain
                                    
                                    break
                                    
                              }
                              else{
                                    newcohort_dlt = sum(newcohort)
                                    y[d] = y[d] + newcohort_dlt
                                    
                                    n[d] = n[d] + cohortsize
                                    
                              }
                              
                             
                              
                        }
                        #calculate enrollment duration per cohort
                        new_cohort_duration = calculate_duration(newcohort,newcohort_dlt, monthly_pts_enroll,follow_up_time )
                        trial_duration = trial_duration + new_cohort_duration
                        
                        #check if dose should be eliminated
                        if (!is.na(b.elim[n[d]])) {
                              if (y[d] >= b.elim[n[d]]) {
                                    elimi[d:ndose] = 1
                                    #if dose one is too toxic
                                    if (d == 1) {
                                          earlystop = 1
                                          break
                                    }
                              }
                              #if extrasafe is used, apply for stringent stopping rule
                              if (extrasafe) {
                                    #if dose level is 1 and num of pats >=3
                                    if (d == 1 && n[1] >= 3) {
                                          if (1 - pbeta(target, y[1] + 1, n[1] - y[1] +
                                                        1) > cutoff.eli - offset) {
                                                earlystop = 1
                                                break
                                          }
                                    }
                              }
                        }
                        
                        #if num of pats in current dose > n.earlystop
                        if (n[d] >= n.earlystop &&
                            #keep current dose level
                            ((y[d] > b.e[n[d]] && y[d] < b.d[n[d]]) ||
                             #dose level 1 and need de-escalate
                             (d == 1 && y[d] >= b.d[n[d]]) ||
                             #dose level ndose or next dose is eliminated 
                             #and need escalate
                             ((d == ndose ||
                               elimi[d + 1] == 1) && y[d] <= b.e[n[d]])))
                              break
                        
                        #if decision is to escalate dose level
                        if (y[d] <= b.e[n[d]] && d != ndose) {
                              if (elimi[d + 1] == 0)
                                    d = d + 1
                        }
                        #if decision is to de-escalate dose level
                        else if (y[d] >= b.d[n[d]] && d != 1) {
                              d = d - 1
                        }
                        #keep dose level
                        else {
                              d = d
                        }
                        
                  }
                  #####end trial
        
                  #select mtd
                  Y[trial,] = y
                  N[trial,] = n
                  if (earlystop == 1) {
                        dselect[trial] = 99
                  }
                  else {
                        dselect[trial] = select.mtd(
                              target,
                              n,
                              y,
                              cutoff.eli,
                              extrasafe,
                              offset,
                              boundMTD = boundMTD,
                              p.tox = p.tox
                        )$MTD
                  }
                  all_trials_duration[trial] = trial_duration
                  
                  #number of patients allocated to MTD
                  if(dselect[trial] != 0){
                        pat_in_mtd[trial] = n[dselect[trial]]
                  }
                  
                  #total number of toxicity > number of patients* target
                  if(sum(y) > ncohort*cohortsize*target){
                        check_toxicity = check_toxicity+1
                  }
            }
            
            #calculate risk of poor allocation: percetage of pats allocated to MTD < that of standard non-sequential design
            pats_non_seq = ncohort*cohortsize/ndose
            risk_poor_allocation = sum(pat_in_mtd < pats_non_seq, na.rm=T )/ntrial
            if(is.na(risk_poor_allocation)){
                  print(sum(pat_in_mtd < pats_non_seq ))
            }
            #percentage of total number of toxicity > num_of_pts * target rate
            risk_high_toxicity = check_toxicity/ntrial
            
            #selection percentage at each dose level
            selpercent = rep(0, ndose)
            #average dlt percentage per dose
            perct_DLT_per_dose = apply(Y/N,2,mean,na.rm=TRUE)
            #average percentage of patients per dose
            perct_pts_per_dose = t(apply(N,1, function(x) x/sum(x)))
            #average number of pats per dose
            nptsdose = apply(N, 2, mean)
            #average number of toxicity observed per dose
            ntoxdose = apply(Y, 2, mean)
            for (i in 1:ndose) {
                  selpercent[i] = sum(dselect == i) / ntrial * 100
            }
            #average trial_duration(month)
            average_trial_duration = sum(all_trials_duration)/ntrial/30.4375
            # 25% quantile trail_duration
            twentyFive_Perct_trial_duration = quantile(all_trials_duration, prob = 0.25)/30.4375
            # 75% quantile trail_duration
            senventyFive_Perct_trial_duration = quantile(all_trials_duration, prob = 0.75)/30.4375
            
            
            #if target is within the p.true, print overdosing60 and overdosing80
            if (length(which(p.true == target)) > 0) {
                  if (which(p.true == target) == ndose - 1) {
                        overdosing60 = mean(N[, p.true > target] > 0.6 *
                                                  npts) * 100
                        overdosing80 = mean(N[, p.true > target] > 0.8 *
                                                  npts) * 100
                  }
                  else {
                        overdosing60 = mean(rowSums(N[, p.true > target]) >
                                                  0.6 * npts) * 100
                        overdosing80 = mean(rowSums(N[, p.true > target]) >
                                                  0.8 * npts) * 100
                  }
                  out = list(
                        p.true = p.true,
                        selpercent = selpercent,
                        perct_DLT_per_dose=perct_DLT_per_dose,
                        perct_pts_per_dose=perct_pts_per_dose,
                        average_trial_duration=average_trial_duration,
                        twentyFive_Perct_trial_duration=twentyFive_Perct_trial_duration,
                        senventyFive_Perct_trial_duration=senventyFive_Perct_trial_duration,
                        risk_poor_allocation=risk_poor_allocation,
                        risk_high_toxicity=risk_high_toxicity,
                        npatients = nptsdose,
                        ntox = ntoxdose,
                        totaltox = sum(Y) / ntrial,
                        totaln = sum(N) / ntrial,
                        percentstop = sum(dselect == 99) / ntrial * 100,
                        overdose60 = overdosing60,
                        overdose80 = overdosing80,
                        simu.setup = data.frame(
                              target = target,
                              p.true = p.true,
                              ncohort = ncohort,
                              cohortsize = cohortsize,
                              startdose = startdose,
                              p.saf = p.saf,
                              p.tox = p.tox,
                              cutoff.eli = cutoff.eli,
                              extrasafe = extrasafe,
                              offset = offset,
                              ntrial = ntrial,
                              dose = 1:ndose
                        ),
                        flowchart = TRUE,
                        lambda_e = lambda_e,
                        lambda_d = lambda_d
                  )
            }
          
            #if the target is not in the p.true, no overdosing60 and overdosing80
            else {
                  out = list(
                        p.true=p.true,
                        selpercent = selpercent,
                        perct_DLT_per_dose=perct_DLT_per_dose,
                        perct_pts_per_dose=perct_pts_per_dose,
                        average_trial_duration=average_trial_duration,
                        twentyFive_Perct_trial_duration = twentyFive_Perct_trial_duration,
                        senventyFive_Perct_trial_duration=senventyFive_Perct_trial_duration,
                        risk_poor_allocation=risk_poor_allocation,
                        risk_high_toxicity=risk_high_toxicity,
                        npatients = nptsdose,
                        ntox = ntoxdose,
                        totaltox = sum(Y) / ntrial,
                        totaln = sum(N) / ntrial,
                        percentstop = sum(dselect == 99) / ntrial * 100,
                        simu.setup = data.frame(
                              target = target,
                              p.true = p.true,
                              ncohort = ncohort,
                              cohortsize = cohortsize,
                              startdose = startdose,
                              p.saf = p.saf,
                              p.tox = p.tox,
                              cutoff.eli = cutoff.eli,
                              extrasafe = extrasafe,
                              offset = offset,
                              ntrial = ntrial,
                              dose = 1:ndose
                        ),
                        flowchart = TRUE,
                        lambda_e = lambda_e,
                        lambda_d = lambda_d
                  )
            }
            class(out) <- "boin"
            return(out)
      }





























########### Run Simulation with different scenario######################

boin_summary<-function(scenarios, file_dir, 
                      target,
                       ncohort,
                       cohortsize,
                       n.earlystop = 100,
                       startdose = 1,
                       titration = FALSE,
                       p.saf = 0.6 * target,
                       p.tox = 1.4 *target,
                       cutoff.eli = 0.95,
                       extrasafe = FALSE,
                       offset = 0.05,
                       boundMTD = FALSE,
                       ntrial = 1000,
                       follow_up_time = 30, 
                       monthly_pts_enroll = 3,
                       seed = 6){
      
      parameters = c("target", "ncohort", "cohortsize","n.earlystop",
                     "startdose", "titration", "p.saf", "p.tox",
                     "cutoff.eli","extrasafe","offset","boundMTD",
                     "ntrial", "follow_up_time","monthly_pts_enroll","seed")
      
      meaning = c(
                  "the target toxicity rate",
                  "the total number of times to add patients",
                  "number of patients in each cohort",
                  "the max number of patients assigned to each dose",
                  "The starting dose level for treating the first cohort of patients",
                  "If titration=TRUE dose titration is performed to accelerate dose escalation",
                  "The highest toxicity probability that is deemed sub therapeutic",
                  "The lowest toxicity probability that is deemed overly toxic",
                  "The cutoff to eliminate the overly toxic dose for safety",
                  "Set extrasafe=TRUE to impose a stricter stopping rule",
                  "A small positive number(between0and0.5) to control how strict the stopping rule is",
                  "set boundMTD=TRUE to impose the special condition",
                  "The number of trials to be simulated",
                  "follow up time for non DLT patients",
                  "number of patients enrolled per month",
                  "Set a seed for the random number generator")
      values = c(target, ncohort, cohortsize,n.earlystop,
                 startdose, titration, p.saf, p.tox,
                 cutoff.eli, extrasafe, offset, boundMTD,
                 ntrial, follow_up_time, monthly_pts_enroll, seed)
      df = data.frame(parameters, meaning, values)
      
      write.table(df, file_dir, na = "NA", append = TRUE, 
                  col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
      write.table('\n', file_dir, append = TRUE, 
                  col.names = FALSE, row.names = FALSE, sep = ",", quote = FALSE)
      
      #scenarios: 
      for(i in 1: ncol(scenarios)){
            true_toxic = scenarios[, i]
            sim1 = boin_simulation(
                  target = target,
                  p.true = true_toxic,
                  ncohort = ncohort,
                  cohortsize = cohortsize,
                  n.earlystop = n.earlystop,
                  startdose = startdose,
                  titration = titration,
                  p.saf = 0.6 * target,
                  p.tox = 1.4 * target,
                  cutoff.eli = cutoff.eli,
                  extrasafe = extrasafe,
                  offset = offset,
                  boundMTD = boundMTD,
                  ntrial = ntrial,
                  follow_up_time = follow_up_time, 
                  monthly_pts_enroll = monthly_pts_enroll,
                  seed = seed
            )
            summary.boin_csv(file_dir = file_dir, sim1)
      }
      
}
