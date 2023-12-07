library("BOIN")
source("boin_helper_func_1.R")
#BOIN Package
##################get decision boundary and decision table ################################################
#     file_dir:   where to save the decision boundary & table
#    
#     target:     The target DLT rate
#     
#     ncohort:    The total number of cohorts
#     
#     cohortsize: The cohort size
#     
#     n.earlystop:the early stopiing parameter. If the number of patients treated at the current 
#                 dose reaches n.earlystop, stop the trial early and select the
#                 MTD based on the observed data. The default valu of is 100 essentially
#                 turns off this type of early stopping
#     
#     cutoff.eli: The cutoff to eliminate the overly toxic dose for safety.
#                 We recommend the default value cutoff.eli=0.95 for general use.
#     
#     extrasafe:  Set extrasafe=TRUE to impose a stricter stopping rule.
#     
#     offset:     A small positive number(between0and0.5) to control how strict the stopping 
#                 rule is when extrasafe=TRUE.A larger value leads to a stricter stopping rule.
#                 The default value offset=0.05 generally works well.
##############################################################################
target = 0.3
ncohort = 12
cohortsize = 3
n.earlystop=12
cutoff.eli=0.95
extrasafe=FALSE
offset=0.05
file_dir = "boin_decision_boundary.csv"
get_boin_decision(file_dir = file_dir, target = target,
                               ncohort = ncohort,
             cohortsize=cohortsize,
             n.earlystop=n.earlystop, 
             cutoff.eli=cutoff.eli, 
             extrasafe=extrasafe,
             offset=offset)




##################Run Simulation################################################
#     scenarios:  different P.true 
#  
#     file_dir:   where to save the simulation results
#
#     target:     The target DLT rate
#     
#     ncohort:    The total number of cohorts
#     
#     cohortsize: The cohort size
#     
#     n.earlystop:he early stopiing parameter. If the number of patients treated at the current 
#                 dose reaches n.earlystop, stop the trial early and select the
#                 MTD based on the observed data. The default valu of is 100 essentially
#                 turns off this type of early stopping
#     
#     startdose:  The starting dose level for treating the first cohort of patients.
#                 The default value is start dose=1,i.e.,starting from the lowest dose.
#     
#     titration:  If titration=TRUE,dose titration is performed to accelerate dose escalation at the
#                 beginning of the trial,where patients are treated one by one(i.e.,cohort size=1),
#                 starting from start dose.If no DLT is observed,escalate the dose;otherwise 
#                 switch to the specified cohortsize = cohortsize
#     
#     p.saf:      The highest toxicity probability that is deemed sub therapeutic(i.e.,below the MTD)
#                 such that dose escalation should be made.The default value of p.saf= 0.6*target.
#     
#     p.tox:      The lowest toxicity probability that is deemed overly toxic such that 
#                 dose de-escalation is required.
#                 The default value of p.tox=1.4*target.
#     
#     cutoff.eli: The cutoff to eliminate the overly toxic dose for safety.
#                 We recommend the default value cutoff.eli=0.95 for general use.
#     
#     extrasafe:  Set extrasafe=TRUE to impose a stricter stopping rule.
#    
#     offset:     A small positive number(between0and0.5) to control how strict the stopping 
#                 rule is when extrasafe=TRUE.A larger value leads to a stricter stopping rule.
#                 The default value offset=0.05 generally works well.
#
#     boundMTD:   set boundMTD=TRUE to impose the condition: 
#                 the isotonic estimate of toxicity
#                 probability for the selected MTD must be less than de-escalation boundary.
#
#     ntrial:     The number of trials to be simulated.
#
#     follow_up_time : follow up time for non DLT patients
#
#     monthly_pts_enroll: number of patients enrolled per month
#
#     seed:       Set a seed for the random number generator.
##############################################################################
######################################################################################
target = 0.3
ncohort = 12
cohortsize = 3
n.earlystop = 18
startdose = 1
titration = FALSE
cutoff.eli = 0.95
extrasafe = FALSE
offset = 0.05
boundMTD = FALSE
ntrial = 1000
follow_up_time = 30
monthly_pts_enroll = 3
seed = 6
file_dir = "boin_simulation_result.csv"
#####Scenarios###################
scenarios = data.frame(s1 = c(0.3,0.46,0.50,0.54,0.58),
                       s2 = c(0.16,0.3,0.47,0.54,0.60),
                       s3 = c(0.04,0.15,0.3,0.48,0.68),
                       s4 = c(0.02,0.07,0.12,0.3,0.45),
                       s5 = c(0.02,0.06,0.1,0.13,0.3))
boin_summary(
      scenarios=scenarios, file_dir = file_dir, target = target,
      ncohort = ncohort, cohortsize = cohortsize, n.earlystop = n.earlystop,
      startdose = startdose, titration = titration, p.saf = 0.6 * target,
      p.tox = 1.4 * target, cutoff.eli = cutoff.eli, extrasafe = extrasafe,
      offset = offset, boundMTD = boundMTD, ntrial = ntrial,
      follow_up_time = follow_up_time, monthly_pts_enroll = monthly_pts_enroll,
      seed = seed
)



