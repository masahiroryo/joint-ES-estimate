#####################################################################
#
#  Rillig et al.(20XX) 
#
#  Appendix: joint effect size estimation and visualization
#
#
#####################################################################


#####################################################################
#####                    Library                                #####
#####################################################################
library(ggplot2)
library(ggepi)
library(ggridges)
library(patchwork)

#devtools::install_github("lwjohnst86/ggepi")
#devtools::install_github("thomasp85/patchwork")
#ggridges: https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html


#####################################################################
#####                    Data info                              #####
#####################################################################

df = data.frame(read.csv("../data/treatment_NA_interpolated.csv", sep=","))

levels = c("Level 2", "Level 5", "Level 8", "Level 10")
stressors = c("drought", "Ndep", "temp", "microplastic", "glyphosate", 
              "antibiotics", "fungicide", "copper", "salinity", "insecticide")
treatment = unique(df$remark)

responses = c("WSA", "wdpt.sec.", "Percent.decomposition..corr.", "CO2_ppm_ThirdWeek",
              "Rich.ASV", "Comm.Comp1", "Comm.Disp")

n_iter   = 10000      # permutaion size

response = responses[1]



#####################################################################
#####                    Functions                              #####
#####################################################################

#----------------------------------------------
# Estimating mean and its 95% confidence interval
#----------------------------------------------

BootStrap_mean = function(response, data=df, target = treatment, n_perm = n_iter){
  summary = list()

  for(treatment in target){
      bs = numeric(0)
      population = df[df$remark==treatment, response]
      size = length(population)
      
      for(id in c(1:n_perm)){
        k = mean(sample(population, size, replace = T))
        bs = append(bs, k)
      }
      summary[[treatment]] = c(quantile(bs, .025), mean(bs), quantile(bs, .975))
      names(summary[[treatment]]) = c("2.5%", "mean", "97.5%")
  }
  summary = t(data.frame(summary))
  summary = data.frame("target" = target, summary); row.names(summary) = c()
  return(summary)
}


#----------------------------------------------
# Estimating unstandardized effect size and its 95% confidence interval
#----------------------------------------------

BootStrap_ES_rep = function(response, data=df, target = treatment, n_perm = n_iter){
  resampled = list()

  for(treatment in target){
    bs = numeric(0)
    population_CT = df[df$remark=="control", response]
    population_TR = df[df$remark==treatment, response]
    size_CT = length(population_CT)
    size_TR = length(population_TR)
    
    for(id in c(1:n_perm)){
      k_CT = mean(sample(population_CT, size_CT, replace = T))
      k_TR = mean(sample(population_TR, size_TR, replace = T))
      bs = append(bs, k_TR - k_CT)
    }
    resampled[[treatment]] = bs
  }
  resampled[["control"]] = rep(0, n_perm)
  return(resampled)
}

BootStrap_ES_summary = function(data){
  summary = list()
  p = 0
  summary[["control"]] = c(0,0,0,1)
  target = names(data)
  
  for(treatment in target[-1]){
    bs = data[[treatment]]
    p = length(which(bs>0))/length(bs)
    p = min(p, 1-p)
    summary[[treatment]] = c(quantile(bs, .025), mean(bs), quantile(bs, .975), p)
  }
  summary = t(data.frame(summary))
  colnames(summary) = c("2.5%", "mean", "97.5%", "p_value")
  color_p = rep("#999999", nrow(summary))
  for (s in c(0.05,0.01,0.001)) color_p[which(summary[,"p_value"]<s)] = c("#FFD35C", "#FF8201", "#FF0000")[which(s==c(0.05,0.01,0.001))]
  summary = data.frame(target, summary, color_p); row.names(summary) = c()
  
  return(summary)
}

  
#-------------------------------------------
#  Estimating unstandardized ES of joint stressors and its 95% confidence interval
#-------------------------------------------

Null_distribution_rep = function(response, data=df, n_perm=n_iter){

  output = list()
  for(Lv in levels){
    
    resampled = list()
    
    # Checking which stressor combinations were jointly tested
    combination = data[data["remark"]==Lv,c(2:11)]
    level = sum(combination[1,])
    
    # Null distributions can be taken based on three different assumptions
    for(type in c("Additive", "Multiplicative", "Dominative")){
      bs = numeric(0)
      population_CT = df[df$remark=="control", response]
      size_CT = length(population_CT)
      
      # For each combination, bootstrap resampling is conducted
      for(j in c(1:nrow(combination))){
        selected_stressors = stressors[which(combination[j,]==1)]
        sub_n_perm = ceiling(n_perm/nrow(combination))
        
        # bootstrap resampling
        for(id in c(1:sub_n_perm)){
          each_effect = numeric(0)
          k_CT = mean(sample(population_CT, size_CT, replace = T))
          
          for(treatment in selected_stressors){
            population_TR = df[df$remark==treatment, response]
            size_TR = length(population_TR)
            k_TR = mean(sample(population_TR, size_TR, replace = T))
            
            # ES estimate depending on the type of null hypotheses
            if(type=="Additive")       each_effect = append(each_effect, (k_TR - k_CT))
            if(type=="Multiplicative") each_effect = append(each_effect, (k_TR - k_CT)/k_CT)
            if(type=="Dominative")      each_effect = append(each_effect, (k_TR - k_CT))
          }
          
          # Calculating an expected ES after collecting the ESs of all relevant single stressors
          if(type=="Additive")       joint_effect = sum(each_effect)
          if(type=="Multiplicative"){
            z = 1
            for(m in c(1:level)) z = z * (1 + each_effect[m])
            joint_effect = (z - 1)*k_CT
          }
          if(type=="Dominative")      joint_effect = each_effect[which(max(abs(each_effect))==abs(each_effect))]
  
          bs = append(bs, joint_effect)
        }
      }
      resampled[[type]] = bs[1:n_perm]
    }
    output[[Lv]] = resampled
  }  
  return(output)
} 

NHST_summary = function(null_data, Actual_data){
  output = list()
  for(Lv in levels){
    summary = list()
    summary[["Actual"]] = c(quantile(Actual_data[[Lv]], .025), mean(Actual_data[[Lv]]), quantile(Actual_data[[Lv]], .975), 1)
    p = 0
    assumptions = names(null_data[[Lv]])
    
    for(i_assumption in assumptions){
      bs   = (Actual_data[[Lv]] - null_data[[Lv]][[i_assumption]])
      p = length(which(bs>0))/length(bs)
      p = min(p, 1-p)
      summary[[i_assumption]] = c(quantile(null_data[[Lv]][[i_assumption]], .025), mean(null_data[[Lv]][[i_assumption]]), quantile(null_data[[Lv]][[i_assumption]], .975), p)
    }
    summary = t(data.frame(summary))
    colnames(summary) = c("2.5%", "mean", "97.5%", "p_value")
    summary = data.frame(ES = c("Actual", "Additive","Multiplicative","Dominative"), summary); row.names(summary) = c()
    
    output[[Lv]] = summary
  }
  
  return(output)
}


NHST_summary_transform = function(data){
  output = list()
  for(i in 1:4){
    summary = rbind(data[["Level 2"]][i, 2:4], data[["Level 5"]][i, 2:4],
                    data[["Level 8"]][i, 2:4], data[["Level 10"]][i, 2:4])
    summary = cbind(levels, summary)
    colnames(summary) = c("Lv", "Low", "Mean", "High")
    output[[c("Actual", "Additive", "Multiplicative", "Dominative")[i]]] = summary
  }
  return(output)
}
    
#####################################################################
#####                  Main analysis                            #####
#####################################################################

# Bootstrap estimate: single stressors
response_mean = BootStrap_mean(response)
response_ES_bs  = BootStrap_ES_rep(response)
response_ES     = BootStrap_ES_summary(response_ES_bs)

# Bootstrap estimate: joint stressors
Null_ES_bs      = Null_distribution_rep(response)
joint_ES_null   = NHST_summary(Null_ES_bs, response_ES_bs)
ES_plot         = NHST_summary_transform(joint_ES_null)

# Output as table
write.csv(response_mean,   "mean estimate for each condition.csv")
write.csv(response_ES[,1:5], "Unstandardized effect size (difference) for each treatment.csv")
write.csv(joint_ES_null, "Test if joint effect size can be predicted.csv")


# Ploting the raw data with the mean & its 95% confidence intervals
g_rawdata = 
ggplot() +
  ylab("Treatment") +
  xlab(response) + 
  scale_y_discrete(limits = rev(treatment)) +
  stat_density_ridges(data=df, aes_string(x = response, y = "remark"),
                      geom = "density_ridges_gradient",
                      rel_min_height = 0.01, 
                      jittered_points = TRUE, color="#00000000",
                      alpha = .5,
                      position = position_points_jitter(height = .1, yoffset = .15),
                      point_size = 1, point_alpha = .3, 
                      scale = .5) +
  theme(legend.position = 'none') +
  geom_estci(data=response_mean, aes(x = mean, y = target, xmin=X2.5., xmax=X97.5., 
                                       xintercept=response_mean[1,"mean"]), center.linecolour = "black",
             colour = response_ES$color_p,
             size=0.7, ci.linesize = 0.9, position=position_nudge(y = -0.15))

  
# Ploting the number of stressors and ES relationship
data_plot = rbind(
  data.frame(
    Lv = rep(levels, each = n_iter),
    Assumption = rep("Actual", each = n_iter),
    Values = c(unlist(response_ES_bs$`Level 2`),unlist(response_ES_bs$`Level 5`),
               unlist(response_ES_bs$`Level 8`),unlist(response_ES_bs$`Level 10`)),
    Pcolor = rep(response_ES$color_p[12:15], each = n_iter),
    row.names = c()
  ),
  data.frame(
    Lv = rep(levels, each = n_iter*3),
    Assumption = rep(c("Additive", "Multiplicative", "Dominative"), each = n_iter),
    Values = unlist(Null_ES_bs),
    Pcolor = rep(NA, n_iter*3),
    row.names = c()
  )
)

g_ms = 
  ggplot()+
  xlab("Effect size") +
  ylab("Level: the number of stressors")+ 
  theme(legend.position = "none")  +
  geom_vline(xintercept=0,colour="#00000040") +
  xlim(low = min(unlist(ES_plot)*1.1), high=0.5) +
  #coord_flip() +
  scale_y_discrete(limits = rev(unique(data_plot$Lv))) +
  
  ### Distributions ###
  # 3 assumptions
    geom_density_ridges(data=data_plot[data_plot$Assumption=="Additive",],  
                        aes(x = Values, y = Lv, group = Lv),
                        rel_min_height = 0.025,
                        scale = 0.3, alpha = .30, fill="#9e8900",color= "#00000000", position=position_nudge(y = -0.2)) +
    geom_density_ridges(data=data_plot[data_plot$Assumption=="Multiplicative",],  
                        aes(x = Values, y = Lv, group = Lv),
                        rel_min_height = 0.025,
                        scale = 0.3, alpha = .3, fill="#ff0083",color= "#00000000", position=position_nudge(y = -0.35)) +
    geom_density_ridges(data=data_plot[data_plot$Assumption=="Dominative",],  
                        aes(x = Values, y = Lv, group = Lv),
                        rel_min_height = 0.025,
                        scale = 0.3, alpha = .3, fill="#A000A0",color= "#00000000", position=position_nudge(y = -0.5)) +
  ### Mean & CI ###
  # 3 assumptions
    geom_estci(data=ES_plot[["Additive"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=0), 
               color="#9e8900AA", linetype = "dashed",size=0.5, ci.linesize = 0.8, position=position_nudge(y = -0.21)) +
    geom_estci(data=ES_plot[["Multiplicative"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=0), 
               color="#ff0083AA", linetype = "dashed",size=0.5, ci.linesize = 0.8, position=position_nudge(y = -0.36)) +
    geom_estci(data=ES_plot[["Dominative"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=0), 
               color="#A000A0AA", linetype = "dashed",size=0.5, ci.linesize = 0.8, position=position_nudge(y = -0.51)) +
    
  # Actual ES
    geom_estci(data=ES_plot[["Actual"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=0), 
    color=response_ES$color_p[12:15], size=0.9, ci.linesize = 0.9, position=position_nudge(y = 0))
  


g_rawdata   | g_ms



#####################################################################
#####         Supplementary analysis                            #####
#####################################################################

#-------------------------------------------
#  LINEAR REGRESSION
#-------------------------------------------
lm_model = c(); lm_R2 = c()
for(response_i in responses){
  i = which(response_i==responses)
  eval(parse(text=(paste("model = summary(lm(", response_i, "~Lv, data = df))")))) 
  lm_model = rbind(lm_model, model$coefficients)
  lm_R2 = rbind(lm_R2, model$r.squared)
}



#-------------------------------------------
#  BIMODALITY
#-------------------------------------------

rescale = cbind(df[df$Lv>1,"Lv"], apply(df[df$Lv>1, responses], 2, scale),df[df$Lv>1,2:11])
colnames(rescale)[1] = "Lv"

df.bimodality = data.frame(
  Lv = rep(rescale[,1], 7),
  Attr = rep(colnames(rescale[,c(2:8)]), each=40),
  Values = unlist(list(rescale[,c(2:8)]))
)
df.bimodality = cbind(
  df.bimodality,
  rbind(rescale[,c(9:18)], rescale[,c(9:18)],rescale[,c(9:18)],rescale[,c(9:18)])
)


g_bimodal_factorwise = 
  ggplot()+
  xlab("Standardized score") +
  theme(legend.position = "none")  +
  scale_y_discrete(limits = rev(unique(df.bimodality$Attr))) +
  

  ### Distributions ###
  geom_density_ridges(data=df.bimodality,
                      aes(x = Values, y = Attr, group = Attr),
                      rel_min_height = 0.05,
                      jittered_points = TRUE,
                      scale = 0.9, alpha = .30, 
                      fill="#9e8900",color= "#00000000")



g_bimodal_levelwise = 
  ggplot()+
  xlab("Standardized score") +
  theme(legend.position = "none")  +
  scale_y_discrete(limits = rev(unique(as.factor(df.bimodality$Lv)))) +

  
  ### Distributions ###
  geom_density_ridges(data=df.bimodality,
                      aes(x = Values, y = as.factor(Lv), group = as.factor(Lv)),
                      rel_min_height = 0.05,
                      jittered_points = TRUE,
                      scale = 0.9, alpha = .30, 
                      fill="#9e8900",color= "#00000000")

g_bimodal_factorwise | g_bimodal_levelwise


