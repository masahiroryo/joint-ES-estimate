#####################################################################
#
#  Rillig et al.(2019) the role of multiple global change factors in driving soil functions and microbial biodiversity 
#
#  Appendix: Data analysis in R
#
#  Author of this script: Masahiro Ryo (masahiroryo@gmail.com)
#
#####################################################################


#####################################################################
#####                    Library                                #####
#####################################################################
# NOTE) ggepi & patchwork can be installed via devtools 
#devtools::install_github("lwjohnst86/ggepi")
#devtools::install_github("thomasp85/patchwork")


library(ggplot2)
library(ggepi)
library(ggridges)
library(patchwork)
library(party)
library(caret)
library(dplyr)



#####################################################################
#####                    Data info                              #####
#####################################################################

df = data.frame(read.csv("Rillig_data_R.csv", sep=","))
df$remark = factor(df$remark, levels = unique(df$remark))

levels = c("1", "2", "5", "8", "10")
stressors = c("D", "N", "T", "M", "G", 
              "A", "F", "C", "S", "I")

treatment = as.vector(unique(df$remark))
responses = c("WSA", "wdpt.sec.", "Percent.decomposition..corr.", "CO2_ppm_ThirdWeek",
               "Rich.ASV", "Comm.Comp1", "Comm.Disp")

n_iter   = 100      # permutaion size




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
      if(treatment=="1") population = data[data$remark%in%stressors, response]
      if(treatment!="1") population = data[data$remark==treatment, response]
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

  population_CT = data[data$remark=="CT", response]

  for(treatment in target){
    bs = numeric(0)
    if(treatment=="1") population_TR = data[data$remark%in%stressors, response]
    if(treatment!="1") population_TR = data[data$remark==treatment, response]
    size_CT = length(population_CT)
    size_TR = length(population_TR)
    
    for(id in c(1:n_perm)){
      k_CT = mean(sample(population_CT, size_CT, replace = T))
      k_TR = mean(sample(population_TR, size_TR, replace = T))
      bs = append(bs, k_TR - k_CT)
    }
    resampled[[treatment]] = bs
  }
  resampled[["CT"]] = rep(0, n_perm)
  return(resampled)
}

BootStrap_ES_summary = function(data){
  summary = list()
  p = 0
  summary[["CT"]] = c(0,0,0,1)
  target = names(data)
  
  for(treatment in target[-1]){
    bs = data[[treatment]]
    p = length(which(bs>0))/length(bs)
    p = min(p, 1-p)
    summary[[treatment]] = c(quantile(bs, .025), mean(bs), quantile(bs, .975), p)
  }
  summary = t(data.frame(summary))
  colnames(summary) = c("2.5%", "mean", "97.5%", "p_value")
  summary = data.frame(target, summary); row.names(summary) = c()
  
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
    if(Lv=="1") combination = data[data$remark%in%stressors,c(2:11)]
    if(Lv!="1") combination = data[data["remark"]==Lv,c(2:11)]
    Level = sum(combination[1,])
    
      
    # Null distributions can be taken based on three different assumptions
    for(type in c("Additive", "Multiplicative", "Dominative")){

      population_CT = df[df$remark=="CT", response]
      size_CT = length(population_CT)
      
      # For each combination, bootstrap resampling is conducted
      for(j in c(1:nrow(combination))){
        bs = numeric(0)
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
            for(m in c(1:Level)) z = z * (1 + each_effect[m])
            joint_effect = (z - 1)*k_CT
          }
          if(type=="Dominative")      joint_effect = each_effect[which(max(abs(each_effect))==abs(each_effect))]
  
          bs = append(bs, joint_effect)
        }
        resampled[[type]][[j]] = bs
      }
      
    }
    output[[Lv]] = resampled
  }  
  return(output)
} 

Null_distribution_rep_transform = function(data){
  output = list()
  for(Lv in levels){
    for(type in c("Additive", "Multiplicative", "Dominative")){
      output[[Lv]][[type]] = sample(unlist(data[[Lv]][[type]]), n_iter, replace=F)
    }
  }
  return(output)
}

NHST_summary = function(null_data, Actual_data){
  output = list()
  for(Lv in levels){
    summary = list()
    summary[["Actual"]] = c(quantile(Actual_data[[Lv]], .025), mean(Actual_data[[Lv]]), quantile(Actual_data[[Lv]], .975), 1)
    p = 0
    assumptions = c("Additive", "Multiplicative", "Dominative")
    
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
    summary = rbind(data[["1"]][i, 2:4], data[["2"]][i, 2:4], data[["5"]][i, 2:4],
                    data[["8"]][i, 2:4], data[["10"]][i, 2:4])
    summary = cbind(levels, summary)
    colnames(summary) = c("Lv", "Low", "Mean", "High")
    output[[c("Actual", "Additive", "Multiplicative", "Dominative")[i]]] = summary
  }
  return(output)
}

Expected_ES_for_each = function(data){
  output = numeric(0)
  for(type in c("Additive", "Multiplicative", "Dominative")){
    tmp = numeric(0)
    for(Lv in levels){
      n_len = length(data[[Lv]][[type]])
      for(i in 1:n_len){
         tmp = append(tmp, mean(data[[Lv]][[type]][[i]]))
      }
    }
    output = cbind(output,tmp)
  }
  colnames(output)= c("E1", "E2", "E3")
  return(output)
}
    
#####################################################################
#####                  Main analysis                            #####
#####################################################################

#----------------------------------------------
# Step 1: Effect size estimate and null hypothesis significance testing
#----------------------------------------------

response_mean_all = list()
response_ES_all   = list()
joint_ES_null_all = list()
ES_for_each_all   = list()
g_rawdata_all     = list()
g_ms_all          = list()


for(i_response in responses){
  # Bootstrap estimate: single stressors
  response_mean   = BootStrap_mean(i_response)
  response_ES_bs  = BootStrap_ES_rep(i_response)
  response_ES     = BootStrap_ES_summary(response_ES_bs)
  
  # Bootstrap estimate: joint stressors
  Null_ES_bs0     = Null_distribution_rep(i_response)
  Null_ES_bs      = Null_distribution_rep_transform(Null_ES_bs0)
  joint_ES_null   = NHST_summary(Null_ES_bs, response_ES_bs)
  ES_plot         = NHST_summary_transform(joint_ES_null)
  ES_for_each     = Expected_ES_for_each(Null_ES_bs0)

  
  #################################################################
  #------ Store the information for each response ----------------
  # Summary table combining the information from all response variables
  response_mean_all[[i_response]] = response_mean
  response_ES_all[[i_response]]   = response_ES
  joint_ES_null_all[[i_response]] = bind_rows(joint_ES_null, .id = "column_label")
  ES_for_each_all[[i_response]]   = ES_for_each
  #################################################################
  
  # Ploting the raw data with the mean & its 95% confidence intervals

  g_rawdata_all[[i_response]] =  local({
   i_response = i_response
   response_mean = response_mean

   ggplot() +
     theme_bw()+
     theme(legend.position = 'none', axis.title.x=element_blank())+
     xlab(i_response) + 
     coord_flip() +
     stat_density_ridges(data=df, aes_string(x = i_response, y = "remark"),
                         geom = "density_ridges_gradient",
                         rel_min_height = 0.01, 
                         jittered_points = TRUE, color="#00000000",
                         alpha = .5,
                         position = position_points_jitter(height = .2, yoffset = .15),
                         point_size = 1, point_alpha = .3, 
                         scale = .5) +
     geom_estci(data=response_mean, aes(x = mean, y = target, xmin=X2.5., xmax=X97.5., 
                                                          xintercept=response_mean[1,"mean"]), center.linecolour = "black",
                size=0.6, ci.linesize = 0.5, position=position_nudge(y = -0.15)) 
  })

  # Ploting the number of stressors and ES relationship

  g_ms_all[[i_response]] = local({
  i_response = i_response
  ct_value = response_mean[1,"mean"]
  ES_plot = ES_plot
  for(i in 1:4) ES_plot[[i]][,2:4] = ES_plot[[i]][,2:4] + ct_value

  ggplot()+
    theme_bw()+
    theme(legend.position = 'none', axis.title.x=element_blank(), axis.title.y=element_blank())+
    theme(legend.position = "none")  +
    coord_flip() +
    scale_y_discrete(limits = factor(c("1","2","5","8","10"), levels=c("1","2","5","8","10"))) +
    
    ### Mean & CI ###
    # 3 assumptions
      geom_estci(data=ES_plot[["Additive"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=ct_value), 
                 color="#008900AA", size=0.4, ci.linesize = 0.4, position=position_nudge(y = +0.1)) +
      geom_estci(data=ES_plot[["Multiplicative"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=ct_value), 
                 color="#ff0083AA", size=0.4, ci.linesize = 0.4, position=position_nudge(y = +0.2)) +
      geom_estci(data=ES_plot[["Dominative"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=ct_value), 
                 color="#0000A0AA", size=0.4, ci.linesize = 0.4, position=position_nudge(y = +0.3)) +
      
    # Actual ES
      geom_estci(data=ES_plot[["Actual"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=ct_value), 
                                    size=0.6, ci.linesize = 0.6, position=position_nudge(y = 0))
  })

  
}


response_mean_all = bind_rows(response_mean_all, .id = "column_label")
response_ES_all   = bind_rows(response_ES_all,   .id = "column_label")
joint_ES_null_all = bind_rows(joint_ES_null_all, .id = "column_label")
colnames(joint_ES_null_all) = c("Response", "Lv", colnames(joint_ES_null_all)[3:ncol(joint_ES_null_all)]) 


# Output as csv
# write.csv(response_mean_all,   "mean estimate for each condition.csv")
# write.csv(response_ES_all[,1:6], "Unstandardized effect size (difference) for each treatment.csv")
# write.csv(joint_ES_null_all, "Test if joint effect size can be predicted.csv")



#----------------------------------------------
# Step 2: Predictability comparison among random forest models
#----------------------------------------------

df.rf = df[df[, "remark"] %in% levels,]

lv_list  = unique(df.rf[,"Lv"])
id_lv1   = which(df.rf[,"Lv"]==1)
id_lvh   = which(df.rf[,"Lv"]>1)

n_data   = nrow(df.rf)
n_lv1    = sum(df.rf[,"Lv"]==10)
n_lvh    = sum(df.rf[,"Lv"]!=1)
n_eachlv = 10
n_tree   = 100
n_iter2  = 100


rf.r2 = data.frame(matrix(NA, ncol=3, nrow=3*n_iter2*length(responses)))
rf.r2[,1] = rep(responses,each=3*n_iter2)
rf.r2[,2] = rep(c("Lv", "Lv+ID", "All"), n_iter2*length(responses))
rf.r2[,2] = factor(rf.r2[,2], levels=c("Lv", "Lv+ID", "All"))
colnames(rf.r2) = c("Response", "Model", "R2")

rf.pred = data.frame(matrix(NA, ncol=6, nrow=n_iter2*length(responses)))
rf.pred[,1] = rep(responses,each=n_iter2)
colnames(rf.pred) = c("Response", lv_list)

rf.vimp = data.frame(matrix(NA, ncol=15, nrow=n_iter2*length(responses)))
rf.vimp[,1] = rep(responses,each=n_iter2) 
colnames(rf.vimp) = c("Response", "Lv", stressors, "E1", "E2", "E3")

rf.prediction.all = list()

j = 0
for(i_response in responses){
  
  df.rf.tmp = cbind(df.rf, ES_for_each_all[[i_response]])
  eval(parse(text=(paste("fml      = formula(",i_response,"~", paste(c("Lv", stressors, colnames(ES_for_each_all[[i_response]])), collapse=" + "),")", sep=""))))
  eval(parse(text=(paste("fml.Lv   = formula(",i_response,"~Lv)", sep=""))))
  eval(parse(text=(paste("fml.LvID = formula(",i_response,"~", paste(c("Lv", stressors), collapse=" + "),")", sep=""))))
  
  for(i in 1:n_iter2){
    j = j + 1
    # bootstrap resampling
    set.seed(j)
    # take 10 sample from Lv1, and take 40 sample from the other levels for balanced resampling
    rid = c(sample(id_lv1, n_eachlv, replace=T), sample(id_lvh, n_lvh, replace=T))
    rdf = df.rf.tmp[rid,]
    rf_model      = tryCatch({
                    cforest(fml,      data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))},
                    error = function(e) {cforest(fml.Lv,   data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))})
    rf_model.Lv   = cforest(fml.Lv,   data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
    rf_model.LvID = cforest(fml.LvID, data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
    
    # shuffling corresponding variables for evaluating R2 reduction
    # rdf.test.LvID = rdf
    # rdf.test.LvID[,colnames(ES_for_each_all[[i_response]])] = apply(rdf.test.LvID[,colnames(ES_for_each_all[[i_response]])],2,sample)  
    
    # evaluating fitting performance
    rf_prdct.Lv   = tryCatch({predict(rf_model.Lv, OOB=T)}, error=function(e){rep(0,length(rid))})
    rf_prdct.LvID = tryCatch({predict(rf_model.LvID, OOB=T)}, error=function(e){rep(0,length(rid))})
    rf_prdct.All  = tryCatch({predict(rf_model, OOB=T)}, error=function(e){rep(0,length(rid))})

    rf.r2[3*j-2,3]    = postResample(rf_prdct.Lv,rdf[,i_response])[2]
    rf.r2[3*j-1,3]    = postResample(rf_prdct.LvID,rdf[,i_response])[2]
    rf.r2[3*j-0,3]    = postResample(rf_prdct.All,rdf[,i_response])[2]
    
    # variable importance
    set.seed(j)
    tmp.vimp = numeric(0)
    for(itmp in 1:5) tmp.vimp = rbind(tmp.vimp,varimp(rf_model))
    tmp.vimp = apply(tmp.vimp,2,mean)
    tmp.vimp = tmp.vimp/sum(tmp.vimp)
    rf.vimp[j, 2:(length(tmp.vimp)+1)] = tmp.vimp*rf.r2[2*j-0,3]*100
    
    # fitting curve
    tmp.curve = c()
    for(i_lv in lv_list){
      tmp.curve = append(tmp.curve, mean(rf_prdct.All[rdf[,"Lv"]==i_lv]))
    }
    rf.pred[j,2:6]    =  tmp.curve
  }
  
  rf_model      = cforest(fml,      data = df.rf.tmp, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
  rf_model.Lv   = cforest(fml.Lv,   data = df.rf.tmp, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
  rf_model.LvID = cforest(fml.LvID, data = df.rf.tmp, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
  
  rf.prediction.all[[i_response]] = data.frame(
    Model     = rep(c("Lv", "LvID", "All"), each = nrow(df.rf.tmp)),
    Predicted = c(predict(rf_model.Lv, OOB=F), predict(rf_model.LvID, OOB=F), predict(rf_model, OOB=F)),
    Observed  = rep(df.rf.tmp[,i_response], 3))
  
}
rf.r2.summary = data.frame(matrix(NA,ncol=5,nrow=3*length(responses)))
colnames(rf.r2.summary) = c("Response","Model", "CI.low", "Mean", "CI.high")
rf.r2.summary[,1] = rep(responses,each=3)
rf.r2.summary[,2] = rep(c("Lv", "Lv+ID", "All"),length(responses))
rf.r2.summary[,2] = factor(rf.r2.summary[,2],levels=c("Lv", "Lv+ID", "All"))


rf.pred.summary = data.frame(matrix(NA,ncol=5,nrow=5*length(responses)))
colnames(rf.pred.summary) = c("Response","Lv","CI.low", "Mean", "CI.high")
rf.pred.summary[,1] = rep(responses,each=length(lv_list))
rf.pred.summary[,2] = rep(lv_list, length(responses))

rf.vimp.summary = data.frame(matrix(NA,ncol=5,nrow=14*length(responses)))
colnames(rf.vimp.summary) = c("Response","Variable","CI.low", "Mean", "CI.high")
rf.vimp.summary[,1] = rep(responses,each=14)
rf.vimp.summary[,2] = rep(colnames(rf.vimp)[2:15], length(responses))
rf.vimp.summary[,2] = factor(rf.vimp.summary[,2], levels = unique(rf.vimp.summary[,2]))

j  = 0
for(i_response in responses){
  jj = 0
  jjj = 0
  rf.r2.summary[3*j+1,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="Lv"),   3],c(.025,.50,.975), na.rm=T)
  rf.r2.summary[3*j+2,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="Lv+ID"),3],c(.025,.50,.975), na.rm=T)
  rf.r2.summary[3*j+3,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="All"),  3],c(.025,.50,.975), na.rm=T)
  
  for(i_lv in lv_list){
    jj = jj + 1
    rf.pred.summary[5*j+jj,3:5]     = quantile(rf.pred[which(rf.pred[,1]==i_response),as.character(i_lv)],c(.025,.50,.975), na.rm=T)
  }
  
  for(i_var in colnames(rf.vimp)[2:15]){
    jjj = jjj + 1
    rf.vimp.summary[14*j+jjj,3:5]     = quantile(rf.vimp[which(rf.vimp[,1]==i_response),as.character(i_var)],c(.025,.50,.975), na.rm=T)
  }
  
  j = j + 1
  
}


# write.csv(rf.r2.summary, "randomforest_r2_summary.csv")

g_rf_all     = list()
g_vimp_all   = list()
g_r2_all     = list()
g_correl_all = list()
for(i in 1:(length(responses)+7)){
  g_rf_all[[i]] =  local({
    i = i
    ggrf_select = rf.pred.summary[rf.pred.summary[,1]==responses[i],]
    ggdf_select = df.rf[,c("Lv",responses[i])]
    
    ggplot(data=ggrf_select)+
      geom_point(data=ggdf_select, aes_string(x="Lv",y=responses[i]))+
      geom_line(aes(x=Lv, y=Mean)) +
      geom_ribbon(aes(x=Lv,  ymax=CI.high, ymin=CI.low), colour = NA, fill="#999999",alpha=.5)+
      theme_bw()+
      theme(legend.position = 'none', axis.title.x=element_blank(),axis.title.y=element_blank())
  })  
  g_vimp_all[[i]] = local({
    i = i
    ggplot(data=rf.vimp.summary[rf.vimp.summary[,1]==responses[i],],aes(x=Variable, y=Mean))+
      geom_bar(stat="identity", fill="#999999", alpha=0.5) +
      xlab("Variability explained [%]")+
      geom_errorbar(aes(ymax=CI.high, ymin=CI.low), width=.2, position=position_dodge(width=0.0)) +
      theme_bw()+
      theme(legend.position = 'none', axis.title.x=element_blank(),axis.title.y=element_blank())
  })
  g_r2_all[[i]] = local({
    i = i
    rf.r2 = rf.r2
    rf.r2.summary = rf.r2.summary 
    ggplot(data=rf.r2[rf.r2[,1]==responses[i],],aes(x=Model,y=R2, fill=Model))+
      geom_violin( color="#00000000",alpha=.5,position=position_dodge(width=0.3),trim=F)+
      geom_pointrange(data=rf.r2.summary[rf.r2.summary[,1]==responses[i],], aes(y=Mean, ymax=CI.high, ymin=CI.low,color=Model), position=position_dodge(width=0.2)) +
      scale_fill_manual(values  = c("#999999",  "#999999",  "#999999"))+ 
      scale_color_manual(values = c("#505050", "#505050", "#505050"))+ 
      theme_bw() + ylim(c(0,1.0))+
      theme(legend.position = 'none', axis.title.x=element_blank())
  })
  
  g_correl_all[[i]] = local({
    i = i
    rf.prediction.all = rf.prediction.all
    ggplot(data = rf.prediction.all[[i]], aes(y=Predicted, x=Observed, group=Model, color=Model))+
    theme_bw()+
    geom_abline(slope=1, intercept=0) + geom_point() + geom_smooth(method="lm", fullrange=F,size= 1.5) +
    xlim(range(rf.prediction.all[[i]][,2:3])) + ylim(range(rf.prediction.all[[i]][,2:3])) +
    theme(legend.position = 'none', axis.title.x=element_blank(),axis.title.y=element_blank()) +
    scale_color_manual(values = c("#599a9780", "#1e3a6180", "#c9b75f80"))
  }) 
  
}


#-------------------------------------------
#  Step 3: Visualization
#-------------------------------------------

g_rawdata_all[[1]]+g_ms_all[[1]]+g_rf_all[[1]]+g_r2_all[[1]]+
  g_rawdata_all[[2]]+g_ms_all[[2]]+g_rf_all[[2]]+g_r2_all[[2]]+
  g_rawdata_all[[3]]+g_ms_all[[3]]+g_rf_all[[3]]+g_r2_all[[3]]+
  g_rawdata_all[[4]]+g_ms_all[[4]]+g_rf_all[[4]]+g_r2_all[[4]]+
  plot_layout(ncol=4, widths=c(4,2,2,1))

g_rawdata_all[[5]]+g_ms_all[[5]]+g_rf_all[[5]]+g_r2_all[[5]]+
  g_rawdata_all[[6]]+g_ms_all[[6]]+g_rf_all[[6]]+g_r2_all[[6]]+
  g_rawdata_all[[7]]+g_ms_all[[7]]+g_rf_all[[7]]+g_r2_all[[7]]+
  plot_spacer()+
  plot_layout(ncol=4, widths=c(4,2,2,1))

g_correl_all[[1]] + g_correl_all[[5]]  + 
  g_correl_all[[2]] + g_correl_all[[6]]  + 
  g_correl_all[[3]] + g_correl_all[[7]]  +
  g_correl_all[[4]] + plot_spacer()  +
  plot_layout(ncol=2, widths = c(1,1))

