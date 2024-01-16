# Server inputs: Dose, dosing regimen, dosing frequency,
# dosing cycle definition, number of dosing cycles

# read objects from "rx_shiny_data.rda" in the  AppDir folder,
# objects include: pkMod1, pkMod1_params, inits

load("rx_shiny_data.rda")
require(Hmisc)

######################Simulate with variability######################
# Perform simulation 100 times
sim_params<-function(nsub,ppkparams){ 
  sigma<-matrix(c(0.0798,0.0628,0.0628,0.0661),2,2)
  mv<-mvrnorm(n=nsub,rep(0,2),sigma)
  CL<-ppkparams["CL"]*exp(mv[,1])
  V1<-ppkparams["V1"]*exp(mv[,2])
  
  sigma1<-matrix(c(0.000179),1,1)
  mv1<-mvrnorm(n=nsub,rep(0,1),sigma1)
  V2<-ppkparams["V2"]*exp(mv1[,1])
  Q<-ppkparams["Q"]
  
  sigma2<-matrix(c(0.0629,0.0505,0.0505,0.0521),2,2)
  mv2<-mvrnorm(n=nsub,rep(0,2),sigma2)
  alpha<-ppkparams["alpha"]*exp(mv2[,1])
  beta <-ppkparams["beta"]*exp(mv2[,2])+1   
  
  params.all<- data.frame(CL=CL,V1=V1,Q=Q,V2=V2,alpha=alpha,beta=beta)
  return(params.all)
}

get.PK <- function(ds1, ds2, ds3, params.mod) {    
  time<-sort(c(0:(14*24))) 
  res <- NULL #Create an empty matrix for storing results
  dose1<-ds1*1000  #convert mg to mcg
  dose2<-ds2*1000  #convert mg to mcg
  dose3<-ds3*1000  #convert mg to mcg
  
  #define time-varying parameters
  params.mod<-params.mod %>%
    mutate(CL_1=CL*(ds1/300)^0.401,
           CL_2=CL*(ds2/300)^0.401,
           CL_3=CL*(ds3/300)^0.401,
           V1_1=V1*(ds1/300)^0.252,
           V1_2=V1*(ds2/300)^0.252,
           V1_3=V1*(ds3/300)^0.252)

  for (i in 1:nrow(params.mod)){
    ## loading doses
    day1_1 <- et(timeUnits="hr") %>%
      et(amt=dose1)
    day1_2 <- et(timeUnits="hr") %>%
      et(amt=dose2,ii=reg1,addl=2)
    ## maintenance dose
    day3 <- et(timeUnits="hr") %>%
      et(amt=dose3,ii=reg1,addl=(ncyc1-1))
    ##Event table
    et <- seq(day1_1,day1_2,day3)  %>%
      add.sampling(set_units(time,hours))
    
    #add time-varying parameters
    et<- et %>% 
      mutate(CL=ifelse(time <= set_units(12, h), params.mod[i,"CL_1"], 
                    ifelse(time <= set_units(36, h),params.mod[i,"CL_2"],
                           params.mod[i,"CL_3"])),
             V1=ifelse(time <= set_units(12, h), params.mod[i,"V1_1"], 
                       ifelse(time <= set_units(36, h),params.mod[i,"V1_2"],
                              params.mod[i,"V1_3"]))
      )
    
    x<-rxSolve(pkMod1,params.mod[i,] %>% dplyr::select(-c(CL,V1,CL_1,CL_2,CL_3,V1_1,V1_2,V1_3)),
               et,inits,method = "lsoda")
    res<-cbind(res,x[,"C1"])
  }
  
  res.q.t <- apply(res, 1, quantile, prob = c(.05, .5, .95))
  res.q.t<-data.frame(time=time+24,t(res.q.t))   #Correct time to start from Day 1
  
  write.csv(res.q.t,"simPK.csv",row.names=F)
  return(res.q.t)
}

calc_auc<-function(data)
{
  auccalc<-data
  auccalc$h=auccalc$time*auccalc$PK
  auccalc$lagH = Lag(auccalc$h,shift=1)
  auccalc$lagTime = Lag(auccalc$time,shift=1)
  auccalc$lagConc = Lag(auccalc$PK,shift=1)
  auccalc$trapAuc=(auccalc$time-auccalc$lagTime)*(auccalc$PK+auccalc$lagConc)/2
  return(sum(auccalc$trapAuc,na.rm=T))
}

get.PKparams<-function(ds,start,end){
  res.pk <- ds %>%
    gather(quantile,PK,-time) %>%
    mutate(quantile=factor(quantile,levels=c("X5.","X50.","X95."),labels=c("5th","Median","95th")))
  doseint<-res.pk %>% filter(time>=start,time<=end)
  cavg1<-NULL
  for (i in unique(doseint$quantile)){
    cavg1<-rbind(cavg1,data.frame(Quantile=i,
                                  AUC=calc_auc(doseint %>%
                                                 filter(quantile==i))))
  }
  sum<-inner_join(res.pk %>%
                    group_by(quantile) %>%
                    filter(time==end) %>%
                    dplyr::select(-time),
                  doseint %>%
                    group_by(quantile) %>%
                    dplyr::summarise(Cmax=max(PK)))
  names(sum)<-c("Quantile","Cmin, ng/mL","Cmax, ng/mL")
  sum<-inner_join(sum,cavg1)
  names(sum)[4]<-c("AUC(0-24h), ng/mL")
  return(sum)
}

# Define server logic 
shinyServer(function(input, output) {
  ds1 <- reactive({input$Dose1}) 
  ds2 <- reactive({input$Dose2})
  ds3 <- reactive({input$Dose3}) 
  semilog<-reactive({switch(input$semilog, "linear"=F, "semi-log"=T)})
  
  sim.params <- eventReactive(input$submit_mod, {
    pkparams<-sim_params(200,pkMod1_params) 
    return(pkparams)
  })
  
  get.cp <- eventReactive(input$submit_mod, {   
    get.PK(ds1(), ds2(), ds3(), sim.params())
  })
    
  output$CpPlot <- renderPlot({
    res.pk <- get.cp() %>%
      mutate(time=time/(24))                         #convert hr to day
    
    p1<-ggplot(res.pk, aes(time,X50.))+
      labs(x="Time, d",y="Conc, ng/mL") +
      geom_line() + 
      geom_ribbon(aes(ymin=X5.,ymax=X95.),alpha=.2)+
      scale_x_continuous(breaks=c(1,3,6,10,14))+
      theme_bw()+
      theme(panel.grid.minor = element_blank())
    
    if (semilog()) {
      p1<-p1+
        scale_y_log10()
        # coord_cartesian(ylim = c(0.01,NA))
    }
    # if (overlayArm2()) {
    #   res.pk_2 <- sim2.PK()
    #   p1<-p1 +
    #     geom_line(data=res.pk_2,aes(time,X50.),col="blue") +
    #     geom_ribbon(data=res.pk_2,aes(ymin=X5.,ymax=X95.),alpha=.2,fill="blue")
    # }
    print(p1)
  })
  
  #Tables of PK parameters for first dosing interval
  output$day3params <- renderTable({
    get.PKparams(get.cp(),3*24,4*24)
  }, caption = "PK parameters after the Day 3 dose:",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  output$ssparams <- renderTable({
    get.PKparams(get.cp(),10*24,11*24)
  }, caption = "PK parameters at Steady State (after Day 10 doses):",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
})

