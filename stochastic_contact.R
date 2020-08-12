# SEIR model for COVID19 # use colocation data
rm(list = ls())
path = "~/2019-ncov/"
library(deSolve)
library(reshape2)
library(dplyr)
library(magrittr)
library(data.table)
library(foreach)
library(doParallel)
options(stringsAsFactors = FALSE)

# Contact matrix 
location_list= c("Keelung", "NewTpi", "Taipei", "Taoyuan", "HschuCnty", "Hsinchu", "Miaoli", 
                 "Taichung", "Changhua", "Yunlin", "ChiayiCnty", "Chiayi", "Nantou", "Tainan", 
                 "Kaohsiung", "Pingtung", "Taitung", "Hualien", "Yilan")
names(location_list) = c("基隆市", "新北市", "台北市", "桃園縣", "新竹縣", "新竹市", "苗栗縣", 
                         "台中市", "彰化縣", "雲林縣", "嘉義縣", "嘉義市", "南投縣", "台南市", 
                         "高雄市", "屏東縣", "台東縣", "花蓮縣", "宜蘭縣")
numpatch= length(location_list)

normal_date = dir(paste0(path, "Taiwan-Coronavirus-Disease-Prevention-Map-Jan-27-2020/Taiwan Coronavirus Disease Prevention Map Jan 27 2020 Colocation Map/")) %>% 
              strsplit(split = "_") %>% sapply(function(x) x[[2]]) %>% sub(".csv", "", .) %>% as.Date() %>% 
              .[. > as.Date("2020-03-03")] %>% as.character() %>% .[-grep("2020-04-07|2020-06-30", .)]
Cmatrix = array(NA, dim = c(length(location_list), length(location_list), length(normal_date)), 
                dimnames = list(location_list, location_list, normal_date))
for (date in normal_date) {
  Cmatrix[, , date] = read.csv(paste0(path, "Taiwan-Coronavirus-Disease-Prevention-Map-Jan-27-2020/Taiwan Coronavirus Disease Prevention Map Jan 27 2020 Colocation Map/", 
                                      "Taiwan Coronavirus Disease Prevention Map Jan 27 2020 Colocation Map_", date, ".csv"), fileEncoding = "UTF-8") %>% 
                      filter(country == "TWN") %>% acast(polygon1_name~polygon2_name, value.var = "link_value") %>% 
                      .[names(location_list), names(location_list)] %>% set_colnames(location_list) %>% set_rownames(location_list)
}
Cmatrix = apply(Cmatrix, MARGIN = c(1,2), mean, na.rm = TRUE)
remove(date, normal_date)

# get density D, population size N
density = read.csv(paste0(path,"data/population_density_TW_city_English.csv"), fileEncoding="UTF-8") %>%
          data.table(key = "City")
loc_density = density[location_list, "population_density"] %>% unlist %>% set_names(location_list)
loc_N = density[location_list, "population"] %>% unlist %>% set_names(location_list)
remove(density)

# parameters (need to draw from distribution instead with these as means?)
Di <- 3         
De <- 3.5       
R0_wuhan= 2.4   #sensitivity analysis, 0.9, 1.2 or 2.4
adjust= 1       #difference between Wuhan and Taipei #sensitivity analysis, 0.5-2
loc_N_matrix = matrix(rep(loc_N, length(loc_N)), byrow=TRUE, nrow=length(loc_N))
diag(loc_N_matrix) = diag(loc_N_matrix) - 1
R0_matrix = (R0_wuhan*adjust*t(Cmatrix)*loc_N_matrix / (Cmatrix[3, 3]*(loc_N[3]-1)))
paras= list (Di= Di, De= De, numpatch= numpatch, N= loc_N, R0_matrix= R0_matrix) 
remove(R0_wuhan, Di, De, loc_N, R0_matrix, adjust, loc_N_matrix)

# initial conditions 
E1= rep(0, numpatch) %>% set_names(location_list)
I1= rep(0, numpatch) %>% set_names(location_list)
# I1["Keelung"]= 10 #can be changed
R1= rep(0, numpatch) %>% set_names(location_list)
S1= (paras$N - I1 - E1 - R1) %>% set_names(location_list)
  
# function to calculate event rates (dS, dE, dI, dR) given current S, E, I, R
event_rate_coloc = function (SEIR_vector, paras){
  N= paras$N
  R0_matrix= paras$R0_matrix
  Di= paras$Di
  De= paras$De
  numpatch= paras$numpatch

  S= SEIR_vector[1:numpatch]
  E= SEIR_vector[(numpatch+1):(numpatch*2)]
  I= SEIR_vector[(numpatch*2+1):(numpatch*3)]
  R= SEIR_vector[(numpatch*3+1):(numpatch*4)]

  # things that can happen, S->E, E->I, I->R
  eventSE = colSums(S * I * R0_matrix / N / Di)
  eventEI = E / De
  eventIR = I / Di
  return(c(eventSE, eventEI, eventIR))
}

# function to calculate changes in SEIR from the "event"
delta_SEIR= function(event){
  numpatch= length(event) / 3
  eventSE= event[1:numpatch]
  eventEI= event[(1+numpatch): (numpatch*2)]
  eventIR= event[(1+numpatch*2): (numpatch*3)]
  
  changeS= -eventSE
  changeE=  eventSE - eventEI
  changeI=  eventEI - eventIR
  changeR=  eventIR
  
  return(c(changeS, changeE, changeI, changeR))
}

# stochastic simulation
time_steps=10000 #can make it bigger
# SEIR_matrix= matrix(NA, ncol= numpatch*4+1, nrow=time_steps) 

nsim <- 1000 # number of replicates

decline_ratio = seq(0.9, 0, by = -0.1)

prob_loc = time_loc = array(NA, dim = c(length(decline_ratio), 10, numpatch),
                            dimnames = list(decline_ratio, 1:10, location_list))
accum_loc = array(NA, dim = c(length(decline_ratio), 10, numpatch, numpatch),
                  dimnames = list(decline_ratio, 1:10, location_list, location_list))  # decline ratio, initial I number, initial I location, accumlated I+R city 
sus_month_list = c(1, 2, 3, 4, 5, 6) 

case_threshold_list = c(50, 100, 200, 500, 1000)

reduceMethod = c("changeNotDiag", "changeOnlyDiag", "changeAll")

R0_matrix_original = paras$R0_matrix

cl = makeCluster(8)
registerDoParallel(cl)

for(method in reduceMethod){

  for (case_threshold in case_threshold_list) {

    for (sus_month in sus_month_list) {

      for (decline in decline_ratio) {
        ## Reduce all R0 matrix except diagonal
        if(method == "changeNotDiag"){       
          paras$R0_matrix = R0_matrix_original * decline * (diag(-1, numpatch, numpatch) + 1) + diag(diag(R0_matrix_original))
        } 
        ## Reduce just diag of R0 matrix
        else if(method == "changeOnlyDiag"){  
          paras$R0_matrix = R0_matrix_original * (diag(-1, numpatch, numpatch) + 1) + diag(diag(R0_matrix_original)) * decline
        } 
        ## Reduce all of R0 matrix
        else if (method == "changeAll"){      
          paras$R0_matrix = R0_matrix_original * decline
        }
        
        for (iniI in c(1:10)) {
        
          for (loc_name in location_list) {
            I1= rep(0, numpatch) %>% set_names(location_list)
            I1[loc_name] = iniI
            
            
            SEIR_matrix_master = foreach (sim = 1:nsim, .combine = "rbind", .packages = "deSolve") %dopar% {
              paras_temp = paras
              switchFlag = TRUE

              # SEIR_matrix= matrix(NA, ncol= numpatch*5+1, nrow=time_steps) 
              S1_alter = S1
              S1_alter[loc_name] = S1_alter[loc_name] - iniI
              SEIR_vector = c(S1_alter, E1, I1, R1, I1) #C1= I1
              # SEIR_matrix[1, ]= c(0, SEIR_vector)
              
              ## For examination
              # change %>% {data.frame(S = .[1:numpatch], E = .[(numpatch+1):(numpatch*2)], I = .[(numpatch*2+1):(numpatch*3)], R = .[(numpatch*3+1):(numpatch*4)])} %>% View
              # SEIR_vector %>% {data.frame(S = .[1:numpatch], E = .[(numpatch+1):(numpatch*2)], I = .[(numpatch*2+1):(numpatch*3)], R = .[(numpatch*3+1):(numpatch*4)], C = .[(numpatch*4+1):(numpatch*5)])} %>% View
              
              day_previous = 0
              step = 2
              # for (step in (2:time_steps)){
              while (TRUE) {
                e_rate= event_rate_coloc(SEIR_vector, paras_temp)    # update event rate based on SEIR vector
                delta_t= rexp(1, rate= sum(e_rate))             # sample time, exponential distribution
                event= as.vector(rmultinom(1, 1, prob= e_rate)) # sample event, multinomial distribution
            
                change= delta_SEIR(event)
                SEIR_vector[1: (numpatch*4)]= SEIR_vector[1: (numpatch*4)] + change     # update SEIR vector
                
                ## if E + I = 0, then stop this simulation (no one can infect others)
                if((sum(SEIR_vector[(numpatch+1): (numpatch*3)]) <= 0)) return(NULL)
                
                if (sum(change[(numpatch*2+1): (numpatch*3)]) >0) {
                  SEIR_vector[(numpatch*4+1): (numpatch*5)]= SEIR_vector[(numpatch*4+1): (numpatch*5)] +
                                                             change[(numpatch*2+1): (numpatch*3)]
                }
                # SEIR_matrix[step, ]= c((SEIR_matrix[(step-1),1]+ delta_t), SEIR_vector) # save data to matrix
                
                ## if the accumulated number of people who had been infected is more than set number, stop the simulation
                if(sum(SEIR_vector[(numpatch*4+1): (numpatch*5)]) >= case_threshold)
                  return(c(sim, (day_previous + delta_t), SEIR_vector))
                
                if(((day_previous + delta_t) >= (sus_month * 30.5)) & switchFlag){
                  paras_temp$R0_matrix = R0_matrix_original
                  switchFlag = FALSE
                }
                
                day_previous = day_previous + delta_t
                step = step + 1
              }  
              # return(cbind(rep(sim, time_steps), SEIR_matrix))
              return(NULL)
            }
          
            if(is.null(nrow(SEIR_matrix_master))) 
              prob_loc[as.character(decline), as.character(iniI), loc_name] = 0
            else{
              prob_loc[as.character(decline), as.character(iniI), loc_name] = nrow(SEIR_matrix_master) / nsim
              time_loc[as.character(decline), as.character(iniI), loc_name] = paste(SEIR_matrix_master[,2], collapse = ",")
              if(nrow(SEIR_matrix_master) == 1)
                accum_loc[as.character(decline), as.character(iniI), loc_name, ] = paste(SEIR_matrix_master[,(numpatch*4+1+2): (numpatch*5+2)], collapse = ",")
              else
                accum_loc[as.character(decline), as.character(iniI), loc_name, ] = apply(SEIR_matrix_master[,(numpatch*4+1+2): (numpatch*5+2)], MARGIN = 2, paste, collapse = ",")
              
            }
          }
        }
        paras$R0_matrix = R0_matrix_original
      }
      save(prob_loc, time_loc, accum_loc, file=paste0(path, "data/R0_2.4/stochasticModel/contact/stocModel_", sus_month, "month_", method, "_", case_threshold, "case.RData"))
      remove(loc_name, decline, iniI)
    }

  }
}
stopCluster(cl)
