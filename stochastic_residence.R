# SEIR model for COVID19 # use proportion of time people living in i spending in j
rm(list = ls())
path = "~/2019-ncov/"
library(deSolve)
library(reshape2)
library(magrittr)
library(data.table)
library(foreach)
library(doParallel)
options(stringsAsFactors = FALSE)

##### Pij matrix #####
Pmatrix <- read.csv(paste0(path, "data/movingData_Fi_Sij/Sij_weekday.csv"), row.names = 1) %>% as.matrix()
location_list = row.names(Pmatrix)
location_list = c("Keelung", "NewTpi", "Taipei", "Taoyuan", "HschuCnty", "Hsinchu", "Miaoli", 
                  "Taichung", "Changhua", "Yunlin", "ChiayiCnty", "Chiayi", "Nantou", "Tainan", 
                  "Kaohsiung", "Pingtung", "Taitung", "Hualien", "Yilan") %>% 
                {location_list[match(., location_list)]}
numpatch = length(location_list)
Pmatrix = Pmatrix[location_list, location_list] %>% as.matrix()

##### get density D, population size N #####
density = read.csv(paste0(path,"data/population_density_TW_city_English.csv"), fileEncoding="UTF-8") %>%
          data.table(key = "City")
loc_density = density[location_list, "population_density"] %>% unlist %>% set_names(location_list)
loc_N = density[location_list, "population"] %>% unlist %>% set_names(location_list)
remove(density)


##### parameters (need to draw from distribution instead with these as means?) #####
Di <- 3       
De <- 3.5     
R0_wuhan= 2.4 #sensitivity analysis, 0.9, 1.2 or 2.4
Density_wuhan= 1271.176254
Model= "Frequency" # Density or Frequency

##### make list for density, population and R0 of each location #####
Density <- as.vector(loc_density)
N <- as.vector(loc_N)
if (Model== "Density")
  R0 = R0_wuhan / Density_wuhan * Density
if (Model== "Frequency") 
  R0 = rep(R0_wuhan, numpatch)
  
paras= list (Di= Di, De= De, numpatch= numpatch, N= N, R0= R0, Pmatrix= Pmatrix) #fill in later
remove(R0_wuhan, Density, Density_wuhan, Di, De, N, R0, Pmatrix, loc_N)

##### initial conditions #####
E1= rep(0, numpatch) %>% set_names(location_list)
I1= rep(0, numpatch) %>% set_names(location_list)
# I1["Keelung"]= 10 #can be changed
R1= rep(0, numpatch) %>% set_names(location_list)
S1= (paras$N - I1 - E1 - R1) %>% set_names(location_list)
  
##### function to calculate event rates (dS, dE, dI, dR) given current S, E, I, R #####
event_rate_move = function (SEIR_vector, paras){
  N= paras$N
  R0= paras$R0
  Di= paras$Di
  De= paras$De
  P= paras$Pmatrix
  numpatch= paras$numpatch
  
  S= SEIR_vector[1:numpatch]
  E= SEIR_vector[(numpatch+1):(numpatch*2)] 
  I= SEIR_vector[(numpatch*2+1):(numpatch*3)]
  R= SEIR_vector[(numpatch*3+1):(numpatch*4)]
  
  # things that can happen, S->E, E->I, I->R
  eventSE = (S * diag(P) * R0 / Di * colSums(I*P) / colSums(N*P)) + 
            (S * colSums((t(P) * I * R0 / N) * (diag(-1, numpatch, numpatch) + 1)) / Di)
  eventEI = E / De
  eventIR = I / Di
  return(c(eventSE, eventEI, eventIR))
}

##### function to calculate changes in SEIR from the "event" #####
delta_SEIR= function(event){
  numpatch= length(event)/3
  eventSE= event[1:numpatch]
  eventEI= event[(1+numpatch): (numpatch*2)]
  eventIR= event[(1+numpatch*2): (numpatch*3)]
  
  changeS= -eventSE
  changeE=  eventSE - eventEI
  changeI=  eventEI - eventIR
  changeR=  eventIR
  
  return(c(changeS, changeE, changeI, changeR))
}

##### Stochastic simulation #####

nsim <- 1000 # number of replicates

decline_ratio = seq(0.9, 0, by = -0.1)

prob_loc = time_loc = array(NA, dim = c(length(decline_ratio), 10, numpatch),
                            dimnames = list(decline_ratio, 1:10, location_list))
accum_loc = array(NA, dim = c(length(decline_ratio), 10, numpatch, numpatch),
                  dimnames = list(decline_ratio, 1:10, location_list, location_list))  # decline ratio, initial I number, initial I location, accumlated I+R city 
sus_month_list = c(1, 2, 3, 4, 5, 6) 

reduceMethod = c("changeNotDiag", "changeR0", "changeAll")

case_threshold_list = c(50, 100, 200, 500, 1000)

Pij_original = paras$Pmatrix
R0_backup = paras$R0

cl = makeCluster(8)
registerDoParallel(cl)

for(method in reduceMethod){

  for (case_threshold in case_threshold_list) {

    for (sus_month in sus_month_list) {
      
      for (decline in decline_ratio) {  
        ## Reduce all Pij matrix except diagonal
        if (method == "changeNotDiag"){
          paras$Pmatrix = Pij_original * decline * (diag(-1, numpatch, numpatch) + 1) + diag(diag(Pij_original))
          diag(paras$Pmatrix) = diag(paras$Pmatrix) + 1-rowSums(paras$Pmatrix)
        } 
        ## Reduce R0
        else if (method == "changeR0"){
          paras$R0 = R0_backup * decline
        } 
        ## Reduce all Pij matrix except diagonal and R0
        else if (method == "changeAll"){
          paras$Pmatrix = Pij_original * decline * (diag(-1, numpatch, numpatch) + 1) + diag(diag(Pij_original))
          diag(paras$Pmatrix) = diag(paras$Pmatrix) + 1-rowSums(paras$Pmatrix)
          paras$R0 = R0_backup * decline
        }
        
        
        for (iniI in c(1:10)) {
          
          for (loc_name in location_list) {
            I1= rep(0, numpatch) %>% set_names(location_list)
            I1[loc_name] = iniI
            
            
            SEIR_matrix_master = foreach (sim = 1:nsim, .combine = "rbind", .packages = "deSolve") %dopar% {
              paras_temp = paras
              switchFlag = TRUE
              
              S1_alter = S1
              S1_alter[loc_name] = S1_alter[loc_name] - iniI
              SEIR_vector = c(S1_alter, E1, I1, R1, I1) #C1= I1
              
              ## For examination
              # change %>% {data.frame(S = .[1:numpatch], E = .[(numpatch+1):(numpatch*2)], I = .[(numpatch*2+1):(numpatch*3)], R = .[(numpatch*3+1):(numpatch*4)])} %>% View
              # SEIR_vector %>% {data.frame(S = .[1:numpatch], E = .[(numpatch+1):(numpatch*2)], I = .[(numpatch*2+1):(numpatch*3)], R = .[(numpatch*3+1):(numpatch*4)], C = .[(numpatch*4+1):(numpatch*5)])} %>% View
              
              day_previous = 0
              step = 2
              while (TRUE) {
                e_rate= event_rate_move(SEIR_vector, paras_temp)    # update event rate based on SEIR vector
                delta_t= rexp(1, rate= sum(e_rate))                 # sample time, exponential distribution
                event= as.vector(rmultinom(1, 1, prob= e_rate))     # sample event, multinomial distribution
                
                change= delta_SEIR(event)
                SEIR_vector[1: (numpatch*4)]= SEIR_vector[1: (numpatch*4)] + change     # update SEIR vector
                
                ## if E + I = 0, then stop this simulation (no one can infect others)
                if((sum(SEIR_vector[(numpatch+1): (numpatch*3)]) <= 0)) return(NULL)
                
                if (sum(change[(numpatch*2+1): (numpatch*3)]) >0) {
                  SEIR_vector[(numpatch*4+1): (numpatch*5)]= SEIR_vector[(numpatch*4+1): (numpatch*5)] +
                    change[(numpatch*2+1): (numpatch*3)]
                }
                
                ## if the accumulated number of people who had been infected is more than set number, stop the simulation
                if(sum(SEIR_vector[(numpatch*4+1): (numpatch*5)]) >= case_threshold)
                  return(c(sim, (day_previous + delta_t), SEIR_vector))
                
                if(((day_previous + delta_t) >= (sus_month * 30.5)) & switchFlag){
                  paras_temp$Pmatrix = Pij_original
                  paras_temp$R0 = R0_backup
                  switchFlag = FALSE
                }
                
                day_previous = day_previous + delta_t
                step = step + 1
              }  
              return(NULL)
            }
            
            if(is.null(nrow(SEIR_matrix_master))) 
              prob_loc[as.character(decline), as.character(iniI), loc_name] = 0
            else{
              if(nrow(SEIR_matrix_master) == 1){
                accum_loc[as.character(decline), as.character(iniI), loc_name, ] = SEIR_matrix_master[,(numpatch*4+1+2): (numpatch*5+2)]
              }
              else{
                accum_loc[as.character(decline), as.character(iniI), loc_name, ] = SEIR_matrix_master[,(numpatch*4+1+2): (numpatch*5+2)] %>%
                                                                                   apply(MARGIN = 2, paste, collapse = ",")
              }
              prob_loc[as.character(decline), as.character(iniI), loc_name] = nrow(SEIR_matrix_master) / nsim
              time_loc[as.character(decline), as.character(iniI), loc_name] = paste(SEIR_matrix_master[,2], collapse = ",")
            }
          }
        }
        paras$Pmatrix = Pij_original
        paras$R0 = R0_backup
      }
      save(prob_loc, time_loc, accum_loc, file=paste0(path, "data/R0_2.4/stochasticModel/residence/stocModel_", sus_month, "month_", method, "_", case_threshold, "case.RData"))
      remove(loc_name, iniI, decline)
    }

  }
}


stopCluster(cl)
remove(decline_ratio, Pij_original)