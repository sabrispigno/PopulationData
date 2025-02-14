require(GenSA)

# Percorsi file
metadata_file <- "C:/Users/spigno/Documents/GitHub/PopulationData/MetagenomicData/filtered_metadata.csv"
metagenomic_data_file <- "C:/Users/spigno/Documents/GitHub/PopulationData/MetagenomicData/metagenomic_data_filtered.csv"

# Caricare i dataset
meta_microbiota <- read.csv(metadata_file, stringsAsFactors = FALSE)
microbiota <- read.csv(metagenomic_data_file, stringsAsFactors = FALSE)
microbiota <- microbiota[,-1]
save_distances = list()

### change the parameters
nbins=12
mqd_treshold=100
n_iter=3
value_xmax= 70 # la massima abbondanza relativa 
Rmin=0
Rmax=100
lammin=0
lammax=2

### to create the pdf action here
#   pdf("~/#fITT10.pdf", width=7, height=11)
#   the name of the dataframe can be changed depending on the dataset used

for (sample_counter in 1){ #:length(microbiota)
  
  ###################################################################################### EXTRACTING SAMPLE
  
  ## finding sample and size corresponding to loop counter
  station = colnames(microbiota)[1]
  diet= meta_microbiota[1,"diet"]
  print('---------------------------------------------')
  print(sample_counter)
  print(station)
  
  ## creating vector of data from a single station and size
  vec = microbiota[,station]
  
  ## removing zero-abundances from the vector
  vec = vec[vec != 0]
  
  ## maximum abundance value
  abund_min = min(vec)
  
  #### DATA HISTOGRAM #####
  ## creating a sequence of bins depending on the last one
  ## j-min calculation
  j_min = log10(abund_min)
  
  # is the position of the start and the end of each bin
  breaks_vec = 10^seq(j_min, 2, length.out = nbins) # mettere la lunghezza di breaks vec maggiore di uno 
  
  # Norm vec is the distance between each bin
  norm_vec=numeric()
  
  for (i in 1:(length(breaks_vec)-1)){
    norm_vec[i]= breaks_vec[i+1]-breaks_vec[i]
  }
  
  norm_vec= head(breaks_vec, -1)
  
  # ## creating histogram of each sample at each time 
  histo_data = hist(vec, breaks = breaks_vec, plot=FALSE)
  print(histo_data)
  
  ###### DECLARATIONS AND PARAMETERS #####
  
  ## setting starting value for x_max
  start_x_max = as.integer(10) # da cambiare in base al dato
  # da dove il modello inizia a considerare la percentuale di abbondanza
  
  ## setting end value for x_max
  end_x_max = as.integer(value_xmax)
  
  ## initialize x_max
  x_max = as.integer(0)
  
  ## initialize counter
  counter = as.integer(0)
  
  ## estimated parameters vectors, r = beta, l = lambda
  vec_estimated_r = numeric()
  vec_estimated_l = numeric()

  
  ## x_max vector
  vec_x_max = numeric()
  
  ## n vector
  vec_n = numeric()
  
  ## norm_const vector
  vec_norm_const = numeric()

  ################ LOOP OVER x_max #############
  
  exp_value = numeric()
  distances = numeric()
  vec_mqd = numeric()
  
  while (x_max <= end_x_max){
    print('*********')
    ## stepping
    x_max = start_x_max + (n_iter*counter) # parte da 0 e con counter puoi cambiare il numero di iterazioni
    counter = counter + 1
    vec_x_max[counter] = x_max
    
    ## calculating n == n(x_max)
    n = sum(vec <= x_max) # how much species are below the x_max
    
    vec_n[counter] = n # every interation it counts the number of species (datapoints) that are minor of the current xmax
    
    print('NUMBER OF DATA POINTS')
    print(n)
    
    #####-----Defining functions of Power Law----
    ## defining distribution parameters
    param_rlk = numeric(2) # i parametri del modello power law
    
    ## defining function A(r) # è la funzione che nel power law definisce il cutoff
    # Nei data Tara oceans va da 0 a 10^-2
    # it is the function to model the parameter r 
    A_func = function(param_rlk){
      out = 0
      for (i in 1 : length(vec)){
        if (vec[i] <= x_max){
          out = out + vec[i]
        }
      }
      out = -(out*param_rlk[1]) ## qui dovrebbe modellare -rx
      return(out)
    }

    # It is the funciont 
    ## defining function N(r,lam)
    N_func = function(param_rlk){
      out = 0
      for (x in 1 : x_max){
        out = out + (x^(-param_rlk[2]))*exp(-param_rlk[1]*x)
      }
      return(out)
    } 
    #Power Law 
    
    ## defining function C(l) # Dovrebbe essere la funzione del parametro 2 quindi del nuovo lambda
    C_func = function(param_rlk){
      out = 0
      for (i in 1 : length(vec)){
        if (vec[i] <= x_max){
          out = out + (vec[i]^-param_rlk[2]) # visto che dovrebbe essere gamma vec x e par alfa, qua fa gamma x+alfa io faccio direttamente x^-lam
        }
      }
      #out = -out
      return(out)
    }
    
    ## defining likelihood ### Questa è l'unica cosa da cambiare che non capisco per fare il modello con solo esponenziale 
    L_func = function(param_rlk){
      out = 0
      out = A_func(param_rlk) + - n * log(N_func(-param_rlk)) + C_func(-param_rlk) 
      return((-out))
    }
    ######
    
    #lower and upper bounds determined based on previous work estimates
    lower = c(Rmin,lammin) # [R, alpha (l), beta (k)] Le references per R erano da 0 a 10^-5
    upper = c(Rmax, lammax)
    parametri = (lower+upper)/2

    estimated_parameters = GenSA(parametri, L_func, lower = lower, upper = upper, control = list(verbose = T, maxit = 300))$par
    
    vec_estimated_r[counter] = estimated_parameters[1]
    vec_estimated_l[counter] = estimated_parameters[2]

    print("************")
    print(paste("x-max: ", x_max))
    print(estimated_parameters)
    
    ### Misura la distanza quadratica media !!
    
    #Estimated
    #expected values from the model for all points lower than x_max
    exp_densities_threshold = histo_data$mids[! histo_data$mids >= x_max] #cosa significa? Toglie i mids maggiori di x max
    # exp densities è il modello da fittare Power bend (se è bend dipende dalla r che stima)
    exp_densities = (exp_densities_threshold^(-vec_estimated_l[counter]))*exp(-(vec_estimated_r[counter]*exp_densities_threshold))
    #calcola i valori che la funzione stima in corrispondenza dei bins 
    #observed
    #obtained values from the histogram up to x_max
    obs_densities = (histo_data$counts/norm_vec)/sum(histo_data$counts/norm_vec)
    length(obs_densities) = length(exp_densities_threshold)
    
    #Costante di normalizzazione
    
    # è la somma delle densità attese e delle densità osservate
    #calculating norm_const
    norm_const = sum(exp_densities)/sum(obs_densities)
    vec_norm_const[counter] = norm_const
    exp_densities = exp_densities/norm_const
    
    #calculating the mean quadratic square distance for each x_max
    # se il midpoint fosse 0 ovvero che in quel bin hai 0 osservazioni,
    # la distanza log di 0 andrebbe all'infinito
    distances = (log(exp_densities) - log(obs_densities))^2 #distanza quadratica media
    mean_quadratic_distance = sum(distances)/length(exp_densities_threshold)
    vec_mqd[counter] = mean_quadratic_distance
    
  }
  
  final_r = NaN
  final_l = NaN
  final_x_max = NaN
  final_n = NaN
  final_mqd = NaN
  final_norm_const = NaN
  
  for (i in 1 : length(vec_mqd)) {
    if (is.finite(vec_mqd[i]) & vec_mqd[i] <= mqd_treshold) { # controlla che la distanza quadratica media sia finita e minore di un treshold
      final_x_max = vec_x_max[i]
      
      final_r = vec_estimated_r[i]
      final_l = vec_estimated_l[i]

      final_n = vec_n[i]
      final_mqd = vec_mqd[i]
      final_norm_const = vec_norm_const[i]
    }
  }
  ######---
  
  ## plotting in log-log and saving as .pdf
  
  par(mfrow=c(2,1))
  
  ## fit and log-bin pdf
  plot(histo_data$mids, ((histo_data$counts/norm_vec)/sum(histo_data$counts/norm_vec)), main=c(station, diet), log = "xy", xlab = "log(abundance)", ylab = "log(density)")
  
  
  if (is.finite(final_x_max) & is.finite(final_l)){
    #equazione power law bend
    curve((exp(-final_l))*exp(-(final_r*x))/final_norm_const, from=abund_min, to=final_x_max, add=TRUE)
    abline(v=final_x_max, col=4) # dove printa x_max  ## che cos'e final norm_const???abline(v=final_x_max, col=4) # dove printa x_max  ## che cos'e final norm_const???abline(v=final_x_max, col=4) # dove printa x_max  ## che cos'e final norm_const???
    legend("topright", legend=c(paste("λ:", final_l), paste("x_max:", final_x_max), paste("final_mqd:", final_mqd)), bty="o")
  } # il primo è lambda, il secondo è xmax e il terzo è la distanza quadratica media
  
  if (is.nan(final_x_max) & is.nan(final_l)){
    legend("topleft", legend="NO FIT", bty="n", cex=3)
  }
  
  ## rank-plot
  vec = sort(vec, decreasing = TRUE)
  plot(vec, main=station, log = "xy", xlab = "log(rank)", ylab = "log(abundance)")
  if (is.finite(final_x_max) & is.finite(final_l)){
    abline(h=final_x_max, col=4)
  }
  
  save_distances[[sample_counter]] = vec_mqd
  ################################### FINAL VALUES OUPUT FILE WRITING
  lambda = final_l

} ## loop all samples
#dev.off()
name_file = "~/#INSERT.csv"
#write.csv(nonpolardcm, file = name_file)