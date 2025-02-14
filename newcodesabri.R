require(GenSA)

# Percorsi file
metadata_file <- "C:/Users/spigno/Documents/GitHub/PopulationData/MetagenomicData/MetadataMetgendata.csv"
metagenomic_data_file <- "C:/Users/spigno/Documents/GitHub/PopulationData/MetagenomicData/RelabMetgendata.csv"

# Caricare i dataset
meta_microbiota <- read.csv(metadata_file, stringsAsFactors = FALSE)
microbiota <- read.csv(metagenomic_data_file, stringsAsFactors = FALSE)

save_distances = list()

pdf("~/#INSERT.pdf", width=7, height=11)
#the name of the dataframe can be changed depending on the dataset used
nonpolardcm = data.frame(matrix(NA, nrow = 0, ncol = 9))
colnames(nonpolardcm) = c("station","lambda","final_r","final_l","final_k","final_x_max","final_mqd","final_norm_const","final_n")
for (sample_counter in 2 : length(colnames(microbiota))){
  
  ###################################################################################### EXTRACTING SAMPLE
  
  ## finding station and size corresponding to loop counter
  station = colnames(microbiota)[sample_counter]
  print('---------------------------------------------')
  print(sample_counter)
  print(station)
  
  ## creating vector of data from a single station and size
  vec = microbiota[,station]
  
  ## removing zero-abundances from the vector
  vec = vec[vec != 0]
  
  ## maximum abundance value
  abund_max = max(vec)
  
  ###################################################################################### DECLARATIONS AND PARAMETERS
  
  ## setting starting value for x_max
  start_x_max = as.integer(5)
  
  ## setting end value for x_max
  end_x_max = as.integer(60)
  
  ## initialize x_max
  x_max = as.integer(0)
  
  ## initialize counter
  counter = as.integer(0)
  
  ## estimated parameters vectors, l = alpha, k = beta
  vec_estimated_r = numeric()
  vec_estimated_l = numeric()
  vec_estimated_k = numeric()
  
  ## x_max vector
  vec_x_max = numeric()
  
  ## n vector
  vec_n = numeric()
  
  ## norm_const vector
  vec_norm_const = numeric()
  
  ## j-max calculation
  j_max = log2(abund_max) + 2
  
  ## breaks and normalizing vector declaration
  breaks_vec = numeric(j_max)
  norm_vec = numeric(j_max)
  
  ## breaks vector construction
  for (j in 1 : j_max){
    breaks_vec[j] = 2^(j-1) - 1
    norm_vec[j] = 2^(j-1)
  }
  
  ## removing last element from norm_vec
  norm_vec = head(norm_vec, -1)
  
  ## creating histogram vector
  histo_data = hist(vec, breaks = breaks_vec, plot=FALSE)
  #print(histo_data)
  
  exp_value = numeric()
  distances = numeric()
  vec_mqd = numeric()
  
  ###################################################################################### LOOP OVER x_max
  
  while (x_max <= end_x_max){
    print('*********')
    ## stepping
    x_max = start_x_max + (5*counter)
    counter = counter + 1
    vec_x_max[counter] = x_max
    
    ## calculating n == n(x_max)
    n = sum(vec <= x_max)
    vec_n[counter] = n
    print('NUMBER OF DATA POINTS')
    print(n)
    
    ## defining distribution parameters
    param_rlk = numeric(3)
    
    ## defining function A(r)
    A_func = function(param_rlk){
      out = 0
      for (i in 1 : length(vec)){
        if (vec[i] <= x_max){
          out = out + vec[i]
        }
      }
      out = -(out*param_rlk[1])
      return(out)
    }
    
    ## defining function N(r,l,k)
    N_func = function(param_rlk){
      out = 0
      for (x in 1 : x_max){
        out = out + exp(-param_rlk[1]*x)*exp(lgamma(x+param_rlk[2])-lgamma(x+param_rlk[3]+1))
      }
      return(out)
    }
    
    ## defining function C(l)
    C_func = function(param_rlk){
      out = 0
      for (i in 1 : length(vec)){
        if (vec[i] <= x_max){
          out = out + lgamma(vec[i]+param_rlk[2])
        }
      }
      return(out)
    }   
    
    
    ## defining function D(k)
    D_func = function(param_rlk){
      out = 0
      for (i in 1 : length(vec)){
        if (vec[i] <= x_max){
          out = out + lgamma(vec[i]+param_rlk[3]+1)
        }
      }
      out = -out
      return(out)
    }   
    
    ## defining likelihood
    L_func = function(param_rlk){
      out = 0
      out = A_func(param_rlk) + -n*log(N_func(param_rlk)) + C_func(param_rlk) + D_func(param_rlk)
      return(-out)
    }
    
    #lower and upper bounds determined based on previous work estimates
    lower = c(0,-5,-5)
    upper = c(0.1, 8, 8)
    parametri = (lower+upper)/2
    estimated_parameters = GenSA(parametri, L_func, lower = lower, upper = upper, control = list(verbose = T, max.time = 10))$par
    
    vec_estimated_r[counter] = estimated_parameters[1]
    vec_estimated_l[counter] = estimated_parameters[2]
    vec_estimated_k[counter] = estimated_parameters[3]
    
    
    print("************")
    print(paste("x-max: ", x_max))
    print(estimated_parameters)
    
    #expected values from the model for all points lower than x_max
    exp_densities_threshold = histo_data$mids[! histo_data$mids >= x_max]
    exp_densities = exp(lgamma(exp_densities_threshold+vec_estimated_l[counter])-lgamma(exp_densities_threshold+1+vec_estimated_k[counter]))*exp(-(vec_estimated_r[counter]*exp_densities_threshold))
    
    #obtained values from the histogram up to x_max
    obs_densities = (histo_data$counts/norm_vec)/sum(histo_data$counts/norm_vec)
    length(obs_densities) = length(exp_densities_threshold)
    
    #calculating norm_const
    norm_const = sum(exp_densities)/sum(obs_densities)
    vec_norm_const[counter] = norm_const
    exp_densities = exp_densities/norm_const
    
    #calculating the mean quadratic square distance for each x_max
    distances = (log(exp_densities) - log(obs_densities))^2
    mean_quadratic_distance = sum(distances)/length(exp_densities_threshold)
    vec_mqd[counter] = mean_quadratic_distance
    
  }
  
  final_r = NaN
  final_l = NaN
  final_k = NaN
  final_x_max = NaN
  final_n = NaN
  final_mqd = NaN
  final_norm_const = NaN
  
  for (i in 1 : length(vec_mqd)) {
    if (is.finite(vec_mqd[i]) & vec_mqd[i] <= 5) {
      final_x_max = vec_x_max[i]
      
      final_r = vec_estimated_r[i]
      final_l = vec_estimated_l[i]
      final_k = vec_estimated_k[i]
      
      final_n = vec_n[i]
      final_mqd = vec_mqd[i]
      final_norm_const = vec_norm_const[i]
    }
  }

  # for (i in 1 : 1341) {
  #   if (context_general$sample_material[i] == station) {
  #     biome = context_general$biomeplot[i]
  #     size = context_general$sizeplot[i]
  #   }
  # }

  ## plotting in log-log and saving as .pdf

  par(mfrow=c(2,1))

  ## fit and log-bin pdf
  plot(histo_data$mids, ((histo_data$counts/norm_vec)/sum(histo_data$counts/norm_vec)), main=station, log = "xy", xlab = "log(abundance)", ylab = "log(density)", ylim = c(0.0001,1))

  if (is.finite(final_x_max) & is.finite(final_l)){
    curve(exp(lgamma(x+final_l)-lgamma(x+1+final_k))*exp(-(final_r*x))/final_norm_const, from=0.5, to=final_x_max, add=TRUE)
    abline(v=final_x_max, col=4)
    legend("topright", legend=c((-1+final_l-final_k),final_x_max,final_mqd), bty="o")
  }

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
  ####################################################################################### FINAL VALUES OUPUT FILE WRITING
  
  #station_number = substring(station, 6, 8)
  #station_size = substring(station, 14, 50)
  lambda = -1+final_l-final_k
  # estimated_data = data.frame(station,lambda,final_r,final_l,final_k,final_x_max,final_mqd,final_norm_const,final_n)
  # nonpolardcm = rbind(nonpolardcm, estimated_data)
} ## loop all samples
dev.off()
name_file = "~/#INSERT.csv"
write.csv(nonpolardcm, file = name_file)