require(GenSA)

### 📌 1. Caricamento dati
metadata_file <- "C:/Users/spigno/Documents/GitHub/PopulationData/MetagenomicData/filtered_metadata.csv"
metagenomic_data_file <- "C:/Users/spigno/Documents/GitHub/PopulationData/MetagenomicData/metagenomic_data_filtered.csv"

meta_microbiota <- read.csv(metadata_file, stringsAsFactors = FALSE)
microbiota <- read.csv(metagenomic_data_file, stringsAsFactors = FALSE)

### 📌 2. Selezione di UNA colonna per il debug
station = colnames(microbiota)[2]  # Usa la seconda colonna se la prima è un indice
vec = microbiota[, station]

# Rimuove abbondanze zero
vec = vec[vec != 0]
abund_min = min(vec)

print(paste("📍 Analizzando la stazione:", station))
print(paste("📉 Minima abbondanza:", abund_min))

### 📌 3. Parametri iniziali e intervalli
nbins = 12
mqd_treshold = 100
n_iter = 5
value_xmax = 20

Rmin = 0
Rmax = 0.001
lammin = -2
lammax = 100

start_x_max = 0.1
end_x_max = 10  # Range ridotto per debug
x_max = 0
counter = 0

vec_estimated_r = numeric()
vec_estimated_l = numeric()
vec_x_max = numeric()
vec_mqd = numeric()
vec_norm_const = numeric()

### 📌 4. Creazione dell'istogramma
j_min = log10(abund_min)
breaks_vec = 10^seq(j_min, 2, length.out = nbins)
norm_vec = head(breaks_vec, -1)

histo_data = hist(vec, breaks = breaks_vec, plot = FALSE)
print(histo_data)

### 📌 5. Funzioni di supporto
A_func = function(param_rlk) {
  out = sum(vec[vec <= x_max]) * param_rlk[1]
  return(out)
}

N_func = function(param_rlk) {
  out = sum((1:x_max)^(-param_rlk[2]) * exp(-param_rlk[1] * (1:x_max)))
  return(out)
}

C_func = function(param_rlk) {
  out = sum(vec[vec <= x_max]^-param_rlk[2])
  return(out)
}

L_func = function(param_rlk) {
  out = A_func(param_rlk) - length(vec[vec <= x_max]) * log(N_func(param_rlk)) + C_func(param_rlk)
  print(paste("🔎 Testing lambda:", param_rlk[2], "Loss:", abs(out)))
  return(abs(out))
}

### 📌 6. Ciclo su x_max
while (x_max <= end_x_max) {
  print(paste('🔄 Iterazione per x_max =', x_max))
  
  x_max = start_x_max + (n_iter * counter)
  counter = counter + 1
  vec_x_max[counter] = x_max
  
  n = sum(vec <= x_max)
  
  # Impostazione parametri iniziali
  lower = c(Rmin, lammin)
  upper = c(Rmax, lammax)
  parametri = c(0.01, 1)  # Cambiato da (0.01,5) per non partire troppo vicino al bordo
  
  # 🔍 Debug: ottimizzazione con GenSA
  estimated_parameters = GenSA(parametri, L_func, lower = lower, upper = upper, control = list(verbose = TRUE, maxit = 100))$par
  vec_estimated_r[counter] = estimated_parameters[1]
  vec_estimated_l[counter] = estimated_parameters[2]
  
  print(paste("✅ Fit completato per x_max =", x_max))
  print(paste("🔹 r stimato:", estimated_parameters[1], "🔹 lambda stimato:", estimated_parameters[2]))
  
  ### 📌 7. Misura della distanza quadratica media
  exp_densities_threshold = histo_data$mids[!histo_data$mids >= x_max]
  exp_densities = (exp_densities_threshold^(-vec_estimated_l[counter])) * exp(-(vec_estimated_r[counter] * exp_densities_threshold))
  
  obs_densities = (histo_data$counts/norm_vec)/sum(histo_data$counts/norm_vec)
  length(obs_densities) = length(exp_densities_threshold)
  
  norm_const = sum(exp_densities)/sum(obs_densities)
  vec_norm_const[counter] = norm_const
  exp_densities = exp_densities / norm_const
  
  distances = (log(exp_densities) - log(obs_densities))^2
  mean_quadratic_distance = sum(distances) / length(exp_densities_threshold)
  vec_mqd[counter] = mean_quadratic_distance
}

### 📌 8. Selezione del miglior fit
best_idx = which.min(vec_mqd[is.finite(vec_mqd)])
final_x_max = vec_x_max[best_idx]
final_r = vec_estimated_r[best_idx]
final_l = vec_estimated_l[best_idx]
final_norm_const = vec_norm_const[best_idx]

print(paste("🏆 Miglior Fit: x_max =", final_x_max, "r =", final_r, "lambda =", final_l))

### 📌 9. Plot Log-Likelihood vs Lambda
lambda_range = seq(lammin, lammax, length.out = 100)
log_likelihood_values = sapply(lambda_range, function(l) {
  L_func(c(0.01, l))
})

plot(lambda_range, log_likelihood_values, type = "l", main = "Log-Likelihood vs Lambda", xlab = "Lambda", ylab = "Neg Log-Likelihood")

### 📌 10. Plot Distribuzione
#pdf("Debug_PowerLaw_Fit.pdf", width = 7, height = 11)
par(mfrow = c(2, 1))

plot(histo_data$mids, ((histo_data$counts/norm_vec)/sum(histo_data$counts/norm_vec)), log = "xy", xlab = "log(abundance)", ylab = "log(density)", main = c(station))
if (is.finite(final_x_max) & is.finite(final_l)){
  curve((x^(-final_l)) * exp((-final_r * x)) / final_norm_const, from = abund_min, to = final_x_max, add = TRUE)
  abline(v = final_x_max, col = 4)
  legend("topright", legend = c(paste("λ:", final_l), paste("x_max:", final_x_max)), bty = "o")
}

plot(sort(vec, decreasing = TRUE), log = "xy", xlab = "log(rank)", ylab = "log(abundance)", main = "Rank-Abundance Plot")
if (is.finite(final_x_max)) abline(h = final_x_max, col = 4)


lambda_range = seq(lammin, lammax, length.out = 10)  # Ridurre per debug
log_likelihood_values = sapply(lambda_range, function(l) {
  value = L_func(c(0.01, l))
  print(paste("Lambda:", l, "Loss:", value))  # Debug
  return(value)
})

valid_idx = which(is.finite(log_likelihood_values))  # Indici validi

if (length(valid_idx) > 0) {
  plot(lambda_range[valid_idx], log_likelihood_values[valid_idx], type = "l", 
       main = "Log-Likelihood vs Lambda", xlab = "Lambda", ylab = "Negative Log-Likelihood")
} else {
  print("⚠️ Nessun valore finito trovato per la Loss. Controllare L_func()!")
}
