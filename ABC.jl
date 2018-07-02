using DataFrames

#Summary statistics
#Calculated from the first two surveys from the Sri Lankan dataset
#I use logs of raw egg counts, so that the data is more normal and
#using SD in measuring the distance is more valid.
#The order is always N. americanus, Ascaris, Trichuris
ec = [2.278871, 4.793951, 2.769747]
ec_sd = [1.365808, 2.003533, 1.551476]
p = [0.2851153, 0.5262055, 0.6750524]
p_sd = [0.02067136, 0.02286196, 0.02144451]

#Prior particles
beta_pri_N_a = rand(Uniform(1e-12, 2e-12), 33)
beta_pri_A_l = rand(Uniform(1e-13, 2e-13), 33)
beta_pri_T_t = rand(Uniform(3e-14, 4e-14), 33)

#Map model over priors
fit_model_run = function(beta_priors)
  main(10000, 1000, SpPars, 18.2, ts, halflife, pc_dr, stool_samp, beta_priors)
end

prep_data_set = function(data, sp)
  ds = DataFrame(ec = 0., p = 0.)
  for i in data
    push!(ds, hcat(i[1][sp], i[2][sp]))
  end
  deleterows!(ds, 1)
  return ds
end

N_a = [fit_model_run([x, mean(beta_pri_A_l), mean(beta_pri_T_t)]) for x in beta_pri_N_a]
N_a = prep_data_set(N_a, 1)

A_l = [fit_model_run([mean(beta_pri_N_a), x, mean(beta_pri_T_t)]) for x in beta_pri_A_l]
A_l = prep_data_set(A_l, 1)

T_t = [fit_model_run([mean(beta_pri_N_a), mean(beta_pri_A_l), x]) for x in beta_pri_T_t]
T_t = prep_data_set(T_t, 1)


#Calculate euclidean distance
ds[:wb_d] = repeat(wb, inner = 5) .- ds[:wb]
ds[:p_d] = repeat(p, inner = 5) .- ds[:p]

euc_dist = function(dists, sds)
    ss(d, sd) = ( d / sd ) ^ 2
    sum([ss(dists[i], sds[i]) for i in 1:length(dists)]) ^ 1/2
end

wb_sd = repeat(wb_sd, inner = 5)
p_sd = repeat(p_sd, inner = 5)

euc_dist_data = Float64[]
for i in 1:nrow(ds)
    push!(euc_dist_data, euc_dist([ds[:wb_d][i], ds[:p_d][i]],
                                  [wb_sd[i], p_sd[i]] ) )
end

ds[:euc] = euc_dist_data

ds

#Get posteriors
#Squeeze posteriors

Pkg.update()
