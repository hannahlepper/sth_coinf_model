using Distributions

#Summary statistics - need changing, plan to use mean and k from Sri Lanka dataset.
#These are currently dummy numbers.
#The order is always N. americanus, Ascaris, Trichuris
mean_obs = [14, 4.793951, 2.769747]
k_obs = [0.2851153, 0.5262055, 0.6750524]

sd_mean = [14, 14, 14]
sd_k = [0.1, 0.1, 0.1]

#Prior particles - using uniform distribution
a = rand(Uniform(0, 1), 333)
b = rand(Uniform(0, 1), 333)
c = rand(Uniform(0, 1), 333)

#Map model over priors.
#Use toy model for experimenting
function fit_model_run(pars)
  rand(NegativeBinomial(0.3, 0.3/(10 * sum(pars) + 0.3)), 100)
end

#fit_model_run = function(beta_priors)
  #main(10000, 1000, SpPars, 18.2, ts, halflife, pc_dr, stool_samp, beta_priors)
#end

#Generate simulated dataset
sim_ds = zeros(Int64, 333, 100)

for s in 1:333
  sim_ds[s,:] = fit_model_run([a[s], b[s], c[s]])
end

#Extract summary statistics

#Need an R function; returns the k and mean of the distribution
using RCall
R"source('fitNB.R')"
@rget fit_nb

fitted_vals = zeros(Float64, 333, 2)
for s in 1:333
  fitted_vals[s,:] = fit_nb(sim_ds[s,:])
end

#Calculate basic distances
dists = zeros(Float64, 333, 2)
dists[:,1] = fitted_vals[:,1] .- k_obs[1]
dists[:,2] = fitted_vals[:,2] .- mean_obs[1]

#Calculate standardised euclidean distance - from Mahalanobis page Wikipedia
function euc(dists, sds)
    ss(d, sd) = ( d^2 / sd^2 )
    sqrt(sum([ss(dists[i], sds[i]) for i in 1:length(dists)]))
end

euc_dists = zeros(Float64, 333, 1)
for s in 1:333
  euc_dists[s,:] = euc(dists[s,:], [sd_mean[1], sd_k[1]])
end

using Plots
histogram(euc_dists[:,1])

#Get posteriors
#Accept small percent

accept = find(i -> i < 6, euc_dists) #35 kept
a_accept = a[accept]
b_accept = b[accept]
c_accept = c[accept]

#Print un-squeezed posteriors. Look about right.
using StatPlots
density(a_accept)
density(b_accept)
density(c_accept)


#Squeeze posteriors
#1. weight according to distance

R"abc_ <- abc::abc"
@rget abc_

Na_a = abc_([0.3, 15], a, fitted_vals, 0.1, "loclinear")
Na_b = abc_([0.3, 15], b, fitted_vals, 0.1, "loclinear")
Na_c = abc_([0.3, 15], c, fitted_vals, 0.1, "loclinear")

Al_a = abc_([0.3, 15], a, fitted_vals, 0.1, "loclinear")
Al_b = abc_([0.3, 15], b, fitted_vals, 0.1, "loclinear")
Al_c = abc_([0.3, 15], c, fitted_vals, 0.1, "loclinear")

density(Na_a[:adj_values])
density(Na_b[:adj_values])
density(Na_c[:adj_values])
