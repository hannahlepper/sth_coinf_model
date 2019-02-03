using Distributions

# Sumary statistics from the data - based on Sri Lanka data, obtained using
# fitdistr function from MASS
mean_obs = [6.9, 249.6, 38.8]
k_obs = [0.08, 0.1, 0.21]

# Prior particles - using uniform distribution
Imme_activation = rand(Uniform(0, 1),33)
Immf_activation = rand(Uniform(0, 1),33)
est_modulation = rand(Uniform(0, 1),33)
fec_modulation = rand(Uniform(0, 1),33)


#Map model over priors.

function run_particle(sp, particle)
  eggs = zeros(Float64, 2000, 3)

  SpPars[sp].Imme_activation = particle[1] #need to change to be the parameters I am actually fitting to
  SpPars[sp].Immf_activation = particle[2]
  SpPars[sp].est_modulation = particle[3]
  SpPars[sp].fec_modulation = particle[4]

  #rr is 0 (if set run_record = 0 get record of whole run); a is ages
  rr, final_record, a = main(1000, 2000)
  eggs[:,:] = final_record[:EOut][:,:]
  return eggs
end

#Generate simulated dataset
sim_ds = zeros(Float64, 2000, 3, 33) #empty to contain the simulated values

@time for s in 1:33
  sim_ds[:,:,s] = run_particle(1, [Imme_activation[s],
                                   Immf_activation[s],
                                   est_modulation[s],
                                   fec_modulation[s]])
end #fill with simulated values

#Extract summary statistics
#Need to simulate sampling process
function stool_sample(inds_eggs, ts=1/73, stool_size=100, stool_samp_size=0.0054)
  eggs_per_day = inds_eggs .* (ts * 365)
  eggs_per_gram_stool = eggs_per_day./stool_size
  eggs_in_samp = eggs_per_gram_stool .* stool_samp_size
  map(i -> rand(Poisson(i)), eggs_in_samp)
end

sim_ds_sampled = stool_sample(sim_ds)

#Need an R function; returns the k and mean of the distribution
using RCall
R"source('fitNB.R')"
@rget fit_nb

fitted_vals = zeros(Float64, 2, 3, 33)
for s in 1:33
  for sp in 1:3
    fitted_vals[:,sp,s] = fit_nb(sim_ds_sampled[:,sp,s])
  end
end

fitted_vals[1,1,:]
# check for correlation
plot(fitted_vals[1,1,:], fitted_vals[1,2,:], st = :scatter)

#Calculate basic distances between sim vals and target vals
dists = zeros(Float64, 33, 2)
dists[:,1] = fitted_vals[1,1,:] .- k_obs[1]
dists[:,2] = fitted_vals[2,1,:] .- mean_obs[1]

#Calculate standardised euclidean distance
function euc(dists)
    sqrt(sum(dists.^2))
end

euc_dists = zeros(Float64, 33, 1)
for s in 1:33
  euc_dists[s] = euc(dists[s,:])
end

using Plots
histogram(euc_dists[:,1])

#Get posteriors
#Accept small percent

accept = findall(i -> i < 8, euc_dists) #change 6 after looking at histogram

#Print un-squeezed posteriors
using StatsPlots
density(Imme_activation[accept])
density(Immf_activation[accept])
density(fec_modulation[accept])
density(est_modulation[accept])


#Squeeze posteriors

R"abc_ <- abc::abc" #use R abc package
@rget abc_ #port function back into Julia

Na_a = abc_([0.3, 15], a, fitted_vals, 0.1, "loclinear")
Na_b = abc_([0.3, 15], b, fitted_vals, 0.1, "loclinear")
Na_c = abc_([0.3, 15], c, fitted_vals, 0.1, "loclinear")

Al_a = abc_([0.3, 15], a, fitted_vals, 0.1, "loclinear")
Al_b = abc_([0.3, 15], b, fitted_vals, 0.1, "loclinear")
Al_c = abc_([0.3, 15], c, fitted_vals, 0.1, "loclinear")

density(Na_a[:adj_values])
density(Na_b[:adj_values])
density(Na_c[:adj_values])

# To Do
# Need to use a better distance measure; maybe Mahalanobis?
# Need to check the 'squeeze' procedure
# Need to make runable for all three species
# Need to have the fitting data being put straight into a CSV during fitting process to avoid data loss


# Notes:
# Using the RCall#master branch
