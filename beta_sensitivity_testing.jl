# #Beta sensitivity testing
# In this script I do a basic sensitivity analysis of the parameters beta and k.

include("coinf_model.jl")

# Values of beta for each species, based on trial and error:
b_N_a_full = linspace(5e-11, 5e-9, 10)
b_A_l_full = linspace(5e-11, 5e-8, 10)
b_T_t_full = linspace(1e-9, 5e-8, 10)

# Values of k for each species, based on fitted values of k from the Sri Lanka data:
k_N_a = [0.01, 0.2, 0.4]
k_A_l = [0.01, 0.2, 0.4]
k_T_t = [0.11, 0.21, 0.31]

# Get all combinations of the beta and k vaules for testing
N_a_tup = collect(Base.product(b_N_a_full,k_N_a))
A_l_tup = collect(Base.product(b_A_l_full,k_A_l))
T_t_tup = collect(Base.product(b_T_t_full,k_T_t))

# Function to do the testing
sp_test = function(p_var, p, sp)

  eggs = zeros(Float64, 2000, length(p_var))

  for i in 1:length(p_var)

    p[sp].b = p_var[i][1]
    p[sp].k = p_var[i][2]

    #rr is 0 (if set run_record = 0 get record of whole run); a is ages
    rr, final_record, a = main(1000, 2000)
    eggs[:,i] = final_record[:EOut][:,sp]
  end

  return eggs
end

# Each set of testing takes just under 2 minutes
@time N_a_b_sensitivity = sp_test(N_a_tup, SpPars, 1)
@time A_l_b_sensitivity = sp_test(A_l_tup, SpPars, 1)
@time T_t_b_sensitivity = sp_test(T_t_tup, SpPars, 1)

# Need to scale the eggs to match the raw data. Egg output from model is all eggs that can be output
# for 5 days (timestep = 1/73). Stool sample size if 54mg (0.054g). I assume 100g of stool per day,
# therefore scale to per day (i/5) and then to correct concentration (* 0.00054)
scale_N_a_sensitivity = map(i -> i/5 * 0.00054, N_a_b_sensitivity)
scale_A_l_sensitivity = map(i -> i/5 * 0.00054, A_l_b_sensitivity)
scale_T_t_sensitivity = map(i -> i/5 * 0.00054, T_t_b_sensitivity)

# Using the fitdistr function from the MASS package in R - it's just easier. I have adapted the
# function slightly, to handle all or almost all 0 data, see fitNB.R file.
using RCall
R"source('fitNB.R')"
@rget fit_nb

# Get mu and k (or mu and size) from the simulated data. Firstly I fit the data using the function
# above, I then prepare it for plotting and plot away.

# The first plots are model sensitivity to beta and k for N. americanus.
# For the simulated k - it's not completely surprising that the k that is simulated is a bit
# different to the k that is given as a parameter, because the parameter k is fitted from EPG
# and used in defining infective larvae acquisition, and the simulated k is from the egg output.
# Target k for N. americanus is 0.08.
using Plots
N_a_size = hcat([fit_nb(round(scale_N_a_sensitivity[:,x]))[1] for x in 1:30])
N_a_size = hcat(N_a_size[1:10], N_a_size[11:20], N_a_size[21:30])
plot(b_N_a_full, N_a_size,
      title = "Simulated N. americanus EPG k ~ beta", label = ["k = 0.01", "k = 0.2", "k = 0.4"])
hline!([0.08], line = (1, :dash, 0.6, :black), label = :target)
#

# Simulated mu for N. americanus look sensible, target N. a. mu is 6.9.
N_a_mu = hcat([fit_nb(round(scale_N_a_sensitivity[:,x]))[2] for x in 1:30])#src
N_a_mu = hcat(N_a_mu[1:10], N_a_mu[11:20], N_a_mu[21:30])#src
plot(b_N_a_full, N_a_mu,
      title = "Simulated N. americanus EPG mu ~ beta", label = ["k = 0.01", "k = 0.2", "k = 0.4"])
hline!([6.9], line = (1, :dash, 0.6, :black), label = :target)
#

# Simulated Ascaris k. Target is 0.1.
A_l_size = hcat([fit_nb(round(scale_A_l_sensitivity[:,x]))[1] for x in 1:30])#src
A_l_size = hcat(A_l_size[1:10], A_l_size[11:20], A_l_size[21:30])#src
plot(b_A_l_full, A_l_size,
      title = "Simulated Ascaris EPG k ~ beta", label = ["k = 0.01", "k = 0.2", "k = 0.4"])
hline!([0.1], line = (1, :dash, 0.6, :black), label = :target)
#

# Simulated Ascaris mu. Target is 249.6
A_l_mu = hcat([fit_nb(round(scale_A_l_sensitivity[:,x]))[2] for x in 1:30])#src
A_l_mu = hcat(A_l_mu[1:10], A_l_mu[11:20], A_l_mu[21:30])#src
plot(b_A_l_full, A_l_mu,
      title = "Simulated Ascaris EPG mu ~ beta", label = ["k = 0.01", "k = 0.2", "k = 0.4"])
hline!([249.6], line = (1, :dash, 0.6, :black), label = :target)
#

# Simulated Trichuris k; target is 0.21
T_t_size = hcat([fit_nb(round(scale_T_t_sensitivity[:,x]))[1] for x in 1:30])#src
T_t_size = hcat(T_t_size[1:10], T_t_size[11:20], T_t_size[21:30])#src
plot(b_T_t_full, T_t_size,
      title = "Simulated Trichuris EPG k ~ beta", label = ["k = 0.11", "k = 0.21", "k = 0.31"])
hline!([0.21], line = (1, :dash, 0.6, :black), label = :target)
#

# Simulated Trichuris mu; target 38.8
T_t_mu = hcat([fit_nb(round(scale_T_t_sensitivity[:,x]))[2] for x in 1:30])#src
T_t_mu = hcat(T_t_mu[1:10], T_t_mu[11:20], T_t_mu[21:30])#src
plot(b_T_t_full, T_t_mu,
      title = "Simulated Trichuris EPG mu ~ beta", label = ["k = 0.11", "k = 0.21", "k = 0.31"])
hline!([38.8], line = (1, :dash, 0.6, :black), label = :target)
