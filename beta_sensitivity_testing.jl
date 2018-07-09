#Beta sensitivity testing


#aim data:
|sp|bna|bal|btt|kna|kal|ktt|
|ECs...|
Head to R to use pscl and plotting, plot mu and k of zinbs
----------------------------------------
#Choose 10 values for each beta
#For each species, run through 20 beta of interest (BOI)
#for each of the 3 other values of betas and the kOI.
#For a total of: 20 * 9 = 180 runs
#For each of the three species, 180 * 3 = 540 runs
#(540 * 15)/60/60
#And expect it to take around 2.5 hours.

using CSV
using DataFrames

#Beta run function
b_vary_run = function(bs, p, sp)
  eggs = zeros(Float64, 5000, length(bs))

  for i in 1:length(bs)
    p[sp].b = bs[i]
    eggs_all = main(2000, 5000, p, 18.2, ts,
                          halflife, pc_dr, stool_samp)
    eggs[:,i] = eggs_all[:,sp]
  end

  return eggs
end

b_N_a_full = linspace(5e-12, 5e-10, 10)
b_A_l_full = linspace(1e-13, 1e-11, 10)
b_T_t_full = linspace(1e-14, 1e-12, 10)

b_N_a_short = linspace(5e-12, 5e-10, 3)
b_A_l_short = linspace(1e-13, 1e-11, 3)
b_T_t_short = linspace(1e-14, 1e-12, 3)

k_N_a = [0.17, 0.27, 0.37]
k_A_l = [0.24, 0.34, 0.44]
k_T_t = [0.11, 0.21, 0.31]

Na_par = hcat(repeat(b_A_l_short, inner = 9),
              repeat(b_T_t_short, inner = 3, outer = 3),
              repeat(k_N_a, outer = 9))
Al_par = hcat(repeat(b_N_a_short, inner = 9),
              repeat(b_T_t_short, inner = 3, outer = 3),
              repeat(k_A_l, outer = 9))
Tt_par = hcat(repeat(b_N_a_short, inner = 9),
              repeat(b_A_l_short, inner = 3, outer = 3),
              repeat(k_T_t, outer = 9))


@time test = b_vary_run(b_N_a_short, SpPars, 1)

writecsv("b_long.csv", hcat(b_N_a_full, b_A_l_full, b_T_t_full))
writecsv("b_short.csv", hcat(b_N_a_short, b_A_l_short, b_T_t_short))
writecsv("k.csv", hcat(k_N_a, k_A_l, k_T_t))

sp_test = function(bs, vp, p, sp_order)
  n = length(bs)
  m = length(vp[:,1])

  eggs = zeros(Float64, 5000, n * m)

  for i in 1:m

    p[sp_order[2]].b = vp[i,1]
    p[sp_order[3]].b = vp[i,2]
    p[sp_order[1]].k = vp[i,3]

    eggs[:, ((i * n) - (n-1)) : (i * n)] = b_vary_run(bs, p, sp_order[1])
  end

  return eggs
end

@time sp_test(b_N_a_short[1], Na_par, SpPars, [1,2,3])

@time sensitivity_data_N_a = sp_test(b_N_a_full, Na_par, SpPars, [1,2,3])
open(file -> serialize(file, sensitivity_data_N_a), "N_a.jls", "w")

@time sensitivity_data_A_l = sp_test(b_A_l_full, Al_par, SpPars, [2,1,3])
open(file -> serialize(file, sensitivity_data_A_l), "A_l.jls", "w")

@time sensitivity_data_T_t = sp_test(b_T_t_full, Tt_par, SpPars, [3,1,2])
open(file -> serialize(file, sensitivity_data_T_t), "T_t.jls", "w")

writecsv("N_a.csv", sensitivity_data_N_a)
writecsv("A_l.csv", sensitivity_data_A_l)
writecsv("T_t.csv", sensitivity_data_T_t)

#Shashi's Text Pass might be faster for this
# Tune betas to figure out realistic ranges for them

#----------
b_N_a = linspace(5.2e-11, 5.5e-11, 3)
EC_Na_data = zeros(Float64, 3)
prev_Na_data = zeros(Float64, 3)

@time for i in 1:3
  N_a.b = b_N_a[i]
  EC, prevs = main(2000, 5000, SpPars, 18.2, ts, halflife, pc_dr, stool_samp)

  EC_Na_data[i] = EC[2000, 1]
  prev_Na_data[i] = prevs[2000,1]
end

plot(b_N_a, EC_Na_data)
plot(1:2000, EC)
#5.5e-11
EC_Na_data
#-------------
b_A_l = linspace(9e-13,2e-12, 10)
EC_Al_data = zeros(Float64, 10)
prev_Al_data = zeros(Float64, 10)

for i in 2
  A_l.b = b_A_l[i]
  EC, prevs = main(2000, 5000, SpPars, 18.2, ts, halflife, pc_dr, stool_samp)

  EC_Al_data[i] = EC[2000,2]
  prev_Al_data[i] = prevs[2000,2]
end

plot(b_A_l, EC_Al_data)
EC_Al_data
plot(1:2000, EC)
#8.05e-12

#-------------------
b_T_t = linspace(8e-14, 2e-13, 3)
EC_Tt_data = zeros(Float64, 3)
prev_Tt_data = zeros(Float64, 3)

for i in 1:3
  T_t.b = b_T_t[i]
  EC, prevs = main(2000, 5000, SpPars, 18.2, ts,
                        halflife, pc_dr, stool_samp)
  if isnan(EC[2000, 3])
    EC_Tt_data[i] = 0
    prev_Tt_data[i] = 0
  else
    EC_Tt_data[i] = EC[2000, 3]
    prev_Tt_data[i] = prevs[2000, 3]
  end
end

plot(float(b_T_t), EC_Tt_data)
EC
plot(1:2000, EC)
