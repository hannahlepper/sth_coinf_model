###### STH co-infection model ######

using Distributions #Package contaning negative binomial distribution

###### Making NegativeBinomial faster ######
#Courtesy of my partner Mike who does programming all the time

#struct NegativeBinomial{T}
#  p0::T
#  cdf::Vector{T}
#end

#function NegativeBinomial{T}(r::Real,p::Real) where T
#  ps = cdf.(Distributions.NegativeBinomial(r,p), 0:40)
#  NegativeBinomial{T}(ps[1], ps)
#end

###### global parameters ######

#It is faster to define them as constants but we can't overwrite them
const ts = 1/73 #time step in years; 5 days
const halflife = 8.1 #halflife of immunity in hosts in days
const av_age = 18.2 #Average age in years in population - used very roughly
const pc_dr = 8/1000 * ts #per capita death rate
const stool_samp = 0.054 #Stool sample used in measuring egg deposition


###### Worm species specific parameters ######

#names
sp = [:Na, :Al, :Tt]

#I am defining a data structure with the elements b to init_pars, which will all be of the same type (class)
#We will be able to assign a set of values to this structure (see below), and will be able to access
#elements in the structure using '.', e.g N_a.b
#T means type, e.g. Int, Float64.
#The struct is mutable so that parameters can vary if necessary

mutable struct Par{T}
  b                    ::T   # 1.rate of exposure per day - has to be v v low
  Imme_activation      ::T   # 2.to fit - activation of immunity per worm
  Immf_activation      ::T   # 3.as above but for anti fecundity immunity
  WfN                  ::T   # 4.max eggs per day per female worm
  mu_le                ::T   # 5.death rate of early larvae per day
  mu_ll                ::T   # 6.death rate of late larvae per day
  mu_adults            ::T   # 7.death rate of adult worms per day
  M_le                 ::T   # 8.rate of maturation of early larvae per day
  M_ll                 ::T   # 9.rate of maturation of late larvae per day
  pool_egg_loss        ::T   # 10.rate of loss of eggs from field per day
  pool_egg_maturation  ::T   # 11.rate of maturation of eggs in field oer day
  pool_infectives_loss ::T   # 12.rate of loss of infective larvae from field per day
  dens_effect          ::T   # 13.density dependent effect on per worm fecundity
  k                    ::T   # 14.k of negative binomial distribution (mean, k)
  weightings           ::T   # 15.relative weight of worms
  est_modulation       ::T   # 16.modulation of antiestablishment immunity
  fec_modulation       ::T   # 17.modulation of antifecundity immunity
  init_mean            ::T   # 18.mean of initial distribution of worms in pop.
end

#Inputting parameters for each species
N_a = Par{Float64}(
    1.9e-10 * 365 * ts,
    0.001, 0.001, #1,2,3
    20000 * 365 * ts, 0.067 * 365 * ts, 0.00182 * 365 * ts, #4,5,6
    0.00182 * 365 * ts, 0.07 * 365 * ts, 0.0467 * 365 * ts, #7,8,9
    0.011 * 365 * ts, 0.11 * 365 * ts, 0.15 * 365 * ts, #10,11,12
    0.019, 0.23, 0.037, 0, 0, 2) #13,14,15,16,17,18


A_l = Par{Float64}(
    1.65e-11 * 365 * ts,
    0.001, 0.001, #1,2,3
    200000 * 365 * ts, 0.067 * 365 * ts, 0.00183 * 365 * ts, #4,5,6
    0.00183 * 365 * ts, 0.10 * 365 * ts, 0.0714 * 365 * ts, #7,8,9
    0.0085 * 365 * ts, 0.0286 * 365 * ts, 0.03 * 365 * ts, #10,11,12
    0.00425, 0.57, 1, 0, 0, 0.5) #13,14,15,16,17,18

T_t = Par{Float64}(
    3.8e-12 * 365 * ts,
    0.001, 0.001, #1,2,3
    20000 * 365 * ts, 1, 0.00182 * 365 * ts, #4,5,6
    0.00182 * 365 * ts, 0.4 * 365 * ts, 0.0133 * 365 * ts, #7,8,9
    0.00192 * 365 * ts, 0.0286 * 365 * ts, 0.05 * 365 * ts, #10,11,12
    0.001, 0.35, 0.0148, 0, 0, 0.5) #13,14,15,16,17,18

#We keep these three Pars structs in an array, accessing as SpPars[1] etc
SpPars = [N_a, A_l, T_t]

##### Host attributes data structures #####

#Each host has infections attributes per species in a data structure:

struct Infection{T}
    Imme    ::T #Anti-establishment immunity strength
    Immf    ::T #Anti-fecundity immunity strength
    PEL     ::T #Pre establishment larvae
    EL      ::T #Established larvae
    AW      ::T #Adult worms
    EOut    ::T #Egg output
end

#We define a method so that we can easily make an empty structure.
Infection{T}() where T = Infection{T}(0,0,0,0,0,0)

##### Pool attributes data structre #####

struct Soil{T}
    PIS ::T #Pre-infective stages
    IS  ::T #Infective stages
end

Soil{T}() where T = Soil{T}(0,0)

##### model functions ######

#arguements are: individual, worm burden, infective stages in pool, halflife, parameters
function update_Infection(i, WB, IS, halflife, p)
    #New infections
    mean = p.b * IS
    if mean > 0
        #exposure = float(rand(NegativeBinomial{Float64}(mean + float(1e-6), p.k)))
        exposure = float(rand(NegativeBinomial(mean, p.k)))
    else
        exposure = 0.0
    end
    newPEL = ((1-p.mu_le) * i.PEL) + exposure
    PEL = newPEL * (1-p.M_le)

    #Anti establishment immunity. We assume this is happening sequentially, using WB from the
    #previous time step - we are not taking the AW into account yet as we have not yet calculated them
    activation = PEL * p.Imme_activation
    modulation = exp(-(p.est_modulation * WB))
    Imme = ((0.5^(1/halflife) * i.Imme) + activation) * modulation

    #New established larvae
    newEL = ((1-p.mu_ll) * i.EL) + (p.M_le * newPEL * exp(-Imme))
    EL = newEL * (1-p.M_ll)

    #Anti fecundity immunity
    activation = EL * p.Immf_activation
    modulation = exp(-(p.fec_modulation * WB))
    Immf = ((0.5^(1/halflife) * i.Immf) + activation) * modulation

    #New adults
    AW = ((1-p.mu_adults) * i.AW) + (p.M_le * newEL)

    #New eggs - we are using the current time step's AW but the previous timestep's
    #WB in Imme and Immf as we assume they were altered only once per timestep
    modulation = exp(-(Immf + (p.dens_effect * AW/2)))
    EOut = AW/2 * p.WfN * modulation

    Infection{Float64}(Imme, Immf, PEL, EL, AW, EOut)
end

function update_WBs(pop, pars, n_hosts)
  [sum([x.AW for x in pop[i,:]] .* [p.weightings for p in pars]) for i in 1:n_hosts]
end

#Arguments: soil, population of infections, parameters
function update_pool(S, pop, p)
    Eggs = sum([x.EOut for x in pop])
    PIS = ((1 - (p.pool_egg_loss + p.pool_egg_maturation)) * S.PIS) + Eggs
    IS = ((1-p.pool_infectives_loss) * S.IS) + (p.pool_egg_maturation * S.PIS)
    Soil{Float64}(PIS, IS)
end

#Birth and death process - select all individuals over 80 years old and randomly select from the rest
#of the population - reset ages and infections so population size remains constant
function reset_inds_sys(population, ages, pc_dr)
  for i = 1:length(ages)
    if ages[i] > 80 || rand(Binomial(1, pc_dr)) == 1
      population[i, 1:3] .= Infection{Float64}()
      ages[i] = 0
    end
  end
  return population, ages
end

function update_ages(ages, ts)
    [a += ts for a in ages]
end


#Set up initial arrays etc. Can return multiple objects.
function SystemSetUp(n_hosts, pars, av_age)
    #Initialise ages, roughly scattered around average age, avoiding negative ages using abs
    ages = [abs(rand(Normal(av_age, av_age))) for i in 1:n_hosts]

    #Initialise pool with eggs and infective stages
    Pool = [Soil{Float64}(100, 1000000) for sp in 1:3]

    #Initialise worms with n. binom draw for adults
    init_worms(r, p) = Infection{Float64}(0, 0, 0, 0, float(rand(NegativeBinomial(r,p))), 0)
    pop_infections = [init_worms(p.init_mean, p.k) for i in 1:n_hosts, p in pars]

    #Initialise worms burdens based on initial adult burdens
    WBs = update_WBs(pop_infections, pars, n_hosts)
    return ages, Pool, pop_infections, WBs
end

function run_mod(n_hosts, pars, ts, halflife, pop_infections, Pool, ages, WBs, pc_dr)

    #Update ages and remove some individuals
    ages = update_ages(ages, ts)
    pop_infections, ages = reset_inds_sys(pop_infections, ages, pc_dr)

    #Species specific calculations
    for sp in 1:3
        p = SpPars[sp]
        IS = Pool[sp].IS

        #Per host calculations
        for i in 1:n_hosts
            pop_infections[i, sp] = update_Infection(pop_infections[i, sp], WBs[i], IS, halflife, p)
        end

        #Update pool
        Pool[sp] = update_pool(Pool[sp], pop_infections[:, sp], p)
    end

    #Update worm burdens
    WBs = update_WBs(pop_infections, SpPars, n_hosts)

    return ages, pop_infections, Pool, WBs
end

function main(n_runs, n_hosts, pars, av_age, ts, halflife, pc_dr, stool_samp)

  #Set arrays up
  ages, Pool, pop_infections, WBs = SystemSetUp(n_hosts, pars, av_age)

  #For storing summary statistics
  EC = zeros(Float64, n_runs, 3)
  prevs = zeros(Float64, n_runs, 3)
  rec = 0#n_runs - 6

  #Loop through the runs
  for r in 1:n_runs
    ages, pop_infections, Pool, WBs = run_mod(n_hosts, pars, ts, halflife, pop_infections, Pool, ages, WBs, pc_dr)
    if r > rec
      for sp in 1:3
        #Take mean egg output from only infected individuals
        eggs = filter(!iszero, [x.EOut for x in pop_infections[:,sp]])
        EC[r - rec,sp] = mean(eggs)/(365 * ts) * stool_samp
        prevs[r - rec,sp] = count(i -> i > 0.0, x.AW for x in pop_infections[:,sp])/n_hosts
      end
    end
  end
#  EC = [log(mean(EC[:,sp])) for sp in 1:3]
#  prevs = [mean(prevs[:,sp]) for sp in 1:3]
  return EC, prevs, pop_infections
end

#Currently taking about 10 secs, give or take
@time EC, prevs, pop = main(10000, 1000, SpPars, 18.2, ts, halflife, pc_dr, stool_samp)

using Plots

plot(1:10000, EC)
plot(1:10000, prevs)
