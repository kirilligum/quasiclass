----------------------------------------
----------------------------------------
------ Atomic units (from http://en.wikipedia.org/wiki/Atomic_units )
------ electron mass = elementary charge = reduced Planck's constant = Coulomb's constant = 1
----------------------------------------
---- Fundamental atomic units in SI
-- mass  
electron_rest_mass =   9.10938291e-31; -- kg; 
-- charge  
elementary_charge =  1.602176565e-19; -- C; 
-- angular momentum  
reduced_Plancks_constant =  1.054571726e-34; --  J*s; 
-- electric constant 
Coulomb_force_constant =   8.9875517873681e9; -- kg*m^3*s^(-2)*C-^(-2); 
----------------------------------------
---- Physical constants in atomic units
fine_structure_constant = 1/137;
classical_electron_radius = 5.32e-5;
proton_mass = 1836;
hbar = 1.0;
boltzmann_constant = 1.0; -- atomic units
--boltzmann_constant = 3.166811429e-6; -- E_h/K ; from http://en.wikipedia.org/wiki/Boltzmann_constant
----------------------------------------
---- Physical constants in SI
si_boltzmann_constant = 1.3806488e-23; -- J/K ; from http://en.wikipedia.org/wiki/Boltzmann_constant
----------------------------------------
---- Derived atomic units in SI
-- length  
--bohr_radius =    5.2917720859e-11; -- m  ; 
--bohr_radius = reduced_Plancks_constant/(electron_rest_mass*speed_of_light*fine_structure_constant); -- m  ; 
bohr_radius = reduced_Plancks_constant^2/(Coulomb_force_constant*electron_rest_mass*elementary_charge^2)
-- energy 
--hartree_energy =     4.35974417e-18; -- J ;
--hartree_energy = fine_structure_constant^2*electron_rest_mass*speed_of_light*speed_of_light; -- J ;
hartree_energy = electron_rest_mass*elementary_charge^4*Coulomb_force_constant^2/reduced_Plancks_constant^2; -- J ;
-- time
--atomic_time =    2.418884326505e-17; -- s; 
atomic_time = reduced_Plancks_constant/hartree_energy; -- s; 
-- velocity
--atomic_velocity = 2.1876912633e6; -- m*s^(-1); 
--atomic_velocity = fine_structure_constant*speed_of_light; -- m*s^(-1); 
atomic_velocity = bohr_radius*hartree_energy/reduced_Plancks_constant; -- m*s^(-1); 
-- force
--atomic_force = 8.2387225e-8; --  N;
atomic_force = hartree_energy/bohr_radius; --  N;
-- temperature
atomic_temperature = 3.1577464e05; -- K  ;
--atomic_temperature = hartree_energy/si_boltzmann_constant; -- K  ;
-- pressure
--atomic_pressure = 2.9421912e13; -- Pa  ;
atomic_pressure = hartree_energy/(bohr_radius^3); -- Pa  ;
-- electric field
--atomic_electric_field = 5.14220652e11; --  V*m^(-1);
atomic_electric_field = hartree_energy/(elementary_charge*bohr_radius); --  V*m^(-1);
-- electric dipole moment
--atomic_electric_dipole_moment = 8.47835326e-30; --  C*m;
atomic_electric_dipole_moment = elementary_charge*bohr_radius; --  C*m;
----------------------------------------
---- Conversions to atomic units
--ev = 0.03674932; -- (*"hartrees"*)
ev = elementary_charge/hartree_energy; -- J/C*e (*"hartrees"*)
--speed_of_light = 299792458; -- m*s^(-1)
speed_of_light = elementary_charge^2*Coulomb_force_constant/(reduced_Plancks_constant*fine_structure_constant); -- m*s^(-1) in SI units
--i_cm = atomic_time/speed_of_light; -- reciprocal centimetre to Hz to atomic_time^(-1)
i_cm = 1.23984e-4*ev; -- reciprocal centimetre to Hz to atomic_time^(-1)

---------------------------------------- model
-- total system parameters
m_p = 1.0;
m=m_p;
--e1 = 0.0*ev; -- energy of donor
e1 = 0.5*ev; -- energy of donor
e2 = -0.5*ev; -- energy of acceptor
--v12 = 0.17*ev; -- coupling beteen donor and accetor
--v12 = 0.03*ev; -- coupling beteen donor and accetor
v12 = 1.0; -- for spin bosson
rp0 = -0.5; -- minimum of ground state potential
--rp1 = 0; -- minimum of donor potential
rp1 = -0.5; -- minimum of donor potential
rp2 = -0.5; -- minimum of acceptor potential
w0 = 1.0;
w1 = 1.0;
w2 = 1.0;
--w0 = 0000*i_cm;
--w1 = 0000*i_cm;
--w2 = 0000*i_cm;
w_c = 2.5;
w_max_factor = 5;
w_max = w_max_factor*w_c;
--n_modes = 10;
n_modes = 100;
--dw = w_c;
dw =w_max/n_modes;
--dw = 1.0e-1*w0;
--dw = 2.0e-3*w0;
--eta = 5e0; -- system-bath coupling (dimensionless)
--eta = 0.0;
alpha =0.09;-- Makri
--eta = 0.0;
eta = alpha*3.1415926*0.5;
--eta = 12*math.pi; -- system-bath coupling
--eta = alpha*M_PI_2;

----------------------------------------
-- total system parameters
--m_p = 1.0*proton_mass; -- mass of hydrogen
--m=m_p;
----e1 = 0.0*ev; -- energy of donor
--e1 = 0.0*ev; -- energy of donor
--e2 = 0.0*ev; -- energy of acceptor
----v12 = 0.17*ev; -- coupling beteen donor and accetor
--v12 = 0.03*ev; -- coupling beteen donor and accetor
--rp0 = -0.5; -- minimum of ground state potential
----rp1 = 0; -- minimum of donor potential
----rp1 = -0.5; -- minimum of donor potential
--rp1 = -0.0; -- minimum of donor potential
--rp2 = -0.5; -- minimum of acceptor potential
--w0 = 3000*i_cm;
--w1 = 3000*i_cm;
--w2 = 3000*i_cm;
----w0 = 0000*i_cm;
----w1 = 0000*i_cm;
----w2 = 0000*i_cm;
--w_c = 600*i_cm;  -- bath cutoff frequency
--w_max_factor = 1;
--w_max = w_max_factor*w_c;
--n_modes = 0;
----n_modes = 100;
--dw = w_c;
----dw =w_max/n_modes;
----dw = 1.0e-1*w0;
----dw = 2.0e-3*w0;
----eta = 5e0; -- system-bath coupling (dimensionless)
--eta = 0.0;
----eta = 12*math.pi; -- system-bath coupling
----eta = alpha*M_PI_2;

--alpha = 0.09;
--wc = 2.5;
--beta = 0.1 ;
--beta = 5 ;
--delta = 1;
--eps = 0;
--end_time = 20 ;
--pp = 0.0;
--qp = 1.0;
--n1 = 1.0;
--n2 = 0.0;

----------------------------------------
-- initial conditions
--pp=0.0;
--qp=rp0;
--n1=1.0;
--n2=1.0-n1;
--temperature = 300/atomic_temperature; -- K to atomic_temperature
----temperature = 300;
--beta=1.0/boltzmann_constant/temperature;

---------------------------------------- model
-- initial conditions
pp=0.0;
qp=rp0;
n1=1.0;
n2=1.0-n1;
beta=5.0;
--temperature = 300/atomic_temperature; -- K to atomic_temperature
temperature = 1/(boltzmann_constant*beta);
--temperature = 300;

----------------------------------------
-- propagation parameters
end_time = 15e0; -- sec to atomic_time
--end_time = 10; -- sec to atomic_time
--end_time = 1e1*1e-15/atomic_time; -- sec to atomic_time
--end_time = 2e1*1e-15/atomic_time; -- sec to atomic_time
--n_times = 5; 
n_times = 100; 
time_step = end_time/n_times;

----------------------------------------
-- method parameters
seed = 1;
--n_bin_width = 0.1;
n_shift = 0.35;
bin_start = 0.35;
bin_end = 0.35;
n_trajs =1;
--n_trajs =2;
--n_trajs = 20;
--n_trajs = 500;
--n_trajs = 2000;
print_int_states = 0;
print_energies = 0;


-- dofile "config.lua"
-- print(sys_k2,"\n",sys_ki,"\n",sys_k4,"\n",sys_m,"\n",sys_hb,"\n",s0)
-- sed -i "s@./data/ki_........e...@$(find ./data/ki* -type d)@g" *.plt
-- rm -rf data; time ./cmqhd nano.lua; sed -i "s@./data/ki_........e...@$(find ./data/ki* -type d)@g" *.plt;killall gnuplot; gnuplot -p -e 'load "pqpss.plt"'
--lua -e 'dofile "config.lua"; print( "\n bohr_radius = ",bohr_radius, "\n hartree_energy = ",hartree_energy, "\n ev = ", ev, "\n atomic_velocity = ", atomic_velocity, "\n speed_of_light = ",speed_of_light, "\n w0 = ", w0, "\n end_time = ", end_time,"\n temperature = ",temperature )'
