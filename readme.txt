Grossowicz et al. 2016 MatLab scripts
25 Jul 2016

Content:
(1) Basic model
(2) Excretion model
(3) Sensetivity analysis
(4) Monte Carlo
(5) Genetic Algorithem


BASIC MODEL
In the "Basic model" folder, run: Pro_LowN_model.m. 
The other files are functions used by the .m file.

EXCRETION MODEL
In the "Excretion model" folder, run: Pro_LowN_model_excretion.m. 
The other files are functions used by the .m file.

SENSETIVITY ANALYSIS
In the "Sensetivity analysis" folder, run: multiplot_fig4.m. 
The other files are functions used by the .m file, the .mat file 
has the values of the parameter's ranges, and the model 
forPCA9312_LowN_axenic.m which is operated directly from 
multiplot_fig4.m file.

MONTE CARLO
In the "Monte Carlo" folder, run: MonteCarlo_Csat_h.m. 
The other files are functions used by the .m file, and the model 
forMonteCarlo_LowN_Csat_h.m which is operated directly from 
multiplot_fig4.m file.

GENETIC ALGORITHM
In the "Genetic Algorithem" folder, run: run.m. calc_csat.m and 
Pro_Csat.m are the model itself. The other functions: ga.m, 
fitness.m, gene_to_values.m, and squared_error.m are related 
to the GA analysis and were written by Dr. George Bezerra, 
taken from: http://matlabgeeks.com/tips-tutorials/ode-tips-tutorials/modeling-with-odes-in-matlab-part-4b/. 