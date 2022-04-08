# Validation scripts
This folder contains the script that was used to validate the transmission network reconstruction
done by [Phyloscanner](https://github.com/BDI-pathogens/phyloscanner) and [TNet](https://github.com/sauravdhr/tnet_python).

to run the script you should have the following folders structure:
```
~home_folder
- phyloscanner
- projects
-- tnet
-- network_favites
```
phyloscanner and tnet for for the tools, and network_favites for the script and FAVITES input files. Inside network_favites
there should be Favites_inputs folder with following structure:

![Favites_inputs content](/validation_scripts/FAVITIS_inputs.jpg)

There are three steps to run the tool:
- data preparation(converting the names from input tree to run the tools)
- running tools
- running analyses of the tools based on outputs from the previous step
Example of commands that need to be executed in order:
```
python3 favites_validation.py Favites_inputs/ data
python3 favites_validation.py Favites_inputs/ tools
python3 favites_validation.py Favites_inputs/ analyze
```
After running the last step the output to console will be:
```
Category log_t10SI_contemp. Total samples 100
tnet sensitivity = 0.296,specificity = 0.296,f1 = 0.296
phylo sensitivity = 0.016,specificity = 0.295,f1 = 0.030
Category origSIR_diff. Total samples 71
tnet sensitivity = 0.689,specificity = 0.589,f1 = 0.629
phylo sensitivity = 0.654,specificity = 0.830,f1 = 0.712
Category expSI_contemp. Total samples 100
tnet sensitivity = 0.502,specificity = 0.502,f1 = 0.502
phylo sensitivity = 0.414,specificity = 0.625,f1 = 0.481
Category origExpSIR_diff. Total samples 67
tnet sensitivity = 0.673,specificity = 0.520,f1 = 0.577
phylo sensitivity = 0.708,specificity = 0.823,f1 = 0.753
```

For some samples phyloscanner might not produce the output since it has only one host. The script skips such cases.
In 'results_time.txt' there should be detailed statistics for each sample for each category:

```
log_t10SI_contemp
sample_name,tnet_sensitivity,tnet_specificity,tnet_f1,phylo_sensitivity,phylo_specificity,phylo_f1
FAVITES_output_log_t10SI_contemp_T200_N100_E1_98,0.021,0.021,0.021,0.021,1.000,0.042
FAVITES_output_log_t10SI_contemp_T200_N100_E1_74,0.034,0.034,0.034,0.000,0.000,0.000
FAVITES_output_log_t10SI_contemp_T200_N100_E1_47,0.500,0.500,0.500,0.000,0.000,0.000
FAVITES_output_log_t10SI_contemp_T200_N100_E1_71,0.200,0.200,0.200,0.000,0.000,0.000
FAVITES_output_log_t10SI_contemp_T200_N100_E1_80,1.000,1.000,1.000,0.000,0.000,0.000
FAVITES_output_log_t10SI_contemp_T200_N100_E1_81,1.000,1.000,1.000,0.000,0.000,0.000
FAVITES_output_log_t10SI_contemp_T200_N100_E1_87,0.714,0.714,0.714,0.000,0.000,0.000
FAVITES_output_log_t10SI_contemp_T200_N100_E1_59,0.355,0.355,0.355,0.032,1.000,0.062
FAVITES_output_log_t10SI_contemp_T200_N100_E1_100,0.077,0.077,0.077,0.038,1.000,0.074
FAVITES_output_log_t10SI_contemp_T200_N100_E1_93,0.571,0.571,0.571,0.000,0.000,0.000
FAVITES_output_log_t10SI_contemp_T200_N100_E1_72,0.000,0.000,0.000,0.000,0.000,0.000
FAVITES_output_log_t10SI_contemp_T200_N100_E1_46,0.031,0.031,0.031,0.015,1.000,0.030
FAVITES_output_log_t10SI_contemp_T200_N100_E1_69,0.034,0.034,0.034,0.000,0.000,0.000
...
```



