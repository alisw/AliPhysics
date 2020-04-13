# PostProcessing repository: PWGLF-RESONANCE #

This directory is intended for the post-processing code. Basically all the macros that people use after they get the results of their analysis tasks (fitting peaks, calculating systematic uncertainties, fitting spectra, …). 


## Folder structure of PostProcessing directory##

— 'RESONANCE/PostProcessing/Name_Of_Particle(_and_Mass)' for your analysis task and AddTask macro
— 'RESONANCE/PostProcessing/Name_Of_Particle(_and_Mass)/macros_Collisions_System,Energy' for post-processing macros that can reproduce your analysis.


e.g. 'RESONANCE/PostProcessing/Xi1530’ for your analysis task and AddTask macro
     'RESONANCE/PostProcessing/Xi1530/macros_pp13TeV’ for post-processing macros that can reproduce your analysis.

Together with your codes a description on how to run it should be put inside a `README.md` file in that directory.

