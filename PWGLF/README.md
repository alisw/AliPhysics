# PWGLF analysis code repository

Welcome! This folder is hosting all the analysis code of the PWGLF.
In order to have analyses that are as much as possible reproducible,
in this folder we host both the analysis tasks and the post processing
macros. In the following you will find a couple of guidelines on how
to add your code in the right place here.
Remember: if you are not sure about where to put the code send an email
to your PAG coordinators or to the LF Offline responsible.

## Folder structure
As you can see the top level of the PWGLF folder contains the natural
division in PAGs and some additional service directories.
The service directories are intended for QA and common analysis related
code, typically your analysis code will live inside your PAG directory.

Inside your PAG directory you can find different ways of organising the
code. The typical strategy is grouping analyses on similar topics in
a common folder (e.g. all the nuclei analyses code go to the `NUCLEX/Nuclei`
directory, while the exotica ones go to the `NUCLEX/Exotica` folder).

> When you start your analysis is strongly suggested that you check if a
> task covering your topic is already present. If it is already there talk
> to the people responsible for the task! Usually joining forces is better
> than re-inventing something.

The best practice is to create a new folder in your PAG directory with the
name of the analysed particle/observable. For instance: the analysis on the 
Xi1530 production sits in the `RESONANCES/Xi1530`.
Inside the folder you should put your analysis task and your `AddTask` macro.

## Post-processing macros
The post-processing macros are vital to reproduce your analysis! You have to
make them public when your results goes for approval (preliminary and/or paper).
In the same folder where the task you used sits, create a subfolder with a name
reporting the system and the energy of the system you analysed. A couple of example:
- the nuclei production in pp at 13 TeV post-processing macros are in the folder `NUCLEX/Nuclei/NucleiPbPb/macros_pp13TeV`;
- the Xi1530 production in pp at 13 TeV post-processing macros are in the folder `RESONANCES/Xi1530/macros_pp13TeV`.

Together with your code a description on how to run it should be put inside a
`README.md` file in that directory.