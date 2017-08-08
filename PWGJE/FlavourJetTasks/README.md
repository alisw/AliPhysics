# Pico V0 package

- Package for V0-jet analysis in pp, p-Pb/Pb-p and Pb-Pb


## output objects in delta AOD

### Utils
- AliPicoBase: nanespace for constants, such as the V0 mass, which will be used in the analysis

### Event level info
- AliPicoHeaderV0: delta AOD header which strores event level info (eg. vtx, EP, trg...)
- AliPicoHeaderJE: AliPicoHeaderV0 + jet background density

### V0 info
- AliPicoV0: basic class which stores V0 info
- AliPicoV0RD: AliPicoV0 + info used for data analysis
- AliPicoV0MC: AliPicoV0 + info used for MC analysis

### Jet info
- AliPicoJet: stores info of jet kinematics, area and leading particle pT


## Analysis tasks

### V0 maker
- AliAnalysisTaskSEPicoV0Maker
  - Inherits from AliAnalysisTaskSE
  - V0 candidate selection for both data and MC
  - Add task: macros/AddTaskPicoV0Maker.C

- AliAnalysisTaskSEPicoV0MakerMC
  - Inherits from AliAnalysisTaskSE
  - V0 candidate selection for only MC + std cuts
  - Add task: macros/AddTaskPicoV0MakerMC.C

### Filters for generating delta AOD
- AliAnalysisTaskSEPicoV0Filter
  - Inherits from AliAnalysisTaskSE
  - Fills delta AOD with V0 branch alone
  - Add task: macros/AddTaskPicoV0Filter.C

- AliAnalysisTaskEmcalJetV0Filter
  - Inherits from AliAnalysisTaskEmcalJet
  - Fills delta AOD with jets and/or V0s
  - Add task: macros/AddTaskEmcalJetV0Filter.C

### Efficiency tool
- AliAnalysisTaskEmcalJetV0CF
  - Inherits from AliAnalysisTaskEmcalJet
  - Tool for getting efficiency of V0s matched with jets
  - Add task: macros/AddTaskEmcalJetV0CF.C


## Cfg example

### Local analysis
- macros/AnalysisTrainEMCalJetV0Local.C
  - Example for local cfgs of AliAnalysisTaskEmcalJetV0Filter
