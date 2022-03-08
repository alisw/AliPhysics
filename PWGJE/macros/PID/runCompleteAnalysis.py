import os
import sys
from helperFunctions import *

foldersToCreate = ['SummedSystematics', 'SingleSystematicResults', 'UEsubtractedJetResults', 'Efficiencycorrected', 'Efficiencycorrected/PDF']
folderEffCorrection = 'Efficiencycorrected'

if len(sys.argv) < 2:
  print('Provide name of analysis folder')
  exit()

analysisFolder=sys.argv[1]
systematicDay=sys.argv[2]

#### Read config ###
config = getConfig(analysisFolder)
config['analysisFolder'] = analysisFolder

#### Produce correction factor parameters ###
parameters = list(())
for spec in ('Electron', 'Muon', 'Pion', 'Kaon', 'Proton'):
  for charge in ('neg', 'pos'):
    configName = 'FS_Parameters_' + spec + '_' + charge
    parameters.append('-'.join(config[configName]))

arguments = {
  'effFile' : config['analysisFolder'] + "/Data/MCData/" + config['efficiencyFileNamePattern'].format('Jets_Inclusive'),
  'outputfile' : config['analysisFolder'] + "/Data/MCData/fastSimulationParameters" + "_" + config['MCRunName'] + ".root",
  'parameters' : parametersString
}
callRootMacro("FitCorrFactors", arguments)

#### Produce Corrections factors of full MC run
arguments = {
  'pathNameEfficiency' : config['analysisFolder'] + "/Data/MCData/" + config['efficiencyFileNamePattern'].format('Jets'),
  'outfileName' : config['pathMCCorrectionFile'],
}
callRootMacro("createFileForBbBCorrections", arguments) #TODO: Give observables, jetPts as argument, check if it can be merged with writeOutCorrectionFiles

#### Write correction files ###
for eff in config['MCVariationsEff']:
  varFolder="Eff" + eff + "_Res100/"
  arguments = {
    'effFilePath': config['analysisFolder'] + "/Data/MCData/" + varFolder + config['efficiencyFileNamePattern'].format('Jets'),
    'outFilePath' : config['analysisFolder'] + "/Data/MCData/",
    'addToFile' : varFolder
  }
  callRootMacro("writeOutCorrectionFiles", arguments) #TODO: Give observables, jetPts as argument (also below). Check if this triple call can be simplified
  
for res in config['MCVariationsRes']:
  varFolder="Eff100" + "_Res" + res + "/"
  arguments = {
    'effFilePath': config['analysisFolder'] + "/Data/MCData/" + varFolder + config['efficiencyFileNamePattern'].format('Jets'),
    'outFilePath' : config['analysisFolder'] + "/Data/MCData/",
    'addToFile' : varFolder
  }
  callRootMacro("writeOutCorrectionFiles", arguments)
  
for varFolder in config['MCVariationsLowPt']:
  arguments = {
    'effFilePath': config['analysisFolder'] + "/Data/MCData/" + varFolder + config['efficiencyFileNamePattern'].format('Jets'),
    'outFilePath' : config['analysisFolder'] + "/Data/MCData/",
    'addToFile' : varFolder
  }
  callRootMacro("writeOutCorrectionFiles", arguments)
    
#### Produce correction files for fast simulation ###  
arguments = {
  'effFilePath': config['analysisFolder'],
  'outFilePath' : config['analysisFolder'] + "/MCSystematicsFiles",
  'addToFile' : varFolder
}
callRootMacro("sysErrorsPythiaFastJet_new", arguments) #TODO: Give jetPtString etc.
    
#### Read systematic variables ###
systematics = getSystematicConfig(analysisFolder)
      
#### Make result folders ###      
for resFolder in (foldersToCreate):
  os.makedirs(analysisFolder +  '/' + resFolder, exist_ok = True)

#### Run individual systematic Error Estimation
if config['doInclusive'] == 1:
  print('Do inclusive analysis')
  jetString = 'Jets_Inclusive'
  runSystematicProcess(config['systematicsInclusive'], systematics, config, jetString, systematicDay) 
  runCalculateEfficiency(jetString, config, systematicDay, systematicDay)

if config['doJets'] == 1:
  print('Do jet analysis')
  jetString = 'Jets'
  runSystematicProcess(config['systematicsJets'], systematics, config, jetString, systematicDay)
  runCalculateEfficiency(jetString, config, systematicDay, systematicDay)
  
if config['doUE'] == 1:
  print('Do UE analysis')
  jetString = 'Jets_UE'
  runSystematicProcess(config['systematicsUE'], systematics, config, jetString, systematicDay)
  runCalculateEfficiency(jetString, config, systematicDay, systematicDay)
  
#### Subtract Underlying event ###
arguments = {
  'jetFilePattern' : analysisFolder + "/Efficiencycorrected/output_EfficiencyCorrection_outputSystematicsTotal_SummedSystematicErrors_Jets_%s__%s%s__" + systematicDay + "__" + systematicDay + ".root",
  'ueFilePattern' : analysisFolder + "/Efficiencycorrected/output_EfficiencyCorrection_outputSystematicsTotal_SummedSystematicErrors_Jets_UE_%s__%s%s__" + systematicDay + "__" + systematicDay + ".root",
  'jetPtStepsString' : config['jetPtString'],
  'centStepsString' : config['centString'],
  'modesInputString' : config['modesJetsString'],
  'outputFilePattern' : analysisFolder + "/UEsubtractedJetResults/output_EfficiencyCorrection_outputSystematicsTotal_SummedSystematicErrors_Jets_%s__%s%s__" + systematicDay + "__" + systematicDay + "_UEsubtractedJetResults.root",
}

callRootMacro("SubtractUnderlyingEvent", arguments)
