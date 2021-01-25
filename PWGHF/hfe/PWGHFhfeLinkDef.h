#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class  AliHFEtools+;
#pragma link C++ class  AliHFEparamBag+;
#pragma link C++ class  AliHFEcollection+;
#pragma link C++ class  AliHFEvarManager+;
#pragma link C++ class  AliHFEvarManager::AliHFEvariable+;
#pragma link C++ class  AliHFEcontainer+;
#pragma link C++ class  AliHFEcontainer::AliHFEvarInfo+;
#pragma link C++ class  AliHFEmcQA+;
#pragma link C++ class  AliHFEpairs+;
#pragma link C++ class  AliHFEsecVtxs+;
#pragma link C++ class  AliHFEsecVtx+;
#pragma link C++ class  AliHFEpriVtx+;
#pragma link C++ class  AliHFEelecbackground+;
#pragma link C++ class  AliHFEspectrum+;
#pragma link C++ class  AliHFEtaggedTrackAnalysis+;

#pragma link C++ class  AliHFEV0info+;
#pragma link C++ class  AliHFEV0pid+;
#pragma link C++ class  AliHFEV0cuts+;
#pragma link C++ class  AliHFEV0pidMC+;
#pragma link C++ class  AliHFEpidQA+;
#pragma link C++ class  AliHFEV0taginfo+;
#pragma link C++ class  AliHFEV0taginfo::AliHFEV0tag+;
#pragma link C++ class  AliHFEitsPIDqa+;
#pragma link C++ class  AliHFEtrdPIDqa+;
#pragma link C++ class  AliHFEdetPIDqa+;
#pragma link C++ class  AliHFEtpcPIDqa+;
#pragma link C++ class  AliHFEtofPIDqa+;
#pragma link C++ class  AliHFEemcalPIDqa+;
#pragma link C++ class  AliHFEpidEMCAL+;
#pragma link C++ class  AliHFEtrdPIDqaV1+;
#pragma link C++ class  AliHFEpidQAmanager+;
#pragma link C++ class  AliHFEpid+;
#pragma link C++ class  AliHFEpidBase+;
#pragma link C++ class  AliHFEpidITS+;
#pragma link C++ class  AliHFEpidTPC+;
#pragma link C++ class  AliHFEpidTRD+;
#pragma link C++ class  AliHFEpidTOF+;
#pragma link C++ class  AliHFEpidMC+;

#pragma link C++ class  AliHFEcuts+;
#pragma link C++ class  AliHFEcutStep+;
#pragma link C++ class  AliHFEtrackFilter+;
#pragma link C++ class  AliHFEextraCuts+;
#pragma link C++ class  AliHFEextraEventCuts+;
#pragma link C++ class  AliHFEsignalCuts+;

#pragma link C++ class  AliHFEpostAnalysis+;
#pragma link C++ class  AliAnalysisTaskHFE+;
#pragma link C++ class  AliAnalysisTaskHFEpidQA+;
#pragma link C++ class  AliHFEefficiency+;

#pragma link C++ class  AliHFEOADBThresholdsTRD+;
#pragma link C++ class  AliHFEbayesPIDqa+;
#pragma link C++ class  AliHFEpidBayes+;

#pragma link C++ class  AliHFEpidObject+;
#pragma link C++ class  AliAnalysisTaskEHCorrel+;
#pragma link C++ class  AliAnalysisTaskFlowTPCEMCalEP+;
#pragma link C++ class  AliAnalysisTaskFlowTPCEMCalRun2+;
#pragma link C++ class  AliAnalysisTaskHFECal+;
#pragma link C++ class  AliAnalysisTaskEMCalHFEpA+;

#pragma link C++ class  AliAnalysisTaskHFEpACorrelation+;
#pragma link C++ class  AliEHCParticle+;
#pragma link C++ class  AliHFEHCParticle+;

#pragma link C++ class  AliHFEdebugTreeTask+;

#pragma link C++ class  AliHFEVZEROEventPlane+;
#pragma link C++ class  AliAnalysisTaskFlowTPCTOFEPSP+;

#pragma link C++ class  AliSelectNonHFE+;
#pragma link C++ class  AliHFENonPhotonicElectron+;
#pragma link C++ class  AliHFEdebugTreeTaskAOD+;

#pragma link C++ class  AliHFECorrectSpectrumBase+;
#pragma link C++ class  AliHFEInclusiveSpectrum+;
#pragma link C++ class  AliHFEInclusiveSpectrumQA+;
#pragma link C++ class  AliHFEBeautySpectrum+;
#pragma link C++ class  AliHFEBeautySpectrumQA+;
#pragma link C++ class  AliHFEsmearDCA+;
#pragma link C++ class  AliAnalysisTaskFlowTPCEMCalQCSP+;

#pragma link C++ class  AliHFEreducedEventCreatorAOD+;
#pragma link C++ class  AliHFEreducedEventCreatorESD+;
#pragma link C++ class  AliHFEreducedEvent+;
#pragma link C++ class  AliHFEreducedTrack+;
#pragma link C++ class  AliHFEreducedMCParticle+;
#pragma link C++ class  AliHFEminiTrack+;
#pragma link C++ class  AliHFEminiEvent+;
#pragma link C++ class  AliHFEminiEventCreator+;
#pragma link C++ class  AliAnalysisTaskHFEQA+;
#pragma link C++ class  AliAnalysisTaskHFEemcQA+;
#pragma link C++ class  AliAnalysisTaskBeautyCal+;
#pragma link C++ class  AliAnalysisTaskFlowITSTPCTOFQCSP+;
#pragma link C++ class  AliAnalysisTaskHFEIPDistribution+;
#pragma link C++ class  AliAnalysisTaskHFEBeautyMCTemplates+;
#pragma link C++ class  AliAnalysisTaskHFEBeautyMCTemplatesRun2+;
#pragma link C++ class  AliAnalysisTaskHFEIPCorrection+;
#pragma link C++ class  AliAnalysisTaskHFEMulti+;
#pragma link C++ class  AliAnalysisTaskHFEMultiplicity+;
#pragma link C++ class  AliAnalysisTaskHFEtemplate+;
#pragma link C++ class  AliHFEhadronicbackground+;
#pragma link C++ class  AliHFEUnfolding+;
#pragma link C++ class  AliAnalysisTaskHFEEfficiency+;
#pragma link C++ class  AliAnalysisTaskHaHFECorrel+;
#pragma link C++ class  AliAnalysisTaskCaloHFEpPbRun2+;
#pragma link C++ class  AliAnalysisTaskCaloHFEpp+;

#pragma link C++ class  AliBasicParticleHaHFE+;
#pragma link C++ class  AliehDPhiBasicParticle+;
#pragma link C++ class  AliehDPhiBasicParticlepPbRun2+;

#pragma link C++ class  AliAnalysisTaskTPCCalBeauty+;
#pragma link C++ class  AliAnalysisHFETPCTOFNew+;
#pragma link C++ class  AliAnalysisHFETPCTOFBeauty+;
#pragma link C++ class  AliAnalysisHFEppTPCTOFBeauty+;
#pragma link C++ class  AliAnalysisHFEppTPCTOFBeauty5TeVNew+;

#pragma link C++ class  AliAnalysisTaskHFEmultTPCTOF+;
#pragma link C++ class  AliAnalysisTaskHFEBESpectraEMC+;

#pragma link C++ class  AliAnalysisTaskHFETPCTOFMultiplicity+;
#pragma link C++ class  AliAnalysisTaskQAHFE+;
#pragma link C++ class  AliAnalysisTask_QA_EMCALElectrons+;
#pragma link C++ class  AliAnalysisTaskHFEmultTPCEMCAL+;
#pragma link C++ class  AliAnalysisTaskHFEBeautyMultiplicity+;
#pragma link C++ class  AliAnalysisHFEppEMCalBeauty+;
#pragma link C++ class  AliAnalysisTaskIPResolBeautyppCal+;
#endif
