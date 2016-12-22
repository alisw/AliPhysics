// -*- mode: c++ -*-
/* Copyright (C) 2007 Christian Holm Christensen <cholm@nbi.dk>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */
#ifdef __CINT__
/**
 * @file   PWGLFforward2LinkDef.h
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri May 24 09:24:36 2013
 *
 * @brief  Link specifications
 */
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclasses;
// ZDC tasks
#pragma link C++ class AliAnalysisTaskZDCTreeMaker+;
#pragma link C++ class AliAnalysisTaskZDCpAcalib+;

// PMD tasks
#pragma link C++ class AliAnalysisTaskPMD+;
#pragma link C++ class AliAnalysisTaskPMDSim+;

// AOD objects
#pragma link C++ class AliAODCentralMult+;
#pragma link C++ class AliAODForwardEP+;
#pragma link C++ class AliAODForwardMult+;
#pragma link C++ class AliAODMultEventClass+;


// Interface to OADB backed by a TTree
#pragma link C++ class AliOADBForward+;
#pragma link C++ class AliOADBForward::Entry+;
#pragma link C++ class AliOADBForward::Table+;

// Correction managers
#pragma link C++ class AliCorrectionManagerBase+;
#pragma link C++ class AliCorrectionManagerBase::Correction+;
// Note: custom streamer to ensure singleton consistency!
#pragma link C++ class AliForwardCorrectionManager-;
#pragma link C++ class AliCentralCorrectionManager-;

// Base tasks etc.
#pragma link C++ class AliBaseESDTask+;
#pragma link C++ class AliBaseAODTask+;
#pragma link C++ class AliBaseMCCorrectionsTask+;
#if ROOT_VERSION_CODE < 0x56300/* ROOT_VERSION(5,99,0)*/ && __GNUC__ < 5
// Used internally and never streamed.
#pragma link C++ class AliBaseMCCorrectionsTask::VtxBin+;
#endif
#pragma link C++ class AliBasedNdetaTask+;
#if ROOT_VERSION_CODE < 0x56300/* ROOT_VERSION(5,99,0)*/ && __GNUC__ < 5
// Used internally and never streamed.
#pragma link C++ class AliBasedNdetaTask::CentralityBin+;
// Used internally and never streamed.
#pragma link C++ class AliBasedNdetaTask::Sum+;
#endif
#pragma link C++ class AliBaseMCTrackDensity+;

// Central (SPD) code
#pragma link C++ class AliCentralCorrAcceptance+;
#pragma link C++ class AliCentralCorrSecondaryMap+;
#pragma link C++ class AliCentraldNdetaTask+;
#pragma link C++ class AliCentralMCCorrectionsTask+;
#if ROOT_VERSION_CODE < 0x56300/* ROOT_VERSION(5,99,0)*/ && __GNUC__ < 5
// Used internally and never streamed.
#pragma link C++ class AliCentralMCCorrectionsTask::VtxBin+;
#endif
#pragma link C++ class AliCentralMCMultiplicityTask+;
#pragma link C++ class AliCentralMultiplicityTask+;
#if ROOT_VERSION_CODE < 0x56300/* ROOT_VERSION(5,99,0)*/ && __GNUC__ < 5
// Used internally and never streamed
#pragma link C++ class AliCentralMultiplicityTask::VtxBin+;
#endif
#pragma link C++ class AliSPDMCTrackDensity+;

// Aux tasks and code
#pragma link C++ class AliCopyHeaderTask+;
#pragma link C++ class AliMCTruthdNdetaTask+;
#if ROOT_VERSION_CODE < 0x56300/* ROOT_VERSION(5,99,0)*/ && __GNUC__ < 5
// Used internally and never streamed.
#pragma link C++ class AliMCTruthdNdetaTask::CentralityBin+;
#endif
#pragma link C++ class AliDisplacedVertexSelection+;
#pragma link C++ class AliDisplacedVertexSelectionAD+;
#pragma link C++ class AliPoissonCalculator+;
#pragma link C++ class AliMCAuxHandler+;
#pragma link C++ class AliTestAD+;


// Forward AUX (Cuts, etc.)
#pragma link C++ class AliFMDMultCuts+;
#pragma link C++ class AliForwardUtil+;
#pragma link C++ class AliForwardUtil::Histos+;
#pragma link C++ class AliForwardUtil::RingHistos+;

// FMD corrections
#pragma link C++ class AliFMDCorrAcceptance+;
#pragma link C++ class AliFMDCorrDoubleHit+;
#pragma link C++ class AliFMDCorrELossFit+;
#pragma link C++ class AliFMDCorrELossFit::ELossFit+;
#pragma link C++ class AliFMDCorrMergingEfficiency+;
#pragma link C++ class AliFMDCorrSecondaryMap+;
#pragma link C++ class AliFMDCorrVertexBias+;
#pragma link C++ class AliFMDCorrNoiseGain+;

// FMD algorithms
#pragma link C++ class AliFMDEnergyFitter+;
#pragma link C++ class AliFMDEnergyFitter::RingHistos+;
#pragma link C++ class AliFMDEventInspector+;
#pragma link C++ class AliFMDEventPlaneFinder+;
#pragma link C++ class AliFMDHistCollector+;
#pragma link C++ class AliFMDESDFixer+;
#pragma link C++ class AliFMDSharingFilter+;
#if ROOT_VERSION_CODE < 0x56300/* ROOT_VERSION(5,99,0)*/ && __GNUC__ < 5
// Used internally and never streamed
#pragma link C++ class AliFMDSharingFilter::RingHistos+;
#endif
#pragma link C++ class AliFMDDensityCalculator+;
#if ROOT_VERSION_CODE < 0x56300/* ROOT_VERSION(5,99,0)*/ && __GNUC__ < 5
// Used internally and never streamed
#pragma link C++ class AliFMDDensityCalculator::RingHistos+;
#endif
#pragma link C++ class AliFMDCorrector+;
#if ROOT_VERSION_CODE < 0x56300/* ROOT_VERSION(5,99,0)*/ && __GNUC__ < 5
// Used internally and never streamed
#pragma link C++ class AliFMDCorrector::RingHistos+;
#endif
#pragma link C++ class AliMultEventClassifier+;

// FMD MC algorithms
#pragma link C++ class AliFMDMCCorrector+;
#pragma link C++ class AliFMDMCDensityCalculator+;
#pragma link C++ class AliFMDMCEventInspector+;
#pragma link C++ class AliFMDMCSharingFilter+;
#pragma link C++ class AliFMDMCTrackDensity+;

// FMD MC investigations
#pragma link C++ class AliFMDMCTrackELoss+;
#pragma link C++ class AliFMDMCTrackELoss::Hit+;
#pragma link C++ class AliFMDMCTrackInspector+;
#pragma link C++ class AliFMDMCTrackInspector::RingHistos+;
#pragma link C++ class AliFMDMCTrackInspectorTask+;

// Forward (FMD) tasks
#pragma link C++ class AliFMDEnergyFitterTask+;
#pragma link C++ class AliFMDEventPlaneTask+;
#pragma link C++ class AliForwarddNdetaTask+;
#if ROOT_VERSION_CODE < 0x56300/* ROOT_VERSION(5,99,0)*/ && __GNUC__ < 5
// Used internally and never streamed.
#pragma link C++ class AliForwarddNdetaTask::CentralityBin+;
#endif
#pragma link C++ class AliForwardFlowTaskQC+;
#if ROOT_VERSION_CODE < 0x56300/* ROOT_VERSION(5,99,0)*/ && __GNUC__ < 5
// Used internally and never streamed.
#pragma link C++ class AliForwardFlowTaskQC::CumuHistos+;
#pragma link C++ class AliForwardFlowTaskQC::VertexBin+;
#endif
#pragma link C++ class AliForwardMCCorrectionsTask+;
#if ROOT_VERSION_CODE < 0x56300/* ROOT_VERSION(5,99,0)*/ && __GNUC__ < 5
// Used internally and never streamed.
#pragma link C++ class AliForwardMCCorrectionsTask::VtxBin+;
#endif
#pragma link C++ class AliForwardMCFlowTaskQC+;
#pragma link C++ class AliForwardMCMultiplicityTask+;
#pragma link C++ class AliForwardMultiplicityBase+;
#pragma link C++ class AliForwardMultiplicityTask+;
#pragma link C++ class AliForwardQATask+;
#pragma link C++ class AliBaseMultTask+;
#pragma link C++ class AliBaseMultTask::Bin+;
#pragma link C++ class AliForwardTriggerBiasCorrection+;
#pragma link C++ class AliForwardTriggerBiasCorrection::Bin+;
#pragma link C++ class AliForwardCreateResponseMatrices+;
#pragma link C++ class AliForwardCreateResponseMatrices::Bin+;
#pragma link C++ class AliForwardMultiplicityDistribution+;
#pragma link C++ class AliForwardMultiplicityDistribution::Bin+;
#pragma link C++ class AliForwardMultDists+;
#pragma link C++ class AliForwardMultDists::EtaBin+;
#pragma link C++ class AliForwardMultDists::BinSpec+;

//  MC Weights
#pragma link C++ class AliForwardFlowWeights+;
#pragma link C++ class AliBaseMCWeights+;
#pragma link C++ class AliSimplePidWeights+;
#pragma link C++ class AliPtEtaPidWeights+;

#pragma link C++ class AliMultEventClassifierTask;
#endif
//
// EOF
//
