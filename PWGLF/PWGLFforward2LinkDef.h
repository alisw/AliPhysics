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
#pragma link C++ class AliAnalysisTaskZDCPbPb+;
#pragma link C++ class AliAnalysisTaskZDCTreeMaker+;
#pragma link C++ class AliAnalysisTaskZDCpAcalib+;

// PMD tasks
#pragma link C++ class AliAnalysisTaskPMD+;
#pragma link C++ class AliAnalysisTaskPMDSim+;

// AOD objects
#pragma link C++ class AliAODCentralMult+;
#pragma link C++ class AliAODForwardEP+;
#pragma link C++ class AliAODForwardMult+; 


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

// Central (SPD) code 
#pragma link C++ class AliCentralCorrAcceptance+;
#pragma link C++ class AliCentralCorrSecondaryMap+;
#pragma link C++ class AliCentraldNdetaTask+;
#pragma link C++ class AliCentralMCCorrectionsTask+;
#pragma link C++ class AliCentralMCCorrectionsTask::VtxBin+;
#pragma link C++ class AliCentralMCMultiplicityTask+;
#pragma link C++ class AliCentralMultiplicityTask+;
#pragma link C++ class AliCentralMultiplicityTask::VtxBin+;
#pragma link C++ class AliSPDMCTrackDensity+;

// Aux tasks and code 
#pragma link C++ class AliCopyHeaderTask+;
#pragma link C++ class AliBasedNdetaTask+;
#pragma link C++ class AliBasedNdetaTask::CentralityBin+;
#pragma link C++ class AliBasedNdetaTask::Sum+;
#pragma link C++ class AliBaseMCTrackDensity+;
#pragma link C++ class AliMCTruthdNdetaTask+;
#pragma link C++ class AliMCTruthdNdetaTask::CentralityBin+;
#pragma link C++ class AliDisplacedVertexSelection+;
#pragma link C++ class AliPoissonCalculator+;
#pragma link C++ class AliMCAuxHandler+;

// Forward AUX (Cuts, etc.)
#pragma link C++ class AliFMDMultCuts+;
#pragma link C++ class AliForwardFlowWeights+;
#pragma link C++ class AliForwardUtil+;
#pragma link C++ class AliForwardUtil::Histos+;
#pragma link C++ class AliForwardUtil::RingHistos+;

// FMD corrections
#pragma link C++ class AliFMDCorrAcceptance+;
#pragma link C++ class AliFMDCorrDoubleHit+;
#pragma link C++ class AliFMDCorrector+;
#pragma link C++ class AliFMDCorrector::RingHistos+;
#pragma link C++ class AliFMDCorrELossFit+;
#pragma link C++ class AliFMDCorrELossFit::ELossFit+;
#pragma link C++ class AliFMDCorrMergingEfficiency+;
#pragma link C++ class AliFMDCorrSecondaryMap+;
#pragma link C++ class AliFMDCorrVertexBias+;

// FMD algorithms 
#pragma link C++ class AliFMDDensityCalculator+;
#pragma link C++ class AliFMDDensityCalculator::RingHistos+;
#pragma link C++ class AliFMDEnergyFitter+;
#pragma link C++ class AliFMDEnergyFitter::RingHistos+;
#pragma link C++ class AliFMDEventInspector+;
#pragma link C++ class AliFMDEventPlaneFinder+;
#pragma link C++ class AliFMDHistCollector+;
#pragma link C++ class AliFMDSharingFilter+;
#pragma link C++ class AliFMDSharingFilter::RingHistos+;

// FMD MC algorithms
#pragma link C++ class AliFMDMCCorrector+;
#pragma link C++ class AliFMDMCDensityCalculator+;
#pragma link C++ class AliFMDMCEventInspector+;
#pragma link C++ class AliFMDMCSharingFilter+;
#pragma link C++ class AliFMDMCTrackDensity+;

// Forward (FMD) tasks 
#pragma link C++ class AliBaseESDTask+;
#pragma link C++ class AliBaseAODTask+;
#pragma link C++ class AliBaseMCCorrectionsTask+;
#pragma link C++ class AliBaseMCCorrectionsTask::VtxBin+;
#pragma link C++ class AliFMDEnergyFitterTask+;
#pragma link C++ class AliFMDEventPlaneTask+;
#pragma link C++ class AliForwarddNdetaTask+;
#pragma link C++ class AliForwarddNdetaTask::CentralityBin+;
#pragma link C++ class AliForwardFlowTaskQC+;
#pragma link C++ class AliForwardFlowTaskQC::VertexBin+;
#pragma link C++ class AliForwardMCCorrectionsTask+;
#pragma link C++ class AliForwardMCCorrectionsTask::VtxBin+;
#pragma link C++ class AliForwardMCFlowTaskQC+;
#pragma link C++ class AliForwardMCMultiplicityTask+;
#pragma link C++ class AliForwardMultiplicityBase+;
#pragma link C++ class AliForwardMultiplicityTask+;
#pragma link C++ class AliForwardQATask+;
#pragma link C++ class AliForwardCreateResponseMatrices+;
#pragma link C++ class AliForwardCreateResponseMatrices::Bin+;
#pragma link C++ class AliForwardMultiplicityDistribution+;
#pragma link C++ class AliForwardMultiplicityDistribution::Bin+;
#pragma link C++ class AliForwardMultDists+;
#pragma link C++ class AliForwardMultDists::EtaBin+;
#pragma link C++ class AliForwardMultDists::BinSpec+;

#else
# error Not for compilation 
#endif
//
// EOF
//
