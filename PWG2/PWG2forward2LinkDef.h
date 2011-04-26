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
/* $Id: FMDflowLinkDef.h 23165 2007-12-19 01:36:20Z cholm $ */
/** @file    FMDbaseLinkDef.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 14:18:46 2006
    @brief   Link specifications for base library 
*/
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclasses;

#pragma link C++ class AliForwardUtil+;
#pragma link C++ class AliForwardUtil::Histos+;
#pragma link C++ class AliForwardUtil::RingHistos+;
#pragma link C++ class AliFMDEventInspector+;
#pragma link C++ class AliFMDMCEventInspector+;
#pragma link C++ class AliFMDSharingFilter+;
#pragma link C++ class AliFMDSharingFilter::RingHistos+;
#pragma link C++ class AliFMDMCSharingFilter+;
//#pragma link C++ class AliFMDMCSharingFilter::RingHistos+;
#pragma link C++ class AliFMDMCTrackDensity+;
#pragma link C++ class AliFMDEnergyFitter+;
#pragma link C++ class AliFMDEnergyFitter::RingHistos+;
#pragma link C++ class AliFMDEnergyFitterTask+;
#pragma link C++ class AliFMDDensityCalculator+;
#pragma link C++ class AliFMDDensityCalculator::RingHistos+;
#pragma link C++ class AliFMDMCDensityCalculator+;
#pragma link C++ class AliFMDCorrector+;
#pragma link C++ class AliFMDCorrector::RingHistos+;
#pragma link C++ class AliFMDMCCorrector+;
//#pragma link C++ class AliFMDMCCorrections::RingHistos+;
#pragma link C++ class AliFMDHistCollector+;
#pragma link C++ class AliFMDCorrAcceptance+;
#pragma link C++ class AliFMDCorrELossFit+;
#pragma link C++ class AliFMDCorrELossFit::ELossFit+;
#pragma link C++ class AliFMDCorrSecondaryMap+;
#pragma link C++ class AliFMDCorrDoubleHit+;
#pragma link C++ class AliFMDCorrVertexBias+;
#pragma link C++ class AliFMDCorrMergingEfficiency+;
#pragma link C++ class AliAODForwardMult+;
#pragma link C++ class AliForwardMultiplicityBase+;
#pragma link C++ class AliForwardMultiplicityTask+;
#pragma link C++ class AliForwardMCMultiplicityTask+;
// Note: custom streamer to ensure singleton consistency!
#pragma link C++ class AliForwardCorrectionManager-;
#pragma link C++ class AliForwardMCCorrectionsTask+;
#pragma link C++ class AliForwardMCCorrectionsTask::VtxBin+;
#pragma link C++ class AliForwarddNdetaTask+;
#pragma link C++ class AliForwarddNdetaTask::CentralityBin+;
#pragma link C++ class AliBasedNdetaTask+;
#pragma link C++ class AliBasedNdetaTask::CentralityBin+;

#pragma link C++ class AliCentralMultiplicityTask+;
#pragma link C++ class AliCentralMultiplicityTask::Manager+;
#pragma link C++ class AliCentralMCMultiplicityTask+;
#pragma link C++ class AliCentralMCCorrectionsTask+;
#pragma link C++ class AliCentralMCCorrectionsTask::VtxBin+;
#pragma link C++ class AliAODCentralMult+;
#pragma link C++ class AliCentralCorrSecondaryMap+;
#pragma link C++ class AliCentralCorrAcceptance+;
#pragma link C++ class AliCentraldNdetaTask+;
#pragma link C++ class AliForwardFlowUtil+;
#pragma link C++ class AliForwardFlowTaskQC+;
#pragma link C++ class AliSPDMCTrackDensity+;

#else
# error Not for compilation 
#endif
//
// EOF
//
