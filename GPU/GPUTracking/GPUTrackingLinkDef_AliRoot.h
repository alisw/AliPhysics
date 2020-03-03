//**************************************************************************\
//* This file is property of and copyright by the ALICE Project            *\
//* ALICE Experiment at CERN, All rights reserved.                         *\
//*                                                                        *\
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *\
//*                  for The ALICE HLT Project.                            *\
//*                                                                        *\
//* Permission to use, copy, modify and distribute this software and its   *\
//* documentation strictly for non-commercial purposes is hereby granted   *\
//* without fee, provided that the above copyright notice appears in all   *\
//* copies and that both the copyright notice and this permission notice   *\
//* appear in the supporting documentation. The authors make no claims     *\
//* about the suitability of this software for any purpose. It is          *\
//* provided "as is" without express or implied warranty.                  *\
//**************************************************************************

/// \file GPUTrackingLinkDef_AliRoot.h
/// \author David Rohr

#if defined(__CINT__) || defined(__CLING__)

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliGPU::gpu::GPUTPCTrack + ;
#pragma link C++ class AliGPU::gpu::GPUTPCTracklet + ;
#pragma link C++ class AliGPU::gpu::GPUTPCBaseTrackParam + ;
#pragma link C++ class AliGPU::gpu::GPUTPCTrackParam + ;
#pragma link C++ class AliGPU::gpu::GPUTPCRow + ;
#pragma link C++ class AliGPU::gpu::GPUTPCGrid + ;
#pragma link C++ class GPUTPCTrackerComponent + ;
#pragma link C++ class AliGPU::gpu::GPUTPCNeighboursFinder + ;
#pragma link C++ class AliGPU::gpu::GPUTPCNeighboursCleaner + ;
#pragma link C++ class AliGPU::gpu::GPUTPCStartHitsFinder + ;
#pragma link C++ class AliGPU::gpu::GPUTPCTrackletConstructor + ;
#pragma link C++ class AliGPU::gpu::GPUTPCTrackletSelector + ;
#pragma link C++ class GPUTPCGlobalMergerComponent + ;
#pragma link C++ class AliGPU::gpu::GPUTPCSliceOutput + ;
#pragma link C++ class AliGPU::gpu::GPUTPCGMTrackParam + ;
#pragma link C++ class AliGPU::gpu::GPUTPCGMSliceTrack + ;
#pragma link C++ class AliGPU::gpu::GPUTPCGMPolynomialField + ;
#pragma link C++ class AliGPU::gpu::GPUTPCGMPropagator + ;
#pragma link C++ class AliGPU::gpu::GPUTPCGMPhysicalTrackModel + ;
#pragma link C++ class GPUTPCGMPolynomialFieldManager + ;
#pragma link C++ class AliHLTTPCClusterStatComponent + ;

//#pragma link C++ class AliGPU::gpu::GPUTRDTrack+; //Templated, should add linkdef for specialization, but with an ifdef for ROOT >= 6 only
//#pragma link C++ class AliGPU::gpu::GPUTRDTracker+;
#pragma link C++ class GPUTRDTrackerComponent + ;
//#pragma link C++ class AliGPU::gpu::GPUTRDTrackletWord+;
#pragma link C++ class GPUTRDTrackletReaderComponent + ;

#endif
