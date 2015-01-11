/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \file MUONrawLinkDef.h
/// \brief The CINT link definitions for \ref raw 

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// raw data
#pragma link C++ class AliMUONDDLTrigger+;
#pragma link C++ class AliMUONDarcHeader+;
#pragma link C++ class AliMUONRegHeader+;
#pragma link C++ class AliMUONLocalStruct+;
#pragma link C++ class AliMUONPayloadTrigger+;
#pragma link C++ class AliMUONDDLTracker+;
#pragma link C++ class AliMUONDspHeader+;
#pragma link C++ class AliMUONBlockHeader+;
#pragma link C++ class AliMUONBusStruct+;
#pragma link C++ class AliMUONPayloadTracker+;
#pragma link C++ class AliMUONVRawStreamTracker+;
#pragma link C++ class AliMUONRawStreamTracker+;
#pragma link C++ class AliMUONRawStreamTrackerHP+;
#pragma link C++ class AliMUONRawStreamTrackerHP::AliBlockHeader+;
#pragma link C++ class AliMUONRawStreamTrackerHP::AliDspHeader+;
#pragma link C++ class AliMUONRawStreamTrackerHP::AliBusPatch+;
#pragma link C++ class AliMUONVRawStreamTrigger+;
#pragma link C++ class AliMUONRawStreamTrigger+;
#pragma link C++ class AliMUONRawStreamTriggerHP+;
#pragma link C++ class AliMUONRawStreamTriggerHP::AliHeader+;
#pragma link C++ class AliMUONRawStreamTriggerHP::AliRegionalHeader+;
#pragma link C++ class AliMUONRawStreamTriggerHP::AliLocalStruct+;
#pragma link C++ class AliMUONRawStream+;

#endif


