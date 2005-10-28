#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ global gAlice;
#pragma link C++ global gMC;
 
// AliDCSClient classes ...
#pragma link C++ class AliSimpleValue;
#pragma link C++ class AliSimpleValue::BoolHolder;
#pragma link C++ class AliSimpleValue::ByteHolder;
#pragma link C++ class AliSimpleValue::IntHolder;
#pragma link C++ class AliSimpleValue::UIntHolder;
#pragma link C++ class AliSimpleValue::FloatHolder;
#pragma link C++ class AliSimpleValue::DynHolder;
#pragma link C++ class AliSimpleValue::DynBoolHolder;
#pragma link C++ class AliSimpleValue::DynByteHolder;
#pragma link C++ class AliSimpleValue::DynIntHolder;
#pragma link C++ class AliSimpleValue::DynUIntHolder;
#pragma link C++ class AliSimpleValue::DynFloatHolder;
#pragma link C++ class AliDCSValue;

#pragma link C++ class AliDCSMessage;
#pragma link C++ class AliDCSClient;

// Shuttle classes ...
#pragma link C++ class  AliShuttleConfig;
#pragma link C++ class  AliShuttleConfig::ConfigHolder;
#pragma link C++ class  AliShuttle;
#pragma link C++ class  AliCDBPreProcessor;

#endif
