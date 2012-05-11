//-*- Mode: C++ -*-
// $Id$

/// @file   LoadLibraries.C
/// @author Matthias.Richter@ift.uib.no
/// @date   2012-05-11
/// @brief  Load library dependencies for PWGHF/correlationHF
///

void LoadLibraries()
{
  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWGHF/correlationHF");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWGflowBase.so");
  gSystem->Load("libPWGflowTasks.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGHFvertexingHF");
}
