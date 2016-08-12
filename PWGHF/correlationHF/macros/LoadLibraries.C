//-*- Mode: C++ -*-
// $Id$

/// @file   LoadLibraries.C
/// @author Matthias.Richter@ift.uib.no
/// @date   2012-05-11
/// @brief  Load library dependencies for PWGHF/correlationHF
///

void LoadLibraries()
{
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/CDB -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWGHF/hfe -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWGHF/correlationHF -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -g");
 // gSystem->AddIncludePath("-I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWGHF/correlationHF");
/*  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGHFvertexingHF");
  */
	//Gui,Proof,Minuit,XMLParser,RAWDatabase,RAWDatarec,CDB,STEER,TOFbase,TOFrec,
	//CORRFW,PWGflowBase,PWGflowTasks,PWGHFbase,PWGHFvertexingHF


  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libPhysics");
  gSystem->Load("libVMC");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGHFbase");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
  gSystem->Load("libPWGHFvertexingHF");
  gSystem->Load("libPWGHFcorrelationHF");
	
	return;
}
