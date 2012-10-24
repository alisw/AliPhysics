//-*- Mode: C++ -*-
// $Id$

/// @file   LoadLibraries.C
/// @author Matthias.Richter@ift.uib.no
/// @date   2012-05-11
/// @brief  Load library dependencies for PWGHF/correlationHF
///

void LoadLibraries()
{
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/CDB -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGHF -I$ALICE_ROOT/PWGHF/base -I$ALICE_ROOT/PWGHF/hfe -I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWGHF/correlationHF -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_ROOT/PWG/FLOW/Tasks -g");
 // gSystem->AddIncludePath("-I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWGHF/correlationHF");
/*  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWGflowBase.so");
  gSystem->Load("libPWGflowTasks.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGHFvertexingHF.so");
  */
	//Gui,Proof,Minuit,XMLParser,RAWDatabase,RAWDatarec,CDB,STEER,TOFbase,TOFrec,
	//CORRFW,PWGflowBase,PWGflowTasks,PWGHFbase,PWGHFvertexingHF
	
	gSystem->Load("libTree.so");
	//gSystem->Load("libGeom.so");
	gSystem->Load("libPhysics.so");
//	gSystem->Load("libVMC.so");
//	gSystem->Load("libMinuit.so");
	gSystem->Load("libSTEERBase.so");
	gSystem->Load("libESD.so");
	gSystem->Load("libAOD.so"); 
	gSystem->Load("libANALYSIS.so");
	gSystem->Load("libOADB.so");
	gSystem->Load("libANALYSISalice.so");
	gSystem->Load("libCORRFW.so");
	gSystem->Load("libPWGHFbase.so");
	gSystem->Load("libPWGflowBase.so");
	gSystem->Load("libPWGflowTasks.so");
	gSystem->Load("libPWGHFvertexingHF.so");
	
	gSystem->Load("libPWGHFcorrelationHF.so");
	
	return;
}
