//-*- Mode: C++ -*-
// $Id$

/// @file   setupDxHFE.C
/// @author Matthias.Richter@cern.ch
/// @date   2012-09-12
/// @brief  Setup the environment specifically for running DxHFE analysis
///
/// Helper macro to initialize the environment for DxHFE. The macro can just
/// prepend other macros like run-single-task in the command line.
/// Usage:
/// aliroot -b -q -l setupDxHFE.C'("localAodDir")' run-single-task.C'(...)'
///
/// Example:
/// aliroot -b -q -l setupDxHFE.C run-single-task.C'(...)'
///
/// The macro has the following tasks:
/// - load the necessary libraries, in order to have those library names also
///   available for the alien handler initialization, a specific configuration
///   object is created
/// - setting a default analysis name via a configuration object
/// - the optional parameter 'localAodDir' allows to create an input chain from
///   local AODs; either a single AliAOD.root, or a folder containing directories
///   named "1, 2, ..."

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// environment specific for DxHFE
//
const char* includePath="-I$ALICE_ROOT/PWGHF/correlationHF -I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWGHF/hfe";
const char* libraryDependencies=
  "libSTEERBase "
  "libESD "
  "libAOD "
  "libANALYSIS "
  "libANALYSISalice "
  "libPWGflowBase.so "
  "libPWGflowTasks.so "
  "libCORRFW.so "
  "libPWGHFvertexingHF.so "
  "libPWGHFhfe.so "
  "libPWGHFcorrelationHF.so "
  ;

void setupDxHFE(const char* aodDirectory=NULL)
{
  gSystem->AddIncludePath(includePath);
  TString libraries=libraryDependencies;
  TObjArray* pTokens=libraries.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntriesFast(); i++) {
      if (gSystem->Load(pTokens->At(i)->GetName())==0) {
	cout << "loading " << pTokens->At(i)->GetName() << endl;
      }
    }
    delete pTokens;
  }
  libraries="";

  // allow run-single-task to fetch the analysis name and library names
  if (gDirectory) gDirectory->Add(new TNamed("analysis_name", "DxHFECorrelation"));
  if (gDirectory) gDirectory->Add(new TNamed("analysis_libraries", libraryDependencies));

  if (aodDirectory) {
    // create AOD tree from local files
    // the created object is added automatically to gDirectory and can be fetched
    // from there later
    gROOT->LoadMacro("$ALICE_ROOT/PWGHF/vertexingHF/MakeAODInputChain.C");
    TString aodPathName(aodDirectory);
    if (!aodPathName.EndsWith("/")) aodPathName+="/";
    aodPathName+="AliAOD.root";
    if (gSystem->AccessPathName(aodPathName)==0) {
      // Create a chain with one set of AliAOD.root and AliAOD.VertexingHF.root. The set needs 
      // to be located in the same folder as you run from (physically or linked)
      ::Info("setupDxHFE.C", Form("make chain from single chunk %s", aodPathName));
      TChain* chain = MakeAODInputChain(aodDirectory ,1, -1);
    } else {
      // Assume several folders containing different AODs. 
      // The AODs need to be in folders named 1, 2,...
      ::Info("setupDxHFE.C", Form("make chain from directory %s", aodDirectory));
      chain=MakeAODInputChain(aodDirectory, 1, 10);
    }
    ::Info("setupDxHFE.C", Form("local AOD chain: %d entries", chain->GetEntries()));
  }
}
