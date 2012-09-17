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
///   locakl AODs

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// environment specific for DxHFE
//
const char* includePath="-I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWGHF/hfe";
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
    // TODO: decide depending on running mode and input to be used
    gROOT->LoadMacro("$ALICE_ROOT/PWGHF/vertexingHF/MakeAODInputChain.C");
    // Create a chain with one set of AliAOD.root and AliAOD.VertexingHF.root. The set needs 
    // to be located in the same folder as you run from (physically or linked)
    TChain* chain = MakeAODInputChain(aodDirectory ,1, -1);
    // If you have several folders containing different AODs, use below. 
    // From the MakeAODInputChain.C: The AODs need to be in folders named 1, 2,...
    //if(useMC)chain =MakeAODInputChain("/scratch/Data/2010/MC/LHC10f7a/130375/AOD051/",1,2);//73
    //chain =MakeAODInputChain("/scratch/Data/2010/MC/LHC10f6a/126437/AOD041/",1,10);//84
    //else chain =MakeAODInputChain("/scratch/Data/2010/LHC10d/000126437/ESDs/pass2/AOD057/",1,10);//LHC10f7a/130375/AOD051/",1,20);
    cout << "local AOD chain: " << chain->GetEntries() << " entries" << endl;
  }
}
