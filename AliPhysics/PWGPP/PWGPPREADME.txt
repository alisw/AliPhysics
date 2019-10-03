

How to make a comparison using MC information:


0. Build the MC info tree
gSystem->Load("libANALYSIS");
gSystem->Load("libPWGPP")
AliGenInfoMaker *t = new AliGenInfoMaker("galice.root","genTracks.root",0,0)
t->Exec();


1. Build the reconstructed info tree

gSystem->Load("libANALYSIS");
gSystem->Load("libPWGPP");
//
AliRecInfoMaker *t2 = new AliRecInfoMaker("genTracks.root","cmpESDTracks.root","galice.root",0,0);
t2->Exec();



2. Make a chain of the information tree
gSystem->Load("libANALYSIS");

gSystem->Load("libPWGPP");

//GSI example
.x ~/rootlogon.C
 gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros")
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")

 AliXRDPROOFtoolkit tool;
 TChain * chain = tool.MakeChain("cmp.txt","ESDcmpTracks",0,1000)
 chain->Lookup();
 .L $ALICE_PHYSICS/PWGPP/AliComparisonSelector.cxx+

 
3. 

a.) Use AliTreeDraw for fast prototyping the queries - analysis type:

gSystem->Load("libPWGPP");
AliTreeDraw comp;
comp.SetTree(tree)


b.) Or use Selector

chain->Process("$ALICE_PHYSICS/PWGPP/AliComparisonSelector.cxx+")


TFile f("Output.root");
