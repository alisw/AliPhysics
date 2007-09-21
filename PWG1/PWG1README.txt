

How to make a comparison using MC information:


0. Build the MC info tree

gSystem->Load("libPWG1.so")
AliGenInfoMaker *t = new AliGenInfoMaker("galice.root","genTracks.root",0,0)
t->Exec();


1. Build the reconstructed info tree


gSystem->Load("libPWG1.so");
//
AliRecInfoMaker *t2 = new AliRecInfoMaker("genTracks.root","cmpESDTracks.root","galice.root",0,0);
t2->Exec();


2. Make a chain of the information tree

gSystem->Load("libPWG1.so");

//GSI example
.x ~/rootlogon.C
.L /u/miranov/macroxrdproof64/AliXRDPROOFtoolkit.cxx+
 AliXRDPROOFtoolkit tool;
 TChain * chain = tool.MakeChain("comp.txt","ESDcmpTracks",0,1000)
 chain->Lookup();
 .L $ALICE_ROOT/PWG1/AliComparisonSelector.cxx+

 
3. 

a.) Use AliTreeDraw for fast prototyping the queries - analysis type:

gSystem->Load("libPWG1.so");
AliTreeDraw comp;
comp.SetTree(tree)


b.) Or use Selector

chain->Process("$ALICE_ROOT/PWG1/AliComparisonSelector.cxx+")


TFile f("Output.root");
