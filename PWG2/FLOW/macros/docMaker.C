{
 gSystem->Load("libPhysics.so") ;
 gSystem->Load("libPWG2flow.so") ;
 THtml gd ;  
 char *htmlPath = gSystem->ExpandPathName("$ALICE_ROOT/PWG2/FLOW/htmldoc/.") ;
 gd.SetOutputDir(htmlPath) ;
 gd.MakeClass("AliFlowConstants") ;
 gd.MakeClass("AliFlowTrack") ;
 gd.MakeClass("AliFlowEvent") ;
 gd.MakeClass("AliFlowV0") ;
 gd.MakeClass("AliFlowSelection") ;
 gd.MakeClass("AliFlowMaker") ;
 gd.MakeClass("AliFlowKineMaker") ;
 gd.MakeClass("AliFlowWeighter") ;
 gd.MakeClass("AliFlowAnalyser") ;
// gd.MakeClass("AliFlow...") ;
// gd.MakeIndex() ;
}
