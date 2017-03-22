// this macro compiles and runs AnalyzeTree.cxx, which takes as input the 
// starlight.root file produced by convertStarlightAsciiToTree.cxx
// output histograms are stored in starlight_histos.root 
//
void AnaTree(){
gROOT->ProcessLine(".L AnalyzeTree.cxx");
AnalyzeTree* l = new AnalyzeTree();
l->Loop();
}
