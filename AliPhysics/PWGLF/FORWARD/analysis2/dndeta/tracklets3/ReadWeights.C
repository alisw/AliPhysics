/**
 * @file   ReadWeights.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Sep 20 17:18:29 2016
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_tracklets
 * 
 */
/** 
 * @{ 
 * @name Read in weights 
 *
 * @ingroup pwglf_forward_tracklets
 */
/** 
 * Read in weights and plot
 * 
 * @param filename 
 * @relates AliTrackletBaseWeights
 */
void ReadWeights(const char* filename="weights.root")
{
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include -I${ALICE_PHYSICS}/include");
  gROOT->LoadMacro("AliAODTracklet.C+g");
  gROOT->LoadMacro("AliTrackletWeights.C+g");

  TFile* file = TFile::Open(filename,"READ");
  file->ls();
  AliTrackletBaseWeights* weights =
    static_cast<AliTrackletBaseWeights*>(file->Get("weights"));
  TCanvas* c = new TCanvas("c","c");
  weights->Draw();
}

/* @} */
//
// EOF
// 
