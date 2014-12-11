/**
 * @file   RunUnfoldMult.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Nov 12 09:25:07 2013
 * 
 * @brief  Run the unfolding.
 * 
 * This wrapper is here to load RooUnfold first
 * @ingroup pwglf_forward_multdist
 */

/** 
 * Run the unfolding 
 * 
 * @ingroup pwglf_forward_multdist
 */
void
RunUnfoldMult()
{
  TString rooUnfold = gSystem->Getenv("ROOUNFOLD");
  if (!rooUnfold.IsNull()) {
    gSystem->AddIncludePath(Form("-I%s/src", rooUnfold.Data()));
    gSystem->AddDynamicPath(rooUnfold);
  }
  gSystem->Load("libRooUnfold");
  gROOT->Macro("UnfoldMult.C++");
}
/*
 * EOF
 */
