/**
 * @file   RunToyModel.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Sep  1 10:54:23 2016
 * 
 * @brief  Run the toy model of the SPD trackleting 
 * 
 * @ingroup pwglf_forward_tracklets_toy
 */
/** 
 * Run the tracklet model 
 * 
 * @param opts    Options  
 * @param n       Number of tracks, steps, or events 
 * @param sigma   Cluster variance 
 * @param y1      Location of second layer 
 *
 * @relates ToyModel
 *
 * @ingroup pwglf_forward_tracklets_toy
 */
ToyModel* RunToyModel(const char* opts="help",
		   Int_t       n=20,
		   Double_t    sigma=0.0008,
		   Double_t    y1=.75)
{
  TString dir("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2");
  if (gSystem->Getenv("ANA_SRC")) dir="$ANA_SRC";

  if (!gROOT->GetClass("ToyModel"))
    gROOT->LoadMacro(Form("%s/dndeta/tracklets3/toymodel/ToyModel.C+g",
			  dir.Data()));

  return ToyModel::Run(opts, n, sigma, y1);
}
//
// EOF
//
