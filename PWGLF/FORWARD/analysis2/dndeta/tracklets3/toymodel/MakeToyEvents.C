/**
 * @file   MakeToyEvents.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Aug 24 22:25:45 2016
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_tracklets_toy
 */
/** 
 * Make some events 
 * 
 * @param full If true, make 500 tracks per event, otherwise 476
 * tracks per event
 *
 * @ingroup pwglf_forward_tracklets_toy
 */
void
MakeToyEvents(Bool_t full=true)
{
  if (!gROOT->GetClass("Model"))
    gROOT->LoadMacro("$ANA_SRC/dndeta/tracklets3/toymodel/TrackletModel.C+g");


  Model* m = new Model(0.75);
  m->SetOut2In(true);     // Track outside-in
  m->SetGlobal(true);     // Select tracklets globally
  m->SetNoReuse(true);    // Do not reuse clusters 
  m->SetSigma(0.0008);    // Cluster size
  m->SetNoise(false);     // Do not add noise
  m->SetVary(false);      // Do not vary size of clusters
  m->SetMerge(false);     // Do not merge close-adjecent clusters
  m->SetSecondary(false); // Do not make secondaries 

  Int_t nEvents = 500;
  Int_t nTracks = (full ? 500 : 476);
  m->RunDist(nEvents, nTracks);

  gSystem->Rename("dist.root", full ? "real.root" : "reduced.root");
}
//
// EOF
//
