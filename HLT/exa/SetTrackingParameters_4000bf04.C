void SetTrackingParameters(AliLevel3 *tracker)
{
  Int_t phi_segments,eta_segments,trackletlength,tracklength;
  Int_t rowscopetracklet,rowscopetrack;
  Double_t min_pt_fit,maxangle,goodDist,hitChi2Cut,xyerror,zerror;
  Double_t goodHitChi2,trackChi2Cut,maxphi,maxeta;
  
  phi_segments = 50;//50;
  eta_segments = 100;//100;
  trackletlength = 3;
  tracklength = 10;
  rowscopetracklet = 2;
  rowscopetrack = 10;
  min_pt_fit = 0;
  maxangle = 0.1745; //AliHLTTransform::Deg2Rad(10);
  goodDist = 5;
  maxphi = 0.1;
  maxeta = 0.1;
  hitChi2Cut = 20;
  goodHitChi2 = 5;
  trackChi2Cut = 10;
      //77.8
  xyerror = -1;
  zerror =  -1;
  
  tracker->SetClusterFinderParam(xyerror,zerror,kTRUE);
  tracker->SetTrackerParam(phi_segments,eta_segments,trackletlength,tracklength,
			   rowscopetracklet,rowscopetrack,
			   min_pt_fit,maxangle,goodDist,hitChi2Cut,
			   goodHitChi2,trackChi2Cut,50,maxphi,maxeta,kTRUE);
  tracker->SetMergerParameters(2,3,0.003,0.1,0.05);
}
