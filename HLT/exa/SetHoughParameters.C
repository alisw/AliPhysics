//$Id$

void SetHoughParameters(AliL3Hough *hough,Char_t *path)
{
  
  bool binary = kTRUE;    //binary files input
  int n_eta_segments=100;
  double histptmin = 0.5; //mininum pt to find (controls the histogram range)

  int threshold=6000;  //peak threshold
  int nxbins = 140;
  int nybins = 150;
  
  int patch=-1; //-1 -> Hough transform on slices (means adding histograms)
  
  hough->SetThreshold(4); //noise threshold on single digits
  hough->SetTransformerParams(nxbins,nybins,histptmin,patch);

/*
  hough->SetTransformerParams(0.068,histptmin,histptmax,20,0);
  hough->SetTransformerParams(0.043,histptmin,histptmax,23,1);
  hough->SetTransformerParams(0.03,histptmin,histptmax,25,2);
  hough->SetTransformerParams(0.02,histptmin,histptmax,25,3);
  hough->SetTransformerParams(0.016,histptmin,histptmax,30,4);
  hough->SetTransformerParams(0.012,histptmin,histptmax,30,5);
*/
  hough->SetPeakThreshold(threshold,patch);
  hough->Init(path,binary,n_eta_segments);
}
