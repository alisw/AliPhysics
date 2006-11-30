//$Id$

void SetHoughParameters(AliHLTHough *hough,Char_t *path,Int_t tversion=1)
{
  bool binary = kTRUE;    //binary files input
  int n_eta_segments=100;
  double histptmin = 0.5; //mininum pt to find (controls the histogram range)

  int threshold=6000;  //peak threshold
  int nxbins = 140;
  int nybins = 150;
  //Int_t threshold=5000;  //peak threshold
  //Int_t nxbins = 190;
  //Int_t nybins = 200;
  
  int patch=-1; //-1 -> Hough transform on slices (means adding histograms)
  
  hough->SetThreshold(4); //noise threshold on single digits
  hough->SetTransformerParams(nxbins,nybins,histptmin,patch);
  hough->SetPeakThreshold(threshold,patch);
  hough->Init(path,binary,n_eta_segments,kFALSE,tversion);
}
