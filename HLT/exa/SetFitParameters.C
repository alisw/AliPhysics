//$Id$

void SetFitParameters(AliHLTClusterFitter *fitter)
{
  fitter->SetInnerWidthFactor(1,1);
  fitter->SetOuterWidthFactor(1,1);
  fitter->SetNmaxOverlaps(5);

  fitter->SetChiSqMax(5,0);
  fitter->SetChiSqMax(5,1);
  fitter->SetChiSqMax(5,2);
}
