//$Id$

void SetFitParameters(AliL3ClusterFitter *fitter)
{

  fitter->SetInnerWidthFactor(1,1.5);
  fitter->SetOuterWidthFactor(1,1.5);
  fitter->SetNmaxOverlaps(5);
  
  //fitter->SetChiSqMax(5,kFALSE); //isolated clusters
  fitter->SetChiSqMax(5,kTRUE);  //overlapping clusters

}
