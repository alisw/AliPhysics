#include "AliPHOSClusterSelection.h"

AliPHOSClusterSelection::AliPHOSClusterSelection()
  : fMinChargedParticleTrackDistance(-1.), 
    fNotUnfolded(false),
    fMaxDispR2(-1.),
    fMaxDispCoreR2(-1.),
    fMaxTOF(-1.)
{
  // Defaults to the most lenient selection allowable
	return;
}

AliPHOSClusterSelection::~AliPHOSClusterSelection()
{
}

Bool_t AliPHOSClusterSelection::IsSelected(AliVCluster* cluster) const
{
  return IsSelectedCPV(cluster)
    && IsSelectedUnfolded(cluster)
    && IsSelectedDisp(cluster)
    && IsSelectedDispCore(cluster)
    && IsSelectedTOF(cluster);
}

Bool_t IsSelectedCPV(AliVCluster* cluster) const
{
  if( 0 > SetMinChargedParticleTrackDistance )
    return true; 
  //TODO: implement positive case
}

AliPHOSClusterSelection* AliPHOSClusterSelection::SetMinChargedParticleTrackDistance(Float_t distance)
{
  // 'distance' set the minimal allowable distance between the cluster
  // and the nearest extrapolated track.
  // If 'distance' is negative, then all clusters are sellected, the selection
  // being "not applied" or "disabled".
  
  fMinChargedParticleTrackDistance = distance;
}

TString AliPHOSClusterSelection::ToString() const
{
  // returns a string an quasi-unique string for whatever selection 
  // parameters the instance contains. The uniqueness of the string
  // is limited by the precision given in the formatting of the string.
  // Take care that the precision is sufficient for your needs.

  return TString::Format("%f_%i_%f_%f_%f",
			 fMinChargedParticleTrackDistance,
			 fNotUnfolded,
			 fMaxDispR2,
			 fMaxDispCoreR2,
			 fMaxTOF
			 );
}


Float_t AliPHOSClusterSelection::SetMinChargedParticleTrackDistance(const TString& string)
{
  TObjArray * objarray = string.Tokenize("_");
  Float_t flt = objarray->At(0)->Atof();
  delete objarray;
  return flt;
}
