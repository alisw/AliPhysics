#include "AliFlatESDVertex.h"

void AliFlatESDVertex::Set(const AliESDVertex &v )
{
  fPosition[0] = v.GetX();
  fPosition[1] = v.GetY();
  fPosition[2] = v.GetZ();
  Double_t c[6]; 
  v.GetCovarianceMatrix( c ); 
  for( int i=0; i<6; i++) fCov[i] = c[i];
  fNContributors = v.GetNContributors();
  fChi2 = v.GetChi2();
}
