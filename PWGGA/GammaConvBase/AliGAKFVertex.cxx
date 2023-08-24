#include "AliLog.h"
#include "AliGAKFVertex.h"
#ifdef PWGGAUSEKFPARTICLE
ClassImp(AliGAKFVertex)

AliGAKFVertex::AliGAKFVertex(const AliVVertex&vertex ) : KFVertex()  {
  //    AliError("THIS CLASS MEMBER NEEDS TO BE IMPLEMENTED");
  double myfP_t[3];
  double myfC_t[6]; 

  vertex.GetXYZ( myfP_t );
  vertex.GetCovarianceMatrix( myfC_t );  

  fChi2 = vertex.GetChi2();  

  fNDF = 2*vertex.GetNContributors() - 3;


  for (int i = 0; i < 3; ++i) fP[i] = float(myfP_t[i]);
  for (int i = 0; i < 6; ++i) fC[i] = float(myfC_t[i]);

  fQ = 0;
  fAtProductionVertex = 0;
  fSFromDecay = 0;

}
#endif // PWGGAUSEKFPARTICLE
