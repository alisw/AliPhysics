#include "AliHBTMonDistributionFctns.h"

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp(AliHBTMonPxDistributionFctn)

AliHBTMonPxDistributionFctn::AliHBTMonPxDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  Rename("Px","Px");
}

/******************************************************************/

ClassImp(AliHBTMonPyDistributionFctn)

AliHBTMonPyDistributionFctn::
AliHBTMonPyDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  Rename("Py","Py");
}

/******************************************************************/

ClassImp(AliHBTMonPzDistributionFctn)

AliHBTMonPzDistributionFctn::AliHBTMonPzDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  Rename("Pz","Pz");
}

/******************************************************************/

ClassImp(AliHBTMonPDistributionFctn)

AliHBTMonPDistributionFctn::AliHBTMonPDistributionFctn
  (Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  Rename("P","P");
}


/******************************************************************/

ClassImp(AliHBTMonPtDistributionFctn)

AliHBTMonPtDistributionFctn::AliHBTMonPtDistributionFctn
       (Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  Rename("Pt","Pt");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTMonPxDistributionVsPtFctn )

AliHBTMonPxDistributionVsPtFctn::AliHBTMonPxDistributionVsPtFctn
       (Int_t nXbins, Double_t maxXval, Double_t minXval, 
        Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTMonOneParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("PxDistVsPt","Px vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonPyDistributionVsPtFctn )

AliHBTMonPyDistributionVsPtFctn::AliHBTMonPyDistributionVsPtFctn
       (Int_t nXbins, Double_t maxXval, Double_t minXval, 
        Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTMonOneParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("PyDistVsPt","Py vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonPzDistributionVsPtFctn )

AliHBTMonPzDistributionVsPtFctn::AliHBTMonPzDistributionVsPtFctn
        (Int_t nXbins, Double_t maxXval, Double_t minXval, 
         Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTMonOneParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("PzDistVsPt","Pz vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonPDistributionVsPtFctn )

AliHBTMonPDistributionVsPtFctn::AliHBTMonPDistributionVsPtFctn
     (Int_t nXbins, Double_t maxXval, Double_t minXval, 
      Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTMonOneParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("PDistVsPt","P vs. Pt");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp(AliHBTMonPhiDistributionFctn)

AliHBTMonPhiDistributionFctn::AliHBTMonPhiDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
  AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  Rename("Phi","Phi");
}
/******************************************************************/
ClassImp(AliHBTMonThetaDistributionFctn)

AliHBTMonThetaDistributionFctn::AliHBTMonThetaDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
  AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  Rename("Theta","Theta");
}
/******************************************************************/
ClassImp( AliHBTMonPhiDistributionVsPtFctn )

AliHBTMonPhiDistributionVsPtFctn::AliHBTMonPhiDistributionVsPtFctn
     (Int_t nXbins, Double_t maxXval, Double_t minXval, 
      Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTMonOneParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("PhiDistVsPt","Phi vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonThetaDistributionVsPtFctn )

AliHBTMonThetaDistributionVsPtFctn::AliHBTMonThetaDistributionVsPtFctn
     (Int_t nXbins, Double_t maxXval, Double_t minXval, 
      Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTMonOneParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("ThetaDistVsPt","Theta vs. Pt");
}
/******************************************************************/

ClassImp(AliHBTMonVxDistributionFctn)

AliHBTMonVxDistributionFctn::
AliHBTMonVxDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  Rename("Vx","X of Vertex");
}
/******************************************************************/

ClassImp(AliHBTMonVyDistributionFctn)

AliHBTMonVyDistributionFctn::
AliHBTMonVyDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  Rename("Vy","Y of Vertex");
}
/******************************************************************/

ClassImp(AliHBTMonVzDistributionFctn)

AliHBTMonVzDistributionFctn::
AliHBTMonVzDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  Rename("Vz","Z of Vertex");
}
/******************************************************************/
ClassImp(AliHBTMonRDistributionFctn)

AliHBTMonRDistributionFctn::
AliHBTMonRDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  Rename("VertexDistanceFromCenter","Distance of Particle Vertex From Center");
}
/******************************************************************/

ClassImp(AliHBTMonVyDistributionVsVxFctn)
AliHBTMonVyDistributionVsVxFctn::AliHBTMonVyDistributionVsVxFctn
     (Int_t nXbins, Double_t maxXval, Double_t minXval, 
      Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTMonOneParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("VyDistVsVx","Vrtex Y position versus X vertex position");
}
/******************************************************************/

ClassImp(AliHBTMonRtDistributionVsVzFctn)
AliHBTMonRtDistributionVsVzFctn::AliHBTMonRtDistributionVsVzFctn
     (Int_t nXbins, Double_t maxXval, Double_t minXval, 
      Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTMonOneParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("RDistVsVz","Distance of vertex position from center in trensverse plane versus Z vertex position");
}

/******************************************************************/
/******************************************************************/




