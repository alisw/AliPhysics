#include "AliHBTMonDistributionFctns.h"
//______________________________________________________________
////////////////////////////////////////////////////////////////
//
// class AliHBTMonPxDistributionFctn;
// class AliHBTMonPxDistributionVsPtFctn;
// class AliHBTMonPyDistributionFctn;
// class AliHBTMonPyDistributionVsPtFctn;
// class AliHBTMonPzDistributionFctn;
// class AliHBTMonPzDistributionVsPtFctn;
// class AliHBTMonPDistributionFctn;
// class AliHBTMonPDistributionVsPtFctn;
// class AliHBTMonPtDistributionFctn;
// class AliHBTMonVxDistributionFctn;
// class AliHBTMonVyDistributionFctn;
// class AliHBTMonVzDistributionFctn;
// class AliHBTMonRDistributionFctn;
// class AliHBTMonVyDistributionVsVxFctn;
// class AliHBTMonRtDistributionVsVzFctn;
//
// added by Zbigniew.Chajecki@cern.ch
// this classes create distribution functions of particle momentum
//
/////////////////////////////////////////////////////////////////
/******************************************************************/

ClassImp(AliHBTMonPxDistributionFctn)

AliHBTMonPxDistributionFctn::AliHBTMonPxDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor
  Rename("Px","Px");
}

Double_t AliHBTMonPxDistributionFctn::GetValue(AliVAODParticle * particle) const
{ 
 //returns value for that function
 return particle->Px();
}

/******************************************************************/

ClassImp(AliHBTMonPyDistributionFctn)

AliHBTMonPyDistributionFctn::AliHBTMonPyDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor
  Rename("Py","Py");
}

/******************************************************************/

ClassImp(AliHBTMonPzDistributionFctn)

AliHBTMonPzDistributionFctn::AliHBTMonPzDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor  
  Rename("Pz","Pz");
}

/******************************************************************/

ClassImp(AliHBTMonPDistributionFctn)

AliHBTMonPDistributionFctn::AliHBTMonPDistributionFctn
  (Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor
  Rename("P","P");
}


/******************************************************************/

ClassImp(AliHBTMonPtDistributionFctn)

AliHBTMonPtDistributionFctn::AliHBTMonPtDistributionFctn
       (Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor
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
  //constructor
 Rename("PxDistVsPt","Px vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonPyDistributionVsPtFctn )

AliHBTMonPyDistributionVsPtFctn::AliHBTMonPyDistributionVsPtFctn
       (Int_t nXbins, Double_t maxXval, Double_t minXval, 
        Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTMonOneParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
  //constructor
  Rename("PyDistVsPt","Py vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonPzDistributionVsPtFctn )

AliHBTMonPzDistributionVsPtFctn::AliHBTMonPzDistributionVsPtFctn
        (Int_t nXbins, Double_t maxXval, Double_t minXval, 
         Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTMonOneParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 //constructor  
 Rename("PzDistVsPt","Pz vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonPDistributionVsPtFctn )

AliHBTMonPDistributionVsPtFctn::AliHBTMonPDistributionVsPtFctn
     (Int_t nXbins, Double_t maxXval, Double_t minXval, 
      Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTMonOneParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
  //constructor
  Rename("PDistVsPt","P vs. Pt");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp(AliHBTMonPhiDistributionFctn)

AliHBTMonPhiDistributionFctn::AliHBTMonPhiDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
  AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor  
  Rename("Phi","Phi");
}
/******************************************************************/
ClassImp(AliHBTMonThetaDistributionFctn)

AliHBTMonThetaDistributionFctn::AliHBTMonThetaDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
  AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor  
  Rename("Theta","Theta");
}
/******************************************************************/
ClassImp( AliHBTMonPhiDistributionVsPtFctn )

AliHBTMonPhiDistributionVsPtFctn::AliHBTMonPhiDistributionVsPtFctn
     (Int_t nXbins, Double_t maxXval, Double_t minXval, 
      Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTMonOneParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
  //constructor  
 Rename("PhiDistVsPt","Phi vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonThetaDistributionVsPtFctn )

AliHBTMonThetaDistributionVsPtFctn::AliHBTMonThetaDistributionVsPtFctn
     (Int_t nXbins, Double_t maxXval, Double_t minXval, 
      Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTMonOneParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
  //constructor
   Rename("ThetaDistVsPt","Theta vs. Pt");
}
/******************************************************************/

ClassImp(AliHBTMonVxDistributionFctn)

AliHBTMonVxDistributionFctn::
AliHBTMonVxDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor  
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
  //constructor  
  Rename("Vz","Z of Vertex");
}
/******************************************************************/
ClassImp(AliHBTMonRDistributionFctn)

AliHBTMonRDistributionFctn::
AliHBTMonRDistributionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTMonOneParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor
    Rename("VertexDistanceFromCenter","Distance of Particle Vertex From Center");
}
/******************************************************************/

ClassImp(AliHBTMonVyDistributionVsVxFctn)
AliHBTMonVyDistributionVsVxFctn::AliHBTMonVyDistributionVsVxFctn
     (Int_t nXbins, Double_t maxXval, Double_t minXval, 
      Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTMonOneParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
  //constructor
  Rename("VyDistVsVx","Vrtex Y position versus X vertex position");
}
/******************************************************************/

ClassImp(AliHBTMonRtDistributionVsVzFctn)
AliHBTMonRtDistributionVsVzFctn::AliHBTMonRtDistributionVsVzFctn
     (Int_t nXbins, Double_t maxXval, Double_t minXval, 
      Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTMonOneParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
  //constructor
  Rename("RDistVsVz","Distance of vertex position from center in trensverse plane versus Z vertex position");
}

/******************************************************************/
/******************************************************************/




