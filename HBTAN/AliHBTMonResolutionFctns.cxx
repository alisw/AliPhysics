#include "AliHBTMonResolutionFctns.h"
//_______________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////////
//
// class AliHBTMonPxResolutionFctn;
// class AliHBTMonPyResolutionFctn;
// class AliHBTMonPzResolutionFctn;
// class AliHBTMonPResolutionFctn;
// class AliHBTMonPtResolutionFctn;
// class AliHBTMonPhiResolutionFctn;
// class AliHBTMonThetaResolutionFctn;
// class AliHBTMonPxResolutionVsPtFctn;
// class AliHBTMonPyResolutionVsPtFctn;
// class AliHBTMonPzResolutionVsPtFctn;
// class AliHBTMonPResolutionVsPtFctn;
// class AliHBTMonPtResolutionVsPtFctn;
// class AliHBTMonPhiResolutionVsPtFctn;
// class AliHBTMonThetaResolutionVsPtFctn;
//
// Caution: On 2D plots on X axis in simulated values
// That is contrary to two-particle resolutions where it is reconstructed one
//
// added by Zbigniew.Chajecki@cern.ch
// this classes create resolution functions of particle momentum 
//
//////////////////////////////////////////////////////////////////////////////////

ClassImp(AliHBTMonPxResolutionFctn)

AliHBTMonPxResolutionFctn::
AliHBTMonPxResolutionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTMonTwoParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor
  Rename("PxResolution","PxResolution");
}
/******************************************************************/

ClassImp(AliHBTMonPyResolutionFctn)

AliHBTMonPyResolutionFctn::
AliHBTMonPyResolutionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTMonTwoParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor
  Rename("PyResolution","PyResolution");
}
/******************************************************************/

ClassImp(AliHBTMonPzResolutionFctn)

AliHBTMonPzResolutionFctn::
AliHBTMonPzResolutionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTMonTwoParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor
  Rename("PzResolution","PzResolution");
}
/******************************************************************/

ClassImp(AliHBTMonPResolutionFctn)

AliHBTMonPResolutionFctn::
AliHBTMonPResolutionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTMonTwoParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor
  Rename("PResolution","PResolution");
}
/******************************************************************/

ClassImp(AliHBTMonPtResolutionFctn)

AliHBTMonPtResolutionFctn::
AliHBTMonPtResolutionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTMonTwoParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor
  Rename("PtResolution","PtResolution");
}
/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTMonPxResolutionVsPtFctn )

AliHBTMonPxResolutionVsPtFctn::
AliHBTMonPxResolutionVsPtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTMonTwoParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
  //constructor
 Rename("PxResolVsPt","Px resolution vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonPyResolutionVsPtFctn )

AliHBTMonPyResolutionVsPtFctn::
AliHBTMonPyResolutionVsPtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTMonTwoParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
  //constructor
 Rename("PyResolVsPt","Py resolution vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonPzResolutionVsPtFctn )

AliHBTMonPzResolutionVsPtFctn::
AliHBTMonPzResolutionVsPtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTMonTwoParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
  //constructor
 Rename("PzResolVsPt","Pz resolution vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonPResolutionVsPtFctn )

AliHBTMonPResolutionVsPtFctn::
AliHBTMonPResolutionVsPtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTMonTwoParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
  //constructor
 Rename("PResolVsPt","P resolution vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonPtResolutionVsPtFctn )

AliHBTMonPtResolutionVsPtFctn::
AliHBTMonPtResolutionVsPtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTMonTwoParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
  //constructor
 Rename("PtResolVsPt","Pt resolution vs. Pt");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp(AliHBTMonPhiResolutionFctn)

AliHBTMonPhiResolutionFctn::
AliHBTMonPhiResolutionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTMonTwoParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor
  Rename("PhiResolution","PhiResolution");
}
/******************************************************************/
ClassImp(AliHBTMonThetaResolutionFctn)

AliHBTMonThetaResolutionFctn::
AliHBTMonThetaResolutionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTMonTwoParticleFctn1D(nbins,maxXval,minXval)
{
  //constructor
  Rename("ThetaResolution","ThetaResolution");
}
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTMonPhiResolutionVsPtFctn )

AliHBTMonPhiResolutionVsPtFctn::
AliHBTMonPhiResolutionVsPtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTMonTwoParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
  //constructor
 Rename("PhiResolVsPt","Phi resolution vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonThetaResolutionVsPtFctn )

AliHBTMonThetaResolutionVsPtFctn::
AliHBTMonThetaResolutionVsPtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTMonTwoParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 //constructor
 Rename("ThetaResolVsPt","Theta resolution vs. Pt");
}
/******************************************************************/


ClassImp( AliHBTMonPhiResolutionVsPhiFctn )

AliHBTMonPhiResolutionVsPhiFctn::
AliHBTMonPhiResolutionVsPhiFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTMonTwoParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
  //constructor
 Rename("PhiResolVsPhi","Phi resolution vs. Phi");
}
/******************************************************************/
ClassImp( AliHBTMonThetaResolutionVsThetaFctn )

AliHBTMonThetaResolutionVsThetaFctn::
AliHBTMonThetaResolutionVsThetaFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTMonTwoParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
  //constructor
 Rename("ThetaResolVsTheta","Theta resolution vs. Theta");
}
/******************************************************************/


