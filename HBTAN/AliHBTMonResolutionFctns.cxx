#include "AliHBTMonResolutionFctns.h"

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp(AliHBTMonPxResolutionFctn)

AliHBTMonPxResolutionFctn::
AliHBTMonPxResolutionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTMonTwoParticleFctn1D(nbins,maxXval,minXval)
{
  Rename("PxResolution","PxResolution");
}
/******************************************************************/

ClassImp(AliHBTMonPyResolutionFctn)

AliHBTMonPyResolutionFctn::
AliHBTMonPyResolutionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTMonTwoParticleFctn1D(nbins,maxXval,minXval)
{
  Rename("PyResolution","PyResolution");
}
/******************************************************************/

ClassImp(AliHBTMonPzResolutionFctn)

AliHBTMonPzResolutionFctn::
AliHBTMonPzResolutionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTMonTwoParticleFctn1D(nbins,maxXval,minXval)
{
  Rename("PzResolution","PzResolution");
}
/******************************************************************/

ClassImp(AliHBTMonPResolutionFctn)

AliHBTMonPResolutionFctn::
AliHBTMonPResolutionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTMonTwoParticleFctn1D(nbins,maxXval,minXval)
{
  Rename("PResolution","PResolution");
}
/******************************************************************/

ClassImp(AliHBTMonPtResolutionFctn)

AliHBTMonPtResolutionFctn::
AliHBTMonPtResolutionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTMonTwoParticleFctn1D(nbins,maxXval,minXval)
{
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
 Rename("PxResolVsPt","Px resolution vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonPyResolutionVsPtFctn )

AliHBTMonPyResolutionVsPtFctn::
AliHBTMonPyResolutionVsPtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTMonTwoParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("PyResolVsPt","Py resolution vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonPzResolutionVsPtFctn )

AliHBTMonPzResolutionVsPtFctn::
AliHBTMonPzResolutionVsPtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTMonTwoParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("PzResolVsPt","Pz resolution vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonPResolutionVsPtFctn )

AliHBTMonPResolutionVsPtFctn::
AliHBTMonPResolutionVsPtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTMonTwoParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("PResolVsPt","P resolution vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonPtResolutionVsPtFctn )

AliHBTMonPtResolutionVsPtFctn::
AliHBTMonPtResolutionVsPtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTMonTwoParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
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
  Rename("PhiResolution","PhiResolution");
}
/******************************************************************/
ClassImp(AliHBTMonThetaResolutionFctn)

AliHBTMonThetaResolutionFctn::
AliHBTMonThetaResolutionFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
                        AliHBTMonTwoParticleFctn1D(nbins,maxXval,minXval)
{
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
 Rename("PhiResolVsPt","Phi resolution vs. Pt");
}
/******************************************************************/
ClassImp( AliHBTMonThetaResolutionVsPtFctn )

AliHBTMonThetaResolutionVsPtFctn::
AliHBTMonThetaResolutionVsPtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTMonTwoParticleFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("ThetaResolVsPt","Theta resolution vs. Pt");
}
/******************************************************************/




