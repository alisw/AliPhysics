#include "AliHBTQDistributionFctns.h"

/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQInvDistributionVsKtFctn )

AliHBTQInvDistributionVsKtFctn::
AliHBTQInvDistributionVsKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTOnePairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QInvDistributionVsKt","Q_{Inv} Distribution vs. K_{t}");
}

/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQOutDistributionVsKtFctn )

AliHBTQOutDistributionVsKtFctn::
AliHBTQOutDistributionVsKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTOnePairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QOutDistributionVsKt","Q_{Out} Distribution vs. K_{t}");
}

/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQSideDistributionVsKtFctn )

AliHBTQSideDistributionVsKtFctn::
AliHBTQSideDistributionVsKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTOnePairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QSideDistributionVsKt","Q_{Side} Distribution vs. K_{t}");
}

/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQLongDistributionVsKtFctn )

AliHBTQLongDistributionVsKtFctn::
AliHBTQLongDistributionVsKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTOnePairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QLongDistributionVsKt","Q_{Long} Distribution vs. K_{t}");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/




