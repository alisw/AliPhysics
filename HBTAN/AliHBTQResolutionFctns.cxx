#include "AliHBTQResolutionFctns.h"

AliHBTQOutResolVSQInvFctn* xxx = new AliHBTQOutResolVSQInvFctn();
/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQOutResolVSQInvFctn )
AliHBTQOutResolVSQInvFctn::
AliHBTQOutResolVSQInvFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                          Int_t nYbins, Double_t maxYval, Double_t minYval):
	         AliHBTFourPartFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QOutResolVSQInv","Q_{Out} Resolution vs. Q_{Inv}");
}
/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQSideResolVSQInvFctn )

AliHBTQSideResolVSQInvFctn::
AliHBTQSideResolVSQInvFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                          Int_t nYbins, Double_t maxYval, Double_t minYval):
	         AliHBTFourPartFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QSideResolVSQInv","Q_{Side} Resolution vs. Q_{Inv}");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQLongResolVSQInvFctn )

AliHBTQLongResolVSQInvFctn::
AliHBTQLongResolVSQInvFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTFourPartFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QLongResolVSQInv","Q_{Long} Resolution vs. Q_{Inv}");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQInvResolVSKtFctn )

AliHBTQInvResolVSKtFctn::
AliHBTQInvResolVSKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTFourPartFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QInvResolVSKt","Q_{Inv} Resolution vs. K_{t}");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTQOutResolVSKtFctn )

AliHBTQOutResolVSKtFctn::
AliHBTQOutResolVSKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTFourPartFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QOutResolVSKt","Q_{Out} Resolution vs. K_{t} ");
}


/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTQSideResolVSKtFctn )

AliHBTQSideResolVSKtFctn::
AliHBTQSideResolVSKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTFourPartFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QSideResolVSKt","Q_{Side} Resolution vs. K_{t} ");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTQLongResolVSKtFctn )

AliHBTQLongResolVSKtFctn::
AliHBTQLongResolVSKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTFourPartFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QLongResolVSKt","Q_{Long} Resolution vs. K_{t} ");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQOutResolVSQOutFctn)

AliHBTQOutResolVSQOutFctn::
AliHBTQOutResolVSQOutFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTFourPartFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QOutResolVSQOut","Q_{Out} Resolution vs. Q_{Out} ");
}
 
/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTQSideResolVSQSideFctn )

AliHBTQSideResolVSQSideFctn::
AliHBTQSideResolVSQSideFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTFourPartFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QSideResolVSQSide","Q_{Side} Resolution vs. Q_{Side} ");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTQLongResolVSQLongFctn )

AliHBTQLongResolVSQLongFctn::
AliHBTQLongResolVSQLongFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTFourPartFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QLongResolVSQLong","Q_{Long} Resolution vs. Q_{Long} ");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/




