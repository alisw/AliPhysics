#include "AliHBTQResolutionFctns.h"


/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTQInvResolVsQInvFctn )
AliHBTQInvResolVsQInvFctn::
AliHBTQInvResolVsQInvFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                          Int_t nYbins, Double_t maxYval, Double_t minYval):
	         AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QInvResolVsQInv","Q_{Inv} Resolution vs. Q_{Inv}");
}
/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQOutResolVsQInvFctn )
AliHBTQOutResolVsQInvFctn::
AliHBTQOutResolVsQInvFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                          Int_t nYbins, Double_t maxYval, Double_t minYval):
	         AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QOutResolVsQInv","Q_{Out} Resolution vs. Q_{Inv}");
}
/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQSideResolVsQInvFctn )

AliHBTQSideResolVsQInvFctn::
AliHBTQSideResolVsQInvFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                          Int_t nYbins, Double_t maxYval, Double_t minYval):
	         AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QSideResolVsQInv","Q_{Side} Resolution vs. Q_{Inv}");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQLongResolVsQInvFctn )

AliHBTQLongResolVsQInvFctn::
AliHBTQLongResolVsQInvFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QLongResolVsQInv","Q_{Long} Resolution vs. Q_{Inv}");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQInvResolVsKtFctn )

AliHBTQInvResolVsKtFctn::
AliHBTQInvResolVsKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QInvResolVsKt","Q_{Inv} Resolution vs. K_{t}");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTQOutResolVsKtFctn )

AliHBTQOutResolVsKtFctn::
AliHBTQOutResolVsKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QOutResolVsKt","Q_{Out} Resolution vs. K_{t} ");
}


/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTQSideResolVsKtFctn )

AliHBTQSideResolVsKtFctn::
AliHBTQSideResolVsKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QSideResolVsKt","Q_{Side} Resolution vs. K_{t} ");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTQLongResolVsKtFctn )

AliHBTQLongResolVsKtFctn::
AliHBTQLongResolVsKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                             AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QLongResolVsKt","Q_{Long} Resolution vs. K_{t} ");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTQOutResolVsQOutFctn)

AliHBTQOutResolVsQOutFctn::
AliHBTQOutResolVsQOutFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QOutResolVsQOut","Q_{Out} Resolution vs. Q_{Out} ");
}
 
/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTQSideResolVsQSideFctn )

AliHBTQSideResolVsQSideFctn::
AliHBTQSideResolVsQSideFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QSideResolVsQSide","Q_{Side} Resolution vs. Q_{Side} ");
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTQLongResolVsQLongFctn )

AliHBTQLongResolVsQLongFctn::
AliHBTQLongResolVsQLongFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("QLongResolVsQLong","Q_{Long} Resolution vs. Q_{Long} ");
}



/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTPairThetaResolVsQInvFctn )

AliHBTPairThetaResolVsQInvFctn::
AliHBTPairThetaResolVsQInvFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("PairThetaResolVsQInv","Pair Theta Angle Resolution vs. Q_{Inv} ");
}
/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTPairPhiResolVsQInvFctn )

AliHBTPairPhiResolVsQInvFctn::
AliHBTPairPhiResolVsQInvFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("PairPhiResolVsQInv","Pair Phi Angle Resolution vs. Q_{Inv} ");
}
/******************************************************************/
/******************************************************************/
/******************************************************************/


ClassImp( AliHBTPairThetaResolVsKtFctn )

AliHBTPairThetaResolVsKtFctn::
AliHBTPairThetaResolVsKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("PairThetaResolVsKt","Pair Theta Angle Resolution vs. K_{t} ");
}
/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTPairPhiResolVsKtFctn )

AliHBTPairPhiResolVsKtFctn::
AliHBTPairPhiResolVsKtFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
 Rename("PairPhiResolVsKt","Pair Phi Angle Resolution vs. K_{t} ");
}
/******************************************************************/
/******************************************************************/
/******************************************************************/


