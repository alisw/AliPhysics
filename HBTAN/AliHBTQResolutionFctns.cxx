//__________________________________________________________________
////////////////////////////////////////////////////////////////////
//                                                                //
// General Remark:                                                //
// CMSLC means                                                    //
// Center of Mass System Longitudially Co-moving                  //
//                                                                //
//                                                                //
// This class creates resolution function of Qout                 //
// (difference of simulated pair Qout and recontructed pair)      //
// in function of QInv                                            //
// it inherits from AliHBTTwoPairFctn2D                           //
//  it needs two pairs to compare                                 //
//  and is two dimentional: numerator and denominator are TH2D    //
//                                                                //
////////////////////////////////////////////////////////////////////

#include "AliHBTQResolutionFctns.h"


/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTKtResolVsQInvFctn )
AliHBTKtResolVsQInvFctn::
AliHBTKtResolVsQInvFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                        Int_t nYbins, Double_t maxYval, Double_t minYval):
 AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
//ctor
 Rename("KtResolVsQInv","K_{t} Resolution vs. Q_{Inv}");
}
/******************************************************************/
/******************************************************************/
/******************************************************************/
ClassImp( AliHBTQInvResolVsQInvFctn )
AliHBTQInvResolVsQInvFctn::
AliHBTQInvResolVsQInvFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                          Int_t nYbins, Double_t maxYval, Double_t minYval):
	         AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
//ctor
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
//ctor
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
//ctor
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
//ctor
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
//ctor
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
//ctor
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
//ctor
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
//ctor
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
//ctor
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
//ctor
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
//ctor
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
//ctor
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
//ctor
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
//ctor
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
//ctor
 Rename("PairPhiResolVsKt","Pair Phi Angle Resolution vs. K_{t} ");
}
/******************************************************************/
/******************************************************************/
/******************************************************************/


