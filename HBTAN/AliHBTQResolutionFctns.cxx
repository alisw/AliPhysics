#include "AliHBTQResolutionFctns.h"

//__________________________________________________________________
////////////////////////////////////////////////////////////////////
//                                                                //
// General Remark:                                                //
// LCMS means                                                    //
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

void AliHBTKtResolVsQInvFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
{
//returns values of the functiion  
  y = partpair->GetKt() - trackpair->GetKt();
  x = partpair->GetQInv();
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
void AliHBTQInvResolVsQInvFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
{
//returns values of the functiion  
 y = partpair->GetQInv() - trackpair->GetQInv();
 x = partpair->GetQInv();
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

void AliHBTQOutResolVsQInvFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
{
  //returns Qoutsim-Qoutrec for y
  //returns Qinv for x
  Double_t tqout = trackpair->GetQOutLCMS();
  y = partpair->GetQOutLCMS() - tqout;
  if (tqout < 0.0) y = -y;
  x = partpair->GetQInv();
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

void AliHBTQSideResolVsQInvFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair,  Double_t& x, Double_t& y) const
{
  //returns Qsidesim-Qsiderec for y
  //returns Qinv for x
  y = partpair->GetQSideLCMS() - trackpair->GetQSideLCMS();
  if (trackpair->GetQSideLCMS() < 0.0) y = -y;
  x = partpair->GetQInv();
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

void AliHBTQLongResolVsQInvFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
{
  //returns Qlongsim-Qlongrec for y
  //returns Qinv for x
  y = partpair->GetQLongLCMS() - trackpair->GetQLongLCMS();
  if (trackpair->GetQLongLCMS() < 0.0) y = -y;
  x = partpair->GetQInv();
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

void AliHBTQInvResolVsKtFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
{
 //returns values of the function
 y = partpair->GetQInv() - trackpair->GetQInv();
 x = partpair->GetKt();
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

void AliHBTQOutResolVsKtFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
{
  //returns Qoutsim-Qoutrec for y
  //returns Kt for x
  y = partpair->GetQOutLCMS() - trackpair->GetQOutLCMS();
  if (trackpair->GetQOutLCMS() < 0.0) y = -y;
  x = partpair->GetKt();
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

void AliHBTQSideResolVsKtFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
{
  //returns Qsidesim-Qsiderec for y
  //returns Kt for x
  y = partpair->GetQSideLCMS() - trackpair->GetQSideLCMS();
  if (trackpair->GetQSideLCMS() < 0.0) y = -y;
  x = partpair->GetKt();
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

void AliHBTQLongResolVsKtFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
{
  //returns Qlongsim-Qlongrec for y
  //returns Kt for x
  y = partpair->GetQLongLCMS() - trackpair->GetQLongLCMS();
  if (trackpair->GetQLongLCMS() < 0.0) y = -y;
  x = partpair->GetKt();
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

void AliHBTQOutResolVsQOutFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
{
//returns values of the function
  x = partpair->GetQOutLCMS();
  y = x - trackpair->GetQOutLCMS();
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

void AliHBTQSideResolVsQSideFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
{
//returns values of the function
  x = partpair->GetQSideLCMS(); 
  y = x - trackpair->GetQSideLCMS();
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

void AliHBTQLongResolVsQLongFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
{
//returns values of the function
 x = partpair->GetQLongLCMS(); 
 y = x - trackpair->GetQLongLCMS();
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

void AliHBTPairThetaResolVsQInvFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
 {
  //returns Pair Theta sim - Pair Theta rec for y
  //returns Qinv for x
   Double_t partTheta = partpair->Particle1()->Theta() - partpair->Particle2()->Theta();
   Double_t trackTheta = trackpair->Particle1()->Theta() - trackpair->Particle2()->Theta();
   y = partTheta - trackTheta;
   x = partpair->GetQInv();
 }
/******************************************************************/
/******************************************************************/
/******************************************************************/

ClassImp( AliHBTPairThetaResolVsPairThetaFctn )

AliHBTPairThetaResolVsPairThetaFctn::
AliHBTPairThetaResolVsPairThetaFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
//ctor
 Rename("PairThetaResolVsPairTheta","Pair Theta Angle Resolution vs. Pair Theta ");
}
/******************************************************************/
void AliHBTPairThetaResolVsPairThetaFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
{
  //returns Pair Theta sim - Pair Theta rec for y
  //returns Pair Theta sim for x
  Double_t partTheta = partpair->Particle1()->Theta() - partpair->Particle2()->Theta();
  Double_t trackTheta = trackpair->Particle1()->Theta() - trackpair->Particle2()->Theta();
  y = partTheta - trackTheta;
  if (trackTheta < 0.0) y = -y;
  x = trackTheta;
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

void AliHBTPairPhiResolVsQInvFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
 {
  //returns Pair Phi sim - Pair Phi rec for y
  //returns QInv sim for x
  Double_t partPhi = partpair->Particle1()->Phi() - partpair->Particle2()->Phi();
  Double_t trackPhi = trackpair->Particle1()->Phi() - trackpair->Particle2()->Phi();
  y = partPhi - trackPhi;
  x = partpair->GetQInv();
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

void AliHBTPairThetaResolVsKtFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
{
  //returns Pair Theta sim - Pair Theta rec for y
  //returns Kt sim for x
  Double_t partTheta = partpair->Particle1()->Theta() - partpair->Particle2()->Theta();
  Double_t trackTheta = trackpair->Particle1()->Theta() - trackpair->Particle2()->Theta();
  y = partTheta - trackTheta;
  x = partpair->GetKt();
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

void AliHBTPairPhiResolVsKtFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
{
  //returns Pair Phi sim - Pair Phi rec for y
  //returns Kt sim for x
  Double_t partPhi = partpair->Particle1()->Phi() - partpair->Particle2()->Phi();
  Double_t trackPhi = trackpair->Particle1()->Phi() - trackpair->Particle2()->Phi();
  y = partPhi - trackPhi;
  x = partpair->GetKt();
}

/******************************************************************/
/******************************************************************/
/******************************************************************/


ClassImp( AliHBTPairPhiResolVsPairPhiFctn )

AliHBTPairPhiResolVsPairPhiFctn::
AliHBTPairPhiResolVsPairPhiFctn(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                           Int_t nYbins, Double_t maxYval, Double_t minYval):
                           AliHBTTwoPairFctn2D(nXbins,maxXval,minXval,nYbins,maxYval,minYval)
{
//ctor
 Rename("PairPhiResolVsPairPhi","Pair Phi Angle Resolution vs. Pair Phi ");
}
/******************************************************************/

void AliHBTPairPhiResolVsPairPhiFctn::GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const
{
  //returns Pair Phi sim - Pair Phi rec for y
  //returns Pair Phi sim for x
  Double_t partPhi = partpair->Particle1()->Phi() - partpair->Particle2()->Phi();
  Double_t trackPhi = trackpair->Particle1()->Phi() - trackpair->Particle2()->Phi();
  y = partPhi - trackPhi;
  x = trackPhi;
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
