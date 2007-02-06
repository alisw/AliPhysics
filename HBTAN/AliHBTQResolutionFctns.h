#ifndef ALIHBTQRESOLUTIONFCTNS_H
#define ALIHBTQRESOLUTIONFCTNS_H
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

#include <TMath.h>

class AliHBTKtResolVsQInvFctn;  //Kt Res   Vs   QInvLCMS 

class AliHBTQOutResolVsQInvFctn;  //QOutLCMS  Res   Vs   QInvLCMS 
class AliHBTQSideResolVsQInvFctn; //QSideLCMS Res   Vs   QInvLCMS 
class AliHBTQLongResolVsQInvFctn; //QLongLCMS Res   Vs   QInvLCMS 
class AliHBTQInvResolVsQInvFctn;  //QInvLCMS Res   Vs   QInvLCMS 

class AliHBTPairThetaResolVsQInvFctn;//PairTheta Res   Vs   QInvLCMS
class AliHBTPairPhiResolVsQInvFctn;  //PairPhi   Res   Vs   QInvLCMS

class AliHBTQInvResolVsKtFctn;    //QInvLCMS  Res   Vs   Kt
class AliHBTQOutResolVsKtFctn;    //QOutLCMS  Res   Vs   Kt
class AliHBTQSideResolVsKtFctn;   //QSideLCMS Res   Vs   Kt
class AliHBTQLongResolVsKtFctn;   //QLongLCMS Res   Vs   Kt

class AliHBTPairThetaResolVsKtFctn;   //PairTheta Res   Vs   Kt
class AliHBTPairPhiResolVsKtFctn;   //PairPhi Res   Vs   Kt


class AliHBTQOutResolVsQOutFctn;  //QOutLCMS  Res   Vs   QOut
class AliHBTQSideResolVsQSideFctn;//QSideLCMS Res   Vs   QSide
class AliHBTQLongResolVsQLongFctn;//QLongLCMS Res   Vs   QLong


 
#include "AliHBTFunction.h"
/***********************************************************************/
/***********************************************************************/
class AliHBTKtResolVsQInvFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTKtResolVsQInvFctn(Int_t nXbins = 200, Double_t maxXval = 0.2, Double_t minXval = 0.0, 
                           Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   
   virtual ~AliHBTKtResolVsQInvFctn(){}
   
   TH1* GetResult(){return this->fNumerator;}  
   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   ClassDef(AliHBTKtResolVsQInvFctn,1)
 };


/***********************************************************************/
/***********************************************************************/
class AliHBTQInvResolVsQInvFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTQInvResolVsQInvFctn(Int_t nXbins = 200, Double_t maxXval = 0.2, Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   
   virtual ~AliHBTQInvResolVsQInvFctn(){}
   
   TH1* GetResult(){return this->fNumerator;}  
   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   ClassDef(AliHBTQInvResolVsQInvFctn,1)
 };


/***********************************************************************/
/***********************************************************************/
class AliHBTQOutResolVsQInvFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTQOutResolVsQInvFctn(Int_t nXbins = 200, Double_t maxXval = 0.2, Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   
   virtual ~AliHBTQOutResolVsQInvFctn(){}
   
   TH1* GetResult(){return this->fNumerator;}  
   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   ClassDef(AliHBTQOutResolVsQInvFctn,1)
 };

/***********************************************************************/
/***********************************************************************/
class AliHBTQSideResolVsQInvFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTQSideResolVsQInvFctn(Int_t nXbins = 200, Double_t maxXval = 0.2, Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .05, Double_t minYval =-0.05);
   virtual ~AliHBTQSideResolVsQInvFctn(){}

   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair,  Double_t& x, Double_t& y) const;
   TH1* GetResult(){return this->fNumerator;} 
   ClassDef(AliHBTQSideResolVsQInvFctn,1)
 };

/***********************************************************************/
/***********************************************************************/
class AliHBTQLongResolVsQInvFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTQLongResolVsQInvFctn(Int_t nXbins = 200, Double_t maxXval = 0.2, Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .05, Double_t minYval =-0.05);
   virtual ~AliHBTQLongResolVsQInvFctn(){}

   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   TH1* GetResult(){return this->fNumerator;} 
   ClassDef(AliHBTQLongResolVsQInvFctn,1)
 };

/***********************************************************************/
/***********************************************************************/
class AliHBTQInvResolVsKtFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTQInvResolVsKtFctn(Int_t nXbins = 200, Double_t maxXval = 1., Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .05, Double_t minYval =-0.05);
   virtual ~AliHBTQInvResolVsKtFctn(){};

   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   TH1* GetResult(){return this->fNumerator;} 
   ClassDef(AliHBTQInvResolVsKtFctn,1)
 };
/***********************************************************************/
/***********************************************************************/
class AliHBTQOutResolVsKtFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTQOutResolVsKtFctn(Int_t nXbins = 200, Double_t maxXval = 1., Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   virtual ~AliHBTQOutResolVsKtFctn(){}
   TH1* GetResult(){return this->GetNumerator();}
   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   ClassDef(AliHBTQOutResolVsKtFctn,1)
 };
/***********************************************************************/
/***********************************************************************/
class AliHBTQSideResolVsKtFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTQSideResolVsKtFctn(Int_t nXbins = 200, Double_t maxXval = 1., Double_t minXval = 0.0, 
                            Int_t nYbins = 500, Double_t maxYval = .05, Double_t minYval =-0.05);
   virtual ~AliHBTQSideResolVsKtFctn(){}
   TH1* GetResult(){return this->GetNumerator();}
   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   ClassDef(AliHBTQSideResolVsKtFctn,1)
 };
/***********************************************************************/
/***********************************************************************/
class AliHBTQLongResolVsKtFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTQLongResolVsKtFctn(Int_t nXbins = 200, Double_t maxXval = 1., Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .05, Double_t minYval =-0.05);
   virtual ~AliHBTQLongResolVsKtFctn(){}

   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   TH1* GetResult(){return this->GetNumerator();}
   ClassDef(AliHBTQLongResolVsKtFctn,1)
 };
/***********************************************************************/
/***********************************************************************/
class AliHBTQOutResolVsQOutFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTQOutResolVsQOutFctn(Int_t nXbins = 200, Double_t maxXval = 0.2, Double_t minXval = -0.2, 
                             Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   virtual ~AliHBTQOutResolVsQOutFctn(){}

   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   TH1* GetResult(){return this->fNumerator;}  
   ClassDef(AliHBTQOutResolVsQOutFctn,1)
 };

/***********************************************************************/
/***********************************************************************/

class AliHBTQSideResolVsQSideFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTQSideResolVsQSideFctn(Int_t nXbins = 200, Double_t maxXval = 0.2, Double_t minXval = -0.2, 
                             Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   virtual ~AliHBTQSideResolVsQSideFctn(){}

   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   TH1* GetResult(){return this->fNumerator;}  
   ClassDef(AliHBTQSideResolVsQSideFctn,1)
 };

/***********************************************************************/
/***********************************************************************/

class AliHBTQLongResolVsQLongFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTQLongResolVsQLongFctn(Int_t nXbins = 200, Double_t maxXval = 0.2, Double_t minXval = -0.2,
                             Int_t nYbins = 500, Double_t maxYval = .05, Double_t minYval =-0.05);
   virtual ~AliHBTQLongResolVsQLongFctn(){}

   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   TH1* GetResult(){return this->fNumerator;}  
   ClassDef(AliHBTQLongResolVsQLongFctn,1)
 };
/***********************************************************************/
/***********************************************************************/

class AliHBTPairThetaResolVsQInvFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTPairThetaResolVsQInvFctn(Int_t nXbins = 200, Double_t maxXval = 0.2, Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   
   virtual ~AliHBTPairThetaResolVsQInvFctn(){}
   
   TH1* GetResult(){return this->fNumerator;}  
   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   ClassDef(AliHBTPairThetaResolVsQInvFctn,1)
 };

/***********************************************************************/
/***********************************************************************/
class AliHBTPairThetaResolVsPairThetaFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTPairThetaResolVsPairThetaFctn(Int_t nXbins = 200, Double_t maxXval = TMath::PiOver2(), Double_t minXval = -TMath::PiOver2(), 
                             Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   
   virtual ~AliHBTPairThetaResolVsPairThetaFctn(){}
   
   TH1* GetResult(){return this->fNumerator;}  
   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   ClassDef(AliHBTPairThetaResolVsPairThetaFctn,1)
 };


/***********************************************************************/
/***********************************************************************/
class AliHBTPairPhiResolVsQInvFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTPairPhiResolVsQInvFctn(Int_t nXbins = 200, Double_t maxXval = 0.2, Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   
   virtual ~AliHBTPairPhiResolVsQInvFctn(){}
   
   TH1* GetResult(){return this->fNumerator;}  
   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   ClassDef(AliHBTPairPhiResolVsQInvFctn,1)
 };

/***********************************************************************/
/***********************************************************************/
class AliHBTPairThetaResolVsKtFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTPairThetaResolVsKtFctn(Int_t nXbins = 200, Double_t maxXval = 0.2, Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   
   virtual ~AliHBTPairThetaResolVsKtFctn(){}
   
   TH1* GetResult(){return this->fNumerator;}  
   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   ClassDef(AliHBTPairThetaResolVsKtFctn,1)
 };

/***********************************************************************/
/***********************************************************************/
class AliHBTPairPhiResolVsKtFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTPairPhiResolVsKtFctn(Int_t nXbins = 200, Double_t maxXval = 0.2, Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   
   virtual ~AliHBTPairPhiResolVsKtFctn(){}
   
   TH1* GetResult(){return this->fNumerator;}  
   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   ClassDef(AliHBTPairPhiResolVsKtFctn,1)
 };

/***********************************************************************/
/***********************************************************************/
class AliHBTPairPhiResolVsPairPhiFctn: public AliHBTTwoPairFctn2D
 {
  public: 
   AliHBTPairPhiResolVsPairPhiFctn(Int_t nXbins = 200, Double_t maxXval = TMath::TwoPi(), Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   
   virtual ~AliHBTPairPhiResolVsPairPhiFctn(){}
   
   TH1* GetResult(){return this->fNumerator;}  
   void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const;
   ClassDef(AliHBTPairPhiResolVsPairPhiFctn,1)
 };


#endif
