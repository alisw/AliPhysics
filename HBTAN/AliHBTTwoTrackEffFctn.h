#ifndef AliHBTTwoTrackEffFctn_H
#define AliHBTTwoTrackEffFctn_H
//____________________________________________________________________
//////////////////////////////////////////////////////////////////////
//                                                                  //
//  class AliHBTTwoTrackEffFctn                                     //
//                                                                  //
//  classes for calculating two track efficiency of the tracking    //
//  binning is done using value of simulated pair montum difference // 
//  pair must be recontructed, that is why we need both pairs       //
//  (simulated and recontructed), thus functions are "two pair"     //
//  Piotr.Skowronski@cern.ch                                        //
//                                                                  //
//////////////////////////////////////////////////////////////////////
#include "AliHBTPair.h"
#include "AliHBTFunction.h"

class AliHBTTwoTrackEffFctn: public AliHBTOnePairFctn1D, public AliHBTCorrelFunction
 {
  public:
    AliHBTTwoTrackEffFctn();
    AliHBTTwoTrackEffFctn(Int_t nbins, Double_t maxval, Double_t minval);
    virtual ~AliHBTTwoTrackEffFctn(){}
    TH1* GetResult();
  protected:
    Double_t GetValue(AliHBTPair* pair) const {return pair->GetDeltaP();}
  private:
    ClassDef(AliHBTTwoTrackEffFctn,2)
 };
/******************************************************************/

class AliHBTTwoTrackEffFctnPxPyPz: public AliHBTOnePairFctn3D, public AliHBTCorrelFunction
 {
  public:
    AliHBTTwoTrackEffFctnPxPyPz(Int_t nXbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0,
                                Int_t nYbins = 100, Double_t maxYval = 0.15, Double_t minYval = 0.0,
	            Int_t nZbins = 100, Double_t maxZval = 0.15, Double_t minZval = 0.0);
    virtual ~AliHBTTwoTrackEffFctnPxPyPz(){}
    TH1* GetResult();
  protected:
    void GetValues(AliHBTPair* pair,Double_t& x, Double_t& y,Double_t& z) const;
  private:
    ClassDef(AliHBTTwoTrackEffFctnPxPyPz,2)
 };
/******************************************************************/

class AliHBTTwoTrackEffFctnPtThetaPhi: public AliHBTOnePairFctn3D, public AliHBTCorrelFunction
 {
  public:
    AliHBTTwoTrackEffFctnPtThetaPhi(Int_t nXbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0,
                                    Int_t nYbins = 100, Double_t maxYval = 0.3, Double_t minYval = 0.0,
	                Int_t nZbins = 100, Double_t maxZval = 0.3, Double_t minZval = 0.0);
    virtual ~AliHBTTwoTrackEffFctnPtThetaPhi(){}
    TH1* GetResult();
  protected:
    void GetValues(AliHBTPair* pair,Double_t& x, Double_t& y,Double_t& z) const;
  private:
    ClassDef(AliHBTTwoTrackEffFctnPtThetaPhi,1)
 };

#endif
