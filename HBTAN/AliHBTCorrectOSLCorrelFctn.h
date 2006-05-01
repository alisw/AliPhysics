#ifndef ALIHBTCORRECTOSLCORRELFCTN_H
#define ALIHBTCORRECTOSLCORRELFCTN_H
//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTCorrectQ3DCorrelFctn                        //
//                                                   //
// Class for calculating Q Invariant correlation     //
// taking to the account resolution of the           //
// detector and coulomb effects.                     //
//                                                   //
//          N[meas]   N[ideal]/D[ideal]
//  C(Q) =  ------- * -----------------
//          D[meas]   N[smear]/D[smear]
//
// if smeared is eqal to the measured than  we get ideal.
///////////////////////////////////////////////////////

#include "AliHBTCorrectQInvCorrelFctn.h"


class AliHBTCorrectOSLCorrelFctn: public AliHBTOnePairFctn3D, public AliHBTCorrectedCorrelFctn
{
  public:
   AliHBTCorrectOSLCorrelFctn(const char* name = "qoslcorrectedCF", 
                              const char* title= "Corrected Q_{out}-Q_{side}-Q_{long} Correlation Fonction");

  AliHBTCorrectOSLCorrelFctn(const Char_t *name, const Char_t *title,
	         Int_t nXbins, Double_t maxXval, Double_t minXval, 
	         Int_t nYbins, Double_t maxYval, Double_t minYval, 
	         Int_t nZbins, Double_t maxZval, Double_t minZval);
	           
   AliHBTCorrectOSLCorrelFctn(const AliHBTCorrectOSLCorrelFctn& in);
   virtual ~AliHBTCorrectOSLCorrelFctn();

   void     ProcessSameEventParticles(AliHBTPair* pair);//process particles from same event (real pair)
   void     ProcessDiffEventParticles(AliHBTPair* pair);//process particles coming from different events (mixed pairs)

   void     SetInitialValues(Double_t lambda, Double_t rout, Double_t rside, Double_t rlong);
   void     Init();
   Int_t    WriteFunction();//overloaded 

   TH1*     GetResult();//returns the result histogram
   void     GetValues(AliHBTPair* pair, Double_t& x, Double_t& y, Double_t& z) const ;
   
   Double_t GetModelValue(Double_t qout, Double_t qside, Double_t qlong) const;
   
   
  protected:

    void BuildHistos(Int_t nxbins, Float_t xmax, Float_t xmin,
                     Int_t nybins, Float_t ymax, Float_t ymin,
	 Int_t nzbins, Float_t zmax, Float_t zmin);
    virtual void BuildHistos() {AliHBTFunction3D::BuildHistos();}

    TH3F* fMeasCorrelFctn; //!Measured correlation function
    
    TH3F* fSmearedNumer; //! Numerator of smeard q
    TH3F* fSmearedDenom; //! Denominator of smeard q
    TH3F* fMeasNumer;  //! Numerator of ideal q calculated on basis of model equation
    TH3F* fMeasDenom;  //! Denominator of ideal q calculated on basis of model equation

    Double_t fLambda;
    Double_t fROutSq;
    Double_t fRSideSq;
    Double_t fRLongSq;
    
  private:
  
    ClassDef(AliHBTCorrectOSLCorrelFctn,1)
};

inline Double_t AliHBTCorrectOSLCorrelFctn::GetModelValue(Double_t qout, Double_t qside, Double_t qlong) const
{
 //returns model value of the cf
  return 1.0 + fLambda*TMath::Exp(( fROutSq*qout*qout + fRSideSq*qside*qside + fRLongSq*qlong*qlong) / (-0.038936366329));
}

#endif
