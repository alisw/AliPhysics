#ifndef ALIHBTCORRECTQINVCORRELFCTN_H
#define ALIHBTCORRECTQINVCORRELFCTN_H
//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTCorrectQInvCorrelFctn                       //
//                                                   //
// Class for calculating Q Invariant correlation     //
// taking to the account resolution of the           //
// detector and coulomb effects.                     //
//                                                   //
///////////////////////////////////////////////////////

#include "AliHBTFunction.h"

class AliHBTCorrectQInvCorrelFctn: public AliHBTOnePairFctn1D
{
  public:
    AliHBTCorrectQInvCorrelFctn(const char* name = "qinvcorrectedCF", 
                                const char* title= "Corrected Q_{inv} Correlation Fonction");

    AliHBTCorrectQInvCorrelFctn(const char* name, const char* title,
	            Int_t nbins, Float_t maxXval, Float_t minXval);

    AliHBTCorrectQInvCorrelFctn(TH1D* measqinv, 
                                const char* name = "qinvcorrectedCF", 
                                const char* title= "Corrected Q_{inv} Correlation Fonction");
    AliHBTCorrectQInvCorrelFctn(const AliHBTCorrectQInvCorrelFctn& in);
    
    virtual ~AliHBTCorrectQInvCorrelFctn();
    
    void     SetInitialValues(Double_t lambda, Double_t r);
    void     Init();
    void     ProcessSameEventParticles(AliHBTPair* pair);//process particles from same event (real pair)
    void     ProcessDiffEventParticles(AliHBTPair* pair);//process particles coming from different events (mixed pairs)
    void     SetMeasuredHistogram(TH1D* meas){fMeasCorrelFctn = meas;}
    TH1*     GetResult();//returns the result histogram
    Double_t GetRadius()const{ return TMath::Sqrt(fR2);}//returns assumed radius
    Double_t GetLambda()const{ return fLambda;}//retutrns assumed intercept parameter
    void     SetRadiusConvergenceTreshold(Double_t ct){fRConvergenceTreshold=ct;}//if fitted and assumed R us different less then that number con
    void     SetLambdaConvergenceTreshold(Double_t ct){fLambdaConvergenceTreshold=ct;}
    Bool_t   IsConverged();
    void     Fit();
    Double_t GetFittedRadius();
    Double_t GetFittedLambda();
    void     WriteAll();
    void     SetMeasNum(TH1D* measnum){fMeasNumer = measnum;}
    void     SetMeasDen(TH1D* h){fMeasDenom = h;}
    void     MakeMeasCF();
    
  protected:
    virtual void BuildHistos(Int_t nbins, Float_t max, Float_t min);
    Double_t GetCoulombCorrection(AliHBTPair* /*pair*/){return 1.0;}
    Double_t GetValue(AliHBTPair * pair) const {return pair->GetQInv();}
    void Smear(AliHBTPair* pair,AliHBTPair& smeared);
    void Smear(AliHBTParticle* part, AliHBTParticle* smeared);
    Double_t GetModelValue(Double_t qinv);

    //Our ideal numerator 
    TH1D* fMeasCorrelFctn; //Measured correlation function
    TH1D* fMeasNumer;//Measured numerator correlation function
    TH1D* fMeasDenom;  //Measured denominator correlation function
    
    TH1D* fSmearedNumer; //! Numerator of smeard q
    TH1D* fSmearedDenom; //! Denominator of smeard q
    
    //Parameters of Pt RMS
    //linear dependence dPt/Pt from Pt itself 
    Float_t fDPtOverPtRMS; //RMS of dPt/Pt 
    
    //We assume that RMS of Theta and Phisangle depends on Pt Like A+B*(Pt)^Alpha
    //Idea copied from Star HBT Maker (Fabrice Retiere)
    //Parameters comes from Monte Carlo Resolution Analysis

    Float_t fThetaA; //"A" parameter of theta RMS dependence
    Float_t fThetaB; //"B" parameter of theta RMS dependence
    Float_t fThetaAlpha; //"Alpha" parameter (power) of theta RMS dependence

    Float_t fPhiA;//"A" parameter of phi RMS dependence
    Float_t fPhiB;//"B" parameter of phi RMS dependence
    Float_t fPhiAlpha;//"Alpha" parameter (power) of phi RMS dependence
    
    Double_t fR2;//square of radius
    Double_t fLambda;//Interception parameter

    Double_t fFittedR;//fitted radius
    Double_t fFittedLambda;//fitted Interception parameter
        
    Float_t  fRConvergenceTreshold;//fRConvergenceTreshold
    Float_t  fLambdaConvergenceTreshold;//fLambdaConvergenceTreshold
    
  private:
    ClassDef(AliHBTCorrectQInvCorrelFctn,1)
};

inline Double_t AliHBTCorrectQInvCorrelFctn::GetModelValue(Double_t qinv)
{
  //factor 0.038936366329 conected with units change GeV<->SI
  return 1.0 + fLambda*TMath::Exp(-fR2*qinv*qinv/0.038936366329);
}

#endif
