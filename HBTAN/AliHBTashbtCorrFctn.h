#ifndef ALIHBTASHBTCORRFCTN_H
#define ALIHBTASHBTCORRFCTN_H

///////////////////////////////////////////////////////
//                                                   //
// AliHBTashbtCorrFctn.h                           //
//                                                   //
// Class for calculating 3D ashbt correlation       //
// functions                                         //
//                                                   //
///////////////////////////////////////////////////////

#include "AliHBTFunction.h"


class AliHBTashbtCorrFctn: public AliHBTOnePairFctn1D
{
  public:
   AliHBTashbtCorrFctn(const char* name = "asejdzbitiCF", 
                         const char* title= "asHBT Correlation Function");

   AliHBTashbtCorrFctn(const char* name, const char* title,
			 Int_t nbins, Float_t maxXval, Float_t minXval);
   AliHBTashbtCorrFctn(const AliHBTashbtCorrFctn& in);

   virtual ~AliHBTashbtCorrFctn();

   void Init();
   void ProcessSameEventParticles(AliHBTPair* pair);
   void ProcessDiffEventParticles(AliHBTPair* pair);

   void WriteFunction();
   
   TH1*     GetResult();
   
 protected:

   Double_t GetValue(AliHBTPair* pair) const {return pair->GetQInv();}
   void BuildHistos(Int_t nbins, Float_t max, Float_t min);
   
   TH1D* fNumOut1;
   TH1D* fNumOut2;
   TH1D* fNumOut3;
   TH1D* fNumOut4;
   TH1D* fNumOut5;
   TH1D* fNumOut6;
   TH1D* fNumOut7;
   TH1D* fNumOut8;

   TH1D* fDenOut1;
   TH1D* fDenOut2;
   TH1D* fDenOut3;
   TH1D* fDenOut4;
   TH1D* fDenOut5;
   TH1D* fDenOut6;
   TH1D* fDenOut7;
   TH1D* fDenOut8;
   
   TH1D* fRatOut1;
   TH1D* fRatOut2;
   TH1D* fRatOut3;
   TH1D* fRatOut4;
   TH1D* fRatOut5;
   TH1D* fRatOut6;
   TH1D* fRatOut7;
   TH1D* fRatOut8;
   

   TH1D* fNumSide1;
   TH1D* fNumSide2;
   TH1D* fNumSide3;
   TH1D* fNumSide4;
   TH1D* fNumSide5;
   TH1D* fNumSide6;
   TH1D* fNumSide7;
   TH1D* fNumSide8;

   TH1D* fDenSide1;
   TH1D* fDenSide2;
   TH1D* fDenSide3;
   TH1D* fDenSide4;
   TH1D* fDenSide5;
   TH1D* fDenSide6;
   TH1D* fDenSide7;
   TH1D* fDenSide8;
   
   TH1D* fRatSide1;
   TH1D* fRatSide2;
   TH1D* fRatSide3;
   TH1D* fRatSide4;
   TH1D* fRatSide5;
   TH1D* fRatSide6;
   TH1D* fRatSide7;
   TH1D* fRatSide8;
   
   TH1D* fNumLong1;
   TH1D* fNumLong2;
   TH1D* fNumLong3;
   TH1D* fNumLong4;
   TH1D* fNumLong5;
   TH1D* fNumLong6;
   TH1D* fNumLong7;
   TH1D* fNumLong8;

   TH1D* fDenLong1;
   TH1D* fDenLong2;
   TH1D* fDenLong3;
   TH1D* fDenLong4;
   TH1D* fDenLong5;
   TH1D* fDenLong6;
   TH1D* fDenLong7;
   TH1D* fDenLong8;
   
   TH1D* fRatLong1;
   TH1D* fRatLong2;
   TH1D* fRatLong3;
   TH1D* fRatLong4;
   TH1D* fRatLong5;
   TH1D* fRatLong6;
   TH1D* fRatLong7;
   TH1D* fRatLong8;
   

  private:
  
    ClassDef(AliHBTashbtCorrFctn,1)
};

#endif
