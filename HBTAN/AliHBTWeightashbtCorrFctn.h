#ifndef ALIHBTWEIGHTASHBTCORRFCTN_H
#define ALIHBTWEIGHTASHBTCORRFCTN_H

///////////////////////////////////////////////////////
//                                                   //
// AliHBTashbtCorrFctn.h                           //
//                                                   //
// Class for calculating 3D ashbt correlation       //
// functions                                         //
//                                                   //
///////////////////////////////////////////////////////

#include "AliHBTFunction.h"


class AliHBTWeightashbtCorrFctn: public AliHBTOnePairFctn1D
{
  public:
   AliHBTWeightashbtCorrFctn(const char* name = "asejdzbitiCF", 
                         const char* title= "asHBT Correlation Function");

   AliHBTWeightashbtCorrFctn(const char* name, const char* title,
			 Int_t nbins, Float_t maxXval, Float_t minXval);
   AliHBTWeightashbtCorrFctn(const AliHBTWeightashbtCorrFctn& in);

   virtual ~AliHBTWeightashbtCorrFctn();

   void Init();
   void ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
   void ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);

   void WriteFunction();
   
   TH1*     GetResult();
   
 protected:

   Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair) 
     {
	return trackpair->GetQInv()-partpair->GetQInv();
     }
   
   void BuildHistos(Int_t nbins, Float_t max, Float_t min);
   
   TH1D* fWeightNumOut1;
   TH1D* fWeightNumOut2;
   TH1D* fWeightNumOut3;
   TH1D* fWeightNumOut4;
   TH1D* fWeightNumOut5;
   TH1D* fWeightNumOut6;
   TH1D* fWeightNumOut7;
   TH1D* fWeightNumOut8;

   TH1D* fWeightDenOut1;
   TH1D* fWeightDenOut2;
   TH1D* fWeightDenOut3;
   TH1D* fWeightDenOut4;
   TH1D* fWeightDenOut5;
   TH1D* fWeightDenOut6;
   TH1D* fWeightDenOut7;
   TH1D* fWeightDenOut8;
   
   TH1D* fWeightRatOut1;
   TH1D* fWeightRatOut2;
   TH1D* fWeightRatOut3;
   TH1D* fWeightRatOut4;
   TH1D* fWeightRatOut5;
   TH1D* fWeightRatOut6;
   TH1D* fWeightRatOut7;
   TH1D* fWeightRatOut8;
   

   TH1D* fWeightNumSide1;
   TH1D* fWeightNumSide2;
   TH1D* fWeightNumSide3;
   TH1D* fWeightNumSide4;
   TH1D* fWeightNumSide5;
   TH1D* fWeightNumSide6;
   TH1D* fWeightNumSide7;
   TH1D* fWeightNumSide8;

   TH1D* fWeightDenSide1;
   TH1D* fWeightDenSide2;
   TH1D* fWeightDenSide3;
   TH1D* fWeightDenSide4;
   TH1D* fWeightDenSide5;
   TH1D* fWeightDenSide6;
   TH1D* fWeightDenSide7;
   TH1D* fWeightDenSide8;
   
   TH1D* fWeightRatSide1;
   TH1D* fWeightRatSide2;
   TH1D* fWeightRatSide3;
   TH1D* fWeightRatSide4;
   TH1D* fWeightRatSide5;
   TH1D* fWeightRatSide6;
   TH1D* fWeightRatSide7;
   TH1D* fWeightRatSide8;
   
   TH1D* fWeightNumLong1;
   TH1D* fWeightNumLong2;
   TH1D* fWeightNumLong3;
   TH1D* fWeightNumLong4;
   TH1D* fWeightNumLong5;
   TH1D* fWeightNumLong6;
   TH1D* fWeightNumLong7;
   TH1D* fWeightNumLong8;

   TH1D* fWeightDenLong1;
   TH1D* fWeightDenLong2;
   TH1D* fWeightDenLong3;
   TH1D* fWeightDenLong4;
   TH1D* fWeightDenLong5;
   TH1D* fWeightDenLong6;
   TH1D* fWeightDenLong7;
   TH1D* fWeightDenLong8;
   
   TH1D* fWeightRatLong1;
   TH1D* fWeightRatLong2;
   TH1D* fWeightRatLong3;
   TH1D* fWeightRatLong4;
   TH1D* fWeightRatLong5;
   TH1D* fWeightRatLong6;
   TH1D* fWeightRatLong7;
   TH1D* fWeightRatLong8;
   

  private:
  
    ClassDef(AliHBTWeightashbtCorrFctn,1)
};

#endif
