#ifndef ALIHBTWEIGHTNONID3DCORRFCTN_H
#define ALIHBTWEIGHTNONID3DCORRFCTN_H

///////////////////////////////////////////////////////
//                                                   //
// AliHBTWeightNonId3DCorrFctn.h                     //
//                                                   //
// Class for calculating 3D non-id correlation       //
// functions using method of weights                 //
//                                                   //
///////////////////////////////////////////////////////

#include "AliHBTFunction.h"


class AliHBTWeights;

class AliHBTWeightNonId3DCorrFctn: public AliHBTTwoPairFctn1D
{
  public:
   AliHBTWeightNonId3DCorrFctn(const char* name = "nonid3DCF", 
			       const char* title= "3D Non-Id Correlation Function");

   AliHBTWeightNonId3DCorrFctn(const char* name, const char* title,
			       Int_t nbinsX, Float_t maxXval, Float_t minXval);
   AliHBTWeightNonId3DCorrFctn(const AliHBTWeightNonId3DCorrFctn& in);

   virtual ~AliHBTWeightNonId3DCorrFctn();

   void Init(); // InitFunction();
   void ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
   void ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);

   void WriteFunction();
   
   TH1*     GetResult();
   
 protected:

   Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair) {return trackpair->GetQInv()-partpair->GetQInv();}
   void BuildHistos(Int_t nbins, Float_t max, Float_t min);
   
   TH1D* fWeightNumOutP;
   TH1D* fWeightDenOutP;
   TH1D* fWeightRatOutP;
   TH1D* fWeightNumOutN;
   TH1D* fWeightDenOutN;
   TH1D* fWeightRatOutN;
   TH1D* fWeightRatOut;
   TH1D* fWeightRatOutNOverP;

   TH1D* fWeightNumSideP;
   TH1D* fWeightDenSideP;
   TH1D* fWeightRatSideP;
   TH1D* fWeightNumSideN;
   TH1D* fWeightDenSideN;
   TH1D* fWeightRatSideN;
   TH1D* fWeightRatSide;
   TH1D* fWeightRatSideNOverP;

   TH1D* fWeightNumLongP;
   TH1D* fWeightDenLongP;
   TH1D* fWeightRatLongP;
   TH1D* fWeightNumLongN;
   TH1D* fWeightDenLongN;
   TH1D* fWeightRatLongN;
   TH1D* fWeightRatLong;
   TH1D* fWeightRatLongNOverP;
   
    
  private:
  
    ClassDef(AliHBTWeightNonId3DCorrFctn,1)
};

#endif
