#ifndef ALIHBTNONID3DCORRFCTN_H
#define ALIHBTNONID3DCORRFCTN_H

///////////////////////////////////////////////////////
//                                                   //
// AliHBTNonId3DCorrFctn.h                           //
//                                                   //
// Class for calculating 3D non-id correlation       //
// functions                                         //
//                                                   //
///////////////////////////////////////////////////////

#include "AliHBTFunction.h"


class AliHBTNonId3DCorrFctn: public AliHBTOnePairFctn1D
{
  public:
   AliHBTNonId3DCorrFctn(const char* name = "nonid3DCF", 
                         const char* title= "3D Non-Id Correlation Function");

   AliHBTNonId3DCorrFctn(const char* name, const char* title,
			 Int_t nbins, Float_t maxXval, Float_t minXval);
   AliHBTNonId3DCorrFctn(const AliHBTNonId3DCorrFctn& in);

   virtual ~AliHBTNonId3DCorrFctn();

   void Init();
   void ProcessSameEventParticles(AliHBTPair* pair);
   void ProcessDiffEventParticles(AliHBTPair* pair);

   void WriteFunction();
   
      TH1*     GetResult();
   
 protected:

   Double_t GetValue(AliHBTPair* pair) {return pair->GetQInv();}
   void BuildHistos(Int_t nbins, Float_t max, Float_t min);
   
   TH1D* fNumOutP;
   TH1D* fDenOutP;
   TH1D* fRatOutP;
   TH1D* fNumOutN;
   TH1D* fDenOutN;
   TH1D* fRatOutN;
   TH1D* fRatOut;
   TH1D* fRatOutNOverP;

   TH1D* fNumSideP;
   TH1D* fDenSideP;
   TH1D* fRatSideP;
   TH1D* fNumSideN;
   TH1D* fDenSideN;
   TH1D* fRatSideN;
   TH1D* fRatSide;
   TH1D* fRatSideNOverP;

   TH1D* fNumLongP;
   TH1D* fDenLongP;
   TH1D* fRatLongP;
   TH1D* fNumLongN;
   TH1D* fDenLongN;
   TH1D* fRatLongN;
   TH1D* fRatLong;
   TH1D* fRatLongNOverP;
   
    
  private:
  
    ClassDef(AliHBTNonId3DCorrFctn,1)
};

#endif
