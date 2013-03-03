#ifndef ALI_TPC_EX_BB_SHAPE_H
#define ALI_TPC_EX_BB_SHAPE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// _________________________________________________________________
//
// Begin_Html
//   <h2>AliExBBShape class </h2>                                                    
//   The class calculates the space point distortions due to the B field    
//   shape imperfections using a second order technique based on integrals  
//   over Bz (e.g. int By/Bz) obtained via the AliMagF class. The essential 
//   input for this class is the magnetic field maps which can be set via the function SetBField.  
//   <p>
//   The class allows "effective Omega Tau" corrections.                    
//   End_Html
//   Begin_Macro(source)
//   {
//   gROOT->SetStyle("Plain"); gStyle->SetPalette(1);
//   TCanvas *c2 = new TCanvas("cAliTPCExBBShape","cAliTPCExBBShape",500,300); 
//   AliTPCExBBShape exb;                                                 
//   AliMagF mag("mag","mag");        // 0.5 Tesla (solenoid)
//   exb.SetBField(&mag);             // use Bfield from AliMagF        
//   exb.SetOmegaTauT1T2(-0.32,1.,1.); // values ideally from OCDB     
//   exb.CreateHistoDRPhiinZR(0,100,100)->Draw("surf2"); 
//   return c2;
//   } 
// End_Macro
// Begin_Html
//   <p>
//   Date: 27/04/2010 <br>                                                      
//   Authors: Magnus Mager, Jim Thomas, Stefan Rossegger                    
// End_Html 
// _________________________________________________________________

#include "AliTPCCorrection.h"

class AliMagF;

class AliTPCExBBShape : public AliTPCCorrection {
public:
  AliTPCExBBShape();
  virtual ~AliTPCExBBShape();

 // initialization and update functions
  virtual void Init();
  virtual void Update(const TTimeStamp &timeStamp);

  // common setters and getters for ExB
  virtual void SetOmegaTauT1T2(Float_t omegaTau,Float_t t1,Float_t t2) {
    fT1=t1; fT2=t2;
    const Float_t wt1=t1*omegaTau;    fC1=wt1/(1.+wt1*wt1);
    const Float_t wt2=t2*omegaTau;    fC2=wt2*wt2/(1.+wt2*wt2);
  };
  void SetC1C2(Float_t c1,Float_t c2) {fC1=c1;fC2=c2;} // CAUTION: USE WITH CARE
  Float_t GetC1() const {return fC1;}
  Float_t GetC2() const {return fC2;}

  // setters and getters for the magentic field map
  void SetBField(const AliMagF *bField) {fBField=(AliMagF*)bField;}
  AliMagF* GetBField() const {return fBField;}

  virtual void GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]);
  void GetBxAndByOverBz(const Float_t x[],const Short_t roc,Float_t BxByOverBz[]);

  virtual void Print(Option_t* option="") const;
  static Double_t GetBFieldXYZ(Double_t gx, Double_t gy, Double_t gz, Int_t axisType);

private:
  Float_t fC1; // coefficient C1          (compare Jim Thomas's notes for definitions)
  Float_t fC2; // coefficient C2          (compare Jim Thomas's notes for definitions)

  AliMagF *fBField;       // pointer to magnetic field

  AliTPCExBBShape & operator =(const AliTPCExBBShape &); // dummy assignment operator
  AliTPCExBBShape(const AliTPCExBBShape&); //dummy copy contructor

  ClassDef(AliTPCExBBShape,2);
};

#endif
