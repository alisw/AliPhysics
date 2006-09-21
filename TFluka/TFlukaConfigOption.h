#ifndef TFLUKACONFIGOPTION
#define TFLUKACONFIGOPTION

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
// Class to store FLUKA and VMC configuration options:                       // 
// Cuts, Processes, User Scoring                                             // 
//                                                                           //
//                                                                           //
// Author: andreas.morsch@cern.ch                                            // 
//                                                                           //  
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
//


//

#include <TObject.h>

typedef enum {kDCAY, kPAIR, kCOMP, kPHOT, kPFIS, kDRAY, kANNI, kBREM,
	      kMUNU, kCKOV, kHADR, kLOSS, kMULS, kRAYL, kSTRA} FlukaProcessOption_t;
typedef enum {kCUTGAM, kCUTELE, kCUTNEU, kCUTHAD, kCUTMUO, kBCUTE, kBCUTM, kDCUTE, kDCUTM, kPPCUTM, kTOFMAX}  FlukaCutOption_t;
typedef enum {kPRIMION, kPRIMIOE}  FlukaModelOption_t;

class TFlukaMCGeometry;
class TGeoMaterial;

class TFlukaConfigOption : public TObject
{
public:
    // Constructors
    TFlukaConfigOption();
    TFlukaConfigOption(Int_t imed);
    // Getters
    Double_t Cut(FlukaCutOption_t i)              const {return fCutValue[i];}
    Int_t    Flag(FlukaProcessOption_t i)         const {return fProcessFlag[i];}
    Double_t  ModelParameter(FlukaModelOption_t i) const {return fModelParameter[i];}
    Int_t    Medium()                             const {return fMedium;}    
    // Setters
    void     SetCut(const char* flagname, Double_t val);
    void     SetModelParameter(const char* flagname, Double_t val);
    void     SetProcess(const char* flagname, Int_t flagValue);
    void     SetMedium(Int_t imed)       {fMedium = imed;}
    //
    void     WriteFlukaInputCards();
    void     ProcessDCAY();
    void     ProcessPAIR();
    void     ProcessBREM();
    void     ProcessCOMP();
    void     ProcessPHOT();
    void     ProcessANNI();
    void     ProcessPFIS();
    void     ProcessMUNU();
    void     ProcessRAYL();
    void     ProcessCKOV();
    void     ProcessHADR();
    void     ProcessMULS();
    void     ProcessLOSS();
    void     ProcessSTRA();
    
    
    void     ProcessCUTGAM();
    void     ProcessCUTELE();
    void     ProcessCUTNEU();
    void     ProcessCUTHAD();
    void     ProcessCUTMUO();
    void     ProcessTOFMAX();
    void     ProcessSensitiveMedium();
    
    //
    static void SetStaticInfo(FILE* file, Float_t matMin, Float_t matMax, TFlukaMCGeometry* geom)
	{fgFile = file; fgMatMin = matMin; fgMatMax = matMax; fgGeom = geom;}
    static Double_t DefaultCut(FlukaCutOption_t i) {return fgDCutValue[i];}
    static Int_t    DefaultProcessFlag(FlukaProcessOption_t i) {return fgDProcessFlag[i];}
 
    
 protected:
    Double_t fCutValue[11];            // User cut
    Int_t    fProcessFlag[15];         // User flag assigned to processes
    Double_t fModelParameter[15];      // User model parameter
    Int_t    fMedium;                  // Material assigned to user settings
    Float_t  fCMatMin;                 // Minimum material number used for current card 
    Float_t  fCMatMax;                 // Maximum material number used for current card
    TGeoMaterial* fCMaterial;          // Current material
    
    // static
    static Double_t  fgDCutValue[11];     // User default cut
    static Int_t     fgDProcessFlag[15];  // User default flag assigned to processes
    static Float_t   fgMatMin;            // Minimum material number 
    static Float_t   fgMatMax;            // Maximum meterial number
    static FILE*     fgFile;              // Output file
    static TFlukaMCGeometry* fgGeom;      // Pointer to geometry     

 private:
    // Copy constructor and operator= declared but not implemented (-Weff++ flag)
    TFlukaConfigOption(const TFlukaConfigOption&);
    TFlukaConfigOption& operator=(const TFlukaConfigOption&);
    
    ClassDef(TFlukaConfigOption, 1)    // Fluka Configuration Option
};
	
#endif
	
