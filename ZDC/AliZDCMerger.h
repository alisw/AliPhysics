#ifndef ALIZDCMERGER_H
#define ALIZDCMERGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  		Merger class for ZDC	      //
////////////////////////////////////////////////

class AliZDC;
class AliZDCHit;
class AliZDCMergedHit;
class AliZDCDigit;

typedef enum {kDigitize=0, kMerge=1} MergeMode_t;

class AliZDCMerger{

public:
   AliZDCMerger();
   virtual ~AliZDCMerger();
   
   void    InitMerging();
   void    Background(Float_t &b, Int_t &nspecn, Int_t &nspecp);
   void    Fragmentation(Float_t b, Int_t nspecn, Int_t nspecp,
                         Int_t &nfreespn, Int_t &nfreespp);
   void    Mixing();
   void    ExtractSignal(Int_t SpecType);
   
   // Inline functions to return TCA of MergerHits to Hits2SDigits()
   TClonesArray *MergedHits()   const {return fMHits;}
   int   GetNMhits()   		const {return fNMhits;}

   // Inline function to return background file Hits2SDigits()
   TFile *BgrFile()		const {return fBgrFile;}

   // Inline function to return event number
   Int_t EvNum()		const {return fNEvBgr;}
   
   // Setters 
   void SetMode(MergeMode_t mode)         {fMerge = mode;}
   void SetBackgroundFileName(char* file) {fFnBgr = file;}        
   void SetBackgroundEventNum(Int_t nev)  {fNEvBgr = nev;}        

   //Open the background file 
   TFile *OpenBgrFile(); 

protected:
   MergeMode_t  fMerge;		// Merging type kDigitize, kMerge
   
   // Background event
   char         *fFnBgr;	// Background file name
   TFile 	*fBgrFile;	// Pointer to background file
   Int_t	fNEvBgr;	// Number of events in background file
   TTree        *fTrHBgr;	// Hits tree for background event
   TClonesArray *fHitsBgr;	// TClonesArray of background hits
   Float_t 	fImpPar;	// Impact Parameter of the collision
   Int_t        fSpecn;		// Number of spectator n
   Int_t        fSpecp;		// Number of spectator p

   // Signal events
   Int_t        fFreeSpn;       // Signal event number x spectator n
   Int_t        fFreeSpp;       // Signal event number x spectator p

   char      	*fFnSpecn;      // Spectator n file name
   TFile      	*fSpecnFile;    // Pointer to signal file -> spectator n
   char       	*fFnSpecp;      // Spectator p file name
   TFile      	*fSpecpFile;    // Pointer to signal file -> spectator p
   
   Int_t	fNMhits;	// Number of Merged hits for background
   TClonesArray *fMHits;	// TCA for "merged" hits  

public:   
  // *** Digits
  // --- Parameters for conversion of light yield in ADC channels
  Float_t fPMGain[3][5];      // PM gain
  Float_t fADCRes;	      // ADC conversion factor
  // --- Digitization parameters setters and getters
  //  PM gain
  void SetPMGain(Int_t Det, Int_t PMDet, Int_t PMGain)
       {fPMGain[Det][PMDet] = PMGain;}
  Float_t GetPMGain(Int_t Det, Int_t PMDet)
       {return fPMGain[Det][PMDet];}
  //  Conversion factor from charge to ADC channels
  //	      F = 1.6E-19 / Resolution [Coulomb/ch]
  void SetADCRes(Int_t ADCRes) {fADCRes =  ADCRes;}
  Float_t GetADCRes() {return fADCRes;}

  
       
    ClassDef(AliZDCMerger,1)
};    
#endif
