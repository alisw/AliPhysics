#ifndef ALIZDCMERGER_H
#define ALIZDCMERGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  		Merger class for ZDC	      //
////////////////////////////////////////////////

class AliZDC;
class AliZDCHit;
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
   void    Mixing(Int_t nfreespn, Int_t nfreespp);
   void    ExtractSignal(Int_t detector, Int_t quadrant, Float_t &quadLight, 
                         Float_t &comLight);
//   void    MergeHit(Int_t track, Int_t *vol, Float_t *hits);

   void    Digitize();
   Int_t   Phe2ADCch(Int_t Detector, Int_t Quadrant, Int_t Light);
   Int_t   AddPedestal();
   
   // Setters 
   void SetBackgroundFileName(char* file) {fFnBgr = file;}        
   void SetMode(MergeMode_t mode)         {fMerge = mode;}
     
private:
   //Open the background file 
   TFile *OpenBgrFile(); 

protected:
   MergeMode_t  fMerge;		// Merging type kDigitize, kMerge
   
   // Background event
   Int_t	fNEvBgr;	// Number of events in the background file
   char         *fFnBgr;	// Background file name
   TFile 	*fBgrFile;	// Pointer to background file
   TTree        *fTrHBgr;	// Hits tree for background event
   TTree        *fTrSDBgr;	// SDigits tree for background event
   Float_t 	fImpPar;	// Impact Parameter of the collision
   Int_t        fSpecn;		// Number of spectator n
   Int_t        fSpecp;		// Number of spectator p
   TClonesArray *fHitsBgr;	// TClonesArray of background hits

   // Signal events
   Int_t        fFreeSpn;       // Signal event number x spectator n
   Int_t        fFreeSpp;       // Signal event number x spectator p
      
   // File containing hit histograms for spectators
   char       	*fFnSpecn;      // Spectator n file name
   TFile      	*fSpecnFile;    // Pointer to signal file -> spectator n
   char       	*fFnSpecp;      // Spectator p file name
   TFile      	*fSpecpFile;    // Pointer to signal file -> spectator p
   Float_t	fQuadLight;     // Light produced in tower PM
   Float_t	fComLight;      // Light produced in common PM
   
   Int_t 	fTrack;

//  // *** Digits
//  // --- Digitization parameters setters and getters
//  // 	PM gain
//  void SetPMGain(Int_t Det, Int_t PMDet, Int_t PMGain)
//       {fPMGain[Det][PMDet] = PMGain;}
//  Float_t GetPMGain(Int_t Det, Int_t PMDet)
//       {return fPMGain[Det][PMDet];}
//  // 	Conversion factor from charge to ADC channels
//  //   	F = 1.6E-19 / Resolution [Coulomb/ch]
//  void SetADCRes(Int_t ADCRes) {fADCRes =  ADCRes;}
//  Float_t GetADCRes() {return fADCRes;}
//
//  // --- Parameters for conversion of light yield in ADC channels
//  Float_t fPMGain[3][5];      // PM gain
//  Float_t fADCRes;            // ADC conversion factor
  
       
    ClassDef(AliZDCMerger,0)
};    
#endif
