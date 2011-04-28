#ifndef AliOADBPWG2Spectra_H
#define AliOADBPWG2Spectra_H

//-------------------------------------------------------------------------
//     OADB interface for the PWG2 spectra
//     Author: Michele Floris, CERN
//    
// -------------------------------------------------------------------------
#include "TNamed.h"
#include "TList.h"

class TList;
class TH1D;
class TH1;
class AliOADBPWG2Spectra : public TNamed {

 public :

  enum EPWG2SpectraDetector { kITSsa=0, kITSTPC, kTPC, kTOF, kTOFTPC, kNDetectors, kDetDummy };
  enum EPWG2SpectraPIDType  { kGaussFit=0, kNSigma, kBayes, kKinks, kNPIDTypes };
  enum EPWG2SpectraCharge   { kPos=0, kNeg, kNCharge };
  enum EPWG2SpectraParticle { kPion = 0, kKaon, kProton, kNParticle };
  AliOADBPWG2Spectra();
  AliOADBPWG2Spectra(const char* name);
  virtual ~AliOADBPWG2Spectra();

  void Init();

  static const char * GetOADBPWG2SpectraFileName();
  const char * GetHistoName(Int_t det, Int_t pidType, Int_t part, 
				   Int_t charge, const char * centrTag = 0, Int_t centrBin = -1) ;
  void AddHisto(TH1D * h, Int_t det, Int_t pidType, Int_t part, 
		  Int_t charge, const char * centrTag = 0, Int_t centrBin = -1) ;
  TH1D * GetHisto(Int_t det, Int_t pidType, Int_t part, 
		  Int_t charge, const char * centrTag = 0, Int_t centrBin = -1);
  TH1D * BookHisto(Int_t det, Int_t pidType, Int_t part, 
		   Int_t charge, const char * centrTag = 0, Int_t centrBin = -1) ;
    
  TH1D * GetHistoStandardBinning(const TH1D* h, Int_t det, Int_t pidType, Int_t part, 
				 Int_t charge, const char * centrTag, Int_t centrBin) ;

  static const char * GetDetectorName(Int_t i) {return fgkDetectorNames[i];}
  static const char * GetPIDName(Int_t i)      {return fgkPidTypeNames[i];}
  static const char * GetChargeName(Int_t i)   {return fgkChargeTags[i];}
  static const char * GetParticleName(Int_t i) {return fgkParticleNames[i];}

  Bool_t CompareBinning(TH1 * h1, TH1 * h2) ;

  virtual void Print (const Option_t * opt = "") const  { fHistos->Print(opt); }

// Browsable
  virtual Bool_t	IsFolder() const { return kTRUE; }
  void Browse(TBrowser *b);

  static const char * fgkDetectorNames[] ; // Detector tags
  static const char * fgkPidTypeNames[]  ; // Name of the PID technique
  static const char * fgkChargeTags[]    ; // tags for charges
  static const char * fgkParticleNames[] ; // Particle tags
  

 private:
  AliOADBPWG2Spectra(const AliOADBPWG2Spectra& cont); 
  AliOADBPWG2Spectra& operator=(const AliOADBPWG2Spectra& cont);

 private :
  
  TList * fHistos; // List containtig all the histograms

  ClassDef(AliOADBPWG2Spectra, 2);
};

#endif
