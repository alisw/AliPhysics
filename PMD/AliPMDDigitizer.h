#ifndef PMDDigitizer_H
#define PMDDigitizer_H
//-----------------------------------------------------//
//                                                     //
//  Header File : PMDDigitization.h, Version 00        //
//                                                     //
//  Date   : September 20 2002                         //
//                                                     //
//-----------------------------------------------------//

#include <Riostream.h>
#include <stdlib.h>
#include <math.h>
#include <TMath.h>

class TClonesArray;
class TFile;

class TObjArray;
class TParticle;
class TTree;
class TNtuple;

class AliLoader;
class AliRunLoader;
class AliRun;
class AliDetector;
class AliPMDhit;
class AliHit;
class AliHeader;

class AliPMDcell;
class AliPMDsdigit;
class AliPMDdigit;
class AliPMDClustering;

class AliPMDDigitizer
{
 protected:
  AliRunLoader *fRunLoader;
  AliRun       *gAlice;
  AliPMDhit    *pmdHit;   /* Pointer to specific detector hits. */
  AliDetector  *PMD;      /* Get pointers to Alice detectors 
			     and Hits containers */
  AliLoader    *pmdloader;

  TClonesArray *PMDhits;
  TObjArray    *Particles;
  TParticle    *particle;

  TTree        *treeH;
  TTree        *treeS;
  TTree        *treeD;

  TClonesArray *fSDigits;
  TClonesArray *fDigits;

  TObjArray    *fCell;
  AliPMDcell   *pmdcell;

  Int_t fNsdigit;
  Int_t fNdigit;
  Int_t fDetNo;
  Float_t fZPos;

  static const Int_t fTotUM = 24;
  static const Int_t fRow   = 48;
  static const Int_t fCol   = 96;
  Float_t fCPV[fTotUM][fRow][fCol];
  Float_t fPMD[fTotUM][fRow][fCol];
  Int_t   fPMDCounter[fTotUM][fRow][fCol];
  Int_t   fPMDTrackNo[fTotUM][fRow][fCol];
  Int_t   fCPVTrackNo[fTotUM][fRow][fCol];

 public:

  AliPMDDigitizer();
  virtual ~AliPMDDigitizer();

  void OpengAliceFile(char * /* galice.root */, Option_t * /* option */);

  void Hits2SDigits(Int_t /* ievt */);
  void Hits2Digits(Int_t /* ievt */);
  void SDigits2Digits(Int_t /* ievt */);
  void TrackAssignment2Cell();
  void MeV2ADC(Float_t /* mev */, Float_t & /* adc */);
  void AddSDigit(Int_t /* trnumber */, Int_t /* det */, Int_t /* smnumber */, 
		 Int_t /* cellnumber */, Float_t /* adc */);
  void AddDigit(Int_t /* trnumber */, Int_t /* det */, Int_t /* smnumber */, 
		Int_t /* cellnumber */, Float_t /* adc */);
  Int_t Convert2RealSMNumber(Int_t /* smnumber1 */ );
  void SetZPosition(Float_t /* zpos */);
  Float_t GetZPosition() const;
  void ResetCell();
  void ResetSDigit();
  void ResetDigit();
  void ResetCellADC();
  void UnLoad(Option_t * /* option */);

  ClassDef(AliPMDDigitizer,2)
};
#endif

