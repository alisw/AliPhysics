#ifndef ALITOFRECONSTRUCTIONERV2_H
#define ALITOFRECONSTRUCTIONERV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
//  Task Class for Reconstruction V2 in TOF      
//                  
//-- Author: F. Pierella


#include "TTask.h"

class AliTOFDigitMap;
class AliTOFHitMap;
class TClonesArray;
class TString;
class TTree;
class TBranch;

class AliTOFReconstructionerV2: public TTask {

public:
  AliTOFReconstructionerV2() ;          // ctor
  AliTOFReconstructionerV2(char* tpcBackTracks, char* tofDigits="digits.root");
  AliTOFReconstructionerV2(const AliTOFReconstructionerV2 & rec);
  virtual ~AliTOFReconstructionerV2() ; // dtor   
  Bool_t        BackPropagation(){return kTRUE;};
  //void          CreateNTuple();
  void          Comparison(Int_t* rtIndex); // for MC comparison
  virtual void  Exec(Option_t* option); // do the main work
  Int_t         GetDbgFlag()           const   {return fdbg;}
  const char*   GetTOFDigitsFile()     const {return (char*)fTOFDigitsFile.Data();}
  const char*   GetTPCBackTracksFile() const {return (char*)fTPCBackTracksFile.Data();}
  void          Init(Option_t* opt);
  Int_t         LoadTPCTracks();
  Int_t         LoadTOFDigits();
  Int_t         LoadTRDTracks();
  Int_t         SaveTracks(const Char_t* outname="tofTracks.root", Int_t split=0);
  void          SetDbg(Int_t dbgflag)        {fdbg=dbgflag;}
  void          SetTOFDigitsFile(char * tofDigitsFile ) {fTOFDigitsFile=tofDigitsFile;}
  void          SetTPCBackTracksFile(char * tpcBackTracksFile ){fTPCBackTracksFile=tpcBackTracksFile;}
  void          SetField(Float_t field)      {fField=field;}
  void          SetNDummyTracks(Int_t nDummyTracks){fNDummyTracks=nDummyTracks;}
  void          SetScaleSigmaFactor(Float_t factor){fScaleSigmaFactor=factor;}
  void          SetStep(Float_t step)        {fStep=step;}
  Bool_t   operator == (const AliTOFReconstructionerV2 & tofrecv2) const ;

private:

  Int_t   fdbg;              //! Flag for debug, 0 no debug, 1 debug
  AliTOFDigitMap* fDigitsMap;//! pointer to the map of TOF digits
  //AliTOFHitMap* fDigitsMap;  //! pointer to the map of TOF digits
  Float_t fField;            //! mag field value [Tesla]
  Int_t   fNDummyTracks;     //  number of test tracks used to search
                             //  the signal on TOF
  Float_t fScaleSigmaFactor; //  scale factor for sigma (common value for sigmaY and sigmaZ)
  Float_t fStep;             //! step inside the TOF volumes during 
                             //  back propagation 
  TClonesArray* fTOFDigits;  //! pointer to the TClonesArray with TOF digits
  TClonesArray* fTOFTracks;  //! pointer to the TClonesArray with TOF tracks
  TString fTOFDigitsFile;    //! file with TOF digits 
  TString fTPCBackTracksFile;//! seed TPC to TOF file name
  TTree*  fKalmanTree       ;//! tree with reconstructed tracks in TPC
  TBranch* fBranchWithTracks;//! branch with backpropagated tracks in TPC

  virtual void  IsInsideThePad(Float_t x, Float_t y, Float_t z, Int_t *nGeom, Float_t& zPad, Float_t& xPad);
  virtual void  GetGlobalXYZ(Double_t alpha, Double_t& x, Double_t& y, Double_t& z);
  Bool_t        DigitFinder(TArrayI *secArray, TArrayI *plaArray, TArrayI *strArray, TArrayI *pdzArray, TArrayI *pdxArray, Int_t* assignedVol, Int_t* digitTrackArray, Float_t& tdc);

 protected:
  // the class is assumed to have a streamer
  ClassDef(AliTOFReconstructionerV2,1)  // Task class for TOF reconstruction V2

};

#endif // AliTOFRECONSTRUCTIONERV2_H
