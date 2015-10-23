#ifndef ALIITSUSUZE02_H
#define ALIITSUSUZE02_H

#include "TMatrixFfwd.h"
class TH1F;

//*******************************************************************
//
//  Simulation of the SUZE02 readout
//  Origin: Serhiy.Senuykov@cern.ch
//  (see macros/ScanDigitsSuze02*.C for an example of use)
//
//*******************************************************************

class AliITSUSuze02{
  public:
  AliITSUSuze02(Int_t Nrows, Int_t Ncols);
  AliITSUSuze02(const AliITSUSuze02& suze);
  AliITSUSuze02 &operator=(const AliITSUSuze02& suze);
  virtual ~AliITSUSuze02();
  
  void SetEncodingWindowSize(Int_t Wrows, Int_t Wcols);
  void SetQuotas(Int_t q32, Int_t qHalfFSBB, Int_t qFSBB);
  //void InitHistos(); 
  void AddDigit(Int_t row, Int_t col);
  void Process(TH1F* OverflowCodes, TH1F* NDigitsPerEncodingWindowDist, Bool_t Verbose=kFALSE); 
  //void Process(Bool_t Verbose=kFALSE); 
       
  void GetResults();  
  Int_t GetNDigitsEncoded() {return fNDigitsEncoded;}
  Int_t GetNEncodedWindows() {return fNEncodedWindows;} 
  Int_t GetNDigitsLost() {return fNDigitsLost;}
  Int_t GetNLostWindows() {return fNLostWindows;}    
  
  Int_t GetDataSize() {return fDataSizePerChip;}  
  
  Int_t GetNWindowsPer32colsMin() {return fNWindowsPer32colsMin;}
  Int_t GetNWindowsPerHalfFSBBMin() {return fNWindowsPerHalfFSBBMin;}
  Int_t GetNWindowsPerFSBBMin() {return fNWindowsPerFSBBMin;}
  
  void ResetChip();
  
  private:         
  static const Int_t kNumberOfFSBB=3; 
  static const Int_t kNumberOfHalfFSBB=2;
  
  //matrix to be processed by Suze02   
  Int_t fNRowsChip;
  Int_t fNColsChip;
  TMatrixF* fChip; 
               
  //Suze02 parameters
  Int_t fTestColumnSize;   //Number of rows of the encoding window
  Int_t fTestRowSize;   //Number of columns of the encoding window
  
  //Suze02 quotas
  Int_t fNWindowsPer32colsMax;
  Int_t fNWindowsPerHalfFSBBMax;
  Int_t fNWindowsPerFSBBMax;
  
  //results    
  Int_t fNDigitsEncoded; 
  Int_t fNEncodedWindows;
  
  Int_t fNDigitsLost;
  Int_t fNLostWindows;
  
  //TH1F* fOverflowCodes;
  //TH1F* fNDigitsPerEncodingWindowDist;
  
  Int_t fDataSizePerChip;
  
  Int_t fNWindowsPer32colsMin;
  Int_t fNWindowsPerHalfFSBBMin;
  Int_t fNWindowsPerFSBBMin;

  ClassDef(AliITSUSuze02,1)
};

#endif
