#ifndef AliHMPIDCalib_h
#define AliHMPIDCalib_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Class of HMPID to manage digits ---> pads
//.
//.
//.

//#include "TTreePlayer.h"
//#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
//#include <TH1S.h>
#include <TMath.h>
#include <TFile.h>
#include <TString.h>
#include "AliHMPIDParam.h"
#include "AliHMPIDRawStream.h"

class TFile;
class AliHMPIDCalib: public TObject { 


public:
  AliHMPIDCalib();
  virtual ~AliHMPIDCalib();
          void Init();
          void FillPedestal(Int_t pad,Int_t q);                             //absolute pad number and the charge of the pad
          void FillErrors(Int_t nDDL,Int_t nErrType, Int_t nErr);           //Fill the errors from RawStream
          void FillDDLCnt(Int_t iddl,Int_t inDDL, Int_t outDDL);            //Fill the errors from RawStream
        Bool_t CalcPedestal(Int_t nDDL, Char_t* name, Int_t nEv);           //number of the DDL, name of the output file and the number of events processed
        Bool_t CalcPedestalPaolo(Int_t nDDL, Char_t* name, Int_t nEv);      //number of the DDL, name of the output file and the number of events processed
        
        Bool_t WriteErrors(Int_t nDDL, Char_t* name, Int_t nEv);            //number of the DDL, name of the output file and the number of events processed
         void InitHisto(Int_t q,Int_t histocnt,Char_t* name);               //Init the pad histograms
         void FillHisto(Int_t histocnt,Int_t q);                            //Fill the ADC histograms
         void InitFile(Int_t inVal);                                        //Init the ADC histo output file (one per LDC or one per DDL)
         void CloseFile();                                                  //Close the file
         void SetRunParams(ULong_t runNum,Int_t timeStamp, Int_t ldcId);    //Set Run Parameters such as Run Number, TimeStamp, LDCid 
         void SetSigCut(Int_t nSigCut) { fSigCut=nSigCut;}                  //Set Sigma Cuts from Setter
         void SetSigCutFromFile(Char_t* name);                              //Set Sigma Cuts from File
         void SetSigCutFromShell(Char_t* name);                             //Set Sigma Cuts from Bash Shell
         void SetDaOutFromShell(Char_t* name);                              //Set out dir. of DA from Bash Shell
         void SetFeeInFromShell(Char_t* name);                              //Set out dir. for Fe2C from Bash Shell
             
         TString GetDaOutFromShell() {return  fDaOut;}                      //Get out dir. of DA from Bash Shell
         TString GetFeeInFromShell() {return  fFeeIn;}                      //Get out dir. for Fe2C from Bash Shell
         void SetWriteHistoPads(Bool_t isOn) {fWritePads=isOn;}             //Set wether ADC histos of pads are written or not
         void SetWriteHistoPads(Bool_t isOn,Bool_t isLarge,Int_t nDDL) {fWritePads=isOn;fLargeHisto=isLarge;fSelectDDL=nDDL;}             //Set wether ADC histos of pads are written or not
         Bool_t GetWritePads()       const{return fWritePads;}              //Set wether ADC histos of pads are written or not
         Bool_t GetLargePads()       const{return fLargeHisto;}             //Set wether ADC histos of pads are written or not
         Bool_t GetSelectedDDL()     const{return fSelectDDL;}              //Set wether ADC histos of pads are written or not
protected: 

    Bool_t     *faddl;                                                         //check is ddl is filled
    Float_t ****fsq;                                                           //Sum of pad Q
    Float_t ****fsq2;                                                          //Sum of pad Q^2
    Int_t   ****fnpc;                                                          //# of the pad was called with non zero charge
    Int_t   ****fpedQ0;                                                        //Check how many times a pad gives 0 charge in pedestal runs
    Int_t     **fErr;                                                          // Store the numner of errors for a given error type and a given DDL
    TH1I      **fPadAdc;                                                       //Charge distribution for pads    
    Bool_t     *fIsPad;                                                        //Check if the ADC histo for the pad is booked or not
    TFile      *fFile;                                                         //ADC histo output file (one per LDC)      
    UInt_t      fLdcId;                                                        //Ldc ID 
    UInt_t      fTimeStamp;                                                    //Time Stamp
    Int_t       fRunNum;                                                       //Run Number
    Int_t       fSigCut;                                                       //n. of pedestal distribution sigmas used to create zero suppresion table                          
    Bool_t      fWritePads;                                                    //Select wether to write ADC pad histograms or not
    Int_t      *fnDDLInStream;                                                 // if the DDL is in the raw data
    Int_t      *fnDDLOutStream;                                                // if the DDL is in the raw data
    Bool_t      fLargeHisto;                                                   //Default is kFALSE.if kTRUE then write large pad histograms with 4093 bins!!!! Only if you have 2GB of RAM!!!   
    Int_t       fSelectDDL;                                                    //Select the DDL to write for the in the large histograms. Only ONE at one time!

    TString     fDaOut;                                                        //Store the DA output files in this directory
    TString     fFeeIn;                                                        //
private:
  AliHMPIDCalib(const AliHMPIDCalib& c);              //dummy copy constructor
  AliHMPIDCalib &operator=(const AliHMPIDCalib& c);   //dummy assignment operator
    
    ClassDef(AliHMPIDCalib,3)                                                  //HMPID calibration and pedestal class        
};
#endif
