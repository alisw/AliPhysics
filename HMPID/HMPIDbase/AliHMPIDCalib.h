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
#include <THnSparse.h>
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
        Bool_t CalcPedestal(Int_t nDDL, Char_t* name, Char_t *name2,Int_t nEv);           //number of the DDL, name of the output file and the number of events processed
        
        Bool_t WriteErrors(Int_t nDDL, Char_t* name, Int_t nEv);            //number of the DDL, name of the output file and the number of events processed
         void SetRunParams(ULong_t runNum,Int_t timeStamp, Int_t ldcId);    //Set Run Parameters such as Run Number, TimeStamp, LDCid 
         void SetSigCut(Int_t nSigCut) { fSigCut=nSigCut;}                  //Set Sigma Cuts from Setter
         void SetSigCutFromFile(TString name);                              //Set Sigma Cuts from File
         void SetDeadChannelMapFromFile(TString name);
         Bool_t GetSelectedDDL()     const  {return fSelectDDL;}              //Set wether ADC histos of pads are written or not
         TH2F   *GetPedMeanMap(Int_t iDDL)  {return fPedMeanMap[iDDL];}       //Get the pedestal mean map for a DDL to send to AMORE
         TH2F   *GetPedSigMap(Int_t iDDL)   {return fPedSigMap[iDDL];}        //Get the pedestal sigma map for a DDL to send to AMORE
         TH1F   *GetPedMean(Int_t iChFee)   {return f1DPedMean[iChFee];}      //Get the pedestal mean map for a FEE channel to send to AMORE
         TH1F   *GetPedSigma(Int_t iChFee)  {return f1DPedSigma[iChFee];}     //Get the pedestal Sigma map for a FEE channel to send to AMORE
     THnSparse  *GetDeadMap()               {return fDeadMap;}                //Get the masked channel map from the DAQ database
         Int_t   GetNumMaskedPads()         {return fNumMaskedPads;}          //Get the number of masked channels
         Int_t   GetNumDeadPads()           {return fNumDeadPads;}            //Get the number of masked channels
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
    Int_t      *fnDDLInStream;                                                 // if the DDL is in the raw data
    Int_t      *fnDDLOutStream;                                                // if the DDL is in the raw data
    Bool_t      fLargeHisto;                                                   //Default is kFALSE.if kTRUE then write large pad histograms with 4093 bins!!!! Only if you have 2GB of RAM!!!   
    Int_t       fSelectDDL;                                                    //Select the DDL to write for the in the large histograms. Only ONE at one time!
    THnSparse  *fDeadMap;                                                      //Dead Channel Map
    TH2F       **fPedMeanMap;                                                  //2D mean pedestal map to export to AMORE
    TH2F       **fPedSigMap;                                                   //2D pedestal sigma map to export to AMORE
    TH1F      **f1DPedMean;                                                    //1D mean pedestal map to export to AMORE
    TH1F      **f1DPedSigma;                                                   //1D pedestal sigma map to export to AMORE
    Int_t       fNumMaskedPads;                                                //Number of masked pads     
    Int_t       fNumDeadPads;                                                  //Number of currently dead channels   
            
  private:
                                           
  AliHMPIDCalib(const AliHMPIDCalib& c);                                       //dummy copy constructor
  AliHMPIDCalib &operator=(const AliHMPIDCalib& c);                            //dummy assignment operator
     
    ClassDef(AliHMPIDCalib,5)                                                  //HMPID calibration and pedestal class        
};
#endif

