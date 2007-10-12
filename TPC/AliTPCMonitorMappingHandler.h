#ifndef ALITPCMONITORMAPPINGHANDLER_H
#define ALITPCMONITORMAPPINGHANDLER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
// AliTPCMonitorMappingHandler class
//
// Class for handling mapping information TPC  
//  
// Authors: Roland Bramm, 
//          Stefan Kniege, IKF, Frankfurt
//       
/////////////////////////////////////////////////////////////////////////

#include "TNamed.h" 

class AliTPCMonitorMappingHandler:   public TNamed {
 public:
    
    AliTPCMonitorMappingHandler(Char_t* name, Char_t* title);
    AliTPCMonitorMappingHandler(const  AliTPCMonitorMappingHandler &maphand);
    AliTPCMonitorMappingHandler& operator= (const AliTPCMonitorMappingHandler& maphand);

    ~AliTPCMonitorMappingHandler();
    
    void     ReadMapping(char* mapfile);
    void     ReadRowMappingGlob(char* fpathtoMappingRowfile) ;
    Int_t    GetNumOfChannels();
    Int_t    GetSizeofArray();
    Short_t* GetLine(          Int_t channel);
    Int_t    GetIndex(         Int_t channel);
    Int_t    GetPadRow(        Int_t channel);
    Int_t    GetPad(           Int_t channel);
    Int_t    GetConnector(     Int_t channel);
    Int_t    GetPin(           Int_t channel);
    Int_t    GetFEC(           Int_t channel);
    Int_t    GetFECchannel(    Int_t channel);
    Int_t    GetFECconnector(  Int_t channel);
    Int_t    GetAltroChannel(  Int_t channel);
    Int_t    GetAltro(         Int_t channel);
    Int_t    GetPadAddInRow(   Int_t row, Int_t pad);
    Int_t    GetNumofPads(     Int_t row);
    
    Int_t    ReadFECMapping(   char* u2ftestfile);
    void     ReadfecHwMap(     Int_t sector);
    void     ReadfecGainMap(   char* fecgainmap);
     
    Int_t    U2fGetBranch(     Int_t fecnr);
    Int_t    U2fGetRCU(        Int_t fecnr);
    Int_t    U2fGetFECinRCU(   Int_t fecnr);
    Int_t    U2fGetFECinBranch(Int_t fecnr);
    Int_t    U2fGetSide(       Int_t fecnr);
    Int_t    U2fGetSector(     Int_t fecnr); 
    Int_t    U2fGetFECnr(      Int_t index);
    
    Int_t    GetFECfromHw(Int_t hw)            { return fMapHwFECglobal[hw][0];}
    Int_t    GetFECChfromHw(Int_t hw)          { return fMapHwFECglobal[hw][1];}
    Float_t  GetFECchGain(Int_t fec, Int_t ch) { return fecGainMap[fec][ch];}
    
 private:
    
 
    Int_t     fnumofChannels;            // Max number of channels
    Int_t     fmaxHWAdress;              // Max value of hardware addresses
    Int_t     fsizeofArray;              // Set to max value of hardware addresses
    Short_t** fmapping;                  // global  mapping array
    Int_t**   fmappingChannelinRow;      // mapping of hardware addresses in one pad row
    
    Short_t** fu2ftestmapping;           // mapping of global FEC numbers in sectors (determined during installation with U2F card)  
    Int_t**   fMapHwFECglobal;           // mapping of global FEC numbers to hardware addresses in these FECs              
    Float_t** fecGainMap;                // global gain calibration map
	
    ClassDef(AliTPCMonitorMappingHandler,1);
};

#endif
