#ifndef ALIPMDDDLDATA_H
#define ALIPMDDDLDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author : B.K. Nandi

#include <TObject.h>

class AliPMDddldata : public TObject {

public:
   AliPMDddldata();
   AliPMDddldata(const AliPMDddldata &ddl);
   AliPMDddldata& operator=(const AliPMDddldata &ddl);

   virtual ~AliPMDddldata();

   void  SetDetector(Int_t idet)      {fDetector = idet;}
   void  SetSMN(Int_t ismn)           {fSMN = ismn;}
   void  SetModule(Int_t imod)        {fModule = imod;}
   void  SetPatchBusId(Int_t ipbusid) {fPatchBus = ipbusid;}
   void  SetMCM(Int_t imcm)           {fMCM = imcm;}
   void  SetChannel(Int_t ich)        {fChannel = ich;}
   void  SetRow(Int_t irow)           {fRow = irow;}
   void  SetColumn(Int_t icol)        {fCol = icol;}
   void  SetSignal(Int_t isig)        {fSignal = isig;}
   void  SetParityBit(Int_t ibit)     {fBit = ibit;}

   Int_t  GetDetector() const   {return fDetector;}
   Int_t  GetSMN() const        {return fSMN;}
   Int_t  GetModule() const     {return fModule;}
   Int_t  GetPatchBusId() const {return fPatchBus;}
   Int_t  GetMCM() const        {return fMCM;}
   Int_t  GetChannel() const    {return fChannel;}
   Int_t  GetRow() const        {return fRow;}
   Int_t  GetColumn() const     {return fCol;}
   Int_t  GetSignal() const     {return fSignal;}
   Int_t  GetParityBit() const  {return fBit;}



 private:

   Int_t     fDetector;        // Detector  (PRE:0, CPV:47)
   Int_t     fSMN;             // Serial Module number (0 to 23 for each det)
   Int_t     fModule;          // Module number (0 to 47)
   Int_t     fPatchBus;        // Patch bus number
   Int_t     fMCM;             // MCM number
   Int_t     fChannel;         // Channel number
   Int_t     fRow;             // Row number
   Int_t     fCol;             // Column number
   Int_t     fSignal;          // ADC of the cell
   Int_t     fBit;             // Parity Bit

   ClassDef(AliPMDddldata,0)  // PMD DDL container
};
#endif
