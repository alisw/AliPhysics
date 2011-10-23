#ifndef ALIITSUPGRADEDIGITIZER_H
#define ALIITSUPGRADEDIGITIZER_H

/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////////////
//								//
//                  Digitizer class for ITS Upgrade		//
//								//
//								//
//    Authors: A.Mastroserio, C.Terrevoli   			//
//	       annalisa.mastroserio@cern.ch     		//
//	       cristina.terrevoli@ba.infn.it			//
//                                                              //
//////////////////////////////////////////////////////////////////


#include <AliDigitizer.h>
#include <TArrayD.h>

class TClonesArray;
class TObjArray;

class AliITSupgradeDigitizer : public AliDigitizer //TObject-TNamed-TTask-AliDigitizer-AliITSupgradeDigitizer
{
 public:
  AliITSupgradeDigitizer():AliDigitizer(),fNxCells(0),fNzCells(0),fNlayers(0)              {;}
    AliITSupgradeDigitizer(AliDigitizationInput *pRunDig):AliDigitizer(pRunDig),fNxCells(0),fNzCells(0),fNlayers(0){;}
      virtual ~AliITSupgradeDigitizer()                                                {;}
      void   SetConfiguration(TArrayD xcell, TArrayD zcell);
      void   Digitize(Option_t* option=0);                //virtual
   
      void    Sdigits2Digits(TClonesArray *pSDigitList,TObjArray *pDigitLst);
 
 protected:
 
      enum {maxLab=12}; // maximum number of MC labels associated to the digit (4 times as much as can be stored in the "mother class")

      TArrayD fNxCells;
      TArrayD fNzCells;
      Short_t fNlayers;    
 
      ClassDef(AliITSupgradeDigitizer,0)
	};    

#endif


