#ifndef ALIITSMODLUESPD_H
#define ALIITSMODLUESPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliITS.h"
#include "AliITSmodule.h"



//____________________________________________________________________
//
//  Class AliITSmoduleSPD
//  describes one SPD module
//  The main function of modules is to simulate DIGITS from  
//  GEANT HITS and produce POINTS from DIGITS
//                   
//  ver. 0.0    CERN, 17.09.1999  
// 
//___________________________________________________________________
//



class AliITSmoduleSPD: public AliITSmodule {

    
public:       
    
    //________________________________________________________________
    //
    // Constructors and deconstructor
    //________________________________________________________________
    //
    
    AliITSmoduleSPD();
    AliITSmoduleSPD(Int_t index);
    ~AliITSmoduleSPD();
  
    //________________________________________________________________
    //
    // Data process methods
    //________________________________________________________________
    //
    
    void HitToDigit() {};              // Not impemented yet
    void DigitToPoint() {};            // Not impemented yet
    void HitToPoint() {};              // Not impemented yet
    
       
  
protected:

public:
    ClassDef(AliITSmoduleSPD, 1)
};

#endif
