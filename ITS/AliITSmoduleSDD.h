#ifndef ALIITSMODLUESDD_H
#define ALIITSMODLUESDD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliITS.h"
#include "AliITSmodule.h"



//____________________________________________________________________
//
//  Class AliITSmoduleSDD
//  describes one SDD module
//  The main function of modules is to simulate DIGITS from  
//  GEANT HITS and produce POINTS from DIGITS
//                   
//  ver. 0.0    CERN, 17.09.1999  
// 
//___________________________________________________________________
//



class AliITSmoduleSDD: public AliITSmodule {

    
public:       
    
    //________________________________________________________________
    //
    // Constructors and deconstructor
    //________________________________________________________________
    //
    
    AliITSmoduleSDD();
    AliITSmoduleSDD(Int_t index);
    ~AliITSmoduleSDD();
  
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
    ClassDef(AliITSmoduleSDD, 1)
};

#endif
