#ifndef ALIITSMODULE_H
#define ALIITSMODULE_H

//#include "AliITSpoint.h"
#include "AliITS.h"
#include "AliITSdigit.h"
#include "AliITShit.h"
#include "TObjArray.h"
#include "TArray.h"


//________________________________________________________________
//
//  Class AliITSmodule 
//  is a superclass for AliITSmoduleSSD, SPD and SDD.
//  The main function of modules is to simulate DIGITS from  
//  GEANT HITS and produce POINTS from DIGITS
//  It also make fast simulation without use of DIGITS
//
//  created by: A.Boucham, W.Peryt, S.Radomski, P.Skowronski
//              R.Barbera, B. Batynia, B. Nilsen
//  ver. 1.0    CERN, 16.09.1999  
// 
//________________________________________________________________
//


class AliITSmodule : public TObject {

public:

    //________________________________________________________________
    //
    // Constructors and deconstructor
    //________________________________________________________________
    //
    
    AliITSmodule();             // default constructor
    AliITSmodule(Int_t index);   // index in TObjectArray in ITS object

    virtual ~AliITSmodule() ;

    //________________________________________________________________
    //
    // Position menagenent
    //________________________________________________________________
    //
    
    inline Int_t GetIndex()  { return fIndex;}
    //inline Int_t GetLayer()  { return fLayer;}
    //inline Int_t GetLadder() { return fLadder;}
    //inline Int_t GetDet()    { return fDet;}
    
    
    //________________________________________________________________
    //
    // Hits menagement
    //________________________________________________________________
    //
    
    inline Int_t GetNhits() { return fNhitsM;} 
    // returns number of hits in this module
                  
    inline TObjArray *GetHits() { return fHitsM; }
    // returns pointer to array (TClonesArray) of pointers to hits
    
    Int_t  AddHit(AliITShit *hit); 
    // Adds pointer of hit belonging to this module
    // and returns number of hits in this module
    
    //________________________________________________________________ 
    //
    // Full Simulation
    //________________________________________________________________
    //
    
    virtual void HitToDigit() {};
    // this functon is virtual, becouse each type of a detector 
    // make simulation in its own specific methods
    
    //________________________________________________________________
    //
    // Reconstruction
    //________________________________________________________________
    //
    
    virtual void DigitToPoint() {};
    // this functon is virtual, becouse each type of a detector 
    // make reconstruction in its own specific methods
    // classes DIGIT and POINT are specyfic to detector type
    
    //________________________________________________________________
    //
    // Fast Simulation
    //________________________________________________________________
    //
	 
    virtual void HitToPoint() {};
    // this functon is virtual, becouse each type of a detector 
    // make fast simulation in its own specific methods
    // this simulation do not use DIGITS 
    
    
protected:
    
    //________________________________________________________________
    //
    // Data members
    //________________________________________________________________
    //
    
    AliDetector *fITS;      // Pointer to ITS detector
    Int_t       fIndex;     // Index of this module in ITSmodules TObjectArray
    
    TObjArray   *fHitsM;    // Pointer to list of hits on this module
    Int_t       fNhitsM;    // Number of hits 
    
    TArrayI     *fIdigits;  // Indexes of DIGITS belonging to this module
                            // in array in ITS
    Int_t       fNdigits;   // Number of DIGITS
                               
    TArrayI     *fIpoints;  // Indexes of POINTS belonging to this module
                            // in array in ITS  
    Int_t       fNpoints;   // Number of POINTS
      
public:

    //________________________________________________________________
    //
    // ROOT compatibility
    //________________________________________________________________
    //
    ClassDef(AliITSmodule,1)  
};

#endif
