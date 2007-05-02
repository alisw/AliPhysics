#ifndef ALIMUONSIMDATA_H
#define ALIMUONSIMDATA_H
//
/*  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004
//
/// \ingroup sim
/// \class AliMUONSimData
/// \brief Class containing MUON data from simulation: hits, digits, globaltrigger, 
/// localtrigger, etc ...
///
//  Author: Gines Martinez, Subatech,  September 2003

#include "AliLoader.h"

#include "AliMUONData.h"

class AliRunLoader;

class TClonesArray;
class TObjArray;
class TTree;


//__________________________________________________________________
/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliMUONSimData                                              //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliMUONSimData : public AliMUONData
{
  public:
  
  //  enum EChamberIteration { kAllChambers, kTrackingChambers, kTriggerChambers };
  
    AliMUONSimData();
    AliMUONSimData(AliLoader * loader, const char* name, const char* title);
    AliMUONSimData(const char* galiceFile);
    virtual ~AliMUONSimData();  
    //virtual void   AddSDigit(Int_t id, const AliMUONDigit& digit); // use copy constructor
    virtual void   AddHit(Int_t fIshunt, Int_t track, Int_t detElemId, 
			  Int_t idpart, Float_t X, Float_t Y, Float_t Z, 
			  Float_t tof, Float_t momentum, Float_t theta, 
			  Float_t phi, Float_t length, Float_t destep, 
			  Float_t Xref,Float_t Yref,Float_t Zref);
    
    TClonesArray*  Hits() {return fHits;} ///< Return hits
    //TClonesArray*  SDigits(Int_t DetectionPlane) const;

    void           GetTrack(Int_t it) const;
    Int_t          GetNtracks() const;
    //void           GetSDigits() const;

    
    virtual void   Fill(Option_t* opt=" ");
    virtual void   MakeBranch(Option_t *opt=" ");
    virtual void   SetDataContainer(Option_t *opt=" ");
    virtual void   SetTreeAddress(Option_t *opt=" ");
    
    virtual void   ResetHits();
    //virtual void   ResetSDigits();
  
                   /// Return tree with hits
    TTree*         TreeH() {return fLoader->TreeH(); }
                   /// Return tree with summable digits
    //TTree*         TreeS() {return fLoader->TreeS(); }
                   /// Return tree with particles
    TTree*         TreeP() {return fLoader->TreeP(); }

                   // Methods to dump data
    void DumpKine(Int_t event2Check=0);
    void DumpHits(Int_t event2Check=0, Option_t* opt="full");
    //void DumpSDigits(Int_t event2Check=0, Option_t* opt="tracks");
    
  protected: 
    /// Not implemented
    AliMUONSimData(const AliMUONSimData& rhs);
    /// Not implemented
    AliMUONSimData& operator=(const AliMUONSimData& rhs);

    TClonesArray*   fHits;    ///< One event in treeH per primary track
    //TObjArray*      fSDigits; ///< One event in treeS and one branch per detection plane

    Int_t           fNhits;   //!< Number of Hits
    //Int_t*          fNSdigits;//!< Number of Digits

    mutable Int_t fCurrentEvent; ///< Current event we're dealing with
    
private:  
    void   FillOwn(Option_t* opt=" ");
    void   MakeOwnBranch(Option_t *opt=" ");
    void   SetOwnDataContainer(Option_t *opt=" ");
    void   SetOwnTreeAddress(Option_t *opt=" ");

  ClassDef(AliMUONSimData,3) // Data accessor for MUON module
      
};
// inline functions


/// Load hits for \a i th entry in hits three
inline void AliMUONSimData::GetTrack(Int_t it) const  {
  if (fLoader && fLoader->TreeH())
    fLoader->TreeH()->GetEvent(it);
}
/*
/// Load sdigits tree
inline void AliMUONSimData::GetSDigits() const {
  if (fLoader && fLoader->TreeS())
    fLoader->TreeS()->GetEvent(0);
}
*/


#endif

