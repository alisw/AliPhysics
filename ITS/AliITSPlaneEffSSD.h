#ifndef ALIITSPLANEEFFSSD_H
#define ALIITSPLANEEFFSSD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliITSPlaneEff.h"

///////////////////////////////////////////
//                                       //
// ITS Plane Efficiency class            //
//       for SSD                         //
// Origin: Giuseppe.Bruno@ba.infn.it     //
///////////////////////////////////////////

/* $Id$ */
  
class AliITSPlaneEffSSD :  public AliITSPlaneEff {
 public:
    AliITSPlaneEffSSD(); // default constructor
    virtual ~AliITSPlaneEffSSD(); // destructror
    // copy constructor
    AliITSPlaneEffSSD(const AliITSPlaneEffSSD &source);
    // ass. operator
    AliITSPlaneEffSSD& operator=(const AliITSPlaneEffSSD &s);
    virtual AliITSPlaneEff& operator=(const AliITSPlaneEff &source);
    // Simple way to add another class (i.e. statistics). 
    AliITSPlaneEffSSD& operator +=( const AliITSPlaneEffSSD &add);
    // Getters for average Plane efficiency (icluding dead/noisy)
    Double_t PlaneEff(const UInt_t mod) const;
    Double_t ErrPlaneEff(const UInt_t mod) const;
    // Methods to update the Plane efficiency (specific of the SSD segmentation) 
    Bool_t UpDatePlaneEff(const Bool_t Kfound, const UInt_t mod);
    //
    enum {kNModule = 1698}; // The number of modules
    enum {kNChip = 6}; // The number of chips per side of a module (2 sides: 12 chips)
    enum {kNSide = 2}; // The number of sides of a module (p and n side)
    enum {kNStrip = 128}; // The number of strips per chip (in a module 2*768 strips)
//
//  UInt_t GetChip(const UInt_t col) const; // get the chip number (from 0 to kNChip)
//  Plane efficiency for active  detector (excluding dead/noisy channels)
//  access to DB is needed
    virtual Double_t LivePlaneEff(UInt_t mod) const;
    virtual Double_t ErrLivePlaneEff(UInt_t mod) const;
    // Compute the fraction of Live area (of the module for the SSD)
    virtual Double_t GetFracLive(const UInt_t mod) const;
    // Compute the fraction of bad (i.e. dead and noisy) area (of the module for the SSD)
    virtual Double_t GetFracBad(const UInt_t mod) const;
    virtual Bool_t WriteIntoCDB() const;
    virtual Bool_t ReadFromCDB(); // this method reads Data Members (statistics) from DataBase
    virtual Bool_t AddFromCDB()   // this method updates Data Members (statistics) from DataBase
      {AliError("AddFromCDB: Still To be implemented"); return kFALSE;}
   // method to locate a basic block from Detector Local coordinate (to be used in tracking)
   // see file cxx for numbering convention.
   // here idet runs from 0 to 747 for layer 4 and from 0 to 949 for layer 5
    UInt_t GetKeyFromDetLocCoord(Int_t ilay,Int_t idet, Float_t, Float_t locz) const;
    UInt_t Nblock() const; // return the number of basic blocks

 protected:
    virtual void Copy(TObject &obj) const;
    Int_t GetMissingTracksForGivenEff(Double_t eff, Double_t RelErr, UInt_t im) const;

// 
    Int_t fFound[kNModule];  // number of associated clusters in a given module
    Int_t fTried[kNModule];  // number of tracks used for module efficiency evaluation
 private:
    UInt_t GetKey(const UInt_t mod) const; // unique key to locate the basic 
                                           // block of the SSD (the module itself)
    UInt_t GetModFromKey(const UInt_t key) const;
    void GetBadInModule(const UInt_t mod, UInt_t& bad) const;

    ClassDef(AliITSPlaneEffSSD,1) // SSD Plane Efficiency class
};
//
inline UInt_t AliITSPlaneEffSSD::Nblock() const {return kNModule;}
//
#endif

