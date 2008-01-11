#ifndef ALIITSPLANEEFFSDD_H
#define ALIITSPLANEEFFSDD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliITSPlaneEff.h"

///////////////////////////////////////////
//                                       //
// ITS Plane Efficiency class            //
//       for SDD                         //
// Origin: Giuseppe.Bruno@ba.infn.it     //
///////////////////////////////////////////

/* $Id$ */
  
class AliITSPlaneEffSDD :  public AliITSPlaneEff {
 public:
    AliITSPlaneEffSDD(); // default constructor
    virtual ~AliITSPlaneEffSDD(); // destructror
    // copy constructor
    AliITSPlaneEffSDD(const AliITSPlaneEffSDD &source);
    // ass. operator
    AliITSPlaneEffSDD& operator=(const AliITSPlaneEffSDD &s);
    virtual AliITSPlaneEff& operator=(const AliITSPlaneEff &source);
    // Simple way to add another class (i.e. statistics). 
    AliITSPlaneEffSDD& operator +=( const AliITSPlaneEffSDD &add);
    // Getters for average Plane efficiency (icluding dead/noisy)
    Double_t PlaneEff(const UInt_t mod, const UInt_t chip, 
                      const UInt_t wing, const UInt_t subw=0) const;
    Double_t ErrPlaneEff(const UInt_t mod, const UInt_t chip,
                         const UInt_t wing, const UInt_t subw=0) const;
    Double_t PlaneEff(const UInt_t key) const 
       {return PlaneEff(GetModFromKey(key),GetChipFromKey(key),
                        GetWingFromKey(key),GetSubWingFromKey(key));};
    Double_t ErrPlaneEff(const UInt_t key) const 
       {return ErrPlaneEff(GetModFromKey(key),GetChipFromKey(key),
                           GetWingFromKey(key),GetSubWingFromKey(key));};
    // Methods to update the Plane efficiency (specific of the SDD segmentation) 
    Bool_t UpDatePlaneEff(const Bool_t Kfound, const UInt_t mod, 
                          const UInt_t chip, const UInt_t wing, const UInt_t subw=0);
    virtual Bool_t UpDatePlaneEff(const Bool_t Kfound, const UInt_t key)
      {return UpDatePlaneEff(Kfound,GetModFromKey(key),GetChipFromKey(key),
                             GetWingFromKey(key),GetSubWingFromKey(key));};
    //
    enum {kNModule = 260}; // The number of modules (i.e. detector 7.25*7.53 cm^2)
    enum {kNChip = 4}; // The number of chips per half module (i.e. per wing, in total 4+4 chips)
    enum {kNWing = 2}; // The number of wings (this is hardware division of the module) 
    enum {kNSubWing = 1}; // Eventually sub-divide each wing (by 2 ?) to account for different 
                          // efficiencies due to different drift times.  
    enum {kNAnode = 64};   // Number of channels/chip (i.e. anodes per chip)
    enum {kNTimeBin = 72};   // granularity along drift direction (i.e. segmentation in r-phi)
//
//  Plane efficiency for active  detector (excluding dead/noisy channels)
//  access to DB is needed
    virtual Double_t LivePlaneEff(UInt_t key) const;
    Double_t LivePlaneEff(const UInt_t mod, const UInt_t chip, 
                          const UInt_t wing, const UInt_t subw=0) const
       {return LivePlaneEff(GetKey(mod,chip,wing,subw));};
    virtual Double_t ErrLivePlaneEff(UInt_t key) const;
    Double_t ErrLivePlaneEff(const UInt_t mod, const UInt_t chip,
                             const UInt_t wing, const UInt_t subw=0) const
       {return ErrLivePlaneEff(GetKey(mod,chip,wing,subw));};
    // Compute the fraction of Live area (of the CHIP/SubWing for the SDD)
    virtual Double_t GetFracLive(const UInt_t key) const;
    // Compute the fraction of bad (i.e. dead and noisy) area (of the CHIP/SubWing for the SDD)
    virtual Double_t GetFracBad(const UInt_t key) const;
    virtual Bool_t WriteIntoCDB() const;
    virtual Bool_t ReadFromCDB(); // this method reads Data Members (statistics) from DataBase
    virtual Bool_t AddFromCDB()   // this method updates Data Members (statistics) from DataBase
      {AliError("AddFromCDB: Still To be implemented"); return kFALSE;}

 protected:
    virtual void Copy(TObject &obj) const;
    Int_t GetMissingTracksForGivenEff(Double_t eff, Double_t RelErr, 
                                      UInt_t im, UInt_t ic, UInt_t iw,  UInt_t isw=0) const;
// 
    Int_t fFound[kNModule*kNChip*kNWing*kNSubWing]; // number of associated clusters in a given block
    Int_t fTried[kNModule*kNChip*kNWing*kNSubWing]; // number of tracks used for efficiency evaluation
 private:
    UInt_t GetKey(const UInt_t mod, const UInt_t chip,           // unique key to locate the 
                  const UInt_t wing, const UInt_t subw=0) const; // basic block of the SDD
    UInt_t GetModFromKey(const UInt_t key) const;
    UInt_t GetChipFromKey(const UInt_t key) const;
    UInt_t GetWingFromKey(const UInt_t key) const;
    UInt_t GetSubWingFromKey(const UInt_t key) const;
    // getters for chip and wing numbers, given the anode number [0,511]
    UInt_t ChipFromAnode(const UInt_t anode) const; // return the chip number (from 0 to kNChip-1)
    UInt_t WingFromAnode(const UInt_t anode) const; // return the wing number (from 0 to kNWing-1)
    void   ChipAndWingFromAnode(const UInt_t anode,UInt_t& chip,UInt_t& wing) const; 
    //
    void GetAllFromKey(const UInt_t key, UInt_t& mod, UInt_t& chip, 
                       UInt_t& wing, UInt_t& subw) const;
    void GetBadInBlock(const UInt_t key, UInt_t& bad) const;

    ClassDef(AliITSPlaneEffSDD,1) // SDD Plane Efficiency class
};
#endif

