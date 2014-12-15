#ifndef ALIITSPLANEEFFSSD_H
#define ALIITSPLANEEFFSSD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TH1F.h>
#include <TH2I.h>
#include "AliITSPlaneEff.h"

class AliCDBId;

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
    //    virtual AliITSPlaneEff& operator=(const AliITSPlaneEff &source);
    // Simple way to add another class (i.e. statistics). 
    AliITSPlaneEffSSD& operator +=( const AliITSPlaneEffSSD &add);
    // Getters for average Plane efficiency (icluding dead/noisy)
    Double_t PlaneEff(const UInt_t mod) const;
    Double_t ErrPlaneEff(const UInt_t mod) const;
    // Getters for fFound[] and fTried[]
    Int_t GetFound(const UInt_t key) const;
    Int_t GetTried(const UInt_t key) const;
    // Methods to update the Plane efficiency (specific of the SSD segmentation) 
    Bool_t UpDatePlaneEff(const Bool_t Kfound, const UInt_t mod);
    //
    enum {kNModule = 1698}; // The number of modules
    enum {kNChip = 6}; // The number of chips per side of a module (2 sides: 12 chips)
    enum {kNSide = 2}; // The number of sides of a module (p and n side)
    enum {kNStrip = 128}; // The number of strips per chip (in a module 2*768 strips)
//
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
    Bool_t AddFromCDB(AliCDBId *cdbId);   // this method updates Data Members (statistics) from DataBase
    virtual Bool_t AddFromCDB() {AliCDBId *cdbId=0; return  AddFromCDB(cdbId);}
   // method to locate a basic block from Detector Local coordinate (to be used in tracking)
   // see file cxx for numbering convention.
   // here idet runs from 0 to 747 for layer 4 and from 0 to 949 for layer 5
    UInt_t GetKeyFromDetLocCoord(Int_t ilay,Int_t idet, Float_t, Float_t locz) const;
    UInt_t Nblock() const; // return the number of basic blocks
   // compute the geometrical limit of a basic block (chip) in detector local coordinate system
    Bool_t GetBlockBoundaries(const UInt_t key,Float_t& xmn,Float_t& xmx,Float_t& zmn,Float_t& zmx) const;
  // Methods for dealing with auxiliary histograms
    // method to set on/off the creation/updates of histograms (Histos are created/destroyed)
    virtual void   SetCreateHistos(Bool_t his=kFALSE)
         {fHis=his; if(fHis) {DeleteHistos(); InitHistos();} else DeleteHistos(); return; }
    virtual Bool_t FillHistos(UInt_t key, Bool_t found, Float_t *track, Float_t *cluster, Int_t *ctype,Float_t*);
    virtual Bool_t WriteHistosToFile(TString filename="PlaneEffSSDHistos.root",Option_t* option = "RECREATE");
    virtual Bool_t ReadHistosFromFile(TString filename="PlaneEffSSDHistos.root"); // histos must exist already !
                                                                          // This method increases the
                                                                          // statistics of histos by adding
                                                                          // those of the input file.
    UInt_t GetKey(const UInt_t mod) const; // unique key to locate the basic
                                           // block of the SSD (the module itself)
 protected:
    virtual void Copy(TObject &obj) const;
    Int_t GetMissingTracksForGivenEff(Double_t eff, Double_t RelErr, UInt_t im) const;
    UInt_t GetModFromKey(const UInt_t key) const;
    void GetBadInModule(const UInt_t mod, UInt_t& bad) const;

// 
    Int_t fFound[kNModule];  // number of associated clusters in a given module
    Int_t fTried[kNModule];  // number of tracks used for module efficiency evaluation
//
 private:
    enum {kNHisto = kNModule}; // The number of histograms: module by module.
    //enum {kNclu = 3};          // Build specific histos of residuals up to cluster size kNclu.
                               // If you change them, then you must change implementation of
                               // the method FillHistos.
    virtual void InitHistos();    // create histos by allocating memory for them
    virtual void DeleteHistos();  // deletete histos (memory is freed)
    void CopyHistos(AliITSPlaneEffSSD& target) const; // copy only histograms to target
//
    TH1F **fHisResX; //! histos with residual distribution (track-cluster) along local X (r-phi)
    TH1F **fHisResZ; //! histos with residual distribution (track-cluster) along local Z
    TH2F **fHisResXZ; //! 2-d histos with residual distribution (track-cluster) along local X and Z
    TH2I **fHisClusterSize; //! histos with cluster-size distribution
    TH1F **fHisTrackErrX; //! histos with track prediction error on Local X
    TH1F **fHisTrackErrZ; //! histos with track prediction error on Local Z
    TH1F **fHisClusErrX; //! histos with Local_X cluster error
    TH1F **fHisClusErrZ; //! histos with Local_Z cluster error
//
    ClassDef(AliITSPlaneEffSSD,3) // SSD Plane Efficiency class
};
//
inline UInt_t AliITSPlaneEffSSD::Nblock() const {return kNModule;}
inline Bool_t AliITSPlaneEffSSD::GetBlockBoundaries(const UInt_t key,Float_t& xmn,Float_t& xmx,
                                                    Float_t& zmn,Float_t& zmx) const {
//  This method return the geometrical boundaries of the active volume of a given
//  basic block, in the detector reference system.
//
if(key>=kNModule)
  {AliWarning("GetBlockBoundaries: you asked for a non existing key"); return kFALSE;}
const Float_t kDxDefault = 72960.; // For Plane Eff. purpouses, default values 
const Float_t kDzDefault = 40000.; // are precise enough !!!
const Float_t kconv = 1.0E-04;  //converts microns to cm.
xmn=-kconv*kDxDefault/2.; xmx=kconv*kDxDefault/2.;
zmn=-kconv*kDzDefault/2.; zmx=kconv*kDzDefault/2.;
return kTRUE;
}
//
inline Int_t AliITSPlaneEffSSD::GetFound(const UInt_t key) const {
 if(key>=kNModule) {AliWarning("GetFound: you asked for a non existing key"); return -1;}
 return fFound[key];
}
inline Int_t AliITSPlaneEffSSD::GetTried(const UInt_t key) const {
 if(key>=kNModule) {AliWarning("GetTried: you asked for a non existing key"); return -1;}
 return fTried[key];
}
//
#endif
