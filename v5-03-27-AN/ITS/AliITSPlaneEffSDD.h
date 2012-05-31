#ifndef ALIITSPLANEEFFSDD_H
#define ALIITSPLANEEFFSDD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TH2I.h>
#include <TProfile.h>
#include "AliITSPlaneEff.h"

class AliCDBId;

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
    // Getters for fFound[] and fTried[]
    Int_t GetFound(const UInt_t key) const;
    Int_t GetTried(const UInt_t key) const;
    // Methods to update the Plane efficiency (specific of the SDD segmentation) 
    Bool_t UpDatePlaneEff(const Bool_t Kfound, const UInt_t mod, 
                          const UInt_t chip, const UInt_t wing, const UInt_t subw=0);
    Bool_t UpDatePlaneEff(const Bool_t Kfound, const UInt_t key)
      {return UpDatePlaneEff(Kfound,GetModFromKey(key),GetChipFromKey(key),
                             GetWingFromKey(key),GetSubWingFromKey(key));};
    //
    enum {kNModule = 260}; // The number of modules (i.e. detector 7.25*7.53 cm^2)
    enum {kNChip = 4}; // The number of chips per half module (i.e. per wing, in total 4+4 chips)
    enum {kNWing = 2}; // The number of wings (this is hardware division of the module) 
    enum {kNSubWing = 1}; // Eventually sub-divide each wing (by 2 ?) to account for different 
                          // efficiencies due to different drift times.  
    enum {kNAnode = 64};   // Number of channels/chip (i.e. anodes per chip)
    //enum {kNTimeBin = 174};   // granularity along drift direction (i.e. segmentation in r-phi)

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
    Bool_t AddFromCDB(AliCDBId *cdbId);   // this method updates Data Members (statistics) from DataBase
    virtual Bool_t AddFromCDB()  {AliCDBId *cdbId=0; return  AddFromCDB(cdbId);}
    // method to locate a basic block from Detector Local coordinate (to be used in tracking)
    // see file cxx for numbering convention.
    // here idet runs from 0 to 83 for layer 2 and from 0 to 175 for layer 3
    UInt_t GetKeyFromDetLocCoord(Int_t ilay,Int_t idet, Float_t locx, Float_t locz) const;
    UInt_t Nblock() const; // return the number of basic blocks
    // compute the geometrical limit of a basic block in detector local coordinate system
    Bool_t GetBlockBoundaries(const UInt_t key,Float_t& xmn,Float_t& xmx,Float_t& zmn,Float_t& zmx) const;
  // Methods for dealing with auxiliary histograms
    // method to set on/off the creation/updates of histograms (Histos are created/destroyed)
    void   SetCreateHistos(Bool_t his=kFALSE)
         {fHis=his; if(fHis) {DeleteHistos(); InitHistos();} else DeleteHistos(); return; }
    virtual Bool_t FillHistos(UInt_t key, Bool_t found, Float_t *track, Float_t *cluster, Int_t *ctype, Float_t *);
    virtual Bool_t WriteHistosToFile(TString filename="PlaneEffSDDHistos.root",Option_t* option = "RECREATE");
    virtual Bool_t ReadHistosFromFile(TString filename="PlaneEffSDDHistos.root"); // histos must exist already !
                                                                          // This method increases the
                                                                          // statistics of histos by adding
                                                                          // those of the input file.
    UInt_t GetKey(const UInt_t mod, const UInt_t chip,           // unique key to locate the
                  const UInt_t wing, const UInt_t subw=0) const; // basic block of the SDD
    // return chip [0,3] and wing [0,1] from the "absolute" chip number [0,7] as defined in AliITSsegmentationSDD
    void ChipAndWingFromChip07(const Int_t chip07, UInt_t& chip,  UInt_t& wing) const;
 protected:
    virtual void Copy(TObject &obj) const;
    Int_t GetMissingTracksForGivenEff(Double_t eff, Double_t RelErr, 
                                      UInt_t im, UInt_t ic, UInt_t iw,  UInt_t isw=0) const;
    UInt_t GetModFromKey(const UInt_t key) const;
    UInt_t GetChipFromKey(const UInt_t key) const;
    UInt_t GetWingFromKey(const UInt_t key) const;
    UInt_t GetSubWingFromKey(const UInt_t key) const;
    // getters for chip and wing numbers, given the anode number [0,511]
    UInt_t ChipFromAnode(const UInt_t anode) const; // return the chip number (from 0 to kNChip-1)
    UInt_t WingFromAnode(const UInt_t anode) const; // return the wing number (from 0 to kNWing-1)
    void   ChipAndWingFromAnode(const UInt_t anode,UInt_t& chip,UInt_t& wing) const;
    // return the Subwing  (from 0 to kNSubWing-1) from the cell time bin in the range
    // [0,ntb] and from the number of time bins
    UInt_t SubWingFromTimeBin(const Int_t tb, const Int_t ntb) const;

    void   ChipAndWingAndSubWingFromLocCoor(Float_t locx, Float_t locz,
                                 UInt_t& chip, UInt_t& wing, UInt_t& subw) const;
    //
    void GetAllFromKey(const UInt_t key, UInt_t& mod, UInt_t& chip,
                       UInt_t& wing, UInt_t& subw) const;
    void GetBadInBlock(const UInt_t key, UInt_t& bad) const;
// 
    Int_t fFound[kNModule*kNChip*kNWing*kNSubWing]; // number of associated clusters in a given block
    Int_t fTried[kNModule*kNChip*kNWing*kNSubWing]; // number of tracks used for efficiency evaluation

 private:
    enum {kNHisto = kNModule}; // The number of histograms: module by module.
    enum {kNclu = 3};          // Build specific histos of residuals up to cluster size kNclu.
                               // If you change them, then you must change implementation of
                               // the method FillHistos.

    virtual void InitHistos();     // create histos by allocating memory for them
    virtual void DeleteHistos();   // deletete histos (memory is freed)
    virtual void CopyHistos(AliITSPlaneEffSDD& target) const; // copy only histograms to target

    TH1F **fHisResX; //! histos with residual distribution (track-cluster) along local X (r-phi)
    TH1F **fHisResZ; //! histos with residual distribution (track-cluster) along local Z
    TH2F **fHisResXZ; //! 2-d histos with residual distribution (track-cluster) along local X and Z
    TH2I **fHisClusterSize; //! histos with cluster-size distribution
    //TH1F ***fHisResXclu; //! histos with residual distribution along local X (r-phi) for cluster type
    TProfile **fProfResXvsCluSizeX; //! TProfile of X Residuals vs. cluster size in X 
    TH1F ***fHisResZclu; //! histos with residual distribution along local Z for cluster type
    TProfile **fProfResXvsX; //! TProfile of X Residuals vs. X (of the cluster)
    TProfile **fProfResZvsX; //! TProfile of Z Residuals vs. X (of the cluster)
    TProfile **fProfClustSizeXvsX; //! TProfile of cluster_size_X vs. X (of the cluster)
    TProfile **fProfClustSizeZvsX; //! TProfile of cluster_size_X vs. X (of the cluster)
    TH1F **fHisTrackErrX; //! histos with track prediction error on Local X
    TH1F **fHisTrackErrZ; //! histos with track prediction error on Local Z
    TH1F **fHisClusErrX; //! histos with Local_X cluster error
    TH1F **fHisClusErrZ; //! histos with Local_Z cluster error

    ClassDef(AliITSPlaneEffSDD,2) // SDD Plane Efficiency class
};
//
inline UInt_t AliITSPlaneEffSDD::Nblock() const {return kNModule*kNChip*kNWing*kNSubWing;}

inline Int_t AliITSPlaneEffSDD::GetFound(const UInt_t key) const {
 if(key>=kNModule*kNChip*kNWing*kNSubWing) {AliWarning("GetFound: you asked for a non existing key"); return -1;}
 return fFound[key];
}
inline Int_t AliITSPlaneEffSDD::GetTried(const UInt_t key) const {
 if(key>=kNModule*kNChip*kNWing*kNSubWing) {AliWarning("GetTried: you asked for a non existing key"); return -1;}
 return fTried[key];
}
inline void AliITSPlaneEffSDD::ChipAndWingFromChip07(const Int_t chip07, UInt_t& chip, 
                                                           UInt_t& wing) const {
if(chip07<0 || chip07>7) 
  {AliWarning("ChipAndWingFromChip07:  you asked for a non existing chip"); return;}
else if(chip07<=3) { chip=chip07; wing=0;}
else {chip=chip07-kNChip; wing=1;} 
return;
}
//
#endif
