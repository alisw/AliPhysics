#ifndef ALIITSPLANEEFF_H
#define ALIITSPLANEEFF_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TString.h>
#include "AliPlaneEff.h"

class AliITSsegmentation;
class TF1;
class AliITSgeom;
class AliLog;

////////////////////////////////////////////////////
//                                                //
// ITS virtual base class for Plane Efficiency    //
// Origin: Giuseppe.Bruno@ba.infn.it              //
//                                                //
////////////////////////////////////////////////////

/* $Id$ */

class AliITSPlaneEff : public AliPlaneEff {
 public:
 
    AliITSPlaneEff();// Default constructor
    // Standard constructor
    virtual ~AliITSPlaneEff(){;};
    // copy constructor. See detector specific implementation.
    AliITSPlaneEff(const AliITSPlaneEff &source);
    // Assignment operator. See detector specific implementation.
    AliITSPlaneEff& operator=(const AliITSPlaneEff &source);
    // Simple way to add another class (i.e. statistics). 
    //AliITSPlaneEff& operator +=( const AliITSPlaneEff &){return *this};
    // Average Plane efficiency (including dead/noisy)
    Int_t   GetRunNumber() const {return fRunNumber;}
    void    SetRunNumber(Int_t n) {fRunNumber=n;}
    //
    Double_t PlaneEff(Int_t nfound,Int_t ntried) const;     
    Double_t ErrPlaneEff(Int_t nfound,Int_t ntried) const; 
    virtual void GetPlaneEff(Int_t nfound,Int_t ntried,Double_t &eff, Double_t &err) const
        {eff=PlaneEff(nfound,ntried); err=ErrPlaneEff(nfound,ntried); return;};
    //
    virtual Double_t PlaneEff(const UInt_t key) const=0;
    // Plane efficiency for active  detector (excluding dead/noisy channels)
    virtual Double_t LivePlaneEff(UInt_t) const
       {AliWarning("This method gives just a rough estimate of the live-Det Efficiency!"); 
        return -1.;};
    virtual Double_t ErrLivePlaneEff(UInt_t) const
       {AliError("This method must be implemented in a derived class"); return -1.;};
    // Compute the fraction of Live detector
    virtual Double_t GetFracLive(const UInt_t) const
       {AliError("This method must be implemented in a derived class"); return -1.;}; 
    // Compute the fraction of Bad (i.e. dead + noisy) detector 
    virtual Double_t GetFracBad(const UInt_t) const
       {AliError("This method must be implemented in a derived class"); return -1.;}; 
    // Update the Counting of the plane efficiency
    virtual Bool_t UpDatePlaneEff(const Bool_t, const UInt_t) 
       {AliError("This method must be implemented in a derived class"); return kFALSE;};
    // Estimate of the number of tracks needed for measuring efficiency within RelErr
    virtual Int_t GetNTracksForGivenEff(Double_t eff, Double_t RelErr) const;
    void SetDefaultStorage(const char* uri);
    // Write into the data base 
    virtual Bool_t WriteIntoCDB() const
       {AliError("This method must be implemented in a derived class"); return kFALSE;};
    virtual Bool_t ReadFromCDB()
       {AliError("This method must be implemented in a derived class"); return kFALSE;};
    virtual Bool_t AddFromCDB()
       {AliError("This method must be implemented in a derived class"); return kFALSE;};
    // method to locate a basic block from Detector Local coordinate 
    virtual UInt_t GetKeyFromDetLocCoord(Int_t, Int_t, Float_t, Float_t) const
      {AliError("This method must be implemented in a derived class"); return 999999;};
    virtual UInt_t Nblock() const // return the number of basic blocks
      {AliError("This method must be implemented in a derived class"); return 999999;};
    virtual Bool_t GetBlockBoundaries(const UInt_t,Float_t&,Float_t&,Float_t&,Float_t&) const
      {AliError("This method must be implemented in a derived class"); return kFALSE;};
  // Methods for dealing with auxiliary histograms
    // method to set on/off the creation/updates of histograms (Histos are created/destroyed)
    virtual void   SetCreateHistos(Bool_t)
      {AliError("This method must be implemented in a derived class"); return; }
    virtual Bool_t GetCreateHistos() const {return fHis;};
    virtual Bool_t FillHistos(UInt_t, Bool_t, Float_t*, Float_t*, Int_t*, Float_t*)
      {AliError("This method must be implemented in a derived class"); return kFALSE; }
    virtual Bool_t WriteHistosToFile(TString ,Option_t*)
      {AliError("This method must be implemented in a derived class"); return kFALSE; }
    virtual Bool_t ReadHistosFromFile(TString )
      {AliError("This method must be implemented in a derived class"); return kFALSE; }
    void InitCDB();

 protected:

    virtual void Copy(TObject &obj) const;
    void NotImplemented(const char *method) const {if(gDebug>0)
         Warning(method,"This method is not implemented for this sub-class");}
    Int_t	fRunNumber;	//! run number (to access CDB)
    TString	fCDBUri;	//! Uri of the default CDB storage
    Bool_t	fInitCDBCalled;	//! flag to check if CDB storages are already initialized
    Bool_t      fHis;           //! if true, then histograms are created and filled 
   
 private:
    //Int_t*	fFound;		// number of associated clusters into a given block (e.g. SPD 1200 chip)
    //Int_t*	fTries;		// number of exspected  clusters into a given block (e.g. SPD 1200 chip)
    //Int_t	fRunNumber;	// run number (to access CDB)

    ClassDef(AliITSPlaneEff,2) // ITS Plane Efficiency virtual base class 
};
#endif
