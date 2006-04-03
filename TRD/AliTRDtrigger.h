#ifndef ALITRDTRIGGER_H
#define ALITRDTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD trigger class                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObjArray.h>
#include <TClonesArray.h>
#include <TNamed.h>
#include <TFile.h>

#include "AliTRDzmaps.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

const Int_t kNplan     =   6;
const Int_t kNmaxZchan = 100;       // max number of tracklets per subchannel
const Int_t kNsubZchan =  16;       // total number of subchannels
const Int_t kNmaxTrk   =  12;       // max number of tracklets in one track (6*2)

class AliTRDltuTracklet : public TObject {
  
 public:

  AliTRDltuTracklet(Int_t det, 
		    Int_t row, 
		    Float_t rowz,
		    Float_t slope, 
		    Float_t offset, 
		    Float_t time, 
		    Int_t ncl,
		    Int_t label,
		    Float_t q);
  ~AliTRDltuTracklet(){};

  Bool_t  IsSortable() const { return kTRUE; }
  virtual Int_t   Compare(const TObject *o) const;
  Int_t   GetDetector() { return fDetector; };
  Int_t   GetPlane(Int_t det) { return ((Int_t) (det % kNplan)); };
  Int_t   GetRow()      { return fRow; };
  Int_t   GetNclusters(){ return fNclusters; };
  Float_t GetSlope()    { return fSlope; };
  Float_t GetOffset()   { return fY; };
  Float_t GetTime0()    { return fX; };
  Float_t GetRowz()     { return fRowz; };
  virtual Float_t GetYproj(Float_t xpl);
  virtual Float_t GetZproj(Float_t xpl);
  Int_t   GetLabel()    { return fLabel; };
  Float_t GetPt(Float_t field);
  Float_t GetQ() { return fQ; };

 protected:

  Float_t fX;                              // distance vertex to entrance window
  Float_t fY;                              // tracklet offset at entrance window
  Float_t fSlope;
  Float_t fRowz;
  Int_t   fDetector;
  Int_t   fRow;
  Int_t   fNclusters;
  Int_t   fLabel;
  Float_t fQ;                              // charge sum divided by number of clusters

  ClassDef(AliTRDltuTracklet,1)

};

class AliTRDltuTracklet;

class AliTRDgtuTrack : public TObject {

 public:

  AliTRDgtuTrack();
  AliTRDgtuTrack(const AliTRDgtuTrack& track);
  ~AliTRDgtuTrack(){};

  Bool_t  IsSortable() const { return kTRUE; }
  virtual Int_t   Compare(const TObject *o) const;

  virtual void Reset();
  void  ResetTracklets() { if(fTracklets) fTracklets->Delete(); };
  virtual void AddTracklet(AliTRDltuTracklet *trk);
  virtual AliTRDltuTracklet *GetTracklet(Int_t pos);
  TObjArray     *Tracklets() { 
    if(!fTracklets) fTracklets = new TObjArray(400); return fTracklets; 
  };
  Int_t          GetNtracklets() {
    if (fTracklets) return fTracklets->GetEntriesFast();
    return 0;
  };
  Float_t GetYproj()     { return fYproj; };
  Float_t GetZproj()     { return fZproj; };
  Float_t GetSlope()     { return fSlope; };
  Int_t   GetTracklets() { return fNtracklets; };
  Int_t   GetPlanes()    { return fNplanes; };
  Int_t   GetClusters()  { return fNclusters; };
  Float_t GetPt()        { return fPt; };
  Float_t GetPhi()       { return fPhi; };
  Float_t GetEta()       { return fEta; };
  Int_t   GetLabel()     { return fLabel; };

  virtual void Track(Float_t xpl, Float_t field);

  virtual void CookLabel();

  void  SetDetector(Int_t det) { fDetector = det; };
  Int_t GetDetector() { return fDetector; };

  virtual void MakePID();
  Float_t GetPID() { return fPID; };

  Bool_t  IsElectron() { return fIsElectron; };

 protected:

  TObjArray          *fTracklets;                   //! Array of LTU tracklets

  Float_t fYproj;                                   // Average values calculated
  Float_t fZproj;                                   // from the tracklets 
  Float_t fSlope;

  Int_t   fDetector;                                // First detector in the module

  Int_t   fNtracklets;                              // Number of tracklets
  Int_t   fNplanes;                                 // Number of TRD planes
  Int_t   fNclusters;                               // Total number of clusters

  Float_t fPt;                                      // Transverse momentum
  Float_t fPhi;                                     // Phi angle at the vertex
  Float_t fEta;                                     // Eta at the vertex
  Int_t   fLabel;                                   // Track label
  Float_t fPID;                                     // PID electron likelihood
  Bool_t  fIsElectron;                              // Electron flag

  ClassDef(AliTRDgtuTrack,1)

};

class AliTRDgtuTrack;
class AliTRDtrigParam;

class AliTRDmodule : public TObject {

 public:

  AliTRDmodule(AliTRDtrigParam *trigp);

  virtual void   Reset();

  virtual void   AddTracklet(Int_t det, 
			     Int_t row, 
			     Float_t rowz,
			     Float_t slope, 
			     Float_t offset, 
			     Float_t time, 
			     Int_t ncl,
			     Int_t label,
			     Float_t q);

  TObjArray     *Tracklets() { 
    if(!fTracklets) fTracklets = new TObjArray(400); return fTracklets; 
  };

  void           ResetTracklets() { if(fTracklets) fTracklets->Delete(); };
  void           SortTracklets()  { if(fTracklets) fTracklets->Sort(); };
  virtual AliTRDltuTracklet *GetTracklet(Int_t pos);
  virtual void   RemoveMultipleTracklets();
  virtual void   RemoveTracklet(Int_t pos);
  Int_t          GetNtracklets() {
    if (fTracklets) return fTracklets->GetEntriesFast();
    return 0;
  };

  virtual void   AddTrack();

  TObjArray     *Tracks() { 
    if(!fTracks) fTracks = new TObjArray(400); return fTracks; 
  };

  virtual void   ResetTracks();
  void           SortTracks()  { if(fTracks) fTracks->Sort(); };
  virtual AliTRDgtuTrack *GetTrack(Int_t pos);
  virtual void   RemoveMultipleTracks();
  virtual void   RemoveTrack(Int_t pos);
  Int_t          GetNtracks() {
    if (fTracks) return fTracks->GetEntriesFast();
    return 0;
  };

  virtual void   SortZ(Int_t cha);
  virtual void   InitZLUT();
  virtual void   FindTracks();
  virtual void   FindTracksCombi(Int_t zchan);

 protected:

  Float_t             fXprojPlane;                  //! X (time) coordinate of the
                                                    //  projection plane
  Float_t             fField;                       //! Magnetic field
  TObjArray          *fTracklets;                   //! Array of LTU tracklets
  TObjArray          *fTracks;                      //! Array of GTU tracks

  Int_t fZnchan[kNplan][kNsubZchan];                //! number of LTU tracklets in each 
                                                    //  subchannel
  Int_t fZtrkid[kNplan][kNmaxZchan][kNsubZchan];    //! list of LTU tracklet id's for 
                                                    //  each subchannel

  Float_t  fDeltaY;                        // Y (offset) matching window in the GTU
  Float_t  fDeltaS;                        // Slope matching window in the GTU

  AliTRDltuTracklet  *fLTUtrk;                      //! Current LTU tracklet
  AliTRDgtuTrack     *fGTUtrk;                      //! Current GTU track

  ClassDef(AliTRDmodule,1)

};

class AliTRDdigitsManager;
class AliTRDdataArrayI;

class AliRunLoader;
class AliRawReader;

class AliTRDmcmTracklet;
class AliTRDmcm;
class AliTRDmodule;

class TTree;

class AliTRDtrigger : public TNamed {

 public:  

  enum { kNMCM = 16 };

  AliTRDtrigger();
  AliTRDtrigger(const Text_t* name, const Text_t* title);
  AliTRDtrigger(const AliTRDtrigger &p);   
  virtual ~AliTRDtrigger();

  AliTRDtrigger &operator=(const AliTRDtrigger &p); 
  virtual void Copy(TObject &p) const;

  virtual void Init();

  void SetRunLoader(AliRunLoader *rl) { fRunLoader = rl; };
  virtual Bool_t Open(const Char_t *name, Int_t nEvent = 0);
  virtual Bool_t ReadDigits();
  virtual Bool_t ReadDigits(AliRawReader* rawReader);
  virtual Bool_t MakeTracklets();
  virtual Bool_t WriteTracklets(Int_t det);
  virtual Bool_t ReadTracklets(AliRunLoader *rl);

  virtual void   AddTracklet(Int_t det, Int_t row, Int_t seed, Int_t n);
  TObjArray     *Tracklets() { 
    if(!fTracklets) fTracklets = new TObjArray(400); return fTracklets; 
  };
  void           ResetTracklets() { if(fTracklets) fTracklets->Delete(); };
  virtual void   SetMCMcoordinates(Int_t imcm);
  virtual void   SetParameter(AliTRDtrigParam *trigp) { fTrigParam = trigp; };
  AliTRDtrigParam *GetParameter() { return fTrigParam; };

  virtual void   MakeTracks(Int_t det);

  AliTRDgtuTrack *GetTrack(Int_t i) const {
    return (AliTRDgtuTrack *)fTracks.UncheckedAt(i);
  }
  void AddTrack(const AliTRDgtuTrack *t, Int_t det) {
    AliTRDgtuTrack * track = new(fTracks[fTracks.GetEntriesFast()]) AliTRDgtuTrack(*t);
    track->SetDetector(det);
  }
  Int_t GetNumberOfTracks() const {return fTracks.GetEntriesFast();}

  Int_t GetNPrimary() { return fNPrimary; };

 protected:

  AliTRDtrigParam       *fTrigParam;                   //! Trigger class parameters
  AliRunLoader          *fRunLoader;                   //! Run Loader
  AliTRDdigitsManager   *fDigitsManager;               //! TRD digits manager
  TTree                 *fTrackletTree;                //! Tree with tracklets
  TObjArray             *fTracklets;                   //! Array of tracklets

  Int_t                  fNROB;                        //! Number of ROBs in the current chamber
  AliTRDmcm             *fMCM;                         //! Current MCM
  AliTRDmcmTracklet     *fTrk;                         //! Current tracklet
  AliTRDmodule          *fModule;                      //! Current module
  AliTRDgtuTrack        *fGTUtrk;                      //! Current GTU track

  Int_t                  fNtracklets;                  //! Tracklets counter

  AliTRDdataArrayI *fDigits;                           //! Array with digits
  AliTRDdataArrayI *fTrack0;                           //! Track dictionary 0
  AliTRDdataArrayI *fTrack1;                           //! Track dictionary 1
  AliTRDdataArrayI *fTrack2;                           //! Track dictionary 2

  Int_t fNPrimary;                                     //! Number of primary tracks

  TClonesArray           fTracks;                      //! Array of GTU tracks

  ClassDef(AliTRDtrigger,1)                            //  TRD trigger class

};

#endif
