#ifndef ALI_TRANSPORTMONITOR__H
#define ALI_TRANSPORTMONITOR__H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

#include <map>

#ifndef ROOT_TStopwatch
#include "TStopwatch.h"
#endif

// Class that can be plugged in the simulation to monitor transport timing per 
// particle for each geometry volume.
//
//  andrei.gheata@cern.ch 

class TH2F;


//______________________________________________________________________________
class AliTransportMonitor : public TObject {
public:
  //________________________________________________________________
  class AliTransportMonitorVol : public TNamed {
  public:
    //___________________________________________________
    class AliPMonData {
    public:
      Int_t         fPDG;        // particle PDG
      Double_t      fEdt;        // Energy * dt integral
      Double_t      fTime;       // Total transport time for the particle in this volume
      AliPMonData() : fPDG(0), fEdt(0), fTime(0) {}
      virtual ~AliPMonData() {}
      ClassDef(AliPMonData, 1)     // Basic monitoring info structure
    };   
    //___________________________________________________   
    AliTransportMonitorVol();
    virtual ~AliTransportMonitorVol();
    
    Int_t           GetNtypes() const {return fNtypes;}
    void            StepInfo(Int_t    pdg,
                             Double_t energy, 
                             Double_t dt,
                             Double_t x, Double_t y, Double_t z);
    Double_t        GetTotalTime() const       {return fTotalTime;}
    Double_t        GetTime(Int_t itype) const {return fPData[itype].fTime;}
    Double_t        GetEmed(Int_t itype) const {return (fTotalTime>0)?fPData[itype].fEdt/fTotalTime : 0.;}
    Int_t           GetPDG(Int_t itype)  const {return fPData[itype].fPDG;}
    TH2F           *GetHistogram() const {return fTimeRZ;}
    private:
      AliPMonData  &GetPMonData(Int_t pdg);
      AliTransportMonitorVol(const AliTransportMonitorVol& other) : TNamed(other), fNtypes(0), fTotalTime(0), fPData(0), fTimeRZ(0), fParticles() {}
      AliTransportMonitorVol &operator=(const AliTransportMonitorVol&) {return *this;}
    private:
      Int_t         fNtypes;     // Number of different particle types
      Double_t      fTotalTime;  // Total time spent in this volume
      AliPMonData  *fPData;      //[fNtypes] Array of particle data
      TH2F         *fTimeRZ;     // Timing R-Z histogram per volume
      typedef std::map<Int_t, Int_t>         ParticleMap_t;
      typedef ParticleMap_t::iterator        ParticleMapIt_t;
      ParticleMap_t fParticles;  //! Map of stored particla data  
  
  ClassDef(AliTransportMonitorVol,1)  // Helper to hold particle info per volume
  };
  //________________________________________________________________
private:
  AliTransportMonitor(const AliTransportMonitor&other) : TObject(other), fTotalTime(0), fTimer(), fVolumeMon(0) {}
  AliTransportMonitor &operator=(const AliTransportMonitor&) {return *this;}
public:
  AliTransportMonitor();
  AliTransportMonitor(Int_t nvolumes);
  virtual ~AliTransportMonitor();
  
  void              StepInfo(Int_t    volId,
                             Int_t    pdg,
                             Double_t energy, 
                             Double_t x, Double_t y, Double_t z);
  void              Print(Option_t *volName="") const;
  void              DummyStep();
  void              Start();
  void              Stop();
  void              Export(const char *fname);
  static AliTransportMonitor *Import(const char *fname);
private:
  Double_t          fTotalTime;  // Total simulation time
  TStopwatch        fTimer;      //! Global timer
  TObjArray        *fVolumeMon;  // Array of monitoring objects per volume
  
ClassDef(AliTransportMonitor,1)  // Class to monitor timing per volume 
};
//______________________________________________________________________________


#endif //ALI_TRANSPORTMONITOR__H
