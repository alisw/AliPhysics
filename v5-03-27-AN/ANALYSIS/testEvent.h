#ifndef Ana_Event
#define Ana_Event

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Simple event used for a simple analysis. The event has just a        //
//    TClonesArray of tracks that can be either pi+/-, protons or       //
//    gammas. Some ratio of all particles are pi0's that decay in       //
//    2 gammas.                                                         //
//                                                                      //
// The analysis classes after the definition of AnaEvent/AnaTrack are   //
//   just some simple tasks:                                            //
//                                                                      //
//   TaskGenerate - task to generate 10k events in a tree "T" split in  //
//                  5 files: input_n.root                               //
//   TaskFilter   - task to filter just gammas from the inputs and      //
//                  write them to a single output file                  //
//   TaskReco     - task to reconstruct pi0's from the mixed background //
//                  +signal of gammas                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TGeoMatrix.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

class AnaTrack;
class AnaEvent;

//
// This task will work without using TSelector functionality
//  No Init() or Terminate() need to be defined. No input/output
//  slots are defined; th task will simply create several output files
//
//______________________________________________________________________________
class TaskGenerate : public AliAnalysisTask {

public:
   TaskGenerate(const char *name) : AliAnalysisTask(name,"") {}
   virtual ~TaskGenerate() {}
   
   virtual void   Exec(Option_t *option);
   ClassDef(TaskGenerate,1)  // Simple generator
};

// The next task is filtering the input events coming from a chain and having
//  more than 100 gammas. In Terminate() some histogram will be drawn

//______________________________________________________________________________
class TaskFilter : public AliAnalysisTask {
private:
   AnaEvent      *fEvent;        // Current event
   TTree         *fOutput;       // Output tree
   TList         *fList;         // List containing the output data 
   TH1I          *fHist1;        // Number of gammas per event
   TH1I          *fHist2;        // Number of gammas per event
public:
   TaskFilter() : AliAnalysisTask(), fEvent(0), fOutput(0), fList(0), fHist1(0), fHist2(0) {}
   TaskFilter(const char *name);
   virtual ~TaskFilter() {}
   
   virtual void   ConnectInputData(Option_t *);
   virtual void   CreateOutputObjects();
   virtual void   Exec(Option_t *option);
   virtual void   Terminate(Option_t *);
   ClassDef(TaskFilter,1)  // Event filter
};   

// This task reconstructs pi0's for gammas coming from the same vertex   
//______________________________________________________________________________
class TaskRecoPi0 : public AliAnalysisTask {
private:
   AnaEvent      *fEvent;        // Current event
   TObjArray     *fGammas;       // Array of gammas
   TObjArray     *fPions;        // Array of pi0's
   TH1F          *fHist;         // Pt distrib. of reconstructed pions
public:
   TaskRecoPi0() : AliAnalysisTask(), fEvent(0), fGammas(0), fPions(0), fHist(0) {}
   TaskRecoPi0(const char *name);
   virtual ~TaskRecoPi0();

   virtual void   ConnectInputData(Option_t *);
   virtual void   CreateOutputObjects();
   virtual void   Exec(Option_t *option);
   virtual void   Terminate(Option_t *);
   ClassDef(TaskRecoPi0,1)  // Pi0 reconstructor
};   
      
//______________________________________________________________________________
class AnaTrack : public TObject {

private:
   Double_t       fPx;           // X component of the momentum
   Double_t       fPy;           // Y component of the momentum
   Double_t       fPz;           // Z component of the momentum
   Float_t        fMass;         // The mass of this particle
   Int_t          fCharge;       // Charge
   Double_t       fVertex[3];    // Track vertex position

public:
   AnaTrack() {}
   AnaTrack(Double_t random, Double_t *vertex=0);
   AnaTrack(Double_t px, Double_t py, Double_t pz, Float_t mass, Int_t charge, 
            Double_t vx, Double_t vy, Double_t vz) : TObject(), fPx(px), fPy(py), fPz(pz), fMass(mass), fCharge(charge)
            {fVertex[0]=vx; fVertex[1]=vy; fVertex[2]=vz;}
   virtual ~AnaTrack() {}
   
   Bool_t         Decay(Double_t &px1, Double_t &py1, Double_t &pz1, 
                        Double_t &px2, Double_t &py2, Double_t &pz2);
   Double_t       GetPx() const {return fPx;}
   Double_t       GetPy() const {return fPy;}
   Double_t       GetPz() const {return fPz;}
   void           SetPx(Double_t px) {fPx = px;}
   void           SetPy(Double_t py) {fPy = py;}
   void           SetPz(Double_t pz) {fPz = pz;}
   Double_t       GetPt() const {return TMath::Sqrt(fPx*fPx + fPy*fPy);}
   Double_t       GetP()  const {return TMath::Sqrt(fPx*fPx + fPy*fPy + fPz*fPz);}
   Double_t       GetMass() const {return fMass;}
   Double_t       GetCharge() const {return fCharge;}
   Double_t       GetVertex(Int_t i) const {return (i<3)?fVertex[i]:0.;}
   
   ClassDef(AnaTrack,1)  // A simple track
};   
 
//______________________________________________________________________________
class AnaEvent : public TObject {

private:
   Int_t          fEventNumber;  // Event number
   Int_t          fNtracks;      // Number of tracks
   TClonesArray  *fTracks;       //-> Array of all tracks

   static TClonesArray *fgTracks;

public:
   AnaEvent();
   virtual ~AnaEvent() {Clear();}
   
   AnaTrack      *AddTrack(Double_t rnd, Double_t *vert=0);
   Int_t          Build(Int_t ev);
   static void    CreateEvents(Int_t nevents, const char *filename);
   void           Clear(Option_t *option="");
   
   Int_t          GetEvtNumber() const {return fEventNumber;}
   Int_t          GetNtracks() const   {return fNtracks;}
   AnaTrack      *GetTrack(Int_t i)    {return (AnaTrack*)fTracks->At(i);}

ClassDef(AnaEvent,1)  // An event with tracks   
};

#endif
