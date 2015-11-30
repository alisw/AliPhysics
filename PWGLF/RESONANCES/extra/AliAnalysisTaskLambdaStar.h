#ifndef ALIANALYSISTASKLAMBDASTAR_H
#define ALIANALYSISTASKLAMBDASTAR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Task to study Lambda* resonance (formalism same as committed proton -lambda task of femtoscopy)
// Author: Sadhana Dash


class TH1F;
class TH2F;
class TH3F;
class TH3D;
class TList;
class TAxis;
#include <THnSparse.h>
#include <AliAnalysisTaskSE.h>
class AliTPCPIDResponse;
class AliPIDResponse;

#include <AliAODEvent.h>
#include <AliAODVertex.h>
#include <AliAODv0.h>
#include <AliAODTrack.h>

class AliAnalysisTaskLambdaStar : public AliAnalysisTaskSE {
  // Forward declaration for nested classes
 private:
  class ResoBufferTrack;
  class ResoBufferEvent;
  class ResoBuffer;
 public:
  AliAnalysisTaskLambdaStar();                 // Dummy constructor
  AliAnalysisTaskLambdaStar(const char *name); // Constructor
  virtual ~AliAnalysisTaskLambdaStar();                // Destructor
  AliAnalysisTaskLambdaStar& operator=(const AliAnalysisTaskLambdaStar& atpl); //Not implemented
  
  void   UserCreateOutputObjects();  // Called once at the beginning
  void   UserExec(Option_t *option); // Main loop
  void   Terminate(Option_t *);      // Called once at the end
  void   UseCircPID(Bool_t cpid){fCirc = cpid;}
  void   SetCentrality(Int_t centmin, Int_t centmax)
  {fCentMin = centmin;
  fCentMax = centmax;}
  void  SetNSigma(Double_t nsigma){fNSigma = nsigma;}
    void  SetRejNSigma(Double_t rejnsigma){fRejNSigma = rejnsigma;}
  void  SetNMix(Int_t nmix){fNMix = nmix;}
  void  SetCentPatch(Int_t centp) {fCentPerPatch = centp ;}
  void  SetClusterTPC(Int_t clusterTPC){fClusterTPC = clusterTPC;}
  void  SetDCAxy(Float_t dcaxy){fDCAxy = dcaxy;}
  void  SetFilterBit(Int_t filterbitbit){fFilterBit = filterbitbit;}
  
 private:


  void   ProcessTPC(AliAODTrack *track,Double_t nsig);         
  void   ProcessHybrid(AliAODTrack *track, Bool_t circ,Double_t nsig );         // Use TPC for p < 0.75 GeV/c);
  void   ProcessHybridPro(AliAODTrack *track, Bool_t circ,Double_t nsig );         // Use TPC for p < 0.75 GeV/c);
  void   ProcessReal();               //  minv calc for real events
  void   ProcessMixed();              // minv calc for mixed events
  void   ProcessLikeSignBkg();     // Check procedure with like sign
  Double_t ApplyCentralityPatchPbPb2011(AliCentrality *central);
  


  Float_t MInv(ResoBufferTrack track1,	  ResoBufferTrack track2);               // Calcs Qinv for proton proton
  Float_t Costheta(ResoBufferTrack track1,ResoBufferTrack track2);               // Calcs Qinv for proton proton
  Float_t Costheta1(ResoBufferTrack track1, ResoBufferTrack track2);               // Calcs Qinv for proton proton
  Float_t Pt(ResoBufferTrack track1, ResoBufferTrack track2);   // pt of pair
  Float_t Rapidity(ResoBufferTrack track1, ResoBufferTrack track2);   // pt of pair
  Double_t OpeningAngle(ResoBufferTrack track1, ResoBufferTrack track2);   // pt of pair
  Bool_t  goodDCA(AliAODTrack *track);                  //! Does a cut on the DCAxy value and fills a histogram
  Bool_t  goodDCAKaon(AliAODTrack *track);                  //! Does a cut on the DCAxy value and fills a histogram
  Float_t RapidityProton(AliAODTrack *track);  //! Rapidity assuming proton mass
  Float_t RapidityKaon(AliAODTrack *track);  //! Rapidity assuming proton mass


  // For AOD tracks, need a function to get the dca
  static Float_t DCAxy(const AliAODTrack *track, const AliVEvent *evt); 
  void StoreGlobalTrackReference(AliAODTrack *track);      // Store pointer to global track
  void ResetGlobalTrackReference();                        // Reset all pointers to global tracks
  Bool_t AcceptTrack(const AliAODTrack *track);            // Famous crossed rows / findable clusters
  Bool_t GoodTPCFitMapSharedMap(const AliAODTrack *track); // Rejects shared clusters, primaries
  void FillDedxHist(const AliVTrack *track);               // Fill dE/dx histograms

  const Float_t  fkAbsZvertexCut;        //! are needed for event mixing
  const Float_t  fkCentCut;                //! are needed for event mixing
  const Float_t  fkKaonMass;                //! PDG-Mass of the lambda particle
  const Float_t  fkProMass;                //! PDG-Mass of the lambda particle

  AliPIDResponse  *fPIDResponse;           //! PID response object

  Bool_t  fCirc;           //! PID response object
  Int_t fCentMin;
  Int_t fCentMax;
  Double_t fNSigma;
  Double_t fRejNSigma;
  Int_t fNMix;
  Bool_t  fPatch;           //
  Int_t  fCentPerPatch;
  ULong64_t fTriggerMask; // trigger mask
  Int_t fClusterTPC;
  Float_t fDCAxy;
  Int_t fFilterBit;

  ResoBuffer     *fResoBuffer;           //! Event mixing: event collection
  AliAODEvent     *fAOD;                   //! AOD event
  AliAODVertex    *fPrimaryVtx;            //! AOD vertex
  Double_t       fPrimaryVtxPosition[3];   //! Position of prim. vertex
  TList           *fOutputList;            //! V0 output list 
  TList           *fOutputPrimaries;       //! Primaries output list

  // Store pointers to global tracks for pid and dca
  AliAODTrack    **fGTI;                  //! Array of pointers, just nicely sorted according to the id
  const UShort_t  fTrackBuffSize;          //! Size fo the above array, ~12000 for PbPb
  //
  //! ---------        Histograms        ---------
  //
  TH1F        *fHistGoodEvent;                  //! Control hist for event cuts
  TH1F        *fHistEvent;    // event after trigger selection
  TH2D	      *fHistZVertexCent;	 //	Z coordinate of primary vertex
  TH1F        *fPriHistShare;                   //! Primaries: number of shared clusters
  TH2D        *fHistMassPtPKMin;   
  TH2D        *fHistMassPtPbarKPlus;               
  TH2D        *fHistMassPtPKMinMix;   
  TH2D        *fHistMassPtPbarKPlusMix;               
  TH2D        *fHistMassPtPKPlusLS;   
  TH2D        *fHistMassPtPbarKMinLS;               
  


  // Primary particles
 
    
  TH2F        *fPriHistTPCsignal;            //! Pri: TPC dE/dx signal vs p
  TH2F        *fPriHistTPCsignalproton;
  TH2F        *fPriHistTPCsignalkaon;
  TH3F        *fPriHistDCAxyYPtPro;             //! Pri: DCAxy vs (y,pt) for protons
  TH3F        *fPriHistDCAxyYPtAPro;            //! Pri: DCAxy vs (y,pt) for anti-protons
  TH3F        *fPriHistDCAxyYPtKPlus;             //! Pri: DCAxy vs (y,pt) for protons
  TH3F        *fPriHistDCAxyYPtKMinus;            //! Pri: DCAxy vs (y,pt) for anti-protons

    
    
    // Primary particles
    TH2D        *fPriHistTPCnsigmakaon_nocut;         // nsigma TPC vs pt for kaon
    TH2D        *fPriHistTOFnsigmakaon_nocut;         // nsigma TOF vs pt for kaon
    TH2D        *fPriHistTPCnsigmaproton_nocut;       // nsigma TPC vs pt for proton
    TH2D        *fPriHistTOFnsigmaproton_nocut;       // nsigma TOF vs pt for proton
    TH2D        *fPriHistTPCTOFnsigmakaon_nocut;      // nsigma TPC & TOFvs pt for kaon
    TH2D        *fPriHistTPCTOFnsigmaproton_nocut;    // nsigma TPC & TOFvs pt for proton
    
    
    
    // Primary particles
    TH2D        *fPriHistTPCnsigmakaon;         // nsigma TPC vs pt for kaon
    TH2D        *fPriHistTOFnsigmakaon;         // nsigma TOF vs pt for kaon
    TH2D        *fPriHistTPCnsigmaproton;       // nsigma TPC vs pt for proton
    TH2D        *fPriHistTOFnsigmaproton;       // nsigma TOF vs pt for proton
    TH2D        *fPriHistTPCTOFnsigmakaon;      // nsigma TPC & TOFvs pt for kaon
    TH2D        *fPriHistTPCTOFnsigmaproton;    // nsigma TPC & TOFvs pt for proton
    
    
    
    
    
 
 private:
  //  ------------------------------------
  //
  //    Nested classes for event mixing
  //    Features are 
  //    * low memory usage as only needed variables
  //      are stored and not, e.g., the full AliAODTrack
  //    * it's fast, as memory allocation occours
  //      only once and no deep copying is done
  //      when doing the fifo shift and reseting
  //      the event. Resetting and shifting could
  //      be reduced to only a few integer assignments.
  //      
  //  ------------------------------------


  class ResoBufferTrack{// Charged tracks
  public:
    ResoBufferTrack();                               // Constructor
    ResoBufferTrack(const AliAODTrack *track,
		     const Float_t bfield,
		     const Float_t priVtx[3]);           // Constructor
    ~ResoBufferTrack(){;}                            // Destructor, nothing to do.
    ResoBufferTrack(const ResoBufferTrack& fbt);    // Copy constructor
    ResoBufferTrack& operator=(const ResoBufferTrack& fbt); // Assignment operator
    void Set(const AliAODTrack *track,  const Float_t bfield,const Float_t priVtx[3]);    
    void Set(const AliAODTrack *track,   const Float_t bfield, const Double_t priVtx[3]);

    void SetP(const AliAODTrack *track);
    // Two functions to get the global and shifted positions. Both is not trivial as
    // tracks in ALICE get propagated in the local system
    // Shifted for studying all events shifted to (0,0,0)
    void GetShiftedPositionAtShiftedRadii(const AliAODTrack *track,
					  const Float_t bfield, const Float_t priVtx[3]); 


    Bool_t UseIt() const {return fID!=65535;}          // Use entry? eg. set to false by cleaning procedure
    void   SetBadFlag(){fID=65535;}                   // Can only set 'bad track' flag
    // Interpretation of ALICE coding rule RC14 is that ResoBufferTrack is a private member
    // of AliAnalysisTaskLambdaStar, so the following data member are private
    UShort_t fID;               //! Unique track id (->AliAODTrack.h), UShort_t goes to 65000
    Double_t fP[3];             //! Momentum of track
    Double_t fPx;             //! Momentum of track
    Double_t fPy;             //! Momentum of track
    /* Float_t  fXglobal[9][3];    //! Global positions at different global radii */
    Float_t  fXshifted[9][3];   //! Shifted positions at different shifted radii
  };
  
 
  
  class ResoBufferEvent{// Event 
  public:
    ResoBufferEvent();                       // Constructor
    ResoBufferEvent(const UShort_t priTrackBuff, const Double_t bfield,const Double_t priVtxPos[3]); // Constructor
    ResoBufferEvent(const UShort_t priTrackBuff); // Constructor
    ResoBufferEvent(const ResoBufferEvent &fbe);           // Copy constructor
    // Assignment operator won't change the size of the arrays!
    ResoBufferEvent& operator=(const ResoBufferEvent &fbe); // Assignment operator
    ~ResoBufferEvent();                       // Destructor
    void Reset(const Double_t bfield, const Double_t priVtxPos[3]);// Resets the event with new variables given
    
    // Functions to add particles to the event
    void AddPro(const AliAODTrack *track);      // Add a proton to this event  
    void AddAPro(const AliAODTrack *track);     // Add a anti-proton this actual event  

    void AddKPlus(const AliAODTrack *track);      // Add a Kplus to this event  
    void AddKMin(const AliAODTrack *track);     // Add a Kminus 

    
    // Getters for the event properties, no setters as you're not supposed to change them :)
    // The fBfield and fPriVtxPos get set on each 'reset' for a new event and the limits get set 
    // on creation of the event
    Double_t GetBfield()const{return fBfield;}            // Getter for magnetic field of event
    void GetVtxPos(Double_t xyz[3])const{for(Int_t i=0;i<3;i++)xyz[i]=fPriVtxPos[i];} // Get the xyz of the vertex.
    UShort_t GetPriTrackLim()const{return fPriTrackLim;}  // Get the size of the array for primary tracks

    // The number of tracks stored in the event
    UShort_t GetNPro()const{return fNProTracks;}       // Number of stored protons
    UShort_t GetNAPro()const{return fNAProTracks;}     // .. anti-protons
    UShort_t GetNKPlus()const{return fNKPlusTracks;}       // kplus
    UShort_t GetNKMin()const{return fNKMinTracks;}     // .. kminus

    // Data member, sorry for this private public private mixture,

  private:

    const UShort_t fPriTrackLim;   // Limit for primary tracks
   
  public:
    // Pointer to array of ...
    ResoBufferTrack *fProTracks;          //! Proton tracks
    ResoBufferTrack *fAProTracks;         //! Anti-proton tracks
    
    ResoBufferTrack *fKPlusTracks;          //! kplus tracks
    ResoBufferTrack *fKMinTracks;         //!  kminus tracks
    

  private:
    // Number of stored tracks in the event
    UShort_t fNProTracks;    // Number of stored protons
    UShort_t fNAProTracks;   // Number of stored anti-protons
    UShort_t fNKPlusTracks;    // Number of stored kaons
    UShort_t fNKMinTracks;   // Number of stored kaons


    // Double_t needed??? magnetic field probably is like 5.0 and not 5.0000120047
    Double_t fBfield;               // Magnetic field in ALICE unit [kG]
    Double_t fPriVtxPos[3];         // Primary vtx position
  };
  
  class ResoBuffer { // Holds the events
  public:
    ResoBuffer();           // Dummy constructor
    ResoBuffer(const UChar_t ZvertexBins,const UChar_t CentBins,const UChar_t MixBuff,const UShort_t PriTrackLim,const Float_t AbsZvertexCut,const Float_t CentCut); // Constructor
    ResoBuffer(const ResoBuffer &fb); //Ctor
    ResoBuffer& operator=(const AliAnalysisTaskLambdaStar::ResoBuffer&); // Assignment
    ~ResoBuffer();          // Destructor
    void ShiftAndAdd(AliAODEvent *evt); // Discard last event, shift all, set first one
    void ShiftAndAdd(const Double_t bfield,const Double_t priVtxPos[3],const Float_t centrality); // Discard last event, shift all, set first one
    ResoBufferEvent *GetEvt(const UChar_t i)const{return fCurEvt[i];}; // Returns a pointer to the i'th event of the current event mixing class
    UChar_t GetMixBuffSize()const{return fkMixBuffSize;}// Returns the number of events held in every mixing bin
                                 
  private:
    const UChar_t fkZvertexBins;           // Number of bins in Zvertex
    const UChar_t fkCentBins;              // Number of bins in centrality
    const UChar_t fkMixBuffSize;           // Number of stored events
    const UShort_t fkPriTrackLim;          // Buffer size protons per event
    const TAxis *fZvertexAxis;          //! To find Zvertex bin
    const TAxis *fCentAxis;             //! To find centrality bin

    ResoBufferEvent **fCurEvt;//! Array of pointer to the set of the current events
                                     //  Note that the pointer won't be constant
    ResoBufferEvent ****fEC;  //! The internal thing where the events
                                     //  are stored. 
                                     //  fEC stands for event collection.
  };
  
  ClassDef(AliAnalysisTaskLambdaStar, 1);
};

#endif
