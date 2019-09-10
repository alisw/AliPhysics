// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Modified version of AliAnalysisTaskCheckCascade.h
// Used bits of code from AliAnalysisTaskCheckPerformanceStrange
//
// --- Georgijs Skorodumovs
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef ALIANALYSISTASKSTRANGECASCADESDISCRETE_H
#define ALIANALYSISTASKSTRANGECASCADESDISCRETE_H

class AliRunningCascadeTrack : public TObject//public TObject
{
private:
    // Track properties
    TLorentzVector TLV_pos; // TLorentzVector of pion from Lambda decay
    TLorentzVector TLV_neg; // TLorentzVector of proton from Lambda decay
    TLorentzVector TLV_bach; // TLorentzVector bachelor track, either kaon or pion
    Short_t dca_pos_to_prim[2];  // xy, z, range is 0..180 cm
    Short_t dca_neg_to_prim[2]; // xy, z, range is 0..180 cm
    Short_t dca_bach_to_prim[2]; // xy, z, range is 0..180 cm
    Short_t dca_V0_to_prim;  // 0..120
    Short_t dca_Omega_to_prim[2]; // xy, z; 0..40, pm40
    Short_t dca_pos_to_neg; //0..1.5 cm
    Short_t dca_bach_to_Lambda; //0...2.0 cm
    Short_t dca_bach_to_pos; // 0..500.0 cm
    
    Short_t nSigma_dEdx_pos[2]; //pos is either [0]=proton(if Omega-) or [1]=pion(if Omega+)
    Short_t nSigma_dEdx_neg[2]; //neg is either [0]=pion(Omega-) or [1]=proton
    Short_t nSigma_dEdx_bach[2]; //bach is either [0]=kaon(Omega) or [1]=pion(Xi)
    
    Short_t nSigma_TOF_pos[2];
    Short_t nSigma_TOF_neg[2];
    Short_t nSigma_TOF_bach[2];
    
    Short_t CosPointingAngle; //also stores the charge of the bachelor track! (cos(THeta)*e_bach)
    Int_t   LeastNumberOfTPCclus;
    Short_t MaxChi2perclus;
    Short_t MinTrackLength;
    Float_t CascadeDecayPos[3];
    Float_t V0fromCascadePos[3];
    
    TClonesArray* defaultarray;
    
public:
    //if use declarations only -> results in "undefined reference"
    /*   AliRunningCascadeTrack(); //default constructor
     AliRunningCascadeTrack(const char* name);
     AliRunningCascadeTrack(const AliRunningCascadeTrack&); //copy operator
     AliRunningCascadeTrack& operator=(const AliRunningCascadeTrack&); //copy assignment
     virtual ~AliRunningCascadeTrack(); //default destructor*/
    
    AliRunningCascadeTrack():
    TLV_pos(),
    TLV_neg(),
    TLV_bach(),
    dca_V0_to_prim(-500),
    dca_pos_to_neg(-500),
    dca_bach_to_Lambda(-500),
    dca_bach_to_pos(-500),
    CosPointingAngle(-2),
    LeastNumberOfTPCclus(-1),
    MaxChi2perclus(-1),
    MinTrackLength(-1),
    defaultarray(0)
    {
        //  defaultarray = new TClonesArray("default array", 10);
    }
    
    AliRunningCascadeTrack(const char* name):
    TLV_pos(),
    TLV_neg(),
    TLV_bach(),
    dca_V0_to_prim(-500),
    dca_pos_to_neg(-500),
    dca_bach_to_Lambda(-500),
    dca_bach_to_pos(-500),
    CosPointingAngle(-2),
    LeastNumberOfTPCclus(-1),
    MaxChi2perclus(-1),
    MinTrackLength(-1),
    defaultarray(0)
    {
        //  defaultarray = new TClonesArray("default array", 10);
    }
    
    AliRunningCascadeTrack(const AliRunningCascadeTrack&): //copy constructor
    // TObject(),
    TObject(),
    TLV_pos(),
    TLV_neg(),
    TLV_bach(),
    dca_V0_to_prim(-500),
    dca_pos_to_neg(-500),
    dca_bach_to_Lambda(-500),
    dca_bach_to_pos(-500),
    CosPointingAngle(-2),
    LeastNumberOfTPCclus(-1),
    MaxChi2perclus(-1),
    MinTrackLength(-1),
    defaultarray(0)
    {
    }
    
    AliRunningCascadeTrack& operator=(const AliRunningCascadeTrack&)
    {
        return *this;
        
    }
    
    virtual ~AliRunningCascadeTrack()
    {
        //   delete[] defaultarray;
        // defaultarray = NULL;
    }
    
    // setters
    void set_TLV_pos(TLorentzVector tlv)   { TLV_pos   = tlv; }
    void set_TLV_neg(TLorentzVector tlv)   { TLV_neg   = tlv; }
    void set_TLV_bach(TLorentzVector tlv)   { TLV_bach   = tlv; }
    
    void set_dca_pos_to_prim(Float_t f1, Float_t f2)   { dca_pos_to_prim[0]   = (Short_t)(f1*100.0); dca_pos_to_prim[1]   = (Short_t)(f2*100.0); }
    void set_dca_neg_to_prim(Float_t f1, Float_t f2)   { dca_neg_to_prim[0]   = (Short_t)(f1*100.0); dca_neg_to_prim[1]   = (Short_t)(f2*100.0); }
    void set_dca_bach_to_prim(Float_t f1, Float_t f2)   { dca_bach_to_prim[0]   = (Short_t)(f1*100.0); dca_bach_to_prim[1]   = (Short_t)(f2*100.0); }
    void set_dca_V0_to_prim(Float_t f)  { dca_V0_to_prim  = (Short_t)(f*100.0); }
    void set_dca_Omega_to_prim(Float_t f1, Float_t f2)   { dca_Omega_to_prim[0]   = (Short_t)(f1*100.0); dca_Omega_to_prim[1]   = (Short_t)(f2*100.0); }
    void set_dca_pos_to_neg(Float_t f)  { dca_pos_to_neg = (Short_t)(f*100.0); }
    void set_dca_bach_to_Lambda(Float_t f)  { dca_bach_to_Lambda = (Short_t)(f*100.0); }
    void set_dca_bach_to_pos(Float_t f)  { dca_bach_to_pos = (Short_t)(f*10.0); }
    void set_nSigma_dEdx_pos(Float_t f1, Float_t f2)   { nSigma_dEdx_pos[0]   = (Short_t)(f1*10); nSigma_dEdx_pos[1]   = (Short_t)(f2*10); }
    void set_nSigma_dEdx_neg(Float_t f1, Float_t f2)   { nSigma_dEdx_neg[0]   = (Short_t)(f1*10); nSigma_dEdx_neg[1]   = (Short_t)(f2*10); }
    void set_nSigma_dEdx_bach(Float_t f1, Float_t f2)   { nSigma_dEdx_bach[0]   = (Short_t)(f1*10); nSigma_dEdx_bach[1]   = (Short_t)(f2*10); }
    void set_nSigma_TOF_pos(Float_t f1, Float_t f2)   { nSigma_TOF_pos[0]   = (Short_t)(f1*10); nSigma_TOF_pos[1]   = (Short_t)(f2*10); }
    void set_nSigma_TOF_neg(Float_t f1, Float_t f2)   { nSigma_TOF_neg[0]   = (Short_t)(f1*10); nSigma_TOF_neg[1]   = (Short_t)(f2*10); }
    void set_nSigma_TOF_bach(Float_t f1, Float_t f2)   { nSigma_TOF_bach[0]   = (Short_t)(f1*10); nSigma_TOF_bach[1]   = (Short_t)(f2*10); }
    void set_CosPointingAngle(Float_t f)    {CosPointingAngle = (Short_t)(f*1000.);} //pay attention to this!
    void set_LeastNumberOfTPCclus(Int_t i)      {LeastNumberOfTPCclus = i; }
    void set_MaxChi2perclus(Float_t f)   {MaxChi2perclus = (Short_t)(f*100.); }
    void set_MinTrackLength(Float_t f)    {MinTrackLength = (Short_t)(f*100.); }
    void set_CascadeDecayPos(Float_t f1, Float_t f2, Float_t f3)
    {
        CascadeDecayPos[0] = f1;
        CascadeDecayPos[1] = f2;
        CascadeDecayPos[2] = f3;
    }
    void set_V0fromCascadePos(Float_t f1, Float_t f2, Float_t f3)
    {
        V0fromCascadePos[0] = f1;
        V0fromCascadePos[1] = f2;
        V0fromCascadePos[2] = f3;
    }
    
    // getters
    TLorentzVector get_TLV_pos() const     { return TLV_pos; }
    TLorentzVector get_TLV_neg() const     { return TLV_neg; }
    TLorentzVector get_TLV_bach() const     { return TLV_bach; }
    
    Float_t get_dca_pos_to_prim(Int_t i) const   { return ((Float_t)dca_pos_to_prim[i])/100.0; }
    Float_t get_dca_neg_to_prim(Int_t i) const   { return ((Float_t)dca_neg_to_prim[i])/100.0; }
    Float_t get_dca_bach_to_prim(Int_t i)  const { return ((Float_t)dca_bach_to_prim[i])/100.0; }
    Float_t get_dca_V0_to_prim() const  { return ((Float_t)dca_V0_to_prim)/100.0; }
    Float_t get_dca_Omega_to_prim(Int_t i) const  { return ((Float_t)dca_pos_to_prim[i])/100.0; }
    Float_t get_dca_pos_to_neg() const { return ((Float_t)dca_pos_to_neg)/100.0; }
    Float_t get_dca_bach_to_Lambda()  const { return ((Float_t)dca_bach_to_Lambda)/100.0; }
    Float_t get_dca_bach_to_pos() const { return ((Float_t)dca_bach_to_pos)/10.0; }
    Float_t get_nSigma_dEdx_pos(Int_t i)  const { return ((Float_t)nSigma_dEdx_pos[i])/10.; }
    Float_t get_nSigma_dEdx_neg(Int_t i)  const {return ((Float_t)nSigma_dEdx_neg[i])/10.;}
    Float_t get_nSigma_dEdx_bach(Int_t i)const   { return ((Float_t)nSigma_dEdx_bach[i])/10.; }
    Float_t get_nSigma_TOF_pos(Int_t i) const  { return ((Float_t)nSigma_TOF_pos[i])/10.;}
    Float_t get_nSigma_TOF_neg(Int_t i) const  { return ((Float_t)nSigma_TOF_neg[i])/10.;}
    Float_t get_nSigma_TOF_bach(Int_t i) const  {return ((Float_t)nSigma_TOF_bach[i])/10.;}
    
    Float_t get_CosPointingAngle()  const  {return ((Float_t)CosPointingAngle)/1000.;} //pay attention to cos! should not be negative
    Int_t get_LeastNumberOfTPCclus()   const   {return LeastNumberOfTPCclus; }
    Float_t get_MaxChi2perclus() const  {return ((Float_t)MaxChi2perclus)/100.; }
    Float_t get_MinTrackLength() const   {return ((Float_t)MinTrackLength)/100.; }
    
    Float_t  get_CascadeDecayPos(Int_t i) const {return CascadeDecayPos[i];}
    Float_t  get_V0fromCascadePos(Int_t i) const {return V0fromCascadePos[i];}
    
    //   ClassDef(AliRunningCascadeTrack,1);  // A simple track of a particle
};




class AliRunningCascadeEvent : public TObject //public TObject //public TObject   //public AliRunningCascadeTrack
{
private:
    Float_t x;
    Float_t y;
    Float_t z;
    Int_t     id; // Run id
    Int_t     N_tracks; // total number of tracks
    Float_t   centrality; //
    Bool_t    MVPPileUpFlag;
    Int_t     multiplicity;
    Long64_t  trigger_word;
    Short_t   magfield; //magnetic field
    UShort_t      fNumTracks; // number of tracks in event
    TClonesArray* fTracks;      //->
    
public:
    
    /*  AliRunningCascadeEvent();
     AliRunningCascadeEvent(const char *name);
     virtual ~AliRunningCascadeEvent();
     AliRunningCascadeEvent(const AliRunningCascadeEvent&);
     AliRunningCascadeEvent& operator =(const AliRunningCascadeEvent&) //copy assignment
     {
     return *this;
     }*/
    
    AliRunningCascadeEvent() :
    x(-1),y(-1),z(-1),id(-1),N_tracks(0),centrality(0),
    MVPPileUpFlag(0),multiplicity(0),trigger_word(0),magfield(0),fNumTracks(0), fTracks(0x0)
    {
        fTracks = new TClonesArray( "AliRunningCascadeTrack", 10 );
    }
    
    
    AliRunningCascadeEvent(const char* name) :
    x(-1),y(-1),z(-1),id(-1),N_tracks(0),centrality(0),
    MVPPileUpFlag(0),multiplicity(0),trigger_word(0),magfield(0),fNumTracks(0), fTracks(0x0)
    {
        fTracks      = new TClonesArray( "AliRunningCascadeTrack", 10 );
    }
    
    
    AliRunningCascadeEvent(const AliRunningCascadeEvent&): TObject(), // copy constructor
    x(-1),y(-1),z(-1),id(-1),N_tracks(0),centrality(0),
    MVPPileUpFlag(0),multiplicity(0),trigger_word(0),magfield(0),fNumTracks(0), fTracks(0x0)
    {
        
    }
    
    AliRunningCascadeEvent& operator =(const AliRunningCascadeEvent&) //copy assignment
    {
        return *this;
    }
    
    virtual ~AliRunningCascadeEvent()
    {
        delete fTracks;
        fTracks = NULL;
    }
    
    
    //----setters and getters-------------------------------------------------------------------
    void       setx(Float_t r)                    { x = r;                         }
    Float_t    getx() const                       { return x;                      }
    
    void       sety(Float_t r)                    { y = r;                         }
    Float_t    gety() const                       { return y;                      }
    
    void       setz(Float_t r)                    { z = r;                         }
    Float_t    getz() const                       { return z;                      }
    
    void       setid(Int_t  r)                    { id = r;                        }
    Int_t      getid() const                      { return id;                     }
    
    void       setN_tracks(Int_t r)                 { N_tracks = r;                    }
    Int_t      getN_tracks() const                    { return N_tracks;                 }
    
    void       setcentrality(Float_t r)             {centrality  = r;                }
    Float_t    getcentrality() const              { return centrality;             }
    
    void       setMVPPileUpFlag(Bool_t r)             {MVPPileUpFlag  = r;                }
    Bool_t      getMVPPileUpFlag() const              { return MVPPileUpFlag;             }
    
    void       setmultiplicity(Int_t r)             {multiplicity  = r;                }
    Int_t      getmultiplicity() const              { return multiplicity;             }
    
    void       settrigger_word(Long64_t r)             { trigger_word = r;                }
    Long64_t      gettrigger_word() const              { return trigger_word;             }
    
    void       setmagfield(Short_t r)             { magfield = r;                }
    Long64_t      getmagfield() const              { return magfield;             }
    
    UShort_t getNumTracks() const        {return fNumTracks;}
    //---------------------------------------------------------------------------------------------------------------------
    
    AliRunningCascadeTrack* createTrack()
    {
        if (fNumTracks == fTracks->GetSize())
            fTracks->Expand( fNumTracks + 10 );
        if (fNumTracks >= 10000)
        {
            Fatal( "AliRunningCascadeEvent::createTrack()", "ERROR: Too many tracks (>10000)!" );
            exit( 2 );
        }
        
        //  new ((*fTracks)[fNumTracks++]) AliRunningCascadeTrack; //original by Alex
        //  return (AliRunningCascadeTrack*)((*fTracks)[fNumTracks - 1]); //original by Alex
        AliRunningCascadeTrack* track = new ((*fTracks)[fNumTracks++]) AliRunningCascadeTrack;
        return track;
    }
    
    void ClearTrackList()
    {
        fNumTracks   = 0;
        fTracks      ->Clear();
    }
    
    AliRunningCascadeTrack* getTrack(UShort_t i) const
    {
        return i < fNumTracks ? (AliRunningCascadeTrack*)((*fTracks)[i]) : NULL;
    }
    
    //  ClassDef(AliRunningCascadeEvent,1);  // A simple event compiled of tracks
};



class AliAnalysisTaskStrangeCascadesDiscrete : public AliAnalysisTaskSE {
private:
    
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;      //! List of Cascade histograms
    TList  *fListK0Short;        // List of Cascade histograms
    TList  *fListLambda;        // List of Cascade histograms
    TList  *fListAntiLambda;        // List of Cascade histograms
    TList  *fListXiMinus;   // List of XiMinus outputs
    TList  *fListXiPlus;   // List of XiPlus outputs
    TList  *fListOmegaMinus;   // List of XiMinus outputs
    TList  *fListOmegaPlus;   // List of XiPlus outputs
    TList *fEventList;
    TTree  *fTreeEvent;              //! Output Tree, Events
    // TTree  *fTreeV0;              //! Output Tree, V0s
    TTree  *fTreeCascade;              //! Output Tree, Cascades
    // TTree *fTreeV0Cascade;       //! Combined output
    
    //--------------------------------------------------------------------------------------
    AliRunningCascadeEvent* Cascade_Event;
    AliRunningCascadeTrack* Cascade_Track;
    
    TTree *fTreeCascadeAsEvent;     //!
    //   AliRunningCascadeEvent* Cascade_Event;
    //  AliRunningCascadeTrack* Cascade_Track;
    
    //private variables used for tests of programm running (no physical importance)-------
    // TTree *fNumberEventTree;  //! a test tree
    
    //----------------------------------------
    
    AliPIDResponse *fPIDResponse;     //! PID response object
    AliESDtrackCuts *fESDtrackCuts;   //! ESD track cuts used for primary track definition
    AliESDtrackCuts *fESDtrackCutsITSsa2010;  //! ESD track cuts used for ITSsa track definition
    AliESDtrackCuts *fESDtrackCutsGlobal2015; //! ESD track cuts used for global track definition
    AliAnalysisUtils *fUtils;         //! analysis utils (for MV pileup selection)
    Int_t number;
    
    AliEventCuts fEventCuts;                 /// Event cuts class
    AliEventCuts fEventCutsStrictAntipileup; /// Event cuts class
    
    TRandom3 *fRand; //!
    
    //Objects Controlling Task Behaviour
    Bool_t fkSaveEventTree;           //if true, save Event TTree
    Bool_t fkSaveV0Tree;              //if true, save TTree
    Bool_t fkDownScaleV0;
    Double_t fDownScaleFactorV0;
    Bool_t fkPreselectDedx;
    Bool_t fkUseOnTheFlyV0Cascading;
    Bool_t fkDebugWrongPIDForTracking; //if true, add extra information to TTrees for debugging
    Bool_t fkDebugBump; //if true, add extra information to TTrees for debugging
    Bool_t fkDebugOOBPileup; // if true, add extra information to TTrees for pileup study
    Bool_t fkDoExtraEvSels; //if true, rely on AliEventCuts
    Int_t fkPileupRejectionMode; //pileup rejection mode (0=none, 1=ionut, 2=anti-ionut)
    Bool_t fkUseOldCentrality; //if true, use AliCentrality instead of AliMultSelection
    
    Bool_t fkSaveCascadeTree;         //if true, save TTree
    Bool_t fkDownScaleCascade;
    Double_t fDownScaleFactorCascade;
    
    //Objects Controlling Task Behaviour: has to be streamed!
    Bool_t    fkRunVertexers;           // if true, re-run vertexer with loose cuts *** only for CASCADES! ***
    Bool_t    fkUseLightVertexer;       // if true, use AliLightVertexers instead of regular ones
    Bool_t    fkDoV0Refit;              // if true, will invoke AliESDv0::Refit in the vertexing procedure
    Bool_t    fkExtraCleanup;           //if true, perform pre-rejection of useless candidates before going through configs
    
    AliVEvent::EOfflineTriggerTypes fTrigType; // trigger type
    
    Double_t  fV0VertexerSels[7];        // Array to store the 7 values for the different selections V0 related
    Double_t  fCascadeVertexerSels[8];   // Array to store the 8 values for the different selections Casc. related
    
    Double_t fLambdaMassMean[5]; //Array to store the lambda mass mean parametrization
    //[0]+[1]*TMath::Exp([2]*x)+[3]*TMath::Exp([4]*x)
    
    Double_t fLambdaMassSigma[4]; //Array to store the lambda mass sigma parametrization
    //[0]+[1]*x+[2]*TMath::Exp([3]*x)
    
    Float_t fMinPtToSave; //minimum pt above which we keep candidates in TTree output
    Float_t fMaxPtToSave; //maximum pt below which we keep candidates in TTree output
    
    //if true, save sandbox mode info (beware large files!)
    Bool_t fkSandboxMode;
    
    //===========================================================================================
    //   Variables for Event Tree
    //===========================================================================================
    Float_t fCentrality; //!
    Bool_t fMVPileupFlag; //!
    Bool_t fOOBPileupFlag; //!
    
    //TOF info for OOB pileuo study
    Int_t  fNTOFClusters;  //!
    Int_t  fNTOFMatches;   //!
    Int_t  fNTracksITSsa2010; //!
    Int_t  fNTracksGlobal2015; //!
    Int_t  fNTracksGlobal2015TriggerPP; //!
    
    //V0 info for OOB pileup study
    Float_t fAmplitudeV0A; //!
    Float_t fAmplitudeV0C; //!
    
    //IR info for OOB pileup study
    Int_t fClosestNonEmptyBC; //!
    
    //===========================================================================================
    //   Variables for V0 Tree
    //===========================================================================================
    Float_t fTreeVariableChi2V0;         //!
    Float_t fTreeVariableDcaV0Daughters; //!
    Float_t fTreeVariableDcaV0ToPrimVertex; //!
    Float_t fTreeVariableDcaPosToPrimVertex; //!
    Float_t fTreeVariableDcaNegToPrimVertex; //!
    Float_t fTreeVariableV0CosineOfPointingAngle; //!
    Float_t fTreeVariableV0Radius; //!
    Float_t fTreeVariablePt; //!
    
    Float_t fTreeVariableRapK0Short; //!
    Float_t fTreeVariableRapLambda; //!
    Float_t fTreeVariableInvMassK0s; //!
    Float_t fTreeVariableInvMassLambda; //!
    Float_t fTreeVariableInvMassAntiLambda; //!
    Float_t fTreeVariableAlphaV0; //!
    Float_t fTreeVariablePtArmV0;//!
    Float_t fTreeVariableNegEta; //!
    Float_t fTreeVariablePosEta; //!
    
    Float_t fTreeVariableNSigmasPosProton; //!
    Float_t fTreeVariableNSigmasPosPion; //!
    Float_t fTreeVariableNSigmasNegProton; //!
    Float_t fTreeVariableNSigmasNegPion; //!
    
    Float_t fTreeVariableDistOverTotMom;//!
    Int_t   fTreeVariableLeastNbrCrossedRows;//!
    Float_t fTreeVariableLeastRatioCrossedRowsOverFindable;//!
    Float_t fTreeVariableMaxChi2PerCluster; //!
    Float_t fTreeVariableMinTrackLength; //!
    
    //Variables for debugging Wrong PID hypothesis in tracking bug
    // more info at: https://alice.its.cern.ch/jira/browse/PWGPP-218
    Int_t fTreeVariablePosPIDForTracking; //! uses AliPID EParticleType code (0=electron, 1=muon, 2=pion, etc)
    Int_t fTreeVariableNegPIDForTracking; //!
    Float_t fTreeVariablePosdEdx; //!
    Float_t fTreeVariableNegdEdx; //!
    Float_t fTreeVariablePosInnerP; //!
    Float_t fTreeVariableNegInnerP; //!
    //Decay Length issue debugging: ULong_t with track status
    ULong64_t fTreeVariableNegTrackStatus; //!
    ULong64_t fTreeVariablePosTrackStatus; //!
    Float_t fTreeVariableNegDCAz; //!
    Float_t fTreeVariablePosDCAz; //!
    
    //Cluster information for all daughter tracks
    Bool_t fTreeVariablePosITSClusters0;
    Bool_t fTreeVariablePosITSClusters1;
    Bool_t fTreeVariablePosITSClusters2;
    Bool_t fTreeVariablePosITSClusters3;
    Bool_t fTreeVariablePosITSClusters4;
    Bool_t fTreeVariablePosITSClusters5;
    
    Bool_t fTreeVariableNegITSClusters0;
    Bool_t fTreeVariableNegITSClusters1;
    Bool_t fTreeVariableNegITSClusters2;
    Bool_t fTreeVariableNegITSClusters3;
    Bool_t fTreeVariableNegITSClusters4;
    Bool_t fTreeVariableNegITSClusters5;
    
    //Cluster information for all daughter tracks
    Bool_t fTreeVariablePosITSSharedClusters0;
    Bool_t fTreeVariablePosITSSharedClusters1;
    Bool_t fTreeVariablePosITSSharedClusters2;
    Bool_t fTreeVariablePosITSSharedClusters3;
    Bool_t fTreeVariablePosITSSharedClusters4;
    Bool_t fTreeVariablePosITSSharedClusters5;
    
    Bool_t fTreeVariableNegITSSharedClusters0;
    Bool_t fTreeVariableNegITSSharedClusters1;
    Bool_t fTreeVariableNegITSSharedClusters2;
    Bool_t fTreeVariableNegITSSharedClusters3;
    Bool_t fTreeVariableNegITSSharedClusters4;
    Bool_t fTreeVariableNegITSSharedClusters5;
    
    Bool_t fTreeVariableIsCowboy; //store if V0 is cowboy-like or sailor-like in XY plane
    
    //Variables for OOB pileup study (high-multiplicity triggers pp 13 TeV - 2016 data)
    Float_t fTreeVariableNegTOFExpTDiff; //!
    Float_t fTreeVariablePosTOFExpTDiff; //!
    Float_t fTreeVariableNegTOFSignal; //!
    Float_t fTreeVariablePosTOFSignal; //!
    Int_t   fTreeVariableNegTOFBCid; //!
    Int_t   fTreeVariablePosTOFBCid; //!
    //Event info
    Float_t fTreeVariableAmplitudeV0A; //!
    Float_t fTreeVariableAmplitudeV0C; //!
    Int_t   fTreeVariableClosestNonEmptyBC; //!
    
    //Event Multiplicity Variables
    Float_t fTreeVariableCentrality; //!
    Bool_t fTreeVariableMVPileupFlag; //!
    Bool_t fTreeVariableOOBPileupFlag; //!
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //Sandbox V0
    Float_t fTreeVariablePrimVertexX;
    Float_t fTreeVariablePrimVertexY;
    Float_t fTreeVariablePrimVertexZ;
    
    AliExternalTrackParam *fTreeVariablePosTrack; //!
    AliExternalTrackParam *fTreeVariableNegTrack; //!
    
    Float_t fTreeVariableMagneticField;
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    
    //===========================================================================================
    //   Variables for Cascade Candidate Tree
    //===========================================================================================
    Int_t fTreeCascVarCharge;         //!
    Float_t fTreeCascVarMassAsXi;     //!
    Float_t fTreeCascVarMassAsOmega;  //!
    Float_t fTreeCascVarPt;           //!
    Float_t fTreeCascVarP;           //!
    Float_t fTreeCascVarPz;          //!
    Float_t fTreeCascVarRapXi;        //!
    Float_t fTreeCascVarRapOmega;     //!
    Float_t fTreeCascVarNegEta;       //!
    Float_t fTreeCascVarPosEta;       //!
    Float_t fTreeCascVarBachEta;      //!
    Float_t fTreeCascVarDCACascDaughters; //!
    Float_t fTreeCascVarDCABachToPrimVtx; //!
    Float_t fTreeCascVarDCABachToPrimVtxZ; //!
    Float_t fTreeCascVarDCAV0Daughters;   //!
    Float_t fTreeCascVarDCAV0ToPrimVtx;   //!
    Float_t fTreeCascVarDCAPosToPrimVtx;  //!
    Float_t fTreeCascVarDCAPosToPrimVtxZ; //!
    Float_t fTreeCascVarDCANegToPrimVtx;  //!
    Float_t fTreeCascVarDCANegToPrimVtxZ;  //!
    Float_t fTreeCascVarCascCosPointingAngle;         //!
    Float_t fTreeCascVarCascDCAtoPVxy;         //!
    Float_t fTreeCascVarCascDCAtoPVz;         //!
    Float_t fTreeCascVarCascRadius;                   //!
    Float_t fTreeCascVarV0Mass;                       //!
    Float_t fTreeCascVarV0MassLambda;                       //!
    Float_t fTreeCascVarV0MassAntiLambda;                       //!
    Float_t fTreeCascVarV0CosPointingAngle;           //!
    Float_t fTreeCascVarV0CosPointingAngleSpecial;    //!
    Float_t fTreeCascVarV0Radius;                     //!
    Float_t fTreeCascVarDCABachToBaryon;                     //!
    Float_t fTreeCascVarWrongCosPA;                   //!
    Int_t   fTreeCascVarLeastNbrClusters;             //!
    Float_t fTreeCascVarDistOverTotMom;               //!
    Float_t fTreeCascVarMaxChi2PerCluster; //!
    Float_t fTreeCascVarMinTrackLength; //!
    
    //TPC dEdx
    Float_t fTreeCascVarNegNSigmaPion;   //!
    Float_t fTreeCascVarNegNSigmaProton; //!
    Float_t fTreeCascVarPosNSigmaPion;   //!
    Float_t fTreeCascVarPosNSigmaProton; //!
    Float_t fTreeCascVarBachNSigmaPion;  //!
    Float_t fTreeCascVarBachNSigmaKaon;  //!
    
    //TOF (experimental, not corrected for weak decay traj)
    Float_t fTreeCascVarNegTOFNSigmaPion;   //!
    Float_t fTreeCascVarNegTOFNSigmaProton; //!
    Float_t fTreeCascVarPosTOFNSigmaPion;   //!
    Float_t fTreeCascVarPosTOFNSigmaProton; //!
    Float_t fTreeCascVarBachTOFNSigmaPion;  //!
    Float_t fTreeCascVarBachTOFNSigmaKaon;  //!
    
    Float_t fTreeCascVarNegITSNSigmaPion;   //!
    Float_t fTreeCascVarNegITSNSigmaProton; //!
    Float_t fTreeCascVarPosITSNSigmaPion;   //!
    Float_t fTreeCascVarPosITSNSigmaProton; //!
    Float_t fTreeCascVarBachITSNSigmaPion;  //!
    Float_t fTreeCascVarBachITSNSigmaKaon;  //!
    
    //ChiSquares
    Float_t fTreeCascVarChiSquareV0;
    Float_t fTreeCascVarChiSquareCascade;
    
    //Extended information: uncertainties at point closest to pV
    Float_t fTreeCascVarBachDCAPVSigmaX2; //
    Float_t fTreeCascVarBachDCAPVSigmaY2; //
    Float_t fTreeCascVarBachDCAPVSigmaZ2; //
    Float_t fTreeCascVarPosDCAPVSigmaX2; //
    Float_t fTreeCascVarPosDCAPVSigmaY2; //
    Float_t fTreeCascVarPosDCAPVSigmaZ2; //
    Float_t fTreeCascVarNegDCAPVSigmaX2; //
    Float_t fTreeCascVarNegDCAPVSigmaY2; //
    Float_t fTreeCascVarNegDCAPVSigmaZ2; //
    
    //Variables for debugging Wrong PID hypothesis in tracking bug
    // more info at: https://alice.its.cern.ch/jira/browse/PWGPP-218
    Int_t fTreeCascVarPosPIDForTracking; //! uses AliPID EParticleType code (0=electron, 1=muon, 2=pion, etc)
    Int_t fTreeCascVarNegPIDForTracking; //!
    Int_t fTreeCascVarBachPIDForTracking; //!
    Float_t fTreeCascVarNegInnerP; //!
    Float_t fTreeCascVarPosInnerP; //!
    Float_t fTreeCascVarBachInnerP; //!
    Float_t fTreeCascVarNegdEdx; //!
    Float_t fTreeCascVarPosdEdx; //!
    Float_t fTreeCascVarBachdEdx; //!
    
    //Decay Length issue debugging: ULong_t with track status
    ULong64_t fTreeCascVarNegTrackStatus; //!
    ULong64_t fTreeCascVarPosTrackStatus; //!
    ULong64_t fTreeCascVarBachTrackStatus; //!
    
    Float_t fTreeCascVarNegDCAz; //!
    Float_t fTreeCascVarPosDCAz; //!
    Float_t fTreeCascVarBachDCAz; //!
    
    //Variables for debugging the invariant mass bump
    //Full momentum information
    Float_t fTreeCascVarNegPx; //!
    Float_t fTreeCascVarNegPy; //!
    Float_t fTreeCascVarNegPz; //!
    Float_t fTreeCascVarPosPx; //!
    Float_t fTreeCascVarPosPy; //!
    Float_t fTreeCascVarPosPz; //!
    Float_t fTreeCascVarBachPx; //!
    Float_t fTreeCascVarBachPy; //!
    Float_t fTreeCascVarBachPz; //!
    
    Float_t fTreeCascVarV0DecayX; //!
    Float_t fTreeCascVarV0DecayY; //!
    Float_t fTreeCascVarV0DecayZ; //!
    Float_t fTreeCascVarCascadeDecayX; //!
    Float_t fTreeCascVarCascadeDecayY; //!
    Float_t fTreeCascVarCascadeDecayZ; //!
    Float_t fTreeCascVarV0Lifetime; //! //V0 lifetime (actually, mL/p)
    //Track Labels (check for duplicates, etc)
    Int_t fTreeCascVarNegIndex; //!
    Int_t fTreeCascVarPosIndex; //!
    Int_t fTreeCascVarBachIndex; //!
    //Event Number (check same-event index mixups)
    ULong64_t fTreeCascVarEventNumber; //!
    
    //Cluster information for all daughter tracks
    Bool_t fTreeCascVarPosITSClusters0;
    Bool_t fTreeCascVarPosITSClusters1;
    Bool_t fTreeCascVarPosITSClusters2;
    Bool_t fTreeCascVarPosITSClusters3;
    Bool_t fTreeCascVarPosITSClusters4;
    Bool_t fTreeCascVarPosITSClusters5;
    
    Bool_t fTreeCascVarNegITSClusters0;
    Bool_t fTreeCascVarNegITSClusters1;
    Bool_t fTreeCascVarNegITSClusters2;
    Bool_t fTreeCascVarNegITSClusters3;
    Bool_t fTreeCascVarNegITSClusters4;
    Bool_t fTreeCascVarNegITSClusters5;
    
    Bool_t fTreeCascVarBachITSClusters0;
    Bool_t fTreeCascVarBachITSClusters1;
    Bool_t fTreeCascVarBachITSClusters2;
    Bool_t fTreeCascVarBachITSClusters3;
    Bool_t fTreeCascVarBachITSClusters4;
    Bool_t fTreeCascVarBachITSClusters5;
    
    //Cluster information for all daughter tracks
    Bool_t fTreeCascVarPosITSSharedClusters0;
    Bool_t fTreeCascVarPosITSSharedClusters1;
    Bool_t fTreeCascVarPosITSSharedClusters2;
    Bool_t fTreeCascVarPosITSSharedClusters3;
    Bool_t fTreeCascVarPosITSSharedClusters4;
    Bool_t fTreeCascVarPosITSSharedClusters5;
    
    Bool_t fTreeCascVarNegITSSharedClusters0;
    Bool_t fTreeCascVarNegITSSharedClusters1;
    Bool_t fTreeCascVarNegITSSharedClusters2;
    Bool_t fTreeCascVarNegITSSharedClusters3;
    Bool_t fTreeCascVarNegITSSharedClusters4;
    Bool_t fTreeCascVarNegITSSharedClusters5;
    
    Bool_t fTreeCascVarBachITSSharedClusters0;
    Bool_t fTreeCascVarBachITSSharedClusters1;
    Bool_t fTreeCascVarBachITSSharedClusters2;
    Bool_t fTreeCascVarBachITSSharedClusters3;
    Bool_t fTreeCascVarBachITSSharedClusters4;
    Bool_t fTreeCascVarBachITSSharedClusters5;
    
    //Variables for OOB pileup study (high-multiplicity triggers pp 13 TeV - 2016 data)
    Float_t fTreeCascVarNegTOFExpTDiff; //!
    Float_t fTreeCascVarPosTOFExpTDiff; //!
    Float_t fTreeCascVarBachTOFExpTDiff; //!
    Float_t fTreeCascVarNegTOFSignal; //!
    Float_t fTreeCascVarPosTOFSignal; //!
    Float_t fTreeCascVarBachTOFSignal; //!
    Int_t   fTreeCascVarNegTOFBCid; //!
    Int_t   fTreeCascVarPosTOFBCid; //!
    Int_t   fTreeCascVarBachTOFBCid; //!
    //Event info
    Float_t fTreeCascVarAmplitudeV0A; //!
    Float_t fTreeCascVarAmplitudeV0C; //!
    Int_t   fTreeCascVarClosestNonEmptyBC; //!
    
    //Event Multiplicity Variables
    Float_t fTreeCascVarCentrality; //!
    Bool_t fTreeCascVarMVPileupFlag; //!
    Bool_t fTreeCascVarOOBPileupFlag; //!
    
    //Kink tagging
    Bool_t fTreeCascVarBachIsKink;
    Bool_t fTreeCascVarPosIsKink;
    Bool_t fTreeCascVarNegIsKink;
    
    //Cowboy/sailor studies
    Bool_t  fTreeCascVarIsCowboy;   //store if V0 is cowboy-like or sailor-like in XY plane
    Float_t fTreeCascVarCowboyness; //negative -> cowboy, positive -> sailor
    Bool_t  fTreeCascVarIsCascadeCowboy;   //store if V0 is cowboy-like or sailor-like in XY plane
    Float_t fTreeCascVarCascadeCowboyness; //negative -> cowboy, positive -> sailor
    
    //Select charge (testing / checks)
    Int_t fkSelectCharge;
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //Sandbox mode for cascades
    Float_t fTreeCascVarPrimVertexX;
    Float_t fTreeCascVarPrimVertexY;
    Float_t fTreeCascVarPrimVertexZ;
    
    AliExternalTrackParam *fTreeCascVarBachTrack;//!
    AliExternalTrackParam *fTreeCascVarPosTrack; //!
    AliExternalTrackParam *fTreeCascVarNegTrack; //!
    
    //symmetry parameters (alpha etc)
    Double_t fTreeCascVarCosineProtonRestV0RestOmegaMinus; //!
    Double_t fTreeCascVarCosineProtonRestV0RestOmegaPlus; //!
    Double_t fTreeCascVarDzetaFromMomentaOmegaMinus; //!
    
    Float_t fTreeCascVarMagneticField;
    
    //extra
    Double_t lEffMassXi;
    Int_t lNegTrackSign;
    Int_t lPosTrackSign;
    Int_t lBachTrackSign;
    Int_t lEvSelCode;
    Double_t lPtbinlimitsCascade[40];
    Double_t lMinV0PA;
    Float_t lMaxV0CosPA;
    Float_t lMinV0CosPA;
    Float_t lMaxCascCosPA;
    Float_t lMinCascPA;
    Float_t lMaxBBCosPA;
    Float_t lMinBBPA;
    Float_t lLifetimeCut[4];
    Float_t lMass[4];
    Double_t dca;
    
    Float_t fTreeCascVarMassXiScratch;
    Float_t fTreeCascVarMassOmegaScratch;
    Float_t fTreeCascVarBackgroundMassOmega;
    Float_t fTreeCascVarBackgroundMassXi;
    Double_t fTreeVariableAlphaXi;
    Double_t fTreeVariablePtArmXi;
    Int_t fkEsdTrackMultiplicity; //!
    ULong64_t fkEsdEventTriggerWord; //!
    Float_t fTreeCascVarBachP;//!
    Float_t fTreeCascVarBachPt;//!
    Float_t fTreeCascVarNegP;//!
    Float_t fTreeCascVarNegPt;//!
    Float_t fTreeCascVarPosP;//!
    Float_t fTreeCascVarPosPt;//!
    Bool_t fkMomentaNeeded;
    Bool_t fkITSPID;
    Bool_t fkMassesAndBackgroundFromScratch;
    Bool_t fkDiscreteSymmetryInfo;
    Bool_t fkMultiplicityInfo;
    Bool_t fkExtraInfoCascadeKinematics;
    Int_t fRunNum;
    Int_t fNtracks;
    Float_t lMWindow[3];
    
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    
    //===========================================================================================
    //   Histograms
    //===========================================================================================
    
    TH1D *fHistEventCounter; //!
    TH1D *fHistEventCounterDifferential; //!
    TH1D *fHistCentrality; //!
    TH1D *fHistCentralityFast; //!
    TH1D *fHistZVertexEvent; //!
    
public:
    AliAnalysisTaskStrangeCascadesDiscrete();
    AliAnalysisTaskStrangeCascadesDiscrete(Bool_t lSaveEventTree,
                                                Bool_t lSaveV0Tree,
                                                Bool_t lSaveCascadeTree,
                                                Bool_t lRunVertexers,
                                                Bool_t lUseLightVertexer,
                                                Bool_t lMomentaNeeded,
                                                Bool_t lExtraInfoCascadeKinematics,
                                                Bool_t lDiscreteSymmetryInfo,
                                                Bool_t lMassesAndBackgroundFromScratch,
                                                Double_t lCascaderMaxChi2,
                                                Double_t lCascaderV0MinImpactParam,
                                                Double_t lCascaderV0MassWindow,
                                                Double_t lCascaderBachMinImpactParam,
                                                Double_t lCascaderMaxDCAV0andBach,
                                                Double_t lCascaderMinCosAngle,
                                                Double_t lCascaderMinRadius,
                                                Double_t lCascaderMaxRadius,
                                                const char *name, TString lExtraOptions = "");
    
    
    
    AliAnalysisTaskStrangeCascadesDiscrete(const AliAnalysisTaskStrangeCascadesDiscrete&); //copy constructor
    AliAnalysisTaskStrangeCascadesDiscrete& operator=(const AliAnalysisTaskStrangeCascadesDiscrete&);// copy assignment
    virtual ~AliAnalysisTaskStrangeCascadesDiscrete();
    
    virtual void   UserCreateOutputObjects();
    
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    Double_t MyRapidity(Double_t rE, Double_t rPz) const;
    Double_t GetCosOfProtonLambdaRestOmegaRest(Double_t pXi[3],Double_t mXi,
                                               Double_t pL[3], Double_t mL,
                                               Double_t pp[3], Double_t mp);
    Double_t DzetaFromMomenta(Double_t p1[3], Double_t p2[3], Double_t pp[3]);
    
    //Fix on-the-fly v0s
    void CheckChargeV0(AliESDv0 *v0);
    
    void SetSaveV0s                (Bool_t lSaveV0s        = kTRUE ) {
        fkSaveV0Tree        = lSaveV0s;
    }
    void SetSaveCascades           (Bool_t lSaveCascades   = kTRUE ) {
        fkSaveCascadeTree   = lSaveCascades;
    }
    void SetPreselectDedx (Bool_t lPreselectDedx= kTRUE ) {
        fkPreselectDedx   = lPreselectDedx;
    }
    void SetUseOnTheFlyV0Cascading( Bool_t lUseOnTheFlyV0Cascading = kTRUE ){
        //Highly experimental, use with care!
        fkUseOnTheFlyV0Cascading = lUseOnTheFlyV0Cascading;
    }
    
    //---------------------------------------------------------------------------------------
    //Task Configuration: trigger selection
    void SetSelectedTriggerClass(AliVEvent::EOfflineTriggerTypes trigType) { fTrigType = trigType;}
    //---------------------------------------------------------------------------------------
    //Task Configuration: Meant to enable quick re-execution of vertexer if needed
    void SetRunVertexers ( Bool_t lRunVertexers = kTRUE) {
        fkRunVertexers = lRunVertexers;
    }
    void SetUseLightVertexers ( Bool_t lUseLightVertexers = kTRUE) {
        fkUseLightVertexer = lUseLightVertexers;
    }
    void SetDoV0Refit ( Bool_t lDoV0Refit = kTRUE) {
        fkDoV0Refit = lDoV0Refit;
    }
    void SetExtraCleanup ( Bool_t lExtraCleanup = kTRUE) {
        fkExtraCleanup = lExtraCleanup;
    }
    //---------------------------------------------------------------------------------------
    void SetUseExtraEvSels ( Bool_t lUseExtraEvSels = kTRUE) {
        fkDoExtraEvSels = lUseExtraEvSels;
    }
    void SetPileupRejectionMode ( Int_t lMode = 1 ){
        //mode switch
        // 0 -> no rejection
        // 1 -> Ionut
        // 2 -> Anti-Ionut
        fkPileupRejectionMode = lMode;
    }
    void SetUseOldCentrality ( Bool_t lUseOldCent = kTRUE) {
        fkUseOldCentrality = lUseOldCent;
    }
    //---------------------------------------------------------------------------------------
    void SetSelectCharge ( Int_t lCharge = -1) {
        fkSelectCharge = lCharge;
    }
    //---------------------------------------------------------------------------------------
    //Task Configuration: Skip Event Selections after trigger (VZERO test)
    void SetDownScaleV0 ( Bool_t lOpt = kTRUE, Float_t lVal = 0.001) {
        fkDownScaleV0 = lOpt;
        fDownScaleFactorV0 = lVal;
    }
    void SetDownScaleCascade ( Bool_t lOpt = kTRUE, Float_t lVal = 0.001 ) {
        fkDownScaleCascade = lOpt;
        fDownScaleFactorCascade = lVal;
    }
    //---------------------------------------------------------------------------------------
    //Setters for the V0 Vertexer Parameters
    void SetV0VertexerMaxChisquare   ( Double_t lParameter ) {
        fV0VertexerSels[0] = lParameter;
    }
    void SetV0VertexerDCAFirstToPV   ( Double_t lParameter ) {
        fV0VertexerSels[1] = lParameter;
    }
    void SetV0VertexerDCASecondtoPV  ( Double_t lParameter ) {
        fV0VertexerSels[2] = lParameter;
    }
    void SetV0VertexerDCAV0Daughters ( Double_t lParameter ) {
        fV0VertexerSels[3] = lParameter;
    }
    void SetV0VertexerCosinePA       ( Double_t lParameter ) {
        fV0VertexerSels[4] = lParameter;
    }
    void SetV0VertexerMinRadius      ( Double_t lParameter ) {
        fV0VertexerSels[5] = lParameter;
    }
    void SetV0VertexerMaxRadius      ( Double_t lParameter ) {
        fV0VertexerSels[6] = lParameter;
    }
    //---------------------------------------------------------------------------------------
    //Setters for the Cascade Vertexer Parameters
    void SetCascVertexerMaxChisquare         ( Double_t lParameter ) {
        fCascadeVertexerSels[0] = lParameter;
    }
    void SetCascVertexerMinV0ImpactParameter ( Double_t lParameter ) {
        fCascadeVertexerSels[1] = lParameter;
    }
    void SetCascVertexerV0MassWindow         ( Double_t lParameter ) {
        fCascadeVertexerSels[2] = lParameter;
    }
    void SetCascVertexerDCABachToPV          ( Double_t lParameter ) {
        fCascadeVertexerSels[3] = lParameter;
    }
    void SetCascVertexerDCACascadeDaughters  ( Double_t lParameter ) {
        fCascadeVertexerSels[4] = lParameter;
    }
    void SetCascVertexerCascadeCosinePA      ( Double_t lParameter ) {
        fCascadeVertexerSels[5] = lParameter;
    }
    void SetCascVertexerCascadeMinRadius     ( Double_t lParameter ) {
        fCascadeVertexerSels[6] = lParameter;
    }
    void SetCascVertexerCascadeMaxRadius     ( Double_t lParameter ) {
        fCascadeVertexerSels[7] = lParameter;
    }
    //---------------------------------------------------------------------------------------
    void SetMinPt     ( Float_t lMinPt ) {
        fMinPtToSave = lMinPt;
    }
    void SetMaxPt     ( Float_t lMaxPt ) {
        fMaxPtToSave = lMaxPt;
    }
    void SetLambdaWindowParameters     ( Double_t *fMeanPars, Double_t *fSigmaPars ) {
        for(Int_t ipar=0; ipar<5; ipar++) fLambdaMassMean[ipar]  = fMeanPars[ipar];
        for(Int_t ipar=0; ipar<4; ipar++) fLambdaMassSigma[ipar] = fSigmaPars[ipar];
    }
    void SetLambdaWindowParametersStandard (){
        fLambdaMassMean[0] =  1.15768e+00;
        fLambdaMassMean[1] = -4.15945e-02;
        fLambdaMassMean[2] = -7.14294e-04;
        fLambdaMassMean[3] = -1.62793e-02;
        fLambdaMassMean[4] = -7.84067e+00;
        fLambdaMassSigma[0] = 1.30345e-03;
        fLambdaMassSigma[1] = 2.89679e-04;
        fLambdaMassSigma[2] = 1.52661e-03;
        fLambdaMassSigma[3] =-2.58251e+00;
    }
    //---------------------------------------------------------------------------------------
    //Superlight mode: add another configuration, please
    void AddConfiguration( AliV0Result      *lV0Result      );
    void AddConfiguration( AliCascadeResult *lCascadeResult );
    //---------------------------------------------------------------------------------------
    //Functions for analysis Bookkeepinp
    // 1- Configure standard vertexing
    void SetupStandardVertexing();
    void SetupLooseVertexing();
    // 2- Standard Topological Selection QA Sweeps
    void AddTopologicalQAV0(Int_t lRecNumberOfSteps = 100);
    void AddTopologicalQACascade(Int_t lRecNumberOfSteps = 100);
    // 3 - Standard analysis configurations + systematics
    void AddStandardV0Configuration(Bool_t lUseFull = kFALSE, Bool_t lDoSweepLooseTight = kFALSE, Int_t lSweepFullNumb = 0);
    void AddStandardV0RadiusSweep();
    void AddStandardCascadeConfiguration(Bool_t lUseFull = kFALSE, Bool_t lDoSystematics = kTRUE);
    void AddCascadeConfiguration276TeV(); //Adds old 2.76 PbPb cut level analyses
    void AddCascadeConfigurationPreliminaryCrosscheck(); //
    //---------------------------------------------------------------------------------------
    Float_t GetDCAz(AliESDtrack *lTrack);
    Float_t GetCosPA(AliESDtrack *lPosTrack, AliESDtrack *lNegTrack, AliESDEvent *lEvent);
    Bool_t GoodTrack(AliESDtrack* trackESD);
    Int_t GetITSstatus(const AliVTrack * const track, Int_t layer) const;
    Bool_t CheckTOFstatus(AliESDtrack* trackESD);
    Int_t EventTriggerWord(AliESDEvent* lESDevent);
    
    
    
    
    ClassDef(AliAnalysisTaskStrangeCascadesDiscrete, 1);
};

#endif



