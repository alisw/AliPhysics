#ifndef AliAnalysisTaskHe3TriTree_H
#define AliAnalysisTaskHe3TriTree_H
#include "AliTRDonlineTrackMatching.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAnalysisAlien.h"

class AliPIDResponse;

class AliAnalysisTaskHe3TriTree : public AliAnalysisTaskSE  
{
  public:
                            AliAnalysisTaskHe3TriTree();
                            AliAnalysisTaskHe3TriTree(const char *name);
    virtual                 ~AliAnalysisTaskHe3TriTree();
    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);
    void                    SetParamsHe(Double_t params[6]) { for(Int_t i=0; i < 6; i++) fBetheParamsHe[i] = params[i];};
    void                    SetParamsT(Double_t params[6]) {for (Int_t i = 0; i < 6; i++) fBetheParamsT[i] = params[i];};
    void                    SetUseExternalSplines(Bool_t extspl = kFALSE) {fUseExternalSplines = extspl;}; 

  private:
    AliESDInputHandler      *fInputHandler;      //!<! Input handler
    AliMCEvent*             mcEvent;               //! MC event
    AliPIDResponse*         fPIDResponse;          //! pid response object
    AliStack                *fStack;               //!<! MC stack
    Bool_t                  fMCtrue;               //< flag for MC events
    AliESDEvent*            fESDevent;             //! input event
    AliEventCuts            fEventCuts;            //< event cuts as advised by PDG (AliEventCuts)
    TList*                  fOutputList;           //! output list
    TH1F*                   fHistEvents;           //! number of events and trigger info
    TH2D*           fHistEventsTRD[2];     //! histograms for TRD trigger
    TH2F*           fHistdEdx;             //! TPC dEdx histogram
    Double_t        fBetheParamsHe[6];     //< Bethe Aleph He3 Parameter + TPC sigma
    Double_t        fBetheParamsT[6];  
    TTree           *fTree;                //< tree containing He3 information
    Bool_t          fUseExternalSplines;   //< Use Splines given in Add file   
    AliTRDonlineTrackMatching *matching;    //< match TRD with ESD tracks 
    Int_t           fYear;  
    Int_t           tRunNumber;
    ULong64_t       tEventId;               //<  event id
    Int_t           tTrigMB;               // trigger info
    Int_t           tTrigHMV0;
    Int_t           tTrigHMSPD;
    Int_t           tTrigHNU;
    Int_t           tTrigHQU;
    
    Double_t        tMagField;
    
    Int_t           tTRDvalid;             // has valid TRD track
    Int_t           tTRDtrigHNU;  
    Int_t           tTRDtrigHQU;  
    Int_t           tTRDPid;
    Int_t           tTRDnTracklets;
    Int_t           tTRDPt;
    Int_t           tTRDLayerMask;
    Double_t        tTRDSagitta;
    Int_t           tTRDStack;
    Int_t           tTRDSector;
    UInt_t          tTRDPID0;
    UInt_t          tTRDPID1;
    UInt_t          tTRDPID2;
    UInt_t          tTRDPID3;
    UInt_t          tTRDPID4;
    UInt_t          tTRDPID5;  
      
    Double_t        tnTPCcluster;
    Double_t        tnITScluster;
    Double_t        tTPCchi2;

    Double_t        tSPDFiredChips0;  // multiplicity triggers
    Double_t        tSPDFiredChips1;
    Double_t        tSPDTracklets;
    Double_t        tSPDCluster;
    Double_t        tV0Multiplicity;
    Double_t        tNTracks;

    Double_t        tMultV0M;      // multiplicity estimators
    Double_t        tMultOfV0M;      
    Double_t        tMultSPDTracklet;  
    Double_t        tMultSPDCluster;  
    Double_t        tMultRef05;      
    Double_t        tMultRef08;       

    Double_t        tCharge;

    Double_t        tPt;        // He3 track parameter
    Double_t        tY;
    Double_t        tEta;
    Double_t        tPhi;
    Double_t        tPx;        
    Double_t        tPy;
    Double_t        tPz;
    Double_t        tE;
    Int_t           tKink;
    Int_t           tTPCrefit;
    Int_t           tITSrefit;  

    Double_t        tP;          // PID parameter 
    Double_t        tHeDEdx;
    Double_t        tHeSigma;
    Double_t        tTOFSignalHe;

    Double_t        tDcaXY;        // impact parameters
    Double_t        tDcaZ;
    Double_t        tSigmaYX;
    Double_t        tSigmaXYZ;
    Double_t        tSigmaZ;

    Int_t           tMCtrue;      // MC info
    Int_t           tPrimary;
    Int_t           tWeak;
    Int_t           tMaterial;
      
    TTree           *fTreeGenHe;           //< tree containing generated He3/3H information
    Double_t        tGenCharge;
    Double_t        tGenPt;        // He3 track parameter
    Double_t        tGenY;
    Int_t           tGenPrimary;

    TTree           *fTreeGenTri;
    Double_t        tGenTriCharge;
    Double_t        tGenTriPt;        // 3H track parameter
    Double_t        tGenTriY;
    Int_t           tGenTriPrimary;
    
    TTree	          *fTreeTri;
    Double_t		    tTriP;
    Double_t	      tChargeTri;
    Double_t	      tTriPt;
    Double_t	      tTriY;
    Double_t	      tTriEta;
    Double_t        tTriPhi;
    Double_t        tTriPx;        
    Double_t        tTriPy;
    Double_t        tTriPz;
    Double_t        tTriE;
    Int_t           tTriKink;
    Int_t           tTriTPCrefit;
    Int_t           tTriITSrefit;
    Double_t        tTriSigma;
    Double_t        tTriDedx;
    Double_t        tTOFSignalTri;
    Double_t        tTriDcaXY;
    Double_t        tTriDcaZ; 
    Double_t        tTriSigmaYX;
    Double_t        tTriSigmaXYZ;
    Double_t        tTriSigmaZ;

    Long_t          pdgCode;
    Long_t          tPdg;

    void            ProcessMCParticles();
    Double_t        GetInvPtDevFromBC(Int_t b, Int_t c);
    Double_t        Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params);
    Double_t        GetTOFSignalHe3(AliESDtrack& trackHe, Double_t tP);
    Double_t        GetTOFSignalTri(AliESDtrack& trackTri, Double_t tP);
    void            SetBetheBlochParams(Int_t runNumber);
    Double_t        TRDtrack(AliESDtrack* esdTrack);
    void            TRDnEvents();
    AliAnalysisTaskHe3TriTree(const AliAnalysisTaskHe3TriTree&); // not implemented
    AliAnalysisTaskHe3TriTree& operator=(const AliAnalysisTaskHe3TriTree&); // not implemented

    ClassDef(AliAnalysisTaskHe3TriTree, 5);
};

#endif
