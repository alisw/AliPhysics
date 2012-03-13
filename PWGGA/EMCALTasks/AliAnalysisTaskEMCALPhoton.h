#ifndef AliAnalysisTaskEMCALPhoton_cxx
#define AliAnalysisTaskEMCALPhoton_cxx

class TH1;
class TH2;
class TObjArray;
class AliESDEvent;
class AliMCEvent;
class AliStack;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDCaloCells;
class AliEMCALGeometry;
class AliVCluster;
class AliAnalysisTaskEMCALClusterizeFast;
class TParticle;
class MyHeaderObj;
class MyConversionObj;
class MyClusterObj;
class MyCellObj;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALPhoton : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEMCALPhoton() : 
  AliAnalysisTaskSE(), 
  
    fTrCuts(0),
    fPrTrCuts(0),
    fSelTracks(0),
    fSelPrimTracks(0),
    fPhotConvArray(0),
    fMyClusts(0),
    fMyAltClusts(0),
    fMyCells(0),
    fMyTracks(0),
    fMyMcParts(0),
    fHeader(0x0),
    fCaloClusters(0),
    fCaloClustersNew(0),
    fEMCalCells(0),
    fGeom(0x0),
    fTimeResTOF(0),
    fMipResponseTPC(0),
    fGeoName("EMCAL_COMPLETEV1"),
    fPeriod("LHC11d"),
    fIsTrain(0),
    fIsMC(0),
    fIsGrid(0),
    fClusThresh(2.0),

    fClusterizer(0),
    fCaloClustersName("EmcalClusterv2"),

    fESD(0),
    fMCEvent(0),
    fStack(0x0),
    
    fOutputList(0),
    fTree(0),
    
    fNV0sBefAndAftRerun(0),
    fConversionVtxXY(0),
    fInvMassV0(0),
    fInvMassV0KF(0),
    fInvMassV0SS(0),
    fDedxPAll(0)
  
  
  {}
  AliAnalysisTaskEMCALPhoton(const char *name);
  virtual ~AliAnalysisTaskEMCALPhoton() {}

  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *);

  void         SetTrackCuts(AliESDtrackCuts *c)                                  { fTrCuts             =    c;    }
  void         SetPrimTrackCuts(AliESDtrackCuts *c)                              { fPrTrCuts           =    c;    }
  void         SetTimeResTOF(Float_t tr = 130.)                                  { fTimeResTOF         =    tr;   }
  void         SetMipResponseTPC(Float_t mr = 47.9)                              { fMipResponseTPC     =    mr;   }
  void         SetGeoName(const char *n)                                         { fGeoName            =    n;    }
  void         SetPeriod(const char *n)                                          { fPeriod             =    n;    }
  void         SetTrainMode(Bool_t t)                                            { fIsTrain            =    t;    }
  void         SetGridMode(Bool_t g)                                             { fIsGrid             =    g;    }
  void         SetClusThreshold(Double_t et)                                     { fClusThresh         =    et;   }
  void         SetClusterizer(AliAnalysisTaskEMCALClusterizeFast *c)             { fClusterizer        =    c;    }
  void         SetMcMode(Bool_t mc)                                              { fIsMC               =    mc;   }
  void         FindConversions();
  void         FillMyCells();
  void         FillMyClusters();
  void         FillMyAltClusters();
  void         FillHighPtTracks();
  void         FillMcPart(TParticle *mcP, Int_t ipart, Int_t itrack);
  void         GetMcParts();
  Double_t     GetTrackIsolation(Double_t cEta, Double_t cPhi, Double_t radius=0.2, Double_t pt=0.)       const;
  Double_t     GetPhiBandEt(Double_t cEta, Double_t cPhi, Double_t radius=0.2, Double_t pt=0.)            const;
 // Double_t     GetPhiBandEt(Double_t cEta, Double_t cPhi, Double_t radius=0.2, Double_t pt=0.)            const;
  Double_t     GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax);
  Double_t     GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const; 
  
 protected:
  AliESDtrackCuts                       *fTrCuts;                 // track cuts
  AliESDtrackCuts                       *fPrTrCuts;               // primary track cuts
  TObjArray                             *fSelTracks;             //!pointer to selected inclusive tracks
  TObjArray                             *fSelPrimTracks;         //!pointer to selected primary tracks
  TClonesArray                          *fPhotConvArray;         //!array of MyConversionObj
  TClonesArray                          *fMyClusts;              //!array of MyClusterObj
  TClonesArray                          *fMyAltClusts;           //!array of MyClusterObj from the alternative clusterizer
  TClonesArray                          *fMyCells;               //!array of MyCellObj
  TClonesArray                          *fMyTracks;              //!array of MyTrackObj
  TClonesArray                          *fMyMcParts;             //!array of MyMcPartObj
  MyHeaderObj                           *fHeader;                //!
  TRefArray                             *fCaloClusters;          //!pointer to EMCal clusters
  TClonesArray                          *fCaloClustersNew;       //!pointer to EMCal clusters v2
  AliESDCaloCells                       *fEMCalCells;            //!pointer to EMCal cells
  AliEMCALGeometry                      *fGeom;                   // geometry utils
  Float_t                                fTimeResTOF;            //TOF time resolution for track PID
  Float_t                                fMipResponseTPC;        //TPC mip response for track pid
  TString                                fGeoName;                // geometry name (def = EMCAL_FIRSTYEARV1)
  TString                                fPeriod;                 // string to the LHC period
  Bool_t                                 fIsTrain;                //variable to set train mode
  Bool_t                                 fIsMC;                   //variable to switch mcparts branch on/off
  Bool_t                                 fIsGrid;                //variable to set grid mode
  Double_t                               fClusThresh;            //!energy threshold for cluster be saved
  AliAnalysisTaskEMCALClusterizeFast    *fClusterizer;           //!pointer for alternative clusterizer
  TString                                fCaloClustersName;      //alternative clusterizer name
  
  
 private:
  AliESDEvent                           *fESD;      //! ESD object
  AliMCEvent                            *fMCEvent;    //! MC event object
  AliStack                              *fStack;     //!MC particles stack object
  
  TList                                 *fOutputList;            //! Output list
  TTree                                 *fTree;                  //!output tree
  
  //conversion histograms
  TH2F                                  *fNV0sBefAndAftRerun;      //!check the number of V0s before and after rerun
  TH2F                                  *fConversionVtxXY;         //! X x Y of found conversion vertices
  TH1F                                  *fInvMassV0;               //!invariant mass from v0->GetEffMass()
  TH1F                                  *fInvMassV0KF;             //!invariant mass from the v0 tracks
  TH1F                                  *fInvMassV0SS;             //!invariant mass from the tracks in the "dirty" finder
  TH2F                                  *fDedxPAll;               //!dE/dx vs p of all selected tracks
  

   
  AliAnalysisTaskEMCALPhoton(const AliAnalysisTaskEMCALPhoton&); // not implemented
  AliAnalysisTaskEMCALPhoton& operator=(const AliAnalysisTaskEMCALPhoton&); // not implemented
  
  ClassDef(AliAnalysisTaskEMCALPhoton, 1); // example of analysis
};

#endif

#ifndef MyObjs_h
#define MyObjs_h

class MyHeaderObj : public TObject
{
  public: MyHeaderObj() :
        TObject(), fTrClassMask(0), fTrCluster(0), fV0Cent(0), fV0(0), fCl1Cent(0), 
        fCl1(0), fTrCent(0), fTr(0), fNClus(0), fNCells(0)  {;}
  public:
  ULong64_t     fTrClassMask;    //         trigger class mask
  UChar_t       fTrCluster;      //         trigger cluster mask
  Double32_t    fV0Cent;         //[0,0,16] v0 cent
  Double32_t    fV0;             //[0,0,16] v0 result used for cent 
  Double32_t    fCl1Cent;        //[0,0,16] cl1 cent
  Double32_t    fCl1;            //[0,0,16] cl1 result used for cent 
  Double32_t    fTrCent;         //[0,0,16] tr cent
  Double32_t    fTr;             //[0,0,16] tr result used for cent 
  Int_t         fNClus;
  Int_t         fNCells;

  ClassDef(MyHeaderObj,2)
};

class MyConversionObj : public TObject
{
  public: MyConversionObj() : 
        TObject(), fPt(0), fEta(0), fPhi(0), fVR(0), fVEta(0), fVPhi(0), fMass(0), fMcLabel(-1),
               fNegPt(0), fNegEta(0), fNegPhi(0), fNegDedx(0), fNegMcLabel(-1),
               fPosPt(0), fPosEta(0), fPosPhi(0), fPosDedx(0), fPosMcLabel(-1) {;}
 public:
  Double32_t    fPt;               //[0,0,16] pt
  Double32_t    fEta;              //[0,0,16] eta
  Double32_t    fPhi;              //[0,0,16] phi
  Double32_t    fVR;               //[0,0,16] prod r (cylinder)
  Double32_t    fVEta;             //[0,0,16] prod eta
  Double32_t    fVPhi;             //[0,0,16] prod phi
  Double32_t    fMass;             //[0,0,16] if correctly filled, should be <50 MeV
  Short_t       fMcLabel;          //corresponding MC label

  //negative daughter
  Double32_t    fNegPt;               //[0,0,16] pt
  Double32_t    fNegEta;              //[0,0,16] eta
  Double32_t    fNegPhi;              //[0,0,16] phi
  Double32_t    fNegDedx;             //[0,0,16] if correctly filled, should be <50 MeV
  Short_t       fNegMcLabel;          //corresponding MC label

  //positive daughter
  Double32_t    fPosPt;               //[0,0,16] pt
  Double32_t    fPosEta;              //[0,0,16] eta
  Double32_t    fPosPhi;              //[0,0,16] phi
  Double32_t    fPosDedx;             //[0,0,16] if correctly filled, should be <50 MeV
  Short_t       fPosMcLabel;          //corresponding MC label

  ClassDef(MyConversionObj,1) // conversion class

};
class MyClusterObj : public TObject
{
  public: MyClusterObj() : 
  TObject(), fE(0), fEt(0), fR(0), fEta(0), fPhi(0), fN(0),fEmax(0),fTmax(0), fIdmax(0), fEcross(0),fDisp(-1), 
        fM20(-1), fM02(-1),fTrDEta(0), fTrDPhi(0), fTrEp(0), fTrDedx(0), fTrIso01(0), fTrIso02(0), fTrIso03(0), fTrIso04(0), 
        fTrPhiBand01(0), fTrPhiBand02(0), fTrPhiBand03(0), fTrPhiBand04(0), fCellsAbsId(""),fMcLabel(-1)
        {;}
 public:
  Double32_t   fE;
  Double32_t   fEt;
  Double32_t   fR;
  Double32_t   fEta;
  Double32_t   fPhi;
  UShort_t     fN;
  Double_t     fEmax;
  Double_t     fTmax;
  Short_t      fIdmax;
  Double_t     fEcross;
  Double32_t   fDisp;
  Double32_t   fM20;
  Double32_t   fM02;
  Double32_t   fTrDEta;
  Double32_t   fTrDPhi;
  Double32_t   fTrEp;
  Double32_t   fTrDedx;
  Double32_t   fTrIso01;
  Double32_t   fTrIso02;
  Double32_t   fTrIso03;
  Double32_t   fTrIso04;
  Double32_t   fTrPhiBand01;
  Double32_t   fTrPhiBand02;
  Double32_t   fTrPhiBand03;
  Double32_t   fTrPhiBand04;
  TString      fCellsAbsId;           //cluster cells absid
  Short_t      fMcLabel;
  
  
  
  ClassDef(MyClusterObj,6) // cluster class

};

class MyCellObj : public TObject
{
  public: MyCellObj() : 
        TObject(), fAbsID(-1), fE(0), fEt(0), fEta(0), fPhi(0), fTime(0)
        {;}
 public:
  Short_t      fAbsID;
  Double32_t   fE;
  Double32_t   fEt;
  Double32_t   fEta;
  Double32_t   fPhi;
  Double32_t   fTime;
  
  
  
  ClassDef(MyCellObj,1) // cell class

};

class MyTrackObj : public TObject
{
  public: MyTrackObj() :
        TObject(), fPt(0), fEta(0), fPhi(0), fDedx(0), fCharge(0), fMcLabel(-1) {;}
  public:
  Double32_t fPt;
  Double32_t fEta;
  Double32_t fPhi;
  Double32_t fDedx;
  Short_t    fCharge;
  Short_t    fMcLabel;

  ClassDef(MyTrackObj,3)
};

class MyMcPartObj : public TObject
{
  public: MyMcPartObj() :
        TObject(), fLabel(-1), fPdg(0), fPt(0), fEta(0), fPhi(0), 
        fVR(0), fVEta(0), fVPhi(0), fMother(-1)  {;}
  public:
  Short_t    fLabel;
  Short_t    fPdg;
  Double32_t fPt;
  Double32_t fEta;
  Double32_t fPhi;
  Double32_t fVR;
  Double32_t fVEta;
  Double32_t fVPhi;
  Short_t    fMother;

  ClassDef(MyMcPartObj,1)
};

#endif
