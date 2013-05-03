#ifndef AliAnalysisTaskEMCALPhoton_h
#define AliAnalysisTaskEMCALPhoton_h

// $Id$

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
class AliPhotonHeaderObj;
class AliPhotonConvObj;
class AliPhotonClusterObj;
class AliPhotonCellObj;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALPhoton : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEMCALPhoton();
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
  void         SetDebugMode(Bool_t d)                                            { fDebug              =    d;    }
  void         FindConversions();
  void         FillMyCells();
  void         FillMyClusters();
  void         FillMyAltClusters();
  void         FillIsoTracks();
  void         FillMcPart(TParticle *mcP,  Int_t itrack);
  void         GetMcParts();
  Double_t     GetMcIsolation(TParticle *mcP, Int_t itrack, Double_t radius, Double_t pt)                 const;
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
  TClonesArray                          *fPhotConvArray;         //!array of AliPhotonConvObj
  TClonesArray                          *fMyClusts;              //!array of AliPhotonClusterObj
  TClonesArray                          *fMyAltClusts;           //!array of AliPhotonClusterObj from the alternative clusterizer
  TClonesArray                          *fMyCells;               //!array of AliPhotonCellObj
  TClonesArray                          *fMyTracks;              //!array of AliPhotonTrackObj
  TClonesArray                          *fMyMcParts;             //!array of AliPhotonMcPartObj
  AliPhotonHeaderObj                    *fHeader;                //!
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
  Bool_t                                 fDebug;                 //variable to switch debug on/off
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

#ifndef AliPhotonObjs_h
#define AliPhotonObjs_h

class AliPhotonHeaderObj : public TObject
{
  public: AliPhotonHeaderObj() :
  TObject(), fInputFileName(""), fTrClassMask(0), fTrCluster(0), fV0Cent(0), fV0(0), fCl1Cent(0), 
    fCl1(0), fTrCent(0), fTr(0), fNClus(0), fNCells(0), fTrackMult(0), fNMcParts(0)  {;}
  public:
  TString       fInputFileName;  // used for normalization purposes in MC productions
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
  Int_t         fTrackMult;
  Int_t         fNMcParts;

  ClassDef(AliPhotonHeaderObj,5)
};

class AliPhotonConvObj : public TObject
{
  public: AliPhotonConvObj() : 
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

  ClassDef(AliPhotonConvObj,1) // conversion class

};
class AliPhotonClusterObj : public TObject
{
  public: AliPhotonClusterObj() : 
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
  
  
  
  ClassDef(AliPhotonClusterObj,6) // cluster class

};

class AliPhotonCellObj : public TObject
{
  public: AliPhotonCellObj() : 
        TObject(), fAbsID(-1), fE(0), fEt(0), fEta(0), fPhi(0), fTime(0)
        {;}
 public:
  Short_t      fAbsID;
  Double32_t   fE;
  Double32_t   fEt;
  Double32_t   fEta;
  Double32_t   fPhi;
  Double32_t   fTime;
  
  
  
  ClassDef(AliPhotonCellObj,1) // cell class

};

class AliPhotonTrackObj : public TObject
{
  public: AliPhotonTrackObj() :
        TObject(), fPt(0), fEta(0), fPhi(0), fDedx(0), fCharge(0), fMcLabel(-1) {;}
  public:
  Double32_t fPt;
  Double32_t fEta;
  Double32_t fPhi;
  Double32_t fDedx;
  Short_t    fCharge;
  Short_t    fMcLabel;

  ClassDef(AliPhotonTrackObj,3)
};

class AliPhotonMcPartObj : public TObject
{
  public: AliPhotonMcPartObj() :
  TObject(), fLabel(-1), fPdg(0), fPt(0), fEta(0), fPhi(0), 
    fVR(0), fVEta(0), fVPhi(0), fMother(-1), fFirstD(-1),
    fLastD(-1), fStatus(-1), fIso(-1) {;}
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
  Short_t    fFirstD;
  Short_t    fLastD;
  Short_t    fStatus;
  Double32_t fIso;

  ClassDef(AliPhotonMcPartObj,2)
};

#endif
