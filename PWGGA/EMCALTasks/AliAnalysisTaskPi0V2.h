#ifndef ALIANALYSISTASKPI0V2_H
#define ALIANALYSISTASKPI0V2_H

class TH1F;
class TH2F;
class TH3F;
class TList;
class AliESDCaloCluster;
class AliVCluster;
class AliESDtrackCuts;
class AliESDEvent;
class THnSparse;
class TClonesArray;
class TString;
class TProfile;
class TProfile2D;
class AliOADBContainer;
class AliEPFlattener;
class AliEMCALGeometry;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTaskPi0V2 : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskPi0V2();
    AliAnalysisTaskPi0V2(const char *name);
    virtual ~AliAnalysisTaskPi0V2();
    
    virtual void           UserCreateOutputObjects();
    virtual void           UserExec(Option_t *option);
    virtual void           Terminate(Option_t *);

    void                   SetTrigClass(const char *n)      {fTrigClass  = n;}
    void		           SetVzCut(Double_t v )	        {fVzCut     = v;}
    void		           SetClustNCell(Double_t c )	    {fNCellCut   = c;}
    void		           SetClustE(Double_t e )	        {fECut       = e;}
    void		           SetClustEta(Double_t e )	        {fEtaCut     = e;}
    void		           SetV2M02Cut(Double_t m )	        {fV2M02Cut   = m;}
    void                   SetV1M02Cut(Double_t m )         {fV1M02Cut   = m;}
    void		           SetDrCut(Double_t m )	        {fDrCut      = m;}
    void		           SetPi0Asy(Double_t a )	        {fPi0AsyCut  = a;}
    void                   UseV2Clust(Bool_t e)             {fUseV2Clust  = e;}
    void                   UseV1Clust(Bool_t e)             {fUseV1Clust  = e;}
    void                   UseTrack(Bool_t e)               {fUseTrk = e;}
    void                   SetV2ClustName(TString n)        {fV2ClustName = n;} 
    void                   SetV1ClustName(TString n)        {fV1ClustName = n;} 
    void                   SetTrackName(TString n)          {fTrackName = n;}
    void		           UsePhosEPCali(Bool_t e)		    {fUsePhosEPCali = e;}
    void                   SetPhosEPCaliFileName(TString n) {fPhosEPCaliFileName = n;}
    void		           FlattenMostCent(Bool_t e)		{fFlattenMostCent = e;}
    void                   FlattenSemiCent(Bool_t e)        {fFlattenSemiCent = e;}    

 private:
    Int_t                  ConvertToInternalRunNum(Int_t n);
    Bool_t                 IsCentAccepted();
    void                   VZEROEventPlane(Bool_t flattenEP);
    Double_t               FlattenV0A(Double_t phi, Double_t c);
    Double_t               FlattenV0C(Double_t phi, Double_t c);
    Double_t               FlattenTPC(Double_t phi, Double_t c);
    Bool_t                 IsGoodV2Cluster(const AliVCluster *c) const;
    Double_t               GetMaxCellEnergy(const AliVCluster *c, Short_t &id) const;
    Double_t               GetCrossEnergy(const AliVCluster *c, Short_t &idmax) const;
    Bool_t                 IsWithinFiducialVolume(Short_t id) const;
    void                   FillPion(const TLorentzVector& p1, const TLorentzVector& p2);
    Bool_t                 IsGoodV1Cluster(const AliVCluster *c) const;
    Bool_t                 IsPi0Candidate(const AliVCluster *c);
    void                   FillPion(const TLorentzVector& p1, AliVCluster* c);
    void                   GetMom(TLorentzVector& p, const AliVCluster* c, Double_t* vertex);      

    TList*                 fOutput;	        //!Output list
    AliESDEvent*		   fESD;		        //!ESD object
    AliAODEvent*		   fAOD;		        //!AOD object
    AliEMCALGeometry*      fGeom;                 // geometry utils
    TString                fGeomName;
    AliOADBContainer*      fPhosEPCaliContainer;
    TString                fPhosEPCaliFileName;      // Name for calibration
    AliEventplane*         fEventPlane;
    TString                fTrackName;	        // name of track collection
    TString                fV2ClustName;           // name of V1 Clus collection
    TString                fV1ClustName;	        // name of V1 Clus collection
    TString                fTrigClass;	        // trigger class name for event selection
    TClonesArray*          fTrack;		//! pico tracks specific for Skim ESD
    TClonesArray*          fV2Clust;      //! Cluster Array for V2
    TClonesArray*          fV1Clust;		//! Cluster Array for V1
    Int_t			       fRunNum;		//! Run numbers
    Int_t			       fInterRunNum;	//! Run numbers
    Double_t		       fVzCut;		// vertex cut
    Double_t			   fNCellCut;	        // N cells Cut
    Double_t			   fECut;			// Cluster E cut
    Double_t			   fEtaCut;		// Cluster Eta Cut
    Double_t			   fV2M02Cut;		// Cluster long axis cut
    Double_t               fV1M02Cut;      // Cluster long axis cut
    Double_t			   fDrCut;		        // Cluster long axis cut
    Bool_t			       fPi0AsyCut;		// pion Asymetry cut 0=off 1=on
    Bool_t			       fUseV2Clust;		// pion Asymetry cut 0=off 1=on
    Bool_t                 fUseV1Clust;       // pion Asymetry cut 0=off 1=on
    Bool_t                 fUseTrk;       // pion Asymetry cut 0=off 1=on
    Bool_t			       fUsePhosEPCali;		// use Phos flattening

    Double_t			   fCentrality;	  	//! Centrality
    Bool_t                 fFlattenMostCent;
    Bool_t                 fFlattenSemiCent;

    Double_t			   fEPTPC;			//! Evt plane TPC
    Double_t			   fEPTPCReso;		//! resolution of TPC method
    Double_t			   fEPV0;	    		//! EP V0
    Double_t			   fEPV0A;			//! EP V0A
    Double_t			   fEPV0C;	  		//! EP V0C
    Double_t			   fEPV0AR;		//! EP V0A reduced
    Double_t			   fEPV0CR;		//! EP V0C reduced
    Double_t			   fEPV0R;	  		//! EP V0 reduced
    Double_t			   fEPV0AR4;		//! EP V0A ring4 only
    Double_t			   fEPV0AR5;		//! EP V0A ring5 only
    Double_t			   fEPV0AR6;		//! EP V0A ring6 only
    Double_t			   fEPV0AR7;		//! EP V0A ring7 only
    Double_t			   fEPV0CR0;		//! EP V0C ring0 only	
    Double_t			   fEPV0CR1;		//! EP V0C ring1 only	
    Double_t			   fEPV0CR2;		//! EP V0C ring2 only	
    Double_t			   fEPV0CR3;		//! EP V0C ring3 only	
    AliEPFlattener*        fEPTPCFlat;      //! Object for flattening of TPC
    AliEPFlattener*        fEPV0AFlat;      //! Object for flattening of V0A
    AliEPFlattener*        fEPV0CFlat;      //! Object for flattening of V0C

    TH1F*                  hEvtCount;		//!
    TH1F*                  hCentA;			//!
    TH1F*                  hCentB;

    TH2F*                  hEPTPC;		//! 2-D histo EPTPC  vs cent
    TH2F*                  hEPTPCReso;	//! 2-D histo TPC resolution vs cent
    TH2F*                  hEPV0A;		//! 2-D histo EPV0A  vs cent
    TH2F*                  hEPV0C;		//! 2-D histo EPV0C  vs cent
    TH2F*                  hEPTPCFlat;        //! 2-D histo EPTPC  vs cent
    TH2F*                  hEPV0AFlat;        //! 2-D histo EPV0A  vs cent
    TH2F*                  hEPV0CFlat;        //! 2-D histo EPV0C  vs cent
    TH2F*                  hEPV0;         //! 2-D histo EPV0   vs cent
    TH2F*                  hEPV0AR;		//! 2-D histo EPV0Ar vs cent
    TH2F*                  hEPV0CR;		//! 2-D histo EPV0Cr vs cent
    TH2F*                  hEPV0R;		//! 2-D histo EPV0r  vs cent
    TH2F*                  hEPV0AR4;		//! 2-D histo EPV0AR4 vs cent
    TH2F*                  hEPV0AR7;		//! 2-D histo EPV0AR7 vs cent
    TH2F*                  hEPV0CR0;		//! 2-D histo EPV0AR0 vs cent
    TH2F*                  hEPV0CR3;		//! 2-D histo EPV0AR3 vs cent
    TH2F*                  hEPDiffV0A_V0CR0;
    TH2F*                  hEPDiffV0A_V0CR3;
    TH2F*                  hEPDiffV0CR0_V0CR3;
    TH2F*                  hEPDiffV0C_V0AR4;
    TH2F*                  hEPDiffV0C_V0AR7;
    TH2F*                  hEPDiffV0AR4_V0AR7;
    TProfile2D*            hEPRbrCosV0A;
    TProfile2D*            hEPRbrSinV0A;
    TProfile2D*            hEPRbrCosV0C;
    TProfile2D*            hEPRbrSinV0C;
    TProfile2D*            hEPRbrCosTPC;
    TProfile2D*            hEPRbrSinTPC;

    TH2F*                  hV2ClustDxDzA;
    TH2F*                  hV2ClustDxDzB;
    TH3F*                  hV2ClustDphiV0A;
    TH3F*                  hV2ClustDphiV0C;
    TH3F*                  hV2ClustCos2phiV0A;
    TH3F*                  hV2ClustCos2phiV0C;
    TH2F*                  hV1ClustDxDzA;
    TH2F*                  hV1ClustDxDzB;
    TH2F*                  hM02EA;
    TH2F*                  hM02EB;
    TH1F*                  hV1ClustNlmA;
    TH1F*                  hV1ClustNlmB;

    TH2F*                  hTrkPhiEta;
    TH1F*                  hTrkPt;
    TH3F*                  hTrkDphiEmcV0A;
    TH3F*                  hTrkDphiEmcV0C;
    TH3F*                  hTrkCos2phiEmcV0A;
    TH3F*                  hTrkCos2phiEmcV0C;
    TH3F*                  hTrkDphiOutEmcV0A;
    TH3F*                  hTrkDphiOutEmcV0C;
    TH3F*                  hTrkCos2phiOutEmcV0A;
    TH3F*                  hTrkCos2phiOutEmcV0C;

    THnSparse*             fV2ClusterV0A;
    THnSparse*             fV2ClusterV0C;
    THnSparse*             fV2ClusterTPC;
    THnSparse*             fV1ClusterV0A;
    THnSparse*             fV1ClusterV0C;
    THnSparse*             fV1ClusterTPC;

    AliAnalysisTaskPi0V2(const AliAnalysisTaskPi0V2&); // not implemented
    AliAnalysisTaskPi0V2& operator=(const AliAnalysisTaskPi0V2&); // not implemented
    
    ClassDef(AliAnalysisTaskPi0V2, 6); // example of analysis
};
#endif
