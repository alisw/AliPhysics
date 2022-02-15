#include "AliTrackInfoContainer.h"

class AliTrackInfoContainer;    // your analysis class

ClassImp(AliTrackInfoContainer) // classimp: necessary for root

//_____________________________________________________________________________
AliTrackInfoContainer::AliTrackInfoContainer() :
TObject                     (), 
//===========================================================================================
//   Variables for Primary tracks Tree
//===========================================================================================
//-----------BASIC-INFO---------------------------
fTreePrimVarCharge(0),
fTreePrimVarRapPion(0), 
fTreePrimVarRapProton(0), 
fTreePrimVarRapKaon(0),
fTreePrimVarEta(0),
fTreePrimVarPtot(0), 
fTreePrimVarPt(0), 

fTreePrimVarDCAxyToPV(0),
fTreePrimVarDCAzToPV(0),

fTreePrimVarNbrCrossedRows(0),
fTreePrimVarITSNbrClusters(0),
fTreePrimVarTPCNbrClusters(0),
fTreePrimVarRatioCrossedRowsOverFindable(0),
fTreePrimVarNbrCrossedRowsOverLength(0),
fTreePrimVarFractionSharedTPCClusters(0),
fTreePrimVarChi2perNDF(0),
fTreePrimVarITSChi2PerCluster(0),
fTreePrimVarTPCChi2PerCluster(0),
fTreePrimVarChi2TPCConstrainedVsGlobal(0),
fTreePrimVarTrackLength(0),
//---------------PID-TPC-INFO---------------------
fTreePrimVarNSigmaPion(0),
fTreePrimVarNSigmaKaon(0),
fTreePrimVarNSigmaProton(0),
//---------------PID-TOF-INFO---------------------
fTreePrimVarTOFNSigmaPion(0),
fTreePrimVarTOFNSigmaKaon(0),
fTreePrimVarTOFNSigmaProton(0),
//---------------PID-ITS-INFO---------------------
fTreePrimVarITSNSigmaPion(0),
fTreePrimVarITSNSigmaKaon(0),
fTreePrimVarITSNSigmaProton(0),
//---Raw TPC dEdx + PIDForTracking information----
fTreePrimVardEdx(0),
fTreePrimVardEdxN(0),
fTreePrimVarPIDForTracking(0),
//------------------------------------------------
fTreePrimVarTrackStatus(0),
fTreePrimVarFilterMap(0),
//------------------------------------------------
fTreePrimVarLabel(-1),
//------------FULL-MOMENTUM-INFO------------------
fTreePrimVarPx(0),
fTreePrimVarPy(0),
fTreePrimVarPz(0),
//---------------CLUSTER-INFO---------------------
fTreePrimVarITSClusters0(0),
fTreePrimVarITSClusters1(0),
fTreePrimVarITSClusters2(0),
fTreePrimVarITSClusters3(0),
fTreePrimVarITSClusters4(0),
fTreePrimVarITSClusters5(0),
//------------------------------------------------
fTreePrimVarITSSharedClusters0(0),
fTreePrimVarITSSharedClusters1(0),
fTreePrimVarITSSharedClusters2(0),
fTreePrimVarITSSharedClusters3(0),
fTreePrimVarITSSharedClusters4(0),
fTreePrimVarITSSharedClusters5(0),
//---------------OOB-PILEUP-INFO---------------------
fTreePrimVarTOFExpTDiff(0),
fTreePrimVarTOFSignal(0),
fTreePrimVarTOFBCid(0),

fTreePrimVarTOFLength(0),
fTreePrimVarTOFDeltaX(0),
fTreePrimVarTOFDeltaZ(0),
//Kink tagging
fTreePrimVarIsKink(kFALSE),

fTreePrimVarPassesTrackCuts2010(kFALSE),
fTreePrimVarPassesTrackCuts2011(kFALSE),
fTreePrimVarPassesTrackCutsTPCRefit(kFALSE),
fTreePrimVarPassesTrackCuts2011Sys(kFALSE),
fTreePrimVarPassesTrackCutsV0(kFALSE),
fTreePrimVarPassesTrackCutsHybrid_kOff(kFALSE),
fTreePrimVarPassesTrackCutsHybrid_kNone(kFALSE),

fTreePrimVarPosX(-999),
fTreePrimVarPosY(-999),
fTreePrimVarPosZ(-999)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliTrackInfoContainer::~AliTrackInfoContainer() // destructor
{
    // destructor
}
//_____________________________________________________________________________
void AliTrackInfoContainer::Reset()
{
    //===========================================================================================
    //   Variables for Primary tracks Tree
    //===========================================================================================
    //-----------BASIC-INFO---------------------------
    fTreePrimVarCharge = 0;
    fTreePrimVarRapPion = 0; 
    fTreePrimVarRapProton = 0; 
    fTreePrimVarRapKaon = 0;
    fTreePrimVarEta = 0;
    fTreePrimVarPtot = 0; 
    fTreePrimVarPt = 0; 

    fTreePrimVarDCAxyToPV = 0;
    fTreePrimVarDCAzToPV = 0;

    fTreePrimVarNbrCrossedRows = 0;
    fTreePrimVarITSNbrClusters = 0;
    fTreePrimVarTPCNbrClusters = 0;
    fTreePrimVarRatioCrossedRowsOverFindable = 0;
    fTreePrimVarNbrCrossedRowsOverLength = 0;
    fTreePrimVarFractionSharedTPCClusters = 0;
    fTreePrimVarChi2perNDF  = 0 ;
    fTreePrimVarITSChi2PerCluster = 0;
    fTreePrimVarTPCChi2PerCluster = 0;
    fTreePrimVarChi2TPCConstrainedVsGlobal = 0;
    fTreePrimVarTrackLength = 0;
    //---------------PID-TPC-INFO---------------------
    fTreePrimVarNSigmaPion = 0;
    fTreePrimVarNSigmaKaon = 0;
    fTreePrimVarNSigmaProton = 0;
    //---------------PID-TOF-INFO---------------------
    fTreePrimVarTOFNSigmaPion = 0;
    fTreePrimVarTOFNSigmaKaon = 0;
    fTreePrimVarTOFNSigmaProton = 0;
    //---------------PID-ITS-INFO---------------------
    fTreePrimVarITSNSigmaPion = 0;
    fTreePrimVarITSNSigmaKaon = 0;
    fTreePrimVarITSNSigmaProton = 0;
    //---Raw TPC dEdx + PIDForTracking information----
    fTreePrimVardEdx = 0;
    fTreePrimVardEdxN = 0;
    fTreePrimVarPIDForTracking = 0;
    //------------------------------------------------
    fTreePrimVarTrackStatus = 0;
	fTreePrimVarFilterMap = 0;
    //------------------------------------------------
    fTreePrimVarLabel = -1;
    //------------FULL-MOMENTUM-INFO------------------
    fTreePrimVarPx = 0;
    fTreePrimVarPy = 0;
    fTreePrimVarPz = 0;
    //---------------CLUSTER-INFO---------------------
    fTreePrimVarITSClusters0 = 0;
    fTreePrimVarITSClusters1 = 0;
    fTreePrimVarITSClusters2 = 0;
    fTreePrimVarITSClusters3 = 0;
    fTreePrimVarITSClusters4 = 0;
    fTreePrimVarITSClusters5 = 0;
    //------------------------------------------------
    fTreePrimVarITSSharedClusters0 = 0;
    fTreePrimVarITSSharedClusters1 = 0;
    fTreePrimVarITSSharedClusters2 = 0;
    fTreePrimVarITSSharedClusters3 = 0;
    fTreePrimVarITSSharedClusters4 = 0;
    fTreePrimVarITSSharedClusters5 = 0;
    //---------------OOB-PILEUP-INFO---------------------
    fTreePrimVarTOFExpTDiff = 0;
    fTreePrimVarTOFSignal = 0;
    fTreePrimVarTOFBCid = 0;

    fTreePrimVarTOFLength = 0;
    fTreePrimVarTOFDeltaX = 0;
    fTreePrimVarTOFDeltaZ = 0;
    //Kink tagging
    fTreePrimVarIsKink = kFALSE;
    
    fTreePrimVarPassesTrackCuts2010         = kFALSE;
    fTreePrimVarPassesTrackCuts2011         = kFALSE;
    fTreePrimVarPassesTrackCutsTPCRefit     = kFALSE;
    fTreePrimVarPassesTrackCuts2011Sys      = kFALSE;
    fTreePrimVarPassesTrackCutsV0           = kFALSE;
    fTreePrimVarPassesTrackCutsHybrid_kOff  = kFALSE;
    fTreePrimVarPassesTrackCutsHybrid_kNone = kFALSE;
    
    fTreePrimVarPosX = -999;
    fTreePrimVarPosY = -999;
    fTreePrimVarPosZ = -999;
}
//_____________________________________________________________________________
void AliTrackInfoContainer::Fill(const AliAODTrack* primTrack, const AliPIDResponse* fPIDResponse, Double_t lBestPrimaryVtxPos[3], Double_t lMagField) 
{
    //------------------------------------------------
    // Initializations
    //------------------------------------------------
    if (!primTrack ) 
    {
        AliWarning("ERROR: Could not retrieve AOD primary track ...");
        return;
    }
    
    //________________________________________________________________________
    // Primary track cuts
    
    Double_t lPx = 0.; 
    Double_t lPy = 0.; 
    Double_t lPz = 0.;
    
    lPx  = primTrack->Px() ; 
    lPy  = primTrack->Py() ; 
    lPz  = primTrack->Pz() ;
    
    //[0] = DCAxy ; [1] = DCAz
    Float_t lDCA[2]    ; Float_t lCovDCA[3]; 
    primTrack->GetImpactParameters(lDCA, lCovDCA);
    
    //------------------------------------------------
    // Calculation of the variables related to primaries
    //------------------------------------------------
    fTreePrimVarCharge          = primTrack->Charge();
    fTreePrimVarRapPion         = primTrack->Y( AliAODTrack::kPion );
    fTreePrimVarRapProton       = primTrack->Y( AliAODTrack::kProton );
    fTreePrimVarRapKaon         = primTrack->Y( AliAODTrack::kKaon );
    fTreePrimVarEta             = primTrack->Eta();
    fTreePrimVarPtot            = primTrack->P();
    fTreePrimVarPt              = primTrack->Pt();
    
    fTreePrimVarDCAxyToPV       = lDCA[0];
    fTreePrimVarDCAzToPV        = lDCA[1];

    //________________________________________________________________________
    // Track quality cuts
    //Least Nbr of Crossed Rows
	fTreePrimVarNbrCrossedRows = primTrack->GetTPCNCrossedRows();
    
    // Least Nbr of clusters      
    fTreePrimVarITSNbrClusters   = primTrack->GetITSNcls();
    fTreePrimVarTPCNbrClusters   = primTrack->GetTPCNcls();
    
    //Compute ratio Crossed Rows / Findable clusters
    //Note: this test avoids division by zero!
	if( primTrack->GetTPCNclsF() > 0 )
		fTreePrimVarRatioCrossedRowsOverFindable = fTreePrimVarNbrCrossedRows / ((Double_t)primTrack->GetTPCNclsF());
    
    //Min track length
	fTreePrimVarTrackLength   = GetLengthInActiveZone(primTrack, 2.0, 220.0, lMagField);
    
    //Nbr Of Crossed Rows Over Length
    fTreePrimVarNbrCrossedRowsOverLength = fTreePrimVarNbrCrossedRows/(fTreePrimVarTrackLength + 1e-5);
    
    // Fraction of shared TPC clusters
    fTreePrimVarFractionSharedTPCClusters = Float_t(primTrack->GetTPCnclsS())/Float_t(fTreePrimVarTPCNbrClusters + 1e-5);
    
    //Track chi2
    fTreePrimVarChi2perNDF = primTrack->Chi2perNDF();
    
    //ITS Chi2 per cluster
	if( fTreePrimVarITSNbrClusters > 0 )
		fTreePrimVarITSChi2PerCluster = primTrack->GetITSchi2() / ((Float_t) fTreePrimVarITSNbrClusters);
    
    //TPC Chi2 per cluster
	if( fTreePrimVarTPCNbrClusters > 0 )
		fTreePrimVarTPCChi2PerCluster = primTrack->GetTPCchi2() / ((Float_t) fTreePrimVarTPCNbrClusters);
    
    //Chi2 TPC Constrained Vs Global
    fTreePrimVarChi2TPCConstrainedVsGlobal = primTrack->GetChi2TPCConstrainedVsGlobal();
    
    //------------------------------------------------
    // TPC dEdx information
    //------------------------------------------------
    fTreePrimVarNSigmaPion      = fPIDResponse->NumberOfSigmasTPC( primTrack, AliPID::kPion );
    fTreePrimVarNSigmaKaon      = fPIDResponse->NumberOfSigmasTPC( primTrack, AliPID::kKaon );
    fTreePrimVarNSigmaProton    = fPIDResponse->NumberOfSigmasTPC( primTrack, AliPID::kProton );
    
    //------------------------------------------------
    // ITS info 
    //------------------------------------------------
    fTreePrimVarITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( primTrack, AliPID::kPion );
    fTreePrimVarITSNSigmaKaon   = fPIDResponse->NumberOfSigmasITS( primTrack, AliPID::kKaon );
    fTreePrimVarITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( primTrack, AliPID::kProton );
    
    //------------------------------------------------
    // TOF info 
    //------------------------------------------------
    fTreePrimVarTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( primTrack, AliPID::kPion );
    fTreePrimVarTOFNSigmaKaon   = fPIDResponse->NumberOfSigmasTOF( primTrack, AliPID::kKaon );
    fTreePrimVarTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( primTrack, AliPID::kProton );
    
    //------------------------------------------------
    // Raw TPC dEdx + PIDForTracking information
    //------------------------------------------------
    //Acquire TPC Signals
    fTreePrimVardEdx    = primTrack->GetTPCsignal();
    fTreePrimVardEdxN   = primTrack->GetTPCsignalN();
    //Acquire PID For Tracking
    fTreePrimVarPIDForTracking = primTrack->GetPIDForTracking();
    
    //________________________________________________________________________
    // Track status
    fTreePrimVarTrackStatus     = primTrack->GetStatus();
	fTreePrimVarFilterMap		= primTrack->GetFilterMap();
    
    //________________________________________________________________________
    // Momentum info
    fTreePrimVarPx = lPx ;  fTreePrimVarPy = lPy  ;  fTreePrimVarPz = lPz  ; 
    
    //________________________________________________________________________
    // Particle label
    fTreePrimVarLabel = primTrack->GetLabel();
    
    //________________________________________________________________________
    //Check clusters
    fTreePrimVarITSClusters0 = primTrack->HasPointOnITSLayer(0);
    fTreePrimVarITSClusters1 = primTrack->HasPointOnITSLayer(1);
    fTreePrimVarITSClusters2 = primTrack->HasPointOnITSLayer(2);
    fTreePrimVarITSClusters3 = primTrack->HasPointOnITSLayer(3);
    fTreePrimVarITSClusters4 = primTrack->HasPointOnITSLayer(4);
    fTreePrimVarITSClusters5 = primTrack->HasPointOnITSLayer(5);

    //________________________________________________________________________
    //Check its clusters, shared
    fTreePrimVarITSSharedClusters0 = primTrack->HasSharedPointOnITSLayer(0);
    fTreePrimVarITSSharedClusters1 = primTrack->HasSharedPointOnITSLayer(1);
    fTreePrimVarITSSharedClusters2 = primTrack->HasSharedPointOnITSLayer(2);
    fTreePrimVarITSSharedClusters3 = primTrack->HasSharedPointOnITSLayer(3);
    fTreePrimVarITSSharedClusters4 = primTrack->HasSharedPointOnITSLayer(4);
    fTreePrimVarITSSharedClusters5 = primTrack->HasSharedPointOnITSLayer(5);
    
    //________________________________________________________________________
    //GetKinkIndex condition
    fTreePrimVarIsKink = kFALSE;
    if( primTrack->GetProdVertex()->GetType() == AliAODVertex::kKink ) fTreePrimVarIsKink = kTRUE;
    
    //________________________________________________________________________
    //TOF info
    fTreePrimVarTOFExpTDiff     = primTrack->GetTOFExpTDiff( lMagField , kTRUE);
    fTreePrimVarTOFSignal       = primTrack->GetTOFsignal() * 1.e-3; // in ns
    fTreePrimVarTOFBCid         = primTrack->GetTOFBunchCrossing( lMagField);
    
    fTreePrimVarTOFLength       = primTrack->GetIntegratedLength();
    fTreePrimVarTOFDeltaX       = primTrack->GetTOFsignalDx();
    fTreePrimVarTOFDeltaZ       = primTrack->GetTOFsignalDz();
    
    Double_t *p = new Double_t[3];
    primTrack->GetXYZ(p);
    fTreePrimVarPosX = p[0];
    fTreePrimVarPosY = p[1];
    fTreePrimVarPosZ = p[2];
    delete[] p;
}

//_____________________________________________________________________________
Float_t AliTrackInfoContainer::GetLengthInActiveZone(const AliAODTrack *gt, Float_t deltaY, Float_t deltaZ, Float_t b ){
    AliESDtrack esdTrack( gt );
    esdTrack.SetESDEvent((AliESDEvent*) gt->GetEvent() );
    AliExternalTrackParam etp;
    etp.CopyFromVTrack(gt);
    esdTrack.ResetTrackParamIp(&etp);
    return esdTrack.GetLengthInActiveZone(1, deltaY, deltaZ, b);
}
//_____________________________________________________________________________
Float_t AliTrackInfoContainer::GetDCAz(AliAODTrack *lTrack)
//Encapsulation of DCAz calculation
{
    Float_t b[2];
    Float_t bCov[3];
    lTrack->GetImpactParameters(b,bCov);
    if (bCov[0]<=0 || bCov[2]<=0) {
        AliDebug(1, "Estimated b resolution lower or equal to zero!");
        bCov[0]=0; bCov[2]=0;
    }
    //Float_t dcaToVertexXY = b[0];
    Float_t dcaToVertexZ = b[1];
    
    return dcaToVertexZ;
}
//_____________________________________________________________________________
Bool_t AliTrackInfoContainer::IsSelected( const AliESDtrackCuts* lPrimTrackCuts )
{
    Double_t lCutMinDCAToVertexXY   = 0.;
    Double_t lCutMinDCAToVertexZ    = 0.;   
    Double_t lCutMaxDCAToVertexXY   = 0.;
    Double_t lCutMaxDCAToVertexZ    = 0.;
    
    lCutMinDCAToVertexXY    = lPrimTrackCuts->GetMinDCAToVertexXY();
    lCutMinDCAToVertexZ     = lPrimTrackCuts->GetMinDCAToVertexZ();
    lCutMaxDCAToVertexXY    = lPrimTrackCuts->GetMaxDCAToVertexXY();
    lCutMaxDCAToVertexZ     = lPrimTrackCuts->GetMaxDCAToVertexZ();
    
    TString lMinDCAXYFormula    ( lPrimTrackCuts->GetMinDCAToVertexXYPtDep()    );
    TString lMinDCAZFormula     ( lPrimTrackCuts->GetMinDCAToVertexZPtDep()     );
    TString lMaxDCAXYFormula    ( lPrimTrackCuts->GetMaxDCAToVertexXYPtDep()    );
    TString lMaxDCAZFormula     ( lPrimTrackCuts->GetMaxDCAToVertexZPtDep()     );
    
    if( lMinDCAXYFormula.CompareTo("") )
    {
        lMinDCAXYFormula.ReplaceAll("pt","x");
        TFormula lCutMinDCAToVertexXYPtDep("lCutMinDCAToVertexXYPtDep", lMinDCAXYFormula.Data());
        
        lCutMinDCAToVertexXY = lCutMinDCAToVertexXYPtDep.Eval(fTreePrimVarPt);
    }
    if( lMinDCAZFormula.CompareTo("") )
    {
        lMinDCAZFormula.ReplaceAll("pt","x");
        TFormula lCutMinDCAToVertexZPtDep("lCutMinDCAToVertexZPtDep", lMinDCAZFormula.Data());
        
        lCutMinDCAToVertexZ = lCutMinDCAToVertexZPtDep.Eval(fTreePrimVarPt);
    }
    if( lMaxDCAXYFormula.CompareTo(""))
    {
        lMaxDCAXYFormula.ReplaceAll("pt","x");
        TFormula lCutMaxDCAToVertexXYPtDep("lCutMaxDCAToVertexXYPtDep", lMaxDCAXYFormula.Data());
        
        lCutMaxDCAToVertexXY = lCutMaxDCAToVertexXYPtDep.Eval(fTreePrimVarPt);
    }
    if( lMaxDCAZFormula.CompareTo("") )
    {
        lMaxDCAZFormula.ReplaceAll("pt","x");
        TFormula lCutMaxDCAToVertexZPtDep("lCutMaxDCAToVertexZPtDep", lMaxDCAZFormula.Data());
        
        lCutMaxDCAToVertexZ = lCutMaxDCAToVertexZPtDep.Eval(fTreePrimVarPt);
    }
    
    Float_t lCutMinPtTracks     = 0.;
    Float_t lCutMaxPtTracks     = 0.;
    
    Float_t lCutMinEtaTracks    = 0.;
    Float_t lCutMaxEtaTracks    = 0.;
    
    Float_t lCutMinRapTracks    = 0.;
    Float_t lCutMaxRapTracks    = 0.;
    
    lPrimTrackCuts->GetPtRange  (lCutMinPtTracks    , lCutMaxPtTracks   );
    lPrimTrackCuts->GetEtaRange (lCutMinEtaTracks   , lCutMaxEtaTracks  );
    lPrimTrackCuts->GetRapRange (lCutMinRapTracks   , lCutMaxRapTracks  );
    
    if (
        //Check 1: Charge selection
        fTreePrimVarCharge                            != 0                                                          &&
        
        //Check 2: Transverse momentum selection
        lCutMinPtTracks     < fTreePrimVarPt    && fTreePrimVarPt   < lCutMaxPtTracks                               && 
        
        //Check 3: Basic Acceptance cuts
        // Pseudo-rapidity
        lCutMinEtaTracks    < fTreePrimVarEta   && fTreePrimVarEta  < lCutMaxEtaTracks                              && 
        // Rapidity
        (
            ( lCutMinRapTracks  < fTreePrimVarRapPion   && fTreePrimVarRapPion      < lCutMaxRapTracks ) ||
            ( lCutMinRapTracks  < fTreePrimVarRapKaon   && fTreePrimVarRapKaon      < lCutMaxRapTracks ) ||
            ( lCutMinRapTracks  < fTreePrimVarRapProton && fTreePrimVarRapProton    < lCutMaxRapTracks ) 
        )                                                                                                           &&

        //Check 4: Topological Variables
        (   
            (// if DCAToVertex2D is off (default)
                !lPrimTrackCuts->GetDCAToVertex2D() && 
                (  // Classical cut on the DCAxy 
                    TMath::Abs( fTreePrimVarDCAxyToPV ) > lCutMinDCAToVertexXY   && 
                    TMath::Abs( fTreePrimVarDCAxyToPV ) < lCutMaxDCAToVertexXY   &&
                    // Classical cut on the DCAz 
                    TMath::Abs( fTreePrimVarDCAzToPV )  > lCutMinDCAToVertexZ    && 
                    TMath::Abs( fTreePrimVarDCAzToPV )  < lCutMaxDCAToVertexZ
                )
            ) 
            ||
            (// if DCAToVertex2D is on 
                lPrimTrackCuts->GetDCAToVertex2D()  &&
                (  // sqrt( (DCAxy/MaxCutDCAxy)**2 + (DCAz/MaxCutDCAz)**2 ) > 1 IF MaxCutDCAxy & MaxCutDCAz > 0
                    ( 
                        ( lCutMaxDCAToVertexXY <= 0 || lCutMaxDCAToVertexZ <= 0 ) ||
                        TMath::Sqrt( TMath::Power(fTreePrimVarDCAxyToPV/lCutMaxDCAToVertexXY, 2) + 
                                        TMath::Power(fTreePrimVarDCAzToPV /lCutMaxDCAToVertexZ , 2) ) > 1 
                    ) &&
                    // sqrt( (DCAxy/MinCutDCAxy)**2 + (DCAz/MinCutDCAz)**2 ) < 1 IF MinCutDCAxy & MinCutDCAz > 0
                    ( 
                        ( lCutMinDCAToVertexXY <= 0 || lCutMinDCAToVertexZ <= 0 ) ||
                        TMath::Sqrt( TMath::Power(fTreePrimVarDCAxyToPV/lCutMinDCAToVertexXY, 2) + 
                                        TMath::Power(fTreePrimVarDCAzToPV /lCutMinDCAToVertexZ , 2) ) < 1 
                    ) 
                )
            )                                                                                                       
        )                                                                                                            &&
        
        // - Miscellaneous
        fTreePrimVarNbrCrossedRows                >= lPrimTrackCuts->GetMinNCrossedRowsTPC()                        &&
        //fTreePrimVarRatioCrossedRowsOverFindable  >= lPrimTrackCuts->GetMinRatioCrossedRowsOverFindableClustersTPC() &&
        fTreePrimVarTrackLength                   > lPrimTrackCuts->GetMinLengthActiveVolumeTPC()                   &&
        fTreePrimVarFractionSharedTPCClusters     < lPrimTrackCuts->GetMaxFractionSharedTPCClusters()               &&
        
        //Check 5: TPC dEdx selections
        (
			TMath::Abs(fTreePrimVarNSigmaPion)        < 5.  ||
			TMath::Abs(fTreePrimVarNSigmaKaon)        < 5.  ||
			TMath::Abs(fTreePrimVarNSigmaProton)      < 5.                                                               
		)																											&&

        //Check 6: Max Chi2/Clusters 
        fTreePrimVarITSChi2PerCluster             < lPrimTrackCuts->GetMaxChi2PerClusterITS()                       &&
        fTreePrimVarTPCChi2PerCluster             < lPrimTrackCuts->GetMaxChi2PerClusterTPC()                       &&
        
        //Check 7 : Chi2 TPC Constrained Vs Global cut
        fTreePrimVarChi2TPCConstrainedVsGlobal    <= lPrimTrackCuts->GetMaxChi2TPCConstrainedGlobal()               &&
        
        //Check 8: kITSrefit track selection if requested
        (   //do not need ITS refit tracks --> skip this cut
            !( lPrimTrackCuts->GetRequireITSRefit() )                           ||
            //otherwise : need to have ITS refit tracks
            ( fTreePrimVarTrackStatus    & AliAODTrack::kITSrefit ) 
        )                                                                                                           &&

        //Check 9: Kink rejection
        (   //do not reject kink topology            OR  reject them if required
            lPrimTrackCuts->GetAcceptKinkDaughters() || fTreePrimVarIsKink == kFALSE 
        )                                                                                                           &&
        
        //Check 10 : ITS Cluster requirements
        ( 
            CheckITSClusterRequirement( lPrimTrackCuts->GetClusterRequirementITS( AliESDtrackCuts::kSPD ),      
                                        fTreePrimVarITSClusters0, 
                                        fTreePrimVarITSClusters1)  &&
            
            CheckITSClusterRequirement( lPrimTrackCuts->GetClusterRequirementITS( AliESDtrackCuts::kSDD ),
                                        fTreePrimVarITSClusters2, 
                                        fTreePrimVarITSClusters3)  &&
            
            CheckITSClusterRequirement( lPrimTrackCuts->GetClusterRequirementITS( AliESDtrackCuts::kSSD ),
                                        fTreePrimVarITSClusters4, 
                                        fTreePrimVarITSClusters5) 
        )                                                                                                          
        
    )//end major if
    {
        return kTRUE;
    }
    
    return kFALSE;
}
//_____________________________________________________________________________
Bool_t AliTrackInfoContainer::CheckITSClusterRequirement(AliESDtrackCuts::ITSClusterRequirement req, Bool_t clusterL1, Bool_t clusterL2)
{
    //
    // Encapsulation of ITS cluster requirement 
    // (from https://github.com/alisw/AliRoot/blob/master/ANALYSIS/ANALYSISalice/AliESDtrackCuts.cxx#L2268)
    //
    // checks if the cluster requirement is fullfilled (in this case: return kTRUE)
    //
    switch (req)
    {
        case AliESDtrackCuts::kOff          :  return kTRUE;
        case AliESDtrackCuts::kNone         :  return !clusterL1 && !clusterL2;
        case AliESDtrackCuts::kAny          :  return clusterL1 || clusterL2;
        case AliESDtrackCuts::kFirst        :  return clusterL1;
        case AliESDtrackCuts::kOnlyFirst    :  return clusterL1 && !clusterL2;
        case AliESDtrackCuts::kSecond       :  return clusterL2;
        case AliESDtrackCuts::kOnlySecond   :  return clusterL2 && !clusterL1;
        case AliESDtrackCuts::kBoth         :  return clusterL1 && clusterL2;
    }

    return kFALSE;
}
//_____________________________________________________________________________
