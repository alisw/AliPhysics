#include "AliMCV0Container.h"

class AliMCV0Container;    // your analysis class

ClassImp(AliMCV0Container) // classimp: necessary for root

//_____________________________________________________________________________
AliMCV0Container::AliMCV0Container() :
TObject                     (), 
//===========================================================================================
//   Variables for V0 Tree
//===========================================================================================
//-----------BASIC-INFO---------------------------
fTreeV0VarOnFlyStatus(kFALSE),
fTreeV0VarMassAsK0s(0), 
fTreeV0VarMassAsLambda(0), 
fTreeV0VarMassAsAntiLambda(0), 
fTreeV0VarRapK0Short(0), 
fTreeV0VarRapLambda(0), 
fTreeV0VarV0Eta(0), 
fTreeV0VarPosEta(0), 
fTreeV0VarNegEta(0), 
fTreeV0VarPtot(0), 
fTreeV0VarPt(0), 
//-----------INFO-FOR-CUTS------------------------
fTreeV0VarAlpha(0), 
fTreeV0VarPtArm(0),
fTreeV0VarDCAV0Dau(0), 
fTreeV0VarDCAV0ToPV(0), 
fTreeV0VarDCAPosToPV(0), 
fTreeV0VarDCANegToPV(0), 
fTreeV0VarCosPA(0), 
fTreeV0VarRadius(0), 

fTreeV0VarLeastNbrCrossedRows(0),
fTreeV0VarLeastNbrClusters(0),
fTreeV0VarLeastRatioCrossedRowsOverFindable(0),
fTreeV0VarMaxChi2PerCluster(0), 
fTreeV0VarMinTrackLength(0), 
//-----------DECAY-LENGTH-INFO--------------------
fTreeV0VarDistOverTotMom(0),
//---------------PID-TPC-INFO---------------------
fTreeV0VarPosNSigmaProton(0),
fTreeV0VarPosNSigmaPion(0),
fTreeV0VarNegNSigmaProton(0),
fTreeV0VarNegNSigmaPion(0),
//---------------PID-TOF-INFO---------------------
fTreeV0VarPosTOFNSigmaProton(0),
fTreeV0VarPosTOFNSigmaPion(0),
fTreeV0VarNegTOFNSigmaProton(0),
fTreeV0VarNegTOFNSigmaPion(0),
//---------------PID-ITS-INFO---------------------
fTreeV0VarPosITSNSigmaProton(0),
fTreeV0VarPosITSNSigmaPion(0),
fTreeV0VarNegITSNSigmaProton(0),
fTreeV0VarNegITSNSigmaPion(0),
//---Raw TPC dEdx + PIDForTracking information----
fTreeV0VarPosdEdx(0),
fTreeV0VarNegdEdx(0),
fTreeV0VarPosdEdxN(0),
fTreeV0VarNegdEdxN(0),
fTreeV0VarPosPIDForTracking(0), 
fTreeV0VarNegPIDForTracking(0),
//------------------------------------------------
fTreeV0VarChi2V0(0),         
//------------------------------------------------
fTreeV0VarPosTrackStatus(0), 
fTreeV0VarNegTrackStatus(0), 
//------------------------------------------------
fTreeV0VarPosDCAz(0), 
fTreeV0VarNegDCAz(0), 
//------------FULL-MOMENTUM-INFO------------------
fTreeV0VarPosPx(0), 
fTreeV0VarPosPy(0), 
fTreeV0VarPosPz(0), 
fTreeV0VarNegPx(0), 
fTreeV0VarNegPy(0), 
fTreeV0VarNegPz(0), 
//------------------------------------------------
fTreeV0VarV0DecayX(0), 
fTreeV0VarV0DecayY(0), 
fTreeV0VarV0DecayZ(0), 
//------------------------------------------------
fTreeV0VarPosLabel(-1), 
fTreeV0VarNegLabel(-1), 
//---------------CLUSTER-INFO---------------------
fTreeV0VarPosITSClusters0(kFALSE),
fTreeV0VarPosITSClusters1(kFALSE),
fTreeV0VarPosITSClusters2(kFALSE),
fTreeV0VarPosITSClusters3(kFALSE),
fTreeV0VarPosITSClusters4(kFALSE),
fTreeV0VarPosITSClusters5(kFALSE),

fTreeV0VarNegITSClusters0(kFALSE),
fTreeV0VarNegITSClusters1(kFALSE),
fTreeV0VarNegITSClusters2(kFALSE),
fTreeV0VarNegITSClusters3(kFALSE),
fTreeV0VarNegITSClusters4(kFALSE),
fTreeV0VarNegITSClusters5(kFALSE),
//------------------------------------------------
fTreeV0VarPosITSSharedClusters0(kFALSE),
fTreeV0VarPosITSSharedClusters1(kFALSE),
fTreeV0VarPosITSSharedClusters2(kFALSE),
fTreeV0VarPosITSSharedClusters3(kFALSE),
fTreeV0VarPosITSSharedClusters4(kFALSE),
fTreeV0VarPosITSSharedClusters5(kFALSE),

fTreeV0VarNegITSSharedClusters0(kFALSE),
fTreeV0VarNegITSSharedClusters1(kFALSE),
fTreeV0VarNegITSSharedClusters2(kFALSE),
fTreeV0VarNegITSSharedClusters3(kFALSE),
fTreeV0VarNegITSSharedClusters4(kFALSE),
fTreeV0VarNegITSSharedClusters5(kFALSE),
//---------------OOB-PILEUP-INFO---------------------
fTreeV0VarNegTOFExpTDiff(0), 
fTreeV0VarPosTOFExpTDiff(0), 
fTreeV0VarNegTOFSignal(0), 
fTreeV0VarPosTOFSignal(0), 
fTreeV0VarNegTOFBCid(0), 
fTreeV0VarPosTOFBCid(0),  

fTreeV0VarPosTOFLength(0),
fTreeV0VarNegTOFLength(0),

fTreeV0VarPosTOFDeltaX(0),
fTreeV0VarNegTOFDeltaX(0),

fTreeV0VarPosTOFDeltaZ(0),
fTreeV0VarNegTOFDeltaZ(0),

//Kink tagging
fTreeV0VarPosIsKink(kFALSE),
fTreeV0VarNegIsKink(kFALSE),

//Cowboy/sailor studies
fTreeV0VarIsCowboy(kFALSE), 

//------------MONTE-CARLO-INFO------------------
fTreeV0VarPosEtaMC(-100),
fTreeV0VarNegEtaMC(-100),

fTreeV0VarPosPxMC(-100),
fTreeV0VarPosPyMC(-100),
fTreeV0VarPosPzMC(-100),
fTreeV0VarNegPxMC(-100),   
fTreeV0VarNegPyMC(-100),
fTreeV0VarNegPzMC(-100),

fTreeV0VarPosMotherLabel(-1),
fTreeV0VarNegMotherLabel(-1),

fTreeV0VarPosGrandMotherLabel(-1),
fTreeV0VarNegGrandMotherLabel(-1),

fTreeV0VarPosPdgCode(-9999),
fTreeV0VarNegPdgCode(-9999),

fTreeV0VarPosMotherPdgCode(-9999),
fTreeV0VarNegMotherPdgCode(-9999),

fTreeV0VarPosGrandMotherPdgCode(-9999),
fTreeV0VarNegGrandMotherPdgCode(-9999),

fTreeV0VarPdgCode(-9999),
fTreeV0VarRapMC(-100),
fTreeV0VarPtotMC(-100),        
fTreeV0VarPtMC(-100),  

fTreeV0VarMotherPdgCode(-9999),
fTreeV0VarMotherPt(-1),
fTreeV0VarMotherRapMC(-100),

fTreeV0VarDecayXMC(-100),
fTreeV0VarDecayYMC(-100),
fTreeV0VarDecayZMC(-100),

fTreeV0VarIsPhysicalPrimary(0),
fTreeV0VarMotherIsPhysicalPrimary(0),

fTreeV0VarIsPhysicalPrimaryPos(kFALSE),
fTreeV0VarIsPhysicalPrimaryNeg(kFALSE),

fTreeV0VarIsPhysicalPrimaryPosMother(kFALSE),
fTreeV0VarIsPhysicalPrimaryNegMother(kFALSE),

fTreeV0VarIsPhysicalPrimaryPosGrandMother(kFALSE),
fTreeV0VarIsPhysicalPrimaryNegGrandMother(kFALSE)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliMCV0Container::~AliMCV0Container() // destructor
{
    // destructor
}
//_____________________________________________________________________________
void AliMCV0Container::Reset()
{
    //-----------BASIC-INFO---------------------------
    fTreeV0VarOnFlyStatus = kFALSE;
    fTreeV0VarMassAsK0s = 0 ; 
    fTreeV0VarMassAsLambda = 0 ; 
    fTreeV0VarMassAsAntiLambda = 0 ; 
    fTreeV0VarRapK0Short = 0 ; 
    fTreeV0VarRapLambda = 0 ; 
    fTreeV0VarV0Eta = 0 ; 
    fTreeV0VarPosEta = 0 ; 
    fTreeV0VarNegEta = 0 ; 
    fTreeV0VarPtot = 0 ; 
    fTreeV0VarPt = 0 ; 
    //-----------INFO-FOR-CUTS------------------------
    fTreeV0VarAlpha = 0 ; 
    fTreeV0VarPtArm = 0 ;
    fTreeV0VarDCAV0Dau = 0 ; 
    fTreeV0VarDCAV0ToPV = 0 ; 
    fTreeV0VarDCAPosToPV = 0 ; 
    fTreeV0VarDCANegToPV = 0 ; 
    fTreeV0VarCosPA = 0 ; 
    fTreeV0VarRadius = 0 ; 

    fTreeV0VarLeastNbrCrossedRows = 0 ;
    fTreeV0VarLeastNbrClusters = 0 ;
    fTreeV0VarLeastRatioCrossedRowsOverFindable = 0 ;
    fTreeV0VarMaxChi2PerCluster = 0 ; 
    fTreeV0VarMinTrackLength = 0 ; 
    //-----------DECAY-LENGTH-INFO--------------------
    fTreeV0VarDistOverTotMom = 0 ;
    //---------------PID-TPC-INFO---------------------
    fTreeV0VarPosNSigmaProton = 0 ;
    fTreeV0VarPosNSigmaPion = 0 ;
    fTreeV0VarNegNSigmaProton = 0 ;
    fTreeV0VarNegNSigmaPion = 0 ;
    //---------------PID-TOF-INFO---------------------
    fTreeV0VarPosTOFNSigmaProton = 0 ;
    fTreeV0VarPosTOFNSigmaPion = 0 ;
    fTreeV0VarNegTOFNSigmaProton = 0 ;
    fTreeV0VarNegTOFNSigmaPion = 0 ;
    //---------------PID-ITS-INFO---------------------
    fTreeV0VarPosITSNSigmaProton = 0 ;
    fTreeV0VarPosITSNSigmaPion = 0 ;
    fTreeV0VarNegITSNSigmaProton = 0 ;
    fTreeV0VarNegITSNSigmaPion = 0 ;
    //---Raw TPC dEdx + PIDForTracking information----
    fTreeV0VarPosdEdx = 0 ;
    fTreeV0VarNegdEdx = 0 ;
    fTreeV0VarPosdEdxN = 0 ;
    fTreeV0VarNegdEdxN = 0 ;
    fTreeV0VarPosPIDForTracking = 0 ; 
    fTreeV0VarNegPIDForTracking = 0 ;
    //------------------------------------------------
    fTreeV0VarChi2V0 = 0 ;         
    //------------------------------------------------
    fTreeV0VarPosTrackStatus = 0 ; 
    fTreeV0VarNegTrackStatus = 0 ; 
    //------------------------------------------------
    fTreeV0VarPosDCAz = 0 ; 
    fTreeV0VarNegDCAz = 0 ; 
    //------------FULL-MOMENTUM-INFO------------------
    fTreeV0VarPosPx = 0 ; 
    fTreeV0VarPosPy = 0 ; 
    fTreeV0VarPosPz = 0 ; 
    fTreeV0VarNegPx = 0 ; 
    fTreeV0VarNegPy = 0 ; 
    fTreeV0VarNegPz = 0 ; 
    //------------------------------------------------
    fTreeV0VarV0DecayX = 0 ; 
    fTreeV0VarV0DecayY = 0 ; 
    fTreeV0VarV0DecayZ = 0 ; 
    //------------------------------------------------
    fTreeV0VarPosLabel = -1 ; 
    fTreeV0VarNegLabel = -1 ; 
    //---------------CLUSTER-INFO---------------------
    fTreeV0VarPosITSClusters0 = kFALSE ;
    fTreeV0VarPosITSClusters1 = kFALSE ;
    fTreeV0VarPosITSClusters2 = kFALSE ;
    fTreeV0VarPosITSClusters3 = kFALSE ;
    fTreeV0VarPosITSClusters4 = kFALSE ;
    fTreeV0VarPosITSClusters5 = kFALSE ;

    fTreeV0VarNegITSClusters0 = kFALSE ;
    fTreeV0VarNegITSClusters1 = kFALSE ;
    fTreeV0VarNegITSClusters2 = kFALSE ;
    fTreeV0VarNegITSClusters3 = kFALSE ;
    fTreeV0VarNegITSClusters4 = kFALSE ;
    fTreeV0VarNegITSClusters5 = kFALSE ;
    //------------------------------------------------
    fTreeV0VarPosITSSharedClusters0 = kFALSE ;
    fTreeV0VarPosITSSharedClusters1 = kFALSE ;
    fTreeV0VarPosITSSharedClusters2 = kFALSE ;
    fTreeV0VarPosITSSharedClusters3 = kFALSE ;
    fTreeV0VarPosITSSharedClusters4 = kFALSE ;
    fTreeV0VarPosITSSharedClusters5 = kFALSE ;

    fTreeV0VarNegITSSharedClusters0 = kFALSE ;
    fTreeV0VarNegITSSharedClusters1 = kFALSE ;
    fTreeV0VarNegITSSharedClusters2 = kFALSE ;
    fTreeV0VarNegITSSharedClusters3 = kFALSE ;
    fTreeV0VarNegITSSharedClusters4 = kFALSE ;
    fTreeV0VarNegITSSharedClusters5 = kFALSE ;
    //---------------OOB-PILEUP-INFO---------------------
    fTreeV0VarNegTOFExpTDiff = 0 ; 
    fTreeV0VarPosTOFExpTDiff = 0 ; 
    fTreeV0VarNegTOFSignal = 0 ; 
    fTreeV0VarPosTOFSignal = 0 ; 
    fTreeV0VarNegTOFBCid = 0 ; 
    fTreeV0VarPosTOFBCid = 0 ;  

    fTreeV0VarPosTOFLength = 0 ;
    fTreeV0VarNegTOFLength = 0 ;

    fTreeV0VarPosTOFDeltaX = 0 ;
    fTreeV0VarNegTOFDeltaX = 0 ;

    fTreeV0VarPosTOFDeltaZ = 0 ;
    fTreeV0VarNegTOFDeltaZ = 0 ;

    //Kink tagging
    fTreeV0VarPosIsKink = kFALSE ;
    fTreeV0VarNegIsKink = kFALSE ;

    //Cowboy/sailor studies
    fTreeV0VarIsCowboy = kFALSE ; 

    //------------MONTE-CARLO-INFO------------------
    fTreeV0VarPosEtaMC = -100 ;
    fTreeV0VarNegEtaMC = -100 ;

    fTreeV0VarPosPxMC = -100 ;
    fTreeV0VarPosPyMC = -100 ;
    fTreeV0VarPosPzMC = -100 ;
    fTreeV0VarNegPxMC = -100 ;   
    fTreeV0VarNegPyMC = -100 ;
    fTreeV0VarNegPzMC = -100 ;

    fTreeV0VarPosMotherLabel = -1 ;
    fTreeV0VarNegMotherLabel = -1 ;

    fTreeV0VarPosGrandMotherLabel = -1 ;
    fTreeV0VarNegGrandMotherLabel = -1 ;

    fTreeV0VarPosPdgCode = -9999 ;
    fTreeV0VarNegPdgCode = -9999 ;

    fTreeV0VarPosMotherPdgCode = -9999 ;
    fTreeV0VarNegMotherPdgCode = -9999 ;

    fTreeV0VarPosGrandMotherPdgCode = -9999 ;
    fTreeV0VarNegGrandMotherPdgCode = -9999 ;

    fTreeV0VarPdgCode = -9999 ;
    fTreeV0VarRapMC = -100 ;
    fTreeV0VarPtotMC = -100 ;        
    fTreeV0VarPtMC = -100 ;

    fTreeV0VarMotherPdgCode = -9999 ;
    fTreeV0VarMotherPt = -1 ;
    fTreeV0VarMotherRapMC = -100 ;

    fTreeV0VarDecayXMC = -100 ;
    fTreeV0VarDecayYMC = -100 ;
    fTreeV0VarDecayZMC = -100 ;

    fTreeV0VarIsPhysicalPrimary = 0 ;
    fTreeV0VarMotherIsPhysicalPrimary = 0 ;

    fTreeV0VarIsPhysicalPrimaryPos = kFALSE ;
    fTreeV0VarIsPhysicalPrimaryNeg = kFALSE ;

    fTreeV0VarIsPhysicalPrimaryPosMother = kFALSE ;
    fTreeV0VarIsPhysicalPrimaryNegMother = kFALSE ;

    fTreeV0VarIsPhysicalPrimaryPosGrandMother = kFALSE ;
    fTreeV0VarIsPhysicalPrimaryNegGrandMother = kFALSE ;
}
//_____________________________________________________________________________
void AliMCV0Container::Fill(AliAODv0* v0, AliPIDResponse* fPIDResponse, Double_t lBestPrimaryVtxPos[3], Double_t lMagField, AliMCEvent* lMCEvent) 
{
    //------------------------------------------------
    // Initializations
    //------------------------------------------------
    if(!v0)
    {
        printf("Error in AliMCV0Container::Fill : v0 is a null ptr ! \n");
        return;
    }
    
    Double_t lPosV0[3] = {0.,0.,0.}; // Position of V0 

    Double_t lV0MassAsK0s           = 0.;
    Double_t lV0MassAsLambda        = 0.;
    Double_t lV0MassAsAntiLambda    = 0.;
    
    //CheckChargeV0( v0 ); //FIXME this won't work for AOD as it is
    v0->GetXYZ(lPosV0);
    
    //Gather tracks informations
    AliAODTrack *pTrack = dynamic_cast<AliAODTrack*>( v0->GetDaughter(0) ); //0->Positive Daughter
    AliAODTrack *nTrack = dynamic_cast<AliAODTrack*>( v0->GetDaughter(1) ); //1->Negative Daughter
    
    if(!pTrack || !nTrack) 
    {
        AliWarning("ERROR: Could not retrieve one of the V0 daughter tracks");
        return;
    }
    
    //Test different mass hypothesis : K0s, Lambda, AntiLambda
    lV0MassAsK0s            = v0->MassK0Short();
    lV0MassAsLambda         = v0->MassLambda();
    lV0MassAsAntiLambda     = v0->MassAntiLambda();
    
    Double_t lPosPx = 0.; Double_t lPosPy = 0.; Double_t lPosPz = 0.;
    Double_t lNegPx = 0.; Double_t lNegPy = 0.; Double_t lNegPz = 0.;
    lPosPx  = pTrack->Px() ; lPosPy  = pTrack->Py() ; lPosPz  = pTrack->Pz() ;
    lNegPx  = nTrack->Px() ; lNegPy  = nTrack->Py() ; lNegPz  = nTrack->Pz() ;
    
    //------------------------------------------------
    // Calculation of the variables related to V0
    //------------------------------------------------
    fTreeV0VarOnFlyStatus           = v0->GetOnFlyStatus();
    fTreeV0VarChi2V0                = v0->Chi2V0();
    fTreeV0VarMassAsK0s             = lV0MassAsK0s;
    fTreeV0VarMassAsLambda          = lV0MassAsLambda;
    fTreeV0VarMassAsAntiLambda      = lV0MassAsAntiLambda;

    fTreeV0VarPtot                  = TMath::Sqrt( v0->Ptot2V0() );
    fTreeV0VarPt                    = TMath::Sqrt( v0->Pt2V0() );

    fTreeV0VarRapK0Short            = v0->RapK0Short();
    fTreeV0VarRapLambda             = v0->RapLambda();
    fTreeV0VarV0Eta                 = 0.5*TMath::Log(TMath::Abs((fTreeV0VarPtot + (lPosPz + lNegPz))/(fTreeV0VarPtot - (lPosPz + lNegPz) + 1e-13)));
    fTreeV0VarPosEta                = pTrack->Eta();
    fTreeV0VarNegEta                = nTrack->Eta();
    
    fTreeV0VarAlpha                 = v0->AlphaV0();
    fTreeV0VarPtArm                 = v0->PtArmV0();
    fTreeV0VarDCAV0Dau              = v0->DcaV0Daughters();
    fTreeV0VarDCAV0ToPV             = v0->DcaV0ToPrimVertex();
    fTreeV0VarDCAPosToPV            = v0->DcaPosToPrimVertex();
    fTreeV0VarDCANegToPV            = v0->DcaNegToPrimVertex();

    fTreeV0VarRadius                = TMath::Sqrt( lPosV0[0]*lPosV0[0]  +  lPosV0[1]*lPosV0[1] );
    fTreeV0VarCosPA                 = v0->CosPointingAngle(lBestPrimaryVtxPos);

    //________________________________________________________________________
    // Track quality cuts
    //Least Nbr of Crossed Rows
    Int_t lPosNbrCrossedRows  = pTrack->GetTPCNCrossedRows();
    Int_t lNegNbrCrossedRows  = nTrack->GetTPCNCrossedRows();
    
    Int_t lLeastNbrCrossedRows = (Int_t) lPosNbrCrossedRows;
    if( lNegNbrCrossedRows < lLeastNbrCrossedRows ) lLeastNbrCrossedRows = (Int_t) lNegNbrCrossedRows;
    
    fTreeV0VarLeastNbrCrossedRows = lLeastNbrCrossedRows;
    
    //Compute ratio Crossed Rows / Findable clusters
    Double_t lPosTrackCrossedRowsOverFindable = lPosNbrCrossedRows / ((Double_t)(pTrack->GetTPCNclsF() + 1e-5));
    Double_t lNegTrackCrossedRowsOverFindable = lNegNbrCrossedRows / ((Double_t)(nTrack->GetTPCNclsF() + 1e-5));
    
    fTreeV0VarLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
    if( lNegTrackCrossedRowsOverFindable < fTreeV0VarLeastRatioCrossedRowsOverFindable ) 
        fTreeV0VarLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;
    
    // Least Nbr of clusters        
    Int_t lPosTPCClusters        = 0.;
    Int_t lNegTPCClusters        = 0.;

    lPosTPCClusters   = pTrack->GetTPCNcls();
    lNegTPCClusters   = nTrack->GetTPCNcls();
    
    Int_t lLeastNbrOfClusters    = lPosTPCClusters;
    if( lNegTPCClusters < lLeastNbrOfClusters ) lLeastNbrOfClusters = lNegTPCClusters;
    
    fTreeV0VarLeastNbrClusters = lLeastNbrOfClusters;
    
    //Min track length
    Double_t lPosTrackLength        = 0.;
    Double_t lNegTrackLength        = 0.;
    
    lPosTrackLength    = GetLengthInActiveZone(pTrack, 2.0, 220.0, lMagField);
    lNegTrackLength    = GetLengthInActiveZone(nTrack, 2.0, 220.0, lMagField);
    
    Double_t lSmallestTrackLength = lPosTrackLength;
    if ( lNegTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lNegTrackLength;
    
    fTreeV0VarMinTrackLength = lSmallestTrackLength;
    
    //Max Chi2 per cluster
    Double_t lBiggestChi2PerCluster = 0.;
    Double_t lPosChi2PerCluster = pTrack->GetTPCchi2() / ((Float_t) lPosTPCClusters + 1e-13);
    Double_t lNegChi2PerCluster = nTrack->GetTPCchi2() / ((Float_t) lNegTPCClusters + 1e-13);
    
    if( lPosChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lPosChi2PerCluster;
    if( lNegChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lNegChi2PerCluster;

    fTreeV0VarMaxChi2PerCluster = lBiggestChi2PerCluster;
    
    //________________________________________________________________________
    // Decay length info
    fTreeV0VarDistOverTotMom  = TMath::Sqrt(   TMath::Power( lPosV0[0] - lBestPrimaryVtxPos[0] , 2) 
                                            + TMath::Power( lPosV0[1] - lBestPrimaryVtxPos[1] , 2) 
                                            + TMath::Power( lPosV0[2] - lBestPrimaryVtxPos[2] , 2) );
    fTreeV0VarDistOverTotMom /= fTreeV0VarPtot;
    
    
    //------------------------------------------------
    // TPC dEdx information
    //------------------------------------------------
    fTreeV0VarPosNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
    fTreeV0VarPosNSigmaProton = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
    fTreeV0VarNegNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion   );
    fTreeV0VarNegNSigmaProton = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
    
    //------------------------------------------------
    // ITS info 
    //------------------------------------------------
    fTreeV0VarPosITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( pTrack, AliPID::kPion );
    fTreeV0VarPosITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( pTrack, AliPID::kProton );
    fTreeV0VarNegITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( nTrack, AliPID::kPion   );
    fTreeV0VarNegITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( nTrack, AliPID::kProton );
    
    //------------------------------------------------
    // TOF info 
    //------------------------------------------------
    fTreeV0VarPosTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( pTrack, AliPID::kPion );
    fTreeV0VarPosTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( pTrack, AliPID::kProton );
    fTreeV0VarNegTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( nTrack, AliPID::kPion   );
    fTreeV0VarNegTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( nTrack, AliPID::kProton );
    
    //------------------------------------------------
    // Raw TPC dEdx + PIDForTracking information
    //------------------------------------------------
    
    //Acquire TPC Signals
    fTreeV0VarPosdEdx   = pTrack->GetTPCsignal();
    fTreeV0VarNegdEdx   = nTrack->GetTPCsignal();
    
    fTreeV0VarPosdEdxN  = pTrack->GetTPCsignalN();
    fTreeV0VarNegdEdxN  = nTrack->GetTPCsignalN();
    
    //Acquire PID For Tracking : uses AliPID EParticleType code (0=electron, 1=muon, 2=pion, etc)
    fTreeV0VarPosPIDForTracking = pTrack->GetPIDForTracking();
    fTreeV0VarNegPIDForTracking = nTrack->GetPIDForTracking();
    
    //________________________________________________________________________
    // Track status
    fTreeV0VarPosTrackStatus      = pTrack->GetStatus();
    fTreeV0VarNegTrackStatus      = nTrack->GetStatus();

    //________________________________________________________________________
    // Track DCAz
    fTreeV0VarPosDCAz = GetDCAz(pTrack);
    fTreeV0VarNegDCAz = GetDCAz(nTrack);
    
    //________________________________________________________________________
    // Momentum info
    fTreeV0VarPosPx  = lPosPx  ;   fTreeV0VarPosPy  = lPosPy  ;    fTreeV0VarPosPz  = lPosPz ;
    fTreeV0VarNegPx  = lNegPx  ;   fTreeV0VarNegPy  = lNegPy  ;    fTreeV0VarNegPz  = lNegPz ;
    
    //________________________________________________________________________
    // Decay vtx info
    fTreeV0VarV0DecayX = lPosV0[0] ; 
    fTreeV0VarV0DecayY = lPosV0[1] ;  
    fTreeV0VarV0DecayZ = lPosV0[2] ; 
    
    //________________________________________________________________________
    // Particle label
    fTreeV0VarPosLabel    = pTrack->GetLabel();
    fTreeV0VarNegLabel    = nTrack->GetLabel();
    
    //________________________________________________________________________
    //Check clusters
    fTreeV0VarPosITSClusters0 = pTrack->HasPointOnITSLayer(0);
    fTreeV0VarPosITSClusters1 = pTrack->HasPointOnITSLayer(1);
    fTreeV0VarPosITSClusters2 = pTrack->HasPointOnITSLayer(2);
    fTreeV0VarPosITSClusters3 = pTrack->HasPointOnITSLayer(3);
    fTreeV0VarPosITSClusters4 = pTrack->HasPointOnITSLayer(4);
    fTreeV0VarPosITSClusters5 = pTrack->HasPointOnITSLayer(5);
    
    fTreeV0VarNegITSClusters0 = nTrack->HasPointOnITSLayer(0);
    fTreeV0VarNegITSClusters1 = nTrack->HasPointOnITSLayer(1);
    fTreeV0VarNegITSClusters2 = nTrack->HasPointOnITSLayer(2);
    fTreeV0VarNegITSClusters3 = nTrack->HasPointOnITSLayer(3);
    fTreeV0VarNegITSClusters4 = nTrack->HasPointOnITSLayer(4);
    fTreeV0VarNegITSClusters5 = nTrack->HasPointOnITSLayer(5);
    
    //________________________________________________________________________
    //Check its clusters, shared
    fTreeV0VarPosITSSharedClusters0 = pTrack->HasSharedPointOnITSLayer(0);
    fTreeV0VarPosITSSharedClusters1 = pTrack->HasSharedPointOnITSLayer(1);
    fTreeV0VarPosITSSharedClusters2 = pTrack->HasSharedPointOnITSLayer(2);
    fTreeV0VarPosITSSharedClusters3 = pTrack->HasSharedPointOnITSLayer(3);
    fTreeV0VarPosITSSharedClusters4 = pTrack->HasSharedPointOnITSLayer(4);
    fTreeV0VarPosITSSharedClusters5 = pTrack->HasSharedPointOnITSLayer(5);
    
    fTreeV0VarNegITSSharedClusters0 = nTrack->HasSharedPointOnITSLayer(0);
    fTreeV0VarNegITSSharedClusters1 = nTrack->HasSharedPointOnITSLayer(1);
    fTreeV0VarNegITSSharedClusters2 = nTrack->HasSharedPointOnITSLayer(2);
    fTreeV0VarNegITSSharedClusters3 = nTrack->HasSharedPointOnITSLayer(3);
    fTreeV0VarNegITSSharedClusters4 = nTrack->HasSharedPointOnITSLayer(4);
    fTreeV0VarNegITSSharedClusters5 = nTrack->HasSharedPointOnITSLayer(5);
    
    //________________________________________________________________________
    //GetKinkIndex condition
    fTreeV0VarPosIsKink   = kFALSE;
    fTreeV0VarNegIsKink   = kFALSE;
    if( pTrack->GetProdVertex()->GetType() == AliAODVertex::kKink ) fTreeV0VarPosIsKink   = kTRUE;
    if( nTrack->GetProdVertex()->GetType() == AliAODVertex::kKink ) fTreeV0VarNegIsKink   = kTRUE;
    
    //________________________________________________________________________
    //TOF info
    fTreeV0VarPosTOFExpTDiff      = pTrack->GetTOFExpTDiff( lMagField );
    fTreeV0VarNegTOFExpTDiff      = nTrack->GetTOFExpTDiff( lMagField );
    
    fTreeV0VarPosTOFSignal        = pTrack->GetTOFsignal() * 1.e-3; // in ns
    fTreeV0VarNegTOFSignal        = nTrack->GetTOFsignal() * 1.e-3; // in ns
    
    fTreeV0VarPosTOFBCid          = pTrack->GetTOFBunchCrossing( lMagField );
    fTreeV0VarNegTOFBCid          = nTrack->GetTOFBunchCrossing( lMagField );

    fTreeV0VarPosTOFLength        = pTrack->GetIntegratedLength();
    fTreeV0VarNegTOFLength        = nTrack->GetIntegratedLength();
    
    fTreeV0VarPosTOFDeltaX        = pTrack->GetTOFsignalDx();
    fTreeV0VarNegTOFDeltaX        = nTrack->GetTOFsignalDx();
    
    fTreeV0VarPosTOFDeltaZ        = pTrack->GetTOFsignalDz();
    fTreeV0VarNegTOFDeltaZ        = nTrack->GetTOFsignalDz();
    
    //________________________________________________________________________
    //Cowboy/sailor info regarding V0 
    //Calculate vec prod with momenta projected to xy plane
    //Provisions for cowboy/sailor check
    Double_t lModp1 = TMath::Sqrt( lPosPx*lPosPx + lPosPy*lPosPy );
    Double_t lModp2 = TMath::Sqrt( lNegPx*lNegPx + lNegPy*lNegPy );
    
    //Calculate vec prod with momenta projected to xy plane
    Double_t lVecProd = (lPosPx*lNegPy - lPosPy*lNegPx) / (lModp1*lModp2);
    
    if ( lMagField < 0 ) lVecProd *= -1; //invert sign
    
    fTreeV0VarIsCowboy = kFALSE;
    if (lVecProd < 0) fTreeV0VarIsCowboy = kTRUE;
    
    //------------------------------------------------
    // Associate V0 Candidates to Monte Carlo!
    //------------------------------------------------
    
    // Abs value = needed ! question of quality track association ...
    Int_t lPosLabel     = (Int_t) TMath::Abs( pTrack->GetLabel() );
    Int_t lNegLabel     = (Int_t) TMath::Abs( nTrack->GetLabel() );
    
    AliMCParticle* pTrackMC = (AliMCParticle*)lMCEvent->GetTrack( lPosLabel  );
    AliMCParticle* nTrackMC = (AliMCParticle*)lMCEvent->GetTrack( lNegLabel  );
    
    if (!pTrackMC || !nTrackMC ) 
    {
        AliWarning("ERROR: Could not retrieve one of the V0 daughter MC tracks");
        return;
    }
    
    if( pTrackMC->IsPhysicalPrimary() ) fTreeV0VarIsPhysicalPrimaryPos  = kTRUE;
    if( nTrackMC->IsPhysicalPrimary() ) fTreeV0VarIsPhysicalPrimaryNeg  = kTRUE;
    
    //Get MC information
    //________________________________________________________________________
    // Pseudo-rapidity info
    fTreeV0VarPosEtaMC  = pTrackMC->Eta() ;
    fTreeV0VarNegEtaMC  = nTrackMC->Eta() ;
    
    //________________________________________________________________________
    // Momentum info
    fTreeV0VarPosPxMC  = pTrackMC->Px() ;   fTreeV0VarPosPyMC  = pTrackMC->Py() ;    fTreeV0VarPosPzMC  = pTrackMC->Pz() ;
    fTreeV0VarNegPxMC  = nTrackMC->Px() ;   fTreeV0VarNegPyMC  = nTrackMC->Py() ;    fTreeV0VarNegPzMC  = nTrackMC->Pz() ;
    
    //________________________________________________________________________
    // PID & label info
    fTreeV0VarPosPdgCode  = pTrackMC->PdgCode();
    fTreeV0VarNegPdgCode  = nTrackMC->PdgCode();
    
    Int_t lPosMotherLabel   = pTrackMC->GetMother();
    Int_t lNegMotherLabel   = nTrackMC->GetMother();
    
    // Extra: check mother particle pdg code completely independently
    if ( lPosMotherLabel >= 0 )
    {
        AliMCParticle* lPosMother = (AliMCParticle*) lMCEvent->GetTrack( lPosMotherLabel );
        if( lPosMother->IsPhysicalPrimary() ) fTreeV0VarIsPhysicalPrimaryPosMother = kTRUE;
        // Mother PDG 
        fTreeV0VarPosMotherPdgCode  = lPosMother->PdgCode();
        fTreeV0VarPosMotherLabel    = lPosMotherLabel;
        
        //Go further than that, please
        Int_t lGrandMotherLabel = lPosMother->GetMother();
        if( lGrandMotherLabel >= 0 )
        {
            AliMCParticle* lPosGrandMother = (AliMCParticle*) lMCEvent->GetTrack( lGrandMotherLabel );
            if( lPosGrandMother->IsPhysicalPrimary() ) fTreeV0VarIsPhysicalPrimaryPosGrandMother = kTRUE;
            // Grand Mother PDG
            fTreeV0VarPosGrandMotherPdgCode   = lPosGrandMother->PdgCode();
            fTreeV0VarPosGrandMotherLabel     = lGrandMotherLabel;
        }
    }
    
    if ( lNegMotherLabel >= 0 )
    {
        AliMCParticle* lNegMother = (AliMCParticle*) lMCEvent->GetTrack( lNegMotherLabel );
        if( lNegMother->IsPhysicalPrimary() ) fTreeV0VarIsPhysicalPrimaryNegMother = kTRUE;
        // Mother PID 
        fTreeV0VarNegMotherPdgCode  = lNegMother->PdgCode();
        fTreeV0VarNegMotherLabel    = lNegMotherLabel;
        
        //Go further than that, please
        Int_t lGrandMotherLabel = lNegMother->GetMother();
        if( lGrandMotherLabel >= 0 )
        {
            AliMCParticle* lNegGrandMother = (AliMCParticle*) lMCEvent->GetTrack( lGrandMotherLabel );
            if( lNegGrandMother->IsPhysicalPrimary() ) fTreeV0VarIsPhysicalPrimaryNegGrandMother = kTRUE;
            // Grand Mother PDG
            fTreeV0VarNegGrandMotherPdgCode = lNegGrandMother->PdgCode();
            fTreeV0VarNegGrandMotherLabel   = lGrandMotherLabel;
        }
    }
    
    //same mother && mother != primary (!= -1)
    if( lPosMotherLabel == lNegMotherLabel && lPosMotherLabel >= 0 )
    { 
        //either label is fine, they're equal at this stage
        AliMCParticle* pV0MC = (AliMCParticle*)lMCEvent->GetTrack( lPosMotherLabel );
        
        fTreeV0VarPdgCode       = pV0MC->PdgCode();  //PDG Code
        
        //Only Interested if it's a Lambda, AntiLambda or K0s
        //Avoid the Junction Bug! PYTHIA has particles with Px=Py=Pz=E=0 occasionally,
        //having particle code 88 (unrecognized by PDG), for documentation purposes.
        //Even ROOT's TParticle::Y() is not prepared to deal with that exception!
        //Note that TParticle::Pt() is immune (that would just return 0)...
        //Though granted that that should be extremely rare in this precise condition...
        if( TMath::Abs(fTreeV0VarPdgCode) == 3122 || fTreeV0VarPdgCode == 310 ) 
            fTreeV0VarRapMC = pV0MC->Y(); //Perfect Y
        
        fTreeV0VarPtotMC     = pV0MC->P();                          //Perfect Ptot
        fTreeV0VarPtMC       = pV0MC->Pt();                         //Perfect Pt
        
        if( pV0MC->IsPhysicalPrimary       () ) fTreeV0VarIsPhysicalPrimary = 1; //Is Primary!
        if( pV0MC->IsSecondaryFromWeakDecay() ) fTreeV0VarIsPhysicalPrimary = 2; //Weak Decay!
        if( pV0MC->IsSecondaryFromMaterial () ) fTreeV0VarIsPhysicalPrimary = 3; //Material Int!
        
        
        Int_t lV0MotherLabel = pV0MC->GetMother() ;
        //if it has a parent, get it and store specs
        if( lV0MotherLabel >= 0 ) 
        { 
            AliMCParticle* pV0MotherMC = (AliMCParticle*)lMCEvent->GetTrack( lV0MotherLabel );
            
            fTreeV0VarMotherPdgCode = pV0MotherMC->PdgCode(); //V0 Mother PDG
            fTreeV0VarMotherPt      = pV0MotherMC->Pt();         //V0 Mother Pt
            //NOTE: Fill only for charged Xi or Omega
            if( TMath::Abs(fTreeV0VarMotherPdgCode) == 3312 || TMath::Abs(fTreeV0VarMotherPdgCode) == 3334) 
                fTreeV0VarMotherRapMC   = pV0MotherMC->Y();         //V0 Mother rapidity
                
            //Primary Status for the V0 Mother particle
            if( pV0MotherMC->IsPhysicalPrimary       () ) fTreeV0VarMotherIsPhysicalPrimary = 1; //Is Primary!
            if( pV0MotherMC->IsSecondaryFromWeakDecay() ) fTreeV0VarMotherIsPhysicalPrimary = 2; //Weak Decay!
            if( pV0MotherMC->IsSecondaryFromMaterial () ) fTreeV0VarMotherIsPhysicalPrimary = 3; //Material Int!
        }
        
    }
    
    //________________________________________________________________________
    // Decay position info
    //Be careful: Vx, Vy, Vz: Creation vertex. So decay position is the
    //Creation vertex of any one of the daughters!
    fTreeV0VarDecayXMC = pTrackMC->Xv();
    fTreeV0VarDecayYMC = pTrackMC->Yv();
    fTreeV0VarDecayZMC = pTrackMC->Zv();
}

//_____________________________________________________________________________
Float_t AliMCV0Container::GetLengthInActiveZone( AliAODTrack *gt, Float_t deltaY, Float_t deltaZ, Float_t b ){
    AliESDtrack esdTrack( gt );
    esdTrack.SetESDEvent((AliESDEvent*) gt->GetEvent() );
    AliExternalTrackParam etp;
    etp.CopyFromVTrack(gt);
    esdTrack.ResetTrackParamIp(&etp);
    return esdTrack.GetLengthInActiveZone(1, deltaY, deltaZ, b);
}
//_____________________________________________________________________________
Float_t AliMCV0Container::GetDCAz(AliAODTrack *lTrack)
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
Bool_t AliMCV0Container::IsSelected( const AliV0Result* lV0Result )
{
    Double_t lK0sPDGMass    = 0.493677;
    Double_t lLambdaPDGMass = 1.115683;
    
    Bool_t      lIsK0s;
    Double_t    lInvMass;
    Double_t    lPDGMass;
    Double_t    lV0Mass;
    Double_t    lRap;
    Double_t    lPosTPCNSigma;
    Double_t    lNegTPCNSigma;
    
    //Second Selection: rough 20-sigma band, parametric.
    //K0Short: Enough to parametrize peak broadening with linear function.
    Double_t lUpperLimitK0Short = (5.63707e-01) + (1.14979e-02)*fTreeV0VarPt;
    Double_t lLowerLimitK0Short = (4.30006e-01) - (1.10029e-02)*fTreeV0VarPt;
    //Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
    //[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
    Double_t lUpperLimitLambda = (1.13688e+00) + (5.27838e-03)*fTreeV0VarPt + (8.42220e-02)*TMath::Exp(-(3.80595e+00)*fTreeV0VarPt);
    Double_t lLowerLimitLambda = (1.09501e+00) - (5.23272e-03)*fTreeV0VarPt - (7.52690e-02)*TMath::Exp(-(3.46339e+00)*fTreeV0VarPt);
    
    
    //========================================================================
    if ( lV0Result->GetMassHypothesis() == AliV0Result::kK0Short        )
    {
        lIsK0s                  = kTRUE;
        lInvMass                = fTreeV0VarMassAsK0s;
        lPDGMass                = lK0sPDGMass;
        lRap                    = fTreeV0VarRapK0Short;
        
        lPosTPCNSigma           = fTreeV0VarPosNSigmaPion;
        lNegTPCNSigma           = fTreeV0VarNegNSigmaPion;
    }
    if ( lV0Result->GetMassHypothesis() == AliV0Result::kLambda         )
    {
        lIsK0s                  = kFALSE;
        lInvMass                = fTreeV0VarMassAsLambda;
        lPDGMass                = lLambdaPDGMass;
        lRap                    = fTreeV0VarRapLambda;
        
        lPosTPCNSigma           = fTreeV0VarPosNSigmaProton;
        lNegTPCNSigma           = fTreeV0VarNegNSigmaPion;
    }
    if ( lV0Result->GetMassHypothesis() == AliV0Result::kAntiLambda     )
    {
        lIsK0s                  = kFALSE;
        lInvMass                = fTreeV0VarMassAsAntiLambda;
        lPDGMass                = lLambdaPDGMass;
        lRap                    = fTreeV0VarRapLambda;
        
        lPosTPCNSigma           = fTreeV0VarPosNSigmaPion;
        lNegTPCNSigma           = fTreeV0VarNegNSigmaProton;
    }
    
    if (
        //Check 1: Offline Vertexer
        fTreeV0VarOnFlyStatus == lV0Result->GetUseOnTheFly() &&
        
        //Check 2: Basic Acceptance cuts
        // Positive
        fTreeV0VarPosEta                                > lV0Result->GetCutMinEtaTracks()                           && 
        fTreeV0VarPosEta                                < lV0Result->GetCutMaxEtaTracks()                           &&
        // Negative
        fTreeV0VarNegEta                                > lV0Result->GetCutMinEtaTracks()                           &&
        fTreeV0VarNegEta                                < lV0Result->GetCutMaxEtaTracks()                           &&
        
        // Corresponding rapidity
        lRap                                            > lV0Result->GetCutMinRapidity()                            &&
        lRap                                            < lV0Result->GetCutMaxRapidity()                            &&
        
        //Check 3: Topological Variables
        fTreeV0VarDCAPosToPV                            > lV0Result->GetCutDCAPosToPV()                             &&
        fTreeV0VarDCANegToPV                            > lV0Result->GetCutDCANegToPV()                             &&
        fTreeV0VarDCAV0Dau                              < lV0Result->GetCutDCAV0Daughters()                         &&
        fTreeV0VarRadius                                > lV0Result->GetCutV0Radius()                               &&
        fTreeV0VarCosPA                                 > lV0Result->GetCutV0CosPA()                                &&
        fTreeV0VarLeastNbrCrossedRows                   > lV0Result->GetCutLeastNumberOfCrossedRows()               &&
        lPDGMass*fTreeV0VarDistOverTotMom               < lV0Result->GetCutProperLifetime()                         &&
        // - Miscellaneous
        fTreeV0VarLeastRatioCrossedRowsOverFindable     > lV0Result->GetCutLeastNumberOfCrossedRowsOverFindable()   &&
        
        //Check 4: TPC dEdx selections
        TMath::Abs(lPosTPCNSigma )                      < lV0Result->GetCutTPCdEdx()                                &&
        TMath::Abs(lNegTPCNSigma )                      < lV0Result->GetCutTPCdEdx()                                &&
        
        //Check 5: Min Track Length 
        fTreeV0VarMinTrackLength                        > lV0Result->GetCutMinTrackLength()                         &&

        //Check 6: Max Chi2/Clusters 
        fTreeV0VarMaxChi2PerCluster                     < lV0Result->GetCutMaxChi2PerCluster()                      &&
        
        //Check 7a: K0s/Lambda mass rejection
        (   //Is K0s --> cut on the Lambda mass
            (
                lIsK0s && 
                TMath::Abs( fTreeV0VarMassAsLambda - lLambdaPDGMass ) > lV0Result->GetCutCompetingV0Rejection() 
            ) 
            ||
            //Is Lambda --> cut on the K0s mass
            (
                !lIsK0s && 
                TMath::Abs( fTreeV0VarMassAsK0s - lK0sPDGMass       ) > lV0Result->GetCutCompetingV0Rejection() 
            )
        )                                                                                                          &&
        
        //Check 7b: K0s/Lambda mass cut
        (   //Is K0s
            (
                lIsK0s && 
                lLowerLimitK0Short  < lInvMass  &&  lInvMass    < lUpperLimitK0Short
            ) 
            ||
            //Is Lambda
            (
                !lIsK0s && 
                lLowerLimitLambda   < lInvMass  &&  lInvMass    < lUpperLimitLambda
            )
        )                                                                                                          &&
        
        //Check 8: kITSrefit track selection if requested
        (   //do not need ITS refit tracks --> skip this cut
            !( lV0Result->GetCutUseITSRefitTracks() )                           ||
            //otherwise : need to have ITS refit tracks
            (   ( fTreeV0VarPosTrackStatus    & AliAODTrack::kITSrefit )        &&
                ( fTreeV0VarNegTrackStatus    & AliAODTrack::kITSrefit )  
            )
        )                                                                                                           &&
        
        //Check 9: cowboy/sailor for V0
        ( lV0Result->GetCutIsCowboy()  ==  0                                    ||
        (lV0Result->GetCutIsCowboy()   ==  1 && fTreeV0VarIsCowboy           )  ||
        (lV0Result->GetCutIsCowboy()   == -1 && !fTreeV0VarIsCowboy          ) )                                    &&
        
        //Check 10: ITS or TOF required 
        (
            !( lV0Result->GetCutITSorTOF() )                                      || 
            (   //One of the daughter track has ITSrefit 
                (   (fTreeV0VarPosTrackStatus   & AliAODTrack::kITSrefit)       ||
                    (fTreeV0VarNegTrackStatus   & AliAODTrack::kITSrefit) 
                )||
                //Or one of the daughter track has a TOF signal
                (   (TMath::Abs(fTreeV0VarPosTOFExpTDiff+2500.  )     > 1e-6 )  ||
                    (TMath::Abs(fTreeV0VarNegTOFExpTDiff+2500.  )     > 1e-6 ) 
                )
            )
            
        )
        
    )//end major if
    {
        return kTRUE;
    }
    
    return kFALSE;
}
