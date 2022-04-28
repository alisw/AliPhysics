#include "AliMCCascadeContainer.h"

class AliMCCascadeContainer;    // your analysis class

ClassImp(AliMCCascadeContainer) // classimp: necessary for root

//_____________________________________________________________________________
AliMCCascadeContainer::AliMCCascadeContainer() :
TObject                     (), 
//-----------BASIC-INFO---------------------------
fTreeCascVarCharge(0),
fTreeCascVarMassAsXi(0),
fTreeCascVarMassAsOmega(0),
fTreeCascVarPtot(0),
fTreeCascVarV0Ptot(0),
fTreeCascVarPt(0),
fTreeCascVarV0Pt(0),
fTreeCascVarRapXi(0),
fTreeCascVarRapOmega(0),
fTreeCascVarCascEta(0),
fTreeCascVarV0Eta(0),
fTreeCascVarBachEta(0),
fTreeCascVarPosEta(0),
fTreeCascVarNegEta(0),
//-----------INFO-FOR-CUTS------------------------
fTreeCascVarOnFlyStatus(kFALSE),
fTreeCascVarAlpha(0),
fTreeCascVarPtArm(0),
fTreeCascVarAlphaV0(0),
fTreeCascVarPtArmV0(0),
fTreeCascVarDCACascDau(0),
fTreeCascVarDCABachToPV(0),
fTreeCascVarDCAV0Dau(0),
fTreeCascVarDCAV0ToPV(0),
fTreeCascVarDCAPosToPV(0),
fTreeCascVarDCANegToPV(0),
fTreeCascVarCascCosPA(0),
fTreeCascVarCascCosPASpecial(0),
    
fTreeCascVarCascRadius(0),
fTreeCascVarV0Mass(0),
fTreeCascVarV0MassAsLambda(0),
fTreeCascVarV0MassAsAntiLambda(0),
fTreeCascVarV0Radius(0),
fTreeCascVarV0CosPA(0),
fTreeCascVarWrongCosPA(0),
fTreeCascVarDCABachToBaryon(0),
fTreeCascVarLeastNbrCrossedRows(0),
fTreeCascVarLeastNbrClusters(0),
fTreeCascVarLeastRatioCrossedRowsOverFindable(0),
fTreeCascVarNbrCrossedRowsOverLength(0),
fTreeCascVarMaxChi2PerCluster(0),
fTreeCascVarMinTrackLength(0),
//-----------DECAY-LENGTH-INFO--------------------
fTreeCascVarDistOverTotMom(0),
fTreeCascVarV0DistOverTotMom(0),
//------------------------------------------------
//---------------PID-TPC-INFO---------------------
fTreeCascVarBachNSigmaPion(0),
fTreeCascVarBachNSigmaKaon(0),
fTreeCascVarPosNSigmaProton(0),
fTreeCascVarPosNSigmaPion(0),
fTreeCascVarNegNSigmaProton(0),
fTreeCascVarNegNSigmaPion(0),
//---------------PID-TOF-INFO---------------------
fTreeCascVarBachTOFNSigmaPion(0),
fTreeCascVarBachTOFNSigmaKaon(0),
fTreeCascVarPosTOFNSigmaProton(0),
fTreeCascVarPosTOFNSigmaPion(0),
fTreeCascVarNegTOFNSigmaProton(0),
fTreeCascVarNegTOFNSigmaPion(0),
//---------------PID-ITS-INFO---------------------
fTreeCascVarBachITSNSigmaPion(0),
fTreeCascVarBachITSNSigmaKaon(0),
fTreeCascVarPosITSNSigmaProton(0),
fTreeCascVarPosITSNSigmaPion(0),
fTreeCascVarNegITSNSigmaProton(0),
fTreeCascVarNegITSNSigmaPion(0),
//---Raw TPC dEdx + PIDForTracking information----
fTreeCascVarPosdEdx(0),
fTreeCascVarNegdEdx(0),
fTreeCascVarBachdEdx(0),
fTreeCascVarPosdEdxN(0),
fTreeCascVarNegdEdxN(0),
fTreeCascVarBachdEdxN(0),
fTreeCascVarPosPIDForTracking(0),
fTreeCascVarNegPIDForTracking(0),
fTreeCascVarBachPIDForTracking(0),
//------------------------------------------------
fTreeCascVarChi2Cascade(0),
fTreeCascVarChi2CascadePerNDF(0),
fTreeCascVarChi2V0(0),
//------------------------------------------------
fTreeCascVarBachTrackStatus(0),
fTreeCascVarPosTrackStatus(0),
fTreeCascVarNegTrackStatus(0),
//------------------------------------------------
fTreeCascVarBachDCAz(0),
fTreeCascVarPosDCAz(0),
fTreeCascVarNegDCAz(0),
//------------FULL-MOMENTUM-INFO------------------
fTreeCascVarBachPx(0),
fTreeCascVarBachPy(0),
fTreeCascVarBachPz(0),
fTreeCascVarPosPx(0),
fTreeCascVarPosPy(0),
fTreeCascVarPosPz(0),
fTreeCascVarNegPx(0),
fTreeCascVarNegPy(0),
fTreeCascVarNegPz(0),
//------------------------------------------------
fTreeCascVarCascadeDecayX(0),
fTreeCascVarCascadeDecayY(0),
fTreeCascVarCascadeDecayZ(0),
fTreeCascVarV0DecayX(0),
fTreeCascVarV0DecayY(0),
fTreeCascVarV0DecayZ(0),
//------------------------------------------------
fTreeCascVarBachLabel(0),
fTreeCascVarPosLabel(0),
fTreeCascVarNegLabel(0),
//---------------CLUSTER-INFO---------------------
fTreeCascVarBachITSClusters0(kFALSE),
fTreeCascVarBachITSClusters1(kFALSE),
fTreeCascVarBachITSClusters2(kFALSE),
fTreeCascVarBachITSClusters3(kFALSE),
fTreeCascVarBachITSClusters4(kFALSE),
fTreeCascVarBachITSClusters5(kFALSE),
    
fTreeCascVarPosITSClusters0(kFALSE),
fTreeCascVarPosITSClusters1(kFALSE),
fTreeCascVarPosITSClusters2(kFALSE),
fTreeCascVarPosITSClusters3(kFALSE),
fTreeCascVarPosITSClusters4(kFALSE),
fTreeCascVarPosITSClusters5(kFALSE),
    
fTreeCascVarNegITSClusters0(kFALSE),
fTreeCascVarNegITSClusters1(kFALSE),
fTreeCascVarNegITSClusters2(kFALSE),
fTreeCascVarNegITSClusters3(kFALSE),
fTreeCascVarNegITSClusters4(kFALSE),
fTreeCascVarNegITSClusters5(kFALSE),
    
//------------------------------------------------
fTreeCascVarPosITSSharedClusters0(kFALSE),
fTreeCascVarPosITSSharedClusters1(kFALSE),
fTreeCascVarPosITSSharedClusters2(kFALSE),
fTreeCascVarPosITSSharedClusters3(kFALSE),
fTreeCascVarPosITSSharedClusters4(kFALSE),
fTreeCascVarPosITSSharedClusters5(kFALSE),

fTreeCascVarNegITSSharedClusters0(kFALSE),
fTreeCascVarNegITSSharedClusters1(kFALSE),
fTreeCascVarNegITSSharedClusters2(kFALSE),
fTreeCascVarNegITSSharedClusters3(kFALSE),
fTreeCascVarNegITSSharedClusters4(kFALSE),
fTreeCascVarNegITSSharedClusters5(kFALSE),

fTreeCascVarBachITSSharedClusters0(kFALSE),
fTreeCascVarBachITSSharedClusters1(kFALSE),
fTreeCascVarBachITSSharedClusters2(kFALSE),
fTreeCascVarBachITSSharedClusters3(kFALSE),
fTreeCascVarBachITSSharedClusters4(kFALSE),
fTreeCascVarBachITSSharedClusters5(kFALSE),

//---------------OOB-PILEUP-INFO---------------------
fTreeCascVarBachTOFExpTDiff(0),
fTreeCascVarPosTOFExpTDiff(0),
fTreeCascVarNegTOFExpTDiff(0),

fTreeCascVarBachTOFSignal(0),
fTreeCascVarPosTOFSignal(0),
fTreeCascVarNegTOFSignal(0),
    
fTreeCascVarBachTOFBCid(0),
fTreeCascVarPosTOFBCid(0),
fTreeCascVarNegTOFBCid(0),

fTreeCascVarBachTOFLength(0),
fTreeCascVarPosTOFLength(0),
fTreeCascVarNegTOFLength(0),

fTreeCascVarBachTOFDeltaX(0),
fTreeCascVarPosTOFDeltaX(0),
fTreeCascVarNegTOFDeltaX(0),

fTreeCascVarBachTOFDeltaZ(0),
fTreeCascVarPosTOFDeltaZ(0),
fTreeCascVarNegTOFDeltaZ(0),
    
//Kink tagging
fTreeCascVarBachIsKink(kFALSE),
fTreeCascVarPosIsKink(kFALSE),
fTreeCascVarNegIsKink(kFALSE),
    
//Cowboy/sailor studies
fTreeCascVarIsCowboy(kFALSE),   //store if V0 is cowboy-like or sailor-like in XY plane
fTreeCascVarCowboyness(0), //negative -> cowboy, positive -> sailor
fTreeCascVarIsCascadeCowboy(kFALSE),   //store if V0 is cowboy-like or sailor-like in XY plane
fTreeCascVarCascadeCowboyness(0),//negative -> cowboy, positive -> sailor
//------------MONTE-CARLO-INFO------------------
fTreeCascVarBachEtaMC(-100),
fTreeCascVarPosEtaMC(-100),
fTreeCascVarNegEtaMC(-100),

fTreeCascVarBachPxMC(-100),
fTreeCascVarBachPyMC(-100),
fTreeCascVarBachPzMC(-100),
fTreeCascVarPosPxMC(-100),
fTreeCascVarPosPyMC(-100),
fTreeCascVarPosPzMC(-100),
fTreeCascVarNegPxMC(-100),
fTreeCascVarNegPyMC(-100),
fTreeCascVarNegPzMC(-100),

fTreeCascVarBachMotherLabel(-1),
fTreeCascVarPosMotherLabel(-1),
fTreeCascVarNegMotherLabel(-1),

fTreeCascVarBachGrandMotherLabel(-1),
fTreeCascVarPosGrandMotherLabel(-1),
fTreeCascVarNegGrandMotherLabel(-1),

fTreeCascVarBachPdgCode(-9999),
fTreeCascVarPosPdgCode(-9999),
fTreeCascVarNegPdgCode(-9999),

fTreeCascVarBachMotherPdgCode(-9999),
fTreeCascVarPosMotherPdgCode(-9999),
fTreeCascVarNegMotherPdgCode(-9999),

fTreeCascVarBachGrandMotherPdgCode(-9999),
fTreeCascVarPosGrandMotherPdgCode(-9999),
fTreeCascVarNegGrandMotherPdgCode(-9999),

fTreeCascVarPdgCode(-9999),
fTreeCascVarRapMC(-100),
fTreeCascVarPtotMC(-100),
fTreeCascVarPtMC(-100),

fTreeCascVarV0PdgCode(-9999),
fTreeCascVarV0PtotMC(-100),
fTreeCascVarV0PtMC(-100),

fTreeCascVarCascadeDecayXMC(-100),
fTreeCascVarCascadeDecayYMC(-100),
fTreeCascVarCascadeDecayZMC(-100),

fTreeCascVarV0DecayXMC(-100),
fTreeCascVarV0DecayYMC(-100),
fTreeCascVarV0DecayZMC(-100),

fTreeCascVarIsPhysicalPrimary(0),

fTreeCascVarIsPhysicalPrimaryBach(kFALSE),
fTreeCascVarIsPhysicalPrimaryPos(kFALSE),
fTreeCascVarIsPhysicalPrimaryNeg(kFALSE),

fTreeCascVarIsPhysicalPrimaryBachMother(kFALSE),
fTreeCascVarIsPhysicalPrimaryPosMother(kFALSE),
fTreeCascVarIsPhysicalPrimaryNegMother(kFALSE),

fTreeCascVarIsPhysicalPrimaryBachGrandMother(kFALSE),
fTreeCascVarIsPhysicalPrimaryPosGrandMother(kFALSE),
fTreeCascVarIsPhysicalPrimaryNegGrandMother(kFALSE)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliMCCascadeContainer::~AliMCCascadeContainer() // destructor
{
    // destructor
}
//_____________________________________________________________________________
void AliMCCascadeContainer::Reset()
{
    //-----------BASIC-INFO---------------------------
    fTreeCascVarCharge = 0;
    fTreeCascVarMassAsXi = 0;
    fTreeCascVarMassAsOmega = 0;
    fTreeCascVarPtot = 0;
    fTreeCascVarV0Ptot = 0;
    fTreeCascVarPt = 0;
    fTreeCascVarV0Pt = 0;
    fTreeCascVarRapXi = 0;
    fTreeCascVarRapOmega = 0;
    fTreeCascVarCascEta = 0;
    fTreeCascVarV0Eta   = 0;
    fTreeCascVarBachEta = 0;
    fTreeCascVarPosEta = 0;
    fTreeCascVarNegEta = 0;
    //-----------INFO-FOR-CUTS------------------------
    fTreeCascVarOnFlyStatus = kFALSE;
    fTreeCascVarAlpha = 0;
    fTreeCascVarPtArm = 0;
    fTreeCascVarAlphaV0 = 0;
    fTreeCascVarPtArmV0 = 0;
    fTreeCascVarDCACascDau = 0;
    fTreeCascVarDCABachToPV = 0;
    fTreeCascVarDCAV0Dau = 0;
    fTreeCascVarDCAV0ToPV = 0;
    fTreeCascVarDCAPosToPV = 0;
    fTreeCascVarDCANegToPV = 0;
    fTreeCascVarCascCosPA = 0;
    fTreeCascVarCascCosPASpecial = 0;
        
    fTreeCascVarCascRadius = 0;
    fTreeCascVarV0Mass = 0;
    fTreeCascVarV0MassAsLambda = 0;
    fTreeCascVarV0MassAsAntiLambda = 0;
    fTreeCascVarV0Radius = 0;
    fTreeCascVarV0CosPA = 0;
    fTreeCascVarWrongCosPA = 0;
    fTreeCascVarDCABachToBaryon = 0;
    fTreeCascVarLeastNbrCrossedRows = 0;
    fTreeCascVarLeastNbrClusters = 0;
    fTreeCascVarLeastRatioCrossedRowsOverFindable = 0;
    fTreeCascVarNbrCrossedRowsOverLength = 0;
    fTreeCascVarMaxChi2PerCluster = 0;
    fTreeCascVarMinTrackLength = 0;
    //-----------DECAY-LENGTH-INFO--------------------
    fTreeCascVarDistOverTotMom = 0;
    fTreeCascVarV0DistOverTotMom = 0;
    //------------------------------------------------
    //---------------PID-TPC-INFO---------------------
    fTreeCascVarBachNSigmaPion = 0;
    fTreeCascVarBachNSigmaKaon = 0;
    fTreeCascVarPosNSigmaProton = 0;
    fTreeCascVarPosNSigmaPion = 0;
    fTreeCascVarNegNSigmaProton = 0;
    fTreeCascVarNegNSigmaPion = 0;
    //---------------PID-TOF-INFO---------------------
    fTreeCascVarBachTOFNSigmaPion = 0;
    fTreeCascVarBachTOFNSigmaKaon = 0;
    fTreeCascVarPosTOFNSigmaProton = 0;
    fTreeCascVarPosTOFNSigmaPion = 0;
    fTreeCascVarNegTOFNSigmaProton = 0;
    fTreeCascVarNegTOFNSigmaPion = 0;
    //---------------PID-ITS-INFO---------------------
    fTreeCascVarBachITSNSigmaPion = 0;
    fTreeCascVarBachITSNSigmaKaon = 0;
    fTreeCascVarPosITSNSigmaProton = 0;
    fTreeCascVarPosITSNSigmaPion = 0;
    fTreeCascVarNegITSNSigmaProton = 0;
    fTreeCascVarNegITSNSigmaPion = 0;
    //---Raw TPC dEdx + PIDForTracking information----
    fTreeCascVarPosdEdx             = 0;
    fTreeCascVarNegdEdx             = 0;
    fTreeCascVarBachdEdx            = 0;
    fTreeCascVarPosdEdxN            = 0;
    fTreeCascVarNegdEdxN            = 0;
    fTreeCascVarBachdEdxN           = 0;
    fTreeCascVarPosPIDForTracking   = 0;
    fTreeCascVarNegPIDForTracking   = 0;
    fTreeCascVarBachPIDForTracking  = 0;
    //------------------------------------------------
    fTreeCascVarChi2Cascade = 0;
    fTreeCascVarChi2CascadePerNDF = 0;
    fTreeCascVarChi2V0 = 0;
    //------------------------------------------------
    fTreeCascVarBachTrackStatus = 0;
    fTreeCascVarPosTrackStatus = 0;
    fTreeCascVarNegTrackStatus = 0;
    //------------------------------------------------
    fTreeCascVarBachDCAz = 0;
    fTreeCascVarPosDCAz = 0;
    fTreeCascVarNegDCAz = 0;
    //------------FULL-MOMENTUM-INFO------------------
    fTreeCascVarBachPx = 0;
    fTreeCascVarBachPy = 0;
    fTreeCascVarBachPz = 0;
    fTreeCascVarPosPx = 0;
    fTreeCascVarPosPy = 0;
    fTreeCascVarPosPz = 0;
    fTreeCascVarNegPx = 0;
    fTreeCascVarNegPy = 0;
    fTreeCascVarNegPz = 0;
    //------------------------------------------------
    fTreeCascVarCascadeDecayX = 0;
    fTreeCascVarCascadeDecayY = 0;
    fTreeCascVarCascadeDecayZ = 0;
    fTreeCascVarV0DecayX = 0;
    fTreeCascVarV0DecayY = 0;
    fTreeCascVarV0DecayZ = 0;
    //------------------------------------------------
    fTreeCascVarBachLabel = 0;
    fTreeCascVarPosLabel = 0;
    fTreeCascVarNegLabel = 0;
    //---------------CLUSTER-INFO---------------------
    fTreeCascVarBachITSClusters0 = kFALSE;
    fTreeCascVarBachITSClusters1 = kFALSE;
    fTreeCascVarBachITSClusters2 = kFALSE;
    fTreeCascVarBachITSClusters3 = kFALSE;
    fTreeCascVarBachITSClusters4 = kFALSE;
    fTreeCascVarBachITSClusters5 = kFALSE;
        
    fTreeCascVarPosITSClusters0 = kFALSE;
    fTreeCascVarPosITSClusters1 = kFALSE;
    fTreeCascVarPosITSClusters2 = kFALSE;
    fTreeCascVarPosITSClusters3 = kFALSE;
    fTreeCascVarPosITSClusters4 = kFALSE;
    fTreeCascVarPosITSClusters5 = kFALSE;
        
    fTreeCascVarNegITSClusters0 = kFALSE;
    fTreeCascVarNegITSClusters1 = kFALSE;
    fTreeCascVarNegITSClusters2 = kFALSE;
    fTreeCascVarNegITSClusters3 = kFALSE;
    fTreeCascVarNegITSClusters4 = kFALSE;
    fTreeCascVarNegITSClusters5 = kFALSE;
        
    //------------------------------------------------
    fTreeCascVarPosITSSharedClusters0 = kFALSE;
    fTreeCascVarPosITSSharedClusters1 = kFALSE;
    fTreeCascVarPosITSSharedClusters2 = kFALSE;
    fTreeCascVarPosITSSharedClusters3 = kFALSE;
    fTreeCascVarPosITSSharedClusters4 = kFALSE;
    fTreeCascVarPosITSSharedClusters5 = kFALSE;

    fTreeCascVarNegITSSharedClusters0 = kFALSE;
    fTreeCascVarNegITSSharedClusters1 = kFALSE;
    fTreeCascVarNegITSSharedClusters2 = kFALSE;
    fTreeCascVarNegITSSharedClusters3 = kFALSE;
    fTreeCascVarNegITSSharedClusters4 = kFALSE;
    fTreeCascVarNegITSSharedClusters5 = kFALSE;

    fTreeCascVarBachITSSharedClusters0 = kFALSE;
    fTreeCascVarBachITSSharedClusters1 = kFALSE;
    fTreeCascVarBachITSSharedClusters2 = kFALSE;
    fTreeCascVarBachITSSharedClusters3 = kFALSE;
    fTreeCascVarBachITSSharedClusters4 = kFALSE;
    fTreeCascVarBachITSSharedClusters5 = kFALSE;

    //---------------OOB-PILEUP-INFO---------------------
    fTreeCascVarBachTOFExpTDiff = 0;
    fTreeCascVarPosTOFExpTDiff = 0;
    fTreeCascVarNegTOFExpTDiff = 0;

    fTreeCascVarBachTOFSignal = 0;
    fTreeCascVarPosTOFSignal = 0;
    fTreeCascVarNegTOFSignal = 0;
        
    fTreeCascVarBachTOFBCid = 0;
    fTreeCascVarPosTOFBCid = 0;
    fTreeCascVarNegTOFBCid = 0;
    
    fTreeCascVarBachTOFLength = 0;
    fTreeCascVarPosTOFLength = 0;
    fTreeCascVarNegTOFLength = 0;

    fTreeCascVarBachTOFDeltaX = 0;
    fTreeCascVarPosTOFDeltaX = 0;
    fTreeCascVarNegTOFDeltaX = 0;

    fTreeCascVarBachTOFDeltaZ = 0;
    fTreeCascVarPosTOFDeltaZ = 0;
    fTreeCascVarNegTOFDeltaZ = 0;

    //Kink tagging
    fTreeCascVarBachIsKink = kFALSE;
    fTreeCascVarPosIsKink = kFALSE;
    fTreeCascVarNegIsKink = kFALSE;
        
    //Cowboy/sailor studies
    fTreeCascVarIsCowboy = kFALSE;   //store if V0 is cowboy-like or sailor-like in XY plane
    fTreeCascVarCowboyness = 0; //negative -> cowboy, positive -> sailor
    fTreeCascVarIsCascadeCowboy = kFALSE;   //store if V0 is cowboy-like or sailor-like in XY plane
    fTreeCascVarCascadeCowboyness = 0; //negative -> cowboy, positive -> sailor
    
    //________________________________________________________________________
    //Initializations
    fTreeCascVarBachEtaMC = -100 ;
    fTreeCascVarPosEtaMC  = -100 ;
    fTreeCascVarNegEtaMC  = -100 ;
    
    fTreeCascVarBachPxMC = -100 ;   fTreeCascVarBachPyMC = -100 ;    fTreeCascVarBachPzMC = -100 ;
    fTreeCascVarPosPxMC  = -100 ;   fTreeCascVarPosPyMC  = -100 ;    fTreeCascVarPosPzMC  = -100 ;
    fTreeCascVarNegPxMC  = -100 ;   fTreeCascVarNegPyMC  = -100 ;    fTreeCascVarNegPzMC  = -100 ;
    
    fTreeCascVarBachMotherLabel = -1;
    fTreeCascVarPosMotherLabel  = -1;
    fTreeCascVarNegMotherLabel  = -1;
    
    fTreeCascVarBachGrandMotherLabel = -1;
    fTreeCascVarPosGrandMotherLabel  = -1;
    fTreeCascVarNegGrandMotherLabel  = -1;
    
    fTreeCascVarBachPdgCode = -9999;
    fTreeCascVarPosPdgCode  = -9999;
    fTreeCascVarNegPdgCode  = -9999;
    
    fTreeCascVarBachMotherPdgCode  = -9999;
    fTreeCascVarPosMotherPdgCode   = -9999;
    fTreeCascVarNegMotherPdgCode   = -9999;
    
    fTreeCascVarBachGrandMotherPdgCode = -9999;
    fTreeCascVarPosGrandMotherPdgCode  = -9999;
    fTreeCascVarNegGrandMotherPdgCode  = -9999;
    
    fTreeCascVarPdgCode         = -9999;
    fTreeCascVarRapMC       = -100;
    fTreeCascVarPtotMC      = -100;
    fTreeCascVarPtMC        = -100;
    
    fTreeCascVarV0PdgCode         = -9999;
    fTreeCascVarV0PtotMC      = -100;
    fTreeCascVarV0PtMC        = -100;
    
    fTreeCascVarCascadeDecayXMC = -100;
    fTreeCascVarCascadeDecayYMC = -100;
    fTreeCascVarCascadeDecayZMC = -100;
    
    fTreeCascVarV0DecayXMC      = -100;
    fTreeCascVarV0DecayYMC      = -100;
    fTreeCascVarV0DecayZMC      = -100;
    
    fTreeCascVarIsPhysicalPrimary                   = 0;
    
    fTreeCascVarIsPhysicalPrimaryBach               = kFALSE;
    fTreeCascVarIsPhysicalPrimaryPos                = kFALSE;
    fTreeCascVarIsPhysicalPrimaryNeg                = kFALSE;
    
    fTreeCascVarIsPhysicalPrimaryBachMother         = kFALSE;
    fTreeCascVarIsPhysicalPrimaryPosMother          = kFALSE;
    fTreeCascVarIsPhysicalPrimaryNegMother          = kFALSE;
    
    fTreeCascVarIsPhysicalPrimaryBachGrandMother    = kFALSE;
    fTreeCascVarIsPhysicalPrimaryPosGrandMother     = kFALSE;
    fTreeCascVarIsPhysicalPrimaryNegGrandMother     = kFALSE;
}
//_____________________________________________________________________________
void AliMCCascadeContainer::Fill(AliAODcascade* cascade, AliPIDResponse* fPIDResponse, Double_t lBestPrimaryVtxPos[3], Double_t lMagField, AliMCEvent* lMCEvent) 
{
    //------------------------------------------------
    // Initializations
    //------------------------------------------------
    if(!cascade)
    {
        printf("Error in AliMCCascadeContainer::Fill : cascade is a null ptr ! \n");
        return;
    }
    
    Double_t lPosCasc[3]            = {0., 0., 0.};// Position of cascade decay vtx
    Double_t lPosV0[3]              = {0.,0.,0.}; // Position of V0 coming from cascade
    
    Double_t lMassAsXiMinus         = 0.;
    Double_t lMassAsXiPlus          = 0.;
    Double_t lMassAsOmegaMinus      = 0.;
    Double_t lMassAsOmegaPlus       = 0.;
    Double_t lV0MassAsLambda        = 0.;
    Double_t lV0MassAsAntiLambda    = 0.;
    
    lPosCasc[0] = cascade->DecayVertexXiX();
    lPosCasc[1] = cascade->DecayVertexXiY();
    lPosCasc[2] = cascade->DecayVertexXiZ();
    
    lPosV0[0] = cascade->DecayVertexV0X();
    lPosV0[1] = cascade->DecayVertexV0Y();
    lPosV0[2] = cascade->DecayVertexV0Z();
    
    //Gather tracks informations
    AliAODTrack *bTrack = dynamic_cast<AliAODTrack*>( cascade->GetDecayVertexXi()->GetDaughter(0) );
    AliAODTrack *pTrack = dynamic_cast<AliAODTrack*>( cascade->GetDaughter(0) );
    AliAODTrack *nTrack = dynamic_cast<AliAODTrack*>( cascade->GetDaughter(1) );
    
    if (!bTrack || !pTrack || !nTrack ) {
        AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...");
        return;
    }
    
    Double_t lBachPx = 0. ; Double_t lBachPy = 0.; Double_t lBachPz = 0.;
    Double_t lPosPx = 0.; Double_t lPosPy = 0.; Double_t lPosPz = 0.;
    Double_t lNegPx = 0.; Double_t lNegPy = 0.; Double_t lNegPz = 0.;
    lBachPx = bTrack->Px() ; lBachPy = bTrack->Py() ; lBachPz = bTrack->Pz() ;
    lPosPx  = pTrack->Px() ; lPosPy  = pTrack->Py() ; lPosPz  = pTrack->Pz() ;
    lNegPx  = nTrack->Px() ; lNegPy  = nTrack->Py() ; lNegPz  = nTrack->Pz() ;
    
    Double_t lV0Px = lPosPx + lNegPx;
    Double_t lV0Py = lPosPy + lNegPy;
    Double_t lV0Pz = lPosPz + lNegPz;

    //------------------------------------------------
    // Calculation of the variables related to casc.
    //------------------------------------------------
    //Test different mass hypothesis : Xi-/Xi+ and Omega-/Omega+
    if(bTrack->Charge() < 0) lMassAsXiMinus     = cascade->MassXi();
    if(bTrack->Charge() > 0) lMassAsXiPlus      = cascade->MassXi();
    if(bTrack->Charge() < 0) lMassAsOmegaMinus  = cascade->MassOmega();
    if(bTrack->Charge() > 0) lMassAsOmegaPlus   = cascade->MassOmega();
    
    lV0MassAsLambda         = cascade->MassLambda();
    lV0MassAsAntiLambda     = cascade->MassAntiLambda();
    
    fTreeCascVarOnFlyStatus         = cascade->GetOnFlyStatus();
    fTreeCascVarChi2Cascade         = cascade->Chi2Xi();
    fTreeCascVarChi2CascadePerNDF   = ((AliAODVertex*)cascade->GetDecayVertexXi())->GetChi2perNDF();
    fTreeCascVarChi2V0              = cascade->Chi2V0();
    
    fTreeCascVarCharge      = cascade->ChargeXi();
    if(fTreeCascVarCharge > 0){
        fTreeCascVarV0Mass = lV0MassAsAntiLambda;
        fTreeCascVarMassAsXi = lMassAsXiPlus;
        fTreeCascVarMassAsOmega = lMassAsOmegaPlus;
    }
    if(fTreeCascVarCharge < 0){
        fTreeCascVarV0Mass = lV0MassAsLambda;
        fTreeCascVarMassAsXi = lMassAsXiMinus;
        fTreeCascVarMassAsOmega = lMassAsOmegaMinus;
    }
    fTreeCascVarV0MassAsLambda      = lV0MassAsLambda;
    fTreeCascVarV0MassAsAntiLambda  = lV0MassAsAntiLambda;
    
    fTreeCascVarPtot        = TMath::Sqrt( cascade->Ptot2Xi() );
    fTreeCascVarV0Ptot      = TMath::Sqrt( TMath::Power(lPosPx+lNegPx,2) 
                                            + TMath::Power(lPosPy+lNegPy,2) 
                                            + TMath::Power(lPosPz+lNegPz,2) );
    fTreeCascVarPt          = TMath::Sqrt( cascade->Pt2Xi() );
    fTreeCascVarV0Pt        = TMath::Sqrt( TMath::Power(lPosPx+lNegPx,2) 
                                            + TMath::Power(lPosPy+lNegPy,2) );
    fTreeCascVarRapXi       = cascade->RapXi();
    fTreeCascVarRapOmega    = cascade->RapOmega();
    
    fTreeCascVarCascEta     = 0.5*TMath::Log(TMath::Abs((fTreeCascVarPtot + (lV0Pz + lBachPz) )/(fTreeCascVarPtot - (lV0Pz + lBachPz) + 1e-13)));
    fTreeCascVarV0Eta       = 0.5*TMath::Log(TMath::Abs((fTreeCascVarV0Ptot + lV0Pz)/( fTreeCascVarV0Ptot - lV0Pz + 1e-13 )));
    
    fTreeCascVarBachEta     = bTrack->Eta();
    fTreeCascVarPosEta      = pTrack->Eta();
    fTreeCascVarNegEta      = nTrack->Eta();
    
    fTreeCascVarAlpha   = cascade->AlphaXi();
    fTreeCascVarPtArm   = cascade->PtArmXi();
    fTreeCascVarAlphaV0 = cascade->AlphaV0();
    fTreeCascVarPtArmV0 = cascade->PtArmV0();
    
    fTreeCascVarDCACascDau  = cascade->DcaXiDaughters();
    fTreeCascVarDCAV0Dau    = cascade->DcaV0Daughters();
    fTreeCascVarDCAV0ToPV   = cascade->DcaV0ToPrimVertex();
    
    
    bTrack->GetImpactParameters(fTreeCascVarDCABachToPV,fTreeCascVarBachDCAz);
    pTrack->GetImpactParameters(fTreeCascVarDCAPosToPV,fTreeCascVarPosDCAz);
    nTrack->GetImpactParameters(fTreeCascVarDCANegToPV,fTreeCascVarNegDCAz);
    //fTreeCascVarDCABachToPV = cascade->DcaBachToPrimVertex();
    //fTreeCascVarDCAPosToPV  = cascade->DcaPosToPrimVertex();
    //fTreeCascVarDCANegToPV  = cascade->DcaNegToPrimVertex();
    
    
    fTreeCascVarCascRadius  = TMath::Sqrt( lPosCasc[0]*lPosCasc[0]  +  lPosCasc[1]*lPosCasc[1] );
    fTreeCascVarV0Radius    = TMath::Sqrt( lPosV0[0]*lPosV0[0]  +  lPosV0[1]*lPosV0[1] );
    fTreeCascVarCascCosPA   = cascade->CosPointingAngleXi(lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
    fTreeCascVarV0CosPA     = cascade->CosPointingAngle(lBestPrimaryVtxPos);
    
    //FIX ME
    //if( bTrack->Charge() < 0 ) fTreeCascVarWrongCosPA = GetCosPA( bTrack , pTrack, lAODevent);
    //if( bTrack->Charge() > 0 ) fTreeCascVarWrongCosPA = GetCosPA( bTrack , nTrack, lAODevent);
    fTreeCascVarWrongCosPA         = cascade->BachBaryonCosPA();
    fTreeCascVarCascCosPASpecial   = cascade->CosPointingAngle(lPosCasc);
    
    //________________________________________________________________________
    // Track quality cuts
    //Least Nbr of Crossed Rows
    Int_t lBachNbrCrossedRows = bTrack->GetTPCNCrossedRows();
    Int_t lPosNbrCrossedRows  = pTrack->GetTPCNCrossedRows();
    Int_t lNegNbrCrossedRows  = nTrack->GetTPCNCrossedRows();
    
    Int_t lLeastNbrCrossedRows = (Int_t) lPosNbrCrossedRows;
    if( lNegNbrCrossedRows < lLeastNbrCrossedRows ) lLeastNbrCrossedRows = (Int_t) lNegNbrCrossedRows;
    if( lBachNbrCrossedRows < lLeastNbrCrossedRows ) lLeastNbrCrossedRows = (Int_t) lBachNbrCrossedRows;
    
    fTreeCascVarLeastNbrCrossedRows = lLeastNbrCrossedRows;
    
    //Compute ratio Crossed Rows / Findable clusters
    //Note: above test avoids division by zero!
    Double_t lBachTrackCrossedRowsOverFindable = lBachNbrCrossedRows / ((Double_t)(bTrack->GetTPCNclsF() + 1e-5));
    Double_t lPosTrackCrossedRowsOverFindable = lPosNbrCrossedRows / ((Double_t)(pTrack->GetTPCNclsF() + 1e-5 ));
    Double_t lNegTrackCrossedRowsOverFindable = lNegNbrCrossedRows / ((Double_t)(nTrack->GetTPCNclsF() + 1e-5));
    
    fTreeCascVarLeastRatioCrossedRowsOverFindable = lBachTrackCrossedRowsOverFindable;
    if( lPosTrackCrossedRowsOverFindable < fTreeCascVarLeastRatioCrossedRowsOverFindable ) 
        fTreeCascVarLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
    if( lNegTrackCrossedRowsOverFindable < fTreeCascVarLeastRatioCrossedRowsOverFindable ) 
        fTreeCascVarLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;
    
    // Least Nbr of clusters        
    Int_t lPosTPCClusters        = 0.;
    Int_t lNegTPCClusters        = 0.;
    Int_t lBachTPCClusters       = 0.;

    lPosTPCClusters   = pTrack->GetTPCNcls();
    lNegTPCClusters   = nTrack->GetTPCNcls();
    lBachTPCClusters  = bTrack->GetTPCNcls();
    
    Int_t lLeastNbrOfClusters    = lBachTPCClusters;
    if( lPosTPCClusters < lLeastNbrOfClusters ) lLeastNbrOfClusters = lPosTPCClusters;
    if( lNegTPCClusters < lLeastNbrOfClusters ) lLeastNbrOfClusters = lNegTPCClusters;
    
    fTreeCascVarLeastNbrClusters = lLeastNbrOfClusters;
    
    //Min track length
    Double_t lPosTrackLength        = 0.;
    Double_t lNegTrackLength        = 0.;
    Double_t lBachTrackLength       = 0.;
    
    lBachTrackLength   = GetLengthInActiveZone(bTrack, 2.0, 220.0, lMagField);
    lPosTrackLength    = GetLengthInActiveZone(pTrack, 2.0, 220.0, lMagField);
    lNegTrackLength    = GetLengthInActiveZone(nTrack, 2.0, 220.0, lMagField);
    
    Double_t lSmallestTrackLength = lBachTrackLength;
    if ( lPosTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lPosTrackLength;
    if ( lNegTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lNegTrackLength;
    
    fTreeCascVarMinTrackLength = lSmallestTrackLength;
    
    //Nbr Of Crossed Rows Over Length
    Double_t lBachTrackNbrCrOverLength = bTrack->GetTPCClusterInfo(2,1)/(lBachTrackLength-TMath::Max(fTreeCascVarCascRadius-85.,0.) + 1e-5);
    Double_t lPosTrackNbrCrOverLength = pTrack->GetTPCClusterInfo(2,1)/(lPosTrackLength-TMath::Max(fTreeCascVarV0Radius-85.,0.) + 1e-5);
    Double_t lNegTrackNbrCrOverLength = nTrack->GetTPCClusterInfo(2,1)/(lNegTrackLength-TMath::Max(fTreeCascVarV0Radius-85.,0.) + 1e-5);

    Double_t lLeastNbrCrOverLength = (Double_t) lPosTrackNbrCrOverLength;
    if( lNegTrackNbrCrOverLength < lLeastNbrCrOverLength ) lLeastNbrCrOverLength = (Float_t) lNegTrackNbrCrOverLength;
    if( lBachTrackNbrCrOverLength < lLeastNbrCrOverLength ) lLeastNbrCrOverLength = (Float_t) lBachTrackNbrCrOverLength;
    
    fTreeCascVarNbrCrossedRowsOverLength = lLeastNbrCrOverLength;
    
    //Max Chi2 per cluster
    Double_t lBiggestChi2PerCluster = 0.;
    Double_t lPosChi2PerCluster = pTrack->GetTPCchi2() / ((Float_t) lPosTPCClusters + 1e-13);
    Double_t lNegChi2PerCluster = nTrack->GetTPCchi2() / ((Float_t) lNegTPCClusters + 1e-13);
    Double_t lBachChi2PerCluster = bTrack->GetTPCchi2() / ((Float_t) lBachTPCClusters + 1e-13);
    
    if( lPosChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lPosChi2PerCluster;
    if( lNegChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lNegChi2PerCluster;
    if( lBachChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lBachChi2PerCluster;
    
    fTreeCascVarMaxChi2PerCluster = lBiggestChi2PerCluster;
    
    //________________________________________________________________________
    // Decay length info
    //Cascade
    fTreeCascVarDistOverTotMom  = TMath::Sqrt( TMath::Power( lPosCasc[0] - lBestPrimaryVtxPos[0] , 2) 
                                                + TMath::Power( lPosCasc[1] - lBestPrimaryVtxPos[1] , 2) 
                                                + TMath::Power( lPosCasc[2] - lBestPrimaryVtxPos[2] , 2) );
    fTreeCascVarDistOverTotMom /= fTreeCascVarPtot;
    
    //V0
    fTreeCascVarV0DistOverTotMom    = TMath::Sqrt( TMath::Power( lPosV0[0] - lPosCasc[0] , 2) 
                                        + TMath::Power( lPosV0[1] - lPosCasc[1] , 2) 
                                        + TMath::Power( lPosV0[2] - lPosCasc[2] , 2) );
    fTreeCascVarV0DistOverTotMom /= fTreeCascVarV0Ptot;
    
    
    //------------------------------------------------
    // TPC dEdx information
    //------------------------------------------------
    fTreeCascVarBachNSigmaPion  = fPIDResponse->NumberOfSigmasTPC( bTrack, AliPID::kPion );
    fTreeCascVarBachNSigmaKaon  = fPIDResponse->NumberOfSigmasTPC( bTrack, AliPID::kKaon );
    fTreeCascVarPosNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
    fTreeCascVarPosNSigmaProton = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
    fTreeCascVarNegNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion   );
    fTreeCascVarNegNSigmaProton = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
    
    //------------------------------------------------
    // ITS info 
    //------------------------------------------------
    fTreeCascVarBachITSNSigmaPion  = fPIDResponse->NumberOfSigmasITS( bTrack, AliPID::kPion );
    fTreeCascVarBachITSNSigmaKaon  = fPIDResponse->NumberOfSigmasITS( bTrack, AliPID::kKaon );
    fTreeCascVarPosITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( pTrack, AliPID::kPion );
    fTreeCascVarPosITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( pTrack, AliPID::kProton );
    fTreeCascVarNegITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( nTrack, AliPID::kPion   );
    fTreeCascVarNegITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( nTrack, AliPID::kProton );
    
    //------------------------------------------------
    // TOF info 
    //------------------------------------------------
    fTreeCascVarBachTOFNSigmaPion  = fPIDResponse->NumberOfSigmasTOF( bTrack, AliPID::kPion );
    fTreeCascVarBachTOFNSigmaKaon  = fPIDResponse->NumberOfSigmasTOF( bTrack, AliPID::kKaon );
    fTreeCascVarPosTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( pTrack, AliPID::kPion );
    fTreeCascVarPosTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( pTrack, AliPID::kProton );
    fTreeCascVarNegTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( nTrack, AliPID::kPion   );
    fTreeCascVarNegTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( nTrack, AliPID::kProton );
    
    //------------------------------------------------
    // Raw TPC dEdx + PIDForTracking information
    //------------------------------------------------
    
    //Step 2: Acquire TPC Signals
    fTreeCascVarBachdEdx    = bTrack->GetTPCsignal();
    fTreeCascVarPosdEdx     = pTrack->GetTPCsignal();
    fTreeCascVarNegdEdx     = nTrack->GetTPCsignal();
    
    fTreeCascVarBachdEdxN   = bTrack->GetTPCsignalN();
    fTreeCascVarPosdEdxN    = pTrack->GetTPCsignalN();
    fTreeCascVarNegdEdxN    = nTrack->GetTPCsignalN();
    
    
    //Step 3: Acquire PID For Tracking
    fTreeCascVarBachPIDForTracking  = bTrack->GetPIDForTracking();
    fTreeCascVarPosPIDForTracking   = pTrack->GetPIDForTracking();
    fTreeCascVarNegPIDForTracking   = nTrack->GetPIDForTracking();
    
    
    //________________________________________________________________________
    // Track status
    fTreeCascVarBachTrackStatus     = bTrack->GetStatus();
    fTreeCascVarPosTrackStatus      = pTrack->GetStatus();
    fTreeCascVarNegTrackStatus      = nTrack->GetStatus();
    
    //________________________________________________________________________
    // Track DCAz
	// Done above
    //fTreeCascVarBachDCAz = GetDCAz(bTrack);
    //fTreeCascVarPosDCAz = GetDCAz(pTrack);
    //fTreeCascVarNegDCAz = GetDCAz(nTrack);
    
    //________________________________________________________________________
    // Momentum info
    fTreeCascVarBachPx = lBachPx ;   fTreeCascVarBachPy = lBachPy ;    fTreeCascVarBachPz = lBachPz ;
    fTreeCascVarPosPx  = lPosPx  ;   fTreeCascVarPosPy  = lPosPy  ;    fTreeCascVarPosPz  = lPosPz ;
    fTreeCascVarNegPx  = lNegPx  ;   fTreeCascVarNegPy  = lNegPy  ;    fTreeCascVarNegPz  = lNegPz ;
    
    //________________________________________________________________________
    // Decay vtx info
    fTreeCascVarCascadeDecayX = lPosCasc[0] ; 
    fTreeCascVarCascadeDecayY = lPosCasc[1] ;  
    fTreeCascVarCascadeDecayZ = lPosCasc[2] ; 
    
    fTreeCascVarV0DecayX = lPosV0[0] ; 
    fTreeCascVarV0DecayY = lPosV0[1] ;  
    fTreeCascVarV0DecayZ = lPosV0[2] ; 
    
    //________________________________________________________________________
    // Particle index
    fTreeCascVarBachLabel   = bTrack->GetLabel();
    fTreeCascVarPosLabel    = pTrack->GetLabel();
    fTreeCascVarNegLabel    = nTrack->GetLabel();
    
    //________________________________________________________________________
    //Check clusters
    fTreeCascVarPosITSClusters0 = pTrack->HasPointOnITSLayer(0);
    fTreeCascVarPosITSClusters1 = pTrack->HasPointOnITSLayer(1);
    fTreeCascVarPosITSClusters2 = pTrack->HasPointOnITSLayer(2);
    fTreeCascVarPosITSClusters3 = pTrack->HasPointOnITSLayer(3);
    fTreeCascVarPosITSClusters4 = pTrack->HasPointOnITSLayer(4);
    fTreeCascVarPosITSClusters5 = pTrack->HasPointOnITSLayer(5);
    
    fTreeCascVarNegITSClusters0 = nTrack->HasPointOnITSLayer(0);
    fTreeCascVarNegITSClusters1 = nTrack->HasPointOnITSLayer(1);
    fTreeCascVarNegITSClusters2 = nTrack->HasPointOnITSLayer(2);
    fTreeCascVarNegITSClusters3 = nTrack->HasPointOnITSLayer(3);
    fTreeCascVarNegITSClusters4 = nTrack->HasPointOnITSLayer(4);
    fTreeCascVarNegITSClusters5 = nTrack->HasPointOnITSLayer(5);
    
    fTreeCascVarBachITSClusters0 = bTrack->HasPointOnITSLayer(0);
    fTreeCascVarBachITSClusters1 = bTrack->HasPointOnITSLayer(1);
    fTreeCascVarBachITSClusters2 = bTrack->HasPointOnITSLayer(2);
    fTreeCascVarBachITSClusters3 = bTrack->HasPointOnITSLayer(3);
    fTreeCascVarBachITSClusters4 = bTrack->HasPointOnITSLayer(4);
    fTreeCascVarBachITSClusters5 = bTrack->HasPointOnITSLayer(5);
    
    //________________________________________________________________________
    //Check its clusters, shared
    fTreeCascVarPosITSSharedClusters0 = pTrack->HasSharedPointOnITSLayer(0);
    fTreeCascVarPosITSSharedClusters1 = pTrack->HasSharedPointOnITSLayer(1);
    fTreeCascVarPosITSSharedClusters2 = pTrack->HasSharedPointOnITSLayer(2);
    fTreeCascVarPosITSSharedClusters3 = pTrack->HasSharedPointOnITSLayer(3);
    fTreeCascVarPosITSSharedClusters4 = pTrack->HasSharedPointOnITSLayer(4);
    fTreeCascVarPosITSSharedClusters5 = pTrack->HasSharedPointOnITSLayer(5);
    
    fTreeCascVarNegITSSharedClusters0 = nTrack->HasSharedPointOnITSLayer(0);
    fTreeCascVarNegITSSharedClusters1 = nTrack->HasSharedPointOnITSLayer(1);
    fTreeCascVarNegITSSharedClusters2 = nTrack->HasSharedPointOnITSLayer(2);
    fTreeCascVarNegITSSharedClusters3 = nTrack->HasSharedPointOnITSLayer(3);
    fTreeCascVarNegITSSharedClusters4 = nTrack->HasSharedPointOnITSLayer(4);
    fTreeCascVarNegITSSharedClusters5 = nTrack->HasSharedPointOnITSLayer(5);
    
    fTreeCascVarBachITSSharedClusters0 = bTrack->HasSharedPointOnITSLayer(0);
    fTreeCascVarBachITSSharedClusters1 = bTrack->HasSharedPointOnITSLayer(1);
    fTreeCascVarBachITSSharedClusters2 = bTrack->HasSharedPointOnITSLayer(2);
    fTreeCascVarBachITSSharedClusters3 = bTrack->HasSharedPointOnITSLayer(3);
    fTreeCascVarBachITSSharedClusters4 = bTrack->HasSharedPointOnITSLayer(4);
    fTreeCascVarBachITSSharedClusters5 = bTrack->HasSharedPointOnITSLayer(5);
    
    //________________________________________________________________________
    //GetKinkIndex condition
    fTreeCascVarBachIsKink  = kFALSE;
    fTreeCascVarPosIsKink   = kFALSE;
    fTreeCascVarNegIsKink   = kFALSE;
	if( bTrack->GetProdVertex()->GetType() == AliAODVertex::kKink ) fTreeCascVarBachIsKink  = kTRUE;
    if( pTrack->GetProdVertex()->GetType() == AliAODVertex::kKink ) fTreeCascVarPosIsKink   = kTRUE;
    if( nTrack->GetProdVertex()->GetType() == AliAODVertex::kKink ) fTreeCascVarNegIsKink   = kTRUE;
    
    //________________________________________________________________________
    //TOF info
    fTreeCascVarBachTOFExpTDiff     = bTrack->GetTOFExpTDiff( lMagField );
    fTreeCascVarPosTOFExpTDiff      = pTrack->GetTOFExpTDiff( lMagField );
    fTreeCascVarNegTOFExpTDiff      = nTrack->GetTOFExpTDiff( lMagField );
    
    fTreeCascVarBachTOFSignal       = bTrack->GetTOFsignal() * 1.e-3; // in ns
    fTreeCascVarPosTOFSignal        = pTrack->GetTOFsignal() * 1.e-3; // in ns
    fTreeCascVarNegTOFSignal        = nTrack->GetTOFsignal() * 1.e-3; // in ns
    
    fTreeCascVarBachTOFBCid         = bTrack->GetTOFBunchCrossing( lMagField );
    fTreeCascVarPosTOFBCid          = pTrack->GetTOFBunchCrossing( lMagField );
    fTreeCascVarNegTOFBCid          = nTrack->GetTOFBunchCrossing( lMagField );
    
    fTreeCascVarBachTOFLength       = bTrack->GetIntegratedLength();
    fTreeCascVarPosTOFLength        = pTrack->GetIntegratedLength();
    fTreeCascVarNegTOFLength        = nTrack->GetIntegratedLength();
    
    fTreeCascVarBachTOFDeltaX       = bTrack->GetTOFsignalDx();
    fTreeCascVarPosTOFDeltaX        = pTrack->GetTOFsignalDx();
    fTreeCascVarNegTOFDeltaX        = nTrack->GetTOFsignalDx();
    
    fTreeCascVarBachTOFDeltaZ       = bTrack->GetTOFsignalDz();
    fTreeCascVarPosTOFDeltaZ        = pTrack->GetTOFsignalDz();
    fTreeCascVarNegTOFDeltaZ        = nTrack->GetTOFsignalDz();
    
    //________________________________________________________________________
    //Cowboy/sailor info regarding V0 inside cascade
    //Calculate vec prod with momenta projected to xy plane
    //Provisions for cowboy/sailor check
    Double_t lModp1 = TMath::Sqrt( lPosPx*lPosPx + lPosPy*lPosPy );
    Double_t lModp2 = TMath::Sqrt( lNegPx*lNegPx + lNegPy*lNegPy );
    
    Double_t lVecProd = (lPosPx*lNegPy - lPosPy*lNegPx) / (lModp1*lModp2);
    
    if ( lMagField < 0 ) lVecProd *= -1; //invert sign
    
    fTreeCascVarCowboyness = lVecProd;
    
    fTreeCascVarIsCowboy = kFALSE;
    if (lVecProd < 0) fTreeCascVarIsCowboy = kTRUE;
    
    Double_t lBachMod = TMath::Sqrt(lBachPx*lBachPx + lBachPy*lBachPy);
    Double_t lVecProdCasc = (lV0Px*lBachPy - lV0Py*lBachPx) / (TMath::Sqrt(lV0Px*lV0Px + lV0Py*lV0Py)*lBachMod);
    
    if ( lMagField < 0 ) lVecProdCasc *= -1; //invert sign
    
    fTreeCascVarCascadeCowboyness = lVecProdCasc;
    
    fTreeCascVarIsCascadeCowboy = kFALSE;
    if (lVecProdCasc < 0) fTreeCascVarIsCascadeCowboy = kTRUE;

    //------------------------------------------------
    // Associate Cascade Candidates to Monte Carlo!
    //------------------------------------------------
    
    //________________________________________________________________________
    // Get MC tracks
    // Abs value = needed ! question of quality track association ...
    Int_t lBachLabel    = (Int_t) TMath::Abs( bTrack->GetLabel() );
    Int_t lPosLabel     = (Int_t) TMath::Abs( pTrack->GetLabel() );
    Int_t lNegLabel     = (Int_t) TMath::Abs( nTrack->GetLabel() );
    
    AliMCParticle* bTrackMC = (AliMCParticle*)lMCEvent->GetTrack( lBachLabel );
    AliMCParticle* pTrackMC = (AliMCParticle*)lMCEvent->GetTrack( lPosLabel  );
    AliMCParticle* nTrackMC = (AliMCParticle*)lMCEvent->GetTrack( lNegLabel  );
    
    if (!bTrackMC || !pTrackMC || !nTrackMC ) {
        AliWarning("ERROR: Could not retrieve one of the 3 AOD MC daughter tracks of the cascade ...");
        return;
    }

    if( bTrackMC->IsPhysicalPrimary() ) fTreeCascVarIsPhysicalPrimaryBach = kTRUE;
    if( pTrackMC->IsPhysicalPrimary() ) fTreeCascVarIsPhysicalPrimaryPos  = kTRUE;
    if( nTrackMC->IsPhysicalPrimary() ) fTreeCascVarIsPhysicalPrimaryNeg  = kTRUE;
    
    //________________________________________________________________________
    // Pseudo-rapidity info
    fTreeCascVarBachEtaMC = bTrackMC->Eta() ;
    fTreeCascVarPosEtaMC  = pTrackMC->Eta() ;
    fTreeCascVarNegEtaMC  = nTrackMC->Eta() ;
    
    //________________________________________________________________________
    // Momentum info
    fTreeCascVarBachPxMC = bTrackMC->Px() ;   fTreeCascVarBachPyMC = bTrackMC->Py() ;    fTreeCascVarBachPzMC = bTrackMC->Pz() ;
    fTreeCascVarPosPxMC  = pTrackMC->Px() ;   fTreeCascVarPosPyMC  = pTrackMC->Py() ;    fTreeCascVarPosPzMC  = pTrackMC->Pz() ;
    fTreeCascVarNegPxMC  = nTrackMC->Px() ;   fTreeCascVarNegPyMC  = nTrackMC->Py() ;    fTreeCascVarNegPzMC  = nTrackMC->Pz() ;
    
    //________________________________________________________________________
    // PID & label info
    
    fTreeCascVarBachPdgCode = bTrackMC->PdgCode();
    fTreeCascVarPosPdgCode  = pTrackMC->PdgCode();
    fTreeCascVarNegPdgCode  = nTrackMC->PdgCode();
    
    Int_t lBachMotherLabel  = bTrackMC->GetMother();
    Int_t lPosMotherLabel   = pTrackMC->GetMother();
    Int_t lNegMotherLabel   = nTrackMC->GetMother();
    
    
    // Extra: check mother particle pdg code completely independently
    if ( lBachMotherLabel >= 0 )
    {
        AliMCParticle* lBachMother = (AliMCParticle*) lMCEvent->GetTrack( lBachMotherLabel );
        
        if( lBachMother->IsPhysicalPrimary() ) fTreeCascVarIsPhysicalPrimaryBachMother = kTRUE;
        //Mother PdgCode
        fTreeCascVarBachMotherPdgCode   = lBachMother->PdgCode();
        fTreeCascVarBachMotherLabel     = lBachMotherLabel;
        
        //Go further than that, please
        Int_t lGrandMotherLabel = lBachMother->GetMother();
        if( lGrandMotherLabel >= 0 )
        {
            AliMCParticle* lBachGrandMother = (AliMCParticle*) lMCEvent->GetTrack( lGrandMotherLabel );
            if( lBachGrandMother->IsPhysicalPrimary() ) fTreeCascVarIsPhysicalPrimaryBachGrandMother = kTRUE;
            // Grand Mother PdgCode
            fTreeCascVarBachGrandMotherPdgCode      = lBachGrandMother->PdgCode();
            fTreeCascVarBachGrandMotherLabel        = lGrandMotherLabel;
        }
    }
    
    if ( lPosMotherLabel >= 0 )
    {
        AliMCParticle* lPosMother = (AliMCParticle*) lMCEvent->GetTrack( lPosMotherLabel );
        if( lPosMother->IsPhysicalPrimary() ) fTreeCascVarIsPhysicalPrimaryPosMother = kTRUE;
        //Mother PdgCode
        fTreeCascVarPosMotherPdgCode   = lPosMother->PdgCode();
        fTreeCascVarPosMotherLabel     = lPosMotherLabel;
        
        //Go further than that, please
        Int_t lGrandMotherLabel = lPosMother->GetMother();
        if( lGrandMotherLabel >= 0 )
        {
            AliMCParticle* lPosGrandMother = (AliMCParticle*) lMCEvent->GetTrack( lGrandMotherLabel );
            if( lPosGrandMother->IsPhysicalPrimary() ) fTreeCascVarIsPhysicalPrimaryPosGrandMother = kTRUE;
            // Grand Mother PdgCode
            fTreeCascVarPosGrandMotherPdgCode      = lPosGrandMother->PdgCode();
            fTreeCascVarPosGrandMotherLabel        = lGrandMotherLabel;
        }
    }
    
    if ( lNegMotherLabel >= 0 )
    {
        AliMCParticle* lNegMother = (AliMCParticle*) lMCEvent->GetTrack( lNegMotherLabel );
        if( lNegMother->IsPhysicalPrimary() ) fTreeCascVarIsPhysicalPrimaryNegMother = kTRUE;
        //Mother PdgCode
        fTreeCascVarNegMotherPdgCode   = lNegMother->PdgCode();
        fTreeCascVarNegMotherLabel     = lNegMotherLabel;
        
        //Go further than that, please
        Int_t lGrandMotherLabel = lNegMother->GetMother();
        if( lGrandMotherLabel >= 0 )
        {
            AliMCParticle* lNegGrandMother = (AliMCParticle*) lMCEvent->GetTrack( lGrandMotherLabel );
            if( lNegGrandMother->IsPhysicalPrimary() ) fTreeCascVarIsPhysicalPrimaryNegGrandMother = kTRUE;
            // Grand Mother PdgCode
            fTreeCascVarNegGrandMotherPdgCode      = lNegGrandMother->PdgCode();
            fTreeCascVarNegGrandMotherLabel        = lGrandMotherLabel;
        }
    }
    
    if( lPosMotherLabel == lNegMotherLabel && lPosMotherLabel >= 0) // same mother & mother != primary (!= -1)
    { 
        // mothers = Lambda candidate ... a priori
        AliMCParticle* pMotherTrackMC = (AliMCParticle*)lMCEvent->GetTrack( lPosMotherLabel );
        AliMCParticle* nMotherTrackMC = (AliMCParticle*)lMCEvent->GetTrack( lNegMotherLabel );
        
        // Get V0 MC info
        fTreeCascVarV0PdgCode       = pMotherTrackMC->PdgCode();
        fTreeCascVarV0PtotMC        = pMotherTrackMC->P();
        fTreeCascVarV0PtMC          = pMotherTrackMC->Pt();
        
        Int_t lPosGrandMotherLabel = pMotherTrackMC->GetMother() ;
        Int_t lNegGrandMotherLabel = nMotherTrackMC->GetMother() ;
        
        if( lPosGrandMotherLabel == lNegGrandMotherLabel && lPosGrandMotherLabel >= 0) // primary lambda ...
        {
            // Grand mother = Cascade candidate ... a priori
            Int_t lBachMotherLabel = (Int_t) bTrackMC->GetMother();
            
            if( lBachMotherLabel == lPosGrandMotherLabel ) //same mother for bach and V0 daughters
            { 
                AliMCParticle* bMotherTrackMC        = (AliMCParticle*) lMCEvent->GetTrack( lBachMotherLabel );
                AliMCParticle* pGrandMotherTrackMC   = (AliMCParticle*) lMCEvent->GetTrack( lPosGrandMotherLabel );
                AliMCParticle* nGrandMotherTrackMC   = (AliMCParticle*) lMCEvent->GetTrack( lNegGrandMotherLabel );
                
                Int_t lBachMotherPdgCode        = bMotherTrackMC     ->PdgCode();
                Int_t lPosGrandMotherPdgCode    = pGrandMotherTrackMC->PdgCode();
                Int_t lNegGrandMotherPdgCode    = nGrandMotherTrackMC->PdgCode();
                
                if( lBachMotherPdgCode == lNegGrandMotherPdgCode && lBachMotherPdgCode == lPosGrandMotherPdgCode ) 
                {
                    // Get Cascade MC info
                    fTreeCascVarPdgCode = lBachMotherPdgCode;
                    
                    if( TMath::Abs( fTreeCascVarPdgCode ) == 3312 || TMath::Abs( fTreeCascVarPdgCode ) == 3334 )
                        fTreeCascVarRapMC  = bMotherTrackMC->Y();
                    fTreeCascVarPtotMC     = bMotherTrackMC->P();
                    fTreeCascVarPtMC       = bMotherTrackMC->Pt();
                    
                    if( bMotherTrackMC->IsPhysicalPrimary()         ) fTreeCascVarIsPhysicalPrimary = 1; //Is Primary!
                    if( bMotherTrackMC->IsSecondaryFromWeakDecay()  ) fTreeCascVarIsPhysicalPrimary = 2; //Weak Decay!
                    if( bMotherTrackMC->IsSecondaryFromMaterial ()  ) fTreeCascVarIsPhysicalPrimary = 3; //From Material!
                    
                }
                
            }
            
        }
        
    }
    
    //________________________________________________________________________
    // Decay position info
    //Be careful: Vx, Vy, Vz: Creation vertex. So decay position is the
    //Creation vertex of any one of the daughters!
    fTreeCascVarCascadeDecayXMC = bTrackMC->Xv();
    fTreeCascVarCascadeDecayYMC = bTrackMC->Yv();
    fTreeCascVarCascadeDecayZMC = bTrackMC->Zv();
    
    fTreeCascVarV0DecayXMC = pTrackMC->Xv();
    fTreeCascVarV0DecayYMC = pTrackMC->Yv();
    fTreeCascVarV0DecayZMC = pTrackMC->Zv();
}

//_____________________________________________________________________________
Float_t AliMCCascadeContainer::GetLengthInActiveZone( AliAODTrack *gt, Float_t deltaY, Float_t deltaZ, Float_t b ){
    AliESDtrack esdTrack( gt );
    esdTrack.SetESDEvent((AliESDEvent*) gt->GetEvent() );
    AliExternalTrackParam etp;
    etp.CopyFromVTrack(gt);
    esdTrack.ResetTrackParamIp(&etp);
    return esdTrack.GetLengthInActiveZone(1, deltaY, deltaZ, b);
}
//_____________________________________________________________________________
Float_t AliMCCascadeContainer::GetDCAz(AliAODTrack *lTrack)
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
Bool_t AliMCCascadeContainer::IsSelected( const AliCascadeResult* lCascadeResult )
{
    Double_t lLambdaPDGMass = 1.115683;
    Double_t lXiPDGMass     = 1.32171;
    Double_t lOmegaPDGMass  = 1.67245;
    
    Int_t       lCascResultCharge;
    Bool_t      lCascResultIsOmega;
    
    Double_t    lInvMass;
    Double_t    lPDGMass;
    Double_t    lV0Mass;
    Double_t    lRap;
    Double_t    lBachTPCNSigma;
    Double_t    lPosTPCNSigma;
    Double_t    lNegTPCNSigma;
    
    
    //========================================================================
    if( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kXiMinus     )
    {
        lCascResultCharge       = -1;
        lCascResultIsOmega      = kFALSE;
        
        lPDGMass                = lXiPDGMass;
        lInvMass                = fTreeCascVarMassAsXi;
        lV0Mass                 = fTreeCascVarV0MassAsLambda;
        lRap                    = fTreeCascVarRapXi;
        
        lBachTPCNSigma          = fTreeCascVarBachNSigmaPion;
        lPosTPCNSigma           = fTreeCascVarPosNSigmaProton;
        lNegTPCNSigma           = fTreeCascVarNegNSigmaPion;
    }
    if( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kXiPlus      )
    {
        lCascResultCharge       = +1;
        lCascResultIsOmega      = kFALSE;
        
        lPDGMass                = lXiPDGMass;
        lInvMass                = fTreeCascVarMassAsXi;
        lV0Mass                 = fTreeCascVarV0MassAsAntiLambda;
        lRap                    = fTreeCascVarRapXi;
        
        lBachTPCNSigma          = fTreeCascVarBachNSigmaPion;
        lPosTPCNSigma           = fTreeCascVarPosNSigmaPion;
        lNegTPCNSigma           = fTreeCascVarNegNSigmaProton;
    }
    if( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kOmegaMinus     )
    {
        lCascResultCharge       = -1;
        lCascResultIsOmega      = kTRUE;
        
        lPDGMass                = lOmegaPDGMass;
        lInvMass                = fTreeCascVarMassAsOmega;
        lV0Mass                 = fTreeCascVarV0MassAsLambda;
        lRap                    = fTreeCascVarRapOmega;
        
        lBachTPCNSigma          = fTreeCascVarBachNSigmaKaon;
        lPosTPCNSigma           = fTreeCascVarPosNSigmaProton;
        lNegTPCNSigma           = fTreeCascVarNegNSigmaPion;
    }
    if( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kOmegaPlus      )
    {
        lCascResultCharge       = +1;
        lCascResultIsOmega      = kTRUE;
        
        lPDGMass                = lOmegaPDGMass;
        lInvMass                = fTreeCascVarMassAsOmega;
        lV0Mass                 = fTreeCascVarV0MassAsAntiLambda;
        lRap                    = fTreeCascVarRapOmega;
        
        lBachTPCNSigma          = fTreeCascVarBachNSigmaKaon;
        lPosTPCNSigma           = fTreeCascVarPosNSigmaPion;
        lNegTPCNSigma           = fTreeCascVarNegNSigmaProton;
    }
    
    
    if(
        //Check 1 : Charge
        fTreeCascVarCharge == lCascResultCharge                                                     &&
        
        //Check 2: Basic Acceptance cuts
        // Bachelor
        fTreeCascVarBachEta                     > lCascadeResult->GetCutMinEtaTracks()              && 
        fTreeCascVarBachEta                     < lCascadeResult->GetCutMaxEtaTracks()              &&
        // Positive
        fTreeCascVarPosEta                      > lCascadeResult->GetCutMinEtaTracks()              && 
        fTreeCascVarPosEta                      < lCascadeResult->GetCutMaxEtaTracks()              &&
        // Negative
        fTreeCascVarNegEta                      > lCascadeResult->GetCutMinEtaTracks()              &&
        fTreeCascVarNegEta                      < lCascadeResult->GetCutMaxEtaTracks()              &&
        
        // Corresponding rapidity
        lRap                                    > lCascadeResult->GetCutMinRapidity()               &&
        lRap                                    < lCascadeResult->GetCutMaxRapidity()               &&
        
        //Check 3: Topological Variables
        // - V0 Selections
        TMath::Abs(fTreeCascVarDCAPosToPV)      > lCascadeResult->GetCutDCAPosToPV()                &&
        TMath::Abs(fTreeCascVarDCANegToPV)      > lCascadeResult->GetCutDCANegToPV()                &&
        TMath::Abs(fTreeCascVarDCAV0Dau)        < lCascadeResult->GetCutDCAV0Daughters()            &&
        TMath::Abs(fTreeCascVarDCAV0ToPV)       > lCascadeResult->GetCutDCAV0ToPV()                 &&
        fTreeCascVarV0Radius                    > lCascadeResult->GetCutV0Radius()                  &&
        fTreeCascVarV0CosPA                     > lCascadeResult->GetCutV0CosPA()                   &&
        // - Cascade Selections
        TMath::Abs(lV0Mass-lLambdaPDGMass)      < lCascadeResult->GetCutV0Mass()                    &&
        TMath::Abs(fTreeCascVarDCABachToPV)     > lCascadeResult->GetCutDCABachToPV()               &&
        TMath::Abs(fTreeCascVarDCACascDau)      < lCascadeResult->GetCutDCACascDaughters()          &&
        fTreeCascVarCascRadius                  > lCascadeResult->GetCutCascRadius()                &&
        fTreeCascVarCascCosPA                   > lCascadeResult->GetCutCascCosPA()                 &&
        fTreeCascVarLeastNbrCrossedRows         > lCascadeResult->GetCutLeastNumberOfCrossedRows()  &&
        lPDGMass*fTreeCascVarDistOverTotMom     < lCascadeResult->GetCutProperLifetime()            &&
        // - Miscellaneous
        fTreeCascVarLeastNbrClusters            > lCascadeResult->GetCutLeastNumberOfClusters()     &&
        
        //Check 4: TPC dEdx selections
        TMath::Abs(lBachTPCNSigma)              < lCascadeResult->GetCutTPCdEdx()                   &&
        TMath::Abs(lPosTPCNSigma )              < lCascadeResult->GetCutTPCdEdx()                   &&
        TMath::Abs(lNegTPCNSigma )              < lCascadeResult->GetCutTPCdEdx()                   &&
        
        //Check 5: Min/Max V0 Lifetime cut
        lLambdaPDGMass*fTreeCascVarV0DistOverTotMom > lCascadeResult->GetCutMinV0Lifetime()         &&
        lLambdaPDGMass*fTreeCascVarV0DistOverTotMom < lCascadeResult->GetCutMaxV0Lifetime()         &&
        
        //Check 6: Min Track Length 
        fTreeCascVarMinTrackLength              > lCascadeResult->GetCutMinTrackLength()            &&
        
        //Check 7: modern track quality selections
        fTreeCascVarNbrCrossedRowsOverLength    > lCascadeResult->GetCutMinCrossedRowsOverLength()  &&

        //Check 8: Max Chi2/Clusters 
        fTreeCascVarMaxChi2PerCluster           < lCascadeResult->GetCutMaxChi2PerCluster()         &&
        
        //Check 9: Experimental DCA Bachelor to Baryon cut
        fTreeCascVarDCABachToBaryon             > lCascadeResult->GetCutDCABachToBaryon()           &&
        
        //Check 10: Cut on the Bach Baryon pointing angle (not CosPA)
        //TMath::ACos(fTreeCascVarWrongCosPA)  > TMath::ACos(lCascadeResult->GetCutBachBaryonCosPA()) &&
        
        //Check 11a: Xi rejection for Omega analysis
        (   //Is not Omega, so skip this cut
            !lCascResultIsOmega   ||
            //Is Omega --> cut on the Xi mass
            (lCascResultIsOmega   && 
            TMath::Abs( fTreeCascVarMassAsXi - lXiPDGMass ) > lCascadeResult->GetCutXiRejection() ) 
            
        )                                                                                           &&
        
        //Check 11b: Mass cut rejection
        TMath::Abs( lInvMass - lPDGMass )       < 0.075                                             &&       
        
        //Check 12: kITSrefit track selection if requested
        (   //do not need ITS refit tracks --> skip this cut
            !( lCascadeResult->GetCutUseITSRefitTracks() )                  ||
            //otherwise : need to have ITS refit tracks
            (   ( fTreeCascVarBachTrackStatus   & AliAODTrack::kITSrefit )  &&
                ( fTreeCascVarPosTrackStatus    & AliAODTrack::kITSrefit )  &&
                ( fTreeCascVarNegTrackStatus    & AliAODTrack::kITSrefit )  
            )
        )                                                                                           &&
        
        //Check 13: cowboy/sailor for V0
        ( lCascadeResult->GetCutIsCowboy()  ==  0                                            ||
        (lCascadeResult->GetCutIsCowboy()   ==  1 && fTreeCascVarIsCowboy           )        ||
        (lCascadeResult->GetCutIsCowboy()   == -1 && !fTreeCascVarIsCowboy          )      )        &&//end cowboy/sailor
        
        //Check 14: cowboy/sailor for cascade
        ( lCascadeResult->GetCutIsCascadeCowboy()   ==  0                                      ||
        (lCascadeResult->GetCutIsCascadeCowboy()    ==  1 && fTreeCascVarIsCascadeCowboy   )   ||
        (lCascadeResult->GetCutIsCascadeCowboy()    == -1 && !fTreeCascVarIsCascadeCowboy  ) )      &&//end cowboy/sailor
        
        //Check 15: ITS or TOF required 
        (
			!( lCascadeResult->GetCutITSorTOF() ) || 
            (   //One of the daughter track has ITSrefit 
                (   (fTreeCascVarBachTrackStatus  & AliAODTrack::kITSrefit)  ||
                    (fTreeCascVarPosTrackStatus   & AliAODTrack::kITSrefit)  ||
                    (fTreeCascVarNegTrackStatus   & AliAODTrack::kITSrefit) 
                )||
                //Or one of the daughter track has a TOF signal
                (   (TMath::Abs(fTreeCascVarBachTOFExpTDiff+2500. )     > 1e-6 ) ||
                    (TMath::Abs(fTreeCascVarPosTOFExpTDiff+2500.  )     > 1e-6 ) ||
                    (TMath::Abs(fTreeCascVarNegTOFExpTDiff+2500.  )     > 1e-6 ) 
                )
            )
            
        )
        
    )//end major if
    {
        return kTRUE;
    }
    
    return kFALSE;
}
