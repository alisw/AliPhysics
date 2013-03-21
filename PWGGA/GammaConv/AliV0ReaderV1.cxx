/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *								          *
 * Authors: Svein Lindal, Daniel Lohner					  *
 * Version 1.0								  *
 *									  *
 * Permission to use, copy, modify and distribute this software and its	  *
 * documentation strictly for non-commercial purposes is hereby granted	  *
 * without fee, provided that the above copyright notice appears in all	  *
 * copies and that both the copyright notice and this permission notice	  *
 * appear in the supporting documentation. The authors make no claims	  *
 * about the suitability of this software for any purpose. It is	  *
 * provided "as is" without express or implied warranty.		  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class reconstructing conversion photons from V0s
//---------------------------------------------
////////////////////////////////////////////////

#include "AliV0ReaderV1.h"
#include "AliKFParticle.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliKFParticle.h"
#include "AliKFConversionPhoton.h"
#include "AliAODConversionPhoton.h"
#include "AliConversionPhotonBase.h"
#include "TVector.h"
#include "AliKFVertex.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "TChain.h"
#include "TFile.h"
#include "AliStack.h"
#include "TString.h"
#include "TObjArray.h"

class iostream;


using namespace std;

ClassImp(AliV0ReaderV1)

//________________________________________________________________________
AliV0ReaderV1::AliV0ReaderV1(const char *name) : AliAnalysisTaskSE(name),
    fConversionCuts(NULL),
    fConversionGammas(NULL),
    fUseImprovedVertex(kTRUE),
    fUseOwnXYZCalculation(kTRUE),
    fUseConstructGamma(kFALSE),
    kUseAODConversionPhoton(kTRUE),
    fCreateAOD(kFALSE),
    fDeltaAODBranchName("GammaConv"),
    fDeltaAODFilename("AliAODGammaConversion.root"),
    fEventIsSelected(kFALSE),
    fPeriodName("")
{
    // Default constructor

    DefineInput(0, TChain::Class());
}

//________________________________________________________________________
AliV0ReaderV1::~AliV0ReaderV1()
{
    // default deconstructor

    if(fConversionGammas){
	fConversionGammas->Delete();// Clear Objects
	delete fConversionGammas;
	fConversionGammas=0x0;
    }
}
/*
//________________________________________________________________________
AliV0ReaderV1::AliV0ReaderV1(AliV0ReaderV1 &original) : AliAnalysisTaskSE(original),
    fConversionCuts(NULL),
    fConversionGammas(NULL),
    fUseImprovedVertex(original.fUseImprovedVertex),
    fUseOwnXYZCalculation(original.fUseOwnXYZCalculation),
    fUseConstructGamma(original.fUseConstructGamma),
    kUseAODConversionPhoton(original.kUseAODConversionPhoton),
    fCreateAOD(original.fCreateAOD),
    fDeltaAODBranchName(original.fDeltaAODBranchName),
    fDeltaAODFilename(original.fDeltaAODFilename),
    fEventIsSelected(original.fEventIsSelected)
{
    // Default constructor

    DefineInput(0, TChain::Class());
}

//____________________________________________________________
AliV0ReaderV1 &AliV0ReaderV1::operator=(const AliV0ReaderV1 &ref){
	//
	// Assignment operator
	// Only copies pointers, object is not the owner of the references
	//
    if(this != &ref){
	AliAnalysisTaskSE::operator=(ref);
	fUseImprovedVertex=ref.fUseImprovedVertex;
	fUseOwnXYZCalculation=ref.fUseOwnXYZCalculation;
	fUseConstructGamma=ref.fUseConstructGamma;
	kUseAODConversionPhoton=ref.kUseAODConversionPhoton;
	fCreateAOD=ref.fCreateAOD;
	fDeltaAODBranchName=ref.fDeltaAODBranchName;
	fDeltaAODFilename=ref.fDeltaAODFilename;
	fEventIsSelected=ref.fEventIsSelected;
    }
    return *this;
}
*/
//________________________________________________________________________
void AliV0ReaderV1::Init()
{
    // Initialize function to be called once before analysis

    if(fConversionCuts==NULL){
	if(fConversionCuts==NULL)AliError("No Cut Selection initialized");
    }

    if(fCreateAOD){kUseAODConversionPhoton=kTRUE;}

    if(fConversionGammas != NULL){
	delete fConversionGammas;
	fConversionGammas=NULL;
    }

    if(fConversionGammas == NULL){
	if(kUseAODConversionPhoton){
	    fConversionGammas = new TClonesArray("AliAODConversionPhoton",100);}
	else{
	    fConversionGammas = new TClonesArray("AliKFConversionPhoton",100);}
    }
    fConversionGammas->Delete();//Reset the TClonesArray

}

//________________________________________________________________________
void AliV0ReaderV1::UserCreateOutputObjects()
{
    // Create AODs

    if(fCreateAOD){
	if(fConversionCuts){
	    fDeltaAODBranchName.Append("_");
	    fDeltaAODBranchName.Append(fConversionCuts->GetCutNumber());
	    fDeltaAODBranchName.Append("_gamma");
	}
	fConversionGammas->SetName(fDeltaAODBranchName.Data());

	AddAODBranch("TClonesArray", &fConversionGammas, fDeltaAODFilename.Data());
	AliAnalysisManager::GetAnalysisManager()->RegisterExtraFile(fDeltaAODFilename.Data());
    }

}

//________________________________________________________________________
void AliV0ReaderV1::UserExec(Option_t *){

   if (fPeriodName.CompareTo("") == 0){
      AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
      if(man) {
         AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
         if (inputHandler){
            TTree* tree = (TTree*) inputHandler->GetTree();
            TFile* file = (TFile*) tree->GetCurrentFile();
            TString fileName(file->GetName());
            TObjArray *arr = fileName.Tokenize("/");
            for (Int_t i = 0; i < arr->GetEntriesFast();i++ ){
               TObjString* testObjString = (TObjString*)arr->At(i);
               if (testObjString->GetString().Contains("LHC")){
                  fPeriodName = testObjString->GetString();
                  i = arr->GetEntriesFast();
               }
            }     
//             cout << fileName.Data() << "\t" <<fPeriodName.Data() << endl;
         }
      }
   }

    // Check if correctly initialized
    if(!fConversionGammas)Init();

    // User Exec
    fEventIsSelected=ProcessEvent(fInputEvent,fMCEvent);
}

//________________________________________________________________________
Bool_t AliV0ReaderV1::ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent)
{
    //Reset the TClonesArray
    fConversionGammas->Delete();

    fInputEvent=inputEvent;
    fMCEvent=mcEvent;

    if(!fInputEvent){
	AliError("No Input event");
	return kFALSE;
    }

    if(!fConversionCuts){AliError("No ConversionCuts");return kFALSE;}

    // Event Cuts
    if(!fConversionCuts->EventIsSelected(fInputEvent,fMCEvent))return kFALSE;

    // Set Magnetic Field
    AliKFParticle::SetField(fInputEvent->GetMagneticField());

    if(fInputEvent->IsA()==AliESDEvent::Class()){
	ProcessESDV0s();
    }
    if(fInputEvent->IsA()==AliAODEvent::Class()){
	GetAODConversionGammas();
    }

    // AOD Output
    FillAODOutput();

    return kTRUE;
}
///________________________________________________________________________
void AliV0ReaderV1::FillAODOutput()
{
    // Fill AOD Output with reconstructed Photons

    if(fInputEvent->IsA()==AliESDEvent::Class()){
	///Make sure delta aod is filled if standard aod is filled (for synchronization when reading aod with standard aod)
	if(fCreateAOD) {
	    AliAODHandler * aodhandler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
            if (aodhandler && aodhandler->GetFillAOD()) {
               AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillExtension(kTRUE);
	      //PostData(0, fConversionGammas);

	    }
	}
    }
}

///________________________________________________________________________
const AliExternalTrackParam *AliV0ReaderV1::GetExternalTrackParam(AliESDv0 *fCurrentV0,Int_t &tracklabel,Int_t charge){

    // Get External Track Parameter with given charge

    if(!(charge==1||charge==-1)){AliError("Charge not defined");return 0x0;}

    // Check for sign flip
    if(fCurrentV0){
	if(!fCurrentV0->GetParamN()||!fCurrentV0->GetParamP())return 0x0;
	if(!fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetNindex())||!fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetPindex()))return 0x0;
	if((fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetPindex()))->Charge()==charge){
	    tracklabel=fCurrentV0->GetPindex();
	    return fCurrentV0->GetParamP();}
	if((fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetNindex()))->Charge()==charge){
	    tracklabel=fCurrentV0->GetNindex();
	    return fCurrentV0->GetParamN();}
    }
    return 0x0;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::ProcessESDV0s()
{
    // Process ESD V0s for conversion photon reconstruction

    AliESDEvent *fESDEvent=dynamic_cast<AliESDEvent*>(fInputEvent);

    AliKFConversionPhoton *fCurrentMotherKFCandidate=NULL;

    if(fESDEvent){

	for(Int_t currentV0Index=0;currentV0Index<fESDEvent->GetNumberOfV0s();currentV0Index++){
	    AliESDv0 *fCurrentV0=(AliESDv0*)(fESDEvent->GetV0(currentV0Index));
	    if(!fCurrentV0){
		printf("Requested V0 does not exist");
		continue;}

	    fCurrentMotherKFCandidate=ReconstructV0(fCurrentV0,currentV0Index);

	    if(fCurrentMotherKFCandidate){

		// Add Gamma to the TClonesArray

		if(kUseAODConversionPhoton){
		    new((*fConversionGammas)[fConversionGammas->GetEntriesFast()]) AliAODConversionPhoton(fCurrentMotherKFCandidate);
		}
		else{
		    new((*fConversionGammas)[fConversionGammas->GetEntriesFast()]) AliKFConversionPhoton(*fCurrentMotherKFCandidate);
		}

		delete fCurrentMotherKFCandidate;
		fCurrentMotherKFCandidate=NULL;
	    }
	}
    }
    return kTRUE;
}

///________________________________________________________________________
AliKFConversionPhoton *AliV0ReaderV1::ReconstructV0(AliESDv0 *fCurrentV0,Int_t currentV0Index)
{
	// Reconstruct conversion photon from ESD v0
    fConversionCuts->FillPhotonCutIndex(AliConversionCuts::kPhotonIn);

    //checks if on the fly mode is set
    if(!fConversionCuts->SelectV0Finder(fCurrentV0->GetOnFlyStatus())){
       fConversionCuts->FillPhotonCutIndex(AliConversionCuts::kOnFly);
       return 0x0;
    }

    // TrackLabels
    Int_t currentTrackLabels[2]={-1,-1};

    // Get Daughter KF Particles

    const AliExternalTrackParam *fCurrentExternalTrackParamPositive=GetExternalTrackParamP(fCurrentV0,currentTrackLabels[0]);
    const AliExternalTrackParam *fCurrentExternalTrackParamNegative=GetExternalTrackParamN(fCurrentV0,currentTrackLabels[1]);

    if(!fCurrentExternalTrackParamPositive||!fCurrentExternalTrackParamNegative)return 0x0;

    // Apply some Cuts before Reconstruction

    AliVTrack * posTrack = fConversionCuts->GetTrack(fInputEvent,currentTrackLabels[0]);
    AliVTrack * negTrack = fConversionCuts->GetTrack(fInputEvent,currentTrackLabels[1]);

    if(!negTrack || !posTrack) {
       fConversionCuts->FillPhotonCutIndex(AliConversionCuts::kNoTracks);
       return 0x0;
    }

    // Track Cuts
    if(!fConversionCuts->TracksAreSelected(negTrack, posTrack)){
       fConversionCuts->FillPhotonCutIndex(AliConversionCuts::kTrackCuts);
       return 0x0;
    }

    // PID Cuts
    if(!fConversionCuts->dEdxCuts(negTrack) || !fConversionCuts->dEdxCuts(posTrack)) {
       fConversionCuts->FillPhotonCutIndex(AliConversionCuts::kdEdxCuts);
       return 0x0;
    }

    // Reconstruct Photon

    AliKFConversionPhoton *fCurrentMotherKF=NULL;

    AliKFParticle fCurrentNegativeKFParticle(*(fCurrentExternalTrackParamNegative),11);
    AliKFParticle fCurrentPositiveKFParticle(*(fCurrentExternalTrackParamPositive),-11);

    // Reconstruct Gamma

    if(fUseConstructGamma){
	fCurrentMotherKF = new AliKFConversionPhoton();
	fCurrentMotherKF->ConstructGamma(fCurrentNegativeKFParticle,fCurrentPositiveKFParticle);
    }else{
	fCurrentMotherKF = new AliKFConversionPhoton(fCurrentNegativeKFParticle,fCurrentPositiveKFParticle);
	fCurrentMotherKF->SetMassConstraint(0,0);
    }

    // Set Track Labels

    fCurrentMotherKF->SetTrackLabels(currentTrackLabels[0],currentTrackLabels[1]);

    // Set V0 index

    fCurrentMotherKF->SetV0Index(currentV0Index);

    //Set MC Label

    if(fMCEvent){

	AliStack *fMCStack= fMCEvent->Stack();

	Int_t labelp=TMath::Abs(fConversionCuts->GetTrack(fInputEvent,fCurrentMotherKF->GetTrackLabelPositive())->GetLabel());
	Int_t labeln=TMath::Abs(fConversionCuts->GetTrack(fInputEvent,fCurrentMotherKF->GetTrackLabelNegative())->GetLabel());

	TParticle *fNegativeMCParticle = fMCStack->Particle(labeln);
	TParticle *fPositiveMCParticle = fMCStack->Particle(labelp);

	if(fPositiveMCParticle&&fNegativeMCParticle){
	    fCurrentMotherKF->SetMCLabelPositive(labelp);
	    fCurrentMotherKF->SetMCLabelNegative(labeln);
	}
    }

    // Update Vertex (moved for same eta compared to old)
    if(fUseImprovedVertex == kTRUE){
       AliKFVertex primaryVertexImproved(*fInputEvent->GetPrimaryVertex());
       primaryVertexImproved+=*fCurrentMotherKF;
       fCurrentMotherKF->SetProductionVertex(primaryVertexImproved);
    }

    // SetPsiPair

    Double_t PsiPair=GetPsiPair(fCurrentV0,fCurrentExternalTrackParamPositive,fCurrentExternalTrackParamNegative);
    fCurrentMotherKF->SetPsiPair(PsiPair);


    // Recalculate ConversionPoint
    Double_t dca[2]={0,0};
    if(fUseOwnXYZCalculation){
	Double_t convpos[3]={0,0,0};
	if(!GetConversionPoint(fCurrentExternalTrackParamPositive,fCurrentExternalTrackParamNegative,convpos,dca)){
           fConversionCuts->FillPhotonCutIndex(AliConversionCuts::kConvPointFail);
           delete fCurrentMotherKF;
           fCurrentMotherKF=NULL;
           return 0x0;
	}

	fCurrentMotherKF->SetConversionPoint(convpos);
    }

    if(fCurrentMotherKF->GetNDF() > 0.)
	fCurrentMotherKF->SetChi2perNDF(fCurrentMotherKF->GetChi2()/fCurrentMotherKF->GetNDF());   //->Photon is created before all chi2 relevant changes are performed, set it "by hand"


    // Set Dilepton Mass (moved down for same eta compared to old)
    fCurrentMotherKF->SetMass(fCurrentMotherKF->M());

    // Apply Photon Cuts

    if(!fConversionCuts->PhotonCuts(fCurrentMotherKF,fInputEvent)){
	fConversionCuts->FillPhotonCutIndex(AliConversionCuts::kPhotonCuts);
	delete fCurrentMotherKF;
	fCurrentMotherKF=NULL;
	return 0x0;
    }

    fConversionCuts->FillPhotonCutIndex(AliConversionCuts::kPhotonOut);
    return fCurrentMotherKF;
}

///________________________________________________________________________
Double_t AliV0ReaderV1::GetPsiPair(const AliESDv0* v0, const AliExternalTrackParam *positiveparam,const AliExternalTrackParam *negativeparam) const {
    //
    // Angle between daughter momentum plane and plane
    //

   AliExternalTrackParam nt(*negativeparam);
   AliExternalTrackParam pt(*positiveparam);

   Float_t magField = fInputEvent->GetMagneticField();

   Double_t xyz[3] = {0.,0.,0.};
   v0->GetXYZ(xyz[0],xyz[1],xyz[2]);

   // Double_t pPlus[3]  = {pt.Px(),pt.Py(),pt.Pz()};
   // Double_t pMinus[3] = {nt.Px(),nt.Py(),nt.Pz()};

   // Double_t u[3] = {pPlus[0]+pMinus[0],pPlus[1]+pMinus[1],pPlus[2]+pMinus[2]};
   // Double_t normu = sqrt( (u[0]*u[0]) + (u[1]*u[1]) + (u[2]*u[2]) );
   
   // u[0] = u[0] / normu;
   // u[1] = u[1] / normu;
   // u[2] = u[2] / normu;

   // Double_t normpPlus = sqrt( (pPlus[0]*pPlus[0]) + (pPlus[1]*pPlus[1]) + (pPlus[2]*pPlus[2]) );
   // Double_t normpMinus = sqrt( (pMinus[0]*pMinus[0]) + (pMinus[1]*pMinus[1]) + (pMinus[2]*pMinus[2]) );

   // pPlus[0] = pPlus[0] / normpPlus;
   // pPlus[1] = pPlus[1] / normpPlus;
   // pPlus[2] = pPlus[2] / normpPlus;

   // pMinus[0] = pMinus[0] / normpMinus;
   // pMinus[1] = pMinus[1] / normpMinus;
   // pMinus[2] = pMinus[2] / normpMinus;

   // Double_t v[3] = {0,0,0}; // pPlus X pMinus
   // v[0] = (pPlus[1]*pMinus[2]) - (pPlus[2]*pMinus[1]);
   // v[1] = (pPlus[2]*pMinus[0]) - (pPlus[0]*pMinus[2]);
   // v[2] = (pPlus[0]*pMinus[1]) - (pPlus[1]*pMinus[0]);
   
   // Double_t w[3] = {0,0,0}; // u X v
   // w[0] = (u[1]*v[2]) - (u[2]*v[1]);
   // w[1] = (u[2]*v[0]) - (u[0]*v[2]);
   // w[2] = (u[0]*v[1]) - (u[1]*v[0]);

   // Double_t z[3] = {0,0,1};
   // Double_t wc[3] = {0,0,0}; // u X v
   // wc[0] = (u[1]*z[2]) - (u[2]*z[1]);
   // wc[1] = (u[2]*z[0]) - (u[0]*z[2]);
   // wc[2] = (u[0]*z[1]) - (u[1]*z[0]);

   // Double_t PhiV = TMath::ACos((w[0]*wc[0]) + (w[1]*wc[1]) + (w[2]*wc[2]));
   //return abs(PhiV);


   // TVector3 pPlus(pt.Px(),pt.Py(),pt.Pz());
   // TVector3 pMinus(nt.Px(),nt.Py(),nt.Pz());

   // TVector3 u = pMinus + pPlus;
   // u = u*(1/u.Mag());

   // TVector3 pHPlus = pPlus*(1/pPlus.Mag());
   // TVector3 pHMinus = pMinus*(1/pMinus.Mag());

   // TVector3 v = pHPlus.Cross(pHMinus);
   // TVector3 w = u.Cross(v);
   // TVector3 z(0,0,1);
   // TVector3 wc = u.Cross(z);

   // Double_t PhiV = w * wc;

   Double_t mn[3] = {0,0,0};
   Double_t mp[3] = {0,0,0};

   v0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
   v0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter;

   Double_t deltat = 1.;
   deltat = TMath::ATan(mp[2]/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1])+1.e-13)) - TMath::ATan(mn[2]/(TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1])+1.e-13));//difference of angles of the two daughter tracks with z-axis
   Double_t radiussum = TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]) + 50;//radius to which tracks shall be propagated

   Double_t momPosProp[3] = {0,0,0};
   Double_t momNegProp[3] = {0,0,0};

   Double_t psiPair = 4.;
   if(nt.PropagateTo(radiussum,magField) == 0) return psiPair; //propagate tracks to the outside -> Better Purity and Efficiency
   if(pt.PropagateTo(radiussum,magField) == 0) return psiPair; //propagate tracks to the outside -> Better Purity and Efficiency

   pt.GetPxPyPz(momPosProp);//Get momentum vectors of tracks after propagation
   nt.GetPxPyPz(momNegProp);

   Double_t pEle =
       TMath::Sqrt(momNegProp[0]*momNegProp[0]+momNegProp[1]*momNegProp[1]+momNegProp[2]*momNegProp[2]);//absolute momentum value of negative daughter

   Double_t pPos =
       TMath::Sqrt(momPosProp[0]*momPosProp[0]+momPosProp[1]*momPosProp[1]+momPosProp[2]*momPosProp[2]);//absolute momentum value of positive daughter

   Double_t scalarproduct =
       momPosProp[0]*momNegProp[0]+momPosProp[1]*momNegProp[1]+momPosProp[2]*momNegProp[2];//scalar product of propagated positive and negative daughters' momenta

   Double_t chipair = TMath::ACos(scalarproduct/(pEle*pPos));//Angle between propagated daughter tracks

   psiPair =  TMath::Abs(TMath::ASin(deltat/chipair));

   return psiPair;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::GetHelixCenter(const AliExternalTrackParam *track,Double_t center[2]){

    // Get Center of the helix track parametrization

    Int_t charge=track->Charge();
    Double_t b=fInputEvent->GetMagneticField();

    Double_t	helix[6];
    track->GetHelixParameters(helix,b);

    Double_t xpos =	helix[5];
    Double_t ypos =	helix[0];
    Double_t radius = TMath::Abs(1./helix[4]);
    Double_t phi = helix[2];

    if(phi < 0){
	phi = phi + 2*TMath::Pi();
    }

    phi -= TMath::Pi()/2.;
    Double_t xpoint =	radius * TMath::Cos(phi);
    Double_t ypoint =	radius * TMath::Sin(phi);

    if(b<0){
	if(charge > 0){
	    xpoint = - xpoint;
	    ypoint = - ypoint;
	}

	if(charge < 0){
	    xpoint =	xpoint;
	    ypoint =	ypoint;
	}
    }
    if(b>0){
	if(charge > 0){
	    xpoint =	xpoint;
	    ypoint =	ypoint;
	}

	if(charge < 0){
	    xpoint = - xpoint;
	    ypoint = - ypoint;
	}
    }
    center[0] =	xpos + xpoint;
    center[1] =	ypos + ypoint;

    return 1;
}
///________________________________________________________________________
Bool_t AliV0ReaderV1::GetConversionPoint(const AliExternalTrackParam *pparam,const AliExternalTrackParam *nparam,Double_t convpos[3],Double_t dca[2]){

    // Recalculate Conversion Point

    if(!pparam||!nparam)return kFALSE;

    Double_t helixcenterpos[2];
    GetHelixCenter(pparam,helixcenterpos);

    Double_t helixcenterneg[2];
    GetHelixCenter(nparam,helixcenterneg);

	Double_t helixpos[6];
	pparam->GetHelixParameters(helixpos,fInputEvent->GetMagneticField());
	Double_t posradius = TMath::Abs(1./helixpos[4]);

	Double_t helixneg[6];
	nparam->GetHelixParameters(helixneg,fInputEvent->GetMagneticField());
	Double_t negradius = TMath::Abs(1./helixneg[4]);

        // Calculate xy-position

	Double_t xpos = helixcenterpos[0];
	Double_t ypos = helixcenterpos[1];
	Double_t xneg = helixcenterneg[0];
	Double_t yneg = helixcenterneg[1];

	convpos[0] = (xpos*negradius + xneg*posradius)/(negradius+posradius);
	convpos[1] = (ypos*negradius+ yneg*posradius)/(negradius+posradius);


	// Calculate Track XY vertex position

	Double_t deltaXPos = convpos[0] -	xpos;
	Double_t deltaYPos = convpos[1] -	ypos;

	Double_t deltaXNeg = convpos[0] -	xneg;
	Double_t deltaYNeg = convpos[1] -	yneg;

	Double_t alphaPos =	TMath::Pi() + TMath::ATan2(-deltaYPos,-deltaXPos);
	Double_t alphaNeg =	TMath::Pi() + TMath::ATan2(-deltaYNeg,-deltaXNeg);

	Double_t vertexXNeg =	xneg +	TMath::Abs(negradius)*TMath::Cos(alphaNeg);
	Double_t vertexYNeg =	yneg +	TMath::Abs(negradius)*TMath::Sin(alphaNeg);

	Double_t vertexXPos =	xpos +	TMath::Abs(posradius)*TMath::Cos(alphaPos);
	Double_t vertexYPos =	ypos +	TMath::Abs(posradius)*TMath::Sin(alphaPos);

	AliExternalTrackParam p(*pparam);
        AliExternalTrackParam n(*nparam);

	TVector2 vertexPos(vertexXPos,vertexYPos);
	TVector2 vertexNeg(vertexXNeg,vertexYNeg);

	// Convert to local coordinate system
	TVector2 vertexPosRot=vertexPos.Rotate(-p.GetAlpha());
	TVector2 vertexNegRot=vertexNeg.Rotate(-n.GetAlpha());

	// Propagate Track Params to Vertex

	if(!p.PropagateTo(vertexPosRot.X(),fInputEvent->GetMagneticField()))return kFALSE;
	if(!n.PropagateTo(vertexNegRot.X(),fInputEvent->GetMagneticField()))return kFALSE;

	// Check whether propagation was sucessful

	if(TMath::Abs(vertexPos.Mod()-TMath::Sqrt(p.GetX()*p.GetX()+p.GetY()*p.GetY()))>0.01)return kFALSE;
	if(TMath::Abs(vertexNeg.Mod()-TMath::Sqrt(n.GetX()*n.GetX()+n.GetY()*n.GetY()))>0.01)return kFALSE;

	// Calculate z position

	convpos[2] = (p.GetZ()*negradius+n.GetZ()*posradius)/(negradius+posradius);

	// Calculate DCA
	TVector2 vdca=vertexPos-vertexNeg;
	dca[0]=vdca.Mod();
	dca[1]=TMath::Abs(n.GetZ()-p.GetZ());

	return kTRUE;
}
//________________________________________________________________________
Bool_t AliV0ReaderV1::GetAODConversionGammas(){

    // Get reconstructed conversion photons from satellite AOD file

    AliAODEvent *fAODEvent=dynamic_cast<AliAODEvent*>(fInputEvent);

    if(fAODEvent){

       if(fConversionGammas == NULL){
          fConversionGammas = new TClonesArray("AliAODConversionPhoton",100);
       }
	fConversionGammas->Delete();//Reset the TClonesArray

	//Get Gammas from satellite AOD gamma branch

	AliAODConversionPhoton *gamma=0x0;

	TClonesArray *fInputGammas=dynamic_cast<TClonesArray*>(fAODEvent->FindListObject(fDeltaAODBranchName.Data()));
	if(!fInputGammas){
           FindDeltaAODBranchName();
           fInputGammas=dynamic_cast<TClonesArray*>(fAODEvent->FindListObject(fDeltaAODBranchName.Data()));}
	if(!fInputGammas){AliError("No Gamma Satellites found");return kFALSE;}
	// Apply Selection Cuts to Gammas and create local working copy
	if(fInputGammas){
           for(Int_t i=0;i<fInputGammas->GetEntriesFast();i++){
              gamma=dynamic_cast<AliAODConversionPhoton*>(fInputGammas->At(i));
              if(gamma){
                 if(fConversionCuts->PhotonIsSelected(gamma,fInputEvent)){
                    new((*fConversionGammas)[fConversionGammas->GetEntriesFast()]) AliAODConversionPhoton(*gamma);}
              }
           }
	}
    }

    if(fConversionGammas->GetEntries()){return kTRUE;}

    return kFALSE;
}

//________________________________________________________________________
void AliV0ReaderV1::FindDeltaAODBranchName(){

    // Find delta AOD branchname containing reconstructed photons

    TList *list=fInputEvent->GetList();
    for(Int_t ii=0;ii<list->GetEntries();ii++){
	TString name((list->At(ii))->GetName());
	if(name.BeginsWith(fDeltaAODBranchName)&&name.EndsWith("gamma")){
	    fDeltaAODBranchName=name;
	    AliInfo(Form("Set DeltaAOD BranchName to: %s",fDeltaAODBranchName.Data()));
	    return;
	}
    }
}


//________________________________________________________________________
void AliV0ReaderV1::Terminate(Option_t *)
{

}
