/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskStrangenessML
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH3D.h"
#include "TList.h"
#include "TFile.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDcascade.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskStrangenessML.h"
#include "AliPIDResponse.h"
#include "AliPID.h"
#include "AliMultSelection.h"
#include "AliCascadeVertexer.h"

class AliAnalysisTaskStrangenessML;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskStrangenessML) // classimp: necessary for root

AliAnalysisTaskStrangenessML::AliAnalysisTaskStrangenessML() : AliAnalysisTaskSE(),
    fESD(0), fOutputList(0), fHist_XiMinus_NNpred(0), fHist_XiMinus_BDTpred(0), fHist_XiMinus_NN_IM(0), fHist_XiMinus_BDT_IM(0),
    fHist_XiMinus_Std_IM(0), fPIDResponse(0), fHistCascadeCounts(0), fXiMinusNN(0), fXiMinusBDT(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskStrangenessML::AliAnalysisTaskStrangenessML(const char* name) : AliAnalysisTaskSE(name),
    fESD(0), fOutputList(0), fHist_XiMinus_NNpred(0), fHist_XiMinus_BDTpred(0) , fHist_XiMinus_NN_IM(0), fHist_XiMinus_BDT_IM(0),
    fHist_XiMinus_Std_IM(0), fPIDResponse(0), fHistCascadeCounts(0), fXiMinusNN(0), fXiMinusBDT(0)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskStrangenessML::~AliAnalysisTaskStrangenessML()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskStrangenessML::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)
    
    //ML Models
    /*
    fXiMinusNN = new AliNeuralNetwork("XiMinusNN");
    fXiMinusNN->LoadModel("alien:///alice/cern.ch/user/g/garciagg/myWorkingDir/NN_Recipes/XiMinus_Recipe.root");
    cout << "Test if load model worked" << endl;
    double Xtest[11] = {1,1,1,1,1,1,1,1,1,1,1};
    cout << fXiMinusNN->Predict(Xtest,11) << endl;
    */
    //Set_NN_Recipe("alien:///alice/cern.ch/user/g/garciagg/myWorkingDir/NN_Recipes/XiMinus_Recipe.root");

    //Set_BDT_Recipe("alien:///alice/cern.ch/user/g/garciagg/myWorkingDir/BDT_Recipes/XiMinus_BDT_Recipe.root");

    //Histograms    
    fHistCascadeCounts = new TH1F("fHistCascadeCounts","fHistCascadeCounts",1,-1,1);
    fOutputList->Add(fHistCascadeCounts);
    
    fHist_XiMinus_NNpred = new TH1D("fHist_XiMinus_NNpred","fHist_XiMinus_NNpred", 100, 0.0, 1.0);
    fOutputList->Add(fHist_XiMinus_NNpred);

    fHist_XiMinus_BDTpred = new TH1D("fHist_XiMinus_BDTpred","fHist_XiMinus_BDTpred", 100, 0.0, 1.0);
    fOutputList->Add(fHist_XiMinus_BDTpred);

/*
    for(int i=0; i<10; i++)//threshold
        for(int j=0; j<10; j++)//pT
            for(int k=0; k<10; k++)//centrality
            {
                fHist_XiMinus_NN_IM[i][j][k] = new TH1F(Form("fHist_XiMinus_NN_IM_t%d,p%d,c%d",i,j,k),Form("fHist_XiMinus_NN_IM_t%d,p%d,c%d",i,j,k), 2000, 0.0, 2.0);
                fOutputList->Add(Form("fHist_XiMinus_NN_IM_t%d,p%d,c%d",i,j,k));
            }
*/

    fHist_XiMinus_NN_IM = new TH3D("fHist_XiMinus_NN_IM","fHist_XiMinus_NN_IM", 2000, 0.0, 2.0, 100, 0.0, 10.0, 100, 0.0, 1.0);
    fOutputList->AddFirst(fHist_XiMinus_NN_IM);

    fHist_XiMinus_BDT_IM = new TH3D("fHist_XiMinus_BDT_IM","fHist_XiMinus_BDT_IM", 2000, 0.0, 2.0, 100, 0.0, 10.0, 100, 0.0, 1.0);
    fOutputList->AddFirst(fHist_XiMinus_BDT_IM);

    fHist_XiMinus_Std_IM = new TH2F("fHist_XiMinus_Std_IM","fHist_XiMinus_Std_IM", 2000, 0.0, 2.0, 100, 0.0, 10.0);
    fOutputList->Add(fHist_XiMinus_Std_IM);

    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output

    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if (man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskStrangenessML::UserExec(Option_t *)
{
    // user exec
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you 
    // have access to the current event. 
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());    // get an event (called fESD) from the input file
                                                        // there's another event format (ESD) which works in a similar wya
                                                        // but is more cpu/memory unfriendly. for now, we'll stick with ESD's
    if(!fESD) return;                                   // if the pointer to the event is empty (getting it failed) skip this event

    AliMultSelection *multSelection = static_cast<AliMultSelection*>(fESD->FindListObject("MultSelection"));
    if(!multSelection) return;

    float centrality = multSelection->GetMultiplicityPercentile("V0M");

    int c=-1;
    if(centrality>=0 && centrality<10) c = 0;
    else if(centrality>=10 && centrality<20)  c = 1;
    else if(centrality>=20 && centrality<30)  c = 2;
    else if(centrality>=30 && centrality<40)  c = 3;
    else if(centrality>=40 && centrality<50)  c = 4;
    else if(centrality>=50 && centrality<60)  c = 5;
    else if(centrality>=60 && centrality<70)  c = 6;
    else if(centrality>=70 && centrality<80)  c = 7;
    else if(centrality>=80 && centrality<90)  c = 8;
    else if(centrality>=90 && centrality<100) c = 0;
    else return;

    //Centrality selection
    if(c!=0) return;//00-10%

    const AliESDVertex* PV = fESD->GetPrimaryVertex();
    Double_t PVxyz[3] = {-100.0,-100.0,-100.0};
    PV->GetXYZ( PVxyz );

    //if(PVxyz[2]<-10.||PVxyz[2]>10.) return; //Event selection cut

    //CASCADE LOOP: Xi
    
    //fESD->ResetCascades();

    
    AliCascadeVertexer lCascVtxer;

    Double_t  fCascadeVertexerSels[8];
    fCascadeVertexerSels[0] =  33.   ;  // max allowed chi2 (same as PDC07)
    fCascadeVertexerSels[1] =   0.05 ;  // min allowed V0 impact parameter                    (PDC07 : 0.05   / LHC09a4 : 0.025 )
    fCascadeVertexerSels[2] =   0.010;  // "window" around the Lambda mass                    (PDC07 : 0.008  / LHC09a4 : 0.010 )
    fCascadeVertexerSels[3] =   0.03 ;  // min allowed bachelor's impact parameter            (PDC07 : 0.035  / LHC09a4 : 0.025 )
    fCascadeVertexerSels[4] =   2.0  ;  // max allowed DCA between the V0 and the bachelor    (PDC07 : 0.1    / LHC09a4 : 0.2   )
    fCascadeVertexerSels[5] =   0.95 ;  // min allowed cosine of the cascade pointing angle   (PDC07 : 0.9985 / LHC09a4 : 0.998 )
    fCascadeVertexerSels[6] =   0.4  ;  // min radius of the fiducial volume                  (PDC07 : 0.9    / LHC09a4 : 0.2   )
    fCascadeVertexerSels[7] = 100.   ;  // max radius of the fiducial volume                  (PDC07 : 100    / LHC09a4 : 100   )

    lCascVtxer.SetDefaultCuts(fCascadeVertexerSels);
    lCascVtxer.SetCuts(fCascadeVertexerSels);

    lCascVtxer.V0sTracks2CascadeVertices(fESD);

    Int_t N_V0s = 0;
    N_V0s = Int_t(fESD->GetNumberOfV0s());

    Int_t N_Casc = 0;
    N_Casc = Int_t(fESD->GetNumberOfCascades());

    for(int g=0; g<N_Casc; g++)
    {

	fHistCascadeCounts->Fill(0);

	//Get Xi candidate
	AliESDcascade *Xi = fESD->GetCascade(g);
	Double_t lV0quality = 0.;
	Xi->ChangeMassHypothesis(lV0quality,3312);

	//Just want Xi-
	if( Xi->Charge() > 0 ) continue;

	//Skip candidate if fail to get it 
	if(!Xi)
	{
		AliWarning("ERROR: Could not retrieve the cascade candidate ...");
		continue;
	}

	//Cascade momentum
	Double_t p_Xi[3] = {-100.0,-100.0,-100.0};
	Xi->GetPxPyPz( p_Xi[0], p_Xi[1], p_Xi[2] );

	//Cascade Position
	Double_t Xi_xyz[3] = {-100.0,-100.0,-100.0};
	Xi->GetXYZcascade( Xi_xyz[0], Xi_xyz[1], Xi_xyz[2]);

	//Get Xi candidate tracks
	AliESDtrack *PosTrackXi = fESD->GetTrack( Xi->GetPindex() ); //Positive
	AliESDtrack *NegTrackXi = fESD->GetTrack( Xi->GetNindex() ); //Negative
	AliESDtrack *BachTrackXi = fESD->GetTrack( Xi->GetBindex() ); //Bachelor

	//Skip candidate if fail to get tracks 
	if(!PosTrackXi || !NegTrackXi || !BachTrackXi )
	{
		AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
		continue;
        }

	//Get track momenta
	Double_t p_PosTrackXi[3] = {-100.0,-100.0,-100.0};
	Xi->GetPPxPyPz( p_PosTrackXi[0], p_PosTrackXi[1], p_PosTrackXi[2] );

	Double_t p_NegTrackXi[3] = {-100.0,-100.0,-100.0};
	Xi->GetNPxPyPz( p_NegTrackXi[0], p_NegTrackXi[1], p_NegTrackXi[2] );

	Double_t p_BachTrackXi[3] = {-100.0,-100.0,-100.0};
	Xi->GetBPxPyPz( p_BachTrackXi[0], p_BachTrackXi[1], p_BachTrackXi[2] );

	//Selection cuts---------------------------------------------------------
	
	//Rapidity
	if(TMath::Abs(Xi->RapXi())>0.5) continue;

	//TPC dE/dx
	Float_t PosNSigmaProton = fPIDResponse->NumberOfSigmasTPC(PosTrackXi, AliPID::kProton);
	Float_t NegNSigmaPion = fPIDResponse->NumberOfSigmasTPC(NegTrackXi, AliPID::kPion);
	Float_t BachNSigmaPion = fPIDResponse->NumberOfSigmasTPC(BachTrackXi, AliPID::kPion);

	if(TMath::Abs(PosNSigmaProton)>4.0 || TMath::Abs(NegNSigmaPion)>4.0 || TMath::Abs(BachNSigmaPion)>4.0) continue;

	//Proper lifetime mL/p
	float m_Xi = 1.31486;
	
	float L_Xi = TMath::Sqrt( (Xi_xyz[0]-PVxyz[0])*(Xi_xyz[0]-PVxyz[0]) + (Xi_xyz[1]-PVxyz[1])*(Xi_xyz[1]-PVxyz[1]) + (Xi_xyz[2]-PVxyz[2])*(Xi_xyz[2]-PVxyz[2]) );

	float p_totXi = TMath::Sqrt( p_Xi[0]*p_Xi[0] + p_Xi[1]*p_Xi[1] + p_Xi[2]*p_Xi[2] );

	if( m_Xi*L_Xi/(p_totXi+1e-13) > 15.0 ) continue;

	//Number of Crossed Rows (NCR)
	float PosNCR = PosTrackXi->GetTPCClusterInfo(2,1);
	float NegNCR = NegTrackXi->GetTPCClusterInfo(2,1);
	float BachNCR = BachTrackXi->GetTPCClusterInfo(2,1);

	if( PosNCR < 80.0 || NegNCR < 80.0 || BachNCR < 80.0) continue;

	//Track Lenght (TL)
	float PosTL, NegTL, BachTL;

	if( PosTrackXi->GetInnerParam() )
		PosTL = PosTrackXi->GetLengthInActiveZone(1,2.0,220.0,fESD->GetMagneticField());

	if( NegTrackXi->GetInnerParam() )
		NegTL = NegTrackXi->GetLengthInActiveZone(1,2.0,220.0,fESD->GetMagneticField());

	if( BachTrackXi->GetInnerParam() )
		BachTL = BachTrackXi->GetLengthInActiveZone(1,2.0,220.0,fESD->GetMagneticField());

	if( PosTL < 90.0 || NegTL < 90.0 || BachTL < 90.0 ) continue;

	// NCR/TL
	if( PosNCR/PosTL < 0.8 || NegNCR/NegTL < 0.8 || BachNCR/BachTL < 0.8 ) continue;

	//Invariant Mass "Bump"
	float cosPA_BB	= GetCosPA( PosTrackXi, BachTrackXi, fESD );

	if( cosPA_BB > 0.999928 ) continue;
	
	//-----------------------------------------------------------------------------------------

	//Topological variables--------------------------------------------------------------------

	Bool_t std_cut_flag = kTRUE;

	double X[11] = {0};

	//1.DCA V0 tracks
	double DCA_V0_tracks = Xi->GetDcaV0Daughters(); 
	X[0] = DCA_V0_tracks;
	if(DCA_V0_tracks>1.0) std_cut_flag=kFALSE;

	//2.DCA Neg track PV
	double DCA_Neg_PV = TMath::Abs( NegTrackXi->GetD(PVxyz[0],PVxyz[1],fESD->GetMagneticField()) );
	X[1] = DCA_Neg_PV;
	if(DCA_Neg_PV<0.2) std_cut_flag=kFALSE;

	//3.DCA Pos track PV
	double DCA_Pos_PV = TMath::Abs( PosTrackXi->GetD(PVxyz[0],PVxyz[1],fESD->GetMagneticField()) );
	X[2] = DCA_Pos_PV;
	if(DCA_Neg_PV<0.2) std_cut_flag=kFALSE;

	//4.V0 transverse radius
	Double_t V0_xyz[3] = {-100.0,-100.0,-100.0};
	Xi->GetXYZ( V0_xyz[0], V0_xyz[1], V0_xyz[2]);

	double r_V0 = TMath::Sqrt( V0_xyz[0]*V0_xyz[0] + V0_xyz[1]*V0_xyz[1] );
	X[3] = r_V0;
	if(r_V0<3.0) std_cut_flag=kFALSE;

	//5.V0 cosPA
	double V0cosPA = Xi->GetV0CosineOfPointingAngle(PVxyz[0],PVxyz[1],PVxyz[2]);
	X[4] = V0cosPA;
	if(V0cosPA<0.95) std_cut_flag=kFALSE;

	//6.DCA V0 Bachelor
	double  DCA_V0_Bach = Xi->GetDcaXiDaughters();
	X[5] = DCA_V0_Bach;
	if(DCA_V0_Bach>1.0) std_cut_flag=kFALSE;

	//7.DCA Bach track PV
	double DCA_Bach_PV = TMath::Abs( BachTrackXi->GetD(PVxyz[0],PVxyz[1],fESD->GetMagneticField()) );
	X[6] = DCA_Bach_PV;
	if(DCA_Bach_PV<0.1) std_cut_flag=kFALSE;

	//8.DCA V0 PV
	double DCA_V0_PV = Xi->GetD(PVxyz[0],PVxyz[1],PVxyz[2]);
	X[7] = DCA_V0_PV;
	if(DCA_V0_PV<0.1) std_cut_flag=kFALSE;

	//9.Xi transverse radius
	double r_Xi = TMath::Sqrt( Xi_xyz[0]*Xi_xyz[0] + Xi_xyz[1]*Xi_xyz[1] );
	X[8] = r_Xi;
	if(r_Xi<1.2) std_cut_flag=kFALSE;

	//10.Xi cosPA
	double XicosPA = Xi->GetCascadeCosineOfPointingAngle(PVxyz[0],PVxyz[1],PVxyz[2]);
	X[9] = XicosPA;
	if(XicosPA<0.95) std_cut_flag=kFALSE;

	//11.V0 invariant mass
	double IM_V0 = Xi->GetEffMass();
	X[10] = IM_V0;
	if( IM_V0<1.115683-0.005 || IM_V0>1.115683+0.005) std_cut_flag=kFALSE;

	//-----------------------------------------------------------------------------------------

	//Fill output histograms

	float IM_Xi = Xi->GetEffMassXi();
	float pT_Xi = TMath::Sqrt( p_Xi[0]*p_Xi[0] + p_Xi[1]*p_Xi[1] );

	if(std_cut_flag) fHist_XiMinus_Std_IM->Fill( IM_Xi, pT_Xi );

	double NNpred = 0;
	double BDTpred = 0;

	if(fXiMinusNN) NNpred = fXiMinusNN->Predict(X,11);
	fHist_XiMinus_NNpred->Fill(NNpred);
	fHist_XiMinus_NN_IM->Fill(IM_Xi, pT_Xi, NNpred);

	if(fXiMinusBDT) BDTpred = fXiMinusBDT->Predict(X,100);
	fHist_XiMinus_BDTpred->Fill(BDTpred);
	fHist_XiMinus_BDT_IM->Fill(IM_Xi, pT_Xi, BDTpred); 

    }                           
    //PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskStrangenessML::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
Float_t AliAnalysisTaskStrangenessML::GetCosPA(AliESDtrack *lPosTrack, AliESDtrack *lNegTrack, AliESDEvent *lEvent)
//Encapsulation of CosPA calculation (warning: considers AliESDtrack clones)
{
    Float_t lCosPA = -1;
    
    //Get Magnetic field and primary vertex
    Double_t b=lEvent->GetMagneticField();
    const AliESDVertex *vtxT3D=lEvent->GetPrimaryVertex();
    Double_t xPrimaryVertex=vtxT3D->GetX();
    Double_t yPrimaryVertex=vtxT3D->GetY();
    Double_t zPrimaryVertex=vtxT3D->GetZ();
    
    //Copy AliExternalParam for handling
    AliExternalTrackParam nt(*lNegTrack), pt(*lPosTrack), *lNegClone=&nt, *lPosClone=&pt;
    
    //Find DCA
    Double_t xn, xp, dca=lNegClone->GetDCA(lPosClone,b,xn,xp);
    
    //Propagate to it
    nt.PropagateTo(xn,b); pt.PropagateTo(xp,b);
    
    //Create V0 object to do propagation
    AliESDv0 vertex(nt,1,pt,2); //Never mind indices, won't use
    
    //Get CosPA
    lCosPA=vertex.GetV0CosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
    
    //Return value
    return lCosPA;
}
//_____________________________________________________________________________
void AliAnalysisTaskStrangenessML::Set_NN_Recipe( TString FileName )
{

  fXiMinusNN = new AliNeuralNetwork("XiMinusNN");

  fXiMinusNN->LoadModel(FileName);

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskStrangenessML::Set_BDT_Recipe( TString FileName )
{

  fXiMinusBDT = new AliBDT("XiMinusBDT");

  fXiMinusBDT->LoadModel(FileName);

  return;
}

