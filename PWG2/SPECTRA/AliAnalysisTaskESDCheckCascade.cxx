
/**************************************************************************
 *  Authors : Antonin Maire, Boris Hippolyte                              *
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

//-----------------------------------------------------------------
//                 AliAnalysisTaskESDCheckCascade class
//            This task is for QAing the Cascades from the ESD
//            Origin : Antonin Maire Fev2008, antonin.maire@ires.in2p3.fr
//	      Modified : A.M Juillet 2008
//-----------------------------------------------------------------


#include <iostream>

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TMath.h"


#include "AliLog.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliESDv0.h"
#include "AliESDcascade.h"
#include "AliESDVertex.h"

#include "AliAnalysisTaskESDCheckCascade.h"

ClassImp(AliAnalysisTaskESDCheckCascade)



//________________________________________________________________________
AliAnalysisTaskESDCheckCascade::AliAnalysisTaskESDCheckCascade() 
  : AliAnalysisTask(), 

    fesdH(0), fMCevent(0), fESD(0),

	// - Cascade part initialisation
    fListHist_Cascade(0),
    fHistTrackMultiplicity(0), fHistCascadeMultiplicity(0),
    fHistPos_TrkgPrimaryVtx_x(0), fHistPos_TrkgPrimaryVtx_y(0), fHistPos_TrkgPrimaryVtx_z(0), fHist_TrkgPrimaryVtxRadius(0),
    fHistPos_SPDPrimaryVtx_x(0), fHistPos_SPDPrimaryVtx_y(0), fHistPos_SPDPrimaryVtx_z(0), fHist_SPDPrimaryVtxRadius(0),
    f2dHistTrkgPrimVtxVsSPDPrimVtx(0),

    fHistEffMassXi(0),  fHistChi2Xi(0),  
    fHistDcaXiDaughters(0), fHistDcaBachToPrimVertex(0), fHistXiCosineOfPointingAngle(0), fHistXiRadius(0),

    fHistXiMom_transv(0),    fHistXiMom_tot(0),
    fHistBachMom_transv(0),   fHistBachMom_tot(0),


    fHistMassLambdaAsCascDghter(0),
    fHistV0Chi2_Xi(0),
    fHistDcaV0Daughters_Xi(0),
    fHistDcaV0ToPrimVertex_Xi(0), 
    fHistV0CosineOfPointingAngle_Xi(0),
    fHistV0Radius_Xi(0),
    fHistDcaPosToPrimVertex_Xi(0), fHistDcaNegToPrimVertex_Xi(0), 

    fHistMassXiMinus(0), fHistMassXiPlus(0),
    fHistMassOmegaMinus(0), fHistMassOmegaPlus(0)

{
  // Dummy Constructor
}








//________________________________________________________________________
AliAnalysisTaskESDCheckCascade::AliAnalysisTaskESDCheckCascade(const char *name) 
  : AliAnalysisTask(name, ""), 

    fesdH(0), fMCevent(0), fESD(0),    
    
    	// - Cascade part initialisation
    fListHist_Cascade(0),
    fHistTrackMultiplicity(0), fHistCascadeMultiplicity(0),
    fHistPos_TrkgPrimaryVtx_x(0), fHistPos_TrkgPrimaryVtx_y(0), fHistPos_TrkgPrimaryVtx_z(0), fHist_TrkgPrimaryVtxRadius(0),
    fHistPos_SPDPrimaryVtx_x(0), fHistPos_SPDPrimaryVtx_y(0), fHistPos_SPDPrimaryVtx_z(0), fHist_SPDPrimaryVtxRadius(0),
    f2dHistTrkgPrimVtxVsSPDPrimVtx(0),

    fHistEffMassXi(0),  fHistChi2Xi(0),  
    fHistDcaXiDaughters(0), fHistDcaBachToPrimVertex(0), fHistXiCosineOfPointingAngle(0), fHistXiRadius(0),

    fHistXiMom_transv(0),    fHistXiMom_tot(0),
    fHistBachMom_transv(0),   fHistBachMom_tot(0),


    fHistMassLambdaAsCascDghter(0),
    fHistV0Chi2_Xi(0),
    fHistDcaV0Daughters_Xi(0),
    fHistDcaV0ToPrimVertex_Xi(0), 
    fHistV0CosineOfPointingAngle_Xi(0),
    fHistV0Radius_Xi(0),
    fHistDcaPosToPrimVertex_Xi(0), fHistDcaNegToPrimVertex_Xi(0), 

    fHistMassXiMinus(0), fHistMassXiPlus(0),
    fHistMassOmegaMinus(0), fHistMassOmegaPlus(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container (Cascade)
  DefineOutput(0, TList::Class());
}







//________________________________________________________________________
void AliAnalysisTaskESDCheckCascade::ConnectInputData(Option_t *) 
{

   AliLog::SetGlobalLogLevel(AliLog::kError); // to suppress the extensive info prompted by a run with MC

  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
   	
	fesdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    //AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

	if (!fesdH) {
		Printf("ERROR: Could not get ESDInputHandler");
	} else
		fESD = fesdH->GetEvent();
  }
}







//________________________________________________________________________
void AliAnalysisTaskESDCheckCascade::CreateOutputObjects()
{
  // Create histograms
  // Called once



 fListHist_Cascade = new TList();

  
	// - General histos
  
if(! fHistTrackMultiplicity) {
	fHistTrackMultiplicity = new TH1F("fHistTrackMultiplicity", "Track Multiplicity;Number of tracks;Events", 150, 0, 150);
    	//  fHistTrackMultiplicity = new TH1F("fHistTrackMultiplicity", "Multiplicity distribution;Number of tracks;Events", 200, 0, 40000); //HERE for Pb-Pb
	fListHist_Cascade->Add(fHistTrackMultiplicity);
}

if(! fHistCascadeMultiplicity) {
	fHistCascadeMultiplicity = new TH1F("fHistCascadeMultiplicity", "Cascades per event;Number of Cascades;Events", 10, 0, 10);
	//  fHistCascadeMultiplicity = new TH1F("fHistCascadeMultiplicity", "Multiplicity distribution;Number of Cascades;Events", 200, 0, 4000); //HERE
	fListHist_Cascade->Add(fHistCascadeMultiplicity);
}





	// - Vertex Positions
  
if(! fHistPos_TrkgPrimaryVtx_x ){
	fHistPos_TrkgPrimaryVtx_x   = new TH1F( "fHistPos_TrkgPrimaryVtx_x" , "Prim. Vertex Position in x; x (cm); Events" , 100, -0.5, 0.5 );
	fListHist_Cascade->Add(fHistPos_TrkgPrimaryVtx_x);
}


if(! fHistPos_TrkgPrimaryVtx_y){
	fHistPos_TrkgPrimaryVtx_y   = new TH1F( "fHistPos_TrkgPrimaryVtx_y" , "Prim. Vertex Position in y; y (cm); Events" , 100, -0.5, 0.5 );
	fListHist_Cascade->Add(fHistPos_TrkgPrimaryVtx_y);
}

if(! fHistPos_TrkgPrimaryVtx_z ){
	fHistPos_TrkgPrimaryVtx_z   = new TH1F( "fHistPos_TrkgPrimaryVtx_z" , "Prim. Vertex Position in z; z (cm); Events" , 100, -15.0, 15.0 );
	fListHist_Cascade->Add(fHistPos_TrkgPrimaryVtx_z);
}

if(! fHist_TrkgPrimaryVtxRadius ){
	fHist_TrkgPrimaryVtxRadius  = new TH1F( "fHist_TrkgPrimaryVtxRadius",  "Prim.  vertex radius; r (cm); Events" , 100, 0., 15.0 );
	fListHist_Cascade->Add(fHist_TrkgPrimaryVtxRadius);
}




if(! fHistPos_SPDPrimaryVtx_x ){
	fHistPos_SPDPrimaryVtx_x   = new TH1F( "fHistPos_SPDPrimaryVtx_x" , "SPD Prim. Vertex Position in x; x (cm); Events" , 100, -0.5, 0.5 );
	fListHist_Cascade->Add(fHistPos_SPDPrimaryVtx_x);
}

if(! fHistPos_SPDPrimaryVtx_y){
	fHistPos_SPDPrimaryVtx_y   = new TH1F( "fHistPos_SPDPrimaryVtx_y" , "SPD Prim. Vertex Position in y; y (cm); Events" , 100, -0.5, 0.5 );
	fListHist_Cascade->Add(fHistPos_SPDPrimaryVtx_y);
}

if(! fHistPos_SPDPrimaryVtx_z ){
	fHistPos_SPDPrimaryVtx_z   = new TH1F( "fHistPos_SPDPrimaryVtx_z" , "SPD Prim. Vertex Position in z; z (cm); Events" , 100, -15.0, 15.0 );
	fListHist_Cascade->Add(fHistPos_SPDPrimaryVtx_z);
}

if(! fHist_SPDPrimaryVtxRadius ){
	fHist_SPDPrimaryVtxRadius  = new TH1F( "fHist_SPDPrimaryVtxRadius",  "SPD Prim.  vertex radius; r (cm); Events" , 100, 0., 15.0 );
	fListHist_Cascade->Add(fHist_SPDPrimaryVtxRadius);
}

if(! f2dHistTrkgPrimVtxVsSPDPrimVtx) {
	f2dHistTrkgPrimVtxVsSPDPrimVtx = new TH2F( "f2dHistTrkgPrimVtxVsSPDPrimVtx", "r_{Trck Prim. Vtx} Vs r_{SPD Prim. Vtx}; r_{Track Vtx} (cm); r_{SPD Vtx} (cm)", 100, 0., 15.0, 100, 0., 15.);
	fListHist_Cascade->Add(f2dHistTrkgPrimVtxVsSPDPrimVtx);
}




// - Typical histos for cascades


if(! fHistEffMassXi) {
     fHistEffMassXi = new TH1F("fHistEffMassXi", "Cascade candidates ; Invariant Mass (GeV/c^{2}) ; Counts", 200, 1.2, 2.0);
     fListHist_Cascade->Add(fHistEffMassXi);
}
   
if(! fHistChi2Xi ){
	fHistChi2Xi = new TH1F("fHistChi2Xi", "Cascade #chi^{2}; #chi^{2}; Number of Cascades", 160, 0, 40);
	fListHist_Cascade->Add(fHistChi2Xi);
}
  
if(! fHistDcaXiDaughters ){
	fHistDcaXiDaughters = new TH1F( "fHistDcaXiDaughters",  "DCA between Xi Daughters; DCA (cm) ; Number of Cascades", 100, 0., 0.5);
	fListHist_Cascade->Add(fHistDcaXiDaughters);
}

if(! fHistDcaBachToPrimVertex) {
	fHistDcaBachToPrimVertex = new TH1F("fHistDcaBachToPrimVertex", "DCA of Bach. to Prim. Vertex;DCA (cm);Number of Cascades", 100, 0., 0.25);
	fListHist_Cascade->Add(fHistDcaBachToPrimVertex);
}

if(! fHistXiCosineOfPointingAngle) {
	fHistXiCosineOfPointingAngle = new TH1F("fHistXiCosineOfPointingAngle", "Cosine of Xi Pointing Angle; Cos (Xi Point.Angl);Number of Xis", 200, 0.98, 1.0);
	fListHist_Cascade->Add(fHistXiCosineOfPointingAngle);
}

if(! fHistXiRadius ){
	fHistXiRadius  = new TH1F( "fHistXiRadius",  "Casc. decay transv. radius; r (cm); Counts" , 200, 0., 20.0 );
	fListHist_Cascade->Add(fHistXiRadius);
}



if(! fHistXiMom_transv ){
	fHistXiMom_transv  = new TH1F( "fHistXiMom_transv" , "Xi transverse momentum ; p_{t}(#Xi) (GeV/c); Counts", 80, 0.0, 8.0);
	fListHist_Cascade->Add(fHistXiMom_transv);
}

if(! fHistXiMom_tot ){
	fHistXiMom_tot  = new TH1F( "fHistXiMom_tot" , "Xi momentum norm; p_{tot}(#Xi) (GeV/c); Counts", 80, 0.0, 8.0);
	fListHist_Cascade->Add(fHistXiMom_tot);
}



if(! fHistBachMom_transv ){
	fHistBachMom_transv  = new TH1F( "fHistBachMom_transv" , "Bach. transverse momentum ; p_{t}(Bach.) (GeV/c); Counts", 50, 0.0, 2.5);
	fListHist_Cascade->Add(fHistBachMom_transv);
}

if(! fHistBachMom_tot ){
	fHistBachMom_tot  = new TH1F( "fHistBachMom_tot" , "Bach. momentum norm; p_{tot}(Bach.) (GeV/c); Counts", 60, 0.0, 3.0);
	fListHist_Cascade->Add(fHistBachMom_tot);
}




// - Histos about ~ the "V0 part" of the cascade,  coming by inheritance from AliESDv0



if (! fHistMassLambdaAsCascDghter) {
     fHistMassLambdaAsCascDghter = new TH1F("fHistMassLambdaAsCascDghter","#Lambda associated to Casc. candidates;Eff. Mass (GeV/c^{2});Counts", 160,1.00,1.8);
    fListHist_Cascade->Add(fHistMassLambdaAsCascDghter);
}

if (! fHistV0Chi2_Xi) {
	fHistV0Chi2_Xi = new TH1F("fHistV0Chi2_Xi", "V0 #chi^{2}, in cascade; #chi^{2};Counts", 160, 0, 40);
	fListHist_Cascade->Add(fHistV0Chi2_Xi);
}

if (! fHistDcaV0Daughters_Xi) {
	fHistDcaV0Daughters_Xi = new TH1F("fHistDcaV0Daughters_Xi", "DCA between V0 daughters, in cascade;DCA (cm);Number of V0s", 120, 0., 0.6);
	fListHist_Cascade->Add(fHistDcaV0Daughters_Xi);
}

if (! fHistDcaV0ToPrimVertex_Xi) {
	fHistDcaV0ToPrimVertex_Xi = new TH1F("fHistDcaV0ToPrimVertex_Xi", "DCA of V0 to Prim. Vertex, in cascade;DCA (cm);Number of Cascades", 200, 0., 1.);
	fListHist_Cascade->Add(fHistDcaV0ToPrimVertex_Xi);
}

if (! fHistV0CosineOfPointingAngle_Xi) {
	fHistV0CosineOfPointingAngle_Xi = new TH1F("fHistV0CosineOfPointingAngle_Xi", "Cosine of V0 Pointing Angle, in cascade;Cos(V0 Point. Angl); Counts", 200, 0.98, 1.0);
	fListHist_Cascade->Add(fHistV0CosineOfPointingAngle_Xi);
}

if (! fHistV0Radius_Xi) {
	fHistV0Radius_Xi = new TH1F("fHistV0Radius_Xi", "V0 decay radius, in cascade; radius (cm); Counts", 200, 0, 20);
	fListHist_Cascade->Add(fHistV0Radius_Xi);
}

if (! fHistDcaPosToPrimVertex_Xi) {
	fHistDcaPosToPrimVertex_Xi = new TH1F("fHistDcaPosToPrimVertex_Xi", "DCA of V0 pos daughter to Prim. Vertex;DCA (cm);Counts", 300, 0, 3);
	fListHist_Cascade->Add(fHistDcaPosToPrimVertex_Xi);
}

if (! fHistDcaNegToPrimVertex_Xi) {
	fHistDcaNegToPrimVertex_Xi = new TH1F("fHistDcaNegToPrimVertex_Xi", "DCA of V0 neg daughter to Prim. Vertex;DCA (cm);Counts", 300, 0, 3);
	fListHist_Cascade->Add(fHistDcaNegToPrimVertex_Xi);
}




	// - Effective mass histos for cascades.
  
if (! fHistMassXiMinus) {
    fHistMassXiMinus = new TH1F("fHistMassXiMinus","#Xi^{-} candidates;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 200,1.2,2.0);
    fListHist_Cascade->Add(fHistMassXiMinus);
}
  
if (! fHistMassXiPlus) {
    fHistMassXiPlus = new TH1F("fHistMassXiPlus","#Xi^{+} candidates;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",200,1.2,2.0);
    fListHist_Cascade->Add(fHistMassXiPlus);
}

if (! fHistMassOmegaMinus) {
    fHistMassOmegaMinus = new TH1F("fHistMassOmegaMinus","#Omega^{-} candidates;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 200,1.2,2.0);
    fListHist_Cascade->Add(fHistMassOmegaMinus);
}
 
if (! fHistMassOmegaPlus) {
    fHistMassOmegaPlus = new TH1F("fHistMassOmegaPlus","#Omega^{+} candidates;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts",200,1.2,2.0);
    fListHist_Cascade->Add(fHistMassOmegaPlus);
}



}// end CreateOutputObjects






//________________________________________________________________________
void AliAnalysisTaskESDCheckCascade::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  if (!fESD) {
    Printf("ERROR: fESD not available \n");
    return;
  }
 
  //  Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());





  // ---------------------------------------------------------------
  // I - Initialisation of the part dedicated to cascade vertices
  
  Int_t ncascades = -1;
  ncascades = fESD->GetNumberOfCascades();

	// - General histos (filled for any event)

	fHistTrackMultiplicity->Fill(fESD->GetNumberOfTracks());
	fHistCascadeMultiplicity->Fill( ncascades );


  
  	// - 1st part of initialisation : variables needed to store AliESDCascade data members
	Double_t lEffMassXi  = 0. ;
	Double_t lChi2Xi = 0. ;
	Double_t lDcaXiDaughters = 0. ;
	Double_t lXiCosineOfPointingAngle = 0. ;
	Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };
	Double_t lXiRadius= 0. ;
	
	Double_t lXiMom_x  = 0. ;	Double_t lXiMom_y  = 0. ;	Double_t lXiMom_z   = 0. ;
	Double_t lXiMom_transv  = 0. ;
	Double_t lXiMom_tot  = 0. ;
	
	Double_t lBachMom_x  = 0. ;	Double_t lBachMom_y  = 0. ;	Double_t lBachMom_z   = 0. ;
	Double_t lBachMom_transv  = 0. ;
	Double_t lBachMom_tot  = 0. ;
	
	
	
	// - 2nd part of initialisation : about V0 part in cascades

	Double_t lInvMassLambdaAsCascDghter = 0.;
	Double_t lV0Chi2_Xi = 0. ;
	Double_t lDcaV0Daughters_Xi = 0.;
	
	Double_t lDcaBachToPrimVertex_Xi = 0. , lDcaV0ToPrimVertex_Xi = 0.;
	Double_t lDcaPosToPrimVertex_Xi = 0 ;
	Double_t lDcaNegToPrimVertex_Xi = 0;
	Double_t lV0CosineOfPointingAngle_Xi = 0;
	Double_t lPosV0_Xi[3] = { -1000. , -1000., -1000. };
	Double_t lV0Radius_Xi = -1000.0;

	
	// - 3rd part of initialisation : extra tests
  	Double_t lV0quality = 0.;
		Double_t lInvMassXiMinus = 0.;
		Double_t lInvMassXiPlus = 0.;
		Double_t lInvMassOmegaMinus = 0.;
		Double_t lInvMassOmegaPlus = 0.;
		
	
	
  
  
  // -------------------------------------
  // II - Calcultaion Part dedicated to Xi vertices
    
   
  for (Int_t iXi = 0; iXi < ncascades; iXi++) 
    {// This is the begining of the Cascade loop
      AliESDcascade *xi = fESD->GetCascade(iXi);
          
      		if (!xi) continue;
		
      // Just to know which file is currently open : locate the file containing Xi
     cout << "Name of the file containing Xi candidate(s) :" <<  fesdH->GetTree()->GetCurrentFile()->GetName() << endl;
	

		// - II.Step 0 : Characteristics of the event : prim. Vtx + magentic field
		//-------------
	const AliESDVertex *PrimaryVtx_Tracking = fESD->GetPrimaryVertex();  
		// get the vtx stored in ESD found with tracks
	Double_t lTrkgPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
		PrimaryVtx_Tracking->GetXYZ( lTrkgPrimaryVtxPos );
		Double_t lTrkgPrimaryVtxRadius = TMath::Sqrt( lTrkgPrimaryVtxPos[0]*lTrkgPrimaryVtxPos[0] +
						lTrkgPrimaryVtxPos[1] * lTrkgPrimaryVtxPos[1] +
						lTrkgPrimaryVtxPos[2] * lTrkgPrimaryVtxPos[2] );
									
	const AliESDVertex *PrimaryVtx_SPDonly 	= fESD->GetPrimaryVertexSPD(); 		
		// get the vtx found by exclusive use of SPD
	Double_t lSPDPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
		PrimaryVtx_SPDonly->GetXYZ( lSPDPrimaryVtxPos );
		Double_t lSPDPrimaryVtxRadius = TMath::Sqrt( lSPDPrimaryVtxPos[0]*lSPDPrimaryVtxPos[0] +
						lSPDPrimaryVtxPos[1] * lSPDPrimaryVtxPos[1] +
						lSPDPrimaryVtxPos[2] * lSPDPrimaryVtxPos[2] );
	
	// As done in AliCascadeVertexer, we keep, between both retrieved vertices, the one which is the best one available.
	// This one will be used for next calculations (DCA essentially)
		
		Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
		if( PrimaryVtx_Tracking->GetStatus() ) { // if tracking vtx = ok
			lBestPrimaryVtxPos[0] = lTrkgPrimaryVtxPos[0];
			lBestPrimaryVtxPos[1] = lTrkgPrimaryVtxPos[1];
			lBestPrimaryVtxPos[2] = lTrkgPrimaryVtxPos[2];
		}
		else{
			lBestPrimaryVtxPos[0] = lSPDPrimaryVtxPos[0];
			lBestPrimaryVtxPos[1] = lSPDPrimaryVtxPos[1];
			lBestPrimaryVtxPos[2] = lSPDPrimaryVtxPos[2];
		}

		
	Double_t lMagneticField = fESD->GetMagneticField( );
	
	
	
		// - II.Step 1 : Assigning the necessary variables for specific AliESDcascade data members
		//-------------
	lV0quality = 0.;
	xi->ChangeMassHypothesis(lV0quality , 3312); // default working hypothesis : cascade = Xi- decay

	lEffMassXi  			= xi->GetEffMassXi();
	lChi2Xi 			= xi->GetChi2Xi();
	lDcaXiDaughters			= xi->GetDcaXiDaughters();
	lXiCosineOfPointingAngle 	= xi->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
		// Take care : the best available vertex should be used (like in AliCascadeVertexer)
	
	xi->GetXYZcascade( lPosXi[0],  lPosXi[1], lPosXi[2] ); 
		lXiRadius		= TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
		
	xi->GetPxPyPz( lXiMom_x, lXiMom_y, lXiMom_z );
		lXiMom_transv  	= TMath::Sqrt( lXiMom_x*lXiMom_x   + lXiMom_y*lXiMom_y );
		lXiMom_tot  	= TMath::Sqrt( lXiMom_x*lXiMom_x   + lXiMom_y*lXiMom_y   + lXiMom_z*lXiMom_z );
		
				
		
	xi->GetBPxPyPz(  lBachMom_x,  lBachMom_y,  lBachMom_z );
		lBachMom_transv  	= TMath::Sqrt( lBachMom_x*lBachMom_x   + lBachMom_y*lBachMom_y );
		lBachMom_tot  		= TMath::Sqrt( lBachMom_x*lBachMom_x   + lBachMom_y*lBachMom_y  +  lBachMom_z*lBachMom_z  );
		

	

		// - II.Step 2 : around the tracks : Bach + V0 
		// ~ Necessary variables for ESDcascade data members coming from the ESDv0 part (inheritance)
		//-------------
		
		UInt_t lIdxPos_Xi 	= (UInt_t) TMath::Abs( xi->GetPindex() );
		UInt_t lIdxNeg_Xi 	= (UInt_t) TMath::Abs( xi->GetNindex() );
		UInt_t lBachIdx 	= (UInt_t) TMath::Abs( xi->GetBindex() );
			// Care track label can be negative in MC production (linked with the track quality)

	AliESDtrack *pTrack_Xi		= fESD->GetTrack( lIdxPos_Xi );
	AliESDtrack *nTrack_Xi		= fESD->GetTrack( lIdxNeg_Xi );
	AliESDtrack *bachTrack_Xi	= fESD->GetTrack( lBachIdx );
	if (!pTrack_Xi || !nTrack_Xi || !bachTrack_Xi ) {
		Printf("ERROR: Could not retrieve one of the 3 daughter tracks of the cascade ...");
		continue;
	}
	
	lInvMassLambdaAsCascDghter =	xi->GetEffMass();
		// This value shouldn't change, whatever the working hyp. is : Xi-, Xi+, Omega-, Omega+
	lDcaV0Daughters_Xi =    	xi->GetDcaV0Daughters(); 
	lV0Chi2_Xi =    		xi->GetChi2V0();
	
	lV0CosineOfPointingAngle_Xi = 	xi->GetV0CosineOfPointingAngle( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
			// Maybe a pointing angle towards the Xi Vertex would also be appropriate
	lDcaV0ToPrimVertex_Xi =	xi->GetD( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
	
	
	
	lDcaBachToPrimVertex_Xi  = bachTrack_Xi->GetD( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lMagneticField  ); // return an algebraic value ...
	lDcaBachToPrimVertex_Xi  = TMath::Abs( lDcaBachToPrimVertex_Xi );
	
	
	
	
		xi->GetXYZ( lPosV0_Xi[0],  lPosV0_Xi[1], lPosV0_Xi[2] ); 
		lV0Radius_Xi	= TMath::Sqrt( lPosV0_Xi[0]*lPosV0_Xi[0]  +  lPosV0_Xi[1]*lPosV0_Xi[1] );
	
	Float_t  tDcaPosToPrimVertex_Xi[2];
	if(pTrack_Xi) pTrack_Xi->GetImpactParameters(tDcaPosToPrimVertex_Xi[0],tDcaPosToPrimVertex_Xi[1]);
		// void GetImpactParameters(Float_t &xy,Float_t &z) const {xy=fD; z=fZ;}
	else { 
			tDcaPosToPrimVertex_Xi[0]= -999.;  
			tDcaPosToPrimVertex_Xi[1]= -999.;
		}
	lDcaPosToPrimVertex_Xi = TMath::Sqrt(tDcaPosToPrimVertex_Xi[0]*tDcaPosToPrimVertex_Xi[0]
					+   tDcaPosToPrimVertex_Xi[1]*tDcaPosToPrimVertex_Xi[1]);

	Float_t  tDcaNegToPrimVertex_Xi[2];
	if(nTrack_Xi) nTrack_Xi->GetImpactParameters(tDcaNegToPrimVertex_Xi[0],tDcaNegToPrimVertex_Xi[1]);
	else { 
			tDcaNegToPrimVertex_Xi[0]= -999.;  
			tDcaNegToPrimVertex_Xi[1]= -999.;
		}
	lDcaNegToPrimVertex_Xi = TMath::Sqrt(tDcaNegToPrimVertex_Xi[0]*tDcaNegToPrimVertex_Xi[0]
					+   tDcaNegToPrimVertex_Xi[1]*tDcaNegToPrimVertex_Xi[1]);
	
		

  // -------------------------------------
  // III - Filling the TH1Fs
  

		// - III.Step 1	
		fHistPos_TrkgPrimaryVtx_x->Fill( lTrkgPrimaryVtxPos[0]  );
		fHistPos_TrkgPrimaryVtx_y->Fill( lTrkgPrimaryVtxPos[1]  );	
		fHistPos_TrkgPrimaryVtx_z->Fill( lTrkgPrimaryVtxPos[2]  ); 
		fHist_TrkgPrimaryVtxRadius->Fill( lTrkgPrimaryVtxRadius );
		
		fHistPos_SPDPrimaryVtx_x->Fill( lSPDPrimaryVtxPos[0] );
		fHistPos_SPDPrimaryVtx_y->Fill( lSPDPrimaryVtxPos[1] );	
		fHistPos_SPDPrimaryVtx_z->Fill( lSPDPrimaryVtxPos[2] ); 
		fHist_SPDPrimaryVtxRadius->Fill( lSPDPrimaryVtxRadius  );
		
		f2dHistTrkgPrimVtxVsSPDPrimVtx->Fill( lTrkgPrimaryVtxRadius, lSPDPrimaryVtxRadius );

			
		// - III.Step 2
		fHistEffMassXi->Fill( lEffMassXi   );
		fHistChi2Xi->Fill( lChi2Xi  );			// Flag CascadeVtxer: Cut Variable a
		fHistDcaXiDaughters->Fill( lDcaXiDaughters );	// Flag CascadeVtxer: Cut Variable e 
		fHistDcaBachToPrimVertex->Fill( lDcaBachToPrimVertex_Xi  );	// Flag CascadeVtxer: Cut Variable d
		fHistXiCosineOfPointingAngle->Fill( lXiCosineOfPointingAngle ); // Flag CascadeVtxer: Cut Variable f
		fHistXiRadius->Fill( lXiRadius );		// Flag CascadeVtxer: Cut Variable g+h
		
		fHistXiMom_transv->Fill( lXiMom_transv );
		fHistXiMom_tot->Fill( lXiMom_tot );
		
		fHistBachMom_transv->Fill( lBachMom_transv );
		fHistBachMom_tot->Fill( lBachMom_tot );
		

		// - III.Step 3
		fHistMassLambdaAsCascDghter->Fill(  lInvMassLambdaAsCascDghter  ); // Flag CascadeVtxer: Cut Variable c
		fHistV0Chi2_Xi->Fill( lV0Chi2_Xi  );	
		fHistDcaV0Daughters_Xi->Fill( lDcaV0Daughters_Xi );
		fHistV0CosineOfPointingAngle_Xi->Fill( lV0CosineOfPointingAngle_Xi ); 
		fHistV0Radius_Xi->Fill( lV0Radius_Xi );
		
		
		fHistDcaV0ToPrimVertex_Xi->Fill( lDcaV0ToPrimVertex_Xi );	// Flag CascadeVtxer: Cut Variable b
		fHistDcaPosToPrimVertex_Xi->Fill( lDcaPosToPrimVertex_Xi );
		fHistDcaNegToPrimVertex_Xi->Fill( lDcaNegToPrimVertex_Xi );
		

		
	
		// - III.Step 4
       lV0quality = 0.;
       xi->ChangeMassHypothesis(lV0quality , 3312); 	// Calculate the effective mass of the Xi- candidate. 
							// pdg code 3312 = Xi-
		if( bachTrack_Xi->Charge()  < 0 ){
				lInvMassXiMinus = xi->GetEffMassXi();
				fHistMassXiMinus->Fill( lInvMassXiMinus );
		} // end if bachelor = negatively charged ( "Xi" = Xi-)
	
	
       lV0quality = 0.;
       xi->ChangeMassHypothesis(lV0quality , -3312); 	// Calculate the effective mass of the Xi+ candidate. 
      							// pdg code -3312 = Xi+
		if( bachTrack_Xi->Charge()  >  0 ){
				lInvMassXiPlus = xi->GetEffMassXi();
				fHistMassXiPlus->Fill( lInvMassXiPlus );			
		} // end if bachelor = positively charged ( "Xi" = Xi+)
	
	
	lV0quality = 0.;
        xi->ChangeMassHypothesis(lV0quality , 3334); 	// Calculate the effective mass of the Xi- candidate. 
							// pdg code 3334 = Omega-
		if( bachTrack_Xi->Charge()  < 0 ){
				lInvMassOmegaMinus = xi->GetEffMassXi();
				fHistMassOmegaMinus->Fill( lInvMassOmegaMinus );
		} // end if bachelor = negatively charged ( "Xi" = Omega-)
	
	
	lV0quality = 0.;
        xi->ChangeMassHypothesis(lV0quality , -3334); 	// Calculate the effective mass of the Xi+ candidate. 
      							// pdg code -3334  = Omega+
		if( bachTrack_Xi->Charge()  >  0 ){
				lInvMassOmegaPlus = xi->GetEffMassXi();
				fHistMassOmegaPlus->Fill( lInvMassOmegaPlus );			
		} // end if bachelor = positively charged ( "Xi" = Omega+)		
			

 
    }// This is the end of the Cascade loop

 
  // Post output data.
 PostData(0, fListHist_Cascade);
}










//________________________________________________________________________
void AliAnalysisTaskESDCheckCascade::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fHistTrackMultiplicity = dynamic_cast<TH1F*> (  ((TList*)GetOutputData(0))->FindObject("fHistTrackMultiplicity")  );
  if (!fHistTrackMultiplicity) {
    Printf("ERROR: fHistTrackMultiplicity not available");
    return;
  }
 
  fHistCascadeMultiplicity = dynamic_cast<TH1F*> (  ((TList*)GetOutputData(0))->FindObject("fHistCascadeMultiplicity")  );
  if (!fHistCascadeMultiplicity) {
    Printf("ERROR: fHistCascadeMultiplicity not available");
    return;
  }
   
  TCanvas *c2 = new TCanvas("AliAnalysisTaskESDCheckCascade","Multiplicity",10,10,510,510);
  c2->cd(1)->SetLogy();

  fHistTrackMultiplicity->SetMarkerStyle(22);
  fHistTrackMultiplicity->DrawCopy("E");
  fHistCascadeMultiplicity->SetMarkerStyle(26);
  fHistCascadeMultiplicity->DrawCopy("ESAME");
  
}
