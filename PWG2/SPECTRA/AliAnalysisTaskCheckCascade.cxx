
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
//                 AliAnalysisTaskCheckCascade class
//            This task is for QAing the Cascades from ESD and AOD 
//            Origin : Antonin Maire Fev2008, antonin.maire@ires.in2p3.fr
//	      Modified : A.M Juillet 2008
//-----------------------------------------------------------------




class TTree;
class TParticle;
class TVector3;

//class AliMCEventHandler;
//class AliMCEvent;
//class AliStack;

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;

#include <iostream>

#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"

#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
//#include "AliCascadeVertexer.h"

#include "AliESDcascade.h"
#include "AliAODcascade.h"

#include "AliAnalysisTaskCheckCascade.h"

ClassImp(AliAnalysisTaskCheckCascade)



//________________________________________________________________________
AliAnalysisTaskCheckCascade::AliAnalysisTaskCheckCascade() 
  : AliAnalysisTaskSE(), fAnalysisType("ESD"), fCollidingSystems(0),

    	// - Cascade part initialisation
    fListHistCascade(0),
    fHistTrackMultiplicity(0), fHistCascadeMultiplicity(0),
    fHistVtxStatus(0),

    fHistPosTrkgPrimaryVtxX(0), fHistPosTrkgPrimaryVtxY(0), fHistPosTrkgPrimaryVtxZ(0), fHistTrkgPrimaryVtxRadius(0),
    fHistPosBestPrimaryVtxX(0), fHistPosBestPrimaryVtxY(0), fHistPosBestPrimaryVtxZ(0), fHistBestPrimaryVtxRadius(0),
    f2dHistTrkgPrimVtxVsBestPrimVtx(0),

    fHistEffMassXi(0),  fHistChi2Xi(0),  
    fHistDcaXiDaughters(0), fHistDcaBachToPrimVertex(0), fHistXiCosineOfPointingAngle(0), fHistXiRadius(0),

    fHistMassLambdaAsCascDghter(0),
    fHistV0Chi2Xi(0),
    fHistDcaV0DaughtersXi(0),
    fHistDcaV0ToPrimVertexXi(0), 
    fHistV0CosineOfPointingAngleXi(0),
    fHistV0RadiusXi(0),
    fHistDcaPosToPrimVertexXi(0), fHistDcaNegToPrimVertexXi(0), 

    fHistMassXiMinus(0), fHistMassXiPlus(0),
    fHistMassOmegaMinus(0), fHistMassOmegaPlus(0),

    fHistXiTransvMom(0),    fHistXiTotMom(0),
    fHistBachTransvMom(0),   fHistBachTotMom(0),

    fHistChargeXi(0),
    fHistV0toXiCosineOfPointingAngle(0),

    fHistRapXi(0), fHistRapOmega(0), fHistEta(0),
    fHistTheta(0), fHistPhi(0),

    f2dHistArmenteros(0),			
    f2dHistEffMassLambdaVsEffMassXiMinus(0), f2dHistEffMassXiVsEffMassOmegaMinus(0),
    f2dHistEffMassLambdaVsEffMassXiPlus(0), f2dHistEffMassXiVsEffMassOmegaPlus(0),
    f2dHistXiRadiusVsEffMassXiMinus(0), f2dHistXiRadiusVsEffMassXiPlus(0)

{
  // Dummy Constructor
}








//________________________________________________________________________
AliAnalysisTaskCheckCascade::AliAnalysisTaskCheckCascade(const char *name) 
  : AliAnalysisTaskSE(name), fAnalysisType("ESD"), fCollidingSystems(0),
     
    	// - Cascade part initialisation
    fListHistCascade(0),
    fHistTrackMultiplicity(0), fHistCascadeMultiplicity(0),
    fHistVtxStatus(0),

    fHistPosTrkgPrimaryVtxX(0), fHistPosTrkgPrimaryVtxY(0), fHistPosTrkgPrimaryVtxZ(0), fHistTrkgPrimaryVtxRadius(0),
    fHistPosBestPrimaryVtxX(0), fHistPosBestPrimaryVtxY(0), fHistPosBestPrimaryVtxZ(0), fHistBestPrimaryVtxRadius(0),
    f2dHistTrkgPrimVtxVsBestPrimVtx(0),

    fHistEffMassXi(0),  fHistChi2Xi(0),  
    fHistDcaXiDaughters(0), fHistDcaBachToPrimVertex(0), fHistXiCosineOfPointingAngle(0), fHistXiRadius(0),

    fHistMassLambdaAsCascDghter(0),
    fHistV0Chi2Xi(0),
    fHistDcaV0DaughtersXi(0),
    fHistDcaV0ToPrimVertexXi(0), 
    fHistV0CosineOfPointingAngleXi(0),
    fHistV0RadiusXi(0),
    fHistDcaPosToPrimVertexXi(0), fHistDcaNegToPrimVertexXi(0), 

    fHistMassXiMinus(0), fHistMassXiPlus(0),
    fHistMassOmegaMinus(0), fHistMassOmegaPlus(0),

    fHistXiTransvMom(0),    fHistXiTotMom(0),
    fHistBachTransvMom(0),   fHistBachTotMom(0),

    fHistChargeXi(0),
    fHistV0toXiCosineOfPointingAngle(0),

    fHistRapXi(0), fHistRapOmega(0), fHistEta(0),
    fHistTheta(0), fHistPhi(0),

    f2dHistArmenteros(0),			
    f2dHistEffMassLambdaVsEffMassXiMinus(0), f2dHistEffMassXiVsEffMassOmegaMinus(0),
    f2dHistEffMassLambdaVsEffMassXiPlus(0), f2dHistEffMassXiVsEffMassOmegaPlus(0),
    f2dHistXiRadiusVsEffMassXiMinus(0), f2dHistXiRadiusVsEffMassXiPlus(0)



{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  
  // Output slot #0 writes into a TList container (Cascade)
  DefineOutput(1, TList::Class());
}






//________________________________________________________________________
void AliAnalysisTaskCheckCascade::UserCreateOutputObjects()
{
  // Create histograms
  // Called once



 fListHistCascade = new TList();

  
	// - General histos
  
if(! fHistTrackMultiplicity) {	
	if(fCollidingSystems)// AA collisions	
		fHistTrackMultiplicity = new TH1F("fHistTrackMultiplicity", 
			"Multiplicity distribution;Number of tracks;Events", 
			200, 0, 40000); 		
	else // pp collisions
		fHistTrackMultiplicity = new TH1F("fHistTrackMultiplicity", 
			"Track Multiplicity;Nbr of tracks/Evt;Events", 
			200, 0, 200);
	fListHistCascade->Add(fHistTrackMultiplicity);
}

if(! fHistCascadeMultiplicity) {
	if(fCollidingSystems)// AA collisions
		fHistCascadeMultiplicity = new TH1F("fHistCascadeMultiplicity", 
			"Multiplicity distribution;Number of Cascades;Events", 
			25, 0, 25); 		
	else // pp collisions
	  	fHistCascadeMultiplicity = new TH1F("fHistCascadeMultiplicity", 
			"Cascades per event;Nbr of Cascades/Evt;Events", 
			10, 0, 10);
	fListHistCascade->Add(fHistCascadeMultiplicity);
}



if(! fHistVtxStatus ){
	fHistVtxStatus   = new TH1F( "fHistVtxStatus" , "Does a Trckg Prim.vtx exist ?; true=1 or false=0; Nb of Events" , 4, -1.0, 3.0 );
	fListHistCascade->Add(fHistVtxStatus);
}






	// - Vertex Positions
  
if(! fHistPosTrkgPrimaryVtxX ){
	fHistPosTrkgPrimaryVtxX   = new TH1F( "fHistPosTrkgPrimaryVtxX" , "Trkg Prim. Vertex Position in x; x (cm); Events" , 200, -0.5, 0.5 );
	fListHistCascade->Add(fHistPosTrkgPrimaryVtxX);
}


if(! fHistPosTrkgPrimaryVtxY){
	fHistPosTrkgPrimaryVtxY   = new TH1F( "fHistPosTrkgPrimaryVtxY" , "Trkg Prim. Vertex Position in y; y (cm); Events" , 200, -0.5, 0.5 );
	fListHistCascade->Add(fHistPosTrkgPrimaryVtxY);
}

if(! fHistPosTrkgPrimaryVtxZ ){
	fHistPosTrkgPrimaryVtxZ   = new TH1F( "fHistPosTrkgPrimaryVtxZ" , "Trkg Prim. Vertex Position in z; z (cm); Events" , 100, -15.0, 15.0 );
	fListHistCascade->Add(fHistPosTrkgPrimaryVtxZ);
}

if(! fHistTrkgPrimaryVtxRadius ){
	fHistTrkgPrimaryVtxRadius  = new TH1F( "fHistTrkgPrimaryVtxRadius",  "Trkg Prim. Vertex radius; r (cm); Events" , 150, 0., 15.0 );
	fListHistCascade->Add(fHistTrkgPrimaryVtxRadius);
}




if(! fHistPosBestPrimaryVtxX ){
	fHistPosBestPrimaryVtxX   = new TH1F( "fHistPosBestPrimaryVtxX" , "Best Prim. Vertex Position in x; x (cm); Events" , 200, -0.5, 0.5 );
	fListHistCascade->Add(fHistPosBestPrimaryVtxX);
}

if(! fHistPosBestPrimaryVtxY){
	fHistPosBestPrimaryVtxY   = new TH1F( "fHistPosBestPrimaryVtxY" , "Best Prim. Vertex Position in y; y (cm); Events" , 200, -0.5, 0.5 );
	fListHistCascade->Add(fHistPosBestPrimaryVtxY);
}

if(! fHistPosBestPrimaryVtxZ ){
	fHistPosBestPrimaryVtxZ   = new TH1F( "fHistPosBestPrimaryVtxZ" , "Best Prim. Vertex Position in z; z (cm); Events" , 100, -15.0, 15.0 );
	fListHistCascade->Add(fHistPosBestPrimaryVtxZ);
}

if(! fHistBestPrimaryVtxRadius ){
	fHistBestPrimaryVtxRadius  = new TH1F( "fHistBestPrimaryVtxRadius",  "Best Prim.  vertex radius; r (cm); Events" , 150, 0., 15.0 );
	fListHistCascade->Add(fHistBestPrimaryVtxRadius);
}

if(! f2dHistTrkgPrimVtxVsBestPrimVtx) {
	f2dHistTrkgPrimVtxVsBestPrimVtx = new TH2F( "f2dHistTrkgPrimVtxVsBestPrimVtx", "r_{Trck Prim. Vtx} Vs r_{Best Prim. Vtx}; r_{Track Vtx} (cm); r_{Best Vtx} (cm)", 300, 0., 15.0, 300, 0., 15.);
	fListHistCascade->Add(f2dHistTrkgPrimVtxVsBestPrimVtx);
}




// - Typical histos for cascades


if(! fHistEffMassXi) {
     fHistEffMassXi = new TH1F("fHistEffMassXi", "Cascade candidates ; Invariant Mass (GeV/c^{2}) ; Counts", 200, 1.2, 2.0);
     fListHistCascade->Add(fHistEffMassXi);
}
   
if(! fHistChi2Xi ){
	fHistChi2Xi = new TH1F("fHistChi2Xi", "Cascade #chi^{2}; #chi^{2}; Number of Cascades", 160, 0, 40);
	fListHistCascade->Add(fHistChi2Xi);
}
  
if(! fHistDcaXiDaughters ){
	fHistDcaXiDaughters = new TH1F( "fHistDcaXiDaughters",  "DCA between Xi Daughters; DCA (cm) ; Number of Cascades", 100, 0., 0.5);
	fListHistCascade->Add(fHistDcaXiDaughters);
}

if(! fHistDcaBachToPrimVertex) {
	fHistDcaBachToPrimVertex = new TH1F("fHistDcaBachToPrimVertex", "DCA of Bach. to Prim. Vertex;DCA (cm);Number of Cascades", 250, 0., 0.25);
	fListHistCascade->Add(fHistDcaBachToPrimVertex);
}

if(! fHistXiCosineOfPointingAngle) {
	fHistXiCosineOfPointingAngle = new TH1F("fHistXiCosineOfPointingAngle", "Cosine of Xi Pointing Angle; Cos (Xi Point.Angl);Number of Xis", 200, 0.98, 1.0);
	fListHistCascade->Add(fHistXiCosineOfPointingAngle);
}

if(! fHistXiRadius ){
	fHistXiRadius  = new TH1F( "fHistXiRadius",  "Casc. decay transv. radius; r (cm); Counts" , 200, 0., 20.0 );
	fListHistCascade->Add(fHistXiRadius);
}


// - Histos about ~ the "V0 part" of the cascade,  coming by inheritance from AliESDv0



if (! fHistMassLambdaAsCascDghter) {
     fHistMassLambdaAsCascDghter = new TH1F("fHistMassLambdaAsCascDghter","#Lambda associated to Casc. candidates;Eff. Mass (GeV/c^{2});Counts", 160,1.00,1.8);
    fListHistCascade->Add(fHistMassLambdaAsCascDghter);
}

if (! fHistV0Chi2Xi) {
	fHistV0Chi2Xi = new TH1F("fHistV0Chi2Xi", "V0 #chi^{2}, in cascade; #chi^{2};Counts", 160, 0, 40);
	fListHistCascade->Add(fHistV0Chi2Xi);
}

if (! fHistDcaV0DaughtersXi) {
	fHistDcaV0DaughtersXi = new TH1F("fHistDcaV0DaughtersXi", "DCA between V0 daughters, in cascade;DCA (cm);Number of V0s", 120, 0., 0.6);
	fListHistCascade->Add(fHistDcaV0DaughtersXi);
}

if (! fHistDcaV0ToPrimVertexXi) {
	fHistDcaV0ToPrimVertexXi = new TH1F("fHistDcaV0ToPrimVertexXi", "DCA of V0 to Prim. Vertex, in cascade;DCA (cm);Number of Cascades", 200, 0., 1.);
	fListHistCascade->Add(fHistDcaV0ToPrimVertexXi);
}

if (! fHistV0CosineOfPointingAngleXi) {
	fHistV0CosineOfPointingAngleXi = new TH1F("fHistV0CosineOfPointingAngleXi", "Cosine of V0 Pointing Angle, in cascade;Cos(V0 Point. Angl); Counts", 200, 0.98, 1.0);
	fListHistCascade->Add(fHistV0CosineOfPointingAngleXi);
}

if (! fHistV0RadiusXi) {
	fHistV0RadiusXi = new TH1F("fHistV0RadiusXi", "V0 decay radius, in cascade; radius (cm); Counts", 200, 0, 20);
	fListHistCascade->Add(fHistV0RadiusXi);
}

if (! fHistDcaPosToPrimVertexXi) {
	fHistDcaPosToPrimVertexXi = new TH1F("fHistDcaPosToPrimVertexXi", "DCA of V0 pos daughter to Prim. Vertex;DCA (cm);Counts", 300, 0, 3);
	fListHistCascade->Add(fHistDcaPosToPrimVertexXi);
}

if (! fHistDcaNegToPrimVertexXi) {
	fHistDcaNegToPrimVertexXi = new TH1F("fHistDcaNegToPrimVertexXi", "DCA of V0 neg daughter to Prim. Vertex;DCA (cm);Counts", 300, 0, 3);
	fListHistCascade->Add(fHistDcaNegToPrimVertexXi);
}




	// - Effective mass histos for cascades.
  
if (! fHistMassXiMinus) {
    fHistMassXiMinus = new TH1F("fHistMassXiMinus","#Xi^{-} candidates;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 200,1.2,2.0);
    fListHistCascade->Add(fHistMassXiMinus);
}
  
if (! fHistMassXiPlus) {
    fHistMassXiPlus = new TH1F("fHistMassXiPlus","#Xi^{+} candidates;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",200,1.2,2.0);
    fListHistCascade->Add(fHistMassXiPlus);
}

if (! fHistMassOmegaMinus) {
	fHistMassOmegaMinus = new TH1F("fHistMassOmegaMinus","#Omega^{-} candidates;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 250,1.5,2.5);
    fListHistCascade->Add(fHistMassOmegaMinus);
}
 
if (! fHistMassOmegaPlus) {
	fHistMassOmegaPlus = new TH1F("fHistMassOmegaPlus","#Omega^{+} candidates;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 250,1.5,2.5);
    fListHistCascade->Add(fHistMassOmegaPlus);
}



	// - Complements for QA

if(! fHistXiTransvMom ){
	fHistXiTransvMom  = new TH1F( "fHistXiTransvMom" , "Xi transverse momentum ; p_{t}(#Xi) (GeV/c); Counts", 100, 0.0, 10.0);
	fListHistCascade->Add(fHistXiTransvMom);
}

if(! fHistXiTotMom ){
	fHistXiTotMom  = new TH1F( "fHistXiTotMom" , "Xi momentum norm; p_{tot}(#Xi) (GeV/c); Counts", 150, 0.0, 15.0);
	fListHistCascade->Add(fHistXiTotMom);
}


if(! fHistBachTransvMom ){
	fHistBachTransvMom  = new TH1F( "fHistBachTransvMom" , "Bach. transverse momentum ; p_{t}(Bach.) (GeV/c); Counts", 100, 0.0, 5.0);
	fListHistCascade->Add(fHistBachTransvMom);
}

if(! fHistBachTotMom ){
	fHistBachTotMom  = new TH1F( "fHistBachTotMom" , "Bach. momentum norm; p_{tot}(Bach.) (GeV/c); Counts", 100, 0.0, 5.0);
	fListHistCascade->Add(fHistBachTotMom);
}


if(! fHistChargeXi ){
	fHistChargeXi  = new TH1F( "fHistChargeXi" , "Charge of casc. candidates ; Sign ; Counts", 5, -2.0, 3.0);
	fListHistCascade->Add(fHistChargeXi);
}


if (! fHistV0toXiCosineOfPointingAngle) {
	fHistV0toXiCosineOfPointingAngle = new TH1F("fHistV0toXiCosineOfPointingAngle", "Cos. of V0 Ptng Angl / Xi vtx ;Cos(V0 Point. Angl / Xi vtx); Counts", 100, 0.99, 1.0);
	fListHistCascade->Add(fHistV0toXiCosineOfPointingAngle);
}


if(! fHistRapXi ){
	fHistRapXi  = new TH1F( "fHistRapXi" , "Rapidity of Xi candidates ; y ; Counts", 200, -5.0, 5.0);
	fListHistCascade->Add(fHistRapXi);
}

if(! fHistRapOmega ){
	fHistRapOmega  = new TH1F( "fHistRapOmega" , "Rapidity of Omega candidates ; y ; Counts", 200, -5.0, 5.0);
	fListHistCascade->Add(fHistRapOmega);
}

if(! fHistEta ){
	fHistEta  = new TH1F( "fHistEta" , "Pseudo-rap. of casc. candidates ; #eta ; Counts", 120, -3.0, 3.0);
	fListHistCascade->Add(fHistEta);
}

if(! fHistTheta ){
	fHistTheta  = new TH1F( "fHistTheta" , "#theta of casc. candidates ; #theta (deg) ; Counts", 180, 0., 180.0);
	fListHistCascade->Add(fHistTheta);
}

if(! fHistPhi ){
	fHistPhi  = new TH1F( "fHistPhi" , "#phi of casc. candidates ; #phi (deg) ; Counts", 360, 0., 360.);
	fListHistCascade->Add(fHistPhi);
}


if(! f2dHistArmenteros) {
	f2dHistArmenteros = new TH2F( "f2dHistArmenteros", "#alpha_{Arm}(casc. cand.) Vs Pt_{Arm}(casc. cand.); #alpha_{Arm} ; Pt_{Arm} (GeV/c)", 140, -1.2, 1.2, 300, 0., 0.3);
	fListHistCascade->Add(f2dHistArmenteros);
}

if(! f2dHistEffMassLambdaVsEffMassXiMinus) {
	f2dHistEffMassLambdaVsEffMassXiMinus = new TH2F( "f2dHistEffMassLambdaVsEffMassXiMinus", "M_{#Lambda} Vs M_{#Xi^{-} candidates} ; Inv. M_{#Lambda^{0}} (GeV/c^{2}) ; M( #Lambda , #pi^{-} ) (GeV/c^{2})", 300, 1.1,1.13, 200, 1.2, 2.0);
	fListHistCascade->Add(f2dHistEffMassLambdaVsEffMassXiMinus);
}

if(! f2dHistEffMassXiVsEffMassOmegaMinus) {
	f2dHistEffMassXiVsEffMassOmegaMinus = new TH2F( "f2dHistEffMassXiVsEffMassOmegaMinus", "M_{#Xi^{-} candidates} Vs M_{#Omega^{-} candidates} ; M( #Lambda , #pi^{-} ) (GeV/c^{2}) ; M( #Lambda , K^{-} ) (GeV/c^{2})", 200, 1.2, 2.0, 250, 1.5, 2.5);
	fListHistCascade->Add(f2dHistEffMassXiVsEffMassOmegaMinus);
}

if(! f2dHistEffMassLambdaVsEffMassXiPlus) {
	f2dHistEffMassLambdaVsEffMassXiPlus = new TH2F( "f2dHistEffMassLambdaVsEffMassXiPlus", "M_{#Lambda} Vs M_{#Xi^{+} candidates} ; Inv. M_{#Lambda^{0}} (GeV/c^{2}) ; M( #Lambda , #pi^{+} ) (GeV/c^{2})", 300, 1.1,1.13, 200, 1.2, 2.0);
	fListHistCascade->Add(f2dHistEffMassLambdaVsEffMassXiPlus);
}

if(! f2dHistEffMassXiVsEffMassOmegaPlus) {
	f2dHistEffMassXiVsEffMassOmegaPlus = new TH2F( "f2dHistEffMassXiVsEffMassOmegaPlus", "M_{#Xi^{+} candidates} Vs M_{#Omega^{+} candidates} ; M( #Lambda , #pi^{+} ) (GeV/c^{2}) ; M( #Lambda , K^{+} ) (GeV/c^{2})", 200, 1.2, 2.0, 250, 1.5, 2.5);
	fListHistCascade->Add(f2dHistEffMassXiVsEffMassOmegaPlus);
}

if(! f2dHistXiRadiusVsEffMassXiMinus) {
	f2dHistXiRadiusVsEffMassXiMinus = new TH2F( "f2dHistXiRadiusVsEffMassXiMinus", "Transv. R_{Xi Decay} Vs M_{#Xi^{-} candidates}; r_{Xi} (cm); M( #Lambda , #pi^{-} ) (GeV/c^{2}) ", 450, 0., 45.0, 200, 1.2, 2.0);
	fListHistCascade->Add(f2dHistXiRadiusVsEffMassXiMinus);
}

if(! f2dHistXiRadiusVsEffMassXiPlus) {
	f2dHistXiRadiusVsEffMassXiPlus = new TH2F( "f2dHistXiRadiusVsEffMassXiPlus", "Transv. R_{Xi Decay} Vs M_{#Xi^{+} candidates}; r_{Xi} (cm); M( #Lambda , #pi^{+} ) (GeV/c^{2}) ", 450, 0., 45.0, 200, 1.2, 2.0);
	fListHistCascade->Add(f2dHistXiRadiusVsEffMassXiPlus);
}




}// end UserCreateOutputObjects






//________________________________________________________________________
void AliAnalysisTaskCheckCascade::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
	
	AliESDEvent *lESDevent = 0x0;
	AliAODEvent *lAODevent = 0x0;
	Int_t ncascades = -1;
	
	
  // Connect to the InputEvent	
  // After these lines, we should have an ESD/AOD event + the number of cascades in it.
		
  if(fAnalysisType == "ESD"){
	lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
	if (!lESDevent) {
		Printf("ERROR: lESDevent not available \n");
		return;
	}
	ncascades = lESDevent->GetNumberOfCascades();
  }
  
  if(fAnalysisType == "AOD"){  
	  lAODevent = dynamic_cast<AliAODEvent*>( InputEvent() ); 
  	if (!lAODevent) {
		Printf("ERROR: lAODevent not available \n");
		return;
	}
	ncascades = lAODevent->GetNumberOfCascades();
  }
	
  // ---------------------------------------------------------------
  // I - Initialisation of the local variables that will be needed

  
  	// - 0th part of initialisation : around primary vertex ...
  Short_t  lStatusTrackingPrimVtx  = -2;
  Double_t lTrkgPrimaryVtxPos[3]   = {-100.0, -100.0, -100.0};
  Double_t lTrkgPrimaryVtxRadius3D = -500.0;
  
  Double_t lBestPrimaryVtxPos[3]   = {-100.0, -100.0, -100.0};
  Double_t lBestPrimaryVtxRadius3D = -500.0;
  
  	// - 1st part of initialisation : variables needed to store AliESDCascade data members
  Double_t lEffMassXi      = 0. ;
  Double_t lChi2Xi         = 0. ;
  Double_t lDcaXiDaughters = 0. ;
  Double_t lXiCosineOfPointingAngle = 0. ;
  Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };
  Double_t lXiRadius = 0. ;
	
		
	// - 2nd part of initialisation : about V0 part in cascades

  Double_t lInvMassLambdaAsCascDghter = 0.;
  Double_t lV0Chi2Xi         = 0. ;
  Double_t lDcaV0DaughtersXi = 0.;
	
  Double_t lDcaBachToPrimVertexXi = 0., lDcaV0ToPrimVertexXi = 0.;
  Double_t lDcaPosToPrimVertexXi  = 0.;
  Double_t lDcaNegToPrimVertexXi  = 0.;
  Double_t lV0CosineOfPointingAngleXi = 0. ;
  Double_t lPosV0Xi[3] = { -1000. , -1000., -1000. }; // Position of VO coming from cascade
  Double_t lV0RadiusXi = -1000.0;
  Double_t lV0quality  = 0.;

	
	// - 3rd part of initialisation : Effective masses
  
  Double_t lInvMassXiMinus    = 0.;
  Double_t lInvMassXiPlus     = 0.;
  Double_t lInvMassOmegaMinus = 0.;
  Double_t lInvMassOmegaPlus  = 0.;


	// - 4th part of initialisation : extra info for QA
  
  Double_t lXiMomX       = 0., lXiMomY = 0., lXiMomZ = 0.;
  Double_t lXiTransvMom  = 0. ;
  Double_t lXiTotMom     = 0. ;
	
  Double_t lBachMomX       = 0., lBachMomY  = 0., lBachMomZ   = 0.;
  Double_t lBachTransvMom  = 0.;
  Double_t lBachTotMom     = 0.;

  Short_t  lChargeXi = -2;
  Double_t lV0toXiCosineOfPointingAngle = 0. ;

  Double_t lRapXi   = -20.0, lRapOmega = -20.0,  lEta = -20.0, lTheta = 360., lPhi = 720. ;
  Double_t lAlphaXi = -200., lPtArmXi  = -200.0;


  
  //-------------------------------------------------
  // O - Cascade vertexer (ESD)

	//   if(fAnalysisType == "ESD" ){
	//      lESDevent->ResetCascades();
	//  	AliCascadeVertexer CascVtxer;
	//  	CascVtxer.V0sTracks2CascadeVertices(lESDevent);
	//     }
  
  
  // ---------------------------------------------------------------
  // I - General histos (filled for any event) - (ESD)

  fHistTrackMultiplicity  ->Fill( (InputEvent())->GetNumberOfTracks() );
  fHistCascadeMultiplicity->Fill( ncascades );
  
  
  
  
  for (Int_t iXi = 0; iXi < ncascades; iXi++)
  {// This is the begining of the Cascade loop (ESD or AOD)
	   
	  
  if(fAnalysisType == "ESD"){ 
  
  // -------------------------------------
  // II - Calcultaion Part dedicated to Xi vertices (ESD)
  
	AliESDcascade *xi = lESDevent->GetCascade(iXi);
	if (!xi) continue;

	// Just to know which file is currently open : locate the file containing Xi
	// cout << "Name of the file containing Xi candidate(s) :" 
	//	<< fInputHandler->GetTree()->GetCurrentFile()->GetName() 
	//	<< endl;
	

		// - II.Step 1 : Characteristics of the event : prim. Vtx + magnetic field (ESD)
		//-------------
	
	// For new code (v4-16-Release-Rev06 or trunk)
	
	const AliESDVertex *lPrimaryTrackingVtx = lESDevent->GetPrimaryVertexTracks();  
		// get the vtx stored in ESD found with tracks
	
		lPrimaryTrackingVtx->GetXYZ( lTrkgPrimaryVtxPos );
		lTrkgPrimaryVtxRadius3D = TMath::Sqrt(  lTrkgPrimaryVtxPos[0] * lTrkgPrimaryVtxPos[0] +
							lTrkgPrimaryVtxPos[1] * lTrkgPrimaryVtxPos[1] +
							lTrkgPrimaryVtxPos[2] * lTrkgPrimaryVtxPos[2] );

		lStatusTrackingPrimVtx = lPrimaryTrackingVtx->GetStatus();


	const AliESDVertex *lPrimaryBestVtx = lESDevent->GetPrimaryVertex();	
		// get the best primary vertex available for the event
		// As done in AliCascadeVertexer, we keep the one which is the best one available.
		// between : Tracking vertex > SPD vertex > TPC vertex > default SPD vertex
		// This one will be used for next calculations (DCA essentially)
	
		lPrimaryBestVtx->GetXYZ( lBestPrimaryVtxPos );
		lBestPrimaryVtxRadius3D = TMath::Sqrt(  lBestPrimaryVtxPos[0] * lBestPrimaryVtxPos[0] +
							lBestPrimaryVtxPos[1] * lBestPrimaryVtxPos[1] +
							lBestPrimaryVtxPos[2] * lBestPrimaryVtxPos[2] );
		
	
	// For older evts
	
	//	const AliESDVertex *lPrimaryTrackingVtx = lESDevent->GetPrimaryVertexTracks();  
	//		get the vtx stored in ESD found with tracks
	// 	Double_t lTrkgPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
	// 		lPrimaryTrackingVtx->GetXYZ( lTrkgPrimaryVtxPos );
	// 		Double_t lTrkgPrimaryVtxRadius3D = TMath::Sqrt( lTrkgPrimaryVtxPos[0]*lTrkgPrimaryVtxPos[0] +
	// 						lTrkgPrimaryVtxPos[1] * lTrkgPrimaryVtxPos[1] +
	// 						lTrkgPrimaryVtxPos[2] * lTrkgPrimaryVtxPos[2] );
	// 
	// 	const AliESDVertex *lPrimarySPDVtx = lESDevent->GetPrimaryVertexSPD();	
	// 		get the vtx found by exclusive use of SPD
	// 	Double_t lSPDPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
	// 		lPrimarySPDVtx->GetXYZ( lSPDPrimaryVtxPos );
	//		
	//	// As done in AliCascadeVertexer, we keep, between both retrieved vertices, 
	//	// the one which is the best one available.
	//	// This one will be used for next calculations (DCA essentially)
	//		
	//		Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
	//		if( lPrimaryTrackingVtx->GetStatus() ) { // if tracking vtx = ok
	//			lBestPrimaryVtxPos[0] = lTrkgPrimaryVtxPos[0];
	//			lBestPrimaryVtxPos[1] = lTrkgPrimaryVtxPos[1];
	//			lBestPrimaryVtxPos[2] = lTrkgPrimaryVtxPos[2];
	//		}
	//		else{
	//			lBestPrimaryVtxPos[0] = lSPDPrimaryVtxPos[0];
	//			lBestPrimaryVtxPos[1] = lSPDPrimaryVtxPos[1];
	//			lBestPrimaryVtxPos[2] = lSPDPrimaryVtxPos[2];
	//		}
	//
	//	Double_t lBestPrimaryVtxRadius3D = TMath::Sqrt( lBestPrimaryVtxPos[0]*lBestPrimaryVtxPos[0] +
	//						lBestPrimaryVtxPos[1] * lBestPrimaryVtxPos[1] +
	//						lBestPrimaryVtxPos[2] * lBestPrimaryVtxPos[2] );
	
	// Back to normal
			
	Double_t lMagneticField = lESDevent->GetMagneticField( );
	
	
	
		// - II.Step 2 : Assigning the necessary variables for specific AliESDcascade data members (ESD)	
		//-------------
	lV0quality = 0.;
	xi->ChangeMassHypothesis(lV0quality , 3312); // default working hypothesis : cascade = Xi- decay

	lEffMassXi  			= xi->GetEffMassXi();
	lChi2Xi 			= xi->GetChi2Xi();
	lDcaXiDaughters			= xi->GetDcaXiDaughters();
	lXiCosineOfPointingAngle 	= xi->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0],
			 							lBestPrimaryVtxPos[1],
										lBestPrimaryVtxPos[2] );
		// Take care : the best available vertex should be used (like in AliCascadeVertexer)
	
	xi->GetXYZcascade( lPosXi[0],  lPosXi[1], lPosXi[2] ); 
	lXiRadius			= TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
		
		

		// - II.Step 3 : around the tracks : Bach + V0 (ESD)
		// ~ Necessary variables for ESDcascade data members coming from the ESDv0 part (inheritance)
		//-------------
		
		UInt_t lIdxPosXi 	= (UInt_t) TMath::Abs( xi->GetPindex() );
		UInt_t lIdxNegXi 	= (UInt_t) TMath::Abs( xi->GetNindex() );
		UInt_t lBachIdx 	= (UInt_t) TMath::Abs( xi->GetBindex() );
			// Care track label can be negative in MC production (linked with the track quality)
			// However = normally, not the case for track index ...

	AliESDtrack *pTrackXi		= lESDevent->GetTrack( lIdxPosXi );
	AliESDtrack *nTrackXi		= lESDevent->GetTrack( lIdxNegXi );
	AliESDtrack *bachTrackXi	= lESDevent->GetTrack( lBachIdx );
	if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
		Printf("ERROR: Could not retrieve one of the 3 daughter tracks of the cascade ...");
		continue;
	}
	
	lInvMassLambdaAsCascDghter	= xi->GetEffMass();
		// This value shouldn't change, whatever the working hyp. is : Xi-, Xi+, Omega-, Omega+
	lDcaV0DaughtersXi 		= xi->GetDcaV0Daughters(); 
	lV0Chi2Xi 			= xi->GetChi2V0();
	
	lV0CosineOfPointingAngleXi 	= xi->GetV0CosineOfPointingAngle( lBestPrimaryVtxPos[0],
									  lBestPrimaryVtxPos[1],
									  lBestPrimaryVtxPos[2] );

	lDcaV0ToPrimVertexXi 		= xi->GetD( lBestPrimaryVtxPos[0], 
						    lBestPrimaryVtxPos[1], 
						    lBestPrimaryVtxPos[2] );
		
	lDcaBachToPrimVertexXi = TMath::Abs( bachTrackXi->GetD(	lBestPrimaryVtxPos[0], 
						     		lBestPrimaryVtxPos[1], 
						     		lMagneticField  ) ); 
					// Note : AliExternalTrackParam::GetD returns an algebraic value ...
		
		xi->GetXYZ( lPosV0Xi[0],  lPosV0Xi[1], lPosV0Xi[2] ); 
	lV0RadiusXi		= TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );
	
	lDcaPosToPrimVertexXi 	= TMath::Abs( pTrackXi	->GetD(	lBestPrimaryVtxPos[0], 
						 		lBestPrimaryVtxPos[1], 
						 		lMagneticField  )     ); 
	
	lDcaNegToPrimVertexXi 	= TMath::Abs( nTrackXi	->GetD(	lBestPrimaryVtxPos[0], 
					    			lBestPrimaryVtxPos[1], 
					    			lMagneticField  )     ); 
	
	
		// - II.Step 4 : around effective masses (ESD)
		// ~ change mass hypotheses to cover all the possibilities :  Xi-/+, Omega -/+
		//-------------

	lV0quality = 0.;
	xi->ChangeMassHypothesis(lV0quality , 3312); 	// Calculate the effective mass of the Xi- candidate. 
							// pdg code 3312 = Xi-
		if( bachTrackXi->Charge() < 0 )	
			lInvMassXiMinus = xi->GetEffMassXi();
			
	lV0quality = 0.;
	xi->ChangeMassHypothesis(lV0quality , -3312); 	// Calculate the effective mass of the Xi+ candidate. 
							// pdg code -3312 = Xi+
		if( bachTrackXi->Charge() >  0 )
			lInvMassXiPlus = xi->GetEffMassXi();
			
	lV0quality = 0.;
	xi->ChangeMassHypothesis(lV0quality , 3334); 	// Calculate the effective mass of the Xi- candidate. 
							// pdg code 3334 = Omega-
		if( bachTrackXi->Charge() < 0 )
			lInvMassOmegaMinus = xi->GetEffMassXi();
		
	lV0quality = 0.;
	xi->ChangeMassHypothesis(lV0quality , -3334); 	// Calculate the effective mass of the Xi+ candidate. 
							// pdg code -3334  = Omega+
		if( bachTrackXi->Charge() >  0 )
			lInvMassOmegaPlus = xi->GetEffMassXi();

	lV0quality = 0.;
	xi->ChangeMassHypothesis(lV0quality , 3312); 	// Back to default hyp.
							

		// - II.Step 5 : extra info for QA (ESD)
		// miscellaneous pieces onf info that may help regarding data quality assessment.
		//-------------

	xi->GetPxPyPz( lXiMomX, lXiMomY, lXiMomZ );
		lXiTransvMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );
		lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY   + lXiMomZ*lXiMomZ );
		
	xi->GetBPxPyPz(  lBachMomX,  lBachMomY,  lBachMomZ );
		lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY );
		lBachTotMom  	= TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY  +  lBachMomZ*lBachMomZ  );

	lChargeXi = xi->Charge();

	lV0toXiCosineOfPointingAngle = xi->GetV0CosineOfPointingAngle( lPosXi[0], lPosXi[1], lPosXi[2] );
	
	lRapXi    = xi->RapXi();
	lRapOmega = xi->RapOmega();
	lEta      = xi->Eta();
	lTheta    = xi->Theta() *180.0/TMath::Pi();
	lPhi      = xi->Phi()   *180.0/TMath::Pi();
	lAlphaXi  = xi->AlphaXi();
	lPtArmXi  = xi->PtArmXi();
	
	
  }// end of ESD treatment
  
 
  if(fAnalysisType == "AOD"){
	
	// -------------------------------------
	// II - Calcultaion Part dedicated to Xi vertices (ESD)
	
	const AliAODcascade *xi = lAODevent->GetCascade(iXi);
	if (!xi) continue;
		
	// Just to know which file is currently open : locate the file containing Xi
	// cout << "Name of the file containing Xi candidate(s) :" <<  fesdH->GetTree()->GetCurrentFile()->GetName() << endl;
	
	
		// - II.Step 1 : Characteristics of the event : prim. Vtx + magnetic field (AOD)
		//-------------
	
	lTrkgPrimaryVtxPos[0]   = -100.0;
	lTrkgPrimaryVtxPos[1]   = -100.0;
	lTrkgPrimaryVtxPos[2]   = -100.0;
	lTrkgPrimaryVtxRadius3D = -500. ;
	// We don't have the different prim. vertex at the AOD level -> nothing to do.

	const AliAODVertex *lPrimaryBestVtx = lAODevent->GetPrimaryVertex();	
	// get the best primary vertex available for the event
	// We may keep the one which is the best one available = GetVertex(0)
	// Pb with pile-up to expect
	// This one will be used for next calculations (DCA essentially)

	lPrimaryBestVtx->GetXYZ( lBestPrimaryVtxPos );
	lBestPrimaryVtxRadius3D = TMath::Sqrt(  lBestPrimaryVtxPos[0] * lBestPrimaryVtxPos[0] +
						lBestPrimaryVtxPos[1] * lBestPrimaryVtxPos[1] +
						lBestPrimaryVtxPos[2] * lBestPrimaryVtxPos[2] );
		
	
		// - II.Step 2 : Assigning the necessary variables for specific AliAODcascade data members (AOD)	
		//-------------
	
	lEffMassXi  			= xi->MassXi(); // default working hypothesis : cascade = Xi- decay
	lChi2Xi 			= xi->Chi2Xi();
	lDcaXiDaughters			= xi->DcaXiDaughters();
	lXiCosineOfPointingAngle 	= xi->CosPointingAngleXi( lBestPrimaryVtxPos[0], 
								  lBestPrimaryVtxPos[1], 
								  lBestPrimaryVtxPos[2] );
					// Take care : 
					// the best available vertex should be used (like in AliCascadeVertexer)

		lPosXi[0] = xi->DecayVertexXiX();
		lPosXi[1] = xi->DecayVertexXiY();
		lPosXi[2] = xi->DecayVertexXiZ();
	lXiRadius = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
		

		// - II.Step 3 : around the tracks : Bach + V0 (AOD)
		// ~ Necessary variables for AODcascade data members coming from the AODv0 part (inheritance)
		//-------------
		
	lChargeXi 			= xi->ChargeXi();
	
	if( lChargeXi < 0)  	
	  lInvMassLambdaAsCascDghter	= xi->MassLambda();
	else 			
	  lInvMassLambdaAsCascDghter	= xi->MassAntiLambda();

	lDcaV0DaughtersXi 		= xi->DcaV0Daughters(); 
	lV0Chi2Xi 			= xi->Chi2V0();

	lV0CosineOfPointingAngleXi 	= xi->CosPointingAngle( lBestPrimaryVtxPos );
	lDcaV0ToPrimVertexXi 		= xi->DcaV0ToPrimVertex();
	
	lDcaBachToPrimVertexXi 		= xi->DcaBachToPrimVertex(); 
	
	
		lPosV0Xi[0] = xi->DecayVertexV0X();
		lPosV0Xi[1] = xi->DecayVertexV0Y();
		lPosV0Xi[2] = xi->DecayVertexV0Z(); 
	lV0RadiusXi	= TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );

	lDcaPosToPrimVertexXi		= xi->DcaPosToPrimVertex(); 
	lDcaNegToPrimVertexXi		= xi->DcaNegToPrimVertex(); 


		// - II.Step 4 : around effective masses (AOD)
		// ~ change mass hypotheses to cover all the possibilities :  Xi-/+, Omega -/+
		//-------------

	if( lChargeXi < 0 )		lInvMassXiMinus 	= xi->MassXi();
	if( lChargeXi > 0 )		lInvMassXiPlus 		= xi->MassXi();
	if( lChargeXi < 0 )		lInvMassOmegaMinus 	= xi->MassOmega();
	if( lChargeXi > 0 )		lInvMassOmegaPlus 	= xi->MassOmega();

	
		// - II.Step 5 : extra info for QA (AOD)
		// miscellaneous pieces onf info that may help regarding data quality assessment.
		//-------------

		lXiMomX = xi->MomXiX();
		lXiMomY = xi->MomXiY();
		lXiMomZ = xi->MomXiZ();
	lXiTransvMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );
	lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY   + lXiMomZ*lXiMomZ );
	
	 	lBachMomX = xi->MomBachX();
	 	lBachMomY = xi->MomBachY();
	 	lBachMomZ = xi->MomBachZ();		
	lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY );
	lBachTotMom  	= TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY  +  lBachMomZ*lBachMomZ  );

	
	lV0toXiCosineOfPointingAngle = xi->CosPointingAngle( xi->GetDecayVertexXi() );
	
	lRapXi    = xi->RapXi();
	lRapOmega = xi->RapOmega();
	lEta      = xi->Eta();				// Will not work ! need a method Pz(), Py() Px() 
	lTheta    = xi->Theta() *180.0/TMath::Pi();     // in AODcascade.
	lPhi      = xi->Phi()   *180.0/TMath::Pi();     // Here, we will get eta, theta, phi for the V0 ...
	lAlphaXi  = xi->AlphaXi();
	lPtArmXi  = xi->PtArmXi();

	 
  }// end of AOD treatment


  // -------------------------------------
  // III - Filling the TH1Fs
  

	// - III.Step 1	

	fHistVtxStatus			->Fill( lStatusTrackingPrimVtx   );  // 1 if tracking vtx = ok

	if( lStatusTrackingPrimVtx ){
		fHistPosTrkgPrimaryVtxX	->Fill( lTrkgPrimaryVtxPos[0]    );
		fHistPosTrkgPrimaryVtxY	->Fill( lTrkgPrimaryVtxPos[1]    );	
		fHistPosTrkgPrimaryVtxZ	->Fill( lTrkgPrimaryVtxPos[2]    ); 
		fHistTrkgPrimaryVtxRadius->Fill( lTrkgPrimaryVtxRadius3D );
	}

	fHistPosBestPrimaryVtxX		->Fill( lBestPrimaryVtxPos[0]    );
	fHistPosBestPrimaryVtxY		->Fill( lBestPrimaryVtxPos[1]    );
	fHistPosBestPrimaryVtxZ		->Fill( lBestPrimaryVtxPos[2]    );
	fHistBestPrimaryVtxRadius	->Fill( lBestPrimaryVtxRadius3D  );
	
	f2dHistTrkgPrimVtxVsBestPrimVtx->Fill( lTrkgPrimaryVtxRadius3D, lBestPrimaryVtxRadius3D );

		
	// - III.Step 2
	fHistEffMassXi			->Fill( lEffMassXi               );
	fHistChi2Xi			->Fill( lChi2Xi                  );	// Flag CascadeVtxer: Cut Variable a
	fHistDcaXiDaughters		->Fill( lDcaXiDaughters          );	// Flag CascadeVtxer: Cut Variable e 
	fHistDcaBachToPrimVertex	->Fill( lDcaBachToPrimVertexXi   );	// Flag CascadeVtxer: Cut Variable d
	fHistXiCosineOfPointingAngle	->Fill( lXiCosineOfPointingAngle ); 	// Flag CascadeVtxer: Cut Variable f
	fHistXiRadius			->Fill( lXiRadius                );	// Flag CascadeVtxer: Cut Variable g+h
	
	
	// - III.Step 3
	fHistMassLambdaAsCascDghter	->Fill( lInvMassLambdaAsCascDghter );	// Flag CascadeVtxer: Cut Variable c
	fHistV0Chi2Xi			->Fill( lV0Chi2Xi                  );	
	fHistDcaV0DaughtersXi		->Fill( lDcaV0DaughtersXi          );
	fHistV0CosineOfPointingAngleXi	->Fill( lV0CosineOfPointingAngleXi ); 
	fHistV0RadiusXi			->Fill( lV0RadiusXi                );
	
	fHistDcaV0ToPrimVertexXi	->Fill( lDcaV0ToPrimVertexXi       );	// Flag CascadeVtxer: Cut Variable b
	fHistDcaPosToPrimVertexXi	->Fill( lDcaPosToPrimVertexXi      );
	fHistDcaNegToPrimVertexXi	->Fill( lDcaNegToPrimVertexXi      );
		
	
	// - III.Step 4
	if( lChargeXi < 0 )	fHistMassXiMinus	->Fill( lInvMassXiMinus );
	if( lChargeXi > 0 )	fHistMassXiPlus		->Fill( lInvMassXiPlus );
	if( lChargeXi < 0 )	fHistMassOmegaMinus	->Fill( lInvMassOmegaMinus );
	if( lChargeXi > 0 )	fHistMassOmegaPlus	->Fill( lInvMassOmegaPlus );

//	if( bachTrackXi->Charge() < 0 )	fHistMassXiMinus	->Fill( lInvMassXiMinus    );
// 	if( bachTrackXi->Charge() > 0 )	fHistMassXiPlus		->Fill( lInvMassXiPlus     );
// 	if( bachTrackXi->Charge() < 0 )	fHistMassOmegaMinus	->Fill( lInvMassOmegaMinus );
// 	if( bachTrackXi->Charge() > 0 )	fHistMassOmegaPlus	->Fill( lInvMassOmegaPlus  );


	// - III.Step 5
	fHistXiTransvMom	->Fill( lXiTransvMom   );
	fHistXiTotMom		->Fill( lXiTotMom      );
		
	fHistBachTransvMom	->Fill( lBachTransvMom );
	fHistBachTotMom		->Fill( lBachTotMom    );

	fHistChargeXi		->Fill( lChargeXi      );
	fHistV0toXiCosineOfPointingAngle->Fill( lV0toXiCosineOfPointingAngle );

	fHistRapXi		->Fill( lRapXi               );
	fHistRapOmega		->Fill( lRapOmega            );
	fHistEta		->Fill( lEta                 );
	fHistTheta		->Fill( lTheta               );
	fHistPhi		->Fill( lPhi                 );

	f2dHistArmenteros	->Fill( lAlphaXi, lPtArmXi   );
		
	if( lChargeXi < 0 ) {
		f2dHistEffMassLambdaVsEffMassXiMinus->Fill( lInvMassLambdaAsCascDghter, lInvMassXiMinus ); 
		f2dHistEffMassXiVsEffMassOmegaMinus ->Fill( lInvMassXiMinus, lInvMassOmegaMinus );
		f2dHistXiRadiusVsEffMassXiMinus     ->Fill( lXiRadius, lInvMassXiMinus );
	}
	else{
		f2dHistEffMassLambdaVsEffMassXiPlus ->Fill( lInvMassLambdaAsCascDghter, lInvMassXiPlus );
		f2dHistEffMassXiVsEffMassOmegaPlus  ->Fill( lInvMassXiPlus, lInvMassOmegaPlus );
		f2dHistXiRadiusVsEffMassXiPlus->Fill( lXiRadius, lInvMassXiPlus);
	}
		
   
    }// end of the Cascade loop (ESD or AOD)
    
  
  // Post output data.
 PostData(1, fListHistCascade);
}










//________________________________________________________________________
void AliAnalysisTaskCheckCascade::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fHistTrackMultiplicity = dynamic_cast<TH1F*> (   ((TList*)GetOutputData(1))->FindObject("fHistTrackMultiplicity") );
  if (!fHistTrackMultiplicity) {
    Printf("ERROR: fHistTrackMultiplicity not available");
    return;
  }
 
  fHistCascadeMultiplicity = dynamic_cast<TH1F*> ( ((TList*)GetOutputData(1))->FindObject("fHistCascadeMultiplicity"));
  if (!fHistCascadeMultiplicity) {
    Printf("ERROR: fHistCascadeMultiplicity not available");
    return;
  }
   
  TCanvas *c2 = new TCanvas("AliAnalysisTaskCheckCascade","Multiplicity",10,10,510,510);
  c2->cd(1)->SetLogy();

  fHistTrackMultiplicity->SetMarkerStyle(22);
  fHistTrackMultiplicity->DrawCopy("E");
  
  fHistCascadeMultiplicity->SetMarkerStyle(26);
  fHistCascadeMultiplicity->DrawCopy("ESAME");
  
}
