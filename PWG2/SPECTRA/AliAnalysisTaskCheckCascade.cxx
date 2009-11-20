
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
//            (AliAnalysisTaskCheckCascade)
//            This task has four roles :
//              1. QAing the Cascades from ESD and AOD
//                 Origin:  AliAnalysisTaskESDCheckV0 by B.H. Nov2007, hippolyt@in2p3.fr
//              2. Prepare the plots which stand as raw material for yield extraction (wi/wo PID)
//              3. Supply an AliCFContainer meant to define the optimised topological selections
//              4. Rough azimuthal correlation study (Eta, Phi)
//            Adapted to Cascade : A.Maire Mar2008, antonin.maire@ires.in2p3.fr
//            Modified :           A.Maire Nov2009, antonin.maire@ires.in2p3.fr
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

#include <Riostream.h>

#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"


#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
//#include "AliCascadeVertexer.h"
#include "AliTPCpidESD.h"
#include "AliCFContainer.h"
#include "AliMultiplicity.h"

#include "AliESDcascade.h"
#include "AliAODcascade.h"

#include "AliAnalysisTaskCheckCascade.h"

ClassImp(AliAnalysisTaskCheckCascade)



//________________________________________________________________________
AliAnalysisTaskCheckCascade::AliAnalysisTaskCheckCascade() 
  : AliAnalysisTaskSE(), fAnalysisType("ESD"), fCollidingSystems(0), fTpcPidManager(0),

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
    fHistMassWithCombPIDXiMinus(0), fHistMassWithCombPIDXiPlus(0),
    fHistMassWithCombPIDOmegaMinus(0), fHistMassWithCombPIDOmegaPlus(0),

    fHistXiTransvMom(0),    fHistXiTotMom(0),
    fHistBachTransvMom(0),   fHistBachTotMom(0),

    fHistChargeXi(0),
    fHistV0toXiCosineOfPointingAngle(0),

    fHistRapXi(0), fHistRapOmega(0), fHistEta(0),
    fHistTheta(0), fHistPhi(0),

    f2dHistArmenteros(0),			
    f2dHistEffMassLambdaVsEffMassXiMinus(0), f2dHistEffMassXiVsEffMassOmegaMinus(0),
    f2dHistEffMassLambdaVsEffMassXiPlus(0), f2dHistEffMassXiVsEffMassOmegaPlus(0),
    f2dHistXiRadiusVsEffMassXiMinus(0), f2dHistXiRadiusVsEffMassXiPlus(0),
    f2dHistXiRadiusVsEffMassOmegaMinus(0), f2dHistXiRadiusVsEffMassOmegaPlus(0),
    
    f3dHistXiPtVsEffMassVsYXiMinus(0), f3dHistXiPtVsEffMassVsYXiPlus(0),
    f3dHistXiPtVsEffMassVsYOmegaMinus(0), f3dHistXiPtVsEffMassVsYOmegaPlus(0),
    
    f3dHistXiPtVsEffMassVsYWithCombPIDXiMinus(0), f3dHistXiPtVsEffMassVsYWithCombPIDXiPlus(0),
    f3dHistXiPtVsEffMassVsYWithCombPIDOmegaMinus(0),  f3dHistXiPtVsEffMassVsYWithCombPIDOmegaPlus(0),
    f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaMinus(0), f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaPlus(0),
    f3dHistXiPtVsEffMassVsYWithTpcPIDOmegaMinus(0),
    
    fCFContCascadePIDXiMinus(0),
    fCFContCascadePIDXiPlus(0),
    fCFContCascadePIDOmegaMinus(0),
    fCFContCascadePIDOmegaPlus(0),
    fCFContCascadeCuts(0),
    
    fHnSpAngularCorrXiMinus(0), fHnSpAngularCorrXiPlus(0), 
    fHnSpAngularCorrOmegaMinus(0), fHnSpAngularCorrOmegaPlus(0)

{
  // Dummy Constructor
}








//________________________________________________________________________
AliAnalysisTaskCheckCascade::AliAnalysisTaskCheckCascade(const char *name) 
  : AliAnalysisTaskSE(name), fAnalysisType("ESD"), fCollidingSystems(0), fTpcPidManager(0),
     
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
    fHistMassWithCombPIDXiMinus(0), fHistMassWithCombPIDXiPlus(0),
    fHistMassWithCombPIDOmegaMinus(0), fHistMassWithCombPIDOmegaPlus(0),

    fHistXiTransvMom(0),    fHistXiTotMom(0),
    fHistBachTransvMom(0),   fHistBachTotMom(0),

    fHistChargeXi(0),
    fHistV0toXiCosineOfPointingAngle(0),

    fHistRapXi(0), fHistRapOmega(0), fHistEta(0),
    fHistTheta(0), fHistPhi(0),

    f2dHistArmenteros(0),			
    f2dHistEffMassLambdaVsEffMassXiMinus(0), f2dHistEffMassXiVsEffMassOmegaMinus(0),
    f2dHistEffMassLambdaVsEffMassXiPlus(0), f2dHistEffMassXiVsEffMassOmegaPlus(0),
    f2dHistXiRadiusVsEffMassXiMinus(0), f2dHistXiRadiusVsEffMassXiPlus(0),
    f2dHistXiRadiusVsEffMassOmegaMinus(0), f2dHistXiRadiusVsEffMassOmegaPlus(0),
    
    f3dHistXiPtVsEffMassVsYXiMinus(0), f3dHistXiPtVsEffMassVsYXiPlus(0),
    f3dHistXiPtVsEffMassVsYOmegaMinus(0), f3dHistXiPtVsEffMassVsYOmegaPlus(0),
    
    f3dHistXiPtVsEffMassVsYWithCombPIDXiMinus(0), f3dHistXiPtVsEffMassVsYWithCombPIDXiPlus(0),
    f3dHistXiPtVsEffMassVsYWithCombPIDOmegaMinus(0),  f3dHistXiPtVsEffMassVsYWithCombPIDOmegaPlus(0),
    f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaMinus(0), f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaPlus(0),
    f3dHistXiPtVsEffMassVsYWithTpcPIDOmegaMinus(0),
    
    fCFContCascadePIDXiMinus(0),
    fCFContCascadePIDXiPlus(0),
    fCFContCascadePIDOmegaMinus(0),
    fCFContCascadePIDOmegaPlus(0),
    fCFContCascadeCuts(0),
    
    fHnSpAngularCorrXiMinus(0), fHnSpAngularCorrXiPlus(0), 
    fHnSpAngularCorrOmegaMinus(0), fHnSpAngularCorrOmegaPlus(0)

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
			200, 0, 20000); 		
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
			50, 0, 50); 		
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
	fHistXiCosineOfPointingAngle = new TH1F("fHistXiCosineOfPointingAngle", "Cosine of Xi Pointing Angle; Cos (Xi Point.Angl);Number of Xis", 200, 0.99, 1.0);
	fListHistCascade->Add(fHistXiCosineOfPointingAngle);
}

if(! fHistXiRadius ){
	fHistXiRadius  = new TH1F( "fHistXiRadius",  "Casc. decay transv. radius; r (cm); Counts" , 1050, 0., 105.0 );
	fListHistCascade->Add(fHistXiRadius);
}


// - Histos about ~ the "V0 part" of the cascade,  coming by inheritance from AliESDv0



if (! fHistMassLambdaAsCascDghter) {
     fHistMassLambdaAsCascDghter = new TH1F("fHistMassLambdaAsCascDghter","#Lambda associated to Casc. candidates;Eff. Mass (GeV/c^{2});Counts", 300,1.00,1.3);
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
	fHistV0RadiusXi = new TH1F("fHistV0RadiusXi", "V0 decay radius, in cascade; radius (cm); Counts", 1050, 0., 105.0);
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
// By cascade hyp  
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

// By cascade hyp + bachelor PID
if (! fHistMassWithCombPIDXiMinus) {
    fHistMassWithCombPIDXiMinus = new TH1F("fHistMassWithCombPIDXiMinus","#Xi^{-} candidates, with Bach. comb. PID;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 200,1.2,2.0);
    fListHistCascade->Add(fHistMassWithCombPIDXiMinus);
}
  
if (! fHistMassWithCombPIDXiPlus) {
    fHistMassWithCombPIDXiPlus = new TH1F("fHistMassWithCombPIDXiPlus","#Xi^{+} candidates, with Bach. comb. PID;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",200,1.2,2.0);
    fListHistCascade->Add(fHistMassWithCombPIDXiPlus);
}

if (! fHistMassWithCombPIDOmegaMinus) {
	fHistMassWithCombPIDOmegaMinus = new TH1F("fHistMassWithCombPIDOmegaMinus","#Omega^{-} candidates, with Bach. comb. PID;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 250,1.5,2.5);
    fListHistCascade->Add(fHistMassWithCombPIDOmegaMinus);
}
 
if (! fHistMassWithCombPIDOmegaPlus) {
	fHistMassWithCombPIDOmegaPlus = new TH1F("fHistMassWithCombPIDOmegaPlus","#Omega^{+} candidates, with Bach. comb. PID;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 250,1.5,2.5);
    fListHistCascade->Add(fHistMassWithCombPIDOmegaPlus);
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

//-------

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

//-------

if(! f2dHistXiRadiusVsEffMassXiMinus) {
	f2dHistXiRadiusVsEffMassXiMinus = new TH2F( "f2dHistXiRadiusVsEffMassXiMinus", "Transv. R_{Xi Decay} Vs M_{#Xi^{-} candidates}; r_{cascade} (cm); M( #Lambda , #pi^{-} ) (GeV/c^{2}) ", 450, 0., 45.0, 200, 1.2, 2.0);
	fListHistCascade->Add(f2dHistXiRadiusVsEffMassXiMinus);
}

if(! f2dHistXiRadiusVsEffMassXiPlus) {
	f2dHistXiRadiusVsEffMassXiPlus = new TH2F( "f2dHistXiRadiusVsEffMassXiPlus", "Transv. R_{Xi Decay} Vs M_{#Xi^{+} candidates}; r_{cascade} (cm); M( #Lambda , #pi^{+} ) (GeV/c^{2}) ", 450, 0., 45.0, 200, 1.2, 2.0);
	fListHistCascade->Add(f2dHistXiRadiusVsEffMassXiPlus);
}

if(! f2dHistXiRadiusVsEffMassOmegaMinus) {
	f2dHistXiRadiusVsEffMassOmegaMinus = new TH2F( "f2dHistXiRadiusVsEffMassOmegaMinus", "Transv. R_{Xi Decay} Vs M_{#Omega^{-} candidates}; r_{cascade} (cm); M( #Lambda , K^{-} ) (GeV/c^{2}) ", 450, 0., 45.0, 250, 1.5, 2.5);
	fListHistCascade->Add(f2dHistXiRadiusVsEffMassOmegaMinus);
}

if(! f2dHistXiRadiusVsEffMassOmegaPlus) {
	f2dHistXiRadiusVsEffMassOmegaPlus = new TH2F( "f2dHistXiRadiusVsEffMassOmegaPlus", "Transv. R_{Xi Decay} Vs M_{#Omega^{+} candidates}; r_{cascade} (cm); M( #Lambda , K^{+} ) (GeV/c^{2}) ", 450, 0., 45.0, 250, 1.5, 2.5);
	fListHistCascade->Add(f2dHistXiRadiusVsEffMassOmegaPlus);
}


// Part 2 : Raw material for yield extraction -------

if(! f3dHistXiPtVsEffMassVsYXiMinus) {
	f3dHistXiPtVsEffMassVsYXiMinus = new TH3F( "f3dHistXiPtVsEffMassVsYXiMinus", "Pt_{cascade} Vs M_{#Xi^{-} candidates} Vs Y_{#Xi}; Pt_{cascade} (GeV/c); M( #Lambda , #pi^{-} ) (GeV/c^{2}) ;Y_{#Xi} ", 100, 0., 10.0, 400, 1.2, 2.0, 48, -1.2,1.2);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYXiMinus);
}

if(! f3dHistXiPtVsEffMassVsYXiPlus) {
	f3dHistXiPtVsEffMassVsYXiPlus = new TH3F( "f3dHistXiPtVsEffMassVsYXiPlus", "Pt_{cascade} Vs M_{#Xi^{+} candidates} Vs Y_{#Xi}; Pt_{cascade} (GeV/c); M( #Lambda , #pi^{+} ) (GeV/c^{2}); Y_{#Xi}", 100, 0., 10.0, 400, 1.2, 2.0, 48, -1.2,1.2);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYXiPlus);
}

if(! f3dHistXiPtVsEffMassVsYOmegaMinus) {
	f3dHistXiPtVsEffMassVsYOmegaMinus = new TH3F( "f3dHistXiPtVsEffMassVsYOmegaMinus", "Pt_{cascade} Vs M_{#Omega^{-} candidates} Vs Y_{#Omega}; Pt_{cascade} (GeV/c); M( #Lambda , K^{-} ) (GeV/c^{2}); Y_{#Omega}", 100, 0., 10.0, 500, 1.5, 2.5, 48, -1.2,1.2);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYOmegaMinus);
}

if(! f3dHistXiPtVsEffMassVsYOmegaPlus) {
	f3dHistXiPtVsEffMassVsYOmegaPlus = new TH3F( "f3dHistXiPtVsEffMassVsYOmegaPlus", "Pt_{cascade} Vs M_{#Omega^{+} candidates} Vs Y_{#Omega}; Pt_{cascade} (GeV/c); M( #Lambda , K^{+} ) (GeV/c^{2}); Y_{#Omega}", 100, 0., 10.0, 500, 1.5, 2.5, 48, -1.2,1.2);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYOmegaPlus);
}

//--
if(! f3dHistXiPtVsEffMassVsYWithCombPIDXiMinus) {
	f3dHistXiPtVsEffMassVsYWithCombPIDXiMinus = new TH3F( "f3dHistXiPtVsEffMassVsYWithCombPIDXiMinus", "Pt_{cascade} Vs M_{#Xi^{-} candidates, with PID} Vs Y_{#Xi}; Pt_{cascade} (GeV/c); M( #Lambda , #pi^{-} ) (GeV/c^{2}) ;Y_{#Xi} ", 100, 0., 10.0, 400, 1.2, 2.0, 48, -1.2,1.2);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYWithCombPIDXiMinus);
}

if(! f3dHistXiPtVsEffMassVsYWithCombPIDXiPlus) {
	f3dHistXiPtVsEffMassVsYWithCombPIDXiPlus = new TH3F( "f3dHistXiPtVsEffMassVsYWithCombPIDXiPlus", "Pt_{cascade} Vs M_{#Xi^{+} candidates, with PID} Vs Y_{#Xi}; Pt_{cascade} (GeV/c); M( #Lambda , #pi^{+} ) (GeV/c^{2}); Y_{#Xi}", 100, 0., 10.0, 400, 1.2, 2.0, 48, -1.2,1.2);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYWithCombPIDXiPlus);
}

if(! f3dHistXiPtVsEffMassVsYWithCombPIDOmegaMinus) {
	f3dHistXiPtVsEffMassVsYWithCombPIDOmegaMinus = new TH3F( "f3dHistXiPtVsEffMassVsYWithCombPIDOmegaMinus", "Pt_{cascade} Vs M_{#Omega^{-} candidates, with PID} Vs Y_{#Omega}; Pt_{cascade} (GeV/c); M( #Lambda , K^{-} ) (GeV/c^{2}); Y_{#Omega}", 100, 0., 10.0, 500, 1.5, 2.5, 48, -1.2,1.2);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYWithCombPIDOmegaMinus);
}

if(! f3dHistXiPtVsEffMassVsYWithCombPIDOmegaPlus) {
	f3dHistXiPtVsEffMassVsYWithCombPIDOmegaPlus = new TH3F( "f3dHistXiPtVsEffMassVsYWithCombPIDOmegaPlus", "Pt_{cascade} Vs M_{#Omega^{+} candidates, with PID} Vs Y_{#Omega}; Pt_{cascade} (GeV/c); M( #Lambda , K^{+} ) (GeV/c^{2}); Y_{#Omega}", 100, 0., 10.0, 500, 1.5, 2.5, 48, -1.2,1.2);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYWithCombPIDOmegaPlus);
}

//--
if(! f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaMinus) {
	f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaMinus = new TH3F( "f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaMinus", "Pt_{cascade} Vs M_{#Omega^{-} candidates, with 2 PID} Vs Y_{#Omega}; Pt_{cascade} (GeV/c); M( #Lambda , K^{-} ) (GeV/c^{2}); Y_{#Omega}", 100, 0., 10.0, 500, 1.5, 2.5, 48, -1.2,1.2);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaMinus);
}

if(! f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaPlus) {
	f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaPlus = new TH3F( "f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaPlus", "Pt_{cascade} Vs M_{#Omega^{+} candidates, with 2 PID} Vs Y_{#Omega}; Pt_{cascade} (GeV/c); M( #Lambda , K^{+} ) (GeV/c^{2}); Y_{#Omega}", 100, 0., 10.0, 500, 1.5, 2.5, 48, -1.2,1.2);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaPlus);
}


if(! fTpcPidManager){
		
  Double_t lAlephParameters[5] = {0.};
  	// Reasonable parameters extracted for p-p simulation (LHC09a4) - A.Kalweit
	lAlephParameters[0] = 4.23232575531564326e+00;//50*0.76176e-1;
	lAlephParameters[1] = 8.68482806165147636e+00;//10.632; 
	lAlephParameters[2] = 1.34000000000000005e-05;//0.13279e-4;
	lAlephParameters[3] = 2.30445734159456084e+00;//1.8631;
	lAlephParameters[4] = 2.25624744086878559e+00;//1.9479;

  fTpcPidManager = new AliTPCpidESD();
  fTpcPidManager->SetBetheBlochParameters(lAlephParameters[0]/50.,
					  lAlephParameters[1],
					  lAlephParameters[2],
					  lAlephParameters[3],
					  lAlephParameters[4]);
}


if(! f3dHistXiPtVsEffMassVsYWithTpcPIDOmegaMinus) {
	f3dHistXiPtVsEffMassVsYWithTpcPIDOmegaMinus = new TH3F( "f3dHistXiPtVsEffMassVsYWithTpcPIDOmegaMinus", "Pt_{cascade} Vs M_{#Omega^{-} candidates, with TPC PID} Vs Y_{#Omega}; Pt_{cascade} (GeV/c); M( #Lambda , K^{-} ) (GeV/c^{2}); Y_{#Omega}", 100, 0., 10.0, 500, 1.5, 2.5, 48, -1.2,1.2);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYWithTpcPIDOmegaMinus);
}



if(! fCFContCascadePIDXiMinus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4];
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 400;
  lNbBinsPerVar[2] = 48;
  lNbBinsPerVar[3] = 200;
   
  
  fCFContCascadePIDXiMinus = new AliCFContainer("fCFContCascadePIDXiMinus","Pt_{cascade} Vs M_{#Xi^{-} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDXiMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDXiMinus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
  fCFContCascadePIDXiMinus->SetBinLimits(2,  -1.2  ,   1.2 );	// Rapidity
  if(fCollidingSystems) 
	fCFContCascadePIDXiMinus->SetBinLimits(3, 0.0, 20000.0  );    // TrackMultiplicity
  else
	fCFContCascadePIDXiMinus->SetBinLimits(3, 0.0, 200.0  );     // TrackMultiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDXiMinus->SetStepTitle(0, "No PID");
  fCFContCascadePIDXiMinus->SetStepTitle(1, "TPC PID / 3-#sigma cut on Bachelor track");
  fCFContCascadePIDXiMinus->SetStepTitle(2, "TPC PID / 3-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDXiMinus->SetStepTitle(3, "TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDXiMinus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDXiMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDXiMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDXiMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDXiMinus->SetVarTitle(1, "M( #Lambda , #pi^{-} ) (GeV/c^{2})");
  fCFContCascadePIDXiMinus->SetVarTitle(2, "Y_{#Xi}");
  fCFContCascadePIDXiMinus->SetVarTitle(3, "Track Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDXiMinus);
  
}

if(! fCFContCascadePIDXiPlus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4];
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 400;
  lNbBinsPerVar[2] = 48;
  lNbBinsPerVar[3] = 200;
   
  
  fCFContCascadePIDXiPlus = new AliCFContainer("fCFContCascadePIDXiPlus","Pt_{cascade} Vs M_{#Xi^{+} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDXiPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDXiPlus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
  fCFContCascadePIDXiPlus->SetBinLimits(2,  -1.2  ,   1.2 );	// Rapidity
  if(fCollidingSystems) 
	fCFContCascadePIDXiPlus->SetBinLimits(3, 0.0, 20000.0  );    // TrackMultiplicity
  else
	fCFContCascadePIDXiPlus->SetBinLimits(3, 0.0, 200.0  );     // TrackMultiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDXiPlus->SetStepTitle(0, "No PID");
  fCFContCascadePIDXiPlus->SetStepTitle(1, "TPC PID / 3-#sigma cut on Bachelor track");
  fCFContCascadePIDXiPlus->SetStepTitle(2, "TPC PID / 3-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDXiPlus->SetStepTitle(3, "TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDXiPlus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDXiPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDXiPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDXiPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDXiPlus->SetVarTitle(1, "M( #Lambda , #pi^{+} ) (GeV/c^{2})");
  fCFContCascadePIDXiPlus->SetVarTitle(2, "Y_{#Xi}");
  fCFContCascadePIDXiPlus->SetVarTitle(3, "Track Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDXiPlus);
  
}


if(! fCFContCascadePIDOmegaMinus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4];
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 500;
  lNbBinsPerVar[2] = 48;
  lNbBinsPerVar[3] = 200;
   
  
  fCFContCascadePIDOmegaMinus = new AliCFContainer("fCFContCascadePIDOmegaMinus","Pt_{cascade} Vs M_{#Omega^{-} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDOmegaMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDOmegaMinus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
  fCFContCascadePIDOmegaMinus->SetBinLimits(2,  -1.2  ,   1.2 );	// Rapidity
  if(fCollidingSystems) 
	fCFContCascadePIDOmegaMinus->SetBinLimits(3, 0.0, 20000.0  );    // TrackMultiplicity
  else
	fCFContCascadePIDOmegaMinus->SetBinLimits(3, 0.0, 200.0  );     // TrackMultiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDOmegaMinus->SetStepTitle(0, "No PID");
  fCFContCascadePIDOmegaMinus->SetStepTitle(1, "TPC PID / 3-#sigma cut on Bachelor track");
  fCFContCascadePIDOmegaMinus->SetStepTitle(2, "TPC PID / 3-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDOmegaMinus->SetStepTitle(3, "TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDOmegaMinus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDOmegaMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDOmegaMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDOmegaMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDOmegaMinus->SetVarTitle(1, "M( #Lambda , K^{-} ) (GeV/c^{2})");
  fCFContCascadePIDOmegaMinus->SetVarTitle(2, "Y_{#Omega}");
  fCFContCascadePIDOmegaMinus->SetVarTitle(3, "Track Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDOmegaMinus);
  
}

if(! fCFContCascadePIDOmegaPlus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4];
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 500;
  lNbBinsPerVar[2] = 48;
  lNbBinsPerVar[3] = 200;
   
  
  fCFContCascadePIDOmegaPlus = new AliCFContainer("fCFContCascadePIDOmegaPlus","Pt_{cascade} Vs M_{#Omega^{+} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDOmegaPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDOmegaPlus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
  fCFContCascadePIDOmegaPlus->SetBinLimits(2,  -1.2  ,   1.2 );	// Rapidity
  if(fCollidingSystems) 
	fCFContCascadePIDOmegaPlus->SetBinLimits(3, 0.0, 20000.0  );    // TrackMultiplicity
  else
	fCFContCascadePIDOmegaPlus->SetBinLimits(3, 0.0, 200.0  );     // TrackMultiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDOmegaPlus->SetStepTitle(0, "No PID");
  fCFContCascadePIDOmegaPlus->SetStepTitle(1, "TPC PID / 3-#sigma cut on Bachelor track");
  fCFContCascadePIDOmegaPlus->SetStepTitle(2, "TPC PID / 3-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDOmegaPlus->SetStepTitle(3, "TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDOmegaPlus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDOmegaPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDOmegaPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDOmegaPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDOmegaPlus->SetVarTitle(1, "M( #Lambda , K^{+} ) (GeV/c^{2})");
  fCFContCascadePIDOmegaPlus->SetVarTitle(2, "Y_{#Omega}");
  fCFContCascadePIDOmegaPlus->SetVarTitle(3, "Track Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDOmegaPlus);
  
}











// Part 3 : Towards the optimisation of topological selections -------
if(! fCFContCascadeCuts){
   
	// Container meant to store all the relevant distributions corresponding to the cut variables.
	// So far, 19 variables have been identified.
	// The following will be done in quite an inelegant way.
	// Improvement expected later.
  const	Int_t  lNbSteps      =  2 ;
  const Int_t  lNbVariables  =  18 ;
  
  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[18];
  lNbBinsPerVar[0]  = 25;
  lNbBinsPerVar[1]  = 25;
  lNbBinsPerVar[2]  = 20;
  lNbBinsPerVar[3]  = 40;
  lNbBinsPerVar[4]  = 50;
  lNbBinsPerVar[5]  = 12;
  
  lNbBinsPerVar[6]  = 20;
  lNbBinsPerVar[7]  = 40;
  lNbBinsPerVar[8]  = 40;
  lNbBinsPerVar[9]  = 25;
  lNbBinsPerVar[10] = 25;
  
  lNbBinsPerVar[11] = 60;
  lNbBinsPerVar[12] = 50;
  
  lNbBinsPerVar[13] = 20;
  lNbBinsPerVar[14] = 20;
 
  lNbBinsPerVar[15] = 40;
  lNbBinsPerVar[16] = 40;
  lNbBinsPerVar[17] = 35;
    
 fCFContCascadeCuts = new AliCFContainer("fCFContCascadeCuts","Container for Cascade cuts", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid for v4-18-10-AN on)
  fCFContCascadeCuts->SetBinLimits(0,    0.0  ,  0.25 );	// DcaXiDaughters
  fCFContCascadeCuts->SetBinLimits(1,    0.0  ,  0.25 );	// DcaBachToPrimVertexXi
  fCFContCascadeCuts->SetBinLimits(2,    0.995,  1.0  );	// XiCosineOfPointingAngle
  fCFContCascadeCuts->SetBinLimits(3,    0.0  ,  4.0  );	// XiRadius
  fCFContCascadeCuts->SetBinLimits(4,    1.1  ,  1.15  );	// InvMassLambdaAsCascDghter
  fCFContCascadeCuts->SetBinLimits(5,    0.0  ,  0.6  );	// DcaV0DaughtersXi
  fCFContCascadeCuts->SetBinLimits(6,    0.98 ,  1.0  );	// V0CosineOfPointingAngleXi
  fCFContCascadeCuts->SetBinLimits(7,    0.0  , 20.0  );	// V0RadiusXi
  fCFContCascadeCuts->SetBinLimits(8,    0.0  ,  1.0  );	// DcaV0ToPrimVertexXi
  fCFContCascadeCuts->SetBinLimits(9,    0.0  ,  2.5  );	// DcaPosToPrimVertexXi
  fCFContCascadeCuts->SetBinLimits(10,   0.0  ,  2.5  );	// DcaNegToPrimVertexXi
  fCFContCascadeCuts->SetBinLimits(11,   1.25 ,  1.45  );	// InvMassXi
  fCFContCascadeCuts->SetBinLimits(12,   1.6  ,  1.8  );	// InvMassOmega
  fCFContCascadeCuts->SetBinLimits(13,   0.0  , 10.0  );	// XiTransvMom
  fCFContCascadeCuts->SetBinLimits(14, -10.0  , 10.0  );	// BestPrimaryVtxPosZ
  if(fCollidingSystems){
  	fCFContCascadeCuts->SetBinLimits(15,   0.0, 10000.0  );    // TrackMultiplicity
  	fCFContCascadeCuts->SetBinLimits(16,   0.0, 10000.0  );    // SPDTrackletsMultiplicity
  }
  else{  
  	fCFContCascadeCuts->SetBinLimits(15,   0.0, 200.0  );     // TrackMultiplicity
  	fCFContCascadeCuts->SetBinLimits(16,   0.0, 200.0  );     // SPDTrackletsMultiplicity
  }
  fCFContCascadeCuts->SetBinLimits(17,  25.0  ,165.0  );	// BachTPCClusters
  
  
  
  // Setting the number of steps : one for negative cascades (Xi- and Omega-), another for the positve cascades(Xi+ and Omega+)
  fCFContCascadeCuts->SetStepTitle(0, "Negative Cascades");
  fCFContCascadeCuts->SetStepTitle(1, "Positive Cascades");
  
  // Setting the variable title, per axis
  // fCFContCascadeCuts->SetVarTitle(0,  "Chi2Xi");
  fCFContCascadeCuts->SetVarTitle(0,  "DcaXiDaughters");
  fCFContCascadeCuts->SetVarTitle(1,  "DcaBachToPrimVertexXi");
  fCFContCascadeCuts->SetVarTitle(2,  "XiCosineOfPointingAngle");
  fCFContCascadeCuts->SetVarTitle(3,  "XiRadius");
  fCFContCascadeCuts->SetVarTitle(4,  "InvMassLambdaAsCascDghter");
   // fCFContCascadeCuts->SetVarTitle(0,  "V0Chi2Xi");
  fCFContCascadeCuts->SetVarTitle(5,  "DcaV0DaughtersXi");
  
  fCFContCascadeCuts->SetVarTitle(6,  "V0CosineOfPointingAngleXi");
  fCFContCascadeCuts->SetVarTitle(7,  "V0RadiusXi");
  fCFContCascadeCuts->SetVarTitle(8,  "DcaV0ToPrimVertexXi");
  fCFContCascadeCuts->SetVarTitle(9,  "DcaPosToPrimVertexXi");
  fCFContCascadeCuts->SetVarTitle(10, "DcaNegToPrimVertexXi");
  
  fCFContCascadeCuts->SetVarTitle(11, "InvMassXi");
  fCFContCascadeCuts->SetVarTitle(12, "InvMassOmega");
  
  fCFContCascadeCuts->SetVarTitle(13, "XiTransvMom");
  //fCFContCascadeCuts->SetVarTitle(14, "V0toXiCosineOfPointingAngle");
  fCFContCascadeCuts->SetVarTitle(14, "BestPrimaryVtxPosZ");
  
  fCFContCascadeCuts->SetVarTitle(15, "TrackMultiplicity");
  fCFContCascadeCuts->SetVarTitle(16, "SPDTrackletsMultiplicity");
  fCFContCascadeCuts->SetVarTitle(17, "BachTPCClusters");
  
  fListHistCascade->Add(fCFContCascadeCuts);
}



// Part 4 : Angular correlation study -------

if(! fHnSpAngularCorrXiMinus){
	// Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks
	// Delta Phi  = 360 bins de -180., 180.
	// Delta Eta  = 120 bins de -3.0, 3.0
	// Pt Cascade = 100 bins de 0., 10.0,
	// Pt track = 150 bins de 0., 15.0
	
   Int_t    bins[5] = { 360, 120, 100, 150, 40};
   Double_t xmin[5] = {-50., -3.,  0.,  0., 1.30};
   Double_t xmax[5] = { 310., 3., 10., 15., 1.34};
   fHnSpAngularCorrXiMinus = new THnSparseF("fHnSpAngularCorrXiMinus", "Angular Correlation for #Xi^{-}:", 5, bins, xmin, xmax);
	fHnSpAngularCorrXiMinus->GetAxis(0)->SetTitle(" #Delta#phi(Casc,Track) (deg)");
	fHnSpAngularCorrXiMinus->GetAxis(1)->SetTitle(" #Delta#eta(Casc,Track)");
	fHnSpAngularCorrXiMinus->GetAxis(2)->SetTitle(" Pt_{Casc} (GeV/c)");
	fHnSpAngularCorrXiMinus->GetAxis(3)->SetTitle(" Pt_{any track} (GeV/c)");
	fHnSpAngularCorrXiMinus->GetAxis(4)->SetTitle(" Eff. Inv Mass (GeV/c^{2})");
	fHnSpAngularCorrXiMinus->Sumw2();
   fListHistCascade->Add(fHnSpAngularCorrXiMinus);
}

if(! fHnSpAngularCorrXiPlus){
	// Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks
	// Delta Phi  = 360 bins de -180., 180.
	// Delta Eta  = 120 bins de -3.0, 3.0
	// Pt Cascade = 100 bins de 0., 10.0,
	// Pt track = 150 bins de 0., 15.0
   Int_t    bins[5] = { 360, 120, 100, 150, 40};
   Double_t xmin[5] = {-50., -3.,  0.,  0., 1.30};
   Double_t xmax[5] = { 310., 3., 10., 15., 1.34};
   fHnSpAngularCorrXiPlus = new THnSparseF("fHnSpAngularCorrXiPlus", "Angular Correlation for #Xi^{+}:", 5, bins, xmin, xmax);
	fHnSpAngularCorrXiPlus->GetAxis(0)->SetTitle(" #Delta#phi(Casc,Track) (deg)");
	fHnSpAngularCorrXiPlus->GetAxis(1)->SetTitle(" #Delta#eta(Casc,Track)");
	fHnSpAngularCorrXiPlus->GetAxis(2)->SetTitle(" Pt_{Casc} (GeV/c)");
	fHnSpAngularCorrXiPlus->GetAxis(3)->SetTitle(" Pt_{any track} (GeV/c)");
	fHnSpAngularCorrXiPlus->GetAxis(4)->SetTitle(" Eff. Inv Mass (GeV/c^{2})");
	fHnSpAngularCorrXiPlus->Sumw2();
   fListHistCascade->Add(fHnSpAngularCorrXiPlus);
}

if(! fHnSpAngularCorrOmegaMinus){
	// Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks
	// Delta Phi  = 360 bins de -180., 180.
	// Delta Eta  = 120 bins de -3.0, 3.0
	// Pt Cascade = 100 bins de 0., 10.0,
	// Pt track = 150 bins de 0., 15.0
	
   Int_t    bins[5] = { 360, 120, 100, 150, 40};
   Double_t xmin[5] = {-50., -3.,  0.,  0., 1.65};
   Double_t xmax[5] = { 310., 3., 10., 15., 1.69};
   fHnSpAngularCorrOmegaMinus = new THnSparseF("fHnSpAngularCorrOmegaMinus", "Angular Correlation for #Omega^{-}:", 5, bins, xmin, xmax);
	fHnSpAngularCorrOmegaMinus->GetAxis(0)->SetTitle(" #Delta#phi(Casc,Track) (deg)");
	fHnSpAngularCorrOmegaMinus->GetAxis(1)->SetTitle(" #Delta#eta(Casc,Track)");
	fHnSpAngularCorrOmegaMinus->GetAxis(2)->SetTitle(" Pt_{Casc} (GeV/c)");
	fHnSpAngularCorrOmegaMinus->GetAxis(3)->SetTitle(" Pt_{any track} (GeV/c)");
	fHnSpAngularCorrOmegaMinus->GetAxis(4)->SetTitle(" Eff. Inv Mass (GeV/c^{2})");
	fHnSpAngularCorrOmegaMinus->Sumw2();
   fListHistCascade->Add(fHnSpAngularCorrOmegaMinus);
}

if(! fHnSpAngularCorrOmegaPlus){
	// Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks
	// Delta Phi  = 360 bins de -180., 180.
	// Delta Eta  = 120 bins de -3.0, 3.0
	// Pt Cascade = 100 bins de 0., 10.0,
	// Pt track = 150 bins de 0., 15.0
   Int_t    bins[5] = { 360, 120, 100, 150, 40};
   Double_t xmin[5] = {-50., -3.,  0.,  0., 1.65};
   Double_t xmax[5] = { 310., 3., 10., 15., 1.69};
   fHnSpAngularCorrOmegaPlus = new THnSparseF("fHnSpAngularCorrOmegaPlus", "Angular Correlation for #Omega^{+}:", 5, bins, xmin, xmax);
	fHnSpAngularCorrOmegaPlus->GetAxis(0)->SetTitle(" #Delta#phi(Casc,Track) (deg)");
	fHnSpAngularCorrOmegaPlus->GetAxis(1)->SetTitle(" #Delta#eta(Casc,Track)");
	fHnSpAngularCorrOmegaPlus->GetAxis(2)->SetTitle(" Pt_{Casc} (GeV/c)");
	fHnSpAngularCorrOmegaPlus->GetAxis(3)->SetTitle(" Pt_{any track} (GeV/c)");
	fHnSpAngularCorrOmegaPlus->GetAxis(4)->SetTitle(" Eff. Inv Mass (GeV/c^{2})");
	fHnSpAngularCorrOmegaPlus->Sumw2();
   fListHistCascade->Add(fHnSpAngularCorrOmegaPlus);
}


}// end UserCreateOutputObjects






//________________________________________________________________________
void AliAnalysisTaskCheckCascade::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
	
	AliESDEvent *lESDevent = 0x0;
	AliAODEvent *lAODevent = 0x0;
	Int_t ncascades          = -1;
	Int_t nTrackMultiplicity = -1;
	
	
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
	// Printf("Number of cascade(s) = %d \n", ncascades);
  }
  
        nTrackMultiplicity = (InputEvent())->GetNumberOfTracks();
	
  
  //-------------------------------------------------
  // O - Cascade vertexer (ESD)

	//   if(fAnalysisType == "ESD" ){
	//      lESDevent->ResetCascades();
	//  	AliCascadeVertexer CascVtxer;
	//  	CascVtxer.V0sTracks2CascadeVertices(lESDevent);
	//     }
  
  
  // ---------------------------------------------------------------
  // I - General histos (filled for any event)

	fHistTrackMultiplicity  ->Fill( nTrackMultiplicity );
	fHistCascadeMultiplicity->Fill( ncascades );
  
  
  // ---------------------------------------------------------------
  // II - Calcultaion Part dedicated to Xi vertices
  
  for (Int_t iXi = 0; iXi < ncascades; iXi++)
  {// This is the begining of the Cascade loop (ESD or AOD)
	   
    // -------------------------------------
    // II.Init - Initialisation of the local variables that will be needed for ESD/AOD

  
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
  
	// - 4th part of initialisation : PID treatment
	Bool_t   lIsPosInXiProton      = kFALSE;
	Bool_t   lIsPosInXiPion        = kFALSE;
	Bool_t   lIsPosInOmegaProton   = kFALSE;
	Bool_t   lIsPosInOmegaPion     = kFALSE;
			
	Bool_t   lIsNegInXiProton      = kFALSE;
	Bool_t   lIsNegInXiPion        = kFALSE;
	Bool_t   lIsNegInOmegaProton   = kFALSE;
	Bool_t   lIsNegInOmegaPion     = kFALSE;
	
	Bool_t   lIsBachelorKaon       = kFALSE;
	Bool_t   lIsBachelorPion       = kFALSE; 
	
	Bool_t   lIsBachelorKaonForTPC = kFALSE; // For ESD only ...//FIXME : wait for availability in AOD
	Bool_t   lIsBachelorPionForTPC = kFALSE; // For ESD only ...
	Bool_t   lIsNegPionForTPC      = kFALSE; // For ESD only ...
	Bool_t   lIsPosPionForTPC      = kFALSE; // For ESD only ...
	Bool_t   lIsNegProtonForTPC    = kFALSE; // For ESD only ...
	Bool_t   lIsPosProtonForTPC    = kFALSE; // For ESD only ...

	// - 5th part of initialisation : extra info for QA
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
	
  	// - 6th part of initialisation : variables for the AliCFContainer dedicated to cascade cut optmisiation
	Int_t    lSPDTrackletsMultiplicity = -1;
	Int_t    lBachTPCClusters          = -1;
  
  	// - 7th part of initialisation : variables needed for Angular correlations
	TVector3 lTVect3MomXi(0.,0.,0.);
	Int_t    lArrTrackID[3] = {-1, -1, -1};

	  
  if(fAnalysisType == "ESD"){ 
  
  // -------------------------------------
  // II.ESD - Calcultaion Part dedicated to Xi vertices (ESD)
  
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
		Printf("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
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
	
	
		// - II.Step 3' : extra-selection for cascade candidates
		// Towards optimisation of AA selection
	Bool_t kExtraSelections = kFALSE;
	
	if(kExtraSelections){
		// if(lChi2Xi > 2000) continue;
		// if(lV0Chi2Xi > 2000) continue;
		
		if(lDcaXiDaughters > 0.05) continue; // > 0.1 by default
		//if(lXiCosineOfPointingAngle < 0.999 ) continue;
		if(lXiRadius < 1.0) continue;
		if(lXiRadius > 100) continue;
		if(TMath::Abs(lInvMassLambdaAsCascDghter-1.11568) > 0.007) continue;
		if(lDcaV0DaughtersXi > 0.3) continue;
		
		if(lV0CosineOfPointingAngleXi > 0.9999) continue;
		//if(lDcaV0ToPrimVertexXi < 0.09) continue;
		if(lDcaBachToPrimVertexXi < 0.04) continue;
		
		if(lV0RadiusXi < 1.0) continue;
		if(lV0RadiusXi > 100) continue;
		//if(lDcaPosToPrimVertexXi < 0.6) continue;
		//if(lDcaNegToPrimVertexXi < 0.6) continue;
	}
	
	
	
		// - II.Step 4 : around effective masses (ESD)
		// ~ change mass hypotheses to cover all the possibilities :  Xi-/+, Omega -/+
		//-------------

	
	if( bachTrackXi->Charge() < 0 )	{
		lV0quality = 0.;
		xi->ChangeMassHypothesis(lV0quality , 3312); 	
			// Calculate the effective mass of the Xi- candidate. 
			// pdg code 3312 = Xi-
		lInvMassXiMinus = xi->GetEffMassXi();
		
		lV0quality = 0.;
		xi->ChangeMassHypothesis(lV0quality , 3334); 	
			// Calculate the effective mass of the Xi- candidate. 
			// pdg code 3334 = Omega-
		lInvMassOmegaMinus = xi->GetEffMassXi();
					
		lV0quality = 0.;
		xi->ChangeMassHypothesis(lV0quality , 3312); 	// Back to default hyp.
	}// end if negative bachelor
	
	
	if( bachTrackXi->Charge() >  0 ){
		lV0quality = 0.;
		xi->ChangeMassHypothesis(lV0quality , -3312); 	
			// Calculate the effective mass of the Xi+ candidate. 
			// pdg code -3312 = Xi+
		lInvMassXiPlus = xi->GetEffMassXi();
		
		lV0quality = 0.;
		xi->ChangeMassHypothesis(lV0quality , -3334); 	
			// Calculate the effective mass of the Xi+ candidate. 
			// pdg code -3334  = Omega+
		lInvMassOmegaPlus = xi->GetEffMassXi();
		
		lV0quality = 0.;
		xi->ChangeMassHypothesis(lV0quality , -3312); 	// Back to "default" hyp.
	}// end if positive bachelor
	
	
	
		// - II.Step 5 : PID on the daughter tracks
		//-------------
	
	// A - Combined PID
	// Reasonable guess for the priors for the cascade track sample (e-, mu, pi, K, p)
	Double_t lPriorsGuessXi[5]    = {0, 0, 2, 0, 1};
	Double_t lPriorsGuessOmega[5] = {0, 0, 1, 1, 1};
	
	// Combined VO-positive-daughter PID
	AliPID pPidXi;		pPidXi.SetPriors(    lPriorsGuessXi    );
	AliPID pPidOmega;	pPidOmega.SetPriors( lPriorsGuessOmega );
		
	if( pTrackXi->IsOn(AliESDtrack::kESDpid) ){  // Combined PID exists
		Double_t r[10] = {0.}; pTrackXi->GetESDpid(r);
		pPidXi.SetProbabilities(r);
		pPidOmega.SetProbabilities(r);
		
		// Check if the V0 positive track is a proton (case for Xi-)
		Double_t pproton = pPidXi.GetProbability(AliPID::kProton);
		if (pproton > pPidXi.GetProbability(AliPID::kElectron) &&
		    pproton > pPidXi.GetProbability(AliPID::kMuon)     &&
		    pproton > pPidXi.GetProbability(AliPID::kPion)     &&
		    pproton > pPidXi.GetProbability(AliPID::kKaon)     )     lIsPosInXiProton = kTRUE;
		
		// Check if the V0 positive track is a pi+ (case for Xi+)
		Double_t ppion = pPidXi.GetProbability(AliPID::kPion);
		if (ppion > pPidXi.GetProbability(AliPID::kElectron) &&
		    ppion > pPidXi.GetProbability(AliPID::kMuon)     &&
		    ppion > pPidXi.GetProbability(AliPID::kKaon)     &&
		    ppion > pPidXi.GetProbability(AliPID::kProton)   )     lIsPosInXiPion = kTRUE;
		
		
		// Check if the V0 positive track is a proton (case for Omega-)
		pproton = 0.;
		    pproton = pPidOmega.GetProbability(AliPID::kProton);
		if (pproton > pPidOmega.GetProbability(AliPID::kElectron) &&
		    pproton > pPidOmega.GetProbability(AliPID::kMuon)     &&
		    pproton > pPidOmega.GetProbability(AliPID::kPion)     &&
		    pproton > pPidOmega.GetProbability(AliPID::kKaon)     )  lIsPosInOmegaProton = kTRUE;
		
		// Check if the V0 positive track is a pi+ (case for Omega+)
		ppion = 0.;
		    ppion = pPidOmega.GetProbability(AliPID::kPion);
		if (ppion > pPidOmega.GetProbability(AliPID::kElectron) &&
		    ppion > pPidOmega.GetProbability(AliPID::kMuon)     &&
		    ppion > pPidOmega.GetProbability(AliPID::kKaon)     &&
		    ppion > pPidOmega.GetProbability(AliPID::kProton)   )    lIsPosInOmegaPion = kTRUE;
		
	}// end if V0 positive track with existing combined PID	
	
	
	// Combined VO-negative-daughter PID
	AliPID nPidXi;		nPidXi.SetPriors(    lPriorsGuessXi    );
	AliPID nPidOmega;	nPidOmega.SetPriors( lPriorsGuessOmega );
		
	if( nTrackXi->IsOn(AliESDtrack::kESDpid) ){  // Combined PID exists
		Double_t r[10] = {0.}; nTrackXi->GetESDpid(r);
		nPidXi.SetProbabilities(r);
		nPidOmega.SetProbabilities(r);
		
		// Check if the V0 negative track is a pi- (case for Xi-)
		Double_t ppion = nPidXi.GetProbability(AliPID::kPion);
		if (ppion > nPidXi.GetProbability(AliPID::kElectron) &&
		    ppion > nPidXi.GetProbability(AliPID::kMuon)     &&
		    ppion > nPidXi.GetProbability(AliPID::kKaon)     &&
		    ppion > nPidXi.GetProbability(AliPID::kProton)   )     lIsNegInXiPion = kTRUE;

		// Check if the V0 negative track is an anti-proton (case for Xi+)
		Double_t pproton = nPidXi.GetProbability(AliPID::kProton);
		if (pproton > nPidXi.GetProbability(AliPID::kElectron) &&
		    pproton > nPidXi.GetProbability(AliPID::kMuon)     &&
		    pproton > nPidXi.GetProbability(AliPID::kPion)     &&
		    pproton > nPidXi.GetProbability(AliPID::kKaon)     )     lIsNegInXiProton = kTRUE;
		
		// Check if the V0 negative track is a pi- (case for Omega-)
		ppion = 0.;
		    ppion = nPidOmega.GetProbability(AliPID::kPion);
		if (ppion > nPidOmega.GetProbability(AliPID::kElectron) &&
		    ppion > nPidOmega.GetProbability(AliPID::kMuon)     &&
		    ppion > nPidOmega.GetProbability(AliPID::kKaon)     &&
		    ppion > nPidOmega.GetProbability(AliPID::kProton)   )    lIsNegInOmegaPion = kTRUE;
		
		// Check if the V0 negative track is an anti-proton (case for Omega+)
		pproton = 0.;
		         pproton = nPidOmega.GetProbability(AliPID::kProton);
		if (pproton > nPidOmega.GetProbability(AliPID::kElectron) &&
		    pproton > nPidOmega.GetProbability(AliPID::kMuon)     &&
		    pproton > nPidOmega.GetProbability(AliPID::kPion)     &&
		    pproton > nPidOmega.GetProbability(AliPID::kKaon)     )  lIsNegInOmegaProton = kTRUE;
		
	}// end if V0 negative track with existing combined PID	
	
		
	// Combined bachelor PID
	AliPID bachPidXi;	bachPidXi.SetPriors(    lPriorsGuessXi    );
	AliPID bachPidOmega;	bachPidOmega.SetPriors( lPriorsGuessOmega );
	
	if( bachTrackXi->IsOn(AliESDtrack::kESDpid) ){  // Combined PID exists
		Double_t r[10] = {0.}; bachTrackXi->GetESDpid(r);
		bachPidXi.SetProbabilities(r);
		bachPidOmega.SetProbabilities(r);
		// Check if the bachelor track is a pion
		Double_t ppion = bachPidXi.GetProbability(AliPID::kPion);
		if (ppion > bachPidXi.GetProbability(AliPID::kElectron) &&
		    ppion > bachPidXi.GetProbability(AliPID::kMuon)     &&
		    ppion > bachPidXi.GetProbability(AliPID::kKaon)     &&
		    ppion > bachPidXi.GetProbability(AliPID::kProton)   )     lIsBachelorPion = kTRUE;
		// Check if the bachelor track is a kaon
		Double_t pkaon = bachPidOmega.GetProbability(AliPID::kKaon);
		if (pkaon > bachPidOmega.GetProbability(AliPID::kElectron) &&
		    pkaon > bachPidOmega.GetProbability(AliPID::kMuon)     &&
		    pkaon > bachPidOmega.GetProbability(AliPID::kPion)     &&
		    pkaon > bachPidOmega.GetProbability(AliPID::kProton)   )  lIsBachelorKaon = kTRUE;	
	}// end if bachelor track with existing combined PID
	
	
	// B - TPC PID : 3-sigma bands on Bethe-Bloch curve
	// Bachelor
	if (TMath::Abs(fTpcPidManager->GetNumberOfSigmas(bachTrackXi,AliPID::kKaon)) < 3) lIsBachelorKaonForTPC = kTRUE;
	if (TMath::Abs(fTpcPidManager->GetNumberOfSigmas(bachTrackXi,AliPID::kPion)) < 3) lIsBachelorPionForTPC = kTRUE;
	
	// Negative V0 daughter
	if (TMath::Abs(fTpcPidManager->GetNumberOfSigmas(nTrackXi,AliPID::kPion   )) < 3) lIsNegPionForTPC   = kTRUE;
	if (TMath::Abs(fTpcPidManager->GetNumberOfSigmas(nTrackXi,AliPID::kProton )) < 3) lIsNegProtonForTPC = kTRUE;
	
	// Positive V0 daughter
	if (TMath::Abs(fTpcPidManager->GetNumberOfSigmas(pTrackXi,AliPID::kPion   )) < 3) lIsPosPionForTPC   = kTRUE;
	if (TMath::Abs(fTpcPidManager->GetNumberOfSigmas(pTrackXi,AliPID::kProton )) < 3) lIsPosProtonForTPC = kTRUE;
	
	
		
		// - II.Step 6 : extra info for QA (ESD)
		// miscellaneous pieces of info that may help regarding data quality assessment.
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
	
	
	
		// II.Step 7 - Complementary info for monitoring the cascade cut variables
	
	const AliMultiplicity *lAliMult = lESDevent->GetMultiplicity();
	lSPDTrackletsMultiplicity = lAliMult->GetNumberOfTracklets();
	lBachTPCClusters          = bachTrackXi->GetTPCNcls(); // for the Bachelor at least, (we will see later for the baryon, for the meson)
		
	
        	// II.Step 8 - Azimuthal correlation study
		//-------------
	
	lTVect3MomXi.SetXYZ( lXiMomX, lXiMomY, lXiMomZ );
	lArrTrackID[0] = pTrackXi   ->GetID();
	lArrTrackID[1] = nTrackXi   ->GetID();
	lArrTrackID[2] = bachTrackXi->GetID();
	
	
  }// end of ESD treatment
  
 
  if(fAnalysisType == "AOD"){
	
	// -------------------------------------
	// II.AOD - Calcultaion Part dedicated to Xi vertices (ESD)
	
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

	
		// - II.Step 5 : PID on the daughters (To be developed ...)
		//-------------
	
	// Combined PID
	
	/* 
	// Reasonable guess for the priors for the cascade track sample
	Double_t lPriorsGuessXi[5]    = {0.0, 0.0, 2, 0, 1};
	Double_t lPriorsGuessOmega[5] = {0.0, 0.0, 1, 1, 1};
	AliPID bachPidXi;	bachPidXi.SetPriors(    lPriorsGuessXi    );
	AliPID bachPidOmega;	bachPidOmega.SetPriors( lPriorsGuessOmega );
	
	const AliAODTrack *bachTrackXi = lAODevent->GetTrack( xi->GetBachID() ); // FIXME : GetBachID not implemented ?
	
	if( bachTrackXi->IsOn(AliESDtrack::kESDpid) ){  // Combined PID exists, the AOD flags = a copy of the ESD ones
		Double_t r[10]; bachTrackXi->GetPID(r);
		bachPidXi.SetProbabilities(r);
		bachPidOmega.SetProbabilities(r);
		// Check if the bachelor track is a pion
		Double_t ppion = bachPidXi.GetProbability(AliPID::kPion);
		if (ppion > bachPidXi.GetProbability(AliPID::kElectron) &&
		    ppion > bachPidXi.GetProbability(AliPID::kMuon)     &&
		    ppion > bachPidXi.GetProbability(AliPID::kKaon)     &&
		    ppion > bachPidXi.GetProbability(AliPID::kProton)   )     lIsBachelorPion = kTRUE;
		// Check if the bachelor track is a kaon
		Double_t pkaon = bachPidOmega.GetProbability(AliPID::kKaon);
		if (pkaon > bachPidOmega.GetProbability(AliPID::kElectron) &&
		    pkaon > bachPidOmega.GetProbability(AliPID::kMuon)     &&
		    pkaon > bachPidOmega.GetProbability(AliPID::kPion)     &&
		    pkaon > bachPidOmega.GetProbability(AliPID::kProton)   )  lIsBachelorKaon = kTRUE;
		
	}// end if bachelor track with existing combined PID
	*/
	
	// TPC PID
	
		// - II.Step 6 : extra info for QA (AOD)
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

		// II.Step 7 - Complementary info for monitoring the cascade cut variables
	//FIXME
	
		// II.Step 8 - Azimuthal correlation study
		//-------------
	
	lTVect3MomXi.SetXYZ( lXiMomX, lXiMomY, lXiMomZ );
	
	AliAODTrack *pTrackXi    = dynamic_cast<AliAODTrack*>( xi->GetDaughter(0) );
	AliAODTrack *nTrackXi    = dynamic_cast<AliAODTrack*>( xi->GetDaughter(1) );
	AliAODTrack *bachTrackXi = dynamic_cast<AliAODTrack*>( xi->GetDecayVertexXi()->GetDaughter(0) );	
		if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
			Printf("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...");
			continue;
		}
		
	lArrTrackID[0] = pTrackXi   ->GetID();
	lArrTrackID[1] = nTrackXi   ->GetID();
	lArrTrackID[2] = bachTrackXi->GetID();
	
  }// end of AOD treatment


    // -------------------------------------
    // II.Fill - Filling the TH1,2,3Fs
  

	// - II.Fill.Step 1	 : primary vertex

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

		
	// II.Fill.Step 2
	fHistEffMassXi			->Fill( lEffMassXi               );
	fHistChi2Xi			->Fill( lChi2Xi                  );	// Flag CascadeVtxer: Cut Variable a
	fHistDcaXiDaughters		->Fill( lDcaXiDaughters          );	// Flag CascadeVtxer: Cut Variable e 
	fHistDcaBachToPrimVertex	->Fill( lDcaBachToPrimVertexXi   );	// Flag CascadeVtxer: Cut Variable d
	fHistXiCosineOfPointingAngle	->Fill( lXiCosineOfPointingAngle ); 	// Flag CascadeVtxer: Cut Variable f
	fHistXiRadius			->Fill( lXiRadius                );	// Flag CascadeVtxer: Cut Variable g+h
	
	
	// II.Fill.Step 3
	fHistMassLambdaAsCascDghter	->Fill( lInvMassLambdaAsCascDghter );	// Flag CascadeVtxer: Cut Variable c
	fHistV0Chi2Xi			->Fill( lV0Chi2Xi                  );	
	fHistDcaV0DaughtersXi		->Fill( lDcaV0DaughtersXi          );
	fHistV0CosineOfPointingAngleXi	->Fill( lV0CosineOfPointingAngleXi ); 
	fHistV0RadiusXi			->Fill( lV0RadiusXi                );
	
	fHistDcaV0ToPrimVertexXi	->Fill( lDcaV0ToPrimVertexXi       );	// Flag CascadeVtxer: Cut Variable b
	fHistDcaPosToPrimVertexXi	->Fill( lDcaPosToPrimVertexXi      );
	fHistDcaNegToPrimVertexXi	->Fill( lDcaNegToPrimVertexXi      );
		
	
	// II.Fill.Step 4 : extra QA info
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
		
	
	// II.Fill.Step 5 : inv mass plots 1D
	if( lChargeXi < 0 ){
					fHistMassXiMinus	       ->Fill( lInvMassXiMinus    );
					fHistMassOmegaMinus	       ->Fill( lInvMassOmegaMinus );
		if(lIsBachelorPion)	fHistMassWithCombPIDXiMinus    ->Fill( lInvMassXiMinus    );
		if(lIsBachelorKaon)	fHistMassWithCombPIDOmegaMinus ->Fill( lInvMassOmegaMinus );
	}
	
	if( lChargeXi > 0 ){
					fHistMassXiPlus		       ->Fill( lInvMassXiPlus     );
					fHistMassOmegaPlus	       ->Fill( lInvMassOmegaPlus  );
		if(lIsBachelorPion)	fHistMassWithCombPIDXiPlus     ->Fill( lInvMassXiPlus     );
		if(lIsBachelorKaon)	fHistMassWithCombPIDOmegaPlus  ->Fill( lInvMassOmegaPlus  );
	}
	
	
	// II.Fill.Step 6 : inv mass plots 2D, 3D
	if( lChargeXi < 0 ) {
		f2dHistEffMassLambdaVsEffMassXiMinus->Fill( lInvMassLambdaAsCascDghter, lInvMassXiMinus ); 
		f2dHistEffMassXiVsEffMassOmegaMinus ->Fill( lInvMassXiMinus, lInvMassOmegaMinus );
		f2dHistXiRadiusVsEffMassXiMinus     ->Fill( lXiRadius, lInvMassXiMinus );
		f2dHistXiRadiusVsEffMassOmegaMinus  ->Fill( lXiRadius, lInvMassOmegaMinus );
		f3dHistXiPtVsEffMassVsYXiMinus      ->Fill( lXiTransvMom, lInvMassXiMinus,    lRapXi    );
		f3dHistXiPtVsEffMassVsYOmegaMinus   ->Fill( lXiTransvMom, lInvMassOmegaMinus, lRapOmega );
		if(lIsPosInXiProton)  				f3dHistXiPtVsEffMassVsYWithCombPIDXiMinus     ->Fill( lXiTransvMom, lInvMassXiMinus,    lRapXi    );
		if(lIsBachelorKaon)  				f3dHistXiPtVsEffMassVsYWithCombPIDOmegaMinus  ->Fill( lXiTransvMom, lInvMassOmegaMinus, lRapOmega );
		if(lIsBachelorKaon && lIsPosInOmegaProton)  	f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaMinus ->Fill( lXiTransvMom, lInvMassOmegaMinus, lRapOmega );
		if(lIsBachelorKaonForTPC)			f3dHistXiPtVsEffMassVsYWithTpcPIDOmegaMinus   ->Fill( lXiTransvMom, lInvMassOmegaMinus, lRapOmega );
	}
	else{
		f2dHistEffMassLambdaVsEffMassXiPlus ->Fill( lInvMassLambdaAsCascDghter, lInvMassXiPlus );
		f2dHistEffMassXiVsEffMassOmegaPlus  ->Fill( lInvMassXiPlus, lInvMassOmegaPlus );
		f2dHistXiRadiusVsEffMassXiPlus      ->Fill( lXiRadius, lInvMassXiPlus);
		f2dHistXiRadiusVsEffMassOmegaPlus   ->Fill( lXiRadius, lInvMassOmegaPlus );
		f3dHistXiPtVsEffMassVsYXiPlus       ->Fill( lXiTransvMom, lInvMassXiPlus,    lRapXi    );
		f3dHistXiPtVsEffMassVsYOmegaPlus    ->Fill( lXiTransvMom, lInvMassOmegaPlus, lRapOmega );
		if(lIsNegInXiProton)  				f3dHistXiPtVsEffMassVsYWithCombPIDXiPlus     ->Fill( lXiTransvMom, lInvMassXiPlus,    lRapXi    );
		if(lIsBachelorKaon )			   	f3dHistXiPtVsEffMassVsYWithCombPIDOmegaPlus  ->Fill( lXiTransvMom, lInvMassOmegaPlus, lRapOmega );
		if(lIsBachelorKaon && lIsNegInOmegaProton)   	f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaPlus ->Fill( lXiTransvMom, lInvMassOmegaPlus, lRapOmega );
	}
	
	// - Filling the AliCFContainers related to PID
	
	Double_t lContainerPIDVars[4] = {0.0};
	
	// Xi Minus		
	if( lChargeXi < 0 ) {
		lContainerPIDVars[0] = lXiTransvMom       ;
		lContainerPIDVars[1] = lInvMassXiMinus    ;
		lContainerPIDVars[2] = lRapXi             ;
		lContainerPIDVars[3] = nTrackMultiplicity ;
			
		// No PID
			fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorPionForTPC  )
			fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 1); // TPC PID / 3-#sigma cut on Bachelor track
		
		if( lIsBachelorPionForTPC && 
		    lIsPosProtonForTPC     )
			fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 2); // TPC PID / 3-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorPionForTPC && 
		    lIsPosProtonForTPC    && 
		    lIsNegPionForTPC       )
			fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 3); // TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks
		
		// Combined PID
		if( lIsBachelorPion        )
			fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		
		if( lIsBachelorPion       && 
		    lIsPosInXiProton    )
			fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		
		if(lIsBachelorPion     && 
		   lIsPosInXiProton && 
		   lIsNegInXiPion    )
		 	fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
	}
	
	lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; lContainerPIDVars[3] = 0.;
	
	// Xi Plus		
	if( lChargeXi > 0 ) {
		lContainerPIDVars[0] = lXiTransvMom       ;
		lContainerPIDVars[1] = lInvMassXiPlus     ;
		lContainerPIDVars[2] = lRapXi             ;
		lContainerPIDVars[3] = nTrackMultiplicity ;
			
		// No PID
			fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorPionForTPC  )
			fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 1); // TPC PID / 3-#sigma cut on Bachelor track
		
		if( lIsBachelorPionForTPC && 
		    lIsNegProtonForTPC     )
			fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 2); // TPC PID / 3-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorPionForTPC && 
		    lIsNegProtonForTPC    && 
		    lIsPosPionForTPC       )
			fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 3); // TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks
		
		// Combined PID
		if( lIsBachelorPion        )
			fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		
		if( lIsBachelorPion       && 
		    lIsNegInXiProton    )
			fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		
		if(lIsBachelorPion     && 
		   lIsNegInXiProton && 
		   lIsPosInXiPion    )
		 	fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
	}
	
	lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; lContainerPIDVars[3] = 0.;
	
	// Omega Minus		
	if( lChargeXi < 0 ) {
		lContainerPIDVars[0] = lXiTransvMom       ;
		lContainerPIDVars[1] = lInvMassOmegaMinus ;
		lContainerPIDVars[2] = lRapOmega          ;
		lContainerPIDVars[3] = nTrackMultiplicity ;
			
		// No PID
			fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorKaonForTPC  )
			fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 1); // TPC PID / 3-#sigma cut on Bachelor track
		
		if( lIsBachelorKaonForTPC && 
		    lIsPosProtonForTPC     )
			fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 2); // TPC PID / 3-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorKaonForTPC && 
		    lIsPosProtonForTPC    && 
		    lIsNegPionForTPC       )
			fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 3); // TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks
		
		// Combined PID
		if( lIsBachelorKaon        )
			fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		
		if( lIsBachelorKaon       && 
		    lIsPosInOmegaProton    )
			fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		
		if(lIsBachelorKaon     && 
		   lIsPosInOmegaProton && 
		   lIsNegInOmegaPion    )
		 	fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
	}
	
	lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; lContainerPIDVars[3] = 0.;
	
	// Omega Plus		
	if( lChargeXi > 0 ) {
		lContainerPIDVars[0] = lXiTransvMom       ;
		lContainerPIDVars[1] = lInvMassOmegaPlus ;
		lContainerPIDVars[2] = lRapOmega          ;
		lContainerPIDVars[3] = nTrackMultiplicity ;
			
		// No PID
			fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorKaonForTPC  )
			fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 1); // TPC PID / 3-#sigma cut on Bachelor track
		
		if( lIsBachelorKaonForTPC && 
		    lIsNegProtonForTPC     )
			fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 2); // TPC PID / 3-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorKaonForTPC && 
		    lIsNegProtonForTPC    && 
		    lIsPosPionForTPC       )
			fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 3); // TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks
		
		// Combined PID
		if( lIsBachelorKaon        )
			fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		
		if( lIsBachelorKaon       && 
		    lIsNegInOmegaProton    )
			fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		
		if(lIsBachelorKaon     && 
		   lIsNegInOmegaProton && 
		   lIsPosInOmegaPion    )
		 	fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
	}
	
	

	
	
	// II.Fill.Step 7 : filling the AliCFContainer (optimisation of topological selections)
	Double_t lContainerCutVars[18] = {0.0};
			
	lContainerCutVars[0]  = lDcaXiDaughters;
	lContainerCutVars[1]  = lDcaBachToPrimVertexXi;
	lContainerCutVars[2]  = lXiCosineOfPointingAngle;
	lContainerCutVars[3]  = lXiRadius;
	lContainerCutVars[4]  = lInvMassLambdaAsCascDghter;
	lContainerCutVars[5]  = lDcaV0DaughtersXi;
	lContainerCutVars[6]  = lV0CosineOfPointingAngleXi;
	lContainerCutVars[7]  = lV0RadiusXi;
	lContainerCutVars[8]  = lDcaV0ToPrimVertexXi;	
	lContainerCutVars[9]  = lDcaPosToPrimVertexXi;
	lContainerCutVars[10] = lDcaNegToPrimVertexXi;
	
	if( lChargeXi < 0 ) { 
		lContainerCutVars[11] = lInvMassXiMinus;
	        lContainerCutVars[12] = lInvMassOmegaMinus;
	}
	else {
		lContainerCutVars[11] = lInvMassXiPlus;
		lContainerCutVars[12] = lInvMassOmegaPlus;
	}
	
	lContainerCutVars[13] = lXiTransvMom;	
	lContainerCutVars[14] = lBestPrimaryVtxPos[2];
	lContainerCutVars[15] = nTrackMultiplicity;
	lContainerCutVars[16] = lSPDTrackletsMultiplicity; // FIXME : SPDTrackletsMultiplicity is not available for AOD ... = -1
	lContainerCutVars[17] = lBachTPCClusters;          // FIXME : BachTPCClusters          is not available for AOD ... = -1
	
	if( lChargeXi < 0 ) fCFContCascadeCuts->Fill(lContainerCutVars,0); // for negative cascades = Xi- and Omega-
	else                fCFContCascadeCuts->Fill(lContainerCutVars,1); // for negative cascades = Xi+ and Omega+
  
  			
	// II.Fill.Step 8 :  angular correlations
	
	if( lChargeXi < 0 ){
		DoAngularCorrelation("Xi-",    lInvMassXiMinus,    lArrTrackID, lTVect3MomXi, lEta);
		DoAngularCorrelation("Omega-", lInvMassOmegaMinus, lArrTrackID, lTVect3MomXi, lEta);
	}
	else{
		DoAngularCorrelation("Xi+",    lInvMassXiPlus,    lArrTrackID, lTVect3MomXi, lEta);
		DoAngularCorrelation("Omega+", lInvMassOmegaPlus, lArrTrackID, lTVect3MomXi, lEta);
	}
	
	
    }// end of the Cascade loop (ESD or AOD)
    
  
  // Post output data.
 PostData(1, fListHistCascade);
}


void AliAnalysisTaskCheckCascade::DoAngularCorrelation( const Char_t   *lCascType, 
							      Double_t  lInvMassCascade, 
							const Int_t    *lArrTrackID,
							      TVector3 &lTVect3MomXi, 
							      Double_t  lEtaXi ){
  // Perform the Delta(Phi)Delta(Eta) analysis 
  // by properly filling the THnSparseF 
	
	TString lStrCascType( lCascType );
	
	Double_t lCascPdgMass = 0.0;
	if( lStrCascType.Contains("Xi") )       lCascPdgMass = 1.3217;
	if( lStrCascType.Contains("Omega") )    lCascPdgMass = 1.6724;
	
	if( lInvMassCascade > lCascPdgMass + 0.010) return;
	if( lInvMassCascade < lCascPdgMass - 0.010) return;
	// Check the Xi- candidate is within the proper mass window m0 +- 10 MeV
	
	
	// 1st loop: check there is no track with a higher pt ...
	// = The cascade is meant to be a leading particle : Pt(Casc) > any track in the event
	for(Int_t TrckIdx = 0; TrckIdx < (InputEvent())->GetNumberOfTracks() ; TrckIdx++ )
	{// Loop over all the tracks of the event
	
		AliVTrack *lCurrentTrck = dynamic_cast<AliVTrack*>( (InputEvent())->GetTrack( TrckIdx ) );
			if (!lCurrentTrck ) {
				Printf("ERROR Correl. Study : Could not retrieve a track while looping over the event tracks ...");
				continue;
			}
		if(lTVect3MomXi.Pt() < lCurrentTrck->Pt() ) return;	
		// Room for improvement: //FIXME
		// 1. There is a given resolution on pt : maybe release the cut Pt(casc) < Pt(track)*90% ?
		// 2. Apply this cut only when DeltaPhi(casc, track) > 90 deg = when track is in the away-side ?
		// 3. Anti-splitting cut (like in Femto analysis) ?
			
	}// end control loop
	
	// 2nd loop: filling loop
	for(Int_t TrckIdx = 0; TrckIdx < (InputEvent())->GetNumberOfTracks() ; TrckIdx++ )
	{// Loop over all the tracks of the event
	
		AliVTrack *lCurrentTrck = dynamic_cast<AliVTrack*>( (InputEvent())->GetTrack( TrckIdx ) );
			if (!lCurrentTrck ) {
				Printf("ERROR Correl. Study : Could not retrieve a track while looping over the event tracks ...");
				continue;
			}
				
		// Room for improvement: //FIXME
		// 1. Loop only on primary tracks ?	
		// 2. Exclude the tracks that build the condisdered cascade = the bachelor + the V0 dghters
		//     This may bias the outcome, especially for low multplicity events.
		// Note : For ESD event, track ID == track index.
			if(lCurrentTrck->GetID() == lArrTrackID[0]) continue;
			if(lCurrentTrck->GetID() == lArrTrackID[1]) continue;
			if(lCurrentTrck->GetID() == lArrTrackID[2]) continue;
			
		TVector3 lTVect3MomTrck(lCurrentTrck->Px(), lCurrentTrck->Py(), lCurrentTrck->Pz() );
		
		// 2 hypotheses made here :
		//   - The Xi trajectory is a straight line,
		//   - The Xi doesn't loose any energy by crossing the first layer(s) of ITS, if ever;
		//      So, meaning hyp: vect p(Xi) at the emission = vect p(Xi) at the decay vertex
		//      By doing this, we introduce a systematic error on the cascade Phi ...
		// Room for improvement: take into account the curvature of the Xi trajectory //FIXME
	
   		Double_t lHnSpFillVar[5] = {0.};
		lHnSpFillVar[0] = lTVect3MomXi.DeltaPhi(lTVect3MomTrck) * 180.0/TMath::Pi(); // Delta phi(Casc,Track) (deg)
		if(lHnSpFillVar[0] < -50.0) lHnSpFillVar[0] += 360.0; 		
		lHnSpFillVar[1] = lEtaXi - lCurrentTrck->Eta(); 			   // Delta eta(Casc,Track)
		lHnSpFillVar[2] = lTVect3MomXi.Pt();					   // Pt_{Casc}
		lHnSpFillVar[3] = lCurrentTrck->Pt();					   // Pt_{any track}
		lHnSpFillVar[4] = lInvMassCascade;					   // Eff. Inv Mass (control var)
		
		if(      lStrCascType.Contains("Xi-") )      fHnSpAngularCorrXiMinus    ->Fill( lHnSpFillVar );
		else if( lStrCascType.Contains("Xi+") )      fHnSpAngularCorrXiPlus     ->Fill( lHnSpFillVar );
		else if( lStrCascType.Contains("Omega-") )   fHnSpAngularCorrOmegaMinus ->Fill( lHnSpFillVar );
		else if( lStrCascType.Contains("Omega+") )   fHnSpAngularCorrOmegaPlus  ->Fill( lHnSpFillVar );
	
	}// end - Loop over all the tracks in the event

}







//________________________________________________________________________
void AliAnalysisTaskCheckCascade::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  TList *cRetrievedList = 0x0;
         cRetrievedList = (TList*)GetOutputData(1);
	if(!cRetrievedList){
		Printf("ERROR - AliAnalysisTaskCheckCascade: ouput data container list not available\n"); return;
	}

  fHistTrackMultiplicity = dynamic_cast<TH1F*> (   cRetrievedList->FindObject("fHistTrackMultiplicity") );
  if (!fHistTrackMultiplicity) {
		Printf("ERROR - AliAnalysisTaskCheckCascade: fHistTrackMultiplicity not available\n"); return;
	}
 
  fHistCascadeMultiplicity = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistCascadeMultiplicity") );
	if (!fHistCascadeMultiplicity) {
		Printf("ERROR - AliAnalysisTaskCheckCascade: fHistCascadeMultiplicity not available\n"); return;
	}
	
  fHistMassXiMinus    = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassXiMinus") );	
	if (!fHistMassXiMinus) {
		Printf("ERROR - AliAnalysisTaskCheckCascade: fHistMassXiMinus not available\n"); return;
	}
  fHistMassXiPlus     = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassXiPlus") );
	if (!fHistMassXiPlus) {
		Printf("ERROR - AliAnalysisTaskCheckCascade: fHistMassXiPlus not available\n"); return;
	}	
  fHistMassOmegaMinus = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassOmegaMinus") );
	if (!fHistMassOmegaMinus) {
		Printf("ERROR - AliAnalysisTaskCheckCascade: fHistMassOmegaMinus not available\n"); return;
	}
  fHistMassOmegaPlus  = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassOmegaPlus") );	
	if (!fHistMassOmegaPlus) {
		Printf("ERROR - AliAnalysisTaskCheckCascade: fHistMassOmegaPlus not available\n"); return;
	}
  
  TCanvas *canCheckCascade = new TCanvas("AliAnalysisTaskCheckCascade","CheckCascade overview",10,10,510,510);
  canCheckCascade->Divide(2,2);
  
  canCheckCascade->cd(1);
  canCheckCascade->cd(1)->SetLogy();
  fHistTrackMultiplicity->DrawCopy("H");
  
  canCheckCascade->cd(2);  
  canCheckCascade->cd(2)->SetLogy();
  fHistCascadeMultiplicity->SetMarkerStyle(kOpenSquare);
  fHistCascadeMultiplicity->DrawCopy("E");
  
  canCheckCascade->cd(3);  
  fHistMassXiMinus->SetMarkerStyle(kFullCircle);
  fHistMassXiMinus ->SetMarkerSize(0.5);
  fHistMassXiMinus->DrawCopy("E");
  fHistMassXiPlus ->SetMarkerStyle(kOpenCircle);
  fHistMassXiPlus ->SetMarkerSize(0.5);
  fHistMassXiPlus ->DrawCopy("ESAME");
  
  canCheckCascade->cd(4);  
  fHistMassOmegaMinus->SetMarkerStyle(kFullCircle);
  fHistMassOmegaMinus->SetMarkerSize(0.5);
  fHistMassOmegaMinus->DrawCopy("E");
  fHistMassOmegaPlus ->SetMarkerStyle(kOpenCircle);
  fHistMassOmegaPlus ->SetMarkerSize(0.5);
  fHistMassOmegaPlus ->DrawCopy("ESAME");

}
