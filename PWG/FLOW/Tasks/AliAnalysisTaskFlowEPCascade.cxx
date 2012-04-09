#include <iostream>

#include "TFile.h"
#include "TChain.h"
#include "TList.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"
#include "AliVVertex.h"

#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDVZERO.h"
#include "AliESDUtils.h"

#include "AliAODEvent.h"

#include "AliESDInputHandler.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"

#include "AliFlowTrackSimple.h"
//#include "AliFlowEventCuts.h"
#include "AliFlowTrackCuts.h"

#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliPIDResponse.h"

#include "AliOADBContainer.h"

#include "AliAnalysisTaskFlowEPCascade.h"

#include "TMath.h"

#include "AliLog.h"

using namespace std;

ClassImp(AliAnalysisTaskFlowEPCascade)

//_______________________________________________________________
  AliAnalysisTaskFlowEPCascade::AliAnalysisTaskFlowEPCascade():
    AliAnalysisTaskSE()
    ,fCutsDau (0x0)
    ,fPIDResponse (0x0)
    ,fMinCent(0)
    ,fMaxCent(0)
    ,fVtxCut(10.)
    ,fOADB(0x0)
    ,fRun(-1)
    ,fICent(-1)
    ,fMultV0(0x0)
    ,fV0Cpol(100)
    ,fV0Apol(100)
    ,fHistList       (0x0)
    ,fhEvent (0x0)
    ,fhEPangleVZero (0x0)
    ,fhEPangleV0A(0x0)
    ,fhEPangleV0C(0x0)
    ,fhEPangleTPC(0x0)
    ,fh1Chi2Xi(0x0)
    ,fh1DCAXiDaughters(0x0)
    ,fh1DCABachToPrimVertex(0x0)
    ,fh1XiCosOfPointingAngle(0x0)
    ,fh1XiRadius(0x0)
    ,fh1MassLambda(0x0)
    ,fh1V0Chi2(0x0)
    ,fh1V0CosOfPointingAngle(0x0)
    ,fh1V0Radius(0x0)
    ,fh1DcaV0DaughtersXi(0x0)
    ,fh1DcaV0ToPrimVertex(0x0)
    ,fh1DCAPosToPrimVertex(0x0)
    ,fh1DCANegToPrimVertex(0x0)
    ,fh1MassXiMinus(0x0)
    ,fh1MassXiPlus(0x0)
    ,fh1MassOmegaMinus(0x0)
    ,fh1MassOmegaPlus(0x0)
    ,fh1MassXi(0x0)
    ,fh1MassOmega(0x0)
    ,fh1XiPt(0x0)
    ,fh1XiP(0x0)
    ,fh1XiBachPt(0x0)
    ,fh1XiBachP(0x0)
    ,fh1ChargeXi(0x0)
    ,fh1V0toXiCosOfPointingAngle(0x0)
    ,fh1PhiXi(0x0)
    ,fh2Armenteros(0x0)
    ,fh2MassLambdaVsMassXiMinus(0x0)
    ,fh2MassXiVsMassOmegaMinus(0x0)
    ,fh2MassLambdaVsMassXiPlus(0x0)
    ,fh2MassXiVsMassOmegaPlus(0x0)
    ,fh2XiRadiusVsMassXiMinus(0x0)
    ,fh2XiRadiusVsMassXiPlus(0x0)
    ,fh2XiRadiusVsMassOmegaMinus(0x0)
    ,fh2XiRadiusVsMassOmegaPlus(0x0)
    ,fh2TPCdEdxOfCascDghters(0x0)
    ,fh2MassVsPtXiMinus(0x0)
    ,fh2MassVsPtXiPlus(0x0)
    ,fh2MassVsPtXiAll(0x0)
    ,fh2MassVsPtOmegaMinus(0x0)
    ,fh2MassVsPtOmegaPlus(0x0)
    ,fh2MassVsPtOmegaAll(0x0)
    ,fhXiRapidity (0x0)
    ,fhOmegaRapidity (0x0)
    ,fProfResolution (0x0)
    ,fh1DistToVtxZAfter(0x0)
    ,fh1DistToVtxXYAfter(0x0)
    ,fh2DistToVtxZBeforeVsAfter(0x0)
    ,fh2DistToVtxXYBeforeVsAfter(0x0)
    ,fh2PxBeforeVsAfter(0x0)
    ,fh2PyBeforeVsAfter(0x0)
    ,fh2PhiPosBeforeVsAfter(0x0)
    ,fh2PhiNegBeforeVsAfter(0x0)
{

  for(Int_t i = 0; i < 3; i++){
    fProfXiV2PtV0A[i] = NULL;
    fProfOmegaV2PtV0A[i] = NULL;
    fProfXiSinePtV0A[i] = NULL;
    fProfOmegaSinePtV0A[i] = NULL;

    fProfXiV2PtV0C[i] = NULL;
    fProfOmegaV2PtV0C[i] = NULL;
    fProfXiSinePtV0C[i] = NULL;
    fProfOmegaSinePtV0C[i] = NULL;

    fProfXiV2Pt[i] = NULL;
    fProfOmegaV2Pt[i] = NULL;
    fProfXiSinePt[i] = NULL;
    fProfOmegaSinePt[i] = NULL;
    
    fProf2dXiV2PtV0A[i] = NULL;
    fProf2dOmegaV2PtV0A[i] = NULL;
    fProf2dXiV2PtV0C[i] = NULL;
    fProf2dOmegaV2PtV0C[i] = NULL;
    fProf2dXiV2Pt[i] = NULL;
    fProf2dOmegaV2Pt[i] = NULL;
  }

  for(int i=0; i!=3; ++i)
    for(int j=0; j!=2; ++j) {
      fXiBands[i][j] = 0;
      fOmegaBands[i][j] = 0;
    }

  for(Int_t i = 0; i != 2; ++i)
    for(Int_t j = 0; j != 2; ++j)
      for(Int_t iC = 0; iC < 9; iC++){
	fMeanQ[iC][i][j] = 0.;
	fWidthQ[iC][i][j] = 0.;
    }

}

//_________________________________________________________________________
AliAnalysisTaskFlowEPCascade::
AliAnalysisTaskFlowEPCascade(const char *name, double centMin, double centMax,
			     double xis[3][2],
			     double omegas[3][2]) 
  : AliAnalysisTaskSE(name)
  ,fCutsDau (0x0)
  ,fPIDResponse (0x0)
  ,fMinCent(centMin)
  ,fMaxCent(centMax)
  ,fVtxCut(10.)
  ,fOADB(0x0)
  ,fRun(-1)
  ,fICent(-1)
  ,fMultV0(0x0)
  ,fV0Cpol(100)
  ,fV0Apol(100)
  ,fHistList       (0x0)
  ,fhEvent (0x0)
  ,fhEPangleVZero (0x0)
  ,fhEPangleV0A(0x0)
  ,fhEPangleV0C(0x0)
  ,fhEPangleTPC(0x0)
  ,fh1Chi2Xi(0x0)
  ,fh1DCAXiDaughters(0x0)
  ,fh1DCABachToPrimVertex(0x0)
  ,fh1XiCosOfPointingAngle(0x0)
  ,fh1XiRadius(0x0)
  ,fh1MassLambda(0x0)
  ,fh1V0Chi2(0x0)
  ,fh1V0CosOfPointingAngle(0x0)
  ,fh1V0Radius(0x0)
  ,fh1DcaV0DaughtersXi(0x0)
  ,fh1DcaV0ToPrimVertex(0x0)
  ,fh1DCAPosToPrimVertex(0x0)
  ,fh1DCANegToPrimVertex(0x0)
  ,fh1MassXiMinus(0x0)
  ,fh1MassXiPlus(0x0)
  ,fh1MassOmegaMinus(0x0)
  ,fh1MassOmegaPlus(0x0)
  ,fh1MassXi(0x0)
  ,fh1MassOmega(0x0)
  ,fh1XiPt(0x0)
  ,fh1XiP(0x0)
  ,fh1XiBachPt(0x0)
  ,fh1XiBachP(0x0)
  ,fh1ChargeXi(0x0)
  ,fh1V0toXiCosOfPointingAngle(0x0)
  ,fh1PhiXi(0x0)
  ,fh2Armenteros(0x0)
  ,fh2MassLambdaVsMassXiMinus(0x0)
  ,fh2MassXiVsMassOmegaMinus(0x0)
  ,fh2MassLambdaVsMassXiPlus(0x0)
  ,fh2MassXiVsMassOmegaPlus(0x0)
  ,fh2XiRadiusVsMassXiMinus(0x0)
  ,fh2XiRadiusVsMassXiPlus(0x0)
  ,fh2XiRadiusVsMassOmegaMinus(0x0)
  ,fh2XiRadiusVsMassOmegaPlus(0x0)
  ,fh2TPCdEdxOfCascDghters(0x0)
  ,fh2MassVsPtXiMinus(0x0)
  ,fh2MassVsPtXiPlus(0x0)
  ,fh2MassVsPtXiAll(0x0)
  ,fh2MassVsPtOmegaMinus(0x0)
  ,fh2MassVsPtOmegaPlus(0x0)
  ,fh2MassVsPtOmegaAll(0x0)
  ,fhXiRapidity (0x0)
  ,fhOmegaRapidity (0x0)
  ,fProfResolution (0x0)
  ,fh1DistToVtxZAfter(0x0)
  ,fh1DistToVtxXYAfter(0x0)
  ,fh2DistToVtxZBeforeVsAfter(0x0)
  ,fh2DistToVtxXYBeforeVsAfter(0x0)
  ,fh2PxBeforeVsAfter(0x0)
  ,fh2PyBeforeVsAfter(0x0)
  ,fh2PhiPosBeforeVsAfter(0x0)
  ,fh2PhiNegBeforeVsAfter(0x0)
{

  for(Int_t i = 0; i < 3; i++){
    fProfXiV2PtV0A[i] = NULL;
    fProfOmegaV2PtV0A[i] = NULL;
    fProfXiSinePtV0A[i] = NULL;
    fProfOmegaSinePtV0A[i] = NULL;

    fProfXiV2PtV0C[i] = NULL;
    fProfOmegaV2PtV0C[i] = NULL;
    fProfXiSinePtV0C[i] = NULL;
    fProfOmegaSinePtV0C[i] = NULL;

    fProfXiV2Pt[i] = NULL;
    fProfOmegaV2Pt[i] = NULL;
    fProfXiSinePt[i] = NULL;
    fProfOmegaSinePt[i] = NULL;

    fProf2dXiV2PtV0A[i] = NULL;
    fProf2dOmegaV2PtV0A[i] = NULL;
    fProf2dXiV2PtV0C[i] = NULL;
    fProf2dOmegaV2PtV0C[i] = NULL;
    fProf2dXiV2Pt[i] = NULL;
    fProf2dOmegaV2Pt[i] = NULL;

  }

  for(int i=0; i!=3; ++i)
    for(int j=0; j!=2; ++j) {
      fXiBands[i][j] = xis[i][j];
      fOmegaBands[i][j] = omegas[i][j];
    }

  for(Int_t i = 0; i != 2; ++i)
    for(Int_t j = 0; j != 2; ++j)
      for(Int_t iC = 0; iC <9; iC++){
	fMeanQ[iC][i][j] = 0.;
	fWidthQ[iC][i][j] = 0.;
      }
  
  DefineInput( 0,TChain::Class());
  //DefineInput (1, TList::Class());
  DefineOutput(1, TList::Class());
}

void AliAnalysisTaskFlowEPCascade::UserCreateOutputObjects()
{
  cout<<"AliAnalysisTaskFlowEPCascade::UserCreateOutputObjects()"<<endl;
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler
    = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  TString oadbfilename = "$ALICE_ROOT/OADB/PWGCF/VZERO/VZEROcalibEP.root";
  fOADB = TFile::Open(oadbfilename.Data());

  if(!fOADB){
    printf("OADB file %s cannot be opened\n",oadbfilename.Data());
    return;
  }

  fHistList = new TList();
  fHistList ->SetOwner();
  
  //Add histograms to the List
  fhEvent = new TH1I("Event", "Number of Events", 3, 0, 3);
  fHistList->Add(fhEvent);
  
  fhEPangleVZero = new TH1F("hEPangleVZero", 
		       "EP from VZERO; #Psi; Number of Events",
		       15, 0., TMath::Pi());
  fHistList->Add(fhEPangleVZero);

  fhEPangleV0A = new TH1F("hEPangleV0A",
			    "EP from V0A; #Psi; Number of Events",
			    15, 0., TMath::Pi());
  fHistList->Add(fhEPangleV0A);

  fhEPangleV0C = new TH1F("hEPangleV0C",
			    "EP from V0C; #Psi; Number of Events",
			    15, 0., TMath::Pi());
  fHistList->Add(fhEPangleV0C);
  

  fhEPangleTPC = new TH1F("hEPangleTPC",
			  "EP from TPC; #Psi; Number of Events",
			  15, 0., TMath::Pi());
  fHistList->Add(fhEPangleTPC);

  fh1Chi2Xi = new TH1F("Chi2Xi", 
                       "Cascade #chi^{2}; #chi^{2}; Number of Cascades", 
                       160, 0, 160);
  fHistList->Add(fh1Chi2Xi);
  
  fh1DCAXiDaughters 
    = new TH1F( "DcaXiDaughters",  
                "DCA between Xi Daughters; DCA (cm); Number of Cascades", 
                100, 0., 0.5);
  fHistList->Add(fh1DCAXiDaughters);

  fh1DCABachToPrimVertex
    = new TH1F("DcaBachToPrimVertex", 
               "DCA of Bach. to Prim. Vertex; DCA (cm);Number of Cascades", 
               250, 0., 2.5);
  fHistList->Add(fh1DCABachToPrimVertex);

  fh1XiCosOfPointingAngle
    = new TH1F("XiCosineOfPointingAngle",
	       "Cos of Xi Pointing Angle; Cos (Xi Point.Angl);Number of Xis", 
	       200, 0.99, 1.0);
  fHistList->Add(fh1XiCosOfPointingAngle);
  
  fh1XiRadius = new TH1F("XiRadius",  
			 "Casc. decay transv. radius; r (cm); Counts" , 
			 1050, 0., 105.0 );
  fHistList->Add(fh1XiRadius);
  
  fh1MassLambda 
    = new TH1F("MassLambdaAsCascDghter",
	       "#Lambda assoc. to Casc. candidates; Eff. Mass (GeV/c^{2}); Counts", 300,1.00,1.3);
  fHistList->Add(fh1MassLambda);
  
  fh1V0Chi2 = new TH1F("V0Chi2Xi", 
		       "V0 #chi^{2}, in cascade; #chi^{2};Counts", 
		       160, 0, 40);
  fHistList->Add(fh1V0Chi2);
  
  fh1V0CosOfPointingAngle 
    = new TH1F("V0CosOfPointingAngleXi", 
	       "Cos of V0 Pointing Angle, in cascade;Cos(V0 Point. Angl); Counts", 
	       200, 0.98, 1.0);
  fHistList->Add(fh1V0CosOfPointingAngle);
  
  fh1V0Radius  = new TH1F("V0RadiusXi", 
                          "V0 decay radius, in cascade; radius (cm); Counts", 
                          1050, 0., 105.0);
  fHistList->Add(fh1V0Radius);
  
  fh1DcaV0DaughtersXi = new TH1F("DcaV0DaughtersXi", 
                                 "DCA between V0 daughters, in cascade;DCA (cm);Number of V0s", 120, 0., 0.6);
  fHistList->Add(fh1DcaV0DaughtersXi);
  
  fh1DcaV0ToPrimVertex = new TH1F("DcaV0ToPrimVertexXi", 
				  "DCA of V0 to Prim. Vertex, in cascade;DCA (cm);Number of Cascades", 200, 0., 1.);
  fHistList->Add(fh1DcaV0ToPrimVertex);
    
   fh1DCAPosToPrimVertex =
     new TH1F("DcaPosToPrimVertexXi", 
	      "DCA of V0 pos daughter to Prim. Vertex;DCA (cm);Counts", 
	      300, 0, 3);
   fHistList->Add(fh1DCAPosToPrimVertex);

  fh1DCANegToPrimVertex 
    =  new TH1F("DcaNegToPrimVertexXi", 
                "DCA of V0 neg daughter to Prim. Vertex;DCA (cm);Counts", 
                300, 0, 3);
  fHistList->Add(fh1DCANegToPrimVertex);

  fh1MassXiMinus 
    = new TH1F("MassXiMinus", 
               "#Xi^{-} candidates;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts",
               1600, 1.2, 2.0);
  fHistList->Add(fh1MassXiMinus);


  fh1MassXiPlus 
    = new TH1F("MassXiPlus", 
               "#Xi^{+} candidates;M(#bar{#Lambda}^{0}, #pi^{+}) (GeV/c^{2});Counts",
               1600, 1.2, 2.0);
  fHistList->Add(fh1MassXiPlus);

  fh1MassXi 
    = new TH1F("MassXi",
               "#Xi candidates;M(#bar{#Lambda}, #pi) (GeV/c^{2});Counts",
               1600, 1.2, 2.0);
  fHistList->Add(fh1MassXi);

  fh1MassOmegaMinus 
    = new TH1F("MassOmegaMinus", 
               "#Omega^{-} candidates; M(#Lambda, K^{-})(GeV/c^{2});Counts",
               2000, 1.5, 2.5);
  fHistList->Add(fh1MassOmegaMinus);

  fh1MassOmega 
    = new TH1F("MassOmega",
	       "#Omega candidates; M(#Lambda, K)(GeV/c^{2});Counts",
	       2000, 1.5, 2.5);
  fHistList->Add(fh1MassOmega);
  
  fh1MassOmegaPlus 
    =  new TH1F("MassOmegaPlus", 
                "#Omega^{+} candidates;M(#bar{#Lambda}^{0}, K^{+})(GeV/c^{2});Counts", 2000, 1.5, 2.5);
  fHistList->Add(fh1MassOmegaPlus);
  
  fh1XiPt 
    = new TH1F("XiPt" , 
               "#Xi Pt (cand. around the mass peak);p_{t}(#Xi)(GeV/c);Counts", 
               100, 0.0, 10.0);
  fHistList->Add(fh1XiPt);

  fh1XiP 
    =  new TH1F("XiTotMom", 
                "#Xi momentum (cand. around the mass peak); p_{tot}(#Xi)(GeV/c); Counts", 
                150, 0.0, 15.0);
  fHistList->Add(fh1XiP);

  fh1XiBachPt = new TH1F("XiBachPt", 
                         "#Xi Bach. transverse momentum (cand. around the mass peak) ; p_{t}(Bach.) (GeV/c); Counts", 
                         100, 0.0, 5.0);
  fHistList->Add(fh1XiBachPt);

  fh1XiBachP = new TH1F("BachTotMomXi", 
                        "#Xi Bach. momentum (cand. around the mass peak); p_{tot}(Bach.) (GeV/c); Counts", 100, 0.0, 5.0);
  fHistList->Add(fh1XiBachP);
  
  fh1ChargeXi = new TH1F("ChargeXi", 
			 "Charge of casc. candidates; Sign; Counts", 
			 5, -2.0, 3.0);
  fHistList->Add(fh1ChargeXi);
  
  fh1V0toXiCosOfPointingAngle 
    = new TH1F("V0toXiCosineOfPointingAngle", 
	       "Cos. of V0 Ptng Angl Xi vtx; Cos(V0 Point. Angl / Xi vtx); Counts", 
	       100, 0.99, 1.0);
  fHistList->Add(fh1V0toXiCosOfPointingAngle);
  
  fh1PhiXi 
    = new TH1F("PhiXi", 
               "#phi of #Xi candidates (around the mass peak); #phi (deg); Counts", 64, 0., 6.4);
  fHistList->Add(fh1PhiXi);

  fh2Armenteros 
    = new TH2F("Armenteros", 
               "#alpha_{Arm}(casc. cand.) Vs Pt_{Arm}(casc. cand.); #alpha_{Arm} ; Pt_{Arm} (GeV/c)", 140, -1.2, 1.2, 300, 0., 0.3);
  fHistList->Add(fh2Armenteros);

   fh2XiRadiusVsMassXiMinus 
     = new TH2F( "XiRadiusVsEffMassXiMinus", 
		 "Transv. R_{Xi Decay} Vs M_{#Xi^{-} candidates}; r_{cascade} (cm); M( #Lambda , #pi^{-} ) (GeV/c^{2}) ", 
		 450, 0., 45.0, 1600, 1.2, 2.0);
  fHistList->Add(fh2XiRadiusVsMassXiMinus);

  fh2XiRadiusVsMassXiPlus 
    = new TH2F( "XiRadiusVsEffMassXiPlus",
                "Transv. R_{Xi Decay} Vs M_{#Xi^{+} candidates}; r_{cascade} (cm); M( #Lambda , #pi^{+} ) (GeV/c^{2}) ",
		450, 0., 45.0, 1600, 1.2, 2.0);
  fHistList->Add(fh2XiRadiusVsMassXiPlus);

  fh2XiRadiusVsMassOmegaMinus
    = new TH2F( "XiRadiusVsEffMassOmegaMinus", 
                "Transv. R_{Xi Decay} Vs M_{#Omega^{-} candidates}; r_{cascade} (cm); M( #Lambda , K^{-} ) (GeV/c^{2})",  
                450, 0., 45.0, 2000, 1.5, 2.5);
  fHistList->Add(fh2XiRadiusVsMassOmegaMinus);

  fh2XiRadiusVsMassOmegaPlus 
    = new TH2F( "XiRadiusVsEffMassOmegaPlus",
		"Transv. R_{Xi Decay} Vs M_{#Omega^{+} candidates}; r_{cascade} (cm); M( #Lambda , K^{+} ) (GeV/c^{2}) ",
                450, 0., 45.0, 2000, 1.5, 2.5);
  fHistList->Add(fh2XiRadiusVsMassOmegaPlus);

  fh2MassLambdaVsMassXiMinus 
    = new TH2F( "f2dHistEffMassLambdaVsEffMassXiMinus", 
                "M_{#Lambda} Vs M_{#Xi^{-} candidates} ; Inv. M_{#Lambda^{0}}\
		 (GeV/c^{2}); M(#Lambda, #pi^{-}) (GeV/c^{2})", 
                300, 1.1, 1.13, 1600, 1.2, 2.0);
  fHistList->Add(fh2MassLambdaVsMassXiMinus);

  fh2MassXiVsMassOmegaMinus 
    = new TH2F( "EffMassXiVsEffMassOmegaMinus", 
                "M_{#Xi^{-} candidates} Vs M_{#Omega^{-} candidates} ; M( #Lambda , #pi^{-} ) (GeV/c^{2}) ; M( #Lambda , K^{-} ) (GeV/c^{2})", 
                1600, 1.2, 2.0, 2000, 1.5, 2.5);
  fHistList->Add(fh2MassXiVsMassOmegaMinus);

  fh2MassLambdaVsMassXiPlus
    = new TH2F("EffMassLambdaVsEffMassXiPlus", 
               "M_{#Lambda} Vs M_{#Xi^{+} candidates}; Inv. M_{#Lambda^{0}}(GeV/c^{2}); M( #Lambda , #pi^{+} ) (GeV/c^{2})", 
	       300, 1.1,1.13, 1600, 1.2, 2.0);
  fHistList->Add(fh2MassLambdaVsMassXiPlus);

  fh2MassXiVsMassOmegaPlus 
    = new TH2F("EffMassXiVsEffMassOmegaPlus", 
               "M_{#Xi^{+} candidates} Vs M_{#Omega^{+} candidates} ; M( #Lambda, #pi^{+} ) (GeV/c^{2}) ; M( #Lambda , K^{+} ) (GeV/c^{2})", 
               1600, 1.2, 2.0, 2000, 1.5, 2.5);
  fHistList->Add(fh2MassXiVsMassOmegaPlus);

  fh2TPCdEdxOfCascDghters
    = new TH2F( "TPCdEdxOfCascDghters",
                "TPC dE/dx of the cascade daughters; charge x || #vec{p}_{TPC inner wall}(Casc. daughter) || (GeV/c); TPC signal (ADC) ", 
                2000, -10.0, 10.0, 450, 0., 900.);
  fHistList->Add(fh2TPCdEdxOfCascDghters);

  fh2MassVsPtXiMinus 
    = new TH2F("MassVsPtXiMinus", 
	       "M_{#Xi^{-} candidates} vs Pt; Pt (GeV/c); M(#Lambda, #pi^{-}) (GeV/c^{2})",   
	       100, 0., 10., 1600, 1.2, 2.0);
  fHistList->Add(fh2MassVsPtXiMinus);

  fh2MassVsPtXiPlus
    = new TH2F("MassVsPtXiPlus",
               "M_{#Xi^{+} candidates} vs Pt; Pt (GeV/c); M(#Lambda, #pi^{+})(GeV/c^{2})",
               100, 0., 10., 1600, 1.2, 2.0);
  fHistList->Add(fh2MassVsPtXiPlus);

  fh2MassVsPtXiAll 
    = new TH2F("MassVsPtXiAll", 
               "M_{#Xi candidates} vs Pt; Pt (GeV/c); M(#Lambda, #pi) (GeV/c^{2})", 
               100, 0., 10., 1600, 1.2, 2.0);
  fHistList->Add(fh2MassVsPtXiAll);

  fh2MassVsPtOmegaMinus 
    = new TH2F("MassVsPtOmegaMinus", 
               "M_{#Omega^{-} candidates} vs Pt; Pt (GeV/c); M(#Lambda, K^{-}) (GeV/c^{2})", 
	       100, 0., 10., 2000, 1.5, 2.5);
  fHistList->Add(fh2MassVsPtOmegaMinus);

  fh2MassVsPtOmegaPlus 
    = new TH2F("MassVsPtOmegaPlus",
	       "M_{#Omega^{+} candidates} vs Pt; Pt (GeV/c); M(#Lambda, K^{+}) (GeV/c^{2})", 
	       100, 0., 10., 2000, 1.5, 2.5);
  fHistList->Add(fh2MassVsPtOmegaPlus);

  fh2MassVsPtOmegaAll
    = new TH2F("MassVsPtOmegaAll",
	       "M_{#Omega candidates} vs Pt; Pt (GeV/c); M(#Lambda, K^{+}) (GeV/c^{2})",
	       100, 0., 10., 2000, 1.5, 2.5);
  fHistList->Add(fh2MassVsPtOmegaAll);

  fhXiRapidity = new TH1F("hXiRapidity", 
			  "#Xi rapidity distribution before rap. cut; y; Number of counts",
			  20, -1., 1.);
  fHistList->Add(fhXiRapidity);
  
  fhOmegaRapidity = new TH1F("hOmegaRapidity",
			     "#Omega rapidity distr. before rap. cut; y; Number of counts", 
			     20, -1., 1.);
  fHistList->Add(fhOmegaRapidity);

  for(Int_t i = 0; i < 3; i++){
    fProfXiV2PtV0A[i] = new TProfile(Form("hProfXiV2PtV0A%d", i), 
				     "; p_{T}[GeV/c]; v_{2}", 
				     100, 0., 10.);
    fHistList->Add(fProfXiV2PtV0A[i]);

    fProfOmegaV2PtV0A[i] = new TProfile(Form("hProfOmegaV2PtV0A%d", i),
					"; p_{T}[GeV/c]; v_{2}",
					100, 0., 10.);
    fHistList->Add(fProfOmegaV2PtV0A[i]);

    fProfXiSinePtV0A[i] = new TProfile(Form("hProfXiSinePtV0A%d", i),
				       ";p_{T}[GeV/c]; <sin[2(#phi-#Psi)]>",
				       100, 0., 10.);
    fHistList->Add(fProfXiSinePtV0A[i]);

    fProfOmegaSinePtV0A[i] 
      = new TProfile(Form("hProfOmegaSinePtV0A%d", i),
		     "; p_{T}[GeV/c]; <sin[2(#phi-#Psi)]>",
		     100, 0., 10.);
    fHistList->Add(fProfOmegaSinePtV0A[i]);
    

    fProfXiV2PtV0C[i] = new TProfile(Form("hProfXiV2PtV0C%d", i), 
				     "; p_{T}[GeV/c]; v_{2}", 
				     100, 0., 10.);
    fHistList->Add(fProfXiV2PtV0C[i]);

    fProfOmegaV2PtV0C[i] = new TProfile(Form("hProfOmegaV2PtV0C%d", i),
					"; p_{T}[GeV/c]; v_{2}",
					100, 0., 10.);
    fHistList->Add(fProfOmegaV2PtV0C[i]);

    fProfXiSinePtV0C[i] = new TProfile(Form("hProfXiSinePtV0C%d", i),
				       ";p_{T}[GeV/c]; <sin[2(#phi-#Psi)]>",
				       100, 0., 10.);
    fHistList->Add(fProfXiSinePtV0C[i]);

    fProfOmegaSinePtV0C[i] 
      = new TProfile(Form("hProfOmegaSinePtV0C%d", i),
		     "; p_{T}[GeV/c]; <sin[2(#phi-#Psi)]>",
		     100, 0., 10.);
    fHistList->Add(fProfOmegaSinePtV0C[i]);


    fProfXiV2Pt[i] = new TProfile(Form("hProfXiV2Pt%d", i), 
				  "; p_{T}[GeV/c]; v_{2}", 
				  100, 0., 10.);
    fHistList->Add(fProfXiV2Pt[i]);

    fProfOmegaV2Pt[i] = new TProfile(Form("hProfOmegaV2Pt%d", i),
                                  "; p_{T}[GeV/c]; v_{2}",
                                  100, 0., 10.);
    fHistList->Add(fProfOmegaV2Pt[i]);

    fProfXiSinePt[i] = new TProfile(Form("hProfXiSinePt%d", i),
                                  ";p_{T}[GeV/c]; <sin[2(#phi-#Psi)]>",
                                  100, 0., 10.);
    fHistList->Add(fProfXiSinePt[i]);

    fProfOmegaSinePt[i] = new TProfile(Form("hProfOmegaSinePt%d", i),
				     "; p_{T}[GeV/c]; <sin[2(#phi-#Psi)]>",
				     100, 0., 10.);
    fHistList->Add(fProfOmegaSinePt[i]);

    fProf2dXiV2PtV0A[i] = new TProfile2D(Form("hProf2dXiV2PtV0A%d", i), 
					 "; p_{T}[GeV/c]; #Psi; v_{2}",
					 20, 0., 10., 
					 15, 0., TMath::Pi());
    fHistList->Add(fProf2dXiV2PtV0A[i]);

    fProf2dOmegaV2PtV0A[i] = new TProfile2D(Form("hProf2dOmegaV2PtV0A%d", i),
					    "; p_{T}[GeV/c]; #Psi; v_{2}",
					    20, 0., 10.,
					    15, 0., 
					    TMath::Pi());
    fHistList->Add(fProf2dOmegaV2PtV0A[i]);

    fProf2dXiV2PtV0C[i] = new TProfile2D(Form("hProf2dXiV2PtV0C%d", i),
                                         "; p_{T}[GeV/c]; #Psi; v_{2}",
                                         20, 0., 10.,
                                         15, 0., TMath::Pi());
    fHistList->Add(fProf2dXiV2PtV0C[i]);
  
    fProf2dOmegaV2PtV0C[i] = new TProfile2D(Form("hProf2dOmegaV2PtV0C%d", i),
                                            "; p_{T}[GeV/c]; #Psi; v_{2}",
                                            20, 0., 10.,
                                            15, 0.,
                                            TMath::Pi());
    fHistList->Add(fProf2dOmegaV2PtV0C[i]);
    
    fProf2dXiV2Pt[i] = new TProfile2D(Form("hProf2dXiV2Pt%d", i),
				      "; p_{T}[GeV/c]; #Psi; v_{2}",
				      20, 0., 10.,
				      15, 0., TMath::Pi());
    fHistList->Add(fProf2dXiV2Pt[i]);

    fProf2dOmegaV2Pt[i] = new TProfile2D(Form("hProf2dOmegaV2Pt%d", i),
					 "; p_{T}[GeV/c]; #Psi; v_{2}",
					 20, 0., 10.,
					 15, 0.,
					 TMath::Pi());
    fHistList->Add(fProf2dOmegaV2Pt[i]);
  }
  
  fh1DistToVtxZAfter = new TH1F("DistToVtxZAfter",
				"Distance to vtx z after propagation to vtx; z [cm]",
				100, -5., 5.);
  fHistList->Add(fh1DistToVtxZAfter);
  
  fh1DistToVtxXYAfter = new TH1F("DistToVtxXYAfter",
				 "Distance to vtx xy after propagation to vtx",
				 500, 0., 50.);
  fHistList->Add(fh1DistToVtxXYAfter);
  
  fh2DistToVtxZBeforeVsAfter
    = new TH2F("DistToVtxZBeforeVsAfter",
	       "Distance to vtx z before vs after propagation; Distance before [cm]; Distance after [cm]",
	       500, -50., 50., 100, -5., 5.);
  fHistList->Add(fh2DistToVtxZBeforeVsAfter);
  
  fh2DistToVtxXYBeforeVsAfter
    = new TH2F("DistToVtxXYBeforeVsAfter",
	       "Distance to vtx xy before vs after propagation; Distance before [cm]; Distance after [cm]",
	       500, 0., 50, 500, 0., 50);
  fHistList->Add(fh2DistToVtxXYBeforeVsAfter);
  
  fh2PxBeforeVsAfter
    = new TH2F("PxBeforeVsAfter",
	       "Px before vs after propagation; Px [GeV/c]; Px' [GeV/c]",
	       200, -10., 10, 200, -10., 10.);
  fHistList->Add(fh2PxBeforeVsAfter);
  
  fh2PyBeforeVsAfter
    = new TH2F("PyBeforeVsAfter",
	       "Py before vs after propagation; Py [GeV/c]; Py' [GeV/c]",
	       200, -10., 10, 200, -10., 10.);
  fHistList->Add(fh2PyBeforeVsAfter);
  
  fh2PhiPosBeforeVsAfter
    = new TH2F("PhiPosBeforeVsAfter",
	       "Phi for positively charged candidates before vs after propagation; #phi; #phi'",
	       360, 0., 2.0*TMath::Pi(), 360, 0., 2.0*TMath::Pi());
  fHistList->Add(fh2PhiPosBeforeVsAfter);
  
  fh2PhiNegBeforeVsAfter
    = new TH2F("PhiNegBeforeVsAfter",
	       "Phi for negatively charged candidates before vs after propagation; #phi; #phi'",
	       360, 0., 2.0*TMath::Pi(), 360, 0., 2.0*TMath::Pi());
  fHistList->Add(fh2PhiNegBeforeVsAfter);
  
  fProfResolution = new TProfile("hProfResolution", 
				 "correlations between event planes", 
				 4, 0.5, 4.5);
  (fProfResolution->GetXaxis())
    ->SetBinLabel(1, "<cos[2(#Psi^{V0A})-#Psi^{V0C}]>");
  (fProfResolution->GetXaxis())
    ->SetBinLabel(2, "<cos[2(#Psi^{V0A})-#Psi^{TPC}]>");
  (fProfResolution->GetXaxis())
    ->SetBinLabel(3, "<cos[2(#Psi^{V0C})-#Psi^{TPC}]>");
  (fProfResolution->GetXaxis())
    ->SetBinLabel(4, "<cos[2(#Psi^{ZDCA})-#Psi^{ZDCC}]>");
  fHistList->Add(fProfResolution);
  
  PostData(1, fHistList);
}

void AliAnalysisTaskFlowEPCascade::UserExec(Option_t *) 
{

  AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  AliAODEvent *fAOD=dynamic_cast<AliAODEvent*>(InputEvent());

  if(fESD){
    
    Info("UserExec", "This task doesn't work with ESD!");
    //fhEvent->Fill(0);
    
    // if (!fFlowEventCuts->IsSelected(InputEvent())) return;
    
    // fhEvent->Fill(2);
    //    ReadFromESDv0(fESD);
  }
  else if(fAOD){
    //routine for AOD info

    fhEvent->Fill(0);

    //At the momment the cutting class does not handle AOD event properly
    //so we are doing the cuts explicitly here
    AliAODHeader *aodHeader = fAOD->GetHeader();
    if(!aodHeader) return;
    AliCentrality *centrality = aodHeader->GetCentralityP();
    if(!centrality) return;
    Double_t cent = centrality->GetCentralityPercentile("V0M" );
    if(cent<fMinCent||cent>=fMaxCent) return; //centrality cut

    if(cent >= 0 && cent < 5) fICent = 0;
    else if(cent >= 5 && cent < 10) fICent = 1;
    else if(cent >= 10 && cent < 20) fICent = 2;
    else if(cent >= 20 && cent < 30) fICent = 3;
    else if(cent >= 30 && cent < 40) fICent = 4;
    else if(cent >= 40 && cent < 50) fICent = 5;
    else if(cent >= 50 && cent < 60) fICent = 6;
    else if(cent >= 60 && cent < 70) fICent = 7;
    else if(cent >= 70 && cent < 80) fICent = 8;
    else 
      return;

    Double_t zvtx = fAOD->GetPrimaryVertex()->GetZ();
    if(TMath::Abs(zvtx) > fVtxCut) return; //vertex cut
    
    fhEvent->Fill(2);
    ReadFromAODv0(fAOD);
  }
  
  return;
}

/*
void AliAnalysisTaskFlowEPCascade::ReadFromESDv0(AliESDEvent *fESD){

  //  fEPcalib->CalculateEventPlanes(fESD);

  //Get EP informations
  //Double_t psiVZero = fEPcalib->GetPsi(2, AliAnaEPcalib::kV0AC);
  //fhEPangleVZero->Fill(psiVZero);

  //Double_t psiV0A = fEPcalib->GetPsi(2, AliAnaEPcalib::kV0A);
  //Double_t psiV0C = fEPcalib->GetPsi(2, AliAnaEPcalib::kV0C);
  //Double_t psiTPC = fEPcalib->GetPsi(2, AliAnaEPcalib::kTPC);

  //  Double_t psiZDC = fEPcalib->GetPsi(1, AliAnaEPcalib::kZDCAC);
  //fhEPangleZDC->Fill(psiZDC);
  
  //Double_t psiZDCA = fEPcalib->GetPsi(1, AliAnaEPcalib::kZDCA);
  //Double_t psiZDCC = fEPcalib->GetPsi(1, AliAnaEPcalib::kZDCC);

  AliEventplane * ep = fESD->GetEventplane();
  Double_t psiTPC = ep->GetEventplane("Q", fESD, 2);
  
  Int_t run = fESD->GetRunNumber();  
  if(run != fRun){
    // Load the calibrations run dependent
    OpenInfoCalbration(run);
    fRun=run;
  }

  //fill profile for resolution estimation
  fProfResolution->Fill(1, TMath::Cos(2.*(psiV0A - psiV0C)));
  fProfResolution->Fill(2, TMath::Cos(2.*(psiV0A - psiTPC)));
  fProfResolution->Fill(3, TMath::Cos(2.*(psiV0C - psiTPC)));
  //fProfResolution->Fill(4, TMath::Cos(2.*(psiZDCA - psiZDCC)));

  //check Xi candidates
  Double_t trkPrimaryVtxPos[3] = {-100., -100., -100.};
  Double_t bestPrimaryVtxPos[3] = {-100., -100., -100.};
  int nCascades=fESD->GetNumberOfCascades();
  
  const AliESDVertex *primaryTrackingESDVtx = fESD->GetPrimaryVertexTracks();
  primaryTrackingESDVtx->GetXYZ(trkPrimaryVtxPos);
  
  const AliESDVertex *primaryBestESDVtx = fESD->GetPrimaryVertex();
  primaryBestESDVtx->GetXYZ(bestPrimaryVtxPos);

  Double_t b = fESD->GetMagneticField();
  
  for(int i = 0; i != nCascades; ++i) {
    Double_t effMassXi = 0.;
    Double_t chi2Xi = -1.;
    Double_t dcaXiDaughters = -1.;
    Double_t XiCosOfPointingAngle = -1.;
    Double_t posXi[3] = {-1000., -1000., -1000.};
    Double_t XiRadius = -1000.;
    
    Double_t invMassLambdaAsCascDghter = 0.;
    Double_t V0Chi2Xi = -1.;
    Double_t dcaV0DaughtersXi = -1.;
    
    Double_t dcaBachToPrimaryVtxXi = -1.;
    Double_t dcaV0ToPrimaryVtxXi = -1.;
    Double_t dcaPosToPrimaryVtxXi = -1.;
    Double_t dcaNegToPrimaryVtxXi = -1.;
    Double_t V0CosOfPointingAngleXi = -1.;
    Double_t posV0Xi[3] = {-1000., -1000., -1000.};
    Double_t V0RadiusXi = -1000.;
    Double_t V0quality = 0.;

    Double_t invMassXiMinus = 0.;
    Double_t invMassXiPlus = 0.;
    Double_t invMassOmegaMinus = 0.;
    Double_t invMassOmegaPlus = 0.;
    
    Bool_t isPosInXiProton = kFALSE;
    Bool_t isPosInXiPion = kFALSE;
    Bool_t isPosInOmegaProton = kFALSE;
    Bool_t isPosInOmegaPion = kFALSE;
    
    Bool_t isNegInXiProton = kFALSE;
    Bool_t isNegInXiPion = kFALSE;
    Bool_t isNegInOmegaProton = kFALSE;
    Bool_t isNegInOmegaPion = kFALSE;

    Bool_t isBachelorKaon = kFALSE;
    Bool_t isBachelorPion = kFALSE;
    
    Bool_t isBachelorKaonForTPC = kFALSE;
    Bool_t isBachelorPionForTPC = kFALSE;
    Bool_t isNegPionForTPC = kFALSE;
    Bool_t isPosPionForTPC = kFALSE;
    Bool_t isNegProtonForTPC = kFALSE;
    Bool_t isPosProtonForTPC = kFALSE;

    Double_t XiPx = 0., XiPy = 0., XiPz = 0.;
    Double_t XiPt = 0.;
    Double_t XiPtot = 0.;

    Double_t bachPx = 0., bachPy = 0., bachPz = 0.;
    Double_t bachPt = 0.;
    Double_t bachPtot = 0.;
    
    //Short_t chargeXi = -2;
    Double_t V0toXiCosOfPointingAngle = 0.;
    
    Double_t rapXi = -20.;
    Double_t rapOmega = -20.;
    Double_t phi = 6.3;
    Double_t alphaXi = -200.;
    Double_t ptArmXi = -200.;

    Double_t distToVtxZBefore = -999.;
    Double_t distToVtxZAfter = -999.;
    Double_t distToVtxXYBefore = -999.;
    Double_t distToVtxXYAfter = -999.;
    Double_t XiPAfter[3] = {-999., -999., -999.};
    Double_t phiAfter = -999.;

    AliESDcascade *xi = fESD->GetCascade(i);
    if(!xi) continue;
    
    if(xi->Charge()<0)
      xi->ChangeMassHypothesis(V0quality, 3312); // Xi- hypothesis
    else if(xi->Charge() > 0)
      xi->ChangeMassHypothesis(V0quality, -3312);
    else continue;

    effMassXi = xi->GetEffMassXi();
    chi2Xi = xi->GetChi2Xi();
    dcaXiDaughters = xi->GetDcaXiDaughters();
    XiCosOfPointingAngle 
      = xi->GetCascadeCosineOfPointingAngle(bestPrimaryVtxPos[0],
                                            bestPrimaryVtxPos[1],
                                            bestPrimaryVtxPos[2]);
    xi->GetXYZcascade(posXi[0], posXi[1], posXi[2]);
    XiRadius = TMath::Sqrt(posXi[0]*posXi[0]
                           +posXi[1]*posXi[1]
                           +posXi[2]*posXi[2]);
    
    UInt_t idxPosXi = (UInt_t)TMath::Abs(xi->GetPindex());
    UInt_t idxNegXi = (UInt_t)TMath::Abs(xi->GetNindex());
    UInt_t idxBach = (UInt_t)TMath::Abs(xi->GetBindex());

    if(idxBach == idxPosXi || idxBach == idxNegXi) continue;

    AliESDtrack *pTrkXi = fESD->GetTrack(idxPosXi);
    AliESDtrack *nTrkXi = fESD->GetTrack(idxNegXi);
    AliESDtrack *bTrkXi = fESD->GetTrack(idxBach);

    if( !pTrkXi || !nTrkXi || !bTrkXi ) continue;

    if( !fCutsDau->IsSelected(pTrkXi) 
        || !fCutsDau->IsSelected(nTrkXi)
        || !fCutsDau->IsSelected(bTrkXi) ) continue;

    const AliExternalTrackParam *pExtTrk = pTrkXi->GetInnerParam();
    const AliExternalTrackParam *nExtTrk = nTrkXi->GetInnerParam();
    const AliExternalTrackParam *bExtTrk = bTrkXi->GetInnerParam();
    
    if(pExtTrk && pTrkXi->IsOn(AliESDtrack::kTPCin)){
      fh2TPCdEdxOfCascDghters->Fill(pExtTrk->GetP()*pExtTrk->Charge(), 
                                    pTrkXi->GetTPCsignal());
    }
    if(nExtTrk && nTrkXi->IsOn(AliESDtrack::kTPCin)){
      fh2TPCdEdxOfCascDghters->Fill(nExtTrk->GetP()*nExtTrk->Charge(),
                                    nTrkXi->GetTPCsignal());
    }
    if(bExtTrk && bTrkXi->IsOn(AliESDtrack::kTPCin)){
      fh2TPCdEdxOfCascDghters->Fill(bExtTrk->GetP()*bExtTrk->Charge(),
                                    bTrkXi->GetTPCsignal());
    }
    

    invMassLambdaAsCascDghter = xi->GetEffMass(); // from V0
    dcaV0DaughtersXi = xi->GetDcaV0Daughters();
    V0Chi2Xi = xi->GetChi2V0();
    
    V0CosOfPointingAngleXi 
      = xi->GetV0CosineOfPointingAngle(bestPrimaryVtxPos[0],
                                       bestPrimaryVtxPos[1],
                                       bestPrimaryVtxPos[2]);
    dcaV0ToPrimaryVtxXi = xi->GetD(bestPrimaryVtxPos[0], 
                                   bestPrimaryVtxPos[1],
                                   bestPrimaryVtxPos[2]);
    dcaBachToPrimaryVtxXi = TMath::Abs(bTrkXi->GetD(bestPrimaryVtxPos[0],
                                                    bestPrimaryVtxPos[1],
                                                    bestPrimaryVtxPos[2]));
    
    //V0
    xi->GetXYZ(posV0Xi[0], posV0Xi[1], posV0Xi[2]);
    V0RadiusXi = TMath::Sqrt(posV0Xi[0]*posV0Xi[0]
                             +posV0Xi[1]*posV0Xi[1]
                             +posV0Xi[2]*posV0Xi[2]);
    dcaPosToPrimaryVtxXi = TMath::Abs(pTrkXi->GetD(bestPrimaryVtxPos[0],
                                                   bestPrimaryVtxPos[1],
                                                   bestPrimaryVtxPos[2]));
    dcaNegToPrimaryVtxXi = TMath::Abs(nTrkXi->GetD(bestPrimaryVtxPos[0],
                                                   bestPrimaryVtxPos[1],
                                                   bestPrimaryVtxPos[2]));
    
    //apply cuts
    //    if(XiRadius < 1.5 || XiRadius > 100.) continue;
    //if(dcaXiDaughters > 0.3) continue;
    //if(XiCosOfPointingAngle < 0.999) continue;
    //if(dcaV0ToPrimaryVtxXi < 0.03) continue;
    //if(dcaBachToPrimaryVtxXi < 0.03) continue;
    
    //V0 mass cut?
    //if(TMath::Abs(invMassLambdaAsCascDghter-1.11568) > 0.006) continue;

    //if(dcaV0DaughtersXi > .6) continue;
    //if(V0CosOfPointingAngleXi > 0.9999) continue;
    //if(dcaPosToPrimaryVtxXi < 0.1) continue;
    //if(dcaNegToPrimaryVtxXi < 0.1) continue;

    //if(V0RadiusXi < 6.0 || V0RadiusXi > 100) continue;
    
    //other cuts?

    // change mass hypothesis to cover all the possibilities
    if(bTrkXi->Charge()<0){
      V0quality = 0.;
      xi->ChangeMassHypothesis(V0quality, 3312); //Xi- hyp.
      invMassXiMinus = xi->GetEffMassXi();

      V0quality = 0.;
      xi->ChangeMassHypothesis(V0quality, 3334); //Omega- hyp.
      invMassOmegaMinus = xi->GetEffMassXi();

      V0quality = 0.;
      xi->ChangeMassHypothesis(V0quality, 3312); //back to default hyp.
    }

    if(bTrkXi->Charge() > 0){
      V0quality = 0.;
      xi->ChangeMassHypothesis(V0quality, -3312); //anti-Xi- hyp.
      invMassXiPlus = xi->GetEffMassXi();

      V0quality = 0.;
      xi->ChangeMassHypothesis(V0quality, -3334); //anti-Omega- hyp.
      invMassOmegaPlus = xi->GetEffMassXi();

      V0quality = 0.;
      xi->ChangeMassHypothesis(V0quality, -3312); //back to default hyp.
    }

    //PID on the daughter tracks
    //A - Combined PID
    //Resonable guess the priors for the cascade track sample
    //(e, mu, pi, K, p)
    Double_t priorsGuessXi[5] = {0, 0, 2, 0, 1};
    Double_t priorsGuessOmega[5] = {0, 0, 1, 1, 1};

    //Combined bachelor-daughter PID
    AliPID pidXi;
    pidXi.SetPriors(priorsGuessXi);
    AliPID pidOmega;
    pidOmega.SetPriors(priorsGuessOmega);
    
    if(pTrkXi->IsOn(AliESDtrack::kESDpid)){// combined PID exists
      Double_t r[10] = {0.};
      pTrkXi->GetESDpid(r);
      pidXi.SetProbabilities(r);
      pidOmega.SetProbabilities(r);

      //Check if the V0 postive track is proton (case for Xi-)
      Double_t pProton = pidXi.GetProbability(AliPID::kProton);
      if(pProton > pidXi.GetProbability(AliPID::kElectron)
         && pProton > pidXi.GetProbability(AliPID::kMuon)
         && pProton > pidXi.GetProbability(AliPID::kPion)
         && pProton > pidXi.GetProbability(AliPID::kKaon))
        isPosInXiProton = kTRUE;
      
      //Check if the V0 postive track is a pi+ (case for Xi+)
      Double_t pPion = pidXi.GetProbability(AliPID::kPion);
      if(pPion > pidXi.GetProbability(AliPID::kElectron)
         && pPion > pidXi.GetProbability(AliPID::kMuon)
         && pPion > pidXi.GetProbability(AliPID::kKaon)
         && pPion > pidXi.GetProbability(AliPID::kProton))
        isPosInXiPion = kTRUE;
      // Check if the V0 positive track is a proton (case for Omega-)
      pProton = pidOmega.GetProbability(AliPID::kProton);
      if(pProton > pidOmega.GetProbability(AliPID::kElectron)
         && pProton > pidOmega.GetProbability(AliPID::kMuon)
         && pProton > pidOmega.GetProbability(AliPID::kPion)
         && pProton > pidOmega.GetProbability(AliPID::kKaon))
        isPosInOmegaProton = kTRUE;
    
      // Check if the V0 positive track is a pi+ (case for Omega+)
      pPion =  pidOmega.GetProbability(AliPID::kPion);
      if(pPion > pidOmega.GetProbability(AliPID::kElectron)
         && pPion > pidOmega.GetProbability(AliPID::kMuon)
         && pPion > pidOmega.GetProbability(AliPID::kKaon)
         && pPion > pidOmega.GetProbability(AliPID::kProton))
        isPosInOmegaPion = kTRUE;
    }
    
    //Combined V0-negative-daughter PID
    pidXi.SetPriors(priorsGuessXi);
    pidOmega.SetPriors(priorsGuessOmega);
    if(nTrkXi->IsOn(AliESDtrack::kESDpid)){
      Double_t r[10] = {0.};
      nTrkXi->GetESDpid(r);
      pidXi.SetProbabilities(r);
      pidOmega.SetProbabilities(r);
      
      // Check if the V0 negative track is a pi- (case for Xi-)
      Double_t pPion = pidXi.GetProbability(AliPID::kPion);
      if(pPion > pidXi.GetProbability(AliPID::kElectron)
         && pPion > pidXi.GetProbability(AliPID::kMuon)
         && pPion >  pidXi.GetProbability(AliPID::kKaon)
         && pPion > pidXi.GetProbability(AliPID::kProton))
        isNegInXiPion = kTRUE;

      // Check if the V0 negative track is an anti-proton (case for Xi+)
      Double_t pProton = pidXi.GetProbability(AliPID::kProton);
      if(pProton > pidXi.GetProbability(AliPID::kElectron)
         && pProton > pidXi.GetProbability(AliPID::kMuon)
         && pProton > pidXi.GetProbability(AliPID::kPion) 
         && pProton > pidXi.GetProbability(AliPID::kKaon))
        isNegInXiProton = kTRUE;
      
      // Check if the V0 negative track is a pi- (case for Omega-)
      pPion = pidOmega.GetProbability(AliPID::kPion);
      if(pPion > pidOmega.GetProbability(AliPID::kElectron)
         && pPion > pidOmega.GetProbability(AliPID::kMuon)
         && pPion > pidOmega.GetProbability(AliPID::kKaon)
         && pPion > pidOmega.GetProbability(AliPID::kProton))
        isNegInOmegaPion = kTRUE;
      
      // Check if the V0 negative track is an anti-proton (case for Omega+)   
      pProton =  pidOmega.GetProbability(AliPID::kProton);
      if(pProton > pidOmega.GetProbability(AliPID::kElectron)
         && pProton > pidOmega.GetProbability(AliPID::kMuon)
         && pProton > pidOmega.GetProbability(AliPID::kPion)
         && pProton > pidOmega.GetProbability(AliPID::kKaon))
        isNegInOmegaProton = kTRUE;
    }

    // Combined bachelor PID
    pidXi.SetPriors(priorsGuessXi);
    pidOmega.SetPriors(priorsGuessOmega);
    if(bTrkXi->IsOn(AliESDtrack::kESDpid)){//Combined PID exists
      Double_t r[10] = {0.};
      bTrkXi->GetESDpid(r);
      pidXi.SetProbabilities(r);
      pidOmega.SetProbabilities(r);
      
      //Check if the bachelor track is a pion
      Double_t pPion = pidXi.GetProbability(AliPID::kPion);
      if(pPion >  pidXi.GetProbability(AliPID::kElectron)
         && pPion >  pidXi.GetProbability(AliPID::kMuon)
         && pPion > pidXi.GetProbability(AliPID::kKaon)
         && pPion > pidXi.GetProbability(AliPID::kProton))
        isBachelorPion = kTRUE;
     
      // Check if the bachelor track is a kaon
      Double_t pKaon = pidOmega.GetProbability(AliPID::kKaon);
      if(pKaon > pidOmega.GetProbability(AliPID::kElectron)
         && pKaon > pidOmega.GetProbability(AliPID::kMuon)
         && pKaon > pidOmega.GetProbability(AliPID::kPion)
         && pKaon > pidOmega.GetProbability(AliPID::kProton))
        isBachelorKaon = kTRUE;
    }

    //B - TPC PID: 3-sigma bands on Bethe-Bloch curve
    //Bachelor
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bTrkXi, AliPID::kKaon))<3.)
      isBachelorKaonForTPC = kTRUE;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bTrkXi, AliPID::kPion))<3.)
      isBachelorPionForTPC = kTRUE;

    //Negative V0 daughter
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrkXi, AliPID::kPion))<3.)
      isNegPionForTPC = kTRUE;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrkXi, AliPID::kProton))<3.)
      isNegProtonForTPC = kTRUE;
    
    //Positive V0 daughter
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrkXi, AliPID::kPion))<3.)
      isPosPionForTPC = kTRUE;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrkXi, AliPID::kProton))<3.)
      isPosProtonForTPC = kTRUE;

    //    //Extra QA information
    xi->GetPxPyPz(XiPx, XiPy, XiPz);
    XiPt = TMath::Sqrt(XiPx*XiPx + XiPy*XiPy);
    XiPtot= TMath::Sqrt(XiPx*XiPx + XiPy*XiPy + XiPz*XiPz);
    
    xi->GetBPxPyPz(bachPx, bachPy, bachPz);
    bachPt = TMath::Sqrt(bachPx*bachPx + bachPy*bachPy);
    bachPtot = TMath::Sqrt(bachPx*bachPx + bachPy*bachPy + bachPz*bachPz);

     V0toXiCosOfPointingAngle 
       = xi->GetV0CosineOfPointingAngle(posXi[0], posXi[1], posXi[2]);
    rapXi = xi->RapXi();
    rapOmega = xi->RapOmega();
    phi = xi->Phi();
    alphaXi = xi->AlphaXi();
    ptArmXi = xi->PtArmXi();

    XiPAfter[0] = XiPx;
    XiPAfter[1] = XiPy;
    XiPAfter[2] = XiPz;

    distToVtxZBefore = posXi[2]-bestPrimaryVtxPos[2];
    distToVtxXYBefore
      = TMath::Sqrt((posXi[0] - bestPrimaryVtxPos[0])
                    *(posXi[0] - bestPrimaryVtxPos[0])
                    +(posXi[1] - bestPrimaryVtxPos[1])
                    *(posXi[1] - bestPrimaryVtxPos[1]));

    //propagation to the best primary vertex to determine the momentum   
    Propagate(bestPrimaryVtxPos, posXi, XiPAfter, b, xi->Charge());
    distToVtxZAfter = posXi[2] - bestPrimaryVtxPos[2];
    distToVtxXYAfter = TMath::Sqrt((posXi[0] - bestPrimaryVtxPos[0])
                                   *(posXi[0] - bestPrimaryVtxPos[0])
                                   +(posXi[1] - bestPrimaryVtxPos[1])
                                   *(posXi[1] - bestPrimaryVtxPos[1]));
    phiAfter = TMath::Pi() + TMath::ATan2(-XiPAfter[1],-XiPAfter[0]);

    fh1DistToVtxZAfter->Fill(distToVtxZAfter);
    fh1DistToVtxXYAfter->Fill(distToVtxXYAfter);
    fh2DistToVtxZBeforeVsAfter->Fill(distToVtxZBefore, distToVtxZAfter);
    fh2DistToVtxXYBeforeVsAfter->Fill(distToVtxXYBefore, distToVtxXYAfter);
    fh2PxBeforeVsAfter->Fill(XiPx, XiPAfter[0]);
    fh2PyBeforeVsAfter->Fill(XiPy, XiPAfter[1]);
    if(xi->Charge()>0)
      fh2PhiPosBeforeVsAfter->Fill(phi, phiAfter);
    else if(xi->Charge()<0)
      fh2PhiNegBeforeVsAfter->Fill(phi, phiAfter);


    if( isBachelorPion && 
	( (xi->Charge()<0 && isPosInXiProton && isNegInXiPion)
	  || (xi->Charge()>0 && isNegInXiProton && isPosInXiPion ) ) )
      {//Xi candidate
        //for default hypothesis
	fh1Chi2Xi->Fill(chi2Xi);
	fh1DCAXiDaughters->Fill(dcaXiDaughters);
	fh1DCABachToPrimVertex->Fill(dcaBachToPrimaryVtxXi);
	fh1XiCosOfPointingAngle->Fill(XiCosOfPointingAngle);
	fh1XiRadius->Fill(XiRadius);
	
	//V0
	fh1MassLambda->Fill(invMassLambdaAsCascDghter);
	fh1V0Chi2->Fill(V0Chi2Xi);
	fh1DcaV0DaughtersXi->Fill(dcaV0DaughtersXi);
	fh1V0CosOfPointingAngle->Fill(V0CosOfPointingAngleXi);
	fh1V0Radius->Fill(V0RadiusXi);
	fh1DcaV0ToPrimVertex->Fill(dcaV0ToPrimaryVtxXi);
	fh1DCAPosToPrimVertex->Fill(dcaPosToPrimaryVtxXi);
	fh1DCANegToPrimVertex->Fill(dcaNegToPrimaryVtxXi);
	fh1ChargeXi->Fill(xi->Charge());
	fh1V0toXiCosOfPointingAngle->Fill(V0toXiCosOfPointingAngle);
       
	if ( TMath::Abs( invMassXiMinus-1.3217 ) < 0.012 
	     || TMath::Abs( invMassXiPlus-1.3217 ) < 0.012) 
	  {// One InvMass should be different from 0             
	    fh1XiPt->Fill(XiPt);
	    fh1XiP->Fill(XiPtot);
	    
	    fh1XiBachPt->Fill(bachPt);
	    fh1XiBachP->Fill(bachPtot);
	    fh1PhiXi->Fill( xi->Phi() );
	  }
	fh2Armenteros->Fill(alphaXi, ptArmXi);
      }

    if(xi->Charge()<0){
      fh1MassXiMinus->Fill(invMassXiMinus);
      fh1MassOmegaMinus->Fill(invMassOmegaMinus);
      fh1MassXi->Fill(invMassXiMinus);
      fh1MassOmega->Fill(invMassOmegaMinus);
      
      fh2MassLambdaVsMassXiMinus->Fill(invMassLambdaAsCascDghter,
				       invMassXiMinus);
      fh2MassXiVsMassOmegaMinus->Fill(invMassXiMinus, 
				      invMassOmegaMinus);
      fh2XiRadiusVsMassXiMinus->Fill(XiRadius, invMassXiMinus);
      fh2XiRadiusVsMassOmegaMinus->Fill(XiRadius, invMassOmegaMinus);
    }

    if(xi->Charge() > 0){
      fh1MassXiPlus->Fill(invMassXiPlus);
      fh1MassOmegaPlus->Fill(invMassOmegaPlus);
      fh1MassXi->Fill(invMassXiPlus);
      fh1MassOmega->Fill(invMassOmegaPlus);
      
      fh2MassLambdaVsMassXiPlus->Fill(invMassLambdaAsCascDghter, 
				      invMassXiPlus);
      fh2MassXiVsMassOmegaPlus->Fill(invMassXiPlus,
				     invMassOmegaPlus);
      fh2XiRadiusVsMassXiPlus->Fill(XiRadius, invMassXiPlus);
      fh2XiRadiusVsMassOmegaPlus->Fill(XiRadius, invMassOmegaPlus);
    }
    

    Double_t phiNew 
      = ( XiPAfter[0] == 0. && XiPAfter[1] == 0. ) ? 
      0.0 : TMath::ATan2(XiPAfter[1], XiPAfter[0]);
    Double_t phiV0 = phiNew;
    phiV0 -= psiVZero;
    //    if(phiV0 < 0) phiV0 += 2.*TMath::Pi();
    //if(phiV0 > TMath::Pi()) phiV0 -= TMath::Pi();

    Double_t phiV0A = phiNew;
    phiV0A -= psiV0A;
    //if(phiV0A < 0) phiV0A += 2.*TMath::Pi();
    //if(phiV0A > TMath::Pi()) phiV0A -= TMath::Pi();

    Double_t phiV0C = phiNew;
    phiV0C -= psiV0C;
    //if(phiV0C < 0) phiV0C += 2.*TMath::Pi();
    //if(phiV0C > TMath::Pi()) phiV0C -= TMath::Pi();

    //PID cuts with TPC cuts
    if(xi->Charge() < 0){
      if(isPosProtonForTPC
         && isNegPionForTPC){
	if(isBachelorPionForTPC){
	  //Xi
	  fhXiRapidity->Fill(rapXi);
	  if(TMath::Abs(rapXi) < 0.8){
	    fh2MassVsPtXiMinus->Fill(XiPt, invMassXiMinus);
            fh2MassVsPtXiAll->Fill(XiPt, invMassXiMinus);
	    
	    for(int r=0; r!=3; ++r) {
              if(invMassXiMinus > fXiBands[r][0] 
                 && invMassXiMinus < fXiBands[r][1]){

		fProfXiV2PtV0A[r]->Fill(XiPt, TMath::Cos(2.*phiV0A));
                fProfXiSinePtV0A[r]->Fill(XiPt, TMath::Sin(2.*phiV0A));

		fProfXiV2PtV0C[r]->Fill(XiPt, TMath::Cos(2.*phiV0C));
                fProfXiSinePtV0C[r]->Fill(XiPt, TMath::Sin(2.*phiV0C));

                fProfXiV2Pt[r]->Fill(XiPt, TMath::Cos(2.*phiV0));
                fProfXiSinePt[r]->Fill(XiPt, TMath::Sin(2.*phiV0));
              }
            }
	      
	  }
	}
    
	if(isBachelorKaonForTPC){
	  //Omega
	  fhOmegaRapidity->Fill(rapOmega);
	  if(TMath::Abs(rapOmega) < 0.8){
	    fh2MassVsPtOmegaMinus->Fill(XiPt, invMassOmegaMinus);
            fh2MassVsPtOmegaAll->Fill(XiPt, invMassOmegaMinus);

	    for(int r=0; r!=3; ++r) {
              if(invMassOmegaMinus > fOmegaBands[r][0]
                 && invMassOmegaMinus < fOmegaBands[r][1]){

		fProfOmegaV2PtV0A[r]->Fill(XiPt, TMath::Cos(2.*phiV0A));
                fProfOmegaSinePtV0A[r]->Fill(XiPt, TMath::Sin(2.*phiV0A));
		
		fProfOmegaV2PtV0C[r]->Fill(XiPt, TMath::Cos(2.*phiV0C));
                fProfOmegaSinePtV0C[r]->Fill(XiPt, TMath::Sin(2.*phiV0C));

                fProfOmegaV2Pt[r]->Fill(XiPt, TMath::Cos(2.*phiV0));
                fProfOmegaSinePt[r]->Fill(XiPt, TMath::Sin(2.*phiV0));
              }
            }
	    

	  }
	}
      }
    }
  
    if(xi->Charge() > 0){
      if(isNegProtonForTPC
	 && isPosPionForTPC){
	if(isBachelorPionForTPC){
	  //Xi
	  fhXiRapidity->Fill(rapXi);
	  if(TMath::Abs(rapXi) < 0.8){
	    fh2MassVsPtXiPlus->Fill(XiPt, invMassXiPlus);
            fh2MassVsPtXiAll->Fill(XiPt, invMassXiPlus);
	    for(int r=0; r!=3; ++r) {
              if(invMassXiPlus > fXiBands[r][0]
                 && invMassXiPlus < fXiBands[r][1]){
		fProfXiV2PtV0A[r]->Fill(XiPt, TMath::Cos(2.*phiV0A));
                fProfXiSinePtV0A[r]->Fill(XiPt, TMath::Sin(2.*phiV0A));

                fProfXiV2PtV0C[r]->Fill(XiPt, TMath::Cos(2.*phiV0C));
                fProfXiSinePtV0C[r]->Fill(XiPt, TMath::Sin(2.*phiV0C));

                fProfXiV2Pt[r]->Fill(XiPt, TMath::Cos(2.*phiV0));
                fProfXiSinePt[r]->Fill(XiPt, TMath::Sin(2.*phiV0));
              }
            }

	  }
	}
      
	if(isBachelorKaonForTPC){
	  //Omega
	  fhOmegaRapidity->Fill(rapOmega);
	  if(TMath::Abs(rapOmega) < 0.8){
	    fh2MassVsPtOmegaPlus->Fill(XiPt, invMassOmegaPlus);
            fh2MassVsPtOmegaAll->Fill(XiPt, invMassOmegaPlus);
	    
	    for(int r=0; r!=3; ++r) {
              if(invMassOmegaPlus > fOmegaBands[r][0]
                 && invMassOmegaPlus < fOmegaBands[r][1]){
		fProfOmegaV2PtV0A[r]->Fill(XiPt, TMath::Cos(2.*phiV0A));
                fProfOmegaSinePtV0A[r]->Fill(XiPt, TMath::Sin(2.*phiV0A));

                fProfOmegaV2PtV0C[r]->Fill(XiPt, TMath::Cos(2.*phiV0C));
                fProfOmegaSinePtV0C[r]->Fill(XiPt, TMath::Sin(2.*phiV0C));
		
                fProfOmegaV2Pt[r]->Fill(XiPt, TMath::Cos(2.*phiV0));
                fProfOmegaSinePt[r]->Fill(XiPt, TMath::Sin(2.*phiV0));
              }
            }


	  }
	}
      }
    }
  
  }//end of cascade loop

  //  fNEvent++;

  PostData(1, fHistList);
}
*/

void AliAnalysisTaskFlowEPCascade::ReadFromAODv0(AliAODEvent *fAOD){

  AliEventplane * ep = (fAOD->GetHeader())->GetEventplaneP();
  Double_t psiTPC = ep->GetEventplane("Q", fAOD, 2); // in range of [0, pi]
  //  if(psiTPC > TMath::PiOver2()) 
  //  psiTPC -= TMath::Pi();
  fhEPangleTPC->Fill(psiTPC);

  Int_t run = fAOD->GetRunNumber();
  if(run != fRun){
    // Load the calibrations run dependent
    OpenInfoCalbration(run);                    
    fRun=run;          
  }

  //reset Q vector info       
  Double_t Qxa2 = 0, Qya2 = 0;
  Double_t Qxc2 = 0, Qyc2 = 0;
  //Double_t Qx2 = 0, Qy2 = 0;

  //V0 info   
  AliAODVZERO* aodV0 = fAOD->GetVZEROData();

  for (Int_t iv0 = 0; iv0 < 64; iv0++) {
    Double_t phiV0 = TMath::PiOver4()*(0.5 + iv0 % 8);
    Float_t multv0 = aodV0->GetMultiplicity(iv0);
    if (iv0 < 32){ // V0C 
      Qxc2 += TMath::Cos(2*phiV0)*multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
      Qyc2 += TMath::Sin(2*phiV0)*multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
    } else { // V0A
      Qxa2 += TMath::Cos(2*phiV0)*multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
      Qya2 += TMath::Sin(2*phiV0)*multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
    }
  }
  
  // Qx2 = Qxa2 + Qxc2;
  //Qy2 = Qya2 + Qya2;

  //grab for each centrality the proper histo with the Qx and Qy 
  //to do the recentering                    
  Double_t Qxamean2 = fMeanQ[fICent][1][0];
  Double_t Qxarms2  = fWidthQ[fICent][1][0];
  Double_t Qyamean2 = fMeanQ[fICent][1][1];
  Double_t Qyarms2  = fWidthQ[fICent][1][1];

  Double_t Qxcmean2 = fMeanQ[fICent][0][0];
  Double_t Qxcrms2  = fWidthQ[fICent][0][0];
  Double_t Qycmean2 = fMeanQ[fICent][0][1];
  Double_t Qycrms2  = fWidthQ[fICent][0][1];

  Double_t QxaCor2 = (Qxa2 - Qxamean2)/Qxarms2;
  Double_t QyaCor2 = (Qya2 - Qyamean2)/Qyarms2;
  Double_t QxcCor2 = (Qxc2 - Qxcmean2)/Qxcrms2;
  Double_t QycCor2 = (Qyc2 - Qycmean2)/Qycrms2;

  Double_t QxCor2 = (Qxa2 - Qxamean2 + Qxc2 - Qxcmean2)
    /TMath::Sqrt(Qxarms2*Qxarms2 + Qxcrms2*Qxcrms2);
  Double_t QyCor2 = (Qya2 - Qyamean2 + Qyc2 - Qycmean2)
    /TMath::Sqrt(Qyarms2*Qyarms2 + Qycrms2*Qycrms2);

  Double_t psiV0A =(TMath::Pi() + TMath::ATan2(-QyaCor2, -QxaCor2))/2.;
  Double_t psiV0C = (TMath::Pi() + TMath::ATan2(-QycCor2, -QxcCor2))/2.;
  Double_t psiVZero = (TMath::Pi() + TMath::ATan2(-QyCor2, -QxCor2))/2.;

  fhEPangleVZero->Fill(psiVZero);
  fhEPangleV0A->Fill(psiV0A);
  fhEPangleV0C->Fill(psiV0C);

  //fill profile for resolution estimation
  fProfResolution->Fill(1, TMath::Cos(2.*(psiV0A - psiV0C)));
  fProfResolution->Fill(2, TMath::Cos(2.*(psiV0A - psiTPC)));
  fProfResolution->Fill(3, TMath::Cos(2.*(psiV0C - psiTPC)));
  
  Double_t bestPrimaryVtxPos[3] = {-100., -100., -100.};

  Double_t b = fAOD->GetMagneticField();

  int nCascades=fAOD->GetNumberOfCascades();
  //Info("ReadFromAODv0", "# cascades = %d", nCascades);

  const AliAODVertex *primaryBestAODVtx = fAOD->GetPrimaryVertex();
  primaryBestAODVtx->GetXYZ(bestPrimaryVtxPos);
  
  for(Int_t iXi = 0; iXi < nCascades; iXi++){
    // Double_t effMassXi = 0.;
    Double_t chi2Xi = -1.;
    Double_t dcaXiDaughters = -1.;
    Double_t XiCosOfPointingAngle = -1.;
    Double_t posXi[3] = {-1000., -1000., -1000.};
    Double_t XiRadius = -1000.;
    
    Double_t invMassLambdaAsCascDghter = 0.;
    Double_t V0Chi2Xi = -1.;
    Double_t dcaV0DaughtersXi = -1.;
    
    Double_t dcaBachToPrimaryVtxXi = -1.;
    Double_t dcaV0ToPrimaryVtxXi = -1.;
    Double_t dcaPosToPrimaryVtxXi = -1.;
    Double_t dcaNegToPrimaryVtxXi = -1.;
    Double_t V0CosOfPointingAngleXi = -1.;
    Double_t posV0Xi[3] = {-1000., -1000., -1000.};
    Double_t V0RadiusXi = -1000.;
    
    Double_t invMassXiMinus = 0.;
    Double_t invMassXiPlus = 0.;
    Double_t invMassOmegaMinus = 0.;
    Double_t invMassOmegaPlus = 0.;

    /*
    Bool_t isPosInXiProton = kFALSE;
    Bool_t isPosInXiPion = kFALSE;
    Bool_t isPosInOmegaProton = kFALSE;
    Bool_t isPosInOmegaPion = kFALSE;
    
    Bool_t isNegInXiProton = kFALSE;
    Bool_t isNegInXiPion = kFALSE;
    Bool_t isNegInOmegaProton = kFALSE;
    Bool_t isNegInOmegaPion = kFALSE;

    Bool_t isBachelorKaon = kFALSE;
    Bool_t isBachelorPion = kFALSE;
    */

    Bool_t isBachelorKaonForTPC = kFALSE;
    Bool_t isBachelorPionForTPC = kFALSE;
    Bool_t isNegPionForTPC = kFALSE;
    Bool_t isPosPionForTPC = kFALSE;
    Bool_t isNegProtonForTPC = kFALSE;
    Bool_t isPosProtonForTPC = kFALSE;

    Double_t XiPx = 0., XiPy = 0., XiPz = 0.;
    Double_t XiPt = 0.;
    Double_t XiPtot = 0.;
    
    Double_t bachPx = 0., bachPy = 0., bachPz = 0.;
    Double_t bachPt = 0.;
    Double_t bachPtot = 0.;
    
    //Short_t chargeXi = -2;
    Double_t V0toXiCosOfPointingAngle = 0.;
    
    Double_t rapXi = -20.;
    Double_t rapOmega = -20.;
    Double_t phi = 6.3;
    Double_t alphaXi = -200.;
    Double_t ptArmXi = -200.;

    Double_t distToVtxZBefore = -999.;
    Double_t distToVtxZAfter = -999.;
    Double_t distToVtxXYBefore = -999.;
    Double_t distToVtxXYAfter = -999.;
    Double_t XiPAfter[3] = {-999., -999., -999.};
    Double_t phiAfter = -999.;

    const AliAODcascade *xi = fAOD->GetCascade(iXi);
    if (!xi) continue;
    
    //effMassXi = xi->MassXi(); //default working hypothesis: Xi- decay
    chi2Xi = xi->Chi2Xi();
    dcaXiDaughters = xi->DcaXiDaughters();
    XiCosOfPointingAngle = xi->CosPointingAngleXi(bestPrimaryVtxPos[0],
                                                  bestPrimaryVtxPos[1],
                                                  bestPrimaryVtxPos[2]);
    posXi[0] = xi->DecayVertexXiX();
    posXi[1] = xi->DecayVertexXiY();
    posXi[2] = xi->DecayVertexXiZ();
    XiRadius = TMath::Sqrt(posXi[0]*posXi[0]
                           +posXi[1]*posXi[1]
                           +posXi[2]*posXi[2]);
    
    AliAODTrack *pTrkXi = dynamic_cast<AliAODTrack*>( xi->GetDaughter(0) );
    AliAODTrack *nTrkXi = dynamic_cast<AliAODTrack*>( xi->GetDaughter(1) );
    AliAODTrack *bTrkXi 
      = dynamic_cast<AliAODTrack*>( xi->GetDecayVertexXi()->GetDaughter(0) );

    if(!pTrkXi || !nTrkXi || !bTrkXi) continue;
    
    UInt_t idxPosXi  = (UInt_t) TMath::Abs( pTrkXi->GetID() );
    UInt_t idxNegXi  = (UInt_t) TMath::Abs( nTrkXi->GetID() );
    UInt_t idxBach   = (UInt_t) TMath::Abs( bTrkXi->GetID() );

    if(idxBach == idxNegXi || idxBach == idxPosXi) continue;

    if( !fCutsDau->IsSelected(pTrkXi) 
        || !fCutsDau->IsSelected(nTrkXi)
        || !fCutsDau->IsSelected(bTrkXi) ) continue;
    
    if(pTrkXi->IsOn(AliESDtrack::kTPCin)){
      fh2TPCdEdxOfCascDghters->Fill(pTrkXi->P()*pTrkXi->Charge(), 
                                    pTrkXi->GetTPCsignal());
    }
    if( nTrkXi->IsOn(AliESDtrack::kTPCin) ){
      fh2TPCdEdxOfCascDghters->Fill(nTrkXi->P()*nTrkXi->Charge(),
                                    nTrkXi->GetTPCsignal());
    }
    if(bTrkXi->IsOn(AliESDtrack::kTPCin)){
      fh2TPCdEdxOfCascDghters->Fill(bTrkXi->P()*bTrkXi->Charge(),
                                    bTrkXi->GetTPCsignal());
    }


    if(xi->ChargeXi() < 0)
      invMassLambdaAsCascDghter = xi->MassLambda();
    else
      invMassLambdaAsCascDghter = xi->MassAntiLambda();

    dcaV0DaughtersXi = xi->DcaV0Daughters();
    V0Chi2Xi = xi->Chi2V0();
    
     V0CosOfPointingAngleXi 
       = xi->CosPointingAngle(bestPrimaryVtxPos);
     dcaV0ToPrimaryVtxXi = xi->DcaV0ToPrimVertex();
     dcaBachToPrimaryVtxXi = xi->DcaBachToPrimVertex();

     //V0
     posV0Xi[0] = xi->DecayVertexV0X();
     posV0Xi[1] = xi->DecayVertexV0Y();
     posV0Xi[2] = xi->DecayVertexV0Z();
     V0RadiusXi = TMath::Sqrt(posV0Xi[0]*posV0Xi[0]
                             +posV0Xi[1]*posV0Xi[1]
			      +posV0Xi[2]*posV0Xi[2]);
     dcaPosToPrimaryVtxXi = xi->DcaPosToPrimVertex();
     dcaNegToPrimaryVtxXi = xi->DcaNegToPrimVertex();


     //apply cuts ?
     // if(XiRadius < 1. || XiRadius > 100.) continue;
     //if(dcaXiDaughters > 0.1) continue;
     //if(XiCosOfPointingAngle < 0.999) continue;
     //if(dcaV0ToPrimaryVtxXi < 0.05) continue;
     //if(dcaBachToPrimaryVtxXi < 0.03) continue;
    
     //V0 mass cut?
     if(TMath::Abs(invMassLambdaAsCascDghter-1.11568) > 0.01) continue;

     //    if(dcaV0DaughtersXi > 1.) continue;
     //if(V0CosOfPointingAngleXi > 0.9999) continue;
     //if(dcaPosToPrimaryVtxXi < 0.1) continue;
     //if(dcaNegToPrimaryVtxXi < 0.1) continue;

     //if(V0RadiusXi < 1.0 || V0RadiusXi > 100) continue;
    
     //other cuts?

     if(xi->ChargeXi()<0){
       invMassXiMinus = xi->MassXi();
       invMassOmegaMinus = xi->MassOmega();
     }else{
       invMassXiPlus = xi->MassXi();
       invMassOmegaPlus = xi->MassOmega();
     }

     /*
     if(pTrkXi->GetMostProbablePID() == AliAODTrack::kProton) {
       isPosInXiProton = kTRUE;
       isPosInOmegaProton = kTRUE;
     }
     if(pTrkXi->GetMostProbablePID() == AliAODTrack::kPion){
       isPosInXiPion = kTRUE;
       isPosInOmegaPion = kTRUE;
     }
    
     if(nTrkXi->GetMostProbablePID() == AliAODTrack::kPion){
       isNegInXiPion = kTRUE;
       isNegInOmegaPion = kTRUE;
     }
     if(nTrkXi->GetMostProbablePID() == AliAODTrack::kProton){
       isNegInXiProton = kTRUE;
       isNegInOmegaProton = kTRUE;
     }

     if(bTrkXi->GetMostProbablePID() == AliAODTrack::kPion)
       isBachelorPion = kTRUE;
     if(bTrkXi->GetMostProbablePID() == AliAODTrack::kKaon)
       isBachelorKaon = kTRUE;
     */

     //PID with TPC only: ??? Fix me
     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bTrkXi, AliPID::kKaon))<3.)
       isBachelorKaonForTPC = kTRUE;
     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bTrkXi, AliPID::kPion))<3.)
       isBachelorPionForTPC = kTRUE;

     //Negative V0 daughter
     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrkXi, AliPID::kPion))<3.)
       isNegPionForTPC = kTRUE;
     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrkXi, AliPID::kProton))<3.)
       isNegProtonForTPC = kTRUE;
    
     //Positive V0 daughter
     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrkXi, AliPID::kPion))<3.)
       isPosPionForTPC = kTRUE;
     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrkXi, AliPID::kProton))<3.)
       isPosProtonForTPC = kTRUE;

     //Extra QA information
     XiPx = xi->MomXiX();
     XiPy = xi->MomXiY();
     XiPz = xi->MomXiZ();
     XiPt = TMath::Sqrt(XiPx*XiPx + XiPy*XiPy);
     XiPtot= TMath::Sqrt(XiPx*XiPx + XiPy*XiPy + XiPz*XiPz);
    
     bachPx = xi->MomBachX();
     bachPy = xi->MomBachY();
     bachPz = xi->MomBachZ();

     bachPt = TMath::Sqrt(bachPx*bachPx + bachPy*bachPy);
     bachPtot = TMath::Sqrt(bachPx*bachPx + bachPy*bachPy + bachPz*bachPz);
  
     V0toXiCosOfPointingAngle = xi->CosPointingAngle( xi->GetDecayVertexXi() );
    
     rapXi = xi->RapXi();
     rapOmega = xi->RapOmega();
     phi = xi->Phi();
     alphaXi = xi->AlphaXi();
     ptArmXi = xi->PtArmXi();

     distToVtxZBefore = posXi[2]-bestPrimaryVtxPos[2];
     distToVtxXYBefore
       = TMath::Sqrt((posXi[0] - bestPrimaryVtxPos[0])
		     *(posXi[0] - bestPrimaryVtxPos[0])
		     +(posXi[1] - bestPrimaryVtxPos[1])
		     *(posXi[1] - bestPrimaryVtxPos[1]));

     XiPAfter[0] = XiPx;
     XiPAfter[1] = XiPy;
     XiPAfter[2] = XiPz;
     
     //propagation to the best primary vertex to determine the momentum
     Propagate(bestPrimaryVtxPos, posXi, XiPAfter, b, xi->ChargeXi());
     distToVtxZAfter = posXi[2] - bestPrimaryVtxPos[2];
     distToVtxXYAfter = TMath::Sqrt((posXi[0] - bestPrimaryVtxPos[0])
				    *(posXi[0] - bestPrimaryVtxPos[0])
				    +(posXi[1] - bestPrimaryVtxPos[1])
				    *(posXi[1] - bestPrimaryVtxPos[1]));
     phiAfter = TMath::Pi() + TMath::ATan2(-XiPAfter[1],-XiPAfter[0]);
     
     fh1DistToVtxZAfter->Fill(distToVtxZAfter);
     fh1DistToVtxXYAfter->Fill(distToVtxXYAfter);
     fh2DistToVtxZBeforeVsAfter->Fill(distToVtxZBefore, distToVtxZAfter);
     fh2DistToVtxXYBeforeVsAfter->Fill(distToVtxXYBefore, distToVtxXYAfter);
     fh2PxBeforeVsAfter->Fill(XiPx, XiPAfter[0]);
     fh2PyBeforeVsAfter->Fill(XiPy, XiPAfter[1]);
     if(xi->ChargeXi()>0)
       fh2PhiPosBeforeVsAfter->Fill(phi, phiAfter);
     else if(xi->ChargeXi()<0)
       fh2PhiNegBeforeVsAfter->Fill(phi, phiAfter);

     if( (xi->ChargeXi() < 0 && isBachelorPionForTPC 
	  && isPosProtonForTPC && isNegPionForTPC)
	 || (xi->ChargeXi() > 0 && isBachelorPionForTPC 
	     && isNegProtonForTPC && isPosPionForTPC))
       {//Xi candidate
        //for default hypothesis
	 fh1Chi2Xi->Fill(chi2Xi);
	 fh1DCAXiDaughters->Fill(dcaXiDaughters);
	 fh1DCABachToPrimVertex->Fill(dcaBachToPrimaryVtxXi);
	 fh1XiCosOfPointingAngle->Fill(XiCosOfPointingAngle);
	 fh1XiRadius->Fill(XiRadius);
        
	 //V0
	 fh1MassLambda->Fill(invMassLambdaAsCascDghter);
	 fh1V0Chi2->Fill(V0Chi2Xi);
	 fh1DcaV0DaughtersXi->Fill(dcaV0DaughtersXi);
	 fh1V0CosOfPointingAngle->Fill(V0CosOfPointingAngleXi);
	 fh1V0Radius->Fill(V0RadiusXi);
	 fh1DcaV0ToPrimVertex->Fill(dcaV0ToPrimaryVtxXi);
	 fh1DCAPosToPrimVertex->Fill(dcaPosToPrimaryVtxXi);
	 fh1DCANegToPrimVertex->Fill(dcaNegToPrimaryVtxXi);
	 fh1ChargeXi->Fill(xi->ChargeXi());
	 fh1V0toXiCosOfPointingAngle->Fill(V0toXiCosOfPointingAngle);
       
	 if ( TMath::Abs( invMassXiMinus-1.3217 ) < 0.012 
	      || TMath::Abs( invMassXiPlus-1.3217 ) < 0.012) 
	   {// One InvMass should be different from 0             
	     fh1XiPt->Fill(XiPt);
	     fh1XiP->Fill(XiPtot);
            
	     fh1XiBachPt->Fill(bachPt);
	     fh1XiBachP->Fill(bachPtot);
	     fh1PhiXi->Fill( xi->Phi() );
	   }
	 fh2Armenteros->Fill(alphaXi, ptArmXi);
       }
  
     if(xi->ChargeXi()<0){
       fh1MassXiMinus->Fill(invMassXiMinus);
       fh1MassOmegaMinus->Fill(invMassOmegaMinus);
       fh1MassXi->Fill(invMassXiMinus);
       fh1MassOmega->Fill(invMassOmegaMinus);
      
       fh2MassLambdaVsMassXiMinus->Fill(invMassLambdaAsCascDghter,
					invMassXiMinus);
       fh2MassXiVsMassOmegaMinus->Fill(invMassXiMinus, 
				       invMassOmegaMinus);
       fh2XiRadiusVsMassXiMinus->Fill(XiRadius, invMassXiMinus);
       fh2XiRadiusVsMassOmegaMinus->Fill(XiRadius, invMassOmegaMinus);
     }
     
     if(xi->ChargeXi() > 0){
       fh1MassXiPlus->Fill(invMassXiPlus);
       fh1MassOmegaPlus->Fill(invMassOmegaPlus);
       fh1MassXi->Fill(invMassXiPlus);
       fh1MassOmega->Fill(invMassOmegaPlus);
      
       fh2MassLambdaVsMassXiPlus->Fill(invMassLambdaAsCascDghter, 
				       invMassXiPlus);
       fh2MassXiVsMassOmegaPlus->Fill(invMassXiPlus,
				      invMassOmegaPlus);
       fh2XiRadiusVsMassXiPlus->Fill(XiRadius, invMassXiPlus);
       fh2XiRadiusVsMassOmegaPlus->Fill(XiRadius, invMassOmegaPlus);
     }


     //      Double_t phiNew 
     //	= ( XiPAfter[0] == 0. && XiPAfter[1] == 0. ) 
     //	? 0.0 : TMath::ATan2(XiPAfter[1], XiPAfter[0]);
      Double_t phiV0 = phiAfter;
      phiV0 -= psiVZero;
      //if(phiV0 < 0) phiV0 += 2.*TMath::Pi();
      //if(phiV0 > TMath::Pi()) phiV0 -= TMath::Pi();


      Double_t phiV0A = phiAfter;
      phiV0A -= psiV0A;
      //if(phiV0A < 0) phiV0A += 2.*TMath::Pi();
      //if(phiV0A > TMath::Pi()) phiV0A -= TMath::Pi();

      Double_t phiV0C = phiAfter;
      phiV0C -= psiV0C;
      //if(phiV0C < 0) phiV0C += 2.*TMath::Pi();
      //if(phiV0C > TMath::Pi()) phiV0C -= TMath::Pi();
      
      //PID cuts with TPC cuts
      if(xi->ChargeXi() < 0){
	if(isPosProtonForTPC
	   && isNegPionForTPC){
	  if(isBachelorPionForTPC && TMath::Abs(rapXi) < 0.8){
	    //Xi
	    fh2MassVsPtXiMinus->Fill(XiPt, invMassXiMinus);
	    fh2MassVsPtXiAll->Fill(XiPt, invMassXiMinus);
	    fhXiRapidity->Fill(rapXi);
	    for(int r=0; r!=3; ++r) {
	      if(invMassXiMinus > fXiBands[r][0] 
		 && invMassXiMinus < fXiBands[r][1]){
		
		fProfXiV2PtV0A[r]->Fill(XiPt, TMath::Cos(2.*phiV0A));
                fProfXiSinePtV0A[r]->Fill(XiPt, TMath::Sin(2.*phiV0A));
		fProf2dXiV2PtV0A[r]->Fill(XiPt, psiV0A, TMath::Cos(2.*phiV0A));

                fProfXiV2PtV0C[r]->Fill(XiPt, TMath::Cos(2.*phiV0C));
                fProfXiSinePtV0C[r]->Fill(XiPt, TMath::Sin(2.*phiV0C));
		fProf2dXiV2PtV0C[r]->Fill(XiPt, psiV0C, TMath::Cos(2.*phiV0C));

		fProfXiV2Pt[r]->Fill(XiPt, TMath::Cos(2.*phiV0));
		fProfXiSinePt[r]->Fill(XiPt, TMath::Sin(2.*phiV0));
		fProf2dXiV2Pt[r]->Fill(XiPt, psiVZero, TMath::Cos(2.*phiV0));
	      }
	    }
	  }
	  
	  if(isBachelorKaonForTPC && TMath::Abs(rapOmega) < 0.8){
	    fh2MassVsPtOmegaMinus->Fill(XiPt, invMassOmegaMinus);
	    fh2MassVsPtOmegaAll->Fill(XiPt, invMassOmegaMinus); 
	    fhOmegaRapidity->Fill(rapOmega);
	    for(int r=0; r!=3; ++r) {
	      if(invMassOmegaMinus > fOmegaBands[r][0]
		 && invMassOmegaMinus < fOmegaBands[r][1]){
		fProfOmegaV2PtV0A[r]->Fill(XiPt, TMath::Cos(2.*phiV0A));
                fProfOmegaSinePtV0A[r]->Fill(XiPt, TMath::Sin(2.*phiV0A));
		fProf2dOmegaV2PtV0A[r]->Fill(XiPt, psiV0A, 
					     TMath::Cos(2.*phiV0A));

                fProfOmegaV2PtV0C[r]->Fill(XiPt, TMath::Cos(2.*phiV0C));
                fProfOmegaSinePtV0C[r]->Fill(XiPt, TMath::Sin(2.*phiV0C));
		fProf2dOmegaV2PtV0C[r]->Fill(XiPt, psiV0C, 
					     TMath::Cos(2.*phiV0C));

		fProfOmegaV2Pt[r]->Fill(XiPt, TMath::Cos(2.*phiV0));
                fProfOmegaSinePt[r]->Fill(XiPt, TMath::Sin(2.*phiV0));
		fProf2dOmegaV2Pt[r]->Fill(XiPt, psiVZero, 
					  TMath::Cos(2.*phiV0));
	      }
	    }
	  }
	}
      }//if charge < 0

      if(xi->ChargeXi() > 0){
	if(isNegProtonForTPC
	   && isPosPionForTPC){
	  if(isBachelorPionForTPC && TMath::Abs(rapXi) < 0.8){
	    //Xi
	    fh2MassVsPtXiPlus->Fill(XiPt, invMassXiPlus);
	    fh2MassVsPtXiAll->Fill(XiPt, invMassXiPlus);
	    fhXiRapidity->Fill(rapXi);
	    for(int r=0; r!=3; ++r) {
	      if(invMassXiPlus > fXiBands[r][0]
		 && invMassXiPlus < fXiBands[r][1]){
		fProfXiV2PtV0A[r]->Fill(XiPt, TMath::Cos(2.*phiV0A));
                fProfXiSinePtV0A[r]->Fill(XiPt, TMath::Sin(2.*phiV0A));
		fProf2dXiV2PtV0A[r]->Fill(XiPt, psiV0A, TMath::Cos(2.*phiV0A));

                fProfXiV2PtV0C[r]->Fill(XiPt, TMath::Cos(2.*phiV0C));
                fProfXiSinePtV0C[r]->Fill(XiPt, TMath::Sin(2.*phiV0C));
		fProf2dXiV2PtV0C[r]->Fill(XiPt, psiV0C, TMath::Cos(2.*phiV0C));

		fProfXiV2Pt[r]->Fill(XiPt, TMath::Cos(2.*phiV0));
                fProfXiSinePt[r]->Fill(XiPt, TMath::Sin(2.*phiV0));
		fProf2dXiV2Pt[r]->Fill(XiPt, psiVZero, TMath::Cos(2.*phiV0));
	      }
	    }
	  }

	  if(isBachelorKaonForTPC && TMath::Abs(rapOmega) < 0.8){
	    //Omega
	    fh2MassVsPtOmegaPlus->Fill(XiPt, invMassOmegaPlus);
	    fh2MassVsPtOmegaAll->Fill(XiPt, invMassOmegaPlus);
	    fhOmegaRapidity->Fill(rapOmega);
	    for(int r=0; r!=3; ++r) {
	      if(invMassOmegaPlus > fOmegaBands[r][0]
		 && invMassOmegaPlus < fOmegaBands[r][1]){
		fProfOmegaV2PtV0A[r]->Fill(XiPt, TMath::Cos(2.*phiV0A));
                fProfOmegaSinePtV0A[r]->Fill(XiPt, TMath::Sin(2.*phiV0A));
		fProf2dOmegaV2PtV0A[r]->Fill(XiPt, psiV0A,
                                             TMath::Cos(2.*phiV0A));


                fProfOmegaV2PtV0C[r]->Fill(XiPt, TMath::Cos(2.*phiV0C));
                fProfOmegaSinePtV0C[r]->Fill(XiPt, TMath::Sin(2.*phiV0C));
		fProf2dOmegaV2PtV0C[r]->Fill(XiPt, psiV0C,
                                             TMath::Cos(2.*phiV0C));

		fProfOmegaV2Pt[r]->Fill(XiPt, TMath::Cos(2.*phiV0));
                fProfOmegaSinePt[r]->Fill(XiPt, TMath::Sin(2.*phiV0));
		fProf2dOmegaV2Pt[r]->Fill(XiPt, psiVZero,
                                          TMath::Cos(2.*phiV0));
	      }
	    }
	  }
	}
      }// if charge > 0

  }//end of cascade loop

  PostData(1, fHistList);
      
}


AliAnalysisTaskFlowEPCascade::~AliAnalysisTaskFlowEPCascade(){
  if(fHistList) delete fHistList;
}

void AliAnalysisTaskFlowEPCascade::Terminate(Option_t *) 
{
  /*
    fHistList = dynamic_cast<TList*> (GetOutputData(1));
    
    if (!fHistList) {
    printf("ERROR: Output tree not available\n");
    return;
    }
  */
}      


void AliAnalysisTaskFlowEPCascade::Propagate(Double_t vv[3],
					     Double_t x[3],
					     Double_t p[3],
					     Double_t bz,
					     Short_t sign){
  //Propagation to the primary vertex to determine the px and py 
  //x, p are the position and momentum as input and output
  //bz is the magnetic field along z direction
  //sign is the charge of particle for propagation
  
  Double_t pp = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  Double_t len = (vv[2]-x[2])*pp/p[2];
  Double_t a = -kB2C*bz*sign;

  Double_t rho = a/pp;
  x[0] += p[0]*TMath::Sin(rho*len)/a - p[1]*(1-TMath::Cos(rho*len))/a;
  x[1] += p[1]*TMath::Sin(rho*len)/a + p[0]*(1-TMath::Cos(rho*len))/a;
  x[2] += p[2]*len/pp;

  Double_t p0=p[0];
  p[0] = p0  *TMath::Cos(rho*len) - p[1]*TMath::Sin(rho*len);
  p[1] = p[1]*TMath::Cos(rho*len) + p0  *TMath::Sin(rho*len);
}


void 
AliAnalysisTaskFlowEPCascade::OpenInfoCalbration(Int_t run){

  AliOADBContainer *cont = (AliOADBContainer*) fOADB->Get("hMultV0BefCorr");
  if(!cont){
    printf("OADB object hMultV0BefCorr is not available in the file\n");
    return;
  }

  if(!(cont->GetObject(run))){
    printf("OADB object hMultV0BefCorr is not available for run %i (used run 137366)\n",run);
    run = 137366;
  }
  fMultV0 = ((TH2F *) cont->GetObject(run))->ProfileX();

  TF1 *fpol0 = new TF1("fpol0","pol0");
  fMultV0->Fit(fpol0,"","",0,32);
  fV0Cpol = fpol0->GetParameter(0);
  fMultV0->Fit(fpol0,"","",32,64);
  fV0Apol = fpol0->GetParameter(0);

  for(Int_t iside=0;iside<2;iside++){
    for(Int_t icoord=0;icoord<2;icoord++){
      for(Int_t i=0;i  < 9;i++){
	char namecont[100];
	if(iside==0 && icoord==0)
	  snprintf(namecont,100,"hQxc2_%i", i);
	else if(iside==1 && icoord==0)
	  snprintf(namecont,100,"hQxa2_%i", i);
	else if(iside==0 && icoord==1)
	  snprintf(namecont,100,"hQyc2_%i", i);
	else if(iside==1 && icoord==1)
	  snprintf(namecont,100,"hQya2_%i", i);
      
	cont = (AliOADBContainer*) fOADB->Get(namecont);
	if(!cont){ 
	  printf("OADB object %s is not available in the file\n",namecont);
	  return;
	}
	
	if(!(cont->GetObject(run))){
	  printf("OADB object %s is not available for run %i (used run 137366)\n",namecont,run);
	  run = 137366;
	}
	fMeanQ[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetMean();
	fWidthQ[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetRMS();
      }
    }
  }
}
