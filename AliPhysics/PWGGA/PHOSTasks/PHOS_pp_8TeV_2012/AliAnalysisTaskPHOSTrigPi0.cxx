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
#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TRefArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom.h"
#include "THashList.h"
#include "TRandom3.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliMultiEventInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAnalysisTaskSE.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "TGeoManager.h"
#include "AliPHOSEsdCluster.h"
#include "AliPHOSAodCluster.h"
#include "AliPHOSCalibData.h"
#include "AliESDEvent.h"
#include "AliESDCaloCells.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliAODCaloCluster.h"
#include "AliAODJet.h"
#include "AliESDtrackCuts.h"
#include "AliOADBContainer.h"
#include "AliGenCocktailEventHeader.h"
#include "TGraphAsymmErrors.h"
#include "AliVCaloTrigger.h"
#include "AliAODCaloTrigger.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "AliPHOSEmcCalibData.h"
#include "AliAnalysisTaskPHOSTrigPi0.h"
#include "AliAODVZERO.h"
#include "AliAODTracklets.h"
#include "AliPIDResponse.h"
// #include "AliCaloTriggerSimulator.cxx"
#include "AliCaloTriggerSimulator.h"
ClassImp(AliAnalysisTaskPHOSTrigPi0);
//________________________________________________________________________
AliAnalysisTaskPHOSTrigPi0::AliAnalysisTaskPHOSTrigPi0(const char *name)
  :AliAnalysisTaskSE(name),
   fPHOSGeo(NULL),
   fCalibDataEmc(NULL),
   fPHOSCalibData(NULL),
   fCDBstorage(NULL),
   fPIDResponse(NULL),
   AODMCTrackArray(NULL),
   fSPDMultiCorrUnit(NULL),
   fSPDMultiCorrGap1(NULL),
   fSPDMultiCorrGap2(NULL),
   fV0ACMultiCorr(NULL),
   fVertexVector(),
   fRunNumber(0),
   fPeriod(0),
   fIsFiredTrig(""),
   fIsMC(false),
   fIsPileUp(false),
   fIsVtxOut10cm(false),
   fIsVtxNoCont(false),
   f0PH0Event(false),
   f0PH0Event_HighEnergy(false),
   fTOFcut0PH0Event(false),
   fTOFcut0PH0Event_HighEnergy(false),
   fTOFcut(25E-9),
   fMinNCell(3),
   fMinEene(0.3),
   fDispCut(2.5),
   fMultiBinSPDUnitRap(0),
   fMultiBinSPDGapRap1(0),
   fMultiBinSPDGapRap2(0),
   fMultiBinV0AC(0),
   fZvtxBin(0),
   fTriggerThreshold(384),
   fRefSPDMultiUnit(0),
   fRefSPDMultiGap1(0),
   fRefSPDMultiGap2(0),
   fRefV0ACMulti(0),
   fMCTracksUnit(0),
   fMCTracksGap1(0),
   fMCTracksGap2(0),
   fMCTracksV0AC(0),
   fHistTotalEvents(NULL),
   fHistAnalyzedEvents(NULL),
   fHistAnalyzed0PH0Events(NULL),
   fHistPUEvents(NULL),
   fHistPhysSelectionEvents(NULL),
   fHistSPDTrackletsUnitRap(NULL),
   fHistSPDTrackletsGapRap1(NULL),
   fHistSPDTrackletsGapRap2(NULL),
   fHistCorrectedSPDTrackletsUnitRap(NULL),
   fHistCorrectedSPDTrackletsGapRap1(NULL),
   fHistCorrectedSPDTrackletsGapRap2(NULL),
   fHistCorrelationMCTrackSPDTrackletsUnitRap(NULL),
   fHistCorrelationMCTrackSPDTrackletsGapRap1(NULL),
   fHistCorrelationMCTrackSPDTrackletsGapRap2(NULL),
   fHistCorrelationMCTrackCorrectedSPDTrackletsUnitRap(NULL),
   fHistCorrelationMCTrackCorrectedSPDTrackletsGapRap1(NULL),
   fHistCorrelationMCTrackCorrectedSPDTrackletsGapRap2(NULL),
   fHistCorrelationMCTrackV0AC(NULL),
   fHistCorrelationMCTrackCorrectedV0AC(NULL),
   fHistCorrelationSPDTrackletsV0AC(NULL),
   fHistCorrelationCorrectedSPDTrackletsCorrectedV0AC(NULL),
   fHistV0AMulti(NULL),
   fHistV0CMulti(NULL),
   fHistV0ACMulti(NULL),
   fHistCorrectedV0ACMulti(NULL),
   fEneSmearMeanMC(NULL),
   fEneSmearSigmaMC(NULL),
   fHistTrackCutX(NULL),
   fHistTrackCutZ(NULL),
   fHistGammaUnitRap(NULL),
   fHistGammaPi0UnitRap(NULL),
   fHistGammaEtaUnitRap(NULL),
   fHistGammaOmegaUnitRap(NULL),
   fHistGammaEtaPrimeUnitRap(NULL),
   fHistGammaPhiUnitRap(NULL),
   fHistGammaRhoUnitRap(NULL),
   fHistGammaOtherUnitRap(NULL),
   fHistGammaAccept(NULL),
   fHistGammaPi0Accept(NULL),
   fHistGammaEtaAccept(NULL),
   fHistGammaOmegaAccept(NULL),
   fHistGammaEtaPrimeAccept(NULL),
   fHistGammaPhiAccept(NULL),
   fHistGammaRhoAccept(NULL),
   fHistGammaOtherAccept(NULL),
   fHistPHOSAcceptance(NULL),
   fHistRadiusPi0(NULL),
   fHistRadius2DPi0(NULL),
   fHistRadius2DRecPi0(NULL),
   fHistZvertexPosition(NULL),
   fHistV0Timing(NULL),
   fHistV0TimingPhysSel(NULL),
   fHistdEdxElesigma(NULL),
   fHistdEventClassUnitRap(NULL),
   fHistdEventClassGapRap1(NULL),
   fHistdEventClassGapRap2(NULL),
   fHistdEventClassV0AC(NULL),
   fHistdEventClassPHIUnitRap(NULL),
   fHistdEventClassPHIGapRap1(NULL),
   fHistdEventClassPHIGapRap2(NULL),
   fHistdEventClassPHIV0AC(NULL),
   fMeanSPDUnitRap(0),
   fMeanSPDGapRap1(0),
   fMeanSPDGapRap2(0),
   fMeanV0AC(0),
   fRandom(new TRandom3(0)),
   fHistClustEneAll(0),
   fHistClustEneAllTOFcut1(0),
   fHistClustEneAllTOFcut2(0),
   fHistClustEneAllTOFcut3(0),
   fHistClustEneAllTOFcut4(0),
   fHistClustEneNcellCut(0),
   fHistClustEneNcellCutTOFcut1(0),
   fHistClustEneNcellCutTOFcut2(0),
   fHistClustEneNcellCutTOFcut3(0),
   fHistClustEneNcellCutTOFcut4(0),
   fHistGoodTRUClustEneAll(0),
   fHistGoodTRUClustEneAllTOFcut1(0),
   fHistGoodTRUClustEneAllTOFcut2(0),
   fHistGoodTRUClustEneAllTOFcut3(0),
   fHistGoodTRUClustEneAllTOFcut4(0),
   fHistGoodTRUClustEneNcellCut(0),
   fHistGoodTRUClustEneNcellCutTOFcut1(0),
   fHistGoodTRUClustEneNcellCutTOFcut2(0),
   fHistGoodTRUClustEneNcellCutTOFcut3(0),
   fHistGoodTRUClustEneNcellCutTOFcut4(0),
   fHist0PH0ClustEneAll(0),
   fHist0PH0ClustEneAllTOFcut1(0),
   fHist0PH0ClustEneAllTOFcut2(0),
   fHist0PH0ClustEneAllTOFcut3(0),
   fHist0PH0ClustEneAllTOFcut4(0),
   fHist0PH0ClustEneNcellCut(0),
   fHist0PH0ClustEneNcellCutTOFcut1(0),
   fHist0PH0ClustEneNcellCutTOFcut2(0),
   fHist0PH0ClustEneNcellCutTOFcut3(0),
   fHist0PH0ClustEneNcellCutTOFcut4(0),
   fHistGoodTRU0PH0ClustEneAll(0),
   fHistGoodTRU0PH0ClustEneAllTOFcut1(0),
   fHistGoodTRU0PH0ClustEneAllTOFcut2(0),
   fHistGoodTRU0PH0ClustEneAllTOFcut3(0),
   fHistGoodTRU0PH0ClustEneAllTOFcut4(0),
   fHistGoodTRU0PH0ClustEneNcellCut(0),
   fHistGoodTRU0PH0ClustEneNcellCutTOFcut1(0),
   fHistGoodTRU0PH0ClustEneNcellCutTOFcut2(0),
   fHistGoodTRU0PH0ClustEneNcellCutTOFcut3(0),
   fHistGoodTRU0PH0ClustEneNcellCutTOFcut4(0),
   fHistCellTimeHG(NULL),
   fHistCellTimeLG(NULL),
   fHistCellTimeHG_HE(NULL),
   fHistCellTimeLG_HE(NULL),
   fHist0PH0CellTimeHG(NULL),
   fHist0PH0CellTimeLG(NULL),
   fHist0PH0CellTimeHG_HE(NULL),
   fHist0PH0CellTimeLG_HE(NULL),

   fTrigData(NULL)
{
  for(Int_t iBuffer=0; iBuffer<10; ++iBuffer) fAOD[iBuffer] = NULL;
  
  for(Int_t iMulti=0; iMulti<10; ++iMulti){
    
    fHistPi0UnitRapV0AC[iMulti]= NULL;
    fHistPi0AcceptV0AC[iMulti] = NULL;
    fHistPi0UnitRapSPDUnitRap[iMulti]= NULL;
    fHistPi0AcceptSPDUnitRap[iMulti] = NULL;
    fHistPi0UnitRapSPDGapRap1[iMulti]= NULL;
    fHistPi0AcceptSPDGapRap1[iMulti] = NULL;
    fHistPi0UnitRapSPDGapRap2[iMulti]= NULL;
    fHistPi0AcceptSPDGapRap2[iMulti] = NULL;
    
    fHistSPDMultiUnitRapV0ACEventClass[iMulti]=NULL;
    fHistSPDMultiGapRap1V0ACEventClass[iMulti]=NULL;
    fHistSPDMultiGapRap2V0ACEventClass[iMulti]=NULL;

    fHistV0ACMultiSPDUnitRapEventClass[iMulti]=NULL;
    fHistV0ACMultiSPDGapRap1EventClass[iMulti]=NULL;
    fHistV0ACMultiSPDGapRap2EventClass[iMulti]=NULL;

    for(Int_t iZvtx=0; iZvtx<4; ++iZvtx){
      fPHOSClusterList[iMulti][iZvtx] = NULL;
    }
  }

  for(Int_t iMod=0; iMod<3; ++iMod){
    
    fIsGoodMod[iMod] = true;
    fBadMap[iMod]    = NULL;
    fBadTile[iMod]   = NULL;

    for(Int_t iTRU=0; iTRU<8; ++iTRU){
      fIsGoodTRU[iMod][iTRU] = true;
    }
    
    fTimingCalibMapLG[iMod]  = NULL;
    fTimingCalibMapHG[iMod]  = NULL;
    fHistClustOccupMap[iMod] = NULL;
    fHist0PH0ClustOccupMap[iMod] = NULL;
    
    for(Int_t iX=0; iX<64; ++iX){
      for(Int_t iZ=0; iZ<64; ++iZ){
	
	//fIsGoodCell[iMod][iX][iZ] = true;
	//fIsGoodTile[iMod][iX][iZ] = true;

      }
    }
  }
  
  for(Int_t iMod=0; iMod<4; ++iMod){
    fHistClustTOF[iMod]   = NULL;
    fHistClustTOFv1[iMod] = NULL;
    fHistClustTOFv2[iMod] = NULL;
    fHistClustTOFv3[iMod] = NULL;
    fHistClustTOFv4[iMod] = NULL;
    fHistClustTOFv5[iMod] = NULL;
    fHistClustTOFv6[iMod] = NULL;
  }

  for(Int_t iCand=0; iCand<4; ++iCand){
    for(Int_t iPID=0; iPID<4; ++iPID){
      for(Int_t iMod=0; iMod<4; ++iMod){
	
	fHistClustEne[iCand][iPID][iMod]    = NULL;
	fHistClustEneMIP[iCand][iPID][iMod] = NULL;
        fHistClustEp[iCand][iPID][iMod]     = NULL;
        fHistClustM02[iCand][iPID][iMod]    = NULL;
	fHistClustM20[iCand][iPID][iMod]    = NULL;
	fHistClustM02M20[iCand][iPID][iMod] = NULL;
	fHistClustDx[iCand][iPID][iMod]     = NULL;
	fHistClustDz[iCand][iPID][iMod]     = NULL;
	fHistClustDr[iCand][iPID][iMod]     = NULL;
	fHistClustNcells[iCand][iPID][iMod] = NULL;
	
	fHist0PH0ClustEne[iCand][iPID][iMod]    = NULL;
	fHist0PH0ClustEneMIP[iCand][iPID][iMod] = NULL;
        fHist0PH0ClustEp[iCand][iPID][iMod]     = NULL;
        fHist0PH0ClustM02[iCand][iPID][iMod]    = NULL;
	fHist0PH0ClustM20[iCand][iPID][iMod]    = NULL;
	fHist0PH0ClustM02M20[iCand][iPID][iMod] = NULL;
	fHist0PH0ClustDx[iCand][iPID][iMod]     = NULL;
	fHist0PH0ClustDz[iCand][iPID][iMod]     = NULL;
	fHist0PH0ClustDr[iCand][iPID][iMod]     = NULL;
	fHist0PH0ClustNcells[iCand][iPID][iMod] = NULL;
	
	for(Int_t iTrue=0; iTrue<10; ++iTrue){
	  
	  fHistTrueClustEne[iTrue][iCand][iPID][iMod]    = NULL;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][iMod] = NULL;
	  fHistTrueClustEp[iTrue][iCand][iPID][iMod]     = NULL;
          fHistTrueClustM02[iTrue][iCand][iPID][iMod]    = NULL;
	  fHistTrueClustM02M20[iTrue][iCand][iPID][iMod] = NULL;
	  fHistTrueClustM20[iTrue][iCand][iPID][iMod]    = NULL;
	  fHistTrueClustDx[iTrue][iCand][iPID][iMod]     = NULL;
	  fHistTrueClustDz[iTrue][iCand][iPID][iMod]     = NULL;
	  fHistTrueClustDr[iTrue][iCand][iPID][iMod]     = NULL;
	  fHistTrueClustNcells[iTrue][iCand][iPID][iMod] = NULL;
	  
	  fHistTrue0PH0ClustEne[iTrue][iCand][iPID][iMod]    = NULL;
	  fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][iMod] = NULL;
	  fHistTrue0PH0ClustEp[iTrue][iCand][iPID][iMod]     = NULL;
          fHistTrue0PH0ClustM02[iTrue][iCand][iPID][iMod]    = NULL;
	  fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][iMod] = NULL;
	  fHistTrue0PH0ClustM20[iTrue][iCand][iPID][iMod]    = NULL;
	  fHistTrue0PH0ClustDx[iTrue][iCand][iPID][iMod]     = NULL;
	  fHistTrue0PH0ClustDz[iTrue][iCand][iPID][iMod]     = NULL;
	  fHistTrue0PH0ClustDr[iTrue][iCand][iPID][iMod]     = NULL;
	  fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][iMod] = NULL;
	  
	}
      }
    }
  }

  for(Int_t iMod=0; iMod<4; ++iMod){
    for(Int_t iTrig=0; iTrig<2; ++iTrig){
      for(Int_t iTRU=0; iTRU<2; ++iTRU){
        for(Int_t iPID=0; iPID<4; ++iPID){
          for(Int_t iTOF=0; iTOF<3; ++iTOF){
            fHistMassTwoGammas[iTrig][iTRU][iTOF][iPID][iMod]    = NULL; 
	    fHistMixMassTwoGammas[iTrig][iTRU][iTOF][iPID][iMod] = NULL;
          }
          for(Int_t iTrue=0; iTrue<5; ++iTrue){
            fHistTrueMassTwoGammas[iTrig][iTRU][iPID][iTrue][iMod] = NULL;
          }

	}
      }
    }
  }
  
	  
  for(Int_t iTrig=0; iTrig<2; ++iTrig){
    for(Int_t iTRU=0; iTRU<2; ++iTRU){
      for(Int_t iMulti=0; iMulti<10; ++iMulti){

        for(Int_t iTOF=0; iTOF<3; ++iTOF){
	  fHistMassTwoGammasMultiUnitRap[iTrig][iTRU][iTOF][iMulti] = NULL;
	  fHistMassTwoGammasMultiGapRap1[iTrig][iTRU][iTOF][iMulti] = NULL;
	  fHistMassTwoGammasMultiGapRap2[iTrig][iTRU][iTOF][iMulti] = NULL;
	  fHistMassTwoGammasMultiV0AC[iTrig][iTRU][iTOF][iMulti]    = NULL;

	  fHistMixMassTwoGammasMultiUnitRap[iTrig][iTRU][iTOF][iMulti] = NULL;
	  fHistMixMassTwoGammasMultiGapRap1[iTrig][iTRU][iTOF][iMulti] = NULL;
	  fHistMixMassTwoGammasMultiGapRap2[iTrig][iTRU][iTOF][iMulti] = NULL;
	  fHistMixMassTwoGammasMultiV0AC[iTrig][iTRU][iTOF][iMulti]    = NULL;
        }
	
      }
    }
  }

  
  fEneSmearMeanMC = new TF1("fEneSmearMeanMC","([0]*exp([1]*x)+[2])",0.001,100);
  fEneSmearMeanMC->SetParameters(0.4,-3,0.001);
  fEneSmearSigmaMC= new TF1("fEneSmearSigmaMC","gaus",0,2);;
  fEneSmearSigmaMC->SetParameters(1,1,1);

  fHistTrackCutX = new TF1("fHistTrackCutX","[0]+[1]*exp([2]*x)",0,100);
  fHistTrackCutX->SetParameters(1.29654e+00,1.42330e+01,-7.99373e-01);
  fHistTrackCutZ = new TF1("fHistTrackCutZ","[0]+[1]*exp([2]*x)",0,100);
  fHistTrackCutZ->SetParameters(1.80107e+00,2.57932e+01,-9.72993e-01);
  
  fOutputContainer[0] = NULL;
  fOutputContainer[1] = NULL;
  fOutputContainer[2] = NULL;

  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());
  DefineOutput(3,TList::Class());
}
//___________________________________________________________________________
AliAnalysisTaskPHOSTrigPi0::~AliAnalysisTaskPHOSTrigPi0()
{
  
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::UserCreateOutputObjects()
{
  
  fCDBstorage = AliCDBManager::Instance()->GetStorage("alien://folder=/alice/data/2012/OCDB");
  
  for(Int_t i=0; i<3; ++i){
    if(fOutputContainer[i] != NULL){
      delete fOutputContainer[i];
    }
    fOutputContainer[i] = new THashList();
    fOutputContainer[i]->SetOwner(kTRUE);
  }
  
  fPHOSCalibData = new AliPHOSCalibData();
  fPHOSCalibData->SetName("PHOSCalibData");
  
  TString nameCand[] = {"_NoMatch","_Match","_EleCand","_HadronCand"};
  TString namePID[]  = {"","_CPV","_Disp","_DispCPV"};
  TString nameMod[]  = {"","_M1","_M2","_M3"};
  TString nameTrue[] = {"_Photon","_NotPhoton","_Ele","_Proton","_AntiProton","_Neutron","_AntiNeutron","_Charged","_Neutral","_ChargedPion"};

  fHistClustOccupMap[0] = new TH2F("fHistClustOccupMap_M1","",64,0.5,64.5,56,0.5,56.5);
  fHistClustOccupMap[1] = new TH2F("fHistClustOccupMap_M2","",64,0.5,64.5,56,0.5,56.5);
  fHistClustOccupMap[2] = new TH2F("fHistClustOccupMap_M3","",64,0.5,64.5,56,0.5,56.5);
  fOutputContainer[0]->Add(fHistClustOccupMap[0]);
  fOutputContainer[0]->Add(fHistClustOccupMap[1]);
  fOutputContainer[0]->Add(fHistClustOccupMap[2]);

  fHist0PH0ClustOccupMap[0] = new TH2F("fHist0PH0ClustOccupMap_M1","",64,0.5,64.5,56,0.5,56.5);
  fHist0PH0ClustOccupMap[1] = new TH2F("fHist0PH0ClustOccupMap_M2","",64,0.5,64.5,56,0.5,56.5);
  fHist0PH0ClustOccupMap[2] = new TH2F("fHist0PH0ClustOccupMap_M3","",64,0.5,64.5,56,0.5,56.5);
  fOutputContainer[0]->Add(fHist0PH0ClustOccupMap[0]);
  fOutputContainer[0]->Add(fHist0PH0ClustOccupMap[1]);
  fOutputContainer[0]->Add(fHist0PH0ClustOccupMap[2]);
  
  fHistPHOSAcceptance = new TH2F("fHistPHOSAcceptance","",100,-0.5,0.5,360,0,360);
  fOutputContainer[0]->Add(fHistPHOSAcceptance);

  fHistClustEneAll     = new TH1F("fHistClustEneAll","",400,0,40);
  fHistClustEneAllTOFcut1 = new TH1F("fHistClustEneAllTOFcut1","",400,0,40);
  fHistClustEneAllTOFcut2 = new TH1F("fHistClustEneAllTOFcut2","",400,0,40);
  fHistClustEneAllTOFcut3 = new TH1F("fHistClustEneAllTOFcut3","",400,0,40);
  fHistClustEneAllTOFcut4 = new TH1F("fHistClustEneAllTOFcut4","",400,0,40);

  fHistClustEneNcellCut     = new TH1F("fHistClustEneNcellCut","",400,0,40);
  fHistClustEneNcellCutTOFcut1 = new TH1F("fHistClustEneNcellCutTOFcut1","",400,0,40);
  fHistClustEneNcellCutTOFcut2 = new TH1F("fHistClustEneNcellCutTOFcut2","",400,0,40);
  fHistClustEneNcellCutTOFcut3 = new TH1F("fHistClustEneNcellCutTOFcut3","",400,0,40);
  fHistClustEneNcellCutTOFcut4 = new TH1F("fHistClustEneNcellCutTOFcut4","",400,0,40);

  fHistGoodTRUClustEneAll     = new TH1F("fHistGoodTRUClustEneAll","",400,0,40);
  fHistGoodTRUClustEneAllTOFcut1 = new TH1F("fHistGoodTRUClustEneAllTOFcut1","",400,0,40);
  fHistGoodTRUClustEneAllTOFcut2 = new TH1F("fHistGoodTRUClustEneAllTOFcut2","",400,0,40);
  fHistGoodTRUClustEneAllTOFcut3 = new TH1F("fHistGoodTRUClustEneAllTOFcut3","",400,0,40);
  fHistGoodTRUClustEneAllTOFcut4 = new TH1F("fHistGoodTRUClustEneAllTOFcut4","",400,0,40);

  fHistGoodTRUClustEneNcellCut     = new TH1F("fHistGoodTRUClustEneNcellCut","",400,0,40);
  fHistGoodTRUClustEneNcellCutTOFcut1 = new TH1F("fHistGoodTRUClustEneNcellCutTOFcut1","",400,0,40);
  fHistGoodTRUClustEneNcellCutTOFcut2 = new TH1F("fHistGoodTRUClustEneNcellCutTOFcut2","",400,0,40);
  fHistGoodTRUClustEneNcellCutTOFcut3 = new TH1F("fHistGoodTRUClustEneNcellCutTOFcut3","",400,0,40);
  fHistGoodTRUClustEneNcellCutTOFcut4 = new TH1F("fHistGoodTRUClustEneNcellCutTOFcut4","",400,0,40);

  fHist0PH0ClustEneAll     = new TH1F("fHist0PH0ClustEneAll","",400,0,40);
  fHist0PH0ClustEneAllTOFcut1 = new TH1F("fHist0PH0ClustEneAllTOFcut1","",400,0,40);
  fHist0PH0ClustEneAllTOFcut2 = new TH1F("fHist0PH0ClustEneAllTOFcut2","",400,0,40);
  fHist0PH0ClustEneAllTOFcut3 = new TH1F("fHist0PH0ClustEneAllTOFcut3","",400,0,40);
  fHist0PH0ClustEneAllTOFcut4 = new TH1F("fHist0PH0ClustEneAllTOFcut4","",400,0,40);

  fHist0PH0ClustEneNcellCut     = new TH1F("fHist0PH0ClustEneNcellCut","",400,0,40);
  fHist0PH0ClustEneNcellCutTOFcut1 = new TH1F("fHist0PH0ClustEneNcellCutTOFcut1","",400,0,40);
  fHist0PH0ClustEneNcellCutTOFcut2 = new TH1F("fHist0PH0ClustEneNcellCutTOFcut2","",400,0,40);
  fHist0PH0ClustEneNcellCutTOFcut3 = new TH1F("fHist0PH0ClustEneNcellCutTOFcut3","",400,0,40);
  fHist0PH0ClustEneNcellCutTOFcut4 = new TH1F("fHist0PH0ClustEneNcellCutTOFcut4","",400,0,40);

  fHistGoodTRU0PH0ClustEneAll     = new TH1F("fHistGoodTRU0PH0ClustEneAll","",400,0,40);
  fHistGoodTRU0PH0ClustEneAllTOFcut1 = new TH1F("fHistGoodTRU0PH0ClustEneAllTOFcut1","",400,0,40);
  fHistGoodTRU0PH0ClustEneAllTOFcut2 = new TH1F("fHistGoodTRU0PH0ClustEneAllTOFcut2","",400,0,40);
  fHistGoodTRU0PH0ClustEneAllTOFcut3 = new TH1F("fHistGoodTRU0PH0ClustEneAllTOFcut3","",400,0,40);
  fHistGoodTRU0PH0ClustEneAllTOFcut4 = new TH1F("fHistGoodTRU0PH0ClustEneAllTOFcut4","",400,0,40);

  fHistGoodTRU0PH0ClustEneNcellCut     = new TH1F("fHistGoodTRU0PH0ClustEneNcellCut","",400,0,40);
  fHistGoodTRU0PH0ClustEneNcellCutTOFcut1 = new TH1F("fHistGoodTRU0PH0ClustEneNcellCutTOFcut1","",400,0,40);
  fHistGoodTRU0PH0ClustEneNcellCutTOFcut2 = new TH1F("fHistGoodTRU0PH0ClustEneNcellCutTOFcut2","",400,0,40);
  fHistGoodTRU0PH0ClustEneNcellCutTOFcut3 = new TH1F("fHistGoodTRU0PH0ClustEneNcellCutTOFcut3","",400,0,40);
  fHistGoodTRU0PH0ClustEneNcellCutTOFcut4 = new TH1F("fHistGoodTRU0PH0ClustEneNcellCutTOFcut4","",400,0,40);
  
  fOutputContainer[0]->Add(fHistClustEneAll);
  fOutputContainer[0]->Add(fHistClustEneAllTOFcut1);
  fOutputContainer[0]->Add(fHistClustEneAllTOFcut2);
  fOutputContainer[0]->Add(fHistClustEneAllTOFcut3);
  fOutputContainer[0]->Add(fHistClustEneAllTOFcut4);
  fOutputContainer[0]->Add(fHistClustEneNcellCut);
  fOutputContainer[0]->Add(fHistClustEneNcellCutTOFcut1);
  fOutputContainer[0]->Add(fHistClustEneNcellCutTOFcut2);
  fOutputContainer[0]->Add(fHistClustEneNcellCutTOFcut3);
  fOutputContainer[0]->Add(fHistClustEneNcellCutTOFcut4);
  fOutputContainer[0]->Add(fHistGoodTRUClustEneAll);
  fOutputContainer[0]->Add(fHistGoodTRUClustEneAllTOFcut1);
  fOutputContainer[0]->Add(fHistGoodTRUClustEneAllTOFcut2);
  fOutputContainer[0]->Add(fHistGoodTRUClustEneAllTOFcut3);
  fOutputContainer[0]->Add(fHistGoodTRUClustEneAllTOFcut4);
  fOutputContainer[0]->Add(fHistGoodTRUClustEneNcellCut);
  fOutputContainer[0]->Add(fHistGoodTRUClustEneNcellCutTOFcut1);
  fOutputContainer[0]->Add(fHistGoodTRUClustEneNcellCutTOFcut2);
  fOutputContainer[0]->Add(fHistGoodTRUClustEneNcellCutTOFcut3);
  fOutputContainer[0]->Add(fHistGoodTRUClustEneNcellCutTOFcut4);
  fOutputContainer[0]->Add(fHist0PH0ClustEneAll);
  fOutputContainer[0]->Add(fHist0PH0ClustEneAllTOFcut1);
  fOutputContainer[0]->Add(fHist0PH0ClustEneAllTOFcut2);
  fOutputContainer[0]->Add(fHist0PH0ClustEneAllTOFcut3);
  fOutputContainer[0]->Add(fHist0PH0ClustEneAllTOFcut4);
  fOutputContainer[0]->Add(fHist0PH0ClustEneNcellCut);
  fOutputContainer[0]->Add(fHist0PH0ClustEneNcellCutTOFcut1);
  fOutputContainer[0]->Add(fHist0PH0ClustEneNcellCutTOFcut2);
  fOutputContainer[0]->Add(fHist0PH0ClustEneNcellCutTOFcut3);
  fOutputContainer[0]->Add(fHist0PH0ClustEneNcellCutTOFcut4);
  fOutputContainer[0]->Add(fHistGoodTRU0PH0ClustEneAll);
  fOutputContainer[0]->Add(fHistGoodTRU0PH0ClustEneAllTOFcut1);
  fOutputContainer[0]->Add(fHistGoodTRU0PH0ClustEneAllTOFcut2);
  fOutputContainer[0]->Add(fHistGoodTRU0PH0ClustEneAllTOFcut3);
  fOutputContainer[0]->Add(fHistGoodTRU0PH0ClustEneAllTOFcut4);
  fOutputContainer[0]->Add(fHistGoodTRU0PH0ClustEneNcellCut);
  fOutputContainer[0]->Add(fHistGoodTRU0PH0ClustEneNcellCutTOFcut1);
  fOutputContainer[0]->Add(fHistGoodTRU0PH0ClustEneNcellCutTOFcut2);
  fOutputContainer[0]->Add(fHistGoodTRU0PH0ClustEneNcellCutTOFcut3);
  fOutputContainer[0]->Add(fHistGoodTRU0PH0ClustEneNcellCutTOFcut4);

  fHistCellTimeHG = new TH2F("fHistCellTimeHG","",10752,0.5,10752.5,400,-200E-9,200E-9);
  fHistCellTimeLG = new TH2F("fHistCellTimeLG","",10752,0.5,10752.5,400,-200E-9,200E-9);
  fHistCellTimeHG_HE = new TH2F("fHistCellTimeHG_HE","",10752,0.5,10752.5,400,-200E-9,200E-9);
  fHistCellTimeLG_HE = new TH2F("fHistCellTimeLG_HE","",10752,0.5,10752.5,400,-200E-9,200E-9);
  fOutputContainer[0]->Add(fHistCellTimeHG);
  fOutputContainer[0]->Add(fHistCellTimeLG);
  fOutputContainer[0]->Add(fHistCellTimeHG_HE);
  fOutputContainer[0]->Add(fHistCellTimeLG_HE);

  fHist0PH0CellTimeHG = new TH2F("fHist0PH0CellTimeHG","",10752,0.5,10752.5,400,-200E-9,200E-9);
  fHist0PH0CellTimeLG = new TH2F("fHist0PH0CellTimeLG","",10752,0.5,10752.5,400,-200E-9,200E-9);
  fHist0PH0CellTimeHG_HE = new TH2F("fHist0PH0CellTimeHG_HE","",10752,0.5,10752.5,400,-200E-9,200E-9);
  fHist0PH0CellTimeLG_HE = new TH2F("fHist0PH0CellTimeLG_HE","",10752,0.5,10752.5,400,-200E-9,200E-9);
  fOutputContainer[0]->Add(fHist0PH0CellTimeHG);
  fOutputContainer[0]->Add(fHist0PH0CellTimeLG);
  fOutputContainer[0]->Add(fHist0PH0CellTimeHG_HE);
  fOutputContainer[0]->Add(fHist0PH0CellTimeLG_HE);
  
  for(Int_t iMod=0; iMod<4; ++iMod){
    
    

    if(iMod==2) continue;



    fHistClustTOF[iMod] = new TH2F("fHistClustTOF"+nameMod[iMod],"",400,0,40,400,-200E-9,200E-9);
    fHistClustTOFv1[iMod] = new TH2F("fHistClustTOFv1"+nameMod[iMod],"",400,0,40,400,-200E-9,200E-9);
    fHistClustTOFv2[iMod] = new TH2F("fHistClustTOFv2"+nameMod[iMod],"",400,0,40,400,-200E-9,200E-9);
    fHistClustTOFv3[iMod] = new TH2F("fHistClustTOFv3"+nameMod[iMod],"",400,0,40,400,-200E-9,200E-9);
    fHistClustTOFv4[iMod] = new TH2F("fHistClustTOFv4"+nameMod[iMod],"",400,0,40,400,-200E-9,200E-9);
    fHistClustTOFv5[iMod] = new TH2F("fHistClustTOFv5"+nameMod[iMod],"",400,0,40,400,-200E-9,200E-9);
    fHistClustTOFv6[iMod] = new TH2F("fHistClustTOFv6"+nameMod[iMod],"",400,0,40,400,-200E-9,200E-9);
    fOutputContainer[0]->Add(fHistClustTOF[iMod]);
    fOutputContainer[0]->Add(fHistClustTOFv1[iMod]);
    fOutputContainer[0]->Add(fHistClustTOFv2[iMod]);
    fOutputContainer[0]->Add(fHistClustTOFv3[iMod]);
    fOutputContainer[0]->Add(fHistClustTOFv4[iMod]);
    fOutputContainer[0]->Add(fHistClustTOFv5[iMod]);
    fOutputContainer[0]->Add(fHistClustTOFv6[iMod]);

    for(Int_t iCand=0; iCand<4; ++iCand){
      for(Int_t iPID=0; iPID<4; ++iPID){
        
	fHistClustEne[iCand][iPID][iMod]    = new TH1F("fHistClustEne"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40);
	fHistClustEneMIP[iCand][iPID][iMod] = new TH1F("fHistClustEneMIP"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",2000,0,2);
        fHistClustEp[iCand][iPID][iMod]     = new TH2F("fHistClustEp"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,40,0,2);
        fHistClustM02[iCand][iPID][iMod]    = new TH2F("fHistClustM02"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,5);
	fHistClustM20[iCand][iPID][iMod]    = new TH2F("fHistClustM20"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,5);
	fHistClustM02M20[iCand][iPID][iMod] = new TH2F("fHistClustM02M20"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",500,0,5,500,0,5);
	fHistClustDx[iCand][iPID][iMod]     = new TH2F("fHistClustDx"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,100,-20,20);
	fHistClustDz[iCand][iPID][iMod]     = new TH2F("fHistClustDz"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,100,-20,20);
	fHistClustDr[iCand][iPID][iMod]     = new TH2F("fHistClustDr"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,100);
	fHistClustNcells[iCand][iPID][iMod] = new TH2F("fHistClustNcells"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,50);

	fHist0PH0ClustEne[iCand][iPID][iMod]    = new TH1F("fHist0PH0ClustEne"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40);
	fHist0PH0ClustEneMIP[iCand][iPID][iMod] = new TH1F("fHist0PH0ClustEneMIP"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",2000,0,2);
        fHist0PH0ClustEp[iCand][iPID][iMod]     = new TH2F("fHist0PH0ClustEp"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,40,0,2);
        fHist0PH0ClustM02[iCand][iPID][iMod]    = new TH2F("fHist0PH0ClustM02"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,5);
	fHist0PH0ClustM20[iCand][iPID][iMod]    = new TH2F("fHist0PH0ClustM20"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,5);
	fHist0PH0ClustM02M20[iCand][iPID][iMod] = new TH2F("fHist0PH0ClustM02M20"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",500,0,5,500,0,5);
	fHist0PH0ClustDx[iCand][iPID][iMod]     = new TH2F("fHist0PH0ClustDx"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,100,-20,20);
	fHist0PH0ClustDz[iCand][iPID][iMod]     = new TH2F("fHist0PH0ClustDz"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,100,-20,20);
	fHist0PH0ClustDr[iCand][iPID][iMod]     = new TH2F("fHist0PH0ClustDr"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,100);
	fHist0PH0ClustNcells[iCand][iPID][iMod] = new TH2F("fHist0PH0ClustNcells"+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,50);

	


	if(iMod==0){
	  fOutputContainer[0]->Add(fHistClustEneMIP[iCand][iPID][iMod]);
	  fOutputContainer[0]->Add(fHistClustEne[iCand][iPID][iMod]);
	  fOutputContainer[0]->Add(fHistClustEp[iCand][iPID][iMod]);
	  fOutputContainer[0]->Add(fHistClustM02M20[iCand][iPID][iMod]);
	  //fOutputContainer[0]->Add(fHistClustM02[iCand][iPID][iMod]);
	  //fOutputContainer[0]->Add(fHistClustM20[iCand][iPID][iMod]);
	  fOutputContainer[0]->Add(fHistClustDx[iCand][iPID][iMod]);
	  fOutputContainer[0]->Add(fHistClustDz[iCand][iPID][iMod]);
	  fOutputContainer[0]->Add(fHistClustDr[iCand][iPID][iMod]);
	  fOutputContainer[0]->Add(fHistClustNcells[iCand][iPID][iMod]);
	  
	  //fOutputContainer[0]->Add(fHist0PH0ClustEneMIP[iCand][iPID][iMod]);
	  fOutputContainer[0]->Add(fHist0PH0ClustEne[iCand][iPID][iMod]);
	  fOutputContainer[0]->Add(fHist0PH0ClustEp[iCand][iPID][iMod]);
	  //fOutputContainer[0]->Add(fHist0PH0ClustM02M20[iCand][iPID][iMod]);
	  //fOutputContainer[0]->Add(fHist0PH0ClustM02[iCand][iPID][iMod]);
	  //fOutputContainer[0]->Add(fHist0PH0ClustM20[iCand][iPID][iMod]);
	  //fOutputContainer[0]->Add(fHist0PH0ClustDx[iCand][iPID][iMod]);
	  //fOutputContainer[0]->Add(fHist0PH0ClustDz[iCand][iPID][iMod]);
	  //fOutputContainer[0]->Add(fHist0PH0ClustDr[iCand][iPID][iMod]);
	  //fOutputContainer[0]->Add(fHist0PH0ClustNcells[iCand][iPID][iMod]);

	}	  

	for(Int_t iTrue=0; iTrue<10; ++iTrue){
	  
	  fHistTrueClustEne[iTrue][iCand][iPID][iMod]    = new TH1F("fHistTrueClustEne"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40);
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][iMod] = new TH1F("fHistTrueClustEneMIP"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",2000,0,2);
	  fHistTrueClustEp[iTrue][iCand][iPID][iMod]     = new TH2F("fHistTrueClustEp"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,40,0,2);
	  fHistTrueClustM02[iTrue][iCand][iPID][iMod]    = new TH2F("fHistTrueClustM02"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,5);
	  fHistTrueClustM20[iTrue][iCand][iPID][iMod]    = new TH2F("fHistTrueClustM20"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,5);
	  fHistTrueClustM02M20[iTrue][iCand][iPID][iMod] = new TH2F("fHistTrueClustM02M20"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",500,0,5,500,0,5);
	  fHistTrueClustDx[iTrue][iCand][iPID][iMod]     = new TH2F("fHistTrueClustDx"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,100,-20,20);
	  fHistTrueClustDz[iTrue][iCand][iPID][iMod]     = new TH2F("fHistTrueClustDz"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,100,-20,20);
	  fHistTrueClustDr[iTrue][iCand][iPID][iMod]     = new TH2F("fHistTrueClustDr"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,100);
	  fHistTrueClustNcells[iTrue][iCand][iPID][iMod] = new TH2F("fHistTrueClustNcells"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,50);

	  fHistTrue0PH0ClustEne[iTrue][iCand][iPID][iMod]    = new TH1F("fHistTrue0PH0ClustEne"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40);
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][iMod] = new TH1F("fHistTrue0PH0ClustEneMIP"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",2000,0,2);
	  fHistTrue0PH0ClustEp[iTrue][iCand][iPID][iMod]     = new TH2F("fHistTrue0PH0ClustEp"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,40,0,2);
	  //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][iMod]    = new TH2F("fHistTrue0PH0ClustM02"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,5);
	  //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][iMod]    = new TH2F("fHistTrue0PH0ClustM20"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,5);
	  //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][iMod] = new TH2F("fHistTrue0PH0ClustM02M20"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",500,0,5,500,0,5);
	  //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][iMod]     = new TH2F("fHistTrue0PH0ClustDx"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,100,-20,20);
	  //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][iMod]     = new TH2F("fHistTrue0PH0ClustDz"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,100,-20,20);
	  //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][iMod]     = new TH2F("fHistTrue0PH0ClustDr"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,100);
	  //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][iMod] = new TH2F("fHistTrue0PH0ClustNcells"+nameTrue[iTrue]+nameCand[iCand]+namePID[iPID]+nameMod[iMod],"",400,0,40,50,0,50);

	  if(iMod==0){
	    fOutputContainer[0]->Add(fHistTrueClustEneMIP[iTrue][iCand][iPID][iMod]);
	    fOutputContainer[0]->Add(fHistTrueClustEne[iTrue][iCand][iPID][iMod]);
	    fOutputContainer[0]->Add(fHistTrueClustEp[iTrue][iCand][iPID][iMod]);
	    //fOutputContainer[0]->Add(fHistTrueClustM02[iTrue][iCand][iPID][iMod]);
	    //fOutputContainer[0]->Add(fHistTrueClustM20[iTrue][iCand][iPID][iMod]);
	    fOutputContainer[0]->Add(fHistTrueClustM02M20[iTrue][iCand][iPID][iMod]);
	    fOutputContainer[0]->Add(fHistTrueClustDx[iTrue][iCand][iPID][iMod]);
	    fOutputContainer[0]->Add(fHistTrueClustDz[iTrue][iCand][iPID][iMod]);
	    fOutputContainer[0]->Add(fHistTrueClustDr[iTrue][iCand][iPID][iMod]);
	    fOutputContainer[0]->Add(fHistTrueClustNcells[iTrue][iCand][iPID][iMod]);

	    //fOutputContainer[0]->Add(fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][iMod]);
	    fOutputContainer[0]->Add(fHistTrue0PH0ClustEne[iTrue][iCand][iPID][iMod]);
	    fOutputContainer[0]->Add(fHistTrue0PH0ClustEp[iTrue][iCand][iPID][iMod]);
	    //fOutputContainer[0]->Add(fHistTrue0PH0ClustM02[iTrue][iCand][iPID][iMod]);
	    //fOutputContainer[0]->Add(fHistTrue0PH0ClustM20[iTrue][iCand][iPID][iMod]);
	    //fOutputContainer[0]->Add(fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][iMod]);
	    //fOutputContainer[0]->Add(fHistTrue0PH0ClustDx[iTrue][iCand][iPID][iMod]);
	    //fOutputContainer[0]->Add(fHistTrue0PH0ClustDz[iTrue][iCand][iPID][iMod]);
	    //fOutputContainer[0]->Add(fHistTrue0PH0ClustDr[iTrue][iCand][iPID][iMod]);
	    //fOutputContainer[0]->Add(fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][iMod]);

	  }

	}
	
      }
    }
  }

  
  fHistGammaUnitRap         = new TH1F("fHistGammaUnitRap","",400,0,40);
  fHistGammaPi0UnitRap      = new TH1F("fHistGammaPi0UnitRap","",400,0,40);
  fHistGammaEtaUnitRap      = new TH1F("fHistGammaEtaUnitRap","",400,0,40);
  fHistGammaOmegaUnitRap    = new TH1F("fHistGammaOmegaUnitRap","",400,0,40);
  fHistGammaEtaPrimeUnitRap = new TH1F("fHistGammaEtaPrimeUnitRap","",400,0,40);
  fHistGammaPhiUnitRap      = new TH1F("fHistGammaPhiUnitRap","",400,0,40);
  fHistGammaRhoUnitRap      = new TH1F("fHistGammaRhoUnitRap","",400,0,40);
  fHistGammaOtherUnitRap      = new TH1F("fHistGammaOtherUnitRap","",400,0,40);
  if(fIsMC){
    fOutputContainer[0]->Add(fHistGammaUnitRap);
    fOutputContainer[0]->Add(fHistGammaPi0UnitRap);
    fOutputContainer[0]->Add(fHistGammaEtaUnitRap);
    fOutputContainer[0]->Add(fHistGammaOmegaUnitRap);
    fOutputContainer[0]->Add(fHistGammaEtaPrimeUnitRap);
    fOutputContainer[0]->Add(fHistGammaPhiUnitRap);
    fOutputContainer[0]->Add(fHistGammaRhoUnitRap);
    fOutputContainer[0]->Add(fHistGammaOtherUnitRap);
  }
  fHistGammaAccept         = new TH1F("fHistGammaAccept","",400,0,40);
  fHistGammaPi0Accept      = new TH1F("fHistGammaPi0Accept","",400,0,40);
  fHistGammaEtaAccept      = new TH1F("fHistGammaEtaAccept","",400,0,40);
  fHistGammaOmegaAccept    = new TH1F("fHistGammaOmegaAccept","",400,0,40);
  fHistGammaEtaPrimeAccept = new TH1F("fHistGammaEtaPrimeAccept","",400,0,40);
  fHistGammaPhiAccept      = new TH1F("fHistGammaPhiAccept","",400,0,40);
  fHistGammaRhoAccept      = new TH1F("fHistGammaRhoAccept","",400,0,40);
  fHistGammaOtherAccept      = new TH1F("fHistGammaOtherAccept","",400,0,40);
  if(fIsMC){
    fOutputContainer[0]->Add(fHistGammaAccept);
    fOutputContainer[0]->Add(fHistGammaPi0Accept);
    fOutputContainer[0]->Add(fHistGammaEtaAccept);
    fOutputContainer[0]->Add(fHistGammaOmegaAccept);
    fOutputContainer[0]->Add(fHistGammaEtaPrimeAccept);
    fOutputContainer[0]->Add(fHistGammaPhiAccept);
    fOutputContainer[0]->Add(fHistGammaRhoAccept);
    fOutputContainer[0]->Add(fHistGammaOtherAccept);
  }
  
  /*
    TString namePID[]  = {"","_CPV","_Disp","_DispCPV"};
    TString nameMod[]  = {"","_M1","_M2","_M3"};
    TString nameTrue[] = {"_Photon","_NotPhoton","_Ele","_Proton","_AntiProton","_Neutron","_AntiNeutron","_Charged","_Neutral","_ChargedPion"};
  */
  
  TString nameTrig[]           = {"","_PHI"};
  TString nameTRU[]            = {"","_GoodTRU"};
  TString nameTOF[]            = {"","_TOFcut1","_TOFcut2"};
  TString nameTrue_TwoGammas[] =  {"","_Pi0","_SameMother","_Pi0FromK0s","_Pi0FromMaterial"};

  for(Int_t iTrig=0; iTrig<2; ++iTrig){
    for(Int_t iTRU=0; iTRU<2; ++iTRU){
      for(Int_t iPID=0; iPID<4; ++iPID){
	for(Int_t iMod=0; iMod<4; ++iMod){
	  
	  if(iMod==2) continue;

	  for(Int_t iTOF=0; iTOF<3; ++iTOF){	
	    fHistMassTwoGammas[iTrig][iTRU][iTOF][iPID][iMod]
	      = new TH2F("fHistMassTwoGammas"+nameTrig[iTrig]+nameTRU[iTRU]+nameTOF[iTOF]+namePID[iPID]+nameMod[iMod],"",500,0,0.5,400,0,40);
	    fOutputContainer[1]->Add(fHistMassTwoGammas[iTrig][iTRU][iTOF][iPID][iMod]);
	    
	    fHistMixMassTwoGammas[iTrig][iTRU][iTOF][iPID][iMod]
	      = new TH2F("fHistMixMassTwoGammas"+nameTrig[iTrig]+nameTRU[iTRU]+nameTOF[iTOF]+namePID[iPID]+nameMod[iMod],"",500,0,0.5,400,0,40);
	    fOutputContainer[1]->Add(fHistMixMassTwoGammas[iTrig][iTRU][iTOF][iPID][iMod]);
	  
	  }
	  
	  if(fIsMC){
	    for(Int_t iTrue=0; iTrue<5; ++iTrue){
	      fHistTrueMassTwoGammas[iTrig][iTRU][iPID][iTrue][iMod] 
		= new TH2F("fHistTrueMassTwoGammas"+nameTrig[iTrig]+nameTRU[iTRU]+nameTrue_TwoGammas[iTrue]+namePID[iPID]+nameMod[iMod],"",500,0,0.5,400,0,40);
	      fOutputContainer[1]->Add(fHistTrueMassTwoGammas[iTrig][iTRU][iPID][iTrue][iMod]);
	    }
	  }
	  
	}
      }
    }
  }
  
  TString nameMulti[] = {"_Multi0","_Multi1","_Multi2","_Multi3","_Multi4","_Multi5","_Multi6","_Multi7","_Multi8","_Multi9"};

  for(Int_t iMulti=0; iMulti<10; ++iMulti){
    
    fHistSPDMultiUnitRapV0ACEventClass[iMulti] = new TH1F("fHistSPDMultiUnitRapV0ACEventClass"+nameMulti[iMulti],"",150,0,150);
    fHistSPDMultiGapRap1V0ACEventClass[iMulti] = new TH1F("fHistSPDMultiGapRap1V0ACEventClass"+nameMulti[iMulti],"",150,0,150);
    fHistSPDMultiGapRap2V0ACEventClass[iMulti] = new TH1F("fHistSPDMultiGapRap2V0ACEventClass"+nameMulti[iMulti],"",150,0,150);
    
    fHistV0ACMultiSPDUnitRapEventClass[iMulti] = new TH1F("fHistV0ACMultiSPDUnitRapEventClass"+nameMulti[iMulti],"",1000,0,2000);
    fHistV0ACMultiSPDGapRap1EventClass[iMulti] = new TH1F("fHistV0ACMultiSPDGapRap1EventClass"+nameMulti[iMulti],"",1000,0,2000);
    fHistV0ACMultiSPDGapRap2EventClass[iMulti] = new TH1F("fHistV0ACMultiSPDGapRap2EventClass"+nameMulti[iMulti],"",1000,0,2000);
    
    fOutputContainer[2]->Add(fHistSPDMultiUnitRapV0ACEventClass[iMulti]);
    fOutputContainer[2]->Add(fHistSPDMultiGapRap1V0ACEventClass[iMulti]);
    fOutputContainer[2]->Add(fHistSPDMultiGapRap2V0ACEventClass[iMulti]);

    fOutputContainer[2]->Add(fHistV0ACMultiSPDUnitRapEventClass[iMulti]);
    fOutputContainer[2]->Add(fHistV0ACMultiSPDGapRap1EventClass[iMulti]);
    fOutputContainer[2]->Add(fHistV0ACMultiSPDGapRap2EventClass[iMulti]);
    
    fHistPi0UnitRapV0AC[iMulti]       = new TH1F("fHistPi0UnitRapV0AC"+nameMulti[iMulti],"",300,0,30);
    fHistPi0AcceptV0AC[iMulti]        = new TH1F("fHistPi0AcceptV0AC"+nameMulti[iMulti],"",300,0,30);
    fHistPi0UnitRapSPDUnitRap[iMulti] = new TH1F("fHistPi0UnitRapSPDUnitRap"+nameMulti[iMulti],"",300,0,30);
    fHistPi0AcceptSPDUnitRap[iMulti]  = new TH1F("fHistPi0AcceptSPDUnitRap"+nameMulti[iMulti],"",300,0,30);
    fHistPi0UnitRapSPDGapRap1[iMulti] = new TH1F("fHistPi0UnitRapSPDGapRap1"+nameMulti[iMulti],"",300,0,30);
    fHistPi0AcceptSPDGapRap1[iMulti]  = new TH1F("fHistPi0AcceptSPDGapRap1"+nameMulti[iMulti],"",300,0,30);
    fHistPi0UnitRapSPDGapRap2[iMulti] = new TH1F("fHistPi0UnitRapSPDGapRap2"+nameMulti[iMulti],"",300,0,30);
    fHistPi0AcceptSPDGapRap2[iMulti]  = new TH1F("fHistPi0AcceptSPDGapRap2"+nameMulti[iMulti],"",300,0,30);
    
    if(fIsMC){
      fOutputContainer[1]->Add(fHistPi0UnitRapV0AC[iMulti]);
      fOutputContainer[1]->Add(fHistPi0AcceptV0AC[iMulti]);
      fOutputContainer[1]->Add(fHistPi0UnitRapSPDUnitRap[iMulti]);
      fOutputContainer[1]->Add(fHistPi0AcceptSPDUnitRap[iMulti]);
      fOutputContainer[1]->Add(fHistPi0UnitRapSPDGapRap1[iMulti]);
      fOutputContainer[1]->Add(fHistPi0AcceptSPDGapRap1[iMulti]);
      fOutputContainer[1]->Add(fHistPi0UnitRapSPDGapRap2[iMulti]);
      fOutputContainer[1]->Add(fHistPi0AcceptSPDGapRap2[iMulti]);
    }
    
    for(Int_t iTrig=0; iTrig<2; ++iTrig){
      for(Int_t iTRU=0; iTRU<2; ++iTRU){
	for(Int_t iTOF=0; iTOF<3; ++iTOF){	
	  fHistMassTwoGammasMultiUnitRap[iTrig][iTRU][iTOF][iMulti] = new TH2F("fHistMassTwoGammasUnitRap"+nameTrig[iTrig]+nameTRU[iTRU]+nameTOF[iTOF]+nameMulti[iMulti],"",500,0,0.5,200,0,20);
	  fHistMassTwoGammasMultiGapRap1[iTrig][iTRU][iTOF][iMulti] = new TH2F("fHistMassTwoGammasGapRap1"+nameTrig[iTrig]+nameTRU[iTRU]+nameTOF[iTOF]+nameMulti[iMulti],"",500,0,0.5,200,0,20);
	  fHistMassTwoGammasMultiGapRap2[iTrig][iTRU][iTOF][iMulti] = new TH2F("fHistMassTwoGammasGapRap2"+nameTrig[iTrig]+nameTRU[iTRU]+nameTOF[iTOF]+nameMulti[iMulti],"",500,0,0.5,200,0,20);
	  fHistMassTwoGammasMultiV0AC[iTrig][iTRU][iTOF][iMulti]    = new TH2F("fHistMassTwoGammasV0AC"+nameTrig[iTrig]+nameTRU[iTRU]+nameTOF[iTOF]+nameMulti[iMulti],"",500,0,0.5,200,0,20);

	  fHistMixMassTwoGammasMultiUnitRap[iTrig][iTRU][iTOF][iMulti] = new TH2F("fHistMixMassTwoGammasUnitRap"+nameTrig[iTrig]+nameTRU[iTRU]+nameTOF[iTOF]+nameMulti[iMulti],"",500,0,0.5,200,0,20);
	  fHistMixMassTwoGammasMultiGapRap1[iTrig][iTRU][iTOF][iMulti] = new TH2F("fHistMixMassTwoGammasGapRap1"+nameTrig[iTrig]+nameTRU[iTRU]+nameTOF[iTOF]+nameMulti[iMulti],"",500,0,0.5,200,0,20);
	  fHistMixMassTwoGammasMultiGapRap2[iTrig][iTRU][iTOF][iMulti] = new TH2F("fHistMixMassTwoGammasGapRap2"+nameTrig[iTrig]+nameTRU[iTRU]+nameTOF[iTOF]+nameMulti[iMulti],"",500,0,0.5,200,0,20);
	  fHistMixMassTwoGammasMultiV0AC[iTrig][iTRU][iTOF][iMulti]    = new TH2F("fHistMixMassTwoGammasV0AC"+nameTrig[iTrig]+nameTRU[iTRU]+nameTOF[iTOF]+nameMulti[iMulti],"",500,0,0.5,200,0,20);
	  
	  fOutputContainer[1]->Add(fHistMassTwoGammasMultiUnitRap[iTrig][iTRU][iTOF][iMulti]);
	  fOutputContainer[1]->Add(fHistMassTwoGammasMultiGapRap1[iTrig][iTRU][iTOF][iMulti]);
	  fOutputContainer[1]->Add(fHistMassTwoGammasMultiGapRap2[iTrig][iTRU][iTOF][iMulti]);
	  fOutputContainer[1]->Add(fHistMassTwoGammasMultiV0AC[iTrig][iTRU][iTOF][iMulti]);

	  fOutputContainer[1]->Add(fHistMixMassTwoGammasMultiUnitRap[iTrig][iTRU][iTOF][iMulti]);
	  fOutputContainer[1]->Add(fHistMixMassTwoGammasMultiGapRap1[iTrig][iTRU][iTOF][iMulti]);
	  fOutputContainer[1]->Add(fHistMixMassTwoGammasMultiGapRap2[iTrig][iTRU][iTOF][iMulti]);
	  fOutputContainer[1]->Add(fHistMixMassTwoGammasMultiV0AC[iTrig][iTRU][iTOF][iMulti]);
	  
	}
	
      }
    }
  }

  fHistRadiusPi0      = new TH1F("fHistRadiusPi0","",500,0,500);
  fHistRadius2DPi0    = new TH2F("fHistRadius2DPi0","",1000,-500,500,1000,-500,500);
  fHistRadius2DRecPi0 = new TH2F("fHistRadius2DRecPi0","",1000,-500,500,1000,-500,500);
  if(fIsMC){
    fOutputContainer[1]->Add(fHistRadiusPi0);
    fOutputContainer[1]->Add(fHistRadius2DPi0);
    fOutputContainer[1]->Add(fHistRadius2DRecPi0);
  }
  
  fHistTotalEvents         = new TH1F("fHistTotalEvents","",30000,170000,200000);
  fHistAnalyzedEvents      = new TH1F("fHistAnalyzedEvents","",30000,170000,200000);
  fHistPUEvents            = new TH1F("fHistPUEvents","",30000,170000,200000);
  fHistPhysSelectionEvents = new TH1F("fHistPhysSelectionEvents","",30000,170000,200000);
  fHistAnalyzed0PH0Events  = new TH1F("fHistAnalyzed0PH0Events","",30000,170000,200000);
  
  fOutputContainer[2]->Add(fHistTotalEvents);
  fOutputContainer[2]->Add(fHistAnalyzedEvents);
  fOutputContainer[2]->Add(fHistPUEvents);
  fOutputContainer[2]->Add(fHistPhysSelectionEvents);
  fOutputContainer[2]->Add(fHistAnalyzed0PH0Events);
  
  fHistSPDTrackletsUnitRap = new TH2F("fHistSPDTrackletsUnitRap","",200,-10,10,300,0,300);
  fHistSPDTrackletsGapRap1 = new TH2F("fHistSPDTrackletsGapRap1","",200,-10,10,300,0,300);
  fHistSPDTrackletsGapRap2 = new TH2F("fHistSPDTrackletsGapRap2","",200,-10,10,300,0,300);
  fOutputContainer[2]->Add(fHistSPDTrackletsUnitRap);
  fOutputContainer[2]->Add(fHistSPDTrackletsGapRap1);
  fOutputContainer[2]->Add(fHistSPDTrackletsGapRap2);

  fHistCorrectedSPDTrackletsUnitRap = new TH2F("fHistCorrectedSPDTrackletsUnitRap","",200,-10,10,300,0,300);
  fHistCorrectedSPDTrackletsGapRap1 = new TH2F("fHistCorrectedSPDTrackletsGapRap1","",200,-10,10,300,0,300);
  fHistCorrectedSPDTrackletsGapRap2 = new TH2F("fHistCorrectedSPDTrackletsGapRap2","",200,-10,10,300,0,300);
  fOutputContainer[2]->Add(fHistCorrectedSPDTrackletsUnitRap);
  fOutputContainer[2]->Add(fHistCorrectedSPDTrackletsGapRap1);
  fOutputContainer[2]->Add(fHistCorrectedSPDTrackletsGapRap2);

  fHistV0AMulti  = new TH2F("fHistV0AMulti","",200,-10,10,500,0,1000);
  fHistV0CMulti  = new TH2F("fHistV0CMulti","",200,-10,10,500,0,1000);
  fHistV0ACMulti = new TH2F("fHistV0ACMulti","",200,-10,10,500,0,1000);
  fOutputContainer[2]->Add(fHistV0AMulti);
  fOutputContainer[2]->Add(fHistV0CMulti);
  fOutputContainer[2]->Add(fHistV0ACMulti);
  fHistCorrectedV0ACMulti = new TH2F("fHistCorrectedV0ACMulti","",200,-10,10,500,0,1000);
  fOutputContainer[2]->Add(fHistCorrectedV0ACMulti);
  
  fHistCorrelationMCTrackSPDTrackletsUnitRap = new TH2F("fHistCorrelationMCTrackSPDTrackletsUnitRap","",100,0,100,300,0,300);
  fHistCorrelationMCTrackSPDTrackletsGapRap1 = new TH2F("fHistCorrelationMCTrackSPDTrackletsGapRap1","",100,0,100,300,0,300);
  fHistCorrelationMCTrackSPDTrackletsGapRap2 = new TH2F("fHistCorrelationMCTrackSPDTrackletsGapRap2","",100,0,100,300,0,300);
  //fOutputContainer[2]->Add(fHistCorrelationMCTrackSPDTrackletsUnitRap);
  //fOutputContainer[2]->Add(fHistCorrelationMCTrackSPDTrackletsGapRap1);
  //fOutputContainer[2]->Add(fHistCorrelationMCTrackSPDTrackletsGapRap2);

  fHistCorrelationMCTrackCorrectedSPDTrackletsUnitRap = new TH2F("fHistCorrelationMCTrackCorrectedSPDTrackletsUnitRap","",100,0,100,300,0,300);
  fHistCorrelationMCTrackCorrectedSPDTrackletsGapRap1 = new TH2F("fHistCorrelationMCTrackCorrectedSPDTrackletsGapRap1","",100,0,100,300,0,300);
  fHistCorrelationMCTrackCorrectedSPDTrackletsGapRap2 = new TH2F("fHistCorrelationMCTrackCorrectedSPDTrackletsGapRap2","",100,0,100,300,0,300);
  fOutputContainer[2]->Add(fHistCorrelationMCTrackCorrectedSPDTrackletsUnitRap);
  fOutputContainer[2]->Add(fHistCorrelationMCTrackCorrectedSPDTrackletsGapRap1);
  fOutputContainer[2]->Add(fHistCorrelationMCTrackCorrectedSPDTrackletsGapRap2);
  
  fHistCorrelationMCTrackV0AC                        = new TH2F("fHistCorrelationMCTrackV0AC","",100,0,100,1000,0,2000);
  fHistCorrelationMCTrackCorrectedV0AC               = new TH2F("fHistCorrelationMCTrackCorrectedV0AC","",100,0,100,1000,0,2000);
  fHistCorrelationSPDTrackletsV0AC                   = new TH2F("fHistCorrelationSPDTrackletsV0AC","",100,0,100,1000,0,1000);
  fHistCorrelationCorrectedSPDTrackletsCorrectedV0AC = new TH2F("fHistCorrelationCorrectedSPDTrackletsCorrectedV0AC","",100,0,100,10000,0,2000);
  //fOutputContainer[2]->Add(fHistCorrelationMCTrackV0AC);
  //fOutputContainer[2]->Add(fHistCorrelationMCTrackCorrectedV0AC);
  //fOutputContainer[2]->Add(fHistCorrelationSPDTrackletsV0AC);
  fOutputContainer[2]->Add(fHistCorrelationCorrectedSPDTrackletsCorrectedV0AC);

  fHistZvertexPosition = new TH1F("fHistZvertexPosition","",60,-30,30);
  fOutputContainer[2]->Add(fHistZvertexPosition);
  
  fHistV0Timing        = new TH2F("fHistV0Timing","",700,-30,40,700,-20,50);
  fHistV0TimingPhysSel = new TH2F("fHistV0TimingPhysSel","",700,-30,40,700,-20,50);
  fOutputContainer[2]->Add(fHistV0Timing);
  fOutputContainer[2]->Add(fHistV0TimingPhysSel);

  fHistdEdxElesigma = new TH2F("fHistdEdxElesigma","",1000,0,10,200,-10,10);
  fOutputContainer[2]->Add(fHistdEdxElesigma);

  fHistdEventClassUnitRap = new TH1F("fHistdEventClassUnitRap","",10,0,10);
  fHistdEventClassGapRap1 = new TH1F("fHistdEventClassGapRap1","",10,0,10);
  fHistdEventClassGapRap2 = new TH1F("fHistdEventClassGapRap2","",10,0,10);
  fHistdEventClassV0AC    = new TH1F("fHistdEventClassV0AC","",10,0,10);
  
  fOutputContainer[2]->Add(fHistdEventClassUnitRap);
  fOutputContainer[2]->Add(fHistdEventClassGapRap1);
  fOutputContainer[2]->Add(fHistdEventClassGapRap2);
  fOutputContainer[2]->Add(fHistdEventClassV0AC);

  fHistdEventClassPHIUnitRap = new TH1F("fHistdEventClassPHIUnitRap","",10,0,10);
  fHistdEventClassPHIGapRap1 = new TH1F("fHistdEventClassPHIGapRap1","",10,0,10);
  fHistdEventClassPHIGapRap2 = new TH1F("fHistdEventClassPHIGapRap2","",10,0,10);
  fHistdEventClassPHIV0AC    = new TH1F("fHistdEventClassPHIV0AC","",10,0,10);
  
  fOutputContainer[2]->Add(fHistdEventClassPHIUnitRap);
  fOutputContainer[2]->Add(fHistdEventClassPHIGapRap1);
  fOutputContainer[2]->Add(fHistdEventClassPHIGapRap2);
  fOutputContainer[2]->Add(fHistdEventClassPHIV0AC);


  for(Int_t iMulti=0; iMulti<10; ++iMulti){
    for(Int_t iZvtx=0; iZvtx<4; ++iZvtx){
      
      fPHOSClusterList[iMulti][iZvtx] = new TList();
      fPHOSClusterList[iMulti][iZvtx]->SetOwner(true);
      
    }
  }
  
  PostData(1, fOutputContainer[0]);
  PostData(2, fOutputContainer[1]);
  PostData(3, fOutputContainer[2]);
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::UserExec(Option_t *)
{
  
  this->EventFlagInit();
  this->ConnectInputData();
  
  if(!fAOD[0]) 
    return;
  
  if(fRunNumber != fAOD[0]->GetRunNumber()){
    this->SetPeriod();
    this->SetGeometry();
    this->SetCalibration();
  }

  this->SetVertex();
  
  if(!this->IsGoodEvent(fAOD[0]))
    return ;
  this->SetTriggerInfo();
  this->SetMultiplicity();
  this->ProcessMC();
  this->NeutralPion();
  this->DirectPhotonAnalysis(fAOD[0]);
  
  fHistdEventClassUnitRap->Fill(fMultiBinSPDUnitRap);
  fHistdEventClassGapRap1->Fill(fMultiBinSPDGapRap1);
  fHistdEventClassGapRap2->Fill(fMultiBinSPDGapRap2);
  fHistdEventClassV0AC->Fill(fMultiBinV0AC);
  
  if(f0PH0Event){
    fHistdEventClassPHIUnitRap->Fill(fMultiBinSPDUnitRap);
    fHistdEventClassPHIGapRap1->Fill(fMultiBinSPDGapRap1);
    fHistdEventClassPHIGapRap2->Fill(fMultiBinSPDGapRap2);
    fHistdEventClassPHIV0AC->Fill(fMultiBinV0AC);
    fHistAnalyzed0PH0Events->Fill(fRunNumber);
  }
  
  PostData(1, fOutputContainer[0]);
  PostData(2, fOutputContainer[1]);
  PostData(3, fOutputContainer[2]);
  
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::SetPeriod()
{
  fRunNumber = fAOD[0]->GetRunNumber();
  if(176701<=fRunNumber && fRunNumber<=177295) fPeriod = "LHC12a";
  if(177381<=fRunNumber && fRunNumber<=178220) fPeriod = "LHC12b";
  if(179444<=fRunNumber && fRunNumber<=182750) fPeriod = "LHC12c";
  if(183913<=fRunNumber && fRunNumber<=186320) fPeriod = "LHC12d";
  if(188796<=fRunNumber && fRunNumber<=192824) fPeriod = "LHC12h";
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::SetGeometry()
{
  
  fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP"); 
  AliOADBContainer geomContainer("phosGeo");
  
  if(!fIsMC){
    geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
    TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");
    for (Int_t mod=0; mod<5; mod++) {
      if (!matrixes->At(mod)) continue;
      else fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod);
    }
  }
  
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::SetCalibration()
{
  
  AliOADBContainer calibContainer("phosRecalibration");
  TObjArray *recalib = NULL;

  if(fIsMC){
    calibContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSMCCalibrations.root","phosRecalibration");
    recalib = (TObjArray*)calibContainer.GetObject(fRunNumber,"PHOSRecalibration","LHC14e2b");
    AliCDBEntry* entryEmc = fCDBstorage->Get("PHOS/Calib/EmcGainPedestals",fRunNumber);
    fCalibDataEmc = (AliPHOSEmcCalibData*)entryEmc->GetObject();
  }
  else{
    calibContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSCalibrations.root","phosRecalibration");
    recalib = (TObjArray*)calibContainer.GetObject(fRunNumber,"PHOSRecalibration");
  }
  
  fPHOSCalibData = (AliPHOSCalibData*)recalib->At(0) ;
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if(inputHandler->GetPIDResponse())
    fPIDResponse=inputHandler->GetPIDResponse();
  
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::SetVertex()
{

  const AliVVertex *primaryVertex = fAOD[0]->GetPrimaryVertexSPD();
  Double_t fVertex[3];

  if( primaryVertex ) {
    fVertex[0] = primaryVertex->GetX();
    fVertex[1] = primaryVertex->GetY();
    fVertex[2] = primaryVertex->GetZ();
    fIsVtxNoCont = false;
  }
  else {
    AliError("Event has 0 Primary Vertex, defaulting to vtx0");
    fVertex[0] = -999;
    fVertex[1] = -999;
    fVertex[2] = -999;
    fIsVtxNoCont = true;
  }
  if(primaryVertex->GetNContributors()<1){
    fVertex[0] = -999;
    fVertex[1] = -999;
    fVertex[2] = -999;
    fIsVtxNoCont = true;
  }

  fVertexVector = TVector3(fVertex);



}
//________________________________________________________________________
Bool_t AliAnalysisTaskPHOSTrigPi0::SetPhysicsSelection(AliAODEvent* AOD)
{
  
  AliAODVZERO * aodV0 = AOD->GetVZEROData();
  Float_t tA = (Double_t) aodV0->GetV0ATime()*1E-9;
  Float_t tC = (Double_t) aodV0->GetV0CTime()*1E-9;
  
  Double_t minDiff = 5.0E-9;
  Double_t maxDiff = 11.5E-9;
  Double_t minSum  = 11.5E-9;
  Double_t maxSum  = 17.5E-9;
  Bool_t isSelected = false;
  
  fHistV0Timing->Fill((tA-tC)*1E+9,(tA+tC)*1E+9);

  if(tA-tC > minDiff && tA-tC < maxDiff && tA+tC > minSum  && tA+tC < maxSum){
    fHistV0TimingPhysSel->Fill((tA-tC)*1E+9,(tA+tC)*1E+9);
    return true;
  }
  else
    return false;
  
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPHOSTrigPi0::IsGoodEvent(AliAODEvent* AOD){

  Bool_t isPhysicsSelect = false;

  if(this->SetPhysicsSelection(fAOD[0])){
    isPhysicsSelect = true;
  }
  if(AOD->IsPileupFromSPD()){
    fIsPileUp = true;
  }
  if(fabs(fVertexVector.Z())>10.){
    fIsVtxOut10cm= true;
  }
  
  fHistTotalEvents->Fill(fRunNumber);
  
  if(fIsPileUp){
    fHistPUEvents->Fill(fRunNumber);
  }
  if(isPhysicsSelect){
    fHistPhysSelectionEvents->Fill(fRunNumber);
  }

  if(!isPhysicsSelect){
    return false;
  }
  if(fIsPileUp){
    return false;
  }
  if(fIsVtxOut10cm){
    return false;
  }
  if(fIsVtxNoCont){
    return false;
  }  

  fHistAnalyzedEvents->Fill(fRunNumber);
  fHistZvertexPosition->Fill(fVertexVector.Z());

  if(-10<fVertexVector.Z() && fVertexVector.Z()<=-5){
    fZvtxBin = 0;
  }
  else if(-5<fVertexVector.Z() && fVertexVector.Z()<=0){
    fZvtxBin = 1;
  }
  else if(0<fVertexVector.Z() && fVertexVector.Z()<=5){
    fZvtxBin = 2;
  }
  else if(5<fVertexVector.Z() && fVertexVector.Z()<=10){
    fZvtxBin = 3;
  }
  return true;
  
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::ConnectInputData(){
  
  fAOD[0] = dynamic_cast<AliAODEvent*> (InputEvent());
  
  if(fIsMC)
    this->GetMCStack();
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::DirectPhotonAnalysis(AliAODEvent* fEvent)
{
  
  AliAODCaloCells* cells = dynamic_cast<AliAODCaloCells*> (fEvent->GetPHOSCells());

  for (Int_t i=0; i<fEvent->GetNumberOfCaloClusters(); i++) {

    AliAODCaloCluster *clu = dynamic_cast<AliAODCaloCluster*>(fEvent->GetCaloCluster(i));
    fTrigData->Reset();
    
    if( !clu->IsPHOS()) continue; // reject EMCal                                                                                                                                                                                            

    Float_t  position[3];
    clu->GetPosition(position);
    TVector3 global(position) ;
    Int_t relId[4] = {0,0,0,0};
    fPHOSGeo->GlobalPos2RelId(global,relId) ;
    Int_t    mod = relId[0] ;
    Int_t    col = relId[2] ;
    Int_t    row = relId[3] ;
    Int_t ncells = clu->GetNCells();
    Int_t    tru = GetTRUNum(col,row);

    Double_t vtx5[3] = {fVertexVector.X(),fVertexVector.Y(),fVertexVector.Z()};
    TLorentzVector ClusterVector;
    AliPHOSAodCluster cluPHOS( *(AliAODCaloCluster*) (clu) );
    cluPHOS.SetE(AddEnergyCalib(mod,cluPHOS.E()));
    cluPHOS.SetPosition(position);
    cluPHOS.GetMomentum(ClusterVector,vtx5);
    
    if(!(ClusterVector.E()>0.))                continue;
    if(!fIsGoodMod[mod-1])                     continue;
    //if(fBadMap[mod-1]->GetBinContent(col,row)) continue;    

    Int_t  maxAbsId = 0;
    Int_t  maxrelId[4];
    MaxEnergyCellPos(cells,clu,maxAbsId);
    fPHOSGeo->AbsToRelNumbering(maxAbsId,maxrelId);
    
    Double_t M02, M20, t1, t2, t3, t4, t5, t6;
    CalculatedClusterParameters(clu, cells, AddEnergyCalib(mod,ClusterVector.E())/clu->E(), M02, M20,t1, t2,t3, t4, t5, t6);

    M02 = cluPHOS.GetM02();
    M20 = cluPHOS.GetM20();

    Bool_t isDispOK = false;
    
    if(cluPHOS.Chi2() < fDispCut*fDispCut)
      isDispOK = true;
    else
      isDispOK = false;
    
    if(ClusterVector.E()>fMinEene){
      fHistClustTOF[0]->Fill(ClusterVector.E(),clu->GetTOF());
      fHistClustTOFv1[0]->Fill(ClusterVector.E(),t1);
      fHistClustTOFv2[0]->Fill(ClusterVector.E(),t2);
      fHistClustTOFv3[0]->Fill(ClusterVector.E(),t3);
      fHistClustTOFv4[0]->Fill(ClusterVector.E(),t4);
      fHistClustTOFv5[0]->Fill(ClusterVector.E(),t5);
      fHistClustTOFv6[0]->Fill(ClusterVector.E(),t6);
      fHistClustTOF[mod]->Fill(ClusterVector.E(),t6);
    }
    
    fHistClustOccupMap[mod-1]->Fill(col,row,1);
    
    Double_t track_dx = clu->GetTrackDx();
    Double_t track_dz = clu->GetTrackDz();

    Double_t track_R  = 0;
    if(track_dx>0 || track_dz>0){
      track_R = sqrt(pow(track_dx,2) + pow(track_dz,2));
    }
    
    Bool_t isTrig       = false;
    
    if(IsTrigger(fTrigData,cells,maxAbsId,ClusterVector.E(),t6))
      isTrig = true;
    
    Int_t mulDigit=clu->GetNCells() ;
    
    for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
      
      
      Int_t absId=clu->GetCellAbsId(iDigit) ;
      
      Bool_t isHG=false ;
      if(cells->GetCellHighGain(absId))
	isHG=true;
      
      Double_t tof = cells->GetCellTime(absId);
      Double_t amp = cells->GetCellAmplitude(absId);

      if(isHG){
	if(amp>1.0)
	  fHistCellTimeHG_HE->Fill(absId,tof);
	else
	  fHistCellTimeHG->Fill(absId,tof);
      }
      else{
	if(amp>1.0)
	  fHistCellTimeLG_HE->Fill(absId,tof);
	else
	  fHistCellTimeLG->Fill(absId,tof);
      }
     
      if(isTrig){
	if(isHG){
	  if(amp>1.0)
	    fHist0PH0CellTimeHG_HE->Fill(absId,tof);
	  else
	    fHist0PH0CellTimeHG->Fill(absId,tof);
	}
	else{
	  if(amp>1.0)
	    fHist0PH0CellTimeLG_HE->Fill(absId,tof);
	  else
	    fHist0PH0CellTimeLG->Fill(absId,tof);
	}
      }//isTrig
 
    }
    
    if(isTrig)
      fHist0PH0ClustOccupMap[mod-1]->Fill(col,row,1);
    
    Bool_t isTrackMatch = IsTrackCluster(clu->GetNTracksMatched(),track_dx,track_dz,ClusterVector.E());
    
    if(isTrig){
      f0PH0Event = true;
      if(cluPHOS.E()>6.0)
	f0PH0Event_HighEnergy = true;
    }
    
    if(isTrig && fabs(t6)<fTOFcut){
      fTOFcut0PH0Event= true;
      if(cluPHOS.E()>6.0)
	fTOFcut0PH0Event_HighEnergy = true;
    }

    fHistClustEneAll->Fill(ClusterVector.E());
    if(fabs(t6)<25.E-9){
      fHistClustEneAllTOFcut1->Fill(ClusterVector.E());
    }
    if(fabs(t6)<50.E-9){
      fHistClustEneAllTOFcut2->Fill(ClusterVector.E());
    }
    if(fabs(t6)<75.E-9){
      fHistClustEneAllTOFcut3->Fill(ClusterVector.E());
    }    
    if(fabs(t6)<100.E-9){
      fHistClustEneAllTOFcut4->Fill(ClusterVector.E());
    }   

    if(ncells>=fMinNCell){
      fHistClustEneNcellCut->Fill(ClusterVector.E());
      if(fabs(t6)<25.E-9){
	fHistClustEneNcellCutTOFcut1->Fill(ClusterVector.E());
      }
      if(fabs(t6)<50.E-9){
	fHistClustEneNcellCutTOFcut2->Fill(ClusterVector.E());
      }
      if(fabs(t6)<75.E-9){
	fHistClustEneNcellCutTOFcut3->Fill(ClusterVector.E());
      }    
      if(fabs(t6)<100.E-9){
	fHistClustEneNcellCutTOFcut4->Fill(ClusterVector.E());
      }   
      
    }
    if(fIsGoodTRU[mod-1][tru-1]){
      fHistGoodTRUClustEneAll->Fill(ClusterVector.E());
      if(fabs(t6)<25.E-9){
	fHistGoodTRUClustEneAllTOFcut1->Fill(ClusterVector.E());
      }
      if(fabs(t6)<50.E-9){
	fHistGoodTRUClustEneAllTOFcut2->Fill(ClusterVector.E());
      }
      if(fabs(t6)<75.E-9){
	fHistGoodTRUClustEneAllTOFcut3->Fill(ClusterVector.E());
      }    
      if(fabs(t6)<100.E-9){
	fHistGoodTRUClustEneAllTOFcut4->Fill(ClusterVector.E());
      }   
      
      if(ncells>=fMinNCell){
	fHistGoodTRUClustEneNcellCut->Fill(ClusterVector.E());
	if(fabs(t6)<25.E-9){
	  fHistGoodTRUClustEneNcellCutTOFcut1->Fill(ClusterVector.E());
	}
	if(fabs(t6)<50.E-9){
	  fHistGoodTRUClustEneNcellCutTOFcut2->Fill(ClusterVector.E());
	}
	if(fabs(t6)<75.E-9){
	  fHistGoodTRUClustEneNcellCutTOFcut3->Fill(ClusterVector.E());
	}    
	if(fabs(t6)<100.E-9){
	  fHistGoodTRUClustEneNcellCutTOFcut4->Fill(ClusterVector.E());
	}   
	
      }
    }

    if(isTrig){
      fHist0PH0ClustEneAll->Fill(ClusterVector.E());
      if(fabs(t6)<25.E-9){
	fHist0PH0ClustEneAllTOFcut1->Fill(ClusterVector.E());
      }
      if(fabs(t6)<50.E-9){
	fHist0PH0ClustEneAllTOFcut2->Fill(ClusterVector.E());
      }
      if(fabs(t6)<75.E-9){
	fHist0PH0ClustEneAllTOFcut3->Fill(ClusterVector.E());
      }    
      if(fabs(t6)<100.E-9){
	fHist0PH0ClustEneAllTOFcut4->Fill(ClusterVector.E());
      }   
      
      if(ncells>=fMinNCell){
	fHist0PH0ClustEneNcellCut->Fill(ClusterVector.E());
	if(fabs(t6)<25.E-9){
	  fHist0PH0ClustEneNcellCutTOFcut1->Fill(ClusterVector.E());
	}
	if(fabs(t6)<50.E-9){
	  fHist0PH0ClustEneNcellCutTOFcut2->Fill(ClusterVector.E());
	}
	if(fabs(t6)<75.E-9){
	  fHist0PH0ClustEneNcellCutTOFcut3->Fill(ClusterVector.E());
	}    
	if(fabs(t6)<100.E-9){
	  fHist0PH0ClustEneNcellCutTOFcut4->Fill(ClusterVector.E());
	}   
	
      }
      if(fIsGoodTRU[mod-1][tru-1]){
	fHistGoodTRU0PH0ClustEneAll->Fill(ClusterVector.E());
	if(fabs(t6)<25.E-9){
	  fHistGoodTRU0PH0ClustEneAllTOFcut1->Fill(ClusterVector.E());
	}
	if(fabs(t6)<50.E-9){
	  fHistGoodTRU0PH0ClustEneAllTOFcut2->Fill(ClusterVector.E());
	}
	if(fabs(t6)<75.E-9){
	  fHistGoodTRU0PH0ClustEneAllTOFcut3->Fill(ClusterVector.E());
	}    
	if(fabs(t6)<100.E-9){
	  fHistGoodTRU0PH0ClustEneAllTOFcut4->Fill(ClusterVector.E());
	}   
	
	if(ncells>=fMinNCell){
	  fHistGoodTRU0PH0ClustEneNcellCut->Fill(ClusterVector.E());
	  if(fabs(t6)<25.E-9){
	    fHistGoodTRU0PH0ClustEneNcellCutTOFcut1->Fill(ClusterVector.E());
	  }
	  if(fabs(t6)<50.E-9){
	    fHistGoodTRU0PH0ClustEneNcellCutTOFcut2->Fill(ClusterVector.E());
	  }
	  if(fabs(t6)<75.E-9){
	    fHistGoodTRU0PH0ClustEneNcellCutTOFcut3->Fill(ClusterVector.E());
	  }    
	  if(fabs(t6)<100.E-9){
	    fHistGoodTRU0PH0ClustEneNcellCutTOFcut4->Fill(ClusterVector.E());
	  }   
	  
	}
      }
    }

    
    if(fabs(t6) > fTOFcut)
      continue;
    if(!fIsGoodTRU[mod-1][tru-1])
      isTrig = false;

    Bool_t isCPVOK = false;
    
    TString trackPID = "";
    Double_t Ep = -100;
    Double_t p  = -100; 
    if(fPIDResponse && clu->GetNTracksMatched()>0){
      AliAODTrack   *track     = dynamic_cast<AliAODTrack*>(clu->GetTrackMatched(0));
      
      Double_t nSigmaTPC_Electron = fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType)AliPID::kElectron);
      Double_t nSigmaTPC_Proton   = fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType)AliPID::kProton);
      fHistdEdxElesigma->Fill(track->P(),nSigmaTPC_Electron);

      Ep = ClusterVector.E()/track->P();
      p  = track->P();
      trackPID = GetTrackPID(track);
      isCPVOK = true;
    }

    Int_t iCand = 0;
    Int_t iPID  = 0;
    
    if(!isTrackMatch){
      iCand = 0;
      iPID  = 0;
      fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
      fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
      
      if(ncells>=fMinNCell)
	fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
      
      if(ClusterVector.E()>0.)
	fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
      if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
	fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
      }
      if(isCPVOK){
	iPID = 1;
	fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
	  fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	if(ClusterVector.E()>0.)
	  fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	  fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
	  fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	  fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	  fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	  fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	  fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	  fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	}
      }
      if(isDispOK){
	iPID = 2;
	fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
	  fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	if(ClusterVector.E()>0.)
	  fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	  //fHistClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	  fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
	  fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	  fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	  fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	  fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	  fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	  fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	}
      }
      if(isCPVOK && isDispOK){
	iPID = 3;
	fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
	  fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	if(ClusterVector.E()>0.)
	  fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	  //fHistClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	  fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
	  fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	  fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	  fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	  fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	  fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	  fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	}
      }
    }
    if(isTrackMatch){
      iCand = 1;
      iPID  = 0;
      fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
      fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
      if(ncells>=fMinNCell)
	fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
      if(ClusterVector.E()>0.)
        fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
      if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
        //fHistClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
        fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
        fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
        fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
        fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
        fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
        fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
        fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
      }
      if(isCPVOK){
	iPID = 1;
	fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
	  fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	if(ClusterVector.E()>0.)
	  fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	  //fHistClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	  fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
	  fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	  fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	  fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	  fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	  fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	  fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	}
      }
      if(isDispOK){
	iPID = 2;
	fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
	  fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	if(ClusterVector.E()>0.)
	  fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	  //fHistClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	  fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
	  fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	  fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	  fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	  fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	  fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	  fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	}
      }
      if(isCPVOK && isDispOK){
	iPID = 3;
	fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
	  fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	if(ClusterVector.E()>0.)
	  fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	  //fHistClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	  fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
	  fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	  fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	  fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	  fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	  fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	  fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	}
      }
      if(trackPID.Contains("Electron")){
	iCand = 2;
	iPID  = 0;
	fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
	  fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	if(ClusterVector.E()>0.)
	  fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	  //fHistClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	  fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
	  fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	  fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	  fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	  fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	  fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	  fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	}
	if(isCPVOK){
	  iPID = 1;
	  fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	  fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	    fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
	    fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	    fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isDispOK){
	  iPID = 2;
	  fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	  fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	    fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
	    fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	    fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isCPVOK && isDispOK){
	  iPID = 3;
	  fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	  fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	    fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
	    fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	    fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
      }
      if(trackPID.Contains("ChargedHadron")){
	iCand = 3;
	iPID  = 0;
	fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
	  fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	if(ClusterVector.E()>0.)
	  fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	  //fHistClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	  fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
	  fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	  fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	  fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	  fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	  fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	  fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	}
	if(isCPVOK){
	  iPID = 1;
	  fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	  fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	    fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
	    fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	    fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isDispOK){
	  iPID = 2;
	  fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	  fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	    fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
	    fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	    fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isCPVOK && isDispOK){
	  iPID = 3;
	  fHistClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	  fHistClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	    fHistClustEp[iCand][iPID][0]->Fill(p,Ep);
	    fHistClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	    fHistClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    
	  }
	}
      }
    }

    if(isTrig){
      if(!isTrackMatch){
	iCand = 0;
	iPID  = 0;
	//fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	//fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
	  fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	if(ClusterVector.E()>0.)
	  //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	  //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);
	  fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	  //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	  //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	  //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	  //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	  //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	  //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	}
	if(isCPVOK){
	  iPID = 1;
	  //fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	  //fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	 //if(ClusterVector.E()>0.)
	    //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	    fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	    //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	    //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isDispOK){
	  iPID = 2;
	  //fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	  //fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	 //if(ClusterVector.E()>0.)
	    //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	    fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	    //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	    //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isCPVOK && isDispOK){
	  iPID = 3;
	  //fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	  //fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	 //if(ClusterVector.E()>0.)
	    //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	    fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	    //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	    //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
      }
      if(isTrackMatch){
	iCand = 1;
	iPID  = 0;
	//fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	//fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
	  fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	if(ClusterVector.E()>0.)
	  //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	  //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	  fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	  //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	  //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	  //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	  //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	  //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	  //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	}
	if(isCPVOK){
	  iPID = 1;
	  //fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	  //fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	 //if(ClusterVector.E()>0.)
	    //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	    fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	    //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	    //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isDispOK){
	  iPID = 2;
	  //fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	  //fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	 //if(ClusterVector.E()>0.)
	    //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	    fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	    //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	    //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isCPVOK && isDispOK){
	  iPID = 3;
	  //fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	  //fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	 //if(ClusterVector.E()>0.)
	  //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	    fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	    //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	    //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(trackPID.Contains("Electron")){
	  iCand = 2;
	  iPID  = 0;
	  //fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	  //fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	 //if(ClusterVector.E()>0.)
	    //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	    fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	    //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	    //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	  if(isCPVOK){
	    iPID = 1;
	    //fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	   //if(ClusterVector.E()>0.)
	      //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	      fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	      //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	      //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isDispOK){
	    iPID = 2;
	    //fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	   //if(ClusterVector.E()>0.)
	      //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	      fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	      //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	      //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isCPVOK && isDispOK){
	    iPID = 3;
	    //fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	   //if(ClusterVector.E()>0.)
	      //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	      fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	      //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	      //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	}
	if(trackPID.Contains("ChargedHadron")){
	  iCand = 3;
	  iPID  = 0;
	  //fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	  //fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	 //if(ClusterVector.E()>0.)
	    //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	    fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	    //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	    //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	  if(isCPVOK){
	    iPID = 1;
	    //fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	   //if(ClusterVector.E()>0.)
	    //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	      fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	      //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	      //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isDispOK){
	    iPID = 2;
	    //fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	   //if(ClusterVector.E()>0.)
	      //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	      fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	      //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	      //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isCPVOK && isDispOK){
	    iPID = 3;
	    //fHist0PH0ClustEneMIP[iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHist0PH0ClustEneMIP[iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHist0PH0ClustEne[iCand][iPID][0]->Fill(ClusterVector.E());
	   //if(ClusterVector.E()>0.)
	      //fHist0PH0ClustNcells[iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHist0PH0ClustEp[iCand][iPID][0]->Fill(ClusterVector.E(),Ep);    
	      fHist0PH0ClustEp[iCand][iPID][0]->Fill(p,Ep);
	      //fHist0PH0ClustM02[iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHist0PH0ClustM20[iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHist0PH0ClustM02M20[iCand][iPID][0]->Fill(M02,M20);
	      //fHist0PH0ClustDx[iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHist0PH0ClustDz[iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHist0PH0ClustDr[iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	      
	    }
	  }
	}
      }

    }

    ///////////////////////////////////////////////////////////////////////
    //For Event mixing pool
    ///////////////////////////////////////////////////////////////////////
    
    AliAODCaloCluster *cluPool = (AliAODCaloCluster*)clu->Clone("cluPool");
    cluPool->SetTOF(t6);
    cluPool->SetM02(M02);
    cluPool->SetM02(M20);
    if(isTrig)
      cluPool->SetID(1);
    else
      cluPool->SetID(0);
    fPHOSClusterList[fMultiBinSPDUnitRap][fZvtxBin]->AddFirst(cluPool);
    if(fPHOSClusterList[fMultiBinSPDUnitRap][fZvtxBin]->GetEntries() > 10){
      fPHOSClusterList[fMultiBinSPDUnitRap][fZvtxBin]->RemoveLast() ;
    }

    ///////////////////////////////////////////////////////////////////////
    //
    ///////////////////////////////////////////////////////////////////////
    
    if(!fIsMC)
      continue;

    Int_t  PDGcode     = 0;
    Int_t  clustCharge = 0;
    Int_t  momLabel    = 0;
    Bool_t sure = false;
    
    Bool_t isPhoton                   = false; // largest contribution to cluster is photon
    Bool_t isElectron                 = false; // largest contribution to cluster is electron
    Bool_t isProton                   = false; // largest contribution to cluster is proton
    Bool_t isAntiProton               = false; // largest contribution to cluster is antiproton
    Bool_t isNeutron                  = false; // largest contribution to cluster is neutron
    Bool_t isAntiNeutron              = false; // largest contribution to cluster is antineutron
    Bool_t isConversion               = false; // largest contribution to cluster is converted electron
    Bool_t isChargedParticle          = false; // largest contribution to cluster is charged particle
    Bool_t isNeutralParticle          = false; // largest contribution to cluster is neutral particle
    Bool_t isChargedPion              = false; // largest contribution to cluster is charged pion
    Bool_t isRadius1cmCut             = false;

    if(fIsMC){
      PDGcode = GetPDGCode(FindPrimary(clu,sure),clustCharge,momLabel,isRadius1cmCut);
      this->SetClusterFlag(PDGcode,clustCharge,isPhoton,isElectron,isProton,isAntiProton,isNeutron,isAntiNeutron,isConversion,isChargedParticle,isNeutralParticle,isChargedPion);
    }
    
    Int_t iTrue=0;
    
    if(isPhoton){
      iTrue = 0;
    }
    else if(isElectron){
      iTrue = 2;
    }
    else if(isProton){
      iTrue = 3;
    }
    else if(isAntiProton){
      iTrue = 4;
    }
    else if(isNeutron){
      iTrue = 5;
    }
    else if(isAntiNeutron){
      iTrue = 6;
    }
    else if(isChargedPion){
      iTrue = 9;
    }
    
    if(iTrue==0 || iTrue==2 || iTrue==3 || iTrue==4 || iTrue==5 || iTrue==6 || iTrue==9){
      if(!isTrackMatch){
	iCand = 0;
	iPID  = 0;
	fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
	  fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	if(ClusterVector.E()>0.)
	  fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	  //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);
	  fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	  fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	  fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	  fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	  fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	  fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	  fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	}
	if(isCPVOK){
	  iPID = 1;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isDispOK){
	  iPID = 2;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isCPVOK && isDispOK){
	  iPID = 3;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
      }
      else{
	iCand = 1;
	iPID  = 0;
	fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
          fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	if(ClusterVector.E()>0.)
	  fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
        if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
          //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
          fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
          fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
          fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
          fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
          fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
          fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
          fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
        }
	if(isCPVOK){
	  iPID = 1;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isDispOK){
	  iPID = 2;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isCPVOK && isDispOK){
	  iPID = 3;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(trackPID.Contains("Electron")){
	  iCand = 2;
	  iPID  = 0;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	  if(isCPVOK){
	    iPID = 1;
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    if(ClusterVector.E()>0.)
	      fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isDispOK){
	    iPID = 2;
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    if(ClusterVector.E()>0.)
	      fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isCPVOK && isDispOK){
	    iPID = 3;
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    if(ClusterVector.E()>0.)
	      fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	}
	if(trackPID.Contains("ChargedHadron")){
	  iCand = 3;
	  iPID  = 0;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	  if(isCPVOK){
	    iPID = 1;
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    if(ClusterVector.E()>0.)
	      fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isDispOK){
	    iPID = 2;
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    if(ClusterVector.E()>0.)
	      fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isCPVOK && isDispOK){
	    iPID = 3;
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    if(ClusterVector.E()>0.)
	      fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }

	  }
	}
      }
    }
    if(isChargedParticle){
      iTrue = 7;
    }
    else{
      iTrue = 8;
    }
    if(iTrue==7 || iTrue==8){
      if(!isTrackMatch){
	iCand = 0;
	iPID  = 0;
	fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
          fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	if(ClusterVector.E()>0.)
	  fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
        if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
          //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
          fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
          fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
          fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
          fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
          fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
          fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
          fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
        }
	if(isCPVOK){
	  iPID = 1;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isDispOK){
	  iPID = 2;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isCPVOK && isDispOK){
	  iPID = 3;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
      }
      else{
	iCand = 1;
	iPID  = 0;
	fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
          fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	if(ClusterVector.E()>0.)
	  fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
        if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
          //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
          fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
          fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
          fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
          fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
          fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
          fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
          fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
        }
	if(isCPVOK){
	  iPID = 1;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isDispOK){
	  iPID = 2;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isCPVOK && isDispOK){
	  iPID = 3;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(trackPID.Contains("Electron")){
	  iCand = 2;
	  iPID  = 0;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	  if(isCPVOK){
	    iPID = 1;
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    if(ClusterVector.E()>0.)
	      fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isDispOK){
	    iPID = 2;
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    if(ClusterVector.E()>0.)
	      fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isCPVOK && isDispOK){
	    iPID = 3;
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    if(ClusterVector.E()>0.)
	      fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	}
	if(trackPID.Contains("ChargedHadron")){
	  iCand = 3;
	  iPID  = 0;
	  fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  if(ClusterVector.E()>0.)
	    fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	  if(isCPVOK){
	    iPID = 1;
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    if(ClusterVector.E()>0.)
	      fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isDispOK){
	    iPID = 2;
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    if(ClusterVector.E()>0.)
	      fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isCPVOK && isDispOK){
	    iPID = 3;
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    fHistTrueClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrueClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    if(ClusterVector.E()>0.)
	      fHistTrueClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);
	      fHistTrueClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      fHistTrueClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      fHistTrueClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      fHistTrueClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      fHistTrueClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      fHistTrueClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      fHistTrueClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	}
      }
    }

    if(isTrig){
      if(iTrue==0 || iTrue==2 || iTrue==3 || iTrue==4 || iTrue==5 || iTrue==6 || iTrue==9){
	if(!isTrackMatch){
	  iCand = 0;
	  iPID  = 0;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	  //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	if(isCPVOK){
	  iPID = 1;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	  //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isDispOK){
	  iPID = 2;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	  //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isCPVOK && isDispOK){
	  iPID = 3;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	  //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
      }
      else{
	iCand = 1;
	iPID  = 0;
	//fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	//fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
          fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	//if(ClusterVector.E()>0.)
	//fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
        if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
          //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
          fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
          //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
          //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
          //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
          //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
          //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
          //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
        }
	if(isCPVOK){
	  iPID = 1;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	  //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isDispOK){
	  iPID = 2;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	  //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isCPVOK && isDispOK){
	  iPID = 3;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	  //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(trackPID.Contains("Electron")){
	  iCand = 2;
	  iPID  = 0;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	  //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	  if(isCPVOK){
	    iPID = 1;
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //if(ClusterVector.E()>0.)
	    //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isDispOK){
	    iPID = 2;
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //if(ClusterVector.E()>0.)
	    //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isCPVOK && isDispOK){
	    iPID = 3;
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //if(ClusterVector.E()>0.)
	    //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	}
	if(trackPID.Contains("ChargedHadron")){
	  iCand = 3;
	  iPID  = 0;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	  //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	  if(isCPVOK){
	    iPID = 1;
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //if(ClusterVector.E()>0.)
	    //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isDispOK){
	    iPID = 2;
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //if(ClusterVector.E()>0.)
	    //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isCPVOK && isDispOK){
	    iPID = 3;
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //if(ClusterVector.E()>0.)
	    //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }

	  }
	}
      }
    }
    if(isChargedParticle){
      iTrue = 7;
    }
    else{
      iTrue = 8;
    }
    if(iTrue==7 || iTrue==8){
      if(!isTrackMatch){
	iCand = 0;
	iPID  = 0;
	//fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	//fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
          fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	//if(ClusterVector.E()>0.)
	//fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
        if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
          //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
          fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
          //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
          //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
          //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
          //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
          //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
          //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
        }
	if(isCPVOK){
	  iPID = 1;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	  //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isDispOK){
	  iPID = 2;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	    //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isCPVOK && isDispOK){
	  iPID = 3;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	    //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
      }
      else{
	iCand = 1;
	iPID  = 0;
	//fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	//fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	if(ncells>=fMinNCell)
          fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	//if(ClusterVector.E()>0.)
	  //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
        if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
          //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
          fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
          //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
          //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
          //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
          //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
          //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
          //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
        }
	if(isCPVOK){
	  iPID = 1;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	    //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isDispOK){
	  iPID = 2;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	    //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(isCPVOK && isDispOK){
	  iPID = 3;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	    //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	}
	if(trackPID.Contains("Electron")){
	  iCand = 2;
	  iPID  = 0;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	    //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	  if(isCPVOK){
	    iPID = 1;
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //if(ClusterVector.E()>0.)
	      //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isDispOK){
	    iPID = 2;
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //if(ClusterVector.E()>0.)
	      //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isCPVOK && isDispOK){
	    iPID = 3;
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //if(ClusterVector.E()>0.)
	      //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	}
	if(trackPID.Contains("ChargedHadron")){
	  iCand = 3;
	  iPID  = 0;
	  //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
          //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	  if(ncells>=fMinNCell)
	    fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	  //if(ClusterVector.E()>0.)
	    //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	  if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	    //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	    fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	    //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	    //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	    //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	    //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	    //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	    //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	  }
	  if(isCPVOK){
	    iPID = 1;
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //if(ClusterVector.E()>0.)
	      //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isDispOK){
	    iPID = 2;
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //if(ClusterVector.E()>0.)
	      //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);      
	      fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	  if(isCPVOK && isDispOK){
	    iPID = 3;
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //fHistTrue0PH0ClustEneMIP[iTrue][iCand][iPID][mod]->Fill(ClusterVector.E());
	    if(ncells>=fMinNCell)
	      fHistTrue0PH0ClustEne[iTrue][iCand][iPID][0]->Fill(ClusterVector.E());
	    //if(ClusterVector.E()>0.)
	      //fHistTrue0PH0ClustNcells[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),ncells);
	    if(ncells>=fMinNCell && ClusterVector.E()>fMinEene){
	      //fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),Ep);
	      fHistTrue0PH0ClustEp[iTrue][iCand][iPID][0]->Fill(p,Ep);
	      //fHistTrue0PH0ClustM02[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M02);
	      //fHistTrue0PH0ClustM20[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),M20);
	      //fHistTrue0PH0ClustM02M20[iTrue][iCand][iPID][0]->Fill(M02,M20);
	      //fHistTrue0PH0ClustDx[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dx);
	      //fHistTrue0PH0ClustDz[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_dz);
	      //fHistTrue0PH0ClustDr[iTrue][iCand][iPID][0]->Fill(ClusterVector.E(),track_R);
	    }
	  }
	}
      }
    }

    }
  }//end of loop i 
  
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::NeutralPion()
{

  Double_t true_prim_vtx[3]={};
  if(fIsMC){
    AliAODMCParticle *particle  = (AliAODMCParticle*)AODMCTrackArray->At(0);
    particle->XvYvZv(true_prim_vtx);
  }

  AliAODCaloTrigger *fTrigData1 = (AliAODCaloTrigger*)fTrigData->Clone("fTrigData1");
  AliAODCaloTrigger *fTrigData2 = (AliAODCaloTrigger*)fTrigData->Clone("fTrigData2");
  AliAODCaloCells* cells = dynamic_cast<AliAODCaloCells*> (fAOD[0]->GetPHOSCells());
  
  Double_t vtx5[3] = {fVertexVector.X(),fVertexVector.Y(),fVertexVector.Z()};
  
  for (Int_t i=0; i<fAOD[0]->GetNumberOfCaloClusters(); i++) {
    AliAODCaloCluster *clu1 = dynamic_cast<AliAODCaloCluster*>(fAOD[0]->GetCaloCluster(i));
    fTrigData1->Reset();

    if( !clu1->IsPHOS()) continue;

    Float_t  position1[3];
    clu1->GetPosition(position1);
    TVector3 global1(position1) ;
    Int_t relId1[4] = {0,0,0,0};
    fPHOSGeo->GlobalPos2RelId(global1,relId1) ;
    Int_t    mod1 = relId1[0] ;
    Int_t    col1 = relId1[2] ;
    Int_t    row1 = relId1[3] ;
    Int_t    tru1 = GetTRUNum(col1,row1);

    TLorentzVector ClusterVector1;
    AliPHOSAodCluster cluPHOS1( *(AliAODCaloCluster*) (clu1) );
    cluPHOS1.SetE(AddEnergyCalib(mod1,cluPHOS1.E()));
    cluPHOS1.SetPosition(position1);
    cluPHOS1.GetMomentum(ClusterVector1,vtx5);

    if(!IsGoodCluster(ClusterVector1.E(),clu1->GetNCells(),mod1,col1,row1)) continue;

    Int_t  maxAbsId1 = 0;
    Int_t  maxrelId1[4];
    MaxEnergyCellPos(cells,clu1,maxAbsId1);
    fPHOSGeo->AbsToRelNumbering(maxAbsId1,maxrelId1);

    Double_t M02_1, M20_1, t1_1, t2_1, t3_1, t4_1, t5_1, t6_1;
    CalculatedClusterParameters(clu1, cells, AddEnergyCalib(mod1,ClusterVector1.E())/clu1->E(), M02_1, M20_1,t1_1, t2_1,t3_1, t4_1, t5_1, t6_1);

    M02_1 = cluPHOS1.GetM02();
    M20_1 = cluPHOS1.GetM20();

    Bool_t isDispOK1 = false;
    
    if(cluPHOS1.Chi2() < fDispCut*fDispCut)
      isDispOK1 = true;
    else
      isDispOK1 = false;
    
    Bool_t isTOFOK1 = false;

    if(fabs(t6_1) < fTOFcut)
      isTOFOK1 = true;
    else
      isTOFOK1 = false;

    Double_t track_dx1 = clu1->GetTrackDx();
    Double_t track_dz1 = clu1->GetTrackDz();
    
    Double_t track_R1  = 0;
    if(track_dx1>0 || track_dz1>0){
      track_R1 = sqrt(pow(track_dx1,2) + pow(track_dz1,2));
    }

    Bool_t isTrig1       = IsTrigger(fTrigData1,cells,maxAbsId1,ClusterVector1.E(),t6_1);
    
    Bool_t isTrackMatch1 = IsTrackCluster(clu1->GetNTracksMatched(),track_dx1,track_dz1,ClusterVector1.E());
    
    TString trackPID1 = "";
    Double_t Ep1 = -100;
    if(clu1->GetNTracksMatched()>0){
      AliAODTrack   *track1     = dynamic_cast<AliAODTrack*>(clu1->GetTrackMatched(0));
      Ep1 = ClusterVector1.E()/track1->P();
      trackPID1 = GetTrackPID(track1);
    }

    //TString nameTrue[] = {"_Photon","_NotPhoton","_Ele","_Proton","_AntiProton","_Neutron","_AntiNeutron","_Charged","_Neutral"};
    
    Int_t  PDGcode1     = 0;
    Int_t  clustCharge1 = 0;
    Int_t  momLabel1    = 0;
    Bool_t sure1 = false;
    
    Bool_t isPhoton1                   = false; // largest contribution to cluster is photon
    Bool_t isElectron1                 = false; // largest contribution to cluster is electron
    Bool_t isProton1                   = false; // largest contribution to cluster is proton
    Bool_t isAntiProton1               = false; // largest contribution to cluster is antiproton
    Bool_t isNeutron1                  = false; // largest contribution to cluster is neutron
    Bool_t isAntiNeutron1              = false; // largest contribution to cluster is antineutron
    Bool_t isConversion1               = false; // largest contribution to cluster is converted electron
    Bool_t isChargedParticle1          = false; // largest contribution to cluster is charged particle
    Bool_t isNeutralParticle1          = false; // largest contribution to cluster is neutral particle
    Bool_t isChargedPion1              = false; // largest contribution to cluster is charged pion
    Bool_t isRadius1cmCut1             = false;

    if(fIsMC){
      PDGcode1 = GetPDGCode(FindPrimary(clu1,sure1),clustCharge1,momLabel1,isRadius1cmCut1);
      this->SetClusterFlag(PDGcode1,clustCharge1,isPhoton1,isElectron1,isProton1,isAntiProton1,isNeutron1,isAntiNeutron1,isConversion1,isChargedParticle1,isNeutralParticle1,isChargedPion1);
    }

    for (Int_t j=i+1; j<fAOD[0]->GetNumberOfCaloClusters(); j++) {
      
      AliAODCaloCluster *clu2 = dynamic_cast<AliAODCaloCluster*>(fAOD[0]->GetCaloCluster(j));
      fTrigData2->Reset();
      
      if( !clu2->IsPHOS()) continue;

      Float_t  position2[3];
      clu2->GetPosition(position2);
      TVector3 global2(position2) ;
      Int_t relId2[4] = {0,0,0,0};
      fPHOSGeo->GlobalPos2RelId(global2,relId2) ;
      Int_t    mod2 = relId2[0] ;
      Int_t    col2 = relId2[2] ;
      Int_t    row2 = relId2[3] ;
      Int_t    tru2 = GetTRUNum(col2,row2);
      
      TLorentzVector ClusterVector2;
      AliPHOSAodCluster cluPHOS2( *(AliAODCaloCluster*) (clu2) );
      cluPHOS2.SetE(AddEnergyCalib(mod2,cluPHOS2.E()));
      cluPHOS2.SetPosition(position2);
      cluPHOS2.GetMomentum(ClusterVector2,vtx5);

      if(!IsGoodCluster(ClusterVector2.E(),clu2->GetNCells(),mod2,col2,row2)) continue;

      Int_t  maxAbsId2 = 0;
      Int_t  maxrelId2[4];
      MaxEnergyCellPos(cells,clu2,maxAbsId2);
      fPHOSGeo->AbsToRelNumbering(maxAbsId2,maxrelId2);

      Double_t M02_2, M20_2, t1_2, t2_2, t3_2, t4_2, t5_2, t6_2;
      CalculatedClusterParameters(clu2, cells, AddEnergyCalib(mod2,ClusterVector2.E())/clu2->E(), M02_2, M20_2,t1_2, t2_2,t3_2, t4_2, t5_2, t6_2);
      
      M02_2 = cluPHOS2.GetM02();
      M20_2 = cluPHOS2.GetM20();

      Bool_t isDispOK2 = false;

      if(cluPHOS2.Chi2() < fDispCut*fDispCut)
	isDispOK2 = true;
      else
	isDispOK2 = false;
      
      Bool_t isTOFOK2 = false;

      if(fabs(t6_2) < fTOFcut)
	isTOFOK2 = true;
      else
	isTOFOK2 = false;
      
      Double_t track_dx2 = clu2->GetTrackDx();
      Double_t track_dz2 = clu2->GetTrackDz();

      Double_t track_R2  = 0;
      if(track_dx2>0 || track_dz2>0){
	track_R2 = sqrt(pow(track_dx2,2) + pow(track_dz2,2));
      }

      Bool_t isTrig2       = IsTrigger(fTrigData2,cells,maxAbsId2,ClusterVector2.E(),t6_2);
      
      TLorentzVector ClusterVector12 = ClusterVector1 + ClusterVector2;
      
      Bool_t isTrackMatch2 = IsTrackCluster(clu2->GetNTracksMatched(),track_dx2,track_dz2,ClusterVector2.E());

      TString trackPID2 = "";
      Double_t Ep2 = -100;
      if(clu2->GetNTracksMatched()>0){
	AliAODTrack   *track2     = dynamic_cast<AliAODTrack*>(clu2->GetTrackMatched(0));
	Ep2 = ClusterVector2.E()/track2->P();
	trackPID2 = GetTrackPID(track2);
      }
      
      Int_t  PDGcode2     = 0;
      Int_t  clustCharge2 = 0;
      Int_t  momLabel2    = 0;
      Bool_t sure2 = false;
      
      Bool_t isPhoton2                   = false; // largest contribution to cluster is photon
      Bool_t isElectron2                 = false; // largest contribution to cluster is electron
      Bool_t isProton2                   = false; // largest contribution to cluster is proton
      Bool_t isAntiProton2               = false; // largest contribution to cluster is antiproton
      Bool_t isNeutron2                  = false; // largest contribution to cluster is neutron
      Bool_t isAntiNeutron2              = false; // largest contribution to cluster is antineutron
      Bool_t isConversion2               = false; // largest contribution to cluster is converted electron
      Bool_t isChargedParticle2          = false; // largest contribution to cluster is charged particle
      Bool_t isNeutralParticle2          = false; // largest contribution to cluster is neutral particle
      Bool_t isChargedPion2              = false; // largest contribution to cluster is charged pion
      Bool_t isRadius1cmCut2             = false;
      
      if(fIsMC){
	PDGcode2 = GetPDGCode(FindPrimary(clu2,sure2),clustCharge2,momLabel2,isRadius1cmCut2);
	this->SetClusterFlag(PDGcode2,clustCharge2,isPhoton2,isElectron2,isProton2,isAntiProton2,isNeutron2,isAntiNeutron2,isConversion2,isChargedParticle2,isNeutralParticle2,isChargedPion2);
      }

      fHistMassTwoGammas[0][0][0][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
      fHistMassTwoGammasMultiUnitRap[0][0][0][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
      fHistMassTwoGammasMultiGapRap1[0][0][0][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
      fHistMassTwoGammasMultiGapRap2[0][0][0][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
      fHistMassTwoGammasMultiV0AC[0][0][0][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
      
      if(isTOFOK1 && isTOFOK2){
	fHistMassTwoGammas[0][0][2][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMassTwoGammasMultiUnitRap[0][0][2][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMassTwoGammasMultiGapRap1[0][0][2][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMassTwoGammasMultiGapRap2[0][0][2][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMassTwoGammasMultiV0AC[0][0][2][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
      }
      if(isTOFOK1 || isTOFOK2){
	fHistMassTwoGammas[0][0][1][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMassTwoGammasMultiUnitRap[0][0][1][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMassTwoGammasMultiGapRap1[0][0][1][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMassTwoGammasMultiGapRap2[0][0][1][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMassTwoGammasMultiV0AC[0][0][1][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
      }
      
      if(mod1==mod2){
	fHistMassTwoGammas[0][0][0][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	if(isTOFOK1 && isTOFOK2){
	  fHistMassTwoGammas[0][0][2][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	}
	if(isTOFOK1 || isTOFOK2){
	  fHistMassTwoGammas[0][0][1][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	}
      }
      
      if(isTrig1 || isTrig2){
	
	fHistMassTwoGammas[1][0][0][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMassTwoGammasMultiUnitRap[1][0][0][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMassTwoGammasMultiGapRap1[1][0][0][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMassTwoGammasMultiGapRap2[1][0][0][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMassTwoGammasMultiV0AC[1][0][0][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());

	if(isTOFOK1 && isTOFOK2){
	  fHistMassTwoGammas[1][0][2][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiUnitRap[1][0][2][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiGapRap1[1][0][2][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiGapRap2[1][0][2][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiV0AC[1][0][2][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());

	}
	if(isTOFOK1 || isTOFOK2){
	  fHistMassTwoGammas[1][0][1][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiUnitRap[1][0][1][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiGapRap1[1][0][1][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiGapRap2[1][0][1][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiV0AC[1][0][1][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());

	}
	
	if(mod1==mod2){
	  fHistMassTwoGammas[1][0][0][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  if(isTOFOK1 && isTOFOK2){
	    fHistMassTwoGammas[1][0][2][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  }
	  if(isTOFOK1 || isTOFOK2){
	    fHistMassTwoGammas[1][0][1][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  }
	  
	}
      }
      
      if((isTrig1 && fIsGoodTRU[mod1-1][tru1-1]) || (isTrig2 && fIsGoodTRU[mod2-1][tru2-1])){
	fHistMassTwoGammas[1][1][0][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMassTwoGammasMultiUnitRap[1][1][0][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
        fHistMassTwoGammasMultiGapRap1[1][1][0][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
        fHistMassTwoGammasMultiGapRap2[1][1][0][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
        fHistMassTwoGammasMultiV0AC[1][1][0][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());

	if(isTOFOK1 && isTOFOK2){
	  fHistMassTwoGammas[1][1][2][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiUnitRap[1][1][2][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiGapRap1[1][1][2][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiGapRap2[1][1][2][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiV0AC[1][1][2][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	}
	if(isTOFOK1 || isTOFOK2){
	  fHistMassTwoGammas[1][1][1][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiUnitRap[1][1][1][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiGapRap1[1][1][1][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiGapRap2[1][1][1][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMassTwoGammasMultiV0AC[1][1][1][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	}
	
	if(mod1==mod2){
	  fHistMassTwoGammas[1][1][0][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  if(isTOFOK1 && isTOFOK2){
	    fHistMassTwoGammas[1][1][2][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  }
	  if(isTOFOK1 || isTOFOK2){
	    fHistMassTwoGammas[1][1][1][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  }
	}
      }

      Bool_t isSameMom         = false;
      Bool_t isSamePi0         = false;
      Bool_t isSameCharged     = false;
      Bool_t isSamePi0K0s      = false;
      Bool_t isSamePi0Material = false;

      Int_t iTrue = 0;
      
      if(fIsMC){

	if(momLabel1 == momLabel2){

	  isSameMom = true;
	  
	  AliAODMCParticle *particle  = (AliAODMCParticle*)AODMCTrackArray->At(momLabel1);
	  if(!particle) continue;
	  Double_t pvtx[3]={};
	  particle->XvYvZv(pvtx);
	  Double_t R = sqrt(pow(true_prim_vtx[0]-pvtx[0],2) + pow(true_prim_vtx[1]-pvtx[1],2))*10;
	  
	  if(particle){
	    if(isPhoton1 && isPhoton2 && particle->GetPdgCode()==111){
	      if(R < 1.){
		isSamePi0 = true;
		iTrue     = 1;
	      }
	      else if(particle->IsSecondaryFromWeakDecay() == true){
		Int_t gmamLabel = particle->GetMother();
		AliAODMCParticle *gparticle  = (AliAODMCParticle*)AODMCTrackArray->At(gmamLabel);
		if(gparticle){
		  if(gparticle->GetPdgCode() == 310){
		    isSamePi0K0s = true;
		    iTrue        = 3;
		  }
		}
	      }
	      else if(particle->IsSecondaryFromMaterial() == true){
		isSamePi0Material = true;
		iTrue             = 4;
		fHistRadiusPi0->Fill(sqrt(pow(pvtx[0],2) + pow(pvtx[1],2)));
		fHistRadius2DRecPi0->Fill(pvtx[0],pvtx[1]);
	      }
	      fHistTrueMassTwoGammas[0][0][0][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	    }//isPhoton1 && isPhoton2 && particle->GetPdgCode()==111
	    
	    if(isChargedParticle1 && isChargedParticle2){
	      isSameCharged = true;
	      iTrue         = 2;
	    }

	  }//particle
	  
	  if(iTrue==1 || iTrue==2 || iTrue==3 || iTrue==4){
	    
	    fHistTrueMassTwoGammas[0][0][0][iTrue][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	    if(mod1==mod2){
	      fHistTrueMassTwoGammas[0][0][iTrue][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	    }
	    
	    if(isTrig1 || isTrig2){
	      fHistTrueMassTwoGammas[1][0][0][iTrue][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	      if(mod1==mod2){
		fHistTrueMassTwoGammas[1][0][0][iTrue][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	      }
	    }
	    
	    if((isTrig1 && fIsGoodTRU[mod1-1][tru1-1]) || (isTrig2 && fIsGoodTRU[mod2-1][tru2-1])){
	      fHistTrueMassTwoGammas[1][1][0][iTrue][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	      if(mod1==mod2){
		fHistTrueMassTwoGammas[1][1][0][iTrue][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	      }
	    }
	  }//iTrue==1 || iTrue==2 || iTrue==3 || iTrue==4
	}//momLabel1 == momLabel
      }// fIsMC
    }

    for (Int_t j=0; j<fPHOSClusterList[fMultiBinSPDUnitRap][fZvtxBin]->GetSize(); j++) {

      AliAODCaloCluster *clu2 = (AliAODCaloCluster*)fPHOSClusterList[fMultiBinSPDUnitRap][fZvtxBin]->At(j);

      if( !clu2->IsPHOS()) continue;

      Float_t  position2[3];
      clu2->GetPosition(position2);
      TVector3 global2(position2) ;
      Int_t relId2[4] = {0,0,0,0};
      fPHOSGeo->GlobalPos2RelId(global2,relId2) ;
      Int_t    mod2 = relId2[0] ;
      Int_t    col2 = relId2[2] ;
      Int_t    row2 = relId2[3] ;
      Int_t    tru2 = GetTRUNum(col2,row2);

      TLorentzVector ClusterVector2;
      AliPHOSAodCluster cluPHOS2( *(AliAODCaloCluster*) (clu2) );
      cluPHOS2.SetE(AddEnergyCalib(mod2,cluPHOS2.E()));
      cluPHOS2.SetPosition(position2);
      cluPHOS2.GetMomentum(ClusterVector2,vtx5);

      if(!IsGoodCluster(ClusterVector2.E(),clu2->GetNCells(),mod2,col2,row2)) continue;

      Double_t M02_2, M20_2, t1_2, t2_2, t3_2, t4_2, t5_2, t6_2;
      t6_2  = clu2->GetTOF();
      M02_2 = cluPHOS2.GetM02();
      M20_2 = cluPHOS2.GetM20();

      Bool_t isDispOK2 = false;

      if(cluPHOS2.Chi2() < fDispCut*fDispCut)
	isDispOK2 = true;
      else
	isDispOK2 = false;
      
      Bool_t isTOFOK2 = false;

      if(fabs(t6_2) < fTOFcut)
	isTOFOK2 = true;
      else
	isTOFOK2 = false;
      
      Double_t track_dx2 = clu2->GetTrackDx();
      Double_t track_dz2 = clu2->GetTrackDz();

      Double_t track_R2  = 0;
      if(track_dx2>0 || track_dz2>0){
	track_R2 = sqrt(pow(track_dx2,2) + pow(track_dz2,2));
      }
      
      Bool_t isTrig2       = false;
      if(clu2->GetID() == 1)
	isTrig2 = true;
      
      TLorentzVector ClusterVector12 = ClusterVector1 + ClusterVector2;

      Bool_t isTrackMatch2 = IsTrackCluster(clu2->GetNTracksMatched(),track_dx2,track_dz2,ClusterVector2.E());

      //[iTrig][iTRU][iTOF][iPID][iMod] 
      fHistMixMassTwoGammas[0][0][0][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
      fHistMixMassTwoGammasMultiUnitRap[0][0][0][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
      fHistMixMassTwoGammasMultiGapRap1[0][0][0][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
      fHistMixMassTwoGammasMultiGapRap2[0][0][0][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
      fHistMixMassTwoGammasMultiV0AC[0][0][0][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
      
      if(isTOFOK1 && isTOFOK2){
	fHistMixMassTwoGammas[0][0][2][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiUnitRap[0][0][2][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiGapRap1[0][0][2][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiGapRap2[0][0][2][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiV0AC[0][0][2][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
      }
      if(isTOFOK1 || isTOFOK2){
	fHistMixMassTwoGammas[0][0][1][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiUnitRap[0][0][1][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiGapRap1[0][0][1][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiGapRap2[0][0][1][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiV0AC[0][0][1][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
      }
      
      if(mod1==mod2){
	fHistMixMassTwoGammas[0][0][0][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	if(isTOFOK1 && isTOFOK2){
	  fHistMixMassTwoGammas[0][0][2][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	}
	if(isTOFOK1 || isTOFOK2){
	  fHistMixMassTwoGammas[0][0][1][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	}
      }

      if(isTrig1 || isTrig2){
	fHistMixMassTwoGammas[1][0][0][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiUnitRap[1][0][0][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiGapRap1[1][0][0][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiGapRap2[1][0][0][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiV0AC[1][0][0][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	if(isTOFOK1 && isTOFOK2){
	  fHistMixMassTwoGammas[1][0][2][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiUnitRap[1][0][2][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiGapRap1[1][0][2][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiGapRap2[1][0][2][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiV0AC[1][0][2][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	}
	if(isTOFOK1 || isTOFOK2){
	  fHistMixMassTwoGammas[1][0][1][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiUnitRap[1][0][1][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiGapRap1[1][0][1][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiGapRap2[1][0][1][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiV0AC[1][0][1][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	}
	
	if(mod1==mod2){
	  fHistMixMassTwoGammas[1][0][0][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  if(isTOFOK1 && isTOFOK2){
	    fHistMixMassTwoGammas[1][0][2][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  }
	  if(isTOFOK1 || isTOFOK2){
	    fHistMixMassTwoGammas[1][0][1][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  }

	}
      }

      if((isTrig1 && fIsGoodTRU[mod1-1][tru1-1]) || (isTrig2 && fIsGoodTRU[mod2-1][tru2-1])){
	fHistMixMassTwoGammas[1][1][0][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiUnitRap[1][1][0][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiGapRap1[1][1][0][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiGapRap2[1][1][0][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	fHistMixMassTwoGammasMultiV0AC[1][1][0][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	if(isTOFOK1 && isTOFOK2){
	  fHistMixMassTwoGammas[1][1][2][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiUnitRap[1][1][2][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiGapRap1[1][1][2][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiGapRap2[1][1][2][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiV0AC[1][1][2][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	}
	if(isTOFOK1 || isTOFOK2){
	  fHistMixMassTwoGammas[1][1][1][0][0]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiUnitRap[1][1][1][fMultiBinSPDUnitRap]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiGapRap1[1][1][1][fMultiBinSPDGapRap1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiGapRap2[1][1][1][fMultiBinSPDGapRap2]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  fHistMixMassTwoGammasMultiV0AC[1][1][1][fMultiBinV0AC]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	}
	
	if(mod1==mod2){
	  fHistMixMassTwoGammas[1][1][0][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  if(isTOFOK1 && isTOFOK2){
	    fHistMixMassTwoGammas[1][1][2][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  }
	  if(isTOFOK1 || isTOFOK2){
	    fHistMixMassTwoGammas[1][1][1][0][mod1]->Fill(ClusterVector12.M(),ClusterVector12.Pt());
	  }

	}
      }

    }
  }
  
}
//________________________________________________________________________
Double_t AliAnalysisTaskPHOSTrigPi0::AddEnergyCalib(Int_t mod,Double_t ene)
{
  
  Double_t calib=1.;
  if(!fIsMC){
    
    if(fPeriod=="LHC12c"){
      if(mod==1) calib = 0.99881*0.99697;
      else       calib = 0.99887*0.98259;
    }
    if(fPeriod=="LHC12d"){
      if(mod==1) calib = 1.00835;
      else       calib = 0.99660;
    }
    if(fPeriod=="LHC12h"){
      if(mod==1) calib = 1.00543;
      else       calib = 0.99522;
    }
    return ene * calib;
  }  
  else{
    
    Double_t ParA=1; Double_t ParB=1;
    Double_t iw1=1;  Double_t iw2=1;

    if(mod==1){
      calib=0.979*1.00478145626526683;
      ParA = 0.000 + 6*0.005;
      ParB = 0.1   + 5*0.1;
      iw1 = 15;
      iw2 = 9;
    }
    else if(mod==3){
      calib=0.980*0.996980394768978040;
      ParA = 0.000 + 7*0.005;
      ParB = 0.1   + 5*0.1;
      iw1 = 7;
      iw2 = 8;
    }
  
    fEneSmearMeanMC->SetParameters(0.05+iw1*0.015, -1.0-iw2*0.2, 0.011);
    Double_t mean  = (1+ParA/(1+pow(ene/ParB,2))) * ene *calib;
    fEneSmearSigmaMC->SetParameter(2,fEneSmearMeanMC->Eval(mean));
    Double_t smear = fEneSmearSigmaMC->GetRandom();
    
    return mean*smear;
  }
  
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPHOSTrigPi0::IsGoodCluster(Double_t ene, Int_t ncell, Int_t mod, Int_t ix, Int_t iz)
{
  
  if(ene<fMinEene)                       return false;
  if(ncell<fMinNCell)                    return false;
  if(!IsGoodChannel("PHOS",mod,ix,iz))   return false;
  if(!fIsGoodMod[mod-1])                 return false;

  return true;

}
//________________________________________________________________________
 Bool_t AliAnalysisTaskPHOSTrigPi0::IsTrackCluster(Int_t nTrack,Double_t dx, Double_t dz, Double_t ene)
{
  
  if(nTrack>0){
    if(fabs(dx)<fHistTrackCutX->Eval(ene) && fabs(dz)<fHistTrackCutZ->Eval(ene)){
      return true;
    }
    else
      false;
  }
  else{
    return false;
  }

  
}
//___________________________________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::MaxEnergyCellPos(AliAODCaloCells *cells, AliAODCluster* clu, Int_t& maxId)
{
  Double_t eMax = -111;
  
  for (Int_t iDig=0; iDig< clu->GetNCells(); iDig++) {

    Int_t cellAbsId = clu->GetCellAbsId(iDig);
    Int_t cellRelId[4];

    Double_t eCell  = cells->GetCellAmplitude(cellAbsId)*clu->GetCellAmplitudeFraction(iDig);

    fPHOSGeo->AbsToRelNumbering(cellAbsId,cellRelId) ;

    if(eCell>eMax) {
      eMax = eCell;
      maxId = cellAbsId;
    }
  }

}

//________________________________________________________________________
Bool_t AliAnalysisTaskPHOSTrigPi0::IsTrigger(AliAODCaloTrigger *trigData, AliAODCaloCells *cells, Int_t maxAbsId, Double_t cluE, Double_t tof)
{
  Int_t maxrelId[4] = {};
  fPHOSGeo->AbsToRelNumbering(maxAbsId, maxrelId);
  Int_t tmod   = 0;
  Int_t tAbsId = 0;
  Int_t trelId[4] = {};
  Bool_t fTrig = false;
  
  Int_t relIdAnalogSum[4]={};
  Int_t AbsIdAnalogSum=0;
  
  while(trigData->Next()) {
    trigData->GetPosition(tmod,tAbsId);
    fPHOSGeo->AbsToRelNumbering(tAbsId,trelId);

    if (IsFireTile(trelId,maxrelId)){
      fTrig = true;
    } // Is fire trigger
  }
  
  return fTrig;

}
//________________________________________________________________________
Bool_t AliAnalysisTaskPHOSTrigPi0::IsFireTile(Int_t *trig_relid, Int_t *cluster_relid)
{
  if( trig_relid[0] != cluster_relid[0] )            return kFALSE; // different modules!
  if( cluster_relid[2]-trig_relid[2]>3 || cluster_relid[2]-trig_relid[2]<0) return kFALSE; // X-distance too large!
  if( cluster_relid[3]-trig_relid[3]>3 || cluster_relid[3]-trig_relid[3]<0) return kFALSE; // Z-distance too large!
  return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::CalculatedClusterParameters(AliAODCaloCluster* clu, AliAODCaloCells* cells, Double_t calib,
							     Double_t& CalcM20, Double_t& CalcM02, Double_t& time1, Double_t& time2, Double_t& time3,Double_t& time4, Double_t& time5, Double_t& time6){
  
     const Double_t rCut = 4.5;
     
     TVector3 iDigitGlobPos;
     Int_t    iDigitRelId[4];
     Float_t xi=0;
     Float_t zi=0;
     Double_t *elist     = clu->GetCellsAmplitudeFraction();
     Double_t recal_sumE = 0;
     
     for(Int_t iDigit=0; iDigit<clu->GetNCells(); ++iDigit){
       Int_t absId = clu->GetCellAbsId(iDigit);
       fPHOSGeo->AbsToRelNumbering(absId,iDigitRelId);
       fPHOSGeo->RelPosInModule(iDigitRelId, xi, zi);
       Double_t ie = cells->GetCellAmplitude(absId) * elist[iDigit] * fPHOSCalibData->GetADCchannelEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]) * calib ;//Calculate real icell energy
       recal_sumE += ie;
     }

     Double_t meanX = 0;
     Double_t meanZ = 0;
     Double_t wtot  = 0;

     for(Int_t iDigit=0; iDigit<clu->GetNCells(); ++iDigit){
       Int_t absId = clu->GetCellAbsId(iDigit);
       fPHOSGeo->AbsToRelNumbering(absId,iDigitRelId);
       fPHOSGeo->RelPosInModule(iDigitRelId, xi, zi);
       Double_t ie = cells->GetCellAmplitude(absId) * elist[iDigit] * fPHOSCalibData->GetADCchannelEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]) * calib ;//Calculate real icell energy
       Double_t w = TMath::Max(0.,4.5 + TMath::Log(ie/recal_sumE));
       meanX += xi * w;
       meanZ += zi * w;
       wtot  += w;
     }
     if(wtot>0){
       meanX /= wtot;
       meanZ /= wtot;
     }
     
     wtot = 0;

     Double_t d    = 0.;
     Double_t dxx  = 0.;
     Double_t dzz  = 0.;
     Double_t dxz  = 0.;
     Double_t xCut = 0.;
     Double_t zCut = 0.;
     
     for(Int_t iDigit=0; iDigit < clu->GetNCells(); iDigit++) {
       Float_t xi ;
       Float_t zi ;
       Int_t absId = clu->GetCellAbsId(iDigit);
       fPHOSGeo->AbsToRelNumbering(absId,iDigitRelId) ;
       fPHOSGeo->RelPosInModule(iDigitRelId, xi, zi);

       Double_t ie = cells->GetCellAmplitude(absId) * elist[iDigit] * fPHOSCalibData->GetADCchannelEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]) * calib ;//Calculate real icell energy
       
       if (recal_sumE > 0 && ie > 0) {
	 if((xi-meanX)*(xi-meanX)+(zi-meanZ)*(zi-meanZ)<rCut*rCut){
       
	   Float_t w = TMath::Max( 0., 4.5 + TMath::Log(ie/recal_sumE) ) ;
	   
	   d    += w * ( (xi-meanX)*(xi-meanX) + (zi-meanZ)*(zi-meanZ) ) ;
	   xCut += w * xi ;
	   zCut += w * zi ;
	   dxx  += w * xi * xi ;
	   dzz  += w * zi * zi ;
	   dxz  += w * xi * zi ;
	   wtot += w ;
	 
	 }
       }
     }
     
     Double_t M02=0;
     Double_t M20=0;
     
     if (wtot>0) {
       xCut/= wtot ;
       zCut/= wtot ;
       dxx /= wtot ;
       dzz /= wtot ;
       dxz /= wtot ;
       dxx -= xCut * xCut ;
       dzz -= zCut * zCut ;
       dxz -= xCut * zCut ;
       
       M02 = 0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz );
       M20 = 0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz );
     }
     else{
       d=0. ;
       M02=0;
       M20=0;
     }
     
     CalcM02 = M02;
     CalcM20 = M20;
     
     Int_t mulDigit=clu->GetNCells() ;
     
     Float_t eMax1=0.;
     Float_t eMax2=0.;
     Float_t eMax3=0.;
     Float_t eMax4=0.;
     Float_t tMax1=0.; //Time at the maximum
     Float_t tMax2=0.; //Time at the maximum
     Float_t tMax3=0.; //Time at the maximum
     Float_t tMax4=0.; //Time at the maximum
     
     for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
       Int_t absId=clu->GetCellAbsId(iDigit) ;
       fPHOSGeo->AbsToRelNumbering(absId,iDigitRelId);
       Double_t ie = cells->GetCellAmplitude(absId) * elist[iDigit] * fPHOSCalibData->GetADCchannelEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]);//Calculate real icell energy
       
       Bool_t isHG=false ;
       if(cells->GetCellHighGain(absId))
	 isHG=true;
       
       if( ie>eMax1){
	 Double_t tof1 = cells->GetCellTime(absId);
	 if(isHG){
	   tof1 -= fPHOSCalibData->GetTimeShiftEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]);
	   tof1 -= fTimingCalibMapHG[iDigitRelId[0]-1]->GetBinContent(iDigitRelId[2],iDigitRelId[3]);
	 }
	 else{
	   tof1 -= fPHOSCalibData->GetLGTimeShiftEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]);
	   tof1 -= fTimingCalibMapLG[iDigitRelId[0]-1]->GetBinContent(iDigitRelId[2],iDigitRelId[3]);
	 }
	 tMax1 = tof1;
	 eMax1 = ie ;
       }
       if( ie>eMax2){
	 Double_t tof2 = cells->GetCellTime(absId);
	 if(isHG){
	   tof2 += fPHOSCalibData->GetTimeShiftEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]);
	   tof2 -= fTimingCalibMapHG[iDigitRelId[0]-1]->GetBinContent(iDigitRelId[2],iDigitRelId[3]);
	 }
	 else{
	   tof2 += fPHOSCalibData->GetLGTimeShiftEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]);
	   tof2 -= fTimingCalibMapLG[iDigitRelId[0]-1]->GetBinContent(iDigitRelId[2],iDigitRelId[3]);
	 }
	 tMax2 = tof2;
	 eMax2 = ie ;
       }
       if( ie>eMax3){
	 Double_t tof3 = 0;
	 if(isHG){
	   tof3 = cells->GetCellTime(absId);
	   tof3 -= fPHOSCalibData->GetTimeShiftEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]);
	   tof3 -= fTimingCalibMapHG[iDigitRelId[0]-1]->GetBinContent(iDigitRelId[2],iDigitRelId[3]);
	   eMax3 = ie ;
	   tMax3 = tof3;
	 }
       }
       if( ie>eMax4){
	 Double_t tof4 = 0;
	 if(isHG){
	   tof4  = cells->GetCellTime(absId);
	   tof4 += fPHOSCalibData->GetTimeShiftEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]);
	   tof4 -= fTimingCalibMapHG[iDigitRelId[0]-1]->GetBinContent(iDigitRelId[2],iDigitRelId[3]);
	   eMax4 = ie ;
	   tMax4 = tof4;
	 }
       }

     }
     
     Double_t eMin1=TMath::Min(0.5,0.2*eMax1) ;
     Double_t eMin2=TMath::Min(0.5,0.2*eMax2) ;
     Double_t eMin3=TMath::Min(0.5,0.2*eMax3) ;
     Double_t eMin4=TMath::Min(0.5,0.2*eMax4) ;
     Float_t  wtot1 = 0.;
     Float_t  wtot2 = 0.;
     Float_t  wtot3 = 0.;
     Float_t  wtot4 = 0.;
     Float_t  wtot5 = 0.;
     Float_t  wtot6 = 0.;
     Double_t t1 = 0. ;
     Double_t t2 = 0. ;
     Double_t t3 = 0. ;
     Double_t t4 = 0. ;
     Double_t t5 = tMax3;
     Double_t t6 = tMax4 ;
     
     Int_t nt1 =0;
     Int_t nt2 =0;
     Int_t nt3 =0;
     Int_t nt4 =0;
     Int_t nt5 =0;
     Int_t nt6 =0;

     for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
       Int_t absId=clu->GetCellAbsId(iDigit) ;
       fPHOSGeo->AbsToRelNumbering(absId,iDigitRelId);
       Double_t ie = cells->GetCellAmplitude(absId) * elist[iDigit] * fPHOSCalibData->GetADCchannelEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]);//Calculate real icell energy
       Bool_t isHG=false ;
       if(cells->GetCellHighGain(absId))
	 isHG=true;
       Double_t tof1 = cells->GetCellTime(absId);
       Double_t tof2 = cells->GetCellTime(absId);
       Double_t tof3 = 0;
       Double_t tof4 = 0;

       if(isHG){
	 tof3 = cells->GetCellTime(absId);
	 tof4 = cells->GetCellTime(absId);
	 tof1 -= fPHOSCalibData->GetTimeShiftEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]);
	 tof1 -= fTimingCalibMapHG[iDigitRelId[0]-1]->GetBinContent(iDigitRelId[2],iDigitRelId[3]);
	 tof2 += fPHOSCalibData->GetTimeShiftEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]);
	 tof2 -= fTimingCalibMapHG[iDigitRelId[0]-1]->GetBinContent(iDigitRelId[2],iDigitRelId[3]);
	 tof3 -= fPHOSCalibData->GetTimeShiftEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]);
	 tof3 -= fTimingCalibMapHG[iDigitRelId[0]-1]->GetBinContent(iDigitRelId[2],iDigitRelId[3]);
	 tof4 += fPHOSCalibData->GetTimeShiftEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]);
	 tof4 -= fTimingCalibMapHG[iDigitRelId[0]-1]->GetBinContent(iDigitRelId[2],iDigitRelId[3]);
       }
       else{
	 tof1 -= fPHOSCalibData->GetLGTimeShiftEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]);
	 tof1 -= fTimingCalibMapLG[iDigitRelId[0]-1]->GetBinContent(iDigitRelId[2],iDigitRelId[3]);
	 tof2 += fPHOSCalibData->GetLGTimeShiftEmc(iDigitRelId[0],iDigitRelId[3],iDigitRelId[2]);
	 tof2 -= fTimingCalibMapLG[iDigitRelId[0]-1]->GetBinContent(iDigitRelId[2],iDigitRelId[3]);
       }
       
       if(TMath::Abs(tof1-tMax1)<50.e-9){
	 //Remove too soft cells
	 if(ie<eMin1)
	   continue ;
	 if(ie>0){
	   //weight = 1./sigma^2
	   //Sigma is parameterization of TOF resolution 16.05.2013
	   Double_t wi1=0.;
	   if(isHG)
	     wi1=1./(2.4e-9 + 3.9e-9/ie) ;
	   else
	     wi1=1./(2.4e-9 + 3.9e-9/(0.1*ie)) ; //E of LG digit is 1/16 of correcponding HG
	   t1    += tof1*wi1 ;
	   wtot1 += wi1 ;
	 }
       }
       
       if(TMath::Abs(tof2-tMax2)<50.e-9){
	 //Remove too soft cells
	 if(ie<eMin2)
	   continue ;
	 if(ie>0){
	   //weight = 1./sigma^2
	   //Sigma is parameterization of TOF resolution 16.05.2013
	   Double_t wi2=0.;
	   if(isHG)
	     wi2=1./(2.4e-9 + 3.9e-9/ie) ;
	   else
	     wi2=1./(2.4e-9 + 3.9e-9/(0.1*ie)) ; //E of LG digit is 1/16 of correcponding HG
	   t2    += tof2*wi2 ;
	   wtot2 += wi2 ;
	 }
       }
       
       if(TMath::Abs(tof3-tMax3)<50.e-9){
	 //Remove too soft cells
	 if(ie<eMin3)
	   continue ;
	 if(ie>0){
	   //weight = 1./sigma^2
	   //Sigma is parameterization of TOF resolution 16.05.2013
	   Double_t wi3=0.;
	   if(isHG){
	     wi3=1./(2.4e-9 + 3.9e-9/ie) ;
	   }
	   else{
	     continue;
	     //wi2=1./(2.4e-9 + 3.9e-9/(0.1*ie)) ; //E of LG digit is 1/16 of correcponding HG
	   }	   
	   
	   t3    += tof3*wi3 ;
	   wtot3 += wi3 ;
	 }
       }

       if(isHG){
	 t4    += tof4*ie;
	 wtot4 += ie;
	 ++nt4;
       }

     }     
     
     if(wtot1>0){
       t1=t1/wtot1 ;
     }
     if(wtot2>0){
       t2=t2/wtot2 ;
     }
     if(wtot3>0){
       t3=t3/wtot3 ;
     }
     if(wtot4>0){
       t4=t4/wtot4/nt4 ;
     }
     
     time1 = t1;
     time2 = t2;
     time3 = t3;
     time4 = t4;
     time5 = t5;
     time6 = t6;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPHOSTrigPi0::IsGoodChannel(const char * det, Int_t mod, Int_t ix, Int_t iz)
{
  
  if(fBadMap[mod-1]->GetBinContent(ix,iz))
    return false;
  else
    return true;
  
}
//________________________________________________________________________                                                                                                                                                                   
Bool_t AliAnalysisTaskPHOSTrigPi0::IsGoodTile(const char * det, Int_t mod, Int_t ix, Int_t iz)
{
  
  //return fIsGoodTile[mod-1][ix-1][iz-1];

}
//_______________________________________________________________________________
Int_t AliAnalysisTaskPHOSTrigPi0::GetTRUNum(Int_t cellX, Int_t cellZ)
{
  //Return TRU region number for given cell.
  //cellX: [1-64], cellZ: [1-56]
  Int_t iTRU=-111;
  //RCU0: TRU 1,2  
  if(1<=cellX&&cellX<=16) {
    if(1<=cellZ&&cellZ<=28) iTRU=2;
    else iTRU=1;
  }
  //RCU1: TRU 3,4
  if(17<=cellX&&cellX<=32) {
    if(1<=cellZ&&cellZ<=28) iTRU=4;
    else iTRU=3;
  }
  //RCU2: TRU 5,6
  if(33<=cellX&&cellX<=48) {
    if(1<=cellZ&&cellZ<=28) iTRU=6;
    else iTRU=5;
  }
  //RCU3: TRU 7,8
  if(49<=cellX&&cellX<=64) {
    if(1<=cellZ&&cellZ<=28) iTRU=8;
    else iTRU=7;
  }

  return iTRU;
}
//_______________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::EventFlagInit()
{
  
  fIsPileUp                   = false;
  fIsVtxOut10cm               = false;
  f0PH0Event                  = false;
  f0PH0Event_HighEnergy       = false;
  fTOFcut0PH0Event            = false;
  fTOFcut0PH0Event_HighEnergy = false;
  fMCTracksUnit = 0;
  fMCTracksGap1 = 0;
  fMCTracksGap2 = 0;
  fMCTracksV0AC = 0;
  fMultiBinSPDUnitRap = 0;
  fMultiBinSPDGapRap1 = 0;
  fMultiBinSPDGapRap2 = 0;
  fMultiBinV0AC       = 0;
  fZvtxBin            = 0;
  
}
//_______________________________________________________________________________
TString AliAnalysisTaskPHOSTrigPi0::GetTrackPID(AliAODTrack *track)
{
  
  if(!fPIDResponse) return "";

  Double_t CrossRowsTPC = track->GetTPCNCrossedRows();
  Double_t nSignalTPC   = track->GetTPCsignalN();

  //if(nSignalTPC<80) return "";
  //if(nSignalTPC/CrossRowsTPC<0.6) return "";

  Double_t nSigmaTPC_Electron = fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType)AliPID::kElectron);
  Double_t nSigmaTPC_Proton = fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType)AliPID::kProton);

  if(-1<nSigmaTPC_Electron && nSigmaTPC_Electron<3)
    return "Electron";
  else if(nSigmaTPC_Electron < -3.5)
    return "ChargedHadron";
  else 
    return ""; 
  
}
//________________________________________________________________________
Int_t AliAnalysisTaskPHOSTrigPi0::FindPrimary(AliAODCluster*clu,  Bool_t&sure){
  
  //Finds primary and estimates if it unique one?                                                                                                                                
  //First check can it be photon/electron                                                                                                                                          
  const Double_t emFraction=0.9; //part of energy of cluster to be assigned to EM particle                                                                                          
  Int_t nLabel=clu->GetNLabels() ;
    
  for(Int_t i=0;  i<nLabel;  i++){
    AliAODMCParticle *particle  = (AliAODMCParticle*)AODMCTrackArray->At(clu->GetLabelAt(i));
    if(!particle) continue;
    Int_t pdg = particle->GetPdgCode() ;
    if(pdg==22  ||  pdg==11 || pdg == -11){
      if(particle->E() > emFraction*clu->E()){
	sure=kTRUE ;
	return clu->GetLabelAt(i);
      }
    }
  }
  
  Double_t*  Ekin=  new  Double_t[nLabel] ;
  for(Int_t i=0;  i<nLabel;  i++){
    AliAODMCParticle *particle  = (AliAODMCParticle*)AODMCTrackArray->At(clu->GetLabelAt(i));
    if(!particle) continue;
    Ekin[i]=particle->P() ;  // estimate of kinetic energy                                                                                                                                  
    if(particle->GetPdgCode()==-2212  ||  particle->GetPdgCode()==-2112){
      Ekin[i]+=1.8  ;  //due to annihilation                                                                                                                                         
    }
  }

  Int_t iMax=0;
  Double_t eMax=0., eSubMax=0. ;
  for(Int_t i=0;  i<nLabel;  i++){
    if(Ekin[i]>eMax){
      eSubMax=eMax;
      eMax=Ekin[i];
      iMax=i;
    }
  }

  if(eSubMax>0.8*eMax)//not obvious primary                                                                                                                                          
    sure=kFALSE;
  else
    sure=kTRUE;
  delete[]  Ekin;
  return  clu->GetLabelAt(iMax) ;

  
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::GetMCStack()
{
  AODMCTrackArray = dynamic_cast<TClonesArray*>(fAOD[0]->FindListObject(AliAODMCParticle::StdBranchName()));
}
//________________________________________________________________________
Int_t AliAnalysisTaskPHOSTrigPi0::GetPDGCode(Int_t iLabel, Int_t &clustCharge, Int_t &momLabel, Bool_t &isRadius1cmCut)
{
  
  AliAODMCParticle *primary  = (AliAODMCParticle*)AODMCTrackArray->At(0);
  Double_t prim_vtx[3] = {primary->Xv(),primary->Yv(),primary->Zv()};

  AliAODMCParticle *particle  = (AliAODMCParticle*)AODMCTrackArray->At(iLabel);
  if(!particle) return 0;
  clustCharge = particle->Charge();
  momLabel    = particle->GetMother();
  
  Double_t Vx = prim_vtx[0] - particle->Xv();
  Double_t Vy = prim_vtx[1] - particle->Yv();
  Double_t R = 0;
  if( (Vx*Vx + Vy*Vy) > 0 )
    R = sqrt(Vx*Vx + Vy*Vy);
  else
    R = 0;
  if(R<1.)
    isRadius1cmCut = true;
  else
    isRadius1cmCut = false;

  return particle->GetPdgCode();

}
//________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::SetClusterFlag(Int_t PDGcode,Int_t clustCharge,Bool_t &isPhoton,Bool_t &isElectron,Bool_t &isProton,
						Bool_t &isAntiProton,Bool_t &isNeutron,Bool_t &isAntiNeutron,Bool_t &isConversion,Bool_t &isChargedParticle,Bool_t &isNeutralParticle, Bool_t &isChargedPion)
{
  if(PDGcode == 22){
    isPhoton = true;
  }
  else if(fabs(PDGcode) == 11){
    isElectron = true;
  }
  else if(PDGcode == 2212){
    isProton = true;
  }
  else if(PDGcode == -2212){
    isAntiProton = true;
  }
  else if(PDGcode == 2112){
    isNeutron = true;
  }
  else if(PDGcode == -2112){
    isAntiNeutron = true;
  }
  else if(PDGcode == 211){
    isChargedPion = true;
  }

  if(clustCharge != 0){
    isChargedParticle = true;
  }
  else{
    isNeutralParticle = true;
  }
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::SetMultiplicity()
{

  AliAODTracklets *trackletAOD = (AliAODTracklets*)fAOD[0]->GetTracklets();
  
  Int_t tot_Tracklet_unitRap    = 0;
  Int_t tot_Tracklet_GapRap1 = 0;
  Int_t tot_Tracklet_GapRap2 = 0;
  
  for(Int_t iTracklet=0; iTracklet<trackletAOD->GetNumberOfTracklets(); ++iTracklet){

    Double_t eta = -1*TMath::Log(TMath::Tan(trackletAOD->GetTheta(iTracklet)/2.));
    if(fabs(eta)<1.0){
      ++tot_Tracklet_unitRap;
    }
    if(0.5<fabs(eta) && fabs(eta)<1.0){
      ++tot_Tracklet_GapRap1;
    }
    if(0.13<fabs(eta) && fabs(eta)<1.0){
      ++tot_Tracklet_GapRap2;
    }

  }

  fHistSPDTrackletsUnitRap->Fill(fVertexVector.Z(),tot_Tracklet_unitRap);
  fHistSPDTrackletsGapRap1->Fill(fVertexVector.Z(),tot_Tracklet_GapRap1);
  fHistSPDTrackletsGapRap2->Fill(fVertexVector.Z(),tot_Tracklet_GapRap2);
  
  AliAODVZERO * aodV0 = fAOD[0]->GetVZEROData();
  Double_t multiV0A   = aodV0->GetMTotV0A();
  Double_t multiV0C   = aodV0->GetMTotV0C();
  Double_t multiV0AC  = multiV0A+multiV0C;
  
  fHistV0AMulti->Fill(fVertexVector.Z(),multiV0A);
  fHistV0CMulti->Fill(fVertexVector.Z(),multiV0C);
  fHistV0ACMulti->Fill(fVertexVector.Z(),multiV0AC);
  
  if(!fSPDMultiCorrUnit) return;

  Double_t MeanUnitRap = fSPDMultiCorrUnit->GetBinContent(fSPDMultiCorrUnit->GetXaxis()->FindBin(fVertexVector.Z()));
  Double_t delta_N_UnitRap = tot_Tracklet_unitRap * fRefSPDMultiUnit/MeanUnitRap - tot_Tracklet_unitRap;
  Int_t signUnitRap = (delta_N_UnitRap > 0.) ? 1 : -1;
  Int_t corrected_Tracklet_unitRap = TMath::Max(tot_Tracklet_unitRap + signUnitRap*fRandom->Poisson(TMath::Abs(delta_N_UnitRap)), 0);

  Double_t MeanGapRap1 = fSPDMultiCorrGap1->GetBinContent(fSPDMultiCorrGap1->GetXaxis()->FindBin(fVertexVector.Z()));
  Double_t delta_N_GapRap1 = tot_Tracklet_GapRap1 * fRefSPDMultiGap1/MeanGapRap1 - tot_Tracklet_GapRap1;
  Int_t signGapRap1 = (delta_N_GapRap1 > 0.) ? 1 : -1;
  Int_t corrected_Tracklet_GapRap1 = TMath::Max(tot_Tracklet_GapRap1 + signGapRap1*fRandom->Poisson(TMath::Abs(delta_N_GapRap1)), 0);

  Double_t MeanGapRap2 = fSPDMultiCorrGap2->GetBinContent(fSPDMultiCorrGap2->GetXaxis()->FindBin(fVertexVector.Z()));
  Double_t delta_N_GapRap2 = tot_Tracklet_GapRap2 * fRefSPDMultiGap2/MeanGapRap2 - tot_Tracklet_GapRap2;
  Int_t signGapRap2 = (delta_N_GapRap2 > 0.) ? 1 : -1;
  Int_t corrected_Tracklet_GapRap2 = TMath::Max(tot_Tracklet_GapRap2 + signGapRap2*fRandom->Poisson(TMath::Abs(delta_N_GapRap2)), 0);
  
  fHistCorrectedSPDTrackletsUnitRap->Fill(fVertexVector.Z(),corrected_Tracklet_unitRap);
  fHistCorrectedSPDTrackletsGapRap1->Fill(fVertexVector.Z(),corrected_Tracklet_GapRap1);
  fHistCorrectedSPDTrackletsGapRap2->Fill(fVertexVector.Z(),corrected_Tracklet_GapRap2);
  
  fHistCorrelationMCTrackSPDTrackletsUnitRap->Fill(fMCTracksUnit,tot_Tracklet_unitRap);
  fHistCorrelationMCTrackSPDTrackletsGapRap1->Fill(fMCTracksGap1,tot_Tracklet_GapRap1);
  fHistCorrelationMCTrackSPDTrackletsGapRap2->Fill(fMCTracksGap2,tot_Tracklet_GapRap2);

  fHistCorrelationMCTrackCorrectedSPDTrackletsUnitRap->Fill(fMCTracksUnit,corrected_Tracklet_unitRap);
  fHistCorrelationMCTrackCorrectedSPDTrackletsGapRap1->Fill(fMCTracksGap1,corrected_Tracklet_GapRap1);
  fHistCorrelationMCTrackCorrectedSPDTrackletsGapRap2->Fill(fMCTracksGap2,corrected_Tracklet_GapRap2);
  
  Double_t MeanV0AC     = fV0ACMultiCorr->GetBinContent(fV0ACMultiCorr->GetXaxis()->FindBin(fVertexVector.Z()));
  Double_t delta_N_V0AC = multiV0AC * fRefV0ACMulti/MeanV0AC - multiV0AC;
  Int_t signV0AC        = (delta_N_V0AC > 0.) ? 1 : -1;
  Double_t corrected_multiV0AC = TMath::Max(multiV0AC + signV0AC*fRandom->Poisson(TMath::Abs(delta_N_V0AC)), 0.);

  fHistCorrectedV0ACMulti->Fill(fVertexVector.Z(),corrected_multiV0AC);
  
  fHistCorrelationMCTrackV0AC->Fill(fMCTracksV0AC,multiV0AC);
  fHistCorrelationMCTrackCorrectedV0AC->Fill(fMCTracksV0AC,corrected_multiV0AC);

  fHistCorrelationSPDTrackletsV0AC->Fill(tot_Tracklet_unitRap, multiV0AC);
  fHistCorrelationCorrectedSPDTrackletsCorrectedV0AC->Fill(corrected_Tracklet_unitRap,corrected_multiV0AC);

  Double_t MultiFractionSPDUnitRap = corrected_Tracklet_unitRap/fMeanSPDUnitRap;
  Double_t MultiFractionSPDGapRap1 = corrected_Tracklet_GapRap1/fMeanSPDGapRap1;
  Double_t MultiFractionSPDGapRap2 = corrected_Tracklet_GapRap2/fMeanSPDGapRap2;
  Double_t MultiFractionV0AC       = corrected_multiV0AC/fMeanV0AC;
  //70 -100%
  if(0<corrected_multiV0AC && corrected_multiV0AC<=36.){
    fMultiBinV0AC=0;
  }
  //50 - 70%
  else if(36.<corrected_multiV0AC && corrected_multiV0AC<=58.){
    fMultiBinV0AC=1;
  }
  //30 - 50%
  else if(58.<corrected_multiV0AC && corrected_multiV0AC<=92.){
    fMultiBinV0AC=2;
  }
  //20 - 30%
  else if(92.<corrected_multiV0AC && corrected_multiV0AC<=118.){
    fMultiBinV0AC=3;
  }
  //10 - 20 %
  else if(118.<corrected_multiV0AC && corrected_multiV0AC<=160.){
    fMultiBinV0AC=4;
  }
  //5 - 10%
  else if(160.<corrected_multiV0AC && corrected_multiV0AC<=196.){
    fMultiBinV0AC=5;
  }
  //1 - 5%
  else if(196.<corrected_multiV0AC && corrected_multiV0AC<=270){
    fMultiBinV0AC=6;
  }
  //0.1 - 1.0%
  else if(270.<corrected_multiV0AC && corrected_multiV0AC<=360){
    fMultiBinV0AC=7;
  }
  //0.0 - 0.1%
  else if(360.<corrected_multiV0AC && corrected_multiV0AC<=504){
    fMultiBinV0AC=8;
  }
  //other 
  else{
    fMultiBinV0AC=9;
  }

  if(0<MultiFractionSPDUnitRap && MultiFractionSPDUnitRap<=0.5){
    fMultiBinSPDUnitRap=0;
  }
  else if(0.5<MultiFractionSPDUnitRap && MultiFractionSPDUnitRap<=1.5){
    fMultiBinSPDUnitRap=1;
  }
  else if(1.5<MultiFractionSPDUnitRap && MultiFractionSPDUnitRap<=2.5){
    fMultiBinSPDUnitRap=2;
  }
  else if(2.5<MultiFractionSPDUnitRap && MultiFractionSPDUnitRap<=3.5){
    fMultiBinSPDUnitRap=3;
  }
  else if(3.5<MultiFractionSPDUnitRap && MultiFractionSPDUnitRap<=4.5){
    fMultiBinSPDUnitRap=4;
  }
  else if(4.5<MultiFractionSPDUnitRap && MultiFractionSPDUnitRap<=5.5){
    fMultiBinSPDUnitRap=5;
  }
  else if(5.5<MultiFractionSPDUnitRap && MultiFractionSPDUnitRap<=6.5){
    fMultiBinSPDUnitRap=6;
  }
  else if(6.5<MultiFractionSPDUnitRap && MultiFractionSPDUnitRap<=7.5){
    fMultiBinSPDUnitRap=7;
  }
  else if(7.5<MultiFractionSPDUnitRap && MultiFractionSPDUnitRap<=8.5){
    fMultiBinSPDUnitRap=8;
  }
  else{
    fMultiBinSPDUnitRap=9;
  }

  if(0<MultiFractionSPDGapRap1 && MultiFractionSPDGapRap1<=0.5){
    fMultiBinSPDGapRap1=0;
  }
  else if(0.5<MultiFractionSPDGapRap1 && MultiFractionSPDGapRap1<=1.5){
    fMultiBinSPDGapRap1=1;
  }
  else if(1.5<MultiFractionSPDGapRap1 && MultiFractionSPDGapRap1<=2.5){
    fMultiBinSPDGapRap1=2;
  }
  else if(2.5<MultiFractionSPDGapRap1 && MultiFractionSPDGapRap1<=3.5){
    fMultiBinSPDGapRap1=3;
  }
  else if(3.5<MultiFractionSPDGapRap1 && MultiFractionSPDGapRap1<=4.5){
    fMultiBinSPDGapRap1=4;
  }
  else if(4.5<MultiFractionSPDGapRap1 && MultiFractionSPDGapRap1<=5.5){
    fMultiBinSPDGapRap1=5;
  }
  else if(5.5<MultiFractionSPDGapRap1 && MultiFractionSPDGapRap1<=6.5){
    fMultiBinSPDGapRap1=6;
  }
  else if(6.5<MultiFractionSPDGapRap1 && MultiFractionSPDGapRap1<=7.5){
    fMultiBinSPDGapRap1=7;
  }
  else if(7.5<MultiFractionSPDGapRap1 && MultiFractionSPDGapRap1<=8.5){
    fMultiBinSPDGapRap1=8;
  }
  else{
    fMultiBinSPDGapRap1=9;
  }

  if(0<MultiFractionSPDGapRap2 && MultiFractionSPDGapRap2<=0.5){
    fMultiBinSPDGapRap2=0;
  }
  else if(0.5<MultiFractionSPDGapRap2 && MultiFractionSPDGapRap2<=1.5){
    fMultiBinSPDGapRap2=1;
  }
  else if(1.5<MultiFractionSPDGapRap2 && MultiFractionSPDGapRap2<=2.5){
    fMultiBinSPDGapRap2=2;
  }
  else if(2.5<MultiFractionSPDGapRap2 && MultiFractionSPDGapRap2<=3.5){
    fMultiBinSPDGapRap2=3;
  }
  else if(3.5<MultiFractionSPDGapRap2 && MultiFractionSPDGapRap2<=4.5){
    fMultiBinSPDGapRap2=4;
  }
  else if(4.5<MultiFractionSPDGapRap2 && MultiFractionSPDGapRap2<=5.5){
    fMultiBinSPDGapRap2=5;
  }
  else if(5.5<MultiFractionSPDGapRap2 && MultiFractionSPDGapRap2<=6.5){
    fMultiBinSPDGapRap2=6;
  }
  else if(6.5<MultiFractionSPDGapRap2 && MultiFractionSPDGapRap2<=7.5){
    fMultiBinSPDGapRap2=7;
  }
  else if(7.5<MultiFractionSPDGapRap2 && MultiFractionSPDGapRap2<=8.5){
    fMultiBinSPDGapRap2=8;
  }
  else{
    fMultiBinSPDGapRap2=9;
  }

  fHistSPDMultiUnitRapV0ACEventClass[fMultiBinV0AC]->Fill(corrected_Tracklet_unitRap);
  fHistSPDMultiGapRap1V0ACEventClass[fMultiBinV0AC]->Fill(corrected_Tracklet_GapRap1);
  fHistSPDMultiGapRap2V0ACEventClass[fMultiBinV0AC]->Fill(corrected_Tracklet_GapRap2);

  fHistV0ACMultiSPDUnitRapEventClass[fMultiBinSPDUnitRap]->Fill(corrected_multiV0AC);
  fHistV0ACMultiSPDGapRap1EventClass[fMultiBinSPDGapRap1]->Fill(corrected_multiV0AC);
  fHistV0ACMultiSPDGapRap2EventClass[fMultiBinSPDGapRap2]->Fill(corrected_multiV0AC);
  
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::UpdataPHOSClusterPool()
{
  /*
  for (Int_t i=0; i<fAOD[0]->GetNumberOfCaloClusters(); i++) {

    AliAODCaloCluster *clu = (AliAODCaloCluster*)fAOD[0]->GetCaloCluster(i)->Clone();
    
    if( !clu->IsPHOS()) continue; // reject EMCal                                                                                                                                                                                            

    fPHOSClusterList[fMultiBin][fZvtxBin]->AddFirst(clu);
    
    if(fPHOSClusterList[fMultiBin][fZvtxBin]->GetEntries() > 50){
      fPHOSClusterList[fMultiBin][fZvtxBin]->RemoveLast() ;
    }
  }
  */
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTrigPi0::ProcessMC()
{
  
  if(!fIsMC) return;

  AliAODMCParticle *primary  = (AliAODMCParticle*)AODMCTrackArray->At(0);
  Double_t prim_vtx[3] = {primary->Xv(),primary->Yv(),primary->Zv()};

  Int_t nTrack = AODMCTrackArray->GetSize();

  for(Int_t iTrack=0;  iTrack<nTrack;  iTrack++){
    
    AliAODMCParticle *particle  = (AliAODMCParticle*)AODMCTrackArray->At(iTrack);
    if(!particle)                  continue;

    Double_t vtx[3]={};
    particle->XvYvZv(vtx);
    
    Double_t ene      = particle->E();
    Double_t Pt       = particle->Pt();
    Double_t rapidity = particle->Y();
    Double_t eta      = particle->Eta();
    Double_t phi      = particle->Phi()*180./TMath::Pi();
    
    Double_t Vx = prim_vtx[0] - vtx[0];
    Double_t Vy = prim_vtx[1] - vtx[1];
    Double_t R = 0;
    
    if( (Vx*Vx + Vy*Vy) > 0 )
      R = sqrt(Vx*Vx + Vy*Vy);
    else
      R = 0;
    
    if(fabs(eta)<0.5 && particle->IsPhysicalPrimary() && particle->Charge()!=0){
      ++fMCTracksUnit;
    }
    if(0.5<fabs(eta) && fabs(eta)<1.0 && particle->IsPhysicalPrimary() && particle->Charge()!=0){
      ++fMCTracksGap1;
    }
    if(0.13<fabs(eta) && fabs(eta)<1.0 && particle->IsPhysicalPrimary() && particle->Charge()!=0){
      ++fMCTracksGap2;
    }
    if( 2.8<eta && eta<5.1 && particle->IsPhysicalPrimary() && particle->Charge()!=0){
      ++fMCTracksV0AC;
    }
    if( -3.7<eta && eta<-1.7 && particle->IsPhysicalPrimary() && particle->Charge()!=0){
      ++fMCTracksV0AC;
    }
    
    if(R<1. && particle->GetPdgCode()==111){
      fHistPi0UnitRapV0AC[fMultiBinV0AC]->Fill(Pt);
      fHistPi0UnitRapSPDUnitRap[fMultiBinSPDUnitRap]->Fill(Pt);
      fHistPi0UnitRapSPDGapRap1[fMultiBinSPDGapRap1]->Fill(Pt);
      fHistPi0UnitRapSPDGapRap2[fMultiBinSPDGapRap2]->Fill(Pt);
      if(fabs(eta)<0.13){
	fHistPi0AcceptV0AC[fMultiBinV0AC]->Fill(Pt);
	fHistPi0AcceptSPDUnitRap[fMultiBinSPDUnitRap]->Fill(Pt);
	fHistPi0AcceptSPDGapRap1[fMultiBinSPDGapRap1]->Fill(Pt);
	fHistPi0AcceptSPDGapRap2[fMultiBinSPDGapRap2]->Fill(Pt);
      }
    }

    if(particle->GetPdgCode()!=22) continue;
    
    Int_t momTrack = particle->GetMother();
    AliAODMCParticle *mom_particle = NULL;
    mom_particle = (AliAODMCParticle*)AODMCTrackArray->At(momTrack);
    if(!mom_particle) continue;
    Int_t momPID = mom_particle->GetPdgCode();
    
    Bool_t isUnitRap = false;
    if(fabs(rapidity)<0.5)
      isUnitRap = true;
    Bool_t isAcceptPHOS = false;
    if(fabs(eta)<0.13 && 260<phi && phi<320)
      isAcceptPHOS = true;
    
    if(isUnitRap==true && momPID==111 && mom_particle->IsSecondaryFromMaterial()==true){
      fHistRadius2DPi0->Fill(mom_particle->Xv(),mom_particle->Yv());
    }

    if(R<1.){
      if(isAcceptPHOS) fHistGammaAccept->Fill(ene);
      if(isUnitRap)    fHistGammaUnitRap->Fill(ene);
      //pi0
      if(momPID == 111){
	if(isAcceptPHOS) fHistGammaPi0Accept->Fill(ene);
	if(isUnitRap)    fHistGammaPi0UnitRap->Fill(ene);
      }
      //eta
      else if(momPID==221){
	if(isAcceptPHOS) fHistGammaEtaAccept->Fill(ene);
	if(isUnitRap)    fHistGammaEtaUnitRap->Fill(ene);
      }
      //omega
      else if(momPID==223){
	if(isAcceptPHOS) fHistGammaOmegaAccept->Fill(ene);
	if(isUnitRap)    fHistGammaOmegaUnitRap->Fill(ene);
      }
      //eta'
      else if(momPID==331){
	if(isAcceptPHOS) fHistGammaEtaPrimeAccept->Fill(ene);
	if(isUnitRap)    fHistGammaEtaPrimeUnitRap->Fill(ene);
      }
      //phi
      else if(momPID==333){
	if(isAcceptPHOS) fHistGammaPhiAccept->Fill(ene);
	if(isUnitRap)    fHistGammaPhiUnitRap->Fill(ene);
      }
      //rho
      else if(momPID==113){
	if(isAcceptPHOS) fHistGammaRhoAccept->Fill(ene);
	if(isUnitRap)    fHistGammaRhoUnitRap->Fill(ene);
      }
      else{
	if(isAcceptPHOS) fHistGammaOtherAccept->Fill(ene);
	if(isUnitRap)    fHistGammaOtherUnitRap->Fill(ene);
      }

    }
    else{


    }
  
  }//iTrack

}

void AliAnalysisTaskPHOSTrigPi0::SetTriggerInfo(){
  
  if(fIsMC){
    AliCaloTriggerSimulator* TrigMC = new AliCaloTriggerSimulator(fAOD[0]);
    TrigMC->SetThreshold(fTriggerThreshold);
    TrigMC->SetDB(fPHOSGeo,fCalibDataEmc,fPHOSCalibData,fCDBstorage);
    fTrigData= (AliAODCaloTrigger*)TrigMC->CreateTriggerMap();
  }
  else{
    fTrigData = fAOD[0]->GetCaloTrigger("PHOS");
    //cout<<"Real data trigger info is used"<<endl;
  }

}
