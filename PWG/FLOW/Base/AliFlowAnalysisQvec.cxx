/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/*************************************************
 * Save Qvec                                     *
 *                                               *
 * author: Shi Qiu                               *
 *         (s.qiu@nikhef.nl)                     *
 *************************************************/

#define AliFlowAnalysisQvec_cxx

#include "Riostream.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "TChain.h"

#include "TFile.h"
#include "TList.h"
#include "TGraph.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TMath.h"
#include "TArrow.h"
#include "TPaveLabel.h"
#include "TCanvas.h"
#include "AliFlowEventSimple.h"
#include "AliFlowVector.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowAnalysisQvec.h"
#include "AliLog.h"
#include "TRandom.h"
#include "TF1.h"
#include "TNtuple.h"
#include "THnSparse.h"
#include <complex>
#include <cmath>

class TH1;
class TH2;
class TGraph;
class TPave;
class TLatex;
class TMarker;
class TRandom3;
class TObjArray;
class TList;
class TCanvas;
class TSystem;
class TROOT;
class AliFlowVector;
class TVector;

//==============================================================================================================

using std::endl;
using std::cout;
using std::flush;
ClassImp(AliFlowAnalysisQvec)

AliFlowAnalysisQvec::AliFlowAnalysisQvec(const char* name,
                                           Int_t nCen,
                                           Double_t CenWidth):
TNamed(name,name),
// 0.) base:
// 0.1) event list:
fEventList(NULL),
// 1.) common:
// 2a.) particle weights:
fUsePtWeights(kFALSE),
fCutMultiplicityOutliers(kFALSE),
fUseZDCESEMulWeights(kFALSE),
fUseZDCESESpecWeights(kFALSE),
// 2b.) event weights:
fMultiplicityWeight(NULL),
// 3.) integrated flow:
fExactNoRPs(0),
fNumberOfRPsEBE(0.),
fNumberOfPOIsEBE(0.),
fReferenceMultiplicityEBE(0.),
fReferenceMultiplicityRecEBE(0.),
fCentralityEBE(0.),
fZDCESEclEbE(0),
fReQGF(NULL),
fImQGF(NULL),
fNewMetricLEBE(0.),
fNewMetricDEBE(0.),
fNewMetricL2EBE(0.),
fNewMetricD2EBE(0.),
fCentralityCL1EBE(0.),
fNITSCL1EBE(0.),
fCentralityTRKEBE(0.),
// 4.) differential flow:
// 5.) other differential correlators:
// 6.) distributions:
// 8.) debugging and cross-checking:
// 9.) mixed harmonics:
// 10.) Control histograms:
// 11.) Bootstrap:
fRandom(NULL),
// 12.) Charge-Eta Asymmetry:
fCalculateCRC(kTRUE),
fCalculateCME(kFALSE),
fUseZDC(kFALSE),
fRunNum(0),
fCachedRunNum(0),
fRunBin(0),
fCenBin(0),
fDivSigma(kTRUE),
fVariousList(NULL),
fCRCQVecList(NULL),
fZDCESEList(NULL),
fCenWeightEbE(0.),
fCRCVZEROCalibList(NULL),
fStoreZDCQVecVtxPos(kFALSE),
fInvertZDC(kFALSE),
fPOIExtraWeights(kNone),
fSelectCharge(kAllCh),
fUsePhiEtaCuts(kFALSE),
fCRCQVecListTPC(NULL),
fZDCCalibListFinalCommonPart(NULL),
fZDCCalibListFinalRunByRun(NULL),
fCRCZDC2DCutList(NULL),
fCRCZDCCalibList(NULL),
fCRCZDCResList(NULL),
fbFlagIsPosMagField(kFALSE),
fbFlagIsBadRunForC34(kFALSE),
fCRCEtaMin(0.),
fCRCEtaMax(0.),
fCRCnRun(211),
fDataSet(kAny),
fInteractionRate(kAll),
fCRCnCen(nCen),
fZNCQ0(0.),
fZNAQ0(0.),
fZNCen(0.),
fZNAen(0.),
fZPCen(0.),
fZPAen(0.),
fEnNucl(1.),
fQAZDCCuts(kFALSE),
fUseTracklets(kFALSE),
fQAZDCCutsFlag(kTRUE),
fMinMulZN(0)
{
  this->InitializeArraysForVarious();

  // constructor
  // base list to hold all output objects:
  fHistList = new TList();
  fHistList->SetName("cobjQC");
  fHistList->SetOwner(kTRUE);

  // multiplicity weight:
  fMultiplicityWeight = new TString("combinations");
  
  // base list to hold all temp objects:
  fTempList = new TList();
  fTempList->SetName("temp");
  fTempList->SetOwner(kTRUE);
  
  // initialize all arrays:
  this->InitializeArraysForVarious();
  
  fRunList = TArrayI();
  fAvVtxPosX = TArrayD();
  fAvVtxPosY = TArrayD();
  fAvVtxPosZ = TArrayD();
  //@Shi add ave vtx IR split
  fAvVtxPosX15oIRSplit = TArrayD();
  fAvVtxPosY15oIRSplit = TArrayD();
  fAvVtxPosZ15oIRSplit = TArrayD();

  this->InitializeCostantsForCRC();
  this->InitializeArraysForParticleWeights();
  this->InitializeArraysForCRC();
  this->InitializeArraysForCRCVZ();
  this->InitializeArraysForCRCZDC();
  this->InitializeArraysForQVec();
  
  
  for(Int_t i=0; i<fkGFPtB; i++) {
    fReQGFPt[i] = NULL;
    fImQGFPt[i] = NULL;
  }

} // end of constructor

//================================================================================================================

AliFlowAnalysisQvec::~AliFlowAnalysisQvec()
{
  // destructor
  delete fHistList;
  delete fTempList;
  if(fCRCZDCCalibList)    delete fCRCZDCCalibList;
  if(fZDCCalibListFinalCommonPart)       delete fZDCCalibListFinalCommonPart;
  if(fZDCCalibListFinalRunByRun)       delete fZDCCalibListFinalRunByRun;
  if(fCRCZDC2DCutList)    delete fCRCZDC2DCutList;
  if(fCRCZDCResList)      delete fCRCZDCResList;
  if(fCRCVZEROCalibList)  delete fCRCVZEROCalibList;
  if(fZDCESEList)         delete fZDCESEList;
  delete[] fCorrMap;
  
} // end of AliFlowAnalysisQvec::~AliFlowAnalysisQvec()

//================================================================================================================

void AliFlowAnalysisQvec::Init()
{
  // a) Cross check if the settings make sense before starting the QC adventure;
  // b) Access all common constants;
  // c) Book all objects;
  // d) Store flags for integrated and differential flow;
  // e) Store flags for distributions of corelations;
  // f) Store harmonic which will be estimated;
  // g) Store flags for mixed harmonics;
  // h) Store flags for control histograms;
  // i) Store bootstrap flags.

  //save old value and prevent histograms from being added to directory
  //to avoid name clashes in case multiple analaysis objects are used
  //in an analysis
  Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fRandom = new TRandom3(0); // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID

  // a) Cross check if the settings make sense before starting the QC adventure;
  this->CrossCheckSettings();
  // b) Access all common constants and book a profile to hold them:
  this->CommonConstants("Init");
  // c) Book all objects:
  this->BookAndFillWeightsHistograms();
  this->BookAndNestAllLists();
  this->BookEverythingForIntegratedFlow();
  this->SetRunList();
  if(fCalculateCRC) {
    this->BookEverythingForQVec();
    this->BookEverythingForCME();
  }
  this->BookEverythingForVarious();
  this->SetCentralityWeights();

  fpQvecEvent = new AliFlowAnalysisQvecEvent();
  
  fEventList = new TList();
  fEventList->SetName("Event List");
  fEventList->SetOwner(kTRUE);
  fHistList->Add(fEventList);
  
  treeEvent = new TTree("events", "event");
  treeEvent->Branch("event", &fpQvecEvent);
  
  fEventList->Add(treeEvent);
  
  TH1::AddDirectory(oldHistAddStatus);
  
  // printf("Stuff booked \n");

} // end of void AliFlowAnalysisQvec::Init()

//================================================================================================================

void AliFlowAnalysisQvec::Make(AliFlowEventSimple* anEvent)
{
  // Running over data only in this method.
  // a) Check all pointers used in this method;
  // b) Define local variables;
  // c) Fill the common control histograms and call the method to fill fAvMultiplicity;
  // d) Loop over data and calculate e-b-e quantities Q_{n,k}, S_{p,k} and s_{p,k};
  // d.1) initialize particle weights
  // e) Calculate the final expressions for S_{p,k} and s_{p,k} (important !!!!);
  // f) Call the methods which calculate correlations for reference flow;
  // g) Call the methods which calculate correlations for differential flow;
  // h) Call the methods which calculate correlations for 2D differential flow;
  // i) Call the methods which calculate other differential correlators;
  // j) Distributions of correlations;
  // k) Store various;
  // l) Cross-check with nested loops correlators for reference flow;
  // m) Cross-check with nested loops correlators for differential flow;
  // n) Reset all event-by-event quantities (very important !!!!).

  // a) Check all pointers used in this method:
  this->CheckPointersUsedInMake();

  // b) Define local variables:
  Double_t dPhi = 0.; // azimuthal angle in the laboratory frame
  Double_t dPt  = 0.; // transverse momentum
  Double_t dEta = 0.; // pseudorapidity
  Double_t wPhi = 1.; // phi weight
  Double_t wPt  = 1.; // pt weight
  Double_t wEta = 1.; // eta weight
  Double_t wTrack = 1.; // track weight
  Double_t wPhiEta = 1.;
  Double_t wt = 1.;
  Int_t nCounterNoRPs = 0; // needed only for shuffling
  fNumberOfRPsEBE = anEvent->GetNumberOfRPs(); // number of RPs (i.e. number of reference particles)
  fNumberOfPOIsEBE = anEvent->GetNumberOfPOIs(); // number of POIs (i.e. number of particles of interest)
  fReferenceMultiplicityEBE = anEvent->GetReferenceMultiplicity(); // reference multiplicity for current event
  fCentralityEBE = anEvent->GetCentrality(); // centrality percentile for current event
  fCentralityCL1EBE = anEvent->GetCentralityCL1(); // centrality percentile for current event (alternative estimation)
  fCentralityTRKEBE = anEvent->GetCentralityTRK(); // centrality percentile for current event (alternative estimation)
  fNITSCL1EBE = anEvent->GetNITSCL1();

//  printf("begin AliFlowAnalysisQvec::Make \n");
  
  //if(fExactNoRPs > 0 && fNumberOfRPsEBE<fExactNoRPs){return;} // no need to have fExactNoRPs
  if(!fCentralityEBE){return;}
  if(!fNumberOfRPsEBE || !fNumberOfPOIsEBE){return;}
  fhCenvsMul[0]->Fill(fCentralityEBE,fReferenceMultiplicityEBE); // some QA hists
  if(fCutMultiplicityOutliers) {
    if((fDataSet==k2015 || fDataSet==k2015v6 || fDataSet==k2015pidfix) && !MultCut2015o()){return;}
  }

  if(!fRefMultRbRPro) { // not useful for Q vector
    fReferenceMultiplicityRecEBE = fReferenceMultiplicityEBE-fMultCutAv->GetBinContent(fMultCutAv->FindBin(fCentralityEBE));
  } else {
    Int_t runbin = fRefMultRbRPro->GetXaxis()->FindBin(Form("%d",fRunNum));
    Int_t cenbin = fRefMultRbRPro->GetYaxis()->FindBin(fCentralityEBE);
    fReferenceMultiplicityRecEBE = fReferenceMultiplicityEBE-fRefMultRbRPro->GetBinContent(runbin,cenbin);
  }

  // centrality flattening with weights
  fCenWeightEbE = 1.;
  if(fCenWeigCalHist) fCenWeightEbE = fCenWeigCalHist->GetBinContent(fCenWeigCalHist->FindBin(fCentralityEBE));

  // primary vertex position (x,y,z)
  anEvent->GetVertexPosition(fVtxPos);
  // re-centered around zer (implemented only for run2)

  if(fDataSet==k2015 || fDataSet==k2015v6 || fDataSet==k2015pidfix) {
    if(fAvVtxPosX[fRunBin]) fVtxPosCor[0] = fVtxPos[0]-fAvVtxPosX[fRunBin];
    if(fAvVtxPosY[fRunBin]) fVtxPosCor[1] = fVtxPos[1]-fAvVtxPosY[fRunBin];
    if(fAvVtxPosZ[fRunBin]) fVtxPosCor[2] = fVtxPos[2]-fAvVtxPosZ[fRunBin];
    //cout<<"===> fRunBin === "<<fRunBin<<endl;
    //cout<<"===> fVtxPosCor15oIRSplit[0] = "<<fVtxPosCor15oIRSplit[0]<<" = fVtxPos[0]-fAvVtxPosX15oIRSplit[fRunBin] = "<<fVtxPos[0]<<" - "<<fAvVtxPosX15oIRSplit[fRunBin]<<endl;
    //cout<<"===> fVtxPosCor15oIRSplit[1] = "<<fVtxPosCor15oIRSplit[1]<<" = fVtxPos[1]-fAvVtxPosY15oIRSplit[fRunBin] = "<<fVtxPos[1]<<" - "<<fAvVtxPosY15oIRSplit[fRunBin]<<endl;
    //cout<<"===> fVtxPosCor15oIRSplit[2] = "<<fVtxPosCor15oIRSplit[2]<<" = fVtxPos[2]-fAvVtxPosY15oIRSplit[fRunBin] = "<<fVtxPos[2]<<" - "<<fAvVtxPosZ15oIRSplit[fRunBin]<<endl;
    if(fAvVtxPosX15oIRSplit[fRunBin]) fVtxPosCor15oIRSplit[0] = fVtxPos[0]-fAvVtxPosX15oIRSplit[fRunBin];
    if(fAvVtxPosY15oIRSplit[fRunBin]) fVtxPosCor15oIRSplit[1] = fVtxPos[1]-fAvVtxPosY15oIRSplit[fRunBin];
    if(fAvVtxPosZ15oIRSplit[fRunBin]) fVtxPosCor15oIRSplit[2] = fVtxPos[2]-fAvVtxPosZ15oIRSplit[fRunBin];
  } else {
    fVtxPosCor[0] = fVtxPos[0];
    fVtxPosCor[1] = fVtxPos[1];
    fVtxPosCor[2] = fVtxPos[2];
    
    fVtxPosCor15oIRSplit[0] = fVtxPos[0];
    fVtxPosCor15oIRSplit[1] = fVtxPos[1];
    fVtxPosCor15oIRSplit[2] = fVtxPos[2];
  }

  Double_t ptEta[2] = {0.,0.}; // 0 = dPt, 1 = dEta
  Int_t dCharge = 0; // charge

  // d) Loop over data and calculate e-b-e quantities Q_{n,k}, S_{p,k} and s_{p,k}:
  Int_t nPrim = anEvent->NumberOfTracks();  // nPrim = total number of primary tracks
  AliFlowTrackSimple *aftsTrack = NULL;

  // d.1) Initialize particle weights
  Int_t cw = 0;

  if (fDataSet==kAny) {fRunBin = 0;}
  else {fRunBin = GetCRCRunBin(fRunNum);}
  if(fRunBin<0 || fRunBin>=fCRCnRun) {return;}
  fCenBin = GetCRCCenBin(fCentralityEBE);
  if(fCenBin<0 || fCenBin>=fCRCnCen) {return;}

  fhCenvsMul[2]->Fill(fNumberOfRPsEBE,fNumberOfPOIsEBE); //defined
  fhCenvsMul[3]->Fill(fReferenceMultiplicityEBE,fNumberOfRPsEBE);
  fhCenvsMul[4]->Fill(fReferenceMultiplicityEBE,fNumberOfPOIsEBE);
  if(fhAvAbsOrbit) {
    UInt_t TimeStamp = (UInt_t)anEvent->GetAbsOrbit();
    fhAvAbsOrbit->Fill(fRunBin+0.5,(Double_t)TimeStamp); // defined 
  }

  if(fDataSet==k2015 || fDataSet==k2015v6) {
    if(fRunNum!=fCachedRunNum) { // fRunNum and fCachedRunNum defined
      fbFlagIsPosMagField = kFALSE;
      Int_t dRun15hPos[] = {246390, 246391, 246392, 246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424};
      for (Int_t i=0; i<40; i++) {
        if(fRunNum==dRun15hPos[i]) fbFlagIsPosMagField = kTRUE;
      }
      fbFlagIsBadRunForC34 = kFALSE;
      Int_t BadRunList[] = {245705, 246042, 246049, 246087, 246151, 246181, 246217, 246222, 246272, 246275};
      for (Int_t i=0; i<10; i++) {
        if(fRunNum==BadRunList[i]) fbFlagIsBadRunForC34 = kTRUE;
      }
    }
  }

  // VZERO *********************************************************************************************************

  if(fUseVZERO) {
    //cout<<"========> fCRCnHar ====== "<<fCRCnHar<<endl;
    for(Int_t h=0; h<fCRCnHar; h++) {
		//cout<<"========> fCRCnHar hhhhh ====== "<<h<<endl;
      // Get Q vectors for the subevents
      AliFlowVector vQarray[2];
      anEvent->GetV02Qsub(vQarray,h+1);
      fVZFlowVect[0][h] = vQarray[0]; // defined
      fVZFlowVect[1][h] = vQarray[1];

	  if (h==1) {
		fpQvecEvent->setVZCRe(fVZFlowVect[0][h].X());
		fpQvecEvent->setVZCIm(fVZFlowVect[0][h].Y());
		fpQvecEvent->setVZCM(fVZFlowVect[0][h].GetMult());
		fpQvecEvent->setVZARe(fVZFlowVect[1][h].X());
		fpQvecEvent->setVZAIm(fVZFlowVect[1][h].Y());
		fpQvecEvent->setVZAM(fVZFlowVect[1][h].GetMult()); 
	  }
	  
      // re-center VZERO Q-vectors
      if(fCRCVZEROCalibList) this->RecenterCRCQVecVZERO();

      // fill Q-vector RbR
//      if(fCRCVZEROQVec[fRunBin][h]) {
//        fCRCVZEROQVec[fRunBin][h]->Fill(0.5,fCentralityEBE,fVZFlowVect[0][h].X());
//        fCRCVZEROQVec[fRunBin][h]->Fill(1.5,fCentralityEBE,fVZFlowVect[0][h].Y());
//        fCRCVZEROQVec[fRunBin][h]->Fill(2.5,fCentralityEBE,fVZFlowVect[1][h].X());
//        fCRCVZEROQVec[fRunBin][h]->Fill(3.5,fCentralityEBE,fVZFlowVect[1][h].Y());
//        fCRCVZEROQVec[fRunBin][h]->Fill(4.5,fCentralityEBE,fVZFlowVect[0][h].X()*fVZFlowVect[1][h].X());
//        fCRCVZEROQVec[fRunBin][h]->Fill(5.5,fCentralityEBE,fVZFlowVect[0][h].Y()*fVZFlowVect[1][h].Y());
//        fCRCVZEROQVec[fRunBin][h]->Fill(6.5,fCentralityEBE,fVZFlowVect[0][h].X()*fVZFlowVect[1][h].Y());
//        fCRCVZEROQVec[fRunBin][h]->Fill(7.5,fCentralityEBE,fVZFlowVect[0][h].Y()*fVZFlowVect[1][h].X());
//      }
    } // end of for(Int_t h=0; h<fCRCnHar; h++)
  } // end of if(fUseVZERO)

  // ZDC *********************************************************************************************************

  if(fUseZDC) {
    // Get Q vectors for the subevents
    AliFlowVector vQarray[2];
    anEvent->GetZDC2Qsub(vQarray);
    fZDCFlowVect[0] = vQarray[0];
    fZDCFlowVect[1] = vQarray[1];
    fZNCQ0 = anEvent->GetZNCQ0()/fEnNucl;
    fZNAQ0 = anEvent->GetZNAQ0()/fEnNucl;
    fZNCen = anEvent->GetZNCEnergy()/fEnNucl;
    fZNAen = anEvent->GetZNAEnergy()/fEnNucl;
    fZPCen = anEvent->GetZPCEnergy()/fEnNucl;
    fZPAen = anEvent->GetZPAEnergy()/fEnNucl;
  } // end of if(fUseZDC)
  this->CalculateCRCQVec();

  fpQvecEvent->setZCRe(fZDCFlowVect[0].X());
  fpQvecEvent->setZCIm(fZDCFlowVect[0].Y());
  fpQvecEvent->setZCM(fZDCFlowVect[0].GetMult());
  fpQvecEvent->setZARe(fZDCFlowVect[1].X());
  fpQvecEvent->setZAIm(fZDCFlowVect[1].Y());
  fpQvecEvent->setZAM(fZDCFlowVect[1].GetMult());
  
  if(fUseZDC) {
    this->PassQAZDCCuts();
    if(fRecenterZDC) {
      //this->RecenterCRCQVecZDC();
      this->RecenterCRCQVecZDC2();

    }
  }
  
  // ZDC-C (eta < -8.8)
  Double_t ZCRe = fZDCFlowVect[0].X();
  Double_t ZCIm = fZDCFlowVect[0].Y();
  Double_t ZCM  = fZDCFlowVect[0].GetMult();
  // ZDC-A (eta > 8.8)
  Double_t ZARe = fZDCFlowVect[1].X();
  Double_t ZAIm = fZDCFlowVect[1].Y();
  Double_t ZAM  = fZDCFlowVect[1].GetMult();
  if( fInvertZDC ) ZARe = -ZARe;

  ////////////////////////////////////////////////////////////////// test ///////////////////////////////////////////////
  // VZ eta < 0
  Double_t VZCRe = fVZFlowVect[0][1].X();
  Double_t VZCIm = fVZFlowVect[0][1].Y();
  Double_t VZCM  = fVZFlowVect[0][1].GetMult();
  Double_t EvPlVZC = TMath::ATan2(VZCIm,VZCRe)/2.;
  // VZ eta > 0
  Double_t VZARe = fVZFlowVect[1][1].X();
  Double_t VZAIm = fVZFlowVect[1][1].Y();
  Double_t VZAM  = fVZFlowVect[1][1].GetMult();
  Double_t EvPlVZA = TMath::ATan2(VZAIm,VZARe)/2.;
  ////////////////////////////////////////////////////////////////// test ///////////////////////////////////////////////
  
  // ZDC QA cuts
  Bool_t bPassZDCcuts = kTRUE;
  if( ZCM<=0. || ZAM<=0. || sqrt(ZCRe*ZCRe+ZCIm*ZCIm)<1.E-6 || sqrt(ZARe*ZARe+ZAIm*ZAIm)<1.E-6 ) bPassZDCcuts=kFALSE;
  if( !std::isfinite(fZDCFlowVect[0].Mod()) || !std::isfinite(fZDCFlowVect[1].Mod())) bPassZDCcuts=kFALSE;
  if(fQAZDCCuts && !fQAZDCCutsFlag) bPassZDCcuts=kFALSE;
  if(bPassZDCcuts) fEventCounter->Fill(1.5);

  // EbE flow *********************************************************************************************************

  // run-by-run corrections ********************************************************************************************

  if(fRunNum!=fCachedRunNum) {
    for (Int_t cb=0; cb<fCRCnCen; cb++) {
      if(fPtWeightsHist[cb]){
        for (Int_t bx=1; bx<=fPtWeightsHist[cb]->GetNbinsX(); bx++) {
          fPtWeightsCent->SetBinContent(bx,cb+1,fPtWeightsHist[cb]->GetBinContent(bx));
        }
      }
    }

    if(fPOIExtraWeights==AliFlowAnalysisQvec::kEtaPhiRbR) {
      if(fWeightsList->FindObject(Form("fCRCQVecPhiHistRbR[%d]",fRunNum))) {
        fPhiEtaRbRWeights = (TH3D*)(fWeightsList->FindObject(Form("fCRCQVecPhiRbRHist[%d]",fRunNum)));
      } else {
        AliWarning("WARNING: POIExtraWeights (kEtaPhiRbR) not found ! \n");
      }
    }
    if(fPOIExtraWeights==AliFlowAnalysisQvec::kEtaPhiChRbR) {
      for (Int_t i=0; i<2; i++) {
        if(fWeightsList->FindObject(Form("fCRCQVecPhiRbRHistCh[%d][%d]",fRunNum,i))) {
          fPhiEtaRbRWeightsCh[i] = (TH3D*)(fWeightsList->FindObject(Form("fCRCQVecPhiRbRHistCh[%d][%d]",fRunNum,i)));
        } else {
          AliWarning("WARNING: POIExtraWeights (kEtaPhiChRbR) not found ! \n");
        }
      }
    }
    if(fPOIExtraWeights==AliFlowAnalysisQvec::kEtaPhiVtxRbR) {
      for (Int_t cb=0; cb<fCRCnCen; cb++) {
        if(fWeightsList->FindObject(Form("CRCQVecPhiHistVtx[%d][%d]",cb,fRunNum))) {
          fPhiEtaWeightsVtx[cb] = (TH3D*)(fWeightsList->FindObject(Form("CRCQVecPhiHistVtx[%d][%d]",cb,fRunNum)));
        } else {
          AliWarning("WARNING: POIExtraWeights (kEtaPhiVtxRbR) not found ! \n");
        }
      }
    }
  }

  // loop over particles **********************************************************************************************
  for(Int_t i=0;i<nPrim;i++) {
    if(fExactNoRPs > 0 && nCounterNoRPs>fExactNoRPs){continue;}
    aftsTrack=anEvent->GetTrack(i);
    if(aftsTrack) {

      if(!(aftsTrack->InRPSelection() || aftsTrack->InPOISelection() || aftsTrack->InPOISelection(2))){continue;} // safety measure: consider only tracks which are RPs or POIs

      // RPs *********************************************************************************************************

      /*if(aftsTrack->InRPSelection()) {
        nCounterNoRPs++;
        dPhi = aftsTrack->Phi();
        dPt  = aftsTrack->Pt();
        dEta = aftsTrack->Eta();
        dCharge = aftsTrack->Charge();

        if(fSelectCharge==kPosCh && dCharge<0.) continue;
        if(fSelectCharge==kNegCh && dCharge>0.) continue;

        cw = (dCharge > 0. ? 0 : 1);
        wPhi = 1.;
        wPt  = 1.;
        wEta = 1.;
        wTrack = 1.;
        wPhiEta = 1.;

        // pT weights
        if(fUsePtWeights && fPtWeightsHist[fCenBin]) {
          if(dPt>fPtWeightsCent->GetXaxis()->GetXmin() && dPt<fPtWeightsCent->GetXaxis()->GetXmax()) wt = fPtWeightsCent->Interpolate(dPt,fCentralityEBE);
          else if(dPt<fPtWeightsCent->GetXaxis()->GetXmin())  wt = fPtWeightsCent->Interpolate(fPtWeightsCent->GetXaxis()->GetXmin(),fCentralityEBE);
          else if(dPt>fPtWeightsCent->GetXaxis()->GetXmax())  wt = fPtWeightsCent->Interpolate(fPtWeightsCent->GetXaxis()->GetXmax(),fCentralityEBE);
          if(std::isfinite(1./wt)) wPhiEta *= 1./wt;
        }

        // extra weights: eta, phi, ch, vtx
        if(fPOIExtraWeights==kEtaPhi && fPhiEtaWeights) // determine phieta weight for POI:
        {
          wt = fPhiEtaWeights->GetBinContent(fPhiEtaWeights->FindBin(fCentralityEBE,dPhi,dEta));
          if(std::isfinite(1./wt)) wPhiEta *= 1./wt;
        }
        if(fPOIExtraWeights==kEtaPhiCh && fPhiEtaWeightsCh[cw]) // determine phieta weight for POI, ch dep:
        {
          wt = fPhiEtaWeightsCh[cw]->GetBinContent(fPhiEtaWeightsCh[cw]->FindBin(fCentralityEBE,dPhi,dEta));
          if(std::isfinite(1./wt)) wPhiEta *= 1./wt;
        }
        if((fPOIExtraWeights==kEtaPhiVtx || fPOIExtraWeights==kEtaPhiVtxRbR) && fPhiEtaWeightsVtx[fCenBin]) // determine phieta weight for POI:
        {
          wt = fPhiEtaWeightsVtx[fCenBin]->GetBinContent(fPhiEtaWeightsVtx[fCenBin]->FindBin(fVtxPosCor[2],dPhi,dEta));
          if(wt==0.) continue;
          if(std::isfinite(1./wt)) wPhiEta *= 1./wt;
        }
        Int_t ptbebe = (dPt>1.? 2 : (dPt>0.5 ? 1 : 0)); // hardcoded
        if(fPOIExtraWeights==kEtaPhiChPt && fPhiEtaWeightsChPt[cw][ptbebe]) // determine phieta weight for POI, ch dep:
        {
          wt = fPhiEtaWeightsChPt[cw][ptbebe]->GetBinContent(fPhiEtaWeightsChPt[cw][ptbebe]->FindBin(fCentralityEBE,dPhi,dEta));
          if(std::isfinite(1./wt)) wPhiEta *= 1./wt;
        }
        // run-by-run
        if(fPOIExtraWeights==kEtaPhiRbR && fPhiEtaRbRWeights) // determine phieta weight for POI:
        {
          wt = fPhiEtaRbRWeights->GetBinContent(fPhiEtaRbRWeights->FindBin(fCentralityEBE,dPhi,dEta));
          if(std::isfinite(1./wt)) wPhiEta *= 1./wt;
        }
        if(fPOIExtraWeights==kEtaPhiChRbR && fPhiEtaRbRWeightsCh[cw]) // determine phieta weight for POI, ch dep:
        {
          wt = fPhiEtaRbRWeightsCh[cw]->GetBinContent(fPhiEtaRbRWeightsCh[cw]->FindBin(fCentralityEBE,dPhi,dEta));
          if(std::isfinite(1./wt)) wPhiEta *= 1./wt;
        }

        if(fPhiExclZoneHist) {
          if(fPhiExclZoneHist->GetBinContent(fPhiExclZoneHist->FindBin(dEta,dPhi))<0.5) continue;
        }

        // Calculate Re[Q_{m*n,k}] and Im[Q_{m*n,k}] for this event (m = 1,2,...,12, k = 0,1,...,8):
        for(Int_t m=0;m<12;m++) // to be improved - hardwired 6
        {
          for(Int_t k=0;k<9;k++) // to be improved - hardwired 9
          {
            (*fReQ)(m,k)+=pow(wPhiEta*wPhi*wPt*wEta*wTrack,k)*TMath::Cos((m+1)*n*dPhi);
            (*fImQ)(m,k)+=pow(wPhiEta*wPhi*wPt*wEta*wTrack,k)*TMath::Sin((m+1)*n*dPhi);
          }
        }
        // Calculate S_{p,k} for this event (Remark: final calculation of S_{p,k} follows after the loop over data bellow):
        for(Int_t p=0;p<8;p++)
        {
          for(Int_t k=0;k<9;k++)
          {
            (*fSpk)(p,k)+=pow(wPhiEta*wPhi*wPt*wEta*wTrack,k);
          }
        }

      } // end of if(pTrack->InRPSelection())*/

      // POIs ********************************************************************************************************

      if(aftsTrack->InPOISelection() || aftsTrack->InPOISelection(2)) {

        if(!fUseTracklets && aftsTrack->InPOISelection(2)) continue;
        if(fUseTracklets && !aftsTrack->InPOISelection(2)) continue;

        dPhi = aftsTrack->Phi();
        dPt  = aftsTrack->Pt();
        dEta = aftsTrack->Eta();
        dCharge = aftsTrack->Charge();
        // Int_t ITStype = aftsTrack->ITStype();

        if(fSelectCharge==kPosCh && dCharge<0.) continue;
        if(fSelectCharge==kNegCh && dCharge>0.) continue;

        Bool_t IsSplitMergedTracks = kFALSE;
        if(fRemoveSplitMergedTracks) {
          IsSplitMergedTracks = EvaulateIfSplitMergedTracks(anEvent,aftsTrack,i);
        }
        if(IsSplitMergedTracks) continue;

        if(fCRCTestSin) {
          if(dCharge > 0.) dPt += 1.E-2;
          else dPt -= 1.E-2;
        }

        cw = (dCharge > 0. ? 0 : 1);
        wPhi = 1.;
        wPt  = 1.;
        wEta = 1.;
        wTrack = 1.;
        wPhiEta = 1.;

        if(fMinMulZN==0 && dPhi>3.141593e-01 && dPhi<1.256637e+00) {
          if(dEta>0.) continue;
          if(dEta<0.) wPhiEta *= 2.;
        }

        // pT weights
        if(fUsePtWeights && fPtWeightsHist[fCenBin]) {
          if(dPt>fPtWeightsCent->GetXaxis()->GetXmin() && dPt<fPtWeightsCent->GetXaxis()->GetXmax()) wt = fPtWeightsCent->Interpolate(dPt,fCentralityEBE);
          else if(dPt<fPtWeightsCent->GetXaxis()->GetXmin())  wt = fPtWeightsCent->Interpolate(fPtWeightsCent->GetXaxis()->GetXmin(),fCentralityEBE);
          else if(dPt>fPtWeightsCent->GetXaxis()->GetXmax())  wt = fPtWeightsCent->Interpolate(fPtWeightsCent->GetXaxis()->GetXmax(),fCentralityEBE);
          if(std::isfinite(1./wt)) wPhiEta *= 1./wt;
        }

        // extra weights: eta, phi, ch, vtx
        if(fPOIExtraWeights==kEtaPhi && fPhiEtaWeights) // determine phieta weight for POI:
        {
          wt = fPhiEtaWeights->GetBinContent(fPhiEtaWeights->FindBin(fCentralityEBE,dPhi,dEta));
          if(std::isfinite(1./wt)) wPhiEta *= 1./wt;
        }
        if(fPOIExtraWeights==kEtaPhiCh && fPhiEtaWeightsCh[cw]) // determine phieta weight for POI, ch dep:
        {
          wt = fPhiEtaWeightsCh[cw]->GetBinContent(fPhiEtaWeightsCh[cw]->FindBin(fCentralityEBE,dPhi,dEta));
          if(std::isfinite(1./wt)) wPhiEta *= 1./wt;
        }
        if((fPOIExtraWeights==kEtaPhiVtx || fPOIExtraWeights==kEtaPhiVtxRbR) && fPhiEtaWeightsVtx[fCenBin]) // determine phieta weight for POI:
        {
          wt = fPhiEtaWeightsVtx[fCenBin]->GetBinContent(fPhiEtaWeightsVtx[fCenBin]->FindBin(fVtxPosCor[2],dPhi,dEta));
          if(wt==0.) continue;
          if(std::isfinite(1./wt)) wPhiEta *= 1./wt;
        }
        Int_t ptbebe = (dPt>1.? 2 : (dPt>0.5 ? 1 : 0)); // hardcoded
        if(fPOIExtraWeights==kEtaPhiChPt && fPhiEtaWeightsChPt[cw][ptbebe]) // determine phieta weight for POI, ch dep:
        {
          wt = fPhiEtaWeightsChPt[cw][ptbebe]->GetBinContent(fPhiEtaWeightsChPt[cw][ptbebe]->FindBin(fCentralityEBE,dPhi,dEta));
          if(std::isfinite(1./wt)) wPhiEta *= 1./wt;
        }
        // run-by-run
        if(fPOIExtraWeights==kEtaPhiRbR && fPhiEtaRbRWeights) // determine phieta weight for POI:
        {
          wt = fPhiEtaRbRWeights->GetBinContent(fPhiEtaRbRWeights->FindBin(fCentralityEBE,dPhi,dEta));
          if(std::isfinite(1./wt)) wPhiEta *= 1./wt;
        }
        if(fPOIExtraWeights==kEtaPhiChRbR && fPhiEtaRbRWeightsCh[cw]) // determine phieta weight for POI, ch dep:
        {
          wt = fPhiEtaRbRWeightsCh[cw]->GetBinContent(fPhiEtaRbRWeightsCh[cw]->FindBin(fCentralityEBE,dPhi,dEta));
          if(std::isfinite(1./wt)) wPhiEta *= 1./wt;
        }

        if(fUsePhiEtaCuts)
        {
          // test: remove region with low SPD efficiency
          if(dPhi>2.136283 && dPhi<2.324779) continue;
        }

        // Generic Framework: Calculate Re[Q_{m*n,k}] and Im[Q_{m*n,k}] for this event (m = 1,2,...,12, k = 0,1,...,8):
        Double_t MaxPtCut = 3.;
        if(fMinMulZN==99) MaxPtCut = 1.;
        if(dPt<MaxPtCut) {
          for(Int_t m=0;m<21;m++) // to be improved - hardwired 6
          {
            for(Int_t k=0;k<9;k++) // to be improved - hardwired 9
            {
              (*fReQGF)(m,k) += pow(wPhiEta*wPhi*wPt*wEta*wTrack,k)*TMath::Cos(m*dPhi);
              (*fImQGF)(m,k) += pow(wPhiEta*wPhi*wPt*wEta*wTrack,k)*TMath::Sin(m*dPhi);
            }
          }
        }

        for(Int_t ptb=0; ptb<fkGFPtB; ptb++) {
          if(ptb==0 && dPt>0.5) continue;
          if(ptb==1 && (dPt<0.5 || dPt>1.)) continue;
          if(ptb==2 && (dPt<1. || dPt>2.)) continue;
          if(ptb==3 && dPt<2.) continue;
          if(ptb==4 && (dPt<1. || dPt>2.5)) continue;
          if(ptb==5 && dPt<2.5) continue;
          if(ptb==6 && (dPt<1. || dPt>3.)) continue;
          if(ptb==7 && dPt<3.) continue;
          for(Int_t m=0;m<21;m++) // to be improved - hardwired 6
          {
            for(Int_t k=0;k<9;k++) // to be improved - hardwired 9
            {
              (*fReQGFPt[ptb])(m,k) += pow(wPhiEta*wPhi*wPt*wEta*wTrack,k)*TMath::Cos(m*dPhi);
              (*fImQGFPt[ptb])(m,k) += pow(wPhiEta*wPhi*wPt*wEta*wTrack,k)*TMath::Sin(m*dPhi);
            }
          }
        }

        // Charge-Rapidity Correlations
        for (Int_t h=0;h<fCRCnHar;h++) {
            
          if(fCalculateCME) {
            Double_t SpecWeig = 1.;
            if(fUseZDCESESpecWeights && fZDCESESpecWeightsHist[fZDCESEclEbE] && dPt>0.2 && dPt<20.2) {
              Double_t weraw = fZDCESESpecWeightsHist[fZDCESEclEbE]->GetBinContent(fZDCESESpecWeightsHist[fZDCESEclEbE]->FindBin(fCentralityEBE,dPt));
              if(weraw > 0.) SpecWeig = 1./weraw;
            }
            
            fCMEQRe[cw][h]->Fill(dEta,SpecWeig*wPhiEta*TMath::Cos((h+1.)*dPhi));
            fCMEQIm[cw][h]->Fill(dEta,SpecWeig*wPhiEta*TMath::Sin((h+1.)*dPhi));
            fCMEMult[cw][h]->Fill(dEta,SpecWeig*wPhiEta);
            fCMEQRe[2+cw][h]->Fill(dEta,pow(SpecWeig*wPhiEta,2.)*TMath::Cos((h+1.)*dPhi));
            fCMEQIm[2+cw][h]->Fill(dEta,pow(SpecWeig*wPhiEta,2.)*TMath::Sin((h+1.)*dPhi));
            fCMEMult[2+cw][h]->Fill(dEta,pow(SpecWeig*wPhiEta,2.));
            
            //@shi add histogram for both charges
            fCMEQReBothCharge[0][h]->Fill(dEta,SpecWeig*wPhiEta*TMath::Cos((h+1.)*dPhi));
            fCMEQImBothCharge[0][h]->Fill(dEta,SpecWeig*wPhiEta*TMath::Sin((h+1.)*dPhi));
            fCMEMultBothCharge[0][h]->Fill(dEta,SpecWeig*wPhiEta);
            fCMEQReBothCharge[1][h]->Fill(dEta,pow(SpecWeig*wPhiEta,2.)*TMath::Cos((h+1.)*dPhi));
            fCMEQImBothCharge[1][h]->Fill(dEta,pow(SpecWeig*wPhiEta,2.)*TMath::Sin((h+1.)*dPhi));
            fCMEMultBothCharge[1][h]->Fill(dEta,pow(SpecWeig*wPhiEta,2.));
			
          }

        } // end of for (Int_t h=0;h<fCRCnHar;h++)


      } // end of if(pTrack->InPOISelection())
    } else // to if(aftsTrack)
    {
      printf("\n WARNING (QC): No particle (i.e. aftsTrack is a NULL pointer in AFAWQC::Make())!!!!\n\n");
    }
  } // end of for(Int_t i=0;i<nPrim;i++)

  // ************************************************************************************************************
  
  //QvecEvent = new AliFlowAnalysisQvecEvent();
  fpQvecEvent->setRunNum(fRunNum);
  fpQvecEvent->setCentrality(fCentralityEBE);
  fpQvecEvent->setVtxPosX(fVtxPos[0]);
  fpQvecEvent->setVtxPosY(fVtxPos[1]);
  fpQvecEvent->setVtxPosZ(fVtxPos[2]);

  
  this->CalculateCMESPPP();
  
  treeEvent->Fill();
  
  // o) Reset all event-by-event quantities (very important !!!!):
  this->ResetEventByEventQuantities();
  
//  printf("end AliFlowAnalysisQvec::Make \n");

} // end of AliFlowAnalysisQvec::Make(AliFlowEventSimple* anEvent)

//=======================================================================================================================

void AliFlowAnalysisQvec::CalculateCMESPPP()
{
  //************************************************ Weights **************************************
  Double_t MulWeig = 1.;
  //@shi no file supplied
  if(fUseZDCESEMulWeights && fZDCESEMultWeightsHist[fZDCESEclEbE]) {
    Double_t weraw = fZDCESEMultWeightsHist[fZDCESEclEbE]->GetBinContent(fZDCESEMultWeightsHist[fZDCESEclEbE]->FindBin(fCentralityEBE,fNumberOfPOIsEBE));
    if(weraw > 0.) MulWeig = 1./weraw;
  }
  
  //************************************************ Get all variables ****************************
  //*********************************************** TPC part *************************************
  Int_t h = 0; //@Shi used for TPC and v0 part. For ZDCpart, it is set to 1
  Double_t e = 1E-5;
  
  Double_t uPNReTPC=0., uPNImTPC=0., uPN2ReTPC=0., uPN2ImTPC=0., uPN2Re2TPC=0., uPN2Im2TPC=0., uPNMTPC=0., uPN2MTPC=0.;
  Double_t uPNReTPCPosEta=0., uPNImTPCPosEta=0., uPN2ReTPCPosEta=0., uPN2ImTPCPosEta=0., uPN2Re2TPCPosEta=0., uPN2Im2TPCPosEta=0., uPNMTPCPosEta=0., uPN2MTPCPosEta=0.;
  Double_t uPNReTPCNegEta=0., uPNImTPCNegEta=0., uPN2ReTPCNegEta=0., uPN2ImTPCNegEta=0., uPN2Re2TPCNegEta=0., uPN2Im2TPCNegEta=0., uPNMTPCNegEta=0., uPN2MTPCNegEta=0.;
  
  Double_t uPReTPCPosEta=0., uPImTPCPosEta=0., uP2ReTPCPosEta=0., uP2ImTPCPosEta=0., uP2Re2TPCPosEta=0., uP2Im2TPCPosEta=0., uPMTPCPosEta=0., uP2MTPCPosEta=0.;
  Double_t uNReTPCPosEta=0., uNImTPCPosEta=0., uN2ReTPCPosEta=0., uN2ImTPCPosEta=0., uN2Re2TPCPosEta=0., uN2Im2TPCPosEta=0., uNMTPCPosEta=0., uN2MTPCPosEta=0.;
  
  Double_t uPReTPCNegEta=0., uPImTPCNegEta=0., uP2ReTPCNegEta=0., uP2ImTPCNegEta=0., uP2Re2TPCNegEta=0., uP2Im2TPCNegEta=0., uPMTPCNegEta=0., uP2MTPCNegEta=0.;
  Double_t uNReTPCNegEta=0., uNImTPCNegEta=0., uN2ReTPCNegEta=0., uN2ImTPCNegEta=0., uN2Re2TPCNegEta=0., uN2Im2TPCNegEta=0., uNMTPCNegEta=0., uN2MTPCNegEta=0.;

  for(Int_t EBin=1; EBin<=fCMEQRe[0][0]->GetNbinsX(); EBin++) {
    // both charge all region [0][h]: first index is the power of weight, h is cos((h+1)phi)
    uPNReTPC += fCMEQReBothCharge[0][h]->GetBinContent(EBin); // SpecWeig*wPhiEta*TMath::Cos(dPhi)
    uPNImTPC += fCMEQImBothCharge[0][h]->GetBinContent(EBin);
    uPN2ReTPC += fCMEQReBothCharge[0][h+1]->GetBinContent(EBin); // SpecWeig*wPhiEta*TMath::Cos(2dPhi)
    uPN2ImTPC += fCMEQImBothCharge[0][h+1]->GetBinContent(EBin);
    uPN2Re2TPC += fCMEQReBothCharge[1][h+1]->GetBinContent(EBin); // pow(SpecWeig*wPhiEta,2)*TMath::Cos(2dPhi)
    uPN2Im2TPC += fCMEQImBothCharge[1][h+1]->GetBinContent(EBin);
    uPNMTPC += fCMEMultBothCharge[0][h]->GetBinContent(EBin);
    uPN2MTPC += fCMEMultBothCharge[1][h]->GetBinContent(EBin);
    
    if (EBin == fCMEQRe[0][h]->FindBin(0.4)) { // positive eta region
		// both charge pos eta region
		uPNReTPCPosEta += fCMEQReBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCPosEta += fCMEQImBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCPosEta += fCMEQReBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCPosEta += fCMEQImBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCPosEta += fCMEQReBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCPosEta += fCMEQImBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCPosEta += fCMEMultBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCPosEta += fCMEMultBothCharge[1][h]->GetBinContent(EBin);
		
		// positive charge positive eta region
		uPReTPCPosEta += fCMEQRe[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCPosEta += fCMEQIm[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCPosEta += fCMEQRe[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCPosEta += fCMEQIm[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCPosEta += fCMEQRe[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCPosEta += fCMEQIm[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCPosEta += fCMEMult[0][h]->GetBinContent(EBin); // w
		uP2MTPCPosEta += fCMEMult[2][h]->GetBinContent(EBin); // w^2
		
		// negative charge positive eta region
		uNReTPCPosEta += fCMEQRe[1][h]->GetBinContent(EBin);
		uNImTPCPosEta += fCMEQIm[1][h]->GetBinContent(EBin);
		uN2ReTPCPosEta += fCMEQRe[1][h+1]->GetBinContent(EBin);
		uN2ImTPCPosEta += fCMEQIm[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCPosEta += fCMEQRe[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCPosEta += fCMEQIm[3][h+1]->GetBinContent(EBin);
		uNMTPCPosEta += fCMEMult[1][h]->GetBinContent(EBin);
		uN2MTPCPosEta += fCMEMult[3][h]->GetBinContent(EBin);
	} else if (EBin == fCMEQRe[0][h]->FindBin(-0.4)) { // negative eta region
		// both charge neg eta region
		uPNReTPCNegEta += fCMEQReBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCNegEta += fCMEQImBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCNegEta += fCMEQReBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCNegEta += fCMEQImBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCNegEta += fCMEQReBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCNegEta += fCMEQImBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCNegEta += fCMEMultBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCNegEta += fCMEMultBothCharge[1][h]->GetBinContent(EBin);
		
		// positive charge negative eta region
		uPReTPCNegEta += fCMEQRe[0][h]->GetBinContent(EBin);
		uPImTPCNegEta += fCMEQIm[0][h]->GetBinContent(EBin);
		uP2ReTPCNegEta += fCMEQRe[0][h+1]->GetBinContent(EBin);
		uP2ImTPCNegEta += fCMEQIm[0][h+1]->GetBinContent(EBin);
		uP2Re2TPCNegEta += fCMEQRe[2][h+1]->GetBinContent(EBin);
		uP2Im2TPCNegEta += fCMEQIm[2][h+1]->GetBinContent(EBin);
		uPMTPCNegEta += fCMEMult[0][h]->GetBinContent(EBin);
		uP2MTPCNegEta += fCMEMult[2][h]->GetBinContent(EBin);
		
		// negative charge negative eta region
		uNReTPCNegEta += fCMEQRe[1][h]->GetBinContent(EBin);
		uNImTPCNegEta += fCMEQIm[1][h]->GetBinContent(EBin);
		uN2ReTPCNegEta += fCMEQRe[1][h+1]->GetBinContent(EBin);
		uN2ImTPCNegEta += fCMEQIm[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCNegEta += fCMEQRe[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCNegEta += fCMEQIm[3][h+1]->GetBinContent(EBin);
		uNMTPCNegEta += fCMEMult[1][h]->GetBinContent(EBin);
		uN2MTPCNegEta += fCMEMult[3][h]->GetBinContent(EBin);
	}
  }
  
  // set fpQvecEvent
  fpQvecEvent->setTPCRePosChPosEta( uPReTPCPosEta ); // w * cos(theta+) eta+
  fpQvecEvent->setTPCImPosChPosEta( uPImTPCPosEta ); // w * sin(theta+) eta+
  fpQvecEvent->setTPC2RePosChPosEta( uP2ReTPCPosEta ); // w * cos(2theta+) eta+
  fpQvecEvent->setTPC2ImPosChPosEta( uP2ImTPCPosEta ); // w * sin(2theta+) eta+
  fpQvecEvent->setTPC2Re2PosChPosEta( uP2Re2TPCPosEta ); // w^2 * cos(2theta+) eta+
  fpQvecEvent->setTPC2Im2PosChPosEta( uP2Im2TPCPosEta ); // w^2 * sin(2theta+) eta+
  fpQvecEvent->setTPCMPosChPosEta( uPMTPCPosEta );   // w ch+ eta+
  fpQvecEvent->setTPC2MPosChPosEta( uP2MTPCPosEta );   // w^2 ch+ eta+
  
  fpQvecEvent->setTPCRePosChNegEta( uPReTPCNegEta ); // w * cos(theta+) eta-
  fpQvecEvent->setTPCImPosChNegEta( uPImTPCNegEta ); // w * sin(theta+) eta-
  fpQvecEvent->setTPC2RePosChNegEta( uP2ReTPCNegEta ); // w * cos(2theta+) eta-
  fpQvecEvent->setTPC2ImPosChNegEta( uP2ImTPCNegEta ); // w * sin(2theta+) eta-
  fpQvecEvent->setTPC2Re2PosChNegEta( uP2Re2TPCNegEta ); // w^2 * cos(2theta+) eta-
  fpQvecEvent->setTPC2Im2PosChNegEta( uP2Im2TPCNegEta ); // w^2 * sin(2theta+) eta-
  fpQvecEvent->setTPCMPosChNegEta( uPMTPCNegEta );   // w ch+ eta-
  fpQvecEvent->setTPC2MPosChNegEta( uP2MTPCNegEta );   // w^2 ch+ eta-
  
  fpQvecEvent->setTPCReNegChPosEta( uNReTPCPosEta ); // w * cos(theta-) eta+
  fpQvecEvent->setTPCImNegChPosEta( uNImTPCPosEta ); // w * sin(theta-) eta+
  fpQvecEvent->setTPC2ReNegChPosEta( uN2ReTPCPosEta ); // w * cos(2theta-) eta+
  fpQvecEvent->setTPC2ImNegChPosEta( uN2ImTPCPosEta ); // w * sin(2theta-) eta+
  fpQvecEvent->setTPC2Re2NegChPosEta( uN2Re2TPCPosEta ); // w^2 * cos(2theta-) eta+
  fpQvecEvent->setTPC2Im2NegChPosEta( uN2Im2TPCPosEta ); // w^2 * sin(2theta-) eta+
  fpQvecEvent->setTPCMNegChPosEta( uNMTPCPosEta );   // w ch- eta+
  fpQvecEvent->setTPC2MNegChPosEta( uN2MTPCPosEta );   // w^2  h- eta+
  
  fpQvecEvent->setTPCReNegChNegEta( uNReTPCNegEta ); // w * cos(theta-) eta-
  fpQvecEvent->setTPCImNegChNegEta( uNImTPCNegEta ); // w * sin(theta-) eta-
  fpQvecEvent->setTPC2ReNegChNegEta( uN2ReTPCNegEta ); // w * cos(2theta-) eta-
  fpQvecEvent->setTPC2ImNegChNegEta( uN2ImTPCNegEta ); // w * sin(2theta-) eta-
  fpQvecEvent->setTPC2Re2NegChNegEta( uN2Re2TPCNegEta ); // w^2 * cos(2theta-) eta-
  fpQvecEvent->setTPC2Im2NegChNegEta( uN2Im2TPCNegEta ); // w^2 * sin(2theta-) eta-
  fpQvecEvent->setTPCMNegChNegEta( uNMTPCNegEta );   // w ch- eta-
  fpQvecEvent->setTPC2MNegChNegEta( uN2MTPCNegEta );   // w^2 ch- eta-

}

//=======================================================================================================================

void AliFlowAnalysisQvec::Finish()
{
  // Calculate the final results.

  // a) Check all pointers used in this method;
  // b) Access the constants;
  // c) Access the flags;
  // d) Calculate reference cumulants (not corrected for detector effects);
  // e) Correct reference cumulants for detector effects;
  // f) Calculate reference flow;
  // g) Store results for reference flow in AliFlowCommonHistResults and print them on the screen;
  // h) Calculate the final results for differential flow (without/with weights);
  // i) Correct the results for differential flow (without/with weights) for effects of non-uniform acceptance (NUA);
  // j) Calculate the final results for integrated flow (RP/POI) and store in AliFlowCommonHistResults;
  // k) Store results for differential flow in AliFlowCommonHistResults;
  // l) Print the final results for integrated flow (RP/POI) on the screen;
  // m) Cross-checking: Results from Q-vectors vs results from nested loops;
  // n) Calculate cumulants for mixed harmonics;
  // o) Calculate charge-rapidity correlations;
  // p) Calculate cumulants for bootstrap;
  // q) Finalize various;

  // a) Check all pointers used in this method:
  this->CheckPointersUsedInFinish();

} // end of AliFlowAnalysisQvec::Finish()

//=======================================================================================================================

void AliFlowAnalysisQvec::CrossCheckSettings()
{
  // a) Cross-check if the choice for multiplicity weights make sense;
  // b) Cross-check if the choice for multiplicity itself make sense.

  // a) Cross-check if the choice for multiplicity weights make sense:
  if((!fMultiplicityWeight->Contains("combinations")) &&
     (!fMultiplicityWeight->Contains("unit")) &&
     (!fMultiplicityWeight->Contains("multiplicity")) )
  {
    cout<<"WARNING (QC): Multiplicity weight can be either \"combinations\", \"unit\""<<endl;
    cout<<"              or \"multiplicity\". Certainly not \""<<fMultiplicityWeight->Data()<<"\"."<<endl;
    exit(0);
  }

} // end of void AliFlowAnalysisQvec::CrossCheckSettings()

//=======================================================================================================================

void AliFlowAnalysisQvec::CommonConstants(TString method)
{
  // Access and store common constants.

  // a) If this method was called in Init() access common constants from AliFlowCommonConstants;
  // b) If this method was called in Init() book and fill TProfile to hold constants accessed in a);
  // c) If this method was called in Finish() access common constants from TProfile booked and filled in b).

  if(method == "Init") {
    // a) If this method was called in Init() access common constants from AliFlowCommonConstants:

    // b) If this method was called in Init() book and fill TProfile to hold constants accessed in a):

  } // end of if(method == "Init")

  else if(method == "Finish") {
    // c) If this method was called in Finish() access common constants from TProfile booked and filled in b):

  } // end of else if(method == "Finish")

} // end of void AliFlowAnalysisQvec::CommonConstants(TString method)

//=======================================================================================================================

void AliFlowAnalysisQvec::GetOutputHistograms(TList *outputListHistos)
{
  // a) Get pointers for common control and common result histograms;
  // b) Get pointers for histograms holding particle weights;
  // c) Get pointers for reference flow histograms;
  // d) Get pointers for differential flow histograms;
  // e) Get pointers for 2D differential flow histograms;
  // f) Get pointers for other differential correlators;
  // g) Get pointers for mixed harmonics histograms;
  // h) Get pointers for nested loops' histograms;
  // i) Get pointers for control histograms;
  // j) Get pointers for bootstrap.
  // k) Get pointers for CRC histograms;

  if(outputListHistos)
  {
    this->SetHistList(outputListHistos);
    if(!fHistList)
    {
      printf("\n WARNING (QC): fHistList is NULL in AFAWQC::GOH() !!!!\n\n");
      exit(0);
    }
    this->GetPointersForCommonHistograms();
    this->GetPointersForParticleWeightsHistograms();
    this->GetPointersForQVec();
    this->GetPointersForVarious();
  } else
  {
    printf("\n WARNING (QC): outputListHistos is NULL in AFAWQC::GOH() !!!!\n\n");
    exit(0);
  }

} // end of void AliFlowAnalysisQvec::GetOutputHistograms(TList *outputListHistos)

//=======================================================================================================================

void AliFlowAnalysisQvec::SetCentralityWeights()
{
  if(!fCenWeightsHist) return;
  fCenWeigCalHist = (TH1D*)(fCenWeightsHist->Clone("fCenWeigCalHist"));
  TF1 *CenFit = new TF1("CenFit","pol0", 0., 100.);
  if(fDataSet!=AliFlowAnalysisQvec::k2011) {
    fCenWeigCalHist->Fit("CenFit","QNR","",0.,50.);
    Double_t CenAv = CenFit->GetParameter(0);
    for(Int_t b=1; b<=fCenWeigCalHist->GetNbinsX(); b++) {
      Double_t newbin = fCenWeigCalHist->GetBinContent(b);
      if(newbin) {
        fCenWeigCalHist->SetBinContent(b,CenAv/newbin);
      } else {
        fCenWeigCalHist->SetBinContent(b,1.);
      }
    }
  }
  else {
    fCenWeigCalHist->Fit("CenFit","QNR","",0.,8.);
    Double_t CenAv = CenFit->GetParameter(0);
    fCenWeigCalHist->Fit("CenFit","QNR","",12.,50.);
    Double_t SemiCenAv = CenFit->GetParameter(0);
    for(Int_t b=1; b<=fCenWeigCalHist->GetNbinsX(); b++) {
      Double_t newbin = fCenWeigCalHist->GetBinContent(b);
      if(newbin) {
        if(b<=10) fCenWeigCalHist->SetBinContent(b,CenAv/newbin);
        if(b>10 && b<=50) fCenWeigCalHist->SetBinContent(b,SemiCenAv/newbin);
        if(b>50) fCenWeigCalHist->SetBinContent(b,1.);
      } else {
        fCenWeigCalHist->SetBinContent(b,1.);
      }
    }
  }
  fCenWeigCalHist->SetName("CenWeights");
  fVariousList->Add(fCenWeigCalHist);
} // end of AliFlowAnalysisQvec::SetCentralityWeights()

//=======================================================================================================================

void AliFlowAnalysisQvec::ResetEventByEventQuantities()
{
  // Reset all event by event quantities.

  // Reference flow:
  //fReQ->Zero();
  //fImQ->Zero();
  //fSpk->Zero();
  fReQGF->Zero();
  fImQGF->Zero();
  for(Int_t i=0; i<fkGFPtB; i++) {
    fReQGFPt[i]->Zero();
    fImQGFPt[i]->Zero();
  }

  // Flow Vector

  for(Int_t i=0;i<2;i++) {
    fZDCFlowVect[i].Clear();
    for (Int_t h=0;h<fCRCnHar;h++) {
      fVZFlowVect[i][h].Clear();
    }
  }
  // CME
  if (fCalculateCME) {
    for(Int_t c=0;c<4;c++) {
      for (Int_t h=0;h<fCRCnHar;h++) {
        if(fCMEQRe[c][h]) fCMEQRe[c][h]->Reset();
        if(fCMEQIm[c][h]) fCMEQIm[c][h]->Reset();
        if(fCMEMult[c][h]) fCMEMult[c][h]->Reset();
      }
    }
    //@shi reset CME Qvector for spectator plane participant plane method
    for(Int_t c=0;c<2;c++) {
      for (Int_t h=0;h<fCRCnHar;h++) {
        if(fCMEQReBothCharge[c][h]) fCMEQReBothCharge[c][h]->Reset();
        if(fCMEQImBothCharge[c][h]) fCMEQImBothCharge[c][h]->Reset();
        if(fCMEMultBothCharge[c][h]) fCMEMultBothCharge[c][h]->Reset();
      }
    }
  }
  
  //delete QvecEvent;
} // end of void AliFlowAnalysisQvec::ResetEventByEventQuantities();

//=======================================================================================================================

Bool_t AliFlowAnalysisQvec::MultCut2015o()
{
  Bool_t PassMultCut = kTRUE;
  if(fCentralityEBE<90.) {
    if(fReferenceMultiplicityEBE<fMultCutMin->GetBinContent(fMultCutMin->FindBin(fCentralityEBE))) PassMultCut = kFALSE;
    if(fReferenceMultiplicityEBE>fMultCutMax->GetBinContent(fMultCutMax->FindBin(fCentralityEBE))) PassMultCut = kFALSE;
  }
  if(PassMultCut) fhCenvsMul[1]->Fill(fCentralityEBE,fReferenceMultiplicityEBE);
  return PassMultCut;
}

//=======================================================================================================================

Bool_t AliFlowAnalysisQvec::EvaulateIfSplitMergedTracks(AliFlowEventSimple* anEvent, AliFlowTrackSimple* aftsTrack, Int_t it1)
{
  const Float_t kLimit1 = 0.02 * 3;
  Float_t bSign = (fbFlagIsPosMagField? -1 : 1);

  Int_t nTracks = anEvent->NumberOfTracks();  // nPrim = total number of primary tracks

  Bool_t isNoSplit = kFALSE;

  //your cuts
  if (it1 < nTracks - 1) {
    for (Int_t itll2 = it1 + 1; itll2 < nTracks; itll2++) {

      AliFlowTrackSimple* aftsTrack2 = (AliFlowTrackSimple*)anEvent->GetTrack(itll2);
      if (!aftsTrack2) {
        delete aftsTrack2;
        continue;
      }

      if(!aftsTrack2->InPOISelection()) continue;

      Double_t deta1 = aftsTrack->Eta() - aftsTrack2->Eta();
      // phi in rad
      Float_t phi1rad1 = aftsTrack->Phi();
      Float_t phi2rad1 = aftsTrack2->Phi();
      Double_t dphi1 = TMath::ASin(TMath::Sin(phi1rad1-phi2rad1));
      Float_t dphistarminabs1 = 1e5;
      Bool_t IsNoSpliTrack = kFALSE;

      if (TMath::Abs(deta1) < 0.1 && aftsTrack->Charge()==aftsTrack2->Charge()) {

        // check first boundaries to see if is worth to loop and find the minimum
        Float_t dphistar11 = GetDPhiStar(phi1rad1, aftsTrack->Pt(), aftsTrack->Charge(), phi2rad1, aftsTrack2->Pt(), aftsTrack2->Charge(), 0.8, bSign);
        Float_t dphistar21 = GetDPhiStar(phi1rad1, aftsTrack->Pt(), aftsTrack->Charge(), phi2rad1, aftsTrack2->Pt(), aftsTrack2->Charge(), 2.5, bSign);

        if (TMath::Abs(dphistar11) < kLimit1 || TMath::Abs(dphistar21) < kLimit1 || dphistar11 * dphistar21 < 0 ) {

          for (Double_t rad1 = 0.8; rad1 < 2.51; rad1 += 0.01) {
            Float_t dphistar1 = GetDPhiStar(phi1rad1, aftsTrack->Pt(), aftsTrack->Charge(), phi2rad1, aftsTrack2->Pt(), aftsTrack2->Charge(), rad1, bSign);
            Float_t dphistarabs1 = TMath::Abs(dphistar1);
            if (dphistarabs1 < dphistarminabs1) {
              dphistarminabs1 = dphistarabs1;
            }
          }

          if (dphistarminabs1 < 0.017 && TMath::Abs(deta1) < 0.012) {
            // printf("HBT: Removed track pair %d %d with [[%f %f]] %f | %f %f %d %f %f %d %f \n", it1, itll2, TMath::Abs(deta1), TMath::Abs(phi1rad1-phi2rad1), dphistarminabs1, phi1rad1, aftsTrack->Pt(), aftsTrack->Charge(), phi2rad1, aftsTrack2->Pt(), aftsTrack2->Charge(), bSign);
            // isNoSplit = kTRUE;
            // IsNoSpliTrack = kTRUE;
          }

        }

        if (TMath::Abs(dphi1) < TMath::TwoPi()/100. && TMath::Abs(deta1) < 0.006) {
          isNoSplit = kTRUE;
          IsNoSpliTrack = kTRUE;
        }

        fTwoTrackDistanceLS[0]->Fill(deta1, dphi1, 0.5*TMath::Abs(aftsTrack->Pt()+aftsTrack2->Pt()));
        if(!IsNoSpliTrack) fTwoTrackDistanceLS[1]->Fill(deta1, dphi1, 0.5*TMath::Abs(aftsTrack->Pt()+aftsTrack2->Pt()));
      }

      IsNoSpliTrack = kFALSE;
      if (TMath::Abs(deta1) < 0.1 && aftsTrack->Charge()!=aftsTrack2->Charge()) {

        Double_t dphi1 = TMath::ASin(TMath::Sin(phi1rad1-phi2rad1));
        if (TMath::Abs(dphi1) < TMath::TwoPi()/100. && TMath::Abs(deta1) < 0.006) {
          IsNoSpliTrack = kTRUE;
        }

        fTwoTrackDistanceUS[0]->Fill(deta1, dphi1, 0.5*TMath::Abs(aftsTrack->Pt()+aftsTrack2->Pt()));
        if(!IsNoSpliTrack) fTwoTrackDistanceUS[1]->Fill(deta1, dphi1, 0.5*TMath::Abs(aftsTrack->Pt()+aftsTrack2->Pt()));
      }

    }
  }

  return isNoSplit;
}

//=======================================================================================================================

Double_t AliFlowAnalysisQvec::GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign)
{
  // calculates dphistar

  Double_t dphistar = phi1 - phi2 - charge1 * bSign * TMath::ASin(0.07510020733 * radius / pt1) + charge2 * bSign * TMath::ASin(0.07510020733 * radius / pt2);

  // circularity
  if (dphistar > TMath::Pi()) dphistar = TMath::Pi() * 2. - dphistar;
  if (dphistar < -TMath::Pi()) dphistar = -TMath::Pi() * 2. - dphistar;
  if (dphistar > TMath::Pi()) dphistar = TMath::Pi() * 2. - dphistar;

  return dphistar;
}

//=======================================================================================================================

void AliFlowAnalysisQvec::InitializeCostantsForCRC()
{
  fCorrMap = new Double_t[90];
  Double_t CorrMap[90] = {3.560467e-01, 4.546940e-01, 5.517416e-01, 6.353035e-01, 7.213249e-01, 8.144493e-01, 8.919082e-01, 9.569955e-01, 1.030360e+00, 1.102777e+00, 1.156254e+00, 1.229825e+00, 1.281096e+00, 1.348774e+00, 1.394619e+00, 1.457919e+00, 1.491964e+00, 1.535202e+00, 1.590917e+00, 1.622521e+00, 1.630109e+00, 1.664196e+00, 1.690919e+00, 1.739245e+00, 1.775229e+00, 1.787735e+00, 1.829441e+00, 1.825088e+00, 1.844297e+00, 1.848143e+00, 1.858541e+00, 1.860540e+00, 1.908495e+00, 1.898985e+00, 1.894655e+00, 1.908905e+00, 1.885567e+00, 1.891935e+00, 1.883991e+00, 1.874601e+00, 1.878468e+00, 1.863729e+00, 1.854691e+00, 1.853367e+00, 1.834473e+00, 1.818142e+00, 1.803534e+00, 1.787149e+00, 1.790569e+00, 1.754662e+00, 1.728960e+00, 1.709801e+00, 1.668947e+00, 1.664165e+00, 1.619436e+00, 1.615641e+00, 1.555998e+00, 1.541695e+00, 1.520922e+00, 1.505968e+00, 1.452864e+00, 1.407158e+00, 1.380039e+00, 1.334520e+00, 1.281019e+00, 1.263069e+00, 1.231441e+00, 1.155737e+00, 1.138802e+00, 1.090618e+00, 1.060680e+00, 1.015759e+00, 9.835138e-01, 9.273367e-01, 9.013176e-01, 8.226471e-01, 7.811623e-01, 7.459602e-01, 6.798878e-01, 6.366282e-01, 6.528206e-01, 6.415005e-01, 5.730953e-01, 5.413170e-01, 5.330439e-01, 5.183582e-01, 4.493151e-01, 4.687033e-01, 3.770695e-01, 3.543272e-01};
  for(Int_t r=0; r<90; r++) {
    fCorrMap[r] = CorrMap[r];
  }
}

//=======================================================================================================================

void AliFlowAnalysisQvec::InitializeArraysForQVec()
{
  for (Int_t i=0; i<3; i++) {
    fVZEROCenHist[i] = NULL;
  }
  
  for(Int_t r=0;r<fCRCnRun;r++) {
    fCRCQVecListRun[r] = NULL;
    
    for(Int_t i=0;i<2;i++) {
      fCRCZDCQVecA[r][i] = NULL;
      fCRCZDCQVecC[r][i] = NULL;
      fCRCZDCQVecACorr[r][i] = NULL;
      fCRCZDCQVecCCorr[r][i] = NULL;
//      fCRCVZQVecA[r][i] = NULL;
//      fCRCVZQVecC[r][i] = NULL;
    }
    
    for (Int_t i=0;i<4;i++) {
      fCRCZDCQVecRes[r][i] = NULL;
    }
  }
  
  for(Int_t k=0; k<2; k++) {
    fZDCESEAvHist[k] =  NULL;
  }
  
  for(Int_t k=0; k<fZDCESEnPol; k++) {
    fZDCESECutsHist[k] =  NULL;
  }
  
  for(Int_t k=0; k<fkNZDCResHist; k++) {
    fZDCResHist[k] = NULL;
  }
  
  for (Int_t i=0; i<4; i++) {
	  fAvr_Run_CentQ[i] = NULL;
	  fAvr_Run_VtxXYZQ[i] = NULL;
  }
  
  for (Int_t c=0; c<20; c++) {
	  for (Int_t i=0; i<4; i++) {
		  fAvr_Cent_VtxXYZQ[c][i] = NULL;
	  }
  }
  
  for(Int_t k=0; k<4; k++) {
    fZDCVtxHist[k] = NULL;
    fZDCEcomHist[k] = NULL;
    fZDCEcomTotHist[k] = NULL;
    fZDCVtxFitHist[k] = NULL;
    for(Int_t i=0; i<3; i++) {
      fZDCVtxFitCenProjHist[k][i] = NULL;
    }
    fZDCVtxFitHist2[k] = NULL;
    for(Int_t i=0; i<3; i++) {
      fZDCVtxFitCenProjHist2[k][i] = NULL;
    }
  }
  
  for(Int_t c=0; c<10; c++) {
    fZDCBinsCenRefMult[c] = NULL;
  }
  
  for(Int_t k=0; k<4; k++) {
    fZDCBinsCenRefMultRbR[k] = NULL;
    fZDCBinsCenRefMultTot[k] = NULL;
    for(Int_t c=0; c<10; c++) {
      fZDCBinsCenRefMultRbRProf[c][k] = NULL;
      fZDCBinsCenRefMultTotProf[c][k] = NULL;
      fZDCBinsCenRefMultRbRProj[c][k] = NULL;
      fZDCBinsCenRefMultTotProj[c][k] = NULL;
    }
  }
  for(Int_t i=0; i<3; i++) {
    for(Int_t k=0; k<4; k++) {
      fZDCBinsVtxCenEZDC[i][k] = NULL;
    }
  }
  for(Int_t i=0; i<10; i++) {
    for(Int_t z=0; z<10; z++) {
      for(Int_t k=0; k<4; k++) {
        fZDCQVecVtxCenEZDC3D[i][z][k] = NULL;
      }
    }
  }
  for(Int_t i=0; i<2; i++) {
    for(Int_t z=0; z<10; z++) {
      fCRCZDC2DCutZDCC[i][z] = NULL;
      fCRCZDC2DCutZDCA[i][z] = NULL;
    }
  }
  fZDCQVecVtxCenEZDCFit0 = NULL;
  fZDCQVecVtxCenEZDCFit1 = NULL;
  for(Int_t k=0; k<12; k++) {
    fZDCEcomTotvsVtxHist[k] = NULL;
  }
  
  for(Int_t c=0; c<10; c++) {
    for(Int_t k=0; k<4; k++) {
      fZDCVtxCenHist[c][k] = NULL;
    }
    for(Int_t k=0; k<8; k++) {
      fZDCVtxCenHistMagPol[c][k] = NULL;
    }
  }
  
  for(Int_t k=0; k<12; k++) {
    fZDCQHist[k] = NULL;
  }
  
  for(Int_t k=0; k<2; k++) {
    fZDCEPHist[k] = NULL;
  }
  
  for(Int_t i=0;i<2;i++) {
    for(Int_t c=0;c<fCRCnCen;c++) {
      for (Int_t j=0;j<2;j++) {
        fCRCZDCQVecDis[i][c][j] = NULL;
      }
    }
  }
  
  for(Int_t c=0;c<4;c++) {
    fCRCZDCQVecCorSteps[c] = NULL;
  }
  
  for(Int_t i=0; i<2; i++) {
    for(Int_t z=0; z<10; z++) {
      fCRCZDC2DCutZDCC[i][z] = NULL;
      fCRCZDC2DCutZDCA[i][z] = NULL;
    }
  }
}

//=======================================================================================================================

void AliFlowAnalysisQvec::InitializeArraysForVarious()
{
  fEventCounter = NULL;
  for (Int_t c=0; c<fZDCESEnCl+1; c++) {
    fhCenvsMul[c] = NULL;
  }
  fCenWeightsHist = NULL;
  fhAvRefMulRbR = NULL;
  fhAvQMCRbR = NULL;
  fhAvQMARbR = NULL;
  fhAvAbsOrbit = NULL;
  fRefMultRbRPro = NULL;
  fAvEZDCCRbRPro = NULL;
  fAvEZDCARbRPro = NULL;
  fVtxPos[0]=0.;
  fVtxPos[1]=0.;
  fVtxPos[2]=0.;
  fVtxPosCor[0]=0.;
  fVtxPosCor[1]=0.;
  fVtxPosCor[2]=0.;
  fVtxPosCor15oIRSplit[0]=0.;
  fVtxPosCor15oIRSplit[1]=0.;
  fVtxPosCor15oIRSplit[2]=0.;
  for(Int_t i=0; i<2; i++) {
    fTwoTrackDistanceLS[i] = NULL;
    fTwoTrackDistanceUS[i] = NULL;
  }
  for(Int_t k=0; k<2; k++) {
    fPolMin[k] = NULL;
    fPolMax[k] = NULL;
    fPolAv[k] = NULL;
    //    fPolDer[k] = NULL;
    //    fPolInt[k] = NULL;
    fPolDist[k] = NULL;
    fPolSlope[k] = NULL;
  }
  for (Int_t c=0; c<2; c++) {
    fhZNvsCen[c] = NULL;
  }
  for(Int_t k=0; k<fZDCESEnPol; k++) {
    fPolCuts[k] = NULL;
  }
  for(Int_t h=0; h<fCRCnCen; h++) {
    fhZNCvsZNA[h] = NULL;
  }
  for(Int_t c=0; c<10; c++) {
    fPtWeightsHist[c] = NULL;
  }
  fhZNvsMul = NULL;
  fMultCutMin = NULL;
  fMultCutMax = NULL;
  fMultCutAv = NULL;
  fEZNCutMin = NULL;
  fEZNCutMax = NULL;
  
  fCenMetric = NULL;
  
  for(Int_t k=0; k<fZDCESEnCl; k++) {
    fZDCESEMultWeightsHist[k] = NULL;
  }
} // end of InitializeArraysForVarious()

//======================================================================================================================
void AliFlowAnalysisQvec::InitializeArraysForCRC()
{
	
}
//=======================================================================================================================

void AliFlowAnalysisQvec::InitializeArraysForCRCVZ()
{
  for(Int_t c=0; c<2; c++) {
    for (Int_t h=0;h<fCRCnHar;h++) {
      fVZFlowVect[c][h] = AliFlowVector();
    }
  }

} // end of AliFlowAnalysisQvec::InitializeArraysForCRCVZ()

//=======================================================================================================================

void AliFlowAnalysisQvec::InitializeArraysForCRCZDC()
{
  for(Int_t c=0;c<2;c++) {
    fZDCFlowVect[c] = AliFlowVector();
  }
}

//=======================================================================================================================

void AliFlowAnalysisQvec::GetPointersForCommonHistograms()
{
	
}

//=======================================================================================================================

void AliFlowAnalysisQvec::GetPointersForQVec()
{
  // Q-vectors

  TList *CRCQVecList = dynamic_cast<TList*>(fHistList->FindObject("Q Vectors"));
  if (CRCQVecList) {
    this->SetCRCQVecList(CRCQVecList);
  } else {
    cout<<"WARNING: CRCQVecList is NULL in AFAWQC::GPFCRC() !!!!"<<endl;
    exit(0);
  }

  for(Int_t r=0;r<fCRCnRun;r++) {
    TList *CRCQVecListRun = dynamic_cast<TList*>(fCRCQVecList->FindObject(Form("Run %d",fRunList[r])));
    if (CRCQVecListRun) {
      this->SetCRCQVecListRun(CRCQVecListRun,r);
    } else { cout<<"WARNING: CRCQVecRunList is NULL in AFAWQC::GPFCRC() !!!!"<<endl; }

//    if(fUseCRCRecenter) {
//      for(Int_t h=0;h<fCRCnHar;h++) {
//        TProfile *CRCQnRe = dynamic_cast<TProfile*>(fCRCQVecListRun[r]->FindObject(Form("fCRCQnRe[%d][%d]",fRunList[r],h)));
//        if(CRCQnRe) { this->SetCRCQnReHist(CRCQnRe,r,h); }
//        else { cout<<"WARNING: CRCQn is NULL in AFAWQC::GPFCRC() !!!!"<<endl; }
//        TProfile *CRCQnIm = dynamic_cast<TProfile*>(fCRCQVecListRun[r]->FindObject(Form("fCRCQnIm[%d][%d]",fRunList[r],h)));
//        if(CRCQnIm) { this->SetCRCQnImHist(CRCQnIm,r,h); }
//        else { cout<<"WARNING: CRCQn is NULL in AFAWQC::GPFCRC() !!!!"<<endl; }
//        TProfile *CRCQnReCorr = dynamic_cast<TProfile*>(fCRCQVecListRun[r]->FindObject(Form("fCRCQnReCorr[%d][%d]",fRunList[r],h)));
//        if(CRCQnReCorr) { this->SetCRCQnReCorrHist(CRCQnReCorr,r,h); }
//        else { cout<<"WARNING: CRCQnCorr is NULL in AFAWQC::GPFCRC() !!!!"<<endl; }
//        TProfile *CRCQnImCorr = dynamic_cast<TProfile*>(fCRCQVecListRun[r]->FindObject(Form("fCRCQnImCorr[%d][%d]",fRunList[r],h)));
//        if(CRCQnImCorr) { this->SetCRCQnImCorrHist(CRCQnImCorr,r,h); }
//        else { cout<<"WARNING: CRCQnCorr is NULL in AFAWQC::GPFCRC() !!!!"<<endl; }
//      }
//    }
  }
  
  if (fUseZDC && fRecenterZDC) {
    for(Int_t r=0;r<fCRCnRun;r++) {
      for(Int_t i=0;i<2;i++) {
        TProfile *CRCZDCQVecA = dynamic_cast<TProfile*>(fCRCQVecListRun[r]->FindObject(Form("fCRCZDCQVecA[%d][%d]",fRunList[r],i)));
        if(CRCZDCQVecA) { this->SetCRCZDCQVecAHist(CRCZDCQVecA,r,i); }
        else { cout<<"WARNING: CRCZDCQVecA is NULL in AFAWQC::GPFCRC() !!!!"<<endl; }
        TProfile *CRCZDCQVecC = dynamic_cast<TProfile*>(fCRCQVecListRun[r]->FindObject(Form("fCRCZDCQVecC[%d][%d]",fRunList[r],i)));
        if(CRCZDCQVecC) { this->SetCRCZDCQVecCHist(CRCZDCQVecC,r,i); }
        else { cout<<"WARNING: CRCZDCQVecC is NULL in AFAWQC::GPFCRC() !!!!"<<endl; }
        TProfile *CRCZDCQVecACorr = dynamic_cast<TProfile*>(fCRCQVecListRun[r]->FindObject(Form("fCRCZDCQVecACorr[%d][%d]",fRunList[r],i)));
        if(CRCZDCQVecACorr) { this->SetCRCZDCQVecACorrHist(CRCZDCQVecACorr,r,i); }
        else { cout<<"WARNING: CRCZDCQVecACorr is NULL in AFAWQC::GPFCRC() !!!!"<<endl; }
        TProfile *CRCZDCQVecCCorr = dynamic_cast<TProfile*>(fCRCQVecListRun[r]->FindObject(Form("fCRCZDCQVecCCorr[%d][%d]",fRunList[r],i)));
        if(CRCZDCQVecCCorr) { this->SetCRCZDCQVecCCorrHist(CRCZDCQVecCCorr,r,i); }
        else { cout<<"WARNING: CRCZDCQVecCCorr is NULL in AFAWQC::GPFCRC() !!!!"<<endl; }
      }
//      for(Int_t i=0;i<4;i++) {
//        TH2D *CRCZDCQVecEP = dynamic_cast<TH2D*>(fCRCQVecListRun[r]->FindObject(Form("fCRCZDCQVecEP[%d][%d]",fRunList[r],i)));
//        if(CRCZDCQVecEP) { this->SetCRCZDCQVecEP(CRCZDCQVecEP,r,i); }
//        else { cout<<"WARNING: CRCZDCQVecEP is NULL in AFAWQC::GPFCRC() !!!!"<<endl; }
//      }
      for(Int_t i=0;i<4;i++) {
        TProfile *CRCZDCQVecRes = dynamic_cast<TProfile*>(fCRCQVecListRun[r]->FindObject(Form("fCRCZDCQVecRes[%d][%d]",fRunList[r],i)));
        if(CRCZDCQVecRes) { this->SetCRCZDCQVecRes(CRCZDCQVecRes,r,i); }
        else { cout<<"WARNING: CRCZDCQVecRes is NULL in AFAWQC::GPFCRC() !!!!"<<endl; }
      }
//      for(Int_t i=0;i<fCRCQVecnCov;i++) {
//        TProfile2D *CRCZDCQVecCov = dynamic_cast<TProfile2D*>(fCRCQVecListRun[r]->FindObject(Form("fCRCZDCQVecCov[%d][%d]",fRunList[r],i)));
//        if(CRCZDCQVecCov) { this->SetCRCZDCQVecCov(CRCZDCQVecCov,r,i); }
//        else { cout<<"WARNING: CRCZDCQVecCov is NULL in AFAWQC::GPFCRC() !!!!"<<endl; }
//      }
    }
  } // end of if (fUseZDC)
  
} // end void AliFlowAnalysisQvec::GetPointersForQVec()

//=======================================================================================================================

void AliFlowAnalysisQvec::GetPointersForVarious()
{
  // Get VariousList
  TList *VariousList = dynamic_cast<TList*>(fHistList->FindObject("Various"));
  if (VariousList) {
    this->SetVariousList(VariousList);
  } else {
    cout<<"WARNING: VariousList is NULL in AFAWQC::GPFV() !!!!"<<endl;
    exit(0);
  }
  for (Int_t c=0; c<2; c++) {
    TH2F* ZNvsCen = dynamic_cast<TH2F*>(fVariousList->FindObject(Form("fhZNvsCen[%d]",c)));
    if(ZNvsCen) this->SetZNvsCen(ZNvsCen,c);
  }
  for (Int_t c=0; c<fZDCESEnCl+1; c++) {
    TH2F* CenvsMul = dynamic_cast<TH2F*>(fVariousList->FindObject(Form("fhCenvsMul[%d]",c)));
    if(CenvsMul) this->SetCenvsMul(CenvsMul,c);
  }
  for(Int_t h=0; h<fCRCnCen; h++) {
    TH2F* ZNCvsZNA = dynamic_cast<TH2F*>(fVariousList->FindObject(Form("fhZNCvsZNA[%d]",h)));
    if(ZNCvsZNA) this->SetZNCvsZNA(ZNCvsZNA,h);
  }
  TH2F* ZNvsMul = dynamic_cast<TH2F*>(fVariousList->FindObject("fhZNvsMul"));
  if(ZNvsMul) this->SetZNvsMul(ZNvsMul);
}

//=======================================================================================================================

void AliFlowAnalysisQvec::WriteHistograms(TString outputFileName)
{
  //store the final results in output .root file
  TFile *output = new TFile(outputFileName.Data(),"RECREATE");
  //output->WriteObject(fHistList, "cobjQC","SingleKey");
  fHistList->Write(fHistList->GetName(), TObject::kSingleKey);
  delete output;
}


//=======================================================================================================================


void AliFlowAnalysisQvec::WriteHistograms(TDirectoryFile *outputFileName)
{
  //store the final results in output .root file
  fHistList->SetName("cobjQC");
  fHistList->SetOwner(kTRUE);
  outputFileName->Add(fHistList);
  outputFileName->Write(outputFileName->GetName(), TObject::kSingleKey);
}

//=======================================================================================================================

void AliFlowAnalysisQvec::BookEverythingForVarious()
{
  Double_t MulMax = ((fDataSet==k2015 || fDataSet==k2015v6 || fDataSet==k2015pidfix) ? 40.E3 : 25.E3);
  fEventCounter = new TH1D("EventCounter","EventCounter",10,0.,10.);
  fEventCounter->GetXaxis()->SetBinLabel(1,"input events");
  fEventCounter->GetXaxis()->SetBinLabel(2,"vn{SPZDC} & SC");
  fEventCounter->GetXaxis()->SetBinLabel(3,"vn{QC}");
  fVariousList->Add(fEventCounter);
  for (Int_t c=0; c<fZDCESEnCl+1; c++) {
    if(c<2) fhCenvsMul[c] = new TH2F(Form("fhCenvsMul[%d]",c), Form("fhCenvsMul[%d]",c), 100, 0., 100., 150, 0., 3000.);
    else    fhCenvsMul[c] = new TH2F(Form("fhCenvsMul[%d]",c), Form("fhCenvsMul[%d]",c), 150, 0., 3000., 150, 0., 3000.);
    fhCenvsMul[c]->Sumw2();
    fVariousList->Add(fhCenvsMul[c]);
  }
  if(fRemoveSplitMergedTracks) {
    for (Int_t c=0; c<2; c++) {
      fTwoTrackDistanceLS[c] = new TH3F(Form("fTwoTrackDistanceLS[%d]",c), ";#Delta#eta;#Delta#varphi;#Delta p_{T}/2.", 100, -0.1, 0.1, 100, -TMath::Pi()/2., TMath::Pi()/2., 24, 0.2, 5.);
      fVariousList->Add(fTwoTrackDistanceLS[c]);
      fTwoTrackDistanceUS[c] = new TH3F(Form("fTwoTrackDistanceUS[%d]",c), ";#Delta#eta;#Delta#varphi;#Delta p_{T}/2.", 100, -0.1, 0.1, 100, -TMath::Pi()/2., TMath::Pi()/2., 24, 0.2, 5.);
      fVariousList->Add(fTwoTrackDistanceUS[c]);
    }
  }
  for(Int_t h=0; h<fCRCnCen; h++) {
    fhZNCvsZNA[h] = new TH2F(Form("hZNCvsZNA[%d]",h),Form("hZNCvsZNA[%d]",h), 100, 0., 100., 100, 0., 100.);
    fVariousList->Add(fhZNCvsZNA[h]);
  }

  for (Int_t c=0; c<2; c++) {
    fhZNvsCen[c] = new TH2F(Form("fhZNvsCen[%d]",c), Form("fhZNvsCen[%d]",c), 100, 0., 100., 500, 0., 200.);
    fhZNvsCen[c]->Sumw2();
    fVariousList->Add(fhZNvsCen[c]);
  }

  fhZNvsMul = new TH2F("fhZNvsMul","fhZNvsMul", 300, 0., MulMax, 250, 0., 200.);
  fhZNvsMul->Sumw2();
  fVariousList->Add(fhZNvsMul);
  
  fMultCutMin = new TH1F("fMultCutMin","fMultCutMin",90,0.,90.);
  Double_t xcutmin[] = {1.702092e+03, 1.659690e+03, 1.594381e+03, 1.530953e+03, 1.469878e+03, 1.411837e+03, 1.354687e+03, 1.300595e+03, 1.248377e+03, 1.198430e+03, 1.149708e+03, 1.102213e+03, 1.057167e+03, 1.013113e+03, 9.712141e+02, 9.295018e+02, 8.904428e+02, 8.518189e+02, 8.149603e+02, 7.785151e+02, 7.446684e+02, 7.110022e+02, 6.784082e+02, 6.467663e+02, 6.168426e+02, 5.874845e+02, 5.590562e+02, 5.316805e+02, 5.050987e+02, 4.798622e+02, 4.549363e+02, 4.310033e+02, 4.081221e+02, 3.862957e+02, 3.646275e+02, 3.439611e+02, 3.247634e+02, 3.051531e+02, 2.865213e+02, 2.690371e+02, 2.520364e+02, 2.356551e+02, 2.201636e+02, 2.052666e+02, 1.909701e+02, 1.770509e+02, 1.641352e+02, 1.515636e+02, 1.394562e+02, 1.282475e+02, 1.172621e+02, 1.072948e+02, 9.723249e+01, 8.826328e+01, 7.950206e+01, 7.095379e+01, 6.349188e+01, 5.601530e+01, 4.910669e+01, 4.276826e+01, 3.667785e+01, 3.096218e+01, 2.581305e+01, 2.112812e+01, 1.660213e+01, 1.255342e+01, 8.751017e+00, 5.289824e+00, 2.208584e+00, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  for (Int_t i=1; i<=90; i++) { fMultCutMin->SetBinContent(i,xcutmin[i-1]); }
  fVariousList->Add(fMultCutMin);
  fMultCutMax = new TH1F("fMultCutMax","fMultCutMax",90,0.,90.);
  Double_t xcutmax[] = {2.977831e+03, 2.789981e+03, 2.672021e+03, 2.563939e+03, 2.463466e+03, 2.368142e+03, 2.279347e+03, 2.194917e+03, 2.114559e+03, 2.037801e+03, 1.964346e+03, 1.894544e+03, 1.827056e+03, 1.762559e+03, 1.699773e+03, 1.640634e+03, 1.582244e+03, 1.526511e+03, 1.472669e+03, 1.420538e+03, 1.369823e+03, 1.321126e+03, 1.274265e+03, 1.228004e+03, 1.183782e+03, 1.140569e+03, 1.098898e+03, 1.058714e+03, 1.019339e+03, 9.809862e+02, 9.442490e+02, 9.083504e+02, 8.735706e+02, 8.394027e+02, 8.070248e+02, 7.755483e+02, 7.442536e+02, 7.149877e+02, 6.864865e+02, 6.583045e+02, 6.312841e+02, 6.048965e+02, 5.792457e+02, 5.545794e+02, 5.306581e+02, 5.079286e+02, 4.853815e+02, 4.638178e+02, 4.430323e+02, 4.228056e+02, 4.034712e+02, 3.844499e+02, 3.667303e+02, 3.490117e+02, 3.322493e+02, 3.164227e+02, 3.005651e+02, 2.856931e+02, 2.713716e+02, 2.575504e+02, 2.443261e+02, 2.315942e+02, 2.193395e+02, 2.075676e+02, 1.964631e+02, 1.857374e+02, 1.755505e+02, 1.657988e+02, 1.563373e+02, 1.475873e+02, 1.390231e+02, 1.310338e+02, 1.233011e+02, 1.159695e+02, 1.090223e+02, 1.024656e+02, 9.602331e+01, 9.007991e+01, 8.425771e+01, 7.879926e+01, 7.353029e+01, 6.864368e+01, 6.395624e+01, 5.952758e+01, 5.516085e+01, 5.105927e+01, 4.718699e+01, 4.349916e+01, 4.015761e+01, 3.710833e+01};
  for (Int_t i=1; i<=90; i++) { fMultCutMax->SetBinContent(i,xcutmax[i-1]); }
  fVariousList->Add(fMultCutMax);
  fMultCutAv = new TH1F("fMultCutAv","fMultCutAv",90,0.,90.);
  Double_t xav[] = {2.324784e+03, 2.212053e+03, 2.120706e+03, 2.034993e+03, 1.954277e+03, 1.877610e+03, 1.804846e+03, 1.735610e+03, 1.669782e+03, 1.606641e+03, 1.545644e+03, 1.487097e+03, 1.431140e+03, 1.377076e+03, 1.324829e+03, 1.274444e+03, 1.226002e+03, 1.178903e+03, 1.133614e+03, 1.089623e+03, 1.047356e+03, 1.006281e+03, 9.667492e+02, 9.278635e+02, 8.907950e+02, 8.548160e+02, 8.197275e+02, 7.861789e+02, 7.532557e+02, 7.216491e+02, 6.908557e+02, 6.611987e+02, 6.324036e+02, 6.047124e+02, 5.777042e+02, 5.517775e+02, 5.266738e+02, 5.023742e+02, 4.790001e+02, 4.563775e+02, 4.344788e+02, 4.133227e+02, 3.928839e+02, 3.732591e+02, 3.542738e+02, 3.362057e+02, 3.186165e+02, 3.017719e+02, 2.855503e+02, 2.700526e+02, 2.550622e+02, 2.407724e+02, 2.271138e+02, 2.139625e+02, 2.014145e+02, 1.894186e+02, 1.779280e+02, 1.670108e+02, 1.565556e+02, 1.466505e+02, 1.371798e+02, 1.281117e+02, 1.195846e+02, 1.115543e+02, 1.038702e+02, 9.662303e+01, 8.979202e+01, 8.334000e+01, 7.720447e+01, 7.142981e+01, 6.605105e+01, 6.101889e+01, 5.625787e+01, 5.185071e+01, 4.768692e+01, 4.382575e+01, 4.016922e+01, 3.681930e+01, 3.370155e+01, 3.076244e+01, 2.804824e+01, 2.552284e+01, 2.315275e+01, 2.098654e+01, 1.904220e+01, 1.729321e+01, 1.574531e+01, 1.437631e+01, 1.320221e+01, 1.218386e+01};
  for (Int_t i=1; i<=90; i++) { fMultCutAv->SetBinContent(i,xav[i-1]); }
  fVariousList->Add(fMultCutAv);
  
  fEZNCutMin = new TH1F("fEZNCutMin","fEZNCutMin",90,0.,90.);
  Double_t ecutmin[] = {5.679458e+00, 7.607178e+00, 1.068261e+01, 1.404126e+01, 1.738144e+01, 2.053910e+01, 2.374938e+01, 2.678917e+01, 2.969818e+01, 3.224316e+01, 3.484165e+01, 3.706608e+01, 3.916280e+01, 4.150452e+01, 4.349076e+01, 4.506416e+01, 4.706689e+01, 4.854908e+01, 4.999278e+01, 5.126242e+01, 5.268097e+01, 5.404030e+01, 5.528360e+01, 5.611058e+01, 5.719087e+01, 5.779483e+01, 5.885472e+01, 5.962037e+01, 5.999710e+01, 6.095082e+01, 6.131383e+01, 6.178263e+01, 6.213514e+01, 6.273090e+01, 6.296995e+01, 6.341473e+01, 6.347308e+01, 6.371836e+01, 6.378971e+01, 6.374707e+01, 6.363084e+01, 6.344150e+01, 6.365763e+01, 6.374055e+01, 6.311442e+01, 6.304190e+01, 6.244914e+01, 6.215758e+01, 6.182734e+01, 6.109873e+01, 6.081862e+01, 6.007977e+01, 5.950013e+01, 5.840487e+01, 5.762088e+01, 5.683201e+01, 5.579488e+01, 5.453968e+01, 5.337461e+01, 5.212841e+01, 5.088454e+01, 4.925251e+01, 4.765313e+01, 4.580473e+01, 4.431336e+01, 4.231223e+01, 4.012934e+01, 3.824709e+01, 3.588116e+01, 3.379921e+01, 3.105872e+01, 2.873682e+01, 2.618550e+01, 2.325735e+01, 2.037027e+01, 1.757194e+01, 1.411580e+01, 1.130037e+01, 8.275737e+00, 5.348862e+00, 2.500213e+00, -3.652928e-01, -2.176934e+00, -5.381661e+00, -7.183125e+00, -9.352529e+00, -1.137439e+01, -1.380040e+01, -1.547430e+01, -1.662124e+01};
  for (Int_t i=1; i<=90; i++) { fEZNCutMin->SetBinContent(i,ecutmin[i-1]); }
  fVariousList->Add(fEZNCutMin);
  fEZNCutMax = new TH1F("fEZNCutMax","fEZNCutMax",90,0.,90.);
  Double_t ecutmax[] = {3.054919e+01, 3.692352e+01, 4.254322e+01, 4.777653e+01, 5.259061e+01, 5.719378e+01, 6.133480e+01, 6.523429e+01, 6.885126e+01, 7.244460e+01, 7.568765e+01, 7.893172e+01, 8.196279e+01, 8.447118e+01, 8.702351e+01, 8.971290e+01, 9.168324e+01, 9.397844e+01, 9.609081e+01, 9.811322e+01, 9.979812e+01, 1.013239e+02, 1.028019e+02, 1.045299e+02, 1.058101e+02, 1.073907e+02, 1.083592e+02, 1.094939e+02, 1.108315e+02, 1.114612e+02, 1.125904e+02, 1.134947e+02, 1.143157e+02, 1.148556e+02, 1.155402e+02, 1.159961e+02, 1.166812e+02, 1.170944e+02, 1.175560e+02, 1.180785e+02, 1.185336e+02, 1.189505e+02, 1.188809e+02, 1.188482e+02, 1.194567e+02, 1.193798e+02, 1.197406e+02, 1.197367e+02, 1.196424e+02, 1.198297e+02, 1.195133e+02, 1.195347e+02, 1.192998e+02, 1.194903e+02, 1.192825e+02, 1.190254e+02, 1.188784e+02, 1.188011e+02, 1.185301e+02, 1.182328e+02, 1.179031e+02, 1.177838e+02, 1.175172e+02, 1.174881e+02, 1.168439e+02, 1.167271e+02, 1.166109e+02, 1.161081e+02, 1.159685e+02, 1.153956e+02, 1.154842e+02, 1.149880e+02, 1.146008e+02, 1.145904e+02, 1.143930e+02, 1.140997e+02, 1.142903e+02, 1.139589e+02, 1.136466e+02, 1.133392e+02, 1.128740e+02, 1.125602e+02, 1.110557e+02, 1.110887e+02, 1.099240e+02, 1.092136e+02, 1.084315e+02, 1.081597e+02, 1.075144e+02, 1.062025e+02};
  for (Int_t i=1; i<=90; i++) { fEZNCutMax->SetBinContent(i,ecutmax[i-1]); }
  fVariousList->Add(fEZNCutMax);
}

void AliFlowAnalysisQvec::BookAndNestAllLists()
{
  // c) Book and nest list for particle weights:
  fWeightsList->SetName("Weights");
  fWeightsList->SetOwner(kTRUE);
  fTempList->Add(fWeightsList);
  
  // e) Book and nest list for various unclassified objects:
  fVariousList = new TList();
  fVariousList->SetName("Various");
  fVariousList->SetOwner(kTRUE);
  fHistList->Add(fVariousList);

  // k) Book and nest list for QEA:
  fCRCQVecList = new TList();
  fCRCQVecList->SetName("Q Vectors");
  fCRCQVecList->SetOwner(kTRUE);
  fHistList->Add(fCRCQVecList);
}

//=======================================================================================================================

void AliFlowAnalysisQvec::BookAndFillWeightsHistograms()
{
  // Book and fill histograms which hold phi, pt and eta weights.
  if(!fWeightsList)
  {
    printf("\n WARNING (QC): fWeightsList is NULL in AFAWQC::BAFWH() !!!! \n\n");
    exit(0);
  }
  
  Double_t ptbinsforweights[200001] = {0.};
  for (Int_t phib=0; phib<200001; phib++) {
    ptbinsforweights[phib] = 0.2 + phib*(50.-0.2)/200000.;
  }
  Double_t cenbinsforweights[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
  if(fPtWeightsHist[0]) {
    fPtWeightsCent = new TH2F("fPtWeightsCent","fPtWeightsCent",200000,ptbinsforweights,10,cenbinsforweights);
    fTempList->Add(fPtWeightsCent);
  }
}

//=======================================================================================================================

void AliFlowAnalysisQvec::BookEverythingForIntegratedFlow()
{
  // Book all objects for integrated flow:
  //  a) Book profile to hold all flags for integrated flow;
  //  b) Book event-by-event quantities;
  //  c) Book profiles; // to be improved (comment)
  //  d) Book histograms holding the final results.
  
  // b) Book event-by-event quantities:
  // Re[Q_{m*n,k}], Im[Q_{m*n,k}] and S_{p,k}^M:
  fReQGF = new TMatrixD(21,9);
  fImQGF = new TMatrixD(21,9);
  for(Int_t i=0; i<fkGFPtB; i++) {
    fReQGFPt[i] = new TMatrixD(21,9);
    fImQGFPt[i] = new TMatrixD(21,9);
  }
  
}

//=======================================================================================================================

void AliFlowAnalysisQvec::BookEverythingForQVec()
{
  for(Int_t k=0; k<12; k++) {
    fZDCQHist[k] = new TProfile();
    fTempList->Add(fZDCQHist[k]);
  }
  
  for(Int_t k=0; k<fkNZDCResHist; k++) {
    fZDCResHist[k] = new TH1D();
    fTempList->Add(fZDCResHist[k]);
  }
  
  for(Int_t k=0; k<4; k++) {
    fZDCVtxHist[k] = new TProfile3D();
    fTempList->Add(fZDCVtxHist[k]);
    fZDCEcomHist[k] = new TProfile2D();
    fTempList->Add(fZDCEcomHist[k]);
    fZDCEcomTotHist[k] = new TProfile2D();
    fTempList->Add(fZDCEcomTotHist[k]);
    fZDCVtxFitHist[k] = new TH3D();
    fTempList->Add(fZDCVtxFitHist[k]);
    for(Int_t i=0; i<3; i++) {
      fZDCVtxFitCenProjHist[k][i] = new TH1D();
      fTempList->Add(fZDCVtxFitCenProjHist[k][i]);
    }
    fZDCVtxFitHist2[k] = new TH3D();
    fTempList->Add(fZDCVtxFitHist2[k]);
    for(Int_t i=0; i<3; i++) {
      fZDCVtxFitCenProjHist2[k][i] = new TH1D();
      fTempList->Add(fZDCVtxFitCenProjHist2[k][i]);
    }
  }
  
  for(Int_t c=0; c<10; c++) {
    fZDCBinsCenRefMult[c] = new TH3D();
    fTempList->Add(fZDCBinsCenRefMult[c]);
  }
  
  for(Int_t k=0; k<4; k++) {
    fZDCBinsCenRefMultRbR[k] = new TProfile2D();
    fTempList->Add(fZDCBinsCenRefMultRbR[k]);
    fZDCBinsCenRefMultTot[k] = new TProfile2D();
    fTempList->Add(fZDCBinsCenRefMultTot[k]);
    for(Int_t c=0; c<10; c++) {
      fZDCBinsCenRefMultRbRProf[c][k] = new TProfile();
      fTempList->Add(fZDCBinsCenRefMultRbRProf[c][k]);
      fZDCBinsCenRefMultTotProf[c][k] = new TProfile();
      fTempList->Add(fZDCBinsCenRefMultTotProf[c][k]);
      fZDCBinsCenRefMultRbRProj[c][k] = new TH1D();
      fTempList->Add(fZDCBinsCenRefMultRbRProj[c][k]);
      fZDCBinsCenRefMultTotProj[c][k] = new TH1D();
      fTempList->Add(fZDCBinsCenRefMultTotProj[c][k]);
    }
  }
  for(Int_t i=0; i<3; i++) {
    for(Int_t k=0; k<4; k++) {
      fZDCBinsVtxCenEZDC[i][k] = new TH3D();
      fTempList->Add(fZDCBinsVtxCenEZDC[i][k]);
    }
  }
  for(Int_t i=0; i<10; i++) {
    for(Int_t z=0; z<10; z++) {
      for(Int_t k=0; k<4; k++) {
        fZDCQVecVtxCenEZDC3D[i][z][k] = new TH3D();
        fTempList->Add(fZDCQVecVtxCenEZDC3D[i][z][k]);

      }
    }
  }
  fZDCQVecVtxCenEZDCFit0 = new TH3D();
  fTempList->Add(fZDCQVecVtxCenEZDCFit0);
  fZDCQVecVtxCenEZDCFit1 = new TH3D();
  fTempList->Add(fZDCQVecVtxCenEZDCFit1);
  
  for(Int_t k=0; k<12; k++) {
    fZDCEcomTotvsVtxHist[k] = new TProfile3D();
    fTempList->Add(fZDCEcomTotvsVtxHist[k]);
  }
  
  for(Int_t c=0; c<10; c++) {
    for(Int_t k=0; k<4; k++) {
      fZDCVtxCenHist[c][k] = new TProfile3D();
      fTempList->Add(fZDCVtxCenHist[c][k]);
    }
    for(Int_t k=0; k<8; k++) {
      fZDCVtxCenHistMagPol[c][k] = new TProfile3D();
      fTempList->Add(fZDCVtxCenHistMagPol[c][k]);
    }
  }
  
  for(Int_t r=0;r<fCRCnRun;r++) {
    fCRCQVecListRun[r] = new TList();
    fCRCQVecListRun[r]->SetName(Form("Run %d",fRunList[r]));
    fCRCQVecListRun[r]->SetOwner(kTRUE);
    fCRCQVecList->Add(fCRCQVecListRun[r]);

  }

  //@Shi add my ZDC calib hists
  for (Int_t i=0; i<4; i++) {
	  fAvr_Run_CentQ[i] = new TProfile();
	  fAvr_Run_VtxXYZQ[i] = new TProfile3D();
	  fTempList->Add(fAvr_Run_CentQ[i]);
	  fTempList->Add(fAvr_Run_VtxXYZQ[i]);
  }
  for (Int_t c=0; c<20; c++) {
	  for (Int_t i=0; i<4; i++) {
		  fAvr_Cent_VtxXYZQ[c][i] = new TProfile3D();
		  fTempList->Add(fAvr_Cent_VtxXYZQ[c][i]);
	  }
  }
  
  if(fCRCZDC2DCutList) {
    for(Int_t i=0; i<2; i++) {
      for(Int_t z=0; z<10; z++) {
        fCRCZDC2DCutZDCC[i][z] = (TH2D*)(fCRCZDC2DCutList->FindObject(Form("CutZDCC_cen%d_magfielpol%d",z,i)));
        fCRCZDC2DCutZDCA[i][z] = (TH2D*)(fCRCZDC2DCutList->FindObject(Form("CutZDCA_cen%d_magfielpol%d",z,i)));
      }
    }
  }
  
  if(fZDCESEList) { //@Shi many histograms are unnecessary
    fZDCESEAvHist[0] = (TH1D*)(fZDCESEList->FindObject("AvMulHis[0]"));
    fZDCESEAvHist[1] = (TH1D*)(fZDCESEList->FindObject("AvMulHis[1]"));
    if(fZDCESEAvHist[0]) fTempList->Add(fZDCESEAvHist[0]);
    else printf("WARNING: fZDCESEAvHist not found! \n");
    if(fZDCESEAvHist[1]) fTempList->Add(fZDCESEAvHist[1]);
    else printf("WARNING: fZDCESEAvHist not found! \n");

    // functions
    for (Int_t k=0; k<2; k++) {
      fPolMin[k] = new TF1(Form("fPolMin[%d]",k),"pol9",0.,100.);
      fPolMax[k] = new TF1(Form("fPolMax[%d]",k),"pol9",0.,100.);
      fTempList->Add(fPolMin[k]);
      fTempList->Add(fPolMax[k]);

      // to calculate new metric
      fPolAv[k] = new TF1(Form("fPolAv[%d]",k),"pol4",5.E3,2.5E4);
      fPolDist[k] = new TF1(Form("fPolDist[%d]",k),"TMath::Abs(pol1(0)-pol1(2))",5.E3,2.5E4);
      if(fZDCESEAvHist[k]) {
        fPolAv[k]->SetParameter(0,4.92958e+01);
        fPolAv[k]->SetParameter(1,-2.28486e-03);
        for (Int_t p=0; p<=fPolAv[k]->GetNpar(); p++) {
          fPolDist[k]->SetParameter(p+2,fPolAv[k]->GetParameter(p));
        }
      }
      fTempList->Add(fPolAv[k]);
      fTempList->Add(fPolDist[k]);

      // new cuts on ZDC-ESE
      for(Int_t k=0; k<fZDCESEnPol; k++) {
        fZDCESECutsHist[k] = (TH1D*)(fZDCESEList->FindObject(Form("CutsHis[%d]",k)));
        if(fZDCESECutsHist[k]) fTempList->Add(fZDCESECutsHist[k]);
        else printf("WARNING: fZDCESECutsHist not found! \n");

        fPolCuts[k] = new TF1(Form("fPolCuts[%d]",k),"pol9",0.,100.);
        if(fZDCESECutsHist[k]) {
          fZDCESECutsHist[k]->Fit(fPolCuts[k],"QRN","",0.,90.);
        }
        fTempList->Add(fPolCuts[k]);
      }

      const Int_t nbins = 3410;

      Double_t xd[nbins] = {1.800005e+00, 1.800005e+00, 1.800005e+00, 5.400014e+00, 1.080003e+01, 1.440004e+01, 1.800005e+01, 2.160006e+01, 2.520007e+01, 2.700007e+01, 3.060008e+01, 3.240008e+01, 3.600009e+01, 3.780010e+01, 4.140011e+01, 4.320011e+01, 4.500012e+01, 4.860013e+01, 5.040013e+01, 5.220014e+01, 5.580015e+01, 5.760015e+01, 5.940016e+01, 6.120016e+01, 6.480017e+01, 6.660017e+01, 6.840018e+01, 7.020018e+01, 7.380019e+01, 7.560020e+01, 7.740020e+01, 7.920021e+01, 8.280022e+01, 8.460022e+01, 8.640023e+01, 8.820023e+01, 9.180024e+01, 9.360024e+01, 9.540025e+01, 9.900026e+01, 1.008003e+02, 1.026003e+02, 1.044003e+02, 1.080003e+02, 1.098003e+02, 1.116003e+02, 1.134003e+02, 1.152003e+02, 1.188003e+02, 1.206003e+02, 1.224003e+02, 1.242003e+02, 1.278003e+02, 1.296003e+02, 1.314003e+02, 1.332003e+02, 1.368004e+02, 1.386004e+02, 1.404004e+02, 1.422004e+02, 1.458004e+02, 1.476004e+02, 1.494004e+02, 1.530004e+02, 1.548004e+02, 1.566004e+02, 1.584004e+02, 1.620004e+02, 1.638004e+02, 1.656004e+02, 1.674004e+02, 1.692004e+02, 1.728005e+02, 1.746005e+02, 1.764005e+02, 1.800005e+02, 1.818005e+02, 1.836005e+02, 1.854005e+02, 1.890005e+02, 1.908005e+02, 1.926005e+02, 1.944005e+02, 1.980005e+02, 1.998005e+02, 2.016005e+02, 2.034005e+02, 2.070005e+02, 2.088005e+02, 2.106005e+02, 2.124006e+02, 2.160006e+02, 2.178006e+02, 2.196006e+02, 2.214006e+02, 2.250006e+02, 2.268006e+02, 2.286006e+02, 2.304006e+02, 2.340006e+02, 2.358006e+02, 2.376006e+02, 2.412006e+02, 2.430006e+02, 2.448006e+02, 2.466006e+02, 2.502007e+02, 2.520007e+02, 2.538007e+02, 2.556007e+02, 2.592007e+02, 2.610007e+02, 2.628007e+02, 2.646007e+02, 2.682007e+02, 2.700007e+02, 2.718007e+02, 2.754007e+02, 2.772007e+02, 2.790007e+02, 2.808007e+02, 2.844007e+02, 2.862007e+02, 2.880008e+02, 2.916008e+02, 2.934008e+02, 2.952008e+02, 2.970008e+02, 3.006008e+02, 3.024008e+02, 3.042008e+02, 3.078008e+02, 3.096008e+02, 3.114008e+02, 3.132008e+02, 3.168008e+02, 3.186008e+02, 3.204008e+02, 3.240008e+02, 3.258009e+02, 3.276009e+02, 3.294009e+02, 3.330009e+02, 3.348009e+02, 3.366009e+02, 3.402009e+02, 3.420009e+02, 3.438009e+02, 3.456009e+02, 3.492009e+02, 3.510009e+02, 3.528009e+02, 3.564009e+02, 3.582009e+02, 3.600009e+02, 3.618009e+02, 3.654010e+02, 3.672010e+02, 3.690010e+02, 3.726010e+02, 3.744010e+02, 3.762010e+02, 3.780010e+02, 3.816010e+02, 3.834010e+02, 3.852010e+02, 3.888010e+02, 3.906010e+02, 3.924010e+02, 3.960010e+02, 3.978010e+02, 3.996010e+02, 4.014010e+02, 4.050011e+02, 4.068011e+02, 4.086011e+02, 4.122011e+02, 4.140011e+02, 4.158011e+02, 4.194011e+02, 4.212011e+02, 4.230011e+02, 4.266011e+02, 4.284011e+02, 4.302011e+02, 4.320011e+02, 4.356011e+02, 4.374011e+02, 4.392011e+02, 4.428012e+02, 4.446012e+02, 4.464012e+02, 4.500012e+02, 4.518012e+02, 4.536012e+02, 4.554012e+02, 4.590012e+02, 4.608012e+02, 4.626012e+02, 4.662012e+02, 4.680012e+02, 4.698012e+02, 4.734012e+02, 4.752012e+02, 4.770012e+02, 4.806013e+02, 4.824013e+02, 4.842013e+02, 4.860013e+02, 4.896013e+02, 4.914013e+02, 4.932013e+02, 4.968013e+02, 4.986013e+02, 5.004013e+02, 5.040013e+02, 5.058013e+02, 5.076013e+02, 5.112013e+02, 5.130013e+02, 5.148013e+02, 5.184014e+02, 5.202014e+02, 5.220014e+02, 5.256014e+02, 5.274014e+02, 5.292014e+02, 5.328014e+02, 5.346014e+02, 5.364014e+02, 5.400014e+02, 5.418014e+02, 5.436014e+02, 5.472014e+02, 5.490014e+02, 5.508014e+02, 5.544014e+02, 5.562015e+02, 5.580015e+02, 5.616015e+02, 5.634015e+02, 5.652015e+02, 5.688015e+02, 5.706015e+02, 5.724015e+02, 5.760015e+02, 5.778015e+02, 5.796015e+02, 5.832015e+02, 5.850015e+02, 5.868015e+02, 5.886015e+02, 5.922015e+02, 5.940016e+02, 5.958016e+02, 5.994016e+02, 6.012016e+02, 6.030016e+02, 6.066016e+02, 6.084016e+02, 6.102016e+02, 6.138016e+02, 6.156016e+02, 6.174016e+02, 6.210016e+02, 6.228016e+02, 6.246016e+02, 6.282016e+02, 6.300016e+02, 6.318016e+02, 6.354017e+02, 6.372017e+02, 6.390017e+02, 6.426017e+02, 6.444017e+02, 6.462017e+02, 6.498017e+02, 6.516017e+02, 6.534017e+02, 6.570017e+02, 6.588017e+02, 6.606017e+02, 6.642017e+02, 6.660017e+02, 6.696017e+02, 6.714018e+02, 6.732018e+02, 6.768018e+02, 6.786018e+02, 6.804018e+02, 6.840018e+02, 6.858018e+02, 6.876018e+02, 6.912018e+02, 6.930018e+02, 6.948018e+02, 6.984018e+02, 7.002018e+02, 7.020018e+02, 7.056018e+02, 7.074018e+02, 7.092019e+02, 7.128019e+02, 7.146019e+02, 7.164019e+02, 7.200019e+02, 7.218019e+02, 7.236019e+02, 7.272019e+02, 7.290019e+02, 7.326019e+02, 7.344019e+02, 7.362019e+02, 7.380019e+02, 7.416019e+02, 7.434019e+02, 7.470019e+02, 7.488020e+02, 7.506020e+02, 7.542020e+02, 7.560020e+02, 7.578020e+02, 7.614020e+02, 7.632020e+02, 7.650020e+02, 7.686020e+02, 7.704020e+02, 7.740020e+02, 7.758020e+02, 7.776020e+02, 7.812020e+02, 7.830020e+02, 7.848020e+02, 7.884021e+02, 7.902021e+02, 7.920021e+02, 7.956021e+02, 7.974021e+02, 8.010021e+02, 8.028021e+02, 8.046021e+02, 8.082021e+02, 8.100021e+02, 8.118021e+02, 8.154021e+02, 8.172021e+02, 8.190021e+02, 8.226021e+02, 8.244022e+02, 8.280022e+02, 8.298022e+02, 8.316022e+02, 8.352022e+02, 8.370022e+02, 8.406022e+02, 8.424022e+02, 8.442022e+02, 8.478022e+02, 8.496022e+02, 8.514022e+02, 8.550022e+02, 8.568022e+02, 8.604022e+02, 8.622023e+02, 8.640023e+02, 8.676023e+02, 8.694023e+02, 8.712023e+02, 8.748023e+02, 8.766023e+02, 8.802023e+02, 8.820023e+02, 8.838023e+02, 8.874023e+02, 8.892023e+02, 8.910023e+02, 8.946023e+02, 8.964023e+02, 9.000023e+02, 9.018024e+02, 9.036024e+02, 9.072024e+02, 9.090024e+02, 9.108024e+02, 9.144024e+02, 9.162024e+02, 9.198024e+02, 9.216024e+02, 9.234024e+02, 9.270024e+02, 9.288024e+02, 9.324024e+02, 9.342024e+02, 9.360024e+02, 9.396025e+02, 9.414025e+02, 9.432025e+02, 9.468025e+02, 9.486025e+02, 9.522025e+02, 9.540025e+02, 9.576025e+02, 9.594025e+02, 9.612025e+02, 9.648025e+02, 9.666025e+02, 9.702025e+02, 9.720025e+02, 9.738025e+02, 9.774026e+02, 9.792026e+02, 9.828026e+02, 9.846026e+02, 9.864026e+02, 9.900026e+02, 9.918026e+02, 9.954026e+02, 9.972026e+02, 9.990026e+02, 1.002603e+03, 1.004403e+03, 1.006203e+03, 1.009803e+03, 1.011603e+03, 1.015203e+03, 1.017003e+03, 1.018803e+03, 1.022403e+03, 1.024203e+03, 1.027803e+03, 1.029603e+03, 1.031403e+03, 1.035003e+03, 1.036803e+03, 1.040403e+03, 1.042203e+03, 1.044003e+03, 1.047603e+03, 1.049403e+03, 1.053003e+03, 1.054803e+03, 1.058403e+03, 1.060203e+03, 1.062003e+03, 1.065603e+03, 1.067403e+03, 1.071003e+03, 1.072803e+03, 1.076403e+03, 1.078203e+03, 1.080003e+03, 1.083603e+03, 1.085403e+03, 1.089003e+03, 1.090803e+03, 1.094403e+03, 1.096203e+03, 1.098003e+03, 1.101603e+03, 1.103403e+03, 1.105203e+03, 1.108803e+03, 1.110603e+03, 1.114203e+03, 1.116003e+03, 1.117803e+03, 1.121403e+03, 1.123203e+03, 1.126803e+03, 1.128603e+03, 1.132203e+03, 1.134003e+03, 1.135803e+03, 1.139403e+03, 1.141203e+03, 1.144803e+03, 1.146603e+03, 1.150203e+03, 1.152003e+03, 1.153803e+03, 1.157403e+03, 1.159203e+03, 1.162803e+03, 1.164603e+03, 1.168203e+03, 1.170003e+03, 1.173603e+03, 1.175403e+03, 1.177203e+03, 1.180803e+03, 1.182603e+03, 1.186203e+03, 1.188003e+03, 1.191603e+03, 1.193403e+03, 1.195203e+03, 1.198803e+03, 1.200603e+03, 1.204203e+03, 1.206003e+03, 1.209603e+03, 1.211403e+03, 1.215003e+03, 1.216803e+03, 1.218603e+03, 1.222203e+03, 1.224003e+03, 1.227603e+03, 1.229403e+03, 1.233003e+03, 1.234803e+03, 1.238403e+03, 1.240203e+03, 1.243803e+03, 1.245603e+03, 1.247403e+03, 1.251003e+03, 1.252803e+03, 1.256403e+03, 1.258203e+03, 1.261803e+03, 1.263603e+03, 1.265403e+03, 1.269003e+03, 1.270803e+03, 1.274403e+03, 1.276203e+03, 1.279803e+03, 1.281603e+03, 1.285203e+03, 1.287003e+03, 1.290603e+03, 1.292403e+03, 1.296003e+03, 1.297803e+03, 1.299603e+03, 1.303203e+03, 1.305003e+03, 1.308603e+03, 1.310403e+03, 1.314003e+03, 1.315803e+03, 1.319403e+03, 1.321203e+03, 1.324803e+03, 1.326603e+03, 1.328403e+03, 1.332003e+03, 1.333803e+03, 1.337403e+03, 1.339203e+03, 1.342804e+03, 1.344604e+03, 1.348204e+03, 1.350004e+03, 1.353604e+03, 1.355404e+03, 1.359004e+03, 1.360804e+03, 1.364404e+03, 1.366204e+03, 1.368004e+03, 1.371604e+03, 1.373404e+03, 1.377004e+03, 1.378804e+03, 1.382404e+03, 1.384204e+03, 1.387804e+03, 1.389604e+03, 1.393204e+03, 1.395004e+03, 1.398604e+03, 1.400404e+03, 1.404004e+03, 1.405804e+03, 1.407604e+03, 1.411204e+03, 1.413004e+03, 1.416604e+03, 1.418404e+03, 1.422004e+03, 1.423804e+03, 1.427404e+03, 1.429204e+03, 1.432804e+03, 1.434604e+03, 1.438204e+03, 1.440004e+03, 1.443604e+03, 1.445404e+03, 1.449004e+03, 1.450804e+03, 1.454404e+03, 1.456204e+03, 1.458004e+03, 1.461604e+03, 1.463404e+03, 1.467004e+03, 1.468804e+03, 1.472404e+03, 1.474204e+03, 1.477804e+03, 1.479604e+03, 1.483204e+03, 1.485004e+03, 1.488604e+03, 1.490404e+03, 1.494004e+03, 1.495804e+03, 1.499404e+03, 1.501204e+03, 1.504804e+03, 1.506604e+03, 1.510204e+03, 1.512004e+03, 1.515604e+03, 1.517404e+03, 1.521004e+03, 1.522804e+03, 1.526404e+03, 1.528204e+03, 1.531804e+03, 1.533604e+03, 1.537204e+03, 1.539004e+03, 1.542604e+03, 1.544404e+03, 1.548004e+03, 1.549804e+03, 1.553404e+03, 1.555204e+03, 1.558804e+03, 1.560604e+03, 1.564204e+03, 1.566004e+03, 1.569604e+03, 1.571404e+03, 1.575004e+03, 1.576804e+03, 1.580404e+03, 1.582204e+03, 1.585804e+03, 1.587604e+03, 1.591204e+03, 1.593004e+03, 1.596604e+03, 1.598404e+03, 1.602004e+03, 1.603804e+03, 1.607404e+03, 1.609204e+03, 1.612804e+03, 1.616404e+03, 1.618204e+03, 1.621804e+03, 1.623604e+03, 1.627204e+03, 1.629004e+03, 1.632604e+03, 1.634404e+03, 1.638004e+03, 1.639804e+03, 1.641604e+03, 1.645204e+03, 1.647004e+03, 1.650604e+03, 1.654204e+03, 1.656004e+03, 1.659604e+03, 1.661404e+03, 1.665004e+03, 1.666804e+03, 1.670404e+03, 1.672204e+03, 1.675804e+03, 1.677604e+03, 1.681204e+03, 1.683004e+03, 1.686604e+03, 1.688404e+03, 1.692004e+03, 1.693804e+03, 1.697404e+03, 1.699204e+03, 1.702804e+03, 1.706404e+03, 1.708204e+03, 1.710004e+03, 1.713604e+03, 1.715404e+03, 1.719004e+03, 1.722604e+03, 1.724405e+03, 1.728005e+03, 1.729805e+03, 1.733405e+03, 1.735205e+03, 1.738805e+03, 1.740605e+03, 1.744205e+03, 1.746005e+03, 1.749605e+03, 1.751405e+03, 1.755005e+03, 1.756805e+03, 1.760405e+03, 1.762205e+03, 1.765805e+03, 1.767605e+03, 1.771205e+03, 1.773005e+03, 1.776605e+03, 1.778405e+03, 1.782005e+03, 1.785605e+03, 1.787405e+03, 1.791005e+03, 1.792805e+03, 1.796405e+03, 1.798205e+03, 1.801805e+03, 1.803605e+03, 1.807205e+03, 1.809005e+03, 1.812605e+03, 1.814405e+03, 1.818005e+03, 1.819805e+03, 1.823405e+03, 1.827005e+03, 1.828805e+03, 1.832405e+03, 1.834205e+03, 1.837805e+03, 1.839605e+03, 1.843205e+03, 1.845005e+03, 1.848605e+03, 1.850405e+03, 1.854005e+03, 1.857605e+03, 1.859405e+03, 1.861205e+03, 1.864805e+03, 1.868405e+03, 1.870205e+03, 1.873805e+03, 1.875605e+03, 1.879205e+03, 1.881005e+03, 1.884605e+03, 1.888205e+03, 1.890005e+03, 1.893605e+03, 1.895405e+03, 1.899005e+03, 1.902605e+03, 1.904405e+03, 1.908005e+03, 1.909805e+03, 1.913405e+03, 1.915205e+03, 1.918805e+03, 1.920605e+03, 1.924205e+03, 1.926005e+03, 1.929605e+03, 1.933205e+03, 1.935005e+03, 1.938605e+03, 1.940405e+03, 1.944005e+03, 1.945805e+03, 1.949405e+03, 1.953005e+03, 1.954805e+03, 1.958405e+03, 1.960205e+03, 1.963805e+03, 1.965605e+03, 1.969205e+03, 1.972805e+03, 1.974605e+03, 1.978205e+03, 1.980005e+03, 1.983605e+03, 1.985405e+03, 1.989005e+03, 1.992605e+03, 1.994405e+03, 1.998005e+03, 1.999805e+03, 2.003405e+03, 2.005205e+03, 2.008805e+03, 2.010605e+03, 2.014205e+03, 2.017805e+03, 2.019605e+03, 2.023205e+03, 2.025005e+03, 2.028605e+03, 2.030405e+03, 2.034005e+03, 2.037605e+03, 2.039405e+03, 2.043005e+03, 2.046605e+03, 2.048405e+03, 2.052005e+03, 2.053805e+03, 2.057405e+03, 2.059205e+03, 2.062805e+03, 2.066405e+03, 2.068205e+03, 2.071805e+03, 2.073605e+03, 2.077205e+03, 2.080805e+03, 2.082605e+03, 2.086205e+03, 2.088005e+03, 2.091605e+03, 2.095205e+03, 2.097005e+03, 2.100605e+03, 2.102405e+03, 2.106005e+03, 2.109606e+03, 2.111406e+03, 2.115006e+03, 2.116806e+03, 2.120406e+03, 2.122206e+03, 2.125806e+03, 2.129406e+03, 2.131206e+03, 2.134806e+03, 2.136606e+03, 2.140206e+03, 2.143806e+03, 2.145606e+03, 2.149206e+03, 2.152806e+03, 2.154606e+03, 2.158206e+03, 2.160006e+03, 2.163606e+03, 2.167206e+03, 2.169006e+03, 2.172606e+03, 2.174406e+03, 2.178006e+03, 2.181606e+03, 2.183406e+03, 2.187006e+03, 2.188806e+03, 2.192406e+03, 2.196006e+03, 2.197806e+03, 2.201406e+03, 2.205006e+03, 2.206806e+03, 2.210406e+03, 2.212206e+03, 2.215806e+03, 2.219406e+03, 2.221206e+03, 2.224806e+03, 2.226606e+03, 2.230206e+03, 2.233806e+03, 2.235606e+03, 2.239206e+03, 2.241006e+03, 2.244606e+03, 2.248206e+03, 2.250006e+03, 2.253606e+03, 2.257206e+03, 2.259006e+03, 2.262606e+03, 2.264406e+03, 2.268006e+03, 2.271606e+03, 2.273406e+03, 2.277006e+03, 2.280606e+03, 2.282406e+03, 2.286006e+03, 2.287806e+03, 2.291406e+03, 2.295006e+03, 2.296806e+03, 2.300406e+03, 2.304006e+03, 2.305806e+03, 2.309406e+03, 2.313006e+03, 2.314806e+03, 2.318406e+03, 2.320206e+03, 2.323806e+03, 2.327406e+03, 2.329206e+03, 2.332806e+03, 2.336406e+03, 2.338206e+03, 2.341806e+03, 2.345406e+03, 2.347206e+03, 2.350806e+03, 2.354406e+03, 2.356206e+03, 2.359806e+03, 2.363406e+03, 2.365206e+03, 2.368806e+03, 2.372406e+03, 2.374206e+03, 2.377806e+03, 2.381406e+03, 2.383206e+03, 2.386806e+03, 2.390406e+03, 2.392206e+03, 2.395806e+03, 2.397606e+03, 2.401206e+03, 2.404806e+03, 2.406606e+03, 2.410206e+03, 2.413806e+03, 2.415606e+03, 2.419206e+03, 2.422806e+03, 2.424606e+03, 2.428206e+03, 2.431806e+03, 2.433606e+03, 2.437206e+03, 2.439006e+03, 2.442606e+03, 2.446206e+03, 2.448006e+03, 2.451606e+03, 2.455206e+03, 2.457006e+03, 2.460606e+03, 2.464206e+03, 2.466006e+03, 2.469606e+03, 2.473206e+03, 2.475006e+03, 2.478606e+03, 2.480406e+03, 2.484006e+03, 2.487606e+03, 2.489406e+03, 2.493007e+03, 2.496607e+03, 2.498407e+03, 2.502007e+03, 2.505607e+03, 2.507407e+03, 2.511007e+03, 2.514607e+03, 2.516407e+03, 2.520007e+03, 2.523607e+03, 2.527207e+03, 2.529007e+03, 2.532607e+03, 2.536207e+03, 2.538007e+03, 2.541607e+03, 2.545207e+03, 2.547007e+03, 2.550607e+03, 2.554207e+03, 2.557807e+03, 2.559607e+03, 2.563207e+03, 2.566807e+03, 2.568607e+03, 2.572207e+03, 2.575807e+03, 2.577607e+03, 2.581207e+03, 2.584807e+03, 2.588407e+03, 2.590207e+03, 2.593807e+03, 2.597407e+03, 2.599207e+03, 2.602807e+03, 2.606407e+03, 2.610007e+03, 2.611807e+03, 2.615407e+03, 2.619007e+03, 2.620807e+03, 2.624407e+03, 2.628007e+03, 2.629807e+03, 2.633407e+03, 2.637007e+03, 2.640607e+03, 2.642407e+03, 2.646007e+03, 2.649607e+03, 2.651407e+03, 2.655007e+03, 2.658607e+03, 2.660407e+03, 2.664007e+03, 2.667607e+03, 2.669407e+03, 2.673007e+03, 2.676607e+03, 2.680207e+03, 2.682007e+03, 2.685607e+03, 2.689207e+03, 2.692807e+03, 2.694607e+03, 2.698207e+03, 2.701807e+03, 2.705407e+03, 2.707207e+03, 2.710807e+03, 2.714407e+03, 2.716207e+03, 2.719807e+03, 2.723407e+03, 2.727007e+03, 2.728807e+03, 2.732407e+03, 2.736007e+03, 2.737807e+03, 2.741407e+03, 2.745007e+03, 2.748607e+03, 2.750407e+03, 2.754007e+03, 2.757607e+03, 2.759407e+03, 2.763007e+03, 2.766607e+03, 2.770207e+03, 2.772007e+03, 2.775607e+03, 2.779207e+03, 2.781007e+03, 2.784607e+03, 2.788207e+03, 2.790007e+03, 2.793607e+03, 2.797207e+03, 2.800807e+03, 2.802607e+03, 2.806207e+03, 2.809807e+03, 2.813407e+03, 2.815207e+03, 2.818807e+03, 2.822407e+03, 2.824207e+03, 2.827807e+03, 2.831407e+03, 2.833207e+03, 2.836807e+03, 2.840407e+03, 2.844007e+03, 2.845807e+03, 2.849407e+03, 2.853007e+03, 2.856607e+03, 2.858407e+03, 2.862007e+03, 2.865607e+03, 2.867407e+03, 2.871007e+03, 2.874608e+03, 2.878208e+03, 2.880008e+03, 2.883608e+03, 2.887208e+03, 2.890808e+03, 2.892608e+03, 2.896208e+03, 2.899808e+03, 2.903408e+03, 2.905208e+03, 2.908808e+03, 2.912408e+03, 2.916008e+03, 2.917808e+03, 2.921408e+03, 2.925008e+03, 2.928608e+03, 2.930408e+03, 2.934008e+03, 2.937608e+03, 2.941208e+03, 2.943008e+03, 2.946608e+03, 2.950208e+03, 2.953808e+03, 2.955608e+03, 2.959208e+03, 2.962808e+03, 2.966408e+03, 2.968208e+03, 2.971808e+03, 2.975408e+03, 2.979008e+03, 2.980808e+03, 2.984408e+03, 2.988008e+03, 2.991608e+03, 2.993408e+03, 2.997008e+03, 3.000608e+03, 3.004208e+03, 3.006008e+03, 3.009608e+03, 3.013208e+03, 3.016808e+03, 3.018608e+03, 3.022208e+03, 3.025808e+03, 3.029408e+03, 3.033008e+03, 3.034808e+03, 3.038408e+03, 3.042008e+03, 3.045608e+03, 3.047408e+03, 3.051008e+03, 3.054608e+03, 3.058208e+03, 3.061808e+03, 3.063608e+03, 3.067208e+03, 3.070808e+03, 3.074408e+03, 3.076208e+03, 3.079808e+03, 3.083408e+03, 3.087008e+03, 3.088808e+03, 3.092408e+03, 3.096008e+03, 3.099608e+03, 3.101408e+03, 3.105008e+03, 3.108608e+03, 3.112208e+03, 3.115808e+03, 3.117608e+03, 3.121208e+03, 3.124808e+03, 3.128408e+03, 3.130208e+03, 3.133808e+03, 3.137408e+03, 3.141008e+03, 3.142808e+03, 3.146408e+03, 3.150008e+03, 3.153608e+03, 3.155408e+03, 3.159008e+03, 3.162608e+03, 3.166208e+03, 3.168008e+03, 3.171608e+03, 3.175208e+03, 3.178808e+03, 3.182408e+03, 3.184208e+03, 3.187808e+03, 3.191408e+03, 3.195008e+03, 3.198608e+03, 3.200408e+03, 3.204008e+03, 3.207608e+03, 3.211208e+03, 3.214808e+03, 3.216608e+03, 3.220208e+03, 3.223808e+03, 3.227408e+03, 3.231008e+03, 3.234608e+03, 3.236408e+03, 3.240008e+03, 3.243608e+03, 3.247208e+03, 3.250808e+03, 3.252608e+03, 3.256208e+03, 3.259809e+03, 3.263409e+03, 3.267009e+03, 3.268809e+03, 3.272409e+03, 3.276009e+03, 3.279609e+03, 3.283209e+03, 3.285009e+03, 3.288609e+03, 3.292209e+03, 3.295809e+03, 3.297609e+03, 3.301209e+03, 3.304809e+03, 3.308409e+03, 3.312009e+03, 3.313809e+03, 3.317409e+03, 3.321009e+03, 3.324609e+03, 3.328209e+03, 3.330009e+03, 3.333609e+03, 3.337209e+03, 3.340809e+03, 3.344409e+03, 3.346209e+03, 3.349809e+03, 3.353409e+03, 3.357009e+03, 3.360609e+03, 3.362409e+03, 3.366009e+03, 3.369609e+03, 3.373209e+03, 3.376809e+03, 3.380409e+03, 3.382209e+03, 3.385809e+03, 3.389409e+03, 3.393009e+03, 3.396609e+03, 3.398409e+03, 3.402009e+03, 3.405609e+03, 3.409209e+03, 3.412809e+03, 3.416409e+03, 3.418209e+03, 3.421809e+03, 3.425409e+03, 3.429009e+03, 3.432609e+03, 3.434409e+03, 3.438009e+03, 3.441609e+03, 3.445209e+03, 3.447009e+03, 3.450609e+03, 3.454209e+03, 3.457809e+03, 3.461409e+03, 3.465009e+03, 3.468609e+03, 3.470409e+03, 3.474009e+03, 3.477609e+03, 3.481209e+03, 3.484809e+03, 3.488409e+03, 3.490209e+03, 3.493809e+03, 3.497409e+03, 3.501009e+03, 3.504609e+03, 3.508209e+03, 3.511809e+03, 3.513609e+03, 3.517209e+03, 3.520809e+03, 3.524409e+03, 3.528009e+03, 3.531609e+03, 3.535209e+03, 3.537009e+03, 3.540609e+03, 3.544209e+03, 3.547809e+03, 3.551409e+03, 3.555009e+03, 3.556809e+03, 3.560409e+03, 3.564009e+03, 3.567609e+03, 3.571209e+03, 3.574809e+03, 3.576609e+03, 3.580209e+03, 3.583809e+03, 3.587409e+03, 3.591009e+03, 3.594609e+03, 3.598209e+03, 3.601809e+03, 3.605409e+03, 3.609009e+03, 3.610809e+03, 3.614409e+03, 3.618009e+03, 3.621609e+03, 3.625209e+03, 3.628809e+03, 3.632409e+03, 3.634209e+03, 3.637809e+03, 3.641410e+03, 3.645010e+03, 3.648610e+03, 3.652210e+03, 3.655810e+03, 3.659410e+03, 3.661210e+03, 3.664810e+03, 3.668410e+03, 3.672010e+03, 3.675610e+03, 3.679210e+03, 3.682810e+03, 3.686410e+03, 3.688210e+03, 3.691810e+03, 3.695410e+03, 3.699010e+03, 3.702610e+03, 3.706210e+03, 3.709810e+03, 3.713410e+03, 3.715210e+03, 3.718810e+03, 3.722410e+03, 3.726010e+03, 3.729610e+03, 3.733210e+03, 3.736810e+03, 3.740410e+03, 3.744010e+03, 3.745810e+03, 3.749410e+03, 3.753010e+03, 3.756610e+03, 3.760210e+03, 3.763810e+03, 3.767410e+03, 3.771010e+03, 3.774610e+03, 3.778210e+03, 3.780010e+03, 3.783610e+03, 3.787210e+03, 3.790810e+03, 3.794410e+03, 3.798010e+03, 3.801610e+03, 3.805210e+03, 3.807010e+03, 3.810610e+03, 3.814210e+03, 3.817810e+03, 3.821410e+03, 3.825010e+03, 3.828610e+03, 3.832210e+03, 3.834010e+03, 3.837610e+03, 3.841210e+03, 3.844810e+03, 3.848410e+03, 3.852010e+03, 3.855610e+03, 3.859210e+03, 3.862810e+03, 3.866410e+03, 3.870010e+03, 3.871810e+03, 3.875410e+03, 3.879010e+03, 3.882610e+03, 3.886210e+03, 3.889810e+03, 3.893410e+03, 3.897010e+03, 3.900610e+03, 3.904210e+03, 3.907810e+03, 3.909610e+03, 3.913210e+03, 3.916810e+03, 3.920410e+03, 3.924010e+03, 3.927610e+03, 3.931210e+03, 3.934810e+03, 3.938410e+03, 3.940210e+03, 3.943810e+03, 3.947410e+03, 3.951010e+03, 3.954610e+03, 3.958210e+03, 3.961810e+03, 3.965410e+03, 3.969010e+03, 3.972610e+03, 3.976210e+03, 3.979810e+03, 3.983410e+03, 3.985210e+03, 3.988810e+03, 3.992410e+03, 3.996010e+03, 3.999610e+03, 4.003210e+03, 4.006810e+03, 4.010410e+03, 4.014010e+03, 4.017610e+03, 4.021210e+03, 4.024811e+03, 4.028411e+03, 4.032011e+03, 4.033811e+03, 4.037411e+03, 4.041011e+03, 4.044611e+03, 4.048211e+03, 4.051811e+03, 4.055411e+03, 4.059011e+03, 4.062611e+03, 4.066211e+03, 4.069811e+03, 4.073411e+03, 4.077011e+03, 4.080611e+03, 4.084211e+03, 4.087811e+03, 4.091411e+03, 4.095011e+03, 4.098611e+03, 4.102211e+03, 4.105811e+03, 4.109411e+03, 4.113011e+03, 4.114811e+03, 4.118411e+03, 4.122011e+03, 4.125611e+03, 4.129211e+03, 4.132811e+03, 4.136411e+03, 4.140011e+03, 4.143611e+03, 4.147211e+03, 4.150811e+03, 4.154411e+03, 4.158011e+03, 4.161611e+03, 4.165211e+03, 4.168811e+03, 4.172411e+03, 4.176011e+03, 4.179611e+03, 4.183211e+03, 4.186811e+03, 4.190411e+03, 4.194011e+03, 4.197611e+03, 4.199411e+03, 4.203011e+03, 4.206611e+03, 4.210211e+03, 4.213811e+03, 4.217411e+03, 4.221011e+03, 4.224611e+03, 4.228211e+03, 4.231811e+03, 4.235411e+03, 4.239011e+03, 4.242611e+03, 4.246211e+03, 4.249811e+03, 4.253411e+03, 4.257011e+03, 4.260611e+03, 4.264211e+03, 4.267811e+03, 4.271411e+03, 4.275011e+03, 4.278611e+03, 4.282211e+03, 4.285811e+03, 4.289411e+03, 4.293011e+03, 4.296611e+03, 4.300211e+03, 4.303811e+03, 4.307411e+03, 4.311011e+03, 4.314611e+03, 4.318211e+03, 4.321811e+03, 4.325411e+03, 4.329011e+03, 4.332611e+03, 4.336211e+03, 4.339811e+03, 4.343411e+03, 4.347011e+03, 4.350611e+03, 4.354211e+03, 4.357811e+03, 4.361411e+03, 4.365011e+03, 4.368611e+03, 4.372211e+03, 4.375811e+03, 4.379411e+03, 4.383011e+03, 4.386611e+03, 4.390211e+03, 4.393811e+03, 4.397411e+03, 4.401011e+03, 4.404611e+03, 4.408212e+03, 4.411812e+03, 4.415412e+03, 4.419012e+03, 4.422612e+03, 4.426212e+03, 4.429812e+03, 4.433412e+03, 4.437012e+03, 4.440612e+03, 4.444212e+03, 4.447812e+03, 4.451412e+03, 4.455012e+03, 4.460412e+03, 4.462212e+03, 4.467612e+03, 4.471212e+03, 4.474812e+03, 4.478412e+03, 4.482012e+03, 4.485612e+03, 4.489212e+03, 4.492812e+03, 4.496412e+03, 4.500012e+03, 4.503612e+03, 4.507212e+03, 4.510812e+03, 4.514412e+03, 4.518012e+03, 4.521612e+03, 4.525212e+03, 4.528812e+03, 4.532412e+03, 4.536012e+03, 4.539612e+03, 4.545012e+03, 4.548612e+03, 4.552212e+03, 4.555812e+03, 4.559412e+03, 4.563012e+03, 4.566612e+03, 4.570212e+03, 4.573812e+03, 4.577412e+03, 4.581012e+03, 4.584612e+03, 4.590012e+03, 4.593612e+03, 4.597212e+03, 4.600812e+03, 4.604412e+03, 4.608012e+03, 4.611612e+03, 4.615212e+03, 4.618812e+03, 4.622412e+03, 4.626012e+03, 4.629612e+03, 4.633212e+03, 4.638612e+03, 4.642212e+03, 4.645812e+03, 4.649412e+03, 4.653012e+03, 4.656612e+03, 4.660212e+03, 4.663812e+03, 4.667412e+03, 4.672812e+03, 4.676412e+03, 4.680012e+03, 4.683612e+03, 4.687212e+03, 4.690812e+03, 4.694412e+03, 4.698012e+03, 4.703412e+03, 4.707012e+03, 4.710612e+03, 4.714212e+03, 4.717812e+03, 4.721412e+03, 4.725012e+03, 4.730412e+03, 4.734012e+03, 4.737612e+03, 4.741212e+03, 4.744812e+03, 4.748412e+03, 4.752012e+03, 4.755612e+03, 4.761012e+03, 4.764612e+03, 4.768212e+03, 4.771812e+03, 4.775412e+03, 4.780812e+03, 4.784412e+03, 4.788012e+03, 4.791613e+03, 4.795213e+03, 4.798813e+03, 4.802413e+03, 4.807813e+03, 4.811413e+03, 4.815013e+03, 4.818613e+03, 4.822213e+03, 4.825813e+03, 4.829413e+03, 4.834813e+03, 4.838413e+03, 4.842013e+03, 4.845613e+03, 4.849213e+03, 4.852813e+03, 4.858213e+03, 4.861813e+03, 4.865413e+03, 4.869013e+03, 4.872613e+03, 4.876213e+03, 4.881613e+03, 4.885213e+03, 4.888813e+03, 4.892413e+03, 4.897813e+03, 4.901413e+03, 4.905013e+03, 4.908613e+03, 4.912213e+03, 4.915813e+03, 4.921213e+03, 4.924813e+03, 4.928413e+03, 4.932013e+03, 4.937413e+03, 4.941013e+03, 4.944613e+03, 4.948213e+03, 4.951813e+03, 4.957213e+03, 4.960813e+03, 4.964413e+03, 4.968013e+03, 4.973413e+03, 4.977013e+03, 4.980613e+03, 4.984213e+03, 4.989613e+03, 4.993213e+03, 4.996813e+03, 5.000413e+03, 5.005813e+03, 5.009413e+03, 5.013013e+03, 5.016613e+03, 5.022013e+03, 5.025613e+03, 5.029213e+03, 5.032813e+03, 5.036413e+03, 5.041813e+03, 5.045413e+03, 5.049013e+03, 5.052613e+03, 5.058013e+03, 5.061613e+03, 5.065213e+03, 5.068813e+03, 5.074213e+03, 5.077813e+03, 5.081413e+03, 5.086813e+03, 5.090413e+03, 5.094013e+03, 5.099413e+03, 5.103013e+03, 5.106613e+03, 5.110213e+03, 5.113813e+03, 5.119213e+03, 5.122813e+03, 5.126413e+03, 5.130013e+03, 5.133613e+03, 5.139013e+03, 5.142613e+03, 5.146213e+03, 5.151613e+03, 5.155213e+03, 5.158813e+03, 5.162413e+03, 5.167813e+03, 5.171413e+03, 5.175014e+03, 5.178614e+03, 5.184014e+03, 5.187614e+03, 5.191214e+03, 5.194814e+03, 5.200214e+03, 5.203814e+03, 5.207414e+03, 5.211014e+03, 5.216414e+03, 5.220014e+03, 5.223614e+03, 5.229014e+03, 5.232614e+03, 5.236214e+03, 5.239814e+03, 5.245214e+03, 5.248814e+03, 5.252414e+03, 5.256014e+03, 5.261414e+03, 5.265014e+03, 5.268614e+03, 5.272214e+03, 5.277614e+03, 5.281214e+03, 5.284814e+03, 5.290214e+03, 5.293814e+03, 5.297414e+03, 5.301014e+03, 5.306414e+03, 5.310014e+03, 5.313614e+03, 5.317214e+03, 5.322614e+03, 5.326214e+03, 5.329814e+03, 5.335214e+03, 5.338814e+03, 5.342414e+03, 5.346014e+03, 5.351414e+03, 5.355014e+03, 5.358614e+03, 5.362214e+03, 5.367614e+03, 5.371214e+03, 5.374814e+03, 5.380214e+03, 5.383814e+03, 5.387414e+03, 5.391014e+03, 5.396414e+03, 5.400014e+03, 5.403614e+03, 5.409014e+03, 5.412614e+03, 5.416214e+03, 5.419814e+03, 5.425214e+03, 5.428814e+03, 5.432414e+03, 5.437814e+03, 5.441414e+03, 5.445014e+03, 5.448614e+03, 5.454014e+03, 5.457614e+03, 5.461214e+03, 5.464814e+03, 5.470214e+03, 5.473814e+03, 5.477414e+03, 5.481014e+03, 5.486414e+03, 5.490014e+03, 5.493614e+03, 5.499014e+03, 5.502614e+03, 5.506214e+03, 5.509814e+03, 5.515214e+03, 5.518814e+03, 5.522414e+03, 5.526014e+03, 5.531414e+03, 5.535014e+03, 5.538614e+03, 5.544014e+03, 5.547614e+03, 5.551214e+03, 5.554814e+03, 5.560215e+03, 5.563815e+03, 5.567415e+03, 5.571015e+03, 5.576415e+03, 5.580015e+03, 5.583615e+03, 5.589015e+03, 5.592615e+03, 5.596215e+03, 5.601615e+03, 5.605215e+03, 5.608815e+03, 5.612415e+03, 5.617815e+03, 5.621415e+03, 5.625015e+03, 5.630415e+03, 5.634015e+03, 5.637615e+03, 5.641215e+03, 5.646615e+03, 5.650215e+03, 5.653815e+03, 5.659215e+03, 5.662815e+03, 5.666415e+03, 5.670015e+03, 5.675415e+03, 5.679015e+03, 5.682615e+03, 5.688015e+03, 5.691615e+03, 5.695215e+03, 5.700615e+03, 5.704215e+03, 5.707815e+03, 5.713215e+03, 5.716815e+03, 5.720415e+03, 5.725815e+03, 5.729415e+03, 5.733015e+03, 5.736615e+03, 5.742015e+03, 5.745615e+03, 5.749215e+03, 5.754615e+03, 5.758215e+03, 5.761815e+03, 5.767215e+03, 5.770815e+03, 5.774415e+03, 5.778015e+03, 5.783415e+03, 5.787015e+03, 5.790615e+03, 5.796015e+03, 5.799615e+03, 5.803215e+03, 5.808615e+03, 5.812215e+03, 5.815815e+03, 5.821215e+03, 5.824815e+03, 5.828415e+03, 5.833815e+03, 5.837415e+03, 5.841015e+03, 5.846415e+03, 5.850015e+03, 5.853615e+03, 5.859015e+03, 5.862615e+03, 5.866215e+03, 5.871615e+03, 5.875215e+03, 5.878815e+03, 5.884215e+03, 5.887815e+03, 5.891415e+03, 5.896815e+03, 5.900415e+03, 5.904015e+03, 5.909415e+03, 5.913015e+03, 5.916615e+03, 5.922015e+03, 5.925615e+03, 5.929215e+03, 5.934615e+03, 5.938216e+03, 5.943616e+03, 5.947216e+03, 5.950816e+03, 5.956216e+03, 5.959816e+03, 5.963416e+03, 5.968816e+03, 5.972416e+03, 5.976016e+03, 5.981416e+03, 5.985016e+03, 5.988616e+03, 5.994016e+03, 5.997616e+03, 6.003016e+03, 6.006616e+03, 6.010216e+03, 6.015616e+03, 6.019216e+03, 6.024616e+03, 6.028216e+03, 6.031816e+03, 6.037216e+03, 6.040816e+03, 6.044416e+03, 6.049816e+03, 6.053416e+03, 6.058816e+03, 6.062416e+03, 6.066016e+03, 6.071416e+03, 6.075016e+03, 6.078616e+03, 6.084016e+03, 6.087616e+03, 6.093016e+03, 6.096616e+03, 6.100216e+03, 6.105616e+03, 6.109216e+03, 6.112816e+03, 6.118216e+03, 6.121816e+03, 6.125416e+03, 6.130816e+03, 6.134416e+03, 6.139816e+03, 6.143416e+03, 6.147016e+03, 6.152416e+03, 6.156016e+03, 6.159616e+03, 6.165016e+03, 6.168616e+03, 6.174016e+03, 6.177616e+03, 6.181216e+03, 6.186616e+03, 6.190216e+03, 6.193816e+03, 6.199216e+03, 6.202816e+03, 6.208216e+03, 6.211816e+03, 6.215416e+03, 6.220816e+03, 6.224416e+03, 6.229816e+03, 6.233416e+03, 6.238816e+03, 6.242416e+03, 6.246016e+03, 6.251416e+03, 6.255016e+03, 6.260416e+03, 6.264016e+03, 6.267616e+03, 6.273016e+03, 6.276616e+03, 6.282016e+03, 6.285616e+03, 6.289216e+03, 6.294616e+03, 6.298216e+03, 6.303616e+03, 6.307216e+03, 6.310816e+03, 6.316216e+03, 6.319816e+03, 6.325217e+03, 6.328817e+03, 6.334217e+03, 6.337817e+03, 6.341417e+03, 6.346817e+03, 6.350417e+03, 6.355817e+03, 6.359417e+03, 6.364817e+03, 6.368417e+03, 6.372017e+03, 6.377417e+03, 6.381017e+03, 6.386417e+03, 6.390017e+03, 6.395417e+03, 6.399017e+03, 6.404417e+03, 6.408017e+03, 6.413417e+03, 6.417017e+03, 6.420617e+03, 6.426017e+03, 6.431417e+03, 6.435017e+03, 6.438617e+03, 6.444017e+03, 6.447617e+03, 6.453017e+03, 6.456617e+03, 6.462017e+03, 6.465617e+03, 6.469217e+03, 6.474617e+03, 6.478217e+03, 6.483617e+03, 6.487217e+03, 6.492617e+03, 6.496217e+03, 6.501617e+03, 6.505217e+03, 6.510617e+03, 6.514217e+03, 6.519617e+03, 6.523217e+03, 6.528617e+03, 6.532217e+03, 6.537617e+03, 6.541217e+03, 6.546617e+03, 6.550217e+03, 6.555617e+03, 6.559217e+03, 6.564617e+03, 6.568217e+03, 6.573617e+03, 6.577217e+03, 6.582617e+03, 6.586217e+03, 6.591617e+03, 6.595217e+03, 6.600617e+03, 6.604217e+03, 6.609617e+03, 6.613217e+03, 6.618617e+03, 6.622217e+03, 6.627617e+03, 6.631217e+03, 6.636617e+03, 6.640217e+03, 6.645617e+03, 6.649217e+03, 6.654617e+03, 6.658217e+03, 6.663617e+03, 6.667217e+03, 6.672617e+03, 6.676217e+03, 6.681617e+03, 6.687017e+03, 6.690617e+03, 6.696017e+03, 6.699617e+03, 6.705018e+03, 6.708618e+03, 6.714018e+03, 6.717618e+03, 6.723018e+03, 6.726618e+03, 6.732018e+03, 6.735618e+03, 6.741018e+03, 6.744618e+03, 6.750018e+03, 6.753618e+03, 6.757218e+03, 6.762618e+03, 6.768018e+03, 6.771618e+03, 6.777018e+03, 6.780618e+03, 6.786018e+03, 6.789618e+03, 6.795018e+03, 6.798618e+03, 6.804018e+03, 6.809418e+03, 6.813018e+03, 6.818418e+03, 6.822018e+03, 6.827418e+03, 6.831018e+03, 6.836418e+03, 6.840018e+03, 6.845418e+03, 6.849018e+03, 6.854418e+03, 6.858018e+03, 6.863418e+03, 6.867018e+03, 6.872418e+03, 6.877818e+03, 6.881418e+03, 6.886818e+03, 6.892218e+03, 6.895818e+03, 6.901218e+03, 6.904818e+03, 6.910218e+03, 6.913818e+03, 6.919218e+03, 6.922818e+03, 6.928218e+03, 6.931818e+03, 6.937218e+03, 6.942618e+03, 6.946218e+03, 6.951618e+03, 6.955218e+03, 6.960618e+03, 6.966018e+03, 6.969618e+03, 6.975018e+03, 6.978618e+03, 6.984018e+03, 6.989418e+03, 6.993018e+03, 6.998418e+03, 7.002018e+03, 7.007418e+03, 7.012818e+03, 7.016418e+03, 7.021818e+03, 7.025418e+03, 7.030818e+03, 7.036218e+03, 7.039818e+03, 7.045218e+03, 7.048818e+03, 7.054218e+03, 7.057818e+03, 7.063218e+03, 7.068618e+03, 7.072218e+03, 7.077618e+03, 7.081218e+03, 7.086618e+03, 7.092019e+03, 7.095619e+03, 7.101019e+03, 7.104619e+03, 7.110019e+03, 7.115419e+03, 7.119019e+03, 7.124419e+03, 7.129819e+03, 7.133419e+03, 7.138819e+03, 7.142419e+03, 7.147819e+03, 7.153219e+03, 7.156819e+03, 7.162219e+03, 7.167619e+03, 7.171219e+03, 7.176619e+03, 7.180219e+03, 7.185619e+03, 7.191019e+03, 7.194619e+03, 7.200019e+03, 7.205419e+03, 7.209019e+03, 7.214419e+03, 7.219819e+03, 7.223419e+03, 7.228819e+03, 7.234219e+03, 7.237819e+03, 7.243219e+03, 7.248619e+03, 7.252219e+03, 7.257619e+03, 7.261219e+03, 7.266619e+03, 7.272019e+03, 7.275619e+03, 7.281019e+03, 7.286419e+03, 7.290019e+03, 7.295419e+03, 7.300819e+03, 7.304419e+03, 7.309819e+03, 7.315219e+03, 7.320619e+03, 7.324219e+03, 7.329619e+03, 7.335019e+03, 7.338619e+03, 7.344019e+03, 7.347619e+03, 7.353019e+03, 7.358419e+03, 7.362019e+03, 7.367419e+03, 7.372819e+03, 7.376419e+03, 7.381819e+03, 7.387219e+03, 7.390819e+03, 7.396219e+03, 7.401619e+03, 7.405219e+03, 7.410619e+03, 7.416019e+03, 7.419619e+03, 7.425019e+03, 7.430419e+03, 7.434019e+03, 7.439419e+03, 7.444819e+03, 7.448419e+03, 7.453819e+03, 7.459219e+03, 7.462819e+03, 7.468219e+03, 7.473620e+03, 7.477220e+03, 7.482620e+03, 7.488020e+03, 7.493420e+03, 7.497020e+03, 7.502420e+03, 7.507820e+03, 7.511420e+03, 7.516820e+03, 7.522220e+03, 7.525820e+03, 7.531220e+03, 7.536620e+03, 7.540220e+03, 7.545620e+03, 7.551020e+03, 7.556420e+03, 7.561820e+03, 7.565420e+03, 7.570820e+03, 7.576220e+03, 7.579820e+03, 7.585220e+03, 7.590620e+03, 7.596020e+03, 7.599620e+03, 7.605020e+03, 7.608620e+03, 7.614020e+03, 7.619420e+03, 7.624820e+03, 7.628420e+03, 7.633820e+03, 7.639220e+03, 7.644620e+03, 7.650020e+03, 7.653620e+03, 7.659020e+03, 7.664420e+03, 7.668020e+03, 7.673420e+03, 7.678820e+03, 7.684220e+03, 7.689620e+03, 7.693220e+03, 7.698620e+03, 7.704020e+03, 7.709420e+03, 7.713020e+03, 7.718420e+03, 7.723820e+03, 7.729220e+03, 7.732820e+03, 7.738220e+03, 7.743620e+03, 7.747220e+03, 7.752620e+03, 7.758020e+03, 7.763420e+03, 7.767020e+03, 7.772420e+03, 7.777820e+03, 7.783220e+03, 7.786820e+03, 7.792220e+03, 7.797620e+03, 7.803020e+03, 7.806620e+03, 7.812020e+03, 7.817420e+03, 7.822820e+03, 7.828220e+03, 7.831820e+03, 7.837220e+03, 7.842620e+03, 7.846220e+03, 7.851620e+03, 7.857021e+03, 7.862421e+03, 7.867821e+03, 7.871421e+03, 7.876821e+03, 7.882221e+03, 7.887621e+03, 7.893021e+03, 7.896621e+03, 7.902021e+03, 7.907421e+03, 7.912821e+03, 7.918221e+03, 7.921821e+03, 7.927221e+03, 7.932621e+03, 7.938021e+03, 7.943421e+03, 7.947021e+03, 7.952421e+03, 7.957821e+03, 7.963221e+03, 7.968621e+03, 7.972221e+03, 7.977621e+03, 7.983021e+03, 7.988421e+03, 7.993821e+03, 7.999221e+03, 8.002821e+03, 8.008221e+03, 8.013621e+03, 8.019021e+03, 8.024421e+03, 8.028021e+03, 8.033421e+03, 8.038821e+03, 8.044221e+03, 8.049621e+03, 8.055021e+03, 8.058621e+03, 8.064021e+03, 8.069421e+03, 8.074821e+03, 8.080221e+03, 8.085621e+03, 8.089221e+03, 8.094621e+03, 8.100021e+03, 8.105421e+03, 8.110821e+03, 8.116221e+03, 8.119821e+03, 8.125221e+03, 8.130621e+03, 8.136021e+03, 8.141421e+03, 8.145021e+03, 8.150421e+03, 8.155821e+03, 8.161221e+03, 8.166621e+03, 8.172021e+03, 8.175621e+03, 8.181021e+03, 8.186421e+03, 8.191821e+03, 8.197221e+03, 8.202621e+03, 8.208021e+03, 8.213421e+03, 8.218821e+03, 8.224221e+03, 8.227821e+03, 8.233221e+03, 8.238622e+03, 8.244022e+03, 8.249422e+03, 8.254822e+03, 8.260222e+03, 8.265622e+03, 8.271022e+03, 8.276422e+03, 8.281822e+03, 8.287222e+03, 8.290822e+03, 8.296222e+03, 8.301622e+03, 8.307022e+03, 8.312422e+03, 8.317822e+03, 8.323222e+03, 8.328622e+03, 8.334022e+03, 8.339422e+03, 8.344822e+03, 8.348422e+03, 8.353822e+03, 8.359222e+03, 8.364622e+03, 8.370022e+03, 8.375422e+03, 8.380822e+03, 8.386222e+03, 8.391622e+03, 8.397022e+03, 8.400622e+03, 8.406022e+03, 8.411422e+03, 8.416822e+03, 8.422222e+03, 8.427622e+03, 8.433022e+03, 8.438422e+03, 8.443822e+03, 8.449222e+03, 8.454622e+03, 8.460022e+03, 8.465422e+03, 8.470822e+03, 8.476222e+03, 8.481622e+03, 8.487022e+03, 8.490622e+03, 8.496022e+03, 8.501422e+03, 8.506822e+03, 8.512222e+03, 8.517622e+03, 8.523022e+03, 8.528422e+03, 8.533822e+03, 8.539222e+03, 8.544622e+03, 8.550022e+03, 8.555422e+03, 8.560822e+03, 8.566222e+03, 8.569822e+03, 8.575222e+03, 8.580622e+03, 8.586022e+03, 8.591422e+03, 8.596822e+03, 8.602222e+03, 8.607622e+03, 8.613022e+03, 8.618422e+03, 8.623823e+03, 8.629223e+03, 8.634623e+03, 8.640023e+03, 8.645423e+03, 8.650823e+03, 8.656223e+03, 8.661623e+03, 8.667023e+03, 8.672423e+03, 8.677823e+03, 8.683223e+03, 8.688623e+03, 8.694023e+03, 8.699423e+03, 8.704823e+03, 8.710223e+03, 8.715623e+03, 8.721023e+03, 8.726423e+03, 8.731823e+03, 8.737223e+03, 8.742623e+03, 8.748023e+03, 8.753423e+03, 8.758823e+03, 8.764223e+03, 8.769623e+03, 8.775023e+03, 8.780423e+03, 8.785823e+03, 8.791223e+03, 8.796623e+03, 8.802023e+03, 8.807423e+03, 8.812823e+03, 8.818223e+03, 8.823623e+03, 8.829023e+03, 8.834423e+03, 8.839823e+03, 8.845223e+03, 8.850623e+03, 8.856023e+03, 8.861423e+03, 8.866823e+03, 8.872223e+03, 8.877623e+03, 8.883023e+03, 8.888423e+03, 8.893823e+03, 8.901023e+03, 8.906423e+03, 8.911823e+03, 8.917223e+03, 8.922623e+03, 8.928023e+03, 8.933423e+03, 8.938823e+03, 8.944223e+03, 8.949623e+03, 8.955023e+03, 8.962223e+03, 8.967623e+03, 8.973023e+03, 8.978423e+03, 8.983823e+03, 8.989223e+03, 8.994623e+03, 9.000023e+03, 9.005424e+03, 9.010824e+03, 9.016224e+03, 9.021624e+03, 9.028824e+03, 9.034224e+03, 9.039624e+03, 9.045024e+03, 9.050424e+03, 9.055824e+03, 9.061224e+03, 9.066624e+03, 9.072024e+03, 9.079224e+03, 9.084624e+03, 9.090024e+03, 9.095424e+03, 9.100824e+03, 9.106224e+03, 9.111624e+03, 9.117024e+03, 9.124224e+03, 9.129624e+03, 9.135024e+03, 9.140424e+03, 9.145824e+03, 9.151224e+03, 9.156624e+03, 9.163824e+03, 9.169224e+03, 9.174624e+03, 9.180024e+03, 9.185424e+03, 9.190824e+03, 9.196224e+03, 9.201624e+03, 9.207024e+03, 9.212424e+03, 9.217824e+03, 9.225024e+03, 9.230424e+03, 9.235824e+03, 9.241224e+03, 9.248424e+03, 9.253824e+03, 9.259224e+03, 9.264624e+03, 9.270024e+03, 9.275424e+03, 9.282624e+03, 9.288024e+03, 9.293424e+03, 9.298824e+03, 9.304224e+03, 9.309624e+03, 9.315024e+03, 9.320424e+03, 9.327624e+03, 9.333024e+03, 9.338424e+03, 9.343824e+03, 9.351024e+03, 9.356424e+03, 9.361824e+03, 9.367224e+03, 9.372624e+03, 9.378024e+03, 9.383424e+03, 9.390625e+03, 9.396025e+03, 9.401425e+03, 9.406825e+03, 9.412225e+03, 9.417625e+03, 9.424825e+03, 9.430225e+03, 9.435625e+03, 9.441025e+03, 9.448225e+03, 9.453625e+03, 9.459025e+03, 9.464425e+03, 9.469825e+03, 9.475225e+03, 9.482425e+03, 9.487825e+03, 9.493225e+03, 9.498625e+03, 9.504025e+03, 9.511225e+03, 9.516625e+03, 9.522025e+03, 9.527425e+03, 9.534625e+03, 9.540025e+03, 9.545425e+03, 9.550825e+03, 9.556225e+03, 9.563425e+03, 9.568825e+03, 9.574225e+03, 9.579625e+03, 9.586825e+03, 9.592225e+03, 9.597625e+03, 9.603025e+03, 9.608425e+03, 9.613825e+03, 9.621025e+03, 9.626425e+03, 9.631825e+03, 9.637225e+03, 9.644425e+03, 9.649825e+03, 9.655225e+03, 9.660625e+03, 9.666025e+03, 9.673225e+03, 9.678625e+03, 9.684025e+03, 9.689425e+03, 9.694825e+03, 9.702025e+03, 9.707425e+03, 9.712825e+03, 9.720025e+03, 9.725425e+03, 9.730825e+03, 9.738025e+03, 9.743425e+03, 9.748825e+03, 9.754225e+03, 9.761425e+03, 9.766825e+03, 9.772226e+03, 9.777626e+03, 9.784826e+03, 9.790226e+03, 9.795626e+03, 9.802826e+03, 9.808226e+03, 9.813626e+03, 9.819026e+03, 9.826226e+03, 9.831626e+03, 9.837026e+03, 9.844226e+03, 9.849626e+03, 9.856826e+03, 9.862226e+03, 9.867626e+03, 9.873026e+03, 9.880226e+03, 9.885626e+03, 9.892826e+03, 9.898226e+03, 9.903626e+03, 9.909026e+03, 9.916226e+03, 9.921626e+03, 9.927026e+03, 9.934226e+03, 9.939626e+03, 9.946826e+03, 9.952226e+03, 9.957626e+03, 9.964826e+03, 9.970226e+03, 9.977426e+03, 9.982826e+03, 9.988226e+03, 9.995426e+03, 1.000083e+04, 1.000623e+04, 1.001343e+04, 1.001883e+04, 1.002423e+04, 1.003143e+04, 1.003683e+04, 1.004223e+04, 1.004763e+04, 1.005483e+04, 1.006023e+04, 1.006743e+04, 1.007283e+04, 1.007823e+04, 1.008543e+04, 1.009083e+04, 1.009803e+04, 1.010343e+04, 1.011063e+04, 1.011603e+04, 1.012143e+04, 1.012863e+04, 1.013403e+04, 1.014123e+04, 1.014663e+04, 1.015203e+04, 1.015923e+04, 1.016463e+04, 1.017003e+04, 1.017723e+04, 1.018263e+04, 1.018803e+04, 1.019523e+04, 1.020063e+04, 1.020783e+04, 1.021323e+04, 1.021863e+04, 1.022583e+04, 1.023123e+04, 1.023843e+04, 1.024383e+04, 1.025103e+04, 1.025643e+04, 1.026363e+04, 1.026903e+04, 1.027443e+04, 1.028163e+04, 1.028703e+04, 1.029423e+04, 1.029963e+04, 1.030683e+04, 1.031223e+04, 1.031763e+04, 1.032483e+04, 1.033023e+04, 1.033743e+04, 1.034283e+04, 1.035003e+04, 1.035543e+04, 1.036083e+04, 1.036803e+04, 1.037343e+04, 1.038063e+04, 1.038603e+04, 1.039323e+04, 1.039863e+04, 1.040583e+04, 1.041123e+04, 1.041843e+04, 1.042383e+04, 1.043103e+04, 1.043643e+04, 1.044363e+04, 1.044903e+04, 1.045623e+04, 1.046163e+04, 1.046703e+04, 1.047423e+04, 1.047963e+04, 1.048683e+04, 1.049223e+04, 1.049943e+04, 1.050483e+04, 1.051203e+04, 1.051923e+04, 1.052463e+04, 1.053003e+04, 1.053723e+04, 1.054443e+04, 1.054983e+04, 1.055703e+04, 1.056243e+04, 1.056783e+04, 1.057503e+04, 1.058043e+04, 1.058763e+04, 1.059483e+04, 1.060023e+04, 1.060563e+04, 1.061283e+04, 1.061823e+04, 1.062543e+04, 1.063263e+04, 1.063803e+04, 1.064523e+04, 1.065063e+04, 1.065783e+04, 1.066503e+04, 1.067043e+04, 1.067763e+04, 1.068303e+04, 1.069023e+04, 1.069563e+04, 1.070283e+04, 1.070823e+04, 1.071543e+04, 1.072083e+04, 1.072803e+04, 1.073523e+04, 1.074063e+04, 1.074783e+04, 1.075323e+04, 1.076043e+04, 1.076763e+04, 1.077303e+04, 1.078023e+04, 1.078563e+04, 1.079283e+04, 1.079823e+04, 1.080543e+04, 1.081263e+04, 1.081803e+04, 1.082523e+04, 1.083063e+04, 1.083783e+04, 1.084323e+04, 1.085043e+04, 1.085583e+04, 1.086303e+04, 1.087023e+04, 1.087563e+04, 1.088283e+04, 1.088823e+04, 1.089543e+04, 1.090263e+04, 1.090803e+04, 1.091523e+04, 1.092063e+04, 1.092783e+04, 1.093323e+04, 1.094043e+04, 1.094583e+04, 1.095303e+04, 1.095843e+04, 1.096563e+04, 1.097103e+04, 1.097823e+04, 1.098543e+04, 1.099083e+04, 1.099803e+04, 1.100343e+04, 1.101063e+04, 1.101783e+04, 1.102323e+04, 1.103043e+04, 1.103763e+04, 1.104303e+04, 1.105023e+04, 1.105563e+04, 1.106283e+04, 1.107003e+04, 1.107723e+04, 1.108263e+04, 1.108983e+04, 1.109703e+04, 1.110243e+04, 1.110963e+04, 1.111683e+04, 1.112403e+04, 1.112943e+04, 1.113663e+04, 1.114383e+04, 1.114923e+04, 1.115643e+04, 1.116363e+04, 1.116903e+04, 1.117623e+04, 1.118343e+04, 1.119063e+04, 1.119603e+04, 1.120323e+04, 1.121043e+04, 1.121763e+04, 1.122303e+04, 1.123023e+04, 1.123563e+04, 1.124283e+04, 1.125003e+04, 1.125543e+04, 1.126263e+04, 1.126983e+04, 1.127523e+04, 1.128243e+04, 1.128963e+04, 1.129503e+04, 1.130223e+04, 1.130943e+04, 1.131663e+04, 1.132203e+04, 1.132923e+04, 1.133643e+04, 1.134363e+04, 1.135083e+04, 1.135623e+04, 1.136343e+04, 1.137063e+04, 1.137783e+04, 1.138323e+04, 1.139043e+04, 1.139763e+04, 1.140483e+04, 1.141023e+04, 1.141743e+04, 1.142463e+04, 1.143003e+04, 1.143723e+04, 1.144443e+04, 1.144983e+04, 1.145703e+04, 1.146423e+04, 1.147143e+04, 1.147863e+04, 1.148583e+04, 1.149123e+04, 1.149843e+04, 1.150563e+04, 1.151283e+04, 1.152003e+04, 1.152543e+04, 1.153263e+04, 1.153983e+04, 1.154703e+04, 1.155423e+04, 1.155963e+04, 1.156683e+04, 1.157403e+04, 1.158123e+04, 1.158663e+04, 1.159383e+04, 1.160103e+04, 1.160823e+04, 1.161363e+04, 1.162083e+04, 1.162803e+04, 1.163523e+04, 1.164243e+04, 1.164963e+04, 1.165683e+04, 1.166403e+04, 1.167123e+04, 1.167663e+04, 1.168383e+04, 1.169103e+04, 1.169823e+04, 1.170543e+04, 1.171263e+04, 1.171983e+04, 1.172523e+04, 1.173243e+04, 1.173963e+04, 1.174683e+04, 1.175403e+04, 1.175943e+04, 1.176663e+04, 1.177383e+04, 1.178103e+04, 1.178823e+04, 1.179543e+04, 1.180263e+04, 1.180803e+04, 1.181523e+04, 1.182243e+04, 1.182963e+04, 1.183683e+04, 1.184223e+04, 1.184943e+04, 1.185663e+04, 1.186383e+04, 1.187103e+04, 1.187823e+04, 1.188543e+04, 1.189263e+04, 1.189803e+04, 1.190523e+04, 1.191243e+04, 1.191963e+04, 1.192683e+04, 1.193403e+04, 1.194123e+04, 1.194843e+04, 1.195563e+04, 1.196283e+04, 1.196823e+04, 1.197543e+04, 1.198263e+04, 1.198983e+04, 1.199703e+04, 1.200423e+04, 1.201143e+04, 1.201863e+04, 1.202403e+04, 1.203123e+04, 1.203843e+04, 1.204563e+04, 1.205283e+04, 1.206003e+04, 1.206723e+04, 1.207443e+04, 1.207983e+04, 1.208703e+04, 1.209423e+04, 1.210143e+04, 1.210863e+04, 1.211583e+04, 1.212303e+04, 1.213023e+04, 1.213743e+04, 1.214463e+04, 1.215003e+04, 1.215723e+04, 1.216443e+04, 1.217163e+04, 1.217883e+04, 1.218603e+04, 1.219323e+04, 1.220043e+04, 1.220763e+04, 1.221483e+04, 1.222203e+04, 1.222923e+04, 1.223643e+04, 1.224363e+04, 1.225083e+04, 1.225803e+04, 1.226523e+04, 1.227243e+04, 1.227963e+04, 1.228683e+04, 1.229403e+04, 1.230123e+04, 1.230843e+04, 1.231563e+04, 1.232283e+04, 1.233003e+04, 1.233723e+04, 1.234443e+04, 1.235163e+04, 1.235883e+04, 1.236603e+04, 1.237323e+04, 1.238043e+04, 1.238763e+04, 1.239483e+04, 1.240203e+04, 1.240923e+04, 1.241643e+04, 1.242363e+04, 1.243083e+04, 1.243803e+04, 1.244523e+04, 1.245243e+04, 1.245963e+04, 1.246683e+04, 1.247583e+04, 1.248303e+04, 1.249023e+04, 1.249743e+04, 1.250463e+04, 1.251183e+04, 1.251903e+04, 1.252623e+04, 1.253343e+04, 1.254063e+04, 1.254783e+04, 1.255503e+04, 1.256223e+04, 1.256943e+04, 1.257843e+04, 1.258563e+04, 1.259283e+04, 1.260003e+04, 1.260723e+04, 1.261443e+04, 1.262163e+04, 1.262883e+04, 1.263603e+04, 1.264323e+04, 1.265043e+04, 1.265763e+04, 1.266663e+04, 1.267383e+04, 1.268103e+04, 1.268823e+04, 1.269543e+04, 1.270263e+04, 1.270983e+04, 1.271703e+04, 1.272603e+04, 1.273323e+04, 1.274043e+04, 1.274763e+04, 1.275483e+04, 1.276203e+04, 1.276923e+04, 1.277643e+04, 1.278363e+04, 1.279263e+04, 1.279983e+04, 1.280703e+04, 1.281423e+04, 1.282143e+04, 1.282863e+04, 1.283763e+04, 1.284483e+04, 1.285203e+04, 1.285923e+04, 1.286643e+04, 1.287363e+04, 1.288083e+04, 1.288983e+04, 1.289703e+04, 1.290423e+04, 1.291143e+04, 1.291863e+04, 1.292763e+04, 1.293483e+04, 1.294203e+04, 1.294923e+04, 1.295643e+04, 1.296543e+04, 1.297263e+04, 1.297983e+04, 1.298703e+04, 1.299423e+04, 1.300143e+04, 1.300863e+04, 1.301583e+04, 1.302483e+04, 1.303203e+04, 1.303923e+04, 1.304643e+04, 1.305363e+04, 1.306083e+04, 1.306803e+04, 1.307703e+04, 1.308423e+04, 1.309143e+04, 1.309863e+04, 1.310583e+04, 1.311483e+04, 1.312203e+04, 1.312923e+04, 1.313823e+04, 1.314543e+04, 1.315263e+04, 1.315983e+04, 1.316703e+04, 1.317603e+04, 1.318323e+04, 1.319043e+04, 1.319943e+04, 1.320663e+04, 1.321383e+04, 1.322103e+04, 1.323003e+04, 1.323723e+04, 1.324443e+04, 1.325343e+04, 1.326063e+04, 1.326963e+04, 1.327683e+04, 1.328403e+04, 1.329303e+04, 1.330023e+04, 1.330743e+04, 1.331463e+04, 1.332363e+04, 1.333083e+04, 1.333803e+04, 1.334523e+04, 1.335423e+04, 1.336143e+04, 1.337043e+04, 1.337763e+04, 1.338483e+04, 1.339203e+04, 1.340103e+04, 1.340823e+04, 1.341724e+04, 1.342444e+04, 1.343164e+04, 1.344064e+04, 1.344784e+04, 1.345504e+04, 1.346404e+04, 1.347124e+04, 1.347844e+04, 1.348744e+04, 1.349464e+04, 1.350364e+04, 1.351084e+04, 1.351804e+04, 1.352704e+04, 1.353424e+04, 1.354324e+04, 1.355044e+04, 1.355764e+04, 1.356664e+04, 1.357384e+04, 1.358284e+04, 1.359004e+04, 1.359904e+04, 1.360624e+04, 1.361524e+04, 1.362244e+04, 1.362964e+04, 1.363864e+04, 1.364584e+04, 1.365484e+04, 1.366204e+04, 1.366924e+04, 1.367824e+04, 1.368544e+04, 1.369444e+04, 1.370344e+04, 1.371064e+04, 1.371784e+04, 1.372684e+04, 1.373404e+04, 1.374304e+04, 1.375024e+04, 1.375924e+04, 1.376824e+04, 1.377544e+04, 1.378444e+04, 1.379164e+04, 1.380064e+04, 1.380784e+04, 1.381684e+04, 1.382404e+04, 1.383304e+04, 1.384024e+04, 1.384924e+04, 1.385644e+04, 1.386544e+04, 1.387444e+04, 1.388164e+04, 1.389064e+04, 1.389784e+04, 1.390684e+04, 1.391584e+04, 1.392304e+04, 1.393204e+04, 1.394104e+04, 1.394824e+04, 1.395724e+04, 1.396624e+04, 1.397524e+04, 1.398244e+04, 1.399144e+04, 1.399864e+04, 1.400764e+04, 1.401484e+04, 1.402384e+04, 1.403284e+04, 1.404004e+04, 1.404904e+04, 1.405804e+04, 1.406704e+04, 1.407424e+04, 1.408324e+04, 1.409224e+04, 1.409944e+04, 1.410844e+04, 1.411744e+04, 1.412644e+04, 1.413364e+04, 1.414444e+04, 1.415164e+04, 1.416064e+04, 1.416964e+04, 1.417864e+04, 1.418764e+04, 1.419484e+04, 1.420384e+04, 1.421284e+04, 1.422184e+04, 1.423084e+04, 1.423984e+04, 1.424884e+04, 1.425784e+04, 1.426684e+04, 1.427584e+04, 1.428484e+04, 1.429384e+04, 1.430464e+04, 1.431364e+04, 1.432264e+04, 1.433164e+04, 1.434064e+04, 1.435144e+04, 1.436044e+04, 1.436944e+04, 1.437844e+04, 1.438924e+04, 1.439824e+04, 1.440904e+04, 1.441804e+04, 1.442884e+04, 1.443784e+04, 1.444684e+04, 1.445764e+04, 1.446664e+04, 1.447744e+04, 1.448824e+04, 1.449724e+04, 1.450804e+04, 1.451884e+04, 1.452784e+04, 1.453864e+04, 1.454944e+04, 1.456024e+04, 1.457284e+04, 1.458364e+04, 1.459444e+04, 1.460524e+04, 1.461784e+04, 1.462864e+04, 1.464124e+04, 1.465204e+04, 1.466464e+04, 1.467724e+04, 1.468984e+04, 1.470244e+04, 1.471324e+04, 1.472764e+04, 1.474024e+04, 1.475464e+04, 1.476724e+04, 1.478164e+04, 1.479604e+04, 1.481044e+04, 1.482484e+04, 1.483924e+04, 1.485544e+04, 1.487164e+04, 1.488784e+04, 1.490404e+04, 1.492024e+04, 1.493824e+04, 1.495444e+04, 1.497424e+04, 1.499404e+04, 1.501384e+04, 1.503544e+04, 1.505884e+04, 1.508224e+04, 1.510744e+04, 1.513444e+04, 1.516324e+04, 1.519384e+04, 1.522804e+04, 1.526404e+04, 1.530544e+04, 1.535404e+04, 1.541344e+04, 1.549444e+04, 1.561144e+04};
      Double_t yd[nbins] = {3.409000e+01, 3.408000e+01, 3.407000e+01, 3.406000e+01, 3.405000e+01, 3.404000e+01, 3.403000e+01, 3.402000e+01, 3.401000e+01, 3.400000e+01, 3.399000e+01, 3.398000e+01, 3.397000e+01, 3.396000e+01, 3.395000e+01, 3.394000e+01, 3.393000e+01, 3.392000e+01, 3.391000e+01, 3.390000e+01, 3.389000e+01, 3.388000e+01, 3.387000e+01, 3.386000e+01, 3.385000e+01, 3.384000e+01, 3.383000e+01, 3.382000e+01, 3.381000e+01, 3.380000e+01, 3.379000e+01, 3.378000e+01, 3.377000e+01, 3.376000e+01, 3.375000e+01, 3.374000e+01, 3.373000e+01, 3.372000e+01, 3.371000e+01, 3.370000e+01, 3.369000e+01, 3.368000e+01, 3.367000e+01, 3.366000e+01, 3.365000e+01, 3.364000e+01, 3.363000e+01, 3.362000e+01, 3.361000e+01, 3.360000e+01, 3.359000e+01, 3.358000e+01, 3.357000e+01, 3.356000e+01, 3.355000e+01, 3.354000e+01, 3.353000e+01, 3.352000e+01, 3.351000e+01, 3.350000e+01, 3.349000e+01, 3.348000e+01, 3.347000e+01, 3.346000e+01, 3.345000e+01, 3.344000e+01, 3.343000e+01, 3.342000e+01, 3.341000e+01, 3.340000e+01, 3.339000e+01, 3.338000e+01, 3.337000e+01, 3.336000e+01, 3.335000e+01, 3.334000e+01, 3.333000e+01, 3.332000e+01, 3.331000e+01, 3.330000e+01, 3.329000e+01, 3.328000e+01, 3.327000e+01, 3.326000e+01, 3.325000e+01, 3.324000e+01, 3.323000e+01, 3.322000e+01, 3.321000e+01, 3.320000e+01, 3.319000e+01, 3.318000e+01, 3.317000e+01, 3.316000e+01, 3.315000e+01, 3.314000e+01, 3.313000e+01, 3.312000e+01, 3.311000e+01, 3.310000e+01, 3.309000e+01, 3.308000e+01, 3.307000e+01, 3.306000e+01, 3.305000e+01, 3.304000e+01, 3.303000e+01, 3.302000e+01, 3.301000e+01, 3.300000e+01, 3.299000e+01, 3.298000e+01, 3.297000e+01, 3.296000e+01, 3.295000e+01, 3.294000e+01, 3.293000e+01, 3.292000e+01, 3.291000e+01, 3.290000e+01, 3.289000e+01, 3.288000e+01, 3.287000e+01, 3.286000e+01, 3.285000e+01, 3.284000e+01, 3.283000e+01, 3.282000e+01, 3.281000e+01, 3.280000e+01, 3.279000e+01, 3.278000e+01, 3.277000e+01, 3.276000e+01, 3.275000e+01, 3.274000e+01, 3.273000e+01, 3.272000e+01, 3.271000e+01, 3.270000e+01, 3.269000e+01, 3.268000e+01, 3.267000e+01, 3.266000e+01, 3.265000e+01, 3.264000e+01, 3.263000e+01, 3.262000e+01, 3.261000e+01, 3.260000e+01, 3.259000e+01, 3.258000e+01, 3.257000e+01, 3.256000e+01, 3.255000e+01, 3.254000e+01, 3.253000e+01, 3.252000e+01, 3.251000e+01, 3.250000e+01, 3.249000e+01, 3.248000e+01, 3.247000e+01, 3.246000e+01, 3.245000e+01, 3.244000e+01, 3.243000e+01, 3.242000e+01, 3.241000e+01, 3.240000e+01, 3.239000e+01, 3.238000e+01, 3.237000e+01, 3.236000e+01, 3.235000e+01, 3.234000e+01, 3.233000e+01, 3.232000e+01, 3.231000e+01, 3.230000e+01, 3.229000e+01, 3.228000e+01, 3.227000e+01, 3.226000e+01, 3.225000e+01, 3.224000e+01, 3.223000e+01, 3.222000e+01, 3.221000e+01, 3.220000e+01, 3.219000e+01, 3.218000e+01, 3.217000e+01, 3.216000e+01, 3.215000e+01, 3.214000e+01, 3.213000e+01, 3.212000e+01, 3.211000e+01, 3.210000e+01, 3.209000e+01, 3.208000e+01, 3.207000e+01, 3.206000e+01, 3.205000e+01, 3.204000e+01, 3.203000e+01, 3.202000e+01, 3.201000e+01, 3.200000e+01, 3.199000e+01, 3.198000e+01, 3.197000e+01, 3.196000e+01, 3.195000e+01, 3.194000e+01, 3.193000e+01, 3.192000e+01, 3.191000e+01, 3.190000e+01, 3.189000e+01, 3.188000e+01, 3.187000e+01, 3.186000e+01, 3.185000e+01, 3.184000e+01, 3.183000e+01, 3.182000e+01, 3.181000e+01, 3.180000e+01, 3.179000e+01, 3.178000e+01, 3.177000e+01, 3.176000e+01, 3.175000e+01, 3.174000e+01, 3.173000e+01, 3.172000e+01, 3.171000e+01, 3.170000e+01, 3.169000e+01, 3.168000e+01, 3.167000e+01, 3.166000e+01, 3.165000e+01, 3.164000e+01, 3.163000e+01, 3.162000e+01, 3.161000e+01, 3.160000e+01, 3.159000e+01, 3.158000e+01, 3.157000e+01, 3.156000e+01, 3.155000e+01, 3.154000e+01, 3.153000e+01, 3.152000e+01, 3.151000e+01, 3.150000e+01, 3.149000e+01, 3.148000e+01, 3.147000e+01, 3.146000e+01, 3.145000e+01, 3.144000e+01, 3.143000e+01, 3.142000e+01, 3.141000e+01, 3.140000e+01, 3.139000e+01, 3.138000e+01, 3.137000e+01, 3.136000e+01, 3.135000e+01, 3.134000e+01, 3.133000e+01, 3.132000e+01, 3.131000e+01, 3.130000e+01, 3.129000e+01, 3.128000e+01, 3.127000e+01, 3.126000e+01, 3.125000e+01, 3.124000e+01, 3.123000e+01, 3.122000e+01, 3.121000e+01, 3.120000e+01, 3.119000e+01, 3.118000e+01, 3.117000e+01, 3.116000e+01, 3.115000e+01, 3.114000e+01, 3.113000e+01, 3.112000e+01, 3.111000e+01, 3.110000e+01, 3.109000e+01, 3.108000e+01, 3.107000e+01, 3.106000e+01, 3.105000e+01, 3.104000e+01, 3.103000e+01, 3.102000e+01, 3.101000e+01, 3.100000e+01, 3.099000e+01, 3.098000e+01, 3.097000e+01, 3.096000e+01, 3.095000e+01, 3.094000e+01, 3.093000e+01, 3.092000e+01, 3.091000e+01, 3.090000e+01, 3.089000e+01, 3.088000e+01, 3.087000e+01, 3.086000e+01, 3.085000e+01, 3.084000e+01, 3.083000e+01, 3.082000e+01, 3.081000e+01, 3.080000e+01, 3.079000e+01, 3.078000e+01, 3.077000e+01, 3.076000e+01, 3.075000e+01, 3.074000e+01, 3.073000e+01, 3.072000e+01, 3.071000e+01, 3.070000e+01, 3.069000e+01, 3.068000e+01, 3.067000e+01, 3.066000e+01, 3.065000e+01, 3.064000e+01, 3.063000e+01, 3.062000e+01, 3.061000e+01, 3.060000e+01, 3.059000e+01, 3.058000e+01, 3.057000e+01, 3.056000e+01, 3.055000e+01, 3.054000e+01, 3.053000e+01, 3.052000e+01, 3.051000e+01, 3.050000e+01, 3.049000e+01, 3.048000e+01, 3.047000e+01, 3.046000e+01, 3.045000e+01, 3.044000e+01, 3.043000e+01, 3.042000e+01, 3.041000e+01, 3.040000e+01, 3.039000e+01, 3.038000e+01, 3.037000e+01, 3.036000e+01, 3.035000e+01, 3.034000e+01, 3.033000e+01, 3.032000e+01, 3.031000e+01, 3.030000e+01, 3.029000e+01, 3.028000e+01, 3.027000e+01, 3.026000e+01, 3.025000e+01, 3.024000e+01, 3.023000e+01, 3.022000e+01, 3.021000e+01, 3.020000e+01, 3.019000e+01, 3.018000e+01, 3.017000e+01, 3.016000e+01, 3.015000e+01, 3.014000e+01, 3.013000e+01, 3.012000e+01, 3.011000e+01, 3.010000e+01, 3.009000e+01, 3.008000e+01, 3.007000e+01, 3.006000e+01, 3.005000e+01, 3.004000e+01, 3.003000e+01, 3.002000e+01, 3.001000e+01, 3.000000e+01, 2.999000e+01, 2.998000e+01, 2.997000e+01, 2.996000e+01, 2.995000e+01, 2.994000e+01, 2.993000e+01, 2.992000e+01, 2.991000e+01, 2.990000e+01, 2.989000e+01, 2.988000e+01, 2.987000e+01, 2.986000e+01, 2.985000e+01, 2.984000e+01, 2.983000e+01, 2.982000e+01, 2.981000e+01, 2.980000e+01, 2.979000e+01, 2.978000e+01, 2.977000e+01, 2.976000e+01, 2.975000e+01, 2.974000e+01, 2.973000e+01, 2.972000e+01, 2.971000e+01, 2.970000e+01, 2.969000e+01, 2.968000e+01, 2.967000e+01, 2.966000e+01, 2.965000e+01, 2.964000e+01, 2.963000e+01, 2.962000e+01, 2.961000e+01, 2.960000e+01, 2.959000e+01, 2.958000e+01, 2.957000e+01, 2.956000e+01, 2.955000e+01, 2.954000e+01, 2.953000e+01, 2.952000e+01, 2.951000e+01, 2.950000e+01, 2.949000e+01, 2.948000e+01, 2.947000e+01, 2.946000e+01, 2.945000e+01, 2.944000e+01, 2.943000e+01, 2.942000e+01, 2.941000e+01, 2.940000e+01, 2.939000e+01, 2.938000e+01, 2.937000e+01, 2.936000e+01, 2.935000e+01, 2.934000e+01, 2.933000e+01, 2.932000e+01, 2.931000e+01, 2.930000e+01, 2.929000e+01, 2.928000e+01, 2.927000e+01, 2.926000e+01, 2.925000e+01, 2.924000e+01, 2.923000e+01, 2.922000e+01, 2.921000e+01, 2.920000e+01, 2.919000e+01, 2.918000e+01, 2.917000e+01, 2.916000e+01, 2.915000e+01, 2.914000e+01, 2.913000e+01, 2.912000e+01, 2.911000e+01, 2.910000e+01, 2.909000e+01, 2.908000e+01, 2.907000e+01, 2.906000e+01, 2.905000e+01, 2.904000e+01, 2.903000e+01, 2.902000e+01, 2.901000e+01, 2.900000e+01, 2.899000e+01, 2.898000e+01, 2.897000e+01, 2.896000e+01, 2.895000e+01, 2.894000e+01, 2.893000e+01, 2.892000e+01, 2.891000e+01, 2.890000e+01, 2.889000e+01, 2.888000e+01, 2.887000e+01, 2.886000e+01, 2.885000e+01, 2.884000e+01, 2.883000e+01, 2.882000e+01, 2.881000e+01, 2.880000e+01, 2.879000e+01, 2.878000e+01, 2.877000e+01, 2.876000e+01, 2.875000e+01, 2.874000e+01, 2.873000e+01, 2.872000e+01, 2.871000e+01, 2.870000e+01, 2.869000e+01, 2.868000e+01, 2.867000e+01, 2.866000e+01, 2.865000e+01, 2.864000e+01, 2.863000e+01, 2.862000e+01, 2.861000e+01, 2.860000e+01, 2.859000e+01, 2.858000e+01, 2.857000e+01, 2.856000e+01, 2.855000e+01, 2.854000e+01, 2.853000e+01, 2.852000e+01, 2.851000e+01, 2.850000e+01, 2.849000e+01, 2.848000e+01, 2.847000e+01, 2.846000e+01, 2.845000e+01, 2.844000e+01, 2.843000e+01, 2.842000e+01, 2.841000e+01, 2.840000e+01, 2.839000e+01, 2.838000e+01, 2.837000e+01, 2.836000e+01, 2.835000e+01, 2.834000e+01, 2.833000e+01, 2.832000e+01, 2.831000e+01, 2.830000e+01, 2.829000e+01, 2.828000e+01, 2.827000e+01, 2.826000e+01, 2.825000e+01, 2.824000e+01, 2.823000e+01, 2.822000e+01, 2.821000e+01, 2.820000e+01, 2.819000e+01, 2.818000e+01, 2.817000e+01, 2.816000e+01, 2.815000e+01, 2.814000e+01, 2.813000e+01, 2.812000e+01, 2.811000e+01, 2.810000e+01, 2.809000e+01, 2.808000e+01, 2.807000e+01, 2.806000e+01, 2.805000e+01, 2.804000e+01, 2.803000e+01, 2.802000e+01, 2.801000e+01, 2.800000e+01, 2.799000e+01, 2.798000e+01, 2.797000e+01, 2.796000e+01, 2.795000e+01, 2.794000e+01, 2.793000e+01, 2.792000e+01, 2.791000e+01, 2.790000e+01, 2.789000e+01, 2.788000e+01, 2.787000e+01, 2.786000e+01, 2.785000e+01, 2.784000e+01, 2.783000e+01, 2.782000e+01, 2.781000e+01, 2.780000e+01, 2.779000e+01, 2.778000e+01, 2.777000e+01, 2.776000e+01, 2.775000e+01, 2.774000e+01, 2.773000e+01, 2.772000e+01, 2.771000e+01, 2.770000e+01, 2.769000e+01, 2.768000e+01, 2.767000e+01, 2.766000e+01, 2.765000e+01, 2.764000e+01, 2.763000e+01, 2.762000e+01, 2.761000e+01, 2.760000e+01, 2.759000e+01, 2.758000e+01, 2.757000e+01, 2.756000e+01, 2.755000e+01, 2.754000e+01, 2.753000e+01, 2.752000e+01, 2.751000e+01, 2.750000e+01, 2.749000e+01, 2.748000e+01, 2.747000e+01, 2.746000e+01, 2.745000e+01, 2.744000e+01, 2.743000e+01, 2.742000e+01, 2.741000e+01, 2.740000e+01, 2.739000e+01, 2.738000e+01, 2.737000e+01, 2.736000e+01, 2.735000e+01, 2.734000e+01, 2.733000e+01, 2.732000e+01, 2.731000e+01, 2.730000e+01, 2.729000e+01, 2.728000e+01, 2.727000e+01, 2.726000e+01, 2.725000e+01, 2.724000e+01, 2.723000e+01, 2.722000e+01, 2.721000e+01, 2.720000e+01, 2.719000e+01, 2.718000e+01, 2.717000e+01, 2.716000e+01, 2.715000e+01, 2.714000e+01, 2.713000e+01, 2.712000e+01, 2.711000e+01, 2.710000e+01, 2.709000e+01, 2.708000e+01, 2.707000e+01, 2.706000e+01, 2.705000e+01, 2.704000e+01, 2.703000e+01, 2.702000e+01, 2.701000e+01, 2.700000e+01, 2.699000e+01, 2.698000e+01, 2.697000e+01, 2.696000e+01, 2.695000e+01, 2.694000e+01, 2.693000e+01, 2.692000e+01, 2.691000e+01, 2.690000e+01, 2.689000e+01, 2.688000e+01, 2.687000e+01, 2.686000e+01, 2.685000e+01, 2.684000e+01, 2.683000e+01, 2.682000e+01, 2.681000e+01, 2.680000e+01, 2.679000e+01, 2.678000e+01, 2.677000e+01, 2.676000e+01, 2.675000e+01, 2.674000e+01, 2.673000e+01, 2.672000e+01, 2.671000e+01, 2.670000e+01, 2.669000e+01, 2.668000e+01, 2.667000e+01, 2.666000e+01, 2.665000e+01, 2.664000e+01, 2.663000e+01, 2.662000e+01, 2.661000e+01, 2.660000e+01, 2.659000e+01, 2.658000e+01, 2.657000e+01, 2.656000e+01, 2.655000e+01, 2.654000e+01, 2.653000e+01, 2.652000e+01, 2.651000e+01, 2.650000e+01, 2.649000e+01, 2.648000e+01, 2.647000e+01, 2.646000e+01, 2.645000e+01, 2.644000e+01, 2.643000e+01, 2.642000e+01, 2.641000e+01, 2.640000e+01, 2.639000e+01, 2.638000e+01, 2.637000e+01, 2.636000e+01, 2.635000e+01, 2.634000e+01, 2.633000e+01, 2.632000e+01, 2.631000e+01, 2.630000e+01, 2.629000e+01, 2.628000e+01, 2.627000e+01, 2.626000e+01, 2.625000e+01, 2.624000e+01, 2.623000e+01, 2.622000e+01, 2.621000e+01, 2.620000e+01, 2.619000e+01, 2.618000e+01, 2.617000e+01, 2.616000e+01, 2.615000e+01, 2.614000e+01, 2.613000e+01, 2.612000e+01, 2.611000e+01, 2.610000e+01, 2.609000e+01, 2.608000e+01, 2.607000e+01, 2.606000e+01, 2.605000e+01, 2.604000e+01, 2.603000e+01, 2.602000e+01, 2.601000e+01, 2.600000e+01, 2.599000e+01, 2.598000e+01, 2.597000e+01, 2.596000e+01, 2.595000e+01, 2.594000e+01, 2.593000e+01, 2.592000e+01, 2.591000e+01, 2.590000e+01, 2.589000e+01, 2.588000e+01, 2.587000e+01, 2.586000e+01, 2.585000e+01, 2.584000e+01, 2.583000e+01, 2.582000e+01, 2.581000e+01, 2.580000e+01, 2.579000e+01, 2.578000e+01, 2.577000e+01, 2.576000e+01, 2.575000e+01, 2.574000e+01, 2.573000e+01, 2.572000e+01, 2.571000e+01, 2.570000e+01, 2.569000e+01, 2.568000e+01, 2.567000e+01, 2.566000e+01, 2.565000e+01, 2.564000e+01, 2.563000e+01, 2.562000e+01, 2.561000e+01, 2.560000e+01, 2.559000e+01, 2.558000e+01, 2.557000e+01, 2.556000e+01, 2.555000e+01, 2.554000e+01, 2.553000e+01, 2.552000e+01, 2.551000e+01, 2.550000e+01, 2.549000e+01, 2.548000e+01, 2.547000e+01, 2.546000e+01, 2.545000e+01, 2.544000e+01, 2.543000e+01, 2.542000e+01, 2.541000e+01, 2.540000e+01, 2.539000e+01, 2.538000e+01, 2.537000e+01, 2.536000e+01, 2.535000e+01, 2.534000e+01, 2.533000e+01, 2.532000e+01, 2.531000e+01, 2.530000e+01, 2.529000e+01, 2.528000e+01, 2.527000e+01, 2.526000e+01, 2.525000e+01, 2.524000e+01, 2.523000e+01, 2.522000e+01, 2.521000e+01, 2.520000e+01, 2.519000e+01, 2.518000e+01, 2.517000e+01, 2.516000e+01, 2.515000e+01, 2.514000e+01, 2.513000e+01, 2.512000e+01, 2.511000e+01, 2.510000e+01, 2.509000e+01, 2.508000e+01, 2.507000e+01, 2.506000e+01, 2.505000e+01, 2.504000e+01, 2.503000e+01, 2.502000e+01, 2.501000e+01, 2.500000e+01, 2.499000e+01, 2.498000e+01, 2.497000e+01, 2.496000e+01, 2.495000e+01, 2.494000e+01, 2.493000e+01, 2.492000e+01, 2.491000e+01, 2.490000e+01, 2.489000e+01, 2.488000e+01, 2.487000e+01, 2.486000e+01, 2.485000e+01, 2.484000e+01, 2.483000e+01, 2.482000e+01, 2.481000e+01, 2.480000e+01, 2.479000e+01, 2.478000e+01, 2.477000e+01, 2.476000e+01, 2.475000e+01, 2.474000e+01, 2.473000e+01, 2.472000e+01, 2.471000e+01, 2.470000e+01, 2.469000e+01, 2.468000e+01, 2.467000e+01, 2.466000e+01, 2.465000e+01, 2.464000e+01, 2.463000e+01, 2.462000e+01, 2.461000e+01, 2.460000e+01, 2.459000e+01, 2.458000e+01, 2.457000e+01, 2.456000e+01, 2.455000e+01, 2.454000e+01, 2.453000e+01, 2.452000e+01, 2.451000e+01, 2.450000e+01, 2.449000e+01, 2.448000e+01, 2.447000e+01, 2.446000e+01, 2.445000e+01, 2.444000e+01, 2.443000e+01, 2.442000e+01, 2.441000e+01, 2.440000e+01, 2.439000e+01, 2.438000e+01, 2.437000e+01, 2.436000e+01, 2.435000e+01, 2.434000e+01, 2.433000e+01, 2.432000e+01, 2.431000e+01, 2.430000e+01, 2.429000e+01, 2.428000e+01, 2.427000e+01, 2.426000e+01, 2.425000e+01, 2.424000e+01, 2.423000e+01, 2.422000e+01, 2.421000e+01, 2.420000e+01, 2.419000e+01, 2.418000e+01, 2.417000e+01, 2.416000e+01, 2.415000e+01, 2.414000e+01, 2.413000e+01, 2.412000e+01, 2.411000e+01, 2.410000e+01, 2.409000e+01, 2.408000e+01, 2.407000e+01, 2.406000e+01, 2.405000e+01, 2.404000e+01, 2.403000e+01, 2.402000e+01, 2.401000e+01, 2.400000e+01, 2.399000e+01, 2.398000e+01, 2.397000e+01, 2.396000e+01, 2.395000e+01, 2.394000e+01, 2.393000e+01, 2.392000e+01, 2.391000e+01, 2.390000e+01, 2.389000e+01, 2.388000e+01, 2.387000e+01, 2.386000e+01, 2.385000e+01, 2.384000e+01, 2.383000e+01, 2.382000e+01, 2.381000e+01, 2.380000e+01, 2.379000e+01, 2.378000e+01, 2.377000e+01, 2.376000e+01, 2.375000e+01, 2.374000e+01, 2.373000e+01, 2.372000e+01, 2.371000e+01, 2.370000e+01, 2.369000e+01, 2.368000e+01, 2.367000e+01, 2.366000e+01, 2.365000e+01, 2.364000e+01, 2.363000e+01, 2.362000e+01, 2.361000e+01, 2.360000e+01, 2.359000e+01, 2.358000e+01, 2.357000e+01, 2.356000e+01, 2.355000e+01, 2.354000e+01, 2.353000e+01, 2.352000e+01, 2.351000e+01, 2.350000e+01, 2.349000e+01, 2.348000e+01, 2.347000e+01, 2.346000e+01, 2.345000e+01, 2.344000e+01, 2.343000e+01, 2.342000e+01, 2.341000e+01, 2.340000e+01, 2.339000e+01, 2.338000e+01, 2.337000e+01, 2.336000e+01, 2.335000e+01, 2.334000e+01, 2.333000e+01, 2.332000e+01, 2.331000e+01, 2.330000e+01, 2.329000e+01, 2.328000e+01, 2.327000e+01, 2.326000e+01, 2.325000e+01, 2.324000e+01, 2.323000e+01, 2.322000e+01, 2.321000e+01, 2.320000e+01, 2.319000e+01, 2.318000e+01, 2.317000e+01, 2.316000e+01, 2.315000e+01, 2.314000e+01, 2.313000e+01, 2.312000e+01, 2.311000e+01, 2.310000e+01, 2.309000e+01, 2.308000e+01, 2.307000e+01, 2.306000e+01, 2.305000e+01, 2.304000e+01, 2.303000e+01, 2.302000e+01, 2.301000e+01, 2.300000e+01, 2.299000e+01, 2.298000e+01, 2.297000e+01, 2.296000e+01, 2.295000e+01, 2.294000e+01, 2.293000e+01, 2.292000e+01, 2.291000e+01, 2.290000e+01, 2.289000e+01, 2.288000e+01, 2.287000e+01, 2.286000e+01, 2.285000e+01, 2.284000e+01, 2.283000e+01, 2.282000e+01, 2.281000e+01, 2.280000e+01, 2.279000e+01, 2.278000e+01, 2.277000e+01, 2.276000e+01, 2.275000e+01, 2.274000e+01, 2.273000e+01, 2.272000e+01, 2.271000e+01, 2.270000e+01, 2.269000e+01, 2.268000e+01, 2.267000e+01, 2.266000e+01, 2.265000e+01, 2.264000e+01, 2.263000e+01, 2.262000e+01, 2.261000e+01, 2.260000e+01, 2.259000e+01, 2.258000e+01, 2.257000e+01, 2.256000e+01, 2.255000e+01, 2.254000e+01, 2.253000e+01, 2.252000e+01, 2.251000e+01, 2.250000e+01, 2.249000e+01, 2.248000e+01, 2.247000e+01, 2.246000e+01, 2.245000e+01, 2.244000e+01, 2.243000e+01, 2.242000e+01, 2.241000e+01, 2.240000e+01, 2.239000e+01, 2.238000e+01, 2.237000e+01, 2.236000e+01, 2.235000e+01, 2.234000e+01, 2.233000e+01, 2.232000e+01, 2.231000e+01, 2.230000e+01, 2.229000e+01, 2.228000e+01, 2.227000e+01, 2.226000e+01, 2.225000e+01, 2.224000e+01, 2.223000e+01, 2.222000e+01, 2.221000e+01, 2.220000e+01, 2.219000e+01, 2.218000e+01, 2.217000e+01, 2.216000e+01, 2.215000e+01, 2.214000e+01, 2.213000e+01, 2.212000e+01, 2.211000e+01, 2.210000e+01, 2.209000e+01, 2.208000e+01, 2.207000e+01, 2.206000e+01, 2.205000e+01, 2.204000e+01, 2.203000e+01, 2.202000e+01, 2.201000e+01, 2.200000e+01, 2.199000e+01, 2.198000e+01, 2.197000e+01, 2.196000e+01, 2.195000e+01, 2.194000e+01, 2.193000e+01, 2.192000e+01, 2.191000e+01, 2.190000e+01, 2.189000e+01, 2.188000e+01, 2.187000e+01, 2.186000e+01, 2.185000e+01, 2.184000e+01, 2.183000e+01, 2.182000e+01, 2.181000e+01, 2.180000e+01, 2.179000e+01, 2.178000e+01, 2.177000e+01, 2.176000e+01, 2.175000e+01, 2.174000e+01, 2.173000e+01, 2.172000e+01, 2.171000e+01, 2.170000e+01, 2.169000e+01, 2.168000e+01, 2.167000e+01, 2.166000e+01, 2.165000e+01, 2.164000e+01, 2.163000e+01, 2.162000e+01, 2.161000e+01, 2.160000e+01, 2.159000e+01, 2.158000e+01, 2.157000e+01, 2.156000e+01, 2.155000e+01, 2.154000e+01, 2.153000e+01, 2.152000e+01, 2.151000e+01, 2.150000e+01, 2.149000e+01, 2.148000e+01, 2.147000e+01, 2.146000e+01, 2.145000e+01, 2.144000e+01, 2.143000e+01, 2.142000e+01, 2.141000e+01, 2.140000e+01, 2.139000e+01, 2.138000e+01, 2.137000e+01, 2.136000e+01, 2.135000e+01, 2.134000e+01, 2.133000e+01, 2.132000e+01, 2.131000e+01, 2.130000e+01, 2.129000e+01, 2.128000e+01, 2.127000e+01, 2.126000e+01, 2.125000e+01, 2.124000e+01, 2.123000e+01, 2.122000e+01, 2.121000e+01, 2.120000e+01, 2.119000e+01, 2.118000e+01, 2.117000e+01, 2.116000e+01, 2.115000e+01, 2.114000e+01, 2.113000e+01, 2.112000e+01, 2.111000e+01, 2.110000e+01, 2.109000e+01, 2.108000e+01, 2.107000e+01, 2.106000e+01, 2.105000e+01, 2.104000e+01, 2.103000e+01, 2.102000e+01, 2.101000e+01, 2.100000e+01, 2.099000e+01, 2.098000e+01, 2.097000e+01, 2.096000e+01, 2.095000e+01, 2.094000e+01, 2.093000e+01, 2.092000e+01, 2.091000e+01, 2.090000e+01, 2.089000e+01, 2.088000e+01, 2.087000e+01, 2.086000e+01, 2.085000e+01, 2.084000e+01, 2.083000e+01, 2.082000e+01, 2.081000e+01, 2.080000e+01, 2.079000e+01, 2.078000e+01, 2.077000e+01, 2.076000e+01, 2.075000e+01, 2.074000e+01, 2.073000e+01, 2.072000e+01, 2.071000e+01, 2.070000e+01, 2.069000e+01, 2.068000e+01, 2.067000e+01, 2.066000e+01, 2.065000e+01, 2.064000e+01, 2.063000e+01, 2.062000e+01, 2.061000e+01, 2.060000e+01, 2.059000e+01, 2.058000e+01, 2.057000e+01, 2.056000e+01, 2.055000e+01, 2.054000e+01, 2.053000e+01, 2.052000e+01, 2.051000e+01, 2.050000e+01, 2.049000e+01, 2.048000e+01, 2.047000e+01, 2.046000e+01, 2.045000e+01, 2.044000e+01, 2.043000e+01, 2.042000e+01, 2.041000e+01, 2.040000e+01, 2.039000e+01, 2.038000e+01, 2.037000e+01, 2.036000e+01, 2.035000e+01, 2.034000e+01, 2.033000e+01, 2.032000e+01, 2.031000e+01, 2.030000e+01, 2.029000e+01, 2.028000e+01, 2.027000e+01, 2.026000e+01, 2.025000e+01, 2.024000e+01, 2.023000e+01, 2.022000e+01, 2.021000e+01, 2.020000e+01, 2.019000e+01, 2.018000e+01, 2.017000e+01, 2.016000e+01, 2.015000e+01, 2.014000e+01, 2.013000e+01, 2.012000e+01, 2.011000e+01, 2.010000e+01, 2.009000e+01, 2.008000e+01, 2.007000e+01, 2.006000e+01, 2.005000e+01, 2.004000e+01, 2.003000e+01, 2.002000e+01, 2.001000e+01, 2.000000e+01, 1.999000e+01, 1.998000e+01, 1.997000e+01, 1.996000e+01, 1.995000e+01, 1.994000e+01, 1.993000e+01, 1.992000e+01, 1.991000e+01, 1.990000e+01, 1.989000e+01, 1.988000e+01, 1.987000e+01, 1.986000e+01, 1.985000e+01, 1.984000e+01, 1.983000e+01, 1.982000e+01, 1.981000e+01, 1.980000e+01, 1.979000e+01, 1.978000e+01, 1.977000e+01, 1.976000e+01, 1.975000e+01, 1.974000e+01, 1.973000e+01, 1.972000e+01, 1.971000e+01, 1.970000e+01, 1.969000e+01, 1.968000e+01, 1.967000e+01, 1.966000e+01, 1.965000e+01, 1.964000e+01, 1.963000e+01, 1.962000e+01, 1.961000e+01, 1.960000e+01, 1.959000e+01, 1.958000e+01, 1.957000e+01, 1.956000e+01, 1.955000e+01, 1.954000e+01, 1.953000e+01, 1.952000e+01, 1.951000e+01, 1.950000e+01, 1.949000e+01, 1.948000e+01, 1.947000e+01, 1.946000e+01, 1.945000e+01, 1.944000e+01, 1.943000e+01, 1.942000e+01, 1.941000e+01, 1.940000e+01, 1.939000e+01, 1.938000e+01, 1.937000e+01, 1.936000e+01, 1.935000e+01, 1.934000e+01, 1.933000e+01, 1.932000e+01, 1.931000e+01, 1.930000e+01, 1.929000e+01, 1.928000e+01, 1.927000e+01, 1.926000e+01, 1.925000e+01, 1.924000e+01, 1.923000e+01, 1.922000e+01, 1.921000e+01, 1.920000e+01, 1.919000e+01, 1.918000e+01, 1.917000e+01, 1.916000e+01, 1.915000e+01, 1.914000e+01, 1.913000e+01, 1.912000e+01, 1.911000e+01, 1.910000e+01, 1.909000e+01, 1.908000e+01, 1.907000e+01, 1.906000e+01, 1.905000e+01, 1.904000e+01, 1.903000e+01, 1.902000e+01, 1.901000e+01, 1.900000e+01, 1.899000e+01, 1.898000e+01, 1.897000e+01, 1.896000e+01, 1.895000e+01, 1.894000e+01, 1.893000e+01, 1.892000e+01, 1.891000e+01, 1.890000e+01, 1.889000e+01, 1.888000e+01, 1.887000e+01, 1.886000e+01, 1.885000e+01, 1.884000e+01, 1.883000e+01, 1.882000e+01, 1.881000e+01, 1.880000e+01, 1.879000e+01, 1.878000e+01, 1.877000e+01, 1.876000e+01, 1.875000e+01, 1.874000e+01, 1.873000e+01, 1.872000e+01, 1.871000e+01, 1.870000e+01, 1.869000e+01, 1.868000e+01, 1.867000e+01, 1.866000e+01, 1.865000e+01, 1.864000e+01, 1.863000e+01, 1.862000e+01, 1.861000e+01, 1.860000e+01, 1.859000e+01, 1.858000e+01, 1.857000e+01, 1.856000e+01, 1.855000e+01, 1.854000e+01, 1.853000e+01, 1.852000e+01, 1.851000e+01, 1.850000e+01, 1.849000e+01, 1.848000e+01, 1.847000e+01, 1.846000e+01, 1.845000e+01, 1.844000e+01, 1.843000e+01, 1.842000e+01, 1.841000e+01, 1.840000e+01, 1.839000e+01, 1.838000e+01, 1.837000e+01, 1.836000e+01, 1.835000e+01, 1.834000e+01, 1.833000e+01, 1.832000e+01, 1.831000e+01, 1.830000e+01, 1.829000e+01, 1.828000e+01, 1.827000e+01, 1.826000e+01, 1.825000e+01, 1.824000e+01, 1.823000e+01, 1.822000e+01, 1.821000e+01, 1.820000e+01, 1.819000e+01, 1.818000e+01, 1.817000e+01, 1.816000e+01, 1.815000e+01, 1.814000e+01, 1.813000e+01, 1.812000e+01, 1.811000e+01, 1.810000e+01, 1.809000e+01, 1.808000e+01, 1.807000e+01, 1.806000e+01, 1.805000e+01, 1.804000e+01, 1.803000e+01, 1.802000e+01, 1.801000e+01, 1.800000e+01, 1.799000e+01, 1.798000e+01, 1.797000e+01, 1.796000e+01, 1.795000e+01, 1.794000e+01, 1.793000e+01, 1.792000e+01, 1.791000e+01, 1.790000e+01, 1.789000e+01, 1.788000e+01, 1.787000e+01, 1.786000e+01, 1.785000e+01, 1.784000e+01, 1.783000e+01, 1.782000e+01, 1.781000e+01, 1.780000e+01, 1.779000e+01, 1.778000e+01, 1.777000e+01, 1.776000e+01, 1.775000e+01, 1.774000e+01, 1.773000e+01, 1.772000e+01, 1.771000e+01, 1.770000e+01, 1.769000e+01, 1.768000e+01, 1.767000e+01, 1.766000e+01, 1.765000e+01, 1.764000e+01, 1.763000e+01, 1.762000e+01, 1.761000e+01, 1.760000e+01, 1.759000e+01, 1.758000e+01, 1.757000e+01, 1.756000e+01, 1.755000e+01, 1.754000e+01, 1.753000e+01, 1.752000e+01, 1.751000e+01, 1.750000e+01, 1.749000e+01, 1.748000e+01, 1.747000e+01, 1.746000e+01, 1.745000e+01, 1.744000e+01, 1.743000e+01, 1.742000e+01, 1.741000e+01, 1.740000e+01, 1.739000e+01, 1.738000e+01, 1.737000e+01, 1.736000e+01, 1.735000e+01, 1.734000e+01, 1.733000e+01, 1.732000e+01, 1.731000e+01, 1.730000e+01, 1.729000e+01, 1.728000e+01, 1.727000e+01, 1.726000e+01, 1.725000e+01, 1.724000e+01, 1.723000e+01, 1.722000e+01, 1.721000e+01, 1.720000e+01, 1.719000e+01, 1.718000e+01, 1.717000e+01, 1.716000e+01, 1.715000e+01, 1.714000e+01, 1.713000e+01, 1.712000e+01, 1.711000e+01, 1.710000e+01, 1.709000e+01, 1.708000e+01, 1.707000e+01, 1.706000e+01, 1.705000e+01, 1.704000e+01, 1.703000e+01, 1.702000e+01, 1.701000e+01, 1.700000e+01, 1.699000e+01, 1.698000e+01, 1.697000e+01, 1.696000e+01, 1.695000e+01, 1.694000e+01, 1.693000e+01, 1.692000e+01, 1.691000e+01, 1.690000e+01, 1.689000e+01, 1.688000e+01, 1.687000e+01, 1.686000e+01, 1.685000e+01, 1.684000e+01, 1.683000e+01, 1.682000e+01, 1.681000e+01, 1.680000e+01, 1.679000e+01, 1.678000e+01, 1.677000e+01, 1.676000e+01, 1.675000e+01, 1.674000e+01, 1.673000e+01, 1.672000e+01, 1.671000e+01, 1.670000e+01, 1.669000e+01, 1.668000e+01, 1.667000e+01, 1.666000e+01, 1.665000e+01, 1.664000e+01, 1.663000e+01, 1.662000e+01, 1.661000e+01, 1.660000e+01, 1.659000e+01, 1.658000e+01, 1.657000e+01, 1.656000e+01, 1.655000e+01, 1.654000e+01, 1.653000e+01, 1.652000e+01, 1.651000e+01, 1.650000e+01, 1.649000e+01, 1.648000e+01, 1.647000e+01, 1.646000e+01, 1.645000e+01, 1.644000e+01, 1.643000e+01, 1.642000e+01, 1.641000e+01, 1.640000e+01, 1.639000e+01, 1.638000e+01, 1.637000e+01, 1.636000e+01, 1.635000e+01, 1.634000e+01, 1.633000e+01, 1.632000e+01, 1.631000e+01, 1.630000e+01, 1.629000e+01, 1.628000e+01, 1.627000e+01, 1.626000e+01, 1.625000e+01, 1.624000e+01, 1.623000e+01, 1.622000e+01, 1.621000e+01, 1.620000e+01, 1.619000e+01, 1.618000e+01, 1.617000e+01, 1.616000e+01, 1.615000e+01, 1.614000e+01, 1.613000e+01, 1.612000e+01, 1.611000e+01, 1.610000e+01, 1.609000e+01, 1.608000e+01, 1.607000e+01, 1.606000e+01, 1.605000e+01, 1.604000e+01, 1.603000e+01, 1.602000e+01, 1.601000e+01, 1.600000e+01, 1.599000e+01, 1.598000e+01, 1.597000e+01, 1.596000e+01, 1.595000e+01, 1.594000e+01, 1.593000e+01, 1.592000e+01, 1.591000e+01, 1.590000e+01, 1.589000e+01, 1.588000e+01, 1.587000e+01, 1.586000e+01, 1.585000e+01, 1.584000e+01, 1.583000e+01, 1.582000e+01, 1.581000e+01, 1.580000e+01, 1.579000e+01, 1.578000e+01, 1.577000e+01, 1.576000e+01, 1.575000e+01, 1.574000e+01, 1.573000e+01, 1.572000e+01, 1.571000e+01, 1.570000e+01, 1.569000e+01, 1.568000e+01, 1.567000e+01, 1.566000e+01, 1.565000e+01, 1.564000e+01, 1.563000e+01, 1.562000e+01, 1.561000e+01, 1.560000e+01, 1.559000e+01, 1.558000e+01, 1.557000e+01, 1.556000e+01, 1.555000e+01, 1.554000e+01, 1.553000e+01, 1.552000e+01, 1.551000e+01, 1.550000e+01, 1.549000e+01, 1.548000e+01, 1.547000e+01, 1.546000e+01, 1.545000e+01, 1.544000e+01, 1.543000e+01, 1.542000e+01, 1.541000e+01, 1.540000e+01, 1.539000e+01, 1.538000e+01, 1.537000e+01, 1.536000e+01, 1.535000e+01, 1.534000e+01, 1.533000e+01, 1.532000e+01, 1.531000e+01, 1.530000e+01, 1.529000e+01, 1.528000e+01, 1.527000e+01, 1.526000e+01, 1.525000e+01, 1.524000e+01, 1.523000e+01, 1.522000e+01, 1.521000e+01, 1.520000e+01, 1.519000e+01, 1.518000e+01, 1.517000e+01, 1.516000e+01, 1.515000e+01, 1.514000e+01, 1.513000e+01, 1.512000e+01, 1.511000e+01, 1.510000e+01, 1.509000e+01, 1.508000e+01, 1.507000e+01, 1.506000e+01, 1.505000e+01, 1.504000e+01, 1.503000e+01, 1.502000e+01, 1.501000e+01, 1.500000e+01, 1.499000e+01, 1.498000e+01, 1.497000e+01, 1.496000e+01, 1.495000e+01, 1.494000e+01, 1.493000e+01, 1.492000e+01, 1.491000e+01, 1.490000e+01, 1.489000e+01, 1.488000e+01, 1.487000e+01, 1.486000e+01, 1.485000e+01, 1.484000e+01, 1.483000e+01, 1.482000e+01, 1.481000e+01, 1.480000e+01, 1.479000e+01, 1.478000e+01, 1.477000e+01, 1.476000e+01, 1.475000e+01, 1.474000e+01, 1.473000e+01, 1.472000e+01, 1.471000e+01, 1.470000e+01, 1.469000e+01, 1.468000e+01, 1.467000e+01, 1.466000e+01, 1.465000e+01, 1.464000e+01, 1.463000e+01, 1.462000e+01, 1.461000e+01, 1.460000e+01, 1.459000e+01, 1.458000e+01, 1.457000e+01, 1.456000e+01, 1.455000e+01, 1.454000e+01, 1.453000e+01, 1.452000e+01, 1.451000e+01, 1.450000e+01, 1.449000e+01, 1.448000e+01, 1.447000e+01, 1.446000e+01, 1.445000e+01, 1.444000e+01, 1.443000e+01, 1.442000e+01, 1.441000e+01, 1.440000e+01, 1.439000e+01, 1.438000e+01, 1.437000e+01, 1.436000e+01, 1.435000e+01, 1.434000e+01, 1.433000e+01, 1.432000e+01, 1.431000e+01, 1.430000e+01, 1.429000e+01, 1.428000e+01, 1.427000e+01, 1.426000e+01, 1.425000e+01, 1.424000e+01, 1.423000e+01, 1.422000e+01, 1.421000e+01, 1.420000e+01, 1.419000e+01, 1.418000e+01, 1.417000e+01, 1.416000e+01, 1.415000e+01, 1.414000e+01, 1.413000e+01, 1.412000e+01, 1.411000e+01, 1.410000e+01, 1.409000e+01, 1.408000e+01, 1.407000e+01, 1.406000e+01, 1.405000e+01, 1.404000e+01, 1.403000e+01, 1.402000e+01, 1.401000e+01, 1.400000e+01, 1.399000e+01, 1.398000e+01, 1.397000e+01, 1.396000e+01, 1.395000e+01, 1.394000e+01, 1.393000e+01, 1.392000e+01, 1.391000e+01, 1.390000e+01, 1.389000e+01, 1.388000e+01, 1.387000e+01, 1.386000e+01, 1.385000e+01, 1.384000e+01, 1.383000e+01, 1.382000e+01, 1.381000e+01, 1.380000e+01, 1.379000e+01, 1.378000e+01, 1.377000e+01, 1.376000e+01, 1.375000e+01, 1.374000e+01, 1.373000e+01, 1.372000e+01, 1.371000e+01, 1.370000e+01, 1.369000e+01, 1.368000e+01, 1.367000e+01, 1.366000e+01, 1.365000e+01, 1.364000e+01, 1.363000e+01, 1.362000e+01, 1.361000e+01, 1.360000e+01, 1.359000e+01, 1.358000e+01, 1.357000e+01, 1.356000e+01, 1.355000e+01, 1.354000e+01, 1.353000e+01, 1.352000e+01, 1.351000e+01, 1.350000e+01, 1.349000e+01, 1.348000e+01, 1.347000e+01, 1.346000e+01, 1.345000e+01, 1.344000e+01, 1.343000e+01, 1.342000e+01, 1.341000e+01, 1.340000e+01, 1.339000e+01, 1.338000e+01, 1.337000e+01, 1.336000e+01, 1.335000e+01, 1.334000e+01, 1.333000e+01, 1.332000e+01, 1.331000e+01, 1.330000e+01, 1.329000e+01, 1.328000e+01, 1.327000e+01, 1.326000e+01, 1.325000e+01, 1.324000e+01, 1.323000e+01, 1.322000e+01, 1.321000e+01, 1.320000e+01, 1.319000e+01, 1.318000e+01, 1.317000e+01, 1.316000e+01, 1.315000e+01, 1.314000e+01, 1.313000e+01, 1.312000e+01, 1.311000e+01, 1.310000e+01, 1.309000e+01, 1.308000e+01, 1.307000e+01, 1.306000e+01, 1.305000e+01, 1.304000e+01, 1.303000e+01, 1.302000e+01, 1.301000e+01, 1.300000e+01, 1.299000e+01, 1.298000e+01, 1.297000e+01, 1.296000e+01, 1.295000e+01, 1.294000e+01, 1.293000e+01, 1.292000e+01, 1.291000e+01, 1.290000e+01, 1.289000e+01, 1.288000e+01, 1.287000e+01, 1.286000e+01, 1.285000e+01, 1.284000e+01, 1.283000e+01, 1.282000e+01, 1.281000e+01, 1.280000e+01, 1.279000e+01, 1.278000e+01, 1.277000e+01, 1.276000e+01, 1.275000e+01, 1.274000e+01, 1.273000e+01, 1.272000e+01, 1.271000e+01, 1.270000e+01, 1.269000e+01, 1.268000e+01, 1.267000e+01, 1.266000e+01, 1.265000e+01, 1.264000e+01, 1.263000e+01, 1.262000e+01, 1.261000e+01, 1.260000e+01, 1.259000e+01, 1.258000e+01, 1.257000e+01, 1.256000e+01, 1.255000e+01, 1.254000e+01, 1.253000e+01, 1.252000e+01, 1.251000e+01, 1.250000e+01, 1.249000e+01, 1.248000e+01, 1.247000e+01, 1.246000e+01, 1.245000e+01, 1.244000e+01, 1.243000e+01, 1.242000e+01, 1.241000e+01, 1.240000e+01, 1.239000e+01, 1.238000e+01, 1.237000e+01, 1.236000e+01, 1.235000e+01, 1.234000e+01, 1.233000e+01, 1.232000e+01, 1.231000e+01, 1.230000e+01, 1.229000e+01, 1.228000e+01, 1.227000e+01, 1.226000e+01, 1.225000e+01, 1.224000e+01, 1.223000e+01, 1.222000e+01, 1.221000e+01, 1.220000e+01, 1.219000e+01, 1.218000e+01, 1.217000e+01, 1.216000e+01, 1.215000e+01, 1.214000e+01, 1.213000e+01, 1.212000e+01, 1.211000e+01, 1.210000e+01, 1.209000e+01, 1.208000e+01, 1.207000e+01, 1.206000e+01, 1.205000e+01, 1.204000e+01, 1.203000e+01, 1.202000e+01, 1.201000e+01, 1.200000e+01, 1.199000e+01, 1.198000e+01, 1.197000e+01, 1.196000e+01, 1.195000e+01, 1.194000e+01, 1.193000e+01, 1.192000e+01, 1.191000e+01, 1.190000e+01, 1.189000e+01, 1.188000e+01, 1.187000e+01, 1.186000e+01, 1.185000e+01, 1.184000e+01, 1.183000e+01, 1.182000e+01, 1.181000e+01, 1.180000e+01, 1.179000e+01, 1.178000e+01, 1.177000e+01, 1.176000e+01, 1.175000e+01, 1.174000e+01, 1.173000e+01, 1.172000e+01, 1.171000e+01, 1.170000e+01, 1.169000e+01, 1.168000e+01, 1.167000e+01, 1.166000e+01, 1.165000e+01, 1.164000e+01, 1.163000e+01, 1.162000e+01, 1.161000e+01, 1.160000e+01, 1.159000e+01, 1.158000e+01, 1.157000e+01, 1.156000e+01, 1.155000e+01, 1.154000e+01, 1.153000e+01, 1.152000e+01, 1.151000e+01, 1.150000e+01, 1.149000e+01, 1.148000e+01, 1.147000e+01, 1.146000e+01, 1.145000e+01, 1.144000e+01, 1.143000e+01, 1.142000e+01, 1.141000e+01, 1.140000e+01, 1.139000e+01, 1.138000e+01, 1.137000e+01, 1.136000e+01, 1.135000e+01, 1.134000e+01, 1.133000e+01, 1.132000e+01, 1.131000e+01, 1.130000e+01, 1.129000e+01, 1.128000e+01, 1.127000e+01, 1.126000e+01, 1.125000e+01, 1.124000e+01, 1.123000e+01, 1.122000e+01, 1.121000e+01, 1.120000e+01, 1.119000e+01, 1.118000e+01, 1.117000e+01, 1.116000e+01, 1.115000e+01, 1.114000e+01, 1.113000e+01, 1.112000e+01, 1.111000e+01, 1.110000e+01, 1.109000e+01, 1.108000e+01, 1.107000e+01, 1.106000e+01, 1.105000e+01, 1.104000e+01, 1.103000e+01, 1.102000e+01, 1.101000e+01, 1.100000e+01, 1.099000e+01, 1.098000e+01, 1.097000e+01, 1.096000e+01, 1.095000e+01, 1.094000e+01, 1.093000e+01, 1.092000e+01, 1.091000e+01, 1.090000e+01, 1.089000e+01, 1.088000e+01, 1.087000e+01, 1.086000e+01, 1.085000e+01, 1.084000e+01, 1.083000e+01, 1.082000e+01, 1.081000e+01, 1.080000e+01, 1.079000e+01, 1.078000e+01, 1.077000e+01, 1.076000e+01, 1.075000e+01, 1.074000e+01, 1.073000e+01, 1.072000e+01, 1.071000e+01, 1.070000e+01, 1.069000e+01, 1.068000e+01, 1.067000e+01, 1.066000e+01, 1.065000e+01, 1.064000e+01, 1.063000e+01, 1.062000e+01, 1.061000e+01, 1.060000e+01, 1.059000e+01, 1.058000e+01, 1.057000e+01, 1.056000e+01, 1.055000e+01, 1.054000e+01, 1.053000e+01, 1.052000e+01, 1.051000e+01, 1.050000e+01, 1.049000e+01, 1.048000e+01, 1.047000e+01, 1.046000e+01, 1.045000e+01, 1.044000e+01, 1.043000e+01, 1.042000e+01, 1.041000e+01, 1.040000e+01, 1.039000e+01, 1.038000e+01, 1.037000e+01, 1.036000e+01, 1.035000e+01, 1.034000e+01, 1.033000e+01, 1.032000e+01, 1.031000e+01, 1.030000e+01, 1.029000e+01, 1.028000e+01, 1.027000e+01, 1.026000e+01, 1.025000e+01, 1.024000e+01, 1.023000e+01, 1.022000e+01, 1.021000e+01, 1.020000e+01, 1.019000e+01, 1.018000e+01, 1.017000e+01, 1.016000e+01, 1.015000e+01, 1.014000e+01, 1.013000e+01, 1.012000e+01, 1.011000e+01, 1.010000e+01, 1.009000e+01, 1.008000e+01, 1.007000e+01, 1.006000e+01, 1.005000e+01, 1.004000e+01, 1.003000e+01, 1.002000e+01, 1.001000e+01, 1.000000e+01, 9.990000e+00, 9.980000e+00, 9.970000e+00, 9.960000e+00, 9.950000e+00, 9.940000e+00, 9.930000e+00, 9.920000e+00, 9.910000e+00, 9.900000e+00, 9.890000e+00, 9.880000e+00, 9.870000e+00, 9.860000e+00, 9.850000e+00, 9.840000e+00, 9.830000e+00, 9.820000e+00, 9.810000e+00, 9.800000e+00, 9.790000e+00, 9.780000e+00, 9.770000e+00, 9.760000e+00, 9.750000e+00, 9.740000e+00, 9.730000e+00, 9.720000e+00, 9.710000e+00, 9.700000e+00, 9.690000e+00, 9.680000e+00, 9.670000e+00, 9.660000e+00, 9.650000e+00, 9.640000e+00, 9.630000e+00, 9.620000e+00, 9.610000e+00, 9.600000e+00, 9.590000e+00, 9.580000e+00, 9.570000e+00, 9.560000e+00, 9.550000e+00, 9.540000e+00, 9.530000e+00, 9.520000e+00, 9.510000e+00, 9.500000e+00, 9.490000e+00, 9.480000e+00, 9.470000e+00, 9.460000e+00, 9.450000e+00, 9.440000e+00, 9.430000e+00, 9.420000e+00, 9.410000e+00, 9.400000e+00, 9.390000e+00, 9.380000e+00, 9.370000e+00, 9.360000e+00, 9.350000e+00, 9.340000e+00, 9.330000e+00, 9.320000e+00, 9.310000e+00, 9.300000e+00, 9.290000e+00, 9.280000e+00, 9.270000e+00, 9.260000e+00, 9.250000e+00, 9.240000e+00, 9.230000e+00, 9.220000e+00, 9.210000e+00, 9.200000e+00, 9.190000e+00, 9.180000e+00, 9.170000e+00, 9.160000e+00, 9.150000e+00, 9.140000e+00, 9.130000e+00, 9.120000e+00, 9.110000e+00, 9.100000e+00, 9.090000e+00, 9.080000e+00, 9.070000e+00, 9.060000e+00, 9.050000e+00, 9.040000e+00, 9.030000e+00, 9.020000e+00, 9.010000e+00, 9.000000e+00, 8.990000e+00, 8.980000e+00, 8.970000e+00, 8.960000e+00, 8.950000e+00, 8.940000e+00, 8.930000e+00, 8.920000e+00, 8.910000e+00, 8.900000e+00, 8.890000e+00, 8.880000e+00, 8.870000e+00, 8.860000e+00, 8.850000e+00, 8.840000e+00, 8.830000e+00, 8.820000e+00, 8.810000e+00, 8.800000e+00, 8.790000e+00, 8.780000e+00, 8.770000e+00, 8.760000e+00, 8.750000e+00, 8.740000e+00, 8.730000e+00, 8.720000e+00, 8.710000e+00, 8.700000e+00, 8.690000e+00, 8.680000e+00, 8.670000e+00, 8.660000e+00, 8.650000e+00, 8.640000e+00, 8.630000e+00, 8.620000e+00, 8.610000e+00, 8.600000e+00, 8.590000e+00, 8.580000e+00, 8.570000e+00, 8.560000e+00, 8.550000e+00, 8.540000e+00, 8.530000e+00, 8.520000e+00, 8.510000e+00, 8.500000e+00, 8.490000e+00, 8.480000e+00, 8.470000e+00, 8.460000e+00, 8.450000e+00, 8.440000e+00, 8.430000e+00, 8.420000e+00, 8.410000e+00, 8.400000e+00, 8.390000e+00, 8.380000e+00, 8.370000e+00, 8.360000e+00, 8.350000e+00, 8.340000e+00, 8.330000e+00, 8.320000e+00, 8.310000e+00, 8.300000e+00, 8.290000e+00, 8.280000e+00, 8.270000e+00, 8.260000e+00, 8.250000e+00, 8.240000e+00, 8.230000e+00, 8.220000e+00, 8.210000e+00, 8.200000e+00, 8.190000e+00, 8.180000e+00, 8.170000e+00, 8.160000e+00, 8.150000e+00, 8.140000e+00, 8.130000e+00, 8.120000e+00, 8.110000e+00, 8.100000e+00, 8.090000e+00, 8.080000e+00, 8.070000e+00, 8.060000e+00, 8.050000e+00, 8.040000e+00, 8.030000e+00, 8.020000e+00, 8.010000e+00, 8.000000e+00, 7.990000e+00, 7.980000e+00, 7.970000e+00, 7.960000e+00, 7.950000e+00, 7.940000e+00, 7.930000e+00, 7.920000e+00, 7.910000e+00, 7.900000e+00, 7.890000e+00, 7.880000e+00, 7.870000e+00, 7.860000e+00, 7.850000e+00, 7.840000e+00, 7.830000e+00, 7.820000e+00, 7.810000e+00, 7.800000e+00, 7.790000e+00, 7.780000e+00, 7.770000e+00, 7.760000e+00, 7.750000e+00, 7.740000e+00, 7.730000e+00, 7.720000e+00, 7.710000e+00, 7.700000e+00, 7.690000e+00, 7.680000e+00, 7.670000e+00, 7.660000e+00, 7.650000e+00, 7.640000e+00, 7.630000e+00, 7.620000e+00, 7.610000e+00, 7.600000e+00, 7.590000e+00, 7.580000e+00, 7.570000e+00, 7.560000e+00, 7.550000e+00, 7.540000e+00, 7.530000e+00, 7.520000e+00, 7.510000e+00, 7.500000e+00, 7.490000e+00, 7.480000e+00, 7.470000e+00, 7.460000e+00, 7.450000e+00, 7.440000e+00, 7.430000e+00, 7.420000e+00, 7.410000e+00, 7.400000e+00, 7.390000e+00, 7.380000e+00, 7.370000e+00, 7.360000e+00, 7.350000e+00, 7.340000e+00, 7.330000e+00, 7.320000e+00, 7.310000e+00, 7.300000e+00, 7.290000e+00, 7.280000e+00, 7.270000e+00, 7.260000e+00, 7.250000e+00, 7.240000e+00, 7.230000e+00, 7.220000e+00, 7.210000e+00, 7.200000e+00, 7.190000e+00, 7.180000e+00, 7.170000e+00, 7.160000e+00, 7.150000e+00, 7.140000e+00, 7.130000e+00, 7.120000e+00, 7.110000e+00, 7.100000e+00, 7.090000e+00, 7.080000e+00, 7.070000e+00, 7.060000e+00, 7.050000e+00, 7.040000e+00, 7.030000e+00, 7.020000e+00, 7.010000e+00, 7.000000e+00, 6.990000e+00, 6.980000e+00, 6.970000e+00, 6.960000e+00, 6.950000e+00, 6.940000e+00, 6.930000e+00, 6.920000e+00, 6.910000e+00, 6.900000e+00, 6.890000e+00, 6.880000e+00, 6.870000e+00, 6.860000e+00, 6.850000e+00, 6.840000e+00, 6.830000e+00, 6.820000e+00, 6.810000e+00, 6.800000e+00, 6.790000e+00, 6.780000e+00, 6.770000e+00, 6.760000e+00, 6.750000e+00, 6.740000e+00, 6.730000e+00, 6.720000e+00, 6.710000e+00, 6.700000e+00, 6.690000e+00, 6.680000e+00, 6.670000e+00, 6.660000e+00, 6.650000e+00, 6.640000e+00, 6.630000e+00, 6.620000e+00, 6.610000e+00, 6.600000e+00, 6.590000e+00, 6.580000e+00, 6.570000e+00, 6.560000e+00, 6.550000e+00, 6.540000e+00, 6.530000e+00, 6.520000e+00, 6.510000e+00, 6.500000e+00, 6.490000e+00, 6.480000e+00, 6.470000e+00, 6.460000e+00, 6.450000e+00, 6.440000e+00, 6.430000e+00, 6.420000e+00, 6.410000e+00, 6.400000e+00, 6.390000e+00, 6.380000e+00, 6.370000e+00, 6.360000e+00, 6.350000e+00, 6.340000e+00, 6.330000e+00, 6.320000e+00, 6.310000e+00, 6.300000e+00, 6.290000e+00, 6.280000e+00, 6.270000e+00, 6.260000e+00, 6.250000e+00, 6.240000e+00, 6.230000e+00, 6.220000e+00, 6.210000e+00, 6.200000e+00, 6.190000e+00, 6.180000e+00, 6.170000e+00, 6.160000e+00, 6.150000e+00, 6.140000e+00, 6.130000e+00, 6.120000e+00, 6.110000e+00, 6.100000e+00, 6.090000e+00, 6.080000e+00, 6.070000e+00, 6.060000e+00, 6.050000e+00, 6.040000e+00, 6.030000e+00, 6.020000e+00, 6.010000e+00, 6.000000e+00, 5.990000e+00, 5.980000e+00, 5.970000e+00, 5.960000e+00, 5.950000e+00, 5.940000e+00, 5.930000e+00, 5.920000e+00, 5.910000e+00, 5.900000e+00, 5.890000e+00, 5.880000e+00, 5.870000e+00, 5.860000e+00, 5.850000e+00, 5.840000e+00, 5.830000e+00, 5.820000e+00, 5.810000e+00, 5.800000e+00, 5.790000e+00, 5.780000e+00, 5.770000e+00, 5.760000e+00, 5.750000e+00, 5.740000e+00, 5.730000e+00, 5.720000e+00, 5.710000e+00, 5.700000e+00, 5.690000e+00, 5.680000e+00, 5.670000e+00, 5.660000e+00, 5.650000e+00, 5.640000e+00, 5.630000e+00, 5.620000e+00, 5.610000e+00, 5.600000e+00, 5.590000e+00, 5.580000e+00, 5.570000e+00, 5.560000e+00, 5.550000e+00, 5.540000e+00, 5.530000e+00, 5.520000e+00, 5.510000e+00, 5.500000e+00, 5.490000e+00, 5.480000e+00, 5.470000e+00, 5.460000e+00, 5.450000e+00, 5.440000e+00, 5.430000e+00, 5.420000e+00, 5.410000e+00, 5.400000e+00, 5.390000e+00, 5.380000e+00, 5.370000e+00, 5.360000e+00, 5.350000e+00, 5.340000e+00, 5.330000e+00, 5.320000e+00, 5.310000e+00, 5.300000e+00, 5.290000e+00, 5.280000e+00, 5.270000e+00, 5.260000e+00, 5.250000e+00, 5.240000e+00, 5.230000e+00, 5.220000e+00, 5.210000e+00, 5.200000e+00, 5.190000e+00, 5.180000e+00, 5.170000e+00, 5.160000e+00, 5.150000e+00, 5.140000e+00, 5.130000e+00, 5.120000e+00, 5.110000e+00, 5.100000e+00, 5.090000e+00, 5.080000e+00, 5.070000e+00, 5.060000e+00, 5.050000e+00, 5.040000e+00, 5.030000e+00, 5.020000e+00, 5.010000e+00, 5.000000e+00, 4.990000e+00, 4.980000e+00, 4.970000e+00, 4.960000e+00, 4.950000e+00, 4.940000e+00, 4.930000e+00, 4.920000e+00, 4.910000e+00, 4.900000e+00, 4.890000e+00, 4.880000e+00, 4.870000e+00, 4.860000e+00, 4.850000e+00, 4.840000e+00, 4.830000e+00, 4.820000e+00, 4.810000e+00, 4.800000e+00, 4.790000e+00, 4.780000e+00, 4.770000e+00, 4.760000e+00, 4.750000e+00, 4.740000e+00, 4.730000e+00, 4.720000e+00, 4.710000e+00, 4.700000e+00, 4.690000e+00, 4.680000e+00, 4.670000e+00, 4.660000e+00, 4.650000e+00, 4.640000e+00, 4.630000e+00, 4.620000e+00, 4.610000e+00, 4.600000e+00, 4.590000e+00, 4.580000e+00, 4.570000e+00, 4.560000e+00, 4.550000e+00, 4.540000e+00, 4.530000e+00, 4.520000e+00, 4.510000e+00, 4.500000e+00, 4.490000e+00, 4.480000e+00, 4.470000e+00, 4.460000e+00, 4.450000e+00, 4.440000e+00, 4.430000e+00, 4.420000e+00, 4.410000e+00, 4.400000e+00, 4.390000e+00, 4.380000e+00, 4.370000e+00, 4.360000e+00, 4.350000e+00, 4.340000e+00, 4.330000e+00, 4.320000e+00, 4.310000e+00, 4.300000e+00, 4.290000e+00, 4.280000e+00, 4.270000e+00, 4.260000e+00, 4.250000e+00, 4.240000e+00, 4.230000e+00, 4.220000e+00, 4.210000e+00, 4.200000e+00, 4.190000e+00, 4.180000e+00, 4.170000e+00, 4.160000e+00, 4.150000e+00, 4.140000e+00, 4.130000e+00, 4.120000e+00, 4.110000e+00, 4.100000e+00, 4.090000e+00, 4.080000e+00, 4.070000e+00, 4.060000e+00, 4.050000e+00, 4.040000e+00, 4.030000e+00, 4.020000e+00, 4.010000e+00, 4.000000e+00, 3.990000e+00, 3.980000e+00, 3.970000e+00, 3.960000e+00, 3.950000e+00, 3.940000e+00, 3.930000e+00, 3.920000e+00, 3.910000e+00, 3.900000e+00, 3.890000e+00, 3.880000e+00, 3.870000e+00, 3.860000e+00, 3.850000e+00, 3.840000e+00, 3.830000e+00, 3.820000e+00, 3.810000e+00, 3.800000e+00, 3.790000e+00, 3.780000e+00, 3.770000e+00, 3.760000e+00, 3.750000e+00, 3.740000e+00, 3.730000e+00, 3.720000e+00, 3.710000e+00, 3.700000e+00, 3.690000e+00, 3.680000e+00, 3.670000e+00, 3.660000e+00, 3.650000e+00, 3.640000e+00, 3.630000e+00, 3.620000e+00, 3.610000e+00, 3.600000e+00, 3.590000e+00, 3.580000e+00, 3.570000e+00, 3.560000e+00, 3.550000e+00, 3.540000e+00, 3.530000e+00, 3.520000e+00, 3.510000e+00, 3.500000e+00, 3.490000e+00, 3.480000e+00, 3.470000e+00, 3.460000e+00, 3.450000e+00, 3.440000e+00, 3.430000e+00, 3.420000e+00, 3.410000e+00, 3.400000e+00, 3.390000e+00, 3.380000e+00, 3.370000e+00, 3.360000e+00, 3.350000e+00, 3.340000e+00, 3.330000e+00, 3.320000e+00, 3.310000e+00, 3.300000e+00, 3.290000e+00, 3.280000e+00, 3.270000e+00, 3.260000e+00, 3.250000e+00, 3.240000e+00, 3.230000e+00, 3.220000e+00, 3.210000e+00, 3.200000e+00, 3.190000e+00, 3.180000e+00, 3.170000e+00, 3.160000e+00, 3.150000e+00, 3.140000e+00, 3.130000e+00, 3.120000e+00, 3.110000e+00, 3.100000e+00, 3.090000e+00, 3.080000e+00, 3.070000e+00, 3.060000e+00, 3.050000e+00, 3.040000e+00, 3.030000e+00, 3.020000e+00, 3.010000e+00, 3.000000e+00, 2.990000e+00, 2.980000e+00, 2.970000e+00, 2.960000e+00, 2.950000e+00, 2.940000e+00, 2.930000e+00, 2.920000e+00, 2.910000e+00, 2.900000e+00, 2.890000e+00, 2.880000e+00, 2.870000e+00, 2.860000e+00, 2.850000e+00, 2.840000e+00, 2.830000e+00, 2.820000e+00, 2.810000e+00, 2.800000e+00, 2.790000e+00, 2.780000e+00, 2.770000e+00, 2.760000e+00, 2.750000e+00, 2.740000e+00, 2.730000e+00, 2.720000e+00, 2.710000e+00, 2.700000e+00, 2.690000e+00, 2.680000e+00, 2.670000e+00, 2.660000e+00, 2.650000e+00, 2.640000e+00, 2.630000e+00, 2.620000e+00, 2.610000e+00, 2.600000e+00, 2.590000e+00, 2.580000e+00, 2.570000e+00, 2.560000e+00, 2.550000e+00, 2.540000e+00, 2.530000e+00, 2.520000e+00, 2.510000e+00, 2.500000e+00, 2.490000e+00, 2.480000e+00, 2.470000e+00, 2.460000e+00, 2.450000e+00, 2.440000e+00, 2.430000e+00, 2.420000e+00, 2.410000e+00, 2.400000e+00, 2.390000e+00, 2.380000e+00, 2.370000e+00, 2.360000e+00, 2.350000e+00, 2.340000e+00, 2.330000e+00, 2.320000e+00, 2.310000e+00, 2.300000e+00, 2.290000e+00, 2.280000e+00, 2.270000e+00, 2.260000e+00, 2.250000e+00, 2.240000e+00, 2.230000e+00, 2.220000e+00, 2.210000e+00, 2.200000e+00, 2.190000e+00, 2.180000e+00, 2.170000e+00, 2.160000e+00, 2.150000e+00, 2.140000e+00, 2.130000e+00, 2.120000e+00, 2.110000e+00, 2.100000e+00, 2.090000e+00, 2.080000e+00, 2.070000e+00, 2.060000e+00, 2.050000e+00, 2.040000e+00, 2.030000e+00, 2.020000e+00, 2.010000e+00, 2.000000e+00, 1.990000e+00, 1.980000e+00, 1.970000e+00, 1.960000e+00, 1.950000e+00, 1.940000e+00, 1.930000e+00, 1.920000e+00, 1.910000e+00, 1.900000e+00, 1.890000e+00, 1.880000e+00, 1.870000e+00, 1.860000e+00, 1.850000e+00, 1.840000e+00, 1.830000e+00, 1.820000e+00, 1.810000e+00, 1.800000e+00, 1.790000e+00, 1.780000e+00, 1.770000e+00, 1.760000e+00, 1.750000e+00, 1.740000e+00, 1.730000e+00, 1.720000e+00, 1.710000e+00, 1.700000e+00, 1.690000e+00, 1.680000e+00, 1.670000e+00, 1.660000e+00, 1.650000e+00, 1.640000e+00, 1.630000e+00, 1.620000e+00, 1.610000e+00, 1.600000e+00, 1.590000e+00, 1.580000e+00, 1.570000e+00, 1.560000e+00, 1.550000e+00, 1.540000e+00, 1.530000e+00, 1.520000e+00, 1.510000e+00, 1.500000e+00, 1.490000e+00, 1.480000e+00, 1.470000e+00, 1.460000e+00, 1.450000e+00, 1.440000e+00, 1.430000e+00, 1.420000e+00, 1.410000e+00, 1.400000e+00, 1.390000e+00, 1.380000e+00, 1.370000e+00, 1.360000e+00, 1.350000e+00, 1.340000e+00, 1.330000e+00, 1.320000e+00, 1.310000e+00, 1.300000e+00, 1.290000e+00, 1.280000e+00, 1.270000e+00, 1.260000e+00, 1.250000e+00, 1.240000e+00, 1.230000e+00, 1.220000e+00, 1.210000e+00, 1.200000e+00, 1.190000e+00, 1.180000e+00, 1.170000e+00, 1.160000e+00, 1.150000e+00, 1.140000e+00, 1.130000e+00, 1.120000e+00, 1.110000e+00, 1.100000e+00, 1.090000e+00, 1.080000e+00, 1.070000e+00, 1.060000e+00, 1.050000e+00, 1.040000e+00, 1.030000e+00, 1.020000e+00, 1.010000e+00, 1.000000e+00, 9.900000e-01, 9.800000e-01, 9.700000e-01, 9.600000e-01, 9.500000e-01, 9.400000e-01, 9.300000e-01, 9.200000e-01, 9.100000e-01, 9.000000e-01, 8.900000e-01, 8.800000e-01, 8.700000e-01, 8.600000e-01, 8.500000e-01, 8.400000e-01, 8.300000e-01, 8.200000e-01, 8.100000e-01, 8.000000e-01, 7.900000e-01, 7.800000e-01, 7.700000e-01, 7.600000e-01, 7.500000e-01, 7.400000e-01, 7.300000e-01, 7.200000e-01, 7.100000e-01, 7.000000e-01, 6.900000e-01, 6.800000e-01, 6.700000e-01, 6.600000e-01, 6.500000e-01, 6.400000e-01, 6.300000e-01, 6.200000e-01, 6.100000e-01, 6.000000e-01, 5.900000e-01, 5.800000e-01, 5.700000e-01, 5.600000e-01, 5.500000e-01, 5.400000e-01, 5.300000e-01, 5.200000e-01, 5.100000e-01, 5.000000e-01, 4.900000e-01, 4.800000e-01, 4.700000e-01, 4.600000e-01, 4.500000e-01, 4.400000e-01, 4.300000e-01, 4.200000e-01, 4.100000e-01, 4.000000e-01, 3.900000e-01, 3.800000e-01, 3.700000e-01, 3.600000e-01, 3.500000e-01, 3.400000e-01, 3.300000e-01, 3.200000e-01, 3.100000e-01, 3.000000e-01, 2.900000e-01, 2.800000e-01, 2.700000e-01, 2.600000e-01, 2.500000e-01, 2.400000e-01, 2.300000e-01, 2.200000e-01, 2.100000e-01, 2.000000e-01, 1.900000e-01, 1.800000e-01, 1.700000e-01, 1.600000e-01, 1.500000e-01, 1.400000e-01, 1.300000e-01, 1.200000e-01, 1.100000e-01, 1.000000e-01, 9.000000e-02, 8.000000e-02, 7.000000e-02, 6.000000e-02, 5.000000e-02, 4.000000e-02, 3.000000e-02, 2.000000e-02, 1.000000e-02, 0.000000e+00};

      fCenMetric = new TGraph(nbins,xd,yd);
      fTempList->Add(fCenMetric);

    }

    fPolSlope[0] = new TF1("slopeper","expo",1.E2,1.E4);
    Double_t parper[] = {9.45678e-01,-5.66325e-04};
    fPolSlope[0]->SetParameters(parper);
    fTempList->Add(fPolSlope[0]);

    fPolSlope[1] = new TF1("slopecen","landaun",1.E4,2.5E4);
    Double_t parcen[] = {1.35360e+02,8.35468e+03,2.39137e+03};
    fPolSlope[1]->SetParameters(parcen);
    fTempList->Add(fPolSlope[1]);

  }
  
  if(fZDCESEList) {
	// functions
    for (Int_t k=0; k<2; k++) {
      fPolMin[k] = new TF1(Form("fPolMin[%d]",k),"pol9",0.,100.);
      fPolMax[k] = new TF1(Form("fPolMax[%d]",k),"pol9",0.,100.);
      fTempList->Add(fPolMin[k]);
      fTempList->Add(fPolMax[k]);
    }
  }

  if(fUseZDC && fRecenterZDC && fStoreZDCQVecVtxPos) {
     for(Int_t i=0;i<2;i++) {
        for(Int_t c=0;c<fCRCnCen;c++) {
          for (Int_t j=0;j<2;j++) {
            fCRCZDCQVecDis[i][c][j] = new TH2D(Form("fCRCZDCQVecDis[%d][%d][%d]",i,c,j),Form("fCRCZDCQVecDis[%d][%d][%d]",i,c,j),200,-2.,2.,200,-2.,2.);
            fCRCZDCQVecDis[i][c][j]->Sumw2();
            fCRCQVecList->Add(fCRCZDCQVecDis[i][c][j]);
          }
        }
     }
      Double_t cenbins[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
      Double_t phibinsforZDCEP[101] = {0.};
      for (Int_t phib=0; phib<101; phib++) {
        phibinsforZDCEP[phib] = phib*TMath::TwoPi()/100.;
      }
      for(Int_t k=0; k<2; k++) {
        fZDCEPHist[k] = new TH2D(Form("fZDCEPHist[%d]",k),Form("fZDCEPHist[%d]",k),10,cenbins,100,phibinsforZDCEP);
        fCRCQVecList->Add(fZDCEPHist[k]);
      }
  }
  if (fUseZDC && fRecenterZDC) {
    for(Int_t r=0;r<fCRCnRun;r++) {

      for(Int_t i=0;i<2;i++) {
        fCRCZDCQVecA[r][i] = new TProfile(Form("fCRCZDCQVecA[%d][%d]",fRunList[r],i),
                                          Form("fCRCZDCQVecA[%d][%d]",fRunList[r],i),100,0.,100.,"s");
        fCRCZDCQVecA[r][i]->Sumw2();
        fCRCQVecListRun[r]->Add(fCRCZDCQVecA[r][i]);
        fCRCZDCQVecC[r][i] = new TProfile(Form("fCRCZDCQVecC[%d][%d]",fRunList[r],i),
                                          Form("fCRCZDCQVecC[%d][%d]",fRunList[r],i),100,0.,100.,"s");
        fCRCZDCQVecC[r][i]->Sumw2();
        fCRCQVecListRun[r]->Add(fCRCZDCQVecC[r][i]);

        fCRCZDCQVecACorr[r][i] = new TProfile(Form("fCRCZDCQVecACorr[%d][%d]",fRunList[r],i),
                                              Form("fCRCZDCQVecACorr[%d][%d]",fRunList[r],i),100,0.,100.,"s");
        fCRCZDCQVecACorr[r][i]->Sumw2();
        fCRCQVecListRun[r]->Add(fCRCZDCQVecACorr[r][i]);
        fCRCZDCQVecCCorr[r][i] = new TProfile(Form("fCRCZDCQVecCCorr[%d][%d]",fRunList[r],i),
                                              Form("fCRCZDCQVecCCorr[%d][%d]",fRunList[r],i),100,0.,100.,"s");
        fCRCZDCQVecCCorr[r][i]->Sumw2();
        fCRCQVecListRun[r]->Add(fCRCZDCQVecCCorr[r][i]);
        
	  }
	 
      for(Int_t i=0;i<4;i++) {
        fCRCZDCQVecRes[r][i] = new TProfile(Form("fCRCZDCQVecRes[%d][%d]",fRunList[r],i),
                                            Form("fCRCZDCQVecRes[%d][%d]",fRunList[r],i),100,0.,100.,"s");
        fCRCZDCQVecRes[r][i]->Sumw2();
        fCRCQVecListRun[r]->Add(fCRCZDCQVecRes[r][i]);
      }
    }
    
    for(Int_t c=0;c<4;c++) {
      fCRCZDCQVecCorSteps[c] = new TProfile2D(Form("fCRCZDCQVecCorSteps[%d]",c),Form("fCRCZDCQVecCorSteps[%d]",c),fkNsteps,0.,1.*fkNsteps,20,0.,100.);
      fCRCZDCQVecCorSteps[c]->Sumw2();
      fCRCQVecList->Add(fCRCZDCQVecCorSteps[c]);
    }
    
  }
  
  fCRCQVecListTPC = new TList();
  fCRCQVecListTPC->SetName("various");
  fCRCQVecListTPC->SetOwner(kTRUE);
  fCRCQVecList->Add(fCRCQVecListTPC);
  
  fhAvRefMulRbR = new TProfile2D("fhAvRefMulRbR","AvRefMulRbR;run;centrality(%)",fCRCnRun,0.,1.*fCRCnRun,90,0.,90.);
  for (Int_t r=0; r<fCRCnRun; r++) fhAvRefMulRbR->GetXaxis()->SetBinLabel(r+1,Form("%d",fRunList[r]));
  fCRCQVecListTPC->Add(fhAvRefMulRbR);
  fhAvQMCRbR = new TProfile2D("fhAvQMCRbR","fhAvQMCRbR;run;centrality(%)",fCRCnRun,0.,1.*fCRCnRun,90,0.,90.);
  for (Int_t r=0; r<fCRCnRun; r++) fhAvQMCRbR->GetXaxis()->SetBinLabel(r+1,Form("%d",fRunList[r]));
  fCRCQVecListTPC->Add(fhAvQMCRbR);
  fhAvQMARbR = new TProfile2D("fhAvQMARbR","fhAvQMARbR;run;centrality(%)",fCRCnRun,0.,1.*fCRCnRun,90,0.,90.);
  for (Int_t r=0; r<fCRCnRun; r++) fhAvQMARbR->GetXaxis()->SetBinLabel(r+1,Form("%d",fRunList[r]));
  fCRCQVecListTPC->Add(fhAvQMARbR);

  fhAvAbsOrbit = new TProfile("fhAvAbsOrbit","fhAvAbsOrbit;run",fCRCnRun,0.,1.*fCRCnRun,"s");
  for (Int_t r=0; r<fCRCnRun; r++) fhAvAbsOrbit->GetXaxis()->SetBinLabel(r+1,Form("%d",fRunList[r]));
  fCRCQVecListTPC->Add(fhAvAbsOrbit);
  
  if(fPOIExtraWeights==AliFlowAnalysisQvec::kEtaPhi) {
    if(fWeightsList->FindObject("fCRCQVecPhiHist")) {
      fPhiEtaWeights = (TH3D*)(fWeightsList->FindObject("fCRCQVecPhiHist"));
    } else {
      AliWarning("WARNING: POIExtraWeights (kEtaPhi) not found ! \n");
    }
  }
  if(fPOIExtraWeights==AliFlowAnalysisQvec::kEtaPhiCh) {
    for (Int_t i=0; i<2; i++) {
      if(fWeightsList->FindObject(Form("fCRCQVecPhiHistCh[%d]",i))) {
        fPhiEtaWeightsCh[i] = (TH3D*)(fWeightsList->FindObject(Form("fCRCQVecPhiHistCh[%d]",i)));
      } else {
        AliWarning("WARNING: POIExtraWeights (kEtaPhiCh) not found ! \n");
      }
    }
  }
  if(fPOIExtraWeights==AliFlowAnalysisQvec::kEtaPhiVtx) {
    for (Int_t cb=0; cb<fCRCnCen; cb++) {
      if(fWeightsList->FindObject(Form("fCRCQVecPhiHistVtx[%d]",cb))) {
        fPhiEtaWeightsVtx[cb] = (TH3D*)(fWeightsList->FindObject(Form("fCRCQVecPhiHistVtx[%d]",cb)));
      } else {
        AliWarning("WARNING: POIExtraWeights (kEtaPhiVtx) not found ! \n");
      }
    }
  }
  if(fPOIExtraWeights==AliFlowAnalysisQvec::kEtaPhiChPt) {
    for (Int_t i=0; i<2; i++) {
      for (Int_t j=0; j<3; j++) {
        if(fWeightsList->FindObject(Form("fCRCQVecPhiHistChPt[%d][%d]",i,j))) {
          fPhiEtaWeightsChPt[i][j] = (TH3D*)(fWeightsList->FindObject(Form("fCRCQVecPhiHistChPt[%d][%d]",i,j)));
        } else {
          AliWarning("WARNING: POIExtraWeights (kEtaPhiChPt) not found ! \n");
        }
      }
    }
  }
}

//=======================================================================================================================

void AliFlowAnalysisQvec::BookEverythingForCME()
{
  // EbE quantities
  for(Int_t c=0;c<4;c++) {
    for (Int_t h=0;h<fCRCnHar;h++) {
      fCMEQRe[c][h] = new TH1D(Form("fCMEQRe[%d][%d]",c,h),Form("fCMEQRe[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaMin,fCRCEtaMax);
      fTempList->Add(fCMEQRe[c][h]);
      fCMEQIm[c][h] = new TH1D(Form("fCMEQIm[%d][%d]",c,h),Form("fCMEQIm[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaMin,fCRCEtaMax);
      fTempList->Add(fCMEQIm[c][h]);
      fCMEMult[c][h] = new TH1D(Form("fCMEMult[%d][%d]",c,h),Form("fCMEMult[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaMin,fCRCEtaMax);
      fTempList->Add(fCMEMult[c][h]);
    }
  }
  

  //@shi book CME Qvector for spectator plane participant plane method
  for(Int_t c=0;c<2;c++) { // c is the index for the power of weight: weight^c*cos((h+1)*phi)
    for (Int_t h=0;h<fCRCnHar;h++) {
      fCMEQReBothCharge[c][h] = new TH1D(Form("fCMEQReBothCharge[%d][%d]",c,h),Form("fCMEQReBothCharge[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaMin,fCRCEtaMax);
      fTempList->Add(fCMEQReBothCharge[c][h]);
      fCMEQImBothCharge[c][h] = new TH1D(Form("fCMEQImBothCharge[%d][%d]",c,h),Form("fCMEQImBothCharge[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaMin,fCRCEtaMax);
      fTempList->Add(fCMEQImBothCharge[c][h]);
      fCMEMultBothCharge[c][h] = new TH1D(Form("fCMEMultBothCharge[%d][%d]",c,h),Form("fCMEMultBothCharge[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaMin,fCRCEtaMax);
      fTempList->Add(fCMEMultBothCharge[c][h]);
    }
  }
  
  if(!fCalculateCME){return;}
  if(!fUseZDC){return;}

}

//=======================================================================================================================

void AliFlowAnalysisQvec::GetPointersForParticleWeightsHistograms()
{
  // Get pointers for histograms with particle weights.
//  TList *weightsList = dynamic_cast<TList*>(fHistList->FindObject("Weights"));
//  if(!weightsList){printf("\n WARNING (QC): weightsList is NULL in AFAWQC::GPFPWH() !!!!\n");exit(0);}
//  this->SetWeightsList(weightsList);
//  TString fUseParticleWeightsName = "fUseParticleWeightsQC"; // to be improved (hirdwired label QC)
//  fUseParticleWeightsName += fAnalysisLabel->Data();
//  TProfile *useParticleWeights = dynamic_cast<TProfile*>(weightsList->FindObject(fUseParticleWeightsName.Data()));
//  if(useParticleWeights)
//  {
//    this->SetUseParticleWeights(useParticleWeights);
//    fUsePhiWeights = (Int_t)fUseParticleWeights->GetBinContent(1);
//    fUsePtWeights = (Int_t)fUseParticleWeights->GetBinContent(2);
//    fUseEtaWeights = (Int_t)fUseParticleWeights->GetBinContent(3);
//    fUseTrackWeights = (Int_t)fUseParticleWeights->GetBinContent(4);
//    fUsePhiEtaWeights = (Int_t)fUseParticleWeights->GetBinContent(5);
//    fUsePhiEtaWeightsChDep = (Int_t)fUseParticleWeights->GetBinContent(6);
//    fUsePhiEtaWeightsVtxDep = (Int_t)fUseParticleWeights->GetBinContent(7);
//    fUsePhiEtaWeightsChPtDep = (Int_t)fUseParticleWeights->GetBinContent(8);
//  }
} // end of void AliFlowAnalysisQvec::GetPointersForParticleWeightsHistograms();

//=======================================================================================================================

void AliFlowAnalysisQvec::InitializeArraysForParticleWeights()
{
  fWeightsList = new TList();
  fPtWeightsCent = NULL;
  fPhiEtaWeights = NULL;
  fPhiEtaRbRWeights = NULL;
  for (Int_t i=0; i<2; i++) {
	fPhiEtaWeightsCh[i] = NULL;
    fPhiEtaRbRWeightsCh[i] = NULL;
  }
  for (Int_t i=0; i<2; i++) {
    for (Int_t j=0; j<3; j++) {
      fPhiEtaWeightsChPt[i][j] = NULL;
    }
  }
  for (Int_t cb=0; cb<fCRCnCen; cb++) {
    fPhiEtaWeightsVtx[cb] = NULL;
  }
}

//=======================================================================================================================

void AliFlowAnalysisQvec::RecenterCRCQVecVZERO()
{
  if(!fCRCVZEROCalibList) {
    cout << " WARNING: no weights provided for VZERO recentering !!! " << endl;
    return;
  }
  //cout<<"========> RecenterCRCQVecVZERO() ====== "<<endl;
  if(fRunNum!=fCachedRunNum) {
    for(Int_t k=0; k<3; k++) {
      fVZEROCenHist[k] = (TProfile2D*)(fCRCVZEROCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCVZEROQVec[%d][%d]",fRunNum,k)));
    }
  }

  for(Int_t k=0; k<3; k++) {
	  //cout<<"=======> k ===== "<<k<<endl;
    // ZDCN-C
    Double_t QCRe = fVZFlowVect[0][k].X();
    Double_t QCIm = fVZFlowVect[0][k].Y();
    Double_t QMC  = fVZFlowVect[0][k].GetMult();
    // ZDCN-A
    Double_t QARe = fVZFlowVect[1][k].X();
    Double_t QAIm = fVZFlowVect[1][k].Y();
    Double_t QMA  = fVZFlowVect[1][k].GetMult();
	
	//cout<<"======> QMC ======== "<<QMC<<endl;
	//cout<<"======> QMA ======== "<<QMC<<endl;
    Double_t QCReR=QCRe, QCImR=QCIm, QAReR=QARe, QAImR=QAIm;
    if(QMC<1. || QMA < 1.) return;

    // recenter vs centrality
    if (fVZEROCenHist[k]) {
		//cout<<"======> Okay fVZEROCenHist["<<k<<"] exists! "<<endl;
      Double_t AvQCRe = fVZEROCenHist[k]->GetBinContent(fVZEROCenHist[k]->FindBin((Double_t)k+0.5,fCentralityEBE));
      // Double_t SDQCRe = fVZEROCenHist[k]->GetBinError(fVZEROCenHist[k]->FindBin((Double_t)k+0.5,fCentralityEBE));
      Double_t AvQCIm = fVZEROCenHist[k]->GetBinContent(fVZEROCenHist[k]->FindBin((Double_t)k+0.5,fCentralityEBE));
      // Double_t SDQCIm = fVZEROCenHist[k]->GetBinError(fVZEROCenHist[k]->FindBin((Double_t)k+0.5,fCentralityEBE));

      Double_t AvQARe = fVZEROCenHist[k]->GetBinContent(fVZEROCenHist[k]->FindBin((Double_t)k+0.5,fCentralityEBE));
      // Double_t SDQARe = fVZEROCenHist[k]->GetBinError(fVZEROCenHist[k]->FindBin((Double_t)k+0.5,fCentralityEBE));
      Double_t AvQAIm = fVZEROCenHist[k]->GetBinContent(fVZEROCenHist[k]->FindBin((Double_t)k+0.5,fCentralityEBE));
      // Double_t SDQAIm = fVZEROCenHist[k]->GetBinError(fVZEROCenHist[k]->FindBin((Double_t)k+0.5,fCentralityEBE));

      if(std::isfinite(AvQCRe) && std::isfinite(AvQCIm)) {
		  //cout<<"======> Okay std::isfinite(AvQCRe) && std::isfinite(AvQCIm)! "<<endl;
		  //cout<<"======> QCReR = QCRe-AvQCRe ===== "<<QCRe<<"-"<<AvQCRe<<" == "<<QCReR<<endl;
          //cout<<"======> QCImR = QCIm-AvQCIm ===== "<<QCIm<<"-"<<AvQCIm<<" == "<<QCImR<<endl;
        QCReR = QCRe-AvQCRe;
        QCImR = QCIm-AvQCIm;
//        if(SDQCRe>0. && SDQCIm>0.) {
//          QCReR /= SDQCRe;
//          QCImR /= SDQCIm;
//        }
        fVZFlowVect[0][k].Set(QCReR,QCImR);
        //@Shi VZ eta < 0
        /*Double_t VZCRe = fVZFlowVect[0][k].X();
        Double_t VZCIm = fVZFlowVect[0][k].Y();
        Double_t EvPlVZC = TMath::ATan2(VZCIm,VZCRe)/2.;
        cout<<"=====> VZCRe [0]["<<k<<"] === "<<VZCRe<<endl;
        cout<<"=====> VZCIm [0]["<<k<<"] === "<<VZCIm<<endl;
        cout<<"=====> EvPlVZC ["<<k<<"] === "<<EvPlVZC<<endl;*/
      }

      if(std::isfinite(AvQCRe) && std::isfinite(AvQCIm)) {
        QAReR = QARe-AvQARe;
        QAImR = QAIm-AvQAIm;
//        if(SDQARe>0. && SDQAIm>0.) {
//          QAReR /= SDQARe;
//          QAImR /= SDQAIm;
//        }
        fVZFlowVect[1][k].Set(QAReR,QAImR);
        
        //@shi VZ eta > 0
        /*Double_t VZARe = fVZFlowVect[1][k].X();
        Double_t VZAIm = fVZFlowVect[1][k].Y();
        Double_t EvPlVZA = TMath::ATan2(VZAIm,VZARe)/2.;
        cout<<"=====> VZARe [1]["<<k<<"] === "<<VZARe<<endl;
        cout<<"=====> VZAIm [1]["<<k<<"] === "<<VZAIm<<endl;
        cout<<"=====> EvPlVZA ["<<k<<"] === "<<EvPlVZA<<endl;*/
      }
    }
  }

}

//=======================================================================================================================

void AliFlowAnalysisQvec::CalculateCRCQVec()
{
  // ZDC
  if(fUseZDC && fRecenterZDC) {
    // ZDC-C (eta < -8.8)
    Double_t VCRe = fZDCFlowVect[0].X();
    Double_t VCIm = fZDCFlowVect[0].Y();
    Double_t VCM  = fZDCFlowVect[0].GetMult();
    if(VCM>0. && sqrt(VCRe*VCRe+VCIm*VCIm)>1.E-6) {
      fCRCZDCQVecC[fRunBin][0]->Fill(fCentralityEBE,VCRe);
      fCRCZDCQVecC[fRunBin][1]->Fill(fCentralityEBE,VCIm);
    }
    // ZDC-A (eta > 8.8)
    Double_t VARe = fZDCFlowVect[1].X();
    Double_t VAIm = fZDCFlowVect[1].Y();
    Double_t VAM  = fZDCFlowVect[1].GetMult();
    if(VAM>0. && sqrt(VARe*VARe+VAIm*VAIm)>1.E-6) {
      fCRCZDCQVecA[fRunBin][0]->Fill(fCentralityEBE,VARe);
      fCRCZDCQVecA[fRunBin][1]->Fill(fCentralityEBE,VAIm);
    }
  } // end of if(fUseVZERO)
	
}

//=======================================================================================================================

void AliFlowAnalysisQvec::SetRunList()
{
  Int_t dAny[1] = {1};

  Int_t dRun10h[92] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161};

  Int_t dRun10hPos[92] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364};

  Int_t dRun10hNeg[92] = {138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161};

  Int_t dRun11h[119] = {167902, 167903, 167915, 167920, 167985, 167987, 167988, 168066, 168068, 168069, 168076, 168104, 168105, 168107, 168108, 168115, 168212, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168461, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168984, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169143, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169923, 169956, 169965, 170027, 170036,170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593};

  // 12 low IR: 244917, 244918, 244975, 244980, 244982, 244983, 245064, 245066, 245068, 246390, 246391, 246392
  // 78 high IR ("CentralBarrelTracking" good runs): 246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246871, 246870, 246867, 246865, 246864, 246859, 246858, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246540, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683

  Int_t dRun15h[] = {244917, 244918, 244975, 244980, 244982, 244983, 245064, 245066, 245068, 246390, 246391, 246392, 246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683, 246148}; // @Shi add 246148

  Int_t dRun15ov6[] = {244918, 244975, 244980, 244982, 244983, 245064, 245066, 245068, 246390, 246391, 246392, 246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683};

  Int_t dRun15opidfix[] = {245145, 245146, 245151, 245152, 245231, 245232, 245259, 245343, 245345, 245346, 245347, 245349, 245353, 245396, 245397, 245401, 245407, 245409, 245441, 245446, 245450, 245454, 245496, 245497, 245501, 245504, 245505, 245507, 245535, 245540, 245542, 245543, 245544, 245545, 245554};

  Int_t dRun15hHIR[] = {246982, 246980, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246805, 246804, 246751, 246750, 246488, 246487, 246431, 246428, 246424, 246271, 246217, 246180, 246178, 246115, 246113, 246042, 246037, 246036, 245923, 245683};

  Int_t dRun15hLIR[] = {246994, 246991, 246989, 246810, 246809, 246766, 246765, 246763, 246760, 246495, 246493, 246276, 246275, 246225, 246185, 246153, 246089, 246053, 246052, 246012, 246003, 245954, 245952, 245949, 245833, 245831, 245705, 245702, 245700};

  Int_t dRun15hPos[] = {246390, 246391, 246392, 246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424};

  Int_t dRun15hNeg[] = {244917, 244918, 244975, 244980, 244982, 244983, 245064, 245066, 245068, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683};

  Double_t dVtxPosX15o[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,7.619407e-02, 7.612905e-02, 7.609009e-02, 7.610981e-02, 7.608885e-02, 7.609981e-02, 7.559263e-02, 7.563009e-02, 7.551201e-02, 7.570994e-02, 7.571927e-02, 7.575639e-02, 7.571133e-02, 7.570653e-02, 7.528412e-02, 7.535235e-02, 7.539954e-02, 7.535435e-02, 7.541641e-02, 7.543658e-02, 7.527343e-02, 7.526024e-02, 7.528295e-02, 7.533821e-02, 7.540461e-02, 7.538317e-02, 7.531677e-02, 7.539861e-02, 7.537667e-02, 7.659318e-02, 7.656796e-02, 7.662898e-02, 7.664257e-02, 7.597872e-02, 7.597437e-02, 7.599091e-02, 7.601310e-02, 7.000359e-02, 6.999659e-02, 6.992559e-02, 6.996793e-02, 7.028519e-02, 7.032696e-02, 7.033503e-02, 6.952509e-02, 6.956378e-02, 6.952446e-02, 6.959759e-02, 6.956048e-02, 6.933134e-02, 6.932882e-02, 6.939338e-02, 6.950613e-02, 6.943631e-02, 6.946196e-02, 6.950454e-02, 7.030973e-02, 7.030203e-02, 7.032272e-02, 7.030936e-02, 7.038967e-02, 7.035136e-02, 7.024752e-02, 6.942316e-02, 6.940115e-02, 6.936367e-02, 6.860689e-02, 6.881501e-02, 6.886743e-02, 6.932714e-02, 6.970325e-02, 6.966504e-02, 6.957355e-02, 6.932303e-02, 6.938184e-02, 6.944933e-02, 6.952461e-02, 6.964167e-02, 0.}; //@Shi add 0 for 246148
  Double_t dVtxPosY15o[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,3.361709e-01, 3.361818e-01, 3.362205e-01, 3.363199e-01, 3.363092e-01, 3.362369e-01, 3.374328e-01, 3.374148e-01, 3.375140e-01, 3.361514e-01, 3.361743e-01, 3.362329e-01, 3.361395e-01, 3.361633e-01, 3.367675e-01, 3.366963e-01, 3.366845e-01, 3.366490e-01, 3.366937e-01, 3.366825e-01, 3.373764e-01, 3.373762e-01, 3.373721e-01, 3.373705e-01, 3.373943e-01, 3.373675e-01, 3.374071e-01, 3.373368e-01, 3.373442e-01, 3.375773e-01, 3.375333e-01, 3.377335e-01, 3.378285e-01, 3.362674e-01, 3.362492e-01, 3.362604e-01, 3.363473e-01, 3.295003e-01, 3.295046e-01, 3.295761e-01, 3.296100e-01, 3.291527e-01, 3.292071e-01, 3.290824e-01, 3.299371e-01, 3.300008e-01, 3.300078e-01, 3.300391e-01, 3.300740e-01, 3.300345e-01, 3.300776e-01, 3.301195e-01, 3.289427e-01, 3.289736e-01, 3.296084e-01, 3.297025e-01, 3.297724e-01, 3.298166e-01, 3.298278e-01, 3.298682e-01, 3.297381e-01, 3.296875e-01, 3.297720e-01, 3.298361e-01, 3.298561e-01, 3.299325e-01, 3.300111e-01, 3.301161e-01, 3.302630e-01, 3.289954e-01, 3.292915e-01, 3.293319e-01, 3.294174e-01, 3.314355e-01, 3.314431e-01, 3.316189e-01, 3.318682e-01, 3.323906e-01, 0.}; //@Shi add 0 for 246148
  Double_t dVtxPosZ15o[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,5.559279e-01, 3.535446e-01, 4.846955e-01, 4.525585e-01, 3.684501e-01, 2.485494e-01, 2.372653e-01, 1.707859e-01, 3.314213e-01, 1.709195e-01, 2.209753e-01, 3.125757e-01, 3.422085e-01, 3.868156e-01, 4.859695e-01, 4.780697e-01, 4.400149e-01, 4.014992e-01, 3.049883e-01, 3.708501e-01, 3.883566e-01, 3.940632e-01, 4.197670e-01, 3.938399e-01, 3.814413e-01, 3.335539e-01, 3.181929e-01, 2.300734e-01, 2.722395e-01, 5.241033e-01, 3.225908e-01, 1.925791e-01, 1.892765e-01, 3.384066e-01, 2.026459e-01, 2.495699e-01, 3.569992e-01, 3.891381e-01, 4.603724e-01, 3.696685e-01, 3.002207e-01, 2.929533e-01, 3.095468e-01, 3.517200e-01, 2.784445e-01, 3.866626e-01, 3.058719e-01, 3.336752e-01, 3.226473e-01, 3.222815e-01, 3.428469e-01, 3.728514e-01, 2.858642e-01, 2.832485e-01, 3.378933e-01, 3.547548e-01, 3.799414e-01, 4.043543e-01, 4.314049e-01, 4.141138e-01, 3.888746e-01, 4.103586e-01, 3.871045e-01, 4.614473e-01, 4.023404e-01, 4.203531e-01, 4.401272e-01, 6.450558e-01, 6.819582e-01, 2.588529e-01, 3.693471e-01, 3.990708e-01, 3.813842e-01, 3.471682e-01, 3.356156e-01, 2.550150e-01, 3.830723e-01, 4.293259e-01, 0.}; //@Shi add 0 for 246148
  Double_t dVtxPosX15oPos[] = {0.000000e+00, 0.000000e+00, 0.000000e+00, 7.619407e-02, 7.612905e-02, 7.609009e-02, 7.610981e-02, 7.608885e-02, 7.609981e-02, 7.559263e-02, 7.563009e-02, 7.551201e-02, 7.570994e-02, 7.571927e-02, 7.575639e-02, 7.571133e-02, 7.570653e-02, 7.528412e-02, 7.535235e-02, 7.539954e-02, 7.535435e-02, 7.541641e-02, 7.543658e-02, 7.527343e-02, 7.526024e-02, 7.528295e-02, 7.533821e-02, 7.540461e-02, 7.538317e-02, 7.531677e-02, 7.539861e-02, 7.537667e-02, 7.659318e-02, 7.656796e-02, 7.662898e-02, 7.664257e-02, 7.597872e-02, 7.597437e-02, 7.599091e-02, 7.601310e-02};
  Double_t dVtxPosY15oPos[] = {0.000000e+00, 0.000000e+00, 0.000000e+00, 3.361709e-01, 3.361818e-01, 3.362205e-01, 3.363199e-01, 3.363092e-01, 3.362369e-01, 3.374328e-01, 3.374148e-01, 3.375140e-01, 3.361514e-01, 3.361743e-01, 3.362329e-01, 3.361395e-01, 3.361633e-01, 3.367675e-01, 3.366963e-01, 3.366845e-01, 3.366490e-01, 3.366937e-01, 3.366825e-01, 3.373764e-01, 3.373762e-01, 3.373721e-01, 3.373705e-01, 3.373943e-01, 3.373675e-01, 3.374071e-01, 3.373368e-01, 3.373442e-01, 3.375773e-01, 3.375333e-01, 3.377335e-01, 3.378285e-01, 3.362674e-01, 3.362492e-01, 3.362604e-01, 3.363473e-01};
  Double_t dVtxPosZ15oPos[] = {0.000000e+00, 0.000000e+00, 0.000000e+00, 5.559279e-01, 3.535446e-01, 4.846955e-01, 4.525585e-01, 3.684501e-01, 2.485494e-01, 2.372653e-01, 1.707859e-01, 3.314213e-01, 1.709195e-01, 2.209753e-01, 3.125757e-01, 3.422085e-01, 3.868156e-01, 4.859695e-01, 4.780697e-01, 4.400149e-01, 4.014992e-01, 3.049883e-01, 3.708501e-01, 3.883566e-01, 3.940632e-01, 4.197670e-01, 3.938399e-01, 3.814413e-01, 3.335539e-01, 3.181929e-01, 2.300734e-01, 2.722395e-01, 5.241033e-01, 3.225908e-01, 1.925791e-01, 1.892765e-01, 3.384066e-01, 2.026459e-01, 2.495699e-01, 3.569992e-01};
  Double_t dVtxPosX15oNeg[] = {0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 7.000359e-02, 6.999659e-02, 6.992559e-02, 6.996793e-02, 7.028519e-02, 7.032696e-02, 7.033503e-02, 6.952509e-02, 6.956378e-02, 6.952446e-02, 6.959759e-02, 6.956048e-02, 6.933134e-02, 6.932882e-02, 6.939338e-02, 6.950613e-02, 6.943631e-02, 6.946196e-02, 6.950454e-02, 7.030973e-02, 7.030203e-02, 7.032272e-02, 7.030936e-02, 7.038967e-02, 7.035136e-02, 7.024752e-02, 6.942316e-02, 6.940115e-02, 6.936367e-02, 6.860689e-02, 6.881501e-02, 6.886743e-02, 6.932714e-02, 6.970325e-02, 6.966504e-02, 6.957355e-02, 6.932303e-02, 6.938184e-02, 6.944933e-02, 6.952461e-02, 6.964167e-02};
  Double_t dVtxPosY15oNeg[] = {0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 3.295003e-01, 3.295046e-01, 3.295761e-01, 3.296100e-01, 3.291527e-01, 3.292071e-01, 3.290824e-01, 3.299371e-01, 3.300008e-01, 3.300078e-01, 3.300391e-01, 3.300740e-01, 3.300345e-01, 3.300776e-01, 3.301195e-01, 3.289427e-01, 3.289736e-01, 3.296084e-01, 3.297025e-01, 3.297724e-01, 3.298166e-01, 3.298278e-01, 3.298682e-01, 3.297381e-01, 3.296875e-01, 3.297720e-01, 3.298361e-01, 3.298561e-01, 3.299325e-01, 3.300111e-01, 3.301161e-01, 3.302630e-01, 3.289954e-01, 3.292915e-01, 3.293319e-01, 3.294174e-01, 3.314355e-01, 3.314431e-01, 3.316189e-01, 3.318682e-01, 3.323906e-01};
  Double_t dVtxPosZ15oNeg[] = {0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 3.891381e-01, 4.603724e-01, 3.696685e-01, 3.002207e-01, 2.929533e-01, 3.095468e-01, 3.517200e-01, 2.784445e-01, 3.866626e-01, 3.058719e-01, 3.336752e-01, 3.226473e-01, 3.222815e-01, 3.428469e-01, 3.728514e-01, 2.858642e-01, 2.832485e-01, 3.378933e-01, 3.547548e-01, 3.799414e-01, 4.043543e-01, 4.314049e-01, 4.141138e-01, 3.888746e-01, 4.103586e-01, 3.871045e-01, 4.614473e-01, 4.023404e-01, 4.203531e-01, 4.401272e-01, 6.450558e-01, 6.819582e-01, 2.588529e-01, 3.693471e-01, 3.990708e-01, 3.813842e-01, 3.471682e-01, 3.356156e-01, 2.550150e-01, 3.830723e-01, 4.293259e-01};
  Double_t dVtxPosX15oHI[] = {7.608885e-02, 7.609981e-02, 7.559263e-02, 7.563009e-02, 7.551201e-02, 7.570994e-02, 7.571927e-02, 7.575639e-02, 7.571133e-02, 7.570653e-02, 7.541641e-02, 7.543658e-02, 7.539861e-02, 7.537667e-02, 7.662898e-02, 7.664257e-02, 7.597437e-02, 7.599091e-02, 7.601310e-02, 6.996793e-02, 7.033503e-02, 6.959759e-02, 6.956048e-02, 6.950613e-02, 6.943631e-02, 7.038967e-02, 7.035136e-02, 7.024752e-02, 6.932714e-02, 6.964167e-02};
  Double_t dVtxPosY15oHI[] = {3.363092e-01, 3.362369e-01, 3.374328e-01, 3.374148e-01, 3.375140e-01, 3.361514e-01, 3.361743e-01, 3.362329e-01, 3.361395e-01, 3.361633e-01, 3.366937e-01, 3.366825e-01, 3.373368e-01, 3.373442e-01, 3.377335e-01, 3.378285e-01, 3.362492e-01, 3.362604e-01, 3.363473e-01, 3.296100e-01, 3.290824e-01, 3.300391e-01, 3.300740e-01, 3.289427e-01, 3.289736e-01, 3.297381e-01, 3.296875e-01, 3.297720e-01, 3.289954e-01, 3.323906e-01};
  Double_t dVtxPosZ15oHI[] = {3.684501e-01, 2.485494e-01, 2.372653e-01, 1.707859e-01, 3.314213e-01, 1.709195e-01, 2.209753e-01, 3.125757e-01, 3.422085e-01, 3.868156e-01, 3.049883e-01, 3.708501e-01, 2.300734e-01, 2.722395e-01, 1.925791e-01, 1.892765e-01, 2.026459e-01, 2.495699e-01, 3.569992e-01, 3.002207e-01, 3.517200e-01, 3.336752e-01, 3.226473e-01, 2.858642e-01, 2.832485e-01, 3.888746e-01, 4.103586e-01, 3.871045e-01, 2.588529e-01, 4.293259e-01};
  Double_t dVtxPosX15oLI[] = {7.619407e-02, 7.612905e-02, 7.609009e-02, 7.528412e-02, 7.535235e-02, 7.527343e-02, 7.526024e-02, 7.528295e-02, 7.533821e-02, 7.659318e-02, 7.656796e-02, 7.000359e-02, 6.999659e-02, 7.028519e-02, 6.952509e-02, 6.933134e-02, 6.946196e-02, 7.030973e-02, 7.030203e-02, 6.942316e-02, 6.940115e-02, 6.860689e-02, 6.881501e-02, 6.886743e-02, 6.970325e-02, 6.966504e-02, 6.932303e-02, 6.938184e-02, 6.944933e-02};
  Double_t dVtxPosY15oLI[] = {3.361709e-01, 3.361818e-01, 3.362205e-01, 3.367675e-01, 3.366963e-01, 3.373764e-01, 3.373762e-01, 3.373721e-01, 3.373705e-01, 3.375773e-01, 3.375333e-01, 3.295003e-01, 3.295046e-01, 3.291527e-01, 3.299371e-01, 3.300345e-01, 3.296084e-01, 3.297724e-01, 3.298166e-01, 3.298361e-01, 3.298561e-01, 3.300111e-01, 3.301161e-01, 3.302630e-01, 3.292915e-01, 3.293319e-01, 3.314355e-01, 3.314431e-01, 3.316189e-01};
  Double_t dVtxPosZ15oLI[] = {5.559279e-01, 3.535446e-01, 4.846955e-01, 4.859695e-01, 4.780697e-01, 3.883566e-01, 3.940632e-01, 4.197670e-01, 3.938399e-01, 5.241033e-01, 3.225908e-01, 3.891381e-01, 4.603724e-01, 2.929533e-01, 2.784445e-01, 3.222815e-01, 3.378933e-01, 3.799414e-01, 4.043543e-01, 4.614473e-01, 4.023404e-01, 4.401272e-01, 6.450558e-01, 6.819582e-01, 3.693471e-01, 3.990708e-01, 3.471682e-01, 3.356156e-01, 2.550150e-01};
  Int_t InEvRbR[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1041837, 284865, 2046572, 1371621, 120035, 2122479, 268047, 595924, 50724, 387311, 539797, 229419, 375195, 249722, 122848, 1446667, 55789, 952960, 982841, 406604, 93217, 237924, 398054, 707103, 63508, 462518, 222559, 456809, 720679, 261426, 838890, 2640961, 296971, 1302825, 638981, 81286, 957516, 414835, 745119, 2065488, 795171, 788238, 1497695, 1408541, 135788, 1116152, 1089236, 746357, 349322, 2419892, 241610, 1119258, 144318, 504444, 207053, 5421086, 1133858, 110452, 1655491, 166314, 1439016, 494122, 90017, 716385, 469564, 1625482, 1948001, 2530161, 432644, 100636, 305525, 1332600, 879753, 404590, 123586, 672729, 890654, 1182816, 0}; //@Shi add 0 for 246148
  Double_t dVtxPosX15opidfix[] = {6.793435e-02, 6.802185e-02, 6.801235e-02, 6.804823e-02, 6.842972e-02, 6.839652e-02, 6.851932e-02, 6.976507e-02, 6.989692e-02, 6.994544e-02, 6.994261e-02, 6.997887e-02, 7.001687e-02, 6.934462e-02, 6.958349e-02, 6.907266e-02, 6.905944e-02, 6.895395e-02, 7.006562e-02, 7.008493e-02, 7.012736e-02, 6.964645e-02, 6.960466e-02, 6.962255e-02, 6.979086e-02, 6.985343e-02, 6.983755e-02, 6.957177e-02, 6.875991e-02, 6.871756e-02, 6.871021e-02, 6.871769e-02, 6.869493e-02, 6.874049e-02, 6.860300e-02};
  Double_t dVtxPosY15opidfix[] = {3.315020e-01, 3.312268e-01, 3.310778e-01, 3.310524e-01, 3.314478e-01, 3.312986e-01, 3.311297e-01, 3.324064e-01, 3.322524e-01, 3.322019e-01, 3.321221e-01, 3.321050e-01, 3.319118e-01, 3.317922e-01, 3.314658e-01, 3.315735e-01, 3.316331e-01, 3.316525e-01, 3.308030e-01, 3.308038e-01, 3.306947e-01, 3.305741e-01, 3.316492e-01, 3.316117e-01, 3.314973e-01, 3.314110e-01, 3.313450e-01, 3.313649e-01, 3.325841e-01, 3.324226e-01, 3.323649e-01, 3.323381e-01, 3.322566e-01, 3.322077e-01, 3.320860e-01};
  Double_t dVtxPosZ15opidfix[] = {4.723797e-01, 4.684324e-01, 4.609304e-01, 4.554974e-01, 4.523016e-01, 3.769890e-01, 4.485548e-01, 5.024484e-01, 5.200088e-01, 5.261731e-01, 5.392851e-01, 5.399264e-01, 5.155504e-01, 4.267668e-01, 5.348764e-01, 4.526746e-01, 4.045626e-01, 4.261759e-01, 5.889205e-01, 6.364843e-01, 5.896163e-01, 3.768637e-01, 4.440771e-01, 4.687029e-01, 4.794467e-01, 4.313422e-01, 3.954777e-01, 3.983129e-01, 3.608064e-01, 2.627038e-01, 3.665826e-01, 4.275667e-01, 3.335445e-01, 3.250815e-01, 3.022907e-01};

// @Shi ave vtx 15o IR splitting
  Double_t dVtxPosX15oIRSplit[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0760396, 0.0761476, 0.0754612, 0.0760119, 0.0771416, 0.0758959, 0.0758307, 0.0764312, 0.0712992, 0.0740006, 0.076035, 0.0795941, 0.0758193, 0.0753836, 0.0759469, 0.0753271, 0.0748559, 0.0755779, 0.0747179, 0.0743789, 0.074226, 0.0738555, 0.0741127, 0.0755049, 0.079727, 0.0754529, 0.0747599, 0.0744282, 0.0742795, 0.0750923, 0.0765961, 0.0762358, 0.0765928, 0.0752035, 0.0767834, 0.0759724, 0.0758235, 0.0690952, 0.0693622, 0.0695388, 0.0704506, 0.070026, 0.0703322, 0.0702859, 0.0695319, 0.0684041, 0.0683909, 0.0696078, 0.0699702, 0.0689661, 0.0677066, 0.0689856, 0.0714685, 0.0690362, 0.0703379, 0.0692874, 0.0702451, 0.0693919, 0.0693631, 0.0702106, 0.0703336, 0.0696804, 0.0668393, 0.0696303, 0.0684486, 0.0693902, 0.0682269, 0.0686902, 0.0688619, 0.069442, 0.0705462, 0.0695982, 0.069336, 0.0685833, 0.0677059, 0.0690834, 0.0691257, 0.0690399, 0.0695431};
  
  Double_t dVtxPosY15oIRSplit[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.336358, 0.336082, 0.33601, 0.336299, 0.338571, 0.336182, 0.336955, 0.33598, 0.337921, 0.334241, 0.335462, 0.337443, 0.334584, 0.336416, 0.335418, 0.336588, 0.338577, 0.335504, 0.336177, 0.336241, 0.338079, 0.338119, 0.337332, 0.336716, 0.340298, 0.337025, 0.337512, 0.337696, 0.336138, 0.338704, 0.336543, 0.337053, 0.335586, 0.335519, 0.335771, 0.334203, 0.335871, 0.32961, 0.329341, 0.328825, 0.330096, 0.328709, 0.329233, 0.329063, 0.329943, 0.330227, 0.329343, 0.330058, 0.32979, 0.330226, 0.330673, 0.330379, 0.325801, 0.329745, 0.327493, 0.329334, 0.329097, 0.331733, 0.330179, 0.329786, 0.330113, 0.327863, 0.331576, 0.329589, 0.329758, 0.32966, 0.329914, 0.329771, 0.330217, 0.327307, 0.32939, 0.329085, 0.329112, 0.331714, 0.327878, 0.331697, 0.330765, 0.331914, 0.33046};
  
  Double_t dVtxPosZ15oIRSplit[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.575355, 0.371561, 0.513012, 0.471856, 0.406659, 0.284534, 0.257454, 0.199246, 0.325792, 0.212815, 0.26548, 0.348875, 0.368373, 0.431977, 0.528148, 0.508048, 0.478282, 0.436153, 0.330369, 0.406381, 0.40795, 0.411249, 0.445349, 0.412348, 0.391552, 0.353029, 0.338251, 0.251904, 0.293615, 0.544099, 0.352431, 0.221797, 0.232368, 0.35809, 0.234556, 0.300599, 0.375358, 0.418464, 0.476625, 0.385246, 0.333402, 0.314478, 0.326505, 0.375008, 0.289914, 0.410377, 0.33794, 0.331634, 0.347134, 0.343325, 0.367387, 0.400036, 0.307101, 0.300977, 0.357842, 0.377861, 0.401782, 0.432738, 0.446801, 0.43286, 0.416691, 0.423076, 0.398294, 0.479613, 0.422342, 0.443408, 0.455862, 0.656827, 0.704932, 0.289011, 0.392294, 0.419466, 0.396562, 0.377537, 0.347602, 0.296413, 0.405798, 0.462996, 0.440022};
///////////////////////////////

  switch(fDataSet) {
    case kAny:
      fCRCnRun=1;
      fRunList=TArrayI(1,dAny);
      break;
    case k2010:
      if(fInteractionRate==kAll) {
        fCRCnRun=92;
        fRunList=TArrayI(fCRCnRun,dRun10h);
      } else if (fInteractionRate==kPos) {
        fCRCnRun=45;
        fRunList=TArrayI(fCRCnRun,dRun10hPos);
      } else if (fInteractionRate==kNeg) {
        fCRCnRun=47;
        fRunList=TArrayI(fCRCnRun,dRun10hNeg);
      }
      fEnNucl=1380.;
      break;
    case k2011:
      fCRCnRun=119;
      fRunList=TArrayI(fCRCnRun,dRun11h);
      fEnNucl=1380.;
      break;
    case k2015:
      if(fInteractionRate==kAll) {
        fCRCnRun=91;
        fRunList=TArrayI(fCRCnRun,dRun15h);
        fAvVtxPosX=TArrayD(fCRCnRun,dVtxPosX15o);
        fAvVtxPosY=TArrayD(fCRCnRun,dVtxPosY15o);
        fAvVtxPosZ=TArrayD(fCRCnRun,dVtxPosZ15o);

        // @Shi diff ave vtx for IR split
        fAvVtxPosX15oIRSplit=TArrayD(fCRCnRun,dVtxPosX15oIRSplit);
        fAvVtxPosY15oIRSplit=TArrayD(fCRCnRun,dVtxPosY15oIRSplit);
        fAvVtxPosZ15oIRSplit=TArrayD(fCRCnRun,dVtxPosZ15oIRSplit);
      } else if (fInteractionRate==kHigh) {
        fCRCnRun=30;
        fRunList=TArrayI(fCRCnRun,dRun15hHIR);
        fAvVtxPosX=TArrayD(fCRCnRun,dVtxPosX15oHI);
        fAvVtxPosY=TArrayD(fCRCnRun,dVtxPosY15oHI);
        fAvVtxPosZ=TArrayD(fCRCnRun,dVtxPosZ15oHI);
      } else if (fInteractionRate==kLow) {
        fCRCnRun=29;
        fRunList=TArrayI(fCRCnRun,dRun15hLIR);
        fAvVtxPosX=TArrayD(fCRCnRun,dVtxPosX15oLI);
        fAvVtxPosY=TArrayD(fCRCnRun,dVtxPosY15oLI);
        fAvVtxPosZ=TArrayD(fCRCnRun,dVtxPosZ15oLI);
      } else if (fInteractionRate==kPos) {
        fCRCnRun=40;
        fRunList=TArrayI(fCRCnRun,dRun15hPos);
        fAvVtxPosX=TArrayD(fCRCnRun,dVtxPosX15oPos);
        fAvVtxPosY=TArrayD(fCRCnRun,dVtxPosY15oPos);
        fAvVtxPosZ=TArrayD(fCRCnRun,dVtxPosZ15oPos);
      } else if (fInteractionRate==kNeg) {
        fCRCnRun=50;
        fRunList=TArrayI(fCRCnRun,dRun15hNeg);
        fAvVtxPosX=TArrayD(fCRCnRun,dVtxPosX15oNeg);
        fAvVtxPosY=TArrayD(fCRCnRun,dVtxPosY15oNeg);
        fAvVtxPosZ=TArrayD(fCRCnRun,dVtxPosZ15oNeg);
      }
      fEnNucl=2511.;
      break;
    case k2015v6:
      fCRCnRun=91;
      fRunList=TArrayI(fCRCnRun,dRun15ov6);
      fEnNucl=2511.;
      break;
    case k2015pidfix:
      fCRCnRun=35;
      fRunList=TArrayI(fCRCnRun,dRun15opidfix);
      fAvVtxPosX=TArrayD(fCRCnRun,dVtxPosX15opidfix);
      fAvVtxPosY=TArrayD(fCRCnRun,dVtxPosY15opidfix);
      fAvVtxPosZ=TArrayD(fCRCnRun,dVtxPosZ15opidfix);
      fEnNucl=2511.;
      break;
  }

} // end of AliFlowAnalysisQvec::SetRunList()

//=======================================================================================================================

Int_t AliFlowAnalysisQvec::GetCRCRunBin(Int_t RunNum)
{
  Int_t CRCBin=-1,bin=0;
  for(Int_t c=0;c<fCRCnRun;c++) {
    if(fRunList[c]==RunNum) CRCBin=bin;
    else bin++;
  }
  return CRCBin;
} // end of AliFlowAnalysisQvec::GetCRCRunBin();

//=======================================================================================================================

void AliFlowAnalysisQvec::PassQAZDCCuts()
{
  // VZ eta < 0
  Double_t VZCM = fVZFlowVect[0][1].GetMult();
  // VZ eta > 0
  Double_t VZAM = fVZFlowVect[1][1].GetMult();
  // ZDCN-C
  Double_t QMC  = fZDCFlowVect[0].GetMult();
  // ZDCN-A
  Double_t QMA  = fZDCFlowVect[1].GetMult();

  // get re-centered QM*
  Double_t QMCrec = QMC;
  Double_t QMArec = QMA;
  if(fAvEZDCCRbRPro && fAvEZDCARbRPro) {
    Int_t runbin = fAvEZDCCRbRPro->GetXaxis()->FindBin(Form("%d",fRunNum));
    Int_t cenbin = fAvEZDCCRbRPro->GetYaxis()->FindBin(fCentralityEBE);
    QMCrec -= fAvEZDCCRbRPro->GetBinContent(runbin,cenbin);
    QMArec -= fAvEZDCARbRPro->GetBinContent(runbin,cenbin);
  }

  // cut on multiplicity
  if( fZNCen<=0. || fZNAen<=0. ) fQAZDCCutsFlag = kFALSE;
  // exclude ZDC bad runs
  if(fRunNum==138469 || fRunNum==138870 || fRunNum==139028 || fRunNum==139029 || fRunNum==139036) fQAZDCCutsFlag = kFALSE;
  Int_t CenBin = (Int_t)(fCentralityEBE);

  // cut on #neutrons
  if(fMinMulZN==1) {
    if(fZNCen > fPolMin[0]->Eval(fCentralityEBE) || fZNAen > fPolMin[1]->Eval(fCentralityEBE)) fQAZDCCutsFlag = kFALSE;
  }
  if(fMinMulZN==2) {
    if(fZNCen < fPolMax[0]->Eval(fCentralityEBE) || fZNAen < fPolMax[1]->Eval(fCentralityEBE)) fQAZDCCutsFlag = kFALSE;
  }

  // temporary mapping: Z*M = a*V0M + b, cut at cen. 40
  Double_t ZNM = (fZNCen+fZNAen)/2.;
  Double_t VZM = VZCM+VZAM;

  if(fMinMulZN==3 || fMinMulZN==4) {
    Double_t ZCMd = fZNCen;
    Double_t ZAMd = fZNAen;
    if(fCorrMap[CenBin]) {
      ZAMd += fCorrMap[CenBin];
      ZCMd -= fCorrMap[CenBin];
    }
    if(fMinMulZN==3) {
      if(ZCMd-ZAMd > fPolMin[0]->Eval(fCentralityEBE)) fQAZDCCutsFlag = kFALSE;
    }
    if(fMinMulZN==4) {
      if(ZCMd-ZAMd < fPolMax[0]->Eval(fCentralityEBE)) fQAZDCCutsFlag = kFALSE;
    }

    // new centrality
    if (VZM>5.E3) {
      Double_t par1 = 0.;
      if(VZM<1.E4) par1 = fPolSlope[0]->Eval(VZM);
      else         par1 = fPolSlope[1]->Eval(VZM);
      Double_t par0 = ZNM-par1*VZM;

      fPolDist[0]->SetParameter(0,par0);
      fPolDist[0]->SetParameter(1,par1);
      Double_t VZMint = fPolDist[0]->GetMinimumX();
      fNewMetricLEBE = sqrt(pow(VZMint-5.E3,2.)+pow(fPolAv[0]->Eval(VZMint)-fPolAv[0]->Eval(5.E3),2.));
      fNewMetricL2EBE = sqrt(pow((VZMint-5.E3)/1.E3,2.)+pow(fPolAv[0]->Eval(VZMint)-fPolAv[0]->Eval(5.E3),2.));

      fPolAv[1]->SetParameter(0,par0);
      fPolAv[1]->SetParameter(1,par1);
      fNewMetricDEBE = sqrt(pow(VZM-VZMint,2.)+pow(fPolAv[1]->Eval(VZM)-fPolAv[1]->Eval(VZMint),2.));
      fNewMetricD2EBE = sqrt(pow((VZM-VZMint)/1.E3,2.)+pow(fPolAv[1]->Eval(VZM)-fPolAv[1]->Eval(VZMint),2.));
      if(fPolAv[1]->Eval(VZM)<fPolAv[0]->Eval(VZM)) {
        fNewMetricDEBE *= -1.;
        fNewMetricD2EBE *= -1.;
      }

      fNewCentralityEBE = fCenMetric->Eval(fNewMetricLEBE);
    } else {
      fNewMetricLEBE = -1.;
      fNewMetricDEBE = -1.;
      fNewCentralityEBE = fCentralityEBE;
    }
  }

  if(fMinMulZN==7) {
    if (VZM>5.E3) {
      fQAZDCCutsFlag = kTRUE;
    } else {
      fQAZDCCutsFlag = kFALSE;
    }
  }

  // new ZDC ESE selection
  if(fMinMulZN==5) {
    fZDCESEclEbE=-1;
    Double_t DifMin=1.E3;
    Double_t ZNS = fZNCen+fZNAen;
    for (Int_t k=0; k<fZDCESEnPol; k++) {
      Double_t PolV = fPolCuts[k]->Eval(fCentralityEBE);
      if(ZNS<PolV && PolV-ZNS<DifMin) {
        fZDCESEclEbE=k;
        DifMin=PolV-ZNS;
      }
    }
    if (fZDCESEclEbE==-1) fZDCESEclEbE=fZDCESEnPol;
    if (fCentralityEBE>=90.) fQAZDCCutsFlag = kFALSE;
  }

  fhZNvsCen[0]->Fill(fCentralityEBE,fZNCen+fZNAen);

  // cut on ZDC spectra for LHC15o
  if(fMinMulZN==8) {
    if(fCentralityEBE<90.) {
      if(fZNCen+fZNAen<fEZNCutMin->GetBinContent(fEZNCutMin->FindBin(fCentralityEBE))) fQAZDCCutsFlag = kFALSE;
      if(fZNCen+fZNAen>fEZNCutMax->GetBinContent(fEZNCutMax->FindBin(fCentralityEBE))) fQAZDCCutsFlag = kFALSE;
    }
  }

  if(fMinMulZN==9) {
    if(fabs(fVtxPosCor[0])>4.25e-3) fQAZDCCutsFlag = kFALSE;
    if(fabs(fVtxPosCor[1])>3.9e-3) fQAZDCCutsFlag = kFALSE;
    if(fabs(fVtxPosCor[2])>5.) fQAZDCCutsFlag = kFALSE;
  }

  if(fMinMulZN>=10) {
    if(fRunNum==246087) fQAZDCCutsFlag = kFALSE;
  }

  if(fMinMulZN==11) {
    if(fabs(fVtxPosCor[0])>4.25e-3) fQAZDCCutsFlag = kFALSE;
    if(fabs(fVtxPosCor[1])>3.9e-3) fQAZDCCutsFlag = kFALSE;
    if(fabs(fVtxPosCor[2])>5.) fQAZDCCutsFlag = kFALSE;
  }

  if(fMinMulZN==12) {
    // cut outliers
    if(fabs(QMCrec)>100.) fQAZDCCutsFlag = kFALSE;
    if(fabs(QMArec)>100.) fQAZDCCutsFlag = kFALSE;
  }

  // fill QA plots
  if(fQAZDCCutsFlag) {
    fhZNCvsZNA[fCenBin]->Fill(fZNCen,fZNAen);
    fhZNvsCen[1]->Fill(fCentralityEBE,fZNCen+fZNAen);
    fhZNvsMul->Fill(fNITSCL1EBE,fZNCen+fZNAen,fCenWeightEbE);
  }
}

//=======================================================================================================================
void AliFlowAnalysisQvec::RecenterCRCQVecZDC2()
{
  if(!fZDCCalibListFinalCommonPart && !fZDCCalibListFinalRunByRun) {
    cout << " WARNING: no weights provided for ZDC recentering !!! " << endl;
    return;
  }
  // Read input histograms
  if(fRunNum!=fCachedRunNum) {
	 // Load Calib hist
	 fAvr_Run_CentQ[0] = (TProfile*)(fZDCCalibListFinalRunByRun->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fRun_CentQ[%d][%d]",fRunNum,0)));
	 fAvr_Run_CentQ[1] = (TProfile*)(fZDCCalibListFinalRunByRun->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fRun_CentQ[%d][%d]",fRunNum,1)));
	 fAvr_Run_CentQ[2] = (TProfile*)(fZDCCalibListFinalRunByRun->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fRun_CentQ[%d][%d]",fRunNum,2)));
	 fAvr_Run_CentQ[3] = (TProfile*)(fZDCCalibListFinalRunByRun->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fRun_CentQ[%d][%d]",fRunNum,3)));
     // 
     for(Int_t k=0; k<4; k++) {
        fAvr_Run_VtxXYZQ[k] = (TProfile3D*)(fZDCCalibListFinalRunByRun->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fRun_VtxXYZQ[%d][%d]",fRunNum,k)));
     }
     
     for(Int_t c=0; c<20; c++) {
	    for(Int_t k=0; k<4; k++) {
           fAvr_Cent_VtxXYZQ[c][k] = (TProfile3D*)(fZDCCalibListFinalCommonPart->FindObject(Form("fCent_VtxXYZQ[%d][%d]",c,k)));
        }
	 }
  }
  
  Int_t CentBin = Int_t(fCentralityEBE/(100/20)); // 0 to 19
  if (fCentralityEBE == 100) CentBin = 19;

  // ZDCN-C
  Double_t QCRe = fZDCFlowVect[0].X(); 
  Double_t QCIm = fZDCFlowVect[0].Y();
  Double_t QMC  = fZDCFlowVect[0].GetMult();
  // ZDCN-A
  Double_t QARe = fZDCFlowVect[1].X();
  Double_t QAIm = fZDCFlowVect[1].Y();
  Double_t QMA  = fZDCFlowVect[1].GetMult();
  
  Double_t QCReR=QCRe, QCImR=QCIm, QAReR=QARe, QAImR=QAIm;
  Double_t fillstep=0.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);
    
  // recenter vs centrality
  if (fAvr_Run_CentQ[0]) {
    Double_t AvQCRe = fAvr_Run_CentQ[0]->GetBinContent(fAvr_Run_CentQ[0]->FindBin(fCentralityEBE));
    Double_t SDQCRe = fAvr_Run_CentQ[0]->GetBinError(fAvr_Run_CentQ[0]->FindBin(fCentralityEBE));
    Double_t AvQCIm = fAvr_Run_CentQ[1]->GetBinContent(fAvr_Run_CentQ[1]->FindBin(fCentralityEBE));
    Double_t SDQCIm = fAvr_Run_CentQ[1]->GetBinError(fAvr_Run_CentQ[1]->FindBin(fCentralityEBE));

    Double_t AvQARe = fAvr_Run_CentQ[2]->GetBinContent(fAvr_Run_CentQ[2]->FindBin(fCentralityEBE));
    Double_t SDQARe = fAvr_Run_CentQ[2]->GetBinError(fAvr_Run_CentQ[2]->FindBin(fCentralityEBE));
    Double_t AvQAIm = fAvr_Run_CentQ[3]->GetBinContent(fAvr_Run_CentQ[3]->FindBin(fCentralityEBE));
    Double_t SDQAIm = fAvr_Run_CentQ[3]->GetBinError(fAvr_Run_CentQ[3]->FindBin(fCentralityEBE));
	//cout<<"===>Step 1=="<<endl;
    if(AvQCRe && AvQCIm && QMC>0. && sqrt(QCRe*QCRe+QCIm*QCIm)>1.E-6) {
      QCReR = QCRe-AvQCRe;
      QCImR = QCIm-AvQCIm;
      //cout<<"QCReR = "<<QCReR<<" = "<<QCRe<<" - "<<AvQCRe<<endl;
      //cout<<"QCImR = "<<QCImR<<" = "<<QCIm<<" - "<<AvQCIm<<endl;
      if(fDivSigma && SDQCRe>0. && SDQCIm>0.) {
        QCReR /= SDQCRe;
        QCImR /= SDQCIm;
      }
      fZDCFlowVect[0].Set(QCReR,QCImR);
    }

    if(AvQARe && AvQAIm && QMA>0. && sqrt(QARe*QARe+QAIm*QAIm)>1.E-6) {
      QAReR = QARe-AvQARe;
      QAImR = QAIm-AvQAIm;
      //cout<<"QAReR = "<<QAReR<<" = "<<QARe<<" - "<<AvQARe<<endl;
      //cout<<"QAImR = "<<QAImR<<" = "<<QAIm<<" - "<<AvQAIm<<endl;
      if(fDivSigma && SDQARe>0. && SDQAIm>0.) {
        QAReR /= SDQARe;
        QAImR /= SDQAIm;
      }
      fZDCFlowVect[1].Set(QAReR,QAImR);
    }
  }
  //cout<<"===>Step 2=="<<endl;
  fillstep=1.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);
  
  Bool_t withinvtx = kTRUE;
  if(fVtxPosCor15oIRSplit[0] < fAvr_Cent_VtxXYZQ[CentBin][0]->GetXaxis()->GetXmin() || fVtxPosCor15oIRSplit[0] > fAvr_Cent_VtxXYZQ[CentBin][0]->GetXaxis()->GetXmax()) withinvtx = kFALSE;
  if(fVtxPosCor15oIRSplit[1] < fAvr_Cent_VtxXYZQ[CentBin][0]->GetYaxis()->GetXmin() || fVtxPosCor15oIRSplit[1] > fAvr_Cent_VtxXYZQ[CentBin][0]->GetYaxis()->GetXmax()) withinvtx = kFALSE;
  if(fVtxPosCor15oIRSplit[2] < fAvr_Cent_VtxXYZQ[CentBin][0]->GetZaxis()->GetXmin() || fVtxPosCor15oIRSplit[2] > fAvr_Cent_VtxXYZQ[CentBin][0]->GetZaxis()->GetXmax()) withinvtx = kFALSE;
  
  if(fAvr_Cent_VtxXYZQ[CentBin][0]) {
    if(withinvtx) {
	  QCReR -= fAvr_Cent_VtxXYZQ[CentBin][0]->GetBinContent(fAvr_Cent_VtxXYZQ[CentBin][0]->FindBin(fVtxPosCor15oIRSplit[0],fVtxPosCor15oIRSplit[1],fVtxPosCor15oIRSplit[2]));
	  QCImR -= fAvr_Cent_VtxXYZQ[CentBin][1]->GetBinContent(fAvr_Cent_VtxXYZQ[CentBin][1]->FindBin(fVtxPosCor15oIRSplit[0],fVtxPosCor15oIRSplit[1],fVtxPosCor15oIRSplit[2]));

	  QAReR -= fAvr_Cent_VtxXYZQ[CentBin][2]->GetBinContent(fAvr_Cent_VtxXYZQ[CentBin][2]->FindBin(fVtxPosCor15oIRSplit[0],fVtxPosCor15oIRSplit[1],fVtxPosCor15oIRSplit[2]));
	  QAImR -= fAvr_Cent_VtxXYZQ[CentBin][3]->GetBinContent(fAvr_Cent_VtxXYZQ[CentBin][3]->FindBin(fVtxPosCor15oIRSplit[0],fVtxPosCor15oIRSplit[1],fVtxPosCor15oIRSplit[2]));
    } else {
	  Double_t vx = fVtxPosCor15oIRSplit[0];
	  Double_t vy = fVtxPosCor15oIRSplit[1];
	  Double_t vz = fVtxPosCor15oIRSplit[2];
	  
	  if(fVtxPosCor15oIRSplit[0] < fAvr_Cent_VtxXYZQ[CentBin][0]->GetXaxis()->GetXmin()) vx = fAvr_Cent_VtxXYZQ[CentBin][0]->GetXaxis()->GetBinCenter(1);
	  if(fVtxPosCor15oIRSplit[0] > fAvr_Cent_VtxXYZQ[CentBin][0]->GetXaxis()->GetXmax()) vx = fAvr_Cent_VtxXYZQ[CentBin][0]->GetXaxis()->GetBinCenter(fAvr_Cent_VtxXYZQ[CentBin][0]->GetNbinsX());
	  if(fVtxPosCor15oIRSplit[1] < fAvr_Cent_VtxXYZQ[CentBin][0]->GetYaxis()->GetXmin()) vy = fAvr_Cent_VtxXYZQ[CentBin][0]->GetYaxis()->GetBinCenter(1);
	  if(fVtxPosCor15oIRSplit[1] > fAvr_Cent_VtxXYZQ[CentBin][0]->GetYaxis()->GetXmax()) vy = fAvr_Cent_VtxXYZQ[CentBin][0]->GetYaxis()->GetBinCenter(fAvr_Cent_VtxXYZQ[CentBin][0]->GetNbinsY());
	  if(fVtxPosCor15oIRSplit[2] < fAvr_Cent_VtxXYZQ[CentBin][0]->GetYaxis()->GetXmin()) vz = fAvr_Cent_VtxXYZQ[CentBin][0]->GetZaxis()->GetBinCenter(1);
	  if(fVtxPosCor15oIRSplit[2] > fAvr_Cent_VtxXYZQ[CentBin][0]->GetYaxis()->GetXmax()) vz = fAvr_Cent_VtxXYZQ[CentBin][0]->GetZaxis()->GetBinCenter(fAvr_Cent_VtxXYZQ[CentBin][0]->GetNbinsZ());
  
	  QCReR -= fAvr_Cent_VtxXYZQ[CentBin][0]->GetBinContent(fAvr_Cent_VtxXYZQ[CentBin][0]->FindBin(vx,vy,vz));
	  QCImR -= fAvr_Cent_VtxXYZQ[CentBin][1]->GetBinContent(fAvr_Cent_VtxXYZQ[CentBin][1]->FindBin(vx,vy,vz));
	  QAReR -= fAvr_Cent_VtxXYZQ[CentBin][2]->GetBinContent(fAvr_Cent_VtxXYZQ[CentBin][2]->FindBin(vx,vy,vz));
	  QAImR -= fAvr_Cent_VtxXYZQ[CentBin][3]->GetBinContent(fAvr_Cent_VtxXYZQ[CentBin][3]->FindBin(vx,vy,vz));
    }
  }
  fZDCFlowVect[0].Set(QCReR,QCImR);
  fZDCFlowVect[1].Set(QAReR,QAImR);
  //cout<<"===>Step 3=="<<endl;
  fillstep=2.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);
  
  // correct vs vtx vs run number
  // check if possible to interpolate
  Bool_t bInterp = kTRUE;
  Int_t bx = fAvr_Run_VtxXYZQ[0]->GetXaxis()->FindBin(fVtxPosCor15oIRSplit[0]);
  Int_t by = fAvr_Run_VtxXYZQ[0]->GetYaxis()->FindBin(fVtxPosCor15oIRSplit[1]);
  Int_t bz = fAvr_Run_VtxXYZQ[0]->GetZaxis()->FindBin(fVtxPosCor15oIRSplit[2]);
  if(bx==1 || bx==fAvr_Run_VtxXYZQ[0]->GetXaxis()->GetNbins()) bInterp = kFALSE;
  if(by==1 || by==fAvr_Run_VtxXYZQ[0]->GetYaxis()->GetNbins()) bInterp = kFALSE;
  if(bz==1 || bz==fAvr_Run_VtxXYZQ[0]->GetZaxis()->GetNbins()) bInterp = kFALSE;
  
  if(fAvr_Run_VtxXYZQ[0]) {
    if(bInterp) {
	  QCReR -= fAvr_Run_VtxXYZQ[0]->Interpolate(fVtxPosCor15oIRSplit[0],fVtxPosCor15oIRSplit[1],fVtxPosCor15oIRSplit[2]);
	  QCImR -= fAvr_Run_VtxXYZQ[1]->Interpolate(fVtxPosCor15oIRSplit[0],fVtxPosCor15oIRSplit[1],fVtxPosCor15oIRSplit[2]);
	  QAReR -= fAvr_Run_VtxXYZQ[2]->Interpolate(fVtxPosCor15oIRSplit[0],fVtxPosCor15oIRSplit[1],fVtxPosCor15oIRSplit[2]);
	  QAImR -= fAvr_Run_VtxXYZQ[3]->Interpolate(fVtxPosCor15oIRSplit[0],fVtxPosCor15oIRSplit[1],fVtxPosCor15oIRSplit[2]);
    } else {
	  QCReR -= fAvr_Run_VtxXYZQ[0]->GetBinContent(bx,by,bz);
	  QCImR -= fAvr_Run_VtxXYZQ[1]->GetBinContent(bx,by,bz);
	  QAReR -= fAvr_Run_VtxXYZQ[2]->GetBinContent(bx,by,bz);
	  QAImR -= fAvr_Run_VtxXYZQ[3]->GetBinContent(bx,by,bz);
    }
  }
  fZDCFlowVect[0].Set(QCReR,QCImR);
  fZDCFlowVect[1].Set(QAReR,QAImR);

  fillstep=3.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);
  
  // if correctly set (via fMinMulZN), cut on ZDC Q-vector distribution
  // WARNING: QARe inverted !!!
  Bool_t pass = kTRUE;
  pass = PassCutZDCQVecDis(QCReR,QCImR,-QAReR,QAImR);
  // ***************************************************************************
  // store results after correction
  // ***************************************************************************

  Int_t bw = (fbFlagIsPosMagField==kTRUE?0:1);
  if(QMC>0. && QMA>0. && sqrt(QCReR*QCReR+QCImR*QCImR)>1.E-6 && sqrt(QAReR*QAReR+QAImR*QAImR)>1.E-6 && pass) {
    if( fInvertZDC ) QAReR = -QAReR;
  
    fCRCZDCQVecRes[fRunBin][0]->Fill(fCentralityEBE,QCReR*QAReR);
    fCRCZDCQVecRes[fRunBin][1]->Fill(fCentralityEBE,QCImR*QAImR);
    fCRCZDCQVecRes[fRunBin][2]->Fill(fCentralityEBE,QCReR*QAImR);
    fCRCZDCQVecRes[fRunBin][3]->Fill(fCentralityEBE,QCImR*QAReR);

    if(fStoreZDCQVecVtxPos) {
      fCRCZDCQVecDis[bw][fCenBin][0]->Fill(QCReR*QAImR,QCImR*QAReR);
      fCRCZDCQVecDis[bw][fCenBin][1]->Fill(QCReR*QAReR,QCImR*QAImR);
    }

    Double_t EvPlZNC = TMath::ATan2(QCImR,QCReR);
    if(EvPlZNC<0.) EvPlZNC += TMath::TwoPi();
    Double_t EvPlZNA = TMath::ATan2(QAImR,QAReR);
    if(EvPlZNA<0.) EvPlZNA += TMath::TwoPi();
    fZDCEPHist[0]->Fill(fCentralityEBE,EvPlZNC);
    fZDCEPHist[1]->Fill(fCentralityEBE,EvPlZNA);

  } else {
    fQAZDCCutsFlag = kFALSE;
  }
}
//=======================================================================================================================

void AliFlowAnalysisQvec::RecenterCRCQVecZDC()
{
  if(!fCRCZDCCalibList) {
    cout << " WARNING: no weights provided for ZDC recentering !!! " << endl;
    return;
  }

  Int_t qb[4] = {0};
  if(fbFlagIsPosMagField) { qb[0]=0; qb[1]=1; qb[2]=4; qb[3]=5; }
  else                    { qb[0]=2; qb[1]=3; qb[2]=6; qb[3]=7; }
  Int_t bw = (fbFlagIsPosMagField==kTRUE?0:1);

  if(fRunNum!=fCachedRunNum) {
    fZDCQHist[0] = (TProfile*)(fCRCZDCCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCZDCQVecC[%d][%d]",fRunNum,0)));
    fZDCQHist[1] = (TProfile*)(fCRCZDCCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCZDCQVecC[%d][%d]",fRunNum,1)));
    fZDCQHist[2] = (TProfile*)(fCRCZDCCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCZDCQVecA[%d][%d]",fRunNum,0)));
    fZDCQHist[3] = (TProfile*)(fCRCZDCCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCZDCQVecA[%d][%d]",fRunNum,1)));
    fZDCQHist[4] = (TProfile*)(fCRCZDCCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCZDCQVecCCorr[%d][%d]",fRunNum,0)));
    fZDCQHist[5] = (TProfile*)(fCRCZDCCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCZDCQVecCCorr[%d][%d]",fRunNum,1)));
    fZDCQHist[6] = (TProfile*)(fCRCZDCCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCZDCQVecACorr[%d][%d]",fRunNum,0)));
    fZDCQHist[7] = (TProfile*)(fCRCZDCCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCZDCQVecACorr[%d][%d]",fRunNum,1)));
    fZDCQHist[8] = (TProfile*)(fCRCZDCCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCZDCQVecCAdd[%d][%d]",fRunNum,0))); // @Shi not exist for 15oHi
    fZDCQHist[9] = (TProfile*)(fCRCZDCCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCZDCQVecCAdd[%d][%d]",fRunNum,1))); // @Shi not exist for 15oHi
    fZDCQHist[10] = (TProfile*)(fCRCZDCCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCZDCQVecAAdd[%d][%d]",fRunNum,0))); // @Shi not exist for 15oHi
    fZDCQHist[11] = (TProfile*)(fCRCZDCCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCZDCQVecAAdd[%d][%d]",fRunNum,1))); // @Shi not exist for 15oHi

    for(Int_t k=0; k<4; k++) {
      fZDCVtxHist[k] = (TProfile3D*)(fCRCZDCCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCZDCQVecVtxPos540[%d][%d]",fRunNum,k))); // @Shi not exist for 15oHi

      fZDCEcomHist[k] = (TProfile2D*)(fCRCZDCCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCZDCQVecECom[%d][%d]",fRunNum,k))); // @Shi not exist for 15oHi
      fZDCEcomTotHist[k] = (TProfile2D*)(fCRCZDCCalibList->FindObject(Form("fCRCZDCQVecEComTot[%d]",k))); // @Shi not exist for 15oHi
    }
    for(Int_t k=0; k<12; k++) {
      fZDCEcomTotvsVtxHist[k] = (TProfile3D*)(fCRCZDCCalibList->FindObject(Form("fCRCZDCQVecEComTotVtx[%d]",k))); // @Shi not exist for 15oHi
    }

    for(Int_t c=0; c<10; c++) {
      for(Int_t k=0; k<4; k++) {
        fZDCVtxCenHist[c][k] = (TProfile3D*)(fCRCZDCCalibList->FindObject(Form("fCRCZDCQVecVtxPosCen[%d][%d]",c,k))); // @Shi not exist for 15oHi
      }
    }

    if(fCRCZDCResList) {
      // ZNC
      fZDCResHist[0] = (TH1D*)(fCRCZDCResList->FindObject(Form("fZNPhase[%d][%d]",fRunNum,0))); // @Shi fZDCResHist[] is not set in AddTask
      fZDCResHist[1] = (TH1D*)(fCRCZDCResList->FindObject(Form("fZNRatio[%d][%d]",fRunNum,0)));
      // ZNA
      fZDCResHist[2] = (TH1D*)(fCRCZDCResList->FindObject(Form("fZNPhase[%d][%d]",fRunNum,1)));
      fZDCResHist[3] = (TH1D*)(fCRCZDCResList->FindObject(Form("fZNRatio[%d][%d]",fRunNum,1)));
    }

     for(Int_t k=0; k<4; k++) {
       fZDCVtxFitHist[k] = (TH3D*)fCRCZDCCalibList->FindObject(Form("TH3SlopeRunCenVtx[%d]",k)); // @Shi do exist for 15oHi
       if(fZDCVtxFitHist[k]) {
         fZDCVtxFitHist[k]->Sumw2(kFALSE);
         Int_t runbin = fZDCVtxFitHist[k]->GetXaxis()->FindBin(Form("%d",fRunNum));
         fZDCVtxFitHist[k]->GetXaxis()->SetRange(runbin,runbin);
         for(Int_t i=0; i<3; i++) {
           fZDCVtxFitHist[k]->GetZaxis()->SetRange(i+1,i+1);
           fZDCVtxFitCenProjHist[k][i] = (TH1D*)fZDCVtxFitHist[k]->Project3D("y")->Clone(Form("proj[%d][%d]",k,i));
         }
       }
     }
    for(Int_t k=0; k<4; k++) {
      fZDCVtxFitHist2[k] = (TH3D*)fCRCZDCCalibList->FindObject(Form("TH3SlopePol3RunVtx[%d]",k)); // @Shi do exist for 15oHi
      if(fZDCVtxFitHist2[k]) {
        fZDCVtxFitHist2[k]->Sumw2(kFALSE);
        Int_t runbin = fZDCVtxFitHist2[k]->GetXaxis()->FindBin(Form("%d",fRunNum));
        fZDCVtxFitHist2[k]->GetXaxis()->SetRange(runbin,runbin);
        for(Int_t i=0; i<3; i++) {
          fZDCVtxFitHist2[k]->GetZaxis()->SetRange(i+1,i+1);
          fZDCVtxFitCenProjHist2[k][i] = (TH1D*)fZDCVtxFitHist2[k]->Project3D("y")->Clone(Form("proj2[%d][%d]",k,i));
        }
      }
    }
    for(Int_t c=0; c<10; c++) {
      for(Int_t k=0; k<8; k++) {
        fZDCVtxCenHistMagPol[c][k] = (TProfile3D*)(fCRCZDCCalibList->FindObject(Form("fZDCVtxCenHistMagPol[%d][%d]",c,k))); // @Shi do exist for 15oHi
      }
    }
    for(Int_t c=0; c<10; c++) {
      fZDCBinsCenRefMult[c] = (TH3D*)(fCRCZDCCalibList->FindObject(Form("ZDCQVecCenRefMul[%d]",c))); // @Shi do exist for 15oHi
    }
    for(Int_t k=0; k<4; k++) {
      fZDCBinsCenRefMultRbR[k] = (TProfile2D*)(fCRCZDCCalibList->FindObject(Form("Run %d",fRunNum))->FindObject(Form("fCRCZDCQVecCenRefMul[%d][%d]",fRunNum,k))); // @Shi do exist for 15oHi
      if(fZDCBinsCenRefMultRbR[k]) {
        for(Int_t c=0; c<10; c++) {
          fZDCBinsCenRefMultRbRProf[c][k] = (TProfile*)fZDCBinsCenRefMultRbR[k]->ProfileY(Form("posrbr%d%d",c,k),c+1,c+1)->Clone(Form("posrbr%d%d",c,k));
          Int_t nFullBin=0;
          for (Int_t bx=0; bx<fZDCBinsCenRefMultRbRProf[c][k]->GetNbinsX(); bx++) {
            if(fabs(fZDCBinsCenRefMultRbRProf[c][k]->GetBinContent(bx))>0.) nFullBin++;
          }
          if(nFullBin>10) {
            Int_t rebn = (Int_t)nFullBin/10.;
            fZDCBinsCenRefMultRbRProf[c][k]->Rebin(rebn+1);
          }
          fZDCBinsCenRefMultRbRProj[c][k] = (TH1D*)fZDCBinsCenRefMultRbRProf[c][k]->Clone(Form("poshisrbr%d%d",c,k));
        }
      }
    }
    for(Int_t k=0; k<4; k++) {
      fZDCBinsCenRefMultTot[k] = (TProfile2D*)(fCRCZDCCalibList->FindObject(Form("fCRCZDCQVecCenRefMulTot[%d][%d]",bw,k))); // @Shi do exist for 15oHi
      if(fZDCBinsCenRefMultTot[k]) {
        for(Int_t c=0; c<10; c++) {
          fZDCBinsCenRefMultTotProf[c][k] = (TProfile*)fZDCBinsCenRefMultTot[k]->ProfileY(Form("postot%d%d",c,k),c+1,c+1)->Clone(Form("postot%d%d",c,k));
          Int_t nFullBin=0;
          for (Int_t bx=0; bx<fZDCBinsCenRefMultTotProf[c][k]->GetNbinsX(); bx++) {
            if(fabs(fZDCBinsCenRefMultTotProf[c][k]->GetBinContent(bx))>0.) nFullBin++;
          }
          if(nFullBin>10) {
            Int_t rebn = (Int_t)nFullBin/10.;
            fZDCBinsCenRefMultTotProf[c][k]->Rebin(rebn+1);
          }
          fZDCBinsCenRefMultTotProj[c][k] = (TH1D*)fZDCBinsCenRefMultTotProf[c][k]->Clone(Form("poshistot%d%d",c,k));
        }
      }
    }

    for(Int_t i=0; i<3; i++) {
      for(Int_t k=0; k<4; k++) {
        fZDCBinsVtxCenEZDC[i][k] = (TH3D*)(fCRCZDCCalibList->FindObject(Form("ZDCQVecVtxCenEZDC[%d][%d]",i,qb[k]))); // @Shi not exist for 15oHi
      }
    }
    fZDCQVecVtxCenEZDCFit0 = (TH3D*)(fCRCZDCCalibList->FindObject("ZDCQVecVtxCenEZDCFit0"));
    fZDCQVecVtxCenEZDCFit1 = (TH3D*)(fCRCZDCCalibList->FindObject("ZDCQVecVtxCenEZDCFit1"));

    for(Int_t i=0; i<10; i++) {
      for(Int_t z=0; z<10; z++) {
        for(Int_t k=0; k<4; k++) {
          fZDCQVecVtxCenEZDC3D[i][z][k] = (TH3D*)(fCRCZDCCalibList->FindObject(Form("ZDCQVecVtxCenEZDC3D[%d][%d][%d]",i,z,qb[k])));
        }
      }
    }

  }

  if(!fQAZDCCutsFlag) return;

  // ZDCN-C
  Double_t QCRe = fZDCFlowVect[0].X(); 
  Double_t QCIm = fZDCFlowVect[0].Y();
  Double_t QMC  = fZDCFlowVect[0].GetMult();
  // ZDCN-A
  Double_t QARe = fZDCFlowVect[1].X();
  Double_t QAIm = fZDCFlowVect[1].Y();
  Double_t QMA  = fZDCFlowVect[1].GetMult();

  // get re-centered QM*
  Double_t QMCrec = QMC;
  Double_t QMArec = QMA;
  if(fAvEZDCCRbRPro && fAvEZDCARbRPro) {
    Int_t runbin = fAvEZDCCRbRPro->GetXaxis()->FindBin(Form("%d",fRunNum));
    Int_t cenbin = fAvEZDCCRbRPro->GetYaxis()->FindBin(fCentralityEBE);
    QMCrec -= fAvEZDCCRbRPro->GetBinContent(runbin,cenbin);
    QMArec -= fAvEZDCARbRPro->GetBinContent(runbin,cenbin);
  }

  Double_t QCReR=QCRe, QCImR=QCIm, QAReR=QARe, QAImR=QAIm;

  Double_t fillstep=0.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);

  // recenter vs centrality
  if (fZDCQHist[0]) {
    Double_t AvQCRe = fZDCQHist[0]->GetBinContent(fZDCQHist[0]->FindBin(fCentralityEBE));
    Double_t SDQCRe = fZDCQHist[0]->GetBinError(fZDCQHist[0]->FindBin(fCentralityEBE));
    Double_t AvQCIm = fZDCQHist[1]->GetBinContent(fZDCQHist[1]->FindBin(fCentralityEBE));
    Double_t SDQCIm = fZDCQHist[1]->GetBinError(fZDCQHist[1]->FindBin(fCentralityEBE));

    Double_t AvQARe = fZDCQHist[2]->GetBinContent(fZDCQHist[2]->FindBin(fCentralityEBE));
    Double_t SDQARe = fZDCQHist[2]->GetBinError(fZDCQHist[2]->FindBin(fCentralityEBE));
    Double_t AvQAIm = fZDCQHist[3]->GetBinContent(fZDCQHist[3]->FindBin(fCentralityEBE));
    Double_t SDQAIm = fZDCQHist[3]->GetBinError(fZDCQHist[3]->FindBin(fCentralityEBE));

    if(AvQCRe && AvQCIm && QMC>0. && sqrt(QCRe*QCRe+QCIm*QCIm)>1.E-6) {
      QCReR = QCRe-AvQCRe;
      QCImR = QCIm-AvQCIm;
      if(fDivSigma && SDQCRe>0. && SDQCIm>0.) {
        QCReR /= SDQCRe;
        QCImR /= SDQCIm;
      }
      fZDCFlowVect[0].Set(QCReR,QCImR);
    }

    if(AvQARe && AvQAIm && QMA>0. && sqrt(QARe*QARe+QAIm*QAIm)>1.E-6) {
      QAReR = QARe-AvQARe;
      QAImR = QAIm-AvQAIm;
      if(fDivSigma && SDQARe>0. && SDQAIm>0.) {
        QAReR /= SDQARe;
        QAImR /= SDQAIm;
      }
      fZDCFlowVect[1].Set(QAReR,QAImR);
    }
  }

  fillstep=1.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);

  Bool_t pass = kTRUE;

  if(!fZDCVtxCenHist[fCenBin][0]) {
//    if(fCentralityEBE>5. && fCentralityEBE<60.) {
//      if(QMC>0. && sqrt(QCRe*QCRe+QCIm*QCIm)>1.E-6) {
//        fCRCZDCQVecVtxPos[fRunBin][0]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QCReR);
//        fCRCZDCQVecVtxPos[fRunBin][1]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QCImR);
//      }
//      if(QMA>0. && sqrt(QARe*QARe+QAIm*QAIm)>1.E-6) {
//        fCRCZDCQVecVtxPos[fRunBin][2]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QAReR);
//        fCRCZDCQVecVtxPos[fRunBin][3]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QAImR);
//      }
//    }
  }

  if (fZDCVtxCenHist[fCenBin][0]) {
    Bool_t withinvtx = kTRUE;
      // method #2: recenter vs vtx in centrality bins
    if(fVtxPosCor[0] < fZDCVtxCenHist[fCenBin][0]->GetXaxis()->GetXmin() || fVtxPosCor[0] > fZDCVtxCenHist[fCenBin][0]->GetXaxis()->GetXmax()) withinvtx = kFALSE;
    if(fVtxPosCor[1] < fZDCVtxCenHist[fCenBin][0]->GetYaxis()->GetXmin() || fVtxPosCor[1] > fZDCVtxCenHist[fCenBin][0]->GetYaxis()->GetXmax()) withinvtx = kFALSE;
    if(fVtxPosCor[2] < fZDCVtxCenHist[fCenBin][0]->GetZaxis()->GetXmin() || fVtxPosCor[2] > fZDCVtxCenHist[fCenBin][0]->GetZaxis()->GetXmax()) withinvtx = kFALSE;

    if(withinvtx) {
      QCReR -= fZDCVtxCenHist[fCenBin][0]->GetBinContent(fZDCVtxCenHist[fCenBin][0]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      QCImR -= fZDCVtxCenHist[fCenBin][1]->GetBinContent(fZDCVtxCenHist[fCenBin][1]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      fZDCFlowVect[0].Set(QCReR,QCImR);
      QAReR -= fZDCVtxCenHist[fCenBin][2]->GetBinContent(fZDCVtxCenHist[fCenBin][2]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      QAImR -= fZDCVtxCenHist[fCenBin][3]->GetBinContent(fZDCVtxCenHist[fCenBin][3]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      fZDCFlowVect[1].Set(QAReR,QAImR);
    } else {
      Double_t vx = fVtxPosCor[0];
      Double_t vy = fVtxPosCor[1];
      Double_t vz = fVtxPosCor[2];
      if(fVtxPosCor[0] < fZDCVtxCenHist[fCenBin][0]->GetXaxis()->GetXmin()) vx = fZDCVtxCenHist[fCenBin][0]->GetXaxis()->GetBinCenter(1);
      if(fVtxPosCor[0] > fZDCVtxCenHist[fCenBin][0]->GetXaxis()->GetXmax()) vx = fZDCVtxCenHist[fCenBin][0]->GetXaxis()->GetBinCenter(fZDCVtxCenHist[fCenBin][0]->GetNbinsX());
      if(fVtxPosCor[1] < fZDCVtxCenHist[fCenBin][0]->GetYaxis()->GetXmin()) vy = fZDCVtxCenHist[fCenBin][0]->GetYaxis()->GetBinCenter(1);
      if(fVtxPosCor[1] > fZDCVtxCenHist[fCenBin][0]->GetYaxis()->GetXmax()) vy = fZDCVtxCenHist[fCenBin][0]->GetYaxis()->GetBinCenter(fZDCVtxCenHist[fCenBin][0]->GetNbinsY());
      if(fVtxPosCor[2] < fZDCVtxCenHist[fCenBin][0]->GetYaxis()->GetXmin()) vz = fZDCVtxCenHist[fCenBin][0]->GetZaxis()->GetBinCenter(1);
      if(fVtxPosCor[2] > fZDCVtxCenHist[fCenBin][0]->GetYaxis()->GetXmax()) vz = fZDCVtxCenHist[fCenBin][0]->GetZaxis()->GetBinCenter(fZDCVtxCenHist[fCenBin][0]->GetNbinsZ());
      QCReR -= fZDCVtxCenHist[fCenBin][0]->GetBinContent(fZDCVtxCenHist[fCenBin][0]->FindBin(vx,vy,vz));
      QCImR -= fZDCVtxCenHist[fCenBin][1]->GetBinContent(fZDCVtxCenHist[fCenBin][1]->FindBin(vx,vy,vz));
      fZDCFlowVect[0].Set(QCReR,QCImR);
      QAReR -= fZDCVtxCenHist[fCenBin][2]->GetBinContent(fZDCVtxCenHist[fCenBin][2]->FindBin(vx,vy,vz));
      QAImR -= fZDCVtxCenHist[fCenBin][3]->GetBinContent(fZDCVtxCenHist[fCenBin][3]->FindBin(vx,vy,vz));
      fZDCFlowVect[1].Set(QAReR,QAImR);
    }
  }

  fillstep=2.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);

  // intermediate step: save vtx pos run-by-run
  //    if(fCentralityEBE>5. && fCentralityEBE<60.) {
  //      if(QMC>0. && sqrt(QCRe*QCRe+QCIm*QCIm)>1.E-6) {
  //        fCRCZDCQVecVtxPos[fRunBin][0]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QCReR);
  //        fCRCZDCQVecVtxPos[fRunBin][1]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QCImR);
  //      }
  //      if(QMA>0. && sqrt(QARe*QARe+QAIm*QAIm)>1.E-6) {
  //        fCRCZDCQVecVtxPos[fRunBin][2]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QAReR);
  //        fCRCZDCQVecVtxPos[fRunBin][3]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QAImR);
  //      }
  //    }

  // correct vs vtx vs run number in centrality 5-40% (only for 15opidfix)
  if (fZDCVtxHist[0]) {
    // check if possible to interpolate
    Bool_t bInterp = kTRUE;
    Int_t bx = fZDCVtxHist[0]->GetXaxis()->FindBin(fVtxPosCor[0]);
    Int_t by = fZDCVtxHist[0]->GetYaxis()->FindBin(fVtxPosCor[1]);
    Int_t bz = fZDCVtxHist[0]->GetZaxis()->FindBin(fVtxPosCor[2]);
    if(bx==1 || bx==fZDCVtxHist[0]->GetXaxis()->GetNbins()) bInterp = kFALSE;
    if(by==1 || by==fZDCVtxHist[0]->GetYaxis()->GetNbins()) bInterp = kFALSE;
    if(bz==1 || bz==fZDCVtxHist[0]->GetZaxis()->GetNbins()) bInterp = kFALSE;
    if(bInterp) {
      QCReR -= fZDCVtxHist[0]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
      QCImR -= fZDCVtxHist[1]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
      QAReR -= fZDCVtxHist[2]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
      QAImR -= fZDCVtxHist[3]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
    } else {
      QCReR -= fZDCVtxHist[0]->GetBinContent(bx,by,bz);
      QCImR -= fZDCVtxHist[1]->GetBinContent(bx,by,bz);
      QAReR -= fZDCVtxHist[2]->GetBinContent(bx,by,bz);
      QAImR -= fZDCVtxHist[3]->GetBinContent(bx,by,bz);
    }
    fZDCFlowVect[0].Set(QCReR,QCImR);
    fZDCFlowVect[1].Set(QAReR,QAImR);
  }

  fillstep=3.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);

  //    Int_t EPCenBin = fCRCZDCQVecEP[fRunBin][0]->GetXaxis()->FindBin(fCentralityEBE)-1;
  //    if (fZDCEPHist[EPCenBin][0]) {
  //      Double_t EvPlZDCC = TMath::ATan2(QCImR,QCReR);
  //      fZDCEPweightEbE[0] = fZDCEPHist[EPCenBin][0]->GetBinContent(fZDCEPHist[EPCenBin][0]->FindBin(EvPlZDCC));
  //      Double_t EvPlZDCA = TMath::ATan2(QAImR,QAReR);
  //      fZDCEPweightEbE[1] = fZDCEPHist[EPCenBin][1]->GetBinContent(fZDCEPHist[EPCenBin][1]->FindBin(EvPlZDCA));
  //      Double_t EvPlZDCfull = TMath::ATan2(QAImR-QCImR,-QAReR-QCReR); // WARNING: fInvertZDC implicit
  //      fZDCEPweightEbE[2] = fZDCEPHist[EPCenBin][2]->GetBinContent(fZDCEPHist[EPCenBin][2]->FindBin(EvPlZDCfull));
  //    }
  //    for (Int_t k=0; k<3; k++) {
  //      if(fZDCEPweightEbE[k]==0) fZDCEPweightEbE[k]=1.;
  //    }

  //  // store cos/sin terms of event planes
  //  if(QMC>0. && sqrt(QCRe*QCRe+QCIm*QCIm)>1.E-6) {
  //    Double_t EvPlZDCC = TMath::ATan2(QCImR,QCReR);
  //    fCRCZDCQVecRes[fRunBin][4]->Fill(fCentralityEBE,cos(2.*EvPlZDCC));
  //    fCRCZDCQVecRes[fRunBin][5]->Fill(fCentralityEBE,sin(2.*EvPlZDCC));
  //  }
  //  if(QMA>0. && sqrt(QARe*QARe+QAIm*QAIm)>1.E-6) {
  //    Double_t EvPlZDCA = TMath::ATan2(QAImR,QAReR);
  //    fCRCZDCQVecRes[fRunBin][6]->Fill(fCentralityEBE,cos(2.*EvPlZDCA));
  //    fCRCZDCQVecRes[fRunBin][7]->Fill(fCentralityEBE,sin(2.*EvPlZDCA));
  //  }

  // rotate, correct, rotate back
  if(fZDCResHist[0] && fZDCResHist[1] && fZDCResHist[2] && fZDCResHist[3]) {
    // ZNA
    Double_t thetaC = fZDCResHist[0]->GetBinContent(fZDCResHist[0]->FindBin(fCentralityEBE));
    Double_t corrfC = fZDCResHist[1]->GetBinContent(fZDCResHist[1]->FindBin(fCentralityEBE));
    // ZNC
    Double_t thetaA = fZDCResHist[2]->GetBinContent(fZDCResHist[2]->FindBin(fCentralityEBE));
    Double_t corrfA = fZDCResHist[3]->GetBinContent(fZDCResHist[3]->FindBin(fCentralityEBE));

    if(corrfC>0. && corrfA>0.) {
      Double_t QCReRt = TMath::Cos(thetaC)*QCReR - TMath::Sin(thetaC)*QCImR;
      Double_t QCImRt = TMath::Sin(thetaC)*QCReR + TMath::Cos(thetaC)*QCImR;
      Double_t QAReRt = TMath::Cos(thetaA)*QAReR - TMath::Sin(thetaA)*QAImR;
      Double_t QAImRt = TMath::Sin(thetaA)*QAReR + TMath::Cos(thetaA)*QAImR;

      QCReRt /= corrfC;
      QAReRt /= corrfA;

      QCReR = TMath::Cos(-thetaC)*QCReRt - TMath::Sin(-thetaC)*QCImRt;
      QCImR = TMath::Sin(-thetaC)*QCReRt + TMath::Cos(-thetaC)*QCImRt;
      QAReR = TMath::Cos(-thetaA)*QAReRt - TMath::Sin(-thetaA)*QAImRt;
      QAImR = TMath::Sin(-thetaA)*QAReRt + TMath::Cos(-thetaA)*QAImRt;
    }
  }

  // if possible, correct Q-vectors vs Energy in the common tower
//  if(fZDCEcomHist[0]) {
//    QCReR -= fZDCEcomHist[0]->GetBinContent(fZDCEcomHist[0]->FindBin(fCentralityEBE,fZNCQ0));
//    QCImR -= fZDCEcomHist[1]->GetBinContent(fZDCEcomHist[1]->FindBin(fCentralityEBE,fZNCQ0));
//    fZDCFlowVect[0].Set(QCReR,QCImR);
//    QAReR -= fZDCEcomHist[2]->GetBinContent(fZDCEcomHist[2]->FindBin(fCentralityEBE,fZNAQ0));
//    QAImR -= fZDCEcomHist[3]->GetBinContent(fZDCEcomHist[3]->FindBin(fCentralityEBE,fZNAQ0));
//    fZDCFlowVect[1].Set(QAReR,QAImR);
//  }
  if(fZDCEcomTotHist[0]) {
    QCReR -= fZDCEcomTotHist[0]->GetBinContent(fZDCEcomTotHist[0]->FindBin(fCentralityEBE,fZNCQ0));
    QCImR -= fZDCEcomTotHist[1]->GetBinContent(fZDCEcomTotHist[1]->FindBin(fCentralityEBE,fZNCQ0));
    fZDCFlowVect[0].Set(QCReR,QCImR);
    QAReR -= fZDCEcomTotHist[2]->GetBinContent(fZDCEcomTotHist[2]->FindBin(fCentralityEBE,fZNAQ0));
    QAImR -= fZDCEcomTotHist[3]->GetBinContent(fZDCEcomTotHist[3]->FindBin(fCentralityEBE,fZNAQ0));
    fZDCFlowVect[1].Set(QAReR,QAImR);
  }

  fillstep=4.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);

  if(fZDCEcomTotvsVtxHist[0]) {
    QCReR -= fZDCEcomTotvsVtxHist[0]->GetBinContent(fZDCEcomTotvsVtxHist[0]->FindBin(fCentralityEBE,fZNCQ0,fVtxPosCor[0]));
    QCImR -= fZDCEcomTotvsVtxHist[1]->GetBinContent(fZDCEcomTotvsVtxHist[1]->FindBin(fCentralityEBE,fZNCQ0,fVtxPosCor[0]));
    QCReR -= fZDCEcomTotvsVtxHist[4]->GetBinContent(fZDCEcomTotvsVtxHist[4]->FindBin(fCentralityEBE,fZNCQ0,fVtxPosCor[1]));
    QCImR -= fZDCEcomTotvsVtxHist[5]->GetBinContent(fZDCEcomTotvsVtxHist[5]->FindBin(fCentralityEBE,fZNCQ0,fVtxPosCor[1]));
    QCReR -= fZDCEcomTotvsVtxHist[8]->GetBinContent(fZDCEcomTotvsVtxHist[8]->FindBin(fCentralityEBE,fZNCQ0,fVtxPosCor[2]));
    QCImR -= fZDCEcomTotvsVtxHist[9]->GetBinContent(fZDCEcomTotvsVtxHist[9]->FindBin(fCentralityEBE,fZNCQ0,fVtxPosCor[2]));
    fZDCFlowVect[0].Set(QCReR,QCImR);
    QAReR -= fZDCEcomTotvsVtxHist[2]->GetBinContent(fZDCEcomTotvsVtxHist[2]->FindBin(fCentralityEBE,fZNAQ0,fVtxPosCor[0]));
    QAImR -= fZDCEcomTotvsVtxHist[3]->GetBinContent(fZDCEcomTotvsVtxHist[3]->FindBin(fCentralityEBE,fZNAQ0,fVtxPosCor[0]));
    QAReR -= fZDCEcomTotvsVtxHist[6]->GetBinContent(fZDCEcomTotvsVtxHist[6]->FindBin(fCentralityEBE,fZNAQ0,fVtxPosCor[1]));
    QAImR -= fZDCEcomTotvsVtxHist[7]->GetBinContent(fZDCEcomTotvsVtxHist[7]->FindBin(fCentralityEBE,fZNAQ0,fVtxPosCor[1]));
    QAReR -= fZDCEcomTotvsVtxHist[10]->GetBinContent(fZDCEcomTotvsVtxHist[10]->FindBin(fCentralityEBE,fZNAQ0,fVtxPosCor[2]));
    QAImR -= fZDCEcomTotvsVtxHist[11]->GetBinContent(fZDCEcomTotvsVtxHist[11]->FindBin(fCentralityEBE,fZNAQ0,fVtxPosCor[2]));
    fZDCFlowVect[1].Set(QAReR,QAImR);
  }

  fillstep=5.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);

  // store Q-vectors vs Energy in the common tower
  if(QMC>0. && QMA>0. && sqrt(QCRe*QCRe+QCIm*QCIm)>1.E-6 && sqrt(QARe*QARe+QAIm*QAIm)>1.E-6) {
//    fCRCZDCQVecECom[fRunBin][0]->Fill(fCentralityEBE,fZNCQ0,QCReR);
//    fCRCZDCQVecECom[fRunBin][1]->Fill(fCentralityEBE,fZNCQ0,QCImR);
//    fCRCZDCQVecECom[fRunBin][2]->Fill(fCentralityEBE,fZNAQ0,QAReR);
//    fCRCZDCQVecECom[fRunBin][3]->Fill(fCentralityEBE,fZNAQ0,QAImR);
//    fCRCZDCQVecEComTot[0]->Fill(fCentralityEBE,fZNCQ0,fVtxPosCor[0],QCReR);
//    fCRCZDCQVecEComTot[1]->Fill(fCentralityEBE,fZNCQ0,fVtxPosCor[0],QCImR);
//    fCRCZDCQVecEComTot[2]->Fill(fCentralityEBE,fZNAQ0,fVtxPosCor[0],QAReR);
//    fCRCZDCQVecEComTot[3]->Fill(fCentralityEBE,fZNAQ0,fVtxPosCor[0],QAImR);
//    fCRCZDCQVecEComTot[4]->Fill(fCentralityEBE,fZNCQ0,fVtxPosCor[1],QCReR);
//    fCRCZDCQVecEComTot[5]->Fill(fCentralityEBE,fZNCQ0,fVtxPosCor[1],QCImR);
//    fCRCZDCQVecEComTot[6]->Fill(fCentralityEBE,fZNAQ0,fVtxPosCor[1],QAReR);
//    fCRCZDCQVecEComTot[7]->Fill(fCentralityEBE,fZNAQ0,fVtxPosCor[1],QAImR);
//    fCRCZDCQVecEComTot[8]->Fill(fCentralityEBE,fZNCQ0,fVtxPosCor[2],QCReR);
//    fCRCZDCQVecEComTot[9]->Fill(fCentralityEBE,fZNCQ0,fVtxPosCor[2],QCImR);
//    fCRCZDCQVecEComTot[10]->Fill(fCentralityEBE,fZNAQ0,fVtxPosCor[2],QAReR);
//    fCRCZDCQVecEComTot[11]->Fill(fCentralityEBE,fZNAQ0,fVtxPosCor[2],QAImR);
  }

  // final recentering vs centrality (if needed)
  if (fZDCQHist[4]) {
    //      Double_t EvPlZDCC = TMath::ATan2(QCImR,QCReR);
    //      Double_t AvCos2C = fZDCFitSec[0]->Eval(fCentralityEBE);
    //      Double_t AvSin2C = fZDCFitSec[1]->Eval(fCentralityEBE);
    //      EvPlZDCC += -AvSin2C*cos(2.*EvPlZDCC) + AvCos2C*sin(2.*EvPlZDCC);
    //      if(EvPlZDCC<-TMath::Pi()) EvPlZDCC += TMath::TwoPi();
    //      if(EvPlZDCC> TMath::Pi()) EvPlZDCC -= TMath::TwoPi();
    //      fZDCFlowVect[0].SetMagPhi(sqrt(QCReR*QCReR+QCImR*QCImR),EvPlZDCC,QMC);
    //      QCReR = fZDCFlowVect[0].X();
    //      QCImR = fZDCFlowVect[0].Y();
    //
    //      Double_t EvPlZDCA = TMath::ATan2(QAImR,QAReR);
    //      Double_t AvCos2A = fZDCFitSec[2]->Eval(fCentralityEBE);
    //      Double_t AvSin2A = fZDCFitSec[3]->Eval(fCentralityEBE);
    //      EvPlZDCA += -AvSin2A*cos(2.*EvPlZDCA) + AvCos2A*sin(2.*EvPlZDCA);
    //      if(EvPlZDCA<-TMath::Pi()) EvPlZDCA += TMath::TwoPi();
    //      if(EvPlZDCA> TMath::Pi()) EvPlZDCA -= TMath::TwoPi();
    //      fZDCFlowVect[1].SetMagPhi(sqrt(QAReR*QAReR+QAImR*QAImR),EvPlZDCA,QMA);
    //      QAReR = fZDCFlowVect[1].X();
    //      QAImR = fZDCFlowVect[1].Y();

    Double_t AvQCRe = fZDCQHist[4]->GetBinContent(fZDCQHist[4]->FindBin(fCentralityEBE));
    Double_t AvQCIm = fZDCQHist[5]->GetBinContent(fZDCQHist[5]->FindBin(fCentralityEBE));

    Double_t AvQARe = fZDCQHist[6]->GetBinContent(fZDCQHist[6]->FindBin(fCentralityEBE));
    Double_t AvQAIm = fZDCQHist[7]->GetBinContent(fZDCQHist[7]->FindBin(fCentralityEBE));

    QCReR -= AvQCRe;
    QCImR -= AvQCIm;
    fZDCFlowVect[0].Set(QCReR,QCImR);

    QAReR -= AvQARe;
    QAImR -= AvQAIm;
    fZDCFlowVect[1].Set(QAReR,QAImR);
  }

  fillstep=6.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);

  // recenter vs vtx vs cen vs run number (through fits)

  if(fZDCVtxFitHist[0]) {
    for (Int_t i=0; i<3; i++) {
      QCReR -= fVtxPosCor[i]*fZDCVtxFitCenProjHist[0][i]->Interpolate(fCentralityEBE);
      QCImR -= fVtxPosCor[i]*fZDCVtxFitCenProjHist[1][i]->Interpolate(fCentralityEBE);
      QAReR -= fVtxPosCor[i]*fZDCVtxFitCenProjHist[2][i]->Interpolate(fCentralityEBE);
      QAImR -= fVtxPosCor[i]*fZDCVtxFitCenProjHist[3][i]->Interpolate(fCentralityEBE);
    }
    fZDCFlowVect[0].Set(QCReR,QCImR);
    fZDCFlowVect[1].Set(QAReR,QAImR);
  }

  fillstep=7.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);

  // second iteration (2D)

  if (fZDCVtxCenHistMagPol[fCenBin][0]) {
    if(fbFlagIsPosMagField) {
      QCReR -= fZDCVtxCenHistMagPol[fCenBin][0]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][0]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      QCImR -= fZDCVtxCenHistMagPol[fCenBin][1]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][1]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      QAReR -= fZDCVtxCenHistMagPol[fCenBin][4]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][4]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      QAImR -= fZDCVtxCenHistMagPol[fCenBin][5]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][5]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
    } else {
      QCReR -= fZDCVtxCenHistMagPol[fCenBin][2]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][2]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      QCImR -= fZDCVtxCenHistMagPol[fCenBin][3]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][3]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      QAReR -= fZDCVtxCenHistMagPol[fCenBin][6]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][6]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      QAImR -= fZDCVtxCenHistMagPol[fCenBin][7]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][7]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
    }
    fZDCFlowVect[0].Set(QCReR,QCImR);
    fZDCFlowVect[1].Set(QAReR,QAImR);
  }

  fillstep=8.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);

  if(fZDCVtxFitHist2[0]) {
    for (Int_t i=0; i<3; i++) {
      Double_t c1 = fZDCVtxFitCenProjHist2[0][i]->GetBinContent(1);
      Double_t c2 = fZDCVtxFitCenProjHist2[0][i]->GetBinContent(2);
      Double_t c3 = fZDCVtxFitCenProjHist2[0][i]->GetBinContent(3);
      QCReR -= fVtxPosCor[i]*c1 + fVtxPosCor[i]*fVtxPosCor[i]*c2 + fVtxPosCor[i]*fVtxPosCor[i]*fVtxPosCor[i]*c3;
      c1 = fZDCVtxFitCenProjHist2[1][i]->GetBinContent(1);
      c2 = fZDCVtxFitCenProjHist2[1][i]->GetBinContent(2);
      c3 = fZDCVtxFitCenProjHist2[1][i]->GetBinContent(3);
      QCImR -= fVtxPosCor[i]*c1 + fVtxPosCor[i]*fVtxPosCor[i]*c2 + fVtxPosCor[i]*fVtxPosCor[i]*fVtxPosCor[i]*c3;
      c1 = fZDCVtxFitCenProjHist2[2][i]->GetBinContent(1);
      c2 = fZDCVtxFitCenProjHist2[2][i]->GetBinContent(2);
      c3 = fZDCVtxFitCenProjHist2[2][i]->GetBinContent(3);
      QAReR -= fVtxPosCor[i]*c1 + fVtxPosCor[i]*fVtxPosCor[i]*c2 + fVtxPosCor[i]*fVtxPosCor[i]*fVtxPosCor[i]*c3;
      c1 = fZDCVtxFitCenProjHist2[3][i]->GetBinContent(1);
      c2 = fZDCVtxFitCenProjHist2[3][i]->GetBinContent(2);
      c3 = fZDCVtxFitCenProjHist2[3][i]->GetBinContent(3);
      QAImR -= fVtxPosCor[i]*c1 + fVtxPosCor[i]*fVtxPosCor[i]*c2 + fVtxPosCor[i]*fVtxPosCor[i]*fVtxPosCor[i]*c3;
    }
    fZDCFlowVect[0].Set(QCReR,QCImR);
    fZDCFlowVect[1].Set(QAReR,QAImR);
  }

  fillstep=9.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);

  Double_t recRefMul = fReferenceMultiplicityEBE-fMultCutAv->GetBinContent(fMultCutAv->FindBin(fCentralityEBE));
  if(fZDCBinsCenRefMultRbRProj[fCenBin][0]) {
    QCReR -= fZDCBinsCenRefMultRbRProj[fCenBin][0]->Interpolate(recRefMul);
    QCImR -= fZDCBinsCenRefMultRbRProj[fCenBin][1]->Interpolate(recRefMul);
    QAReR -= fZDCBinsCenRefMultRbRProj[fCenBin][2]->Interpolate(recRefMul);
    QAImR -= fZDCBinsCenRefMultRbRProj[fCenBin][3]->Interpolate(recRefMul);
    fZDCFlowVect[0].Set(QCReR,QCImR);
    fZDCFlowVect[1].Set(QAReR,QAImR);
  }

  fillstep=10.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);

//  Int_t EZDCCBin = fCRCZDCQVecDummyEZDCBins[fCenBin]->GetXaxis()->FindBin(QMCrec)-1;
//  Int_t EZDCABin = fCRCZDCQVecDummyEZDCBins[fCenBin]->GetXaxis()->FindBin(QMArec)-1;
//
//  if(fZDCQVecVtxCenEZDC3D[0][0][0]) {
//    Bool_t pass2=kTRUE;
//    // exclude events with vtx outside of range
//    if(fVtxPosCor[0] < fZDCQVecVtxCenEZDC3D[0][0][0]->GetXaxis()->GetXmin() || fVtxPosCor[0] > fZDCQVecVtxCenEZDC3D[0][0][0]->GetXaxis()->GetXmax()) pass2 = kFALSE;
//    if(fVtxPosCor[1] < fZDCQVecVtxCenEZDC3D[0][0][0]->GetYaxis()->GetXmin() || fVtxPosCor[1] > fZDCQVecVtxCenEZDC3D[0][0][0]->GetYaxis()->GetXmax()) pass2 = kFALSE;
//    if(fVtxPosCor[2] < fZDCQVecVtxCenEZDC3D[0][0][0]->GetZaxis()->GetXmin() || fVtxPosCor[2] > fZDCQVecVtxCenEZDC3D[0][0][0]->GetZaxis()->GetXmax()) pass2 = kFALSE;
//    // exclude events with EZDC outside of range
//    Int_t EZDCCBin = fCRCZDCQVecDummyEZDCBins[fCenBin]->GetXaxis()->FindBin(QMCrec)-1;
//    Int_t EZDCABin = fCRCZDCQVecDummyEZDCBins[fCenBin]->GetXaxis()->FindBin(QMArec)-1;
////    if(EZDCCBin<=0 || EZDCCBin>=10 || EZDCABin<=0 || EZDCABin>=10) pass2 = kFALSE;
//    if(EZDCCBin<0) EZDCCBin=0;
//    if(EZDCCBin>9) EZDCCBin=9;
//    if(EZDCABin<0) EZDCABin=0;
//    if(EZDCABin>9) EZDCABin=9;
//    if(pass2) {
//      // check if possible to interpolate
//      Bool_t bInterp = kTRUE;
//      Int_t bx = fZDCQVecVtxCenEZDC3D[0][0][0]->GetXaxis()->FindBin(fVtxPosCor[0]);
//      Int_t by = fZDCQVecVtxCenEZDC3D[0][0][0]->GetYaxis()->FindBin(fVtxPosCor[1]);
//      Int_t bz = fZDCQVecVtxCenEZDC3D[0][0][0]->GetZaxis()->FindBin(fVtxPosCor[2]);
//      if(bx==1 || bx==fZDCQVecVtxCenEZDC3D[0][0][0]->GetXaxis()->GetNbins()) bInterp = kFALSE;
//      if(by==1 || by==fZDCQVecVtxCenEZDC3D[0][0][0]->GetYaxis()->GetNbins()) bInterp = kFALSE;
//      if(bz==1 || bz==fZDCQVecVtxCenEZDC3D[0][0][0]->GetZaxis()->GetNbins()) bInterp = kFALSE;
//      if(bInterp) {
//        QCReR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][0]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
//        QCImR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][1]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
//        QAReR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][2]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
//        QAImR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][3]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
//      } else {
//        QCReR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][0]->GetBinContent(bx,by,bz);
//        QCImR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][1]->GetBinContent(bx,by,bz);
//        QAReR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][2]->GetBinContent(bx,by,bz);
//        QAImR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][3]->GetBinContent(bx,by,bz);
//      }
//      fZDCFlowVect[0].Set(QCReR,QCImR);
//      fZDCFlowVect[1].Set(QAReR,QAImR);
//    } else {
//      QCReR = 0.; QCImR = 0.; QAReR = 0.; QAImR = 0.; QMC=0.; QMA=0.;
//      fZDCFlowVect[0].Set(QCReR,QCImR);
//      fZDCFlowVect[0].SetMult(0.);
//      fZDCFlowVect[1].Set(QAReR,QAImR);
//      fZDCFlowVect[1].SetMult(0.);
//    }
//  }

  fillstep=11.5;
  fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
  fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
  fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
  fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);

  // if correctly set (via fMinMulZN), cut on ZDC Q-vector distribution
  // WARNING: QARe inverted !!!

  pass = PassCutZDCQVecDis(QCReR,QCImR,-QAReR,QAImR);

  // additional recentering vs centrality (if needed)
  if (fZDCQHist[8]) {
    QCReR -= fZDCQHist[8]->GetBinContent(fZDCQHist[8]->FindBin(fCentralityEBE));
    QCImR -= fZDCQHist[9]->GetBinContent(fZDCQHist[9]->FindBin(fCentralityEBE));
    fZDCFlowVect[0].Set(QCReR,QCImR);

    QAReR -= fZDCQHist[10]->GetBinContent(fZDCQHist[10]->FindBin(fCentralityEBE));
    QAImR -= fZDCQHist[11]->GetBinContent(fZDCQHist[11]->FindBin(fCentralityEBE));
    fZDCFlowVect[1].Set(QAReR,QAImR);
  }

  if(pass) {
    fillstep=12.5;
    fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
    fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
    fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
    fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);
  }

  // ***************************************************************************
  // store results after correction
  // ***************************************************************************

  if(QMC>0. && QMA>0. && sqrt(QCReR*QCReR+QCImR*QCImR)>1.E-6 && sqrt(QAReR*QAReR+QAImR*QAImR)>1.E-6 && pass) {

    fillstep=13.5;
    fCRCZDCQVecCorSteps[0]->Fill(fillstep,fCentralityEBE,-QCReR*QAReR);
    fCRCZDCQVecCorSteps[1]->Fill(fillstep,fCentralityEBE,QCImR*QAImR);
    fCRCZDCQVecCorSteps[2]->Fill(fillstep,fCentralityEBE,QCReR*QAImR);
    fCRCZDCQVecCorSteps[3]->Fill(fillstep,fCentralityEBE,-QCImR*QAReR);

    fCRCZDCQVecCCorr[fRunBin][0]->Fill(fCentralityEBE,QCReR);
    fCRCZDCQVecCCorr[fRunBin][1]->Fill(fCentralityEBE,QCImR);

//    fCRCZDCQVecCov[fRunBin][0]->Fill(fCentralityEBE,fVtxPosCor[0],QCReR);
//    fCRCZDCQVecCov[fRunBin][1]->Fill(fCentralityEBE,fVtxPosCor[1],QCReR);
//    fCRCZDCQVecCov[fRunBin][2]->Fill(fCentralityEBE,fVtxPosCor[2],QCReR);
//    fCRCZDCQVecCov[fRunBin][3]->Fill(fCentralityEBE,fVtxPosCor[0],QCImR);
//    fCRCZDCQVecCov[fRunBin][4]->Fill(fCentralityEBE,fVtxPosCor[1],QCImR);
//    fCRCZDCQVecCov[fRunBin][5]->Fill(fCentralityEBE,fVtxPosCor[2],QCImR);

//    if(fbFlagIsPosMagField) {
//      fCRCZDCQVecVtxPosCen[fCenBin][0]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QCReR);
//      fCRCZDCQVecVtxPosCen[fCenBin][1]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QCImR);
//    } else {
//      fCRCZDCQVecVtxPosCen[fCenBin][2]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QCReR);
//      fCRCZDCQVecVtxPosCen[fCenBin][3]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QCImR);
//    }

    fCRCZDCQVecACorr[fRunBin][0]->Fill(fCentralityEBE,QAReR);
    fCRCZDCQVecACorr[fRunBin][1]->Fill(fCentralityEBE,QAImR);
//    fCRCZDCQVecCov[fRunBin][6]->Fill(fCentralityEBE,fVtxPosCor[0],QAReR);
//    fCRCZDCQVecCov[fRunBin][7]->Fill(fCentralityEBE,fVtxPosCor[1],QAReR);
//    fCRCZDCQVecCov[fRunBin][8]->Fill(fCentralityEBE,fVtxPosCor[2],QAReR);
//    fCRCZDCQVecCov[fRunBin][9]->Fill(fCentralityEBE,fVtxPosCor[0],QAImR);
//    fCRCZDCQVecCov[fRunBin][10]->Fill(fCentralityEBE,fVtxPosCor[1],QAImR);
//    fCRCZDCQVecCov[fRunBin][11]->Fill(fCentralityEBE,fVtxPosCor[2],QAImR);

//    if(fbFlagIsPosMagField) {
//      fCRCZDCQVecVtxPosCen[fCenBin][4]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QAReR);
//      fCRCZDCQVecVtxPosCen[fCenBin][5]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QAImR);
//    } else {
//      fCRCZDCQVecVtxPosCen[fCenBin][6]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QAReR);
//      fCRCZDCQVecVtxPosCen[fCenBin][7]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QAImR);
//    }

    fhAvRefMulRbR->Fill(fRunBin+0.5,fCentralityEBE,fReferenceMultiplicityEBE);
    fhAvQMCRbR->Fill(fRunBin+0.5,fCentralityEBE,QMC);
    fhAvQMARbR->Fill(fRunBin+0.5,fCentralityEBE,QMA);

    //      for(Int_t i=0;i<3;i++) {
//        if(fbFlagIsPosMagField) {
//          fCRCZDCQVecVtxCenEZDC[i][0]->Fill(fCentralityEBE,QMCrec,fVtxPosCor[i],QCReR);
//          fCRCZDCQVecVtxCenEZDC[i][1]->Fill(fCentralityEBE,QMCrec,fVtxPosCor[i],QCImR);
//          fCRCZDCQVecVtxCenEZDC[i][4]->Fill(fCentralityEBE,QMArec,fVtxPosCor[i],QAReR);
//          fCRCZDCQVecVtxCenEZDC[i][5]->Fill(fCentralityEBE,QMArec,fVtxPosCor[i],QAImR);
//        } else {
//          fCRCZDCQVecVtxCenEZDC[i][2]->Fill(fCentralityEBE,QMCrec,fVtxPosCor[i],QCReR);
//          fCRCZDCQVecVtxCenEZDC[i][3]->Fill(fCentralityEBE,QMCrec,fVtxPosCor[i],QCImR);
//          fCRCZDCQVecVtxCenEZDC[i][6]->Fill(fCentralityEBE,QMArec,fVtxPosCor[i],QAReR);
//          fCRCZDCQVecVtxCenEZDC[i][7]->Fill(fCentralityEBE,QMArec,fVtxPosCor[i],QAImR);
//        }
//      }
//      Int_t EZDCCBin = fCRCZDCQVecDummyEZDCBins[fCenBin]->GetXaxis()->FindBin(QMCrec)-1;
//      Int_t EZDCABin = fCRCZDCQVecDummyEZDCBins[fCenBin]->GetXaxis()->FindBin(QMArec)-1;
//      if(fbFlagIsPosMagField) {
//        fCRCZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][0]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QCReR);
//        fCRCZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][1]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QCImR);
//        fCRCZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][4]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QAReR);
//        fCRCZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][5]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QAImR);
//      } else {
//        fCRCZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][2]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QCReR);
//        fCRCZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][3]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QCImR);
//        fCRCZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][6]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QAReR);
    //        fCRCZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][7]->Fill(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2],QAImR);
    //      }
    //    }

    if( fInvertZDC ) QAReR = -QAReR;

    //    fhZNCenDis[0]->Fill(fCentralityEBE,QCReR,QCImR);
    //    fhZNCenDis[1]->Fill(fCentralityEBE,QAReR,QAImR);

    // Double_t EvPlZDCfull = TMath::ATan2(QAImR-QCImR,QAReR-QCReR);
    //    fCRCZDCQVecEP[fRunBin][2]->Fill(fCentralityEBE,EvPlZDCfull);
    // EvPlZDCfull = TMath::ATan2(QAImR+QCImR,QAReR+QCReR);
    //    fCRCZDCQVecEP[fRunBin][3]->Fill(fCentralityEBE,EvPlZDCfull);
    fCRCZDCQVecRes[fRunBin][0]->Fill(fCentralityEBE,QCReR*QAReR);
    fCRCZDCQVecRes[fRunBin][1]->Fill(fCentralityEBE,QCImR*QAImR);
    fCRCZDCQVecRes[fRunBin][2]->Fill(fCentralityEBE,QCReR*QAImR);
    fCRCZDCQVecRes[fRunBin][3]->Fill(fCentralityEBE,QCImR*QAReR);

    if(fStoreZDCQVecVtxPos) {
      fCRCZDCQVecDis[bw][fCenBin][0]->Fill(QCReR*QAImR,QCImR*QAReR);
      fCRCZDCQVecDis[bw][fCenBin][1]->Fill(QCReR*QAReR,QCImR*QAImR);
    }

    Double_t EvPlZNC = TMath::ATan2(QCImR,QCReR);
    if(EvPlZNC<0.) EvPlZNC += TMath::TwoPi();
    Double_t EvPlZNA = TMath::ATan2(QAImR,QAReR);
    if(EvPlZNA<0.) EvPlZNA += TMath::TwoPi();
    fZDCEPHist[0]->Fill(fCentralityEBE,EvPlZNC);
    fZDCEPHist[1]->Fill(fCentralityEBE,EvPlZNA);

  } else {
    fQAZDCCutsFlag = kFALSE;
  }

} // end of AliFlowAnalysisQvec::RecenterCRCQVecZDC()

//=======================================================================================================================

Bool_t AliFlowAnalysisQvec::PassCutZDCQVecDis(Double_t ZCRe, Double_t ZCIm, Double_t ZARe, Double_t ZAIm)
{
  Int_t bw = (fbFlagIsPosMagField==kTRUE?0:1);
  Bool_t passcut = kTRUE;

  if(fMinMulZN==13) {
    if(fCRCZDC2DCutZDCC[bw][fCenBin]) {
      if(fRandom->Uniform(0.,1.)>fCRCZDC2DCutZDCC[bw][fCenBin]->GetBinContent(fCRCZDC2DCutZDCC[bw][fCenBin]->GetXaxis()->FindBin(ZCRe),fCRCZDC2DCutZDCC[bw][fCenBin]->GetYaxis()->FindBin(ZCIm))) passcut = kFALSE;
    }
  }

  if(fMinMulZN==14) {
    if(fCRCZDC2DCutZDCA[bw][fCenBin]) {
      if(fRandom->Uniform(0.,1.)>fCRCZDC2DCutZDCA[bw][fCenBin]->GetBinContent(fCRCZDC2DCutZDCA[bw][fCenBin]->GetXaxis()->FindBin(ZARe),fCRCZDC2DCutZDCA[bw][fCenBin]->GetYaxis()->FindBin(ZAIm))) passcut = kFALSE;
    }
  }

  if(fMinMulZN==15) {
    if(fCRCZDC2DCutZDCC[bw][fCenBin] && fCRCZDC2DCutZDCA[bw][fCenBin]) {
      if(fRandom->Uniform(0.,1.)>fCRCZDC2DCutZDCC[bw][fCenBin]->GetBinContent(fCRCZDC2DCutZDCC[bw][fCenBin]->GetXaxis()->FindBin(ZCRe),fCRCZDC2DCutZDCC[bw][fCenBin]->GetYaxis()->FindBin(ZCIm))) passcut = kFALSE;
      if(fRandom->Uniform(0.,1.)>fCRCZDC2DCutZDCA[bw][fCenBin]->GetBinContent(fCRCZDC2DCutZDCA[bw][fCenBin]->GetXaxis()->FindBin(ZARe),fCRCZDC2DCutZDCA[bw][fCenBin]->GetYaxis()->FindBin(ZAIm))) passcut = kFALSE;
    }
  }

  if(fMinMulZN==16) {
    if(fCRCZDC2DCutZDCC[bw][fCenBin]) {
      if(fRandom->Uniform(0.,1.)>fCRCZDC2DCutZDCC[bw][fCenBin]->GetBinContent(fCRCZDC2DCutZDCC[bw][fCenBin]->GetXaxis()->FindBin(ZCRe*ZAIm),fCRCZDC2DCutZDCC[bw][fCenBin]->GetYaxis()->FindBin(ZCIm*ZARe))) passcut = kFALSE;
      if(ZCRe*ZAIm < fCRCZDC2DCutZDCC[bw][fCenBin]->GetXaxis()->GetXmin() || ZCRe*ZAIm > fCRCZDC2DCutZDCC[bw][fCenBin]->GetXaxis()->GetXmax()) passcut = kFALSE;
      if(ZCIm*ZARe < fCRCZDC2DCutZDCC[bw][fCenBin]->GetYaxis()->GetXmin() || ZCIm*ZARe > fCRCZDC2DCutZDCC[bw][fCenBin]->GetYaxis()->GetXmax()) passcut = kFALSE;
    }
  }

  return passcut;
}

//=======================================================================================================================

void AliFlowAnalysisQvec::InitializeArraysForCME()
{
  for(Int_t c=0;c<4;c++) {
    for (Int_t h=0;h<fCRCnHar;h++) {
      fCMEQRe[c][h] = NULL;
      fCMEQIm[c][h] = NULL;
      fCMEMult[c][h] = NULL;
    }
  }
  
  //@Shi Define CME Qvector for spectator plane participant plane method
  for(Int_t c=0;c<2;c++) { // index c represents the power of weight weight^c*cos(phi)
    for (Int_t h=0;h<fCRCnHar;h++) { // index h represents harmonics cos((h+1)*phi)
      fCMEQReBothCharge[c][h] = NULL;
      fCMEQImBothCharge[c][h] = NULL;
      fCMEMultBothCharge[c][h] = NULL;
    }
  }
  
}

//=======================================================================================================================

void AliFlowAnalysisQvec::CheckPointersUsedInFinish()
{
  // Check all pointers used in method Finish().

} // end of void AliFlowAnalysisQvec::CheckPointersUsedInFinish()

//=======================================================================================================================

void AliFlowAnalysisQvec::CheckPointersUsedInMake()
{
  // Check all pointers used in method Make(). // to be improved - check other pointers as well

} // end of void AliFlowAnalysisQvec::CheckPointersUsedInMake()

//=======================================================================================================================

Int_t AliFlowAnalysisQvec::GetCRCCenBin(Double_t Centrality)
{
  Int_t CenBin=-1;
  if (Centrality>0. && Centrality<5.) CenBin=0;
  if (Centrality>5. && Centrality<10.) CenBin=1;
  if (Centrality>10. && Centrality<20.) CenBin=2;
  if (Centrality>20. && Centrality<30.) CenBin=3;
  if (Centrality>30. && Centrality<40.) CenBin=4;
  if (Centrality>40. && Centrality<50.) CenBin=5;
  if (Centrality>50. && Centrality<60.) CenBin=6;
  if (Centrality>60. && Centrality<70.) CenBin=7;
  if (Centrality>70. && Centrality<80.) CenBin=8;
  if (Centrality>80. && Centrality<90.) CenBin=9;
  if (Centrality>90. && Centrality<100.) CenBin=10;
  if (CenBin>=fCRCnCen) CenBin=-1;
  if (fCRCnCen==1) CenBin=0;
  return CenBin;
} // end of AliFlowAnalysisCRC::GetCRCCenBin(Double_t Centrality)
