/**************************************************************************
 * Contributors are not mentioned at all.                                 *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright noticxse appears in all *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//-----------------------------------------------------------------
//                 AliAnalysisTaskAllPtcv2 class
//-----------------------------------------------------------------

class TTree;
class TParticle;
class TVector3;

#include "AliAnalysisManager.h"
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0; 

#include <iostream>

#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TChain.h"
#include "Riostream.h"
#include "AliLog.h"
#include "AliCascadeVertexer.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliAODEvent.h"
#include "AliInputEventHandler.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include <TRandom3.h>
#include "AliMultSelection.h"
#include "AliQnCorrectionsManager.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliEventCuts.h"

#include "AliAnalysisTaskAllPtcv2.h"

ClassImp(AliAnalysisTaskAllPtcv2)

using std::cout;
using std::endl;
    
//________________________________________________________________________
AliAnalysisTaskAllPtcv2::AliAnalysisTaskAllPtcv2() 
: AliAnalysisTaskSE(),
  fESDevent(0),                         //! 
  fAODevent(0),                           //! 
  fevent(0),   
  fAnalysisType("AOD"),
  fisPrimCut(kFALSE),
  fptc(1),     
  fVzmax(10),
  fCentrality("V0M"),
  fApplyFlatten(kFALSE),
  fYear(2015),
  fHarmonic(2),
  fListHist(0),              
  fHistEventMultiplicity(0), 
  fHistTrackMultiplicity(0),
  fHistTrackMultiplicityCentral(0),    
  fHistTrackMultiplicitySemiCentral(0),
  fHistTrackMultiplicityMB(0),
  fHistTrackMultiplicityINT7(0),
  fhBB(0),
  fhTOF(0),
  fhMassTOF(0),
  EPVzAvsCentrality(0), 
  EPVzCvsCentrality(0), 
  EPTPCvsCentrality(0), 
  EPVzvsCentrality(0), 
  EPTPCpvsCentrality(0), 
  EPTPCnvsCentrality(0),  
  hEvPlaneTPCvsEvPVz0010(0),                      
  hEvPlaneTPCvsEvPVz1020(0), 
  hEvPlaneTPCvsEvPVz2030(0),
  hEvPlaneTPCvsEvPVz3040(0),                      
  hEvPlaneTPCvsEvPVz4050(0),   
  hEvPlaneTPCvsEvPVz5060(0),                      
  hEvPlaneTPCvsEvPVz6080(0),      
  hCos2DeltaTPCVzAvsCentrality(0),
  hCos2DeltaTPCVzCvsCentrality(0),
  hCos2DeltaVzAVzCvsCentrality(0),
  hCos2DeltaVzMVzAvsCentrality(0),
  hCos2DeltaVzMVzCvsCentrality(0),
  hCos2DeltaVzATPCvsCentrality(0),
  hCos2DeltaVzCTPCvsCentrality(0),
  hCos2DeltaVzCVzAvsCentrality(0),
  hCos2DeltaVzMTPCpvsCentrality(0),
  hCos2DeltaVzMTPCnvsCentrality(0),
  hCos2DeltaTPCpTPCnvsCentrality(0),
  hQVzAQVzCvsCentrality(0),
  hQxVzAvsCentrality(0),
  hQyVzAvsCentrality(0),
  hQxVzCvsCentrality(0),
  hQyVzCvsCentrality(0),
  hQxVzMvsCentrality(0),
  hQyVzMvsCentrality(0),
  hCosdeltaphiTPCvsPtAll0010(0),      
  hCosdeltaphiV0AvsPtAll0010(0),      
  hCosdeltaphiV0CvsPtAll0010(0),   
  huqV0AvsPtAll0010(0),
  huqV0CvsPtAll0010(0),
  hCosdeltaphiTPCvsPtAll1020(0),      
  hCosdeltaphiV0AvsPtAll1020(0),      
  hCosdeltaphiV0CvsPtAll1020(0),       
  huqV0AvsPtAll1020(0),
  huqV0CvsPtAll1020(0),
  hCosdeltaphiTPCvsPtAll2030(0),      
  hCosdeltaphiV0AvsPtAll2030(0),      
  hCosdeltaphiV0CvsPtAll2030(0),      
  huqV0AvsPtAll2030(0),
  huqV0CvsPtAll2030(0),
  hCosdeltaphiTPCvsPtAll3040(0),      
  hCosdeltaphiV0AvsPtAll3040(0),      
  hCosdeltaphiV0CvsPtAll3040(0),       
  huqV0AvsPtAll3040(0),
  huqV0CvsPtAll3040(0),
  hCosdeltaphiTPCvsPtAll4050(0),      
  hCosdeltaphiV0AvsPtAll4050(0),      
  hCosdeltaphiV0CvsPtAll4050(0),      
  huqV0AvsPtAll4050(0),
  huqV0CvsPtAll4050(0),
  hCosdeltaphiTPCvsPtAll5060(0),      
  hCosdeltaphiV0AvsPtAll5060(0),      
  hCosdeltaphiV0CvsPtAll5060(0),      
  huqV0AvsPtAll5060(0),
  huqV0CvsPtAll5060(0),
  hCosdeltaphiTPCvsPtAll6080(0),      
  hCosdeltaphiV0AvsPtAll6080(0),      
  hCosdeltaphiV0CvsPtAll6080(0),      
  huqV0AvsPtAll6080(0),
  huqV0CvsPtAll6080(0),
  hphivsPtAll0010(0),
  hphivsPtAll1020(0),
  hphivsPtAll2030(0),
  hphivsPtAll3040(0),
  hphivsPtAll4050(0),
  hphivsPtAll5060(0),
  hphivsPtAll6080(0),
  huqV0AvsuqV0C0010(0),
  huqV0AvsuqV0C1020(0),
  huqV0AvsuqV0C2030(0),
  huqV0AvsuqV0C3040(0),
  huqV0AvsuqV0C4050(0),
  huqV0AvsuqV0C5060(0),
  huqV0AvsuqV0C6080(0),
  huqV0AxuqV0CvsPtAll0010(0),
  huqV0AxuqV0CvsPtAll1020(0),
  huqV0AxuqV0CvsPtAll2030(0),
  huqV0AxuqV0CvsPtAll3040(0),
  huqV0AxuqV0CvsPtAll4050(0),
  huqV0AxuqV0CvsPtAll5060(0),
  huqV0AxuqV0CvsPtAll6080(0),
  fESDtrackCuts(0),
  fESDtrackCutsEP(0),
  fPIDResponse(0),
  fEventCuts(0)
{
  // Dummy Constructor 
  fESDtrackCuts   = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  fESDtrackCutsEP = new AliESDtrackCuts("AliESDtrackCutsEP","AliESDtrackCutsEP");
  Initialize();
  cout<<"Dummy constructor"<<endl;
}

//________________________________________________________________________
AliAnalysisTaskAllPtcv2::AliAnalysisTaskAllPtcv2(const char *name) 
: AliAnalysisTaskSE(name),
  fESDevent(0),                         //! 
  fAODevent(0),                           //! 
  fevent(0),   
  fAnalysisType("AOD"),
  fisPrimCut(kFALSE),
  fptc(1),     
  fVzmax(10),
  fCentrality("V0M"),
  fApplyFlatten(kFALSE),
  fYear(2015),
  fHarmonic(2),
  fListHist(0),              
  fHistEventMultiplicity(0), 
  fHistTrackMultiplicity(0),
  fHistTrackMultiplicityCentral(0),    
  fHistTrackMultiplicitySemiCentral(0),
  fHistTrackMultiplicityMB(0),
  fHistTrackMultiplicityINT7(0),
  fhBB(0),
  fhTOF(0),
  fhMassTOF(0),
  EPVzAvsCentrality(0), 
  EPVzCvsCentrality(0), 
  EPTPCvsCentrality(0), 
  EPVzvsCentrality(0), 
  EPTPCpvsCentrality(0), 
  EPTPCnvsCentrality(0),  
  hEvPlaneTPCvsEvPVz0010(0),                      
  hEvPlaneTPCvsEvPVz1020(0), 
  hEvPlaneTPCvsEvPVz2030(0),
  hEvPlaneTPCvsEvPVz3040(0),                      
  hEvPlaneTPCvsEvPVz4050(0),   
  hEvPlaneTPCvsEvPVz5060(0),                      
  hEvPlaneTPCvsEvPVz6080(0),      
  hCos2DeltaTPCVzAvsCentrality(0),
  hCos2DeltaTPCVzCvsCentrality(0),
  hCos2DeltaVzAVzCvsCentrality(0),
  hCos2DeltaVzMVzAvsCentrality(0),
  hCos2DeltaVzMVzCvsCentrality(0),
  hCos2DeltaVzATPCvsCentrality(0),
  hCos2DeltaVzCTPCvsCentrality(0),
  hCos2DeltaVzCVzAvsCentrality(0),
  hCos2DeltaVzMTPCpvsCentrality(0),
  hCos2DeltaVzMTPCnvsCentrality(0),
  hCos2DeltaTPCpTPCnvsCentrality(0),
  hQVzAQVzCvsCentrality(0),
  hQxVzAvsCentrality(0),
  hQyVzAvsCentrality(0),
  hQxVzCvsCentrality(0),
  hQyVzCvsCentrality(0),
  hQxVzMvsCentrality(0),
  hQyVzMvsCentrality(0),
  hCosdeltaphiTPCvsPtAll0010(0),      
  hCosdeltaphiV0AvsPtAll0010(0),      
  hCosdeltaphiV0CvsPtAll0010(0),   
  huqV0AvsPtAll0010(0),
  huqV0CvsPtAll0010(0),
  hCosdeltaphiTPCvsPtAll1020(0),      
  hCosdeltaphiV0AvsPtAll1020(0),      
  hCosdeltaphiV0CvsPtAll1020(0),       
  huqV0AvsPtAll1020(0),
  huqV0CvsPtAll1020(0),
  hCosdeltaphiTPCvsPtAll2030(0),      
  hCosdeltaphiV0AvsPtAll2030(0),      
  hCosdeltaphiV0CvsPtAll2030(0),      
  huqV0AvsPtAll2030(0),
  huqV0CvsPtAll2030(0),
  hCosdeltaphiTPCvsPtAll3040(0),      
  hCosdeltaphiV0AvsPtAll3040(0),      
  hCosdeltaphiV0CvsPtAll3040(0),       
  huqV0AvsPtAll3040(0),
  huqV0CvsPtAll3040(0),
  hCosdeltaphiTPCvsPtAll4050(0),      
  hCosdeltaphiV0AvsPtAll4050(0),      
  hCosdeltaphiV0CvsPtAll4050(0),      
  huqV0AvsPtAll4050(0),
  huqV0CvsPtAll4050(0),
  hCosdeltaphiTPCvsPtAll5060(0),      
  hCosdeltaphiV0AvsPtAll5060(0),      
  hCosdeltaphiV0CvsPtAll5060(0),      
  huqV0AvsPtAll5060(0),
  huqV0CvsPtAll5060(0),
  hCosdeltaphiTPCvsPtAll6080(0),      
  hCosdeltaphiV0AvsPtAll6080(0),      
  hCosdeltaphiV0CvsPtAll6080(0),      
  huqV0AvsPtAll6080(0),
  huqV0CvsPtAll6080(0),
  hphivsPtAll0010(0),
  hphivsPtAll1020(0),
  hphivsPtAll2030(0),
  hphivsPtAll3040(0),
  hphivsPtAll4050(0),
  hphivsPtAll5060(0),
  hphivsPtAll6080(0),
  huqV0AvsuqV0C0010(0),
  huqV0AvsuqV0C1020(0),
  huqV0AvsuqV0C2030(0),
  huqV0AvsuqV0C3040(0),
  huqV0AvsuqV0C4050(0),
  huqV0AvsuqV0C5060(0),
  huqV0AvsuqV0C6080(0),
  huqV0AxuqV0CvsPtAll0010(0),
  huqV0AxuqV0CvsPtAll1020(0),
  huqV0AxuqV0CvsPtAll2030(0),
  huqV0AxuqV0CvsPtAll3040(0),
  huqV0AxuqV0CvsPtAll4050(0),
  huqV0AxuqV0CvsPtAll5060(0),
  huqV0AxuqV0CvsPtAll6080(0),
  fESDtrackCuts(0),
  fESDtrackCutsEP(0),
  fPIDResponse(0),
  fEventCuts(0)
{
  // Define input and output slots here
  // Input slot #0 works with a TChain
  //DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container ()

  //
  // create track cuts
  //
  fESDtrackCuts   = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  fESDtrackCutsEP = new AliESDtrackCuts("AliESDtrackCutsEP","AliESDtrackCutsEP");
  //
  cout<<"Real constructor"<<endl;
  Initialize();

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  
}

void AliAnalysisTaskAllPtcv2::Initialize()
{
  //
  // updating parameters in case of changes
  //
  // This were used for the analysis note
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(fisPrimCut,kTRUE);
  fESDtrackCuts->SetMaxDCAToVertexXY(3);
  fESDtrackCuts->SetMaxDCAToVertexZ(2);
  fESDtrackCuts->SetEtaRange(-0.8,0.8);
  
  //Test 1 : Check with only TPC tracks also for POIs
  //  fESDtrackCuts =  AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); 
  
  fESDtrackCutsEP = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); 

}

//________________________________________________________________________
Float_t AliAnalysisTaskAllPtcv2::GetPhi0Pi(Float_t phi){
  // Sets the phi angle in the range 0-pi
  Float_t result=phi;
  while(result<0){
    result=result+TMath::Pi();
  }
  while(result>TMath::Pi()){
    result=result-TMath::Pi();
  }
   return result;
}

//________________________________________________________________________
const AliQnCorrectionsQnVector *AliAnalysisTaskAllPtcv2::GetQnVectorFromList(const TList *list,
									     const char *subdetector,
									     const char *expectedstep,
									     const char *altstep)
{    
  AliQnCorrectionsQnVector *theQnVector = NULL;
  
  TList *pQvecList = dynamic_cast<TList*> (list->FindObject(subdetector));
  if (pQvecList != NULL) {
    /* the detector is present */
    if (TString(expectedstep).EqualTo("latest"))
      theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
    else
      theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(expectedstep);
    
    if (theQnVector == NULL || !(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0)) {
      /* the Qn vector for the expected step was not there */
      if (TString(altstep).EqualTo("latest"))
	theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
      else
	theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(altstep);
    }
  }
  if (theQnVector != NULL) {
    /* check the Qn vector quality */
    if (!(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0))
      /* not good quality, discarded */
      theQnVector = NULL;
  }
  return theQnVector;
}

//==================DEFINITION OF OUTPUT OBJECTS==============================

void AliAnalysisTaskAllPtcv2::UserCreateOutputObjects()
{
  // fUtils = new AliAnalysisUtils();
  // fUtils->SetCutOnZVertexSPD(kFALSE);

  //-------------------------------------------------------
  fListHist = new TList();
  fListHist->SetOwner();  // IMPORTANT!
  
  if(! fHistEventMultiplicity ){

    fHistEventMultiplicity   = new TH1F( "fHistEventMultiplicity" , "Nb of Events" , 14 ,0.5,14.5);

    fHistEventMultiplicity->GetXaxis()->SetBinLabel(1,"All Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(2,"Events w/PV & wo pileup");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(3,"Events w/|Vz|<10cm");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(4,"Central Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(5,"Semi-Central Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(6,"MB Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(7,"INT7 Events");
    //from HF
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(8,"nEvSelected");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(9,"nCandidatesSelected");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(10,"out of pt bounds");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(11,"mismatch lab");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(12,"non valid TPC EP");
    fListHist->Add(fHistEventMultiplicity);
  }
  
  if(! fHistTrackMultiplicity ){
    fHistTrackMultiplicity  = new TH2F( "fHistTrackMultiplicity", "Nb of Tracks MB Events |Vz| < 10", 250,0, 25000,105,0,105);
    fHistTrackMultiplicity->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicity->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicity);
  } 

  if(! fHistTrackMultiplicityCentral ){
    fHistTrackMultiplicityCentral  = new TH2F( "fHistTrackMultiplicityCentral", "Nb of Tracks MB Events |Vz| < 10", 250,0, 25000,105,0,105);
    fHistTrackMultiplicityCentral->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityCentral->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicityCentral);
  } 
  if(! fHistTrackMultiplicitySemiCentral ){
    fHistTrackMultiplicitySemiCentral  = new TH2F( "fHistTrackMultiplicitySemiCentral", "Nb of Tracks MB Events |Vz| < 10", 250,0, 25000,105,0,105);
    fHistTrackMultiplicitySemiCentral->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicitySemiCentral->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicitySemiCentral);
  } 
  if(! fHistTrackMultiplicityMB ){
    fHistTrackMultiplicityMB  = new TH2F( "fHistTrackMultiplicityMB", "Nb of Tracks MB Events |Vz| < 10", 250,0, 25000,105,0,105);
    fHistTrackMultiplicityMB->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityMB->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicityMB);
  } 

  if(! fHistTrackMultiplicityINT7 ){
    fHistTrackMultiplicityINT7  = new TH2F( "fHistTrackMultiplicityINT7", "Nb of Tracks INT7 Events |Vz| < 10", 250,0, 25000,105,0,105);
    fHistTrackMultiplicityINT7->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityINT7->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicityINT7);
  } 

  if(! fhBB ){
    fhBB = new TH2F( "fhBB" , "BetheBlochTPC" , 240,-10,10,250,0,1000);
    fListHist->Add(fhBB);
  }
  
  if(! fhTOF ){
    fhTOF = new TH2F( "fhTOF" , "Scatter Plot TOF" , 240,-10,10,500,0,1.2);
    fListHist->Add(fhTOF);
  }
  if(! fhMassTOF){
    fhMassTOF=new TH1F ("fhMassTOF","Particle Mass - TOF", 100,0 ,10);
    fListHist->Add(fhMassTOF);
  }
  
  EPVzAvsCentrality  = new TH2D("EPVzAvsCentrality" , "EPVzAvsCentrality" , 80, -0.5,2*TMath::Pi()+0.5, 105,0,105);
  EPVzCvsCentrality  = new TH2D("EPVzCvsCentrality" , "EPVzCvsCentrality" , 80, -0.5,2*TMath::Pi()+0.5, 105,0,105);
  EPTPCvsCentrality  = new TH2D("EPTPCvsCentrality" , "EPTPCvsCentrality" , 80, -0.5,2*TMath::Pi()+0.5, 105,0,105);
  EPVzvsCentrality   = new TH2D("EPVzvsCentrality"  , "EPVzvsCentrality"  , 80, -0.5,2*TMath::Pi()+0.5, 105,0,105);
  EPTPCpvsCentrality = new TH2D("EPTPCpvsCentrality", "EPTPCpvsCentrality", 80, -0.5,2*TMath::Pi()+0.5, 105,0,105);
  EPTPCnvsCentrality = new TH2D("EPTPCnvsCentrality", "EPTPCnvsCentrality", 80, -0.5,2*TMath::Pi()+0.5, 105,0,105);

  fListHist->Add(EPVzAvsCentrality);
  fListHist->Add(EPVzCvsCentrality);
  fListHist->Add(EPTPCvsCentrality);
  fListHist->Add(EPVzvsCentrality);
  fListHist->Add(EPTPCpvsCentrality);
  fListHist->Add(EPTPCnvsCentrality);
  
  
  hEvPlaneTPCvsEvPVz0010 = new TH2F("hEvPlaneTPCvsEvPVz0010","hEvPlaneTPCvsEvPVz0010",100,-0.5,2*TMath::Pi()+0.5,100,-0.5,2*TMath::Pi()+0.5);                      
  hEvPlaneTPCvsEvPVz1020 = new TH2F("hEvPlaneTPCvsEvPVz1020","hEvPlaneTPCvsEvPVz1020",100,-0.5,2*TMath::Pi()+0.5,100,-0.5,2*TMath::Pi()+0.5); 
  hEvPlaneTPCvsEvPVz2030 = new TH2F("hEvPlaneTPCvsEvPVz2030","hEvPlaneTPCvsEvPVz2030",100,-0.5,2*TMath::Pi()+0.5,100,-0.5,2*TMath::Pi()+0.5);
  hEvPlaneTPCvsEvPVz3040 = new TH2F("hEvPlaneTPCvsEvPVz3040","hEvPlaneTPCvsEvPVz3040",100,-0.5,2*TMath::Pi()+0.5,100,-0.5,2*TMath::Pi()+0.5);                      
  hEvPlaneTPCvsEvPVz4050 = new TH2F("hEvPlaneTPCvsEvPVz4050","hEvPlaneTPCvsEvPVz4050",100,-0.5,2*TMath::Pi()+0.5,100,-0.5,2*TMath::Pi()+0.5);                      
  hEvPlaneTPCvsEvPVz5060 = new TH2F("hEvPlaneTPCvsEvPVz5060","hEvPlaneTPCvsEvPVz5060",100,-0.5,2*TMath::Pi()+0.5,100,-0.5,2*TMath::Pi()+0.5);   
  hEvPlaneTPCvsEvPVz6080 = new TH2F("hEvPlaneTPCvsEvPVz6080","hEvPlaneTPCvsEvPVz6080",100,-0.5,2*TMath::Pi()+0.5,100,-0.5,2*TMath::Pi()+0.5);   
  
  fListHist->Add(  hEvPlaneTPCvsEvPVz0010 );                    
  fListHist->Add(  hEvPlaneTPCvsEvPVz1020 );
  fListHist->Add(  hEvPlaneTPCvsEvPVz2030 );
  fListHist->Add(  hEvPlaneTPCvsEvPVz3040 );                      
  fListHist->Add(  hEvPlaneTPCvsEvPVz4050 );                      
  fListHist->Add(  hEvPlaneTPCvsEvPVz5060 );   
  fListHist->Add(  hEvPlaneTPCvsEvPVz6080 );                      

  hCos2DeltaTPCVzAvsCentrality   = new TH2F("hCos2DeltaTPCVzAvsCentrality"  ,"hCos2DeltaTPCVzAvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaTPCVzCvsCentrality   = new TH2F("hCos2DeltaTPCVzCvsCentrality"  ,"hCos2DeltaTPCVzCvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzAVzCvsCentrality   = new TH2F("hCos2DeltaVzAVzCvsCentrality"  ,"hCos2DeltaVzAVzCvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzMVzAvsCentrality   = new TH2F("hCos2DeltaVzMVzAvsCentrality"  ,"hCos2DeltaVzMVzAvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzMVzCvsCentrality   = new TH2F("hCos2DeltaVzMVzCvsCentrality"  ,"hCos2DeltaVzMVzCvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzATPCvsCentrality   = new TH2F("hCos2DeltaVzATPCvsCentrality"  ,"hCos2DeltaVzATPCvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzCTPCvsCentrality   = new TH2F("hCos2DeltaVzCTPCvsCentrality"  ,"hCos2DeltaVzCTPCvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzCVzAvsCentrality   = new TH2F("hCos2DeltaVzCVzAvsCentrality"  ,"hCos2DeltaVzCVzAvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzMTPCpvsCentrality  = new TH2F("hCos2DeltaVzMTPCpvsCentrality" ,"hCos2DeltaVzMTPCpvsCentrality" ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzMTPCnvsCentrality  = new TH2F("hCos2DeltaVzMTPCnvsCentrality" ,"hCos2DeltaVzMTPCnvsCentrality" ,100,-1.1,1.1,105,0,105);
  hCos2DeltaTPCpTPCnvsCentrality = new TH2F("hCos2DeltaTPCpTPCnvsCentrality","hCos2DeltaTPCpTPCnvsCentrality",100,-1.1,1.1,105,0,105);

  fListHist->Add(hCos2DeltaTPCVzAvsCentrality);
  fListHist->Add(hCos2DeltaTPCVzCvsCentrality);
  fListHist->Add(hCos2DeltaVzAVzCvsCentrality);
  fListHist->Add(hCos2DeltaVzMVzAvsCentrality);
  fListHist->Add(hCos2DeltaVzMVzCvsCentrality);
  fListHist->Add(hCos2DeltaVzATPCvsCentrality);
  fListHist->Add(hCos2DeltaVzCTPCvsCentrality);
  fListHist->Add(hCos2DeltaVzCVzAvsCentrality);
  fListHist->Add(hCos2DeltaVzMTPCpvsCentrality);  
  fListHist->Add(hCos2DeltaVzMTPCnvsCentrality); 
  fListHist->Add(hCos2DeltaTPCpTPCnvsCentrality);

  if(fHarmonic < 3)
    hQVzAQVzCvsCentrality = new TH2F("hQVzAQVzCvsCentrality","hQVzAQVzCvsCentrality",1000,-100,100,105,0,105);
  else
    hQVzAQVzCvsCentrality = new TH2F("hQVzAQVzCvsCentrality","hQVzAQVzCvsCentrality",5000,-1000,1000,105,0,105);
 
  fListHist->Add(hQVzAQVzCvsCentrality);
  if(fHarmonic < 3){
    hQxVzAvsCentrality = new TH2F("hQxVzAvsCentrality","hQxVzAvsCentrality",100,-20,20,105,0,105);
    hQyVzAvsCentrality = new TH2F("hQyVzAvsCentrality","hQyVzAvsCentrality",100,-20,20,105,0,105);
    hQxVzCvsCentrality = new TH2F("hQxVzCvsCentrality","hQxVzCvsCentrality",100,-20,20,105,0,105);
    hQyVzCvsCentrality = new TH2F("hQyVzCvsCentrality","hQyVzCvsCentrality",100,-20,20,105,0,105);
    hQxVzMvsCentrality = new TH2F("hQxVzMvsCentrality","hQxVzMvsCentrality",100,-20,20,105,0,105);
    hQyVzMvsCentrality = new TH2F("hQyVzMvsCentrality","hQyVzMvsCentrality",100,-20,20,105,0,105);
  }
  
  else{
    hQxVzAvsCentrality = new TH2F("hQxVzAvsCentrality","hQxVzAvsCentrality",2000,-500,500,105,0,105);
    hQyVzAvsCentrality = new TH2F("hQyVzAvsCentrality","hQyVzAvsCentrality",2000,-500,500,105,0,105);
    hQxVzCvsCentrality = new TH2F("hQxVzCvsCentrality","hQxVzCvsCentrality",2000,-500,500,105,0,105);
    hQyVzCvsCentrality = new TH2F("hQyVzCvsCentrality","hQyVzCvsCentrality",2000,-500,500,105,0,105);
    hQxVzMvsCentrality = new TH2F("hQxVzMvsCentrality","hQxVzMvsCentrality",2000,-500,500,105,0,105);
    hQyVzMvsCentrality = new TH2F("hQyVzMvsCentrality","hQyVzMvsCentrality",2000,-500,500,105,0,105);
  }

  fListHist->Add(hQxVzAvsCentrality);
  fListHist->Add(hQyVzAvsCentrality);
  fListHist->Add(hQxVzCvsCentrality);
  fListHist->Add(hQyVzCvsCentrality);
  fListHist->Add(hQxVzMvsCentrality);
  fListHist->Add(hQyVzMvsCentrality);

  
  hCosdeltaphiTPCvsPtAll0010  = new TH2F("hCosdeltaphiTPCvsPtAll0010","hCosdeltaphiTPCvsPtAll0010" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiV0AvsPtAll0010  = new TH2F("hCosdeltaphiV0AvsPtAll0010","hCosdeltaphiV0AvsPtAll0010" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiV0CvsPtAll0010  = new TH2F("hCosdeltaphiV0CvsPtAll0010","hCosdeltaphiV0CvsPtAll0010" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiTPCvsPtAll1020  = new TH2F("hCosdeltaphiTPCvsPtAll1020","hCosdeltaphiTPCvsPtAll1020" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiV0AvsPtAll1020  = new TH2F("hCosdeltaphiV0AvsPtAll1020","hCosdeltaphiV0AvsPtAll1020" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiV0CvsPtAll1020  = new TH2F("hCosdeltaphiV0CvsPtAll1020","hCosdeltaphiV0CvsPtAll1020" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiTPCvsPtAll2030  = new TH2F("hCosdeltaphiTPCvsPtAll2030","hCosdeltaphiTPCvsPtAll2030" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiV0AvsPtAll2030  = new TH2F("hCosdeltaphiV0AvsPtAll2030","hCosdeltaphiV0AvsPtAll2030" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiV0CvsPtAll2030  = new TH2F("hCosdeltaphiV0CvsPtAll2030","hCosdeltaphiV0CvsPtAll2030" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiTPCvsPtAll3040  = new TH2F("hCosdeltaphiTPCvsPtAll3040","hCosdeltaphiTPCvsPtAll3040" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiV0AvsPtAll3040  = new TH2F("hCosdeltaphiV0AvsPtAll3040","hCosdeltaphiV0AvsPtAll3040" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiV0CvsPtAll3040  = new TH2F("hCosdeltaphiV0CvsPtAll3040","hCosdeltaphiV0CvsPtAll3040" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiTPCvsPtAll4050  = new TH2F("hCosdeltaphiTPCvsPtAll4050","hCosdeltaphiTPCvsPtAll4050" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiV0AvsPtAll4050  = new TH2F("hCosdeltaphiV0AvsPtAll4050","hCosdeltaphiV0AvsPtAll4050" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiV0CvsPtAll4050  = new TH2F("hCosdeltaphiV0CvsPtAll4050","hCosdeltaphiV0CvsPtAll4050" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiTPCvsPtAll5060  = new TH2F("hCosdeltaphiTPCvsPtAll5060","hCosdeltaphiTPCvsPtAll5060" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiV0AvsPtAll5060  = new TH2F("hCosdeltaphiV0AvsPtAll5060","hCosdeltaphiV0AvsPtAll5060" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiV0CvsPtAll5060  = new TH2F("hCosdeltaphiV0CvsPtAll5060","hCosdeltaphiV0CvsPtAll5060" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiTPCvsPtAll6080  = new TH2F("hCosdeltaphiTPCvsPtAll6080","hCosdeltaphiTPCvsPtAll6080" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiV0AvsPtAll6080  = new TH2F("hCosdeltaphiV0AvsPtAll6080","hCosdeltaphiV0AvsPtAll6080" ,100,-1.1,1.1,400,0,20);      
  hCosdeltaphiV0CvsPtAll6080  = new TH2F("hCosdeltaphiV0CvsPtAll6080","hCosdeltaphiV0CvsPtAll6080" ,100,-1.1,1.1,400,0,20);      


  fListHist->Add(hCosdeltaphiTPCvsPtAll0010);
  fListHist->Add(hCosdeltaphiV0AvsPtAll0010);
  fListHist->Add(hCosdeltaphiV0CvsPtAll0010);
  fListHist->Add(hCosdeltaphiTPCvsPtAll1020);
  fListHist->Add(hCosdeltaphiV0AvsPtAll1020);
  fListHist->Add(hCosdeltaphiV0CvsPtAll1020);
  fListHist->Add(hCosdeltaphiTPCvsPtAll2030);
  fListHist->Add(hCosdeltaphiV0AvsPtAll2030);
  fListHist->Add(hCosdeltaphiV0CvsPtAll2030);
  fListHist->Add(hCosdeltaphiTPCvsPtAll3040);
  fListHist->Add(hCosdeltaphiV0AvsPtAll3040);
  fListHist->Add(hCosdeltaphiV0CvsPtAll3040);
  fListHist->Add(hCosdeltaphiTPCvsPtAll4050);
  fListHist->Add(hCosdeltaphiV0AvsPtAll4050);
  fListHist->Add(hCosdeltaphiV0CvsPtAll4050);
  fListHist->Add(hCosdeltaphiTPCvsPtAll5060);
  fListHist->Add(hCosdeltaphiV0AvsPtAll5060);
  fListHist->Add(hCosdeltaphiV0CvsPtAll5060);
  fListHist->Add(hCosdeltaphiTPCvsPtAll6080);
  fListHist->Add(hCosdeltaphiV0AvsPtAll6080);
  fListHist->Add(hCosdeltaphiV0CvsPtAll6080);
  
  //phi distro charged particle

   
  hphivsPtAll0010  = new TH2F("hphivsPtAll0010","hphivsPtAll0010",100,0, TMath::Pi(),400,0,20);      
  hphivsPtAll1020  = new TH2F("hphivsPtAll1020","hphivsPtAll1020",100,0, TMath::Pi(),400,0,20);      
  hphivsPtAll2030  = new TH2F("hphivsPtAll2030","hphivsPtAll2030",100,0, TMath::Pi(),400,0,20);      
  hphivsPtAll3040  = new TH2F("hphivsPtAll3040","hphivsPtAll3040",100,0, TMath::Pi(),400,0,20);      
  hphivsPtAll4050  = new TH2F("hphivsPtAll4050","hphivsPtAll4050",100,0, TMath::Pi(),400,0,20);      
  hphivsPtAll5060  = new TH2F("hphivsPtAll5060","hphivsPtAll5060",100,0, TMath::Pi(),400,0,20);      
  hphivsPtAll6080  = new TH2F("hphivsPtAll6080","hphivsPtAll6080",100,0, TMath::Pi(),400,0,20);      

  fListHist->Add(hphivsPtAll0010); 
  fListHist->Add(hphivsPtAll1020); 
  fListHist->Add(hphivsPtAll2030); 
  fListHist->Add(hphivsPtAll3040); 
  fListHist->Add(hphivsPtAll4050); 
  fListHist->Add(hphivsPtAll5060); 
  fListHist->Add(hphivsPtAll6080); 

  huqV0AvsPtAll0010 = new TH2F("huqV0AvsPtAll0010","huqV0AvsPtAll0010",100,-20,20,400,0,20);
  huqV0CvsPtAll0010 = new TH2F("huqV0CvsPtAll0010","huqV0CvsPtAll0010",100,-20,20,400,0,20);
  huqV0AvsPtAll1020 = new TH2F("huqV0AvsPtAll1020","huqV0AvsPtAll1020",100,-20,20,400,0,20);
  huqV0CvsPtAll1020 = new TH2F("huqV0CvsPtAll1020","huqV0CvsPtAll1020",100,-20,20,400,0,20);
  huqV0AvsPtAll2030 = new TH2F("huqV0AvsPtAll2030","huqV0AvsPtAll2030",100,-20,20,400,0,20);
  huqV0CvsPtAll2030 = new TH2F("huqV0CvsPtAll2030","huqV0CvsPtAll2030",100,-20,20,400,0,20);
  huqV0AvsPtAll3040 = new TH2F("huqV0AvsPtAll3040","huqV0AvsPtAll3040",100,-20,20,400,0,20);
  huqV0CvsPtAll3040 = new TH2F("huqV0CvsPtAll3040","huqV0CvsPtAll3040",100,-20,20,400,0,20);
  huqV0AvsPtAll4050 = new TH2F("huqV0AvsPtAll4050","huqV0AvsPtAll4050",100,-20,20,400,0,20);
  huqV0CvsPtAll4050 = new TH2F("huqV0CvsPtAll4050","huqV0CvsPtAll4050",100,-20,20,400,0,20);
  huqV0AvsPtAll5060 = new TH2F("huqV0AvsPtAll5060","huqV0AvsPtAll5060",100,-20,20,400,0,20);
  huqV0CvsPtAll5060 = new TH2F("huqV0CvsPtAll5060","huqV0CvsPtAll5060",100,-20,20,400,0,20);
  huqV0AvsPtAll6080 = new TH2F("huqV0AvsPtAll6080","huqV0AvsPtAll6080",100,-20,20,400,0,20);
  huqV0CvsPtAll6080 = new TH2F("huqV0CvsPtAll6080","huqV0CvsPtAll6080",100,-20,20,400,0,20);

  fListHist->Add(huqV0AvsPtAll0010);
  fListHist->Add(huqV0CvsPtAll0010);
  fListHist->Add(huqV0AvsPtAll1020);
  fListHist->Add(huqV0CvsPtAll1020);
  fListHist->Add(huqV0AvsPtAll2030);
  fListHist->Add(huqV0CvsPtAll2030);
  fListHist->Add(huqV0AvsPtAll3040);
  fListHist->Add(huqV0CvsPtAll3040);
  fListHist->Add(huqV0AvsPtAll4050);
  fListHist->Add(huqV0CvsPtAll4050);
  fListHist->Add(huqV0AvsPtAll5060);
  fListHist->Add(huqV0CvsPtAll5060);
  fListHist->Add(huqV0AvsPtAll6080);
  fListHist->Add(huqV0CvsPtAll6080);

  // new hist

  huqV0AvsuqV0C0010  = new TH2F("huqV0AvsuqV0C0010","huqV0AvsuqV0C0010",100,-20,20,100,-20,20);
  huqV0AvsuqV0C1020  = new TH2F("huqV0AvsuqV0C1020","huqV0AvsuqV0C1020",100,-20,20,100,-20,20);
  huqV0AvsuqV0C2030  = new TH2F("huqV0AvsuqV0C2030","huqV0AvsuqV0C2030",100,-20,20,100,-20,20);
  huqV0AvsuqV0C3040  = new TH2F("huqV0AvsuqV0C3040","huqV0AvsuqV0C3040",100,-20,20,100,-20,20);
  huqV0AvsuqV0C4050  = new TH2F("huqV0AvsuqV0C4050","huqV0AvsuqV0C4050",100,-20,20,100,-20,20);
  huqV0AvsuqV0C5060  = new TH2F("huqV0AvsuqV0C5060","huqV0AvsuqV0C5060",100,-20,20,100,-20,20);
  huqV0AvsuqV0C6080  = new TH2F("huqV0AvsuqV0C6080","huqV0AvsuqV0C6080",100,-20,20,100,-20,20);
 
  huqV0AxuqV0CvsPtAll0010     = new TH2F("huqV0AxuqV0CvsPtAll0010","huqV0AxuqV0CvsPtAll0010",100,-20,20,400,0,20);
  huqV0AxuqV0CvsPtAll1020     = new TH2F("huqV0AxuqV0CvsPtAll1020","huqV0AxuqV0CvsPtAll1020",100,-20,20,400,0,20);
  huqV0AxuqV0CvsPtAll2030     = new TH2F("huqV0AxuqV0CvsPtAll2030","huqV0AxuqV0CvsPtAll2030",100,-20,20,400,0,20);
  huqV0AxuqV0CvsPtAll3040     = new TH2F("huqV0AxuqV0CvsPtAll3040","huqV0AxuqV0CvsPtAll3040",100,-20,20,400,0,20);
  huqV0AxuqV0CvsPtAll4050     = new TH2F("huqV0AxuqV0CvsPtAll4050","huqV0AxuqV0CvsPtAll4050",100,-20,20,400,0,20);
  huqV0AxuqV0CvsPtAll5060     = new TH2F("huqV0AxuqV0CvsPtAll5060","huqV0AxuqV0CvsPtAll5060",100,-20,20,400,0,20);
  huqV0AxuqV0CvsPtAll6080     = new TH2F("huqV0AxuqV0CvsPtAll6080","huqV0AxuqV0CvsPtAll6080",100,-20,20,400,0,20);
 
  fListHist->Add(huqV0AvsuqV0C0010);
  fListHist->Add(huqV0AvsuqV0C1020);
  fListHist->Add(huqV0AvsuqV0C2030);
  fListHist->Add(huqV0AvsuqV0C3040);
  fListHist->Add(huqV0AvsuqV0C4050);
  fListHist->Add(huqV0AvsuqV0C5060);
  fListHist->Add(huqV0AvsuqV0C6080);
  
  fListHist->Add(huqV0AxuqV0CvsPtAll0010);
  fListHist->Add(huqV0AxuqV0CvsPtAll1020);
  fListHist->Add(huqV0AxuqV0CvsPtAll2030);
  fListHist->Add(huqV0AxuqV0CvsPtAll3040);
  fListHist->Add(huqV0AxuqV0CvsPtAll4050);
  fListHist->Add(huqV0AxuqV0CvsPtAll5060);
  fListHist->Add(huqV0AxuqV0CvsPtAll6080);

  fEventCuts.AddQAplotsToList(fListHist);
  
  PostData(1,  fListHist);


}// end UserCreateOutputObjects


//====================== USER EXEC ========================

void AliAnalysisTaskAllPtcv2::UserExec(Option_t *) 
{
  // Main loop
  // Called for EACH event
  //  cout<<"AliAnalysisTaskAllPtcv2 Starting UserExec"<<endl;

  Info("AliAnalysisTaskAllPtcv2","Starting UserExec");  
 
  fHistEventMultiplicity->Fill(1);
    
  AliVEvent *event = InputEvent();
  if (!event) { Printf("ERROR: Could not retrieve event"); return; }


  if(fAnalysisType == "ESD"){
    fESDevent = dynamic_cast<AliESDEvent*>(event);
    if (!fESDevent) {
      AliError("Cannot get the ESD event");
      return;
    }  
    
    fevent = fESDevent;
  }

  else if(fAnalysisType == "AOD"){
    fAODevent = dynamic_cast<AliAODEvent*>(event);
    if (!fAODevent) {
      AliError("Cannot get the AOD event");
      return;
    }  
    fevent = fAODevent;
  }
  
  else{
    AliError("Cannot get any event");
    return;
  }

  cout<<"event: "<<fevent<<endl;
  //  cout<<"event cuts: "<<fEventCuts<<endl;
  cout<<"pass cuts: "<<fEventCuts.AcceptEvent(fevent)<<endl;

  /// Use the event cut class to apply the required selections
  if (!fEventCuts.AcceptEvent(fevent)) {
    PostData(1, fListHist);
    
    return;
  }
  
   
  //______________________________________________________
  // PID
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse=inputHandler->GetPIDResponse(); 

  Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
  AliCentrality* centrality = 0x0;
  
  //-----------------------------------------------------
  
  Float_t  impactXY=-999., impactZ=-999.;
  Double_t TPCSignal=0.;
  
  ULong_t  status=0;
 
  Double_t pmax  = 10.;
  Double_t ptmax = 6.2;

  // Primary vertex cut

  const AliVVertex* vertexmain = fevent->GetPrimaryVertex();
  if (!vertexmain){
    AliWarning("No prim. vertex in ESD... return!");
    
    PostData(1, fListHist);
    return;
  }
  
  vertexmain->GetXYZ( lBestPrimaryVtxPos );

  Bool_t isPileUpSpd=kFALSE;
  if(fAnalysisType == "ESD"){
    isPileUpSpd=fESDevent->IsPileupFromSPD();
  }
  else if(fAnalysisType == "AOD"){
    isPileUpSpd=fAODevent->IsPileupFromSPD();
  }
  if(isPileUpSpd){  
    PostData(1, fListHist);
    return;
  }
  
  fHistEventMultiplicity->Fill(2); // analyzed events with PV w/o pile up
  
  if((TMath::Abs(lBestPrimaryVtxPos[2])) > fVzmax) return;
  
  fHistEventMultiplicity->Fill(3);

  Bool_t isSelectedCentral     = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral);
  Bool_t isSelectedSemiCentral = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kSemiCentral);
  Bool_t isSelectedMB          = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  Bool_t isSelectedINT7        = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
  
  //_____________________________________________________
  //   Centrality 

  Float_t percentile = -999;
  // AliMultSelection *fMultSelection = 0x0; 

  if(fCentrality == "V0M")
    percentile = fEventCuts.GetCentrality(0);
  else if(fCentrality == "CL1")
    percentile = fEventCuts.GetCentrality(1);
  
  Int_t TrackNumber = fevent->GetNumberOfTracks();
  fHistTrackMultiplicity->Fill(TrackNumber,percentile); //tracce per evento
  
  Int_t eventtype = -999;
  
  if(isSelectedCentral){
    fHistEventMultiplicity->Fill(4);
    fHistTrackMultiplicityCentral->Fill(TrackNumber,percentile); 
    eventtype =1;
  }
  
  if(isSelectedSemiCentral){
    fHistEventMultiplicity->Fill(5);
    fHistTrackMultiplicitySemiCentral->Fill(TrackNumber,percentile); 
    eventtype =2;
  }
  
  if(isSelectedMB){
    if(percentile<0)return;
    if(percentile>=80)return;
    fHistEventMultiplicity->Fill(6);
    fHistTrackMultiplicityMB->Fill(TrackNumber,percentile); 
    eventtype =3;
  }
 
  if(isSelectedINT7){
    if(percentile<0)return;
    if(percentile>80)return;
    fHistEventMultiplicity->Fill(7);
    fHistTrackMultiplicityINT7->Fill(TrackNumber,percentile); 
    eventtype =4;
  }
  
  if(eventtype!=1 && eventtype!=2 && eventtype!=3 && eventtype!=4)return;
 
  
   // from D2H task
  
  Double_t eventplaneqncorrTPC[3];
  Double_t eventplaneqncorrVZERO[3];
  TList *qnlist = 0x0;

  /////////////////////////////////////////////////////////////
  //////////////          GET Qn vectors         //////////////
  /////////////////////////////////////////////////////////////
  fHistEventMultiplicity->Fill(8);
  
  AliQnCorrectionsManager *flowQnVectorMgr;
  AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections*>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
  if (flowQnVectorTask != NULL) {
    flowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
    fHistEventMultiplicity->Fill(9);
  }
  else {
    AliWarning("This task needs the Flow Qn vector corrections framework and it is not present. Aborting!!!\n");
      //fHistEventMultiplicity->Fill(10);
    return;
  }
  
  qnlist = flowQnVectorMgr->GetQnVectorList();
  if(qnlist)
    fHistEventMultiplicity->Fill(10);
  if (!qnlist) {
    AliWarning("Qn list is not present. Aborting!!!\n");
    return;
  }

  cout<<"flowQnVectorTask "<<flowQnVectorTask<<endl;

  fHistEventMultiplicity->Fill(11);

  const AliQnCorrectionsQnVector* qnVectTPC[3];
  const AliQnCorrectionsQnVector* qnVectV0[3];

  TString fDetTPCConfName[3];
  TString fDetV0ConfName[3];
  TString fNormMethod="QoverQlength";

  fDetTPCConfName[0] = "TPC";
  fDetTPCConfName[1] = "TPCNegEta";
  fDetTPCConfName[2] = "TPCPosEta";
    
  fDetV0ConfName[0]  = "VZERO";
  fDetV0ConfName[1]  = "VZEROA";
  fDetV0ConfName[2]  = "VZEROC";

  for(int iDet = 0; iDet < 3; iDet++) {
    qnVectTPC[iDet]  = GetQnVectorFromList(qnlist, Form("%s%s",fDetTPCConfName[iDet].Data(),fNormMethod.Data()), "latest", "plain");
    qnVectV0[iDet]   = GetQnVectorFromList(qnlist, Form("%s%s",fDetV0ConfName[iDet].Data(),fNormMethod.Data()), "latest", "raw");
    if(!qnVectTPC[iDet] || !qnVectV0[iDet]) return;
    eventplaneqncorrTPC[iDet]   = qnVectTPC[iDet]->EventPlane(fHarmonic);
    eventplaneqncorrVZERO[iDet] = qnVectV0[iDet]->EventPlane(fHarmonic);
    if(eventplaneqncorrTPC[iDet]<0.)   
      eventplaneqncorrTPC[iDet]    += 2.*(TMath::Pi())/fHarmonic;
    if(eventplaneqncorrVZERO[iDet]<0.) 
      eventplaneqncorrVZERO[iDet]  += 2.*(TMath::Pi())/fHarmonic;
   
    //  fHistEvPlaneQncorrTPC[iDet]->Fill(eventplaneqncorrTPC[iDet]);
    //  fHistEvPlaneQncorrVZERO[iDet]->Fill(eventplaneqncorrVZERO[iDet]);
    

  }
  
  fHistEventMultiplicity->Fill(10);

  Double_t qxEPa = 0, qyEPa = 0;
  Double_t qxEPc = 0, qyEPc = 0;
  Double_t qxEP =  0, qyEP  = 0;

  Double_t Qx2  = 0, Qy2  = 0;
  Double_t Qx2p = 0, Qy2p = 0;
  Double_t Qx2n = 0, Qy2n = 0;

  const AliQnCorrectionsQnVector* QnVectTPC  = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetTPCConfName[0].Data()), "latest", "plain");     
  const AliQnCorrectionsQnVector* QnVectTPCn = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetTPCConfName[1].Data()), "latest", "plain");     
  const AliQnCorrectionsQnVector* QnVectTPCp = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetTPCConfName[2].Data()), "latest", "plain");     
  const AliQnCorrectionsQnVector* QnVectV0   = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetV0ConfName[0].Data()), "latest", "raw");	   
  const AliQnCorrectionsQnVector* QnVectV0A  = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetV0ConfName[1].Data()), "latest", "raw");	   
  const AliQnCorrectionsQnVector* QnVectV0C  = GetQnVectorFromList(qnlist, Form("%sQoverSqrtM",fDetV0ConfName[2].Data()), "latest", "raw");        

  Qx2 = QnVectTPC->Qx(fHarmonic);
  Qy2 = QnVectTPC->Qy(fHarmonic);

  Qx2p  = QnVectTPCp->Qx(fHarmonic);
  Qy2p  = QnVectTPCp->Qy(fHarmonic);

  Qx2n  = QnVectTPCn->Qx(fHarmonic);
  Qy2n  = QnVectTPCn->Qy(fHarmonic);

  qxEPa = QnVectV0A->Qx(fHarmonic);
  qyEPa = QnVectV0A->Qy(fHarmonic);
  	  
  qxEPc = QnVectV0C->Qx(fHarmonic);
  qyEPc = QnVectV0C->Qy(fHarmonic);
  	  
  qxEP  = QnVectV0->Qx(fHarmonic);
  qyEP  = QnVectV0->Qy(fHarmonic);


  Double_t evPlAngV0A = eventplaneqncorrVZERO[1];
  Double_t evPlAngV0C = eventplaneqncorrVZERO[2];
  Double_t evPlAngV0  = eventplaneqncorrVZERO[0];

  Double_t evPlAngTPC  = eventplaneqncorrTPC[0];
  Double_t evPlAngTPCn = eventplaneqncorrTPC[1];
  Double_t evPlAngTPCp = eventplaneqncorrTPC[2];

  EPVzAvsCentrality  ->Fill( evPlAngV0A  , percentile); 
  EPVzCvsCentrality  ->Fill( evPlAngV0C  , percentile); 
  EPVzvsCentrality   ->Fill( evPlAngV0   , percentile); 

  EPTPCvsCentrality  ->Fill( evPlAngTPC  , percentile); 
  EPTPCpvsCentrality ->Fill( evPlAngTPCn , percentile); 
  EPTPCnvsCentrality ->Fill( evPlAngTPCp , percentile); 
 
  //--------------------------------------------------
    
  if(percentile>=0 && percentile<10)
    hEvPlaneTPCvsEvPVz0010 ->Fill(evPlAngTPC,evPlAngV0); 
  if(percentile>=10 && percentile<20)
    hEvPlaneTPCvsEvPVz1020 ->Fill(evPlAngTPC,evPlAngV0); 
  if(percentile>=20 && percentile<30)
    hEvPlaneTPCvsEvPVz2030->Fill(evPlAngTPC,evPlAngV0);
  if(percentile>=30 && percentile<40)
    hEvPlaneTPCvsEvPVz3040->Fill(evPlAngTPC,evPlAngV0);  
  if(percentile>=40 && percentile<50)                    
    hEvPlaneTPCvsEvPVz4050->Fill(evPlAngTPC,evPlAngV0);   
  if(percentile>=50 && percentile<60)                   
    hEvPlaneTPCvsEvPVz5060->Fill(evPlAngTPC,evPlAngV0);           
  if(percentile>=60 && percentile<80)                   
    hEvPlaneTPCvsEvPVz6080->Fill(evPlAngTPC,evPlAngV0);           
  
  // For TPC, V0M, V0c and V0A resolution 

  hCos2DeltaTPCVzAvsCentrality  ->Fill(TMath::Cos(fHarmonic*(evPlAngTPC - evPlAngV0A)) , percentile);
  hCos2DeltaTPCVzCvsCentrality  ->Fill(TMath::Cos(fHarmonic*(evPlAngTPC - evPlAngV0C)) , percentile);
  hCos2DeltaVzAVzCvsCentrality  ->Fill(TMath::Cos(fHarmonic*(evPlAngV0A - evPlAngV0C)) , percentile);
  hCos2DeltaVzMVzAvsCentrality  ->Fill(TMath::Cos(fHarmonic*(evPlAngV0  - evPlAngV0A)) , percentile);
  hCos2DeltaVzMVzCvsCentrality  ->Fill(TMath::Cos(fHarmonic*(evPlAngV0  - evPlAngV0C)) , percentile);
  hCos2DeltaVzATPCvsCentrality  ->Fill(TMath::Cos(fHarmonic*(evPlAngV0A - evPlAngTPC)) , percentile);
  hCos2DeltaVzCTPCvsCentrality  ->Fill(TMath::Cos(fHarmonic*(evPlAngV0C - evPlAngTPC)) , percentile);
  hCos2DeltaVzCVzAvsCentrality  ->Fill(TMath::Cos(fHarmonic*(evPlAngV0C - evPlAngV0A)) , percentile);
  hCos2DeltaVzMTPCpvsCentrality ->Fill(TMath::Cos(fHarmonic*(evPlAngV0  - evPlAngTPCp)), percentile);
  hCos2DeltaVzMTPCnvsCentrality ->Fill(TMath::Cos(fHarmonic*(evPlAngV0  - evPlAngTPCn)), percentile);
  hCos2DeltaTPCpTPCnvsCentrality->Fill(TMath::Cos(fHarmonic*(evPlAngTPCp- evPlAngTPCn)), percentile);

  //Scalar Product
  
  Double_t  QV0AQV0C = qxEPa * qxEPc + qyEPa*qyEPc;
  hQVzAQVzCvsCentrality->Fill(QV0AQV0C,percentile);
  
  //NUA correction
 
  hQxVzAvsCentrality->Fill(qxEPa,percentile);
  hQyVzAvsCentrality->Fill(qyEPa,percentile);
  hQxVzCvsCentrality->Fill(qxEPc,percentile);
  hQyVzCvsCentrality->Fill(qyEPc,percentile);
  hQxVzMvsCentrality->Fill(qxEP ,percentile);
  hQyVzMvsCentrality->Fill(qyEP ,percentile);

  //====================================================================================================================
  
  Double_t uqV0A      = -999.;
  Double_t uqV0C      = -999.; 
  Double_t deltaphiTPC= -999.;
  Double_t deltaphiV0 = -999.;
  Double_t deltaphiV0A= -999.;
  Double_t deltaphiV0C= -999.;

  Double_t eta       = -999.;
  Double_t pt        = -999.;
  Double_t tPhi      = -999.;
  Double_t ptpc      = -999.;
  Double_t sign      = -999.;

  Double_t  beta  = 0;
  Double_t  p     = 0;
  Double_t  gamma = 0;
  Double_t  mass  = -99;
  Double_t  tof   = -99;



  for (Int_t j=0; j<TrackNumber; j++) { //loop on tracks
 
    uqV0A     = -999.;
    uqV0C     = -999.; 
    eta       = -999.;
    pt        = -999.;
    TPCSignal = -999.;
    beta      = 0.;
    gamma     = 0.;
    mass      = -9999.;
    p         = -9999.;
    tof       = -9999.;


    AliVTrack* atrack = (AliVTrack*) fevent->GetTrack(j);
    if (!atrack)
      continue;
    
    Bool_t trkFlag = 0;
    if(fAnalysisType == "ESD"){
      trkFlag = fESDtrackCuts->AcceptTrack((AliESDtrack*)atrack);
    }
    else if(fAnalysisType == "AOD"){
      trkFlag = ((AliAODTrack*)atrack)->TestFilterBit(4);
    }

    if(!trkFlag)continue;

    status  = (ULong_t)atrack->GetStatus();
    
    Bool_t hasTOFout  = status&AliVTrack::kTOFout; 
    Bool_t hasTOFtime = status&AliVTrack::kTIME;
    Bool_t hasTOF     = kFALSE;
    if (hasTOFout && hasTOFtime) 
      hasTOF = kTRUE;
    Float_t length = atrack->GetIntegratedLength(); 
    if (length < 350.) hasTOF = kFALSE;
    
    TPCSignal=atrack->GetTPCsignal(); 
    
    if(TPCSignal<10)continue;
    if(TPCSignal>1000)continue;
                     
    ptpc = atrack->GetTPCmomentum(); // momentum for dEdx determination
    pt   = atrack->Pt();
    p    = atrack->P();
    sign = atrack->Charge();
    
    eta       = atrack->Eta();
    pt        = atrack->Pt();
    tPhi      = atrack->Phi();

    if ((TMath::Abs(eta) > 0.5) || (pt < 0.2) || (pt > 20.0))continue;

    fhBB->Fill(ptpc*sign,TPCSignal);

    if(fAnalysisType == "ESD"){
      AliESDtrack *aesdtrack = static_cast<AliESDtrack *>(atrack);
      aesdtrack->GetImpactParameters(impactXY, impactZ);
    }
    else if(fAnalysisType == "AOD"){
      Double_t d[2], covd[3];
      AliAODTrack* track_clone=(AliAODTrack*)atrack->Clone("track_clone"); // need to clone because PropagateToDCA updates the track parameters
      Bool_t isDCA = track_clone->PropagateToDCA(fevent->GetPrimaryVertex(),fevent->GetMagneticField(),9999.,d,covd);
      delete track_clone;
      if(!isDCA)d[0]=-999.;
      impactXY = d[0];
      impactZ  = d[1];
    }
    

    if (hasTOF) {
      
      beta  = length / (2.99792457999999984e-02 * tof);
      gamma = 1/TMath::Sqrt(1 - beta*beta);
      mass  = p/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
      fhMassTOF->Fill(mass);
      fhTOF->Fill(p*sign,beta);
    }
      
    // Event Plane
    // Remove AutoCorrelation
    
    deltaphiTPC=TMath::Cos(fHarmonic*GetPhi0Pi(tPhi-evPlAngTPC));
    deltaphiV0 =TMath::Cos(fHarmonic*GetPhi0Pi(tPhi-evPlAngV0 ));
    deltaphiV0A=TMath::Cos(fHarmonic*GetPhi0Pi(tPhi-evPlAngV0A));
    deltaphiV0C=TMath::Cos(fHarmonic*GetPhi0Pi(tPhi-evPlAngV0C));
    
    // Scalar Product
    
    uqV0A = TMath::Cos(fHarmonic*tPhi)*qxEPa+TMath::Sin(fHarmonic*tPhi)*qyEPa;
    uqV0C = TMath::Cos(fHarmonic*tPhi)*qxEPc+TMath::Sin(fHarmonic*tPhi)*qyEPc;
    
    if(percentile>=0 && percentile<10) {
      
      hCosdeltaphiTPCvsPtAll0010->Fill(deltaphiTPC,pt); 
      hCosdeltaphiV0AvsPtAll0010->Fill(deltaphiV0A,pt); 
      hCosdeltaphiV0CvsPtAll0010->Fill(deltaphiV0C,pt); 
      huqV0AvsPtAll0010         ->Fill(uqV0A,pt); 
      huqV0CvsPtAll0010         ->Fill(uqV0C,pt); 
      hphivsPtAll0010           ->Fill(tPhi,pt);  
      
      huqV0AvsuqV0C0010         ->Fill(uqV0A,uqV0C); 
      huqV0AxuqV0CvsPtAll0010   ->Fill(uqV0A*uqV0C,pt); 
    }
    
    if(percentile>=10 && percentile<20){
      
      hCosdeltaphiTPCvsPtAll1020->Fill(deltaphiTPC,pt); 
      hCosdeltaphiV0AvsPtAll1020->Fill(deltaphiV0A,pt); 
      hCosdeltaphiV0CvsPtAll1020->Fill(deltaphiV0C,pt); 
      huqV0AvsPtAll1020         ->Fill(uqV0A,pt); 
      huqV0CvsPtAll1020         ->Fill(uqV0C,pt); 
      hphivsPtAll1020 ->Fill(tPhi,pt);  
	
      huqV0AvsuqV0C1020        ->Fill(uqV0A,uqV0C); 
      huqV0AxuqV0CvsPtAll1020  ->Fill(uqV0A*uqV0C,pt); 
    }

    if(percentile>=20 && percentile<30){
      
      hCosdeltaphiTPCvsPtAll2030->Fill(deltaphiTPC,pt); 
      hCosdeltaphiV0AvsPtAll2030->Fill(deltaphiV0A,pt); 
      hCosdeltaphiV0CvsPtAll2030->Fill(deltaphiV0C,pt); 
      huqV0AvsPtAll2030         ->Fill(uqV0A,pt); 
      huqV0CvsPtAll2030         ->Fill(uqV0C,pt); 
      hphivsPtAll2030 ->Fill(tPhi,pt);
   
      huqV0AvsuqV0C2030        ->Fill(uqV0A,uqV0C); 
      huqV0AxuqV0CvsPtAll2030  ->Fill(uqV0A*uqV0C,pt); 
    }
    
    if(percentile>=30 && percentile<40){    
      
      hCosdeltaphiTPCvsPtAll3040->Fill(deltaphiTPC,pt); 
      hCosdeltaphiV0AvsPtAll3040->Fill(deltaphiV0A,pt); 
      hCosdeltaphiV0CvsPtAll3040->Fill(deltaphiV0C,pt); 
      huqV0AvsPtAll3040         ->Fill(uqV0A,pt); 
      huqV0CvsPtAll3040         ->Fill(uqV0C,pt); 
      hphivsPtAll3040 ->Fill(tPhi,pt);  
      
      huqV0AvsuqV0C3040        ->Fill(uqV0A,uqV0C); 
      huqV0AxuqV0CvsPtAll3040  ->Fill(uqV0A*uqV0C,pt); 
      
    }
    
    if(percentile>=40 && percentile<50){      
      
      hCosdeltaphiTPCvsPtAll4050->Fill(deltaphiTPC,pt); 
      hCosdeltaphiV0AvsPtAll4050->Fill(deltaphiV0A,pt); 
      hCosdeltaphiV0CvsPtAll4050->Fill(deltaphiV0C,pt); 
      huqV0AvsPtAll4050         ->Fill(uqV0A,pt); 
      huqV0CvsPtAll4050         ->Fill(uqV0C,pt); 
      hphivsPtAll4050 ->Fill(tPhi,pt);
      
      huqV0AvsuqV0C4050        ->Fill(uqV0A,uqV0C); 
      huqV0AxuqV0CvsPtAll4050  ->Fill(uqV0A*uqV0C,pt); 
      
    }

    if(percentile>=50 && percentile<60){      
      
      hCosdeltaphiTPCvsPtAll5060->Fill(deltaphiTPC,pt); 
      hCosdeltaphiV0AvsPtAll5060->Fill(deltaphiV0A,pt); 
      hCosdeltaphiV0CvsPtAll5060->Fill(deltaphiV0C,pt); 
      huqV0AvsPtAll5060         ->Fill(uqV0A,pt); 
      huqV0CvsPtAll5060         ->Fill(uqV0C,pt); 
      hphivsPtAll5060 ->Fill(tPhi,pt);
      
      huqV0AvsuqV0C5060        ->Fill(uqV0A,uqV0C); 
      huqV0AxuqV0CvsPtAll5060  ->Fill(uqV0A*uqV0C,pt); 
      
    }

    if(percentile>=60 && percentile<80){      
      
      hCosdeltaphiTPCvsPtAll6080->Fill(deltaphiTPC,pt); 
      hCosdeltaphiV0AvsPtAll6080->Fill(deltaphiV0A,pt); 
      hCosdeltaphiV0CvsPtAll6080->Fill(deltaphiV0C,pt); 
      huqV0AvsPtAll6080         ->Fill(uqV0A,pt); 
      huqV0CvsPtAll6080         ->Fill(uqV0C,pt); 
      hphivsPtAll6080 ->Fill(tPhi,pt);
      
      huqV0AvsuqV0C6080        ->Fill(uqV0A,uqV0C); 
      huqV0AxuqV0CvsPtAll6080  ->Fill(uqV0A*uqV0C,pt); 
      
    }

    //}//only if tof
  }  //track
  
  PostData(1, fListHist);
} //end userexec


//________________________________________________________________________

void AliAnalysisTaskAllPtcv2::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
}

