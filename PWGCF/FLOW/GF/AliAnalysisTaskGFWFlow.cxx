/*
Author: Vytautas Vislavicius
Extention of Generic Flow (https://arxiv.org/abs/1312.3572)
*/
#include "AliAnalysisTaskGFWFlow.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisFilter.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliVVertex.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliAODMCParticle.h"
#include "AliExternalTrackParam.h"
#include "AliMultSelection.h"
#include "AliInputEventHandler.h"
#include "TList.h"
#include "TProfile.h"
#include "TClonesArray.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "TParticle.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliGFWWeights.h"
#include "AliGFWFlowContainer.h"
#include "TObjArray.h"
#include "TNamed.h"
#include "AliGFW.h"
#include "TRandom.h"
#include <vector>
#include "AliCollisionGeometry.h"
#include "AliGenEventHeader.h"
#include "AliGenHepMCEventHeader.h"

ClassImp(AliAnalysisTaskGFWFlow);

AliAnalysisTaskGFWFlow::AliAnalysisTaskGFWFlow():
  AliAnalysisTaskSE(),
  fTriggerType(AliVEvent::kINT7),
  fProduceWeights(kTRUE),
  fWeightList(0),
  fPeriod("LHC18q_pass3"),
  fCentMap(0),
  fWeights(0),
  fFC(0),
  fMCEvent(0),
  fIsMC(kFALSE),
  fIsTrain(kFALSE),
  fEvNomFlag(1),
  fTrNomFlag(1),
  fPtAxis(0),
  fPOIpTMin(0.2),
  fPOIpTMax(20),
  fRFpTMin(0.2),
  fRFpTMax(3.0),
  fTotFlags(0), //Total number of flags: 1 (nominal) + fTotTrackFlags + N_Event flags
  fTotTrackFlags(0), //Total number of track flags (without nominal)
  fTotEvFlags(0),
  fRunNo(-1),
  fBypassCalculations(kFALSE),
  fMultiDist(0),
  fCollisionsSystem(2)
{
};
AliAnalysisTaskGFWFlow::AliAnalysisTaskGFWFlow(const char *name, Bool_t ProduceWeights, Bool_t IsMC, Bool_t IsTrain):
  AliAnalysisTaskSE(name),
  fTriggerType(AliVEvent::kINT7),
  fProduceWeights(ProduceWeights),
  fWeightList(0),
  fPeriod("LHC18q_pass3"),
  fCentMap(0),
  fWeights(0),
  fFC(0),
  fIsMC(IsMC),
  fIsTrain(IsTrain),
  fEvNomFlag(1),
  fTrNomFlag(1),
  fPtAxis(new TAxis()),
  fPOIpTMin(0.2),
  fPOIpTMax(20),
  fRFpTMin(0.2),
  fRFpTMax(3.0),
  fTotFlags(0),
  fTotTrackFlags(0),
  fTotEvFlags(0),
  fRunNo(-1),
  fBypassCalculations(kFALSE),
  fMultiDist(0),
  fCollisionsSystem(2)
{
  if(!fProduceWeights) {
    if(fIsTrain) DefineInput(1,TH1D::Class());
    else DefineInput(1,TList::Class());
  }
  DefineOutput(1,(fProduceWeights?TList::Class():AliGFWFlowContainer::Class()));
  DefineOutput(2,TH1D::Class());
};
AliAnalysisTaskGFWFlow::~AliAnalysisTaskGFWFlow() {
};
void AliAnalysisTaskGFWFlow::UserCreateOutputObjects(){
  OpenFile(1);
  if(!fTotTrackFlags) fTotTrackFlags = gNTrackFlags;//If not defined explicitly, use the global values
  if(!fTotEvFlags) fTotEvFlags = gNEventFlags;//If not defined explicitly, use the global values
  fTotFlags = fTotEvFlags + fTotTrackFlags;//AliGFWCuts::gNEventFlags+1;
  if(fProduceWeights) {
    //Initialize selection objects
    fWeightList = new TList();
    fWeightList->SetName("WeightList");
    fWeightList->SetOwner(kTRUE);

    //RunNumbers for different periods
    vector<int> temp = {};
    if (fPeriod.EqualTo("LHC15o")) {
      vector<int> temp = {244917, 244918, 244975, 244980, 244982, 244983, 245061, 245064, 245066, 245068, 245145, 245146, 245148, 245151, 245152, 245231, 245232, 245233, 245259, 245343, 245345, 245346, 245347, 245349, 245353, 245396, 245397, 245401, 245407, 245409, 245410, 245411, 245439, 245441, 245446, 245450, 245452, 245453, 245454, 245496, 245497, 245501, 245504, 245505, 245507, 245535, 245540, 245542, 245543, 245544, 245545, 245554, 245683, 245692, 245700, 245702, 245705, 245729, 245731, 245738, 245752, 245759, 245766, 245775, 245785, 245793, 245829, 245831, 245833, 245923, 245949, 245952, 245954, 245963, 246001, 246003, 246012, 246036, 246037, 246042, 246048, 246049, 246052, 246053, 246087, 246089, 246113, 246115, 246148, 246151, 246152, 246153, 246178, 246180, 246181, 246182, 246185, 246217, 246222, 246225, 246271, 246272, 246275, 246276, 246390, 246391, 246392, 246424, 246428, 246431, 246434, 246487, 246488, 246493, 246495, 246540, 246543, 246553, 246567, 246568, 246575, 246583, 246648, 246671, 246675, 246676, 246750, 246751, 246757, 246758, 246759, 246760, 246763, 246765, 246766, 246804, 246805, 246807, 246808, 246809, 246810, 246844, 246845, 246846, 246847, 246851, 246855, 246858, 246859, 246864, 246865, 246867, 246870, 246871, 246928, 246930, 246937, 246942, 246945, 246948, 246949, 246980, 246982, 246984, 246989, 246991, 246994};
      RunNumber = temp;
    } 
    else if (fPeriod.EqualTo("LHC15o_pass2")) {
      vector<int> temp = {246994, 246991, 246989, 246984, 246982, 246948, 246945, 246928, 246871, 246870, 246867, 246865, 246864, 246859, 246858, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246434, 246431, 246424, 246392, 246391, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245793, 245785, 245775, 245766, 245759, 245752, 245731, 245729, 245705, 245702, 245692, 245683, 245554, 245545, 245544, 245543, 245542, 245540, 245535, 245507, 245505, 245504, 245501, 245497, 245496, 245454, 245453, 245450, 245446, 245441, 245411, 245410, 245409, 245407, 245401, 245397, 245396, 245353, 245349, 245347, 245346, 245345, 245343, 245259, 245233, 245232, 245231, 245152, 245151, 245146, 245145, 245068, 245066, 245064, 244983, 244982, 244980, 244975, 244918, 244917};
      RunNumber = temp;
    }
    else if (fPeriod.EqualTo("LHC18r_pass3")) {
      vector<int> temp = {297595, 297590, 297588, 297558, 297544, 297542, 297541, 297540, 297537, 297512, 297483, 297479, 297452, 297451, 297450, 297446, 297442, 297441, 297415, 297414, 297413, 297406, 297405, 297380, 297379, 297372, 297367, 297366, 297363, 297336, 297335, 297333, 297332, 297317, 297311, 297310, 297278, 297222, 297221, 297218, 297196, 297195, 297193, 297133, 297132, 297129, 297128, 297124, 297123, 297119, 297118, 297117, 297085, 297035, 297031, 296966, 296941, 296938, 296935, 296934, 296932, 296931, 296930, 296903, 296900, 296899, 296894, 296852, 296851, 296850, 296848, 296839, 296838, 296836, 296835, 296799, 296794, 296793, 296790, 296787, 296786, 296785, 296784, 296781, 296752, 296694, 296693, 296691, 296690};
      RunNumber = temp;
    }
    else if (fPeriod.EqualTo("LHC18q_pass3")) {
      vector<int> temp = {296623, 296622, 296621, 296619, 296618, 296616, 296615, 296594, 296553, 296552, 296551, 296550, 296548, 296547, 296516, 296512, 296511, 296510, 296509, 296472, 296433, 296424, 296423, 296420, 296419, 296415, 296414, 296383, 296381, 296380, 296379, 296378, 296377, 296376, 296375, 296312, 296309, 296304, 296303, 296280, 296279, 296273, 296270, 296269, 296247, 296246, 296244, 296243, 296242, 296241, 296240, 296198, 296197, 296196, 296195, 296194, 296192, 296191, 296143, 296142, 296135, 296134, 296133, 296132, 296123, 296074, 296066, 296065, 296063, 296062, 296060, 296016, 295942, 295941, 295937, 295936, 295913, 295910, 295909, 295861, 295860, 295859, 295856, 295855, 295854, 295853, 295831, 295829, 295826, 295825, 295822, 295819, 295818, 295816, 295791, 295788, 295786, 295763, 295762, 295759, 295758, 295755, 295754, 295725, 295723, 295721, 295719, 295718, 295717, 295714, 295712, 295676, 295675, 295673, 295668, 295667, 295666, 295615, 295612, 295611, 295610, 295589, 295588, 295586, 295585};
      RunNumber = temp;
    }
    else if (fPeriod.EqualTo("LHC17n")) {
      vector<int> temp = {280234, 280235};
      RunNumber = temp;
    }
    for (auto RunIter = RunNumber.begin(); RunIter != RunNumber.end(); RunIter++)
    {
      for(Int_t i=0;i<fTotFlags;i++) { //One less needed, because otherwise "Nominal" would be recorded twice
        //fEvNomFlag and fTrNomFlag are already in mask-format, so for those, need to fetch corresponding indeces
        Int_t defEvInd = i<fTotEvFlags?i:BitIndex(fEvNomFlag);
        Int_t defTrInd = (i-fTotEvFlags);
        if(defTrInd<0) defTrInd=BitIndex(fTrNomFlag); //if index < 0, then we're still in the event flags. Thus, select the nominal here
        else if(defTrInd==BitIndex(fTrNomFlag)) continue; //Otherwise, check if we are doing the nominal flag. If so, the skip this, b/c it's already done with nominal event selection
        fWeightList->Add(new AliGFWWeights());
        fWeights = (AliGFWWeights*)fWeightList->Last();
        fWeights->SetName(Form("w%s%s",to_string(*RunIter).c_str(),GetSystPF(defEvInd,defTrInd).Data()));
        fWeights->Init(!fIsMC,fIsMC); // AddData = !fIsMC; AddMC = fIsMC
      }
    }
  } else {
    //Setup the structure of FC:
    TObjArray *OAforPt=new TObjArray();
    //No gap:
    OAforPt->Add(new TNamed("MidV22","MidV22"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV22_pt_%i",i+1),"MidV22_pTDiff"));
    OAforPt->Add(new TNamed("MidV24","MidV24"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV24_pt_%i",i+1),"MidV24_pTDiff"));
    OAforPt->Add(new TNamed("MidV26","MidV26"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV26_pt_%i",i+1),"MidV26_pTDiff"));
    OAforPt->Add(new TNamed("MidV28","MidV28"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV28_pt_%i",i+1),"MidV28_pTDiff"));
    OAforPt->Add(new TNamed("MidV32","MidV32"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV32_pt_%i",i+1),"MidV32_pTDiff"));
    OAforPt->Add(new TNamed("MidV34","MidV34"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV34_pt_%i",i+1),"MidV34_pTDiff"));
    OAforPt->Add(new TNamed("MidV42","MidV42"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV42_pt_%i",i+1),"MidV42_pTDiff"));
    OAforPt->Add(new TNamed("MidV52","MidV52"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV52_pt_%i",i+1),"MidV52_pTDiff"));

    //2SENeg:
    OAforPt->Add(new TNamed("Mid2SENV22","Mid2SENV22"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SENV22_pt_%i",i+1),"Mid2SENV22_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SENV24","Mid2SENV24"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SENV24_pt_%i",i+1),"Mid2SENV24_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SENV26","Mid2SENV26"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SENV26_pt_%i",i+1),"Mid2SENV26_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SENV28","Mid2SENV28"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SENV28_pt_%i",i+1),"Mid2SENV28_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SENV32","Mid2SENV32"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SENV32_pt_%i",i+1),"Mid2SENV32_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SENV42","Mid2SENV42"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SENV42_pt_%i",i+1),"Mid2SENV42_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SENV52","Mid2SENV52"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SENV52_pt_%i",i+1),"Mid2SENV52_pTDiff"));

    //2SEPos:
    OAforPt->Add(new TNamed("Mid2SEPV22","Mid2SEPV22"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEPV22_pt_%i",i+1),"Mid2SEPV22_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEPV24","Mid2SEPV24"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEPV24_pt_%i",i+1),"Mid2SEPV24_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEPV26","Mid2SEPV26"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEPV26_pt_%i",i+1),"Mid2SEPV26_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEPV28","Mid2SEPV28"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEPV28_pt_%i",i+1),"Mid2SEPV28_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEPV32","Mid2SEPV32"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEPV32_pt_%i",i+1),"Mid2SEPV32_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEPV42","Mid2SEPV42"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEPV42_pt_%i",i+1),"Mid2SEPV42_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEPV52","Mid2SEPV52"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEPV52_pt_%i",i+1),"Mid2SEPV52_pTDiff"));

    //|eta|>0.5
    OAforPt->Add(new TNamed("MidGapNV22","MidGapNV22"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapNV22_pt_%i",i+1),"MidGapNV22_pTDiff"));
    OAforPt->Add(new TNamed("MidGapNV24","MidGapNV24"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapNV24_pt_%i",i+1),"MidGapNV24_pTDiff"));
    OAforPt->Add(new TNamed("MidGapNV26","MidGapNV26"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapNV26_pt_%i",i+1),"MidGapNV26_pTDiff"));
    OAforPt->Add(new TNamed("MidGapNV28","MidGapNV28"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapNV28_pt_%i",i+1),"MidGapNV28_pTDiff"));
    OAforPt->Add(new TNamed("MidGapNV32","MidGapNV32"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapNV32_pt_%i",i+1),"MidGapNV32_pTDiff"));
    OAforPt->Add(new TNamed("MidGapNV42","MidGapNV42"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapNV42_pt_%i",i+1),"MidGapNV42_pTDiff"));
    OAforPt->Add(new TNamed("MidGapNV52","MidGapNV52"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapNV52_pt_%i",i+1),"MidGapNV52_pTDiff"));

    //|eta|>0.5
    OAforPt->Add(new TNamed("MidGapPV22","MidGapPV22"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapPV22_pt_%i",i+1),"MidGapPV22_pTDiff"));
    OAforPt->Add(new TNamed("MidGapPV24","MidGapPV24"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapPV24_pt_%i",i+1),"MidGapPV24_pTDiff"));
    OAforPt->Add(new TNamed("MidGapPV26","MidGapPV26"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapPV26_pt_%i",i+1),"MidGapPV26_pTDiff"));
    OAforPt->Add(new TNamed("MidGapPV28","MidGapPV28"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapPV28_pt_%i",i+1),"MidGapPV28_pTDiff"));
    OAforPt->Add(new TNamed("MidGapPV32","MidGapPV32"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapPV32_pt_%i",i+1),"MidGapPV32_pTDiff"));
    OAforPt->Add(new TNamed("MidGapPV42","MidGapPV42"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapPV42_pt_%i",i+1),"MidGapPV42_pTDiff"));
    OAforPt->Add(new TNamed("MidGapPV52","MidGapPV52"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapPV52_pt_%i",i+1),"MidGapPV52_pTDiff"));

    //Multi bins:
    Double_t multibins[] = {0,5,10,20,30,40,50,60,70};
    fFC = new AliGFWFlowContainer();
    fFC->SetName(Form("FC%s",GetSystPF(BitIndex(fEvNomFlag),BitIndex(fTrNomFlag)).Data()));
    fFC->SetXAxis(fPtAxis);
    fFC->Initialize(OAforPt,8,multibins,10); //Statistics only required for nominal profiles, so do not create randomized profiles for systematics
    //Powers per harmonic:
    Int_t NoGap[] = {9,0,8,4,7,2,6,0,5};
    Int_t WithGap[] = {5,0,2,2,3,2,4,0,5};
    Int_t POIPowers[] = {2,0,2,2,2,2}; //POIs will always be with weight of 1 (or 0, for what matters)

    fGFW = new AliGFW();
    //Full regions
    fGFW->AddRegion("poiMid",6,POIPowers,-0.8,0.8,1+fPtAxis->GetNbins(),1);
    fGFW->AddRegion("refMid",9,NoGap,-0.8,0.8,1,2);
    //2 subevets:
    fGFW->AddRegion("poiSENeg",6,POIPowers,-0.8,0.,1+fPtAxis->GetNbins(),1);
    fGFW->AddRegion("refSENeg",9,WithGap,-0.8,0.,1,2);
    fGFW->AddRegion("poiSEPos",6,POIPowers,0.,0.8,1+fPtAxis->GetNbins(),1);
    fGFW->AddRegion("refSEPos",9,WithGap,0.,0.8,1,2);
    //With gap
    fGFW->AddRegion("poiGapNeg",6,POIPowers,-0.8,-0.5,1+fPtAxis->GetNbins(),1);
    fGFW->AddRegion("refGapNeg",9,WithGap,-0.8,-0.5,1,2);
    fGFW->AddRegion("poiGapPos",6,POIPowers,0.5,0.8,1+fPtAxis->GetNbins(),1);
    fGFW->AddRegion("refGapPos",9,WithGap,0.5,0.8,1,2);
    //Overlap:
    fGFW->AddRegion("olMid",9,NoGap,-0.8,0.8,1+fPtAxis->GetNbins(),4);
    fGFW->AddRegion("olSENeg",9,WithGap,-0.8,0.,1+fPtAxis->GetNbins(),4);
    fGFW->AddRegion("olSEPos",9,WithGap,0.,0.8,1+fPtAxis->GetNbins(),4);
    fGFW->AddRegion("olGapNeg",9,WithGap,-0.8,-0.5,1+fPtAxis->GetNbins(),4);
    fGFW->AddRegion("olGapPos",9,WithGap,0.5,0.8,1+fPtAxis->GetNbins(),4);

  };
  if(fProduceWeights) PostData(1,fWeightList);
  else PostData(1,fFC);
  fMultiDist = new TH1D(Form("Multiplicity_distribution%s",GetSystPF(BitIndex(fEvNomFlag),BitIndex(fTrNomFlag)).Data()),"Multiplicity distribution",100, 0, 100);
  PostData(2,fMultiDist);
  if(!fProduceWeights) {
    if(fIsTrain) {
      fCentMap = (TH1D*) GetInputData(1);
      if(!fCentMap) AliFatal("Could not fetch centrality map!\n");
    } else {
      fWeightList = (TList*) GetInputData(1);
      if(!fWeightList) { AliFatal("Could not retrieve weight list!\n"); return; };
    };
    CreateCorrConfigs();
  };
};


AliMCEvent *AliAnalysisTaskGFWFlow::FetchMCEvent(Double_t &impactParameter) {
  /*
  //Reverting to the old implementation. The new one fails on trains (but old one fails locally :( )
  if(!fIsTrain) { AliFatal("Snap, Jim! Ain't no train here!\n"); return 0; }
  AliMCEvent* ev = dynamic_cast<AliMCEvent*>(MCEvent());
  if(!ev) { AliFatal("MC event not found!"); return 0; }
  AliGenHepMCEventHeader *header = dynamic_cast<AliGenHepMCEventHeader*>(ev->GenEventHeader());
  if(!header) { AliFatal("MC event not generated!"); return 0; }
  impactParameter = header->impact_parameter();
  return ev;*/

  //The old implementation
  if(!fIsTrain) { AliFatal("Snap, Jim! Ain't no train here! :(\n"); return 0; }
  AliMCEvent* ev = dynamic_cast<AliMCEvent*>(MCEvent());
  if(!ev) { AliFatal("MC event not found!"); return 0; }
  if(fOverrideCentrality>=0) return ev;
  AliGenEventHeader *header = dynamic_cast<AliGenEventHeader*>(ev->GenEventHeader());
  if(!header) { AliFatal("MC event not generated!"); return 0; }
  AliCollisionGeometry* headerH;
  TString genName;
  TList *ltgen = (TList*)ev->GetCocktailList();
  if (ltgen) {
  for(auto&& listObject: *ltgen){
    genName = Form("%s",listObject->GetName());
    if (genName.Contains("Hijing")) {
      headerH = dynamic_cast<AliCollisionGeometry*>(listObject);
      break;
      }
    }
  }
  else
    headerH = dynamic_cast<AliCollisionGeometry*>(ev->GenEventHeader());
  if(headerH){
      impactParameter = headerH->ImpactParameter();
  }
  return ev;
}

void AliAnalysisTaskGFWFlow::UserExec(Option_t*) {
  if(fIsTrain) { //Separate sniplet for running on trains. If it's a train, then do it and be done with it.
    Double_t lImpactParameter = -1;
    AliMCEvent *fEv = FetchMCEvent(lImpactParameter);
    if(!fEv) return;
    Double_t l_Cent;
    if(fOverrideCentrality>0) l_Cent = fOverrideCentrality; else {
      if(lImpactParameter < 0) AliFatal("Impact parameter is negative!\n");
      l_Cent = GetCentFromIP(lImpactParameter);
    };
    fMultiDist->Fill(l_Cent);
    if(l_Cent>70 || l_Cent<5) return; //not considering anything below 5% or above 70%
    Int_t nTracks = fEv->GetNumberOfPrimaries();
    if(nTracks < 1) { return; }
    //GFW stuff:
    fGFW->Clear();
    for(Int_t i=0;i<nTracks;i++) {
      AliMCParticle* lPart = dynamic_cast<AliMCParticle*>(fEv->GetTrack(i));
      if(!lPart) { continue; };
      Double_t l_pT=lPart->Pt();
      Double_t l_phi=lPart->Phi();
      Double_t l_eta=lPart->Eta();
      Bool_t WithinPtPOI = (fPOIpTMin<l_pT) && (l_pT<fPOIpTMax); //within POI pT range
      Bool_t WithinPtRF  = (fRFpTMin <l_pT) && (l_pT<fRFpTMax);  //within RF pT range
      if(!WithinPtPOI && !WithinPtRF) continue; //if the track is not within any pT range, then continue
      Int_t l_pTInd = fPtAxis->FindBin(l_pT)-1;
      if(WithinPtPOI) fGFW->Fill(l_eta,l_pTInd,l_phi,1,1); //Fill POI (mask = 1). Weights are always 1
      if(WithinPtRF)  fGFW->Fill(l_eta,l_pTInd,l_phi,1,2); //Fit RF (mask = 2). Weights are always 1
      if(WithinPtRF && WithinPtPOI) fGFW->Fill(l_eta,l_pTInd,l_phi,1,4); //Filling overlap. Weights are always 1
    };
    Bool_t filled;
    for(Int_t l_ind=0; l_ind<corrconfigs.size(); l_ind++) {
      filled = FillFCs(corrconfigs.at(l_ind),l_Cent,0);//,DisableOL);
    };
    PostData(1,fFC);
    PostData(2,fMultiDist);
    return;
  }

  //This part is on MC production (for weights) or data (for ordinary calculations)
  //First, check if the event passes event cuts
  AliGFWFlags *lFlags = (AliGFWFlags*)fInputEvent->FindListObject("GFWFlags");
  if(!lFlags) { printf("GFWFlags not found!\n"); return; };
  UInt_t gEventFlag = lFlags->GetEventFlags();
  if(fProduceWeights) { if(!gEventFlag) return; }//if producing weights, any flag is good flag
  else if(!(gEventFlag&fEvNomFlag)) return; //otherwise, if not the selected event flag, then move on
  AliMultSelection *lMultSel = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  Double_t cent = lMultSel->GetMultiplicityPercentile(fCollisionsSystem==1?"V0A":"V0M");
  if(fCollisionsSystem==2) if(!CheckTriggerVsCentrality(cent)) return;
  // if(!lFlags->CheckEventFlag(1)) return;
  fMultiDist->Fill(cent);
  if(fIsMC) {
    fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
    if(!fMCEvent)
      return;
  };
  if(fCollisionsSystem==2) { //only for Pb-Pb, because of the PU
    if(cent<0) return; //Do not consider 0-5%
    if(cent>70) return; //Also, peripheral cutoff //devel comm
  };
  if(fBypassCalculations) return;
  //fFlowEvent->SetRunNumber(fAOD->GetRunNumber());
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  Double_t vz = fAOD->GetPrimaryVertex()->GetZ();
  TClonesArray *tca = 0;
  if(fProduceWeights) {
    if(fIsMC) {
      tca = (TClonesArray*)fInputEvent->FindListObject("mcparticles");
      /*for(Int_t i=0;i<tca->GetEntries();i++) {
      	AliVParticle *lPart;
      	lPart = (AliAODMCParticle*)tca->At(i);
      	if(!AcceptParticle(lPart)) continue;
      	Int_t partbit = GetParticleBit(lPart);
      	if(!partbit) continue; //if no particle bit is set, no need to continue
      	Int_t combinedbit = CombineBits(vtxb,partbit);
      	Int_t count=0;
      	for(Int_t i=0;i<fTotFlags;++i)
      	  if(fSelections[i]->NeedsExtraWeight()) {
      	    if(combinedbit&(1<<i))
      	      ((AliGFWWeights*)fWeightList->At(count))->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),cent,2);
      	    count++;
      	  };
      };*/
    };
    AliAODTrack *lTrack;
    Double_t pt;
    UInt_t gTrackFlags=0;

    auto it = find(RunNumber.begin(), RunNumber.end(), fAOD->GetRunNumber());
    Int_t Runindex = -1;
    if (it == RunNumber.end()) { Runindex = -1;}
    else { Runindex = it - RunNumber.begin();}

    if (Runindex > -1)
    {
      for(Int_t lTr=0;lTr<lFlags->GetNFiltered();lTr++) {
        gTrackFlags = lFlags->GetTrackFlag(lTr);
        Int_t trInd = lFlags->GetTrackIndex(lTr);
        lTrack = (AliAODTrack*)fAOD->GetTrack(trInd);
        pt=lTrack->Pt();
        if(pt<fRFpTMin || pt>fRFpTMax) continue; //NUAs just for RF
        //For event cuts, keep track cuts nominal:
        if(gTrackFlags&fTrNomFlag) for(Int_t jEv=0;jEv<fTotEvFlags;jEv++) if(gEventFlag&(1<<jEv))
           ((AliGFWWeights*)fWeightList->At(Runindex*(fTotFlags-1)+jEv))->Fill(lTrack->Phi(),lTrack->Eta(),vz,pt,cent,0); //Nominal one is calculated twice
        //For track cuts, keep event cuts nominal:
        if(gEventFlag&fEvNomFlag) for(Int_t jTr=1;jTr<fTotTrackFlags;jTr++) if(gTrackFlags&(1<<jTr)) //also, start track flag from 1, b/c 0 flag is nominal and already recorded with nominal ev cuts
          ((AliGFWWeights*)fWeightList->At(Runindex*(fTotFlags-1)+fTotEvFlags+jTr-1))->Fill(lTrack->Phi(),lTrack->Eta(),vz,pt,cent,0);
      };
    }
    PostData(1,fWeightList);
    PostData(2,fMultiDist);
    return;
  } else {
    fGFW->Clear();
    AliAODTrack *lTrack;
    UInt_t gTrackFlags=0;
    Double_t l_pT;
    //Event check performed already
    for(Int_t lTr=0;lTr<lFlags->GetNFiltered();lTr++) {
      gTrackFlags = lFlags->GetTrackFlag(lTr);
      if(!(gTrackFlags&fTrNomFlag)) continue; //Check if we want to accept the track
      Int_t trInd = lFlags->GetTrackIndex(lTr);
      lTrack = (AliAODTrack*)fAOD->GetTrack(trInd);
      l_pT=lTrack->Pt();
      Bool_t WithinPtPOI = (fPOIpTMin<l_pT) && (l_pT<fPOIpTMax); //within POI pT range
      Bool_t WithinPtRF  = (fRFpTMin <l_pT) && (l_pT<fRFpTMax);  //within RF pT range
      if(!WithinPtPOI && !WithinPtRF) continue; //if the track is not within any pT range, then continue
      if(!fWeights) printf("Weights do not exist!\n");
      Double_t nua = fWeights->GetNUA(lTrack->Phi(),lTrack->Eta(),vz);
      Double_t nue = 1; //Since doing pT-diff., we can set this to one for speed up.
    	if(WithinPtPOI) fGFW->Fill(lTrack->Eta(),fPtAxis->FindBin(l_pT)-1,lTrack->Phi(),nua*nue,1); //Fill POI (mask = 1)
      if(WithinPtRF)  fGFW->Fill(lTrack->Eta(),fPtAxis->FindBin(l_pT)-1,lTrack->Phi(),nua*nue,2); //Fit RF (mask = 2)
      if(WithinPtRF && WithinPtPOI) fGFW->Fill(lTrack->Eta(),fPtAxis->FindBin(l_pT)-1,lTrack->Phi(),nua*nue,4); //Filling overlap
    };
    TRandom rndm(0);
    Double_t rndmn=rndm.Rndm();
    Bool_t filled;
    for(Int_t l_ind=0; l_ind<corrconfigs.size(); l_ind++) {
      filled = FillFCs(corrconfigs.at(l_ind),cent,rndmn);//,DisableOL);
    };
    PostData(1,fFC);
    PostData(2,fMultiDist);
  };
};
void AliAnalysisTaskGFWFlow::Terminate(Option_t*) {
  // printf("\n********* Time: %f\n**********",mywatch.RealTime());
  // printf("Filling time: %f\n",mywatchFill.RealTime());
  // printf("Storing time: %f\n",mywatchStore.RealTime());
};
void AliAnalysisTaskGFWFlow::SetPtBins(Int_t nBins, Double_t *bins, Double_t RFpTMin, Double_t RFpTMax) {
  fPtAxis->Set(nBins, bins);
  //Set pT range given by the defined pT axis
  fPOIpTMin=bins[0];
  fPOIpTMax=bins[nBins];
  //if RF pT range is not defined, use the same as for POI
  fRFpTMin=(RFpTMin<0)?fPOIpTMin:RFpTMin;
  fRFpTMax=(RFpTMax<0)?fPOIpTMax:RFpTMax;
}

Bool_t AliAnalysisTaskGFWFlow::CheckTriggerVsCentrality(Double_t l_cent) {
  UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(!(fTriggerType&fSelMask)) { return kFALSE; }; //printf("Returning from the generic check\n");
  if(fSelMask&(fTriggerType&(AliVEvent::kINT7+AliVEvent::kMB))) {return kTRUE; }; //printf("Passed by MB trigger!\n");
  if((fSelMask&fTriggerType&AliVEvent::kCentral) && l_cent>10) {return kFALSE; }; //printf("Returnning from kCent case\n");
  if((fSelMask&fTriggerType&AliVEvent::kSemiCentral) && (l_cent<30 || l_cent>50)) {return kFALSE; }; //printf("Returning from kSC case\n");
  return kTRUE;
}
Bool_t AliAnalysisTaskGFWFlow::LoadWeights(Int_t runno) { //Cannot be used when running on the trains
  TString wName=Form("w%i%s",runno,GetSystPF(BitIndex(fEvNomFlag), BitIndex(fTrNomFlag)).Data());
  if(fWeightList) {
    fWeights = (AliGFWWeights*)fWeightList->FindObject("WeightList")->FindObject(wName.Data());
    if(!fWeights) {
      fWeightList->ls();
      AliFatal("Weights could not be found in the list!\n");
      return kFALSE;
    };
    fWeights->CreateNUA();
    fWeights->CreateNUE();
    return kTRUE;
  } else {
    AliFatal("Weight list (for some reason) not set!\n");
    return kFALSE;
  };
};
void AliAnalysisTaskGFWFlow::NotifyRun() {
  if(fIsTrain) return;
  Int_t runno = fInputEvent->GetRunNumber();
  if(!fProduceWeights)
    LoadWeights(runno);
}
Bool_t AliAnalysisTaskGFWFlow::FillFCs(AliGFW::CorrConfig corconf, Double_t cent, Double_t rndmn, Bool_t DisableOverlap) {
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).real();
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).real()/dnx;
    if(TMath::Abs(val)<1)
      fFC->FillProfile(corconf.Head.c_str(),cent,val,dnx,rndmn);
    return kTRUE;
  };
  for(Int_t i=1;i<=fPtAxis->GetNbins();i++) {
    dnx = fGFW->Calculate(corconf,i-1,kTRUE).real();
    if(dnx==0) continue;
    val = fGFW->Calculate(corconf,i-1,kFALSE).real()/dnx;
    if(TMath::Abs(val)<1)
      fFC->FillProfile(Form("%s_pt_%i",corconf.Head.c_str(),i),cent,val,dnx,rndmn);
  };
  return kTRUE;
};
void AliAnalysisTaskGFWFlow::CreateCorrConfigs() {
  corrconfigs.push_back(GetConf("MidV22","refMid {2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidV22","poiMid refMid | olMid {2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidV24","refMid {2 2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidV24","poiMid refMid | olMid {2 2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidV26","refMid {2 2 2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidV26","poiMid refMid | olMid {2 2 2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidV28","refMid {2 2 2 2 -2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidV28","poiMid refMid | olMid {2 2 2 2 -2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidV32","refMid {3 -3}", kFALSE));
  corrconfigs.push_back(GetConf("MidV32","poiMid refMid | olMid {3 -3}", kTRUE));
  corrconfigs.push_back(GetConf("MidV34","refMid {3 3 -3 -3}", kFALSE));
  corrconfigs.push_back(GetConf("MidV34","poiMid refMid | olMid {3 3 -3 -3}", kTRUE));
  corrconfigs.push_back(GetConf("MidV42","refMid {4 -4}", kFALSE));
  corrconfigs.push_back(GetConf("MidV42","poiMid refMid | olMid {4 -4}", kTRUE));
  corrconfigs.push_back(GetConf("MidV52","poiMid refMid | olMid {5 -5}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SENV22","refSENeg {2} refSEPos {-2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SENV24","refSENeg {2 2} refSEPos {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SENV26","refSENeg {2 2 2} refSEPos {-2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SENV28","refSENeg {2 2 2 2} refSEPos {-2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SENV22","poiSENeg refSENeg | olSENeg {2} refSEPos {-2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SENV24","poiSENeg refSENeg | olSENeg {2 2} refSEPos {-2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SENV26","poiSENeg refSENeg | olSENeg {2 2 2} refSEPos {-2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SENV28","poiSENeg refSENeg | olSENeg {2 2 2 2} refSEPos {-2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SEPV22","refSEPos {2} refSENeg {-2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SEPV24","refSEPos {2 2} refSENeg {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SEPV26","refSEPos {2 2 2} refSENeg {-2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SEPV28","refSEPos {2 2 2 2} refSENeg {-2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SEPV22","poiSEPos refSEPos | olSEPos {2} refSENeg {-2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SEPV24","poiSEPos refSEPos | olSEPos {2 2} refSENeg {-2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SEPV26","poiSEPos refSEPos | olSEPos {2 2 2} refSENeg {-2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SEPV28","poiSEPos refSEPos | olSEPos {2 2 2 2} refSENeg {-2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SENV32","refSENeg {3} refSEPos {-3}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SENV32","poiSENeg refSENeg | olSENeg {3} refSEPos {-3}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SEPV32","refSEPos {3} refSENeg {-3}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SEPV32","poiSEPos refSEPos | olSEPos {3} refSENeg {-3}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SENV42","refSENeg {4} refSEPos {-4}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SENV42","poiSENeg refSENeg | olSENeg {4} refSEPos {-4}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SEPV42","refSEPos {4} refSENeg {-4}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SEPV42","poiSEPos refSEPos | olSEPos {4} refSENeg {-4}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SENV52","refSENeg {5} refSEPos {-5}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SENV52","poiSENeg refSENeg | olSENeg {5} refSEPos {-5}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SEPV52","refSEPos {5} refSENeg {-5}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SEPV52","poiSEPos refSEPos | olSEPos {5} refSENeg {-5}", kTRUE));

  corrconfigs.push_back(GetConf("MidGapNV22","refGapNeg {2} refGapPos {-2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapNV24","refGapNeg {2 2} refGapPos {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapNV26","refGapNeg {2 2 2} refGapPos {-2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapNV28","refGapNeg {2 2 2 2} refGapPos {-2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapNV22","poiGapNeg refGapNeg | olGapNeg {2} refGapPos {-2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapNV24","poiGapNeg refGapNeg | olGapNeg {2 2} refGapPos {-2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapNV26","poiGapNeg refGapNeg | olGapNeg {2 2 2} refGapPos {-2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapNV28","poiGapNeg refGapNeg | olGapNeg {2 2 2 2} refGapPos {-2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapPV22","refGapPos {2} refGapNeg {-2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapPV24","refGapPos {2 2} refGapNeg {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapPV26","refGapPos {2 2 2} refGapNeg {-2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapPV28","refGapPos {2 2 2 2} refGapNeg {-2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapPV22","poiGapPos refGapPos | olGapPos {2} refGapNeg {-2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapPV24","poiGapPos refGapPos | olGapPos {2 2} refGapNeg {-2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapPV26","poiGapPos refGapPos | olGapPos {2 2 2} refGapNeg {-2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapPV28","poiGapPos refGapPos | olGapPos {2 2 2 2} refGapNeg {-2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapNV32","refGapNeg {3} refGapPos {-3}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapNV32","poiGapNeg refGapNeg | olGapNeg {3} refGapPos {-3}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapPV32","refGapPos {3} refGapNeg {-3}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapPV32","poiGapPos refGapPos | olGapPos {3} refGapNeg {-3}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapNV42","refGapNeg {4} refGapPos {-4}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapNV42","poiGapNeg refGapNeg | olGapNeg {4} refGapPos {-4}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapPV42","refGapPos {4} refGapNeg {-4}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapPV42","poiGapPos refGapPos | olGapPos {4} refGapNeg {-4}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapNV52","refGapNeg {5} refGapPos {-5}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapNV52","poiGapNeg refGapNeg | olGapNeg {5} refGapPos {-5}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapPV52","refGapPos {5} refGapNeg {-5}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapPV52","poiGapPos refGapPos | olGapPos {5} refGapNeg {-5}", kTRUE));
}
void AliAnalysisTaskGFWFlow::SetupFlagsByIndex(Int_t ind) {
  SetupFlagsByIndex(ind,fEvNomFlag,fTrNomFlag);
}
void AliAnalysisTaskGFWFlow::SetupFlagsByIndex(const Int_t &ind, UInt_t &l_EvFlag, UInt_t &l_TrFlag) {
  l_EvFlag=1<<kNominal;
  l_TrFlag=1<<kFB96;
  switch(ind) {
    default: // also 0
      break;
    //Event flags:
    case 1:
      l_EvFlag = 1<<kVtx9;
      break;
    case 2:
      l_EvFlag = 1<<kVtx7;
      break;
    case 3:
      l_EvFlag = 1<<kVtx5;
      break;
    //Track flags:
    case 4:
      l_TrFlag = 1<<kFB768;
      break;
    case 5:
      l_TrFlag = 1<<kDCAz10;
      break;
    case 6:
      l_TrFlag = 1<<kDCAz05;
      break;
    case 7:
      l_TrFlag = 1<<kDCA4Sigma;
      break;
    case 8:
      l_TrFlag = 1<<kDCA10Sigma;
      break;
    case 9:
      l_TrFlag = 1<<kChiSq2;
      break;
    case 10:
      l_TrFlag = 1<<kChiSq3;
      break;
    case 11:
      l_TrFlag = 1<<kNTPC80;
      break;
    case 12:
      l_TrFlag = 1<<kNTPC90;
      break;
    case 13:
      l_TrFlag = 1<<kNTPC100;
      break;
    case 14:
      l_TrFlag = 1<<kFB768Tuned;
      break;
    case 15:
      l_TrFlag = 1<<kFB96Tuned;
      break;
    case 16:
      l_TrFlag = 1<<kFB768DCAz;
      break;
    case 17:
      l_TrFlag = 1<<kFB768DCAxyLow;
      break;
    case 18:
      l_TrFlag = 1<<kFB768DCAxyHigh;
      break;
    case 19:
      l_TrFlag = 1<<kFB768ChiSq2;
      break;
    case 20:
      l_TrFlag = 1<<kFB768ChiSq3;
      break;
    case 21:
      l_TrFlag = 1<<kFB768nTPC;
      break;
    case 22:
      l_TrFlag = 1<<kFB96MergedDCA;
      break;
    case 23: //This will now be used for rebinned NUA test
      l_TrFlag = 1<<kChiSq25;
      break;
    case 24:
      l_TrFlag = 1<<kRebinnedNUA;
      break;
  }
}
