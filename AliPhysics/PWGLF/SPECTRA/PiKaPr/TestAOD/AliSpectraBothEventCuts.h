#ifndef ALISPECTRABOTHEVENTCUTS_H
#define ALISPECTRABOTHEVENTCUTS_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliSpectraBothEventCuts
//
//
//
//
// Authors: Michele Floris, CERN, Philip Versteeg, UU, Redmer Bertens, UU
//-------------------------------------------------------------------------

class AliVEvent;
class AliSpectraBothTrackCuts;
//class AliSpectraBothHistoManager;

#include "TH1I.h"
#include "TNamed.h"
#include "AliSpectraBothTrackCuts.h"
#include "AliAnalysisUtils.h"	
#include "AliPPVsMultUtils.h"	
#include "AliMCEvent.h"
#include "AliVEvent.h"


class AliSpectraBothEventCuts : public TNamed
{
 public:
  enum {  kProcessedEvents = 0, // all events 
kPhysSelEvents, // phys selected and not pile-up
kAcceptedEvents, // events passing all cuts  
kGoodVtx,// events with good veretx and   
kVtxCentral, 
kVtxNoEvent, 
kQVector,
kTPCasPV,
kZeroCont,
kNVtxCuts};
enum {kDoNotCheckforSDD=0,kwithSDD,kwithoutSDD};	

  // Constructors
 AliSpectraBothEventCuts() : TNamed(), fAOD(0),fAODEvent(AliSpectraBothTrackCuts::kAODobject),fTrackBits(0), fIsMC(0), fCentEstimator(""), fUseCentPatchAOD049(0),fUseSDDPatchforLHC11a(kDoNotCheckforSDD),fTriggerSettings(AliVEvent::kMB),fTrackCuts(0),
fIsSelected(0), fCentralityCutMin(0), fCentralityCutMax(0), fQVectorCutMin(0), fQVectorCutMax(0), fVertexCutMin(0),fVertexCutMax(0), fMultiplicityCutMin(0), fMultiplicityCutMax(0), fMaxChi2perNDFforVertex(0),fMinRun(0),fMaxRun(0),fetarangeofmultiplicitycut(0.0),fUseAliPPVsMultUtils(false),fNMCProcessType(-1),fEventMCProcessType(0),fEventMCProcessTypeIncluded(0),fchecktypeofveretxbytitle(kTRUE),fvertexselection(-1),fDotheeventcutsinmultselection(kFALSE),
fDotheBGRejection(kTRUE),fDothePileUpRejection(kTRUE),
fHistoCuts(0),fHistoVtxBefSel(0),fHistoVtxAftSel(0),fHistoEtaBefSel(0),fHistoEtaAftSel(0),fHistoNChAftSel(0),fHistoQVector(0),fHistoEP(0), fHistoVtxAftSelwithoutZvertexCut(0),fHistoVtxalltriggerEventswithMCz(0),fHistoVtxAftSelwithoutZvertexCutusingMCz(0),fHistoRunNumbers(0),fHistoCentrality(0),fHistoMultiplicty(0),fAnalysisUtils(0),fAliPPVsMultUtils(0)
{

}
  AliSpectraBothEventCuts(const char *name);
  virtual ~AliSpectraBothEventCuts();// {}
  
  void SetIsMC(Bool_t isMC = kFALSE)    {fIsMC = isMC; };
  Bool_t GetIsMC()           const           { return fIsMC;};
  void SetCentEstimator(TString cent = "V0M")    {fCentEstimator = cent; };
  TString GetCentFromV0()           const           { return fCentEstimator;};
  
  void SetUseCentPatchAOD049(Bool_t useCentPatchAOD049 = kFALSE)    {fUseCentPatchAOD049 = useCentPatchAOD049; };
  Bool_t GetUseCentPatchAOD049()           const           { return fUseCentPatchAOD049;};
  void 	SetUseSDDPatchforLHC11a(Int_t useSDDPatchforLHC11a) {fUseSDDPatchforLHC11a=useSDDPatchforLHC11a;} ;  
  Int_t GetUseSDDPatchforLHC11a() {return fUseSDDPatchforLHC11a;};   

   void   SetTriggerSettings(UInt_t triggerSettings = AliVEvent::kMB) {fTriggerSettings = triggerSettings;};
   UInt_t   GetTriggerSettings() {return fTriggerSettings;};



  // Methods
  Bool_t IsSelected(AliVEvent * aod,AliSpectraBothTrackCuts     *trackcuts, Bool_t isMC=kFALSE, Double_t mcZ=-100.0, TH1F* managerhisteventcuts=0x0);
  Bool_t CheckVtx();
  Bool_t CheckCentralityCut();
  Bool_t CheckMultiplicityCut();
  Bool_t CheckQVectorCut();
  Bool_t CheckVtxChi2perNDF();
  Bool_t CheckMCProcessType(AliMCEvent* mcevent);	
  void  SetCentralityCutMin(Float_t cut)  { fCentralityCutMin = cut; }
  void  SetCentralityCutMax(Float_t cut)  { fCentralityCutMax = cut; }
  //void  SetQVectorPosCut(Float_t min,Float_t max)  { fQVectorPosCutMin = min; fQVectorPosCutMax = max; }
  //void  SetQVectorNegCut(Float_t min,Float_t max)  { fQVectorNegCutMin = min; fQVectorNegCutMax = max; }
  void  SetQVectorCut(Float_t min,Float_t max)  { fQVectorCutMin = min; fQVectorCutMax = max; }
  void  SetVertexCut(Float_t min,Float_t max)  { fVertexCutMin = min; fVertexCutMax = max; }
  void  SetMultiplicityCut(Int_t min,Int_t max)  { fMultiplicityCutMin = min; fMultiplicityCutMax = max; }
  void SetTrackBits(UInt_t TrackBits) {fTrackBits=TrackBits;}
  void SetMaxChi2perNDFforVertex(Float_t cut) { fMaxChi2perNDFforVertex=cut;}
  void SetRunNumberRange(Int_t min,Int_t max);	
  void SetEtaRangeforMultiplictyCut(Float_t eta) {fetarangeofmultiplicitycut=eta;}
  void  SetAnalysisUtils(AliAnalysisUtils* inAnalysisUtils){fAnalysisUtils=inAnalysisUtils;}
  void SetUseAliPPVsMultUtils(Bool_t flag){fUseAliPPVsMultUtils=flag;} 
  void  SetEventMCProcessTypeIncluded(Bool_t flag) {fEventMCProcessTypeIncluded=flag;}
  void SetNMCProcessType(Int_t flag); 
   void AddMCProcessType(Int_t type,Int_t index);
  void SetCheckTypeofVeretx(Bool_t flag) {fchecktypeofveretxbytitle=flag;}
 void SetVertexSelection(Int_t flag){fvertexselection=flag;}  
void SetDotheeventcutsinmultselection(Bool_t flag) {fDotheeventcutsinmultselection=flag;}
void SetDotheBGRejection(Bool_t flag) {fDotheBGRejection=flag;}
void SetDothePileUpRejection(Bool_t flag) {fDothePileUpRejection=flag;}

  UInt_t GetTrackType()  const    { return fTrackBits;}
  TH1I * GetHistoCuts()         {  return fHistoCuts; }
  TH1F * GetHistoVtxBefSel()         {  return fHistoVtxBefSel; }
  TH1F * GetHistoVtxAftSel()         {  return fHistoVtxAftSel; }
  TH1F * GetHistoVtxAftSelwithoutZvertexCut()         {  return fHistoVtxAftSelwithoutZvertexCut; }
  TH1F * GetHistoVtxGenerated()         {  return fHistoVtxalltriggerEventswithMCz; }
  TH1F * GetHistoVtxAftSelwithoutZvertexCutusingMCz()         {  return fHistoVtxAftSelwithoutZvertexCutusingMCz; }

  TH1F * GetHistoEtaBefSel()         {  return fHistoEtaBefSel; }
  TH1F * GetHistoEtaAftSel()         {  return fHistoEtaAftSel; }
  TH1F * GetHistoNChAftSel()         {  return fHistoNChAftSel; }
  /* TH1F * GetHistoQVectorPos()         {  return fHistoQVectorPos; } */
  /* TH1F * GetHistoQVectorNeg()         {  return fHistoQVectorNeg; } */
  TH1F * GetHistoQVector()         {  return fHistoQVector; }
  TH1F * GetHistoEP()         {  return fHistoEP; }
  TH1F * GetHistoRunNumbers()         {  return fHistoRunNumbers; }
  TH2F * GetHistoCentrality()	      {	return fHistoCentrality;}
  TH2F * GetHistoMultiplicty()         { return fHistoMultiplicty; }
  Float_t  GetCentralityMin()  const {  return fCentralityCutMin; }
  Float_t  GetCentralityMax()  const {  return fCentralityCutMax; }
  //Float_t  GetQVectorPosCutMin()  const {  return fQVectorPosCutMin; }
  //Float_t  GetQVectorPosCutMax()  const {  return fQVectorPosCutMax; }
  //Float_t  GetQVectorNegCutMin()  const {  return fQVectorNegCutMin; }
  //Float_t  GetQVectorNegCutMax()  const {  return fQVectorNegCutMax; }
  Float_t  GetQVectorCutMin()  const {  return fQVectorCutMin; }
  Float_t  GetQVectorCutMax()  const {  return fQVectorCutMax; }
  Float_t  GetVertexCutMin()  const {  return fVertexCutMin; }
  Float_t  GetVertexCutMax()  const {  return fVertexCutMax; }
  Float_t  GetMultiplicityCutMin()  const {  return fMultiplicityCutMin; }
  Float_t  GetMultiplicityCutMax()  const {  return fMultiplicityCutMax; }
  Float_t GetMaxChi2perNDFforVertex() const {return fMaxChi2perNDFforVertex;}
   AliAnalysisUtils* GetAnalysisUtils() const {return fAnalysisUtils; }
Bool_t GetUseAliPPVsMultUtils() const {return fUseAliPPVsMultUtils;} 
 Int_t  GetNMCProcessType() const { return fNMCProcessType;}
 Bool_t  GetEventMCProcessTypeIncluded() const {return fEventMCProcessTypeIncluded;}
 Bool_t  GetCheckTypeofVeretx() const {return fchecktypeofveretxbytitle;}
 Int_t  GetVertexSelection()const {return fvertexselection;} 
 Bool_t GetDotheeventcutsinmultselection() const {return fDotheeventcutsinmultselection;}


  void InitHisto();	
  void   PrintCuts();
  Double_t ApplyCentralityPatchAOD049();

  Float_t  NumberOfEvents()     { return fHistoCuts->GetBinContent(kAcceptedEvents+1); }
  Float_t  NumberOfProcessedEvents()     { return fHistoCuts->GetBinContent(kProcessedEvents+1); }
  Float_t  NumberOfPhysSelEvents()     { return fHistoCuts->GetBinContent(kPhysSelEvents+1); }
  Bool_t CheckifESDorAODEvent();
  Bool_t CheckSDDPatchforLHC11a();		

  Long64_t Merge(TCollection* list);


 private:
  
  AliVEvent     *fAOD;              //! AOD event
  Int_t 		 fAODEvent; // flag what type event is conected use values form AliSpeatraAODTrackCuts	
  UInt_t           fTrackBits;       // Type of track to be used in the Qvector calculation
  Bool_t          fIsMC;// true if processing MC
  TString          fCentEstimator;// name of centrality estimator
  Bool_t          fUseCentPatchAOD049;// Patch for centrality selection on AOD049
  Int_t          fUseSDDPatchforLHC11a; // if true will check for ALLNOTRD  in fired trigger class 
  UInt_t fTriggerSettings;     // triger configuration 
  AliSpectraBothTrackCuts     *fTrackCuts;             //! track cuts
  Bool_t          fIsSelected;        // True if cuts are selected
  Float_t         fCentralityCutMin;     // minimum centrality percentile
  Float_t         fCentralityCutMax;     // maximum centrality percentile
  /* Float_t         fQVectorPosCutMin;     // minimum qvecPos */
  /* Float_t         fQVectorPosCutMax;     // maximum qvecPos */
  /* Float_t         fQVectorNegCutMin;     // minimum qvecNeg */
  /* Float_t         fQVectorNegCutMax;     // maximum qvecNeg */
  Float_t         fQVectorCutMin;     // minimum qvecPos
  Float_t         fQVectorCutMax;     // maximum qvecPos
  Float_t         fVertexCutMin;     // minimum vertex position
  Float_t         fVertexCutMax;     // maximum vertex position
  Int_t         fMultiplicityCutMin;     // minimum multiplicity position
  Int_t         fMultiplicityCutMax;     // maximum multiplicity position
  Float_t	  fMaxChi2perNDFforVertex; // maximum value of Chi2perNDF of vertex
  Int_t           fMinRun;                //minmum run number 			 
  Int_t 	  fMaxRun;		  //maximum run number 	 
  Float_t         fetarangeofmultiplicitycut; // eta range fot multipilicty cut 
  Bool_t 	  fUseAliPPVsMultUtils;   // use  AliPPVsMultUtils for centrailty 
  Int_t           fNMCProcessType;  // to include or exlude 
  Int_t 	  *fEventMCProcessType;   //[fNMCProcessType] process typ cut 				  
  Bool_t 	  fEventMCProcessTypeIncluded; // if false those process are excluded if true they are included  
  Bool_t          fchecktypeofveretxbytitle; // if 1 type of veretx check using the title
  Int_t  	   fvertexselection; // -1= no check 0=not track not SPD 1= only track 2=only spd 3= track + SPD
  Bool_t 	   fDotheeventcutsinmultselection; // if true it assume that pile-up/veretx/inel cut are done during the centrality multiplicty selection by AliPPVsMultUtils or loader task 
  Bool_t 	   fDotheBGRejection;              // if 1 and the fDotheeventcutsinmultselection 0 BG is reject using AnalysisUtils
  Bool_t 	   fDothePileUpRejection;        // if 1 the  fDotheeventcutsinmultselection 0  PileUp is reject using AnalysisUtils 	 	  					 
  TH1I            *fHistoCuts;        // Cuts statistics
  TH1F            *fHistoVtxBefSel;        // Vtx distr before event selection 	
  TH1F            *fHistoVtxAftSel;        // Vtx distr after event selection
  TH1F            *fHistoEtaBefSel;        // Eta distr before event selection
  TH1F            *fHistoEtaAftSel;        // Eta distr after event selection
  TH1F            *fHistoNChAftSel;        // NCh distr after event selection
  //TH1F            *fHistoQVectorPos;        // QVectorPos
  //TH1F            *fHistoQVectorNeg;        // QVectorNeg
  TH1F            *fHistoQVector;        // QVector
  TH1F            *fHistoEP;        // EP
  TH1F 		  *fHistoVtxAftSelwithoutZvertexCut;        // Vtx distr after event selection but without z vertex cut
  TH1F 		  *fHistoVtxalltriggerEventswithMCz;        // MC z vertex ditribution of generated events
  TH1F 		  *fHistoVtxAftSelwithoutZvertexCutusingMCz;        // Vtx distr after event selection but without z vertex cut
  TH1F            *fHistoRunNumbers;   // histogram of numbers of events per run
  TH2F 		  *fHistoCentrality; // centrality histogram
  TH2F 		  *fHistoMultiplicty; // multiplicty histogram
  AliAnalysisUtils *fAnalysisUtils;// Analysis Utils which have pile-up cut
  AliPPVsMultUtils  *fAliPPVsMultUtils; //AliPPVsMultUtils class for centrailty

  AliSpectraBothEventCuts(const AliSpectraBothEventCuts&);
  AliSpectraBothEventCuts& operator=(const AliSpectraBothEventCuts&);
  
  
  ClassDef(AliSpectraBothEventCuts, 16);
  
};
#endif

