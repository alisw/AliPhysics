#include <Riostream.h>
#include <TMath.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TSpline.h>
#include <TRandom3.h>

#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliTHn.h"
#include "AliAnalysisTaskTriggeredBF.h"

#include "AliPidBFBase.h"
using std::cout;
using std::endl;
using std::cerr;

//-------------------------------------------------------------------------
// Base class for filling the mixed and same events histogram for BF (++, +- , -- ,-+)
// Noor Alam : vecc, Kolkata, noor1989phyalam@gmail.com
// Thanks to Panos : 
//-------------------------------------------------------------------------


ClassImp(AliPidBFBase)

//____________________________________________________________________//
AliPidBFBase::AliPidBFBase() :
  TObject(), 
  fAnalysisLevel("ESD"),
  fAnalyzedEvents(0) ,
  fCentralityId(0) ,
  fHistP(0),
  fHistN(0),
  fHistPN(0),
  fHistNP(0),
  fHistPP(0),
  fHistNN(0),
  fHistHBTbefore(0),
  fHistHBTafter(0),
  fHistPhiStarHBTbefore(0),
  fHistPhiStarHBTafter(0),
  fHistConversionbefore(0),
  fHistConversionafter(0),
  fHistPsiMinusPhi(0),
  fHistResonancesBefore(0),
  fHistResonancesRho(0),
  fHistResonancesK0(0),
  fHistResonancesLambda(0),
  fHistQbefore(0),
  fHistQafter(0),
  fPsiInterval(15.),
  fDeltaEtaMax(2.0),
  fMomentumOrdering(kTRUE),
  fResonancesCut(kFALSE),
  fHBTCut(kFALSE),
  fHBTCutValue(0.02),
  fConversionCut(kFALSE),
  fInvMassCutConversion(0.04),
  fQCut(kFALSE),
  fDeltaPtMin(0.0),
  fVertexBinning(kFALSE),
  fCustomBinning(""),
  fBinningString(""),
  fEventClass("EventPlane"),
  kTrackVariablesSingle(2),
  kTrackVariablesPair(5){
  // Default constructor
}

//____________________________________________________________________//
AliPidBFBase::AliPidBFBase(const AliPidBFBase& balance):
  TObject(balance), 
  fAnalysisLevel(balance.fAnalysisLevel),
  fAnalyzedEvents(balance.fAnalyzedEvents), 
  fCentralityId(balance.fCentralityId),
  fHistP(balance.fHistP),
  fHistN(balance.fHistN),
  fHistPN(balance.fHistPN),
  fHistNP(balance.fHistNP),
  fHistPP(balance.fHistPP),
  fHistNN(balance.fHistNN),
  fHistHBTbefore(balance.fHistHBTbefore),
  fHistHBTafter(balance.fHistHBTafter),
  fHistPhiStarHBTbefore(balance.fHistPhiStarHBTbefore),
  fHistPhiStarHBTafter(balance.fHistPhiStarHBTafter),
  fHistConversionbefore(balance.fHistConversionbefore),
  fHistConversionafter(balance.fHistConversionafter),
  fHistPsiMinusPhi(balance.fHistPsiMinusPhi),
  fHistResonancesBefore(balance.fHistResonancesBefore),
  fHistResonancesRho(balance.fHistResonancesRho),
  fHistResonancesK0(balance.fHistResonancesK0),
  fHistResonancesLambda(balance.fHistResonancesLambda),
  fHistQbefore(balance.fHistQbefore),
  fHistQafter(balance.fHistQafter),
  fPsiInterval(balance.fPsiInterval),
  fDeltaEtaMax(balance.fDeltaEtaMax),
  fMomentumOrdering(balance.fMomentumOrdering),
  fResonancesCut(balance.fResonancesCut),
  fHBTCut(balance.fHBTCut),
  fHBTCutValue(balance.fHBTCutValue),
  fConversionCut(balance.fConversionCut),
  fInvMassCutConversion(balance.fInvMassCutConversion),
  fQCut(balance.fQCut),
  fDeltaPtMin(balance.fDeltaPtMin),
  fVertexBinning(balance.fVertexBinning),
  fCustomBinning(balance.fCustomBinning),
  fBinningString(balance.fBinningString),
  fEventClass("EventPlane"),
  kTrackVariablesSingle(balance.kTrackVariablesSingle),
  kTrackVariablesPair(balance.kTrackVariablesPair){
  //copy constructor
}

//____________________________________________________________________//
AliPidBFBase::~AliPidBFBase() {
  // Destructor
  delete fHistP;
  delete fHistN;
  delete fHistPN;
  delete fHistNP;
  delete fHistPP;
  delete fHistNN;

  delete fHistHBTbefore;
  delete fHistHBTafter;
  delete fHistPhiStarHBTbefore;
  delete fHistPhiStarHBTafter;
  delete fHistConversionbefore;
  delete fHistConversionafter;
  delete fHistPsiMinusPhi;
  delete fHistResonancesBefore;
  delete fHistResonancesRho;
  delete fHistResonancesK0;
  delete fHistResonancesLambda;
  delete fHistQbefore;
  delete fHistQafter;
    
}

//____________________________________________________________________//
void AliPidBFBase::InitHistograms() {
  // single particle histograms

  // global switch disabling the reference 
  // (to avoid "Replacing existing TH1" if several wagons are created in train)
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  Int_t anaSteps   = 1;       // analysis steps
  Int_t iBinSingle[kTrackVariablesSingle];        // binning for track variables
  Double_t* dBinsSingle[kTrackVariablesSingle];   // bins for track variables  
  TString* axisTitleSingle = new TString[kTrackVariablesSingle]; // axis titles for track variables
  
  // two particle histograms
  Int_t iBinPair[kTrackVariablesPair];         // binning for track variables
  Double_t* dBinsPair[kTrackVariablesPair];    // bins for track variables  
  TString* axisTitlePair = new TString[kTrackVariablesPair];  // axis titles for track variables



  // =========================================================
  // The default string (from older versions of AliPidBFBase)
  // =========================================================
  TString defaultBinningStr;
  defaultBinningStr = "multiplicity: 0,10,20,30,40,50,60,70,80,100,100000\n"                // Multiplicity Bins
    "centrality: 0.,5.,10.,20.,30.,40.,50.,60.,70.,80.\n"                                   // Centrality Bins
    "centralityVertex: 0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.\n" // Centrality Bins (Vertex Binning)
    "eventPlane: -0.5,0.5,1.5,2.5,3.5\n"                                                    // Event Plane Bins (Psi: -0.5->0.5 (in plane), 0.5->1.5 (intermediate), 1.5->2.5 (out of plane), 2.5->3.5 (rest))
    "deltaEta: -1.6, -1.56, -1.52, -1.48, -1.44, -1.4, -1.36, -1.32, -1.28, -1.24, -1.2, -1.16, -1.12, -1.08, -1.04, -1, -0.96, -0.92, -0.88, -0.84, -0.8, -0.76, -0.72, -0.68, -0.64, -0.6, -0.56, -0.52, -0.48, -0.44, -0.4, -0.36, -0.32, -0.28, -0.24, -0.2, -0.16, -0.12, -0.08, -0.04, 0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1, 1.04, 1.08, 1.12, 1.16, 1.2, 1.24, 1.28, 1.32, 1.36, 1.4, 1.44, 1.48, 1.52, 1.56, 1.6\n" // Delta Eta Bins
    "deltaEtaVertex: -1.6, -1.52, -1.44, -1.36, -1.28, -1.2, -1.12, -1.04, -0.96, -0.88, -0.8, -0.72, -0.64, -0.56, -0.48, -0.4, -0.32, -0.24,-0.16, -0.08, 0, 0.08, 0.16, 0.24, 0.32, 0.4, 0.48, 0.56, 0.64, 0.72, 0.8, 0.88, 0.96, 1.04, 1.12, 1.2, 1.28, 1.36, 1.44, 1.52, 1.6\n" // Delta Eta Bins (Vertex Binning)
    "deltaPhi: -1.5708, -1.48353, -1.39626, -1.309, -1.22173, -1.13446, -1.0472, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.0872665, 0, 0.0872665, 0.174533, 0.261799, 0.349066, 0.436332, 0.523599, 0.610865, 0.698132, 0.785398, 0.872665, 0.959931, 1.0472, 1.13446, 1.22173, 1.309, 1.39626, 1.48353, 1.5708, 1.65806, 1.74533, 1.8326, 1.91986, 2.00713, 2.0944, 2.18166, 2.26893, 2.35619, 2.44346, 2.53073, 2.61799, 2.70526, 2.79253, 2.87979, 2.96706, 3.05433, 3.14159, 3.22886, 3.31613, 3.40339, 3.49066, 3.57792, 3.66519, 3.75246, 3.83972, 3.92699, 4.01426, 4.10152, 4.18879, 4.27606, 4.36332, 4.45059, 4.53786, 4.62512, 4.71239\n" // Delta Phi Bins
    "pT: 0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0\n"             // pT Bins
    "pTVertex: 0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0\n"                                              // pT Bins (Vertex Binning)
    "vertex: -10., 10.\n"                                                                   // Vertex Bins
    "vertexVertex: -10., -7., -5., -3., -1., 1., 3., 5., 7., 10.\n"                         // Vertex Bins (Vertex Binning)
    ;

  TObjArray* lines = defaultBinningStr.Tokenize("\n");
  for (Int_t i=0; i<lines->GetEntriesFast(); i++)
  {
    TString line(lines->At(i)->GetName());
    TString tag = line(0, line.Index(":")+1);
    if (!fCustomBinning.BeginsWith(tag) && !fCustomBinning.Contains(TString("\n") + tag))
      fBinningString += line + "\n";
    else
      AliInfo(Form("Using custom binning for %s", tag.Data()));
  }
  delete lines;
  fBinningString += fCustomBinning;
  
  AliInfo(Form("Used AliTHn Binning:\n%s",fBinningString.Data()));

  //Depending on fEventClass Variable, do one thing or the other...
  /**********************************************************
   
  ======> Modification: Change Event Classification Scheme
    
  ---> fEventClass == "EventPlane"
   
   Default operation with Event Plane 
   
  ---> fEventClass == "Multiplicity"
   
   Work with reference multiplicity (from GetReferenceMultiplicity, which one is decided in the configuration)
   
  ---> fEventClass == "Centrality" 
   
   Work with Centrality Bins

  ***********************************************************/
  if(fEventClass == "Multiplicity"){
    dBinsSingle[0]     = GetBinning(fBinningString, "multiplicity", iBinSingle[0]);
    dBinsPair[0]       = GetBinning(fBinningString, "multiplicity", iBinPair[0]);
    axisTitleSingle[0] = "reference multiplicity";
    axisTitlePair[0]   = "reference multiplicity";
  }
  if(fEventClass == "Centrality"){
    // fine binning in case of vertex Z binning
    if(fVertexBinning){
      dBinsSingle[0]     = GetBinning(fBinningString, "centralityVertex", iBinSingle[0]);
      dBinsPair[0]       = GetBinning(fBinningString, "centralityVertex", iBinPair[0]);
    }
    else{
      dBinsSingle[0]     = GetBinning(fBinningString, "centrality", iBinSingle[0]);
      dBinsPair[0]       = GetBinning(fBinningString, "centrality", iBinPair[0]);
    }
    axisTitleSingle[0] = "Centrality percentile [%]";
    axisTitlePair[0]   = "Centrality percentile [%]";
  }
  if(fEventClass == "EventPlane"){
    dBinsSingle[0]     = GetBinning(fBinningString, "eventPlane", iBinSingle[0]);
    dBinsPair[0]       = GetBinning(fBinningString, "eventPlane", iBinPair[0]);
    axisTitleSingle[0] = "#varphi - #Psi_{2} (a.u.)";
    axisTitlePair[0]   = "#varphi - #Psi_{2} (a.u.)";
  }
  

  // Delta Eta and Delta Phi
  // (coarse binning in case of vertex Z binning)
  if(fVertexBinning){
    dBinsPair[1]       = GetBinning(fBinningString, "deltaEtaVertex", iBinPair[1]);
  }
  else{
    dBinsPair[1]       = GetBinning(fBinningString, "deltaEta", iBinPair[1]);
  }
  axisTitlePair[1]  = "#Delta#eta"; 
  
  dBinsPair[2]       = GetBinning(fBinningString, "deltaPhi", iBinPair[2]);
  axisTitlePair[2]   = "#Delta#varphi (rad)";  
  

  // pT Trig and pT Assoc
  // (coarse binning in case of vertex Z binning)
  if(fVertexBinning){
    dBinsSingle[1]   = GetBinning(fBinningString, "pTVertex", iBinSingle[1]);
    dBinsPair[3]     = GetBinning(fBinningString, "pTVertex", iBinPair[3]);
    dBinsPair[4]     = GetBinning(fBinningString, "pTVertex", iBinPair[4]);
  }
  else{
    dBinsSingle[1]   = GetBinning(fBinningString, "pT", iBinSingle[1]);
    dBinsPair[3]     = GetBinning(fBinningString, "pT", iBinPair[3]);
    dBinsPair[4]     = GetBinning(fBinningString, "pT", iBinPair[4]);
  }
  
  axisTitleSingle[1]  = "p_{T,trig.} (GeV/c)"; 
  axisTitlePair[3]    = "p_{T,trig.} (GeV/c)"; 
  axisTitlePair[4]    = "p_{T,assoc.} (GeV/c)";  
 

  // vertex Z binning or not
  if(fVertexBinning){
    dBinsSingle[2]   = GetBinning(fBinningString, "vertexVertex", iBinSingle[2]);
    dBinsPair[5]     = GetBinning(fBinningString, "vertexVertex", iBinPair[5]);
  }
  else{
    dBinsSingle[2]   = GetBinning(fBinningString, "vertex", iBinSingle[2]);
    dBinsPair[5]     = GetBinning(fBinningString, "vertex", iBinPair[5]);
  }

  axisTitleSingle[2]  = "v_{Z} (cm)"; 
  axisTitlePair[5]    = "v_{Z} (cm)"; 



  // =========================================================
  // Create the Output objects (AliTHn)
  // =========================================================

  TString histName;
  //+ triggered particles
  histName = "fHistP"; 
  if(fCentralityId) histName += fCentralityId.Data();
  fHistP = new AliTHn(histName.Data(),histName.Data(),anaSteps,kTrackVariablesSingle,iBinSingle);
  for (Int_t j=0; j<kTrackVariablesSingle; j++) {
    fHistP->SetBinLimits(j, dBinsSingle[j]);
    fHistP->SetVarTitle(j, axisTitleSingle[j]);
  }

  //- triggered particles
  histName = "fHistN"; 
  if(fCentralityId) histName += fCentralityId.Data();
  fHistN = new AliTHn(histName.Data(),histName.Data(),anaSteps,kTrackVariablesSingle,iBinSingle);
  for (Int_t j=0; j<kTrackVariablesSingle; j++) {
    fHistN->SetBinLimits(j, dBinsSingle[j]);
    fHistN->SetVarTitle(j, axisTitleSingle[j]);
  }
  
  //+- pairs
  histName = "fHistPN";
  if(fCentralityId) histName += fCentralityId.Data();
  fHistPN = new AliTHn(histName.Data(),histName.Data(),anaSteps, kTrackVariablesPair, iBinPair);
  for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fHistPN->SetBinLimits(j, dBinsPair[j]);
    fHistPN->SetVarTitle(j, axisTitlePair[j]);
  }

  //-+ pairs
  histName = "fHistNP";
  if(fCentralityId) histName += fCentralityId.Data();
  fHistNP = new AliTHn(histName.Data(),histName.Data(),anaSteps, kTrackVariablesPair, iBinPair);
  for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fHistNP->SetBinLimits(j, dBinsPair[j]);
    fHistNP->SetVarTitle(j, axisTitlePair[j]);
  }

  //++ pairs
  histName = "fHistPP";
  if(fCentralityId) histName += fCentralityId.Data();
  fHistPP = new AliTHn(histName.Data(),histName.Data(),anaSteps, kTrackVariablesPair, iBinPair);
  for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fHistPP->SetBinLimits(j, dBinsPair[j]);
    fHistPP->SetVarTitle(j, axisTitlePair[j]);
  }

  //-- pairs
  histName = "fHistNN";
  if(fCentralityId) histName += fCentralityId.Data();
  fHistNN = new AliTHn(histName.Data(),histName.Data(),anaSteps, kTrackVariablesPair, iBinPair);
  for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fHistNN->SetBinLimits(j, dBinsPair[j]);
    fHistNN->SetVarTitle(j, axisTitlePair[j]);
  }

  AliInfo("Finished setting up the AliTHn");

  // QA histograms
  fHistHBTbefore        = new TH2D("fHistHBTbefore","before HBT cut",200,-0.5,0.5,200,-0.5,0.5);
  fHistHBTafter         = new TH2D("fHistHBTafter","after HBT cut",200,-0.5,0.5,200,-0.5,0.5);
  fHistPhiStarHBTbefore = new TH2D("fHistPhiStarHBTbefore","before PhiStarHBT cut",200,-0.5,0.5,200,-0.5,0.5);
  fHistPhiStarHBTafter  = new TH2D("fHistPhiStarHBTafter","after PhiStarHBT cut",200,-0.5,0.5,200,-0.5,0.5);
  fHistConversionbefore = new TH3D("fHistConversionbefore","before Conversion cut;#Delta#eta;#Delta#phi;M_{inv}^{2}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistConversionafter  = new TH3D("fHistConversionafter","after Conversion cut;#Delta#eta;#Delta#phi;M_{inv}^{2}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistPsiMinusPhi      = new TH2D("fHistPsiMinusPhi","",4,-0.5,3.5,100,0,2.*TMath::Pi());
  fHistResonancesBefore = new TH3D("fHistResonancesBefore","before resonance cut;#Delta#eta;#Delta#phi;M_{inv}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistResonancesRho    = new TH3D("fHistResonancesRho","after #rho resonance cut;#Delta#eta;#Delta#phi;M_{inv}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistResonancesK0     = new TH3D("fHistResonancesK0","after #rho, K0 resonance cut;#Delta#eta;#Delta#phi;M_{inv}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistResonancesLambda = new TH3D("fHistResonancesLambda","after #rho, K0, Lambda resonance cut;#Delta#eta;#Delta#phi;M_{inv}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistQbefore          = new TH3D("fHistQbefore","before momentum difference cut;#Delta#eta;#Delta#phi;|#Delta p_{T}| (GeV/c)",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistQafter           = new TH3D("fHistQafter","after momentum difference cut;#Delta#eta;#Delta#phi;|#Delta p_{T}| (GeV/c)",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);

  if(axisTitleSingle) delete [] axisTitleSingle;
  if(axisTitlePair) delete [] axisTitlePair;
  
  TH1::AddDirectory(oldStatus);

}

//____________________________________________________________________//
void AliPidBFBase::CalculateBalance(Double_t gReactionPlane,
				     TObjArray *particles, 
				     TObjArray *particlesMixed,
				     Float_t bSign,
				     Double_t kMultorCent,
				     Double_t vertexZ) {
  // Calculates the balance function
  fAnalyzedEvents++;
    
  // Initialize histograms if not done yet
  if(!fHistPN){
    AliWarning("Histograms not yet initialized! --> Will be done now");
    AliWarning("This works only in local mode --> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
    InitHistograms();
  }

  Double_t trackVariablesSingle[kTrackVariablesSingle];
  Double_t trackVariablesPair[kTrackVariablesPair];

  if (!particles){
    AliWarning("particles TObjArray is NULL pointer --> return");
    return;
  }
  
  // define end of particle loops
  Int_t iMax = particles->GetEntriesFast();
  Int_t jMax = iMax;
  if (particlesMixed)
    jMax = particlesMixed->GetEntriesFast();


  // Eta() is extremely time consuming, therefore cache it for the inner loop here:
  TObjArray* particlesSecond = (particlesMixed) ? particlesMixed : particles;

  TArrayF secondEta(jMax);
  TArrayF secondPhi(jMax);
  TArrayF secondPt(jMax);
  TArrayS secondCharge(jMax);
  TArrayD secondCorrection(jMax);

  for (Int_t i=0; i<jMax; i++){
    secondEta[i] = ((AliVParticle*) particlesSecond->At(i))->Eta();
    secondPhi[i] = ((AliVParticle*) particlesSecond->At(i))->Phi();
    secondPt[i]  = ((AliVParticle*) particlesSecond->At(i))->Pt();
    secondCharge[i]  = (Short_t)((AliVParticle*) particlesSecond->At(i))->Charge();
    secondCorrection[i]  = (Double_t)((AliBFBasicParticle*) particlesSecond->At(i))->Correction();   //==========================correction
  }
  
  //TLorenzVector implementation for resonances
  TLorentzVector vectorMother, vectorDaughter[2];
  TParticle pPion, pProton, pRho0, pK0s, pLambda;
  pPion.SetPdgCode(211); //pion
  pRho0.SetPdgCode(113); //rho0
  pK0s.SetPdgCode(310); //K0s
  pProton.SetPdgCode(2212); //proton
  pLambda.SetPdgCode(3122); //Lambda
  Double_t gWidthForRho0 = 0.01;
  Double_t gWidthForK0s = 0.01;
  Double_t gWidthForLambda = 0.006;
  Double_t nSigmaRejection = 3.0;

  // 1st particle loop
  for (Int_t i = 0; i < iMax; i++) {
    //AliVParticle* firstParticle = (AliVParticle*) particles->At(i);
    AliBFBasicParticle* firstParticle = (AliBFBasicParticle*) particles->At(i); //==========================correction
    
    // some optimization
    Float_t firstEta = firstParticle->Eta();
    Float_t firstPhi = firstParticle->Phi();
    Float_t firstPt  = firstParticle->Pt();
    Float_t firstCorrection  = firstParticle->Correction();//==========================correction

    // Event plane (determine psi bin)
    Double_t gPsiMinusPhi    =   0.;
    Double_t gPsiMinusPhiBin = -10.;
    gPsiMinusPhi   = TMath::Abs(firstPhi - gReactionPlane);
    //in-plane
    if((gPsiMinusPhi <= 7.5*TMath::DegToRad())||
       ((172.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 187.5*TMath::DegToRad())))
      gPsiMinusPhiBin = 0.0;
    //intermediate
    else if(((37.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 52.5*TMath::DegToRad()))||
	    ((127.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 142.5*TMath::DegToRad()))||
	    ((217.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 232.5*TMath::DegToRad()))||
	    ((307.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 322.5*TMath::DegToRad())))
      gPsiMinusPhiBin = 1.0;
    //out of plane
    else if(((82.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 97.5*TMath::DegToRad()))||
	    ((262.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 277.5*TMath::DegToRad())))
      gPsiMinusPhiBin = 2.0;
    //everything else
    else 
      gPsiMinusPhiBin = 3.0;
    
    fHistPsiMinusPhi->Fill(gPsiMinusPhiBin,gPsiMinusPhi);

    Short_t  charge1 = (Short_t) firstParticle->Charge();
    
    trackVariablesSingle[0]    =  gPsiMinusPhiBin;
    trackVariablesSingle[1]    =  firstPt;
      if(fEventClass=="Multiplicity" || fEventClass == "Centrality" ) trackVariablesSingle[0] = kMultorCent;
    trackVariablesSingle[2]    =  vertexZ;

    
    //fill single particle histograms
    if(charge1 > 0)      fHistP->Fill(trackVariablesSingle,0,firstCorrection); //==========================correction
    else if(charge1 < 0) fHistN->Fill(trackVariablesSingle,0,firstCorrection);  //==========================correction
    
    // 2nd particle loop
    for(Int_t j = 0; j < jMax; j++) {   

      if(!particlesMixed && j == i) continue; // no auto correlations (only for non mixing)

      // pT,Assoc < pT,Trig (if momentum ordering is switched ON)
      if(fMomentumOrdering){
	if(firstPt < secondPt[j]) 
	  continue;
      }

      Short_t charge2 = secondCharge[j];
      
      trackVariablesPair[0]    =  trackVariablesSingle[0];
      trackVariablesPair[1]    =  firstEta - secondEta[j];  // delta eta
      trackVariablesPair[2]    =  firstPhi - secondPhi[j];  // delta phi
      //if (trackVariablesPair[2] > 180.)   // delta phi between -180 and 180 
      //trackVariablesPair[2] -= 360.;
      //if (trackVariablesPair[2] <  - 180.) 
      //trackVariablesPair[2] += 360.;
      if (trackVariablesPair[2] > TMath::Pi()) // delta phi between -pi and pi 
	trackVariablesPair[2] -= 2.*TMath::Pi();
      if (trackVariablesPair[2] <  - TMath::Pi()) 
	trackVariablesPair[2] += 2.*TMath::Pi();
      if (trackVariablesPair[2] <  - TMath::Pi()/2.) 
      trackVariablesPair[2] += 2.*TMath::Pi();
      
      trackVariablesPair[3]    =  firstPt;      // pt trigger
      trackVariablesPair[4]    =  secondPt[j];  // pt
      trackVariablesPair[5]    =  vertexZ;      // z of the primary vertex
      
      //Exclude resonances for the calculation of pairs by looking 
      //at the invariant mass and not considering the pairs that 
      //fall within 3sigma from the mass peak of: rho0, K0s, Lambda
      if(fResonancesCut) {
	if (charge1 * charge2 < 0) {

	  //rho0
	  vectorDaughter[0].SetPtEtaPhiM(firstPt,firstEta,firstPhi,pPion.GetMass());
	  vectorDaughter[1].SetPtEtaPhiM(secondPt[j],secondEta[j],secondPhi[j],pPion.GetMass());
	  vectorMother = vectorDaughter[0] + vectorDaughter[1];
	  fHistResonancesBefore->Fill(trackVariablesPair[1],trackVariablesPair[2],vectorMother.M());
	  if(TMath::Abs(vectorMother.M() - pRho0.GetMass()) <= nSigmaRejection*gWidthForRho0)
	    continue;
	  fHistResonancesRho->Fill(trackVariablesPair[1],trackVariablesPair[2],vectorMother.M());
	  
	  //K0s
	  if(TMath::Abs(vectorMother.M() - pK0s.GetMass()) <= nSigmaRejection*gWidthForK0s)
	    continue;
	  fHistResonancesK0->Fill(trackVariablesPair[1],trackVariablesPair[2],vectorMother.M());
	  
	  
	  //Lambda
	  vectorDaughter[0].SetPtEtaPhiM(firstPt,firstEta,firstPhi,pPion.GetMass());
	  vectorDaughter[1].SetPtEtaPhiM(secondPt[j],secondEta[j],secondPhi[j],pProton.GetMass());
	  vectorMother = vectorDaughter[0] + vectorDaughter[1];
	  if(TMath::Abs(vectorMother.M() - pLambda.GetMass()) <= nSigmaRejection*gWidthForLambda)
	    continue;
	  
	  vectorDaughter[0].SetPtEtaPhiM(firstPt,firstEta,firstPhi,pProton.GetMass());
	  vectorDaughter[1].SetPtEtaPhiM(secondPt[j],secondEta[j],secondPhi[j],pPion.GetMass());
	  vectorMother = vectorDaughter[0] + vectorDaughter[1];
	  if(TMath::Abs(vectorMother.M() - pLambda.GetMass()) <= nSigmaRejection*gWidthForLambda)
	    continue;
	  fHistResonancesLambda->Fill(trackVariablesPair[1],trackVariablesPair[2],vectorMother.M());
	
	}//unlike-sign only
      }//resonance cut

      if(fHBTCut && charge1 * charge2 > 0){  // VERSION 2 (only for LS)
	//if( dphi < 3 || deta < 0.01 ){   // VERSION 1
	//  continue;
	
	Double_t deta = firstEta - secondEta[j];
	Double_t dphi = firstPhi - secondPhi[j];
	if(dphi > TMath::Pi())
	  dphi = secondPhi[j] - firstPhi;

	// for QA: get dphistar in the middle of the TPC R = 1.65
	Float_t  dphistarMiddle = GetDPhiStar(firstPhi, firstPt, charge1, secondPhi[j], secondPt[j], charge2, 1.65, bSign);

	// see e.g. https://indico.cern.ch/materialDisplay.py?contribId=36&sessionId=6&materialId=slides&confId=142700
	fHistHBTbefore->Fill(deta,dphi);
	fHistPhiStarHBTbefore->Fill(deta,dphistarMiddle);
	
	// optimization
	if (TMath::Abs(deta) < fHBTCutValue * 2.5 * 3) //fHBTCutValue = 0.02 [default for dphicorrelations]
	  {
	    Float_t phi1rad = firstPhi;
	    Float_t phi2rad = secondPhi[j];
	    
	    // check first boundaries to see if is worth to loop and find the minimum
	    Float_t dphistar1 = GetDPhiStar(phi1rad, firstPt, charge1, phi2rad, secondPt[j], charge2, 0.8, bSign);
	    Float_t dphistar2 = GetDPhiStar(phi1rad, firstPt, charge1, phi2rad, secondPt[j], charge2, 2.5, bSign);
	    
	    const Float_t kLimit = fHBTCutValue * 3;
	    
	    Float_t dphistarminabs = 1e5;
	    //Float_t dphistarmin = 1e5;
	    
	    if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0 ) {
	      for (Double_t rad=0.8; rad<2.51; rad+=0.01) {
		Float_t dphistar = GetDPhiStar(phi1rad, firstPt, charge1, phi2rad, secondPt[j], charge2, rad, bSign);
		Float_t dphistarabs = TMath::Abs(dphistar);
		
		if (dphistarabs < dphistarminabs) {
		  //dphistarmin = dphistar;
		  dphistarminabs = dphistarabs;
		}
	      }
	      
	      if (dphistarminabs < fHBTCutValue && TMath::Abs(deta) < fHBTCutValue) {
		continue;
	      }
	    }
	  }
	fHistHBTafter->Fill(deta,dphi);
	fHistPhiStarHBTafter->Fill(deta,dphistarMiddle);
      }//HBT cut
	
      // conversions
      if(fConversionCut) {
	if (charge1 * charge2 < 0) {
	  Double_t deta = firstEta - secondEta[j];
	  Double_t dphi = firstPhi - secondPhi[j];
	  
	  Float_t m0 = 0.510e-3;
	  Float_t tantheta1 = 1e10;
	  Float_t phi1rad = firstPhi;
	  Float_t phi2rad = secondPhi[j];
	  
	  if (firstEta < -1e-10 || firstEta > 1e-10)
	    tantheta1 = 2 * TMath::Exp(-firstEta) / ( 1 - TMath::Exp(-2*firstEta));
	  
	  Float_t tantheta2 = 1e10;
	  if (secondEta[j] < -1e-10 || secondEta[j] > 1e-10)
	    tantheta2 = 2 * TMath::Exp(-secondEta[j]) / ( 1 - TMath::Exp(-2*secondEta[j]));
	  
	  Float_t e1squ = m0 * m0 + firstPt * firstPt * (1.0 + 1.0 / tantheta1 / tantheta1);
	  Float_t e2squ = m0 * m0 + secondPt[j] * secondPt[j] * (1.0 + 1.0 / tantheta2 / tantheta2);
	  
	  Float_t masssqu = 2 * m0 * m0 + 2 * ( TMath::Sqrt(e1squ * e2squ) - ( firstPt * secondPt[j] * ( TMath::Cos(phi1rad - phi2rad) + 1.0 / tantheta1 / tantheta2 ) ) );

	  fHistConversionbefore->Fill(deta,dphi,masssqu);
	  
	  if (masssqu < fInvMassCutConversion*fInvMassCutConversion){
	    continue;
	  }
	  fHistConversionafter->Fill(deta,dphi,masssqu);
	}
      }//conversion cut

      // momentum difference cut - suppress femtoscopic effects
      if(fQCut){ 

	//Double_t ptMin        = 0.1; //const for the time being (should be changeable later on)
	Double_t ptDifference = TMath::Abs( firstPt - secondPt[j]);

	fHistQbefore->Fill(trackVariablesPair[1],trackVariablesPair[2],ptDifference);
	if(ptDifference < fDeltaPtMin) continue;
	fHistQafter->Fill(trackVariablesPair[1],trackVariablesPair[2],ptDifference);

      }

      if( charge1 > 0 && charge2 < 0)  fHistPN->Fill(trackVariablesPair,0,firstCorrection*secondCorrection[j]); //==========================correction
      else if( charge1 < 0 && charge2 > 0)  fHistNP->Fill(trackVariablesPair,0,firstCorrection*secondCorrection[j]);//==========================correction 
      else if( charge1 > 0 && charge2 > 0)  fHistPP->Fill(trackVariablesPair,0,firstCorrection*secondCorrection[j]);//==========================correction 
      else if( charge1 < 0 && charge2 < 0)  fHistNN->Fill(trackVariablesPair,0,firstCorrection*secondCorrection[j]);//==========================correction 
      else {
	//AliWarning(Form("Wrong charge combination: charge1 = %d and charge2 = %d",charge,charge2));
	continue;
      }
    }//end of 2nd particle loop
  }//end of 1st particle loop
}  

//____________________________________________________________________//
Float_t AliPidBFBase::GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign) { 
  //
  // calculates dphistar
  //
  Float_t dphistar = phi1 - phi2 - charge1 * bSign * TMath::ASin(0.075 * radius / pt1) + charge2 * bSign * TMath::ASin(0.075 * radius / pt2);
  
  static const Double_t kPi = TMath::Pi();
  
  // circularity
//   if (dphistar > 2 * kPi)
//     dphistar -= 2 * kPi;
//   if (dphistar < -2 * kPi)
//     dphistar += 2 * kPi;
  
  if (dphistar > kPi)
    dphistar = kPi * 2 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi * 2 - dphistar;
  if (dphistar > kPi) // might look funny but is needed
    dphistar = kPi * 2 - dphistar;
  
  return dphistar;
}

//____________________________________________________________________//
Double_t* AliPidBFBase::GetBinning(const char* configuration, const char* tag, Int_t& nBins)
{
  // This method is a copy from AliUEHist::GetBinning
  // takes the binning from <configuration> identified by <tag>
  // configuration syntax example:
  // eta: 2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4
  // phi: .....
  //
  // returns bin edges which have to be deleted by the caller
  
  TString config(configuration);
  TObjArray* lines = config.Tokenize("\n");
  for (Int_t i=0; i<lines->GetEntriesFast(); i++)
  {
    TString line(lines->At(i)->GetName());
    if (line.BeginsWith(TString(tag) + ":"))
    {
      line.Remove(0, strlen(tag) + 1);
      line.ReplaceAll(" ", "");
      TObjArray* binning = line.Tokenize(",");
      Double_t* bins = new Double_t[binning->GetEntriesFast()];
      for (Int_t j=0; j<binning->GetEntriesFast(); j++)
	bins[j] = TString(binning->At(j)->GetName()).Atof();
      
      nBins = binning->GetEntriesFast() - 1;

      delete binning;
      delete lines;
      return bins;
    }
  }
  
  delete lines;
  AliFatal(Form("Tag %s not found in %s", tag, configuration));
  return 0;
}

