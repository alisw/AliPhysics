//////////////////////////////////////////////////////////////////////
//
// $Id: AliFlowEvent.cxx 18618 2007-05-16 15:38:22Z snelling $
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//
//_____________________________________________________________
//
// Description: 
//         AliFlowEvent is the basic class to perform flow study in ALICE. 
// The AliFlowEvent object contains data mebers wich summarize the ESD 
// informations that are most useful for flow study, together with a 
// collection of tracks (AliFlowTrackCollection) and reconstructed v0s. 
// The class also implements member functions to calculate Multiplicity, 
// Mean P and Pt, Q (the event plane vector), Psi (the event plane angle), 
// normalized q (Q/sqrt(Mult) and Q/sqrt(sum of weights^2)), moreover, 
// functions to separate random or eta sub-events, to fill the bayesian vector 
// of particle abundance (for bayesian p.id calculation), to plug-in weights
// (to make the tracks' distribution isotropic and to enhance the resolution), 
// and to select the tracks that are used in the event plane calculation (for 
// each harmonic, selection, and subevent. 
// AliFlowEvent supports the standard analysis as well as the Cumulant 
// Analysis: method to extract the generating functions for the old and new 
// cumulant methods is also there.
// The Flow Event structure (AliFlowEvent, AliFlowTrack, AliFlowV0 and 
// AliFlowSelection) is independent from AliRoot, i.e. once the FlowEvents 
// are filled from the AliESD or the MC KineTree (see AliFlowInterface: here 
// AliRoot is needed), the flow analysis itself can be performed just under 
// ROOT (see AliFlowAnalysisMaker).
//
// AliFlowEvent is adapted from the original class StFlowEvent, succesfully 
// employed to study flow in the STAR experiment at RICH.
// Original Authors:                 Raimond Snellings & Art Poskanzer
//

#include "AliFlowEvent.h"
#include "AliFlowTrack.h"
#include "AliFlowV0.h"
#include "AliFlowSelection.h"
#include "AliFlowConstants.h"

#include "TVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include <iostream>
using namespace std; //required for resolving the 'cout' symbol

// - Flags & Sets
Bool_t  AliFlowEvent::fgPtWgt	       = kFALSE ;  // gives pT-based weights
Bool_t  AliFlowEvent::fgEtaWgt	       = kFALSE ;  // gives eta-based weights
Bool_t  AliFlowEvent::fgOnePhiWgt      = kTRUE  ;  // kTRUE --> ENABLEs SINGLE WEIGHT ISTOGRAM , kFALSE --> ENABLEs 3 WEIGHT ISTOGRAMS
Bool_t  AliFlowEvent::fgNoWgt	       = kFALSE ;  // No Weight is used
Bool_t  AliFlowEvent::fgCustomRespFunc = kFALSE ;  // custom "detector response function" is used for P.Id (see AliFlowTrack)
// - Eta Sub-Events (later used to calculate the resolution)
Int_t   AliFlowEvent::fgEtaSubs        = 0 ;       // makes random subevents (0 = random , 1 = eta , -1 = charged)

ClassImp(AliFlowEvent) 
//-----------------------------------------------------------
AliFlowEvent::AliFlowEvent(Int_t lenght):
  fEventID(0), 
  fRunID(0), 
  fOrigMult(0), 
  fL0Trigger(0), 
  fZDCpart(0), 
  fCentrality(-1), 
  fTrackCollection(0x0), 
  fV0Collection(0x0), 
  fDone(kFALSE)	
{
 // not-Default constructor: initializes the ObjArray of FlowTracks and FlowV0s, 
 // cleans the internal variables, sets all the weights to 1, sets default flags.
 // USE THIS WHEN CREATING ALIFLOWEVENTS
 
 for(int zz=0;zz<3;zz++) { fZDCenergy[zz] = 0. ; }
 for(int i=0;i<AliFlowConstants::kHars;i++) { fExtPsi[i] = 0. ; fExtRes[i] = 0.; }
 
 // Make a new track collection
 fTrackCollection =  new TClonesArray("AliFlowTrack",lenght) ; fTrackCollection->BypassStreamer(kTRUE) ;
 fV0Collection    =  new TClonesArray("AliFlowV0",lenght)    ; fV0Collection->BypassStreamer(kTRUE)    ;
 
 // Set Weights Arrays to 1 (default)
 for(int nS=0;nS<AliFlowConstants::kSels;nS++)
 {
  for(int nH=0;nH<AliFlowConstants::kHars;nH++) 
  {
   for(int nP=0;nP<AliFlowConstants::kPhiBins;nP++) 
   { 
    // enable this with: SetOnePhiWgt()  
    fPhiWgt[nS][nH][nP] = 1.	  ; // cout << nS << nH << nP << "  val  " << fPhiWgt[nS][nH][nP] << endl ; 
    // enable these with: SetFirstLastPhiWgt()  
    fPhiWgtPlus[nS][nH][nP]  = 1. ; // cout << nS << nH << nP << "  val  " << fPhiWgtPlus[nS][nH][nP] << endl ; 
    fPhiWgtMinus[nS][nH][nP] = 1. ; // cout << nS << nH << nP << "  val  " << fPhiWgtMinus[nS][nH][nP] << endl ; 
    fPhiWgtCross[nS][nH][nP] = 1. ; // cout << nS << nH << nP << "  val  " << fPhiWgtCross[nS][nH][nP] << endl ; 
   }
  }
 }
 //for(int nH=0;nH<AliFlowConstants::kHars;nH++) { fExtPsi[nH] = 0. ; fExtRes[nH] = 0. ; }
 
 // The Expected particles abundance is taken directly from AliFlowConstants::fgBayesian[] (see Bayesian P.Id.)
 
}
//-----------------------------------------------------------
AliFlowEvent::AliFlowEvent():
  fEventID(0), 
  fRunID(0), 
  fOrigMult(0), 
  fL0Trigger(0), 
  fZDCpart(0), 
  fCentrality(-1), 
  fTrackCollection(0x0), 
  fV0Collection(0x0), 
  fDone(kFALSE)	
{
 // Default constructor :  USE THIS WHEN READING ALIFLOWEVENTS

 // Set Weights Arrays to 1 (default)
 for(int nS=0;nS<AliFlowConstants::kSels;nS++)
 {
  for(int nH=0;nH<AliFlowConstants::kHars;nH++) 
  {
   for(int nP=0;nP<AliFlowConstants::kPhiBins;nP++) 
   { 
    // enable this with: SetOnePhiWgt()  
    fPhiWgt[nS][nH][nP] = 1.	  ; // cout << nS << nH << nP << "  val  " << fPhiWgt[nS][nH][nP] << endl ; 
    // enable these with: SetFirstLastPhiWgt()  
    fPhiWgtPlus[nS][nH][nP]  = 1. ; // cout << nS << nH << nP << "  val  " << fPhiWgtPlus[nS][nH][nP] << endl ; 
    fPhiWgtMinus[nS][nH][nP] = 1. ; // cout << nS << nH << nP << "  val  " << fPhiWgtMinus[nS][nH][nP] << endl ; 
    fPhiWgtCross[nS][nH][nP] = 1. ; // cout << nS << nH << nP << "  val  " << fPhiWgtCross[nS][nH][nP] << endl ; 
   }
  }
 }
}
//-----------------------------------------------------------
AliFlowEvent::~AliFlowEvent() 
{
 // Default distructor: deletes the ObjArrays. 
 
 if(fTrackCollection) { delete fTrackCollection ; } // fTrackCollection->Delete() ; 
 if(fV0Collection)    { delete fV0Collection ; }    // fV0Collection->Delete()    ; 
}
//-------------------------------------------------------------

//-------------------------------------------------------------
Double_t AliFlowEvent::PhiWeightRaw(Int_t selN, Int_t harN, AliFlowTrack* pFlowTrack) const  
{
 // Weight for making the event plane isotropic in the lab.

 float phi = pFlowTrack->Phi() ; if(phi < 0.) { phi += 2*TMath::Pi() ; }
 Double_t  eta = (Double_t)pFlowTrack->Eta() ;
 int n = (int)((phi/(2*TMath::Pi()))*AliFlowConstants::kPhiBins);

 Double_t phiWgt = 1. ;
 if(OnePhiWgt()) 
 {
  phiWgt = (Double_t)fPhiWgt[selN][harN][n]; //cout << "Selection " << selN << " ; Harmonic " << harN << " ; PhiBin " << n << "  -  Wgt = " << phiWgt << endl ; 
 } 
 else if(FirstLastPhiWgt())
 {
  float zFirstPoint = pFlowTrack->ZFirstPoint(); // float zLastPoint = pFlowTrack->ZLastPoint();
  
  if (zFirstPoint > 0. && eta > 0.)     { phiWgt = (Double_t)fPhiWgtPlus[selN][harN][n] ; } 
  else if(zFirstPoint < 0. && eta < 0.) { phiWgt = (Double_t)fPhiWgtMinus[selN][harN][n] ; } 
  else  				{ phiWgt = (Double_t)fPhiWgtCross[selN][harN][n] ; }
 } 

 return phiWgt ;
}
//-------------------------------------------------------------
Double_t AliFlowEvent::Weight(Int_t selN, Int_t harN, AliFlowTrack* pFlowTrack) const  
{ 
 // Weight for enhancing the resolution (eta gives sign +/- for Odd Harmonics)

 if(selN>AliFlowConstants::kSels) { selN = 0 ; }
 bool oddHar = (harN+1) % 2 ;
 Double_t phiWgt = 1. ;
 if(PtWgt()) 
 {
  Double_t pt = (Double_t)pFlowTrack->Pt();
  if(pt < AliFlowConstants::fgPtWgtSaturation) { phiWgt *= pt ; } 
  else 			   	  { phiWgt *= AliFlowConstants::fgPtWgtSaturation ; } // pt weighting going constant
 }
 Double_t eta = (Double_t)pFlowTrack->Eta();
 Double_t etaAbs = TMath::Abs(eta);
 if(EtaWgt() && oddHar)     { phiWgt *= etaAbs ; }
 if(oddHar && eta < 0.)     { phiWgt *= -1. ; }

 return phiWgt ;
}
//-------------------------------------------------------------
Double_t AliFlowEvent::PhiWeight(Int_t selN, Int_t harN, AliFlowTrack* pFlowTrack) const
{
 // Weight for making the event plane isotropic in the lab and enhancing the resolution
 // (it simply rerurns PhiWeightRaw() * Weight()). If fgNoWgt = kTRUE, returns +/-1 ,
 // basing on Sign(eta), for odd harmonics .
 
 if(fgNoWgt) // no weights (but +/- according to eta)
 { 
  bool oddHar = (harN+1) % 2 ;
  if(oddHar) { return TMath::Sign((Double_t)1.,(Double_t)pFlowTrack->Eta()) ; }
  else       { return (Double_t)1. ; }
 } 
 Double_t phiWgtRaw = PhiWeightRaw(selN, harN, pFlowTrack);
 Double_t weight = Weight(selN, harN, pFlowTrack);
 if(AliFlowConstants::fgDebug) { cout << "[PhiWeight]: phiWgtRaw = " << phiWgtRaw << " , weight = " << weight << " , eta = " << pFlowTrack->Eta() << endl ; }

 return phiWgtRaw * weight;
}
//-------------------------------------------------------------
Int_t AliFlowEvent::Mult(AliFlowSelection* pFlowSelect) const	
{
 // Multiplicity of tracks in the specified Selection 
 
 if(fDone) 
 {
  int sub = pFlowSelect->Sub() ;
  if(sub<0) { return fMult[pFlowSelect->Sel()][pFlowSelect->Har()] ; }
  else      { return fMultSub[sub][pFlowSelect->Sel()][pFlowSelect->Har()] ; }
 }
 // -
 Int_t mult = 0;
 Int_t itr ;
 for(itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if(pFlowSelect->Select(pFlowTrack)) { mult++ ; }
 }
 return mult;
}
//-------------------------------------------------------------
Float_t AliFlowEvent::MeanPt(AliFlowSelection* pFlowSelect) const
{
 // Mean pt of tracks in the specified Selection 

 Double_t meanPt = 0. ;
 Float_t sumPt = 0. ;
 UInt_t mult = 0 ;
 Int_t itr ;
 for(itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if (pFlowSelect->Select(pFlowTrack)) 
  {
   sumPt += pFlowTrack->Pt();
   mult++;
  }
 }
 if(mult) { meanPt = sumPt/(float)mult ; }
 
 return meanPt ;
}
//-------------------------------------------------------------
TVector2 AliFlowEvent::Q(AliFlowSelection* pFlowSelect) const
{
 // Event plane vector for the specified Selection 

 if(fDone) 
 { 
  int sub = pFlowSelect->Sub() ;
  if(sub<0) { return fQ[pFlowSelect->Sel()][pFlowSelect->Har()] ; }
  else      { return fQSub[sub][pFlowSelect->Sel()][pFlowSelect->Har()] ; }
 } 
 // -
 TVector2 mQ ;
 Double_t mQx=0. , mQy=0. ;
 int    selN  = pFlowSelect->Sel() ;
 int    harN  = pFlowSelect->Har() ;
 double order = (double)(harN + 1) ;

 Int_t itr ;
 for(itr=0;itr<TrackCollection()->GetEntries();itr++)
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if(pFlowSelect->Select(pFlowTrack)) 
  {
   double phiWgt = PhiWeight(selN, harN, pFlowTrack);  
   float phi = pFlowTrack->Phi();                      
   mQx += phiWgt * cos(phi * order) ;
   mQy += phiWgt * sin(phi * order) ;
   if(AliFlowConstants::fgDebug) { cout << itr << " phi = " << phi << " ,  wgt = " << phiWgt << endl ; }
  }
 }
 mQ.Set(mQx, mQy);

 return mQ;
}
//-------------------------------------------------------------
Float_t AliFlowEvent::Psi(AliFlowSelection* pFlowSelect) const
{
 // Event plane angle for the specified Selection 

 int   harN = pFlowSelect->Har() ;
 float order = (float)(harN + 1) ;
 Float_t psi = 0. ;

 TVector2 mQ = Q(pFlowSelect);
 if(mQ.Mod())  // if vector is not 0
 {
  psi= mQ.Phi() / order ;
  if (psi < 0.) { psi += 2*TMath::Pi() / order ; }
 }
 return psi;
}
//-------------------------------------------------------------
TVector2 AliFlowEvent::NormQ(AliFlowSelection* pFlowSelect) const  	
{
 // Normalized Q = Q/sqrt(sum of weights^2) for the specified Selection 

 if(fDone) 
 {
  TVector2 mQ = fQ[pFlowSelect->Sel()][pFlowSelect->Har()] ;   
  double sumOfWeightSqr = fSumOfWeightSqr[pFlowSelect->Sel()][pFlowSelect->Har()] ;
  if(sumOfWeightSqr) { mQ /= TMath::Sqrt(sumOfWeightSqr) ; }
  else { mQ.Set(0.,0.) ; }
  return mQ ;
 }
 // -
 TVector2 mQ ;
 Double_t mQx=0. , mQy=0. ;
 double sumOfWeightSqr = 0 ;
 int selN     = pFlowSelect->Sel() ;
 int harN     = pFlowSelect->Har() ;
 double order = (double)(harN + 1) ;
 Int_t itr ;
 for(itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if (pFlowSelect->Select(pFlowTrack)) 
  {
   double phiWgt = PhiWeight(selN, harN, pFlowTrack);
   sumOfWeightSqr += phiWgt*phiWgt;

   float phi = pFlowTrack->Phi();
   mQx += phiWgt * cos(phi * order);
   mQy += phiWgt * sin(phi * order);
  }
 }
 if(sumOfWeightSqr) { mQ.Set(mQx/TMath::Sqrt(sumOfWeightSqr), mQy/TMath::Sqrt(sumOfWeightSqr)); }
 else { mQ.Set(0.,0.); }
  
 return mQ;
}
//-------------------------------------------------------------
Float_t AliFlowEvent::OldQ(AliFlowSelection* pFlowSelect) const 
{ 
 // Magnitude of normalized Q vector (without pt or eta weighting) for the specified Selection 

 TVector2 mQ ;
 Double_t mQx = 0. , mQy = 0. ;
 int selN     = pFlowSelect->Sel() ;
 int harN     = pFlowSelect->Har() ;
 double order = (double)(harN + 1) ;
 double sumOfWeightSqr = 0 ;

 Int_t itr ;
 for(itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if(pFlowSelect->Select(pFlowTrack)) 
  {
   double phiWgt = PhiWeightRaw(selN, harN, pFlowTrack); // Raw
   sumOfWeightSqr += phiWgt*phiWgt;
   float phi = pFlowTrack->Phi();
   mQx += phiWgt * cos(phi * order);
   mQy += phiWgt * sin(phi * order);
  }
 }
 if(sumOfWeightSqr) { mQ.Set(mQx/TMath::Sqrt(sumOfWeightSqr), mQy/TMath::Sqrt(sumOfWeightSqr)); }
 else { mQ.Set(0.,0.); }
  
 return mQ.Mod();
}
//-----------------------------------------------------------------------
Double_t AliFlowEvent::NewG(AliFlowSelection* pFlowSelect, Double_t Zx, Double_t Zy) const
{ 
 // Generating function for the new cumulant method. Eq. 3 in the Practical Guide 

 int selN     = pFlowSelect->Sel();
 int harN     = pFlowSelect->Har();
 double order = (double)(harN + 1);

 double mult = (double)Mult(pFlowSelect);
 Double_t theG = 1.;

 Int_t itr ;
 for(itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if (pFlowSelect->Select(pFlowTrack)) 
  {
   double phiWgt = PhiWeight(selN, harN, pFlowTrack);	   
   float phi = pFlowTrack->Phi();
   theG *= (1. + (phiWgt/mult) * (2.* Zx * cos(phi * order) + 2.* Zy * sin(phi * order) ) );		
  }
 }
 return theG;
}
//-----------------------------------------------------------------------
Double_t AliFlowEvent::OldG(AliFlowSelection* pFlowSelect, Double_t Zx, Double_t Zy) const
{ 
 // Generating function for the old cumulant method (if expanded in Taylor
 // series, one recovers NewG() in new new cumulant method)

 TVector2 normQ = NormQ(pFlowSelect) ;

 return exp(2*Zx*normQ.X() + 2*Zy*normQ.Y());
}
//-----------------------------------------------------------------------
Double_t AliFlowEvent::SumWeightSquare(AliFlowSelection* pFlowSelect) const
{
 // Return sum of weights^2 for the specified Selection (used for normalization)

 if(fDone) { return fSumOfWeightSqr[pFlowSelect->Sel()][pFlowSelect->Har()] ; }

 int selN = pFlowSelect->Sel();
 int harN = pFlowSelect->Har();
 Double_t sumOfWeightSqr = 0;
 Int_t itr ;
 for(itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if (pFlowSelect->Select(pFlowTrack)) 
  {
   double phiWgt = PhiWeight(selN, harN, pFlowTrack);
   sumOfWeightSqr += phiWgt*phiWgt;
  }
 }
 if(sumOfWeightSqr < 0.) { return Mult(pFlowSelect) ; }

 return sumOfWeightSqr;
}
//-------------------------------------------------------------
Double_t AliFlowEvent::WgtMultQ4(AliFlowSelection* pFlowSelect) const 
{ 
 // Used only for the old cumulant method, for getting q4 when weight is on.
 // Replace multiplicity in Eq.(74b) by this quantity when weight is on.
 // This is derived based on (A4) in the old cumulant paper.

 int selN = pFlowSelect->Sel();
 int harN = pFlowSelect->Har();
 double theMult        = 0.;
 double theMeanWj4     = 0.;
 double theMeanWj2     = 0.;
 double theSumOfWgtSqr = 0;
 double phiWgtSq;

 Int_t itr ;
 for(itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if (pFlowSelect->Select(pFlowTrack)) 
  {
   double phiWgt   = PhiWeight(selN, harN, pFlowTrack);
   phiWgtSq	   = phiWgt*phiWgt;
   theSumOfWgtSqr += phiWgtSq;
   theMeanWj4	  += phiWgtSq*phiWgtSq;
   theMult	  += 1.;      
  }
 }
 if (theMult <= 0.) return theMult;

 theMeanWj4 /= theMult;
 theMeanWj2  = theSumOfWgtSqr / theMult;

 return (theSumOfWgtSqr*theSumOfWgtSqr)/(theMult*(-theMeanWj4+2*theMeanWj2*theMeanWj2));
}
//-------------------------------------------------------------
Double_t AliFlowEvent::WgtMultQ6(AliFlowSelection* pFlowSelect) const 
{ 
 // Used only for the old cumulant method. For getting q6 when weight is on.
 // Replace multiplicity in Eq.(74c) by this quantity when weight is on.
 // This is derived based on (A4) in the old cumulant paper.

 int selN = pFlowSelect->Sel();
 int harN = pFlowSelect->Har();
 double theMult        = 0.;
 double theMeanWj6     = 0.;
 double theMeanWj4     = 0.;
 double theMeanWj2     = 0.;
 double theSumOfWgtSqr = 0;
 double phiWgtSq;

 Int_t itr ;
 for(itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if (pFlowSelect->Select(pFlowTrack)) 
  {
   double phiWgt   = PhiWeight(selN, harN, pFlowTrack);
   phiWgtSq	   = phiWgt*phiWgt;
   theSumOfWgtSqr += phiWgtSq;
   theMeanWj4	  += phiWgtSq*phiWgtSq;
   theMeanWj6	  += phiWgtSq*phiWgtSq*phiWgtSq;
   theMult	  += 1.;
  }
 } 
 if (theMult <= 0.) return theMult*theMult;

 theMeanWj6 /= theMult;
 theMeanWj4 /= theMult;
 theMeanWj2  = theSumOfWgtSqr / theMult;

 return 4.*(theSumOfWgtSqr*theSumOfWgtSqr*theSumOfWgtSqr)/(theMult*(theMeanWj6-9.*theMeanWj2*theMeanWj4+12.*theMeanWj2*theMeanWj2*theMeanWj2));
}
//-------------------------------------------------------------
void AliFlowEvent::SetSelections(AliFlowSelection* pFlowSelect) 	 
{
 // Sets the selection of tracks used in the Reaction Plane calculation 
 // for the specific Harmonic and Selection - this does not cut trow away 
 // tracks from the event, neither exclude them from flow study. See also 
 // the AliFlowSelection class.
 // Strategy of Selection: 
 //                     For the specific Harmonic and Selection, IF cuts 
 // are defined (such that low<high) and IF the track satisfies them, THEN 
 // the respective track's flag (bool AliFlowTrack::mSelection[har][sel]) 
 // is set kTRUE so that the track, from now on, will be included in the 
 // R.P. determination for that selection.
 // If NO cuts are defined -> ALL Flags are setted kTRUE (all particle 
 // used for all the Reaction Planes -> no difference in Psi[har][sel]).
 // -------------------------------------------------------------------
 // The first selection (all harmonics) is set kTRUE : no conditions.

 Int_t itr ;
 for(itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  pFlowTrack->ResetSelection() ;  		// this re-sets all the mSelection flags to 0

 // * this sets all the selection n.[0] flag kTRUE (all harmonics) *
  for(int harN=0;harN<AliFlowConstants::kHars;harN++) { pFlowTrack->SetSelect(harN,0) ; } 

 // Track need to be Constrainable
  if(pFlowSelect->ConstrainCut() && !pFlowTrack->IsConstrainable()) continue ;

 // PID - gets rid of the track with AliFlowTrack::Pid() != AliFlowSelection::Pid() (if there)
  if(pFlowSelect->Pid()[0] != '\0') 
  {
   if(strstr(pFlowSelect->Pid(), "h")!=0) 
   {
    int charge = pFlowTrack->Charge() ;
    if(strcmp("h+",pFlowSelect->Pid())==0 && charge != 1)  continue;
    if(strcmp("h-",pFlowSelect->Pid())==0 && charge != -1) continue;
   } 
   else 
   {
    Char_t pid[10];
    strcpy(pid, pFlowTrack->Pid());
    if(strstr(pid, pFlowSelect->Pid())==0) continue;
   }
  }
  double eta = (double)(pFlowTrack->Eta());      
  float  pt  = pFlowTrack->Pt();   		 
  float  gDca = pFlowTrack->Dca() ; 			

 // Global DCA - gets rid of the track with DCA outside the range
  if((pFlowSelect->DcaGlobalCutHi()>pFlowSelect->DcaGlobalCutLo()) && (gDca<pFlowSelect->DcaGlobalCutLo() || gDca>pFlowSelect->DcaGlobalCutHi())) continue ;

 // Pt & Eta - this is done differently for different Harmonic & Selection
  for(int selN = 1; selN < AliFlowConstants::kSels; selN++)               // not even consider the 0th selection (no cut applied there)
  {
  // min. TPC hits required
   if(pFlowSelect->NhitsCut(selN) && (pFlowTrack->FitPtsTPC()<pFlowSelect->NhitsCut(selN))) continue ;

   for(int harN = 0; harN < AliFlowConstants::kHars; harN++) 
   {
   // Eta - gets rid of the track with Eta outside the range
    if((pFlowSelect->EtaCutHi(harN%2,selN)>pFlowSelect->EtaCutLo(harN%2,selN)) && (TMath::Abs(eta)<pFlowSelect->EtaCutLo(harN%2,selN) || TMath::Abs(eta)>pFlowSelect->EtaCutHi(harN%2,selN))) continue ; 
   // Pt - gets rid of the track with Pt outside the range
    if((pFlowSelect->PtCutHi(harN%2,selN)>pFlowSelect->PtCutLo(harN%2,selN)) && (pt<pFlowSelect->PtCutLo(harN%2,selN) || pt>pFlowSelect->PtCutHi(harN%2,selN))) continue ; 
  
    pFlowTrack->SetSelect(harN, selN) ;  // if cuts defined (low<high) && track is in the range -> Set [har][sel] Flag ON

    if(AliFlowConstants::fgDebug) 
    {
     if(harN==1)
     {
      cout << " n. " << itr << " , pFlowTrack->PrintSelection() = " ; pFlowTrack->PrintSelection() ;
      if(pFlowSelect->Pid()[0] != '\0') { cout << " track:  pid  " << pFlowTrack->Pid() << " = "<< pFlowSelect->Pid() << endl ; } 
      cout << " track:  dca  " << pFlowSelect->DcaGlobalCutLo() << " < " << gDca << " < " << pFlowSelect->DcaGlobalCutHi() << endl ;
      cout << " track:  eta  " << pFlowSelect->EtaCutLo(harN,selN) << " < |" << eta << "| < " << pFlowSelect->EtaCutHi(harN,selN) << endl ;
      cout << " track:  pT   " << pFlowSelect->PtCutLo(harN,selN) << " < " << pt << " < " << pFlowSelect->PtCutHi(harN,selN) << endl ;
      cout << " pFlowTrack->PrintSelection() = " ; pFlowTrack->PrintSelection() ;
     }
     // cout << " selN " << selN << " ,  harN " << harN%2 << " - si" << endl ;
    }
   }
  }
 }
}
//-------------------------------------------------------------
void AliFlowEvent::SetPids()
{
 // Re-sets the tracks P.id. (using the current AliFlowConstants::fgBayesian[] array).
 // To use the Raw P.Id (just the detector response function), set fgBayesian[] = {1,1,1,1,1,1}.
 
 const Int_t kCode[] = {11,13,211,321,2212,10010020} ;
 for(Int_t itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;

  Float_t bayPid[AliFlowConstants::kPid] ;
  if(fgCustomRespFunc)  { pFlowTrack->PidProbsC(bayPid) ; }
  else 			{ pFlowTrack->PidProbs(bayPid)  ; }

  Int_t maxN = 2 ; 		  // if No id. -> then is a Pi
  Float_t pidMax = bayPid[2] ;    // (if all equal, Pi probability was the first)
  for(Int_t nP=0;nP<AliFlowConstants::kPid;nP++)
  {
   if(bayPid[nP]>pidMax) { maxN = nP ; pidMax = bayPid[nP] ; }
  }
  Int_t pdgCode = TMath::Sign(kCode[maxN],pFlowTrack->Charge()) ;     
  pFlowTrack->SetMostLikelihoodPID(pdgCode) ;			
 }
}
//-------------------------------------------------------------
void AliFlowEvent::MakeSubEvents() const
{
 // Make random or eta sub-events
 
 if(fgEtaSubs == 1)        { MakeEtaSubEvents() ; }
 else if(fgEtaSubs == -1)  { MakeChrSubEvents() ; }
 else if(fgEtaSubs == 0)   { MakeRndSubEvents() ; }
} 
//-------------------------------------------------------------
void AliFlowEvent::MakeRndSubEvents() const 
{
 // Make random subevents
 
 int eventMult[AliFlowConstants::kHars][AliFlowConstants::kSels] = {{0}};
 int harN, selN, subN = 0;
 
 // loop to count the total number of tracks for each selection
 for(Int_t itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  for (selN = 0; selN < AliFlowConstants::kSels; selN++) 
  {
   for (harN = 0; harN < AliFlowConstants::kHars; harN++) 
   {
    if(pFlowTrack->Select(harN, selN)) { eventMult[harN][selN]++ ; }
   }
  }
 }
 // loop to set the SubEvent member
 for (selN = 0; selN < AliFlowConstants::kSels; selN++) 
 {
  for (harN = 0; harN < AliFlowConstants::kHars; harN++) 
  {
   int subEventMult = eventMult[harN][selN] / AliFlowConstants::kSubs;
   if (subEventMult) 
   {
    subN = 0;
    int countN = 0;
    for(Int_t itr=0;itr<TrackCollection()->GetEntries();itr++) 
    {
     AliFlowTrack* pFlowTrack ;
     pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
     if(pFlowTrack->Select(harN, selN)) 
     {
      pFlowTrack->SetSubevent(harN, selN, subN);
      countN++;
      if (countN % subEventMult == 0.) { subN++ ; }
     }
    }
   }
  }
 }
 return ; 
}
//-------------------------------------------------------------
void AliFlowEvent::MakeEtaSubEvents() const
{
 // Make subevents for positive and negative eta 

 int harN, selN = 0;
 // loop to set the SubEvent member
 for (selN = 0; selN < AliFlowConstants::kSels; selN++) 
 {
  for (harN = 0; harN < AliFlowConstants::kHars; harN++) 
  {
   for(Int_t itr=0;itr<TrackCollection()->GetEntries();itr++) 
   {
    AliFlowTrack* pFlowTrack ;
    pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
    if(pFlowTrack->Select(harN, selN)) 
    {
     float eta = pFlowTrack->Eta();
     if (eta > 0.) { pFlowTrack->SetSubevent(harN, selN, 0) ; } 
     else          { pFlowTrack->SetSubevent(harN, selN, 1) ; }
    }
   }
  }
 }
}
//-------------------------------------------------------------
void AliFlowEvent::MakeChrSubEvents() const
{
 // Make subevents for positive and negative charged tracks
 
 int harN, selN = 0;
 // loop to set the SubEvent member
 for (selN = 0; selN < AliFlowConstants::kSels; selN++) 
 {
  for (harN = 0; harN < AliFlowConstants::kHars; harN++) 
  {
   for(Int_t itr=0;itr<TrackCollection()->GetEntries();itr++) 
   {
    AliFlowTrack* pFlowTrack ;
    pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
    if(pFlowTrack->Select(harN, selN)) 
    {
     float charge = pFlowTrack->Charge();
     if (charge > 0.) { pFlowTrack->SetSubevent(harN, selN, 0) ; } 
     else             { pFlowTrack->SetSubevent(harN, selN, 1) ; }
    }
   }
  }
 }
}
//-------------------------------------------------------------
void AliFlowEvent::RandomShuffle() 
{
 // Randomly re-shuffles the tracks in the array; if a track is not
 // primary, the reference carried by the reconstructed mother (in 
 // the v0 array) is changed accordingly.

 return ; // ...at the moment is disabled ... // 
 
//  Int_t tot = 0 ;
//  UInt_t imax = fTrackCollection->GetEntries() ;
//  TRandom* rnd = new TRandom(0) ; // SetSeed(0) ;
//  TClonesArray* newTrackCollection = new TClonesArray("AliFlowTrack",imax) ;
// 
//  // re-arranges the ObjArray (TrackCollection())
//  for(UInt_t itr=0;itr<imax;itr++)
//  {
//   AliFlowTrack* pFlowTrack ;
//   pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
//   
//   UInt_t rndNumber = rnd->Integer(imax) ;
//   Bool_t put = kFALSE ;
//   while(!put)
//   { 
//    if(!newTrackCollection->AddrAt(rndNumber)) 
//    { 
//     ... new(newTrackCollection[rndNumber]) AliFlowTrack(*pFlowTrack)  ; 
//     put = kTRUE ; tot++ ;  		
//     if(AliFlowConstants::fgDebug) { cout << "  " << itr << " --> " << rndNumber << endl ; } 
//    }
//    else 
//    {
//     rndNumber++ ; if(rndNumber>=imax) { rndNumber -= imax ; }
//    }
//   }
//  }
//  if(AliFlowConstants::fgDebug) { cout << "* RandomShuffle() :  " << tot << "/" << imax << " flow tracks have been shuffled " << endl ; }  
//  fTrackCollection = newTrackCollection ;
}
//-----------------------------------------------------------------------
UInt_t AliFlowEvent::Centrality() 
{
 // Returns the Centrality Class as stored

 if(fCentrality < 0)  { SetCentrality() ; } 
 return fCentrality ;
}
//-----------------------------------------------------------------------
void AliFlowEvent::SetCentrality() 
{
 // Sets the Centrality Classes basing on Multiplicity at mid rapidity, 
 // ... (ZDC information can be added) .
 
 Float_t* cent ;
 Int_t  tracks = MultEta() ;

 if(RunID() == -1)
 { 
  cent = AliFlowConstants::fgCentNorm ;
  //if centrality classes are not defined, does it now (with CentNorm & MaxMult)
  if(cent[AliFlowConstants::kCents-1] <= 1) 
  {
   for(Int_t ic=0;ic<AliFlowConstants::kCents;ic++)
   {
    cent[ic] *= AliFlowConstants::fgMaxMult ; 
    if(AliFlowConstants::fgDebug) { cout << "Centrality[" << ic << "] = " << cent[ic] << " . " << endl ; }
   }
  }
 }
 else if((RunID() != -1) && (CenterOfMassEnergy() == 5500.))
 {
  cent = (Float_t*)AliFlowConstants::fgCent0 ;
 } 
 else // other definition of centrality are possible ...
 {
  cent = (Float_t*)AliFlowConstants::fgCent0 ;
 } 
 if      (tracks < cent[0])  { fCentrality = 0; }
 else if (tracks < cent[1])  { fCentrality = 1; }
 else if (tracks < cent[2])  { fCentrality = 2; }
 else if (tracks < cent[3])  { fCentrality = 3; }
 else if (tracks < cent[4])  { fCentrality = 4; }
 else if (tracks < cent[5])  { fCentrality = 5; }
 else if (tracks < cent[6])  { fCentrality = 6; }
 else if (tracks < cent[7])  { fCentrality = 7; }
 else if (tracks < cent[8])  { fCentrality = 8; }
 else                        { fCentrality = 9; }

 if(AliFlowConstants::fgDebug) { cout << " * Centrality Class :  " << fCentrality << " . " << endl ; }
}
//-----------------------------------------------------------------------
TVector AliFlowEvent::Bayesian() const 
{ 
 // Returns bayesian array of particle abundances (from AliFlowConstants::)
 
 TVector bayes(AliFlowConstants::kPid) ; 
 for(Int_t i=0;i<AliFlowConstants::kPid;i++) { bayes[i] = AliFlowConstants::fgBayesian[i] ; }
 return bayes ;
}
//-----------------------------------------------------------------------
void AliFlowEvent::PrintFlagList() const 
{
 // Prints the list of selection cuts ( called in AliFlowInterface::Finish() )
 
 cout << "#######################################################" << endl;
 cout << "# Weighting and Striping:" << endl;
 if(PtWgt()) 
 {
  cout << "#	PtWgt =  kTRUE " << endl ; 	// (also for output of PhiWgt file?)
  cout << "#	PtWgt Saturation =  " << AliFlowConstants::fgPtWgtSaturation << endl;
 } 
 else 
 {
  cout << "#	PtWgt = kFALSE" << endl;
 }
 if(EtaWgt()) 
 {
  cout << "#	EtaWgt = kTRUE " << endl ; 	// (also for output of PhiWgt file for odd harmonics?)
 } 
 else 
 {
  cout << "#	EtaWgt = kFALSE" << endl;
 }
 cout << "#######################################################" << endl;  
}
//-----------------------------------------------------------------------
Int_t AliFlowEvent::MultEta() const
{
 // Returns the multiplicity in the interval |eta|<(AliFlowConstants::fgEetaMid), used 
 // for centrality measurement (see centrality classes in fCentrality) .
 
 Int_t goodtracks = 0 ;
 for(Int_t itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if((pFlowTrack->Charge()) && (TMath::Abs(pFlowTrack->Eta())<AliFlowConstants::fgEtaMid)) { goodtracks++ ; }
 }
 return goodtracks ; 
}
//-----------------------------------------------------------------------
Int_t AliFlowEvent::UncorrNegMult(Float_t eta)  const
{ 
 // Negative multiplicity in the interval (-eta..eta)
 // (default is  AliFlowConstants::fgEetaGood = 0.9)
 
 Int_t negMult = 0 ;
 for(Int_t itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if((pFlowTrack->Charge()<0) && (TMath::Abs(pFlowTrack->Eta())<TMath::Abs(eta))) { negMult++ ; }
  //delete pFlowTrack ;
 }
 return negMult; 
}
//-----------------------------------------------------------------------
Int_t AliFlowEvent::UncorrPosMult(Float_t eta)  const
{ 
 // Positive multiplicity in the interval (-eta..eta)
 // (default is  AliFlowConstants::fgEetaGood = 0.9)
 
 Int_t posMult = 0 ;
 for(Int_t itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if((pFlowTrack->Charge()>0) && (TMath::Abs(pFlowTrack->Eta())<TMath::Abs(eta))) { posMult++ ; }
  //delete pFlowTrack ;
 }
 return posMult; 
}
//-----------------------------------------------------------------------
TVector3 AliFlowEvent::VertexPos() const
{
 // Returns primary vertex position as a TVector3

 Float_t vertex[3] ;
 VertexPos(vertex) ;
 return TVector3(vertex) ;
}
//-----------------------------------------------------------------------
void AliFlowEvent::MakeAll()
{ 
 // calculates all quantities in 1 shoot 

 Double_t mQx[AliFlowConstants::kSels][AliFlowConstants::kHars] ;
 Double_t mQy[AliFlowConstants::kSels][AliFlowConstants::kHars] ;
 Double_t mQxSub[AliFlowConstants::kSubs][AliFlowConstants::kSels][AliFlowConstants::kHars] ;
 Double_t mQySub[AliFlowConstants::kSubs][AliFlowConstants::kSels][AliFlowConstants::kHars] ;
 // -
 int selN, harN, subN ;
 for(selN=0;selN<AliFlowConstants::kSels;selN++) 
 {
  for(harN=0;harN<AliFlowConstants::kHars;harN++) 
  {
   mQx[selN][harN]    = 0. ;  	
   mQy[selN][harN]    = 0. ;  	
   fMult[selN][harN]  = 0 ;	
   fSumOfWeightSqr[selN][harN] = 0. ;
   for(subN=0;subN<AliFlowConstants::kSubs;subN++)
   {
    mQxSub[subN][selN][harN]   = 0. ;
    mQySub[subN][selN][harN]   = 0. ;
    fMultSub[subN][selN][harN] = 0 ;
   }
  }
 }
 
 double order = 0.  ;
 double phiWgt = 0. ;  
 float  phi = 0.    ; 		     
 // -
 int itr ;
 for(itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  phi = pFlowTrack->Phi();
  for(selN=0;selN<AliFlowConstants::kSels;selN++) 
  {
   for(harN=0;harN<AliFlowConstants::kHars;harN++) 
   {
    order = (double)(harN+1) ;
    if(pFlowTrack->Select(harN,selN)) 
    {
     phiWgt = PhiWeight(selN,harN,pFlowTrack) ;  
     fSumOfWeightSqr[selN][harN] += phiWgt*phiWgt ;
     mQx[selN][harN] += phiWgt * cos(phi * order) ;
     mQy[selN][harN] += phiWgt * sin(phi * order) ;
     fMult[selN][harN]++ ;	
     for(subN=0;subN<AliFlowConstants::kSubs;subN++)
     {
      if(pFlowTrack->Select(harN,selN,subN))
      {
       mQxSub[subN][selN][harN] += phiWgt * cos(phi * order) ;
       mQySub[subN][selN][harN] += phiWgt * sin(phi * order) ;
       fMultSub[subN][selN][harN]++  ;
      }
     }  // sub
    }  // if
   }  // har
  }  // sel
 }  //itr

 for(selN=0;selN<AliFlowConstants::kSels;selN++)  
 {
  for(harN=0;harN<AliFlowConstants::kHars;harN++)
  {
   fQ[selN][harN].Set(mQx[selN][harN],mQy[selN][harN]) ;  
   for(subN=0;subN<AliFlowConstants::kSubs;subN++) { fQSub[subN][selN][harN].Set(mQxSub[subN][selN][harN],mQySub[subN][selN][harN]) ; }
  }
 }

 fDone = kTRUE ;
}
//-----------------------------------------------------------------------
