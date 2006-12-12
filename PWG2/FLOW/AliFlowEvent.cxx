//////////////////////////////////////////////////////////////////////
//
// $Id$
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
#include "AliFlowConstants.h"
#include "TVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include <iostream>
using namespace std; //required for resolving the 'cout' symbol

// - Flags & Sets
Bool_t  AliFlowEvent::fPtWgt	       = kFALSE ;  // gives pT-based weights
Bool_t  AliFlowEvent::fEtaWgt	       = kFALSE ;  // gives eta-based weights
Bool_t  AliFlowEvent::fOnePhiWgt       = kTRUE  ;  // kTRUE --> ENABLEs SINGLE WEIGHT ISTOGRAM , kFALSE --> ENABLEs 3 WEIGHT ISTOGRAMS
Bool_t  AliFlowEvent::fNoWgt	       = kFALSE ;  // No Weight is used
// - Eta Sub-Events (later used to calculate the resolution)
Bool_t  AliFlowEvent::fEtaSubs 	       = kFALSE ;  // makes eta subevents

ClassImp(AliFlowEvent) ;
//-----------------------------------------------------------
AliFlowEvent::AliFlowEvent(Int_t lenght)  
{
 // Default constructor: initializes the ObjArray of FlowTracks and FlowV0s, 
 // cleans the internal variables, sets all the weights to 1, sets default flags.
 
 fEventID = 0 ; 
 fRunID = 0 ;		       
 fOrigMult = 0 ;
 fL0TriggerWord = 0 ;			       
 fCentrality = -1 ;
 fZDCpart = 0 ;
 for(int zz=0;zz<3;zz++) { fZDCenergy[zz] = 0. ; }
 
 // Make a new track collection
 fTrackCollection =  new TObjArray(lenght) ;
 fV0Collection =  new TObjArray(0) ;
 
 // Set Weights Arrays to 1 (default)
 for(int nS=0;nS<Flow::nSels;nS++)
 {
  for(int nH=0;nH<Flow::nHars;nH++) 
  {
   for(int nP=0;nP<Flow::nPhiBins;nP++) 
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
 //for(int nH=0;nH<Flow::nHars;nH++) { fExtPsi[nH] = 0. ; fExtRes[nH] = 0. ; }
 
 // The Expected particles abundance is taken directly from Flow::fBayesian[] (see Bayesian P.Id.)
 
 fDone = kFALSE ;
}
//-----------------------------------------------------------
AliFlowEvent::~AliFlowEvent() 
{
 // Default distructor: deletes the ObjArrays. 
 
 fTrackCollection->Delete() ; delete fTrackCollection ;
 fV0Collection->Delete()    ; delete fV0Collection ;
}
//-------------------------------------------------------------

//-------------------------------------------------------------
Double_t AliFlowEvent::PhiWeightRaw(Int_t selN, Int_t harN, AliFlowTrack* pFlowTrack) const  
{
 // Weight for making the event plane isotropic in the lab.

 float phi = pFlowTrack->Phi() ; if(phi < 0.) { phi += 2*TMath::Pi() ; }
 Double_t  eta = (Double_t)pFlowTrack->Eta() ;
 int n = (int)((phi/(2*TMath::Pi()))*Flow::nPhiBins);

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

 // ... weight can be implemented as well for ITS , TRD , PHOS , etc.
 // 
 // if(...) 
 // {
 //  n = (int)((phi/2*TMath::Pi())*Flow::nPhiBins...) ;
 //  ...
 // }

 return phiWgt ;
}
//-------------------------------------------------------------
Double_t AliFlowEvent::Weight(Int_t selN, Int_t harN, AliFlowTrack* pFlowTrack) const  
{ 
 // Weight for enhancing the resolution (eta gives sign +/- for Odd Harmonics)

 bool oddHar = (harN+1) % 2 ;
 Double_t phiWgt = 1. ;
 if(PtWgt()) 
 {
  Double_t pt = (Double_t)pFlowTrack->Pt();
  if(pt < Flow::fPtWgtSaturation) { phiWgt *= pt ; } 
  else 			   	  { phiWgt *= Flow::fPtWgtSaturation ; } // pt weighting going constant
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
 // (it simply rerurns PhiWeightRaw() * Weight()). If fNoWgt = kTRUE, returns +/-1 ,
 // basing on Sign(eta), for odd harmonics .
 
 if(NoWgt()) // no weights (but +/- according to eta)
 { 
  bool oddHar = (harN+1) % 2 ;
  if(oddHar) { return TMath::Sign((Double_t)1.,(Double_t)pFlowTrack->Eta()) ; }
  else       { return (Double_t)1. ; }
 } 
 Double_t phiWgtRaw = PhiWeightRaw(selN, harN, pFlowTrack);
 Double_t weight = Weight(selN, harN, pFlowTrack);
 if(Flow::fDebug) { cout << "[PhiWeight]: phiWgtRaw = " << phiWgtRaw << " , weight = " << weight << " , eta = " << pFlowTrack->Eta() << endl ; }

 return phiWgtRaw * weight;
}
//-------------------------------------------------------------
Int_t AliFlowEvent::Mult(AliFlowSelection* pFlowSelect)  	
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
Float_t AliFlowEvent::MeanPt(AliFlowSelection* pFlowSelect)
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
TVector2 AliFlowEvent::Q(AliFlowSelection* pFlowSelect)
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
   if(Flow::fDebug) { cout << itr << " phi = " << phi << " ,  wgt = " << phiWgt << endl ; }
  }
 }
 mQ.Set(mQx, mQy);

 return mQ;
}
//-------------------------------------------------------------
Float_t AliFlowEvent::Psi(AliFlowSelection* pFlowSelect)
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
TVector2 AliFlowEvent::NormQ(AliFlowSelection* pFlowSelect)  	
{
 // Normalized Q = Q/sqrt(sum of weights^2) for the specified Selection 

 if(fDone) 
 {
  TVector2 mQ = fQ[pFlowSelect->Sel()][pFlowSelect->Har()] ;   
  double SumOfWeightSqr = fSumOfWeightSqr[pFlowSelect->Sel()][pFlowSelect->Har()] ;
  if(SumOfWeightSqr) { mQ /= TMath::Sqrt(SumOfWeightSqr) ; }
  else { mQ.Set(0.,0.) ; }
  return mQ ;
 }
 // -
 TVector2 mQ ;
 Double_t mQx=0. , mQy=0. ;
 double SumOfWeightSqr = 0 ;
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
   SumOfWeightSqr += phiWgt*phiWgt;

   float phi = pFlowTrack->Phi();
   mQx += phiWgt * cos(phi * order);
   mQy += phiWgt * sin(phi * order);
  }
 }
 if(SumOfWeightSqr) { mQ.Set(mQx/TMath::Sqrt(SumOfWeightSqr), mQy/TMath::Sqrt(SumOfWeightSqr)); }
 else { mQ.Set(0.,0.); }
  
 return mQ;
}
//-------------------------------------------------------------
Float_t AliFlowEvent::q(AliFlowSelection* pFlowSelect) 
{ 
 // Magnitude of normalized Q vector (without pt or eta weighting) for the specified Selection 

 TVector2 mQ ;
 Double_t mQx = 0. , mQy = 0. ;
 int selN     = pFlowSelect->Sel() ;
 int harN     = pFlowSelect->Har() ;
 double order = (double)(harN + 1) ;
 double SumOfWeightSqr = 0 ;

 Int_t itr ;
 for(itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if(pFlowSelect->Select(pFlowTrack)) 
  {
   double phiWgt = PhiWeightRaw(selN, harN, pFlowTrack); // Raw
   SumOfWeightSqr += phiWgt*phiWgt;
   float phi = pFlowTrack->Phi();
   mQx += phiWgt * cos(phi * order);
   mQy += phiWgt * sin(phi * order);
  }
 }
 if(SumOfWeightSqr) { mQ.Set(mQx/TMath::Sqrt(SumOfWeightSqr), mQy/TMath::Sqrt(SumOfWeightSqr)); }
 else { mQ.Set(0.,0.); }
  
 return mQ.Mod();
}
//-----------------------------------------------------------------------
Double_t AliFlowEvent::G_New(AliFlowSelection* pFlowSelect, Double_t Zx, Double_t Zy)
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
Double_t AliFlowEvent::G_Old(AliFlowSelection* pFlowSelect, Double_t Zx, Double_t Zy) 
{ 
 // Generating function for the old cumulant method (if expanded in Taylor
 // series, one recovers G_New() in new new cumulant method)

 TVector2 normQ = NormQ(pFlowSelect);

 return exp(2*Zx*normQ.X() + 2*Zy*normQ.Y());
}
//-----------------------------------------------------------------------
Double_t AliFlowEvent::SumWeightSquare(AliFlowSelection* pFlowSelect)
{
 // Return sum of weights^2 for the specified Selection (used for normalization)

 if(fDone) 
 { 
  return fSumOfWeightSqr[pFlowSelect->Sel()][pFlowSelect->Har()] ; 
 }
 // -
 int selN = pFlowSelect->Sel();
 int harN = pFlowSelect->Har();
 Double_t SumOfWeightSqr = 0;

 Int_t itr ;
 for(itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if (pFlowSelect->Select(pFlowTrack)) 
  {
   double phiWgt = PhiWeight(selN, harN, pFlowTrack);
   SumOfWeightSqr += phiWgt*phiWgt;
  }
 }
 if(SumOfWeightSqr < 0.) { return Mult(pFlowSelect) ; }

 return SumOfWeightSqr;
}
//-------------------------------------------------------------
Double_t AliFlowEvent::WgtMult_q4(AliFlowSelection* pFlowSelect)
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
Double_t AliFlowEvent::WgtMult_q6(AliFlowSelection* pFlowSelect)
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
  for(int harN=0;harN<Flow::nHars;harN++) { pFlowTrack->SetSelect(harN,0) ; } 

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
  float  Pt  = pFlowTrack->Pt();   		 
  float  gDca = pFlowTrack->Dca() ; 			

 // Global DCA - gets rid of the track with DCA outside the range
  if((pFlowSelect->DcaGlobalCutHi()>pFlowSelect->DcaGlobalCutLo()) && (gDca<pFlowSelect->DcaGlobalCutLo() || gDca>pFlowSelect->DcaGlobalCutHi())) continue ;

 // Pt & Eta - this is done differently for different Harmonic & Selection
  for(int selN = 1; selN < Flow::nSels; selN++)               // not even consider the 0th selection (no cut applied there)
  {
  // min. TPC hits required
   if(pFlowSelect->NhitsCut(selN) && (pFlowTrack->FitPtsTPC()<pFlowSelect->NhitsCut(selN))) continue ;

   for(int harN = 0; harN < Flow::nHars; harN++) 
   {
   // Eta - gets rid of the track with Eta outside the range
    if((pFlowSelect->EtaCutHi(harN%2,selN)>pFlowSelect->EtaCutLo(harN%2,selN)) && (TMath::Abs(eta)<pFlowSelect->EtaCutLo(harN%2,selN) || TMath::Abs(eta)>pFlowSelect->EtaCutHi(harN%2,selN))) continue ; 
   // Pt - gets rid of the track with Pt outside the range
    if((pFlowSelect->PtCutHi(harN%2,selN)>pFlowSelect->PtCutLo(harN%2,selN)) && (Pt<pFlowSelect->PtCutLo(harN%2,selN) || Pt>pFlowSelect->PtCutHi(harN%2,selN))) continue ; 
  
    pFlowTrack->SetSelect(harN, selN) ;  // if cuts defined (low<high) && track is in the range -> Set [har][sel] Flag ON

    if(Flow::fDebug) 
    {
     cout << " harN " << harN%2 << " ,  selN " << selN << " - si" << endl ;
     if(pFlowSelect->Pid()[0] != '\0') { cout << " track:  pid  " << pFlowTrack->Pid() << " = "<< pFlowSelect->Pid() << endl ; } 
     cout << " track:  dca  " << pFlowSelect->DcaGlobalCutLo() << " < " << gDca << " < " << pFlowSelect->DcaGlobalCutHi() << endl ;
     cout << " track:  eta  " << pFlowSelect->EtaCutLo(harN,selN) << " < |" << eta << "| < " << pFlowSelect->EtaCutHi(harN,selN) << endl ;
     cout << " track:  pT   " << pFlowSelect->PtCutLo(harN,selN) << " < " << Pt << " < " << pFlowSelect->PtCutHi(harN,selN) << endl ;
     pFlowTrack->PrintSelection() ;
    }
   }
  }
 }
}
//-------------------------------------------------------------
void AliFlowEvent::SetPids()
{
 // Re-sets the tracks P.id. (using the current Flow::fBayesian[] array)
 
 const Int_t code[] = {11,13,211,321,2212,10010020} ;
 for(Int_t itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  TVector bayPid = pFlowTrack->PidProbs() ;
  Int_t maxN = 2 ; 		   // if No id. -> then is a Pi
  Float_t pid_max = bayPid[2] ;    // (if all equal, Pi probability get's the advantage to be the first)
  for(Int_t nP=0;nP<Flow::nPid;nP++)
  {
   if(bayPid[nP]>pid_max) { maxN = nP ; pid_max = bayPid[nP] ; }
  }
  Int_t pdg_code = TMath::Sign(code[maxN],pFlowTrack->Charge()) ;     
  pFlowTrack->SetMostLikelihoodPID(pdg_code);			
 }
}
//-------------------------------------------------------------
void AliFlowEvent::MakeSubEvents() 
{
 // Make random or eta sub-events
 
 if(EtaSubs())  { MakeEtaSubEvents() ; }
 else 		{ MakeRndSubEvents() ; }
} 
//-------------------------------------------------------------
void AliFlowEvent::MakeRndSubEvents() 
{
 // Make random subevents
 
 int eventMult[Flow::nHars][Flow::nSels] = {{0}};
 int harN, selN, subN = 0;
 
 // loop to count the total number of tracks for each selection
 for(Int_t itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  for (selN = 0; selN < Flow::nSels; selN++) 
  {
   for (harN = 0; harN < Flow::nHars; harN++) 
   {
    if(pFlowTrack->Select(harN, selN)) { eventMult[harN][selN]++ ; }
   }
  }
 }
 // loop to set the SubEvent member
 for (selN = 0; selN < Flow::nSels; selN++) 
 {
  for (harN = 0; harN < Flow::nHars; harN++) 
  {
   int subEventMult = eventMult[harN][selN] / Flow::nSubs;
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
void AliFlowEvent::MakeEtaSubEvents() 
{
 // Make subevents for positive and negative eta 
 // (when done, fEtaSubs flag setted to kTRUE).
 
 int harN, selN = 0;
 // loop to set the SubEvent member
 for (selN = 0; selN < Flow::nSels; selN++) 
 {
  for (harN = 0; harN < Flow::nHars; harN++) 
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
void AliFlowEvent::RandomShuffle() 
{
 // Randomly re-shuffles the tracks in the array; if a track is not
 // primary, the reference carried by the reconstructed mother (in 
 // the v0 array) is changed accordingly.

 Int_t tot = 0 ;
 UInt_t imax = TrackCollection()->GetEntries() ;
 TRandom* rnd = new TRandom(0) ; // SetSeed(0) ;
 TObjArray* newTrackCollection = new TObjArray(imax) ;

 // re-arranges the ObjArray (TrackCollection())
 for(UInt_t itr=0;itr<imax;itr++)
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  
  UInt_t rndNumber = rnd->Integer(imax) ;
  Bool_t put = kFALSE ;
  while(!put)
  { 
   if(!newTrackCollection->At(rndNumber)) 
   { 
    newTrackCollection->AddAt(pFlowTrack, rndNumber) ; 
    put = kTRUE ; tot++ ;  		if(Flow::fDebug) { cout << "  " << itr << " --> " << rndNumber << endl ; } 
   }
   else 
   {
    rndNumber++ ; if(rndNumber>=imax) { rndNumber -= imax ; }
   }
  }
 }
 if(Flow::fDebug) { cout << "* RandomShuffle() :  " << tot << "/" << imax << " flow tracks have been shuffled " << endl ; }  
 fTrackCollection = newTrackCollection ;
}
//-----------------------------------------------------------------------
UInt_t AliFlowEvent::Centrality() 
{
 // Returns the Centrality Class as stored

 if(fCentrality < 0)  { SetCentrality() ; } 
 return fCentrality ;
}
//-----------------------------------------------------------------------
void AliFlowEvent::SetCentrality(Int_t cent) 
{
 // Set the Centrality Classes to "cent" .

 fCentrality = cent ; 
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
  cent = Flow::fCentNorm ;
  //if centrality classes are not defined, does it now (with CentNorm & MaxMult)
  if(cent[Flow::nCents-1] <= 1) 
  {
   for(Int_t ic=0;ic<Flow::nCents;ic++)
   {
    cent[ic] *= Flow::fMaxMult ; 
    if(Flow::fDebug) { cout << "Centrality[" << ic << "] = " << cent[ic] << " . " << endl ; }
   }
  }
 }
 else if((RunID() != -1) && (CenterOfMassEnergy() == 5500.))
 {
  cent = (Float_t*)Flow::fCent0 ;
 } 
 else // other definition of centrality are possible...
 {
  cent = (Float_t*)Flow::fCent0 ;
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

 if(Flow::fDebug) { cout << " * Centrality Class :  " << fCentrality << " . " << endl ; }
}
//-----------------------------------------------------------------------
void AliFlowEvent::Bayesian(Double_t bayes[Flow::nPid]) 
{
 // Returns bayesian array of particles' abundances (from Flow::)
 
 for(Int_t i=0;i<Flow::nPid;i++) { bayes[i] = Flow::fBayesian[i] ; }
}
//-----------------------------------------------------------------------
TVector AliFlowEvent::Bayesian() 
{ 
 TVector bayes(Flow::nPid) ; 
 for(Int_t i=0;i<Flow::nPid;i++) { bayes[i] = Flow::fBayesian[i] ; }
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
  cout << "#	PtWgt Saturation =  " << Flow::fPtWgtSaturation << endl;
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
void AliFlowEvent::SetEventID(const Int_t& id) 		  
{ 
 // Sets Event ID and the Event name (name = evtNumber_runId)
 
 fEventID = id ; 
 TString name = "" ;
 name += fEventID ;
 if(fRunID) { name += "_" ; name += fRunID ; }
 SetName(name) ;  // from TNamed::SetName
}
//-----------------------------------------------------------------------
Int_t AliFlowEvent::MultEta()
{
 // Returns the multiplicity in the interval |eta|<(Flow::fEetaMid), used 
 // for centrality measurement (see centrality classes in fCentrality) .
 
 Int_t goodtracks = 0 ;
 for(Int_t itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if((pFlowTrack->Charge()) && (TMath::Abs(pFlowTrack->Eta())<Flow::fEtaMid)) { goodtracks++ ; }
 }
 return goodtracks ; 
}
//-----------------------------------------------------------------------
Int_t AliFlowEvent::UncorrNegMult(Float_t eta)  const
{ 
 // Negative multiplicity in the interval (-eta..eta)
 // (default is  Flow::fEetaGood = 0.9)
 
 Int_t negMult = 0 ;
 for(Int_t itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if((pFlowTrack->Charge()<0) && (TMath::Abs(pFlowTrack->Eta())<TMath::Abs(eta))) { negMult++ ; }
  delete pFlowTrack ;
 }
 return negMult; 
}
//-----------------------------------------------------------------------
Int_t AliFlowEvent::UncorrPosMult(Float_t eta)  const
{ 
 // Positive multiplicity in the interval (-eta..eta)
 // (default is  Flow::fEetaGood = 0.9)
 
 Int_t posMult ;
 for(Int_t itr=0;itr<TrackCollection()->GetEntries();itr++) 
 {
  AliFlowTrack* pFlowTrack ;
  pFlowTrack = (AliFlowTrack*)TrackCollection()->At(itr) ;
  if((pFlowTrack->Charge()>0) && (TMath::Abs(pFlowTrack->Eta())<TMath::Abs(eta))) { posMult++ ; }
  delete pFlowTrack ;
 }
 return posMult; 
}
//-----------------------------------------------------------------------
void AliFlowEvent::VertexPos(Float_t vtx[3])  const	 	
{ 
 for(Int_t ii=0;ii<3;ii++) { vtx[ii] = fVertexPos[ii] ; }
}
//-----------------------------------------------------------------------
TVector3 AliFlowEvent::VertexPos() const
{
 Float_t vertex[3] ;
 VertexPos(vertex) ;
 return TVector3(vertex) ;
}
//-----------------------------------------------------------------------
void AliFlowEvent::SetVertexPos(Float_t v1,Float_t v2,Float_t v3)
{ 
 fVertexPos[0] = v1 ; fVertexPos[1] = v2 ; fVertexPos[2] = v3 ; 
}
//-----------------------------------------------------------------------
void AliFlowEvent::MakeAll()
{ 
 // calculates all quantities in 1 shoot ...
 //  ...

 Double_t mQx[Flow::nSels][Flow::nHars] ;
 Double_t mQy[Flow::nSels][Flow::nHars] ;
 Double_t mQxSub[Flow::nSubs][Flow::nSels][Flow::nHars] ;
 Double_t mQySub[Flow::nSubs][Flow::nSels][Flow::nHars] ;
 // -
 int selN, harN, subN ;
 for(selN=0;selN<Flow::nSels;selN++) 
 {
  for(harN=0;harN<Flow::nHars;harN++) 
  {
   mQx[selN][harN]    = 0. ;  	
   mQy[selN][harN]    = 0. ;  	
   fMult[selN][harN]  = 0 ;	
   fSumOfWeightSqr[selN][harN] = 0. ;
   for(subN=0;subN<Flow::nSubs;subN++)
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
  for(selN=0;selN<Flow::nSels;selN++) 
  {
   for(harN=0;harN<Flow::nHars;harN++) 
   {
    order = (double)(harN+1) ;
    if(pFlowTrack->Select(harN,selN)) 
    {
     phiWgt = PhiWeight(selN,harN,pFlowTrack) ;  
     fSumOfWeightSqr[selN][harN] += phiWgt*phiWgt ;
     mQx[selN][harN] += phiWgt * cos(phi * order) ;
     mQy[selN][harN] += phiWgt * sin(phi * order) ;
     fMult[selN][harN]++ ;	
     for(subN=0;subN<Flow::nSubs;subN++)
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

 for(selN=0;selN<Flow::nSels;selN++)  
 {
  for(harN=0;harN<Flow::nHars;harN++)
  {
   fQ[selN][harN].Set(mQx[selN][harN],mQy[selN][harN]) ;  
   for(subN=0;subN<Flow::nSubs;subN++) { fQSub[subN][selN][harN].Set(mQxSub[subN][selN][harN],mQySub[subN][selN][harN]) ; }
  }
 }

 fDone = kTRUE ;
}
// //-----------------------------------------------------------------------
// Float_t  AliFlowEvent::ExtPsi(Int_t harN) const 	
// { 
//  // external R.P. angle (check input source...)
// 
//  if(harN<Flow::nHars) { return fExtPsi[harN] ; }
//  else 
//  { 
//   cout << "AliFlowEvent::ExtPsi(" << harN << ") : harmonic " << harN+1 << " is not there !" << endl ; 
//   return 0. ; 
//  }
// }
// //-----------------------------------------------------------------------
// Float_t  AliFlowEvent::ExtRes(Int_t harN) const
// { 
//  // external R.P. resolution (check input source...)
// 
//  if(harN<Flow::nHars) { return fExtRes[harN] ; }
//  else 
//  { 
//   cout << "AliFlowEvent::ExtRes(" << harN << ") : harmonic " << harN+1 << " is not there !" << endl ; 
//   return 0. ; 
//  }
// }
// //-----------------------------------------------------------------------
// void AliFlowEvent::SetExtPsi(Int_t harN,Float_t psi)  	 
// { 
//  if(harN<Flow::nHars) { fExtPsi[harN] = psi ; }
// }
// //-----------------------------------------------------------------------
// void AliFlowEvent::SetExtRes(Int_t harN,Float_t res)  	 
// { 
//  if(harN<Flow::nHars) { fExtRes[harN] = res ; }
// }
// //-----------------------------------------------------------------------



//-----------------------------------------------------------------------
// - those will go into the .h as inline functions ..................... 
//-----------------------------------------------------------------------
TObjArray* AliFlowEvent::TrackCollection() const 	{ return fTrackCollection; }
TObjArray* AliFlowEvent::V0Collection() const 		{ return fV0Collection; }
//-----------------------------------------------------------------------
Int_t	 AliFlowEvent::EventID() const  	 	{ return fEventID; }
Int_t	 AliFlowEvent::RunID() const		 	{ return fRunID; }
Double_t AliFlowEvent::CenterOfMassEnergy() const 	{ return Flow::fCenterOfMassEnergy ; }
Double_t AliFlowEvent::MagneticField() const 		{ return Flow::fMagneticField ; }
Short_t  AliFlowEvent::BeamMassNumberEast() const 	{ return Flow::fBeamMassNumberEast ; }
Short_t  AliFlowEvent::BeamMassNumberWest() const 	{ return Flow::fBeamMassNumberWest ; }
UInt_t   AliFlowEvent::OrigMult() const 	 	{ return fOrigMult; }
Long_t   AliFlowEvent::L0TriggerWord() const	 	{ return fL0TriggerWord; }
Int_t    AliFlowEvent::V0Mult() const 	  		{ return V0Collection()->GetEntries() ; }
Int_t    AliFlowEvent::FlowEventMult() const	 	{ return TrackCollection()->GetEntries() ; }
Int_t    AliFlowEvent::ZDCpart() const		 	{ return fZDCpart; }
Float_t  AliFlowEvent::ZDCenergy(Int_t npem) const	{ return fZDCenergy[npem]; }
Float_t  AliFlowEvent::PtWgtSaturation() const   	{ return Flow::fPtWgtSaturation; }
Bool_t   AliFlowEvent::PtWgt() const		 	{ return fPtWgt; }
Bool_t   AliFlowEvent::EtaWgt() const		 	{ return fEtaWgt; }
Bool_t   AliFlowEvent::FirstLastPhiWgt() const   	{ return !fOnePhiWgt ; }
Bool_t   AliFlowEvent::OnePhiWgt() const    		{ return fOnePhiWgt ; }
Bool_t   AliFlowEvent::NoWgt() const		 	{ return fNoWgt; }
Bool_t   AliFlowEvent::EtaSubs() const  	 	{ return fEtaSubs ; }
//-----------------------------------------------------------------------
void   	 AliFlowEvent::SetEtaSubs(Bool_t etasub)  		      { fEtaSubs = etasub ; }
void 	 AliFlowEvent::SetRunID(const Int_t& id)		      { fRunID = id; }
void 	 AliFlowEvent::SetMagneticField(const Double_t& mf)	      { Flow::fMagneticField = mf; }
void 	 AliFlowEvent::SetCenterOfMassEnergy(const Double_t& cms)     { Flow::fCenterOfMassEnergy = cms; }
void 	 AliFlowEvent::SetBeamMassNumberEast(const Short_t& bme)      { Flow::fBeamMassNumberEast = bme; }
void 	 AliFlowEvent::SetBeamMassNumberWest(const Short_t& bmw)      { Flow::fBeamMassNumberWest = bmw; }
void 	 AliFlowEvent::SetOrigMult(const UInt_t& tracks)	      { fOrigMult = tracks; }
void 	 AliFlowEvent::SetL0TriggerWord(const Long_t& trigger)        { fL0TriggerWord = trigger; }
void 	 AliFlowEvent::SetZDCpart(Int_t zdcp)			      { fZDCpart = zdcp ; }
void 	 AliFlowEvent::SetZDCenergy(Float_t n, Float_t p, Float_t em) { fZDCenergy[0] = n ; fZDCenergy[1] = p ; fZDCenergy[2] = em ; }
//-----------------------------------------------------------------------
void AliFlowEvent::SetBayesian(Double_t bayes[Flow::nPid]) 	  
{ 
 for(Int_t i=0;i<Flow::nPid;i++) { Flow::fBayesian[i] = bayes[i] ; } 
}
//-----------------------------------------------------------------------
void AliFlowEvent::SetNoWgt(Bool_t nowgt) 
{ 
 // Sets no phi weightening, but Wgt = 1*sign(Eta) for odd harmonics .
 
 fNoWgt = nowgt ; // cout << " Wgt = +1 (positive Eta) or -1 (negative Eta) . " << endl ;
}
//-----------------------------------------------------------------------
void AliFlowEvent::SetOnePhiWgt()				  { fOnePhiWgt = kTRUE ; }
void AliFlowEvent::SetFirstLastPhiWgt() 			  { fOnePhiWgt = kFALSE ; }
void AliFlowEvent::SetPtWgt(Bool_t PtWgt)			  { fPtWgt = PtWgt; }
void AliFlowEvent::SetEtaWgt(Bool_t EtaWgt)			  { fEtaWgt = EtaWgt; }
//-----------------------------------------------------------------------
#ifndef __CINT__
void AliFlowEvent::SetPhiWeight(const Flow::PhiWgt_t& pPhiWgt)  	 { memcpy (fPhiWgt, pPhiWgt, sizeof(Flow::PhiWgt_t)); }
void AliFlowEvent::SetPhiWeightPlus(const Flow::PhiWgt_t& pPhiWgtPlus)   { memcpy (fPhiWgtPlus,  pPhiWgtPlus,  sizeof(Flow::PhiWgt_t)); }
void AliFlowEvent::SetPhiWeightMinus(const Flow::PhiWgt_t& pPhiWgtMinus) { memcpy (fPhiWgtMinus, pPhiWgtMinus, sizeof(Flow::PhiWgt_t)); }
void AliFlowEvent::SetPhiWeightCross(const Flow::PhiWgt_t& pPhiWgtCross) { memcpy (fPhiWgtCross, pPhiWgtCross, sizeof(Flow::PhiWgt_t)); }
#endif
//-----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////
