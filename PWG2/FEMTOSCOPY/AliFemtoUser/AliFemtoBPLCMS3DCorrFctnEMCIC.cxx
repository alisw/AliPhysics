/***************************************************************************
 *
 * $Id: AliFemtoBPLCMS3DCorrFctnEMCIC.cxx  $
 *
 * Author: Nicolas Bock, Ohio State University, bock@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: Calculates of the 3D Correlation Function, and also
 *              produces histograms to calculate Energy Momentum Conservation
 *              Induced Correlations  (EMCICs)
 *
 * This Class produces the following histograms as function of Qinv
 * (for both real and mixed pairs):
 *        1)   E1 + E2
 *        2)   E1 * E2
 *        3)   Pt1*Pt2
 *        4)   Pz1*Pz2
 *  
 * The class is derived from AliFemtoBPLCMS3DCorrFctn, therefore it produces
 * also the histograms in that class. 
 * 
 * NOTE: The EMCIC histograms are not averaged in this class, to obtain 
 * the average, the user needs to divide the real pair histograms by 
 * the numerator, and the mixed pair histograms by the denominator
 *
 ***************************************************************************
 *
 **************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoBPLCMS3DCorrFctn: a class to calculate 3D correlation         //
// for pairs of identical particles.                                     //
// It also stored the weighted qinv per bin histogram for the coulomb    //
// correction.                                                           //
// In analysis the function should be first created in a macro, then     //
// added to the analysis, and at the end of the macro the procedure to   //
// write out histograms should be called.                                //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliFemtoBPLCMS3DCorrFctnEMCIC.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoAnalysisReactionPlane.h"
//#include "AliFemtoHisto.h"
#include <cstdio>
#include <TVector2.h>

#ifdef __ROOT__ 
ClassImp(AliFemtoBPLCMS3DCorrFctnEMCIC)
#endif

//____________________________
AliFemtoBPLCMS3DCorrFctnEMCIC::AliFemtoBPLCMS3DCorrFctnEMCIC(char* title, const int& nbins, const float& QLo, const float& QHi)
:
AliFemtoBPLCMS3DCorrFctn(title,nbins, QLo, QHi),
//fEnergyTotalReal(0),  
// fEnergyMultReal(0),        
//fPzMultReal(0),      
//fPtMultReal(0),            
  fEnergyTotalMix(0),      
  fEnergyMultMix(0),      
  fPzMultMix(0),            
  fPtMultMix(0) 
{
  
  //Setup EnergyTotalReal
  /*char tTitNum1[100] = "ESumReal";
  strcat(tTitNum1,title);
  fEnergyTotalReal = new TH3D(tTitNum1,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
  
  //Setup EnergyMultReal
  char tTitNum2[100] = "EMultReal";
  strcat(tTitNum2,title);
  fEnergyMultReal = new TH3D(tTitNum2,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
  
  //Setup Pz MultReal
  char tTitNum3[100] = "PzMultReal";
  strcat(tTitNum3,title);
  fPzMultReal = new TH3D(tTitNum3,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
  
  //Setup Pt MultReal
  char tTitNum4[100] = "PtMultReal";
  strcat(tTitNum4,title);
  fPtMultReal = new TH3D(tTitNum4,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);  */
  
  //Setup EnergyTotalMix
  char tTitNum5[100] = "ESumMix";
  strcat(tTitNum5,title);
  fEnergyTotalMix = new TH3D(tTitNum5,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
  
  //Setup EnergyMultMix
  char tTitNum6[100] = "EMultMix";
  strcat(tTitNum6,title);
  fEnergyMultMix = new TH3D(tTitNum6,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
  
  //Setup Pz MultMix
  char tTitNum7[100] = "PzMultMix";
  strcat(tTitNum7,title);
  fPzMultMix = new TH3D(tTitNum7,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
  
  //Setup Pt MultMix
  char tTitNum8[100] = "PtMultMix";
  strcat(tTitNum8,title);
  fPtMultMix = new TH3D(tTitNum8,title,nbins,-QHi,QHi,nbins,-QHi,QHi,nbins,-QHi,QHi);
  // To enable error bar calculation
  
  /*fEnergyTotalReal->Sumw2();
  fEnergyMultReal->Sumw2();
  fPzMultReal->Sumw2();
  fPtMultReal->Sumw2();  */
  fEnergyTotalMix->Sumw2();
  fEnergyMultMix->Sumw2();
  fPzMultMix->Sumw2();
  fPtMultMix->Sumw2();
  
}

AliFemtoBPLCMS3DCorrFctnEMCIC::AliFemtoBPLCMS3DCorrFctnEMCIC(const AliFemtoBPLCMS3DCorrFctnEMCIC& aCorrFctn) :
  AliFemtoBPLCMS3DCorrFctn(aCorrFctn),
  /*fEnergyTotalReal(0),  
  fEnergyMultReal(0),        
  fPzMultReal(0),      
  fPtMultReal(0),*/            
  fEnergyTotalMix (0),      
  fEnergyMultMix (0),  
  fPzMultMix(0),          
  fPtMultMix(0) 
{
  // Copy constructor
  
  /*  fEnergyTotalReal = new TH3D(*aCorrFctn.fEnergyTotalReal);
  fEnergyMultReal = new TH3D(*aCorrFctn.fEnergyMultReal);
  fPzMultReal = new TH3D(*aCorrFctn.fPzMultReal);
  fPtMultReal = new TH3D(*aCorrFctn.fPtMultReal);  */
  fEnergyTotalMix = new TH3D(*aCorrFctn.fEnergyTotalMix);
  fEnergyMultMix = new TH3D(*aCorrFctn.fEnergyMultMix);
  fPzMultMix = new TH3D(*aCorrFctn.fPzMultMix);
  fPtMultMix = new TH3D(*aCorrFctn.fPtMultMix);
}
//____________________________
AliFemtoBPLCMS3DCorrFctnEMCIC::~AliFemtoBPLCMS3DCorrFctnEMCIC(){
  // Destructor
  /*  delete fEnergyTotalReal;
  delete fEnergyMultReal;        
  delete fPzMultReal;     
  delete fPtMultReal;  */
  delete fEnergyTotalMix;      
  delete fEnergyMultMix; 
  delete fPzMultMix;   
  delete fPtMultMix;
}
//_________________________
AliFemtoBPLCMS3DCorrFctnEMCIC& AliFemtoBPLCMS3DCorrFctnEMCIC::operator=(const AliFemtoBPLCMS3DCorrFctnEMCIC& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;
  
  //Emcics
  /*  if (fEnergyTotalReal) delete fEnergyTotalReal;
  fEnergyTotalReal = new TH3D(*aCorrFctn.fEnergyTotalReal);
  if (fEnergyMultReal) delete fEnergyMultReal;
  fEnergyMultReal = new TH3D(*aCorrFctn.fEnergyMultReal);
  if (fPzMultReal) delete fPzMultReal;
  fPzMultReal = new TH3D(*aCorrFctn.fPzMultReal);
  if (fPtMultReal) delete fPtMultReal;
  fPtMultReal = new TH3D(*aCorrFctn.fPtMultReal);  */
  if (fEnergyTotalMix) delete fEnergyTotalMix;
  fEnergyTotalMix = new TH3D(*aCorrFctn.fEnergyTotalMix);
  if (fEnergyMultMix) delete fEnergyMultMix;
  fEnergyMultMix = new TH3D(*aCorrFctn.fEnergyMultMix);
  if (fPzMultMix) delete fPzMultMix;
  fPzMultMix = new TH3D(*aCorrFctn.fPzMultMix);
  if (fPtMultMix) delete fPtMultMix;
  fPtMultMix = new TH3D(*aCorrFctn.fPtMultMix);
  
  return *this;
}

//_________________________
void AliFemtoBPLCMS3DCorrFctnEMCIC::WriteOutHistos(){
  TH3D* Num3D = Numerator();
  TH3D* Den3D = Denominator();
  Num3D->Write();
  Den3D->Write();
  //fEnergyTotalReal->Write();
  //fEnergyMultReal->Write();        
  //fPzMultReal->Write();      
  //fPtMultReal->Write();            
  fEnergyTotalMix->Write();      
  fEnergyMultMix->Write();      
  fPzMultMix->Write();            
  fPtMultMix->Write(); 
  cout << "write histos emcics" << endl;
}
//______________________________
TList* AliFemtoBPLCMS3DCorrFctnEMCIC::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();
  cout << "Getting list from CFemcic" << endl;
  tOutputList->Add(Numerator());
  tOutputList->Add(Denominator());
  /*tOutputList->Add(fEnergyTotalReal);
  tOutputList->Add(fEnergyMultReal);        
  tOutputList->Add(fPzMultReal);      
  tOutputList->Add(fPtMultReal);     */       
  tOutputList->Add(fEnergyTotalMix );      
  tOutputList->Add(fEnergyMultMix );      
  tOutputList->Add(fPzMultMix);            
  tOutputList->Add(fPtMultMix);
  return tOutputList;
}



//____________________________
void AliFemtoBPLCMS3DCorrFctnEMCIC::AddRealPair( AliFemtoPair* pair){
  // perform operations on real pairs
  
  //WILL WE NEED THIS PART HERE>?????????? I think so. 
  //  cout << "ADDing real Pair to 3D emcics!" << endl; 

 if (fPairCut){
    if (fUseRPSelection) {
      AliFemtoKTPairCut *ktc = dynamic_cast<AliFemtoKTPairCut *>(fPairCut);
      if (!ktc) { 
	cout << "RP aware cut requested, but not connected to the CF" << endl;
	if (!(fPairCut->Pass(pair))) return;
      }
      else {
	AliFemtoAnalysisReactionPlane *arp = dynamic_cast<AliFemtoAnalysisReactionPlane *> (HbtAnalysis());
	if (!arp) {
	  cout << "RP aware cut requested, but not connected to the CF" << endl;
	  if (!(fPairCut->Pass(pair))) return;
	}
	if (!(ktc->Pass(pair, arp->GetCurrentReactionPlane()))) return;
      }
    }
    else
      if (!(fPairCut->Pass(pair))) return;
      }

  
 double qOut = pair->QOutCMS();  //Removed fabs() to bin in negative Qosl
  double qSide = pair->QSideCMS();
  double qLong = pair->QLongCMS();

  AliFemtoBPLCMS3DCorrFctn::Numerator()->Fill(qOut,qSide,qLong);

  /*AliFemtoLorentzVector tMom1 = pair->Track1()->FourMomentum();
  AliFemtoLorentzVector tMom2 = pair->Track2()->FourMomentum();
  double tE1 = tMom1.e();
  double tE2 = tMom2.e();
  double tPz1 = tMom1.pz();
  double tPz2 = tMom2.pz();
  TVector2 tPt1;  
  TVector2 tPt2; 
  tPt1.Set(tMom1.px(),tMom1.py());
  tPt2.Set(tMom2.px(),tMom2.py());
  
  
  double tPt1DotPt2 = tPt1*tPt2;
  
  fEnergyTotalReal->Fill(qOut,qSide,qLong,tE1+tE2);
  fEnergyMultReal->Fill(qOut,qSide,qLong,tE1*tE2);
  fPzMultReal->Fill(qOut,qSide,qLong,tPz1*tPz2);
  fPtMultReal->Fill(qOut,qSide,qLong,tPt1DotPt2); */
}
//____________________________
void AliFemtoBPLCMS3DCorrFctnEMCIC::AddMixedPair( AliFemtoPair* pair){
  // perform operations on mixed pairs
//   if (fPairCut){
//     if (!(fPairCut->Pass(pair))) return;
//   }
  if (fPairCut){
    if (fUseRPSelection) {
      AliFemtoKTPairCut *ktc = dynamic_cast<AliFemtoKTPairCut *>(fPairCut);
      if (!ktc) { 
	cout << "RP aware cut requested, but not connected to the CF" << endl;
	if (!(fPairCut->Pass(pair))) return;
      }
      else {
	AliFemtoAnalysisReactionPlane *arp = dynamic_cast<AliFemtoAnalysisReactionPlane *> (HbtAnalysis());
	if (!arp) {
	  cout << "RP aware cut requested, but not connected to the CF" << endl;
	  if (!(fPairCut->Pass(pair))) return;
	}
	if (!(ktc->Pass(pair, arp->GetCurrentReactionPlane()))) return;
      }
    }
    else
      if (!(fPairCut->Pass(pair))) return;
  }

  
  double qOut = pair->QOutCMS();
  double qSide = pair->QSideCMS();
  double qLong = pair->QLongCMS();

  AliFemtoBPLCMS3DCorrFctn::Denominator()->Fill(qOut,qSide,qLong);

  AliFemtoLorentzVector tMom1 = pair->Track1()->FourMomentum();
  AliFemtoLorentzVector tMom2 = pair->Track2()->FourMomentum();
  double tE1 = tMom1.e();
  double tE2 = tMom2.e();
  double tPz1 = tMom1.pz();
  double tPz2 = tMom2.pz();
  TVector2 tPt1;  
  TVector2 tPt2; 
  tPt1.Set(tMom1.px(),tMom1.py());
  tPt2.Set(tMom2.px(),tMom2.py());
  double tPt1DotPt2 = tPt1*tPt2;
  
  fEnergyTotalMix->Fill(qOut,qSide,qLong,tE1+tE2);
  fEnergyMultMix->Fill(qOut,qSide,qLong,tE1*tE2);
  fPzMultMix->Fill(qOut,qSide,qLong,tPz1*tPz2);
  fPtMultMix->Fill(qOut,qSide,qLong,tPt1DotPt2);
  
}

