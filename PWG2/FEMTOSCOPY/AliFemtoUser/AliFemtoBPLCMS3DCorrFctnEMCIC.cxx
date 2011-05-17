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
AliFemtoCorrFctn(),
//fEnergyTotalReal(0),  
// fEnergyMultReal(0),        
//fPzMultReal(0),      
//fPtMultReal(0), 
  fNumerator(0),
  fDenominator(0),
  fEnergyTotalMix(0),      
  fEnergyMultMix(0),      
  fPzMultMix(0),            
  fPtMultMix(0),
  fUseRPSelection(0)
{
  
  // set up numerator
  char tTitNum[101] = "Num";
  strncat(tTitNum,title, 100);
  fNumerator = new TH3D(tTitNum,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // set up denominator
  char tTitDen[101] = "Den";
  strncat(tTitDen,title, 100);
  fDenominator = new TH3D(tTitDen,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);

  //Setup EnergyTotalReal
  /*char tTitNum1[100] = "ESumReal";
  strncat(tTitNum1,title, 100);
  fEnergyTotalReal = new TH3D(tTitNum1,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  
  //Setup EnergyMultReal
  char tTitNum2[100] = "EMultReal";
  strncat(tTitNum2,title, 100);
  fEnergyMultReal = new TH3D(tTitNum2,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  
  //Setup Pz MultReal
  char tTitNum3[100] = "PzMultReal";
  strncat(tTitNum3,title, 100);
  fPzMultReal = new TH3D(tTitNum3,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  
  //Setup Pt MultReal
  char tTitNum4[100] = "PtMultReal";
  strncat(tTitNum4,title, 100);
  fPtMultReal = new TH3D(tTitNum4,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);  */
  
  //Setup EnergyTotalMix
  char tTitNum5[101] = "ESumMix";
  strncat(tTitNum5,title, 100);
  fEnergyTotalMix = new TH3D(tTitNum5,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  
  //Setup EnergyMultMix
  char tTitNum6[101] = "EMultMix";
  strncat(tTitNum6,title, 100);
  fEnergyMultMix = new TH3D(tTitNum6,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  
  //Setup Pz MultMix
  char tTitNum7[101] = "PzMultMix";
  strncat(tTitNum7,title, 100);
  fPzMultMix = new TH3D(tTitNum7,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  
  //Setup Pt MultMix
  char tTitNum8[101] = "PtMultMix";
  strncat(tTitNum8,title, 100);
  fPtMultMix = new TH3D(tTitNum8,title,nbins,QLo,QHi,nbins,QLo,QHi,nbins,QLo,QHi);
  // To enable error bar calculation
  
  /*fEnergyTotalReal->Sumw2();
  fEnergyMultReal->Sumw2();
  fPzMultReal->Sumw2();
  fPtMultReal->Sumw2();  */
  fNumerator->Sumw2();
  fDenominator->Sumw2();
  fEnergyTotalMix->Sumw2();
  fEnergyMultMix->Sumw2();
  fPzMultMix->Sumw2();
  fPtMultMix->Sumw2();
  
}


// Variable bin size constructor :
//qBins array of low-edges for each bin. This is an array of size nbins+1
AliFemtoBPLCMS3DCorrFctnEMCIC::AliFemtoBPLCMS3DCorrFctnEMCIC(char* title, const int& nbins, const float* qBins):
AliFemtoCorrFctn(),
  fNumerator(0),
  fDenominator(0),
  fEnergyTotalMix(0),      
  fEnergyMultMix(0),      
  fPzMultMix(0),            
  fPtMultMix(0),
  fUseRPSelection(0)
{
  
  // set up numerator
  char tTitNum[101] = "Num";
  strncat(tTitNum,title, 100);
  fNumerator = new TH3D(tTitNum,title,nbins,qBins,nbins,qBins,nbins,qBins);
  // set up denominator
  char tTitDen[101] = "Den";
  strncat(tTitDen,title, 100);
  fDenominator = new TH3D(tTitDen,title,nbins,qBins,nbins,qBins,nbins,qBins);

  //Setup EnergyTotalReal
  /*char tTitNum1[100] = "ESumReal";
  strncat(tTitNum1,title, 100);
  fEnergyTotalReal = new TH3D(tTitNum1,title,nbins,qBins,nbins,qBins,nbins,qBins);
  
  //Setup EnergyMultReal
  char tTitNum2[100] = "EMultReal";
  strncat(tTitNum2,title, 100);
  fEnergyMultReal = new TH3D(tTitNum2,title,nbins,qBins,nbins,qBins,nbins,qBins);
  
  //Setup Pz MultReal
  char tTitNum3[100] = "PzMultReal";
  strncat(tTitNum3,title, 100);
  fPzMultReal = new TH3D(tTitNum3,title,nbins,qBins,nbins,qBins,nbins,qBins);
  
  //Setup Pt MultReal
  char tTitNum4[100] = "PtMultReal";
  strncat(tTitNum4,title, 100);
  fPtMultReal = new TH3D(tTitNum4,title,nbins,qBins,nbins,qBins,nbins,qBins);  */
  
  //Setup EnergyTotalMix
  char tTitNum5[101] = "ESumMix";
  strncat(tTitNum5,title, 100);
  fEnergyTotalMix = new TH3D(tTitNum5,title,nbins,qBins,nbins,qBins,nbins,qBins);
  
  //Setup EnergyMultMix
  char tTitNum6[101] = "EMultMix";
  strncat(tTitNum6,title, 100);
  fEnergyMultMix = new TH3D(tTitNum6,title,nbins,qBins,nbins,qBins,nbins,qBins);
  
  //Setup Pz MultMix
  char tTitNum7[101] = "PzMultMix";
  strncat(tTitNum7,title, 100);
  fPzMultMix = new TH3D(tTitNum7,title,nbins,qBins,nbins,qBins,nbins,qBins);
  
  //Setup Pt MultMix
  char tTitNum8[101] = "PtMultMix";
  strncat(tTitNum8,title, 100);
  fPtMultMix = new TH3D(tTitNum8,title,nbins,qBins,nbins,qBins,nbins,qBins);
  // To enable error bar calculation
  
  /*fEnergyTotalReal->Sumw2();
  fEnergyMultReal->Sumw2();
  fPzMultReal->Sumw2();
  fPtMultReal->Sumw2();  */
  fNumerator->Sumw2();
  fDenominator->Sumw2();
  fEnergyTotalMix->Sumw2();
  fEnergyMultMix->Sumw2();
  fPzMultMix->Sumw2();
  fPtMultMix->Sumw2();
}
  






AliFemtoBPLCMS3DCorrFctnEMCIC::AliFemtoBPLCMS3DCorrFctnEMCIC(const AliFemtoBPLCMS3DCorrFctnEMCIC& aCorrFctn) :
  AliFemtoCorrFctn(aCorrFctn),
  /*fEnergyTotalReal(0),  
  fEnergyMultReal(0),        
  fPzMultReal(0),      
  fPtMultReal(0),*/     
  fNumerator(0),
  fDenominator(0),
  fEnergyTotalMix (0),      
  fEnergyMultMix (0),  
  fPzMultMix(0),          
  fPtMultMix(0),
  fUseRPSelection(0)
{
  // Copy constructor
  fNumerator = new TH3D(*aCorrFctn.fNumerator);
  fDenominator = new TH3D(*aCorrFctn.fDenominator);
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
  delete fNumerator;
  delete fDenominator;
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
  if (fNumerator) delete fNumerator;
  fNumerator = new TH3D(*aCorrFctn.fNumerator);
  if (fDenominator) delete fDenominator;
  fDenominator = new TH3D(*aCorrFctn.fDenominator);
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
  
  fUseRPSelection = aCorrFctn.fUseRPSelection;

  return *this;
}

//_________________________
void AliFemtoBPLCMS3DCorrFctnEMCIC::WriteOutHistos(){
  
  fNumerator->Write();
  fDenominator->Write();
  //fEnergyTotalReal->Write();
  //fEnergyMultReal->Write();        
  //fPzMultReal->Write();      
  //fPtMultReal->Write();            
  fEnergyTotalMix->Write();      
  fEnergyMultMix->Write();      
  fPzMultMix->Write();            
  fPtMultMix->Write(); 
  //cout << "write histos emcics" << endl;
}
//______________________________
TList* AliFemtoBPLCMS3DCorrFctnEMCIC::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();
  cout << "Getting list from CFemcic" << endl;
  tOutputList->Add(fNumerator);
  tOutputList->Add(fDenominator);
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
	else if (!(ktc->Pass(pair, arp->GetCurrentReactionPlane()))) return;
      }
    }
    else
      if (!(fPairCut->Pass(pair))) return;
  }
  
  double qOut  = fabs(pair->QOutCMS());  
  double qSide = fabs( pair->QSideCMS());
  double qLong = fabs(pair->QLongCMS());

  fNumerator->Fill(qOut,qSide,qLong);
  
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
	else if (!(ktc->Pass(pair, arp->GetCurrentReactionPlane()))) return;
      }
    }
    else
      if (!(fPairCut->Pass(pair))) return;
  }

  
  double qOut = fabs(pair->QOutCMS());
  double qSide = fabs(pair->QSideCMS());
  double qLong = fabs(pair->QLongCMS());

  fDenominator->Fill(qOut,qSide,qLong);

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


void AliFemtoBPLCMS3DCorrFctnEMCIC::SetUseRPSelection(unsigned short aRPSel)
{
  fUseRPSelection = aRPSel;
}
