/***************************************************************************
 *
 * $Id: AliFemtoCorrFctn3DSphercicalEMCIC.h  $
 *
 * Author: Nicolas Bock, Ohio State University, bock@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: Calculates of the 3D Correlation Function in Spherical
 *              coordinates, and also produces histograms to calculate 
 *              Energy Momentum Conservation Induced Correlations  (EMCICs)
 *
 * This Class produces the following histograms as function of Q, theta, phi
 * (for both real and mixed pairs):
 *        1)   E1 + E2
 *        2)   E1 * E2
 *        3)   Pt1*Pt2
 *        4)   Pz1*Pz2
 *  
 * This class is similar to AliFemtoCorrFctn3DSpherical, but it uses Q
 * instead of K to do the binning. 
 * 
 * NOTE: The EMCIC histograms are not averaged in this class, to obtain 
 * the average, the user needs to divide the real pair histograms by 
 * the numerator, and the mixed pair histograms by the denominator.
 *
 ***************************************************************************
 *
 **************************************************************************/

#include "AliFemtoCorrFctn3DSphericalEMCIC.h"
#include <TMath.h>
#include <TVector2.h>
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoCorrFctn3DSphericalEMCIC)
#endif


//____________________________
AliFemtoCorrFctn3DSphericalEMCIC::AliFemtoCorrFctn3DSphericalEMCIC(char* title, const int& nqbins, const float& QLo, const float& QHi, const int& nphibins, const int& ncthetabins):
  AliFemtoCorrFctn(),
  fNumerator(0),
  fDenominator(0),
/*fEnergyTotalReal(0),  
  fEnergyMultReal(0),        
  fPzMultReal(0),      
  fPtMultReal(0),*/            
  fEnergyTotalMix (0),      
  fEnergyMultMix (0),      
  fPzMultMix(0),            
  fPtMultMix(0),
  fPairCut(0x0)
{
  
  // To better sample on phi shift bin low edge by binsize/2 = Pi/numBins.
  Double_t shiftPhi=TMath::Pi()/nphibins;
  // set up numerator
  char tTitNum[100] = "Num";
  strcat(tTitNum,title);
  fNumerator = new TH3D(tTitNum,title,nqbins,QLo,QHi,nphibins,
			-TMath::Pi()-shiftPhi,TMath::Pi()-shiftPhi,ncthetabins,-1.0,1.0);
  // set up denominator
  char tTitDen[100] = "Den";
  strcat(tTitDen,title);
  fDenominator = new TH3D(tTitDen,title,nqbins,QLo,QHi,nphibins,
			  -TMath::Pi()-shiftPhi,TMath::Pi()-shiftPhi,ncthetabins,-1.0,1.0);

  //Added histograms to calculate EMCICs , Nicolas Bock 19.01.2010
  //Setup EnergyTotalReal
  /*char tTitNum1[100] = "ESumReal";
  strcat(tTitNum1,title);
  fEnergyTotalReal = new TH3D(tTitNum1,title,nqbins,QLo,QHi,nphibins,
			      -TMath::Pi()-shiftPhi,TMath::Pi()-shiftPhi,ncthetabins,-1.0,1.0);
 
  //Setup EnergyMultReal
  char tTitNum2[100] = "EMultReal";
  strcat(tTitNum2,title);
  fEnergyMultReal = new TH3D(tTitNum2,title,nqbins,QLo,QHi,nphibins,
			     -TMath::Pi()-shiftPhi,TMath::Pi()-shiftPhi,ncthetabins,-1.0,1.0);
  
  //Setup Pz MultReal
  char tTitNum3[100] = "PzMultReal";
  strcat(tTitNum3,title);
  fPzMultReal = new TH3D(tTitNum3,title,nqbins,QLo,QHi,nphibins,
			 -TMath::Pi()-shiftPhi,TMath::Pi()-shiftPhi,ncthetabins,-1.0,1.0);

  //Setup Pt MultReal
  char tTitNum4[100] = "PtMultReal";
  strcat(tTitNum4,title);
  fPtMultReal = new TH3D(tTitNum4,title,nqbins,QLo,QHi,nphibins,
			 -TMath::Pi()-shiftPhi,TMath::Pi()-shiftPhi,ncthetabins,-1.0,1.0);
  */ 


  //Setup EnergyTotalMix
  char tTitNum5[100] = "ESumMix";
  strcat(tTitNum5,title);
  fEnergyTotalMix = new TH3D(tTitNum5,title,nqbins,QLo,QHi,nphibins,
			     -TMath::Pi()-shiftPhi,TMath::Pi()-shiftPhi,ncthetabins,-1.0,1.0);
  
  //Setup EnergyMultMix
  char tTitNum6[100] = "EMultMix";
  strcat(tTitNum6,title);
  fEnergyMultMix = new TH3D(tTitNum6,title,nqbins,QLo,QHi,nphibins,
			    -TMath::Pi()-shiftPhi,TMath::Pi()-shiftPhi,ncthetabins,-1.0,1.0);
  
  //Setup Pz MultMix
  char tTitNum7[100] = "PzMultMix";
  strcat(tTitNum7,title);
  fPzMultMix = new TH3D(tTitNum7,title,nqbins,QLo,QHi,nphibins,
			-TMath::Pi()-shiftPhi,TMath::Pi()-shiftPhi,ncthetabins,-1.0,1.0);

  //Setup Pt MultMix
  char tTitNum8[100] = "PtMultMix";
  strcat(tTitNum8,title);
  fPtMultMix = new TH3D(tTitNum8,title,nqbins,QLo,QHi,nphibins,
			-TMath::Pi()-shiftPhi,TMath::Pi()-shiftPhi,
			ncthetabins,-1.0,1.0);

  // to enable error bar calculation...
  fNumerator->Sumw2();
  fDenominator->Sumw2();
  /*fEnergyTotalReal->Sumw2();
  fEnergyMultReal->Sumw2();
  fPzMultReal->Sumw2();
  fPtMultReal->Sumw2();  */
  fEnergyTotalMix->Sumw2();
  fEnergyMultMix->Sumw2();
  fPzMultMix->Sumw2();
  fPtMultMix->Sumw2();


}

AliFemtoCorrFctn3DSphericalEMCIC::AliFemtoCorrFctn3DSphericalEMCIC(const AliFemtoCorrFctn3DSphericalEMCIC& aCorrFctn) :
  AliFemtoCorrFctn(),
  fNumerator(0),
  fDenominator(0),
  /*fEnergyTotalReal(0),
  fEnergyMultReal(0),        
  fPzMultReal(0),      
  fPtMultReal(0),            */
  fEnergyTotalMix (0),      
  fEnergyMultMix (0),      
  fPzMultMix(0),            
  fPtMultMix(0),
  fPairCut(0x0)
{
  // Copy constructor
  fNumerator = new TH3D(*aCorrFctn.fNumerator);
  fDenominator = new TH3D(*aCorrFctn.fDenominator);
  /*fEnergyTotalReal = new TH3D(*aCorrFctn.fEnergyTotalReal);
  fEnergyMultReal = new TH3D(*aCorrFctn.fEnergyMultReal);
  fPzMultReal = new TH3D(*aCorrFctn.fPzMultReal);
  fPtMultReal = new TH3D(*aCorrFctn.fPtMultReal);*/
  fEnergyTotalMix = new TH3D(*aCorrFctn.fEnergyTotalMix);
  fEnergyMultMix = new TH3D(*aCorrFctn.fEnergyMultMix);
  fPzMultMix = new TH3D(*aCorrFctn.fPzMultMix);
  fPtMultMix = new TH3D(*aCorrFctn.fPtMultMix);
  fPairCut = aCorrFctn.fPairCut;
}
//____________________________
AliFemtoCorrFctn3DSphericalEMCIC::~AliFemtoCorrFctn3DSphericalEMCIC(){
  // Destructor
  delete fNumerator;
  delete fDenominator;
  /*delete fEnergyTotalReal;
  delete fEnergyMultReal;        
  delete fPzMultReal;     
  delete fPtMultReal;            */
  delete fEnergyTotalMix;      
  delete fEnergyMultMix; 
  delete fPzMultMix;   
  delete fPtMultMix;
}
//_________________________
AliFemtoCorrFctn3DSphericalEMCIC& AliFemtoCorrFctn3DSphericalEMCIC::operator=(const AliFemtoCorrFctn3DSphericalEMCIC& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (fNumerator) delete fNumerator;
  fNumerator = new TH3D(*aCorrFctn.fNumerator);
  if (fDenominator) delete fDenominator;
  fDenominator = new TH3D(*aCorrFctn.fDenominator);
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
  fPairCut = aCorrFctn.fPairCut;
  
  return *this;
}

//_________________________
void AliFemtoCorrFctn3DSphericalEMCIC::WriteOutHistos(){
  // Write out all histograms to file
  fNumerator->Write();
  fDenominator->Write();
  /*fEnergyTotalReal->Write();
  fEnergyMultReal->Write();        
  fPzMultReal->Write();      
  fPtMultReal->Write();            */
  fEnergyTotalMix->Write();      
  fEnergyMultMix->Write();      
  fPzMultMix->Write();            
  fPtMultMix->Write();
}
//______________________________
TList* AliFemtoCorrFctn3DSphericalEMCIC::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumerator); 
  tOutputList->Add(fDenominator);  
  /*  tOutputList->Add(fEnergyTotalReal);
  tOutputList->Add(fEnergyMultReal);        
  tOutputList->Add(fPzMultReal);      
  tOutputList->Add(fPtMultReal);            */
  tOutputList->Add(fEnergyTotalMix );      
  tOutputList->Add(fEnergyMultMix );      
  tOutputList->Add(fPzMultMix);            
  tOutputList->Add(fPtMultMix);
  return tOutputList;
}

//_________________________
void AliFemtoCorrFctn3DSphericalEMCIC::Finish(){
  // here is where we should normalize, fit, etc...
}

//____________________________
AliFemtoString AliFemtoCorrFctn3DSphericalEMCIC::Report(){
  // Construct the report
  string stemp = "PRF Frame SphericalEMCIC 3D Correlation Function Report:\n";
  char ctemp[100];
  sprintf(ctemp,"Number of entries in numerator:\t%E\n",fNumerator->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in denominator:\t%E\n",fDenominator->GetEntries());
  stemp += ctemp;

  if (fPairCut){
    sprintf(ctemp,"Here is the PairCut specific to this CorrFctn\n");
    stemp += ctemp;
    stemp += fPairCut->Report();
  }
  else{
    sprintf(ctemp,"No PairCut specific to this CorrFctn\n");
    stemp += ctemp;
  }

  //  
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoCorrFctn3DSphericalEMCIC::AddRealPair( AliFemtoPair* pair){
  // perform operations on real pairs
  if (fPairCut){
    if (!(fPairCut->Pass(pair))) return;
  }

  //                          
  double tQO = pair->QOutCMS();  
  double tQS = pair->QSideCMS();  
  double tQL = pair->QLongCMS();  

  double tQR = sqrt(tQO*tQO + tQS*tQS + tQL*tQL);
  double tQC = 0;
  if ( fabs(tQR) < 1e-10 ) tQC = 0.0;
  else tQC = tQL/tQR;
  double tQP = atan2(tQS,tQO);

  fNumerator->Fill(tQR,tQP,tQC);
  
  // EMCICs  
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
  
  fEnergyTotalReal->Fill(tQR,tQP,tQC,tE1+tE2);
  fEnergyMultReal->Fill(tQR,tQP,tQC,tE1*tE2);
  fPzMultReal->Fill(tQR,tQP,tQC,tPz1*tPz2);
  fPtMultReal->Fill(tQR,tQP,tQC,tPt1DotPt2);*/
   

}
//____________________________
void AliFemtoCorrFctn3DSphericalEMCIC::AddMixedPair( AliFemtoPair* pair){
  // perform operations on mixed pairs
  if (fPairCut){
    if (!(fPairCut->Pass(pair))) return;
  }
 //                          //Changed K to Q to be in LCMS, N. Bock
  double tQO = pair->QOutCMS();  
  double tQS = pair->QSideCMS();   
  double tQL = pair->QLongCMS();  

  double tQR = sqrt(tQO*tQO + tQS*tQS + tQL*tQL);
  double tQC;
  if ( fabs(tQR) < 1e-10 ) tQC = 0.0;
  else tQC=tQL/tQR;
  double tQP=atan2(tQS,tQO);

  fDenominator->Fill(tQR,tQP,tQC);

  // EMCICs   
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
  
  fEnergyTotalMix->Fill(tQR,tQP,tQC,tE1+tE2);
  fEnergyMultMix->Fill(tQR,tQP,tQC,tE1*tE2);
  fPzMultMix->Fill(tQR,tQP,tQC,tPz1*tPz2);
  fPtMultMix->Fill(tQR,tQP,tQC,tPt1DotPt2);
  
}
