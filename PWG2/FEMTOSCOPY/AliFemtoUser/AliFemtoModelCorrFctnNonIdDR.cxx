////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoModelCorrFctnNonIdDR - correlation function for non-identical      //
// particles which uses k* as a function variable. Stores the correlation     //
// function separately for positive and negative signs of k* projections into //
// out, side and long directions, enabling the calculations of double ratios  //
// Uses pair weight to simulate the model correlation function.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoModelCorrFctnNonIdDR.h"
#include "AliFemtoModelManager.h"
//#include "AliFemtoHisto.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoModelCorrFctnNonIdDR)
#endif

//____________________________
AliFemtoModelCorrFctnNonIdDR::AliFemtoModelCorrFctnNonIdDR(char* title, const int& nbins, const float& QinvLo, const float& QinvHi):
  AliFemtoModelCorrFctn(title, nbins, QinvLo, QinvHi),
  fNumTOutP(0), 
  fNumTOutN(0),  
  fNumTSideP(0), 
  fNumTSideN(0), 
  fNumTLongP(0), 
  fNumTLongN(0), 
  fNumFOutP(0), 
  fNumFOutN(0),  
  fNumFSideP(0), 
  fNumFSideN(0), 
  fNumFLongP(0), 
  fNumFLongN(0), 
  fDenOutP(0),  
  fDenOutN(0),  
  fDenSideP(0), 
  fDenSideN(0), 
  fDenLongP(0), 
  fDenLongN(0)
{
  // Default constructor
  char bufname[200];

  // set up true numerators
  snprintf(bufname, 200, "NumTOutP%s", title);
  fNumTOutP = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "NumTOutN%s", title);
  fNumTOutN = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "NumTSideP%s", title);
  fNumTSideP = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "NumTSideN%s", title);
  fNumTSideN = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "NumTLongP%s", title);
  fNumTLongP = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "NumTLongN%s", title);
  fNumTLongN = new TH1D(bufname,title,nbins,QinvLo,QinvHi);

  // set up fake numerators
  snprintf(bufname, 200, "NumFOutP%s", title);
  fNumFOutP = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "NumFOutN%s", title);
  fNumFOutN = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "NumFSideP%s", title);
  fNumFSideP = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "NumFSideN%s", title);
  fNumFSideN = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "NumFLongP%s", title);
  fNumFLongP = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "NumFLongN%s", title);
  fNumFLongN = new TH1D(bufname,title,nbins,QinvLo,QinvHi);

  // set up denominators
  snprintf(bufname, 200, "DenOutP%s", title);
  fDenOutP = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "DenOutN%s", title);
  fDenOutN = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "DenSideP%s", title);
  fDenSideP = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "DenSideN%s", title);
  fDenSideN = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "DenLongP%s", title);
  fDenLongP = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "DenLongN%s", title);
  fDenLongN = new TH1D(bufname,title,nbins,QinvLo,QinvHi);

  // to enable error bar calculation...
  fNumTOutP->Sumw2(); 
  fNumTOutN->Sumw2();  
  fNumTSideP->Sumw2(); 
  fNumTSideN->Sumw2(); 
  fNumTLongP->Sumw2(); 
  fNumTLongN->Sumw2(); 
  fNumFOutP->Sumw2(); 
  fNumFOutN->Sumw2();  
  fNumFSideP->Sumw2(); 
  fNumFSideN->Sumw2(); 
  fNumFLongP->Sumw2(); 
  fNumFLongN->Sumw2(); 
  fDenOutP->Sumw2();  
  fDenOutN->Sumw2();  
  fDenSideP->Sumw2(); 
  fDenSideN->Sumw2(); 
  fDenLongP->Sumw2(); 
  fDenLongN->Sumw2();
}

//____________________________
AliFemtoModelCorrFctnNonIdDR::AliFemtoModelCorrFctnNonIdDR(const AliFemtoModelCorrFctnNonIdDR& aCorrFctn) :
  AliFemtoModelCorrFctn(),
  fNumTOutP(0), 
  fNumTOutN(0),  
  fNumTSideP(0), 
  fNumTSideN(0), 
  fNumTLongP(0), 
  fNumTLongN(0), 
  fNumFOutP(0), 
  fNumFOutN(0),  
  fNumFSideP(0), 
  fNumFSideN(0), 
  fNumFLongP(0), 
  fNumFLongN(0), 
  fDenOutP(0),  
  fDenOutN(0),  
  fDenSideP(0), 
  fDenSideN(0), 
  fDenLongP(0), 
  fDenLongN(0)
{
  // copy constructor
  if (aCorrFctn.fNumTOutP)
    fNumTOutP = new TH1D(*aCorrFctn.fNumTOutP);
  if (aCorrFctn.fNumTOutN)
    fNumTOutN = new TH1D(*aCorrFctn.fNumTOutN);
  if (aCorrFctn.fNumTSideP)
    fNumTSideP = new TH1D(*aCorrFctn.fNumTSideP);
  if (aCorrFctn.fNumTSideN)
    fNumTSideN = new TH1D(*aCorrFctn.fNumTSideN);
  if (aCorrFctn.fNumTLongP)
    fNumTLongP = new TH1D(*aCorrFctn.fNumTLongP);
  if (aCorrFctn.fNumTLongN)
    fNumTLongN = new TH1D(*aCorrFctn.fNumTLongN);

  if (aCorrFctn.fNumFOutP)
    fNumFOutP = new TH1D(*aCorrFctn.fNumFOutP);
  if (aCorrFctn.fNumFOutN)
    fNumFOutN = new TH1D(*aCorrFctn.fNumFOutN);
  if (aCorrFctn.fNumFSideP)
    fNumFSideP = new TH1D(*aCorrFctn.fNumFSideP);
  if (aCorrFctn.fNumFSideN)
    fNumFSideN = new TH1D(*aCorrFctn.fNumFSideN);
  if (aCorrFctn.fNumFLongP)
    fNumFLongP = new TH1D(*aCorrFctn.fNumFLongP);
  if (aCorrFctn.fNumFLongN)
    fNumFLongN = new TH1D(*aCorrFctn.fNumFLongN);

  if (aCorrFctn.fDenOutP)
    fDenOutP = new TH1D(*aCorrFctn.fDenOutP);
  if (aCorrFctn.fDenOutN)
    fDenOutN = new TH1D(*aCorrFctn.fDenOutN);
  if (aCorrFctn.fDenSideP)
    fDenSideP = new TH1D(*aCorrFctn.fDenSideP);
  if (aCorrFctn.fDenSideN)
    fDenSideN = new TH1D(*aCorrFctn.fDenSideN);
  if (aCorrFctn.fDenLongP)
    fDenLongP = new TH1D(*aCorrFctn.fDenLongP);
  if (aCorrFctn.fDenLongN)
    fDenLongN = new TH1D(*aCorrFctn.fDenLongN);
}
//____________________________
AliFemtoModelCorrFctnNonIdDR::~AliFemtoModelCorrFctnNonIdDR(){
  // Destructor
  delete fNumTOutP; 
  delete fNumTOutN;  
  delete fNumTSideP; 
  delete fNumTSideN; 
  delete fNumTLongP; 
  delete fNumTLongN; 
  delete fNumFOutP; 
  delete fNumFOutN;  
  delete fNumFSideP; 
  delete fNumFSideN; 
  delete fNumFLongP; 
  delete fNumFLongN; 
  delete fDenOutP;  
  delete fDenOutN;  
  delete fDenSideP; 
  delete fDenSideN; 
  delete fDenLongP; 
  delete fDenLongN;
}
//_________________________
AliFemtoModelCorrFctnNonIdDR& AliFemtoModelCorrFctnNonIdDR::operator=(const AliFemtoModelCorrFctnNonIdDR& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (aCorrFctn.fNumTOutP)
    fNumTOutP = new TH1D(*aCorrFctn.fNumTOutP);
  if (aCorrFctn.fNumTOutN)
    fNumTOutN = new TH1D(*aCorrFctn.fNumTOutN);
  if (aCorrFctn.fNumTSideP)
    fNumTSideP = new TH1D(*aCorrFctn.fNumTSideP);
  if (aCorrFctn.fNumTSideN)
    fNumTSideN = new TH1D(*aCorrFctn.fNumTSideN);
  if (aCorrFctn.fNumTLongP)
    fNumTLongP = new TH1D(*aCorrFctn.fNumTLongP);
  if (aCorrFctn.fNumTLongN)
    fNumTLongN = new TH1D(*aCorrFctn.fNumTLongN);

  if (aCorrFctn.fNumFOutP)
    fNumFOutP = new TH1D(*aCorrFctn.fNumFOutP);
  if (aCorrFctn.fNumFOutN)
    fNumFOutN = new TH1D(*aCorrFctn.fNumFOutN);
  if (aCorrFctn.fNumFSideP)
    fNumFSideP = new TH1D(*aCorrFctn.fNumFSideP);
  if (aCorrFctn.fNumFSideN)
    fNumFSideN = new TH1D(*aCorrFctn.fNumFSideN);
  if (aCorrFctn.fNumFLongP)
    fNumFLongP = new TH1D(*aCorrFctn.fNumFLongP);
  if (aCorrFctn.fNumFLongN)
    fNumFLongN = new TH1D(*aCorrFctn.fNumFLongN);

  if (aCorrFctn.fDenOutP)
    fDenOutP = new TH1D(*aCorrFctn.fDenOutP);
  if (aCorrFctn.fDenOutN)
    fDenOutN = new TH1D(*aCorrFctn.fDenOutN);
  if (aCorrFctn.fDenSideP)
    fDenSideP = new TH1D(*aCorrFctn.fDenSideP);
  if (aCorrFctn.fDenSideN)
    fDenSideN = new TH1D(*aCorrFctn.fDenSideN);
  if (aCorrFctn.fDenLongP)
    fDenLongP = new TH1D(*aCorrFctn.fDenLongP);
  if (aCorrFctn.fDenLongN)
    fDenLongN = new TH1D(*aCorrFctn.fDenLongN);

  return *this;
}

//_________________________
void AliFemtoModelCorrFctnNonIdDR::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  fNumerator->Draw();
  //fDenominator->Draw();
  //fRatio->Draw();
  //  fRatio->Divide(fNumerator,fDenominator,1.0,1.0);

}

//____________________________
AliFemtoString AliFemtoModelCorrFctnNonIdDR::Report(){
  // construct report
  string stemp = "Non-identical particles Model Correlation Function Report:\n";
  char ctemp[100];
  sprintf(ctemp,"Number of entries in numerators:\t%E\n",fNumTOutP->GetEntries()+fNumTOutN->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in denominators:\t%E\n",fDenOutP->GetEntries()+fDenOutN->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoModelCorrFctnNonIdDR::AddRealPair(AliFemtoPair* pair){
  // add true pair
  double tKStar = pair->KStar();
  Double_t weight = fManager->GetWeight(pair);

  if (pair->KOut()>0.0)
    fNumTOutP->Fill(tKStar, weight);
  else
    fNumTOutN->Fill(tKStar, weight);

  if (pair->KSide()>0.0)
    fNumTSideP->Fill(tKStar, weight);
  else
    fNumTSideN->Fill(tKStar, weight);

  if (pair->KLong()>0.0)
    fNumTLongP->Fill(tKStar, weight);
  else
    fNumTLongN->Fill(tKStar, weight);

}
//____________________________
void AliFemtoModelCorrFctnNonIdDR::AddMixedPair(AliFemtoPair* pair){
  // add mixed (background) pair
  double tKStar = pair->KStar();
  Double_t weight = fManager->GetWeight(pair);

  if (pair->KOut()>0.0)
    fNumFOutP->Fill(tKStar, weight);
  else
    fNumFOutN->Fill(tKStar, weight);

  if (pair->KSide()>0.0)
    fNumFSideP->Fill(tKStar, weight);
  else
    fNumFSideN->Fill(tKStar, weight);

  if (pair->KLong()>0.0)
    fNumFLongP->Fill(tKStar, weight);
  else
    fNumFLongN->Fill(tKStar, weight);

  if (pair->KOut()>0.0)
    fDenOutP->Fill(tKStar);
  else
    fDenOutN->Fill(tKStar);

  if (pair->KSide()>0.0)
    fDenSideP->Fill(tKStar);
  else
    fDenSideN->Fill(tKStar);

  if (pair->KLong()>0.0)
    fDenLongP->Fill(tKStar);
  else
    fDenLongN->Fill(tKStar);
}
//____________________________
void AliFemtoModelCorrFctnNonIdDR::Write(){
  // Write out histos
  fNumTOutP->Write(); 
  fNumTOutN->Write();  
  fNumTSideP->Write(); 
  fNumTSideN->Write(); 
  fNumTLongP->Write(); 
  fNumTLongN->Write(); 
  fNumFOutP->Write(); 
  fNumFOutN->Write();  
  fNumFSideP->Write(); 
  fNumFSideN->Write(); 
  fNumFLongP->Write(); 
  fNumFLongN->Write(); 
  fDenOutP->Write();  
  fDenOutN->Write();  
  fDenSideP->Write(); 
  fDenSideN->Write(); 
  fDenLongP->Write(); 
  fDenLongN->Write();
}

TList* AliFemtoModelCorrFctnNonIdDR::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumTOutP); 
  tOutputList->Add(fNumTOutN);  
  tOutputList->Add(fNumTSideP); 
  tOutputList->Add(fNumTSideN); 
  tOutputList->Add(fNumTLongP); 
  tOutputList->Add(fNumTLongN); 
  tOutputList->Add(fNumFOutP); 
  tOutputList->Add(fNumFOutN);  
  tOutputList->Add(fNumFSideP); 
  tOutputList->Add(fNumFSideN); 
  tOutputList->Add(fNumFLongP); 
  tOutputList->Add(fNumFLongN); 
  tOutputList->Add(fDenOutP);  
  tOutputList->Add(fDenOutN);  
  tOutputList->Add(fDenSideP); 
  tOutputList->Add(fDenSideN); 
  tOutputList->Add(fDenLongP); 
  tOutputList->Add(fDenLongN);

  return tOutputList;
}

//_______________________
AliFemtoModelCorrFctn* AliFemtoModelCorrFctnNonIdDR::Clone()
{
  // Create clone
  AliFemtoModelCorrFctnNonIdDR *tCopy = new AliFemtoModelCorrFctnNonIdDR(*this);
  
  return tCopy;
}
