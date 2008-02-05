////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnNonIdDR - correlation function for non-identical particles //
// uses k* as a function variable. Stores the correlation function separately //
// for positive and negative signs of k* projections into out, side and long  //
// directions, enabling the calculations of double ratios                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnNonIdDR.h"
//#include "AliFemtoHisto.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoCorrFctnNonIdDR)
#endif

//____________________________
AliFemtoCorrFctnNonIdDR::AliFemtoCorrFctnNonIdDR(char* title, const int& nbins, const float& QinvLo, const float& QinvHi):
  fNumOutP(0), 
  fNumOutN(0),  
  fNumSideP(0), 
  fNumSideN(0), 
  fNumLongP(0), 
  fNumLongN(0), 
  fDenOutP(0),  
  fDenOutN(0),  
  fDenSideP(0), 
  fDenSideN(0), 
  fDenLongP(0), 
  fDenLongN(0)
{
  // Default constructor
  // set up numerators
  char bufname[200];
  snprintf(bufname, 200, "NumOutP%s", title);
  fNumOutP = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "NumOutN%s", title);
  fNumOutN = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "NumSideP%s", title);
  fNumSideP = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "NumSideN%s", title);
  fNumSideN = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "NumLongP%s", title);
  fNumLongP = new TH1D(bufname,title,nbins,QinvLo,QinvHi);
  snprintf(bufname, 200, "NumLongN%s", title);
  fNumLongN = new TH1D(bufname,title,nbins,QinvLo,QinvHi);

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
  fNumOutP->Sumw2(); 
  fNumOutN->Sumw2();  
  fNumSideP->Sumw2(); 
  fNumSideN->Sumw2(); 
  fNumLongP->Sumw2(); 
  fNumLongN->Sumw2(); 
  fDenOutP->Sumw2();  
  fDenOutN->Sumw2();  
  fDenSideP->Sumw2(); 
  fDenSideN->Sumw2(); 
  fDenLongP->Sumw2(); 
  fDenLongN->Sumw2();
}

//____________________________
AliFemtoCorrFctnNonIdDR::AliFemtoCorrFctnNonIdDR(const AliFemtoCorrFctnNonIdDR& aCorrFctn) :
  fNumOutP(0), 
  fNumOutN(0),  
  fNumSideP(0), 
  fNumSideN(0), 
  fNumLongP(0), 
  fNumLongN(0), 
  fDenOutP(0),  
  fDenOutN(0),  
  fDenSideP(0), 
  fDenSideN(0), 
  fDenLongP(0), 
  fDenLongN(0)
{
  // copy constructor
  if (aCorrFctn.fNumOutP)
    fNumOutP = new TH1D(*aCorrFctn.fNumOutP);
  if (aCorrFctn.fNumOutN)
    fNumOutN = new TH1D(*aCorrFctn.fNumOutN);
  if (aCorrFctn.fNumSideP)
    fNumSideP = new TH1D(*aCorrFctn.fNumSideP);
  if (aCorrFctn.fNumSideN)
    fNumSideN = new TH1D(*aCorrFctn.fNumSideN);
  if (aCorrFctn.fNumLongP)
    fNumLongP = new TH1D(*aCorrFctn.fNumLongP);
  if (aCorrFctn.fNumLongN)
    fNumLongN = new TH1D(*aCorrFctn.fNumLongN);

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
AliFemtoCorrFctnNonIdDR::~AliFemtoCorrFctnNonIdDR(){
  delete fNumOutP; 
  delete fNumOutN;  
  delete fNumSideP; 
  delete fNumSideN; 
  delete fNumLongP; 
  delete fNumLongN; 
  delete fDenOutP;  
  delete fDenOutN;  
  delete fDenSideP; 
  delete fDenSideN; 
  delete fDenLongP; 
  delete fDenLongN;
}
//_________________________
AliFemtoCorrFctnNonIdDR& AliFemtoCorrFctnNonIdDR::operator=(const AliFemtoCorrFctnNonIdDR& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  if (aCorrFctn.fNumOutP)
    fNumOutP = new TH1D(*aCorrFctn.fNumOutP);
  if (aCorrFctn.fNumOutN)
    fNumOutN = new TH1D(*aCorrFctn.fNumOutN);
  if (aCorrFctn.fNumSideP)
    fNumSideP = new TH1D(*aCorrFctn.fNumSideP);
  if (aCorrFctn.fNumSideN)
    fNumSideN = new TH1D(*aCorrFctn.fNumSideN);
  if (aCorrFctn.fNumLongP)
    fNumLongP = new TH1D(*aCorrFctn.fNumLongP);
  if (aCorrFctn.fNumLongN)
    fNumLongN = new TH1D(*aCorrFctn.fNumLongN);

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
void AliFemtoCorrFctnNonIdDR::Finish(){
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
AliFemtoString AliFemtoCorrFctnNonIdDR::Report(){
  // construct report
  string stemp = "Non-identical particles Correlation Function Report:\n";
  char ctemp[100];
  sprintf(ctemp,"Number of entries in numerators:\t%E\n",fNumOutP->GetEntries()+fNumOutN->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in denominators:\t%E\n",fDenOutP->GetEntries()+fDenOutN->GetEntries());
  stemp += ctemp;
  //  stemp += mCoulombWeight->Report();
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoCorrFctnNonIdDR::AddRealPair(AliFemtoPair* pair){
  // add true pair
  double tKStar = pair->KStar();
  if (pair->KOut()>0.0)
    fNumOutP->Fill(tKStar);
  else
    fNumOutN->Fill(tKStar);

  if (pair->KSide()>0.0)
    fNumSideP->Fill(tKStar);
  else
    fNumSideN->Fill(tKStar);

  if (pair->KLong()>0.0)
    fNumLongP->Fill(tKStar);
  else
    fNumLongN->Fill(tKStar);

}
//____________________________
void AliFemtoCorrFctnNonIdDR::AddMixedPair(AliFemtoPair* pair){
  // add mixed (background) pair
  double tKStar = pair->KStar();
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
void AliFemtoCorrFctnNonIdDR::Write(){
  fNumOutP->Write(); 
  fNumOutN->Write();  
  fNumSideP->Write(); 
  fNumSideN->Write(); 
  fNumLongP->Write(); 
  fNumLongN->Write(); 
  fDenOutP->Write();  
  fDenOutN->Write();  
  fDenSideP->Write(); 
  fDenSideN->Write(); 
  fDenLongP->Write(); 
  fDenLongN->Write();
}

TList* AliFemtoCorrFctnNonIdDR::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumOutP); 
  tOutputList->Add(fNumOutN);  
  tOutputList->Add(fNumSideP); 
  tOutputList->Add(fNumSideN); 
  tOutputList->Add(fNumLongP); 
  tOutputList->Add(fNumLongN); 
  tOutputList->Add(fDenOutP);  
  tOutputList->Add(fDenOutN);  
  tOutputList->Add(fDenSideP); 
  tOutputList->Add(fDenSideN); 
  tOutputList->Add(fDenLongP); 
  tOutputList->Add(fDenLongN);

  return tOutputList;
}
