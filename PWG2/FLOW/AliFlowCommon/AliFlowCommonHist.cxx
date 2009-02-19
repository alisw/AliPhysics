/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "Riostream.h"              //needed as include
#include "AliFlowCommonConstants.h" //needed as include
#include "AliFlowCommonHist.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"

#include "TString.h" 
#include "TProfile.h"
#include "TMath.h"   //needed as include
#include "TList.h"
#include "AliFlowVector.h"
#include "TBrowser.h"

class TH1F;
class TH1D;

// AliFlowCommonHist:
//
// Description: Class to organise common histograms for Flow Analysis

// authors: N. van der Kolk (kolk@nikhef.nl), A. Bilandzic (anteb@nikhef.nl), RS


ClassImp(AliFlowCommonHist)

//-----------------------------------------------------------------------

AliFlowCommonHist::AliFlowCommonHist():TNamed(),
  fHistMultOrig(NULL),
  fHistMultInt(NULL),
  fHistMultDiff(NULL),
  fHistPtInt(NULL),
  fHistPtDiff(NULL),
  fHistPhiInt(NULL),
  fHistPhiDiff(NULL),
  fHistEtaInt(NULL),
  fHistEtaDiff(NULL),
  fHistProMeanPtperBin(NULL),
  fHistQ(NULL),
  fHistList(NULL)
{
  
  //default constructor
  
}

AliFlowCommonHist::AliFlowCommonHist(const AliFlowCommonHist& a):
  TNamed(),
  fHistMultOrig(new TH1F(*a.fHistMultOrig)),
  fHistMultInt(new TH1F(*a.fHistMultInt)),
  fHistMultDiff(new TH1F(*a.fHistMultDiff)),
  fHistPtInt(new TH1F(*a.fHistPtInt)),
  fHistPtDiff(new TH1F(*a.fHistPtDiff)),
  fHistPhiInt(new TH1F(*a.fHistPhiInt)),
  fHistPhiDiff(new TH1F(*a.fHistPhiDiff)),
  fHistEtaInt(new TH1F(*a.fHistEtaInt)),
  fHistEtaDiff(new TH1F(*a.fHistEtaDiff)),
  fHistProMeanPtperBin(new TProfile(*a.fHistProMeanPtperBin)),
  fHistQ(new TH1F(*a.fHistQ)),
  fHistList(NULL)
{
  // copy constructor

  fHistList = new TList();
  fHistList-> Add(fHistMultOrig);        
  fHistList-> Add(fHistMultInt);        
  fHistList-> Add(fHistMultDiff);       
  fHistList-> Add(fHistPtInt);          
  fHistList-> Add(fHistPtDiff);         
  fHistList-> Add(fHistPhiInt);          
  fHistList-> Add(fHistPhiDiff);         
  fHistList-> Add(fHistEtaInt);          
  fHistList-> Add(fHistEtaDiff);         
  fHistList-> Add(fHistProMeanPtperBin); 
  fHistList-> Add(fHistQ);           
  //  TListIter next = TListIter(a.fHistList);

}

// AliFlowCommonHist& AliFlowCommonHist::operator=(const AliFlowCommonHist& a) 
// {
//   *fHistMultOrig = *a.fHistMultOrig;
//   *fHistMultInt = *a.fHistMultInt;
//   *fHistMultDiff = *a.fHistMultDiff;
//   *fHistPtInt = *a.fHistPtInt;
//   *fHistPtDiff = *a.fHistPtDiff;
//   *fHistPhiInt = *a.fHistPhiInt;
//   *fHistPhiDiff = *a.fHistPhiDiff;
//   *fHistEtaInt = *a.fHistEtaInt;
//   *fHistEtaDiff = *a.fHistEtaDiff;
//   *fHistProMeanPtperBin = *a.fHistProMeanPtperBin;
//   *fHistQ = *a.fHistQ;
//   //  *fHistList = *a.fHistList;
//   fHistList = NULL;
 
//   return *this;
// }

//-----------------------------------------------------------------------

  AliFlowCommonHist::AliFlowCommonHist(const char *anInput,const char *title):TNamed(anInput,title),
   fHistMultOrig(NULL),
   fHistMultInt(NULL),
   fHistMultDiff(NULL),
   fHistPtInt(NULL),
   fHistPtDiff(NULL),
   fHistPhiInt(NULL),
   fHistPhiDiff(NULL),
   fHistEtaInt(NULL),
   fHistEtaDiff(NULL),
   fHistProMeanPtperBin(NULL),
   fHistQ(NULL),
   fHistList(NULL)
 {

  //constructor creating histograms 
  Int_t iNbinsMult = AliFlowCommonConstants::GetNbinsMult();
  Int_t iNbinsPt = AliFlowCommonConstants::GetNbinsPt();
  Int_t iNbinsPhi = AliFlowCommonConstants::GetNbinsPhi();
  Int_t iNbinsEta = AliFlowCommonConstants::GetNbinsEta();
  Int_t iNbinsQ = AliFlowCommonConstants::GetNbinsQ();
  TString sName;

  Double_t  dMultMin = AliFlowCommonConstants::GetMultMin();            
  Double_t  dMultMax = AliFlowCommonConstants::GetMultMax();
  Double_t  dPtMin = AliFlowCommonConstants::GetPtMin();	     
  Double_t  dPtMax = AliFlowCommonConstants::GetPtMax();
  Double_t  dPhiMin = AliFlowCommonConstants::GetPhiMin();	     
  Double_t  dPhiMax = AliFlowCommonConstants::GetPhiMax();
  Double_t  dEtaMin = AliFlowCommonConstants::GetEtaMin();	     
  Double_t  dEtaMax = AliFlowCommonConstants::GetEtaMax();	     
  Double_t  dQMin = AliFlowCommonConstants::GetQMin();	     
  Double_t  dQMax = AliFlowCommonConstants::GetQMax();	
  
  cout<<"The settings for the common histograms are as follows:"<<endl;
  cout<<"Multiplicity: "<<iNbinsMult<<" bins between "<<dMultMin<<" and "<<dMultMax<<endl;
  cout<<"Pt: "<<iNbinsPt<<" bins between "<<dPtMin<<" and "<<dPtMax<<endl;
  cout<<"Phi: "<<iNbinsPhi<<" bins between "<<dPhiMin<<" and "<<dPhiMax<<endl;
  cout<<"Eta: "<<iNbinsEta<<" bins between "<<dEtaMin<<" and "<<dEtaMax<<endl;
  cout<<"Q: "<<iNbinsQ<<" bins between "<<dQMin<<" and "<<dQMax<<endl;

  //Multiplicity
  sName = "Control_Flow_OrigMult_";
  sName +=anInput;
  fHistMultOrig = new TH1F(sName.Data(), sName.Data(),iNbinsMult, dMultMin, dMultMax);
  fHistMultOrig ->SetXTitle("Original Multiplicity");
  fHistMultOrig ->SetYTitle("Counts");

  sName = "Control_Flow_MultInt_";
  sName +=anInput;
  fHistMultInt = new TH1F(sName.Data(), sName.Data(),iNbinsMult, dMultMin, dMultMax);
  fHistMultInt ->SetXTitle("Multiplicity for integrated flow");
  fHistMultInt ->SetYTitle("Counts");

  sName = "Control_Flow_MultDiff_";
  sName +=anInput;
  fHistMultDiff = new TH1F(sName.Data(), sName.Data(),iNbinsMult, dMultMin, dMultMax);
  fHistMultDiff ->SetXTitle("Multiplicity for differential flow");
  fHistMultDiff ->SetYTitle("Counts");

  //Pt
  sName = "Control_Flow_PtInt_";
  sName +=anInput;
  fHistPtInt = new TH1F(sName.Data(), sName.Data(),iNbinsPt, dPtMin, dPtMax); 
  fHistPtInt ->SetXTitle("P_{t} (GeV/c) for integrated flow");
  fHistPtInt ->SetYTitle("Counts");

  sName = "Control_Flow_PtDiff_";
  sName +=anInput;
  fHistPtDiff = new TH1F(sName.Data(), sName.Data(),iNbinsPt, dPtMin, dPtMax); 
  //binning has to be the same as for fHistProVPt! use to get Nprime!
  fHistPtDiff ->SetXTitle("P_{t} (GeV/c) for differential flow");
  fHistPtDiff ->SetYTitle("Counts");

  //Phi
  sName = "Control_Flow_PhiInt_";
  sName +=anInput;
  fHistPhiInt = new TH1F(sName.Data(), sName.Data(),iNbinsPhi, dPhiMin, dPhiMax);
  fHistPhiInt ->SetXTitle("#phi for integrated flow");
  fHistPhiInt ->SetYTitle("Counts");

  sName = "Control_Flow_PhiDiff_";
  sName +=anInput;
  fHistPhiDiff = new TH1F(sName.Data(), sName.Data(),iNbinsPhi, dPhiMin, dPhiMax);
  fHistPhiDiff ->SetXTitle("#phi for differential flow");
  fHistPhiDiff ->SetYTitle("Counts");

  //Eta
  sName = "Control_Flow_EtaInt_";
  sName +=anInput;
  fHistEtaInt = new TH1F(sName.Data(), sName.Data(),iNbinsEta, dEtaMin, dEtaMax);
  fHistEtaInt ->SetXTitle("#eta for integrated flow");
  fHistEtaInt ->SetYTitle("Counts");

  sName = "Control_Flow_EtaDiff_";
  sName +=anInput;
  fHistEtaDiff = new TH1F(sName.Data(), sName.Data(),iNbinsEta, dEtaMin, dEtaMax);
  fHistEtaDiff ->SetXTitle("#eta for differential flow");
  fHistEtaDiff ->SetYTitle("Counts");

  //Mean Pt per pt bin 
  sName = "Control_FlowPro_MeanPtperBin_";
  sName +=anInput;
  fHistProMeanPtperBin = new TProfile(sName.Data(), sName.Data(),iNbinsPt,dPtMin,dPtMax);
  fHistProMeanPtperBin ->SetXTitle("P_{t} (GeV/c) ");
  fHistProMeanPtperBin ->SetYTitle("<P_{t}> (GeV/c) ");

  //Q vector
  sName = "Control_Flow_Q_";
  sName +=anInput;
  fHistQ = new TH1F(sName.Data(), sName.Data(),iNbinsQ, dQMin, dQMax);
  fHistQ ->SetXTitle("Q_{vector}/Mult");
  fHistQ ->SetYTitle("Counts");  

  //list of histograms if added here also add in copy constructor
  fHistList = new TList();
  fHistList-> Add(fHistMultOrig);        
  fHistList-> Add(fHistMultInt);        
  fHistList-> Add(fHistMultDiff);       
  fHistList-> Add(fHistPtInt);          
  fHistList-> Add(fHistPtDiff);         
  fHistList-> Add(fHistPhiInt);          
  fHistList-> Add(fHistPhiDiff);         
  fHistList-> Add(fHistEtaInt);          
  fHistList-> Add(fHistEtaDiff);         
  fHistList-> Add(fHistProMeanPtperBin); 
  fHistList-> Add(fHistQ);           



}


//----------------------------------------------------------------------- 

AliFlowCommonHist::~AliFlowCommonHist()
{
  //deletes histograms
  delete fHistMultOrig;      
  delete fHistMultInt;       
  delete fHistMultDiff;      
  delete fHistPtInt;         
  delete fHistPtDiff;       
  delete fHistPhiInt;        
  delete fHistPhiDiff;       
  delete fHistEtaInt;        
  delete fHistEtaDiff;
  delete fHistProMeanPtperBin;
  delete fHistQ;
  delete fHistList;
}


//----------------------------------------------------------------------- 

Bool_t AliFlowCommonHist::FillControlHistograms(AliFlowEventSimple* anEvent)
{
  //Fills the control histograms
  if (!anEvent){
    cout<<"##### FillControlHistograms: FlowEvent pointer null"<<endl;
    return kFALSE;
  }

  Double_t dPt, dPhi, dEta;


  //fill the histograms
  Int_t iNumberOfTracks = anEvent->NumberOfTracks();
  fHistMultOrig->Fill(iNumberOfTracks);

  AliFlowVector vQ = anEvent->GetQ(); 
  //weight by the Multiplicity
  Double_t dQX = 0.;
  Double_t dQY = 0.;
  if (vQ.GetMult()!=0) {
    dQX = vQ.X()/vQ.GetMult();
    dQY = vQ.Y()/vQ.GetMult();
  }
  vQ.Set(dQX,dQY);
  fHistQ->Fill(vQ.Mod());

  Int_t iMultInt = 0;
  Int_t iMultDiff = 0;
  
  AliFlowTrackSimple* pTrack = NULL;     

  for (Int_t i=0;i<iNumberOfTracks;i++) {
    pTrack = anEvent->GetTrack(i);
    if (pTrack ) {
      if (pTrack->UseForIntegratedFlow()){
	dPt = pTrack->Pt();
	fHistPtInt->Fill(dPt);
	dPhi = pTrack->Phi();
	if (dPhi<0.) dPhi+=2*TMath::Pi();
	fHistPhiInt->Fill(dPhi);
	dEta = pTrack->Eta();
	fHistEtaInt->Fill(dEta);
	iMultInt++;
      }
      if (pTrack->UseForDifferentialFlow()){
	dPt = pTrack->Pt();
	fHistPtDiff->Fill(dPt);
	dPhi = pTrack->Phi();
	if (dPhi<0.) dPhi+=2*TMath::Pi();
	fHistPhiDiff->Fill(dPhi);
	dEta = pTrack->Eta();
	fHistEtaDiff->Fill(dEta);
	fHistProMeanPtperBin->Fill(dPt,dPt);
	iMultDiff++;
      }
    } //track
  } //loop over tracks
  
  fHistMultInt->Fill(iMultInt);
  fHistMultDiff->Fill(iMultDiff);

  return kTRUE; 
}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetEntriesInPtBinRP(Int_t aBin)
{
  //get entries in bin aBin from fHistPtDiff
  Double_t dEntries = fHistPtInt->GetBinContent(aBin);

  return dEntries;

}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetEntriesInPtBinPOI(Int_t aBin)
{
  //get entries in bin aBin from fHistPtDiff
  Double_t dEntries = fHistPtDiff->GetBinContent(aBin);

  return dEntries;

}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetEntriesInEtaBinRP(Int_t aBin)
{
  //get entries in bin aBin from fHistPtDiff
  Double_t dEntries = fHistEtaInt->GetBinContent(aBin);

  return dEntries;

}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetEntriesInEtaBinPOI(Int_t aBin)
{
  //get entries in bin aBin from fHistPtDiff
  Double_t dEntries = fHistEtaDiff->GetBinContent(aBin);

  return dEntries;

}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetMeanPt(Int_t aBin)
{  
  //Get entry from bin aBin from fHistProMeanPtperBin
  Double_t dMeanPt = fHistProMeanPtperBin->GetBinContent(aBin);

  return dMeanPt;
  
}


//----------------------------------------------------------------------- 
 Double_t AliFlowCommonHist::Merge(TCollection *aList)
{
  //merge fuction
  cout<<"entering merge function"<<endl;
  if (!aList) return 0;
  if (aList->IsEmpty()) return 0; //no merging is needed

  Int_t iCount = 0;
  TIter next(aList); // list is supposed to contain only objects of the same type as this
  AliFlowCommonHist *toMerge;
  // make a temporary list
  TList *pTemp = new TList();
  while ((toMerge=(AliFlowCommonHist*)next())) {
    pTemp->Add(toMerge->GetHistList()); 
    iCount++;
  }
  // Now call merge for fHistList providing temp list
  fHistList->Merge(pTemp);
  // Cleanup
  delete pTemp;
    
  cout<<"Merged"<<endl;
  return (double)iCount;
    
}

void AliFlowCommonHist::Print(Option_t *option) const
{
  //   -*-*-*-*-*Print some global quantities for this histogram collection class *-*-*-*-*-*-*-*
  //             ===============================================
  //   printf( "TH1.Print Name  = %s, Entries= %d, Total sum= %g\n",GetName(),Int_t(fEntries),GetSumOfWeights());
  printf( "Class.Print Name = %s, Histogram list:\n",GetName());

  if (fHistList) {  
    fHistList->Print(option);
  }
  else
    {
      printf( "Empty histogram list \n");
    }
}

//----------------------------------------------------------------------- 
 void AliFlowCommonHist::Browse(TBrowser *b)
{

  if (!b) return;
  if (fHistList) b->Add(fHistList,"AliFlowCommonHistList");
}




