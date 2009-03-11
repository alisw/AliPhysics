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
  fHistMultRP(NULL),
  fHistMultPOI(NULL),
  fHistPtRP(NULL),
  fHistPtPOI(NULL),
  fHistPhiRP(NULL),
  fHistPhiPOI(NULL),
  fHistEtaRP(NULL),
  fHistEtaPOI(NULL),
  fHistProMeanPtperBin(NULL),
  fHistQ(NULL),
  fHistList(NULL)
{
  
  //default constructor
  
}

AliFlowCommonHist::AliFlowCommonHist(const AliFlowCommonHist& a):
  TNamed(),
  fHistMultOrig(new TH1F(*a.fHistMultOrig)),
  fHistMultRP(new TH1F(*a.fHistMultRP)),
  fHistMultPOI(new TH1F(*a.fHistMultPOI)),
  fHistPtRP(new TH1F(*a.fHistPtRP)),
  fHistPtPOI(new TH1F(*a.fHistPtPOI)),
  fHistPhiRP(new TH1F(*a.fHistPhiRP)),
  fHistPhiPOI(new TH1F(*a.fHistPhiPOI)),
  fHistEtaRP(new TH1F(*a.fHistEtaRP)),
  fHistEtaPOI(new TH1F(*a.fHistEtaPOI)),
  fHistProMeanPtperBin(new TProfile(*a.fHistProMeanPtperBin)),
  fHistQ(new TH1F(*a.fHistQ)),
  fHistList(NULL)
{
  // copy constructor

  fHistList = new TList();
  fHistList-> Add(fHistMultOrig);        
  fHistList-> Add(fHistMultRP);        
  fHistList-> Add(fHistMultPOI);       
  fHistList-> Add(fHistPtRP);          
  fHistList-> Add(fHistPtPOI);         
  fHistList-> Add(fHistPhiRP);          
  fHistList-> Add(fHistPhiPOI);         
  fHistList-> Add(fHistEtaRP);          
  fHistList-> Add(fHistEtaPOI);         
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
   fHistMultRP(NULL),
   fHistMultPOI(NULL),
   fHistPtRP(NULL),
   fHistPtPOI(NULL),
   fHistPhiRP(NULL),
   fHistPhiPOI(NULL),
   fHistEtaRP(NULL),
   fHistEtaPOI(NULL),
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

  sName = "Control_Flow_MultRP_";
  sName +=anInput;
  fHistMultRP = new TH1F(sName.Data(), sName.Data(),iNbinsMult, dMultMin, dMultMax);
  fHistMultRP ->SetXTitle("Multiplicity for RP selection");
  fHistMultRP ->SetYTitle("Counts");

  sName = "Control_Flow_MultPOI_";
  sName +=anInput;
  fHistMultPOI = new TH1F(sName.Data(), sName.Data(),iNbinsMult, dMultMin, dMultMax);
  fHistMultPOI ->SetXTitle("Multiplicity for POI selection");
  fHistMultPOI ->SetYTitle("Counts");

  //Pt
  sName = "Control_Flow_PtRP_";
  sName +=anInput;
  fHistPtRP = new TH1F(sName.Data(), sName.Data(),iNbinsPt, dPtMin, dPtMax); 
  fHistPtRP ->SetXTitle("P_{t} (GeV/c) for RP selection");
  fHistPtRP ->SetYTitle("Counts");

  sName = "Control_Flow_PtPOI_";
  sName +=anInput;
  fHistPtPOI = new TH1F(sName.Data(), sName.Data(),iNbinsPt, dPtMin, dPtMax); 
  //binning has to be the same as for fHistProVPt! use to get Nprime!
  fHistPtPOI ->SetXTitle("P_{t} (GeV/c) for POI selection");
  fHistPtPOI ->SetYTitle("Counts");

  //Phi
  sName = "Control_Flow_PhiRP_";
  sName +=anInput;
  fHistPhiRP = new TH1F(sName.Data(), sName.Data(),iNbinsPhi, dPhiMin, dPhiMax);
  fHistPhiRP ->SetXTitle("#phi for RP selection");
  fHistPhiRP ->SetYTitle("Counts");

  sName = "Control_Flow_PhiPOI_";
  sName +=anInput;
  fHistPhiPOI = new TH1F(sName.Data(), sName.Data(),iNbinsPhi, dPhiMin, dPhiMax);
  fHistPhiPOI ->SetXTitle("#phi for POI selection");
  fHistPhiPOI ->SetYTitle("Counts");

  //Eta
  sName = "Control_Flow_EtaRP_";
  sName +=anInput;
  fHistEtaRP = new TH1F(sName.Data(), sName.Data(),iNbinsEta, dEtaMin, dEtaMax);
  fHistEtaRP ->SetXTitle("#eta for RP selection");
  fHistEtaRP ->SetYTitle("Counts");

  sName = "Control_Flow_EtaPOI_";
  sName +=anInput;
  fHistEtaPOI = new TH1F(sName.Data(), sName.Data(),iNbinsEta, dEtaMin, dEtaMax);
  fHistEtaPOI ->SetXTitle("#eta for POI selection");
  fHistEtaPOI ->SetYTitle("Counts");

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
  fHistList-> Add(fHistMultRP);        
  fHistList-> Add(fHistMultPOI);       
  fHistList-> Add(fHistPtRP);          
  fHistList-> Add(fHistPtPOI);         
  fHistList-> Add(fHistPhiRP);          
  fHistList-> Add(fHistPhiPOI);         
  fHistList-> Add(fHistEtaRP);          
  fHistList-> Add(fHistEtaPOI);         
  fHistList-> Add(fHistProMeanPtperBin); 
  fHistList-> Add(fHistQ);           



}


//----------------------------------------------------------------------- 

AliFlowCommonHist::~AliFlowCommonHist()
{
  //deletes histograms
  delete fHistMultOrig;      
  delete fHistMultRP;       
  delete fHistMultPOI;      
  delete fHistPtRP;         
  delete fHistPtPOI;       
  delete fHistPhiRP;        
  delete fHistPhiPOI;       
  delete fHistEtaRP;        
  delete fHistEtaPOI;
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

  Int_t iMultRP = 0;
  Int_t iMultPOI = 0;
  
  AliFlowTrackSimple* pTrack = NULL;     

  for (Int_t i=0;i<iNumberOfTracks;i++) {
    pTrack = anEvent->GetTrack(i);
    if (pTrack ) {
      if (pTrack->InRPSelection()){
	dPt = pTrack->Pt();
	fHistPtRP->Fill(dPt);
	dPhi = pTrack->Phi();
	if (dPhi<0.) dPhi+=2*TMath::Pi();
	fHistPhiRP->Fill(dPhi);
	dEta = pTrack->Eta();
	fHistEtaRP->Fill(dEta);
	iMultRP++;
      }
      if (pTrack->InPOISelection()){
	dPt = pTrack->Pt();
	fHistPtPOI->Fill(dPt);
	dPhi = pTrack->Phi();
	if (dPhi<0.) dPhi+=2*TMath::Pi();
	fHistPhiPOI->Fill(dPhi);
	dEta = pTrack->Eta();
	fHistEtaPOI->Fill(dEta);
	fHistProMeanPtperBin->Fill(dPt,dPt);
	iMultPOI++;
      }
    } //track
  } //loop over tracks
  
  fHistMultRP->Fill(iMultRP);
  fHistMultPOI->Fill(iMultPOI);

  return kTRUE; 
}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetEntriesInPtBinRP(Int_t aBin)
{
  //get entries in bin aBin from fHistPtRP
  Double_t dEntries = fHistPtRP->GetBinContent(aBin);

  return dEntries;

}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetEntriesInPtBinPOI(Int_t aBin)
{
  //get entries in bin aBin from fHistPtPOI
  Double_t dEntries = fHistPtPOI->GetBinContent(aBin);

  return dEntries;

}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetEntriesInEtaBinRP(Int_t aBin)
{
  //get entries in bin aBin from fHistPtRP
  Double_t dEntries = fHistEtaRP->GetBinContent(aBin);

  return dEntries;

}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetEntriesInEtaBinPOI(Int_t aBin)
{
  //get entries in bin aBin from fHistPtPOI
  Double_t dEntries = fHistEtaPOI->GetBinContent(aBin);

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




