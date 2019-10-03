/*
***********************************************************
    Histogram manager for event plane framework QA
    Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
    Histogram manager inspired from PWGDQ/dielectron/AliDielectronHistos by J.Wiechula
    and based on work of Ionut-Cristian Arsene
***********************************************************
*/

#include "AliQnCorrectionsHistos.h"

#include <iostream>
#include <fstream>

#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>
#include <THashList.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <THn.h>
#include <TIterator.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TAxis.h>
#include <TVector.h>
#include <TTimeStamp.h>


ClassImp(AliQnCorrectionsHistos)

Int_t gkNVars = 999;

AliQnCorrectionsHistos* AliQnCorrectionsHistos::fgEventPlaneHistos = NULL;
AliQnCorrectionsHistos::AliQnCorrectionsHistos() :
//   TCollection(),
  TNamed("AliQnCorrectionsHistos","EventPlane Histogram Container"),
  //fHistLists(0x0),   // main histogram list for the current running process
  fHistFile(0x0)      // pointer to a TFile opened for reading
{
  //
  // Default constructor
  //
  fHistLists = new TObjArray();
}

//_____________________________________________________________________________
AliQnCorrectionsHistos::AliQnCorrectionsHistos(const char* name, const char* title) :
//   TCollection(),
  TNamed(name, title),
  //fHistLists(0x0),   // main histogram list for the current running process
  fHistFile(0x0)      // pointer to a TFile opened for reading
{
  //
  // TNamed constructor
  //
  fHistLists = new TObjArray();
}

//_____________________________________________________________________________
AliQnCorrectionsHistos::~AliQnCorrectionsHistos()
{
  //
  // Destructor
  //
  fHistLists->Clear("C");
}


//__________________________________________________________________
void AliQnCorrectionsHistos::FillHistClass(const Char_t* className, Float_t* values){// Bool_t* usedvalues) {
  //
  //  fill a class of histograms
  //
  THashList* hList = (THashList*)fHistLists->FindObject(className);
  if(!hList) {
    return;
  }
  
  TIter next(hList);
  TObject* h=0x0;
  Bool_t isProfile;
  Bool_t isTHn;
  Int_t thnDim=0;
  Double_t fillValues[20]={0.0};
  Bool_t allVarsGood = kTRUE;
  Int_t uid = 0;
  Int_t varX=kNothing, varY=kNothing, varZ=kNothing, varT=kNothing, varW=kNothing;
  Int_t dimension=0;
  while((h=next())) {
    uid = h->GetUniqueID();
    isProfile = (uid%10==1 ? kTRUE : kFALSE);   // units digit encodes the isProfile
    isTHn = ((uid%100)>10 ? kTRUE : kFALSE);      
    if(isTHn) thnDim = (uid%100)-10;        // the excess over 10 from the last 2 digits give the dimension of the THn
    dimension = 0;
    if(!isTHn) dimension = ((TH1*)h)->GetDimension();
        
    uid = (uid-(uid%100))/100;
    varX = kNothing;
    varY = kNothing;
    varZ = kNothing;
    varT = kNothing;
    varW = kNothing;
    if(uid>0) {
      varW = uid%(gkNVars+1)-1;
      if(varW==0) varW=kNothing;
      uid = (uid-(uid%(gkNVars+1)))/(gkNVars+1);
      if(uid>0) varT = uid - 1;
      //cout << "Filling " << h->GetName() << " with varT " << varT << endl;
    }
        
    if(!isTHn) {
      varX = ((TH1*)h)->GetXaxis()->GetUniqueID();
      //if(usedvalues[varX]) {
        switch(dimension) {
          case 1:
            if(isProfile) {
              varY = ((TH1*)h)->GetYaxis()->GetUniqueID();
              //if(!usedvalues[varY]) break;
              if(varW!=kNothing) { 
                //if(!usedvalues[varW]) break;
                ((TProfile*)h)->Fill(values[varX],values[varY],values[varW]);
              }
              else 
                ((TProfile*)h)->Fill(values[varX],values[varY]);
            }
            else {
              if(varW!=kNothing) {
                //if(!usedvalues[varW]) break;
                ((TH1F*)h)->Fill(values[varX],values[varW]);
                //if(!usedvalues[varW]) break;
              }
              else
                ((TH1F*)h)->Fill(values[varX]);
            }
          break;
          case 2:
            varY = ((TH1*)h)->GetYaxis()->GetUniqueID();
            //if(!usedvalues[varY]) break;
            if(isProfile) {
              varZ = ((TH1*)h)->GetZaxis()->GetUniqueID();
              //if(!usedvalues[varZ]) break;
              if(varW!=kNothing) {
                //if(!usedvalues[varW]) break;
                ((TProfile2D*)h)->Fill(values[varX],values[varY],values[varZ],values[varW]);
              }
              else
                ((TProfile2D*)h)->Fill(values[varX],values[varY],values[varZ]);        
            }
            else {
              if(varW!=kNothing) {
                //if(!usedvalues[varW]) break;
                ((TH2F*)h)->Fill(values[varX],values[varY], values[varW]);
              }
              else
                ((TH2F*)h)->Fill(values[varX],values[varY]);
            }
          break;
          case 3:
            varY = ((TH1*)h)->GetYaxis()->GetUniqueID();
            //if(!usedvalues[varY]) break;
            varZ = ((TH1*)h)->GetZaxis()->GetUniqueID();
            //if(!usedvalues[varZ]) break;
            if(isProfile) {
              //if(!usedvalues[varT]) break;
              if(varW!=kNothing) {
                //if(!usedvalues[varW]) break;
                ((TProfile3D*)h)->Fill(values[varX],values[varY],values[varZ],values[varT],values[varW]);
              }
              else
                ((TProfile3D*)h)->Fill(values[varX],values[varY],values[varZ],values[varT]);
            }
            else {
              if(varW!=kNothing) {
                //if(!usedvalues[varW]) break;
                ((TH3F*)h)->Fill(values[varX],values[varY],values[varZ],values[varW]);
              }
              else
                ((TH3F*)h)->Fill(values[varX],values[varY],values[varZ]);
            }
          break;
          default:
          break;
        }  // end switch
      //}
    }  // end if(!isTHn)
    else {
      for(Int_t idim=0;idim<thnDim;++idim) {
        //allVarsGood &= usedvalues[((THnF*)h)->GetAxis(idim)->GetUniqueID()];
        fillValues[idim] = values[((THnF*)h)->GetAxis(idim)->GetUniqueID()];
      }
      if(allVarsGood) {
        if(varW!=kNothing) {
          //if(usedvalues[varW])
            ((THnF*)h)->Fill(fillValues,values[varW]);
        }
        else
          ((THnF*)h)->Fill(fillValues);
      }
    }
  }
}


//__________________________________________________________________
void AliQnCorrectionsHistos::AddHistClass(const Char_t* histClass) {
  //
  // Add a new histogram list
  //
  if(!fHistLists) {
    fHistLists = new TObjArray();
    fHistLists->SetOwner();
    fHistLists->SetName("histos");
  }
  
  if(fHistLists->FindObject(histClass)) {
    //cout << "Warning in AddHistClass: Cannot add histogram class " << histClass
    //     << " because it already exists." << endl;
    return;
  }
  THashList* hList=new THashList;
  hList->SetOwner(kTRUE);
  hList->SetName(histClass);
  fHistLists->Add(hList);
}

//_________________________________________________________________
void AliQnCorrectionsHistos::AddHistogram(const Char_t* histClass,
		                                   const Char_t* name, const Char_t* title, Bool_t isProfile,
                                       Int_t nXbins, Double_t xmin, Double_t xmax, Int_t varX,
		                                   Int_t nYbins, Double_t ymin, Double_t ymax, Int_t varY,
		                                   Int_t nZbins, Double_t zmin, Double_t zmax, Int_t varZ,
                                       const Char_t* xLabels, const Char_t* yLabels, const Char_t* zLabels,
                                       Int_t varT, Int_t varW) {
  //
  // add a histogram
  //
  THashList* hList = (THashList*)fHistLists->FindObject(histClass);
  if(hList->FindObject(name)) {
    //cout << "Warning in AddHistogram(): Histogram " << name << " already exists" << endl;
    return;
  }
  TString hname = name;
  
  Int_t dimension = 1;
  if(varY!=kNothing) dimension = 2;
  if(varZ!=kNothing) dimension = 3;
  
  TString titleStr(title);
  TObjArray* arr=titleStr.Tokenize(";");
  //if(varT!=kNothing) usedvalues[varT] = kTRUE;
  //if(varW!=kNothing) usedvalues[varW] = kTRUE;
  
  TH1* h=0x0;
  switch(dimension) {
    case 1:
      h=new TH1F(hname.Data(),arr->At(0)->GetName(),nXbins,xmin,xmax);
      fBinsAllocated+=nXbins+2;
      h->Sumw2();
      h->SetUniqueID(0);
      if(varW!=kNothing) h->SetUniqueID(100*(varW+1)+0); 
      h->GetXaxis()->SetUniqueID(UInt_t(varX));
      if(arr->At(1)) h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      if(xLabels[0]!='\0') MakeAxisLabels(h->GetXaxis(), xLabels);
      //usedvalues[varX] = kTRUE;
      hList->Add(h);
      break;
    case 2:
      if(isProfile) {
	h=new TProfile(hname.Data(),arr->At(0)->GetName(),nXbins,xmin,xmax);
        fBinsAllocated+=nXbins+2;
        h->SetUniqueID(1);
        if(varW!=kNothing) h->SetUniqueID(100*(varW+1)+1);
      }
      else {
	h=new TH2F(hname.Data(),arr->At(0)->GetName(),nXbins,xmin,xmax,nYbins,ymin,ymax);
        fBinsAllocated+=(nXbins+2)*(nYbins+2);
        h->Sumw2();
        h->SetUniqueID(0);
        if(varW!=kNothing) h->SetUniqueID(100*(varW+1)+0); 
      }
      h->GetXaxis()->SetUniqueID(UInt_t(varX));
      h->GetYaxis()->SetUniqueID(UInt_t(varY));
      if(arr->At(1)) h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      if(xLabels[0]!='\0') MakeAxisLabels(h->GetXaxis(), xLabels);
      if(arr->At(2)) h->GetYaxis()->SetTitle(arr->At(2)->GetName());
      if(yLabels[0]!='\0') MakeAxisLabels(h->GetYaxis(), yLabels);
      //usedvalues[varX] = kTRUE;
      //usedvalues[varY] = kTRUE;
      hList->Add(h);
      break;
    case 3:
      if(isProfile) {
        if(varT!=kNothing) {
          h=new TProfile3D(hname.Data(),arr->At(0)->GetName(),nXbins,xmin,xmax,nYbins,ymin,ymax,nZbins,zmin,zmax);
          fBinsAllocated+=(nXbins+2)*(nYbins+2)*(nZbins+2);
          if(varW!=kNothing) h->SetUniqueID(((varW+1)+(gkNVars+1)*(varT+1))*100+1);   // 4th variable "varT" is encoded in the UniqueId of the histogram
          else h->SetUniqueID((gkNVars+1)*(varT+1)*100+1);
        }
        else {
	  h=new TProfile2D(hname.Data(),arr->At(0)->GetName(),nXbins,xmin,xmax,nYbins,ymin,ymax);
          fBinsAllocated+=(nXbins+2)*(nYbins+2);
          h->SetUniqueID(1);
          if(varW!=kNothing) h->SetUniqueID(100*(varW+1)+1); 
        }
      }
      else {
	h=new TH3F(hname.Data(),arr->At(0)->GetName(),nXbins,xmin,xmax,nYbins,ymin,ymax,nZbins,zmin,zmax);
        fBinsAllocated+=(nXbins+2)*(nYbins+2)*(nZbins+2);
        h->Sumw2();
        h->SetUniqueID(0);
        if(varW!=kNothing) h->SetUniqueID(100*(varW+1)+0); 
      }
      h->GetXaxis()->SetUniqueID(UInt_t(varX));
      h->GetYaxis()->SetUniqueID(UInt_t(varY));
      h->GetZaxis()->SetUniqueID(UInt_t(varZ));
      if(arr->At(1)) h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      if(xLabels[0]!='\0') MakeAxisLabels(h->GetXaxis(), xLabels);
      if(arr->At(2)) h->GetYaxis()->SetTitle(arr->At(2)->GetName());
      if(yLabels[0]!='\0') MakeAxisLabels(h->GetYaxis(), yLabels);
      if(arr->At(3)) h->GetZaxis()->SetTitle(arr->At(3)->GetName());
      if(zLabels[0]!='\0') MakeAxisLabels(h->GetZaxis(), zLabels);
      //usedvalues[varX] = kTRUE;
      //usedvalues[varY] = kTRUE;
      //usedvalues[varZ] = kTRUE;
      hList->Add(h);
      break;
  }
  h->SetDirectory(0);
}

//_________________________________________________________________
void AliQnCorrectionsHistos::AddHistogram(const Char_t* histClass,
		                   const Char_t* name, const Char_t* title, Bool_t isProfile,
                                   Int_t nXbins, Double_t* xbins, Int_t varX,
		                   Int_t nYbins, Double_t* ybins, Int_t varY,
		                   Int_t nZbins, Double_t* zbins, Int_t varZ,
		                   const Char_t* xLabels, const Char_t* yLabels, const Char_t* zLabels,
                                   Int_t varT, Int_t varW) {
  //
  // add a histogram
  //
  THashList* hList = (THashList*)fHistLists->FindObject(histClass);
  if(hList->FindObject(name)) {
    //cout << "Warning in AddHistogram(): Histogram " << name << " already exists" << endl;
    return;
  }
  TString hname = name;
  
  Int_t dimension = 1;
  if(varY!=kNothing) dimension = 2;
  if(varZ!=kNothing) dimension = 3;
  
  
  TString titleStr(title);
  TObjArray* arr=titleStr.Tokenize(";");
  
  TH1* h=0x0;
  switch(dimension) {
    case 1:
      h=new TH1F(hname.Data(),arr->At(0)->GetName(),nXbins,xbins);
      fBinsAllocated+=nXbins+2;
      h->Sumw2();
      h->SetUniqueID(0);
      if(varW!=kNothing) h->SetUniqueID(100*(varW+1)+0); 
      h->GetXaxis()->SetUniqueID(UInt_t(varX));
      if(arr->At(1)) h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      if(xLabels[0]!='\0') MakeAxisLabels(h->GetXaxis(), xLabels);
      hList->Add(h);
      break;
    case 2:
      if(isProfile) {
	h=new TProfile(hname.Data(),arr->At(0)->GetName(),nXbins,xbins);
        fBinsAllocated+=nXbins+2;
        h->SetUniqueID(1);
        if(varW!=kNothing) h->SetUniqueID(100*(varW+1)+1); 
      }
      else {
	h=new TH2F(hname.Data(),arr->At(0)->GetName(),nXbins,xbins,nYbins,ybins);
        fBinsAllocated+=(nXbins+2)*(nYbins+2);
        h->Sumw2();
        h->SetUniqueID(0);
        if(varW!=kNothing) h->SetUniqueID(100*(varW+1)+0);
      }
      h->GetXaxis()->SetUniqueID(UInt_t(varX));
      h->GetYaxis()->SetUniqueID(UInt_t(varY));
      if(arr->At(1)) h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      if(xLabels[0]!='\0') MakeAxisLabels(h->GetXaxis(), xLabels);
      if(arr->At(2)) h->GetYaxis()->SetTitle(arr->At(2)->GetName());
      if(yLabels[0]!='\0') MakeAxisLabels(h->GetYaxis(), yLabels);
      hList->Add(h);
      break;
    case 3:
      if(isProfile) {
        if(varT!=kNothing) {
          h=new TProfile3D(hname.Data(),arr->At(0)->GetName(),nXbins,xbins,nYbins,ybins,nZbins,zbins);
          fBinsAllocated+=(nXbins+2)*(nYbins+2)*(nZbins+2);
          if(varW!=kNothing) h->SetUniqueID(((varW+1)+(gkNVars+1)*(varT+1))*100+1);   // 4th variable "varT" is encoded in the UniqueId of the histogram
          else h->SetUniqueID((gkNVars+1)*(varT+1)*100+1);
        }
        else {
	  h=new TProfile2D(hname.Data(),arr->At(0)->GetName(),nXbins,xbins,nYbins,ybins);
          fBinsAllocated+=(nXbins+2)*(nYbins+2);
          h->SetUniqueID(1);
          if(varW!=kNothing) h->SetUniqueID(100*(varW+1)+1);
        }
      }
      else {
	h=new TH3F(hname.Data(),arr->At(0)->GetName(),nXbins,xbins,nYbins,ybins,nZbins,zbins);
        fBinsAllocated+=(nXbins+2)*(nYbins+2)*(nZbins+2);
        h->Sumw2();
        h->SetUniqueID(0);
        if(varW!=kNothing) h->SetUniqueID(100*(varW+1)+0);
      }
      h->GetXaxis()->SetUniqueID(UInt_t(varX));
      h->GetYaxis()->SetUniqueID(UInt_t(varY));
      h->GetZaxis()->SetUniqueID(UInt_t(varZ));
      if(arr->At(1)) h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      if(xLabels[0]!='\0') MakeAxisLabels(h->GetXaxis(), xLabels);
      if(arr->At(2)) h->GetYaxis()->SetTitle(arr->At(2)->GetName());
      if(yLabels[0]!='\0') MakeAxisLabels(h->GetYaxis(), yLabels);
      if(arr->At(3)) h->GetZaxis()->SetTitle(arr->At(3)->GetName());
      if(zLabels[0]!='\0') MakeAxisLabels(h->GetZaxis(), zLabels);
      hList->Add(h);
      break;
  }
  h->SetDirectory(0);
}


//_________________________________________________________________
void AliQnCorrectionsHistos::AddHistogram(const Char_t* histClass,
                                   const Char_t* name, const Char_t* title, 
                                   Int_t nDimensions, Int_t* vars,
                                   Int_t* nBins, Double_t* xmin, Double_t* xmax,
                                   TString* axLabels, Int_t varW,
				                           Int_t axisA, Double_t * newbinsA,
		    		                       Int_t axisB, Double_t * newbinsB,
				                           Int_t axisC, Double_t * newbinsC) {
  //
  // add a multi-dimensional histogram THnF
  //
  THashList* hList = (THashList*)fHistLists->FindObject(histClass);
  if(hList->FindObject(name)) {
    //cout << "Warning in AddHistogram(): Histogram " << name << " already exists" << endl;
    return;
  }
  TString hname = name;
  
  TString titleStr(title);
  TObjArray* arr=titleStr.Tokenize(";");
  
  
  THnF* h=new THnF(hname.Data(),arr->At(0)->GetName(),nDimensions,nBins,xmin,xmax);
  if(axisA!=-1) {TAxis *axisrebinA = h->GetAxis(axisA);
  axisrebinA->Set(nBins[axisA], newbinsA);}
  if(axisB!=-1) {TAxis *axisrebinB = h->GetAxis(axisB);
  axisrebinB->Set(nBins[axisB], newbinsB);}
  if(axisC!=-1) {TAxis *axisrebinC = h->GetAxis(axisC);
  axisrebinC->Set(nBins[axisC], newbinsC);}
  h->Sumw2();
  if(varW!=kNothing) h->SetUniqueID(10+nDimensions+100*(varW+1));
  else h->SetUniqueID(10+nDimensions);
  ULong_t bins = 1;
  for(Int_t idim=0;idim<nDimensions;++idim) {
    bins*=(nBins[idim]+2);
    TAxis* axis = h->GetAxis(idim);
    axis->SetUniqueID(vars[idim]);
    if(arr->At(1+idim)) axis->SetTitle(arr->At(1+idim)->GetName());
    if(axLabels && !axLabels[idim].IsNull()) 
      MakeAxisLabels(axis, axLabels[idim].Data());
  }
  hList->Add(h);
  fBinsAllocated+=bins;
}

//_________________________________________________________________
void AliQnCorrectionsHistos::AddHistogram(const Char_t* histClass,
                                   const Char_t* name, const Char_t* title,
                                   Int_t nDimensions, Int_t* vars,
                                   TArrayD* binLimits,
                                   TString* axLabels,
                                   Int_t varW) {
  //
  // add a multi-dimensional histogram THnF with equal or variable bin widths
  //
  THashList* hList = (THashList*)fHistLists->FindObject(histClass);
  if(hList->FindObject(name)) {
    //cout << "Warning in AddHistogram(): Histogram " << name << " already exists" << endl;
    return;
  }
  TString hname = name;

  TString titleStr(title);
  TObjArray* arr=titleStr.Tokenize(";");


  Double_t* xmin = new Double_t[nDimensions];
  Double_t* xmax = new Double_t[nDimensions];
  Int_t* nBins = new Int_t[nDimensions];
  for(Int_t idim=0;idim<nDimensions;++idim) {
    nBins[idim] = binLimits[idim].GetSize()-1;
    xmin[idim] = binLimits[idim][0];
    xmax[idim] = binLimits[idim][nBins[idim]];
  }

  THnF* h=new THnF(hname.Data(),arr->At(0)->GetName(),nDimensions,nBins,xmin,xmax);
  for(Int_t idim=0;idim<nDimensions;++idim) {
    TAxis* axis=h->GetAxis(idim);
    axis->Set(nBins[idim], binLimits[idim].GetArray());
  }

  h->Sumw2();
  if(varW!=kNothing) h->SetUniqueID(10+nDimensions+100*(varW+1));
  else h->SetUniqueID(10+nDimensions);
  ULong_t bins = 1;
  for(Int_t idim=0;idim<nDimensions;++idim) {
    bins*=(nBins[idim]+2);
    TAxis* axis = h->GetAxis(idim);
    axis->SetUniqueID(vars[idim]);
    //if(VAR::fUseDefaultVariablesName)
    //  axis->SetTitle(Form("%s %s", VAR::fVariableNames[vars[idim]][0],
    //                (VAR::fVariableNames[vars[idim]][1][0]=='\0' ? "" : Form("(%s)", VAR::fVariableNames[vars[idim]][1]))));
    if(arr->At(1+idim)) axis->SetTitle(arr->At(1+idim)->GetName());
    if(axLabels && !axLabels[idim].IsNull())
      MakeAxisLabels(axis, axLabels[idim].Data());
  }
  hList->Add(h);
  fBinsAllocated+=bins;

  delete [] xmin;
  delete [] xmax;
  delete [] nBins;
  //delete [] binLimits;
}


//____________________________________________________________________________________
void AliQnCorrectionsHistos::MakeAxisLabels(TAxis* ax, const Char_t* labels) {
  //
  // add bin labels to an axis
  //
  TString labelsStr(labels);
  TObjArray* arr=labelsStr.Tokenize(";");
  for(Int_t ib=1; ib<=ax->GetNbins(); ++ib) {
    if(ib>=arr->GetEntries()+1) break;
    ax->SetBinLabel(ib, arr->At(ib-1)->GetName());
  }
}


//__________________________________________________________________
void AliQnCorrectionsHistos::WriteOutput(TFile* save) {
  //
  // Write the histogram lists in the output file
  //
  //cout << "Writing the output to " << save->GetName() << " ... " << flush;
  //TFile* save=new TFile(filename,"RECREATE");
  TDirectory* mainDir = save->mkdir(fHistLists->GetName());
  mainDir->cd();
  for(Int_t i=0; i<fHistLists->GetEntries(); ++i) {
    THashList* list = (THashList*)fHistLists->At(i);
    TDirectory* dir = mainDir->mkdir(list->GetName());
    dir->cd();
    list->Write();
    mainDir->cd();
  }
  save->Close();
  //cout << "done" << endl;
}


