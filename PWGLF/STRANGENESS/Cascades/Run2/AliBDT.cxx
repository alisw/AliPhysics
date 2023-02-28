#include <TNamed.h>
#include "AliMachineLearning.h"
#include "AliBDT.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include <iostream>
#include <TROOT.h>
using namespace std;

ClassImp(AliBDT);
//________________________________________________________________
AliBDT::AliBDT() : AliMachineLearning(),
	fFt(0), fSoS(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//________________________________________________________________
AliBDT::AliBDT(const char * name, const char * title) : AliMachineLearning(name,title),
	fFt(0), fSoS(0)

{
    // Named constructor
    
}

//________________________________________________________________
/*
AliBDT::AliBDT(const AliBDT& lCopyMe, TString lNewName) : AliMachineLearning(lCopyMe),
{
    SetName( lNewName.Data() );
}
*/
//________________________________________________________________
AliBDT::~AliBDT()
{

    // Proper destructor: delete pointer data member
    if(fFt){delete fFt; fFt = 0x0;}
    if(fSoS){delete fSoS; fSoS = 0x0;}
}

//________________________________________________________________
/*
AliBDT& AliBDT::operator=(const AliBDT& lCopyMe)
{
    if (&lCopyMe == this) return *this;
    //Careful with names
    SetName(lCopyMe.GetName());
    SetTitle(lCopyMe.GetTitle());

    return *this;
}
*/
//________________________________________________________________
void AliBDT::LoadModel(TString lModelName)
//Function to load sinaptic weights of neural network and store them in the data members
{
  TFile* FileModel = 0x0;
  FileModel = TFile::Open(lModelName.Data(),"READ");

  //Check existence, please
  if(!FileModel) cout << Form("Cannot open requested Model file %s", lModelName.Data());
  if(!FileModel->IsOpen()) cout <<Form("Cannot open Model file %s", lModelName.Data());

  TH1I* lFt = (TH1I*)FileModel->Get("FeatureIndex");
  if(!lFt) cout << Form("File %s does note contain valid data member", lModelName.Data());
  lFt->SetName("FeatureIndex_to_copy");
  fFt = (TH1I*)lFt->Clone("FeatureIndex");

  TH1D* lSoS = (TH1D*)FileModel->Get("Split_or_Score");
  if(!lSoS) cout << Form("File %s does note contain valid data member", lModelName.Data());
  lSoS->SetName("Split_or_Score_to_copy");
  fSoS = (TH1D*)lSoS->Clone("Split_or_Score");
  
  fFt->SetDirectory(nullptr);
  fSoS->SetDirectory(nullptr);
  FileModel->Close();
}

//________________________________________________________________
double AliBDT::Predict(double* X, int Nt)
//Function to return BDT prediction
//Nt = number of trees
{

  if(!fFt||!fSoS) 
  {
 	cout << "Can not Predict: One or more of the data members are empty!" << endl;
	return(0);
  }

  //PERFORM FORWARD FUNCTION THROUGH BDT

  double P = 0.0;

  int nodes = int(fFt->GetNbinsX()/Nt); //nodes per tree = 2^(tree_depth + 1) - 1

  for(int i=0; i<Nt; i++)
  {
    //booster i structure
    int Lidx = 1;			//local index of a node in each tree [1,nodes]
    int Gidx = (nodes*i)+Lidx;	//global index of a node in Nt trees

    while( Lidx < (nodes+1)/2 )
    {
      if( X[int(fFt->GetBinContent(Gidx))] < fSoS->GetBinContent(Gidx) ) Lidx = 2*Lidx;

      else 								 Lidx = 2*Lidx+1;

      Gidx = (nodes*i)+Lidx;
    }

    //Sum leaf score
    P += fSoS->GetBinContent(Gidx);
  }
  
/*
  //100 boosters
  for(int i=0; i<100; i++)
  {
    //Booster variables
    int idx0 = int( fFt->GetBinContent((7*i)+1) );
    double Split0 = fSoS->GetBinContent((7*i)+1);

    int idx1 = int( fFt->GetBinContent((7*i)+1+1) );
    double Split1 = fSoS->GetBinContent((7*i)+1+1);

    int idx2 = int( fFt->GetBinContent((7*i)+2+1) );
    double Split2 = fSoS->GetBinContent((7*i)+2+1);

    //3 to 6 should be leaves
    int idx3 = int( fFt->GetBinContent((7*i)+3+1) );
    double Score3 = fSoS->GetBinContent((7*i)+3+1);

    int idx4 = int( fFt->GetBinContent((7*i)+4+1) );
    double Score4 = fSoS->GetBinContent((7*i)+4+1);

    int idx5 = int( fFt->GetBinContent((7*i)+5+1) );
    double Score5 = fSoS->GetBinContent((7*i)+5+1);

    int idx6 = int( fFt->GetBinContent((7*i)+6+1) );
    double Score6 = fSoS->GetBinContent((7*i)+6+1);

    if(idx3+idx4+idx5+idx6!=-4)
      cout << "somethin is wrong with a leaf!" << endl;

    //booster i structure
    if( X[idx0]<Split0 )
    {
      //0:yes
      if(X[idx1]<Split1)
        //1:yes
        P+=Score3;
      else
        //1:no
        P+=Score4;
    }
    else
    {
      //0:no
      if(X[idx2]<Split2)
        //2:yes
        P+=Score5;
      else
        //2:no
        P+=Score6;
    }
  }
*/
  P = 1.0/(1.0+TMath::Exp(-1.0*P));

  return(P);
}

//________________________________________________________________
/*
void AliBDT::Print(Option_t *option="")
{
    cout<<"========================================"<<endl;
    cout<<"    AliBDT Configuration      "<<endl;
    cout<<"========================================"<<endl;

    return;
}
*/
