//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold V0 configuration + results histogram 
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#include <TNamed.h>
#include "AliMachineLearning.h"
#include "AliNeuralNetwork.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TMath.h"
#include <iostream>
#include <TROOT.h>
using namespace std;

ClassImp(AliNeuralNetwork);
//________________________________________________________________
AliNeuralNetwork::AliNeuralNetwork() : AliMachineLearning(),
	fW1(0), fB1(0), fW2(0), fB2(0), fW3(0), fB3(0), fW4(0), fB4(0), fW5(0), fB5(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//________________________________________________________________
AliNeuralNetwork::AliNeuralNetwork(const char * name, const char * title) : AliMachineLearning(name,title),
	fW1(0), fB1(0), fW2(0), fB2(0), fW3(0), fB3(0), fW4(0), fB4(0), fW5(0), fB5(0)

{
    // Named constructor
    
}

//________________________________________________________________
/*
AliNeuralNetwork::AliNeuralNetwork(const AliNeuralNetwork& lCopyMe, TString lNewName) : AliMachineLearning(lCopyMe),
{
    SetName( lNewName.Data() );
}
*/
//________________________________________________________________
AliNeuralNetwork::~AliNeuralNetwork()
{

    // Proper destructor: delete pointer data member
    if(fW1){delete fW1; fW1 = 0x0;}
    if(fB1){delete fB1; fB1 = 0x0;}
    if(fW2){delete fW2; fW2 = 0x0;}
    if(fB2){delete fB2; fB2 = 0x0;}
    if(fW3){delete fW3; fW3 = 0x0;}
    if(fB3){delete fB3; fB3 = 0x0;}
    if(fW4){delete fW4; fW4 = 0x0;}
    if(fB4){delete fB4; fB4 = 0x0;}
    if(fW5){delete fW5; fW5 = 0x0;}
    if(fB5){delete fB5; fB5 = 0x0;}

}

//________________________________________________________________
/*
AliNeuralNetwork& AliNeuralNetwork::operator=(const AliNeuralNetwork& lCopyMe)
{
    if (&lCopyMe == this) return *this;
    //Careful with names
    SetName(lCopyMe.GetName());
    SetTitle(lCopyMe.GetTitle());

    return *this;
}
*/
//________________________________________________________________
void AliNeuralNetwork::LoadModel(TString lModelName)
//Function to load sinaptic weights of neural network and store them in the data members
{
  TFile* FileModel = 0x0;
  FileModel = TFile::Open(lModelName.Data());

  //Check existence, please
  if(!FileModel) cout << Form("Cannot open requested Model file %s", lModelName.Data());
  if(!FileModel->IsOpen()) cout <<Form("Cannot open Model file %s", lModelName.Data());

  TH2D* lW1 = (TH2D*)FileModel->Get("W1");
  if(!lW1) cout << Form("File %s does note contain valid W1 layer", lModelName.Data());
  lW1->SetName("W1_to_copy");
  fW1 = (TH2D*) lW1->Clone("W1");

  TH1D* lB1 = (TH1D*)FileModel->Get("B1");
  if(!lB1) cout << Form("File %s does note contain valid B1 layer", lModelName.Data());
  lB1->SetName("B1_to_copy");
  fB1 = (TH1D*) lB1->Clone("B1");

  TH2D* lW2 = (TH2D*)FileModel->Get("W2");
  if(!lW2) cout << Form("File %s does note contain valid W2 layer", lModelName.Data());
  lW2->SetName("W2_to_copy");
  fW2 = (TH2D*) lW2->Clone("W2");

  TH1D* lB2 = (TH1D*)FileModel->Get("B2");
  if(!lB2) cout << Form("File %s does note contain valid B2 layer", lModelName.Data());
  lB2->SetName("B2_to_copy");
  fB2 = (TH1D*) lB2->Clone("B2");

  TH2D* lW3 = (TH2D*)FileModel->Get("W3");
  if(!lW3) cout << Form("File %s does note contain valid W3 layer", lModelName.Data());
  lW3->SetName("W3_to_copy");
  fW3 = (TH2D*) lW3->Clone("W3");

  TH1D* lB3 = (TH1D*)FileModel->Get("B3");
  if(!lB3) cout << Form("File %s does note contain valid B3 layer", lModelName.Data());
  lB3->SetName("B3_to_copy");
  fB3 = (TH1D*) lB3->Clone("B3");

  TH2D* lW4 = (TH2D*)FileModel->Get("W4");
  if(!lW4) cout << Form("File %s does note contain valid W4 layer", lModelName.Data());
  lW4->SetName("W4_to_copy");
  fW4 = (TH2D*) lW4->Clone("W4");

  TH1D* lB4 = (TH1D*)FileModel->Get("B4");
  if(!lB4) cout << Form("File %s does note contain valid B4 layer", lModelName.Data());
  lB4->SetName("B4_to_copy");
  fB4 = (TH1D*) lB4->Clone("B4");

  TH2D* lW5 = (TH2D*)FileModel->Get("W5");
  if(!lW5) cout << Form("File %s does note contain valid W5 layer", lModelName.Data());
  lW5->SetName("W5_to_copy");
  fW5 = (TH2D*) lW5->Clone("W5");

  TH1D* lB5 = (TH1D*)FileModel->Get("B5");
  if(!lB5) cout << Form("File %s does note contain valid B5 layer", lModelName.Data());
  lB5->SetName("B5_to_copy");
  fB5 = (TH1D*) lB5->Clone("B5");


  fW1->SetDirectory(nullptr);
  fB1->SetDirectory(nullptr);
  fW2->SetDirectory(nullptr);
  fB2->SetDirectory(nullptr);
  fW3->SetDirectory(nullptr);
  fB3->SetDirectory(nullptr);
  fW4->SetDirectory(nullptr);
  fB4->SetDirectory(nullptr);
  fW5->SetDirectory(nullptr);
  fB5->SetDirectory(nullptr);
  FileModel->Close();
}

//________________________________________________________________
double AliNeuralNetwork::Predict(double* X, int K)
//Function to return neural network prediction
{

  if(!fW1||!fB1||!fW2||!fB2||!fW3||!fB3||!fW4||!fB4||!fW5||!fB5) 
  {
 	cout << "Can not Predict: One or more of the NN layers are empty!" << endl;
	return(0);
  }

  //PERFORM FORWARD FUNCTION THROUGH NEURAL NETWORK

  //1st: [1][K]x[K][256] = [1][256]
  double v1[256] = {0};
  for(Int_t j=0; j<256; j++)
  {
  	//Layer product
	for(Int_t k=0; k<K; k++)
		v1[j] += X[k]*fW1->GetBinContent(k+1,j+1);
	//Bias sum
	v1[j]+=fB1->GetBinContent(j+1);
	//Activation function: Sigmoid
	v1[j] = 1.0/(1.0+TMath::Exp(-1.0*v1[j]));
  }

  //2nd: [1][256]x[256][64] = [1][64]
  double v2[64] = {0};
  for(Int_t j=0; j<64; j++)
  {
	//Layer product
	for(Int_t k=0; k<256; k++)
		v2[j] += v1[k]*fW2->GetBinContent(k+1,j+1);
	//Bias sum
	v2[j]+=fB2->GetBinContent(j+1);
	//Activation function: Sigmoid
	v2[j] = 1.0/(1.0+TMath::Exp(-1.0*v2[j]));
  }

  //3rd: [1][64]x[64][16] = [1][16]
  double v3[16] = {0};
  for(Int_t j=0; j<16; j++)
  {
	//Layer product
	for(Int_t k=0; k<64; k++)
		v3[j] += v2[k]*fW3->GetBinContent(k+1,j+1);
	//Bias sum
	v3[j]+=fB3->GetBinContent(j+1);
	//Activation function: Sigmoid
	v3[j] = 1.0/(1.0+TMath::Exp(-1.0*v3[j]));
  }

  //4th: [1][16]x[16][4] = [1][4]
  double v4[4] = {0};
  for(Int_t j=0; j<4; j++)
  {
	//Layer product
	for(Int_t k=0; k<16; k++)
		v4[j] += v3[k]*fW4->GetBinContent(k+1,j+1);
	//Bias sum
	v4[j]+=fB4->GetBinContent(j+1);
	//Activation function: Sigmoid
	v4[j] = 1.0/(1.0+TMath::Exp(-1.0*v4[j]));
  }

  //5th: [1][4]x[4][1] = [1][1]
  double v5 = 0;
  //Layer product
  for(Int_t k=0; k<4; k++)
	v5 += v4[k]*fW5->GetBinContent(k+1,1);
  //Bias sum
  v5+=fB5->GetBinContent(1);
  //Activation function: Sigmoid
  v5 = 1.0/(1.0+TMath::Exp(-1.0*v5));

  return(v5);
}

//________________________________________________________________
/*
void AliNeuralNetwork::Print(Option_t *option="")
{
    cout<<"========================================"<<endl;
    cout<<"    AliNeuralNetwork Configuration      "<<endl;
    cout<<"========================================"<<endl;

    return;
}
*/
