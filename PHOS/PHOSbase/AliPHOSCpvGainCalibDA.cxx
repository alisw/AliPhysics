#include "AliPHOSCpvGainCalibDA.h" 
#include "AliPHOSCpvParam.h" 
#include "AliPHOSCpvRawDigiProducer.h"
#include "AliPHOSDigit.h"
#include <fstream>
#include <iostream>
#include <TTree.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TSystem.h>
#include <TTimeStamp.h>
#include "AliPHOSGeometry.h"

#include "TFile.h"

ClassImp(AliPHOSCpvGainCalibDA) ;

using namespace std;

using std::ifstream;
using std::ofstream;
 


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliPHOSCpvGainCalibDA::AliPHOSCpvGainCalibDA():
  TObject(),
  fMinClustSize(4),
  fGeom(0)
{
  //
  //constructor
  //
  
  fGeom=AliPHOSGeometry::GetInstance() ;
  if(!fGeom) fGeom = AliPHOSGeometry::GetInstance("IHEP");

  for(Int_t iDDL = 0; iDDL<2*AliPHOSCpvParam::kNDDL; iDDL++){
    fDeadMap[iDDL] = 0x0;
    fEntriesMap[iDDL] = 0x0;
    for(Int_t iX=0; iX<AliPHOSCpvParam::kPadPcX; iX++) 
      for(Int_t iY=1; iY<AliPHOSCpvParam::kPadPcY; iY++) 
	fAmplA0Histo[iDDL][iX][iY] = 0x0;
  }
  CreateQAHistos();
}  //constructor
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliPHOSCpvGainCalibDA::~AliPHOSCpvGainCalibDA()
{
  //
  //destructor
  //
  for(Int_t iDDL=0; iDDL<2*AliPHOSCpvParam::kNDDL; iDDL++) {
    for(Int_t iX=0; iX<AliPHOSCpvParam::kPadPcX; iX++) {
      for(Int_t iY=1; iY<AliPHOSCpvParam::kPadPcY; iY++) {
	if(fAmplA0Histo[iDDL][iX][iY]) delete fAmplA0Histo[iDDL][iX][iY];
      }//iY
    }//iX
    if(fDeadMap[iDDL])delete fDeadMap[iDDL];
    if(fEntriesMap[iDDL]) delete fEntriesMap[iDDL];
  }//iDDL

  //delete fhErrors;
}  //destructor
//***********************************************************************
void AliPHOSCpvGainCalibDA::InitCalibration(TFile *fCalibrSupplyRoot){
  //tries to open CpvCalibrSupply.root for loading previously filled histograms;
  //creates new histos otherwise.
  //TFile *fCalibrSupplyRoot = TFile::Open("CpvCalibrSupply.root");
  for(Int_t iDDL = 0;iDDL<2*AliPHOSCpvParam::kNDDL;iDDL++){
    if(fCalibrSupplyRoot)
      if(fCalibrSupplyRoot->Get(Form("hEntriesMap%d",iDDL))) 
	fEntriesMap[iDDL]=new TH2I(*(TH2I*)(fCalibrSupplyRoot->Get(Form("hEntriesMap%d",iDDL))));
      else fEntriesMap[iDDL]=0x0;
    else fEntriesMap[iDDL]=0x0;
    for(Int_t iX = 0;iX  <AliPHOSCpvParam::kPadPcX;iX++)
      for(Int_t iY = 0;iY  <AliPHOSCpvParam::kPadPcY;iY++)
	if(fCalibrSupplyRoot){
	  if(fCalibrSupplyRoot->Get(Form("hAmplA0_DDL%d_iX%d_iY%d",iDDL,iX,iY)))
	    fAmplA0Histo[iDDL][iX][iY]=new TH1F(*(TH1F*)(fCalibrSupplyRoot->Get(Form("hAmplA0_DDL%d_iX%d_iY%d",iDDL,iX,iY))));
	  else fAmplA0Histo[iDDL][iX][iY]=0x0;
	}
	else fAmplA0Histo[iDDL][iX][iY]=0x0;
  }
}
//***********************************************************************
void AliPHOSCpvGainCalibDA::CreateA0Histos(Int_t iDDL){
  //create 1D histos for particular DDL to fill them with raw amplitudes later
  if(iDDL<0||iDDL>=2*AliPHOSCpvParam::kNDDL) return; //invalid ddl number
  fEntriesMap[iDDL]=new TH2I(Form("hEntriesMap%d",iDDL),Form("A0 entries map, DDL = %d",iDDL),AliPHOSCpvParam::kPadPcX,0,AliPHOSCpvParam::kPadPcX,AliPHOSCpvParam::kPadPcY,0,AliPHOSCpvParam::kPadPcY);
  fHistosList->Add(fEntriesMap[iDDL]);
  for(Int_t iX = 0;iX  <AliPHOSCpvParam::kPadPcX;iX++)
    for(Int_t iY = 0;iY  <AliPHOSCpvParam::kPadPcY;iY++)
      fAmplA0Histo[iDDL][iX][iY]=new TH1F(Form("hAmplA0_DDL%d_iX%d_iY%d",iDDL,iX,iY),Form("Max amplitude in cluster DDL = %d X = %d Y = %d",iDDL,iX,iY),4096,0.,4096.);
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvGainCalibDA::SetDeadChannelMapFromFile(const char * filename){
  //
  //Set Dead Channel Map Cut from the file, if the input file is not present default value is set!
  //Arguments: the name of the Dead Channel Map file 
  //Returns: false if not possible to open file with provided filename

  TFile *fDeadCh = TFile::Open(filename);
  if(!fDeadCh)return 0;
  for(Int_t iDDL = 0; iDDL < 2*AliPHOSCpvParam::kNDDL; iDDL+=2){
    if(fDeadCh->Get(Form("hBadChMap%d",iDDL)))
      fDeadMap[iDDL] = new TH2I(*(TH2I*)(fDeadCh->Get(Form("hBadChMap%d",iDDL))));
  }
  //fDeadCh->Close();
  return 1;
}//SetDeadChannelMapFromFile()
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliPHOSCpvGainCalibDA::WriteA0HistosToFile(const char* filename) const
{
  if(!filename)return;
  //write all A0 amplitude histos and A0 entries maps to CpvCalibrSupply.root
  TFile * rootF = new TFile(filename,"RECREATE");
  rootF->cd();
  
  for(Int_t iDDL=0; iDDL<2*AliPHOSCpvParam::kNDDL; iDDL+=2){
    if(fEntriesMap[iDDL]) fEntriesMap[iDDL]->Write();
    for(Int_t iX=0; iX<AliPHOSCpvParam::kPadPcX; iX++) 
      for(Int_t iY=0; iY<AliPHOSCpvParam::kPadPcY; iY++) 
	if(fAmplA0Histo[iDDL][iX][iY])fAmplA0Histo[iDDL][iX][iY]->Write();
  }
  rootF->Close();
} //WriteAllHistsToFile()
//*************************************************************
Bool_t AliPHOSCpvGainCalibDA::FillAmplA0Histos(TClonesArray *digits){
  // do a clusterization then find cell with max amlitude (so called A0) in every cluster
  // fill corresponding histos with A0 amplitude
  // return true in case of success (found at least 1 cluster).
  if(!digits) return kFALSE;
  Int_t nDig = digits->GetEntriesFast();
  if(nDig < 1) return kFALSE;//we need at least 1 digit
  Bool_t stop = kFALSE;
  Int_t nExcludedPoints = 0;
  Bool_t *excludedPoints = new Bool_t[nDig];//points which already belongs to other clusters
  for(int i=0;i<nDig;i++)excludedPoints[i]=kFALSE;
  Int_t clusterIndex[100][5][5];//100 clusters max; this array contains digit numbers which belongs to particular cluster
  Int_t clusterDDL[100];// DDL number of particular cluster
  Int_t clusterX[100]; //X coordinate of cluster 
  Int_t clusterY[100]; //Y coordinate of cluster
  Float_t clusterAmplitude[100][5][5];
  //=============================================================================
  //========================= C L U S T E R I S A T I O N =======================
  //=============================================================================
  //here we define cluster as bunch of cells with at least 1 common verticies
// x= |_0_|_1_|_2_|_3_|_4_|
// y=0|   |   |   |   |   |
//    |___|___|___|___|___|
// y=1|   |   |   |   |   |
//    |___|___|___|___|___|
// y=2|   |   |   |   |   |
//    |___|___|___|___|___|
// y=3|   |   |   |   |   |
//    |___|___|___|___|___|
// y=4|   |   |   |   |   |
//    |___|___|___|___|___|
  // initialize clusters array
  for(Int_t iClus=0;iClus<100;iClus++)
    for(Int_t ix=0;ix<5;ix++)
      for(Int_t iy=0;iy<5;iy++)
	clusterIndex[iClus][ix][iy]=-1;
  
  Int_t relId[4];
  Int_t cluNumber = 0;
  while(!stop){//we are going to find 100 or less clusters
    Float_t qMax = 0.;//local maximum value
    Int_t indMax = -1;//local maximum index
    for(Int_t iDig = 0; iDig<nDig;iDig++){//find a local maximum
      if(excludedPoints[iDig]==kTRUE)continue;//is this point already excluded?
      AliPHOSDigit *digit = static_cast<AliPHOSDigit*>( digits->At(iDig) ) ;
      fGeom->AbsToRelNumbering(digit->GetId(),relId);
      if(relId[1]!=-1){//exclude this digit
	nExcludedPoints++; 
	excludedPoints[iDig]=kTRUE;
	continue;//this is not a CPV digit
      }
      int DDL = AliPHOSCpvParam::Mod2DDL(relId[0]);
      if(IsBad(DDL,relId[2]-1,relId[3]-1)){//let's see if it is a bad pad
	nExcludedPoints++; 
	excludedPoints[iDig]=kTRUE;
	continue;
      }
      if( digit->GetEnergy()>qMax) {qMax = digit->GetEnergy(); indMax = iDig;}
    }
    if(indMax<0){//did we find something?
      stop=kTRUE;
      continue;//no new local maximum 
    }
    //set found local maximum as center of new cluster
    nExcludedPoints++; excludedPoints[indMax]=kTRUE; //do not forget to exclude found maximum from consideration
    AliPHOSDigit *digit = static_cast<AliPHOSDigit*>( digits->At(indMax) ) ;
    fGeom->AbsToRelNumbering(digit->GetId(),relId);
    clusterIndex[cluNumber][2][2] = indMax;
    clusterAmplitude[cluNumber][2][2] = qMax;
    clusterDDL[cluNumber] = AliPHOSCpvParam::Mod2DDL(relId[0]);
    clusterX[cluNumber]=relId[2]-1;
    clusterY[cluNumber]=relId[3]-1;
    //let us try to find all other digits which belongs to the same cluster
    Int_t pointsFound = 0;
    do{
      pointsFound=0;
      for(Int_t iDig = 0; iDig<nDig;iDig++){
	if(excludedPoints[iDig]==kTRUE)continue;//is this point already excluded?
	AliPHOSDigit *digit = static_cast<AliPHOSDigit*>( digits->At(iDig) ) ;
	fGeom->AbsToRelNumbering(digit->GetId(),relId);
	if(AliPHOSCpvParam::Mod2DDL(relId[0])!=clusterDDL[cluNumber]) continue;
	//see if this current digit has common vertex with 
	for(int ix = 0;ix<5;ix++)
	  for(int iy = 0;iy<5;iy++)
	    if(clusterIndex[cluNumber][ix][iy]>=0&&
	       (TMath::Abs(relId[2]-1-clusterX[cluNumber]-(ix-2))<=1&& // if X neighbor
		TMath::Abs(relId[3]-1-clusterY[cluNumber]-(iy-2))<=1&& //if Y neighbor
		TMath::Abs(relId[2]-1-clusterX[cluNumber])<=2&& //if digit is within 5x5 region
		TMath::Abs(relId[3]-1-clusterY[cluNumber])<=2)){//we found a cell!
	      pointsFound++;
	      
	      clusterAmplitude[cluNumber][relId[2]-1-clusterX[cluNumber]+2][relId[3]-1-clusterY[cluNumber]+2] = digit->GetEnergy();
	      clusterIndex[cluNumber][relId[2]-1-clusterX[cluNumber]+2][relId[3]-1-clusterY[cluNumber]+2] = iDig;
	      //finally, exclude this point from consideration
	      nExcludedPoints++;
	      excludedPoints[iDig]=kTRUE;
	    }
      }
    }while (pointsFound!=0);
    //OK, we have finished with this cluster
    cluNumber++;
    if(cluNumber>=100) stop=kTRUE; //we found enough clusters
    if(nExcludedPoints>=nDig) stop=kTRUE;//we assigned all the digits
  }
  //cout<<"I found " <<cluNumber<<" clusters"<<endl;

  //now we can operate with clusters
  //=====================================================================================
  //===================== F I L L I N G == O F == H I S T O G R A M S====================
  //=====================================================================================
  
  fhClusterMult->Fill(cluNumber);
  for (Int_t iClu = 0; iClu < cluNumber; iClu++){
    //count cluster size
    Int_t clustSize=0;
    for(int i=0;i<5;i++)
      for(int j=0;j<5;j++)
	if(clusterIndex[iClu][i][j]>0)clustSize++;
    if(clustSize<fMinClustSize)continue;//skip small cluster
    if(!fEntriesMap[clusterDDL[iClu]]) CreateA0Histos(clusterDDL[iClu]);
    // cout<<"iClu = "<<iClu<<
    fAmplA0Histo[clusterDDL[iClu]][clusterX[iClu]][clusterY[iClu]]->Fill(clusterAmplitude[iClu][2][2]);
    fEntriesMap[clusterDDL[iClu]]->Fill(clusterX[iClu],clusterY[iClu]);
    fhA0Value->Fill(clusterAmplitude[iClu][2][2]);
    Double_t totAmpl = 0.;
    for(int ix = 0; ix<5; ix++)
      for(int iy = 0; iy<5; iy++)
	if(clusterIndex[iClu][ix][iy]>=0){
	  fhAmplInClust->Fill(clusterAmplitude[iClu][ix][iy],ix*5+iy);
	  fhClusterShape->Fill(ix-2,iy-2);
	  totAmpl+=clusterAmplitude[iClu][ix][iy];
	}
    fhTotalClusterAmplitude->Fill(totAmpl);
  }
  assert(false);
}
//*************************************************************
void AliPHOSCpvGainCalibDA::CreateQAHistos(){
  fHistosList=new TList();
  
  fhClusterMult = new TH1F("fhClusterMult","Cluster Multiplicity in event",100,0,100);
  fHistosList->Add(fhClusterMult);
  
  fhClusterShape=new TH2F("fhClusterShape","Shape of cluster", 5,-2.5,2.5, 5,-2.5,2.5 );
  fHistosList->Add(fhClusterShape);

  fhA0Value = new TH1F("fhA0Value","Max Amplitude in Cluster ",4096,0.,4096);
  fHistosList->Add(fhA0Value);
  
  fhAmplInClust=new TH2F("fhAmplInClust", "amplitude distribution in cluster", 1000,0.,1000., 25,0.,25.);
  fHistosList->Add(fhAmplInClust);

  fhTotalClusterAmplitude = new TH1F("fhTotalClusterAmplitude", "Total Amplitude in Cluster",4096,0.,4096);
  fHistosList->Add(fhTotalClusterAmplitude);
}
