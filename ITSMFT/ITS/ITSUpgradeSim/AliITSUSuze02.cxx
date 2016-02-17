#include <Riostream.h>
#include "TMath.h"
#include "AliITSUSuze02.h"
#include "TMatrixF.h"
#include "TH1F.h"

//*******************************************************************
//
//  Simulation of the SUZE02 readout
//  Origin: Serhiy.Senuykov@cern.ch
//  
//*******************************************************************

using std::cout;
using std::endl;

ClassImp(AliITSUSuze02)

AliITSUSuze02::AliITSUSuze02(Int_t Nrows, Int_t Ncols):
  fNRowsChip(Nrows),
  fNColsChip(Ncols),
  fChip(new TMatrixF(Nrows,Ncols)), 
  fTestColumnSize(2),
  fTestRowSize(2),
  fNWindowsPer32colsMax(0),
  fNWindowsPerHalfFSBBMax(0),
  fNWindowsPerFSBBMax(0),
  fNDigitsEncoded(0), 
  fNEncodedWindows(0),
  fNDigitsLost(0),
  fNLostWindows(0),
  fDataSizePerChip(0),
  fNWindowsPer32colsMin(0),
  fNWindowsPerHalfFSBBMin(0),
  fNWindowsPerFSBBMin(0)
{
  if (Ncols%(kNumberOfFSBB*kNumberOfHalfFSBB) != 0) {
    printf("Number of columns should be multiple of %d. SUZE matrix wasn't created\n",kNumberOfFSBB*kNumberOfHalfFSBB);
  }
}  
  
AliITSUSuze02::AliITSUSuze02(const AliITSUSuze02& suze): 
  fNRowsChip(suze.fNRowsChip),
  fNColsChip(suze.fNColsChip),
  fChip(new TMatrixF(*suze.fChip)), 
  fTestColumnSize(suze.fTestColumnSize),
  fTestRowSize(suze.fTestRowSize),
  fNWindowsPer32colsMax(suze.fNWindowsPer32colsMax),
  fNWindowsPerHalfFSBBMax(suze.fNWindowsPerHalfFSBBMax),
  fNWindowsPerFSBBMax(suze.fNWindowsPerFSBBMax),
  fNDigitsEncoded(suze.fNDigitsEncoded), 
  fNEncodedWindows(suze.fNEncodedWindows),
  fNDigitsLost(suze.fNDigitsLost),
  fNLostWindows(suze.fNLostWindows),
  fDataSizePerChip(suze.fDataSizePerChip),
  fNWindowsPer32colsMin(suze.fNWindowsPer32colsMin),
  fNWindowsPerHalfFSBBMin(suze.fNWindowsPerHalfFSBBMin),
  fNWindowsPerFSBBMin(suze.fNWindowsPerFSBBMin)
{
}

AliITSUSuze02 &AliITSUSuze02::operator=(const AliITSUSuze02& suze) {
  if (&suze == this) return *this;

  fNRowsChip = suze.fNRowsChip;
  fNColsChip = suze.fNColsChip;
  fChip = new TMatrixF(*suze.fChip);  
  fTestColumnSize = suze.fTestColumnSize;
  fTestRowSize = suze.fTestRowSize;
  fNWindowsPer32colsMax = suze.fNWindowsPer32colsMax;
  fNWindowsPerHalfFSBBMax = suze.fNWindowsPerHalfFSBBMax;
  fNWindowsPerFSBBMax = suze.fNWindowsPerFSBBMax;
  fNDigitsEncoded = suze.fNDigitsEncoded;
  fNEncodedWindows = suze.fNEncodedWindows;
  fNDigitsLost = suze.fNDigitsLost;
  fNLostWindows = suze.fNLostWindows;
  fDataSizePerChip = suze.fDataSizePerChip;
  fNWindowsPer32colsMin = suze.fNWindowsPer32colsMin;
  fNWindowsPerHalfFSBBMin = suze.fNWindowsPerHalfFSBBMin;
  fNWindowsPerFSBBMin = suze.fNWindowsPerFSBBMin;

  return *this;
}

AliITSUSuze02::~AliITSUSuze02() {
  if(fChip) delete fChip;
}

void AliITSUSuze02::SetEncodingWindowSize(Int_t Wrows, Int_t Wcols){
  fTestColumnSize=Wrows;
  fTestRowSize=Wcols;
}

void AliITSUSuze02::SetQuotas(Int_t q32, Int_t qHalfFSBB, Int_t qFSBB){
  fNWindowsPer32colsMax=q32;
  fNWindowsPerHalfFSBBMax=qHalfFSBB;
  fNWindowsPerFSBBMax=qFSBB;
}

void AliITSUSuze02::AddDigit(Int_t row, Int_t col){
  (*fChip)(row,col)++;
}

//void AliITSUSuze02::Process(Bool_t Verbose){  
void AliITSUSuze02::Process(TH1F* OverflowCodes, TH1F* NDigitsPerEncodingWindowDist, Bool_t Verbose) {

  //cout<<"Processing"<<endl;
  //fChip->Print(); 
  
  Int_t NRowsFSBB=fNRowsChip;
  Int_t NColsFSBB=fNColsChip/kNumberOfFSBB;
  TMatrixF FSBB(NRowsFSBB,NColsFSBB);
  
  Int_t NRowsSuperLine=fTestColumnSize;
  Int_t NColsSuperLine=NColsFSBB;   
  TMatrixF SuperLineUp(NRowsSuperLine,NColsSuperLine);
  TMatrixF SuperLineDown(NRowsSuperLine,NColsSuperLine);
  
  Int_t NRowsSuperLineX2=NRowsSuperLine*2; 
  Int_t NColsSuperLineX2=NColsSuperLine; //SuperLineX2 and SuperLine size in columns is equal to FSBB.
  TMatrixF SuperLineX2(NRowsSuperLineX2,NColsSuperLineX2); 
  
  TMatrixF TestRow(1,fTestRowSize);
  TMatrixF TestColumn(fTestColumnSize,1);      
  Int_t TestRowSum=0;
  TMatrixF EncodingWindow(fTestColumnSize,fTestRowSize);
  
  Int_t EncodingWindowStartRow=0;
  Int_t EncodingWindowStopRow=0;
  Int_t EncodingWindowStartCol=0;
  Int_t EncodingWindowStopCol=0;
  
  Int_t nMasks=fTestRowSize-1;
  Int_t MaskSize=NRowsSuperLineX2-1;
  TMatrixF Masks(MaskSize,nMasks);

   //parameters for internal data size calculation
  Int_t DataSizePerWindowInRow=8+fTestColumnSize*fTestRowSize+TMath::Ceil(TMath::Log2(fTestColumnSize));
  Int_t DataSizePerSuperLineX2Header=30;   
  
  Int_t DataSizePerSuperLineX2=0;
  
  Bool_t Overflow32=kFALSE;
  Bool_t OverflowHalfFSBB=kFALSE;
  Bool_t OverflowFSBB=kFALSE;
  
  Int_t nWindowsPer32cols=fNWindowsPer32colsMax;
  Int_t nWindowsPerHalfFSBB=fNWindowsPerHalfFSBBMax;
  Int_t nWindowsPerFSBB=fNWindowsPerFSBBMax; 
  
  fNWindowsPer32colsMin=fNWindowsPer32colsMax;
  fNWindowsPerHalfFSBBMin=fNWindowsPerHalfFSBBMax;
  fNWindowsPerFSBBMin=fNWindowsPerFSBBMax;  
  
  fNEncodedWindows=0;
  fNDigitsEncoded=0;
  fNLostWindows=0;
  fNDigitsLost=0;
  fDataSizePerChip=0;    
  
  for(Int_t FSBBindex=0; FSBBindex<kNumberOfFSBB; FSBBindex++){
    FSBB=fChip->GetSub(0,NRowsFSBB-1,FSBBindex*NColsFSBB,(FSBBindex+1)*NColsFSBB-1);
    SuperLineDown=FSBB.GetSub(0,NRowsSuperLine-1,0,NColsSuperLine-1);
    for(Int_t SuperLineX2StartRow=0; SuperLineX2StartRow<NRowsFSBB; SuperLineX2StartRow+=NRowsSuperLine){
      if(nWindowsPerFSBB<fNWindowsPerFSBBMin) {fNWindowsPerFSBBMin=nWindowsPerFSBB;} //saving the lowest number of remaining windows
      nWindowsPerFSBB=fNWindowsPerFSBBMax; //Reset number of available encoding windows for the new double SuperLine
      OverflowFSBB=kFALSE;
      SuperLineUp=SuperLineDown;
      if(SuperLineX2StartRow+NRowsSuperLineX2<=NRowsFSBB){
        SuperLineDown=FSBB.GetSub(SuperLineX2StartRow+NRowsSuperLine,SuperLineX2StartRow+NRowsSuperLineX2-1,0,NColsSuperLine-1);
      }
      else if(SuperLineX2StartRow+NRowsSuperLine<NRowsFSBB){
        SuperLineDown.Zero();
        SuperLineDown.SetSub(0,0,FSBB.GetSub(SuperLineX2StartRow+NRowsSuperLine,NRowsFSBB-1,0,NColsSuperLine-1));
      } 
      else{
        SuperLineDown.Zero();
      }
      if(SuperLineUp.Sum()>0){
        DataSizePerSuperLineX2=0; 
        SuperLineX2.SetSub(0,0,SuperLineUp);
        SuperLineX2.SetSub(NRowsSuperLine,0,SuperLineDown);
        for(Int_t HalfFSBBindex=0; HalfFSBBindex<kNumberOfHalfFSBB; HalfFSBBindex++){
          if(nWindowsPerHalfFSBB<fNWindowsPerHalfFSBBMin) {fNWindowsPerHalfFSBBMin=nWindowsPerHalfFSBB;}  
          nWindowsPerHalfFSBB=fNWindowsPerHalfFSBBMax; //reset counter per HalfFSBB
          OverflowHalfFSBB=kFALSE;
          for(Int_t i=0; i<MaskSize; i++){ //reset masks to 1111111
            for(Int_t j=0; j<nMasks; j++){
              Masks(i,j)=1;
            }
          }
          
          for(Int_t TestRowStartCol=HalfFSBBindex*NColsSuperLineX2/kNumberOfHalfFSBB; TestRowStartCol<(HalfFSBBindex+1)*NColsSuperLineX2/kNumberOfHalfFSBB; TestRowStartCol++){
            if(TestRowStartCol%32==0){
	            if(nWindowsPer32cols<fNWindowsPer32colsMin) fNWindowsPer32colsMin=nWindowsPer32cols;
              nWindowsPer32cols=fNWindowsPer32colsMax; //reset nWindowsPer32cols counter every 32 columns
              Overflow32=kFALSE;
            }
            //apply masks
            for(Int_t RowIndex=0; RowIndex<MaskSize; RowIndex++){
              for(Int_t MaskIndex=0; MaskIndex<nMasks; MaskIndex++){
                if(Masks(RowIndex,MaskIndex)==0){
                  // cout<<"Mask has zero bit at pos:"<<RowIndex<<":"<<MaskIndex<<endl;
                  // cout<<"will clean the pixels at row:"<<RowIndex<<" from pos.:"<<TestRowStartCol<<" till "<< TestRowStartCol+(nMasks-MaskIndex)<<endl;
                  for(Int_t ColIndex=TestRowStartCol; ColIndex<TestRowStartCol+(nMasks-MaskIndex); ColIndex++){
                    if(ColIndex<(HalfFSBBindex+1)*NColsSuperLineX2/kNumberOfHalfFSBB){
                      SuperLineX2(RowIndex,ColIndex)=0;
                      // cout<<"Mask has zero bit. Cleaning at pos:"<<RowIndex<<":"<<ColIndex<<endl;
                      if(RowIndex>=fTestColumnSize){
                        SuperLineDown(RowIndex-fTestColumnSize,ColIndex)=0;
                      }
                    }
                  }
                  break;
                }
              }
            }
            //shift masks
            for(Int_t RowIndex=0; RowIndex<MaskSize; RowIndex++){
              for(Int_t MaskIndex=nMasks-1; MaskIndex>0; MaskIndex--){
                Masks(RowIndex,MaskIndex)=Masks(RowIndex,MaskIndex-1);
              }
              Masks(RowIndex,0)=1;
            }  
            
            for(Int_t TestRowStartRow=0; TestRowStartRow<fTestColumnSize; TestRowStartRow++){
              TestRowSum=0;
              for(Int_t TestRowIndex=0; TestRowIndex<fTestRowSize; TestRowIndex++){
                if(TestRowStartCol+TestRowIndex<(HalfFSBBindex+1)*NColsSuperLineX2/kNumberOfHalfFSBB){
                  TestRowSum+=SuperLineX2(TestRowStartRow,TestRowStartCol+TestRowIndex);
                }
                if(TestRowSum>0){
                  //cout<<"TestR at col n."<<TestRowStartCol<<" and row n."<<TestRowStartRow<<" has a hit"<<endl;
                  break;
                }
              }
              if(TestRowSum>0){
                TestColumn=SuperLineX2.GetSub(TestRowStartRow,TestRowStartRow+fTestColumnSize-1,TestRowStartCol,TestRowStartCol);   
                if(TestColumn.Sum()>0){
                  EncodingWindowStartRow=TestRowStartRow;
                  EncodingWindowStopRow=EncodingWindowStartRow+fTestColumnSize-1;
                  EncodingWindowStartCol=TestRowStartCol; 
                  
                  if(TestRowStartCol+fTestRowSize>(HalfFSBBindex+1)*NColsSuperLineX2/kNumberOfHalfFSBB){
                    EncodingWindowStopCol=((HalfFSBBindex+1)*NColsSuperLineX2/kNumberOfHalfFSBB)-1;
                  }
                  else{
                    EncodingWindowStopCol=EncodingWindowStartCol+fTestRowSize-1;
                  }

                  EncodingWindow.Zero();
                  EncodingWindow.SetSub(0,0,SuperLineX2.GetSub(EncodingWindowStartRow,EncodingWindowStopRow,EncodingWindowStartCol,EncodingWindowStopCol));
                  
                  if(nWindowsPer32cols && nWindowsPerHalfFSBB && nWindowsPerFSBB){
                    //cout<<"Will encode window starting at "<<TestRowStartRow<<":"<<TestRowStartCol<<endl;
                    fNDigitsEncoded+=EncodingWindow.Sum();
                    fNEncodedWindows++;
                    OverflowCodes->Fill(0);
		                NDigitsPerEncodingWindowDist->Fill(EncodingWindow.Sum());
                    nWindowsPerFSBB--;
                    nWindowsPerHalfFSBB--;
                    nWindowsPer32cols--;
                    DataSizePerSuperLineX2+=DataSizePerWindowInRow;
                    //cout<<"Windows left:"<<nWindowsPerFSBB<<":"<<nWindowsPerHalfFSBB<<":"<<nWindowsPer32cols<<endl;
                  }
                  else{
                    fNDigitsLost+=EncodingWindow.Sum();
                    fNLostWindows++;
                    //cout<<"------ No encoding at col.:"<<TestRowStartCol<<" at SuperLineX2 that starts at row:"<<SuperLineX2StartRow<<endl;
                    if(!nWindowsPer32cols) Overflow32=kTRUE;
                    if(!nWindowsPerHalfFSBB) OverflowHalfFSBB=kTRUE;
                    if(!nWindowsPerFSBB) OverflowFSBB=kTRUE;
                    OverflowCodes->Fill(Overflow32+2*OverflowHalfFSBB+4*OverflowFSBB);
                  } 
                  for(Int_t k=TestRowStartRow;k<TestRowStartRow+fTestColumnSize; k++){
                    Masks(k,0)=0;
                    //cout<<"Setting Mask to 0 at "<<k<<endl;
                    //Setting to 0 the first column of the encoding window. This part if probably present in the real SUZE.
                    SuperLineX2(k,TestRowStartCol)=0;
                    if(k>=fTestColumnSize){
                      SuperLineDown(k-fTestColumnSize,TestRowStartCol)=0;
                    }
                  }
                }  
                break;
              }
            }
          }
        }
        fDataSizePerChip+=(DataSizePerSuperLineX2+DataSizePerSuperLineX2Header);
      }
    }
  }
  if(Verbose){
    cout<<fNDigitsEncoded<<" digits encoded in "<<fNEncodedWindows<<" windows"<<endl;
    cout<<fNDigitsLost<<" digits lost in "<<fNLostWindows<<" windows"<<endl;
  }
}

void AliITSUSuze02::GetResults(){
  
}
/*
void AliITSUSuze02::InitHistos(){
  fOverflowCodes = new TH1F("OverflowCodes","Overflow codes",8,0,8);
  if(fNRowsChip*fNColsChip){
    fNDigitsPerEncodingWindowDist = new TH1F("nDigitsPerEncodingWindowPerChip","nDigitsPerEncodingWindowPerChip",fTestColumnSize*fTestRowSize,1,fTestColumnSize*fTestRowSize+1);
  }
  else{
    printf("Run AliITSUSuze02::SetEncodingWindowSize first\n");
  }
}
*/
void AliITSUSuze02::ResetChip(){
  fChip->Zero();
}
