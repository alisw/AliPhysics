#ifndef __CINT__
#include "AliL3FileHandler.h"
#include "AliL3DigitData.h"
#include "AliL3Transform.h"
#include <stdio.h>
#include <iostream.h>

void QSort(AliL3DigitData **a,Int_t first,Int_t last);
Int_t CompareDigits(AliL3DigitData *a,AliL3DigitData *b);
#endif

void Read2(Int_t triggerevent)
{
  Int_t srow[2] = {0,175};
  UInt_t size;
  const Int_t NEvents = 25;
  //Int_t Offset[25] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240};
  //Int_t Offset[25] = {0,20,40,60,80,100,120,140,160,180,200,220,240,-20,-40,-60,-80,-100,-120,-140,-160,-180,-200,-220,-240};
  Int_t Offset[25]   = {30,60,90,120,150,180,210,240,270,300,330,360,-30,-60,-90,-120,-150,-180,-210,-240,-270,-300,-330,-360,-390};
  Offset[triggerevent]=0;
  AliL3FileHandler *hand[NEvents];
  AliL3DigitRowData *data[NEvents];
  AliL3DigitData *rowData[NEvents];
  AliL3DigitData *test;
  Char_t Carry[255];
  //Int_t event = 0;
  Int_t DigitsTot = 0;
  Int_t MeVDigit;
  Int_t DigitsPerRow = 0;
  //AliL3DigitRowData *NData;
  Byte_t *NData;
  Int_t tot_dig=0;
  
  AliL3Transform *transform = new AliL3Transform();
  const Int_t maxdigits=50000;
  AliL3DigitData *temp[maxdigits];
  AliL3DigitData *temp2[maxdigits];
  
  
  for(Int_t slice=0 ; slice<36 ; slice++){
    for(Int_t event = 0; event < NEvents ; event++){
      hand[event] = new AliL3FileHandler();
      sprintf(Carry,"/prog/alice/data/Rawdata/1_patch/pp/event_%d/digits_%d_0.raw",event,slice);
      if(!hand[event]->SetBinaryInput(Carry))
	cerr<<" Error opening file :"<<Carry<<endl;
      cout << "Open File: " << Carry;
      hand[event]->Init(slice,0,srow);
      data[event] = hand[event]->CompBinary2Memory(size);
      cout<<" done"<<endl;
      hand[event]->CloseBinaryInput();
    }
    DigitsTot=0;
    Int_t sign = slice < 18 ? -1 : 1;
    for(Int_t row = 0 ; row < 176 ; row++){
      
      for(Int_t event = 0; event < NEvents ; event++){
	rowData[event] = (AliL3DigitData *)data[event]->fDigitData;
	for(UInt_t dig=0; dig<data[event]->fNDigit; dig++)
	  {
	    Int_t time = rowData[event][dig].fTime + sign*Offset[event];
	    if(time < transform->GetNTimeBins() && time >= 0)
	      DigitsTot++;
	  }
	//DigitsTot += data[event]->fNDigit;
	hand[event]->UpdateRowPointer(data[event]);
      }
      
    }
    //cout << "Try to allocate : " << DigitsTot << " Digits in row : " << row << endl;
    
    
    //NData = (AliL3DigitRowData*) malloc(AllDigitsSize);// Allocate !
    
    Int_t AllDigitsSize = sizeof(AliL3DigitData) * DigitsTot + sizeof(AliL3DigitRowData) * 176;
    NData = new Byte_t[AllDigitsSize];
    memset(NData,0,AllDigitsSize);
    AliL3DigitRowData *AllRowData = (AliL3DigitRowData*)NData;
    
    //cout << "Allocated " << endl;
    
    //Reset the data pointers, because they changed when doing UpdateRowPointer.
    for(Int_t event=0; event<NEvents; event++)
      data[event] = (AliL3DigitRowData*)hand[event]->GetDataPointer(size);
    
    tot_dig=0;
    for(Int_t row = 0 ; row < 176 ; row++){
      DigitsPerRow = 0;
      for(Int_t event = 0; event < NEvents ; event++){
	rowData[event] = (AliL3DigitData *)data[event]->fDigitData;
	for(UInt_t dig=0; dig<data[event]->fNDigit; dig++)
	  {
	    Int_t time = rowData[event][dig].fTime + sign*Offset[event];
	    if(time < transform->GetNTimeBins() && time >= 0)
	      DigitsPerRow++;
	  }
	//DigitsPerRow += data[event]->fNDigit;// + DigitsPerRow; 
      }
      //tot_dig += DigitsPerRow;
      //cout << "Found : " << DigitsPerRow << " in all Events, Row :" << row << endl;
      
      /*
	AllRowData->fRow = row;
	AllRowData->fNDigit = DigitsPerRow;
	test = (AliL3DigitData *) AllRowData->fDigitData;
      */
      
      //cout << "Got Row Pointer " << endl;
      MeVDigit = 0;
      for(Int_t event = 0; event < NEvents ; event++){
	rowData[event] = (AliL3DigitData *)data[event]->fDigitData;
	for(UInt_t digit = 0; digit < data[event]->fNDigit ; digit++){
	  
	  //Copy the digits to a temporary storage, because we have to sort them before storing them :
	  Int_t time = rowData[event][digit].fTime + sign*Offset[event];
	  if(time >= transform->GetNTimeBins() || time < 0)
	    continue;
	  
	  temp[MeVDigit] = &rowData[event][digit];
	  temp[MeVDigit]->fTime = temp[MeVDigit]->fTime + sign*Offset[event];
	  
	  //cout << MeVDigit << " : " << digit << endl;
	  MeVDigit++;
	}
	hand[event]->UpdateRowPointer(data[event]);
      }
      //Sort the digits:
      QSort(temp,0,MeVDigit);
      int final_count=0;
      
      //Now, loop and check for overlaps:
      for(Int_t c=0; c<MeVDigit; c++)
	{
	  
	  if(c>0 && temp[c]->fTime == temp[c-1]->fTime && temp[c]->fPad == temp[c-1]->fPad)
	    {//this one is overlapping with the previous one:
	      temp2[final_count-1]->fCharge += temp[c+1]->fCharge;
	      if(temp2[final_count-1]->fCharge>=1024)
		temp2[final_count-1]->fCharge=1023;//Saturation
	    }
	  else
	    {
	      temp2[final_count] = temp[c];
	      final_count++;
	    }
	}
      tot_dig+=final_count;
      AllRowData->fRow = row;
      AllRowData->fNDigit = final_count;
      test = (AliL3DigitData *) AllRowData->fDigitData;
      
      for(Int_t c=0; c<final_count; c++)//and copy them to the new event:
	{
	  test[c].fCharge = temp2[c]->fCharge;
	  test[c].fPad = temp2[c]->fPad;
	  test[c].fTime = temp2[c]->fTime;
	  if(c>0 && test[c].fTime == test[c-1].fTime && test[c].fPad == test[c-1].fPad)
	    cout<<"Overlap at row "<<row<<" pad "<<(int)test[c].fPad<<" time "<<(int)test[c].fTime<<endl;
	  //  cout<<"Copying back, pad "<<(int)test[c].fPad<<" time "<<(int)test[c].fTime<<" charge "<<(int)test[c].fCharge<<endl;
	}
      
      if(MeVDigit!=DigitsPerRow)
	cerr<<endl<<"Error: "<<MeVDigit<<" "<<DigitsPerRow<<endl;
      if(final_count > MeVDigit)
	cerr<<"Error; final_count "<<final_count<<" MeVDigit "<<MeVDigit<<endl;
      Byte_t *tmp = (Byte_t *)AllRowData;
      Int_t UpdateSize = sizeof(AliL3DigitRowData) + sizeof(AliL3DigitData)*final_count;//DigitsPerRow;
      tmp += UpdateSize;
      AllRowData = (AliL3DigitRowData *) tmp; 
      
    }
    if(tot_dig>DigitsTot)
      cerr<<endl<<"Mismatching digitcount "<<tot_dig<<" "<<DigitsTot<<endl;
    
    for(Int_t event = 0; event < NEvents ; event++)
      {
	data[event]=0;
	delete hand[event];
      }
    
    AliL3FileHandler *Out = new AliL3FileHandler();
    //sprintf(Carry,"/prog/alice/data/Rawdata/PileUp/digits_1_0.raw");
    AllRowData = (AliL3DigitRowData*)NData;
    sprintf(Carry,"/prog/alice/data/Rawdata/1_patch/pp/pileups/event_%d/digits_%d_0.raw",triggerevent,slice);
    Out->SetBinaryOutput(Carry);
    cout << "Open File: " << Carry << " to wite PileUp Event" << endl;
    Out->Memory2CompBinary(176,AllRowData);
    Out->CloseBinaryOutput();
    delete Out;
    delete [] NData;
  }
  
  delete transform;
}

void QSort(AliL3DigitData **a,Int_t first,Int_t last)
{
  //General sorting routine. only sorting the pointers.

  AliL3DigitData *tmp;
  int i; 
  int j;
  
  while (last - first > 1) {
    i = first;
    j = last;
    for (;;) {
      while (++i < last && CompareDigits(a[i], a[first]) < 0)
	;
      while (--j > first && CompareDigits(a[j], a[first]) > 0)
            ;
         if (i >= j)
            break;

         tmp  = a[i];
         a[i] = a[j];
         a[j] = tmp;
    }
      if (j == first) {
         ++first;
         continue;
      }
      tmp = a[first];
      a[first] = a[j];
      a[j] = tmp;
      if (j - first < last - (j + 1)) {
         QSort(a, first, j);
         first = j + 1;   // QSort(j + 1, last);
      } else {
         QSort(a, j + 1, last);
         last = j;        // QSort(first, j);
      }
   }
  
  
}

Int_t CompareDigits(AliL3DigitData *a,AliL3DigitData *b)
{

  if(a->fPad==b->fPad && a->fTime == b->fTime) return 0;
  
  if(a->fPad<b->fPad) return -1;
  if(a->fPad==b->fPad && a->fTime<b->fTime) return -1;
  
  return 1;
  
}
