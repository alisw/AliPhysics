////////////////////////////////////////////////
//  RawData classes for set:ITS               //
////////////////////////////////////////////////

#include <TMath.h>

#include <iostream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "AliITSHuffman.h"
#include "AliITSRawData.h"

ClassImp(AliITSHNode)

//_____________________________________________________________________________

AliITSHNode::AliITSHNode()
{
  // constructor
    fLeft=0;
    fRight=0;
    fFather=0;
}
//_____________________________________________________________________________

AliITSHNode::AliITSHNode(UChar_t sym, ULong_t freq)
{
  // standard constructor
    fSymbol=sym;
    fFrequency=freq;
    fLeft=0;
    fRight=0;
    fFather=0;
}

//__________________________________________________________________________
AliITSHNode::AliITSHNode(const AliITSHNode &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fSymbol = source.fSymbol;
  this->fFrequency = source.fFrequency;
  this->fLeft = source.fLeft;
  this->fRight = source.fRight;
  this->fFather = source.fFather;
  return;
}

//_________________________________________________________________________
AliITSHNode& 
  AliITSHNode::operator=(const AliITSHNode &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fSymbol = source.fSymbol;
  this->fFrequency = source.fFrequency;
  this->fLeft = source.fLeft;
  this->fRight = source.fRight;
  this->fFather = source.fFather;
  return *this;
}

//____________________________________________
Int_t AliITSHNode::Compare(TObject *obj)
{
  // function called by Sort method of TObjArray

         AliITSHNode *node=(AliITSHNode *)obj;
	 ULong_t f=fFrequency;
         ULong_t fo=node->fFrequency;
         if (f>fo) return 1;
         else if (f<fo) return -1;
         else return 0;
}
//_____________________________________________________________________________


ClassImp(AliITSHTable)

//_____________________________________________________________________________

AliITSHTable::AliITSHTable()
{
  // constructor
    fCodeLen=0;
    fCode=0;
    fHNodes=0;
    fNnodes=0;
   
}
//_____________________________________________________________________________

AliITSHTable::AliITSHTable(Int_t size)
{
  //
  // Creates the look-up table for the 1D compression
  //

  //initialise

  fSize=size;
  fCodeLen = new UChar_t[fSize]; 
  fCode = new ULong_t[fSize]; 
  fHNodes = new TObjArray;
  fNnodes=0;
  fSym= new Short_t[fSize];
  for (Short_t i=0;i<fSize;i++) {
       fSym[i]=i;
  }
  ClearTable(); 

}

//__________________________________________________________________________
AliITSHTable::AliITSHTable(const AliITSHTable &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fSize = source.fSize;
  this->fCodeLen = source.fCodeLen;
  this->fCode = source.fCode;
  this->fSym = source.fSym;
  this->fHNodes = source.fHNodes;
  this->fNnodes = source.fNnodes;
  return;
}

//_________________________________________________________________________
AliITSHTable& 
  AliITSHTable::operator=(const AliITSHTable &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fSize = source.fSize;
  this->fCodeLen = source.fCodeLen;
  this->fCode = source.fCode;
  this->fSym = source.fSym;
  this->fHNodes = source.fHNodes;
  this->fNnodes = source.fNnodes;
  return *this;
}

//_____________________________________________________________________________
void AliITSHTable::GetFrequencies(Int_t len, UChar_t *stream)
{
  // get frequencies
  printf("Get Frequencies: sym %p \n",fSym);

  // use temporarily the fCode array to store the frequencies
  for (Int_t i=0; i< len; i++) {
      Int_t idx=TMath::BinarySearch(fSize,fSym,(Short_t)stream[i]);
      if (idx == (Int_t)stream[i]) fCode[idx]++;
      // test
      if(idx==134) cout<< "idx fCode[idx]  "<<idx<<" "<<fCode[idx]<<endl;
      //printf("idx,fCode[i] %d %d\n",idx,(Int_t)fCode[idx]);
  }


}


//_____________________________________________________________________________
void AliITSHTable::BuildHTable()
{
  // build Htable

  for (Int_t i=0; i< fSize; i++) {
    //printf("i,fCode[i] %d %d\n",i,(Int_t)fCode[i]);
     if (fCode[i] > 0) {
        fNnodes++;
        cout<< "i fCode[i] fNnodes "<<i<<" "<<fCode[i]<<" "<<fNnodes<<endl;
	//printf("i, fCode[i] fNnodes %d %d %d\n",i,fCode[i],fNnodes);
        fHNodes->Add(new AliITSHNode((UChar_t)i,fCode[i]));
     }
  }

  Int_t nentries=fHNodes->GetEntriesFast();
  Int_t nindex=nentries-1;  
  printf("nentries fNnodes nindex %d %d %d\n",nentries,fNnodes,nindex);

  while (nindex > 0) 
    {

     fHNodes->Sort(nindex);
     AliITSHNode *aux = new AliITSHNode(0,0);
     AliITSHNode *node= (AliITSHNode*)fHNodes->UncheckedAt(nindex-1);
     AliITSHNode *node1= (AliITSHNode*)fHNodes->UncheckedAt(nindex);
     aux->fLeft = node;
     aux->fRight = node1;
     aux->fFrequency = node->fFrequency + node1->fFrequency;
     printf("symbol symbol1 freq freq1 %d %d %d %d\n",(int)node->fSymbol,(int)node1->fSymbol,(int)node->fFrequency,(int)node1->fFrequency);
     cout << "aux - frequency "<< (Int_t)(aux->fFrequency) <<endl;
     fHNodes->RemoveAt(nindex-1);
     fHNodes->AddAt(aux,nindex-1);
     nindex--;
     printf("nindex, obj at nindex %d %p \n",nindex,(AliITSHNode*)fHNodes->UncheckedAt(nindex));

    }

    ClearTable();

    AliITSHNode *start= (AliITSHNode*)fHNodes->UncheckedAt(0);
    SpanTree(start,0,0);
    
    // check the Huffman table

    cout << "...Done, Huffman Table is: \n";
    for (int c=0; c <= 255; c++) {
      if (fCodeLen[c] > 0) cout << "Symbol " << c << " Coded as " << fCode[c] << " and long " << (int) fCodeLen[c] << " bits.\n"; 
    }

}

//_____________________________________________________________________________
AliITSHTable::~AliITSHTable()
{
  // HTable
    printf("HTable destructor !\n");
    if (fCodeLen) delete[] fCodeLen;
    if (fCode) delete [] fCode;
    delete fHNodes;
}


//____________________________________________
Bool_t AliITSHTable::SpanTree(AliITSHNode *start, ULong_t code, UChar_t len)
{
  // span tree
  AliITSHNode * visited;
  visited = start;

  printf("outside: code, len %d %d\n",(int)code,(int)len);

  Int_t idx=(Int_t)visited->fSymbol;
  if (!visited->fLeft) {
	fCode[idx] = code; 
	fCodeLen[idx] = len;
        printf("idx, fCode[idx], fCodeLen[idx] %d %d %d\n",idx,(int)fCode[idx],
              (int)fCodeLen[idx]);
	return kTRUE;
  }

// reccursive stuff

  if (SpanTree(visited->fLeft, code << 1, len + 1)) {
          printf("code, len %d %d\n",(int)code,(int)len);
	  if (visited->fRight) 
                  SpanTree(visited->fRight, code << 1 | 0x01, len + 1);
  }
  return kTRUE;
}

//____________________________________________
void AliITSHTable::ResetHNodes()
{
    //
    // Reset number of HNodes and the HNodes array 
    //
	if (fHNodes)  fHNodes->Clear();
	if (fNnodes)  fNnodes=0;

}

//_____________________________________________________________________________
void AliITSHTable::ClearTable()
{
  // clear
    memset(fCodeLen,0,sizeof(UChar_t)*fSize);
    memset(fCode,0,sizeof(ULong_t)*fSize);
}

