/**************************************************************************
 * Copyright(c) 2006-2008, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

////////////////////////////////////////////////
//                                            //
//  RawData classes for set:ITS               //
//                                            //
////////////////////////////////////////////////

#include <TMath.h>
#include <TObjArray.h>
#include <Riostream.h>

#include "AliITSHuffman.h"

ClassImp(AliITSHuffman)

//_____________________________________________________________________________

  AliITSHuffman::AliITSHNode::AliITSHNode(): 
TObject(),
fSymbol(),
fFrequency(0),
fLeft(),
fRight(),
fFather() {
  // default constructor
}
//_____________________________________________________________________________

AliITSHuffman::AliITSHNode::AliITSHNode(UChar_t sym, ULong_t freq):
TObject(),
fSymbol(sym),
fFrequency(freq),
fLeft(),
fRight(),
fFather() {
  // standard constructor
}

//__________________________________________________________________________
AliITSHuffman::AliITSHNode::AliITSHNode(const AliITSHNode &source): 
TObject(source),
fSymbol(source.fSymbol),
fFrequency(source.fFrequency),
fLeft(source.fLeft),
fRight(source.fRight),
fFather(source.fFather) {
  //     Copy Constructor 
  return;
}

//_________________________________________________________________________
AliITSHuffman::AliITSHNode& 
  AliITSHuffman::AliITSHNode::operator=(const AliITSHuffman::AliITSHNode &source) {
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
Int_t AliITSHuffman::AliITSHNode::Compare(const TObject *obj) const
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

AliITSHuffman::AliITSHuffman():
TObject(),
fSize(0),
fCodeLen(),
fCode(),
fSym(),
fHNodes(),
fNnodes(0)
{
  // default constructor
   
}
//_____________________________________________________________________________

AliITSHuffman::AliITSHuffman(Int_t size):
TObject(),
fSize(size),
fCodeLen(),
fCode(),
fSym(),
fHNodes(),
fNnodes(0)
{
  //
  // Creates the look-up table for the 1D compression
  //

  //initialise

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
AliITSHuffman::AliITSHuffman(const AliITSHuffman &source) : 
TObject(source),
fSize(source.fSize),
fCodeLen(source.fCodeLen),
fCode(source.fCode),
fSym(source.fSym),
fHNodes(source.fHNodes),
fNnodes(source.fNnodes)
{
  //     Copy Constructor 
}

//_________________________________________________________________________
AliITSHuffman& 
  AliITSHuffman::operator=(const AliITSHuffman &source) {
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
void AliITSHuffman::GetFrequencies(Int_t len, UChar_t *stream)
{
  // get frequencies
  printf("Get Frequencies: sym %p \n",(void*)fSym);

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
void AliITSHuffman::BuildHTable()
{
  // build Htable

  for (Int_t i=0; i< fSize; i++) {
    //printf("i,fCode[i] %d %d\n",i,(Int_t)fCode[i]);
     if (fCode[i] > 0) {
        fNnodes++;
        cout<< "i fCode[i] fNnodes "<<i<<" "<<fCode[i]<<" "<<fNnodes<<endl;
	//printf("i, fCode[i] fNnodes %d %d %d\n",i,fCode[i],fNnodes);
        fHNodes->Add(new AliITSHuffman::AliITSHNode((UChar_t)i,fCode[i]));
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
     aux->SetLeft(node);
     aux->SetRight(node1);
     aux->SetFrequency(node->GetFrequency() + node1->GetFrequency());
     printf("symbol symbol1 freq freq1 %d %d %d %d\n",(int)node->GetSymbol(),(int)node1->GetSymbol(),(int)node->GetFrequency(),(int)node1->GetFrequency());
     cout << "aux - frequency "<< (Int_t)(aux->GetFrequency()) <<endl;
     fHNodes->RemoveAt(nindex-1);
     fHNodes->AddAt(aux,nindex-1);
     nindex--;
     printf("nindex, obj at nindex %d %p \n",nindex,(void*)fHNodes->UncheckedAt(nindex));

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
AliITSHuffman::~AliITSHuffman()
{
  // HTable
    printf("HTable destructor !\n");
    if (fCodeLen) delete[] fCodeLen;
    if (fCode) delete [] fCode;
    if (fHNodes) {
      fHNodes->Delete();
      delete fHNodes;
    }
}


//____________________________________________
Bool_t AliITSHuffman::SpanTree(AliITSHNode *start, ULong_t code, UChar_t len)
{
  // span tree
  AliITSHNode * visited;
  visited = start;

  printf("outside: code, len %d %d\n",(int)code,(int)len);

  Int_t idx=(Int_t)visited->GetSymbol();
  if (!visited->GetLeft()) {
	fCode[idx] = code; 
	fCodeLen[idx] = len;
        printf("idx, fCode[idx], fCodeLen[idx] %d %d %d\n",idx,(int)fCode[idx],
              (int)fCodeLen[idx]);
	return kTRUE;
  }

// reccursive stuff

  if (SpanTree(visited->GetLeft(), code << 1, len + 1)) {
          printf("code, len %d %d\n",(int)code,(int)len);
	  if (visited->GetRight()) 
                  SpanTree(visited->GetRight(), code << 1 | 0x01, len + 1);
  }
  return kTRUE;
}

//____________________________________________
void AliITSHuffman::ResetHNodes()
{
    //
    // Reset number of HNodes and the HNodes array 
    //
	if (fHNodes)  fHNodes->Clear();
	if (fNnodes)  fNnodes=0;

}

//_____________________________________________________________________________
void AliITSHuffman::ClearTable()
{
  // clear
    memset(fCodeLen,0,sizeof(UChar_t)*fSize);
    memset(fCode,0,sizeof(ULong_t)*fSize);
}

