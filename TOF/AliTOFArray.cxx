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
/* $Id$ */

// 
// Class to hold variable size arrays of Float_t
//

#include <TObject.h>
#include <TArrayF.h>
#include <TCollection.h>
#include "AliTOFArray.h"
//#include "AliLog.h"

ClassImp(AliTOFArray)

//-------------------------------------------------------------------
AliTOFArray::AliTOFArray(Int_t size):
	TObject(),
	fSize(size),
	fArray(new TArrayF*[size]){

	// constructor

	for (Int_t i=0;i<size;i++){
		fArray[i] = NULL;
	}
}
//-------------------------------------------------------------------
AliTOFArray::AliTOFArray(const AliTOFArray & source):
	TObject(),fSize(0),fArray(0x0){ 

	//
	// copy constructor
	//

	this->fSize= source.fSize;
	fArray = new TArrayF*[fSize];
	for (Int_t ich = 0; ich<fSize; ich ++){
		fArray[ich] = new TArrayF();
		fArray[ich]->Set(source.fArray[ich]->GetSize());
		for (Int_t j = 0; j < fArray[ich]->GetSize(); j++){
			fArray[ich]->AddAt(fArray[ich]->GetAt(j),j);
		}
	}
}

//-------------------------------------------------------------------
AliTOFArray& AliTOFArray::operator=(const AliTOFArray & source) { 

	//
	// assignment operator
	//

	if (this != &source){
		this->fSize= source.fSize;
		delete [] fArray;
		fArray = new TArrayF*[fSize];
		for (Int_t ich = 0; ich<fSize; ich ++){
			fArray[ich] = new TArrayF();
			fArray[ich]->Set(source.fArray[ich]->GetSize());
			for (Int_t j = 0; j < fArray[ich]->GetSize(); j++){
				fArray[ich]->AddAt(fArray[ich]->GetAt(j),j);
			}
		}
	}
	return *this;
}

//------------------------------------------------------------------
AliTOFArray::~AliTOFArray(){

	//
	// Destructor
	//

	delete [] fArray;
}

//-------------------------------------------------------------------
void AliTOFArray::SetArray(Int_t pos, Int_t size) {

	//
	// adding an array of Float_t with size=size
	//
	
	if (pos>-1 && pos < fSize){
		if (!fArray[pos]) {
			//			printf("Creating array\n");
			fArray[pos] = new TArrayF();			
		}
		fArray[pos]->Reset();
		fArray[pos]->Set(size);
	}
	else printf("Position out of bounds, returning\n");
	return;
}

//-------------------------------------------------------------------
void AliTOFArray::SetAt(Int_t pos, Int_t nelements, Float_t* content) {

	// pos = index of the array to be modified
	// nelements = n. of elements to be modified in array 
	// content = values to be set in array

	if (pos>-1 && pos < fSize){
		if (fArray[pos]){
			Int_t size = fArray[pos]->GetSize();
			if (nelements <= size){
				for (Int_t i=0;i<nelements;i++){
					fArray[pos]->AddAt(content[i],i);
				}
			}
			else printf("Too many elements to be added, returning without adding any\n");
		}
		else printf("Non-existing array, returning\n");
	}
	else printf("Position out of bounds, returning\n");
	return;
}

//-------------------------------------------------------------------
void AliTOFArray::SetAt(Int_t pos, Int_t ielement, Float_t content) {

	// pos = index of the array to be modified
	// ielement = index of the element to be modified in array 
	// content = value to be set in array

	if (pos>-1 && pos < fSize){
		if (fArray[pos]){
			Int_t size = fArray[pos]->GetSize();
			if (ielement < size){
				//printf("Adding %f content in position %d to array %d \n",content, ielement, pos);
				fArray[pos]->AddAt(content,ielement);
			}
			else if (ielement == size) {
				printf ("Increasing the size by 1 and adding a new element to the array\n");
				fArray[pos]->Set(size+1);
				fArray[pos]->AddAt(content,ielement);
			}
			else printf("Not possible to add element %d, size of the array too small, and this would not be the next entry!\n",ielement);
		}
		else printf("Non-existing array, returning\n");
	}	
	else printf("Position out of bounds, returning\n");
	return;
}

//-------------------------------------------------------------------
void AliTOFArray::RemoveArray(Int_t pos) {

	//
	// removing the array at position pos
	//

	if (fArray[pos]) fArray[pos]->Reset();
	else printf("Not possible to remove array, array does not exist\n");
	return;
}

//-------------------------------------------------------------------
Float_t* AliTOFArray::GetArray(Int_t pos) {

	//
	// Getting back array at position pos
	//

	if  (pos>-1 && pos < fSize){
		if (fArray[pos]){
			return fArray[pos]->GetArray();
		}
		else printf("Non-existing array, returning\n");
	}
	else printf("Position out of bounds, returning\n");
	return 0;
}

//-------------------------------------------------------------------
Float_t AliTOFArray::GetArrayAt(Int_t pos, Int_t ielement) {

	//
	// Getting back ielement of array at position pos
	//

	if  (pos>-1 && pos < fSize){
		if (fArray[pos]){
			if (ielement<fArray[pos]->GetSize()){
				return fArray[pos]->GetAt(ielement);
			}
			else printf("Element in array out of bounds, returning\n");
		}
		else printf("Non-existing array, returning\n");
	}
	else printf("Position out of bounds, returning\n");
	return 0;
}

//-------------------------------------------------------------------
void AliTOFArray::ReSetArraySize(Int_t pos, Int_t size) {

	//
	// Changing size of array at position pos, using TArrayF::Set method
	// (without loosing what is already there)
	//

	if  (pos>-1 && pos < fSize){
		if (fArray[pos]){
			fArray[pos]->Set(size);
		}
		else printf("Non-existing array, returning\n");
	}
	else printf("Position out of bounds, returning\n");
	return;
}

//-------------------------------------------------------------------
Int_t AliTOFArray::GetArraySize(Int_t pos) {

	//
	// Getting back size of array at position pos
	//

	if  (pos>-1 && pos < fSize){
		if (fArray[pos]){
			return fArray[pos]->GetSize();
		}
		else printf("Non-existing array, returning\n");
	}
	else printf("Position out of bounds, returning\n");
	return -1;
}

//-------------------------------------------------------------------
Long64_t AliTOFArray::Merge(TCollection *list){

	//
	// Merging method
	//
	
	if (!list) return 0;
	if (list->IsEmpty()) return 1;
	printf("Merging %d AliTOFArrays %s\n", list->GetSize()+1, GetName());
	
	// iterating over the entries in the TList
	TIter next(list);
	AliTOFArray *tofArray;
	Int_t count = 0; // object counter
	while ((tofArray=(AliTOFArray*)next())) {
		//		printf("Count = %d \n",count);
		//if (!tofArray) continue; // dead_code x coverity
		if (tofArray->GetSize() != fSize){
			printf("Merging with current entry in list not possible, AliTOFArray in the list has size different from the current one\n");
			continue;
		}
		for (Int_t i = 0; i<fSize; i++){
			Float_t* tempArray = tofArray->GetArray(i);
			Int_t tempSize = tofArray->GetArraySize(i);
			Int_t currentSize = GetArraySize(i);
			Int_t mergeSize = currentSize+tempSize;
			fArray[i]->Set(mergeSize);
			if (tempSize !=0){
				for (Int_t j = currentSize; j<mergeSize; j++){
					SetAt(i,j,tempArray[j-currentSize]);
				}
			}			
		}
		count++;
		printf("Count = %d \n",count);
		
	}
	return count+1;
	
}

