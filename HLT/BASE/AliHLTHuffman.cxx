// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Thorsten Kollegger <kollegge@ikf.uni-frankfurt.de>    *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTHuffman.cxx
/// @author Thorsten Kollegger, Matthias Richter
/// @date   2011-08-14
/// @brief  Huffman code generator/encoder/decoder

#include "AliHLTHuffman.h"

#include <iostream>
#include <set>
#include <bitset>
#include <algorithm>

AliHLTHuffmanNode::AliHLTHuffmanNode() 
	: TObject()
	, fValue(-1)
	, fWeight(0.)
{
        // nop
}

AliHLTHuffmanNode::AliHLTHuffmanNode(const AliHLTHuffmanNode& other)
	: TObject()
	, fValue(other.GetValue())
	, fWeight(other.GetWeight())
{
}

AliHLTHuffmanNode& AliHLTHuffmanNode::operator =(const AliHLTHuffmanNode& other) {
        /// assignment operator
	this->fValue = other.fValue;
	this->fWeight = other.fWeight;
	return *this;
}

AliHLTHuffmanNode::~AliHLTHuffmanNode() {
}

void AliHLTHuffmanNode::AssignCode() {
        /// assign code to this node loop to right and left nodes
	if (GetLeftChild()) {
		// shift current bit pattern to left and add one
		std::bitset < 64 > v = (this->GetBinaryCode() << 1);
		v.set(0);
		GetLeftChild()->SetBinaryCode(this->GetBinaryCodeLength() + 1, v);
		GetLeftChild()->AssignCode();
	}
	if (GetRightChild()) {
		std::bitset < 64 > v = (this->GetBinaryCode() << 1);
		v.reset(0);
		GetRightChild()->SetBinaryCode(this->GetBinaryCodeLength() + 1, v);
		GetRightChild()->AssignCode();
	}
}

void AliHLTHuffmanNode::Print(Option_t* /*option*/) const {
        /// print info
	std::cout << "value=" << GetValue() << ", weight=" << GetWeight() << ", length="
			<< GetBinaryCodeLength() << ", code=" << GetBinaryCode().to_string()
			<< std::endl;
	if (GetLeftChild()) {
		GetLeftChild()->Print();
	}
	if (GetRightChild()) {
		GetRightChild()->Print();
	}
}

ClassImp(AliHLTHuffmanNode)

///////////////////////////////////////////////////////////////////////////////////////////////

AliHLTHuffmanTreeNode::AliHLTHuffmanTreeNode() 
	: AliHLTHuffmanNode()
	, fBinaryCodeLength(0)
	, fBinaryCode(0)
	, fLeft(NULL)
	, fRight(NULL)
{
        // nop
}

AliHLTHuffmanTreeNode::AliHLTHuffmanTreeNode(const AliHLTHuffmanTreeNode& other)
	: AliHLTHuffmanNode()
	, fBinaryCodeLength(other.fBinaryCodeLength)
	, fBinaryCode(other.fBinaryCode)
	, fLeft(other.GetLeftChild())
	, fRight(other.GetRightChild())
{

}

AliHLTHuffmanTreeNode& AliHLTHuffmanTreeNode::operator =(const AliHLTHuffmanTreeNode& other)
{
        /// assignment operator
        if (&other==this) return *this;
	this->fBinaryCodeLength = other.GetBinaryCodeLength();
	this->fBinaryCode = other.GetBinaryCode();
	this->fLeft = other.GetLeftChild();
	this->fRight = other.GetRightChild();
	AliHLTHuffmanNode::operator=(other);
	return *this;
}

AliHLTHuffmanTreeNode::AliHLTHuffmanTreeNode(AliHLTHuffmanNode* l, AliHLTHuffmanNode* r)
	: AliHLTHuffmanNode()
	, fBinaryCodeLength(0)
	, fBinaryCode(0)
	, fLeft(l)
	, fRight(r) {
	if (l && r) {
		SetWeight(l->GetWeight() + r->GetWeight());
	} else if (l && !r) {
		SetWeight(l->GetWeight());
	} else if (!l && r) {
		SetWeight(r->GetWeight());
	}
}

AliHLTHuffmanTreeNode::~AliHLTHuffmanTreeNode() 
{
        // nop
}

ClassImp(AliHLTHuffmanTreeNode)

///////////////////////////////////////////////////////////////////////////////////////////////

AliHLTHuffmanLeaveNode::AliHLTHuffmanLeaveNode() 
	: AliHLTHuffmanNode()
	, fBinaryCodeLength(0)
	, fBinaryCode(0)
	, fLeft(NULL)
	, fRight(NULL)
{
        // nop
}

AliHLTHuffmanLeaveNode::AliHLTHuffmanLeaveNode(const AliHLTHuffmanLeaveNode& other)
	: AliHLTHuffmanNode()
	, fBinaryCodeLength(other.fBinaryCodeLength)
	, fBinaryCode(other.fBinaryCode)
	, fLeft(other.GetLeftChild())
	, fRight(other.GetRightChild())
{

}

AliHLTHuffmanLeaveNode& AliHLTHuffmanLeaveNode::operator =(const AliHLTHuffmanLeaveNode& other)
{
        /// assignment operator
        if (&other==this) return *this;
	this->fBinaryCodeLength = other.GetBinaryCodeLength();
	this->fBinaryCode = other.GetBinaryCode();
	this->fLeft = other.GetLeftChild();
	this->fRight = other.GetRightChild();
	AliHLTHuffmanNode::operator=(other);
	return *this;
}

AliHLTHuffmanLeaveNode::~AliHLTHuffmanLeaveNode() 
{
        // nop
}

ClassImp(AliHLTHuffmanLeaveNode)

///////////////////////////////////////////////////////////////////////////////////////////////

AliHLTHuffman::AliHLTHuffman()
	: TNamed()
	, fMaxBits(0)
	, fMaxValue(0)
	, fNodes(0)
	, fHuffTopNode(NULL)
{
        /// nop
}

AliHLTHuffman::AliHLTHuffman(const AliHLTHuffman& other)
	: TNamed()
	, AliHLTLogging()
	, fMaxBits(other.fMaxBits)
	, fMaxValue(other.fMaxValue)
	, fNodes(other.fNodes)
	, fHuffTopNode(NULL) {
        /// nop
}

AliHLTHuffman::AliHLTHuffman(const char* name, UInt_t maxBits)
        : TNamed(name, name)
	, fMaxBits(maxBits)
	, fMaxValue((((AliHLTUInt64_t) 1) << maxBits) - 1)
	, fNodes((((AliHLTUInt64_t) 1) << maxBits))
	, fHuffTopNode(NULL)
 {
        /// standard constructor
	for (AliHLTUInt64_t i = 0; i <= fMaxValue; i++) {
		fNodes[i].SetValue(i);
	}
}

AliHLTHuffman::~AliHLTHuffman() {
        /// destructor, nop
}

const std::bitset<64>& AliHLTHuffman::Encode(const AliHLTUInt64_t v, AliHLTUInt64_t& codeLength) const {
        /// encode a value
	codeLength = 0;
	if (v <= fMaxValue) {
		// valid symbol/value
		if (fHuffTopNode) {
			// huffman code has been generated
			codeLength = fNodes[v].GetBinaryCodeLength();
			return fNodes[v].GetBinaryCode();
		} else {
		  HLTError("encoder '%s' does not seem to be initialized", GetName());
		}
	} else {
	  HLTError("encoder %s: value %llu exceeds range of %d bits", GetName(), v, GetMaxBits());
	}

	static const std::bitset<64> dummy;
	return dummy;
}

Bool_t AliHLTHuffman::Decode(std::bitset<64> bits, AliHLTUInt64_t& value) const {
	// TODO: check decoding logic, righ now it is just as written
	AliHLTHuffmanNode* currNode = fHuffTopNode;
	if (!currNode) return kFALSE;
	if (currNode->GetValue() >= 0) {
		// handle case with just one node - also quite unlikely
		value = currNode->GetValue();
		return kTRUE;
	}
	while (currNode) {
		if (bits[0] && currNode->GetLeftChild()) {
			// follow left branch
			currNode = currNode->GetLeftChild();
			bits >>= 1;
			if (currNode->GetValue() >= 0) {
				value = currNode->GetValue();
				return kTRUE;
			}
		}
		if (!bits[0] && currNode->GetRightChild()) {
			currNode = currNode->GetRightChild();
			bits >>= 1;
			if (currNode->GetValue() >= 0) {
				value = currNode->GetValue();
				return kTRUE;
			}
		}
	}
	value = ((AliHLTUInt64_t)1) << 63;
	return kFALSE;
}

Bool_t AliHLTHuffman::AddTrainingValue(const AliHLTUInt64_t value,
		const float_t weight) {
	if (value > fMaxValue) {
		/* TODO: ERROR message */
		return kFALSE;
	}
	fNodes[value].AddWeight(weight);
	return kTRUE;
}

Bool_t AliHLTHuffman::GenerateHuffmanTree() {
	// insert pointer to nodes into ordered structure to build tree
	std::multiset<AliHLTHuffmanNode*, AliHLTHuffmanNode::less> nodeCollection;
	//	std::copy(fNodes.begin(), fNodes.end(),
	//			std::inserter(freq_coll, freq_coll.begin()));
	for (std::vector<AliHLTHuffmanLeaveNode>::iterator i = fNodes.begin(); i
			!= fNodes.end(); ++i) {
		nodeCollection.insert(&(*i));
	}
	while (nodeCollection.size() > 1) {
		// insert new node into structure, combining the two with lowest probability
		nodeCollection.insert(
				new AliHLTHuffmanTreeNode(*nodeCollection.begin(),
						*++nodeCollection.begin()));
		nodeCollection.erase(nodeCollection.begin());
		nodeCollection.erase(nodeCollection.begin());
	}
	//assign value
	fHuffTopNode = *nodeCollection.begin();
	fHuffTopNode->AssignCode();
	return kTRUE;
}

void AliHLTHuffman::Print(Option_t* option) const {
        std::cout << GetName() << endl;
        bool bPrintShort=strcmp(option, "short")==0;
	if (fHuffTopNode && !bPrintShort) {
		std::cout << "Huffman tree:" << endl;
		fHuffTopNode->Print();
	}
	Double_t uncompressedSize = 0;
	Double_t compressedSize = 0;
	Double_t totalWeight = 0;
	if (!bPrintShort)
	  std::cout << std::endl << "Huffman codes:" << std::endl;
	for (AliHLTUInt64_t i = 0; i <= fMaxValue; i++) {
	  if (!bPrintShort) fNodes[i].Print();
		totalWeight += fNodes[i].GetWeight();
		uncompressedSize += fNodes[i].GetWeight() * fMaxBits;
		compressedSize += fNodes[i].GetWeight()
				* fNodes[i].GetBinaryCodeLength();
	}
	if (uncompressedSize > 0) {
		std::cout << "compression ratio: " << compressedSize
				/ uncompressedSize << std::endl;
		std::cout << "<bits> uncompressed: " << uncompressedSize / totalWeight
				<< std::endl;
		std::cout << "<bits> compressed:   " << compressedSize / totalWeight
				<< std::endl;
	}
}

AliHLTHuffman& AliHLTHuffman::operator =(const AliHLTHuffman& other) {
	fMaxValue = other.fMaxValue;
	fNodes = other.fNodes;
	fHuffTopNode = NULL;
	return *this;
}

ClassImp(AliHLTHuffman)

