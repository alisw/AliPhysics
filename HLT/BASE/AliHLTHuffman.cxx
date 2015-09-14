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
#include <iomanip>
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
        if (this==&other) return *this;
	this->fValue = other.fValue;
	this->fWeight = other.fWeight;
	return *this;
}

AliHLTHuffmanNode::~AliHLTHuffmanNode() {
}

void AliHLTHuffmanNode::AssignCode(bool bReverse) {
        /// assign code to this node loop to right and left nodes
        /// the decoding always has to start from the least significant bit since the
        /// code length is variable. Thats why the bit corresponding to the parent node
        /// has to be right of the bit of child nodes, i.e. bits correspond to the
        /// current code length. For storage in a bit stream however, bits are stored
        /// in a stream from MSB to LSB and overwrapping to the MSBs of the next byte.
        /// Here the reverse code is needed and the word of fixed length read from the
        /// stream needs to be reversed before decoding.
        /// Note: by changing the AliHLTDataDeflater interface to write from LSB to MSB
        /// this can be avoided.
	if (GetLeftChild()) {
	  if (bReverse) {
		std::bitset < 64 > v = (this->GetBinaryCode() << 1);
		v.set(0);
		GetLeftChild()->SetBinaryCode(this->GetBinaryCodeLength() + 1, v);
	  } else {
		std::bitset < 64 > v = (this->GetBinaryCode());
		int codelen=this->GetBinaryCodeLength();
		v.set(codelen);
		GetLeftChild()->SetBinaryCode(codelen + 1, v);
	  }
	  GetLeftChild()->AssignCode(bReverse);
	}
	if (GetRightChild()) {
	  if (bReverse) {
		std::bitset < 64 > v = (this->GetBinaryCode() << 1);
		v.reset(0);
		GetRightChild()->SetBinaryCode(this->GetBinaryCodeLength() + 1, v);
	  } else {
		std::bitset < 64 > v = (this->GetBinaryCode());
		int codelen=this->GetBinaryCodeLength();
		v.reset(codelen);
		GetRightChild()->SetBinaryCode(codelen + 1, v);
	  }
	  GetRightChild()->AssignCode(bReverse);
	}
}

void AliHLTHuffmanNode::Print(Option_t* /*option*/) const {
        /// print info
        std::streamsize precbackup = std::cout.precision();
        std::streamsize widthbackup = std::cout.width();
        std::ios::fmtflags flagsbackup = std::cout.flags();
        std::cout.precision(6);
        std::cout << " value="  << std::setw(6) << std::left << GetValue()
                  << " weight=" << std::setw(12) << std::scientific << GetWeight()
                  << " length=" << std::setw(3) << GetBinaryCodeLength()
                  << " code=" << GetBinaryCode().to_string()
                  << std::endl;
	if (GetLeftChild()) {
		GetLeftChild()->Print();
	}
	if (GetRightChild()) {
		GetRightChild()->Print();
	}
        std::cout.precision(precbackup);
        std::cout.width(widthbackup);
        std::cout.flags(flagsbackup);
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
	: AliHLTHuffmanNode(other)
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
	: AliHLTHuffmanNode(other)
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
	, fReverseCode(true)
	, fMaxCodeLength(0)
	, fDecodingNodes()
	, fDecodingTopNode(NULL)
{
        /// nop
}

AliHLTHuffman::AliHLTHuffman(const AliHLTHuffman& other)
	: TNamed()
	, AliHLTLogging()
	, fMaxBits(other.fMaxBits)
	, fMaxValue(other.fMaxValue)
	, fNodes(other.fNodes)
	, fHuffTopNode(NULL)
	, fReverseCode(other.fReverseCode)
	, fMaxCodeLength(other.fMaxCodeLength)
	, fDecodingNodes()
	, fDecodingTopNode(NULL)
{
        /// nop
}

AliHLTHuffman::AliHLTHuffman(const char* name, UInt_t maxBits)
        : TNamed(name, name)
	, fMaxBits(maxBits)
	, fMaxValue((((AliHLTUInt64_t) 1) << maxBits) - 1)
	, fNodes((((AliHLTUInt64_t) 1) << maxBits))
	, fHuffTopNode(NULL)
	, fReverseCode(true)
	, fMaxCodeLength(0)
	, fDecodingNodes()
	, fDecodingTopNode(NULL)
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

Bool_t AliHLTHuffman::DecodeLSB(std::bitset<64> bits, AliHLTUInt64_t& value,
				AliHLTUInt32_t& length, AliHLTUInt32_t& codeLength) const
{
	// huffman decoding
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
			if (currNode && currNode->GetValue() >= 0) {
				value = currNode->GetValue();
				length = fMaxBits;
				codeLength = currNode->GetBinaryCodeLength();
				return kTRUE;
			}
			continue;
		}
		if (!bits[0] && currNode->GetRightChild()) {
			currNode = currNode->GetRightChild();
			bits >>= 1;
			if (currNode && currNode->GetValue() >= 0) {
				value = currNode->GetValue();
				length = fMaxBits;
				codeLength = currNode->GetBinaryCodeLength();
				return kTRUE;
			}
			continue;
		}
		break;
	}
	value = ((AliHLTUInt64_t)1) << 63;
	return kFALSE;
}

Bool_t AliHLTHuffman::DecodeMSB(std::bitset<64> bits, AliHLTUInt64_t& value,
			        AliHLTUInt32_t& length, AliHLTUInt32_t& codeLength) const
{
	// huffman decoding
	AliHLTHuffmanNode* currNode = fHuffTopNode;
	if (!currNode) return kFALSE;
	if (currNode->GetValue() >= 0) {
		// handle case with just one node - also quite unlikely
		value = currNode->GetValue();
		return kTRUE;
	}
	while (currNode) {
		if (bits[63] && currNode->GetLeftChild()) {
			// follow left branch
			currNode = currNode->GetLeftChild();
			bits <<= 1;
			if (currNode && currNode->GetValue() >= 0) {
				value = currNode->GetValue();
				length = fMaxBits;
				codeLength = currNode->GetBinaryCodeLength();
				return kTRUE;
			}
			continue;
		}
		if (!bits[63] && currNode->GetRightChild()) {
			currNode = currNode->GetRightChild();
			bits <<= 1;
			if (currNode && currNode->GetValue() >= 0) {
				value = currNode->GetValue();
				length = fMaxBits;
				codeLength = currNode->GetBinaryCodeLength();
				return kTRUE;
			}
			continue;
		}
		break;
	}
	value = ((AliHLTUInt64_t)1) << 63;
	return kFALSE;
}

Bool_t AliHLTHuffman::FastDecodeMSB(std::bitset<64> bits, AliHLTUInt64_t& value,
				    AliHLTUInt32_t& length, AliHLTUInt32_t& codeLength) const
{
	// huffman decoding using the binary node map
        HLTFatal("this function needs to be tested and commisioned");
	AliHuffmanDecodingNode* currNode = fDecodingTopNode;
	if (!currNode) return kFALSE;
	if (currNode->fValue >= 0) {
		// handle case with just one node - also quite unlikely
		value = currNode->fValue;
		return kTRUE;
	}
	while (currNode) {
		if (bits[63] && currNode->fLeft) {
			// follow left branch
			currNode = currNode->fLeft;
			bits <<= 1;
			if (currNode && currNode->fValue >= 0) {
				value = currNode->fValue;
				length = fMaxBits;
				codeLength = currNode->fBinaryCodeLength;
				return kTRUE;
			}
			continue;
		}
		if (!bits[63] && currNode->fRight) {
			currNode = currNode->fRight;
			bits <<= 1;
			if (currNode && currNode->fValue >= 0) {
				value = currNode->fValue;
				length = fMaxBits;
				codeLength = currNode->fBinaryCodeLength;
				return kTRUE;
			}
			continue;
		}
		break;
	}
	value = ((AliHLTUInt64_t)1) << 63;
	return kFALSE;
}

Bool_t AliHLTHuffman::AddTrainingValue(const AliHLTUInt64_t value,
		const Float_t weight) {
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
	Double_t totalWeight = 0.;
	Double_t leastWeight = -1.;
	for (std::vector<AliHLTHuffmanLeaveNode>::iterator i = fNodes.begin(); i
			!= fNodes.end(); ++i) {
		nodeCollection.insert(&(*i));
		Double_t nodeWeight = i->GetWeight();
		totalWeight += nodeWeight;
		if (nodeWeight > 0. &&
		    (leastWeight < 0. || leastWeight > nodeWeight))
		  leastWeight = nodeWeight;
	}
	Double_t scaleexp = 0.;
	if (leastWeight > 0.) {
	  scaleexp = std::floor(std::log10(leastWeight));
	}
	// float has 23 bit fraction, to be accurate for the occurence counting
	// number has to fit into that fraction (exponent 0)
	if (totalWeight > (0x1 << 23)) {
	  if (scaleexp > 0.) {
	    HLTInfo("scaling weights by exponent %f, total weight %f", scaleexp, totalWeight);
	    for (std::vector<AliHLTHuffmanLeaveNode>::iterator i = fNodes.begin();
		 i != fNodes.end(); ++i) {
	      i->ScaleWeight(std::pow(10., scaleexp));
	    }
	  } else {
	    // we can not scale, through a warning
	    HLTWarning("single float precision does not sufficiently represent "
		       "small differences at total weight > 2^23 (~ 8.3e6), "
		       "total weight %f; huffman table will be generated correctly "
		       "however the weigth of nodes stored in the object has a "
		       "small inaccuracy",
		       totalWeight
		       );
	  }
	}
	while (nodeCollection.size() > 1) {
		// insert new node into structure, combining the two with lowest probability
		AliHLTHuffmanNode* node=new AliHLTHuffmanTreeNode(*nodeCollection.begin(), *++nodeCollection.begin());
		if (!node) return kFALSE;
		nodeCollection.insert(node);
		nodeCollection.erase(nodeCollection.begin());
		nodeCollection.erase(nodeCollection.begin());
	}
	//assign value
	fHuffTopNode = *nodeCollection.begin();
	fHuffTopNode->AssignCode(fReverseCode);
	InitMaxCodeLength();
	return kTRUE;
}

void AliHLTHuffman::Print(Option_t* option) const {
        std::cout << GetName() << endl;
        bool bPrintShort=strcmp(option, "full")!=0;
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
	  Double_t nodeWeight = fNodes[i].GetWeight();
	  if (!bPrintShort) fNodes[i].Print();
	  totalWeight += nodeWeight;
	  uncompressedSize += nodeWeight * fMaxBits;
	  compressedSize += nodeWeight * fNodes[i].GetBinaryCodeLength();
	}
	if (uncompressedSize > 0) {
		std::cout << "compression ratio: " << compressedSize
				/ uncompressedSize << std::endl;
		std::cout << "<bits> uncompressed: " << uncompressedSize / totalWeight
				<< std::endl;
		std::cout << "<bits> compressed:   " << compressedSize / totalWeight
				<< std::endl;
	} else {
	  std::cout << "   no weights available" << std::endl;
	}
}

AliHLTHuffman& AliHLTHuffman::operator =(const AliHLTHuffman& other) {
        if (this==&other) return *this;
	fMaxValue = other.fMaxValue;
	fNodes = other.fNodes;
	fHuffTopNode = NULL;
	fMaxCodeLength = 0;
	return *this;
}

bool AliHLTHuffman::CheckConsistency() const
{
  if (!fHuffTopNode) {
    cout << "huffman table not yet generated" << endl;
  }

  for (AliHLTUInt64_t v=0; v<GetMaxValue(); v++) {
    AliHLTUInt64_t codeLength=0;
    std::bitset<64> code=AliHLTHuffman::Encode(v, codeLength);
    AliHLTUInt64_t readback=0;
    AliHLTUInt32_t readbacklen=0;
    AliHLTUInt32_t readbackcodelen=0;
    if (fReverseCode) {
      code<<=64-codeLength;
      if (!DecodeMSB(code, readback, readbacklen, readbackcodelen)) {
	cout << "Decode failed" << endl;
	return false;
      }
    } else {
    if (!DecodeLSB(code, readback, readbacklen, readbackcodelen)) {
      cout << "Decode failed" << endl;
      return false;
    }
    }
    if (v!=readback) {
      cout << "readback of value " << v << " code length " << codeLength << " failed: got " << readback << " code length " << readbackcodelen << endl;
      return false;
    }
  }
  return true;
}

UInt_t AliHLTHuffman::InitMaxCodeLength()
{
  // loop over leave nodes and set maximum code length
  fMaxCodeLength=0;
  for (std::vector<AliHLTHuffmanLeaveNode>::const_iterator node=fNodes.begin();
       node!=fNodes.end(); node++) {
    if (fMaxCodeLength<node->GetBinaryCodeLength())
      fMaxCodeLength=node->GetBinaryCodeLength();
  }
  return fMaxCodeLength;
}

int AliHLTHuffman::EnableDecodingMap()
{
  // build decoder nodes from node tree
  fDecodingTopNode=NULL;
  fDecodingNodes.clear();
  if (fNodes.size()==0) {
    return 0;
  }
  fDecodingNodes.reserve(2*fNodes.size());
  fDecodingTopNode=BuildDecodingNode(fHuffTopNode, fDecodingNodes);
  if (!fDecodingTopNode) {
    fDecodingNodes.clear();
  } else {
    HLTError("decoder nodes created sucessfully");
  }
  return 0;
}

AliHLTHuffman::AliHuffmanDecodingNode* AliHLTHuffman::BuildDecodingNode(AliHLTHuffmanNode* node, vector<AliHuffmanDecodingNode>& decodernodes) const
{
  // build and add decoder node structure corresponding to huffman node
  if (!node) {
    HLTError("invalid node pointer");
    return NULL;
  }
  for (vector<AliHuffmanDecodingNode>::iterator i=decodernodes.begin();
       i!=decodernodes.end(); i++) {
    if (i->fParent==node) {
      HLTError("node %p already existing in decoder nodes", node);
      return NULL;
    }
  }

  AliHuffmanDecodingNode* decodernode=NULL;
  if (decodernodes.size()+1>decodernodes.capacity()) {
    HLTError("initial size of array too small, can not add more nodes because pointers will become invalid when growing vector");
  } else {
    AliHuffmanDecodingNode dn;
    dn.fParent=node;
    dn.fValue=node->GetValue();
    dn.fLeft=NULL;
    dn.fRight=NULL;
    dn.fBinaryCodeLength=0;
    decodernodes.push_back(dn);
    decodernode=&(decodernodes.back());
  }

  if (decodernode && decodernode->fValue<0) {
    decodernode->fRight=BuildDecodingNode(node->GetRightChild(), decodernodes);
    decodernode->fLeft=BuildDecodingNode(node->GetLeftChild(), decodernodes);
    if (decodernode->fLeft==NULL || decodernode->fRight==0) {
      HLTError("failed to build corresponding decoder node for node %p", node);
      decodernode=NULL;
    }
  }

  return decodernode;
}

ClassImp(AliHLTHuffman)

