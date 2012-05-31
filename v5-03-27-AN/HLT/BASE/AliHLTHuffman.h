//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTHUFFMAN_H
#define ALIHLTHUFFMAN_H
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTDataHuffman.h
/// @author Thorsten Kollegger, Matthias Richter
/// @date   2011-08-14
/// @brief  Huffman code generator/encoder/decode

#include "AliHLTLogging.h"

#include "TNamed.h"

#include <vector>
#include <bitset>
#include <string>

/**
 * @class AliHLTHuffmanNode
 * Huffman code nodes. This is the base class for two types of nodes,
 * AliHLTHuffmanTreeNode and AliHLTHuffmanLeaveNode. The latter nodes
 * represent the index array of values to huffman codes. The former are
 * used in a tree reaching from the value of highest occuremce to all
 * leaves. The two node types are defined in order to optimize the storage
 * format for the huffman table. Leave nodes don't need persistent pointers
 * to childs, tree nodes don't need persistent binary code values.
 *
 * @ingroup alihlt_base
 */
class AliHLTHuffmanNode: public TObject {
public:
	/// default constructor
	AliHLTHuffmanNode();
	/// copy constructor
	AliHLTHuffmanNode(const AliHLTHuffmanNode& other);
	/// assignment operator
	AliHLTHuffmanNode& operator=(const AliHLTHuffmanNode& other);
	/// destructor
	virtual ~AliHLTHuffmanNode();

	/// set symbol value of node
	void SetValue(AliHLTInt64_t v) {
		fValue = v;
	}
	/// add weight to node
	void AddWeight(Float_t w) {
		fWeight += w;
	}
	/// set weight of node
	void SetWeight(Float_t w) {
		fWeight = w;
	}
	/// return symbol value of node
	AliHLTInt64_t GetValue() const {
		return fValue;
	}
	/// return weight of node
	Float_t GetWeight() const {
		return fWeight;
	}

        /// assign huffman code to this node and its children
        /// bReverse = true the bit corresponding to the parent is shifted and
        /// decoding needs to start from the MSB of the code
	void AssignCode(bool bReverse=false);
	/// set binary huffman code and length of node
	virtual void SetBinaryCode(AliHLTUInt64_t length, std::bitset<64> code) = 0;
	/// set pointer to left child
	virtual void SetLeftChild(AliHLTHuffmanNode* n) = 0;
	/// set pointer to right child
	virtual void SetRightChild(AliHLTHuffmanNode* n) = 0;
	/// Return length of huffman code
	virtual AliHLTUInt64_t GetBinaryCodeLength() const = 0;
	/// Return binary huffman code
	virtual const std::bitset<64>& GetBinaryCode() const = 0;
	/// Return pointer to left Child
	virtual AliHLTHuffmanNode* GetLeftChild() const = 0;
	/// Return pointer to right Child
	virtual AliHLTHuffmanNode* GetRightChild() const = 0;

	/// Print information about node
	void Print(Option_t* option = "") const;
	/// Overload less operator, based on node weights
	Bool_t operator <(const AliHLTHuffmanNode& other) const {
		return (this->GetWeight() < other.GetWeight());
	}
	class less {
	public:
		bool operator()(const AliHLTHuffmanNode* i1,
				const AliHLTHuffmanNode* i2) const {
			//reverse sort, less likely to most likely
			return ((*i1) < (*i2));
		}
	};

private:
	AliHLTInt64_t fValue;        // value
	Float_t fWeight;             // weight

ClassDef(AliHLTHuffmanNode, 1)
};

/**
 * @class AliHLTHuffmanTreeNode
 * Tree nodes store the childs persistently. The binary code is
 * a transient member.
 */
class AliHLTHuffmanTreeNode: public AliHLTHuffmanNode {
public:
	/// default constructor
	AliHLTHuffmanTreeNode();
	/// copy constructor
	AliHLTHuffmanTreeNode(const AliHLTHuffmanTreeNode& other);
	/// assignment operator
	AliHLTHuffmanTreeNode& operator=(const AliHLTHuffmanTreeNode& other);
	/// constructor for internal nodes, based on two input nodes (leaves or internal nodes)
	AliHLTHuffmanTreeNode(AliHLTHuffmanNode* l, AliHLTHuffmanNode* r);
	/// desstructor
	~AliHLTHuffmanTreeNode();

	/// set binary huffman code and length of node
	void SetBinaryCode(AliHLTUInt64_t length, std::bitset<64> code) {
	        fBinaryCodeLength = length;
	        fBinaryCode = code;
	}

	/// set pointer to left child
	void SetLeftChild(AliHLTHuffmanNode* n) {
		fLeft = n;
	}
	/// set pointer to right child
	void SetRightChild(AliHLTHuffmanNode* n) {
		fRight = n;
	}
	/// Return length of huffman code
	AliHLTUInt64_t GetBinaryCodeLength() const {
		return fBinaryCodeLength;
	}
	/// Return binary huffman code
	const std::bitset<64>& GetBinaryCode() const {
		return fBinaryCode;
	}
	/// Return pointer to left Child
	AliHLTHuffmanNode* GetLeftChild() const {
		return fLeft;
	}
	/// Return pointer to right Child
	AliHLTHuffmanNode* GetRightChild() const {
		return fRight;
	}

private:
	AliHLTUInt8_t fBinaryCodeLength; //! code length
	std::bitset<64> fBinaryCode; //! WARNING: this fixed the maximum code length to 128
	AliHLTHuffmanNode* fLeft;    // left neighbor
	AliHLTHuffmanNode* fRight;   // right neighbor

ClassDef(AliHLTHuffmanTreeNode, 1)
};

/**
 * @class AliHLTHuffmanLeaveNode
 * Leave nodes store the binary code persistently, while the childs are transient
 */
class AliHLTHuffmanLeaveNode: public AliHLTHuffmanNode {
public:
	/// default constructor
	AliHLTHuffmanLeaveNode();
	/// copy constructor
	AliHLTHuffmanLeaveNode(const AliHLTHuffmanLeaveNode& other);
	/// assignment operator
	AliHLTHuffmanLeaveNode& operator=(const AliHLTHuffmanLeaveNode& other);
	/// destructor
	~AliHLTHuffmanLeaveNode();

	/// set binary huffman code and length of node
	void SetBinaryCode(AliHLTUInt64_t length, std::bitset<64> code) {
	        fBinaryCodeLength = length;
	        fBinaryCode = code;
	}
	/// set pointer to left child
	void SetLeftChild(AliHLTHuffmanNode* n) {
		fLeft = n;
	}
	/// set pointer to right child
	void SetRightChild(AliHLTHuffmanNode* n) {
		fRight = n;
	}
	/// Return length of huffman code
	AliHLTUInt64_t GetBinaryCodeLength() const {
		return fBinaryCodeLength;
	}
	/// Return binary huffman code
	const std::bitset<64>& GetBinaryCode() const {
		return fBinaryCode;
	}
	/// Return pointer to left Child
	AliHLTHuffmanNode* GetLeftChild() const {
		return fLeft;
	}
	/// Return pointer to right Child
	AliHLTHuffmanNode* GetRightChild() const {
		return fRight;
	}

private:
	AliHLTUInt8_t fBinaryCodeLength; // code length
	std::bitset<64> fBinaryCode; // WARNING: this fixed the maximum code length to 128
	AliHLTHuffmanNode* fLeft;    //! left neighbor
	AliHLTHuffmanNode* fRight;   //! right neighbor

ClassDef(AliHLTHuffmanLeaveNode, 1)
};

/**
 * @class AliHLTHuffman
 * Huffman code generator/encoder/decoder
 *
 * @ingroup alihlt_base
 */
class AliHLTHuffman: public TNamed, public AliHLTLogging {
public:
	AliHLTHuffman();
	AliHLTHuffman(const AliHLTHuffman& other);
	AliHLTHuffman(const char* name, UInt_t maxBits);
	~AliHLTHuffman();

	UInt_t GetMaxBits() const {return fMaxBits;}
	UInt_t GetMaxValue() const {return fMaxValue;}

	/// Return huffman code for a value
	const std::bitset<64>& Encode(const AliHLTUInt64_t v, AliHLTUInt64_t& codeLength) const;

	/// Return value for bit pattern
        Bool_t Decode(std::bitset<64> /*bits*/, AliHLTUInt64_t& /*value*/, AliHLTUInt32_t& length, AliHLTUInt32_t& codeLength) const;
	/// Add a new training value (with optional weight) to the training sample
	Bool_t AddTrainingValue(const AliHLTUInt64_t value,
			const Float_t weight = 1.);
	/// Generate huffman tree from training sample
	Bool_t GenerateHuffmanTree();
	/// Print info about huffman en-/decoder
	void Print(Option_t* option = "short") const;
	/// Overload assignment operator
	AliHLTHuffman& operator =(const AliHLTHuffman& other);

        bool CheckConsistency() const;

private:
	UInt_t fMaxBits;    // bit lenght
	UInt_t fMaxValue;   // maximum value
	std::vector<AliHLTHuffmanLeaveNode> fNodes; // array of nodes
	AliHLTHuffmanNode* fHuffTopNode;       // top node
	bool fReverseCode; // indicate the type of the binary code

ClassDef(AliHLTHuffman, 2)
};

#endif
