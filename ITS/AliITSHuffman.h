#ifndef AliITSHUFFMAN_H
#define AliITSHUFFMAN_H

///////////////////////////////////////////////////
//  Huffman Table associated classes for set:ITS //
///////////////////////////////////////////////////

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Attention! Two classes in this file.
// They have to stay in the same file.
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include <TObject.h>

class AliITSInStream;
class TObjectArray;
class AliITSHNode: public TObject  {

 public:
  AliITSHNode();
  AliITSHNode(UChar_t symbol, ULong_t freq);
  virtual   ~AliITSHNode() {
    // destructor
  }
  AliITSHNode(const AliITSHNode &source); // copy constructor
  AliITSHNode& operator=(const AliITSHNode &source); // ass. op.

  Bool_t IsSortable() const {
    // is sortable
    return kTRUE;
  }
  Int_t Compare(const TObject *obj) const;
  
  ClassDef(AliITSHNode,1)     //HuffT node object for set:ITS

 public:

  UChar_t    fSymbol;        // comment to be written
  ULong_t    fFrequency;     // comment to be written
  AliITSHNode     *fLeft;    // comment to be written
  AliITSHNode     *fRight;   // comment to be written
  AliITSHNode     *fFather;  // not used
};

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  Attention! Next class has kept deliberaty in 
//  the same file as the previous one
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
//___________________________________________
class AliITSHTable: public TObject{
  
public:
  AliITSHTable(); 
  AliITSHTable(Int_t size);
  virtual   ~AliITSHTable();
  AliITSHTable(const AliITSHTable &source); // copy constructor
  AliITSHTable& operator=(const AliITSHTable &source); // ass. op.
  
  Int_t  Size() {
    // size
    return fSize;
  }
  UChar_t   *CodeLen() {
    // code len
    return fCodeLen;
  }
  ULong_t *Code() {
    // code
    return fCode;
  }
  TObjArray  *HNodes() {
    // HNodes
    return fHNodes;
  }
  
  
  void GetFrequencies(Int_t len, UChar_t *stream);
  void BuildHTable();   
  Bool_t SpanTree(AliITSHNode*start, ULong_t code, UChar_t len);
  void ResetHNodes();
  void ClearTable();
  
 protected:

  Int_t          fSize;     // size of the arrays
  UChar_t       *fCodeLen;  //![fSize] number of bits array
  ULong_t       *fCode;     //![fSize] coded symbols array
  
  Short_t       *fSym;      //![fSize] array of input symbols
  TObjArray     *fHNodes;   // array of nodes
  Int_t          fNnodes;   // number of nodes

  ClassDef(AliITSHTable,1)     //Huffman Table  object for set:ITS
    };

#endif
