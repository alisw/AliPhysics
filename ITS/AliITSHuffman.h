#ifndef AliITSHUFFMAN_H
#define AliITSHUFFMAN_H

///////////////////////////////////////////////////
//  Huffman Table associated classes for set:ITS //
///////////////////////////////////////////////////


#include <TObject.h>

class AliITSInStream;
class TObjectArray;


//___________________________________________
class AliITSHuffman: public TObject{
  
public:
class AliITSHNode : public TObject {

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
  UChar_t GetSymbol() const {return fSymbol;}
  ULong_t GetFrequency() const {return fFrequency;}
  AliITSHNode *GetLeft() const {return fLeft;}
  AliITSHNode *GetRight() const {return fRight;}
  AliITSHNode *GetFather() const {return fFather;}
  //  void SetSymbol(UChar_r s){fSymbol=s;}
  void SetFrequency(ULong_t fq){fFrequency=fq;}
  void SetLeft(AliITSHNode *n){fLeft = n;}
  void SetRight(AliITSHNode *n){fRight = n;}
  void SetFather(AliITSHNode *n){fFather = n;}


 private:

  UChar_t    fSymbol;        // comment to be written
  ULong_t    fFrequency;     // comment to be written
  AliITSHNode     *fLeft;    // comment to be written
  AliITSHNode     *fRight;   // comment to be written
  AliITSHNode     *fFather;  // not used
};  
  AliITSHuffman(); 
  AliITSHuffman(Int_t size);
  virtual   ~AliITSHuffman();
  AliITSHuffman(const AliITSHuffman &source); // copy constructor
  AliITSHuffman& operator=(const AliITSHuffman &source); // ass. op.
  
  Int_t  Size() const {
    // size
    return fSize;
  }
  UChar_t   *CodeLen() const {
    // code len
    return fCodeLen;
  }
  ULong_t *Code() const {
    // code
    return fCode;
  }
  TObjArray  *HNodes() const {
    // HNodes
    return fHNodes;
  }
  
  
  void GetFrequencies(Int_t len, UChar_t *stream);
  void BuildHTable();   
  Bool_t SpanTree(AliITSHuffman::AliITSHNode*start, ULong_t code, UChar_t len);
  void ResetHNodes();
  void ClearTable();

 protected:

  Int_t          fSize;     // size of the arrays
  UChar_t       *fCodeLen;  //![fSize] number of bits array
  ULong_t       *fCode;     //![fSize] coded symbols array
  
  Short_t       *fSym;      //![fSize] array of input symbols
  TObjArray     *fHNodes;   // array of nodes
  Int_t          fNnodes;   // number of nodes

  ClassDef(AliITSHuffman,1)     //Huffman Table  object for set:ITS
    };

#endif
