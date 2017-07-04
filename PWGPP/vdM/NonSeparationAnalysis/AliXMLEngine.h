// -*- C++ -*-
#ifndef _ALI_XML_ENGINE_H
#define _ALI_XML_ENGINE_H
#include <algorithm>

#include <TXMLEngine.h>
#include <TString.h>

#include "AliLog.h"

// wrapper for TXMLEngine
class AliXMLEngine : private TXMLEngine {
 public:
  AliXMLEngine(const char* fileName)
    : TXMLEngine()
    , fXMLDoc(ParseFile(fileName)) {
    if (!fXMLDoc)
      AliFatalF("TXMLEngine::ParseFile(\"%s\") returned nullptr", fileName);
  }
  virtual ~AliXMLEngine() {
    if (fXMLDoc)
      FreeDoc(fXMLDoc);
    fXMLDoc = nullptr;
  }

  // policy class for adapting XMLNodePointer_t to an iterator
  struct NodePolicy {
    typedef XMLNodePointer_t ptr_type;
    const char* GetName(AliXMLEngine *xml, ptr_type nodePtr) const { return xml->GetNodeName(nodePtr); }
    const char* GetData(AliXMLEngine *xml, ptr_type nodePtr) const { return xml->GetNodeContent(nodePtr); }
    ptr_type    GetNext(AliXMLEngine *xml, ptr_type nodePtr) const { return xml->GetNext(nodePtr); }
  } ;
  // policy class for adapting XMLAttrPointer_t to an iterator
  struct AttrPolicy {
    typedef XMLAttrPointer_t ptr_type;
    const char* GetName(AliXMLEngine *xml, ptr_type attrPtr) const { return xml->GetAttrName(attrPtr); }
    const char* GetData(AliXMLEngine *xml, ptr_type attrPtr) const { return xml->GetAttrValue(attrPtr); }
    ptr_type    GetNext(AliXMLEngine *xml, ptr_type attrPtr) const { return xml->GetNextAttr(attrPtr); }
  } ;
  // generic iterator (for now const) to be used for nodes and attributes
  template<typename Policy>
  class IteratorBase
    : private Policy
    , private std::iterator<std::forward_iterator_tag,  // iterator_category
                            IteratorBase<Policy>,       // value_type
                            std::ptrdiff_t,             // difference_type
                            typename Policy::ptr_type,  // pointer
                            const IteratorBase<Policy>& // reference
                            > {
  public:
    typedef IteratorBase<Policy> base_type;
    typedef typename base_type::reference reference;
    virtual ~IteratorBase() {}

    const char* GetName() const { return Policy::GetName(fXML, fPtr); }
    const char* GetData() const { return Policy::GetData(fXML, fPtr); }

    base_type& operator++() {
      fPtr = Policy::GetNext(fXML, fPtr);
      return *this;
    }
    bool      operator!=(reference other) const { return fPtr != other.fPtr; }
    reference operator*() const { return *this; }
  protected:
    typedef typename base_type::pointer pointer;
    IteratorBase(AliXMLEngine *xml, pointer ptr)
      : fXML(xml)
      , fPtr(ptr) {}
    pointer   operator&() const { return fPtr; }
  private:
    AliXMLEngine *fXML; //!
    pointer       fPtr; //!
  } ;

  // iterator for attributes: typedef is enough
  typedef IteratorBase<AttrPolicy> Attr;

  // iterator for nodes: we need a derived class in order to have Get{Child*,Attr*} methods
  class Node : public IteratorBase<NodePolicy> {
  public:
    typedef IteratorBase<NodePolicy> NodeIterator;

    Node(const NodeIterator& b)
      : NodeIterator(b) {}
    Node(AliXMLEngine *xml, XMLNodePointer_t nodePtr)
      : NodeIterator(xml, nodePtr) {}
    virtual ~Node() {}

    NodeIterator GetChildBegin() const { return {fXML, fXML->GetChild(&*this)}; }
    NodeIterator GetChildEnd()   const { return {fXML, nullptr}; }
    NodeIterator GetChild(TString name) const {
      return std::find_if(GetChildBegin(), GetChildEnd(), [name](const Node& n) { return (name == n.GetName()); });
    }

    Attr GetAttrBegin() const { return {fXML, fXML->GetFirstAttr(&*this)}; }
    Attr GetAttrEnd()   const { return {fXML, nullptr}; }
    Attr GetAttr(TString name) const {
      return std::find_if(GetAttrBegin(), GetAttrEnd(), [name](const Attr& n) { return (name == n.GetName()); });
    }
  protected:
  private:
  } ;

  Node GetRootNode() { return {this, DocGetRootElement(fXMLDoc)}; }

protected:
private:
  AliXMLEngine();
  AliXMLEngine(const AliXMLEngine&);
  AliXMLEngine& operator=(const AliXMLEngine&);

  XMLDocPointer_t fXMLDoc; //!

  ClassDef(AliXMLEngine, 1);
} ;

#endif // _ALI_XML_ENGINE_H
