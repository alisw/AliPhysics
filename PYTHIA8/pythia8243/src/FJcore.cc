// fjcore -- extracted from FastJet v3.2.1 (http://fastjet.fr)
//
// fjcore constitutes a digest of the main FastJet functionality.
// The files fjcore.hh and fjcore.cc are meant to provide easy access to these 
// core functions, in the form of single files and without the need of a full 
// FastJet installation:
//
//     g++ main.cc fjcore.cc
// 
// with main.cc including fjcore.hh.
//
// A fortran interface, fjcorefortran.cc, is also provided. See the example 
// and the Makefile for instructions.
//
// The results are expected to be identical to those obtained by linking to
// the full FastJet distribution.
//
// NOTE THAT, IN ORDER TO MAKE IT POSSIBLE FOR FJCORE AND THE FULL FASTJET
// TO COEXIST, THE FORMER USES THE "fjcore" NAMESPACE INSTEAD OF "fastjet". 
//
// In particular, fjcore provides:
//
//   - access to all native pp and ee algorithms, kt, anti-kt, C/A.
//     For C/A, the NlnN method is available, while anti-kt and kt
//     are limited to the N^2 one (still the fastest for N < 100k particles)
//   - access to selectors, for implementing cuts and selections
//   - access to all functionalities related to pseudojets (e.g. a jet's
//     structure or user-defined information)
//
// Instead, it does NOT provide:
//
//   - jet areas functionality
//   - background estimation
//   - access to other algorithms via plugins
//   - interface to CGAL
//   - fastjet tools, e.g. filters, taggers
//
// If these functionalities are needed, the full FastJet installation must be
// used. The code will be fully compatible, with the sole replacement of the
// header files and of the fjcore namespace with the fastjet one.
//
// fjcore.hh and fjcore.cc are not meant to be human-readable.
// For documentation, see the full FastJet manual and doxygen at http://fastjet.fr
//
// Like FastJet, fjcore is released under the terms of the GNU General Public
// License version 2 (GPLv2). If you use this code as part of work towards a
// scientific publication, whether directly or contained within another program
// (e.g. Delphes, MadGraph, SpartyJet, Rivet, LHC collaboration software frameworks, 
// etc.), you should include a citation to
// 
//   EPJC72(2012)1896 [arXiv:1111.6097] (FastJet User Manual)
//   and, optionally, Phys.Lett.B641 (2006) 57 [arXiv:hep-ph/0512210]
//
// Copyright (c) 2005-2016, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//
//#include "fjcore.hh"
// For inclusion in Pythia8 line above is replaced by line below.
#include "Pythia8/FJcore.h"
#ifndef __FJCORE_VERSION_HH__
#define __FJCORE_VERSION_HH__
#include<string>
FJCORE_BEGIN_NAMESPACE
const char* fastjet_version = FJCORE_PACKAGE_VERSION;
FJCORE_END_NAMESPACE
#endif // __FJCORE_VERSION_HH__
#ifndef __FJCORE_CLUSTERQUENCE_N2_ICC__
#define __FJCORE_CLUSTERQUENCE_N2_ICC__
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
template<class BJ> void ClusterSequence::_simple_N2_cluster() {
  int n = _jets.size();
  BJ * briefjets = new BJ[n];
  BJ * jetA = briefjets, * jetB;
  for (int i = 0; i< n; i++) {
    _bj_set_jetinfo(jetA, i);
    jetA++; // move on to next entry of briefjets
  }
  BJ * tail = jetA; // a semaphore for the end of briefjets
  BJ * head = briefjets; // a nicer way of naming start
  for (jetA = head + 1; jetA != tail; jetA++) {
    _bj_set_NN_crosscheck(jetA, head, jetA);
  }
  double * diJ = new double[n];
  jetA = head;
  for (int i = 0; i < n; i++) {
    diJ[i] = _bj_diJ(jetA);
    jetA++; // have jetA follow i
  }
  int history_location = n-1;
  while (tail != head) {
    double diJ_min = diJ[0];
    int diJ_min_jet = 0;
    for (int i = 1; i < n; i++) {
      if (diJ[i] < diJ_min) {diJ_min_jet = i; diJ_min  = diJ[i];}
    }
    history_location++;
    jetA = & briefjets[diJ_min_jet];
    jetB = static_cast<BJ *>(jetA->NN);
    diJ_min *= _invR2; 
    if (jetB != NULL) {
      if (jetA < jetB) {std::swap(jetA,jetB);}
      int nn; // new jet index
      _do_ij_recombination_step(jetA->_jets_index, jetB->_jets_index, diJ_min, nn);
      _bj_set_jetinfo(jetB, nn);
    } else {
      _do_iB_recombination_step(jetA->_jets_index, diJ_min);
    }
    tail--; n--;
    *jetA = *tail;
    diJ[jetA - head] = diJ[tail-head];
    for (BJ * jetI = head; jetI != tail; jetI++) {
      if (jetI->NN == jetA || jetI->NN == jetB) {
	_bj_set_NN_nocross(jetI, head, tail);
	diJ[jetI-head] = _bj_diJ(jetI); // update diJ 
      } 
      if (jetB != NULL) {
	double dist = _bj_dist(jetI,jetB);
	if (dist < jetI->NN_dist) {
	  if (jetI != jetB) {
	    jetI->NN_dist = dist;
	    jetI->NN = jetB;
	    diJ[jetI-head] = _bj_diJ(jetI); // update diJ...
	  }
	}
	if (dist < jetB->NN_dist) {
	  if (jetI != jetB) {
	    jetB->NN_dist = dist;
	    jetB->NN      = jetI;}
	}
      }
      if (jetI->NN == tail) {jetI->NN = jetA;}
    }
    if (jetB != NULL) {diJ[jetB-head] = _bj_diJ(jetB);}
  }
  delete[] diJ;
  delete[] briefjets;
}
FJCORE_END_NAMESPACE
#endif // __FJCORE_CLUSTERQUENCE_N2_ICC__
#ifndef __FJCORE_DYNAMICNEARESTNEIGHBOURS_HH__
#define __FJCORE_DYNAMICNEARESTNEIGHBOURS_HH__
#include<vector>
#include<string>
#include<iostream>
#include<sstream>
#include<cassert>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
class EtaPhi {
public:
  double first, second;
  EtaPhi() {}
  EtaPhi(double a, double b) {first = a; second = b;}
  void sanitize() {    
    if (second <  0)     second += twopi; 
    if (second >= twopi) second -= twopi;
  }
};
class DnnError : public Error {
public:
  DnnError(const std::string & message_in) : Error(message_in) {}
};
class DynamicNearestNeighbours {
public:
  virtual int NearestNeighbourIndex(const int ii) const = 0;
  virtual double NearestNeighbourDistance(const int ii) const = 0;
  virtual bool Valid(const int index) const = 0;
  virtual void RemoveAndAddPoints(const std::vector<int> & indices_to_remove,
			  const std::vector<EtaPhi> & points_to_add,
			  std::vector<int> & indices_added,
			  std::vector<int> & indices_of_updated_neighbours) = 0;
  inline void RemovePoint (const int index,
			   std::vector<int> & indices_of_updated_neighbours) {
    std::vector<int> indices_added;
    std::vector<EtaPhi> points_to_add;
    std::vector<int> indices_to_remove(1);
    indices_to_remove[0] = index;
    RemoveAndAddPoints(indices_to_remove, points_to_add, indices_added,
		       indices_of_updated_neighbours
		       );};
  inline void RemoveCombinedAddCombination(
			const int index1, const int index2,
			const EtaPhi & newpoint,
			int & index3,
			std::vector<int> & indices_of_updated_neighbours) {
    std::vector<int> indices_added(1);
    std::vector<EtaPhi> points_to_add(1);
    std::vector<int> indices_to_remove(2);
    indices_to_remove[0] = index1;
    indices_to_remove[1] = index2;
    points_to_add[0] = newpoint;
    RemoveAndAddPoints(indices_to_remove, points_to_add, indices_added,
		       indices_of_updated_neighbours
		       );
    index3 = indices_added[0];
  };
  virtual ~DynamicNearestNeighbours () {}
};
FJCORE_END_NAMESPACE
#endif // __FJCORE_DYNAMICNEARESTNEIGHBOURS_HH__
#ifndef __FJCORE_SEARCHTREE_HH__
#define __FJCORE_SEARCHTREE_HH__
#include<vector>
#include<cassert>
#include<cstddef>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
template<class T> class SearchTree {
public:
  class Node;
  class circulator;
  class const_circulator;
  SearchTree(const std::vector<T> & init);
  SearchTree(const std::vector<T> & init, unsigned int max_size);
  void remove(unsigned node_index);
  void remove(typename SearchTree::Node * node);
  void remove(typename SearchTree::circulator & circ);
  circulator insert(const T & value);
  const Node & operator[](int i) const {return _nodes[i];};
  unsigned int size() const {return _nodes.size() - _available_nodes.size();}
  void verify_structure();
  void verify_structure_linear() const;
  void verify_structure_recursive(const Node * , const Node * , const Node * ) const;
  void print_elements();
#ifdef __FJCORE_SEARCHTREE_TRACK_DEPTH
  inline unsigned int max_depth() const {return _max_depth;};
#else
  inline unsigned int max_depth() const {return 0;};
#endif
  int loc(const Node * node) const ;
  Node * _find_predecessor(const Node *);
  Node * _find_successor(const Node *);
  const Node & operator[](unsigned int i) const {return _nodes[i];};
  const_circulator somewhere() const;
  circulator somewhere();
private:
  void _initialize(const std::vector<T> & init);
  std::vector<Node> _nodes;
  std::vector<Node *> _available_nodes;
  Node * _top_node;
  unsigned int _n_removes;
  void _do_initial_connections(unsigned int this_one, unsigned int scale,
			       unsigned int left_edge, unsigned int right_edge,
			       unsigned int depth);
#ifdef __FJCORE_SEARCHTREE_TRACK_DEPTH
  unsigned int _max_depth;
#endif
};
template<class T> class SearchTree<T>::Node{
public:
  Node() {}; /// default constructor
  bool treelinks_null() const {
    return ((parent==0) && (left==0) && (right==0));};
  inline void nullify_treelinks() {
    parent = NULL; 
    left   = NULL; 
    right  = NULL;
  };
  void reset_parents_link_to_me(Node * XX);
  T      value;
  Node * left;
  Node * right;
  Node * parent;
  Node * successor;
  Node * predecessor;
};
template<class T> void SearchTree<T>::Node::reset_parents_link_to_me(typename SearchTree<T>::Node * XX) {
  if (parent == NULL) {return;}
  if (parent->right == this) {parent->right = XX;}
  else {parent->left = XX;}
}
template<class T> class SearchTree<T>::circulator{
public:
  friend class SearchTree<T>::const_circulator;
  friend class SearchTree<T>;
  circulator() : _node(NULL) {}
  circulator(Node * node) : _node(node) {}
  const T * operator->() const {return &(_node->value);}
  T * operator->() {return &(_node->value);}
  const T & operator*() const {return _node->value;}
  T & operator*() {return _node->value;}
  circulator & operator++() {
    _node = _node->successor; 
    return *this;}
  circulator operator++(int) {
    circulator tmp = *this;
    _node = _node->successor; 
    return tmp;}
  circulator & operator--() {
    _node = _node->predecessor; 
    return *this;}
  circulator operator--(int) {
    circulator tmp = *this;
    _node = _node->predecessor; 
    return tmp;}
  circulator next() const {
    return circulator(_node->successor);}
  circulator previous() const {
    return circulator(_node->predecessor);}
  bool operator!=(const circulator & other) const {return other._node != _node;}
  bool operator==(const circulator & other) const {return other._node == _node;}
private:
  Node * _node;
};
template<class T> class SearchTree<T>::const_circulator{
public:
  const_circulator() : _node(NULL) {}
  const_circulator(const Node * node) : _node(node) {}
  const_circulator(const circulator & circ) :_node(circ._node) {}
  const T * operator->() {return &(_node->value);}
  const T & operator*() const {return _node->value;}
  const_circulator & operator++() {
    _node = _node->successor; 
    return *this;}
  const_circulator operator++(int) {
    const_circulator tmp = *this;
    _node = _node->successor; 
    return tmp;}
  const_circulator & operator--() {
    _node = _node->predecessor; 
    return *this;}
  const_circulator operator--(int) {
    const_circulator tmp = *this;
    _node = _node->predecessor; 
    return tmp;}
  const_circulator next() const {
    return const_circulator(_node->successor);}
  const_circulator previous() const {
    return const_circulator(_node->predecessor);}
  bool operator!=(const const_circulator & other) const {return other._node != _node;}
  bool operator==(const const_circulator & other) const {return other._node == _node;}
private:
  const Node * _node;
};
template<class T> SearchTree<T>::SearchTree(const std::vector<T> & init,
					    unsigned int max_size) :
  _nodes(max_size) {
  _available_nodes.reserve(max_size);
  _available_nodes.resize(max_size - init.size());
  for (unsigned int i = init.size(); i < max_size; i++) {
    _available_nodes[i-init.size()] = &(_nodes[i]);
  }
  _initialize(init);
}
template<class T> SearchTree<T>::SearchTree(const std::vector<T> & init) :
  _nodes(init.size()), _available_nodes(0) {
  _available_nodes.reserve(init.size());
  _initialize(init);
}
template<class T> void SearchTree<T>::_initialize(const std::vector<T> & init) {
  _n_removes = 0;
  unsigned n = init.size();
  assert(n>=1);
#ifdef __FJCORE_SEARCHTREE_TRACK_DEPTH
  _max_depth     = 0;
#endif
  for (unsigned int i = 1; i<n; i++) {
    assert(!(init[i] < init[i-1]));
  }
  for(unsigned int i = 0; i < n; i++) {
    _nodes[i].value = init[i];
    _nodes[i].predecessor = (& (_nodes[i])) - 1;
    _nodes[i].successor   = (& (_nodes[i])) + 1;
    _nodes[i].nullify_treelinks();
  }
  _nodes[0].predecessor = (& (_nodes[n-1]));
  _nodes[n-1].successor = (& (_nodes[0]));
  unsigned int scale = (n+1)/2;
  unsigned int top   = std::min(n-1,scale);
  _nodes[top].parent = NULL;
  _top_node = &(_nodes[top]);
  _do_initial_connections(top, scale, 0, n, 0);
}
template<class T> inline  int SearchTree<T>::loc(const Node * node) const {return node == NULL? 
      -999 : node - &(_nodes[0]);}
template<class T> void SearchTree<T>::_do_initial_connections(
                                         unsigned int this_one, 
					 unsigned int scale,
					 unsigned int left_edge,
					 unsigned int right_edge,
					 unsigned int depth
					 ) {
#ifdef __FJCORE_SEARCHTREE_TRACK_DEPTH
  _max_depth = max(depth, _max_depth);
#endif
  unsigned int ref_new_scale = (scale+1)/2;
  unsigned new_scale = ref_new_scale;
  bool     did_child  = false;
  while(true) {
    int left = this_one - new_scale; // be careful here to use signed int...
    if (left >= static_cast<int>(left_edge) 
	                && _nodes[left].treelinks_null() ) {
      _nodes[left].parent = &(_nodes[this_one]);
      _nodes[this_one].left = &(_nodes[left]);
      _do_initial_connections(left, new_scale, left_edge, this_one, depth+1);
      did_child = true;
      break;
    }
    unsigned int old_new_scale = new_scale;
    new_scale = (old_new_scale + 1)/2;
    if (new_scale == old_new_scale) break;
  }
  if (!did_child) {_nodes[this_one].left = NULL;}
  new_scale = ref_new_scale;
  did_child  = false;
  while(true) {
    unsigned int right = this_one + new_scale;
    if (right < right_edge  && _nodes[right].treelinks_null()) {
      _nodes[right].parent = &(_nodes[this_one]);
      _nodes[this_one].right = &(_nodes[right]);
      _do_initial_connections(right, new_scale, this_one+1,right_edge,depth+1);
      did_child = true;
      break;
    }
    unsigned int old_new_scale = new_scale;
    new_scale = (old_new_scale + 1)/2;
    if (new_scale == old_new_scale) break;
  }
  if (!did_child) {_nodes[this_one].right = NULL;}
}
template<class T> void SearchTree<T>::remove(unsigned int node_index) {
  remove(&(_nodes[node_index]));
}
template<class T> void SearchTree<T>::remove(circulator & circ) {
  remove(circ._node);
}
template<class T> void SearchTree<T>::remove(typename SearchTree<T>::Node * node) {
  assert(size() > 1); // switch this to throw...?
  assert(!node->treelinks_null());
  node->predecessor->successor = node->successor;
  node->successor->predecessor = node->predecessor;
  if (node->left == NULL && node->right == NULL) {
    node->reset_parents_link_to_me(NULL); 
  } else if (node->left != NULL && node->right == NULL){
    node->reset_parents_link_to_me(node->left);
    node->left->parent = node->parent;         
    if (_top_node == node) {_top_node = node->left;}
  } else if (node->left == NULL && node->right != NULL){
    node->reset_parents_link_to_me(node->right);
    node->right->parent = node->parent;   
    if (_top_node == node) {_top_node = node->right;}
  } else {
    Node * replacement;
    bool use_predecessor = (_n_removes % 2 == 1);
    if (use_predecessor) {
      replacement = node->predecessor;
      assert(replacement->right == NULL); // guaranteed if it's our predecessor
      if (replacement != node->left) {
	if (replacement->left != NULL) {
	  replacement->left->parent = replacement->parent;}
	replacement->reset_parents_link_to_me(replacement->left);
	replacement->left   = node->left;
      }
      replacement->parent = node->parent;
      replacement->right  = node->right;
    } else {
      replacement = node->successor;
      assert(replacement->left == NULL); // guaranteed if it's our successor
      if (replacement != node->right) {
	if (replacement->right != NULL) {
	  replacement->right->parent = replacement->parent;}
	replacement->reset_parents_link_to_me(replacement->right);
	replacement->right  = node->right;
      }
      replacement->parent = node->parent;
      replacement->left   = node->left;
    }
    node->reset_parents_link_to_me(replacement);
    if (node->left  != replacement) {node->left->parent  = replacement;}
    if (node->right != replacement) {node->right->parent = replacement;}
    if (_top_node == node) {_top_node = replacement;}
  }
  node->nullify_treelinks();
  node->predecessor = NULL;
  node->successor   = NULL;
  _n_removes++;
  _available_nodes.push_back(node);
}
template<class T> typename SearchTree<T>::circulator SearchTree<T>::insert(const T & value) {
  assert(_available_nodes.size() > 0);
  Node * node = _available_nodes.back();
  _available_nodes.pop_back();
  node->value = value;
  Node * location = _top_node;
  Node * old_location = NULL;
  bool             on_left = true; // (init not needed -- but soothes g++4)
#ifdef __FJCORE_SEARCHTREE_TRACK_DEPTH
  unsigned int depth = 0;
#endif
  while(location != NULL) {
#ifdef __FJCORE_SEARCHTREE_TRACK_DEPTH
    depth++;
#endif
    old_location = location;
    on_left = value < location->value;
    if (on_left) {location = location->left;}
    else {location = location->right;}
  }
#ifdef __FJCORE_SEARCHTREE_TRACK_DEPTH
  _max_depth = max(depth, _max_depth);
#endif
  node->parent = old_location;
  if (on_left) {node->parent->left = node;} 
  else {node->parent->right = node;}
  node->left = NULL;
  node->right = NULL;
  node->predecessor = _find_predecessor(node);
  if (node->predecessor != NULL) {
    node->successor = node->predecessor->successor;
    node->predecessor->successor = node;
    node->successor->predecessor = node;
  } else {
    node->successor = _find_successor(node);
    assert(node->successor != NULL); // can only happen if we're sole element 
    node->predecessor = node->successor->predecessor;
    node->successor->predecessor = node;
    node->predecessor->successor = node;
  }
  return circulator(node);
}
template<class T> void SearchTree<T>::verify_structure() {
  verify_structure_linear();
  const Node * left_limit = _top_node;
  while (left_limit->left != NULL) {left_limit = left_limit->left;}
  const Node * right_limit = _top_node;
  while (right_limit->right != NULL) {right_limit = right_limit->right;}
  verify_structure_recursive(_top_node, left_limit, right_limit);
}
template<class T> void SearchTree<T>::verify_structure_recursive(
		      const typename SearchTree<T>::Node * element, 
		      const typename SearchTree<T>::Node * left_limit,
		      const typename SearchTree<T>::Node * right_limit)  const {
  assert(!(element->value < left_limit->value));
  assert(!(right_limit->value < element->value));
  const Node * left = element->left;
  if (left != NULL) {
    assert(!(element->value < left->value));
    if (left != left_limit) {
      verify_structure_recursive(left, left_limit, element);}
  }
  const Node * right = element->right;
  if (right != NULL) {
    assert(!(right->value < element->value));
    if (right != right_limit) {
      verify_structure_recursive(right, element, right_limit);}
  }
}
template<class T> void SearchTree<T>::verify_structure_linear() const {
  unsigned n_top = 0;
  unsigned n_null = 0;
  for(unsigned i = 0; i < _nodes.size(); i++) {
    const typename SearchTree<T>::Node * node = &(_nodes[i]);
    if (node->treelinks_null()) {n_null++; continue;}
    if (node->parent == NULL) {
      n_top++;
    } else {
      assert((node->parent->left == node) ^ (node->parent->right == node));
    }
    if (node->left != NULL) {
      assert(!(node->value < node->left->value ));}
    if (node->right != NULL) {
      assert(!(node->right->value < node->value ));}
  }
  assert(n_top == 1 || (n_top == 0 && size() <= 1) );
  assert(n_null == _available_nodes.size() ||
	 (n_null == _available_nodes.size() + 1 && size() == 1));
}
template<class T> typename SearchTree<T>::Node * SearchTree<T>::_find_predecessor(const typename SearchTree<T>::Node * node) {
  typename SearchTree<T>::Node * newnode;
  if (node->left != NULL) {
    newnode = node->left;
    while(newnode->right != NULL) {newnode = newnode->right;}
    return newnode;
  } else {
    const typename SearchTree<T>::Node * lastnode = node;
    newnode = node->parent;
    while(newnode != NULL) {
      if (newnode->right == lastnode) {return newnode;}
      lastnode = newnode;
      newnode = newnode->parent;
    }
    return newnode;
  }
}
template<class T> typename SearchTree<T>::Node * SearchTree<T>::_find_successor(const typename SearchTree<T>::Node * node) {
  typename SearchTree<T>::Node * newnode;
  if (node->right != NULL) {
    newnode = node->right;
    while(newnode->left != NULL) {newnode = newnode->left;}
    return newnode;
  } else {
    const typename SearchTree<T>::Node * lastnode = node;
    newnode = node->parent;
    while(newnode != NULL) {
      if (newnode->left == lastnode) {return newnode;}
      lastnode = newnode;
      newnode = newnode->parent;
    }
    return newnode;
  }
}
template<class T> void SearchTree<T>::print_elements() {
  typename SearchTree<T>::Node * base_node = &(_nodes[0]);
  typename SearchTree<T>::Node * node = base_node;
  int n = _nodes.size();
  for(; node - base_node < n ; node++) {
    printf("%4d parent:%4d left:%4d right:%4d pred:%4d succ:%4d value:%10.6f\n",loc(node), loc(node->parent), loc(node->left), loc(node->right), loc(node->predecessor),loc(node->successor),node->value);
  }
}
template<class T> typename SearchTree<T>::circulator SearchTree<T>::somewhere() {
  return circulator(_top_node);
}
template<class T> typename SearchTree<T>::const_circulator SearchTree<T>::somewhere() const {
  return const_circulator(_top_node);
}
FJCORE_END_NAMESPACE
#endif // __FJCORE_SEARCHTREE_HH__
#ifndef __FJCORE_MINHEAP__HH__
#define __FJCORE_MINHEAP__HH__
#include<vector>
#include<cassert>
#include<memory>
#include<limits>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
class MinHeap {
public:
  MinHeap (const std::vector<double> & values, unsigned int max_size) :
    _heap(max_size) {initialise(values);}
  MinHeap (unsigned int max_size) : _heap(max_size) {}
  MinHeap (const std::vector<double> & values) :
    _heap(values.size()) {initialise(values);}
  void initialise(const std::vector<double> & values);
  inline unsigned int minloc() const {
    return (_heap[0].minloc) - &(_heap[0]);}
  inline double       minval() const {return _heap[0].minloc->value;}
  inline double operator[](int i) const {return _heap[i].value;}
  void remove(unsigned int loc) {
    update(loc,std::numeric_limits<double>::max());};
  void update(unsigned int, double);
private:
  struct ValueLoc{
    double value;
    ValueLoc * minloc;
  };
  std::vector<ValueLoc> _heap;
};
FJCORE_END_NAMESPACE
#endif // __FJCORE_MINHEAP__HH__
#ifndef __FJCORE_CLOSESTPAIR2DBASE__HH__
#define __FJCORE_CLOSESTPAIR2DBASE__HH__
#include<vector>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
class Coord2D {
public:
  double x, y;
  Coord2D() : x(0.0), y(0.0) {};
  Coord2D(double a, double b): x(a), y(b) {};
  Coord2D operator-(const Coord2D & other) const {
    return Coord2D(x - other.x,  y - other.y);};
  Coord2D operator+(const Coord2D & other) const {
    return Coord2D(x + other.x,  y + other.y);};
  Coord2D operator*(double factor) const {return Coord2D(factor*x,factor*y);};
  friend Coord2D operator*(double factor, const Coord2D & coord) {
    return Coord2D(factor*coord.x,factor*coord.y);
  }
  Coord2D operator/(double divisor) const {
    return Coord2D(x / divisor,  y / divisor);};
  friend double distance2(const Coord2D & a, const Coord2D & b) {
    double dx = a.x - b.x, dy = a.y-b.y;
    return dx*dx+dy*dy;
  };
  double distance2(const Coord2D & b) const {
    double dx = x - b.x, dy = y-b.y;
    return dx*dx+dy*dy;
  };
};
class ClosestPair2DBase {
public:
  virtual void closest_pair(unsigned int & ID1, unsigned int & ID2, 
			    double & distance2) const = 0;
  virtual void remove(unsigned int ID) = 0;
  virtual unsigned int insert(const Coord2D & position) = 0;
  virtual unsigned int replace(unsigned int ID1, unsigned int ID2, 
			       const Coord2D & position) {
    remove(ID1); 
    remove(ID2); 
    unsigned new_ID = insert(position);
    return(new_ID);
  };
  virtual void replace_many(const std::vector<unsigned int> & IDs_to_remove,
		       const std::vector<Coord2D> & new_positions,
		       std::vector<unsigned int> & new_IDs) {
    for(unsigned i = 0; i < IDs_to_remove.size(); i++) {
      remove(IDs_to_remove[i]);}
    new_IDs.resize(0);
    for(unsigned i = 0; i < new_positions.size(); i++) {
      new_IDs.push_back(insert(new_positions[i]));}
  }
  virtual unsigned int size() = 0;
  virtual ~ClosestPair2DBase() {};
};
FJCORE_END_NAMESPACE
#endif // __FJCORE_CLOSESTPAIR2DBASE__HH__
#ifndef __FJCORE_CLOSESTPAIR2D__HH__
#define __FJCORE_CLOSESTPAIR2D__HH__
#include<vector>
#include<stack>
#include<iostream>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
class ClosestPair2D : public ClosestPair2DBase {
public:
  ClosestPair2D(const std::vector<Coord2D> & positions, 
		const Coord2D & left_corner, const Coord2D & right_corner) {
    _initialize(positions, left_corner, right_corner, positions.size());
  };
  ClosestPair2D(const std::vector<Coord2D> & positions, 
		const Coord2D & left_corner, const Coord2D & right_corner,
		const unsigned int max_size) {
    _initialize(positions, left_corner, right_corner, max_size);
  };
  void closest_pair(unsigned int & ID1, unsigned int & ID2, 
		    double & distance2) const;
  void remove(unsigned int ID);
  unsigned int insert(const Coord2D &);
  virtual unsigned int replace(unsigned int ID1, unsigned int ID2, 
			       const Coord2D & position);
  virtual void replace_many(const std::vector<unsigned int> & IDs_to_remove,
			    const std::vector<Coord2D> & new_positions,
			    std::vector<unsigned int> & new_IDs);
  inline void print_tree_depths(std::ostream & outdev) const {
    outdev    << _trees[0]->max_depth() << " "
	      << _trees[1]->max_depth() << " "
	      << _trees[2]->max_depth() << "\n";
  };
  unsigned int size();
private:
  void _initialize(const std::vector<Coord2D> & positions, 
	      const Coord2D & left_corner, const Coord2D & right_corner,
	      const unsigned int max_size);
  static const unsigned int _nshift = 3;
  class Point; // will be defined below
  template<class T> class triplet {
  public:
    inline const T & operator[](unsigned int i) const {return _contents[i];};
    inline       T & operator[](unsigned int i)       {return _contents[i];};
  private:
    T _contents[_nshift];
  };
  class Shuffle {
  public:
    unsigned int x, y;
    Point * point;
    bool operator<(const Shuffle &) const;
    void operator+=(unsigned int shift) {x += shift; y+= shift;};
  };
  typedef SearchTree<Shuffle>     Tree;
  typedef Tree::circulator        circulator;
  typedef Tree::const_circulator  const_circulator;
  triplet<SharedPtr<Tree> >  _trees;
  SharedPtr<MinHeap>     _heap;
  std::vector<Point>     _points;
  std::stack<Point *>    _available_points;
  std::vector<Point *>   _points_under_review;
  static const unsigned int _remove_heap_entry = 1;
  static const unsigned int _review_heap_entry = 2;
  static const unsigned int _review_neighbour  = 4;
  void _add_label(Point * point, unsigned int review_flag);
  void _set_label(Point * point, unsigned int review_flag);
  void _deal_with_points_to_review();
  void _remove_from_search_tree(Point * point_to_remove);
  void _insert_into_search_tree(Point * new_point);
  void _point2shuffle(Point & , Shuffle & , unsigned int shift);
  Coord2D _left_corner;
  double _range;
  int _ID(const Point *) const;
  triplet<unsigned int> _shifts;     // absolute shifts
  triplet<unsigned int> _rel_shifts; // shifts relative to previous shift
  unsigned int _cp_search_range;
};
class ClosestPair2D::Point {
public:
  Coord2D coord;
  Point * neighbour;
  double  neighbour_dist2;
  triplet<circulator> circ;
  unsigned int review_flag;
  double distance2(const Point & other) const {
    return coord.distance2(other.coord);
  };
};
inline bool floor_ln2_less(unsigned x, unsigned y) {
  if (x>y) return false;
  return (x < (x^y)); // beware of operator precedence...
}
inline int ClosestPair2D::_ID(const Point * point) const {
  return point - &(_points[0]);
}
inline unsigned int ClosestPair2D::size() {
  return _points.size() - _available_points.size();
}
FJCORE_END_NAMESPACE
#endif // __FJCORE_CLOSESTPAIR2D__HH__
#ifndef __FJCORE_LAZYTILING9ALT_HH__
#define __FJCORE_LAZYTILING9ALT_HH__
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
const double tile_edge_security_margin=1.0e-7;
class TiledJet {
public:
  double     eta, phi, kt2, NN_dist;
  TiledJet * NN, *previous, * next; 
  int        _jets_index, tile_index;
  bool _minheap_update_needed;
  inline void label_minheap_update_needed() {_minheap_update_needed = true;}
  inline void label_minheap_update_done()   {_minheap_update_needed = false;}
  inline bool minheap_update_needed() const {return _minheap_update_needed;}
};
const int n_tile_neighbours = 9;
class Tile {
public:
  typedef double (Tile::*DistToTileFn)(const TiledJet*) const;
  typedef std::pair<Tile *, DistToTileFn> TileFnPair;
  TileFnPair begin_tiles[n_tile_neighbours]; 
  TileFnPair *  surrounding_tiles; 
  TileFnPair *  RH_tiles;  
  TileFnPair *  end_tiles; 
  TiledJet * head;    
  bool     tagged;    
  bool     use_periodic_delta_phi;
  double max_NN_dist;
  double eta_min, eta_max, phi_min, phi_max;
  double distance_to_centre(const TiledJet *) const {return 0;}
  double distance_to_left(const TiledJet * jet) const {
    double deta = jet->eta - eta_min;
    return deta*deta;
  }
  double distance_to_right(const TiledJet * jet) const {
    double deta = jet->eta - eta_max;
    return deta*deta;
  }
  double distance_to_bottom(const TiledJet * jet) const {
    double dphi = jet->phi - phi_min;
    return dphi*dphi;
  }
  double distance_to_top(const TiledJet * jet) const {
    double dphi = jet->phi - phi_max;
    return dphi*dphi;
  }
  double distance_to_left_top(const TiledJet * jet) const {
    double deta = jet->eta - eta_min;
    double dphi = jet->phi - phi_max;
    return deta*deta + dphi*dphi;
  }
  double distance_to_left_bottom(const TiledJet * jet) const {
    double deta = jet->eta - eta_min;
    double dphi = jet->phi - phi_min;
    return deta*deta + dphi*dphi;
  }
  double distance_to_right_top(const TiledJet * jet) const {
    double deta = jet->eta - eta_max;
    double dphi = jet->phi - phi_max;
    return deta*deta + dphi*dphi;
  }
  double distance_to_right_bottom(const TiledJet * jet) const {
    double deta = jet->eta - eta_max;
    double dphi = jet->phi - phi_min;
    return deta*deta + dphi*dphi;
  }
};
class LazyTiling9Alt {
public:
  LazyTiling9Alt(ClusterSequence & cs);
  void run();
protected:
  ClusterSequence & _cs;
  const std::vector<PseudoJet> & _jets;
  std::vector<Tile> _tiles;
  double _Rparam, _R2, _invR2;
  double _tiles_eta_min, _tiles_eta_max;
  double _tile_size_eta, _tile_size_phi;
  double _tile_half_size_eta, _tile_half_size_phi;
  int    _n_tiles_phi,_tiles_ieta_min,_tiles_ieta_max;
  std::vector<TiledJet *> _jets_for_minheap;
  void _initialise_tiles();
  inline int _tile_index (int ieta, int iphi) const {
    return (ieta-_tiles_ieta_min)*_n_tiles_phi
                  + (iphi+_n_tiles_phi) % _n_tiles_phi;
  }
  void  _bj_remove_from_tiles(TiledJet * const jet);
  int _tile_index(const double eta, const double phi) const;
  void _tj_set_jetinfo(TiledJet * const jet, const int _jets_index);
  void _print_tiles(TiledJet * briefjets ) const;
  void _add_neighbours_to_tile_union(const int tile_index, 
		 std::vector<int> & tile_union, int & n_near_tiles) const;
  void _add_untagged_neighbours_to_tile_union(const int tile_index, 
		 std::vector<int> & tile_union, int & n_near_tiles);
  void _add_untagged_neighbours_to_tile_union_using_max_info(const TiledJet * const jet, 
		 std::vector<int> & tile_union, int & n_near_tiles);
  void _update_jetX_jetI_NN(TiledJet * jetX, TiledJet * jetI, std::vector<TiledJet *> & jets_for_minheap);
  void _set_NN(TiledJet * jetI, std::vector<TiledJet *> & jets_for_minheap);
  template <class J> double _bj_diJ(const J * const jet) const {
    double kt2 = jet->kt2;
    if (jet->NN != NULL) {if (jet->NN->kt2 < kt2) {kt2 = jet->NN->kt2;}}
    return jet->NN_dist * kt2;
  }
  template <class J> inline void _bj_set_jetinfo(
                            J * const jetA, const int _jets_index) const {
    jetA->eta  = _jets[_jets_index].rap();
    jetA->phi  = _jets[_jets_index].phi_02pi();
    jetA->kt2  = _cs.jet_scale_for_algorithm(_jets[_jets_index]);
    jetA->_jets_index = _jets_index;
    jetA->NN_dist = _R2;
    jetA->NN      = NULL;
  }
  template <class J> inline double _bj_dist(
                const J * const jetA, const J * const jetB) const {
    double dphi = std::abs(jetA->phi - jetB->phi);
    double deta = (jetA->eta - jetB->eta);
    if (dphi > pi) {dphi = twopi - dphi;}
    return dphi*dphi + deta*deta;
  }
  template <class J> inline double _bj_dist_not_periodic(
                const J * const jetA, const J * const jetB) const {
    double dphi = jetA->phi - jetB->phi;
    double deta = (jetA->eta - jetB->eta);
    return dphi*dphi + deta*deta;
  }
};
FJCORE_END_NAMESPACE
#endif // __FJCORE_LAZYTILING9ALT_HH__
#ifndef __FJCORE_LAZYTILING9_HH__
#define __FJCORE_LAZYTILING9_HH__
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
template<int NN>
class Tile2Base {
public:
  Tile2Base *   begin_tiles[NN]; 
  Tile2Base **  surrounding_tiles; 
  Tile2Base **  RH_tiles;  
  Tile2Base **  end_tiles; 
  TiledJet * head;    
  bool     tagged;    
  bool     use_periodic_delta_phi;
  double max_NN_dist;
  double eta_centre, phi_centre;
  int jet_count() const {
    int count = 0;
    const TiledJet * jet = head;
    while (jet != 0) {
      count++;
      jet = jet->next;
    }
    return count;
  }
};
typedef Tile2Base<9> Tile2;
class LazyTiling9 {
public:
  LazyTiling9(ClusterSequence & cs);
  void run();
protected:
  ClusterSequence & _cs;
  const std::vector<PseudoJet> & _jets;
  std::vector<Tile2> _tiles;
#ifdef INSTRUMENT2
  int _ncall; // GPS tmp
  int _ncall_dtt; // GPS tmp
#endif // INSTRUMENT2
  double _Rparam, _R2, _invR2;
  double _tiles_eta_min, _tiles_eta_max;
  double _tile_size_eta, _tile_size_phi;
  double _tile_half_size_eta, _tile_half_size_phi;
  int    _n_tiles_phi,_tiles_ieta_min,_tiles_ieta_max;
  std::vector<TiledJet *> _jets_for_minheap;
  void _initialise_tiles();
  inline int _tile_index (int ieta, int iphi) const {
    return (ieta-_tiles_ieta_min)*_n_tiles_phi
                  + (iphi+_n_tiles_phi) % _n_tiles_phi;
  }
  void  _bj_remove_from_tiles(TiledJet * const jet);
  int _tile_index(const double eta, const double phi) const;
  void _tj_set_jetinfo(TiledJet * const jet, const int _jets_index);
  void _print_tiles(TiledJet * briefjets ) const;
  void _add_neighbours_to_tile_union(const int tile_index, 
		 std::vector<int> & tile_union, int & n_near_tiles) const;
  void _add_untagged_neighbours_to_tile_union(const int tile_index, 
		 std::vector<int> & tile_union, int & n_near_tiles);
  void _add_untagged_neighbours_to_tile_union_using_max_info(const TiledJet * const jet, 
		 std::vector<int> & tile_union, int & n_near_tiles);
  double _distance_to_tile(const TiledJet * bj, const Tile2 *) 
#ifdef INSTRUMENT2
    ;
#else
    const;
#endif 
  void _update_jetX_jetI_NN(TiledJet * jetX, TiledJet * jetI, std::vector<TiledJet *> & jets_for_minheap);
  void _set_NN(TiledJet * jetI, std::vector<TiledJet *> & jets_for_minheap);
  template <class J> double _bj_diJ(const J * const jet) const {
    double kt2 = jet->kt2;
    if (jet->NN != NULL) {if (jet->NN->kt2 < kt2) {kt2 = jet->NN->kt2;}}
    return jet->NN_dist * kt2;
  }
  template <class J> inline void _bj_set_jetinfo(
                            J * const jetA, const int _jets_index) const {
    jetA->eta  = _jets[_jets_index].rap();
    jetA->phi  = _jets[_jets_index].phi_02pi();
    jetA->kt2  = _cs.jet_scale_for_algorithm(_jets[_jets_index]);
    jetA->_jets_index = _jets_index;
    jetA->NN_dist = _R2;
    jetA->NN      = NULL;
  }
  template <class J> inline double _bj_dist(
                const J * const jetA, const J * const jetB) 
#ifdef INSTRUMENT2
    {
    _ncall++; // GPS tmp
#else
    const {
#endif 
    double dphi = std::abs(jetA->phi - jetB->phi);
    double deta = (jetA->eta - jetB->eta);
    if (dphi > pi) {dphi = twopi - dphi;}
    return dphi*dphi + deta*deta;
  }
  template <class J> inline double _bj_dist_not_periodic(
                const J * const jetA, const J * const jetB)
#ifdef INSTRUMENT2
    {
    _ncall++; // GPS tmp
#else
    const {
#endif 
    double dphi = jetA->phi - jetB->phi;
    double deta = (jetA->eta - jetB->eta);
    return dphi*dphi + deta*deta;
  }
};
FJCORE_END_NAMESPACE
#endif // __FJCORE_LAZYTILING9_HH__
#ifndef __FJCORE_LAZYTILING25_HH__
#define __FJCORE_LAZYTILING25_HH__
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
typedef Tile2Base<25> Tile25;
class LazyTiling25 {
public:
  LazyTiling25(ClusterSequence & cs);
  void run();
protected:
  ClusterSequence & _cs;
  const std::vector<PseudoJet> & _jets;
  std::vector<Tile25> _tiles;
#ifdef INSTRUMENT2
  int _ncall; // GPS tmp
  int _ncall_dtt; // GPS tmp
#endif // INSTRUMENT2
  double _Rparam, _R2, _invR2;
  double _tiles_eta_min, _tiles_eta_max;
  double _tile_size_eta, _tile_size_phi;
  double _tile_half_size_eta, _tile_half_size_phi;
  int    _n_tiles_phi,_tiles_ieta_min,_tiles_ieta_max;
  std::vector<TiledJet *> _jets_for_minheap;
  void _initialise_tiles();
  inline int _tile_index (int ieta, int iphi) const {
    return (ieta-_tiles_ieta_min)*_n_tiles_phi
                  + (iphi+_n_tiles_phi) % _n_tiles_phi;
  }
  void  _bj_remove_from_tiles(TiledJet * const jet);
  int _tile_index(const double eta, const double phi) const;
  void _tj_set_jetinfo(TiledJet * const jet, const int _jets_index);
  void _print_tiles(TiledJet * briefjets ) const;
  void _add_neighbours_to_tile_union(const int tile_index, 
		 std::vector<int> & tile_union, int & n_near_tiles) const;
  void _add_untagged_neighbours_to_tile_union(const int tile_index, 
		 std::vector<int> & tile_union, int & n_near_tiles);
  void _add_untagged_neighbours_to_tile_union_using_max_info(const TiledJet * const jet, 
		 std::vector<int> & tile_union, int & n_near_tiles);
  double _distance_to_tile(const TiledJet * bj, const Tile25 *) 
#ifdef INSTRUMENT2
    ;
#else
    const;
#endif 
  void _update_jetX_jetI_NN(TiledJet * jetX, TiledJet * jetI, std::vector<TiledJet *> & jets_for_minheap);
  void _set_NN(TiledJet * jetI, std::vector<TiledJet *> & jets_for_minheap);
  template <class J> double _bj_diJ(const J * const jet) const {
    double kt2 = jet->kt2;
    if (jet->NN != NULL) {if (jet->NN->kt2 < kt2) {kt2 = jet->NN->kt2;}}
    return jet->NN_dist * kt2;
  }
  template <class J> inline void _bj_set_jetinfo(
                            J * const jetA, const int _jets_index) const {
    jetA->eta  = _jets[_jets_index].rap();
    jetA->phi  = _jets[_jets_index].phi_02pi();
    jetA->kt2  = _cs.jet_scale_for_algorithm(_jets[_jets_index]);
    jetA->_jets_index = _jets_index;
    jetA->NN_dist = _R2;
    jetA->NN      = NULL;
  }
  template <class J> inline double _bj_dist(
                const J * const jetA, const J * const jetB) 
#ifdef INSTRUMENT2
    {
    _ncall++; // GPS tmp
#else
    const {
#endif 
    double dphi = std::abs(jetA->phi - jetB->phi);
    double deta = (jetA->eta - jetB->eta);
    if (dphi > pi) {dphi = twopi - dphi;}
    return dphi*dphi + deta*deta;
  }
  template <class J> inline double _bj_dist_not_periodic(
                const J * const jetA, const J * const jetB)
#ifdef INSTRUMENT2
    {
    _ncall++; // GPS tmp
#else
    const {
#endif 
    double dphi = jetA->phi - jetB->phi;
    double deta = (jetA->eta - jetB->eta);
    return dphi*dphi + deta*deta;
  }
};
FJCORE_END_NAMESPACE
#endif // __FJCORE_LAZYTILING25_HH__
#ifndef __FJCORE_TILINGEXTENT_HH__
#define __FJCORE_TILINGEXTENT_HH__
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
class TilingExtent {
public:
  TilingExtent(ClusterSequence & cs);
  TilingExtent(const std::vector<PseudoJet> &particles);
  double minrap() const {return _minrap;}
  double maxrap() const {return _maxrap;}
  double sum_of_binned_squared_multiplicity() const {return _cumul2;}
private:
  double _minrap, _maxrap, _cumul2;
  void _determine_rapidity_extent(const std::vector<PseudoJet> & particles);
};
FJCORE_END_NAMESPACE      // defined in fastjet/internal/base.hh
#endif // __FJCORE_TILINGEXTENT_HH__
#include<limits>
#include<iostream>
#include<iomanip>
#include<algorithm>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
const unsigned int twopow31      = 2147483648U;
using namespace std;
void ClosestPair2D::_point2shuffle(Point & point, Shuffle & shuffle, 
				  unsigned int shift) {
  Coord2D renorm_point = (point.coord - _left_corner)/_range;
  assert(renorm_point.x >=0);
  assert(renorm_point.x <=1);
  assert(renorm_point.y >=0);
  assert(renorm_point.y <=1);
  shuffle.x = static_cast<unsigned int>(twopow31 * renorm_point.x) + shift;
  shuffle.y = static_cast<unsigned int>(twopow31 * renorm_point.y) + shift;
  shuffle.point = &point;
}
bool ClosestPair2D::Shuffle::operator<(const Shuffle & q) const {
  if (floor_ln2_less(x ^ q.x, y ^ q.y)) {
    return (y < q.y);
  } else {
    return (x < q.x);
  }
}
void ClosestPair2D::_initialize(const std::vector<Coord2D> & positions, 
			     const Coord2D & left_corner, 
			     const Coord2D & right_corner,
			     unsigned int max_size) {
  unsigned int n_positions = positions.size();
  assert(max_size >= n_positions);
  _points.resize(max_size);
  for (unsigned int i = n_positions; i < max_size; i++) {
    _available_points.push(&(_points[i]));
  }
  _left_corner = left_corner;
  _range       = max((right_corner.x - left_corner.x),
		     (right_corner.y - left_corner.y));
  vector<Shuffle> shuffles(n_positions);
  for (unsigned int i = 0; i < n_positions; i++) {
    _points[i].coord = positions[i];
    _points[i].neighbour_dist2 = numeric_limits<double>::max();
    _points[i].review_flag = 0;
    _point2shuffle(_points[i], shuffles[i], 0);
  }
  for (unsigned ishift = 0; ishift < _nshift; ishift++) {
   _shifts[ishift] = static_cast<unsigned int>(((twopow31*1.0)*ishift)/_nshift);
    if (ishift == 0) {_rel_shifts[ishift] = 0;}
    else {_rel_shifts[ishift] = _shifts[ishift] - _shifts[ishift-1];}
  }
  _cp_search_range = 30;
  _points_under_review.reserve(_nshift * _cp_search_range);
  for (unsigned int ishift = 0; ishift < _nshift; ishift++) {
    if (ishift > 0) {
      unsigned rel_shift = _rel_shifts[ishift];
      for (unsigned int i = 0; i < shuffles.size(); i++) {
	shuffles[i] += rel_shift; }
    }
    sort(shuffles.begin(), shuffles.end());
    _trees[ishift] = SharedPtr<Tree>(new Tree(shuffles, max_size));
    circulator circ = _trees[ishift]->somewhere(), start=circ;
    unsigned int CP_range = min(_cp_search_range, n_positions-1);
    do {
      Point * this_point = circ->point;
      this_point->circ[ishift] = circ;
      circulator other = circ;
      for (unsigned i=0; i < CP_range; i++) {
	++other;
	double dist2 = this_point->distance2(*other->point);
	if (dist2 < this_point->neighbour_dist2) {
	  this_point->neighbour_dist2 = dist2;
	  this_point->neighbour       = other->point;
	}
      }
    } while (++circ != start);
  }
  vector<double> mindists2(n_positions);
  for (unsigned int i = 0; i < n_positions; i++) {
    mindists2[i] = _points[i].neighbour_dist2;}
  _heap = SharedPtr<MinHeap>(new MinHeap(mindists2, max_size));
}
void ClosestPair2D::closest_pair(unsigned int & ID1, unsigned int & ID2, 
				 double & distance2) const {
  ID1 = _heap->minloc();
  ID2 = _ID(_points[ID1].neighbour);
  distance2 = _points[ID1].neighbour_dist2;
  if (ID1 > ID2) std::swap(ID1,ID2);
}
inline void ClosestPair2D::_add_label(Point * point, unsigned int review_flag) {
  if (point->review_flag == 0) _points_under_review.push_back(point);
  point->review_flag |= review_flag;
}
inline void ClosestPair2D::_set_label(Point * point, unsigned int review_flag) {
  if (point->review_flag == 0) _points_under_review.push_back(point);
  point->review_flag = review_flag;
}
void ClosestPair2D::remove(unsigned int ID) {
  Point * point_to_remove = & (_points[ID]);
  _remove_from_search_tree(point_to_remove);
  _deal_with_points_to_review();
}
void ClosestPair2D::_remove_from_search_tree(Point * point_to_remove) {
  _available_points.push(point_to_remove);
  _set_label(point_to_remove, _remove_heap_entry);
  unsigned int CP_range = min(_cp_search_range, size()-1);
  for (unsigned int ishift = 0; ishift < _nshift; ishift++) {
    circulator removed_circ = point_to_remove->circ[ishift];
    circulator right_end = removed_circ.next();
    _trees[ishift]->remove(removed_circ);
    circulator left_end  = right_end, orig_right_end = right_end;
    for (unsigned int i = 0; i < CP_range; i++) {left_end--;}
    if (size()-1 < _cp_search_range) {
      left_end--; right_end--;
    }
    do {
      Point * left_point = left_end->point;
      if (left_point->neighbour == point_to_remove) {
	// we'll deal with it later...
	_add_label(left_point, _review_neighbour);
      } else {
	// check to see if right point has become its closest neighbour
	double dist2 = left_point->distance2(*right_end->point);
	if (dist2 < left_point->neighbour_dist2) {
	  left_point->neighbour = right_end->point;
	  left_point->neighbour_dist2 = dist2;
	  // NB: (LESSER) REVIEW NEEDED HERE TOO...
	  _add_label(left_point, _review_heap_entry);
	}
      }
      ++right_end;
    } while (++left_end != orig_right_end);
  } // ishift...
}
void ClosestPair2D::_deal_with_points_to_review() {
  unsigned int CP_range = min(_cp_search_range, size()-1);
  while(_points_under_review.size() > 0) {
    Point * this_point = _points_under_review.back();
    _points_under_review.pop_back();  
    if (this_point->review_flag & _remove_heap_entry) {
      assert(!(this_point->review_flag ^ _remove_heap_entry));
      _heap->remove(_ID(this_point));
    } 
    else {
      if (this_point->review_flag & _review_neighbour) {
	this_point->neighbour_dist2 = numeric_limits<double>::max();
	// among all three shifts
	for (unsigned int ishift = 0; ishift < _nshift; ishift++) {
	  circulator other = this_point->circ[ishift];
	  // among points within CP_range
	  for (unsigned i=0; i < CP_range; i++) {
	    ++other;
	    double dist2 = this_point->distance2(*other->point);
	    if (dist2 < this_point->neighbour_dist2) {
	      this_point->neighbour_dist2 = dist2;
	      this_point->neighbour       = other->point;
	    }
	  }
	}
      }
      _heap->update(_ID(this_point), this_point->neighbour_dist2);
    }
    this_point->review_flag = 0; 
  }
}
unsigned int ClosestPair2D::insert(const Coord2D & new_coord) {
  assert(_available_points.size() > 0);
  Point * new_point = _available_points.top();
  _available_points.pop();
  new_point->coord = new_coord;
  _insert_into_search_tree(new_point);
  _deal_with_points_to_review();
  return _ID(new_point);
}
unsigned int ClosestPair2D::replace(unsigned int ID1, unsigned int ID2, 
				    const Coord2D & position) {
  Point * point_to_remove = & (_points[ID1]);
  _remove_from_search_tree(point_to_remove);
  point_to_remove = & (_points[ID2]);
  _remove_from_search_tree(point_to_remove);
  Point * new_point = _available_points.top();
  _available_points.pop();
  new_point->coord = position;
  _insert_into_search_tree(new_point);
  _deal_with_points_to_review();
  return _ID(new_point);
}
void ClosestPair2D::replace_many(
                  const std::vector<unsigned int> & IDs_to_remove,
		  const std::vector<Coord2D> & new_positions,
		  std::vector<unsigned int> & new_IDs) {
  for (unsigned int i = 0; i < IDs_to_remove.size(); i++) {
    _remove_from_search_tree(& (_points[IDs_to_remove[i]]));
  }
  new_IDs.resize(0);
  for (unsigned int i = 0; i < new_positions.size(); i++) {
    Point * new_point = _available_points.top();
    _available_points.pop();
    new_point->coord = new_positions[i];
    _insert_into_search_tree(new_point);
    new_IDs.push_back(_ID(new_point));
  }
  _deal_with_points_to_review();
}
void ClosestPair2D::_insert_into_search_tree(Point * new_point) {
  _set_label(new_point, _review_heap_entry);
  new_point->neighbour_dist2 = numeric_limits<double>::max();
  unsigned int CP_range = min(_cp_search_range, size()-1);
  for (unsigned ishift = 0; ishift < _nshift; ishift++) {
    Shuffle new_shuffle;
    _point2shuffle(*new_point, new_shuffle, _shifts[ishift]);
    circulator new_circ = _trees[ishift]->insert(new_shuffle);
    new_point->circ[ishift] = new_circ;
    circulator right_edge = new_circ; right_edge++;
    circulator left_edge  = new_circ;
    for (unsigned int i = 0; i < CP_range; i++) {left_edge--;}
    do {
      Point * left_point  = left_edge->point;
      Point * right_point = right_edge->point;
      double new_dist2 = left_point->distance2(*new_point);
      if (new_dist2 < left_point->neighbour_dist2) {
	left_point->neighbour_dist2 = new_dist2;
	left_point->neighbour       = new_point;
	_add_label(left_point, _review_heap_entry);
      }
      new_dist2 = new_point->distance2(*right_point);
      if (new_dist2 < new_point->neighbour_dist2) {
	new_point->neighbour_dist2 = new_dist2;
	new_point->neighbour = right_point;
      }
      if (left_point->neighbour == right_point) {
	_add_label(left_point, _review_neighbour);
      }
      right_edge++;
    } while (++left_edge != new_circ);
  }
}
FJCORE_END_NAMESPACE
#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<cassert>
#include<string>
#include<set>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
using namespace std;
std::ostream * ClusterSequence::_fastjet_banner_ostr = &cout;
ClusterSequence::~ClusterSequence () {
  if (_structure_shared_ptr){
    ClusterSequenceStructure* csi = dynamic_cast<ClusterSequenceStructure*>(_structure_shared_ptr.get()); 
    assert(csi != NULL);
    csi->set_associated_cs(NULL);
    if (_deletes_self_when_unused) {
      _structure_shared_ptr.set_count(_structure_shared_ptr.use_count() 
				        + _structure_use_count_after_construction);
    }
  }
}
void ClusterSequence::signal_imminent_self_deletion() const {
  assert(_deletes_self_when_unused);
  _deletes_self_when_unused = false;
}
void ClusterSequence::_initialise_and_run (
				  const JetDefinition & jet_def_in,
				  const bool & writeout_combinations) {
  _decant_options(jet_def_in, writeout_combinations);
  _initialise_and_run_no_decant();
}
void ClusterSequence::_initialise_and_run_no_decant () {
  _fill_initial_history();
  if (n_particles() == 0) return;
  if (_jet_algorithm == plugin_algorithm) {
    _plugin_activated = true;
    _jet_def.plugin()->run_clustering( (*this) );
    _plugin_activated = false;
    _update_structure_use_count();
    return;
  } else if (_jet_algorithm == ee_kt_algorithm ||
	     _jet_algorithm == ee_genkt_algorithm) {
    _strategy = N2Plain;
    if (_jet_algorithm == ee_kt_algorithm) {
      assert(_Rparam > 2.0); 
      _invR2 = 1.0;
    } else {
      if (_Rparam > pi) {
	// choose a value that ensures that back-to-back particles will
	// always recombine 
	//_R2 = 4.0000000000001;
	_R2 = 2 * ( 3.0 + cos(_Rparam) );
      } else {
	_R2    = 2 * ( 1.0 - cos(_Rparam) );
      }
      _invR2 = 1.0/_R2;
    }
    _simple_N2_cluster_EEBriefJet();
    return;
  } else if (_jet_algorithm == undefined_jet_algorithm) {
    throw Error("A ClusterSequence cannot be created with an uninitialised JetDefinition");
  }
  if (_strategy == Best) {
    _strategy = _best_strategy();
#ifdef __FJCORE_DROP_CGAL
    if (_strategy == NlnN) _strategy = N2MHTLazy25;
#endif  // __FJCORE_DROP_CGAL
  } else if (_strategy == BestFJ30) {
    int N = _jets.size();
    if (min(1.0,max(0.1,_Rparam)*3.3)*N <= 30) {
      _strategy = N2Plain;
    } else if (N > 6200/pow(_Rparam,2.0) && _jet_def.jet_algorithm() == cambridge_algorithm) {
      _strategy = NlnNCam;
#ifndef __FJCORE_DROP_CGAL
    } else if ((N > 16000/pow(_Rparam,1.15) && _jet_def.jet_algorithm() != antikt_algorithm)
	       || N > 35000/pow(_Rparam,1.15)) {
      _strategy = NlnN;
#endif  // __FJCORE_DROP_CGAL
    } else if (N <= 450) {
      _strategy = N2Tiled;
    } else {                   
      _strategy = N2MinHeapTiled;
    }
  }
  if (_Rparam >= twopi) {
    if (   _strategy == NlnN
	|| _strategy == NlnN3pi
	|| _strategy == NlnNCam
	|| _strategy == NlnNCam2pi2R
	|| _strategy == NlnNCam4pi) {
#ifdef __FJCORE_DROP_CGAL
      _strategy = N2MinHeapTiled;
#else
      _strategy = NlnN4pi;
#endif    
    }
    if (_jet_def.strategy() != Best && _strategy != _jet_def.strategy()) {
      ostringstream oss;
      oss << "Cluster strategy " << strategy_string(_jet_def.strategy())
	  << " automatically changed to " << strategy_string()
	  << " because the former is not supported for R = " << _Rparam
	  << " >= 2pi";
      _changed_strategy_warning.warn(oss.str());
    }
  }
  if (_strategy == N2Plain) {
    this->_simple_N2_cluster_BriefJet();
  } else if (_strategy == N2Tiled) {
    this->_faster_tiled_N2_cluster();
  } else if (_strategy == N2MinHeapTiled) {
    this->_minheap_faster_tiled_N2_cluster();
  } else if (_strategy == N2MHTLazy9Alt) {
    _plugin_activated = true;
    LazyTiling9Alt tiling(*this);
    tiling.run();
    _plugin_activated = false;
  } else if (_strategy == N2MHTLazy25) {
    _plugin_activated = true;
    LazyTiling25 tiling(*this);
    tiling.run();
    _plugin_activated = false;
  } else if (_strategy == N2MHTLazy9) {
    _plugin_activated = true;
    LazyTiling9 tiling(*this);
    tiling.run();
    _plugin_activated = false;
  } else if (_strategy == N2MHTLazy9AntiKtSeparateGhosts) {
    throw Error("N2MHTLazy9AntiKtSeparateGhosts strategy not supported with FJCORE");
  } else if (_strategy == NlnN) {
    this->_delaunay_cluster();
  } else if (_strategy == NlnNCam) {
    this->_CP2DChan_cluster_2piMultD();
  } else if (_strategy == NlnN3pi || _strategy == NlnN4pi ) {
    this->_delaunay_cluster();
  } else if (_strategy ==  N3Dumb ) {
    this->_really_dumb_cluster();
  } else if (_strategy == N2PoorTiled) {
    this->_tiled_N2_cluster();
  } else if (_strategy == NlnNCam4pi) {
    this->_CP2DChan_cluster();
  } else if (_strategy == NlnNCam2pi2R) {
    this->_CP2DChan_cluster_2pi2R();
  } else {
    ostringstream err;
    err << "Unrecognised value for strategy: "<<_strategy;
    throw Error(err.str());
  }
}
bool ClusterSequence::_first_time = true;
LimitedWarning ClusterSequence::_exclusive_warnings;
string fastjet_version_string() {
  return "FastJet version "+string(fastjet_version)+" [fjcore]";
}
void ClusterSequence::print_banner() {
  if (!_first_time) {return;}
  _first_time = false;
  ostream * ostr = _fastjet_banner_ostr;
  if (!ostr) return;  
  (*ostr) << "#--------------------------------------------------------------------------\n";
  (*ostr) << "#                     FastJet release " << fastjet_version << " [fjcore]" << endl;
  (*ostr) << "#                 M. Cacciari, G.P. Salam and G. Soyez                  \n"; 
  (*ostr) << "#     A software package for jet finding and analysis at colliders      \n";
  (*ostr) << "#                           http://fastjet.fr                           \n"; 
  (*ostr) << "#	                                                                      \n";
  (*ostr) << "# Please cite EPJC72(2012)1896 [arXiv:1111.6097] if you use this package\n";
  (*ostr) << "# for scientific work and optionally PLB641(2006)57 [hep-ph/0512210].   \n";
  (*ostr) << "#                                                                       \n";
  (*ostr) << "# FastJet is provided without warranty under the terms of the GNU GPLv2.\n";
  (*ostr) << "# It uses T. Chan's closest pair algorithm, S. Fortune's Voronoi code";
#ifndef __FJCORE_DROP_CGAL
  (*ostr) << ",\n# CGAL ";
#else
  (*ostr) << "\n# ";
#endif  // __FJCORE_DROP_CGAL
  (*ostr) << "and 3rd party plugin jet algorithms. See COPYING file for details.\n";
  (*ostr) << "#--------------------------------------------------------------------------\n";
  ostr->flush();
}
void ClusterSequence::_decant_options(const JetDefinition & jet_def_in,
                                      const bool & writeout_combinations) {
  _jet_def = jet_def_in;
  _writeout_combinations = writeout_combinations;
  _structure_shared_ptr.reset(new ClusterSequenceStructure(this));
  _decant_options_partial();
}
void ClusterSequence::_decant_options_partial() {
  print_banner();
  _jet_algorithm = _jet_def.jet_algorithm();
  _Rparam = _jet_def.R();  _R2 = _Rparam*_Rparam; _invR2 = 1.0/_R2;
  _strategy = _jet_def.strategy();
  _plugin_activated = false;
  _update_structure_use_count(); // make sure it's correct already here
}
void ClusterSequence::_fill_initial_history () {
  _jets.reserve(_jets.size()*2);
  _history.reserve(_jets.size()*2);
  _Qtot = 0;
  for (int i = 0; i < static_cast<int>(_jets.size()) ; i++) {
    history_element element;
    element.parent1 = InexistentParent;
    element.parent2 = InexistentParent;
    element.child   = Invalid;
    element.jetp_index = i;
    element.dij     = 0.0;
    element.max_dij_so_far = 0.0;
    _history.push_back(element);
    _jet_def.recombiner()->preprocess(_jets[i]);
    _jets[i].set_cluster_hist_index(i);
    _set_structure_shared_ptr(_jets[i]);
    _Qtot += _jets[i].E();
  }
  _initial_n = _jets.size();
  _deletes_self_when_unused = false;
}
string ClusterSequence::strategy_string (Strategy strategy_in)  const {
  string strategy;
  switch(strategy_in) {
  case NlnN:
    strategy = "NlnN"; break;
  case NlnN3pi:
    strategy = "NlnN3pi"; break;
  case NlnN4pi:
    strategy = "NlnN4pi"; break;
  case N2Plain:
    strategy = "N2Plain"; break;
  case N2Tiled:
    strategy = "N2Tiled"; break;
  case N2MinHeapTiled:
    strategy = "N2MinHeapTiled"; break;
  case N2PoorTiled:
    strategy = "N2PoorTiled"; break;
  case N2MHTLazy9:
    strategy = "N2MHTLazy9"; break;
  case N2MHTLazy9Alt:
    strategy = "N2MHTLazy9Alt"; break;
  case N2MHTLazy25:
    strategy = "N2MHTLazy25"; break;
  case N2MHTLazy9AntiKtSeparateGhosts:
    strategy = "N2MHTLazy9AntiKtSeparateGhosts"; break;
  case N3Dumb:
    strategy = "N3Dumb"; break;
  case NlnNCam4pi:
    strategy = "NlnNCam4pi"; break;
  case NlnNCam2pi2R:
    strategy = "NlnNCam2pi2R"; break;
  case NlnNCam:
    strategy = "NlnNCam"; break; // 2piMultD
  case plugin_strategy:
    strategy = "plugin strategy"; break;
  default:
    strategy = "Unrecognized";
  }
  return strategy;
}  
double ClusterSequence::jet_scale_for_algorithm(
				  const PseudoJet & jet) const {
  if (_jet_algorithm == kt_algorithm)             {return jet.kt2();}
  else if (_jet_algorithm == cambridge_algorithm) {return 1.0;}
  else if (_jet_algorithm == antikt_algorithm) {
    double kt2=jet.kt2();
    return kt2 > 1e-300 ? 1.0/kt2 : 1e300;
  } else if (_jet_algorithm == genkt_algorithm) {
    double kt2 = jet.kt2();
    double p   = jet_def().extra_param();
    if (p <= 0 && kt2 < 1e-300) kt2 = 1e-300; // dodgy safety check
    return pow(kt2, p);
  } else if (_jet_algorithm == cambridge_for_passive_algorithm) {
    double kt2 = jet.kt2();
    double lim = _jet_def.extra_param();
    if (kt2 < lim*lim && kt2 != 0.0) {
      return 1.0/kt2;
    } else {return 1.0;}
  } else {throw Error("Unrecognised jet algorithm");}
}
Strategy ClusterSequence::_best_strategy() const {
  int N = _jets.size();
  double bounded_R = max(_Rparam, 0.1);
  if (N <= 30 || N <= 39.0/(bounded_R + 0.6)) {
    return N2Plain;
  } 
  const static _Parabola N_Tiled_to_MHT_lowR             (-45.4947,54.3528,44.6283);
  const static _Parabola L_MHT_to_MHTLazy9_lowR          (0.677807,-1.05006,10.6994);
  const static _Parabola L_MHTLazy9_to_MHTLazy25_akt_lowR(0.169967,-0.512589,12.1572);
  const static _Parabola L_MHTLazy9_to_MHTLazy25_kt_lowR (0.16237,-0.484612,12.3373);
  const static _Parabola L_MHTLazy9_to_MHTLazy25_cam_lowR = L_MHTLazy9_to_MHTLazy25_kt_lowR;
  const static _Parabola L_MHTLazy25_to_NlnN_akt_lowR    (0.0472051,-0.22043,15.9196);
  const static _Parabola L_MHTLazy25_to_NlnN_kt_lowR     (0.118609,-0.326811,14.8287);
  const static _Parabola L_MHTLazy25_to_NlnN_cam_lowR    (0.10119,-0.295748,14.3924);
  const static _Line     L_Tiled_to_MHTLazy9_medR         (-1.31304,7.29621);
  const static _Parabola L_MHTLazy9_to_MHTLazy25_akt_medR = L_MHTLazy9_to_MHTLazy25_akt_lowR;
  const static _Parabola L_MHTLazy9_to_MHTLazy25_kt_medR  = L_MHTLazy9_to_MHTLazy25_kt_lowR;
  const static _Parabola L_MHTLazy9_to_MHTLazy25_cam_medR = L_MHTLazy9_to_MHTLazy25_cam_lowR;
  const static _Parabola L_MHTLazy25_to_NlnN_akt_medR     = L_MHTLazy25_to_NlnN_akt_lowR;
  const static _Parabola L_MHTLazy25_to_NlnN_kt_medR      = L_MHTLazy25_to_NlnN_kt_lowR;
  const static _Parabola L_MHTLazy25_to_NlnN_cam_medR     = L_MHTLazy25_to_NlnN_cam_lowR;
  const static double    N_Plain_to_MHTLazy9_largeR         = 75;
  const static double    N_MHTLazy9_to_MHTLazy25_akt_largeR = 700;
  const static double    N_MHTLazy9_to_MHTLazy25_kt_largeR  = 1000;
  const static double    N_MHTLazy9_to_MHTLazy25_cam_largeR = 1000;
  const static double    N_MHTLazy25_to_NlnN_akt_largeR     = 100000;
  const static double    N_MHTLazy25_to_NlnN_kt_largeR      = 40000;
  const static double    N_MHTLazy25_to_NlnN_cam_largeR     = 15000;
  JetAlgorithm jet_algorithm;
  if (_jet_algorithm == genkt_algorithm) {
    double p   = jet_def().extra_param();
    if (p < 0.0) jet_algorithm = antikt_algorithm;
    else         jet_algorithm =     kt_algorithm;
  } else if (_jet_algorithm == cambridge_for_passive_algorithm) {
    jet_algorithm = kt_algorithm;
  } else {
    jet_algorithm = _jet_algorithm;
  }
  if (bounded_R < 0.65) {
    if          (N    < N_Tiled_to_MHT_lowR(bounded_R))              return N2Tiled;
    double logN = log(double(N));
    if          (logN < L_MHT_to_MHTLazy9_lowR(bounded_R))           return N2MinHeapTiled;
    else {
      if (jet_algorithm == antikt_algorithm){
        if      (logN < L_MHTLazy9_to_MHTLazy25_akt_lowR(bounded_R)) return N2MHTLazy9;
        else if (logN < L_MHTLazy25_to_NlnN_akt_lowR(bounded_R))     return N2MHTLazy25;
        else                                                         return NlnN;
      } else if (jet_algorithm == kt_algorithm){
        if      (logN < L_MHTLazy9_to_MHTLazy25_kt_lowR(bounded_R))  return N2MHTLazy9;
        else if (logN < L_MHTLazy25_to_NlnN_kt_lowR(bounded_R))      return N2MHTLazy25;
        else                                                         return NlnN;
      } else if (jet_algorithm == cambridge_algorithm)  {
        if      (logN < L_MHTLazy9_to_MHTLazy25_cam_lowR(bounded_R)) return N2MHTLazy9;
        else if (logN < L_MHTLazy25_to_NlnN_cam_lowR(bounded_R))     return N2MHTLazy25;
        else                                                         return NlnNCam;
      }
    }
  } else if (bounded_R < 0.5*pi) {
    double logN = log(double(N));
    if      (logN < L_Tiled_to_MHTLazy9_medR(bounded_R))             return N2Tiled;
    else {
      if (jet_algorithm == antikt_algorithm){
        if      (logN < L_MHTLazy9_to_MHTLazy25_akt_medR(bounded_R)) return N2MHTLazy9;
        else if (logN < L_MHTLazy25_to_NlnN_akt_medR(bounded_R))     return N2MHTLazy25;
        else                                                         return NlnN;
      } else if (jet_algorithm == kt_algorithm){
        if      (logN < L_MHTLazy9_to_MHTLazy25_kt_medR(bounded_R))  return N2MHTLazy9;
        else if (logN < L_MHTLazy25_to_NlnN_kt_medR(bounded_R))      return N2MHTLazy25;
        else                                                         return NlnN;
      } else if (jet_algorithm == cambridge_algorithm)  {
        if      (logN < L_MHTLazy9_to_MHTLazy25_cam_medR(bounded_R)) return N2MHTLazy9;
        else if (logN < L_MHTLazy25_to_NlnN_cam_medR(bounded_R))     return N2MHTLazy25;
        else                                                         return NlnNCam;
      }
    }
  } else {
    if      (N    < N_Plain_to_MHTLazy9_largeR)                      return N2Plain;
    else {
      if (jet_algorithm == antikt_algorithm){
        if      (N < N_MHTLazy9_to_MHTLazy25_akt_largeR)             return N2MHTLazy9;
        else if (N < N_MHTLazy25_to_NlnN_akt_largeR)                 return N2MHTLazy25;
        else                                                         return NlnN;
      } else if (jet_algorithm == kt_algorithm){
        if      (N < N_MHTLazy9_to_MHTLazy25_kt_largeR)              return N2MHTLazy9;
        else if (N < N_MHTLazy25_to_NlnN_kt_largeR)                  return N2MHTLazy25;
        else                                                         return NlnN;
      } else if (jet_algorithm == cambridge_algorithm)  {
        if      (N < N_MHTLazy9_to_MHTLazy25_cam_largeR)             return N2MHTLazy9;
        else if (N < N_MHTLazy25_to_NlnN_cam_largeR)                 return N2MHTLazy25;
        else                                                         return NlnNCam;
      }
    }
  }
  assert(0 && "Code should never reach here");
  return N2MHTLazy9;
}
ClusterSequence & ClusterSequence::operator=(const ClusterSequence & cs) {
  if (&cs != this) {
    _deletes_self_when_unused = false;
    transfer_from_sequence(cs);
  }
  return *this;
}
void ClusterSequence::transfer_from_sequence(const ClusterSequence & from_seq,
					     const FunctionOfPseudoJet<PseudoJet> * action_on_jets){
  if (will_delete_self_when_unused()) 
    throw(Error("cannot use CS::transfer_from_sequence after a call to delete_self_when_unused()"));
  _jet_def                 = from_seq._jet_def                ;
  _writeout_combinations   = from_seq._writeout_combinations  ;
  _initial_n               = from_seq._initial_n              ;
  _Rparam                  = from_seq._Rparam                 ;
  _R2                      = from_seq._R2                     ;
  _invR2                   = from_seq._invR2                  ;
  _strategy                = from_seq._strategy               ;
  _jet_algorithm           = from_seq._jet_algorithm          ;
  _plugin_activated        = from_seq._plugin_activated       ;
  if (action_on_jets)
    _jets     = (*action_on_jets)(from_seq._jets);
  else
    _jets     = from_seq._jets;
  _history  = from_seq._history;
  _extras   = from_seq._extras;
  if (_structure_shared_ptr) {
    if (_deletes_self_when_unused) throw Error("transfer_from_sequence cannot be used for a cluster sequence that deletes self when unused");
    ClusterSequenceStructure* csi = dynamic_cast<ClusterSequenceStructure*>(_structure_shared_ptr.get()); 
    assert(csi != NULL);
    csi->set_associated_cs(NULL);
  }
  _structure_shared_ptr.reset(new ClusterSequenceStructure(this));
  _update_structure_use_count();
  for (unsigned int i=0; i<_jets.size(); i++){
    _jets[i].set_cluster_hist_index(from_seq._jets[i].cluster_hist_index());
    _set_structure_shared_ptr(_jets[i]);
  }
}
void ClusterSequence::plugin_record_ij_recombination(
	   int jet_i, int jet_j, double dij, 
	   const PseudoJet & newjet, int & newjet_k) {
  plugin_record_ij_recombination(jet_i, jet_j, dij, newjet_k);
  int tmp_index = _jets[newjet_k].cluster_hist_index();
  _jets[newjet_k] = newjet;
  _jets[newjet_k].set_cluster_hist_index(tmp_index);
  _set_structure_shared_ptr(_jets[newjet_k]);
}
vector<PseudoJet> ClusterSequence::inclusive_jets (const double ptmin) const{
  double dcut = ptmin*ptmin;
  int i = _history.size() - 1; // last jet
  vector<PseudoJet> jets_local;
  if (_jet_algorithm == kt_algorithm) {
    while (i >= 0) {
      if (_history[i].max_dij_so_far < dcut) {break;}
      if (_history[i].parent2 == BeamJet && _history[i].dij >= dcut) {
	// for beam jets
	int parent1 = _history[i].parent1;
	jets_local.push_back(_jets[_history[parent1].jetp_index]);}
      i--;
    }
  } else if (_jet_algorithm == cambridge_algorithm) {
    while (i >= 0) {
      if (_history[i].parent2 != BeamJet) {break;}
      int parent1 = _history[i].parent1;
      const PseudoJet & jet = _jets[_history[parent1].jetp_index];
      if (jet.perp2() >= dcut) {jets_local.push_back(jet);}
      i--;
    }
  } else if (_jet_algorithm == plugin_algorithm 
             || _jet_algorithm == ee_kt_algorithm
             || _jet_algorithm == antikt_algorithm
             || _jet_algorithm == genkt_algorithm
             || _jet_algorithm == ee_genkt_algorithm
             || _jet_algorithm == cambridge_for_passive_algorithm) {
    while (i >= 0) {
      if (_history[i].parent2 == BeamJet) {
	int parent1 = _history[i].parent1;
	const PseudoJet & jet = _jets[_history[parent1].jetp_index];
	if (jet.perp2() >= dcut) {jets_local.push_back(jet);}
      }
      i--;
    }
  } else {throw Error("cs::inclusive_jets(...): Unrecognized jet algorithm");}
  return jets_local;
}
int ClusterSequence::n_exclusive_jets (const double dcut) const {
  int i = _history.size() - 1; // last jet
  while (i >= 0) {
    if (_history[i].max_dij_so_far <= dcut) {break;}
    i--;
  }
  int stop_point = i + 1;
  int njets = 2*_initial_n - stop_point;
  return njets;
}
vector<PseudoJet> ClusterSequence::exclusive_jets (const double dcut) const {
  int njets = n_exclusive_jets(dcut);
  return exclusive_jets(njets);
}
vector<PseudoJet> ClusterSequence::exclusive_jets (const int njets) const {
  if (njets > _initial_n) {
    ostringstream err;
    err << "Requested " << njets << " exclusive jets, but there were only " 
	<< _initial_n << " particles in the event";
    throw Error(err.str());
  }
  return exclusive_jets_up_to(njets);
}
vector<PseudoJet> ClusterSequence::exclusive_jets_up_to (const int njets) const {
  if (( _jet_def.jet_algorithm() != kt_algorithm) &&
      ( _jet_def.jet_algorithm() != cambridge_algorithm) &&
      ( _jet_def.jet_algorithm() != ee_kt_algorithm) &&
      (((_jet_def.jet_algorithm() != genkt_algorithm) && 
	(_jet_def.jet_algorithm() != ee_genkt_algorithm)) || 
       (_jet_def.extra_param() <0)) &&
      ((_jet_def.jet_algorithm() != plugin_algorithm) ||
       (!_jet_def.plugin()->exclusive_sequence_meaningful()))) {
    _exclusive_warnings.warn("dcut and exclusive jets for jet-finders other than kt, C/A or genkt with p>=0 should be interpreted with care.");
  }
  int stop_point = 2*_initial_n - njets;
  if (stop_point < _initial_n) stop_point = _initial_n;
  if (2*_initial_n != static_cast<int>(_history.size())) {
    ostringstream err;
    err << "2*_initial_n != _history.size() -- this endangers internal assumptions!\n";
    throw Error(err.str());
  }
  vector<PseudoJet> jets_local;
  for (unsigned int i = stop_point; i < _history.size(); i++) {
    int parent1 = _history[i].parent1;
    if (parent1 < stop_point) {
      jets_local.push_back(_jets[_history[parent1].jetp_index]);
    }
    int parent2 = _history[i].parent2;
    if (parent2 < stop_point && parent2 > 0) {
      jets_local.push_back(_jets[_history[parent2].jetp_index]);
    }
  }
  if (int(jets_local.size()) != min(_initial_n, njets)) {
    ostringstream err;
    err << "ClusterSequence::exclusive_jets: size of returned vector ("
	 <<jets_local.size()<<") does not coincide with requested number of jets ("
	 <<njets<<")";
    throw Error(err.str());
  }
  return jets_local;
}
double ClusterSequence::exclusive_dmerge (const int njets) const {
  assert(njets >= 0);
  if (njets >= _initial_n) {return 0.0;}
  return _history[2*_initial_n-njets-1].dij;
}
double ClusterSequence::exclusive_dmerge_max (const int njets) const {
  assert(njets >= 0);
  if (njets >= _initial_n) {return 0.0;}
  return _history[2*_initial_n-njets-1].max_dij_so_far;
}
std::vector<PseudoJet> ClusterSequence::exclusive_subjets 
   (const PseudoJet & jet, const double dcut) const {
  set<const history_element*> subhist;
  get_subhist_set(subhist, jet, dcut, 0);
  vector<PseudoJet> subjets;
  subjets.reserve(subhist.size());
  for (set<const history_element*>::iterator elem = subhist.begin(); 
       elem != subhist.end(); elem++) {
    subjets.push_back(_jets[(*elem)->jetp_index]);
  }
  return subjets;
}
int ClusterSequence::n_exclusive_subjets(const PseudoJet & jet, 
                        const double dcut) const {
  set<const history_element*> subhist;
  get_subhist_set(subhist, jet, dcut, 0);
  return subhist.size();
}
std::vector<PseudoJet> ClusterSequence::exclusive_subjets
   (const PseudoJet & jet, int nsub) const {
  vector<PseudoJet> subjets = exclusive_subjets_up_to(jet, nsub);
  if (int(subjets.size()) < nsub) {
    ostringstream err;
    err << "Requested " << nsub << " exclusive subjets, but there were only " 
	<< subjets.size() << " particles in the jet";
    throw Error(err.str());
  }
  return subjets;
}
std::vector<PseudoJet> ClusterSequence::exclusive_subjets_up_to
   (const PseudoJet & jet, int nsub) const {
  set<const history_element*> subhist;
  vector<PseudoJet> subjets;
  if (nsub <  0) throw Error("Requested a negative number of subjets. This is nonsensical.");
  if (nsub == 0) return subjets;
  get_subhist_set(subhist, jet, -1.0, nsub);
  subjets.reserve(subhist.size());
  for (set<const history_element*>::iterator elem = subhist.begin(); 
       elem != subhist.end(); elem++) {
    subjets.push_back(_jets[(*elem)->jetp_index]);
  }
  return subjets;
}
double ClusterSequence::exclusive_subdmerge(const PseudoJet & jet, int nsub) const {
  set<const history_element*> subhist;
  get_subhist_set(subhist, jet, -1.0, nsub);
  set<const history_element*>::iterator highest = subhist.end();
  highest--;
  return (*highest)->dij;
}
double ClusterSequence::exclusive_subdmerge_max(const PseudoJet & jet, int nsub) const {
  set<const history_element*> subhist;
  get_subhist_set(subhist, jet, -1.0, nsub);
  set<const history_element*>::iterator highest = subhist.end();
  highest--;
  return (*highest)->max_dij_so_far;
}
void ClusterSequence::get_subhist_set(set<const history_element*> & subhist,
                                     const  PseudoJet & jet, 
                                     double dcut, int maxjet) const {
  assert(contains(jet));
  subhist.clear();
  subhist.insert(&(_history[jet.cluster_hist_index()]));
  int njet = 1;
  while (true) {
    set<const history_element*>::iterator highest = subhist.end();
    assert (highest != subhist.begin()); 
    highest--;
    const history_element* elem = *highest;
    if (njet == maxjet) break;
    if (elem->parent1 < 0)            break;
    if (elem->max_dij_so_far <= dcut) break;
    subhist.erase(highest);
    subhist.insert(&(_history[elem->parent1]));
    subhist.insert(&(_history[elem->parent2]));
    njet++;
  }
}
bool ClusterSequence::object_in_jet(const PseudoJet & object, 
                                    const PseudoJet & jet) const {
  assert(contains(object) && contains(jet));
  const PseudoJet * this_object = &object;
  const PseudoJet * childp;
  while(true) {
    if (this_object->cluster_hist_index() == jet.cluster_hist_index()) {
      return true;
    } else if (has_child(*this_object, childp)) {
      this_object = childp;
    } else {
      return false;
    }
  }
}
bool ClusterSequence::has_parents(const PseudoJet & jet, PseudoJet & parent1, 
                              PseudoJet & parent2) const {
  const history_element & hist = _history[jet.cluster_hist_index()];
  assert ((hist.parent1 >= 0 && hist.parent2 >= 0) || 
          (hist.parent1 < 0 && hist.parent2 < 0));
  if (hist.parent1 < 0) {
    parent1 = PseudoJet(0.0,0.0,0.0,0.0);
    parent2 = parent1;
    return false;
  } else {
    parent1 = _jets[_history[hist.parent1].jetp_index];
    parent2 = _jets[_history[hist.parent2].jetp_index];
    if (parent1.perp2() < parent2.perp2()) std::swap(parent1,parent2);
    return true;
  }
}
bool ClusterSequence::has_child(const PseudoJet & jet, PseudoJet & child) const {
  const PseudoJet * childp;
  bool res = has_child(jet, childp);
  if (res) {
    child = *childp;
    return true;
  } else {
    child = PseudoJet(0.0,0.0,0.0,0.0);
    return false;
  }
}
bool ClusterSequence::has_child(const PseudoJet & jet, const PseudoJet * & childp) const {
  const history_element & hist = _history[jet.cluster_hist_index()];
  if (hist.child >= 0 && _history[hist.child].jetp_index >= 0) {
    childp = &(_jets[_history[hist.child].jetp_index]);
    return true;
  } else {
    childp = NULL;
    return false;
  }
}
bool ClusterSequence::has_partner(const PseudoJet & jet, 
                              PseudoJet & partner) const {
  const history_element & hist = _history[jet.cluster_hist_index()];
  if (hist.child >= 0 && _history[hist.child].parent2 >= 0) {
    const history_element & child_hist = _history[hist.child];
    if (child_hist.parent1 == jet.cluster_hist_index()) {
      partner = _jets[_history[child_hist.parent2].jetp_index];
    } else {
      partner = _jets[_history[child_hist.parent1].jetp_index];
    }
    return true;
  } else {
    partner = PseudoJet(0.0,0.0,0.0,0.0);
    return false;
  }
}
vector<PseudoJet> ClusterSequence::constituents (const PseudoJet & jet) const {
  vector<PseudoJet> subjets;
  add_constituents(jet, subjets);
  return subjets;
}
void ClusterSequence::print_jets_for_root(const std::vector<PseudoJet> & jets_in, 
                                          ostream & ostr) const {
  for (unsigned i = 0; i < jets_in.size(); i++) {
    ostr << i  << " "
         << jets_in[i].px() << " "
         << jets_in[i].py() << " "
         << jets_in[i].pz() << " "
         << jets_in[i].E() << endl;
    vector<PseudoJet> cst = constituents(jets_in[i]);
    for (unsigned j = 0; j < cst.size() ; j++) {
      ostr << " " << j << " "
           << cst[j].rap() << " "
           << cst[j].phi() << " "
           << cst[j].perp() << endl;
    }
    ostr << "#END" << endl;
  }
}
void ClusterSequence::print_jets_for_root(const std::vector<PseudoJet> & jets_in, 
					  const std::string & filename,
					  const std::string & comment ) const {
  std::ofstream ostr(filename.c_str());
  if (comment != "") ostr << "# " << comment << endl;
  print_jets_for_root(jets_in, ostr);
}
vector<int> ClusterSequence::particle_jet_indices(
                        const vector<PseudoJet> & jets_in) const {
  vector<int> indices(n_particles());
  for (unsigned ipart = 0; ipart < n_particles(); ipart++) 
    indices[ipart] = -1;
  for (unsigned ijet = 0; ijet < jets_in.size(); ijet++) {
    vector<PseudoJet> jet_constituents(constituents(jets_in[ijet]));
    for (unsigned ip = 0; ip < jet_constituents.size(); ip++) {
      unsigned iclust = jet_constituents[ip].cluster_hist_index();
      unsigned ipart = history()[iclust].jetp_index;
      indices[ipart] = ijet;
    }
  }
  return indices;
}
void ClusterSequence::add_constituents (
           const PseudoJet & jet, vector<PseudoJet> & subjet_vector) const {
  int i = jet.cluster_hist_index();
  int parent1 = _history[i].parent1;
  int parent2 = _history[i].parent2;
  if (parent1 == InexistentParent) {
    subjet_vector.push_back(_jets[i]);
    return;
  } 
  add_constituents(_jets[_history[parent1].jetp_index], subjet_vector);
  if (parent2 != BeamJet) {
    add_constituents(_jets[_history[parent2].jetp_index], subjet_vector);
  }
}
void ClusterSequence::_add_step_to_history (
               const int parent1, 
	       const int parent2, const int jetp_index,
	       const double dij) {
  history_element element;
  element.parent1 = parent1;
  element.parent2 = parent2;
  element.jetp_index = jetp_index;
  element.child = Invalid;
  element.dij   = dij;
  element.max_dij_so_far = max(dij,_history[_history.size()-1].max_dij_so_far);
  _history.push_back(element);
  int local_step = _history.size()-1;
  assert(parent1 >= 0);
  if (_history[parent1].child != Invalid){
    throw InternalError("trying to recomine an object that has previsously been recombined");
  }
  _history[parent1].child = local_step;
  if (parent2 >= 0) {
    if (_history[parent2].child != Invalid){
      throw InternalError("trying to recomine an object that has previsously been recombined");
    }
    _history[parent2].child = local_step;
  }
  if (jetp_index != Invalid) {
    assert(jetp_index >= 0);
    _jets[jetp_index].set_cluster_hist_index(local_step);
    _set_structure_shared_ptr(_jets[jetp_index]);
  }
  if (_writeout_combinations) {
    cout << local_step << ": " 
	 << parent1 << " with " << parent2
	 << "; y = "<< dij<<endl;
  }
}
vector<int> ClusterSequence::unique_history_order() const {
  valarray<int> lowest_constituent(_history.size());
  int hist_n = _history.size();
  lowest_constituent = hist_n; // give it a large number
  for (int i = 0; i < hist_n; i++) {
    lowest_constituent[i] = min(lowest_constituent[i],i); 
    if (_history[i].child > 0) lowest_constituent[_history[i].child] 
      = min(lowest_constituent[_history[i].child],lowest_constituent[i]);
  }
  valarray<bool> extracted(_history.size()); extracted = false;
  vector<int> unique_tree;
  unique_tree.reserve(_history.size());
  for (unsigned i = 0; i < n_particles(); i++) {
    if (!extracted[i]) {
      unique_tree.push_back(i);
      extracted[i] = true;
      _extract_tree_children(i, extracted, lowest_constituent, unique_tree);
    }
  }
  return unique_tree;
}
void ClusterSequence::_extract_tree_children(
       int position, 
       valarray<bool> & extracted, 
       const valarray<int> & lowest_constituent,
       vector<int> & unique_tree) const {
  if (!extracted[position]) {
    _extract_tree_parents(position,extracted,lowest_constituent,unique_tree);
  } 
  int child = _history[position].child;
  if (child  >= 0) _extract_tree_children(child,extracted,lowest_constituent,unique_tree);
}
vector<PseudoJet> ClusterSequence::unclustered_particles() const {
  vector<PseudoJet> unclustered;
  for (unsigned i = 0; i < n_particles() ; i++) {
    if (_history[i].child == Invalid) 
      unclustered.push_back(_jets[_history[i].jetp_index]);
  }
  return unclustered;
}
vector<PseudoJet> ClusterSequence::childless_pseudojets() const {
  vector<PseudoJet> unclustered;
  for (unsigned i = 0; i < _history.size() ; i++) {
    if ((_history[i].child == Invalid) && (_history[i].parent2 != BeamJet))
      unclustered.push_back(_jets[_history[i].jetp_index]);
  }
  return unclustered;
}
bool ClusterSequence::contains(const PseudoJet & jet) const {
  return jet.cluster_hist_index() >= 0 
    &&   jet.cluster_hist_index() < int(_history.size())
    &&   jet.has_valid_cluster_sequence()
    &&   jet.associated_cluster_sequence() == this;
}
void ClusterSequence::_extract_tree_parents(
       int position, 
       valarray<bool> & extracted, 
       const valarray<int> & lowest_constituent,
       vector<int> & unique_tree) const {
  if (!extracted[position]) {
    int parent1 = _history[position].parent1;
    int parent2 = _history[position].parent2;
    if (parent1 >= 0 && parent2 >= 0) {
      if (lowest_constituent[parent1] > lowest_constituent[parent2]) 
	std::swap(parent1, parent2);
    }
    if (parent1 >= 0 && !extracted[parent1]) 
      _extract_tree_parents(parent1,extracted,lowest_constituent,unique_tree);
    if (parent2 >= 0 && !extracted[parent2]) 
      _extract_tree_parents(parent2,extracted,lowest_constituent,unique_tree);
    unique_tree.push_back(position);
    extracted[position] = true;
  }
}
void ClusterSequence::_do_ij_recombination_step(
                               const int jet_i, const int jet_j, 
			       const double dij, 
			       int & newjet_k) {
  PseudoJet newjet(false); 
  _jet_def.recombiner()->recombine(_jets[jet_i], _jets[jet_j], newjet);
  _jets.push_back(newjet);
  newjet_k = _jets.size()-1;
  int newstep_k = _history.size();
  _jets[newjet_k].set_cluster_hist_index(newstep_k);
  int hist_i = _jets[jet_i].cluster_hist_index();
  int hist_j = _jets[jet_j].cluster_hist_index();
  _add_step_to_history(min(hist_i, hist_j), max(hist_i,hist_j),
		       newjet_k, dij);
}
void ClusterSequence::_do_iB_recombination_step(
				  const int jet_i, const double diB) {
  _add_step_to_history(_jets[jet_i].cluster_hist_index(),BeamJet,
		       Invalid, diB);
}
LimitedWarning ClusterSequence::_changed_strategy_warning;
void ClusterSequence::_set_structure_shared_ptr(PseudoJet & j) {
  j.set_structure_shared_ptr(_structure_shared_ptr);
  _update_structure_use_count();
}
void ClusterSequence::_update_structure_use_count() {
  _structure_use_count_after_construction = _structure_shared_ptr.use_count();
}
void ClusterSequence::delete_self_when_unused() {
  int new_count = _structure_shared_ptr.use_count() - _structure_use_count_after_construction;
  if (new_count <= 0) {
    throw Error("delete_self_when_unused may only be called if at least one object outside the CS (e.g. a jet) is already associated with the CS");
  }
  _structure_shared_ptr.set_count(new_count);
  _deletes_self_when_unused = true;
}
FJCORE_END_NAMESPACE
#include<limits>
#include<vector>
#include<cmath>
#include<iostream>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
using namespace std;
namespace Private {
  class MirrorInfo{
  public:
    int orig, mirror;
    MirrorInfo(int a, int b) : orig(a), mirror(b) {}
    MirrorInfo() : orig(0), mirror(0) {} // set dummy values to keep static code checkers happy
  };
  bool make_mirror(Coord2D & point, double Dlim) {
    if (point.y < Dlim)       {point.y += twopi; return true;}
    if (twopi-point.y < Dlim) {point.y -= twopi; return true;}
    return false;
  }
}
using namespace Private;
void ClusterSequence::_CP2DChan_limited_cluster (double Dlim) {
  unsigned int n = _initial_n;
  vector<MirrorInfo>   coordIDs(2*n); // coord IDs of a given jetID
  vector<int>          jetIDs(2*n);   // jet ID for a given coord ID
  vector<Coord2D>      coords(2*n);   // our coordinates (and copies)
  double Dlim4mirror = min(Dlim,pi);
  double minrap = numeric_limits<double>::max();
  double maxrap = -minrap;
  int coord_index = -1;
  int n_active = 0;
  for (unsigned jet_i = 0; jet_i < _jets.size(); jet_i++) {
    if (_history[_jets[jet_i].cluster_hist_index()].child != Invalid ||
	(_jets[jet_i].E() == abs(_jets[jet_i].pz()) && 
	 _jets[jet_i].perp2() == 0.0)
	) {continue;}
    n_active++;
    coordIDs[jet_i].orig = ++coord_index;
    coords[coord_index]  = Coord2D(_jets[jet_i].rap(), _jets[jet_i].phi_02pi());
    jetIDs[coord_index]  = jet_i;
    minrap = min(coords[coord_index].x,minrap);
    maxrap = max(coords[coord_index].x,maxrap);
    Coord2D mirror_point(coords[coord_index]);
    if (make_mirror(mirror_point, Dlim4mirror)) {
      coordIDs[jet_i].mirror = ++coord_index;
      coords[coord_index] = mirror_point;
      jetIDs[coord_index] = jet_i;
    } else {
      coordIDs[jet_i].mirror = Invalid;
    }
  }
  coords.resize(coord_index+1);
  Coord2D left_edge(minrap-1.0, -3.15); // a security margin below  -pi
  Coord2D right_edge(maxrap+1.0, 9.45); // a security margin above 3*pi
  ClosestPair2D cp(coords, left_edge, right_edge);
  vector<Coord2D> new_points(2);
  vector<unsigned int> cIDs_to_remove(4);
  vector<unsigned int> new_cIDs(2);
  do {
    unsigned int cID1, cID2;
    double distance2;
    cp.closest_pair(cID1,cID2,distance2);
    if (distance2 > Dlim*Dlim) {break;}
    distance2 *= _invR2;
    int jet_i = jetIDs[cID1];
    int jet_j = jetIDs[cID2];
    assert (jet_i != jet_j); // to catch issue of recombining with mirror point
    int newjet_k;
    _do_ij_recombination_step(jet_i, jet_j, distance2, newjet_k);
    if (--n_active == 1) {break;}
    cIDs_to_remove.resize(0);
    cIDs_to_remove.push_back(coordIDs[jet_i].orig);
    cIDs_to_remove.push_back(coordIDs[jet_j].orig);
    if (coordIDs[jet_i].mirror != Invalid) 
      cIDs_to_remove.push_back(coordIDs[jet_i].mirror);
    if (coordIDs[jet_j].mirror != Invalid) 
      cIDs_to_remove.push_back(coordIDs[jet_j].mirror);
    Coord2D new_point(_jets[newjet_k].rap(),_jets[newjet_k].phi_02pi());
    new_points.resize(0);
    new_points.push_back(new_point);
    if (make_mirror(new_point, Dlim4mirror)) new_points.push_back(new_point);  //< same warning as before concerning the mirroring
    cp.replace_many(cIDs_to_remove, new_points, new_cIDs);
    coordIDs[newjet_k].orig = new_cIDs[0];
    jetIDs[new_cIDs[0]]       = newjet_k;
    if (new_cIDs.size() == 2) {
      coordIDs[newjet_k].mirror = new_cIDs[1];
      jetIDs[new_cIDs[1]]         = newjet_k;
    } else {coordIDs[newjet_k].mirror = Invalid;}
  } while(true);
}
void ClusterSequence::_CP2DChan_cluster_2pi2R () {
  if (_jet_algorithm != cambridge_algorithm) throw Error("CP2DChan clustering method called for a jet-finder that is not the cambridge algorithm");
  _CP2DChan_limited_cluster(_Rparam);
  _do_Cambridge_inclusive_jets();
}
void ClusterSequence::_CP2DChan_cluster_2piMultD () {
  if (_Rparam >= 0.39) {
    _CP2DChan_limited_cluster(min(_Rparam/2,0.3));
  }
  _CP2DChan_cluster_2pi2R ();
}
void ClusterSequence::_CP2DChan_cluster () {
  if (_jet_algorithm != cambridge_algorithm) throw Error("_CP2DChan_cluster called for a jet-finder that is not the cambridge algorithm");
  unsigned int n = _jets.size();
  vector<MirrorInfo>   coordIDs(2*n);  // link from original to mirror indices
  vector<int>          jetIDs(2*n);     // link from mirror to original indices
  vector<Coord2D>      coords(2*n);   // our coordinates (and copies)
  double minrap = numeric_limits<double>::max();
  double maxrap = -minrap;
  int coord_index = 0;
  for (unsigned i = 0; i < n; i++) {
    if (_jets[i].E() == abs(_jets[i].pz()) && _jets[i].perp2() == 0.0) {
      coordIDs[i] = MirrorInfo(BeamJet,BeamJet);
    } else {
      coordIDs[i].orig   = coord_index;
      coordIDs[i].mirror = coord_index+1;
      coords[coord_index]   = Coord2D(_jets[i].rap(), _jets[i].phi_02pi());
      coords[coord_index+1] = Coord2D(_jets[i].rap(), _jets[i].phi_02pi()+twopi);
      jetIDs[coord_index]   = i;
      jetIDs[coord_index+1] = i;
      minrap = min(coords[coord_index].x,minrap);
      maxrap = max(coords[coord_index].x,maxrap);
      coord_index += 2;
    }
  }
  for (unsigned i = n; i < 2*n; i++) {coordIDs[i].orig = Invalid;}
  coords.resize(coord_index);
  Coord2D left_edge(minrap-1.0, 0.0);
  Coord2D right_edge(maxrap+1.0, 2*twopi);
  ClosestPair2D cp(coords, left_edge, right_edge);
  vector<Coord2D> new_points(2);
  vector<unsigned int> cIDs_to_remove(4);
  vector<unsigned int> new_cIDs(2);
  do {
    unsigned int cID1, cID2;
    double distance2;
    cp.closest_pair(cID1,cID2,distance2);
    distance2 *= _invR2;
    if (distance2 > 1.0) {break;}
    int jet_i = jetIDs[cID1];
    int jet_j = jetIDs[cID2];
    assert (jet_i != jet_j); // to catch issue of recombining with mirror point
    int newjet_k;
    _do_ij_recombination_step(jet_i, jet_j, distance2, newjet_k);
    cIDs_to_remove[0] = coordIDs[jet_i].orig;
    cIDs_to_remove[1] = coordIDs[jet_i].mirror;
    cIDs_to_remove[2] = coordIDs[jet_j].orig;
    cIDs_to_remove[3] = coordIDs[jet_j].mirror;
    new_points[0] = Coord2D(_jets[newjet_k].rap(),_jets[newjet_k].phi_02pi());
    new_points[1] = Coord2D(_jets[newjet_k].rap(),_jets[newjet_k].phi_02pi()+twopi);
    new_cIDs[0] = cp.replace(cIDs_to_remove[0], cIDs_to_remove[2], new_points[0]);
    new_cIDs[1] = cp.replace(cIDs_to_remove[1], cIDs_to_remove[3], new_points[1]);
    coordIDs[jet_i].orig = Invalid;
    coordIDs[jet_j].orig = Invalid;
    coordIDs[newjet_k] = MirrorInfo(new_cIDs[0], new_cIDs[1]);
    jetIDs[new_cIDs[0]] = newjet_k;
    jetIDs[new_cIDs[1]] = newjet_k;
    n--;
    if (n == 1) {break;}
  } while(true);
  _do_Cambridge_inclusive_jets();
}
void ClusterSequence::_do_Cambridge_inclusive_jets () {
  unsigned int n = _history.size();
  for (unsigned int hist_i = 0; hist_i < n; hist_i++) {
    if (_history[hist_i].child == Invalid) {
      _do_iB_recombination_step(_history[hist_i].jetp_index, 1.0);
    }
  }
}
FJCORE_END_NAMESPACE
#include<iostream>
#include<sstream>
#include<cmath>
#include <cstdlib>
#include<cassert>
#include<memory>
#ifndef __FJCORE_DROP_CGAL // in case we do not have the code for CGAL
#include "fastjet/internal/Dnn4piCylinder.hh"
#include "fastjet/internal/Dnn3piCylinder.hh"
#include "fastjet/internal/Dnn2piCylinder.hh"
#endif //  __FJCORE_DROP_CGAL 
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
using namespace std;
void ClusterSequence::_delaunay_cluster () {
  int n = _jets.size();
  vector<EtaPhi> points(n); // recall EtaPhi is just a typedef'd pair<double>
  for (int i = 0; i < n; i++) {
    points[i] = EtaPhi(_jets[i].rap(),_jets[i].phi_02pi());
    points[i].sanitize(); // make sure things are in the right range
  }
  SharedPtr<DynamicNearestNeighbours> DNN;
  const bool verbose = false;
#ifndef __FJCORE_DROP_CGAL // strategy = NlnN* are not supported if we drop CGAL...
  bool ignore_nearest_is_mirror = (_Rparam < twopi);
  if (_strategy == NlnN4pi) {
    DNN.reset(new Dnn4piCylinder(points,verbose));
  } else if (_strategy == NlnN3pi) {
    DNN.reset(new Dnn3piCylinder(points,ignore_nearest_is_mirror,verbose));
  } else if (_strategy == NlnN) {
    DNN.reset(new Dnn2piCylinder(points,ignore_nearest_is_mirror,verbose));
  } else 
#else
  if (_strategy == NlnN4pi || _strategy == NlnN3pi || _strategy == NlnN) {
    ostringstream err;
    err << "ERROR: Requested strategy "<<strategy_string()<<" but it is not"<<endl;
    err << "       supported because FastJet was compiled without CGAL"<<endl;
    throw Error(err.str());
  } else
#endif // __FJCORE_DROP_CGAL
  {
    assert(false);
  }
  DistMap DijMap;
  for (int ii = 0; ii < n; ii++) {
    _add_ktdistance_to_map(ii, DijMap, DNN.get());
  }
  for (int i=0;i<n;i++) {
    TwoVertices SmallestDijPair;
    int jet_i, jet_j;
    double SmallestDij;
    bool Valid2;
    bool recombine_with_beam;
    do { 
      SmallestDij = DijMap.begin()->first;
      SmallestDijPair = DijMap.begin()->second;
      jet_i = SmallestDijPair.first;
      jet_j = SmallestDijPair.second;
      if (verbose) cout << "CS_Delaunay found recombination candidate: " << jet_i << " " << jet_j << " " << SmallestDij << endl; // GPS debugging
      DijMap.erase(DijMap.begin());
      recombine_with_beam = (jet_j == BeamJet);
      if (!recombine_with_beam) {Valid2 = DNN->Valid(jet_j);} 
      else {Valid2 = true;}
      if (verbose) cout << "CS_Delaunay validities i & j: " << DNN->Valid(jet_i) << " " << Valid2 << endl;
    } while ( !DNN->Valid(jet_i) || !Valid2);
    if (! recombine_with_beam) {
      int nn; // will be index of new jet
      if (verbose) cout << "CS_Delaunay call _do_ij_recomb: " << jet_i << " " << jet_j << " " << SmallestDij << endl; // GPS debug
      _do_ij_recombination_step(jet_i, jet_j, SmallestDij, nn);
      EtaPhi newpoint(_jets[nn].rap(), _jets[nn].phi_02pi());
      newpoint.sanitize(); // make sure it is in correct range
      points.push_back(newpoint);
    } else {
      if (verbose) cout << "CS_Delaunay call _do_iB_recomb: " << jet_i << " " << SmallestDij << endl; // GPS debug
      _do_iB_recombination_step(jet_i, SmallestDij);
    }
    if (i == n-1) {break;}
    vector<int> updated_neighbours;
    if (! recombine_with_beam) {
      int point3;
      DNN->RemoveCombinedAddCombination(jet_i, jet_j, 
				       points[points.size()-1], point3,
				       updated_neighbours);
      if (static_cast<unsigned int> (point3) != points.size()-1) {
	throw Error("INTERNAL ERROR: point3 != points.size()-1");}
    } else {
      DNN->RemovePoint(jet_i, updated_neighbours);
    }
    vector<int>::iterator it = updated_neighbours.begin();
    for (; it != updated_neighbours.end(); ++it) {
      int ii = *it;
      _add_ktdistance_to_map(ii, DijMap, DNN.get());
    }
  } // end clustering loop 
}
void ClusterSequence::_add_ktdistance_to_map(
                          const int ii, 
			  DistMap & DijMap,
			  const DynamicNearestNeighbours * DNN) {
  double yiB = jet_scale_for_algorithm(_jets[ii]);
  if (yiB == 0.0) {
    DijMap.insert(DijEntry(yiB,  TwoVertices(ii,-1)));
  } else {
    double DeltaR2 = DNN->NearestNeighbourDistance(ii) * _invR2;
    if (DeltaR2 > 1.0) {
      DijMap.insert(DijEntry(yiB,  TwoVertices(ii,-1)));
    } else {
      double kt2i = jet_scale_for_algorithm(_jets[ii]);
      int jj = DNN->NearestNeighbourIndex(ii);
      if (kt2i <= jet_scale_for_algorithm(_jets[jj])) {
	double dij = DeltaR2 * kt2i;
	DijMap.insert(DijEntry(dij, TwoVertices(ii,jj)));
      }
    }
  }
}
FJCORE_END_NAMESPACE
#include<iostream>
#include<cmath>
#include <cstdlib>
#include<cassert>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
using namespace std;
void ClusterSequence::_really_dumb_cluster () {
  vector<PseudoJet *> jetsp(_jets.size());
  vector<int>         indices(_jets.size());
  for (size_t i = 0; i<_jets.size(); i++) {
    jetsp[i] = & _jets[i];
    indices[i] = i;
  }
  for (int n = jetsp.size(); n > 0; n--) {
    int ii, jj;
    double ymin = jet_scale_for_algorithm(*(jetsp[0]));
    ii = 0; jj = -2;
    for (int i = 0; i < n; i++) {
      double yiB = jet_scale_for_algorithm(*(jetsp[i]));
      if (yiB < ymin) {
	ymin = yiB; ii = i; jj = -2;}
    }
    for (int i = 0; i < n-1; i++) {
      for (int j = i+1; j < n; j++) {
	//double y = jetsp[i]->kt_distance(*jetsp[j])*_invR2;
	double y = min(jet_scale_for_algorithm(*(jetsp[i])), 
		       jet_scale_for_algorithm(*(jetsp[j])))
	            * jetsp[i]->plain_distance(*jetsp[j])*_invR2;
	if (y < ymin) {ymin = y; ii = i; jj = j;}
      }
    }
    int newn = 2*jetsp.size() - n;
    if (jj >= 0) {
      int nn; // new jet index
      _do_ij_recombination_step(jetsp[ii]-&_jets[0], 
				jetsp[jj]-&_jets[0], ymin, nn);
      jetsp[ii] = &_jets[nn];
      jetsp[jj] = jetsp[n-1];
      indices[ii] = newn;
      indices[jj] = indices[n-1];
    } else {
      _do_iB_recombination_step(jetsp[ii]-&_jets[0], ymin);
      jetsp[ii] = jetsp[n-1];
      indices[ii] = indices[n-1];
    }
  }
}
FJCORE_END_NAMESPACE
#include<iostream>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
using namespace std;
template<> inline void ClusterSequence::_bj_set_jetinfo(
                           EEBriefJet * const jetA, const int _jets_index) const {
  double E = _jets[_jets_index].E();
  double scale = E*E; // the default energy scale for the kt alg
  double p  = jet_def().extra_param(); // in case we're ee_genkt
  switch (_jet_algorithm) {
  case ee_kt_algorithm:
    assert(_Rparam > 2.0); // force this to be true! [not best place, but works]
    break; 
  case ee_genkt_algorithm:
    if (p <= 0 && scale < 1e-300) scale = 1e-300; // same dodgy safety as genkt
    scale = pow(scale,p);
    break;
  default:
    throw Error("Unrecognised jet algorithm");
  }
  jetA->kt2  = scale; // "kt2" might one day be renamed as "scale" or some such
  double norm = _jets[_jets_index].modp2();
  if (norm > 0) {
    norm = 1.0/sqrt(norm);
    jetA->nx = norm * _jets[_jets_index].px();
    jetA->ny = norm * _jets[_jets_index].py();
    jetA->nz = norm * _jets[_jets_index].pz();
  } else {
    jetA->nx = 0.0;
    jetA->ny = 0.0;
    jetA->nz = 1.0;
  }
  jetA->_jets_index = _jets_index;
  jetA->NN_dist = _R2;
  jetA->NN      = NULL;
}
template<> double ClusterSequence::_bj_dist(
                const EEBriefJet * const jeta, 
                const EEBriefJet * const jetb) const {
  double dist = 1.0 
    - jeta->nx*jetb->nx
    - jeta->ny*jetb->ny
    - jeta->nz*jetb->nz;
  dist *= 2; // distance is _2_*min(Ei^2,Ej^2)*(1-cos theta)
  return dist;
}
void ClusterSequence::_simple_N2_cluster_BriefJet() {  
  _simple_N2_cluster<BriefJet>();
}
void ClusterSequence::_simple_N2_cluster_EEBriefJet() {  
  _simple_N2_cluster<EEBriefJet>();
}
FJCORE_END_NAMESPACE
#include <iostream>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
using namespace std;
ClusterSequenceStructure::~ClusterSequenceStructure(){
  if (_associated_cs != NULL 
      && _associated_cs->will_delete_self_when_unused()) {
    _associated_cs->signal_imminent_self_deletion();
    delete _associated_cs;
  }
}
bool ClusterSequenceStructure::has_valid_cluster_sequence() const{
  return (_associated_cs != NULL);
}
const ClusterSequence* ClusterSequenceStructure::associated_cluster_sequence() const{
  return _associated_cs;
}
const ClusterSequence * ClusterSequenceStructure::validated_cs() const {
  if (!_associated_cs) 
    throw Error("you requested information about the internal structure of a jet, but its associated ClusterSequence has gone out of scope.");
  return _associated_cs;
}
bool ClusterSequenceStructure::has_partner(const PseudoJet &reference, PseudoJet &partner) const{
  return validated_cs()->has_partner(reference, partner);
}
bool ClusterSequenceStructure::has_child(const PseudoJet &reference, PseudoJet &child) const{
  return validated_cs()->has_child(reference, child);
}
bool ClusterSequenceStructure::has_parents(const PseudoJet &reference, PseudoJet &parent1, PseudoJet &parent2) const{
  return validated_cs()->has_parents(reference, parent1, parent2);
}
bool ClusterSequenceStructure::object_in_jet(const PseudoJet &reference, const PseudoJet &jet) const{
  if ((!has_associated_cluster_sequence()) || (!jet.has_associated_cluster_sequence()))
    throw Error("you requested information about the internal structure of a jet, but it is not associated with a ClusterSequence or its associated ClusterSequence has gone out of scope."); 
  if (reference.associated_cluster_sequence() != jet.associated_cluster_sequence()) return false;
  return validated_cs()->object_in_jet(reference, jet);
}
bool ClusterSequenceStructure::has_constituents() const{
  if (!has_associated_cluster_sequence())
    throw Error("you requested information about the internal structure of a jet, but it is not associated with a ClusterSequence or its associated ClusterSequence has gone out of scope."); 
  return true;
}
vector<PseudoJet> ClusterSequenceStructure::constituents(const PseudoJet &reference) const{
  return validated_cs()->constituents(reference);
}
bool ClusterSequenceStructure::has_exclusive_subjets() const{
  if (!has_associated_cluster_sequence())
    throw Error("you requested information about the internal structure of a jet, but it is not associated with a ClusterSequence or its associated ClusterSequence has gone out of scope."); 
  return true;
}
std::vector<PseudoJet> ClusterSequenceStructure::exclusive_subjets (const PseudoJet &reference, const double & dcut) const {
  return validated_cs()->exclusive_subjets(reference, dcut);
}
int ClusterSequenceStructure::n_exclusive_subjets(const PseudoJet &reference, const double & dcut) const {
  return validated_cs()->n_exclusive_subjets(reference, dcut);
}
std::vector<PseudoJet> ClusterSequenceStructure::exclusive_subjets_up_to (const PseudoJet &reference, int nsub) const {
  return validated_cs()->exclusive_subjets_up_to(reference, nsub);
}
double ClusterSequenceStructure::exclusive_subdmerge(const PseudoJet &reference, int nsub) const {
  return validated_cs()->exclusive_subdmerge(reference, nsub);
}
double ClusterSequenceStructure::exclusive_subdmerge_max(const PseudoJet &reference, int nsub) const {
  return validated_cs()->exclusive_subdmerge_max(reference, nsub);
}
bool ClusterSequenceStructure::has_pieces(const PseudoJet &reference) const{
  PseudoJet dummy1, dummy2;
  return has_parents(reference, dummy1, dummy2);
}
vector<PseudoJet> ClusterSequenceStructure::pieces(const PseudoJet &reference) const{
  PseudoJet j1, j2;
  vector<PseudoJet> res;
  if (has_parents(reference, j1, j2)){
    res.push_back(j1);
    res.push_back(j2);
  }
  return res;
}
FJCORE_END_NAMESPACE
#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
using namespace std;
void ClusterSequence::_bj_remove_from_tiles(TiledJet * const jet) {
  Tile * tile = & _tiles[jet->tile_index];
  if (jet->previous == NULL) {
    tile->head = jet->next;
  } else {
    jet->previous->next = jet->next;
  }
  if (jet->next != NULL) {
    jet->next->previous = jet->previous;
  }
}
void ClusterSequence::_initialise_tiles() {
  double default_size = max(0.1,_Rparam);
  _tile_size_eta = default_size;
  _n_tiles_phi   = max(3,int(floor(twopi/default_size)));
  _tile_size_phi = twopi / _n_tiles_phi; // >= _Rparam and fits in 2pi
  TilingExtent tiling_analysis(*this);
  _tiles_eta_min = tiling_analysis.minrap();
  _tiles_eta_max = tiling_analysis.maxrap();
  _tiles_ieta_min = int(floor(_tiles_eta_min/_tile_size_eta));
  _tiles_ieta_max = int(floor( _tiles_eta_max/_tile_size_eta));
  _tiles_eta_min = _tiles_ieta_min * _tile_size_eta;
  _tiles_eta_max = _tiles_ieta_max * _tile_size_eta;
  _tiles.resize((_tiles_ieta_max-_tiles_ieta_min+1)*_n_tiles_phi);
  for (int ieta = _tiles_ieta_min; ieta <= _tiles_ieta_max; ieta++) {
    for (int iphi = 0; iphi < _n_tiles_phi; iphi++) {
      Tile * tile = & _tiles[_tile_index(ieta,iphi)];
      tile->head = NULL; // first element of tiles points to itself
      tile->begin_tiles[0] =  tile;
      Tile ** pptile = & (tile->begin_tiles[0]);
      pptile++;
      tile->surrounding_tiles = pptile;
      if (ieta > _tiles_ieta_min) {
	// with the itile subroutine, we can safely run tiles from
	// idphi=-1 to idphi=+1, because it takes care of
	// negative and positive boundaries
	for (int idphi = -1; idphi <=+1; idphi++) {
	  *pptile = & _tiles[_tile_index(ieta-1,iphi+idphi)];
	  pptile++;
	}	
      }
      *pptile = & _tiles[_tile_index(ieta,iphi-1)];
      pptile++;
      tile->RH_tiles = pptile;
      *pptile = & _tiles[_tile_index(ieta,iphi+1)];
      pptile++;
      if (ieta < _tiles_ieta_max) {
	for (int idphi = -1; idphi <= +1; idphi++) {
	  *pptile = & _tiles[_tile_index(ieta+1,iphi+idphi)];
	  pptile++;
	}	
      }
      tile->end_tiles = pptile;
      tile->tagged = false;
    }
  }
}
int ClusterSequence::_tile_index(const double eta, const double phi) const {
  int ieta, iphi;
  if      (eta <= _tiles_eta_min) {ieta = 0;}
  else if (eta >= _tiles_eta_max) {ieta = _tiles_ieta_max-_tiles_ieta_min;}
  else {
    ieta = int(((eta - _tiles_eta_min) / _tile_size_eta));
    if (ieta > _tiles_ieta_max-_tiles_ieta_min) {
      ieta = _tiles_ieta_max-_tiles_ieta_min;} 
  }
  iphi = int((phi+twopi)/_tile_size_phi) % _n_tiles_phi;
  return (iphi + ieta * _n_tiles_phi);
}
inline void ClusterSequence::_tj_set_jetinfo( TiledJet * const jet,
					      const int _jets_index) {
  _bj_set_jetinfo<>(jet, _jets_index);
  jet->tile_index = _tile_index(jet->eta, jet->phi);
  Tile * tile = &_tiles[jet->tile_index];
  jet->previous   = NULL;
  jet->next       = tile->head;
  if (jet->next != NULL) {jet->next->previous = jet;}
  tile->head      = jet;
}
void ClusterSequence::_print_tiles(TiledJet * briefjets ) const {
  for (vector<Tile>::const_iterator tile = _tiles.begin(); 
       tile < _tiles.end(); tile++) {
    cout << "Tile " << tile - _tiles.begin()<<" = ";
    vector<int> list;
    for (TiledJet * jetI = tile->head; jetI != NULL; jetI = jetI->next) {
      list.push_back(jetI-briefjets);
    }
    sort(list.begin(),list.end());
    for (unsigned int i = 0; i < list.size(); i++) {cout <<" "<<list[i];}
    cout <<"\n";
  }
}
void ClusterSequence::_add_neighbours_to_tile_union(const int tile_index, 
	       vector<int> & tile_union, int & n_near_tiles) const {
  for (Tile * const * near_tile = _tiles[tile_index].begin_tiles; 
       near_tile != _tiles[tile_index].end_tiles; near_tile++){
    tile_union[n_near_tiles] = *near_tile - & _tiles[0];
    n_near_tiles++;
  }
}
inline void ClusterSequence::_add_untagged_neighbours_to_tile_union(
               const int tile_index, 
	       vector<int> & tile_union, int & n_near_tiles)  {
  for (Tile ** near_tile = _tiles[tile_index].begin_tiles; 
       near_tile != _tiles[tile_index].end_tiles; near_tile++){
    if (! (*near_tile)->tagged) {
      (*near_tile)->tagged = true;
      tile_union[n_near_tiles] = *near_tile - & _tiles[0];
      n_near_tiles++;
    }
  }
}
void ClusterSequence::_tiled_N2_cluster() {
  _initialise_tiles();
  int n = _jets.size();
  TiledJet * briefjets = new TiledJet[n];
  TiledJet * jetA = briefjets, * jetB;
  TiledJet oldB;
  oldB.tile_index=0; // prevents a gcc warning
  vector<int> tile_union(3*n_tile_neighbours);
  for (int i = 0; i< n; i++) {
    _tj_set_jetinfo(jetA, i);
    jetA++; // move on to next entry of briefjets
  }
  TiledJet * tail = jetA; // a semaphore for the end of briefjets
  TiledJet * head = briefjets; // a nicer way of naming start
  vector<Tile>::const_iterator tile;
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      for (jetB = tile->head; jetB != jetA; jetB = jetB->next) {
	double dist = _bj_dist(jetA,jetB);
	if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
      }
    }
    for (Tile ** RTile = tile->RH_tiles; RTile != tile->end_tiles; RTile++) {
      for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
	for (jetB = (*RTile)->head; jetB != NULL; jetB = jetB->next) {
	  double dist = _bj_dist(jetA,jetB);
	  if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	  if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
	}
      }
    }
  }
  double * diJ = new double[n];
  jetA = head;
  for (int i = 0; i < n; i++) {
    diJ[i] = _bj_diJ(jetA);
    jetA++; // have jetA follow i
  }
  int history_location = n-1;
  while (tail != head) {
    double diJ_min = diJ[0];
    int diJ_min_jet = 0;
    for (int i = 1; i < n; i++) {
      if (diJ[i] < diJ_min) {diJ_min_jet = i; diJ_min  = diJ[i];}
    }
    history_location++;
    jetA = & briefjets[diJ_min_jet];
    jetB = jetA->NN;
    diJ_min *= _invR2; 
    if (jetB != NULL) {
      if (jetA < jetB) {std::swap(jetA,jetB);}
      int nn; // new jet index
      _do_ij_recombination_step(jetA->_jets_index, jetB->_jets_index, diJ_min, nn);
      _bj_remove_from_tiles(jetA);
      oldB = * jetB;  // take a copy because we will need it...
      _bj_remove_from_tiles(jetB);
      _tj_set_jetinfo(jetB, nn); // also registers the jet in the tiling
    } else {
      _do_iB_recombination_step(jetA->_jets_index, diJ_min);
      _bj_remove_from_tiles(jetA);
    }
    int n_near_tiles = 0;
    _add_neighbours_to_tile_union(jetA->tile_index, tile_union, n_near_tiles);
    if (jetB != NULL) {
      bool sort_it = false;
      if (jetB->tile_index != jetA->tile_index) {
	sort_it = true;
	_add_neighbours_to_tile_union(jetB->tile_index,tile_union,n_near_tiles);
      }
      if (oldB.tile_index != jetA->tile_index && 
	  oldB.tile_index != jetB->tile_index) {
	sort_it = true;
	_add_neighbours_to_tile_union(oldB.tile_index,tile_union,n_near_tiles);
      }
      if (sort_it) {
	// sort the tiles before then compressing the list
	sort(tile_union.begin(), tile_union.begin()+n_near_tiles);
	// and now condense the list
	int nnn = 1;
	for (int i = 1; i < n_near_tiles; i++) {
	  if (tile_union[i] != tile_union[nnn-1]) {
	    tile_union[nnn] = tile_union[i]; 
	    nnn++;
	  }
	}
	n_near_tiles = nnn;
      }
    }
    tail--; n--;
    if (jetA == tail) {
    } else {
      *jetA = *tail;
      diJ[jetA - head] = diJ[tail-head];
      if (jetA->previous == NULL) {
	_tiles[jetA->tile_index].head = jetA;
      } else {
	jetA->previous->next = jetA;
      }
      if (jetA->next != NULL) {jetA->next->previous = jetA;}
    }
    for (int itile = 0; itile < n_near_tiles; itile++) {
      Tile * tile_ptr = &_tiles[tile_union[itile]];
      for (TiledJet * jetI = tile_ptr->head; jetI != NULL; jetI = jetI->next) {
	// see if jetI had jetA or jetB as a NN -- if so recalculate the NN
	if (jetI->NN == jetA || (jetI->NN == jetB && jetB != NULL)) {
	  jetI->NN_dist = _R2;
	  jetI->NN      = NULL;
	  // now go over tiles that are neighbours of I (include own tile)
	  for (Tile ** near_tile  = tile_ptr->begin_tiles; 
	               near_tile != tile_ptr->end_tiles; near_tile++) {
	    // and then over the contents of that tile
	    for (TiledJet * jetJ  = (*near_tile)->head; 
                            jetJ != NULL; jetJ = jetJ->next) {
	      double dist = _bj_dist(jetI,jetJ);
	      if (dist < jetI->NN_dist && jetJ != jetI) {
		jetI->NN_dist = dist; jetI->NN = jetJ;
	      }
	    }
	  }
	  diJ[jetI-head] = _bj_diJ(jetI); // update diJ 
	}
	// check whether new jetB is closer than jetI's current NN and
	// if need to update things
	if (jetB != NULL) {
	  double dist = _bj_dist(jetI,jetB);
	  if (dist < jetI->NN_dist) {
	    if (jetI != jetB) {
	      jetI->NN_dist = dist;
	      jetI->NN = jetB;
	      diJ[jetI-head] = _bj_diJ(jetI); // update diJ...
	    }
	  }
	  if (dist < jetB->NN_dist) {
	    if (jetI != jetB) {
	      jetB->NN_dist = dist;
	      jetB->NN      = jetI;}
	  }
	}
      }
    }
    if (jetB != NULL) {diJ[jetB-head] = _bj_diJ(jetB);}
    for (Tile ** near_tile = _tiles[tail->tile_index].begin_tiles; 
	         near_tile!= _tiles[tail->tile_index].end_tiles; near_tile++){
      for (TiledJet * jetJ = (*near_tile)->head; 
	             jetJ != NULL; jetJ = jetJ->next) {
	if (jetJ->NN == tail) {jetJ->NN = jetA;}
      }
    }
    if (jetB != NULL) {diJ[jetB-head] = _bj_diJ(jetB);}
  }
  delete[] diJ;
  delete[] briefjets;
}
void ClusterSequence::_faster_tiled_N2_cluster() {
  _initialise_tiles();
  int n = _jets.size();
  TiledJet * briefjets = new TiledJet[n];
  TiledJet * jetA = briefjets, * jetB;
  TiledJet oldB;
  oldB.tile_index=0; // prevents a gcc warning
  vector<int> tile_union(3*n_tile_neighbours);
  for (int i = 0; i< n; i++) {
    _tj_set_jetinfo(jetA, i);
    jetA++; // move on to next entry of briefjets
  }
  TiledJet * head = briefjets; // a nicer way of naming start
  vector<Tile>::const_iterator tile;
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      for (jetB = tile->head; jetB != jetA; jetB = jetB->next) {
	double dist = _bj_dist(jetA,jetB);
	if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
      }
    }
    for (Tile ** RTile = tile->RH_tiles; RTile != tile->end_tiles; RTile++) {
      for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
	for (jetB = (*RTile)->head; jetB != NULL; jetB = jetB->next) {
	  double dist = _bj_dist(jetA,jetB);
	  if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	  if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
	}
      }
    }
  }
  struct diJ_plus_link {
    double     diJ; // the distance
    TiledJet * jet; // the jet (i) for which we've found this distance
  };
  diJ_plus_link * diJ = new diJ_plus_link[n];
  jetA = head;
  for (int i = 0; i < n; i++) {
    diJ[i].diJ = _bj_diJ(jetA); // kt distance * R^2
    diJ[i].jet = jetA;  // our compact diJ table will not be in	     
    jetA->diJ_posn = i; // one-to-one corresp. with non-compact jets,
    jetA++; // have jetA follow i 
  }
  int history_location = n-1;
  while (n > 0) {
    diJ_plus_link * best, *stop; // pointers a bit faster than indices
    double diJ_min = diJ[0].diJ; // initialise the best one here.
    best = diJ;                  // and here
    stop = diJ+n;
    for (diJ_plus_link * here = diJ+1; here != stop; here++) {
      if (here->diJ < diJ_min) {best = here; diJ_min  = here->diJ;}
    }
    history_location++;
    jetA = best->jet;
    jetB = jetA->NN;
    diJ_min *= _invR2; 
    if (jetB != NULL) {
      if (jetA < jetB) {std::swap(jetA,jetB);}
      int nn; // new jet index
      _do_ij_recombination_step(jetA->_jets_index, jetB->_jets_index, diJ_min, nn);
      _bj_remove_from_tiles(jetA);
      oldB = * jetB;  // take a copy because we will need it...
      _bj_remove_from_tiles(jetB);
      _tj_set_jetinfo(jetB, nn); // cause jetB to become _jets[nn]
    } else {
      _do_iB_recombination_step(jetA->_jets_index, diJ_min);
      _bj_remove_from_tiles(jetA);
    }
    int n_near_tiles = 0;
    _add_untagged_neighbours_to_tile_union(jetA->tile_index, 
					   tile_union, n_near_tiles);
    if (jetB != NULL) {
      if (jetB->tile_index != jetA->tile_index) {
	_add_untagged_neighbours_to_tile_union(jetB->tile_index,
					       tile_union,n_near_tiles);
      }
      if (oldB.tile_index != jetA->tile_index && 
	  oldB.tile_index != jetB->tile_index) {
	_add_untagged_neighbours_to_tile_union(oldB.tile_index,
					       tile_union,n_near_tiles);
      }
    }
    n--;
    diJ[n].jet->diJ_posn = jetA->diJ_posn;
    diJ[jetA->diJ_posn] = diJ[n];
    for (int itile = 0; itile < n_near_tiles; itile++) {
      Tile * tile_ptr = &_tiles[tile_union[itile]];
      tile_ptr->tagged = false; // reset tag, since we're done with unions
      for (TiledJet * jetI = tile_ptr->head; jetI != NULL; jetI = jetI->next) {
	// see if jetI had jetA or jetB as a NN -- if so recalculate the NN
	if (jetI->NN == jetA || (jetI->NN == jetB && jetB != NULL)) {
	  jetI->NN_dist = _R2;
	  jetI->NN      = NULL;
	  // now go over tiles that are neighbours of I (include own tile)
	  for (Tile ** near_tile  = tile_ptr->begin_tiles; 
	               near_tile != tile_ptr->end_tiles; near_tile++) {
	    // and then over the contents of that tile
	    for (TiledJet * jetJ  = (*near_tile)->head; 
                            jetJ != NULL; jetJ = jetJ->next) {
	      double dist = _bj_dist(jetI,jetJ);
	      if (dist < jetI->NN_dist && jetJ != jetI) {
		jetI->NN_dist = dist; jetI->NN = jetJ;
	      }
	    }
	  }
	  diJ[jetI->diJ_posn].diJ = _bj_diJ(jetI); // update diJ kt-dist
	}
	// check whether new jetB is closer than jetI's current NN and
	// if jetI is closer than jetB's current (evolving) nearest
	// neighbour. Where relevant update things
	if (jetB != NULL) {
	  double dist = _bj_dist(jetI,jetB);
	  if (dist < jetI->NN_dist) {
	    if (jetI != jetB) {
	      jetI->NN_dist = dist;
	      jetI->NN = jetB;
	      diJ[jetI->diJ_posn].diJ = _bj_diJ(jetI); // update diJ...
	    }
	  }
	  if (dist < jetB->NN_dist) {
	    if (jetI != jetB) {
	      jetB->NN_dist = dist;
	      jetB->NN      = jetI;}
	  }
	}
      }
    }
    if (jetB != NULL) {diJ[jetB->diJ_posn].diJ = _bj_diJ(jetB);}
  }
  delete[] diJ;
  delete[] briefjets;
}
void ClusterSequence::_minheap_faster_tiled_N2_cluster() {
  _initialise_tiles();
  int n = _jets.size();
  TiledJet * briefjets = new TiledJet[n];
  TiledJet * jetA = briefjets, * jetB;
  TiledJet oldB;
  oldB.tile_index=0; // prevents a gcc warning
  vector<int> tile_union(3*n_tile_neighbours);
  for (int i = 0; i< n; i++) {
    _tj_set_jetinfo(jetA, i);
    jetA++; // move on to next entry of briefjets
  }
  TiledJet * head = briefjets; // a nicer way of naming start
  vector<Tile>::const_iterator tile;
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      for (jetB = tile->head; jetB != jetA; jetB = jetB->next) {
	double dist = _bj_dist(jetA,jetB);
	if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
      }
    }
    for (Tile ** RTile = tile->RH_tiles; RTile != tile->end_tiles; RTile++) {
      for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
	for (jetB = (*RTile)->head; jetB != NULL; jetB = jetB->next) {
	  double dist = _bj_dist(jetA,jetB);
	  if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	  if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
	}
      }
    }
  }
  vector<double> diJs(n);
  for (int i = 0; i < n; i++) {
    diJs[i] = _bj_diJ(&briefjets[i]);
    briefjets[i].label_minheap_update_done();
  }
  MinHeap minheap(diJs);
  vector<TiledJet *> jets_for_minheap;
  jets_for_minheap.reserve(n); 
  int history_location = n-1;
  while (n > 0) {
    double diJ_min = minheap.minval() *_invR2;
    jetA = head + minheap.minloc();
    history_location++;
    jetB = jetA->NN;
    if (jetB != NULL) {
      if (jetA < jetB) {std::swap(jetA,jetB);}
      int nn; // new jet index
      _do_ij_recombination_step(jetA->_jets_index, jetB->_jets_index, diJ_min, nn);
      _bj_remove_from_tiles(jetA);
      oldB = * jetB;  // take a copy because we will need it...
      _bj_remove_from_tiles(jetB);
      _tj_set_jetinfo(jetB, nn); // cause jetB to become _jets[nn]
    } else {
      _do_iB_recombination_step(jetA->_jets_index, diJ_min);
      _bj_remove_from_tiles(jetA);
    }
    minheap.remove(jetA-head);
    int n_near_tiles = 0;
    _add_untagged_neighbours_to_tile_union(jetA->tile_index, 
					   tile_union, n_near_tiles);
    if (jetB != NULL) {
      if (jetB->tile_index != jetA->tile_index) {
	_add_untagged_neighbours_to_tile_union(jetB->tile_index,
					       tile_union,n_near_tiles);
      }
      if (oldB.tile_index != jetA->tile_index && 
	  oldB.tile_index != jetB->tile_index) {
	// GS: the line below generates a warning that oldB.tile_index
	// may be used uninitialised. However, to reach this point, we
	// ned jetB != NULL (see test a few lines above) and is jetB
	// !=NULL, one would have gone through "oldB = *jetB before
	// (see piece of code ~20 line above), so the index is
	// initialised. We do not do anything to avoid the warning to
	// avoid any potential speed impact.
	_add_untagged_neighbours_to_tile_union(oldB.tile_index,
					       tile_union,n_near_tiles);
      }
      jetB->label_minheap_update_needed();
      jets_for_minheap.push_back(jetB);
    }
    for (int itile = 0; itile < n_near_tiles; itile++) {
      Tile * tile_ptr = &_tiles[tile_union[itile]];
      tile_ptr->tagged = false; // reset tag, since we're done with unions
      for (TiledJet * jetI = tile_ptr->head; jetI != NULL; jetI = jetI->next) {
	// see if jetI had jetA or jetB as a NN -- if so recalculate the NN
	if (jetI->NN == jetA || (jetI->NN == jetB && jetB != NULL)) {
	  jetI->NN_dist = _R2;
	  jetI->NN      = NULL;
	  // label jetI as needing heap action...
	  if (!jetI->minheap_update_needed()) {
	    jetI->label_minheap_update_needed();
	    jets_for_minheap.push_back(jetI);}
	  // now go over tiles that are neighbours of I (include own tile)
	  for (Tile ** near_tile  = tile_ptr->begin_tiles; 
	               near_tile != tile_ptr->end_tiles; near_tile++) {
	    // and then over the contents of that tile
	    for (TiledJet * jetJ  = (*near_tile)->head; 
                            jetJ != NULL; jetJ = jetJ->next) {
	      double dist = _bj_dist(jetI,jetJ);
	      if (dist < jetI->NN_dist && jetJ != jetI) {
		jetI->NN_dist = dist; jetI->NN = jetJ;
	      }
	    }
	  }
	}
	// check whether new jetB is closer than jetI's current NN and
	// if jetI is closer than jetB's current (evolving) nearest
	// neighbour. Where relevant update things
	if (jetB != NULL) {
	  double dist = _bj_dist(jetI,jetB);
	  if (dist < jetI->NN_dist) {
	    if (jetI != jetB) {
	      jetI->NN_dist = dist;
	      jetI->NN = jetB;
	      // label jetI as needing heap action...
	      if (!jetI->minheap_update_needed()) {
		jetI->label_minheap_update_needed();
		jets_for_minheap.push_back(jetI);}
	    }
	  }
	  if (dist < jetB->NN_dist) {
	    if (jetI != jetB) {
	      jetB->NN_dist = dist;
	      jetB->NN      = jetI;}
	  }
	}
      }
    }
    while (jets_for_minheap.size() > 0) {
      TiledJet * jetI = jets_for_minheap.back(); 
      jets_for_minheap.pop_back();
      minheap.update(jetI-head, _bj_diJ(jetI));
      jetI->label_minheap_update_done();
    }
    n--;
  }
  delete[] briefjets;
}
FJCORE_END_NAMESPACE
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
using namespace std;
CompositeJetStructure::CompositeJetStructure(const std::vector<PseudoJet> & initial_pieces, 
					     const JetDefinition::Recombiner * recombiner)
  : _pieces(initial_pieces){
  if (recombiner){};  // ugly trick to prevent a gcc warning
  _area_4vector_ptr = 0;
}
std::string CompositeJetStructure::description() const{ 
  string str = "Composite PseudoJet";
  return str; 
}
bool CompositeJetStructure::has_constituents() const{
  return _pieces.size()!=0;
}
std::vector<PseudoJet> CompositeJetStructure::constituents(const PseudoJet & /*jet*/) const{
  vector<PseudoJet> all_constituents;
  for (unsigned i = 0; i < _pieces.size(); i++) {
    if (_pieces[i].has_constituents()){
      vector<PseudoJet> constits = _pieces[i].constituents();
      copy(constits.begin(), constits.end(), back_inserter(all_constituents));
    } else {
      all_constituents.push_back(_pieces[i]);
    }
  }
  return all_constituents;
}
std::vector<PseudoJet> CompositeJetStructure::pieces(const PseudoJet & /*jet*/) const{
  return _pieces;
}
FJCORE_END_NAMESPACE      // defined in fastjet/internal/base.hh
#include <sstream>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
using namespace std;
bool Error::_print_errors = true;
bool Error::_print_backtrace = false;
ostream * Error::_default_ostr = & cerr;
#if (!defined(FJCORE_HAVE_EXECINFO_H)) || defined(__FJCORE__)
  LimitedWarning Error::_execinfo_undefined;
#endif
Error::Error(const std::string & message_in) {
  _message = message_in; 
  if (_print_errors && _default_ostr){
    ostringstream oss;
    oss << "fjcore::Error:  "<< message_in << endl;
    *_default_ostr << oss.str();
    _default_ostr->flush(); 
  }
}
void Error::set_print_backtrace(bool enabled) {
#if (!defined(FJCORE_HAVE_EXECINFO_H)) || defined(__FJCORE__)
  if (enabled) {
    _execinfo_undefined.warn("Error::set_print_backtrace(true) will not work with this build of FastJet");
  }
#endif    
  _print_backtrace = enabled;
}
FJCORE_END_NAMESPACE
#include <string>
#include <sstream>
using namespace std;
FJCORE_BEGIN_NAMESPACE
FJCORE_END_NAMESPACE
#include<sstream>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
using namespace std;
const double JetDefinition::max_allowable_R = 1000.0;
JetDefinition::JetDefinition(JetAlgorithm jet_algorithm_in, 
			     double R_in, 
			     RecombinationScheme recomb_scheme_in,
			     Strategy strategy_in,
                             int nparameters) :
  _jet_algorithm(jet_algorithm_in), _Rparam(R_in), _strategy(strategy_in) {
  if (_jet_algorithm == ee_kt_algorithm) {
    _Rparam = 4.0; // introduce a fictional R that ensures that
  } else {
    if (R_in > max_allowable_R) {
      ostringstream oss;
      oss << "Requested R = " << R_in << " for jet definition is larger than max_allowable_R = " << max_allowable_R;
      throw Error(oss.str());
    }
  }
  unsigned int nparameters_expected = n_parameters_for_algorithm(jet_algorithm_in);
  if (nparameters != (int) nparameters_expected){
    ostringstream oss;
    oss << "The jet algorithm you requested ("
        << jet_algorithm_in << ") should be constructed with " << nparameters_expected 
        << " parameter(s) but was called with " << nparameters << " parameter(s)\n";
    throw Error(oss.str()); 
  }
  assert (_strategy  != plugin_strategy);
  _plugin = NULL;
  set_recombination_scheme(recomb_scheme_in);
  set_extra_param(0.0); // make sure it's defined
}
bool JetDefinition::is_spherical() const {
  if (jet_algorithm() == plugin_algorithm) {
    return plugin()->is_spherical();
  } else {
    return (jet_algorithm() == ee_kt_algorithm ||  // as of 2013-02-14, the two
            jet_algorithm() == ee_genkt_algorithm  // native spherical algorithms
            );
  }
}
string JetDefinition::description() const {
  ostringstream name;
  name << description_no_recombiner();
  if ((jet_algorithm() == plugin_algorithm) || (jet_algorithm() == undefined_jet_algorithm)){
    return name.str();
  }
  if (n_parameters_for_algorithm(jet_algorithm()) == 0)
    name << " with ";
  else 
    name << " and ";
  name << recombiner()->description();
  return name.str();
}
string JetDefinition::description_no_recombiner() const {
  ostringstream name;
  if (jet_algorithm() == plugin_algorithm) {
    return plugin()->description();
  } else if (jet_algorithm() == undefined_jet_algorithm) {
    return "uninitialised JetDefinition (jet_algorithm=undefined_jet_algorithm)" ;
  }
  name << algorithm_description(jet_algorithm());
  switch (n_parameters_for_algorithm(jet_algorithm())){
  case 0: name << " (NB: no R)"; break;
  case 1: name << " with R = " << R(); break; // the parameter is always R
  case 2: 
    name << " with R = " << R();
    if (jet_algorithm() == cambridge_for_passive_algorithm){
      name << "and a special hack whereby particles with kt < " 
           << extra_param() << "are treated as passive ghosts";
    } else {
      name << ", p = " << extra_param();
    }
  };
  return name.str();
}
string JetDefinition::algorithm_description(const JetAlgorithm jet_alg){
  ostringstream name;
  switch (jet_alg){
  case plugin_algorithm:                return "plugin algorithm";
  case kt_algorithm:                    return "Longitudinally invariant kt algorithm";
  case cambridge_algorithm:             return "Longitudinally invariant Cambridge/Aachen algorithm";
  case antikt_algorithm:                return "Longitudinally invariant anti-kt algorithm";
  case genkt_algorithm:                 return "Longitudinally invariant generalised kt algorithm";
  case cambridge_for_passive_algorithm: return "Longitudinally invariant Cambridge/Aachen algorithm";
  case ee_kt_algorithm:                 return "e+e- kt (Durham) algorithm (NB: no R)";
  case ee_genkt_algorithm:              return "e+e- generalised kt algorithm";
  case undefined_jet_algorithm:         return "undefined jet algorithm";
  default:
    throw Error("JetDefinition::algorithm_description(): unrecognized jet_algorithm");
  };
}
unsigned int JetDefinition::n_parameters_for_algorithm(const JetAlgorithm jet_alg){
  switch (jet_alg) {
  case ee_kt_algorithm:    return 0;
  case genkt_algorithm:
  case ee_genkt_algorithm: return 2;
  default:                 return 1;
  };
}
void JetDefinition::set_recombination_scheme(
                               RecombinationScheme recomb_scheme) {
  _default_recombiner = JetDefinition::DefaultRecombiner(recomb_scheme);
  if (_shared_recombiner) _shared_recombiner.reset();
  _recombiner = 0;
}
void JetDefinition::set_recombiner(const JetDefinition &other_jet_def){
  assert(other_jet_def._recombiner || 
         other_jet_def.recombination_scheme() != external_scheme);
  if (other_jet_def._recombiner == 0){
    set_recombination_scheme(other_jet_def.recombination_scheme());
    return;
  }
  _recombiner = other_jet_def._recombiner;
  _default_recombiner = DefaultRecombiner(external_scheme);
  _shared_recombiner.reset(other_jet_def._shared_recombiner);
}
bool JetDefinition::has_same_recombiner(const JetDefinition &other_jd) const{
  const RecombinationScheme & scheme = recombination_scheme();
  if (other_jd.recombination_scheme() != scheme) return false;
  return (scheme != external_scheme) 
    || (recombiner() == other_jd.recombiner());
}
void JetDefinition::delete_recombiner_when_unused(){
  if (_recombiner == 0){
    throw Error("tried to call JetDefinition::delete_recombiner_when_unused() for a JetDefinition without a user-defined recombination scheme");
  } else if (_shared_recombiner.get()) {
    throw Error("Error in JetDefinition::delete_recombiner_when_unused: the recombiner is already scheduled for deletion when unused (or was already set as shared)");
  }
  _shared_recombiner.reset(_recombiner);
}
void JetDefinition::delete_plugin_when_unused(){
  if (_plugin == 0){
    throw Error("tried to call JetDefinition::delete_plugin_when_unused() for a JetDefinition without a plugin");
  }
  _plugin_shared.reset(_plugin);
}
string JetDefinition::DefaultRecombiner::description() const {
  switch(_recomb_scheme) {
  case E_scheme:
    return "E scheme recombination";
  case pt_scheme:
    return "pt scheme recombination";
  case pt2_scheme:
    return "pt2 scheme recombination";
  case Et_scheme:
    return "Et scheme recombination";
  case Et2_scheme:
    return "Et2 scheme recombination";
  case BIpt_scheme:
    return "boost-invariant pt scheme recombination";
  case BIpt2_scheme:
    return "boost-invariant pt2 scheme recombination";
  case WTA_pt_scheme:
    return "pt-ordered Winner-Takes-All recombination";
  case WTA_modp_scheme:
    return "|3-momentum|-ordered Winner-Takes-All recombination";
  default:
    ostringstream err;
    err << "DefaultRecombiner: unrecognized recombination scheme " 
        << _recomb_scheme;
    throw Error(err.str());
  }
}
void JetDefinition::DefaultRecombiner::recombine(
           const PseudoJet & pa, const PseudoJet & pb,
           PseudoJet & pab) const {
  double weighta, weightb;
  switch(_recomb_scheme) {
  case E_scheme:
    pab.reset(pa.px()+pb.px(),
    	      pa.py()+pb.py(),
    	      pa.pz()+pb.pz(),
    	      pa.E ()+pb.E ());
    return;
  case pt_scheme:
  case Et_scheme:
  case BIpt_scheme:
    weighta = pa.perp(); 
    weightb = pb.perp();
    break;
  case pt2_scheme:
  case Et2_scheme:
  case BIpt2_scheme:
    weighta = pa.perp2(); 
    weightb = pb.perp2();
    break;
  case WTA_pt_scheme:{
    const PseudoJet & phard = (pa.pt2() >= pb.pt2()) ? pa : pb;
    pab.reset_PtYPhiM(pa.pt()+pb.pt(), 
                      phard.rap(), phard.phi(), phard.m());
    return;}
  case WTA_modp_scheme:{
    bool a_hardest = (pa.modp2() >= pb.modp2());
    const PseudoJet & phard = a_hardest ? pa : pb;
    const PseudoJet & psoft = a_hardest ? pb : pa;
    double modp_hard = phard.modp();
    double modp_ab = modp_hard + psoft.modp();
    if (phard.modp2()==0.0){
      pab.reset(0.0, 0.0, 0.0, phard.m());
    } else {
      double scale = modp_ab/modp_hard;
      pab.reset(phard.px()*scale, phard.py()*scale, phard.pz()*scale,
                sqrt(modp_ab*modp_ab + phard.m2()));
    }
    return;}
  default:
    ostringstream err;
    err << "DefaultRecombiner: unrecognized recombination scheme " 
        << _recomb_scheme;
    throw Error(err.str());
  }
  double perp_ab = pa.perp() + pb.perp();
  if (perp_ab != 0.0) { // weights also non-zero...
    double y_ab    = (weighta * pa.rap() + weightb * pb.rap())/(weighta+weightb);
    double phi_a = pa.phi(), phi_b = pb.phi();
    if (phi_a - phi_b > pi)  phi_b += twopi;
    if (phi_a - phi_b < -pi) phi_b -= twopi;
    double phi_ab = (weighta * phi_a + weightb * phi_b)/(weighta+weightb);
    pab.reset_PtYPhiM(perp_ab,y_ab,phi_ab);
  } else { // weights are zero
    pab.reset(0.0, 0.0, 0.0, 0.0);
  }
}
void JetDefinition::DefaultRecombiner::preprocess(PseudoJet & p) const {
  switch(_recomb_scheme) {
  case E_scheme:
  case BIpt_scheme:
  case BIpt2_scheme:
  case WTA_pt_scheme:
  case WTA_modp_scheme:
    break;
  case pt_scheme:
  case pt2_scheme:
    {
      double newE = sqrt(p.perp2()+p.pz()*p.pz());
      p.reset_momentum(p.px(), p.py(), p.pz(), newE);
    }
    break;
  case Et_scheme:
  case Et2_scheme:
    {
      double rescale = p.E()/sqrt(p.perp2()+p.pz()*p.pz());
      p.reset_momentum(rescale*p.px(), rescale*p.py(), rescale*p.pz(), p.E());
    }
    break;
  default:
    ostringstream err;
    err << "DefaultRecombiner: unrecognized recombination scheme " 
        << _recomb_scheme;
    throw Error(err.str());
  }
}
void JetDefinition::Plugin::set_ghost_separation_scale(double /*scale*/) const {
  throw Error("set_ghost_separation_scale not supported");
}
PseudoJet join(const vector<PseudoJet> & pieces, const JetDefinition::Recombiner & recombiner){
  PseudoJet result;  // automatically initialised to 0
  if (pieces.size()>0){
    result = pieces[0];
    for (unsigned int i=1; i<pieces.size(); i++)
      recombiner.plus_equal(result, pieces[i]);
  }
  CompositeJetStructure *cj_struct = new CompositeJetStructure(pieces, &recombiner);
  result.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(cj_struct));
  return result;
}
PseudoJet join(const PseudoJet & j1, 
	       const JetDefinition::Recombiner & recombiner){
  return join(vector<PseudoJet>(1,j1), recombiner);
}
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, 
	       const JetDefinition::Recombiner & recombiner){
  vector<PseudoJet> pieces;
  pieces.push_back(j1);
  pieces.push_back(j2);
  return join(pieces, recombiner);
}
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3, 
	       const JetDefinition::Recombiner & recombiner){
  vector<PseudoJet> pieces;
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  return join(pieces, recombiner);
}
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3, const PseudoJet & j4,
	       const JetDefinition::Recombiner & recombiner){
  vector<PseudoJet> pieces;
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  pieces.push_back(j4);
  return join(pieces, recombiner);
}
FJCORE_END_NAMESPACE
#include <sstream>
#include <limits>
using namespace std;
FJCORE_BEGIN_NAMESPACE
ostream * LimitedWarning::_default_ostr = &cerr;
std::list< LimitedWarning::Summary > LimitedWarning::_global_warnings_summary;
int LimitedWarning::_max_warn_default = 5;
void LimitedWarning::warn(const char * warning, std::ostream * ostr) {
  if (_this_warning_summary == 0) {
    _global_warnings_summary.push_back(Summary(warning, 0));
    _this_warning_summary = & (_global_warnings_summary.back());
  }
  if (_n_warn_so_far < _max_warn) {
    ostringstream warnstr;
    warnstr << "WARNING from FastJet: ";
    warnstr << warning;
    _n_warn_so_far++;
    if (_n_warn_so_far == _max_warn) warnstr << " (LAST SUCH WARNING)";
    warnstr << std::endl;
    if (ostr) {
      (*ostr) << warnstr.str();
      ostr->flush(); // get something written to file even if the program aborts
    }
  }
  if (_this_warning_summary->second < numeric_limits<unsigned>::max()) {
    _this_warning_summary->second++;
  }
}
string LimitedWarning::summary() {
  ostringstream str;
  for (list<Summary>::const_iterator it = _global_warnings_summary.begin();
       it != _global_warnings_summary.end(); it++) {
    str << it->second << " times: " << it->first << endl;
  }
  return str.str();
}
FJCORE_END_NAMESPACE
#include<iostream>
#include<cmath>
#include<limits>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
using namespace std;
void MinHeap::initialise(const std::vector<double> & values){
  for (unsigned i = values.size(); i < _heap.size(); i++) {
    _heap[i].value = std::numeric_limits<double>::max();
    _heap[i].minloc = &(_heap[i]);
  }
  for (unsigned i = 0; i < values.size(); i++) {
    _heap[i].value = values[i];
    _heap[i].minloc = &(_heap[i]);
  }
  for (unsigned i = _heap.size()-1; i > 0; i--) {
    ValueLoc * parent = &(_heap[(i-1)/2]);
    ValueLoc * here   = &(_heap[i]);
    if (here->minloc->value < parent->minloc->value) {
      parent->minloc = here->minloc;
    }
  }
}
void MinHeap::update(unsigned int loc, double new_value) {
  assert(loc < _heap.size());
  ValueLoc * start = &(_heap[loc]);
  if (start->minloc != start && !(new_value < start->minloc->value)) {
    start->value = new_value;
    return;
  }
  start->value = new_value;
  start->minloc = start;
  bool change_made = true;
  ValueLoc * heap_end = (&(_heap[0])) + _heap.size();
  while(change_made) {
    ValueLoc * here = &(_heap[loc]);
    change_made     = false;
    if (here->minloc == start) {
      here->minloc = here; change_made = true;
    }
    ValueLoc * child = &(_heap[2*loc+1]);
    if (child < heap_end && child->minloc->value < here->minloc->value ) {
      here->minloc = child->minloc;
      change_made = true;}
    child++;
    if (child < heap_end && child->minloc->value < here->minloc->value ) {
      here->minloc = child->minloc;
      change_made = true;}
    if (loc == 0) {break;}
    loc = (loc-1)/2;
  }
}
FJCORE_END_NAMESPACE
#include<valarray>
#include<iostream>
#include<sstream>
#include<cmath>
#include<algorithm>
#include <cstdarg>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
using namespace std;
PseudoJet::PseudoJet(const double px_in, const double py_in, const double pz_in, const double E_in) {
  _E  = E_in ;
  _px = px_in;
  _py = py_in;
  _pz = pz_in;
  this->_finish_init();
  _reset_indices();
}
void PseudoJet::_finish_init () {
  _kt2 = this->px()*this->px() + this->py()*this->py();
  _phi = pseudojet_invalid_phi;
  _rap = pseudojet_invalid_rap;
}
void PseudoJet::_set_rap_phi() const {
  if (_kt2 == 0.0) {
    _phi = 0.0; } 
  else {
    _phi = atan2(this->py(),this->px());
  }
  if (_phi < 0.0) {_phi += twopi;}
  if (_phi >= twopi) {_phi -= twopi;} // can happen if phi=-|eps<1e-15|?
  if (this->E() == abs(this->pz()) && _kt2 == 0) {
    double MaxRapHere = MaxRap + abs(this->pz());
    if (this->pz() >= 0.0) {_rap = MaxRapHere;} else {_rap = -MaxRapHere;}
  } else {
    double effective_m2 = max(0.0,m2()); // force non tachyonic mass
    double E_plus_pz    = _E + abs(_pz); // the safer of p+, p-
    _rap = 0.5*log((_kt2 + effective_m2)/(E_plus_pz*E_plus_pz));
    if (_pz > 0) {_rap = - _rap;}
  }
}
valarray<double> PseudoJet::four_mom() const {
  valarray<double> mom(4);
  mom[0] = _px;
  mom[1] = _py;
  mom[2] = _pz;
  mom[3] = _E ;
  return mom;
}
double PseudoJet::operator () (int i) const {
  switch(i) {
  case X:
    return px();
  case Y:
    return py();
  case Z:
    return pz();
  case T:
    return e();
  default:
    ostringstream err;
    err << "PseudoJet subscripting: bad index (" << i << ")";
    throw Error(err.str());
  }
  return 0.;
}  
double PseudoJet::pseudorapidity() const {
  if (px() == 0.0 && py() ==0.0) return MaxRap;
  if (pz() == 0.0) return 0.0;
  double theta = atan(perp()/pz());
  if (theta < 0) theta += pi;
  return -log(tan(theta/2));
}
PseudoJet operator+ (const PseudoJet & jet1, const PseudoJet & jet2) {
  return PseudoJet(jet1.px()+jet2.px(),
		   jet1.py()+jet2.py(),
		   jet1.pz()+jet2.pz(),
		   jet1.E() +jet2.E()  );
} 
PseudoJet operator- (const PseudoJet & jet1, const PseudoJet & jet2) {
  return PseudoJet(jet1.px()-jet2.px(),
		   jet1.py()-jet2.py(),
		   jet1.pz()-jet2.pz(),
		   jet1.E() -jet2.E()  );
} 
PseudoJet operator* (double coeff, const PseudoJet & jet) {
  jet._ensure_valid_rap_phi(); 
  PseudoJet coeff_times_jet(jet);
  coeff_times_jet *= coeff;
  return coeff_times_jet;
} 
PseudoJet operator* (const PseudoJet & jet, double coeff) {
  return coeff*jet;
} 
PseudoJet operator/ (const PseudoJet & jet, double coeff) {
  return (1.0/coeff)*jet;
} 
void PseudoJet::operator*=(double coeff) {
  _ensure_valid_rap_phi(); 
  _px *= coeff;
  _py *= coeff;
  _pz *= coeff;
  _E  *= coeff;
  _kt2*= coeff*coeff;
}
void PseudoJet::operator/=(double coeff) {
  (*this) *= 1.0/coeff;
}
void PseudoJet::operator+=(const PseudoJet & other_jet) {
  _px += other_jet._px;
  _py += other_jet._py;
  _pz += other_jet._pz;
  _E  += other_jet._E ;
  _finish_init(); // we need to recalculate phi,rap,kt2
}
void PseudoJet::operator-=(const PseudoJet & other_jet) {
  _px -= other_jet._px;
  _py -= other_jet._py;
  _pz -= other_jet._pz;
  _E  -= other_jet._E ;
  _finish_init(); // we need to recalculate phi,rap,kt2
}
bool operator==(const PseudoJet & a, const PseudoJet & b) {
  if (a.px() != b.px()) return false;
  if (a.py() != b.py()) return false;
  if (a.pz() != b.pz()) return false;
  if (a.E () != b.E ()) return false;
  if (a.user_index()    != b.user_index()) return false;
  if (a.cluster_hist_index() != b.cluster_hist_index()) return false;
  if (a.user_info_ptr() != b.user_info_ptr()) return false;
  if (a.structure_ptr() != b.structure_ptr()) return false;
  return true;
}
bool operator==(const PseudoJet & jet, const double val) {
  if (val != 0) 
    throw Error("comparing a PseudoJet with a non-zero constant (double) is not allowed.");
  return (jet.px() == 0 && jet.py() == 0 && 
	  jet.pz() == 0 && jet.E() == 0);
}
PseudoJet & PseudoJet::boost(const PseudoJet & prest) {
  if (prest.px() == 0.0 && prest.py() == 0.0 && prest.pz() == 0.0) 
    return *this;
  double m_local = prest.m();
  assert(m_local != 0);
  double pf4  = (  px()*prest.px() + py()*prest.py()
                 + pz()*prest.pz() + E()*prest.E() )/m_local;
  double fn   = (pf4 + E()) / (prest.E() + m_local);
  _px +=  fn*prest.px();
  _py +=  fn*prest.py();
  _pz +=  fn*prest.pz();
  _E = pf4;
  _finish_init(); // we need to recalculate phi,rap,kt2
  return *this;
}
PseudoJet & PseudoJet::unboost(const PseudoJet & prest) {
  if (prest.px() == 0.0 && prest.py() == 0.0 && prest.pz() == 0.0) 
    return *this;
  double m_local = prest.m();
  assert(m_local != 0);
  double pf4  = ( -px()*prest.px() - py()*prest.py()
                 - pz()*prest.pz() + E()*prest.E() )/m_local;
  double fn   = (pf4 + E()) / (prest.E() + m_local);
  _px -=  fn*prest.px();
  _py -=  fn*prest.py();
  _pz -=  fn*prest.pz();
  _E = pf4;
  _finish_init(); // we need to recalculate phi,rap,kt2
  return *this;
}
bool have_same_momentum(const PseudoJet & jeta, const PseudoJet & jetb) {
  return jeta.px() == jetb.px()
    &&   jeta.py() == jetb.py()
    &&   jeta.pz() == jetb.pz()
    &&   jeta.E()  == jetb.E();
}
void PseudoJet::set_cached_rap_phi(double rap_in, double phi_in) {
  _rap = rap_in; _phi = phi_in;
  if (_phi >= twopi) _phi -= twopi;
  if (_phi < 0)      _phi += twopi;
}
void PseudoJet::reset_momentum_PtYPhiM(double pt_in, double y_in, double phi_in, double m_in) {
  assert(phi_in < 2*twopi && phi_in > -twopi);
  double ptm = (m_in == 0) ? pt_in : sqrt(pt_in*pt_in+m_in*m_in);
  double exprap = exp(y_in);
  double pminus = ptm/exprap;
  double pplus  = ptm*exprap;
  double px_local = pt_in*cos(phi_in);
  double py_local = pt_in*sin(phi_in);
  reset_momentum(px_local,py_local,0.5*(pplus-pminus),0.5*(pplus+pminus));
  set_cached_rap_phi(y_in,phi_in);
}
PseudoJet PtYPhiM(double pt, double y, double phi, double m) {
  assert(phi < 2*twopi && phi > -twopi);
  double ptm = (m == 0) ? pt : sqrt(pt*pt+m*m);
  double exprap = exp(y);
  double pminus = ptm/exprap;
  double pplus  = ptm*exprap;
  double px = pt*cos(phi);
  double py = pt*sin(phi);
  PseudoJet mom(px,py,0.5*(pplus-pminus),0.5*(pplus+pminus));
  mom.set_cached_rap_phi(y,phi);
  return mom;
}
double PseudoJet::kt_distance(const PseudoJet & other) const {
  double distance = min(_kt2, other._kt2);
  double dphi = abs(phi() - other.phi());
  if (dphi > pi) {dphi = twopi - dphi;}
  double drap = rap() - other.rap();
  distance = distance * (dphi*dphi + drap*drap);
  return distance;
}
double PseudoJet::plain_distance(const PseudoJet & other) const {
  double dphi = abs(phi() - other.phi());
  if (dphi > pi) {dphi = twopi - dphi;}
  double drap = rap() - other.rap();
  return (dphi*dphi + drap*drap);
}
double PseudoJet::delta_phi_to(const PseudoJet & other) const {
  double dphi = other.phi() - phi();
  if (dphi >  pi) dphi -= twopi;
  if (dphi < -pi) dphi += twopi;
  return dphi;
}
string PseudoJet::description() const{
  if (!_structure)
    return "standard PseudoJet (with no associated clustering information)";
  return _structure->description();
}
bool PseudoJet::has_associated_cluster_sequence() const{
  return (_structure) && (_structure->has_associated_cluster_sequence());
}
const ClusterSequence* PseudoJet::associated_cluster_sequence() const{
  if (! has_associated_cluster_sequence()) return NULL;
  return _structure->associated_cluster_sequence();
}
bool PseudoJet::has_valid_cluster_sequence() const{
  return (_structure) && (_structure->has_valid_cluster_sequence());
}
const ClusterSequence * PseudoJet::validated_cs() const {
  return validated_structure_ptr()->validated_cs();
}
void PseudoJet::set_structure_shared_ptr(const SharedPtr<PseudoJetStructureBase> &structure_in){
  _structure = structure_in;
}
bool PseudoJet::has_structure() const{
  return bool(_structure);
}
const PseudoJetStructureBase* PseudoJet::structure_ptr() const {
  return _structure.get();
}
PseudoJetStructureBase* PseudoJet::structure_non_const_ptr(){
  return _structure.get();
}
const PseudoJetStructureBase* PseudoJet::validated_structure_ptr() const {
  if (!_structure) 
    throw Error("Trying to access the structure of a PseudoJet which has no associated structure");
  return _structure.get();
}
const SharedPtr<PseudoJetStructureBase> & PseudoJet::structure_shared_ptr() const {
  return _structure;
}
bool PseudoJet::has_partner(PseudoJet &partner) const{
  return validated_structure_ptr()->has_partner(*this, partner);
}
bool PseudoJet::has_child(PseudoJet &child) const{
  return validated_structure_ptr()->has_child(*this, child);
}
bool PseudoJet::has_parents(PseudoJet &parent1, PseudoJet &parent2) const{
  return validated_structure_ptr()->has_parents(*this, parent1, parent2);
}
bool PseudoJet::contains(const PseudoJet &constituent) const{
  return validated_structure_ptr()->object_in_jet(constituent, *this);
}
bool PseudoJet::is_inside(const PseudoJet &jet) const{
  return validated_structure_ptr()->object_in_jet(*this, jet);
}
bool PseudoJet::has_constituents() const{
  return (_structure) && (_structure->has_constituents());
}
vector<PseudoJet> PseudoJet::constituents() const{
  return validated_structure_ptr()->constituents(*this);
}
bool PseudoJet::has_exclusive_subjets() const{
  return (_structure) && (_structure->has_exclusive_subjets());
}
std::vector<PseudoJet> PseudoJet::exclusive_subjets (const double dcut) const {
  return validated_structure_ptr()->exclusive_subjets(*this, dcut);
}
int PseudoJet::n_exclusive_subjets(const double dcut) const {
  return validated_structure_ptr()->n_exclusive_subjets(*this, dcut);
}
std::vector<PseudoJet> PseudoJet::exclusive_subjets_up_to (int nsub) const {
  return validated_structure_ptr()->exclusive_subjets_up_to(*this, nsub);
}
std::vector<PseudoJet> PseudoJet::exclusive_subjets (int nsub) const {
  vector<PseudoJet> subjets = exclusive_subjets_up_to(nsub);
  if (int(subjets.size()) < nsub) {
    ostringstream err;
    err << "Requested " << nsub << " exclusive subjets, but there were only " 
	<< subjets.size() << " particles in the jet";
    throw Error(err.str());
  }
  return subjets;
}
double PseudoJet::exclusive_subdmerge(int nsub) const {
  return validated_structure_ptr()->exclusive_subdmerge(*this, nsub);
}
double PseudoJet::exclusive_subdmerge_max(int nsub) const {
  return validated_structure_ptr()->exclusive_subdmerge_max(*this, nsub);
}
bool PseudoJet::has_pieces() const{
  return ((_structure) && (_structure->has_pieces(*this)));
}
std::vector<PseudoJet> PseudoJet::pieces() const{
  return validated_structure_ptr()->pieces(*this);
}
PseudoJet::InexistentUserInfo::InexistentUserInfo() : Error("you attempted to perform a dynamic cast of a PseudoJet's extra info, but the extra info pointer was null")
{}
void sort_indices(vector<int> & indices, 
			 const vector<double> & values) {
  IndexedSortHelper index_sort_helper(&values);
  sort(indices.begin(), indices.end(), index_sort_helper);
}
vector<PseudoJet> sorted_by_pt(const vector<PseudoJet> & jets) {
  vector<double> minus_kt2(jets.size());
  for (size_t i = 0; i < jets.size(); i++) {minus_kt2[i] = -jets[i].kt2();}
  return objects_sorted_by_values(jets, minus_kt2);
}
vector<PseudoJet> sorted_by_rapidity(const vector<PseudoJet> & jets) {
  vector<double> rapidities(jets.size());
  for (size_t i = 0; i < jets.size(); i++) {rapidities[i] = jets[i].rap();}
  return objects_sorted_by_values(jets, rapidities);
}
vector<PseudoJet> sorted_by_E(const vector<PseudoJet> & jets) {
  vector<double> energies(jets.size());
  for (size_t i = 0; i < jets.size(); i++) {energies[i] = -jets[i].E();}
  return objects_sorted_by_values(jets, energies);
}
vector<PseudoJet> sorted_by_pz(const vector<PseudoJet> & jets) {
  vector<double> pz(jets.size());
  for (size_t i = 0; i < jets.size(); i++) {pz[i] = jets[i].pz();}
  return objects_sorted_by_values(jets, pz);
}
PseudoJet join(const vector<PseudoJet> & pieces){
  PseudoJet result;  // automatically initialised to 0
  for (unsigned int i=0; i<pieces.size(); i++)
    result += pieces[i];
  CompositeJetStructure *cj_struct = new CompositeJetStructure(pieces);
  result.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(cj_struct));
  return result;
}
PseudoJet join(const PseudoJet & j1){
  return join(vector<PseudoJet>(1,j1));
}
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2){
  vector<PseudoJet> pieces;
  pieces.reserve(2);
  pieces.push_back(j1);
  pieces.push_back(j2);
  return join(pieces);
}
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3){
  vector<PseudoJet> pieces;
  pieces.reserve(3);
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  return join(pieces);
}
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3, const PseudoJet & j4){
  vector<PseudoJet> pieces;
  pieces.reserve(4);
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  pieces.push_back(j4);
  return join(pieces);
}
FJCORE_END_NAMESPACE
using namespace std;
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
const ClusterSequence* PseudoJetStructureBase::associated_cluster_sequence() const{
  return NULL;
}
const ClusterSequence * PseudoJetStructureBase::validated_cs() const{
  throw Error("This PseudoJet structure is not associated with a valid ClusterSequence");
}
bool PseudoJetStructureBase::has_partner(const PseudoJet & /*reference */, PseudoJet & /*partner*/) const{
  throw Error("This PseudoJet structure has no implementation for has_partner");
}
bool PseudoJetStructureBase::has_child(const PseudoJet & /*reference*/, PseudoJet & /*child*/) const{
  throw Error("This PseudoJet structure has no implementation for has_child");
}
bool PseudoJetStructureBase::has_parents(const PseudoJet & /*reference*/, PseudoJet &/*parent1*/, PseudoJet &/*parent2*/) const{
  throw Error("This PseudoJet structure has no implementation for has_parents");
}
bool PseudoJetStructureBase::object_in_jet(const PseudoJet & /*reference*/, const PseudoJet & /*jet*/) const{
  throw Error("This PseudoJet structure has no implementation for is_inside");
}
vector<PseudoJet> PseudoJetStructureBase::constituents(const PseudoJet &/*reference*/) const{
  throw Error("This PseudoJet structure has no implementation for constituents");
}
vector<PseudoJet> PseudoJetStructureBase::exclusive_subjets (const PseudoJet & /*reference*/, const double & /*dcut*/) const{
  throw Error("This PseudoJet structure has no implementation for exclusive_subjets");
}
int PseudoJetStructureBase::n_exclusive_subjets(const PseudoJet & /*reference*/, const double & /*dcut*/) const{
  throw Error("This PseudoJet structure has no implementation for n_exclusive_subjets");
}
vector<PseudoJet> PseudoJetStructureBase::exclusive_subjets_up_to (const PseudoJet & /*reference*/, int /*nsub*/) const{
  throw Error("This PseudoJet structure has no implementation for exclusive_subjets");
}
double PseudoJetStructureBase::exclusive_subdmerge(const PseudoJet & /*reference*/, int /*nsub*/) const{
  throw Error("This PseudoJet structure has no implementation for exclusive_submerge");
}
double PseudoJetStructureBase::exclusive_subdmerge_max(const PseudoJet & /*reference*/, int /*nsub*/) const{
  throw Error("This PseudoJet structure has no implementation for exclusive_submerge_max");
}
std::vector<PseudoJet> PseudoJetStructureBase::pieces(const PseudoJet & /*reference*/) const{
  throw Error("This PseudoJet structure has no implementation for pieces");  
}
FJCORE_END_NAMESPACE
#include <sstream>
#include <algorithm>
using namespace std;
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
std::vector<PseudoJet> Selector::operator()(const std::vector<PseudoJet> & jets) const {
  std::vector<PseudoJet> result;
  const SelectorWorker * worker_local = validated_worker();
  if (worker_local->applies_jet_by_jet()) {
    for (std::vector<PseudoJet>::const_iterator jet = jets.begin(); 
         jet != jets.end(); jet++) {
      if (worker_local->pass(*jet)) result.push_back(*jet);
    }
  } else {
    std::vector<const PseudoJet *> jetptrs(jets.size());
    for (unsigned i = 0; i < jets.size(); i++) {
      jetptrs[i] = & jets[i];
    }
    worker_local->terminator(jetptrs);
    for (unsigned i = 0; i < jetptrs.size(); i++) {
      if (jetptrs[i]) result.push_back(jets[i]);
    }
  }
  return result;
}
unsigned int Selector::count(const std::vector<PseudoJet> & jets) const {
  unsigned n = 0;
  const SelectorWorker * worker_local = validated_worker();
  if (worker_local->applies_jet_by_jet()) {
    for (unsigned i = 0; i < jets.size(); i++) {
      if (worker_local->pass(jets[i])) n++;
    }
  } else {
    std::vector<const PseudoJet *> jetptrs(jets.size());
    for (unsigned i = 0; i < jets.size(); i++) {
      jetptrs[i] = & jets[i];
    }
    worker_local->terminator(jetptrs);
    for (unsigned i = 0; i < jetptrs.size(); i++) {
      if (jetptrs[i]) n++;
    }
  }
  return n;
}
PseudoJet Selector::sum(const std::vector<PseudoJet> & jets) const {
  PseudoJet this_sum(0,0,0,0);
  const SelectorWorker * worker_local = validated_worker();
  if (worker_local->applies_jet_by_jet()) {
    for (unsigned i = 0; i < jets.size(); i++) {
      if (worker_local->pass(jets[i])) this_sum += jets[i];
    }
  } else {
    std::vector<const PseudoJet *> jetptrs(jets.size());
    for (unsigned i = 0; i < jets.size(); i++) {
      jetptrs[i] = & jets[i];
    }
    worker_local->terminator(jetptrs);
    for (unsigned i = 0; i < jetptrs.size(); i++) {
      if (jetptrs[i]) this_sum += jets[i];
    }
  }
  return this_sum;
}
double Selector::scalar_pt_sum(const std::vector<PseudoJet> & jets) const {
  double this_sum = 0.0;
  const SelectorWorker * worker_local = validated_worker();
  if (worker_local->applies_jet_by_jet()) {
    for (unsigned i = 0; i < jets.size(); i++) {
      if (worker_local->pass(jets[i])) this_sum += jets[i].pt();
    }
  } else {
    std::vector<const PseudoJet *> jetptrs(jets.size());
    for (unsigned i = 0; i < jets.size(); i++) {
      jetptrs[i] = & jets[i];
    }
    worker_local->terminator(jetptrs);
    for (unsigned i = 0; i < jetptrs.size(); i++) {
      if (jetptrs[i]) this_sum += jets[i].pt();
    }
  }
  return this_sum;
}
void Selector::sift(const std::vector<PseudoJet> & jets,
		    std::vector<PseudoJet> & jets_that_pass,
		    std::vector<PseudoJet> & jets_that_fail
		    ) const {
  const SelectorWorker * worker_local = validated_worker();
  jets_that_pass.clear();
  jets_that_fail.clear();
  if (worker_local->applies_jet_by_jet()) {
    for (unsigned i = 0; i < jets.size(); i++) {
      if (worker_local->pass(jets[i])) {
	jets_that_pass.push_back(jets[i]);
      } else {
	jets_that_fail.push_back(jets[i]);
      }
    }
  } else {
    std::vector<const PseudoJet *> jetptrs(jets.size());
    for (unsigned i = 0; i < jets.size(); i++) {
      jetptrs[i] = & jets[i];
    }
    worker_local->terminator(jetptrs);
    for (unsigned i = 0; i < jetptrs.size(); i++) {
      if (jetptrs[i]) {
	jets_that_pass.push_back(jets[i]);
      } else {
	jets_that_fail.push_back(jets[i]);
      }
    }
  }
}
bool SelectorWorker::has_finite_area() const { 
  if (! is_geometric()) return false;
  double rapmin, rapmax;
  get_rapidity_extent(rapmin, rapmax);
  return (rapmax != std::numeric_limits<double>::infinity())
    &&  (-rapmin != std::numeric_limits<double>::infinity());
}
class SW_Identity : public SelectorWorker {
public:
  SW_Identity(){}
  virtual bool pass(const PseudoJet &) const {
    return true;
  }
  virtual void terminator(vector<const PseudoJet *> &) const {
    return;
  }
  virtual string description() const { return "Identity";}
  virtual bool is_geometric() const { return true;}
};
Selector SelectorIdentity() {
  return Selector(new SW_Identity);
}
class SW_Not : public SelectorWorker {
public:
  SW_Not(const Selector & s) : _s(s) {}
  virtual SelectorWorker* copy(){ return new SW_Not(*this);}
  virtual bool pass(const PseudoJet & jet) const {
    if (!applies_jet_by_jet())
      throw Error("Cannot apply this selector worker to an individual jet");
    return ! _s.pass(jet);
  } 
  virtual bool applies_jet_by_jet() const {return _s.applies_jet_by_jet();}
  virtual void terminator(vector<const PseudoJet *> & jets) const {
    if (applies_jet_by_jet()){
      SelectorWorker::terminator(jets);
      return;
    }
    vector<const PseudoJet *> s_jets = jets;
    _s.worker()->terminator(s_jets);
    for (unsigned int i=0; i<s_jets.size(); i++){
      if (s_jets[i]) jets[i] = NULL;
    }
  }
  virtual string description() const {
    ostringstream ostr;
    ostr << "!(" << _s.description() << ")";
    return ostr.str();
  }
  virtual bool is_geometric() const { return _s.is_geometric();}
  virtual bool takes_reference() const { return _s.takes_reference();}
  virtual void set_reference(const PseudoJet &ref) { _s.set_reference(ref);}
protected:
  Selector _s;
};
Selector operator!(const Selector & s) {
  return Selector(new SW_Not(s));
}
class SW_BinaryOperator: public SelectorWorker {
public:
  SW_BinaryOperator(const Selector & s1, const Selector & s2) : _s1(s1), _s2(s2) {
    _applies_jet_by_jet = _s1.applies_jet_by_jet() && _s2.applies_jet_by_jet();
    _takes_reference = _s1.takes_reference() || _s2.takes_reference();
    _is_geometric = _s1.is_geometric() && _s2.is_geometric();
  }
  virtual bool applies_jet_by_jet() const {return _applies_jet_by_jet;}
  virtual bool takes_reference() const{ 
    return _takes_reference;
  }
  virtual void set_reference(const PseudoJet &centre){
    _s1.set_reference(centre);
    _s2.set_reference(centre);
  }
  virtual bool is_geometric() const { return _is_geometric;} 
protected:
  Selector _s1, _s2;
  bool _applies_jet_by_jet;
  bool _takes_reference;
  bool _is_geometric;
};
class SW_And: public SW_BinaryOperator {
public:
  SW_And(const Selector & s1, const Selector & s2) : SW_BinaryOperator(s1,s2){}
  virtual SelectorWorker* copy(){ return new SW_And(*this);}
  virtual bool pass(const PseudoJet & jet) const {
    if (!applies_jet_by_jet())
      throw Error("Cannot apply this selector worker to an individual jet");
    return _s1.pass(jet) && _s2.pass(jet);
  }
  virtual void terminator(vector<const PseudoJet *> & jets) const {
    if (applies_jet_by_jet()){
      SelectorWorker::terminator(jets);
      return;
    }
    vector<const PseudoJet *> s1_jets = jets;
    _s1.worker()->terminator(s1_jets);
    _s2.worker()->terminator(jets);
    for (unsigned int i=0; i<jets.size(); i++){
      if (! s1_jets[i]) jets[i] = NULL;
    }
  }
  virtual void get_rapidity_extent(double & rapmin, double & rapmax) const {
    double s1min, s1max, s2min, s2max;
    _s1.get_rapidity_extent(s1min, s1max);
    _s2.get_rapidity_extent(s2min, s2max);
    rapmax = min(s1max, s2max);
    rapmin = max(s1min, s2min);
  }
  virtual string description() const {
    ostringstream ostr;
    ostr << "(" << _s1.description() << " && " << _s2.description() << ")";
    return ostr.str();
  }
};
Selector operator&&(const Selector & s1, const Selector & s2) {
  return Selector(new SW_And(s1,s2));
}
class SW_Or: public SW_BinaryOperator {
public:
  SW_Or(const Selector & s1, const Selector & s2) : SW_BinaryOperator(s1,s2) {}
  virtual SelectorWorker* copy(){ return new SW_Or(*this);}
  virtual bool pass(const PseudoJet & jet) const {
    if (!applies_jet_by_jet())
      throw Error("Cannot apply this selector worker to an individual jet");
    return _s1.pass(jet) || _s2.pass(jet);
  }
  virtual bool applies_jet_by_jet() const {
    return _s1.applies_jet_by_jet() && _s2.applies_jet_by_jet();
  }
  virtual void terminator(vector<const PseudoJet *> & jets) const {
    if (applies_jet_by_jet()){
      SelectorWorker::terminator(jets);
      return;
    }
    vector<const PseudoJet *> s1_jets = jets;
    _s1.worker()->terminator(s1_jets);
    _s2.worker()->terminator(jets);
    for (unsigned int i=0; i<jets.size(); i++){
      if (s1_jets[i]) jets[i] = s1_jets[i];
    }
  }
  virtual string description() const {
    ostringstream ostr;
    ostr << "(" << _s1.description() << " || " << _s2.description() << ")";
    return ostr.str();
  }
  virtual void get_rapidity_extent(double & rapmin, double & rapmax) const {
    double s1min, s1max, s2min, s2max;
    _s1.get_rapidity_extent(s1min, s1max);
    _s2.get_rapidity_extent(s2min, s2max);
    rapmax = max(s1max, s2max);
    rapmin = min(s1min, s2min);
  }
};
Selector operator ||(const Selector & s1, const Selector & s2) {
  return Selector(new SW_Or(s1,s2));
}
class SW_Mult: public SW_And {
public:
  SW_Mult(const Selector & s1, const Selector & s2) : SW_And(s1,s2) {}
  virtual SelectorWorker* copy(){ return new SW_Mult(*this);}
  virtual void terminator(vector<const PseudoJet *> & jets) const {
    if (applies_jet_by_jet()){
      SelectorWorker::terminator(jets);
      return;
    }
    _s2.worker()->terminator(jets);
    _s1.worker()->terminator(jets);
  }
  virtual string description() const {
    ostringstream ostr;
    ostr << "(" << _s1.description() << " * " << _s2.description() << ")";
    return ostr.str();
  }
};
Selector operator*(const Selector & s1, const Selector & s2) {
  return Selector(new SW_Mult(s1,s2));
}
class QuantityBase{
public:
  QuantityBase(double q) : _q(q){}
  virtual ~QuantityBase(){}
  virtual double operator()(const PseudoJet & jet ) const =0;
  virtual string description() const =0;
  virtual bool is_geometric() const { return false;}
  virtual double comparison_value() const {return _q;}
  virtual double description_value() const {return comparison_value();}
protected:
  double _q;
};  
class QuantitySquareBase : public QuantityBase{
public:
  QuantitySquareBase(double sqrtq) : QuantityBase(sqrtq*sqrtq), _sqrtq(sqrtq){}
  virtual double description_value() const {return _sqrtq;}
protected:
  double _sqrtq;
};  
template<typename QuantityType>
class SW_QuantityMin : public SelectorWorker{
public:
  SW_QuantityMin(double qmin) : _qmin(qmin) {}
  virtual bool pass(const PseudoJet & jet) const {return _qmin(jet) >= _qmin.comparison_value();}
  virtual string description() const {
    ostringstream ostr;
    ostr << _qmin.description() << " >= " << _qmin.description_value();
    return ostr.str();
  }
  virtual bool is_geometric() const { return _qmin.is_geometric();}
protected:
  QuantityType _qmin;     ///< the cut
};
template<typename QuantityType>
class SW_QuantityMax : public SelectorWorker {
public:
  SW_QuantityMax(double qmax) : _qmax(qmax) {}
  virtual bool pass(const PseudoJet & jet) const {return _qmax(jet) <= _qmax.comparison_value();}
  virtual string description() const {
    ostringstream ostr;
    ostr << _qmax.description() << " <= " << _qmax.description_value();
    return ostr.str();
  }
  virtual bool is_geometric() const { return _qmax.is_geometric();}
protected:
  QuantityType _qmax;   ///< the cut
};
template<typename QuantityType>
class SW_QuantityRange : public SelectorWorker {
public:
  SW_QuantityRange(double qmin, double qmax) : _qmin(qmin), _qmax(qmax) {}
  virtual bool pass(const PseudoJet & jet) const {
    double q = _qmin(jet); // we could identically use _qmax
    return (q >= _qmin.comparison_value()) && (q <= _qmax.comparison_value());
  }
  virtual string description() const {
    ostringstream ostr;
    ostr << _qmin.description_value() << " <= " << _qmin.description() << " <= " << _qmax.description_value();
    return ostr.str();
  }
  virtual bool is_geometric() const { return _qmin.is_geometric();}
protected:
  QuantityType _qmin;   // the lower cut 
  QuantityType _qmax;   // the upper cut
};
class QuantityPt2 : public QuantitySquareBase{
public:
  QuantityPt2(double pt) : QuantitySquareBase(pt){}
  virtual double operator()(const PseudoJet & jet ) const { return jet.perp2();}
  virtual string description() const {return "pt";}
};  
Selector SelectorPtMin(double ptmin) {
  return Selector(new SW_QuantityMin<QuantityPt2>(ptmin));
}
Selector SelectorPtMax(double ptmax) {
  return Selector(new SW_QuantityMax<QuantityPt2>(ptmax));
}
Selector SelectorPtRange(double ptmin, double ptmax) {
  return Selector(new SW_QuantityRange<QuantityPt2>(ptmin, ptmax));
}
class QuantityEt2 : public QuantitySquareBase{
public:
  QuantityEt2(double Et) : QuantitySquareBase(Et){}
  virtual double operator()(const PseudoJet & jet ) const { return jet.Et2();}
  virtual string description() const {return "Et";}
};  
Selector SelectorEtMin(double Etmin) {
  return Selector(new SW_QuantityMin<QuantityEt2>(Etmin));
}
Selector SelectorEtMax(double Etmax) {
  return Selector(new SW_QuantityMax<QuantityEt2>(Etmax));
}
Selector SelectorEtRange(double Etmin, double Etmax) {
  return Selector(new SW_QuantityRange<QuantityEt2>(Etmin, Etmax));
}
class QuantityE : public QuantityBase{
public:
  QuantityE(double E) : QuantityBase(E){}
  virtual double operator()(const PseudoJet & jet ) const { return jet.E();}
  virtual string description() const {return "E";}
};  
Selector SelectorEMin(double Emin) {
  return Selector(new SW_QuantityMin<QuantityE>(Emin));
}
Selector SelectorEMax(double Emax) {
  return Selector(new SW_QuantityMax<QuantityE>(Emax));
}
Selector SelectorERange(double Emin, double Emax) {
  return Selector(new SW_QuantityRange<QuantityE>(Emin, Emax));
}
class QuantityM2 : public QuantitySquareBase{
public:
  QuantityM2(double m) : QuantitySquareBase(m){}
  virtual double operator()(const PseudoJet & jet ) const { return jet.m2();}
  virtual string description() const {return "mass";}
};  
Selector SelectorMassMin(double mmin) {
  return Selector(new SW_QuantityMin<QuantityM2>(mmin));
}
Selector SelectorMassMax(double mmax) {
  return Selector(new SW_QuantityMax<QuantityM2>(mmax));
}
Selector SelectorMassRange(double mmin, double mmax) {
  return Selector(new SW_QuantityRange<QuantityM2>(mmin, mmax));
}
class QuantityRap : public QuantityBase{
public:
  QuantityRap(double rap) : QuantityBase(rap){}
  virtual double operator()(const PseudoJet & jet ) const { return jet.rap();}
  virtual string description() const {return "rap";}
  virtual bool is_geometric() const { return true;}
};  
class SW_RapMin : public SW_QuantityMin<QuantityRap>{
public:
  SW_RapMin(double rapmin) : SW_QuantityMin<QuantityRap>(rapmin){}
  virtual void get_rapidity_extent(double &rapmin, double & rapmax) const{
    rapmax = std::numeric_limits<double>::max();     
    rapmin = _qmin.comparison_value();
  }
};
class SW_RapMax : public SW_QuantityMax<QuantityRap>{
public:
  SW_RapMax(double rapmax) : SW_QuantityMax<QuantityRap>(rapmax){}
  virtual void get_rapidity_extent(double &rapmin, double & rapmax) const{
    rapmax = _qmax.comparison_value(); 
    rapmin = -std::numeric_limits<double>::max();
  }
};
class SW_RapRange : public SW_QuantityRange<QuantityRap>{
public:
  SW_RapRange(double rapmin, double rapmax) : SW_QuantityRange<QuantityRap>(rapmin, rapmax){
    assert(rapmin<=rapmax);
  }
  virtual void get_rapidity_extent(double &rapmin, double & rapmax) const{
    rapmax = _qmax.comparison_value();      
    rapmin = _qmin.comparison_value(); 
  }
  virtual bool has_known_area() const { return true;} ///< the area is analytically known
  virtual double known_area() const { 
    return twopi * (_qmax.comparison_value()-_qmin.comparison_value());
  }
};
Selector SelectorRapMin(double rapmin) {
  return Selector(new SW_RapMin(rapmin));
}
Selector SelectorRapMax(double rapmax) {
  return Selector(new SW_RapMax(rapmax));
}
Selector SelectorRapRange(double rapmin, double rapmax) {
  return Selector(new SW_RapRange(rapmin, rapmax));
}
class QuantityAbsRap : public QuantityBase{
public:
  QuantityAbsRap(double absrap) : QuantityBase(absrap){}
  virtual double operator()(const PseudoJet & jet ) const { return abs(jet.rap());}
  virtual string description() const {return "|rap|";}
  virtual bool is_geometric() const { return true;}
};  
class SW_AbsRapMax : public SW_QuantityMax<QuantityAbsRap>{
public:
  SW_AbsRapMax(double absrapmax) : SW_QuantityMax<QuantityAbsRap>(absrapmax){}
  virtual void get_rapidity_extent(double &rapmin, double & rapmax) const{
    rapmax =  _qmax.comparison_value(); 
    rapmin = -_qmax.comparison_value();
  }
  virtual bool has_known_area() const { return true;}   ///< the area is analytically known
  virtual double known_area() const { 
    return twopi * 2 * _qmax.comparison_value();
  }
};
class SW_AbsRapRange : public SW_QuantityRange<QuantityAbsRap>{
public:
  SW_AbsRapRange(double absrapmin, double absrapmax) : SW_QuantityRange<QuantityAbsRap>(absrapmin, absrapmax){}
  virtual void get_rapidity_extent(double &rapmin, double & rapmax) const{
    rapmax =  _qmax.comparison_value(); 
    rapmin = -_qmax.comparison_value();
  }
  virtual bool has_known_area() const { return true;} ///< the area is analytically known
  virtual double known_area() const { 
    return twopi * 2 * (_qmax.comparison_value()-max(_qmin.comparison_value(),0.0)); // this should handle properly absrapmin<0
  }
};
Selector SelectorAbsRapMin(double absrapmin) {
  return Selector(new SW_QuantityMin<QuantityAbsRap>(absrapmin));
}
Selector SelectorAbsRapMax(double absrapmax) {
  return Selector(new SW_AbsRapMax(absrapmax));
}
Selector SelectorAbsRapRange(double rapmin, double rapmax) {
  return Selector(new SW_AbsRapRange(rapmin, rapmax));
}
class QuantityEta : public QuantityBase{
public:
  QuantityEta(double eta) : QuantityBase(eta){}
  virtual double operator()(const PseudoJet & jet ) const { return jet.eta();}
  virtual string description() const {return "eta";}
};  
Selector SelectorEtaMin(double etamin) {
  return Selector(new SW_QuantityMin<QuantityEta>(etamin));
}
Selector SelectorEtaMax(double etamax) {
  return Selector(new SW_QuantityMax<QuantityEta>(etamax));
}
Selector SelectorEtaRange(double etamin, double etamax) {
  return Selector(new SW_QuantityRange<QuantityEta>(etamin, etamax));
}
class QuantityAbsEta : public QuantityBase{
public:
  QuantityAbsEta(double abseta) : QuantityBase(abseta){}
  virtual double operator()(const PseudoJet & jet ) const { return abs(jet.eta());}
  virtual string description() const {return "|eta|";}
  virtual bool is_geometric() const { return true;}
};  
Selector SelectorAbsEtaMin(double absetamin) {
  return Selector(new SW_QuantityMin<QuantityAbsEta>(absetamin));
}
Selector SelectorAbsEtaMax(double absetamax) {
  return Selector(new SW_QuantityMax<QuantityAbsEta>(absetamax));
}
Selector SelectorAbsEtaRange(double absetamin, double absetamax) {
  return Selector(new SW_QuantityRange<QuantityAbsEta>(absetamin, absetamax));
}
class SW_PhiRange : public SelectorWorker {
public:
  SW_PhiRange(double phimin, double phimax) : _phimin(phimin), _phimax(phimax){
    assert(_phimin<_phimax);
    assert(_phimin>-twopi);
    assert(_phimax<2*twopi);
    _phispan = _phimax - _phimin;
  }
  virtual bool pass(const PseudoJet & jet) const {
    double dphi=jet.phi()-_phimin;
    if (dphi >= twopi) dphi -= twopi;
    if (dphi < 0)      dphi += twopi;
    return (dphi <= _phispan);
  }
  virtual string description() const {
    ostringstream ostr;
    ostr << _phimin << " <= phi <= " << _phimax;
    return ostr.str();
  }
  virtual bool is_geometric() const { return true;}
protected:
  double _phimin;   // the lower cut 
  double _phimax;   // the upper cut
  double _phispan;  // the span of the range
};
Selector SelectorPhiRange(double phimin, double phimax) {
  return Selector(new SW_PhiRange(phimin, phimax));
}
class SW_RapPhiRange : public SW_And{
public:
  SW_RapPhiRange(double rapmin, double rapmax, double phimin, double phimax)
    : SW_And(SelectorRapRange(rapmin, rapmax), SelectorPhiRange(phimin, phimax)){
    _known_area = ((phimax-phimin > twopi) ? twopi : phimax-phimin) * (rapmax-rapmin);
  }
  virtual double known_area() const{
    return _known_area;
  }
protected:
  double _known_area;
};
Selector SelectorRapPhiRange(double rapmin, double rapmax, double phimin, double phimax) {
  return Selector(new SW_RapPhiRange(rapmin, rapmax, phimin, phimax));
}
class SW_NHardest : public SelectorWorker {
public:
  SW_NHardest(unsigned int n) : _n(n) {};
  virtual bool pass(const PseudoJet &) const {
    if (!applies_jet_by_jet())
      throw Error("Cannot apply this selector worker to an individual jet");
    return false;
  }
  virtual void terminator(vector<const PseudoJet *> & jets) const {
    if (jets.size() < _n) return;
    vector<double> minus_pt2(jets.size());
    vector<unsigned int> indices(jets.size());
    for (unsigned int i=0; i<jets.size(); i++){
      indices[i] = i;
      minus_pt2[i] = jets[i] ? -jets[i]->perp2() : 0.0;
    }
    IndexedSortHelper sort_helper(& minus_pt2);
    partial_sort(indices.begin(), indices.begin()+_n, indices.end(), sort_helper);
    for (unsigned int i=_n; i<jets.size(); i++)
      jets[indices[i]] = NULL;
  }
  virtual bool applies_jet_by_jet() const {return false;}
  virtual string description() const {
    ostringstream ostr;
    ostr << _n << " hardest";
    return ostr.str();
  }
protected:
  unsigned int _n;
};
Selector SelectorNHardest(unsigned int n) {
  return Selector(new SW_NHardest(n));
}
class SW_WithReference : public SelectorWorker{
public:
  SW_WithReference() : _is_initialised(false){};
  virtual bool takes_reference() const { return true;}
  virtual void set_reference(const PseudoJet &centre){
    _is_initialised = true;
    _reference = centre;
  }
protected:
  PseudoJet _reference;
  bool _is_initialised;
};
class SW_Circle : public SW_WithReference {
public:
  SW_Circle(const double radius) : _radius2(radius*radius) {}
  virtual SelectorWorker* copy(){ return new SW_Circle(*this);}
  virtual bool pass(const PseudoJet & jet) const {
    if (! _is_initialised)
      throw Error("To use a SelectorCircle (or any selector that requires a reference), you first have to call set_reference(...)");
    return jet.squared_distance(_reference) <= _radius2;
  } 
  virtual string description() const {
    ostringstream ostr;
    ostr << "distance from the centre <= " << sqrt(_radius2);
    return ostr.str();
  }
  virtual void get_rapidity_extent(double & rapmin, double & rapmax) const{
    if (! _is_initialised)
      throw Error("To use a SelectorCircle (or any selector that requires a reference), you first have to call set_reference(...)");
    rapmax = _reference.rap()+sqrt(_radius2);
    rapmin = _reference.rap()-sqrt(_radius2);
  }
  virtual bool is_geometric() const { return true;}    ///< implies a finite area
  virtual bool has_finite_area() const { return true;} ///< regardless of the reference 
  virtual bool has_known_area() const { return true;}  ///< the area is analytically known
  virtual double known_area() const { 
    return pi * _radius2;
  }
protected:
  double _radius2;
};
Selector SelectorCircle(const double radius) {
  return Selector(new SW_Circle(radius));
}
class SW_Doughnut : public SW_WithReference {
public:
  SW_Doughnut(const double radius_in, const double radius_out)
    : _radius_in2(radius_in*radius_in), _radius_out2(radius_out*radius_out) {}
  virtual SelectorWorker* copy(){ return new SW_Doughnut(*this);}
  virtual bool pass(const PseudoJet & jet) const {
    if (! _is_initialised)
      throw Error("To use a SelectorDoughnut (or any selector that requires a reference), you first have to call set_reference(...)");
    double distance2 = jet.squared_distance(_reference);
    return (distance2 <= _radius_out2) && (distance2 >= _radius_in2);
  } 
  virtual string description() const {
    ostringstream ostr;
    ostr << sqrt(_radius_in2) << " <= distance from the centre <= " << sqrt(_radius_out2);
    return ostr.str();
  }
  virtual void get_rapidity_extent(double & rapmin, double & rapmax) const{
    if (! _is_initialised)
      throw Error("To use a SelectorDoughnut (or any selector that requires a reference), you first have to call set_reference(...)");
    rapmax = _reference.rap()+sqrt(_radius_out2);
    rapmin = _reference.rap()-sqrt(_radius_out2);
  }
  virtual bool is_geometric() const { return true;}    ///< implies a finite area
  virtual bool has_finite_area() const { return true;} ///< regardless of the reference 
  virtual bool has_known_area() const { return true;}  ///< the area is analytically known
  virtual double known_area() const { 
    return pi * (_radius_out2-_radius_in2);
  }
protected:
  double _radius_in2, _radius_out2;
};
Selector SelectorDoughnut(const double radius_in, const double radius_out) {
  return Selector(new SW_Doughnut(radius_in, radius_out));
}
class SW_Strip : public SW_WithReference {
public:
  SW_Strip(const double delta) : _delta(delta) {}
  virtual SelectorWorker* copy(){ return new SW_Strip(*this);}
  virtual bool pass(const PseudoJet & jet) const {
    if (! _is_initialised)
      throw Error("To use a SelectorStrip (or any selector that requires a reference), you first have to call set_reference(...)");
    return abs(jet.rap()-_reference.rap()) <= _delta;
  } 
  virtual string description() const {
    ostringstream ostr;
    ostr << "|rap - rap_reference| <= " << _delta;
    return ostr.str();
  }
  virtual void get_rapidity_extent(double & rapmin, double & rapmax) const{
    if (! _is_initialised)
      throw Error("To use a SelectorStrip (or any selector that requires a reference), you first have to call set_reference(...)");
    rapmax = _reference.rap()+_delta;
    rapmin = _reference.rap()-_delta;
  }
  virtual bool is_geometric() const { return true;}    ///< implies a finite area
  virtual bool has_finite_area() const { return true;} ///< regardless of the reference 
  virtual bool has_known_area() const { return true;}  ///< the area is analytically known
  virtual double known_area() const { 
    return twopi * 2 * _delta;
  }
protected:
  double _delta;
};
Selector SelectorStrip(const double half_width) {
  return Selector(new SW_Strip(half_width));
}
class SW_Rectangle : public SW_WithReference {
public:
  SW_Rectangle(const double delta_rap, const double delta_phi)
    : _delta_rap(delta_rap),  _delta_phi(delta_phi) {}
  virtual SelectorWorker* copy(){ return new SW_Rectangle(*this);}
  virtual bool pass(const PseudoJet & jet) const {
    if (! _is_initialised)
      throw Error("To use a SelectorRectangle (or any selector that requires a reference), you first have to call set_reference(...)");
    return (abs(jet.rap()-_reference.rap()) <= _delta_rap) && (abs(jet.delta_phi_to(_reference)) <= _delta_phi);
  } 
  virtual string description() const {
    ostringstream ostr;
    ostr << "|rap - rap_reference| <= " << _delta_rap << " && |phi - phi_reference| <= " << _delta_phi ;
    return ostr.str();
  }
  virtual void get_rapidity_extent(double & rapmin, double & rapmax) const{
    if (! _is_initialised)
      throw Error("To use a SelectorRectangle (or any selector that requires a reference), you first have to call set_reference(...)");
    rapmax = _reference.rap()+_delta_rap;
    rapmin = _reference.rap()-_delta_rap;
  }
  virtual bool is_geometric() const { return true;}    ///< implies a finite area
  virtual bool has_finite_area() const { return true;} ///< regardless of the reference 
  virtual bool has_known_area() const { return true;}  ///< the area is analytically known
  virtual double known_area() const { 
    return 4 * _delta_rap * _delta_phi;
  }
protected:
  double _delta_rap, _delta_phi;
};
Selector SelectorRectangle(const double half_rap_width, const double half_phi_width) {
  return Selector(new SW_Rectangle(half_rap_width, half_phi_width));
}
class SW_PtFractionMin : public SW_WithReference {
public:
  SW_PtFractionMin(double fraction) : _fraction2(fraction*fraction){}
  virtual SelectorWorker* copy(){ return new SW_PtFractionMin(*this);}
  virtual bool pass(const PseudoJet & jet) const {
    if (! _is_initialised)
      throw Error("To use a SelectorPtFractionMin (or any selector that requires a reference), you first have to call set_reference(...)");
    return (jet.perp2() >= _fraction2*_reference.perp2());
  }
  virtual string description() const {
    ostringstream ostr;
    ostr << "pt >= " << sqrt(_fraction2) << "* pt_ref";
    return ostr.str();
  }
protected:
  double _fraction2;
};
Selector SelectorPtFractionMin(double fraction){
  return Selector(new SW_PtFractionMin(fraction));
}
class SW_IsZero : public SelectorWorker {
public:
  SW_IsZero(){}
  virtual bool pass(const PseudoJet & jet) const {
    return jet==0;
  }
  virtual string description() const { return "zero";}
};
Selector SelectorIsZero(){
  return Selector(new SW_IsZero());
}
Selector & Selector::operator &=(const Selector & b){
  _worker.reset(new SW_And(*this, b));
  return *this;
}
Selector & Selector::operator |=(const Selector & b){
  _worker.reset(new SW_Or(*this, b));
  return *this;
}
FJCORE_END_NAMESPACE      // defined in fastjet/internal/base.hh
#include <iomanip>
using namespace std;
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
LazyTiling25::LazyTiling25(ClusterSequence & cs) :
  _cs(cs), _jets(cs.jets())
{
#ifdef INSTRUMENT2
  _ncall = 0; // gps tmp
  _ncall_dtt = 0; // gps tmp
#endif // INSTRUMENT2
  _Rparam = cs.jet_def().R();
  _R2 = _Rparam * _Rparam;
  _invR2 = 1.0 / _R2;
  _initialise_tiles();
}
void LazyTiling25::_initialise_tiles() {
  double default_size = max(0.1,_Rparam)/2;
  _tile_size_eta = default_size;
  _n_tiles_phi   = max(5,int(floor(twopi/default_size)));
  _tile_size_phi = twopi / _n_tiles_phi; // >= _Rparam and fits in 2pi
#define _FJCORE_TILING25_USE_TILING_ANALYSIS_
#ifdef  _FASTJET_TILING25_USE_TILING_ANALYSIS_
  TilingExtent tiling_analysis(_cs);
  _tiles_eta_min = tiling_analysis.minrap();
  _tiles_eta_max = tiling_analysis.maxrap();
#else // not _FASTJET_TILING25_USE_TILING_ANALYSIS_
  _tiles_eta_min = 0.0;
  _tiles_eta_max = 0.0;
  const double maxrap = 7.0;
  for(unsigned int i = 0; i < _jets.size(); i++) {
    double eta = _jets[i].rap();
    if (abs(eta) < maxrap) {
      if (eta < _tiles_eta_min) {_tiles_eta_min = eta;}
      if (eta > _tiles_eta_max) {_tiles_eta_max = eta;}
    }
  }
#endif // _FASTJET_TILING25_USE_TILING_ANALYSIS_
# define FJCORE_LAZY25_MIN3TILESY
#ifdef FJCORE_LAZY25_MIN3TILESY
   if (_tiles_eta_max - _tiles_eta_min < 3*_tile_size_eta) {
     _tile_size_eta = (_tiles_eta_max - _tiles_eta_min)/3;
     _tiles_ieta_min = 0;
     _tiles_ieta_max = 2;
     _tiles_eta_max -= _tile_size_eta;
   } else {
#endif //FASTJET_LAZY25_MIN3TILESY
    _tiles_ieta_min = int(floor(_tiles_eta_min/_tile_size_eta));
    _tiles_ieta_max = int(floor( _tiles_eta_max/_tile_size_eta));
    _tiles_eta_min = _tiles_ieta_min * _tile_size_eta;
    _tiles_eta_max = _tiles_ieta_max * _tile_size_eta;
#ifdef FJCORE_LAZY25_MIN3TILESY
   }
#endif
  _tile_half_size_eta = _tile_size_eta * 0.5;
  _tile_half_size_phi = _tile_size_phi * 0.5;
  vector<bool> use_periodic_delta_phi(_n_tiles_phi, false);
  if (_n_tiles_phi <= 5) {
    fill(use_periodic_delta_phi.begin(), use_periodic_delta_phi.end(), true);
  } else {
    use_periodic_delta_phi[0] = true;
    use_periodic_delta_phi[1] = true;
    use_periodic_delta_phi[_n_tiles_phi-2] = true;
    use_periodic_delta_phi[_n_tiles_phi-1] = true;
  }
  _tiles.resize((_tiles_ieta_max-_tiles_ieta_min+1)*_n_tiles_phi);
  for (int ieta = _tiles_ieta_min; ieta <= _tiles_ieta_max; ieta++) {
    for (int iphi = 0; iphi < _n_tiles_phi; iphi++) {
      Tile25 * tile = & _tiles[_tile_index(ieta,iphi)];
      tile->head = NULL; // first element of tiles points to itself
      tile->begin_tiles[0] =  tile;
      Tile25 ** pptile = & (tile->begin_tiles[0]);
      pptile++;
      tile->surrounding_tiles = pptile;
      if (ieta > _tiles_ieta_min) {
	// with the itile subroutine, we can safely run tiles from
	// idphi=-1 to idphi=+1, because it takes care of
	// negative and positive boundaries
	for (int idphi = -2; idphi <=+2; idphi++) {
	  *pptile = & _tiles[_tile_index(ieta-1,iphi+idphi)];
	  pptile++;
	}	
      }
      if (ieta > _tiles_ieta_min + 1) {
	// with the itile subroutine, we can safely run tiles from
	// idphi=-1 to idphi=+1, because it takes care of
	// negative and positive boundaries
	for (int idphi = -2; idphi <= +2; idphi++) {
	  *pptile = & _tiles[_tile_index(ieta-2,iphi+idphi)];
	  pptile++;
	}	
      }
      *pptile = & _tiles[_tile_index(ieta,iphi-1)];
      pptile++;
      *pptile = & _tiles[_tile_index(ieta,iphi-2)];
      pptile++;
      tile->RH_tiles = pptile;
      *pptile = & _tiles[_tile_index(ieta,iphi+1)];
      pptile++;
      *pptile = & _tiles[_tile_index(ieta,iphi+2)];
      pptile++;
      if (ieta < _tiles_ieta_max) {
	for (int idphi = -2; idphi <= +2; idphi++) {
	  *pptile = & _tiles[_tile_index(ieta+1,iphi+idphi)];
	  pptile++;
	}	
      }
      if (ieta < _tiles_ieta_max - 1) {
	for (int idphi = -2; idphi <= +2; idphi++) {
	  *pptile = & _tiles[_tile_index(ieta+2,iphi+idphi)];
	  pptile++;
	}	
      }
      tile->end_tiles = pptile;
      tile->tagged = false;
      tile->use_periodic_delta_phi = use_periodic_delta_phi[iphi];
      tile->max_NN_dist = 0;
      tile->eta_centre = (ieta-_tiles_ieta_min+0.5)*_tile_size_eta + _tiles_eta_min;
      tile->phi_centre = (iphi+0.5)*_tile_size_phi;
    }
  }
}
int LazyTiling25::_tile_index(const double eta, const double phi) const {
  int ieta, iphi;
  if      (eta <= _tiles_eta_min) {ieta = 0;}
  else if (eta >= _tiles_eta_max) {ieta = _tiles_ieta_max-_tiles_ieta_min;}
  else {
    ieta = int(((eta - _tiles_eta_min) / _tile_size_eta));
    if (ieta > _tiles_ieta_max-_tiles_ieta_min) {
      ieta = _tiles_ieta_max-_tiles_ieta_min;} 
  }
  iphi = int((phi+twopi)/_tile_size_phi) % _n_tiles_phi;
  return (iphi + ieta * _n_tiles_phi);
}
inline void LazyTiling25::_tj_set_jetinfo( TiledJet * const jet,
					      const int _jets_index) {
  _bj_set_jetinfo<>(jet, _jets_index);
  jet->tile_index = _tile_index(jet->eta, jet->phi);
  Tile25 * tile = &_tiles[jet->tile_index];
  jet->previous   = NULL;
  jet->next       = tile->head;
  if (jet->next != NULL) {jet->next->previous = jet;}
  tile->head      = jet;
}
void LazyTiling25::_bj_remove_from_tiles(TiledJet * const jet) {
  Tile25 * tile = & _tiles[jet->tile_index];
  if (jet->previous == NULL) {
    tile->head = jet->next;
  } else {
    jet->previous->next = jet->next;
  }
  if (jet->next != NULL) {
    jet->next->previous = jet->previous;
  }
}
void LazyTiling25::_print_tiles(TiledJet * briefjets ) const {
  for (vector<Tile25>::const_iterator tile = _tiles.begin(); 
       tile < _tiles.end(); tile++) {
    cout << "Tile " << tile - _tiles.begin()
         << " at " << setw(10) << tile->eta_centre << "," << setw(10) << tile->phi_centre
         << " = ";
    vector<int> list;
    for (TiledJet * jetI = tile->head; jetI != NULL; jetI = jetI->next) {
      list.push_back(jetI-briefjets);
    }
    sort(list.begin(),list.end());
    for (unsigned int i = 0; i < list.size(); i++) {cout <<" "<<list[i];}
    cout <<"\n";
  }
}
void LazyTiling25::_add_neighbours_to_tile_union(const int tile_index, 
	       vector<int> & tile_union, int & n_near_tiles) const {
  for (Tile25 * const * near_tile = _tiles[tile_index].begin_tiles; 
       near_tile != _tiles[tile_index].end_tiles; near_tile++){
    tile_union[n_near_tiles] = *near_tile - & _tiles[0];
    n_near_tiles++;
  }
}
inline void LazyTiling25::_add_untagged_neighbours_to_tile_union(
               const int tile_index, 
	       vector<int> & tile_union, int & n_near_tiles)  {
  for (Tile25 ** near_tile = _tiles[tile_index].begin_tiles; 
       near_tile != _tiles[tile_index].end_tiles; near_tile++){
    if (! (*near_tile)->tagged) {
      (*near_tile)->tagged = true;
      tile_union[n_near_tiles] = *near_tile - & _tiles[0];
      n_near_tiles++;
    }
  }
}
inline void LazyTiling25::_add_untagged_neighbours_to_tile_union_using_max_info(
               const TiledJet * jet, 
	       vector<int> & tile_union, int & n_near_tiles)  {
  Tile25 & tile = _tiles[jet->tile_index];
  for (Tile25 ** near_tile = tile.begin_tiles; near_tile != tile.end_tiles; near_tile++){
    if ((*near_tile)->tagged) continue;
    double dist = _distance_to_tile(jet, *near_tile) - tile_edge_security_margin;
    if (dist > (*near_tile)->max_NN_dist) continue;
    (*near_tile)->tagged = true;
    tile_union[n_near_tiles] = *near_tile - & _tiles[0];
    n_near_tiles++;
  }
}
inline double LazyTiling25::_distance_to_tile(const TiledJet * bj, const Tile25 * tile) 
#ifdef INSTRUMENT2
   {
  _ncall_dtt++; // GPS tmp
#else
  const {
#endif // INSTRUMENT2
  double deta;
  if (_tiles[bj->tile_index].eta_centre == tile->eta_centre) deta = 0;
  else   deta = std::abs(bj->eta - tile->eta_centre) - _tile_half_size_eta;
  double dphi = std::abs(bj->phi - tile->phi_centre);
  if (dphi > pi) dphi = twopi-dphi;
  dphi -= _tile_half_size_phi;
  if (dphi < 0) dphi = 0;
  return dphi*dphi + deta*deta;
}
inline void LazyTiling25::_update_jetX_jetI_NN(TiledJet * jetX, TiledJet * jetI, vector<TiledJet *> & jets_for_minheap) {
  double dist = _bj_dist(jetI,jetX);
  if (dist < jetI->NN_dist) {
    if (jetI != jetX) {
      jetI->NN_dist = dist;
      jetI->NN = jetX;
      if (!jetI->minheap_update_needed()) {
	jetI->label_minheap_update_needed();
	jets_for_minheap.push_back(jetI);
      }
    }
  }
  if (dist < jetX->NN_dist) {
    if (jetI != jetX) {
      jetX->NN_dist = dist;
      jetX->NN      = jetI;}
  }
}
inline void LazyTiling25::_set_NN(TiledJet * jetI, 
                              vector<TiledJet *> & jets_for_minheap) {
  jetI->NN_dist = _R2;
  jetI->NN      = NULL;
  if (!jetI->minheap_update_needed()) {
    jetI->label_minheap_update_needed();
    jets_for_minheap.push_back(jetI);}
  Tile25 * tile_ptr = &_tiles[jetI->tile_index];
    for (Tile25 ** near_tile  = tile_ptr->begin_tiles; 
         near_tile != tile_ptr->end_tiles; near_tile++) {
      if (jetI->NN_dist < _distance_to_tile(jetI, *near_tile)) continue;
      for (TiledJet * jetJ  = (*near_tile)->head; 
           jetJ != NULL; jetJ = jetJ->next) {
        double dist = _bj_dist(jetI,jetJ);
        if (dist < jetI->NN_dist && jetJ != jetI) {
          jetI->NN_dist = dist; jetI->NN = jetJ;
        }
      }
    }
}
void LazyTiling25::run() {
  int n = _jets.size();
  if (n == 0) return; 
  TiledJet * briefjets = new TiledJet[n];
  TiledJet * jetA = briefjets, * jetB;
  TiledJet oldB = briefjets[0]; 
  vector<int> tile_union(3*25);
  for (int i = 0; i< n; i++) {
    _tj_set_jetinfo(jetA, i);
    jetA++; // move on to next entry of briefjets
  }
  TiledJet * head = briefjets; // a nicer way of naming start
  vector<Tile25>::iterator tile;
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      for (jetB = tile->head; jetB != jetA; jetB = jetB->next) {
	double dist = _bj_dist_not_periodic(jetA,jetB);
	if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
      }
    }
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      if (jetA->NN_dist > tile->max_NN_dist) tile->max_NN_dist = jetA->NN_dist;
    }
  }
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    if (tile->use_periodic_delta_phi) {
      for (Tile25 ** RTile = tile->RH_tiles; RTile != tile->end_tiles; RTile++) {
        for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
          double dist_to_tile = _distance_to_tile(jetA, *RTile);
          bool relevant_for_jetA  = dist_to_tile <= jetA->NN_dist;
          bool relevant_for_RTile = dist_to_tile <= (*RTile)->max_NN_dist;
          if (relevant_for_jetA || relevant_for_RTile) {
            for (jetB = (*RTile)->head; jetB != NULL; jetB = jetB->next) {
              double dist = _bj_dist(jetA,jetB);
              if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
              if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
            }
          } 
        }
      }
    } else {
      for (Tile25 ** RTile = tile->RH_tiles; RTile != tile->end_tiles; RTile++) {
        for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
          double dist_to_tile = _distance_to_tile(jetA, *RTile);
          bool relevant_for_jetA  = dist_to_tile <= jetA->NN_dist;
          bool relevant_for_RTile = dist_to_tile <= (*RTile)->max_NN_dist;
          if (relevant_for_jetA || relevant_for_RTile) {
            for (jetB = (*RTile)->head; jetB != NULL; jetB = jetB->next) {
              double dist = _bj_dist_not_periodic(jetA,jetB);
              if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
              if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
            }
          } 
        }
      }
    }
  }
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    tile->max_NN_dist = 0;
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      if (jetA->NN_dist > tile->max_NN_dist) tile->max_NN_dist = jetA->NN_dist;
    }
  }
#ifdef INSTRUMENT2
  cout << "intermediate ncall, dtt = " << _ncall << " " << _ncall_dtt << endl; // GPS tmp
#endif // INSTRUMENT2
  vector<double> diJs(n);
  for (int i = 0; i < n; i++) {
    diJs[i] = _bj_diJ(&briefjets[i]);
    briefjets[i].label_minheap_update_done();
  }
  MinHeap minheap(diJs);
  vector<TiledJet *> jets_for_minheap;
  jets_for_minheap.reserve(n); 
  int history_location = n-1;
  while (n > 0) {
    double diJ_min = minheap.minval() *_invR2;
    jetA = head + minheap.minloc();
    history_location++;
    jetB = jetA->NN;
    if (jetB != NULL) {
      if (jetA < jetB) {std::swap(jetA,jetB);}
      int nn; // new jet index
      _cs.plugin_record_ij_recombination(jetA->_jets_index, jetB->_jets_index, diJ_min, nn);
      _bj_remove_from_tiles(jetA);
      oldB = * jetB;  // take a copy because we will need it...
      _bj_remove_from_tiles(jetB);
      _tj_set_jetinfo(jetB, nn); // cause jetB to become _jets[nn]
    } else {
      _cs.plugin_record_iB_recombination(jetA->_jets_index, diJ_min);
      _bj_remove_from_tiles(jetA);
    }
    minheap.remove(jetA-head);
    int n_near_tiles = 0;
    if (jetB != NULL) {
      Tile25 & jetB_tile = _tiles[jetB->tile_index];
      for (Tile25 ** near_tile  = jetB_tile.begin_tiles; 
	           near_tile != jetB_tile.end_tiles; near_tile++) {
    	double dist_to_tile = _distance_to_tile(jetB, *near_tile);
    	bool relevant_for_jetB  = dist_to_tile <= jetB->NN_dist;
    	bool relevant_for_near_tile = dist_to_tile <= (*near_tile)->max_NN_dist;
        bool relevant = relevant_for_jetB || relevant_for_near_tile;
        if (! relevant) continue;
        tile_union[n_near_tiles] = *near_tile - & _tiles[0];
        (*near_tile)->tagged = true;
        n_near_tiles++;
        for (TiledJet * jetI = (*near_tile)->head; jetI != NULL; jetI = jetI->next) {
          if (jetI->NN == jetA || jetI->NN == jetB) _set_NN(jetI, jets_for_minheap);
          _update_jetX_jetI_NN(jetB, jetI, jets_for_minheap);
        }
      }
    }
    int n_done_tiles = n_near_tiles;
    _add_untagged_neighbours_to_tile_union_using_max_info(jetA, 
       					   tile_union, n_near_tiles);
    if (jetB != NULL) {
	_add_untagged_neighbours_to_tile_union_using_max_info(&oldB,
							      tile_union,n_near_tiles);
      jetB->label_minheap_update_needed();
      jets_for_minheap.push_back(jetB);
    }
    for (int itile = 0; itile < n_done_tiles; itile++) {
      _tiles[tile_union[itile]].tagged = false;
    }
    for (int itile = n_done_tiles; itile < n_near_tiles; itile++) {
      Tile25 * tile_ptr = &_tiles[tile_union[itile]];
      tile_ptr->tagged = false;
      for (TiledJet * jetI = tile_ptr->head; jetI != NULL; jetI = jetI->next) {
        if (jetI->NN == jetA || (jetI->NN == jetB && jetB != NULL)) {
          _set_NN(jetI, jets_for_minheap);
        }
      }
    }
    while (jets_for_minheap.size() > 0) {
      TiledJet * jetI = jets_for_minheap.back(); 
      jets_for_minheap.pop_back();
      minheap.update(jetI-head, _bj_diJ(jetI));
      jetI->label_minheap_update_done();
      Tile25 & tile_I = _tiles[jetI->tile_index];
      if (tile_I.max_NN_dist < jetI->NN_dist) tile_I.max_NN_dist = jetI->NN_dist;
    }
    n--;
  }
  delete[] briefjets;
#ifdef INSTRUMENT2
  cout << "ncall, dtt = " << _ncall << " " << _ncall_dtt << endl; // GPS tmp
#endif // INSTRUMENT2
}
FJCORE_END_NAMESPACE
#include <iomanip>
#include <limits>
#include <cmath>
using namespace std;
#define _FJCORE_TILING2_USE_TILING_ANALYSIS_
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
LazyTiling9::LazyTiling9(ClusterSequence & cs) :
  _cs(cs), _jets(cs.jets())
{
#ifdef INSTRUMENT2
  _ncall = 0; // gps tmp
  _ncall_dtt = 0; // gps tmp
#endif // INSTRUMENT2
  _Rparam = cs.jet_def().R();
  _R2 = _Rparam * _Rparam;
  _invR2 = 1.0 / _R2;
  _initialise_tiles();
}
void LazyTiling9::_initialise_tiles() {
  double default_size = max(0.1,_Rparam);
  _tile_size_eta = default_size;
  _n_tiles_phi   = max(3,int(floor(twopi/default_size)));
  _tile_size_phi = twopi / _n_tiles_phi; // >= _Rparam and fits in 2pi
#ifdef _FJCORE_TILING2_USE_TILING_ANALYSIS_
  TilingExtent tiling_analysis(_cs);
  _tiles_eta_min = tiling_analysis.minrap();
  _tiles_eta_max = tiling_analysis.maxrap();
#else
  _tiles_eta_min = 0.0;
  _tiles_eta_max = 0.0;
  const double maxrap = 7.0;
  for(unsigned int i = 0; i < _jets.size(); i++) {
    double eta = _jets[i].rap();
    if (abs(eta) < maxrap) {
      if (eta < _tiles_eta_min) {_tiles_eta_min = eta;}
      if (eta > _tiles_eta_max) {_tiles_eta_max = eta;}
    }
  }
#endif
# define FJCORE_LAZY9_MIN2TILESY
#ifdef FJCORE_LAZY9_MIN2TILESY
   if (_tiles_eta_max - _tiles_eta_min < 2*_tile_size_eta) {
     _tile_size_eta = (_tiles_eta_max - _tiles_eta_min)/2;
     _tiles_ieta_min = 0;
     _tiles_ieta_max = 1;
     _tiles_eta_max -= _tile_size_eta;
   } else {
#endif //FASTJET_LAZY9_MIN2TILESY
  _tiles_ieta_min = int(floor(_tiles_eta_min/_tile_size_eta));
  _tiles_ieta_max = int(floor( _tiles_eta_max/_tile_size_eta));
  _tiles_eta_min = _tiles_ieta_min * _tile_size_eta;
  _tiles_eta_max = _tiles_ieta_max * _tile_size_eta;
#ifdef FJCORE_LAZY9_MIN2TILESY
   }
#endif
  _tile_half_size_eta = _tile_size_eta * 0.5;
  _tile_half_size_phi = _tile_size_phi * 0.5;
  vector<bool> use_periodic_delta_phi(_n_tiles_phi, false);
  if (_n_tiles_phi <= 3) {
    fill(use_periodic_delta_phi.begin(), use_periodic_delta_phi.end(), true);
  } else {
    use_periodic_delta_phi[0] = true;
    use_periodic_delta_phi[_n_tiles_phi-1] = true;
  }
  _tiles.resize((_tiles_ieta_max-_tiles_ieta_min+1)*_n_tiles_phi);
  for (int ieta = _tiles_ieta_min; ieta <= _tiles_ieta_max; ieta++) {
    for (int iphi = 0; iphi < _n_tiles_phi; iphi++) {
      Tile2 * tile = & _tiles[_tile_index(ieta,iphi)];
      tile->head = NULL; // first element of tiles points to itself
      tile->begin_tiles[0] =  tile;
      Tile2 ** pptile = & (tile->begin_tiles[0]);
      pptile++;
      tile->surrounding_tiles = pptile;
      if (ieta > _tiles_ieta_min) {
	// with the itile subroutine, we can safely run tiles from
	// idphi=-1 to idphi=+1, because it takes care of
	// negative and positive boundaries
	for (int idphi = -1; idphi <=+1; idphi++) {
	  *pptile = & _tiles[_tile_index(ieta-1,iphi+idphi)];
	  pptile++;
	}	
      }
      *pptile = & _tiles[_tile_index(ieta,iphi-1)];
      pptile++;
      tile->RH_tiles = pptile;
      *pptile = & _tiles[_tile_index(ieta,iphi+1)];
      pptile++;
      if (ieta < _tiles_ieta_max) {
	for (int idphi = -1; idphi <= +1; idphi++) {
	  *pptile = & _tiles[_tile_index(ieta+1,iphi+idphi)];
	  pptile++;
	}	
      }
      tile->end_tiles = pptile;
      tile->tagged = false;
      tile->use_periodic_delta_phi = use_periodic_delta_phi[iphi];
      tile->max_NN_dist = 0;
      tile->eta_centre = (ieta-_tiles_ieta_min+0.5)*_tile_size_eta + _tiles_eta_min;
      tile->phi_centre = (iphi+0.5)*_tile_size_phi;
    }
  }
}
int LazyTiling9::_tile_index(const double eta, const double phi) const {
  int ieta, iphi;
  if      (eta <= _tiles_eta_min) {ieta = 0;}
  else if (eta >= _tiles_eta_max) {ieta = _tiles_ieta_max-_tiles_ieta_min;}
  else {
    ieta = int(((eta - _tiles_eta_min) / _tile_size_eta));
    if (ieta > _tiles_ieta_max-_tiles_ieta_min) {
      ieta = _tiles_ieta_max-_tiles_ieta_min;} 
  }
  iphi = int((phi+twopi)/_tile_size_phi) % _n_tiles_phi;
  return (iphi + ieta * _n_tiles_phi);
}
inline void LazyTiling9::_tj_set_jetinfo( TiledJet * const jet,
					      const int _jets_index) {
  _bj_set_jetinfo<>(jet, _jets_index);
  jet->tile_index = _tile_index(jet->eta, jet->phi);
  Tile2 * tile = &_tiles[jet->tile_index];
  jet->previous   = NULL;
  jet->next       = tile->head;
  if (jet->next != NULL) {jet->next->previous = jet;}
  tile->head      = jet;
}
void LazyTiling9::_bj_remove_from_tiles(TiledJet * const jet) {
  Tile2 * tile = & _tiles[jet->tile_index];
  if (jet->previous == NULL) {
    tile->head = jet->next;
  } else {
    jet->previous->next = jet->next;
  }
  if (jet->next != NULL) {
    jet->next->previous = jet->previous;
  }
}
void LazyTiling9::_print_tiles(TiledJet * briefjets ) const {
  for (vector<Tile2>::const_iterator tile = _tiles.begin(); 
       tile < _tiles.end(); tile++) {
    cout << "Tile " << tile - _tiles.begin()<<" = ";
    vector<int> list;
    for (TiledJet * jetI = tile->head; jetI != NULL; jetI = jetI->next) {
      list.push_back(jetI-briefjets);
    }
    sort(list.begin(),list.end());
    for (unsigned int i = 0; i < list.size(); i++) {cout <<" "<<list[i];}
    cout <<"\n";
  }
}
void LazyTiling9::_add_neighbours_to_tile_union(const int tile_index, 
	       vector<int> & tile_union, int & n_near_tiles) const {
  for (Tile2 * const * near_tile = _tiles[tile_index].begin_tiles; 
       near_tile != _tiles[tile_index].end_tiles; near_tile++){
    tile_union[n_near_tiles] = *near_tile - & _tiles[0];
    n_near_tiles++;
  }
}
inline void LazyTiling9::_add_untagged_neighbours_to_tile_union(
               const int tile_index, 
	       vector<int> & tile_union, int & n_near_tiles)  {
  for (Tile2 ** near_tile = _tiles[tile_index].begin_tiles; 
       near_tile != _tiles[tile_index].end_tiles; near_tile++){
    if (! (*near_tile)->tagged) {
      (*near_tile)->tagged = true;
      tile_union[n_near_tiles] = *near_tile - & _tiles[0];
      n_near_tiles++;
    }
  }
}
inline void LazyTiling9::_add_untagged_neighbours_to_tile_union_using_max_info(
               const TiledJet * jet, 
	       vector<int> & tile_union, int & n_near_tiles)  {
  Tile2 & tile = _tiles[jet->tile_index];
  for (Tile2 ** near_tile = tile.begin_tiles; near_tile != tile.end_tiles; near_tile++){
    if ((*near_tile)->tagged) continue;
    double dist = _distance_to_tile(jet, *near_tile) - tile_edge_security_margin;
    if (dist > (*near_tile)->max_NN_dist) continue;
    (*near_tile)->tagged = true;
    tile_union[n_near_tiles] = *near_tile - & _tiles[0];
    n_near_tiles++;
  }
}
inline double LazyTiling9::_distance_to_tile(const TiledJet * bj, const Tile2 * tile) 
#ifdef INSTRUMENT2
   {
  _ncall_dtt++; // GPS tmp
#else
  const {
#endif // INSTRUMENT2
  double deta;
  if (_tiles[bj->tile_index].eta_centre == tile->eta_centre) deta = 0;
  else   deta = std::abs(bj->eta - tile->eta_centre) - _tile_half_size_eta;
  double dphi = std::abs(bj->phi - tile->phi_centre);
  if (dphi > pi) dphi = twopi-dphi;
  dphi -= _tile_half_size_phi;
  if (dphi < 0) dphi = 0;
  return dphi*dphi + deta*deta;
}
inline void LazyTiling9::_update_jetX_jetI_NN(TiledJet * jetX, TiledJet * jetI, vector<TiledJet *> & jets_for_minheap) {
  double dist = _bj_dist(jetI,jetX);
  if (dist < jetI->NN_dist) {
    if (jetI != jetX) {
      jetI->NN_dist = dist;
      jetI->NN = jetX;
      if (!jetI->minheap_update_needed()) {
	jetI->label_minheap_update_needed();
	jets_for_minheap.push_back(jetI);
      }
    }
  }
  if (dist < jetX->NN_dist) {
    if (jetI != jetX) {
      jetX->NN_dist = dist;
      jetX->NN      = jetI;}
  }
}
inline void LazyTiling9::_set_NN(TiledJet * jetI, 
                              vector<TiledJet *> & jets_for_minheap) {
  jetI->NN_dist = _R2;
  jetI->NN      = NULL;
  if (!jetI->minheap_update_needed()) {
    jetI->label_minheap_update_needed();
    jets_for_minheap.push_back(jetI);}
  Tile2 * tile_ptr = &_tiles[jetI->tile_index];
    for (Tile2 ** near_tile  = tile_ptr->begin_tiles; 
         near_tile != tile_ptr->end_tiles; near_tile++) {
      if (jetI->NN_dist < _distance_to_tile(jetI, *near_tile)) continue;
      for (TiledJet * jetJ  = (*near_tile)->head; 
           jetJ != NULL; jetJ = jetJ->next) {
        double dist = _bj_dist(jetI,jetJ);
        if (dist < jetI->NN_dist && jetJ != jetI) {
          jetI->NN_dist = dist; jetI->NN = jetJ;
        }
      }
    }
}
void LazyTiling9::run() {
  int n = _jets.size();
  if (n == 0) return; 
  TiledJet * briefjets = new TiledJet[n];
  TiledJet * jetA = briefjets, * jetB;
  TiledJet oldB = briefjets[0]; 
  vector<int> tile_union(3*n_tile_neighbours);
  for (int i = 0; i< n; i++) {
    _tj_set_jetinfo(jetA, i);
    jetA++; // move on to next entry of briefjets
  }
  TiledJet * head = briefjets; // a nicer way of naming start
  vector<Tile2>::iterator tile;
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      for (jetB = tile->head; jetB != jetA; jetB = jetB->next) {
	double dist = _bj_dist_not_periodic(jetA,jetB);
	if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
      }
    }
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      if (jetA->NN_dist > tile->max_NN_dist) tile->max_NN_dist = jetA->NN_dist;
    }
  }
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    if (tile->use_periodic_delta_phi) {
      for (Tile2 ** RTile = tile->RH_tiles; RTile != tile->end_tiles; RTile++) {
        for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
          double dist_to_tile = _distance_to_tile(jetA, *RTile);
          bool relevant_for_jetA  = dist_to_tile <= jetA->NN_dist;
          bool relevant_for_RTile = dist_to_tile <= (*RTile)->max_NN_dist;
          if (relevant_for_jetA || relevant_for_RTile) {
            for (jetB = (*RTile)->head; jetB != NULL; jetB = jetB->next) {
              double dist = _bj_dist(jetA,jetB);
              if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
              if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
            }
          } 
        }
      }
    } else {
      for (Tile2 ** RTile = tile->RH_tiles; RTile != tile->end_tiles; RTile++) {
        for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
          double dist_to_tile = _distance_to_tile(jetA, *RTile);
          bool relevant_for_jetA  = dist_to_tile <= jetA->NN_dist;
          bool relevant_for_RTile = dist_to_tile <= (*RTile)->max_NN_dist;
          if (relevant_for_jetA || relevant_for_RTile) {
            for (jetB = (*RTile)->head; jetB != NULL; jetB = jetB->next) {
              double dist = _bj_dist_not_periodic(jetA,jetB);
              if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
              if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
            }
          } 
        }
      }
    }
  }
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    tile->max_NN_dist = 0;
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      if (jetA->NN_dist > tile->max_NN_dist) tile->max_NN_dist = jetA->NN_dist;
    }
  }
#ifdef INSTRUMENT2
  cout << "intermediate ncall, dtt = " << _ncall << " " << _ncall_dtt << endl; // GPS tmp
#endif // INSTRUMENT2
  vector<double> diJs(n);
  for (int i = 0; i < n; i++) {
    diJs[i] = _bj_diJ(&briefjets[i]);
    briefjets[i].label_minheap_update_done();
  }
  MinHeap minheap(diJs);
  vector<TiledJet *> jets_for_minheap;
  jets_for_minheap.reserve(n); 
  int history_location = n-1;
  while (n > 0) {
    double diJ_min = minheap.minval() *_invR2;
    jetA = head + minheap.minloc();
    history_location++;
    jetB = jetA->NN;
    if (jetB != NULL) {
      if (jetA < jetB) {std::swap(jetA,jetB);}
      int nn; // new jet index
      _cs.plugin_record_ij_recombination(jetA->_jets_index, jetB->_jets_index, diJ_min, nn);
      _bj_remove_from_tiles(jetA);
      oldB = * jetB;  // take a copy because we will need it...
      _bj_remove_from_tiles(jetB);
      _tj_set_jetinfo(jetB, nn); // cause jetB to become _jets[nn]
    } else {
      _cs.plugin_record_iB_recombination(jetA->_jets_index, diJ_min);
      _bj_remove_from_tiles(jetA);
    }
    minheap.remove(jetA-head);
    int n_near_tiles = 0;
    if (jetB != NULL) {
      Tile2 & jetB_tile = _tiles[jetB->tile_index];
      for (Tile2 ** near_tile  = jetB_tile.begin_tiles; 
	           near_tile != jetB_tile.end_tiles; near_tile++) {
    	double dist_to_tile = _distance_to_tile(jetB, *near_tile);
    	bool relevant_for_jetB  = dist_to_tile <= jetB->NN_dist;
    	bool relevant_for_near_tile = dist_to_tile <= (*near_tile)->max_NN_dist;
        bool relevant = relevant_for_jetB || relevant_for_near_tile;
        if (! relevant) continue;
        tile_union[n_near_tiles] = *near_tile - & _tiles[0];
        (*near_tile)->tagged = true;
        n_near_tiles++;
        for (TiledJet * jetI = (*near_tile)->head; jetI != NULL; jetI = jetI->next) {
          if (jetI->NN == jetA || jetI->NN == jetB) _set_NN(jetI, jets_for_minheap);
          _update_jetX_jetI_NN(jetB, jetI, jets_for_minheap);
        }
      }
    }
    int n_done_tiles = n_near_tiles;
    _add_untagged_neighbours_to_tile_union_using_max_info(jetA, 
       					   tile_union, n_near_tiles);
    if (jetB != NULL) {
	_add_untagged_neighbours_to_tile_union_using_max_info(&oldB,
							      tile_union,n_near_tiles);
      jetB->label_minheap_update_needed();
      jets_for_minheap.push_back(jetB);
    }
    for (int itile = 0; itile < n_done_tiles; itile++) {
      _tiles[tile_union[itile]].tagged = false;
    }
    for (int itile = n_done_tiles; itile < n_near_tiles; itile++) {
      Tile2 * tile_ptr = &_tiles[tile_union[itile]];
      tile_ptr->tagged = false;
      for (TiledJet * jetI = tile_ptr->head; jetI != NULL; jetI = jetI->next) {
        if (jetI->NN == jetA || (jetI->NN == jetB && jetB != NULL)) {
          _set_NN(jetI, jets_for_minheap);
        }
      }
    }
    while (jets_for_minheap.size() > 0) {
      TiledJet * jetI = jets_for_minheap.back(); 
      jets_for_minheap.pop_back();
      minheap.update(jetI-head, _bj_diJ(jetI));
      jetI->label_minheap_update_done();
      Tile2 & tile_I = _tiles[jetI->tile_index];
      if (tile_I.max_NN_dist < jetI->NN_dist) tile_I.max_NN_dist = jetI->NN_dist;
    }
    n--;
  }
  delete[] briefjets;
#ifdef INSTRUMENT2
  cout << "ncall, dtt = " << _ncall << " " << _ncall_dtt << endl; // GPS tmp
#endif // INSTRUMENT2
}
FJCORE_END_NAMESPACE
#include <iomanip>
using namespace std;
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
LazyTiling9Alt::LazyTiling9Alt(ClusterSequence & cs) :
  _cs(cs), _jets(cs.jets())
{
  _Rparam = cs.jet_def().R();
  _R2 = _Rparam * _Rparam;
  _invR2 = 1.0 / _R2;
  _initialise_tiles();
}
void LazyTiling9Alt::_initialise_tiles() {
  double default_size = max(0.1,_Rparam);
  _tile_size_eta = default_size;
  _n_tiles_phi   = max(3,int(floor(twopi/default_size)));
  _tile_size_phi = twopi / _n_tiles_phi; // >= _Rparam and fits in 2pi
  _tiles_eta_min = 0.0;
  _tiles_eta_max = 0.0;
  const double maxrap = 7.0;
  for(unsigned int i = 0; i < _jets.size(); i++) {
    double eta = _jets[i].rap();
    if (abs(eta) < maxrap) {
      if (eta < _tiles_eta_min) {_tiles_eta_min = eta;}
      if (eta > _tiles_eta_max) {_tiles_eta_max = eta;}
    }
  }
  _tiles_ieta_min = int(floor(_tiles_eta_min/_tile_size_eta));
  _tiles_ieta_max = int(floor( _tiles_eta_max/_tile_size_eta));
  _tiles_eta_min = _tiles_ieta_min * _tile_size_eta;
  _tiles_eta_max = _tiles_ieta_max * _tile_size_eta;
  _tile_half_size_eta = _tile_size_eta * 0.5;
  _tile_half_size_phi = _tile_size_phi * 0.5;
  vector<bool> use_periodic_delta_phi(_n_tiles_phi, false);
  if (_n_tiles_phi <= 3) {
    fill(use_periodic_delta_phi.begin(), use_periodic_delta_phi.end(), true);
  } else {
    use_periodic_delta_phi[0] = true;
    use_periodic_delta_phi[_n_tiles_phi-1] = true;
  }
  _tiles.resize((_tiles_ieta_max-_tiles_ieta_min+1)*_n_tiles_phi);
  for (int ieta = _tiles_ieta_min; ieta <= _tiles_ieta_max; ieta++) {
    for (int iphi = 0; iphi < _n_tiles_phi; iphi++) {
      Tile * tile = & _tiles[_tile_index(ieta,iphi)];
      tile->head = NULL; // first element of tiles points to itself
      tile->begin_tiles[0] =  Tile::TileFnPair(tile,&Tile::distance_to_centre);
      Tile::TileFnPair * pptile = & (tile->begin_tiles[0]);
      pptile++;
      tile->surrounding_tiles = pptile;
      if (ieta > _tiles_ieta_min) {
	// with the itile subroutine, we can safely run tiles from
	// idphi=-1 to idphi=+1, because it takes care of
	// negative and positive boundaries
	//for (int idphi = -1; idphi <=+1; idphi++) {
        *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta-1,iphi-1)],
                                   &Tile::distance_to_left_bottom);
        pptile++;
        *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta-1,iphi)],
                                   &Tile::distance_to_left);
        pptile++;
        *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta-1,iphi+1)],
                                   &Tile::distance_to_left_top);
        pptile++;
      }
      *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta,iphi-1)], 
                                 &Tile::distance_to_bottom);
      pptile++;
      tile->RH_tiles = pptile;
      *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta,iphi+1)], 
                                 &Tile::distance_to_top);
      pptile++;
      if (ieta < _tiles_ieta_max) {
	//for (int idphi = -1; idphi <= +1; idphi++) {
	//  *pptile = & _tiles[_tile_index(ieta+1,iphi+idphi)];
	//  pptile++;
	//}	
        *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta+1,iphi-1)],
                                   &Tile::distance_to_right_bottom);
        pptile++;
        *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta+1,iphi)],
                                   &Tile::distance_to_right);
        pptile++;
        *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta+1,iphi+1)],
                                   &Tile::distance_to_right_top);
        pptile++;
      }
      tile->end_tiles = pptile;
      tile->tagged = false;
      tile->use_periodic_delta_phi = use_periodic_delta_phi[iphi];
      tile->max_NN_dist = 0;
      tile->eta_min = ieta*_tile_size_eta;
      tile->eta_max = (ieta+1)*_tile_size_eta;
      tile->phi_min = iphi*_tile_size_phi;
      tile->phi_max = (iphi+1)*_tile_size_phi;
    }
  }
}
int LazyTiling9Alt::_tile_index(const double eta, const double phi) const {
  int ieta, iphi;
  if      (eta <= _tiles_eta_min) {ieta = 0;}
  else if (eta >= _tiles_eta_max) {ieta = _tiles_ieta_max-_tiles_ieta_min;}
  else {
    ieta = int(((eta - _tiles_eta_min) / _tile_size_eta));
    if (ieta > _tiles_ieta_max-_tiles_ieta_min) {
      ieta = _tiles_ieta_max-_tiles_ieta_min;} 
  }
  iphi = int((phi+twopi)/_tile_size_phi) % _n_tiles_phi;
  return (iphi + ieta * _n_tiles_phi);
}
inline void LazyTiling9Alt::_tj_set_jetinfo( TiledJet * const jet,
					      const int _jets_index) {
  _bj_set_jetinfo<>(jet, _jets_index);
  jet->tile_index = _tile_index(jet->eta, jet->phi);
  Tile * tile = &_tiles[jet->tile_index];
  jet->previous   = NULL;
  jet->next       = tile->head;
  if (jet->next != NULL) {jet->next->previous = jet;}
  tile->head      = jet;
}
void LazyTiling9Alt::_bj_remove_from_tiles(TiledJet * const jet) {
  Tile * tile = & _tiles[jet->tile_index];
  if (jet->previous == NULL) {
    tile->head = jet->next;
  } else {
    jet->previous->next = jet->next;
  }
  if (jet->next != NULL) {
    jet->next->previous = jet->previous;
  }
}
void LazyTiling9Alt::_print_tiles(TiledJet * briefjets ) const {
  for (vector<Tile>::const_iterator tile = _tiles.begin(); 
       tile < _tiles.end(); tile++) {
    cout << "Tile " << tile - _tiles.begin()<<" = ";
    vector<int> list;
    for (TiledJet * jetI = tile->head; jetI != NULL; jetI = jetI->next) {
      list.push_back(jetI-briefjets);
    }
    sort(list.begin(),list.end());
    for (unsigned int i = 0; i < list.size(); i++) {cout <<" "<<list[i];}
    cout <<"\n";
  }
}
void LazyTiling9Alt::_add_neighbours_to_tile_union(const int tile_index, 
	       vector<int> & tile_union, int & n_near_tiles) const {
  for (Tile::TileFnPair const * near_tile = _tiles[tile_index].begin_tiles; 
       near_tile != _tiles[tile_index].end_tiles; near_tile++){
    tile_union[n_near_tiles] = near_tile->first - & _tiles[0];
    n_near_tiles++;
  }
}
inline void LazyTiling9Alt::_add_untagged_neighbours_to_tile_union(
               const int tile_index, 
	       vector<int> & tile_union, int & n_near_tiles)  {
  for (Tile::TileFnPair * near_tile = _tiles[tile_index].begin_tiles; 
       near_tile != _tiles[tile_index].end_tiles; near_tile++){
    if (! (near_tile->first)->tagged) {
      (near_tile->first)->tagged = true;
      tile_union[n_near_tiles] = near_tile->first - & _tiles[0];
      n_near_tiles++;
    }
  }
}
inline void LazyTiling9Alt::_add_untagged_neighbours_to_tile_union_using_max_info(
               const TiledJet * jet, 
	       vector<int> & tile_union, int & n_near_tiles)  {
  Tile & tile = _tiles[jet->tile_index];
  for (Tile::TileFnPair * near_tile = tile.begin_tiles; near_tile != tile.end_tiles; near_tile++){
    if ((near_tile->first)->tagged) continue;
    double dist = (tile.*(near_tile->second))(jet) - tile_edge_security_margin;
    if (dist > (near_tile->first)->max_NN_dist) continue;
    (near_tile->first)->tagged = true;
    tile_union[n_near_tiles] = near_tile->first - & _tiles[0];
    n_near_tiles++;
  }
}
ostream & operator<<(ostream & ostr, const TiledJet & jet) {
  ostr << "j" << setw(3) << jet._jets_index << ":pt2,rap,phi=" ; ostr.flush();
  ostr     << jet.kt2 << ","; ostr.flush();
  ostr     << jet.eta << ","; ostr.flush();
  ostr     << jet.phi; ostr.flush();
  ostr     << ", tile=" << jet.tile_index; ostr.flush();
  return ostr;
}
inline void LazyTiling9Alt::_update_jetX_jetI_NN(TiledJet * jetX, TiledJet * jetI, vector<TiledJet *> & jets_for_minheap) {
  double dist = _bj_dist(jetI,jetX);
  if (dist < jetI->NN_dist) {
    if (jetI != jetX) {
      jetI->NN_dist = dist;
      jetI->NN = jetX;
      if (!jetI->minheap_update_needed()) {
	jetI->label_minheap_update_needed();
	jets_for_minheap.push_back(jetI);
      }
    }
  }
  if (dist < jetX->NN_dist) {
    if (jetI != jetX) {
      jetX->NN_dist = dist;
      jetX->NN      = jetI;}
  }
}
inline void LazyTiling9Alt::_set_NN(TiledJet * jetI, 
                            vector<TiledJet *> & jets_for_minheap) {
  jetI->NN_dist = _R2;
  jetI->NN      = NULL;
  if (!jetI->minheap_update_needed()) {
    jetI->label_minheap_update_needed();
    jets_for_minheap.push_back(jetI);}
  Tile * tile_ptr = &_tiles[jetI->tile_index];
    for (Tile::TileFnPair * near_tile  = tile_ptr->begin_tiles; 
         near_tile != tile_ptr->end_tiles; near_tile++) {
      if (jetI->NN_dist < (tile_ptr->*(near_tile->second))(jetI)) continue;
      for (TiledJet * jetJ  = (near_tile->first)->head; 
           jetJ != NULL; jetJ = jetJ->next) {
        double dist = _bj_dist(jetI,jetJ);
        if (dist < jetI->NN_dist && jetJ != jetI) {
          jetI->NN_dist = dist; jetI->NN = jetJ;
        }
      }
    }
}
void LazyTiling9Alt::run() {
  int n = _jets.size();
  TiledJet * briefjets = new TiledJet[n];
  TiledJet * jetA = briefjets, * jetB;
  TiledJet oldB;
  vector<int> tile_union(3*n_tile_neighbours);
  for (int i = 0; i< n; i++) {
    _tj_set_jetinfo(jetA, i);
    jetA++; // move on to next entry of briefjets
  }
  TiledJet * head = briefjets; // a nicer way of naming start
  vector<Tile>::iterator tile;
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      for (jetB = tile->head; jetB != jetA; jetB = jetB->next) {
	double dist = _bj_dist_not_periodic(jetA,jetB);
	if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
      }
    }
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      if (jetA->NN_dist > tile->max_NN_dist) tile->max_NN_dist = jetA->NN_dist;
    }
  }
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    if (tile->use_periodic_delta_phi) {
      for (Tile::TileFnPair * RTileFnPair = tile->RH_tiles; 
           RTileFnPair != tile->end_tiles; RTileFnPair++) {
        Tile *RTile = RTileFnPair->first;
        for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
          double dist_to_tile = ((*tile).*(RTileFnPair->second))(jetA);
          bool relevant_for_jetA  = dist_to_tile <= jetA->NN_dist;
          bool relevant_for_RTile = dist_to_tile <= RTile->max_NN_dist;
          if (relevant_for_jetA || relevant_for_RTile) {
            for (jetB = RTile->head; jetB != NULL; jetB = jetB->next) {
              double dist = _bj_dist(jetA,jetB);
              if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
              if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
            }
          } 
        }
      }
    } else {
      for (Tile::TileFnPair* RTileFnPair = tile->RH_tiles;
           RTileFnPair != tile->end_tiles; RTileFnPair++) {
        Tile *RTile = RTileFnPair->first;
        for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
          double dist_to_tile = ((*tile).*(RTileFnPair->second))(jetA);
          bool relevant_for_jetA  = dist_to_tile <= jetA->NN_dist;
          bool relevant_for_RTile = dist_to_tile <= RTile->max_NN_dist;
          if (relevant_for_jetA || relevant_for_RTile) {
            for (jetB = RTile->head; jetB != NULL; jetB = jetB->next) {
              double dist = _bj_dist_not_periodic(jetA,jetB);
              if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
              if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
            }
          } 
        }
      }
    }
  }
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    tile->max_NN_dist = 0;
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      if (jetA->NN_dist > tile->max_NN_dist) tile->max_NN_dist = jetA->NN_dist;
    }
  }
  vector<double> diJs(n);
  for (int i = 0; i < n; i++) {
    diJs[i] = _bj_diJ(&briefjets[i]);
    briefjets[i].label_minheap_update_done();
  }
  MinHeap minheap(diJs);
  vector<TiledJet *> jets_for_minheap;
  jets_for_minheap.reserve(n); 
  int history_location = n-1;
  while (n > 0) {
    double diJ_min = minheap.minval() *_invR2;
    jetA = head + minheap.minloc();
    history_location++;
    jetB = jetA->NN;
    if (jetB != NULL) {
      if (jetA < jetB) {std::swap(jetA,jetB);}
      int nn; // new jet index
      _cs.plugin_record_ij_recombination(jetA->_jets_index, jetB->_jets_index, diJ_min, nn);
      _bj_remove_from_tiles(jetA);
      oldB = * jetB;  // take a copy because we will need it...
      _bj_remove_from_tiles(jetB);
      _tj_set_jetinfo(jetB, nn); // cause jetB to become _jets[nn]
    } else {
      _cs.plugin_record_iB_recombination(jetA->_jets_index, diJ_min);
      _bj_remove_from_tiles(jetA);
    }
    minheap.remove(jetA-head);
    int n_near_tiles = 0;
    _add_untagged_neighbours_to_tile_union_using_max_info(jetA, 
       					   tile_union, n_near_tiles);
    if (jetB != NULL) {
	_add_untagged_neighbours_to_tile_union_using_max_info(&oldB,
							      tile_union,n_near_tiles);
      jetB->label_minheap_update_needed();
      jets_for_minheap.push_back(jetB);
    }
    if (jetB != NULL) {
      Tile & jetB_tile = _tiles[jetB->tile_index];
      for (Tile::TileFnPair * near_tile_fn_pair  = jetB_tile.begin_tiles; 
	           near_tile_fn_pair != jetB_tile.end_tiles; near_tile_fn_pair++) {
        Tile * near_tile = near_tile_fn_pair->first;
    	double dist_to_tile = (jetB_tile.*(near_tile_fn_pair->second))(jetB);
    	bool relevant_for_jetB  = dist_to_tile <= jetB->NN_dist;
    	bool relevant_for_near_tile = dist_to_tile <= near_tile->max_NN_dist;
        bool relevant = relevant_for_jetB || relevant_for_near_tile;
        if (relevant) {
          if (near_tile->tagged) {
            for (TiledJet * jetI = near_tile->head; jetI != NULL; jetI = jetI->next) {
              if (jetI->NN == jetA || jetI->NN == jetB) _set_NN(jetI, jets_for_minheap);
              _update_jetX_jetI_NN(jetB, jetI, jets_for_minheap);
            }
          near_tile->tagged = false;
          } else {
            for (TiledJet * jetI = near_tile->head; jetI != NULL; jetI = jetI->next) {
              _update_jetX_jetI_NN(jetB, jetI, jets_for_minheap);
            }
          }
        }
	//     // -- Keep this old inline code for later speed tests
      }
    }
    for (int itile = 0; itile < n_near_tiles; itile++) {
      Tile * tile_ptr = &_tiles[tile_union[itile]];
      if (!tile_ptr->tagged) continue; // because earlier loop may have undone the tag
      tile_ptr->tagged = false;
      for (TiledJet * jetI = tile_ptr->head; jetI != NULL; jetI = jetI->next) {
        if (jetI->NN == jetA || (jetI->NN == jetB && jetB != NULL)) {
          _set_NN(jetI, jets_for_minheap);
        }
      }
    }
    while (jets_for_minheap.size() > 0) {
      TiledJet * jetI = jets_for_minheap.back(); 
      jets_for_minheap.pop_back();
      minheap.update(jetI-head, _bj_diJ(jetI));
      jetI->label_minheap_update_done();
      Tile & tile_I = _tiles[jetI->tile_index];
      if (tile_I.max_NN_dist < jetI->NN_dist) tile_I.max_NN_dist = jetI->NN_dist;
    }
    n--;
  }
  delete[] briefjets;
}
FJCORE_END_NAMESPACE
#include <iomanip>
#include <limits>
#include <cmath>
using namespace std;
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
TilingExtent::TilingExtent(ClusterSequence & cs) {
  _determine_rapidity_extent(cs.jets());
}
TilingExtent::TilingExtent(const vector<PseudoJet> &particles) {
  _determine_rapidity_extent(particles);
}
void TilingExtent::_determine_rapidity_extent(const vector<PseudoJet> & particles) {
  int nrap = 20; 
  int nbins = 2*nrap;
  vector<double> counts(nbins, 0);
  _minrap =  numeric_limits<double>::max();
  _maxrap = -numeric_limits<double>::max();
  int ibin;
  for (unsigned i = 0; i < particles.size(); i++) {
    if (particles[i].E() == abs(particles[i].pz())) continue;
    double rap = particles[i].rap();
    if (rap < _minrap) _minrap = rap;
    if (rap > _maxrap) _maxrap = rap;
    ibin = int(rap+nrap); 
    if (ibin < 0) ibin = 0;
    if (ibin >= nbins) ibin = nbins - 1;
    counts[ibin]++;
  }
  double max_in_bin = 0;
  for (ibin = 0; ibin < nbins; ibin++) {
    if (max_in_bin < counts[ibin]) max_in_bin = counts[ibin];
  }
  const double allowed_max_fraction = 0.25;
  const double min_multiplicity = 4;
  double allowed_max_cumul = floor(max(max_in_bin * allowed_max_fraction, min_multiplicity));
  if (allowed_max_cumul > max_in_bin) allowed_max_cumul = max_in_bin;
  double cumul_lo = 0;
  _cumul2 = 0;
  for (ibin = 0; ibin < nbins; ibin++) {
    cumul_lo += counts[ibin];
    if (cumul_lo >= allowed_max_cumul) {
      double y = ibin-nrap;
      if (y > _minrap) _minrap = y;
      break;
    }
  }
  assert(ibin != nbins); // internal consistency check that you found a bin
  _cumul2 += cumul_lo*cumul_lo;
  int ibin_lo = ibin;
  double cumul_hi = 0;
  for (ibin = nbins-1; ibin >= 0; ibin--) {
    cumul_hi += counts[ibin];
    if (cumul_hi >= allowed_max_cumul) {
      double y = ibin-nrap+1; // +1 here is the rapidity bin width
      if (y < _maxrap) _maxrap = y;
      break;
    }
  }
  assert(ibin >= 0); // internal consistency check that you found a bin
  int ibin_hi = ibin;
  assert(ibin_hi >= ibin_lo); 
  if (ibin_hi == ibin_lo) {
    _cumul2 = pow(double(cumul_lo + cumul_hi - counts[ibin_hi]), 2);
  } else {
    _cumul2 += cumul_hi*cumul_hi;
    for (ibin = ibin_lo+1; ibin < ibin_hi; ibin++) {
      _cumul2 += counts[ibin]*counts[ibin];
    }
  }
}
FJCORE_END_NAMESPACE
