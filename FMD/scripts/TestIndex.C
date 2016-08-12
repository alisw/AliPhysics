//____________________________________________________________________
//
// $Id$
//
/** @file    TestIndex.C
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Thu Mar 30 01:20:02 2006
    @brief   Test AliFMDIndex and AliFMDObjIndex
*/
#ifndef __CINT__
#include <AliFMDIndex.h>
#include <iostream>
#include <TFile.h>
#endif
/** @defgroup FMD_index_test Test of AliFMDIndex and AliFMDObjIndex 
    @ingroup FMD_script
*/
/** Write an AliFMDIndex object to output stream
    @ingroup FMD_index_test
    @param o Stream
    @param i Index object 
    @return @a o */
std::ostream&
operator << (std::ostream& o, const AliFMDIndex& i) 
{
  return o << i.Name();
}

/** Do the comparision, and print to standard out
    @ingroup FMD_index_test
    @param lhs Left hand side
    @param rhs Right hand side
    @return @f$ lhs < rhs @f$ */
bool 
cmp(const AliFMDIndex& lhs, const AliFMDIndex& rhs)
{
  bool ret = lhs < rhs;
  std::cout << "    (" << lhs << " <  " << rhs << "): " << ret << std::endl;
  return ret;
}

/** Test if two index objects are equivilant (are the same).
    Equivilance is defined as 
    @f[
    lhs \equiv rhs: \not (lhs < rhs \vee rhs < lhs)
    @f]
    @ingroup FMD_index_test
    @param lhs Left hand side 
    @param rhs Right hand side 
    @return @c true if @a lhs and @a rhs are equivilant */
bool 
equiv(const AliFMDIndex& lhs, const AliFMDIndex& rhs)
{
  bool ret = !(cmp(lhs,rhs) || cmp(rhs,lhs));
  std::cout << "    (" << lhs << " <> " << rhs << "): " << ret << std::endl;
  return ret;
}

/** Check that @f$ \not (x < x)@f$
    @ingroup FMD_index_test
    @param x Object to test
    @return @c true if @a x is not less than itself */
bool
self(const AliFMDIndex& x)
{
  bool ret = !cmp(x,x);
  std::cout << "  !(" << x << " < " << x << "): " << ret << std::endl;
  return ret;
}

/** Check if @f$ a \wedge b \Rightarrow c@f$
    @ingroup FMD_index_test
    @param a First condition
    @param b Second condition
    @param c Implication
    @return @c true if the implication is valid */
bool
imply(bool a, bool b, bool c) 
{
  bool ret = ((a && b) && c) || (!(a && b));
  return ret;
}

/** Check if the comparison operator is transitive, that is 
    @f[
    (x < y \wedge y < z) \Rightarrow x < z
    @f]
    @ingroup FMD_index_test
    @param x First object
    @param y Second object
    @param z Third object
    @return @c true if the implication is met. */
bool
trans(const AliFMDIndex& x, const AliFMDIndex& y, const AliFMDIndex& z)
{
  bool ret  = imply(cmp(x,y), cmp(y,z), cmp(x,z));
  std::cout << "  (" << x << " < " << y << " && " << y << " < " << z 
	    << ") => " << x << " < " << z << " " 
	    << (ret ? "holds" : "violated") << std::endl;
  return ret;
  
}

/** Check that the comparison operator preserves equivilance, that is 
    @f[ 
    (x \equiv y \wedge y \equiv z) \Rightarrow (x \equiv z)
    @f]
    @ingroup FMD_index_test
    @param x First object
    @param y Second argument
    @param z Third object
    @return @c true if the implication holds. */
bool
equiv(const AliFMDIndex& x, const AliFMDIndex& y, const AliFMDIndex& z)
{
  bool ret  = imply(equiv(x,y), equiv(y,z), equiv(x,z));
  std::cout << "  (" << x << " <> " << y << " && " << y << " <> " << z 
	    << ") => " << x << " <> " << z << " " 
	    << (ret ? "holds" : "violated") << std::endl;
  return ret;
}


/** Check if the comparison operator is a @e strictly @e weak @e ordering
    @ingroup FMD_index_test
 */
void
TestIndex() 
{
  AliFMDIndex i1(1, 'I', 5, 63);
  AliFMDIndex i2(i1);  
  AliFMDIndex i3(i1);  
  AliFMDIndex i4(1, 'I', 5, 127);
  AliFMDIndex i5(1, 'I', 15, 63);
  AliFMDIndex i6(2, 'O', 15, 60);
  std::cout << "Is !(x < x): " << std::endl;
  self(i1);
  std::cout << "Does (x < y && y < z) imply x < z: " << std::endl;
  // true, true
  trans(i1, i4, i5);
  // true, false
  trans(i1, i6, i5);
  // false, true
  trans(i4, i1, i5);
  // false, false
  trans(i6, i5, i1);
  std::cout << "Does !(x < y || x > y) && !(y < z || z > y) imply "
	    << "!(x < z || x > z)" << std::endl;
  // true, true
  equiv(i1, i2, i3);
  // true, false
  equiv(i1, i2, i4);
  // false, true
  equiv(i4, i1, i2);
  // false, false
  equiv(i1, i4, i5);
  


  TFile* file = TFile::Open("index.root", "RECREATE");
  file->WriteObject(&i1,"i1");
  file->Write();
  file->Close();
  
  AliFMDIndex* i7 = 0;
  file = TFile::Open("index.root", "READ");
  file->GetObject("i1", i7);
  file->Close();
  std::cout << *i7 << " == " << i1 << ": " << (*i7 == i1) << std::endl;
  
}

/** Check if the comparison operator is a @e strictly @e weak @e ordering
    @ingroup FMD_index_test
 */
void
TestObjIndex() 
{
  AliFMDObjIndex i1(1, 'I', 5, 63);
  AliFMDObjIndex i2(i1);  
  AliFMDObjIndex i3(i1);  
  AliFMDObjIndex i4(1, 'I', 5, 127);
  AliFMDObjIndex i5(1, 'I', 15, 63);
  AliFMDObjIndex i6(2, 'O', 15, 60);
  std::cout << "Is !(x < x): " << std::endl;
  self(i1);
  std::cout << "Does (x < y && y < z) imply x < z: " << std::endl;
  // true, true
  trans(i1, i4, i5);
  // true, false
  trans(i1, i6, i5);
  // false, true
  trans(i4, i1, i5);
  // false, false
  trans(i6, i5, i1);
  std::cout << "Does !(x < y || x > y) && !(y < z || z > y) imply "
	    << "!(x < z || x > z)" << std::endl;
  // true, true
  equiv(i1, i2, i3);
  // true, false
  equiv(i1, i2, i4);
  // false, true
  equiv(i4, i1, i2);
  // false, false
  equiv(i1, i4, i5);
  


  TFile* file = TFile::Open("index.root", "RECREATE");
  i1.Write("i1");
  file->Write();
  file->Close();
  
  file = TFile::Open("index.root", "READ");
  AliFMDObjIndex* i7 = (AliFMDObjIndex*)(file->Get("i1"));
  file->Close();
  std::cout << *i7 << " == " << i1 << ": " << (*i7 == i1) << std::endl;
  
}

/** Check that we can sort an array of index objects
    @ingroup FMD_index_test
 */
void
SortIndex()
{
  TList l;
  for (size_t i = 0; i < 30; i++) {
    UShort_t det  = gRandom->Integer(3)+1;
    Char_t   ring = (gRandom->Uniform() > .5 ? 'O' : 'I');
    UShort_t sec  = gRandom->Integer(ring == 'I' ?  20 :  40);
    UShort_t str  = gRandom->Integer(ring == 'I' ? 512 : 256);
    l.AddAt(new AliFMDObjIndex(det, ring, sec, str), i);
  }
  std::cout << "Before sort" << std::endl;
  l.ls();
  l.Sort();
  std::cout << "After sort" << std::endl;
  l.ls();
}

    
//____________________________________________________________________
//
// EOF
//
