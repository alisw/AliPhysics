#ifndef __EVTSPINAMP_HH__
#define __EVTSPINAMP_HH__

#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtComplex.hh"

#include <vector>
using std::vector;

#include <cstdarg>

class EvtSpinAmp;
EvtSpinAmp operator*( const EvtComplex&, const EvtSpinAmp& );
EvtSpinAmp operator*( const EvtSpinAmp&, const EvtComplex& );
EvtSpinAmp operator/( const EvtSpinAmp&, const EvtComplex& );

class EvtSpinAmp
{

    friend EvtSpinAmp operator*( const EvtComplex&, const EvtSpinAmp& );
    friend EvtSpinAmp operator*( const EvtSpinAmp&, const EvtComplex& );
    friend EvtSpinAmp operator/( const EvtSpinAmp&, const EvtComplex& );
    friend std::ostream& operator<< ( std::ostream&, const EvtSpinAmp& );
    
public:

    EvtSpinAmp( ) {};
    EvtSpinAmp( const vector<EvtSpinType::spintype>& );
    EvtSpinAmp( const vector<EvtSpinType::spintype>& , const EvtComplex& );
    EvtSpinAmp( const vector<EvtSpinType::spintype>& , const vector<EvtComplex>& );
    EvtSpinAmp( const EvtSpinAmp & );

    ~EvtSpinAmp( ) {};
    
    // Input to the index functions are twice the magnetic quantum number
    EvtComplex& operator()( const vector<int>& ); 
    const EvtComplex& operator()( const vector<int>& ) const; 
    EvtComplex& operator()( int, ... );
    const EvtComplex& operator()( int, ... ) const;

    EvtSpinAmp& operator=( const EvtSpinAmp& );

    EvtSpinAmp operator+( const EvtSpinAmp& ) const;
    EvtSpinAmp& operator+=( const EvtSpinAmp& );

    EvtSpinAmp operator-( const EvtSpinAmp& ) const;
    EvtSpinAmp& operator-=( const EvtSpinAmp& );

    // Direct Product
    EvtSpinAmp operator*( const EvtSpinAmp& ) const;
    EvtSpinAmp& operator*=( const EvtSpinAmp& );
    
    EvtSpinAmp& operator*=( const EvtComplex& );
    EvtSpinAmp& operator/=( const EvtComplex& );

    // Contraction of amplitudes
    void intcont( size_t, size_t );
    void extcont( const EvtSpinAmp &, int, int );
    
    // assign this value to every member in the container
    void assign( const EvtComplex & val ) { _elem.assign( _elem.size(), val ); }
    
    // get the order of the container
    size_t rank( ) const { return _twospin.size(); }

    // get the dimension vector of the container
    const vector<unsigned int> & dims( ) const { return _twospin; } 
    
    // set the elements and the dimensions of the vector - useful for something
    // things eventough it is usually not the cleanest solution
    void addspin( int twospin ) { _twospin.push_back( twospin ); }
    void setelem( const vector<EvtComplex> &elem ) { _elem = elem; }

    bool iterate( vector<int>& index ) const;
    vector<int> iterinit() const;

    bool allowed( const vector<int>& index ) const;
    bool iterateallowed( vector<int>& index) const;
    vector<int> iterallowedinit() const;

private:

    void checkindexargs( const vector<int>& index ) const;
    void checktwospin( const vector<unsigned int>& twospin ) const;
    int findtrueindex( const vector<int>& index ) const;
    vector<unsigned int> calctwospin( const vector<EvtSpinType::spintype>& type ) const;

    vector<EvtSpinType::spintype> _type;
    vector<unsigned int> _twospin;
    vector<EvtComplex> _elem;
};

#endif // __EVTSPINAMP__
