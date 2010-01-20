#ifndef __EVTMTREE_HH__
#define __EVTMTREE_HH__

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtSpinAmp.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include "EvtGenBase/EvtMNode.hh"
#include "EvtGenBase/EvtMParticle.hh"
#include "EvtGenBase/EvtMRes.hh"

#include <vector>
using std::vector;

#include <string>
using std::string;

typedef string::const_iterator ptype;

class EvtParticle;

class EvtMTree {
    
    public:

        EvtMTree( const EvtId * , unsigned int ); 
        ~EvtMTree( );

        // return the invariant amplitude of the entire tree
        EvtSpinAmp amplitude( EvtParticle * ) const;

        // add a decay tree to the list of trees that we posess
        void addtree( const string& );

    private:
        
        vector< EvtMNode * > _root;
        vector<string> _lbltbl;
        double _norm;

        bool parsecheck( char , const string& );
        void parseerror( bool, ptype&, ptype&, ptype& );
        
        string parseId( ptype&, ptype&, ptype& );
        string parseKey( ptype&, ptype&, ptype& );
        vector<string> parseArg( ptype&, ptype&, ptype& );
        vector<EvtComplex> parseAmps( ptype&, ptype&, ptype& );
        vector<EvtMNode *> duplicate( const vector<EvtMNode *>& ) const;
        vector<vector<EvtMNode * > > unionChildren( const string&,
                vector<vector<EvtMNode * > >& );
        vector<vector<EvtMNode * > > parseChildren( ptype&, ptype&, ptype& );
        vector<EvtMNode *> parsenode( const string& , bool );
        bool validTree( const EvtMNode * ) const;
        
        vector<EvtMNode *> makeparticles( const string& );
        EvtMRes * makeresonance( const EvtId&, const string &, const
                vector<string>&, const string& , const vector<EvtComplex>& ,
                const vector<EvtMNode * >& ); 

        EvtSpinAmp getrotation( EvtParticle * ) const;
};

#endif
