#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "EvtGenBase/EvtParticle.hh"

#include "EvtGenBase/EvtMTree.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtReport.hh"

// Make sure to include Lineshapes here
#include "EvtGenBase/EvtMTrivialLS.hh"
#include "EvtGenBase/EvtMBreitWigner.hh"

// Make sure to include Parametrizations
#include "EvtGenBase/EvtMHelAmp.hh"

using std::endl;

EvtMTree::EvtMTree( const EvtId * idtbl, unsigned int ndaug )
{
    for( size_t i=0; i<ndaug; ++i ) {
        _lbltbl.push_back( EvtPDL::name( idtbl[i] ) );
    }
}

EvtMTree::~EvtMTree()
{
    for(size_t i=0; i<_root.size(); ++i) delete _root[i];
}

bool EvtMTree::parsecheck( char arg, const string& chars )
{
    bool ret = false;

    for(size_t i=0; i<chars.size(); ++i) {
        ret = ret || (chars[i]==arg);
    }

    return ret;
}

vector<EvtMNode *> EvtMTree::makeparticles( const string& strid ) 
{
    vector<EvtMNode *> particles;
    vector<int> labels;
   
    for( size_t i = 0; i<_lbltbl.size(); ++i ) {
        if( _lbltbl[i] == strid ) labels.push_back( i );
    }
    
    if( labels.size() == 0 ) {
        report(Severity::Error,"EvtGen")<<"Error unknown particle label "<<strid<<endl;
        ::abort();
    }

    for( size_t i = 0; i<labels.size(); ++i )
        particles.push_back( new EvtMParticle( labels[i], EvtPDL::getId( strid ) ) );

    return particles;
}

EvtMRes * EvtMTree::makeresonance( const EvtId& id, const string& ls,
        const vector<string>& lsarg, const string& type,
        const vector<EvtComplex>& amps, const vector<EvtMNode *>& children )
{
    EvtMRes * resonance = NULL;
    EvtMLineShape * lineshape = NULL;

    if( ls=="BREITWIGNER" ) {
        lineshape = new EvtMBreitWigner( id, lsarg );
    } else if( ls=="TRIVIAL" ) {
        lineshape = new EvtMTrivialLS( id, lsarg );
    } else {
        report(Severity::Error,"EvtGen")<<"Lineshape "<<lineshape
                              <<" not recognized."<<endl;
        ::abort();
    }

    if( type=="HELAMP" ) {
        resonance = new EvtMHelAmp( id, lineshape, children, amps );
    } else {
        report(Severity::Error,"EvtGen")<<"Model "<<type<<" not recognized."<<endl;
        ::abort();
    }

    lineshape->setres( resonance );

    return resonance;
}

void EvtMTree::parseerror( bool flag, ptype& c_iter, ptype& c_begin, 
        ptype& c_end )
{ 
    if(!flag) return;

    string error;
    
    while( c_begin != c_end ) {
        if(c_begin == c_iter) {
            error+='_';
            error+=*c_begin;
            error+='_';
        } else 
            error+=*c_begin;

        ++c_begin;
    }

    report(Severity::Error,"EvtGen")<<"Parse error at: "<<error<<endl;
    ::abort();
}

string EvtMTree::parseId( ptype& c_iter, ptype& c_begin, ptype& c_end ) 
{
    string strid;

    while(c_iter != c_end) {
        parseerror(parsecheck(*c_iter, ")[],"), c_iter, c_begin, c_end);
        if( *c_iter == '(' ) {
            ++c_iter;
            return strid;
        }

        strid += *c_iter;
        ++c_iter;
    }

    return strid;
}

string EvtMTree::parseKey( ptype& c_iter, ptype& c_begin, ptype& c_end ) 
{
    string key;

    while( *c_iter != ',' ) {
        parseerror(c_iter==c_end || parsecheck(*c_iter, "()[]"),
            c_iter, c_begin, c_end);
        key += *c_iter;
        ++c_iter;
    }

    ++c_iter;

    parseerror(c_iter == c_end, c_iter, c_begin, c_end);
    
    return key;
}

vector<string> EvtMTree::parseArg( ptype &c_iter, ptype &c_begin, ptype &c_end )
{
    vector<string> arg;

    if( *c_iter != '[' ) return arg;
    ++c_iter;

    string temp;
    while(true) {
        parseerror( c_iter == c_end || parsecheck(*c_iter, "[()"),
                c_iter, c_begin, c_end );

        if( *c_iter == ']' ) {
            ++c_iter;
            if(temp.size() > 0) arg.push_back( temp );
            break;
        }

        if( *c_iter == ',') {
            arg.push_back( temp );
            temp.clear();
            ++c_iter;
            continue;
        }

        temp += *c_iter;
        ++c_iter;
    }
    parseerror(c_iter == c_end || *c_iter != ',', c_iter, c_begin, c_end);
    ++c_iter;

    return arg;
}

vector<EvtComplex> EvtMTree::parseAmps( ptype &c_iter, 
        ptype &c_begin, ptype &c_end )
{
    vector<string> parg = parseArg( c_iter, c_begin, c_end );
    parseerror( parg.size() == 0, c_iter, c_begin, c_end );

    // Get parametrization amplitudes
    vector<string>::iterator amp_iter = parg.begin();
    vector<string>::iterator amp_end = parg.end();
    vector<EvtComplex> amps;

    while( amp_iter != amp_end ) {
        const char * nptr;
        char * endptr = NULL;
        double amp=0.0, phase=0.0;

        nptr = (*amp_iter).c_str();
        amp = strtod(nptr, &endptr);
        parseerror( nptr==endptr, c_iter, c_begin, c_end );

        ++amp_iter;
        parseerror( amp_iter == amp_end, c_iter, c_begin, c_end );

        nptr = (*amp_iter).c_str();
        phase = strtod(nptr, &endptr);
        parseerror( nptr==endptr, c_iter, c_begin, c_end );

        amps.push_back( amp*exp(EvtComplex(0.0, phase)) );

        ++amp_iter;
    }

    return amps;
}

vector<EvtMNode *> EvtMTree::duplicate( const vector<EvtMNode *>& list ) const
{
    vector<EvtMNode *> newlist;

    for(size_t i=0; i<list.size(); ++i)
        newlist.push_back( list[i]->duplicate() );

    return newlist;
}

// XXX Warning it is unsafe to use cl1 after a call to this function XXX
vector< vector<EvtMNode * > > EvtMTree::unionChildren( const string& nodestr,
        vector< vector<EvtMNode * > >& cl1 ) 
{
    vector<EvtMNode *> cl2 = parsenode( nodestr, false );
    vector< vector<EvtMNode * > > cl;
    
    if( cl1.size() == 0 ) {
        for( size_t i=0; i<cl2.size(); ++i ) {
            vector<EvtMNode *> temp(1, cl2[i]);
            cl.push_back( temp );
        }

        return cl;
    }

    for( size_t i=0; i<cl1.size(); ++i ) {
        for( size_t j=0; j<cl2.size(); ++j ) {
            vector<EvtMNode *> temp;
            temp = duplicate( cl1[i] );
            temp.push_back( cl2[j]->duplicate() );

            cl.push_back( temp );
        }
    }
 
    for(size_t i=0; i<cl1.size(); ++i)
        for(size_t j=0; j<cl1[i].size(); ++j)
            delete cl1[i][j];
    for(size_t i=0; i<cl2.size(); ++i)
        delete (cl2[i]);

    return cl;
}

vector< vector<EvtMNode * > > EvtMTree::parseChildren( ptype &c_iter, 
        ptype &c_begin, ptype &c_end ) 
{
    bool test = true;
    int pcount=0;
    string nodestr;
    vector< vector<EvtMNode * > > children;

    parseerror(c_iter == c_end || *c_iter != '[', c_iter, c_begin, c_end );
    ++c_iter;

    while( test ) {
        parseerror( c_iter==c_end || pcount < 0, c_iter, c_begin, c_end );

        switch( *c_iter ) {
            case ')':
                --pcount;
                nodestr += *c_iter;
                break;
            case '(':
                ++pcount;
                nodestr += *c_iter;
                break;
            case ']':
                if( pcount==0 ) {
                    children = unionChildren( nodestr, children );
                    test=false;
                } else {
                    nodestr += *c_iter;
                }
                break;
            case ',':
                if( pcount==0 ) {
                    children = unionChildren( nodestr, children );
                    nodestr.clear();
                } else {
                    nodestr += *c_iter;
                }
                break;
            default:
                nodestr += *c_iter;
                break;
        }

        ++c_iter;
    }

    return children;
}
    
vector<EvtMNode *> EvtMTree::parsenode( const string& args, bool rootnode )
{
    ptype c_iter, c_begin, c_end;

    c_iter=c_begin=args.begin();
    c_end = args.end();

    string strid = parseId( c_iter, c_begin, c_end );

    // Case 1: Particle
    if( c_iter == c_end ) return makeparticles( strid );

    // Case 2: Resonance - parse further
    EvtId id = EvtPDL::getId(strid);
    parseerror(EvtId( -1, -1 )==id, c_iter, c_begin, c_end);
    
    string ls;
    vector<string> lsarg;

    if( rootnode ) {
        ls = "TRIVIAL";
    } else {
        // Get lineshape (e.g. BREITWIGNER)
        ls = parseKey( c_iter, c_begin, c_end );
        lsarg = parseArg( c_iter, c_begin, c_end );
    }

    // Get resonance parametrization type (e.g. HELAMP)
    string type = parseKey( c_iter, c_begin, c_end );
    vector<EvtComplex> amps = parseAmps( c_iter, c_begin, c_end );

    // Children
    vector<vector<EvtMNode * > > children = parseChildren( c_iter, c_begin,
            c_end );

    report(Severity::Error,"EvtGen")<<children.size()<<endl;
    vector<EvtMNode *> resonances;
    for(size_t i=0; i<children.size(); ++i ) {
        resonances.push_back(makeresonance(id,ls,lsarg,type,amps,children[i]));
    }

    parseerror(c_iter == c_end || *c_iter!=')', c_iter, c_begin, c_end);

    return resonances;
}

bool EvtMTree::validTree( const EvtMNode * root ) const
{
    vector<int> res = root->getresonance();
    vector<bool> check(res.size(), false);

    for( size_t i=0; i<res.size(); ++i) {
        check[res[i]] = true;
    }

    bool ret = true;

    for( size_t i=0; i<check.size(); ++i ) {
        ret = ret&&check[i];
    }

    return ret;
}

void EvtMTree::addtree( const string& str )
{
    vector<EvtMNode *> roots = parsenode( str, true );
    _norm = 0;

    for( size_t i=0; i<roots.size(); ++i ) {
        if( validTree( roots[i] ) ) {
            _root.push_back( roots[i] );
            _norm = _norm + 1;
        } else
            delete roots[i];
    }
    
    _norm = 1.0/sqrt(_norm);
}
EvtSpinAmp EvtMTree::getrotation( EvtParticle * p ) const
{
    // Set up the rotation matrix for the root particle (for now)
    EvtSpinDensity sd = p->rotateToHelicityBasis();
    EvtSpinType::spintype type = EvtPDL::getSpinType(_root[0]->getid());
    int twospin = EvtSpinType::getSpin2(type);
    
    vector<EvtSpinType::spintype> types(2, type);
    EvtSpinAmp rot( types, EvtComplex(0.0, 0.0) );
    vector<int> index = rot.iterallowedinit();
    do {
        rot(index) = sd.get((index[0]+twospin)/2,(index[1]+twospin)/2);
    } while( rot.iterateallowed( index ) );

    return rot;
}

EvtSpinAmp EvtMTree::amplitude( EvtParticle * p ) const
{
    vector<EvtVector4R> product;
    for(size_t i=0; i<p->getNDaug(); ++i)
        product.push_back(p->getDaug(i)->getP4Lab());

    if( _root.size() == 0 ) {
        report(Severity::Error, "EvtGen")<<"No decay tree present."<<endl;
        ::abort();
    }
    
    EvtSpinAmp amp = _root[0]->amplitude( product );
    for( size_t i=1; i<_root.size(); ++i ) {
        // Assume that helicity amplitude is returned 
        amp += _root[i]->amplitude( product );
    }
    amp = _norm*amp;

    //ryd
    return amp;

    // Do Rotation to Proper Frame
    EvtSpinAmp newamp = getrotation( p );
    newamp.extcont(amp, 1, 0);

    return newamp;
}
