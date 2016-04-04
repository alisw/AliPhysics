//--------------------------------------------------------------------------
//
// Module: EvtD0gammaDalitz.cc
//
// Modification history:
//
//    JGT     February 13, 2012         Module created
//
//------------------------------------------------------------------------

#include <cstdlib>
#include <string>
#include <cmath>

#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtResonance.hh"
#include "EvtGenBase/EvtResonance2.hh"
#include "EvtGenModels/EvtD0gammaDalitz.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtFlatte.hh"
#include "EvtGenBase/EvtDecayTable.hh"


// Initialize the static variables.
const EvtSpinType::spintype& EvtD0gammaDalitz::_SCALAR = EvtSpinType::SCALAR;
const EvtSpinType::spintype& EvtD0gammaDalitz::_VECTOR = EvtSpinType::VECTOR;
const EvtSpinType::spintype& EvtD0gammaDalitz::_TENSOR = EvtSpinType::TENSOR;

const EvtDalitzReso::CouplingType& EvtD0gammaDalitz::_EtaPic   = EvtDalitzReso::EtaPic;
const EvtDalitzReso::CouplingType& EvtD0gammaDalitz::_PicPicKK = EvtDalitzReso::PicPicKK;

const EvtDalitzReso::NumType& EvtD0gammaDalitz::_RBW   = EvtDalitzReso::RBW_CLEO_ZEMACH;
const EvtDalitzReso::NumType& EvtD0gammaDalitz::_GS    = EvtDalitzReso::GS_CLEO_ZEMACH;
const EvtDalitzReso::NumType& EvtD0gammaDalitz::_KMAT  = EvtDalitzReso::K_MATRIX;

const EvtCyclic3::Pair& EvtD0gammaDalitz::_AB = EvtCyclic3::AB;
const EvtCyclic3::Pair& EvtD0gammaDalitz::_AC = EvtCyclic3::AC;
const EvtCyclic3::Pair& EvtD0gammaDalitz::_BC = EvtCyclic3::BC;


EvtD0gammaDalitz::EvtD0gammaDalitz()
{
  /* Empty constructor. */
}


EvtD0gammaDalitz::~EvtD0gammaDalitz()
{
  /* Empty destructor. */
}


std::string EvtD0gammaDalitz::getName()
{
  return "D0GAMMADALITZ";
}


EvtDecayBase* EvtD0gammaDalitz::clone()
{
  return new EvtD0gammaDalitz;
}


void EvtD0gammaDalitz::init()
{
  // check that there are 0 arguments
  checkNArg(0);

  // Check that this model is valid for the specified decay.
  checkNDaug( 3 );
  checkSpinParent  (    _SCALAR );
  checkSpinDaughter( 0, _SCALAR );
  checkSpinDaughter( 1, _SCALAR );
  checkSpinDaughter( 2, _SCALAR );

  // Get the values of the EvtId objects from the data files.
  readPDGValues();

  // Get the EvtId of the D0 and its 3 daughters.
  getParentId();

  EvtId dau[ 3 ];
  for ( int index = 0; index < 3; index++ )
  {
    dau[ index ] = getDaug( index );
  }

  // Look for K0bar h+ h-. The order will be K[0SL] h+ h-
  for ( int index = 0; index < 3; index++ )
  {
    if      ( ( dau[ index ] == _K0B ) || ( dau[ index ] == _KS ) || ( dau[ index ] == _KL ) )
    {
      _d1 = index;
    }
    else if ( ( dau[ index ] == _PIP ) || ( dau[ index ] == _KP ) )
    {
      _d2 = index;
    }
    else if ( ( dau[ index ] == _PIM ) || ( dau[ index ] == _KM ) )
    {
      _d3 = index;
    }
    else
    {
      reportInvalidAndExit();
    }
  }

  // Check if we're dealing with Ks pi pi or with Ks K K.
  _isKsPiPi = false;
  if ( dau[ _d2 ] == _PIP || dau[ _d2 ] == _PIM )
  {
    _isKsPiPi = true;
  }

}


void EvtD0gammaDalitz::initProbMax()
{
  setProbMax( 5200. );
}


void EvtD0gammaDalitz::decay( EvtParticle* part )
{
  // Check if the D is from a B+- -> D0 K+- decay with the appropriate model.
  EvtParticle* parent = part->getParent(); // If there are no mistakes, should be B+ or B-.
  if (parent != 0 && EvtDecayTable::getInstance()->getDecayFunc( parent )->getName() == "BTODDALITZCPK" )
  {
    EvtId parId = parent->getId();
    if ( ( parId == _BP ) || ( parId == _BM ) ||
         ( parId == _B0 ) || ( parId == _B0B) )
    {
      _bFlavor = parId;
    }
    else
    {
      reportInvalidAndExit();
    }
  }
  else
  {
    reportInvalidAndExit();
  }

  // Read the D decay parameters from the B decay model.
  // Gamma angle in rad.
  double gamma = EvtDecayTable::getInstance()->getDecayFunc( parent )->getArg( 0 );
  // Strong phase in rad.
  double delta = EvtDecayTable::getInstance()->getDecayFunc( parent )->getArg( 1 );
  // Ratio between B->D0K and B->D0barK
  double rB    = EvtDecayTable::getInstance()->getDecayFunc( parent )->getArg( 2 );

  // Same structure for all of these decays.
  part->initializePhaseSpace( getNDaug(), getDaugs() );
  EvtVector4R pA = part->getDaug( _d1 )->getP4();
  EvtVector4R pB = part->getDaug( _d2 )->getP4();
  EvtVector4R pC = part->getDaug( _d3 )->getP4();

  // Squared invariant masses.
  double mSqAB = ( pA + pB ).mass2();
  double mSqAC = ( pA + pC ).mass2();
  double mSqBC = ( pB + pC ).mass2();

  EvtComplex amp( 1.0, 0.0 );

  // Direct and conjugated amplitudes.
  EvtComplex ampDir;
  EvtComplex ampCnj;

  if ( _isKsPiPi )
  {
    // Direct and conjugated Dalitz points.
    EvtDalitzPoint pointDir( _mKs, _mPi, _mPi, mSqAB, mSqBC, mSqAC );
    EvtDalitzPoint pointCnj( _mKs, _mPi, _mPi, mSqAC, mSqBC, mSqAB );

    // Direct and conjugated amplitudes.
    ampDir = dalitzKsPiPi( pointDir );
    ampCnj = dalitzKsPiPi( pointCnj );
  }
  else
  {
    // Direct and conjugated Dalitz points.
    EvtDalitzPoint pointDir( _mKs, _mK, _mK, mSqAB, mSqBC, mSqAC );
    EvtDalitzPoint pointCnj( _mKs, _mK, _mK, mSqAC, mSqBC, mSqAB );

    // Direct and conjugated amplitudes.
    ampDir = dalitzKsKK( pointDir );
    ampCnj = dalitzKsKK( pointCnj );
  }

  if ( _bFlavor == _BP || _bFlavor == _B0 )
  {
    amp = ampCnj + rB * exp( EvtComplex( 0., delta + gamma ) ) * ampDir;
  }
  else
  {
    amp = ampDir + rB * exp( EvtComplex( 0., delta - gamma ) ) * ampCnj;
  }

  vertex( amp );

  return;

}


EvtComplex EvtD0gammaDalitz::dalitzKsPiPi( const EvtDalitzPoint& point ) const
{
  static const EvtDalitzPlot plot( _mKs, _mPi, _mPi, _mD0 );

  EvtComplex amp = 0.;

  // This corresponds to relativistic Breit-Wigner distributions. Not K-matrix.
  // Defining resonances.
  static EvtDalitzReso KStarm      ( plot, _BC, _AC, _VECTOR, 0.893606, 0.0463407, _RBW );
  static EvtDalitzReso KStarp      ( plot, _BC, _AB, _VECTOR, 0.893606, 0.0463407, _RBW );
  static EvtDalitzReso rho0        ( plot, _AC, _BC, _VECTOR, 0.7758  , 0.1464   , _GS  );
  static EvtDalitzReso omega       ( plot, _AC, _BC, _VECTOR, 0.78259 , 0.00849  , _RBW );
  static EvtDalitzReso f0_980      ( plot, _AC, _BC, _SCALAR, 0.975   , 0.044    , _RBW );
  static EvtDalitzReso f0_1370     ( plot, _AC, _BC, _SCALAR, 1.434   , 0.173    , _RBW );
  static EvtDalitzReso f2_1270     ( plot, _AC, _BC, _TENSOR, 1.2754  , 0.1851   , _RBW );
  static EvtDalitzReso K0Starm_1430( plot, _BC, _AC, _SCALAR, 1.459   , 0.175    , _RBW );
  static EvtDalitzReso K0Starp_1430( plot, _BC, _AB, _SCALAR, 1.459   , 0.175    , _RBW );
  static EvtDalitzReso K2Starm_1430( plot, _BC, _AC, _TENSOR, 1.4256  , 0.0985   , _RBW );
  static EvtDalitzReso K2Starp_1430( plot, _BC, _AB, _TENSOR, 1.4256  , 0.0985   , _RBW );
  static EvtDalitzReso sigma       ( plot, _AC, _BC, _SCALAR, 0.527699, 0.511861 , _RBW );
  static EvtDalitzReso sigma2      ( plot, _AC, _BC, _SCALAR, 1.03327 , 0.0987890, _RBW );
  static EvtDalitzReso KStarm_1680 ( plot, _BC, _AC, _VECTOR, 1.677   , 0.205    , _RBW );

  // Adding terms to the amplitude with their corresponding amplitude and phase terms.
  amp += EvtComplex(   .848984 ,   .893618  );
  amp += EvtComplex( -1.16356  ,  1.19933   ) * KStarm      .evaluate( point );
  amp += EvtComplex(   .106051 , - .118513  ) * KStarp      .evaluate( point );
  amp += EvtComplex(  1.0      ,  0.0       ) * rho0        .evaluate( point );
  amp += EvtComplex( - .0249569,   .0388072 ) * omega       .evaluate( point );
  amp += EvtComplex( - .423586 , - .236099  ) * f0_980      .evaluate( point );
  amp += EvtComplex( -2.16486  ,  3.62385   ) * f0_1370     .evaluate( point );
  amp += EvtComplex(   .217748 , - .133327  ) * f2_1270     .evaluate( point );
  amp += EvtComplex(  1.62128  ,  1.06816   ) * K0Starm_1430.evaluate( point );
  amp += EvtComplex(   .148802 ,   .0897144 ) * K0Starp_1430.evaluate( point );
  amp += EvtComplex(  1.15489  , - .773363  ) * K2Starm_1430.evaluate( point );
  amp += EvtComplex(   .140865 , - .165378  ) * K2Starp_1430.evaluate( point );
  amp += EvtComplex( -1.55556  , - .931685  ) * sigma       .evaluate( point );
  amp += EvtComplex( - .273791 , - .0535596 ) * sigma2      .evaluate( point );
  amp += EvtComplex( -1.69720  ,   .128038  ) * KStarm_1680 .evaluate( point );

  return amp;
}


EvtComplex EvtD0gammaDalitz::dalitzKsKK( const EvtDalitzPoint& point ) const
{
  static const EvtDalitzPlot plot( _mKs, _mK, _mK, _mD0 );

  // Defining resonances.
  static EvtDalitzReso a00_980 ( plot, _AC, _BC, _SCALAR, 0.999  , _RBW, .550173, .324, _EtaPic   );
  static EvtDalitzReso phi     ( plot, _AC, _BC, _VECTOR, 1.01943,       .00459319    , _RBW      );
  static EvtDalitzReso a0p_980 ( plot, _AC, _AB, _SCALAR, 0.999  , _RBW, .550173, .324, _EtaPic   );
  static EvtDalitzReso f0_1370 ( plot, _AC, _BC, _SCALAR, 1.350  ,       .265         , _RBW      );
  static EvtDalitzReso a0m_980 ( plot, _AB, _AC, _SCALAR, 0.999  , _RBW, .550173, .324, _EtaPic   );
  static EvtDalitzReso f0_980  ( plot, _AC, _BC, _SCALAR, 0.965  , _RBW, .695   , .165, _PicPicKK );
  static EvtDalitzReso f2_1270 ( plot, _AC, _BC, _TENSOR, 1.2754 ,       .1851        , _RBW      );
  static EvtDalitzReso a00_1450( plot, _AC, _BC, _SCALAR, 1.474  ,       .265         , _RBW      );
  static EvtDalitzReso a0p_1450( plot, _AC, _AB, _SCALAR, 1.474  ,       .265         , _RBW      );
  static EvtDalitzReso a0m_1450( plot, _AB, _AC, _SCALAR, 1.474  ,       .265         , _RBW      );

  // Adding terms to the amplitude with their corresponding amplitude and phase terms.
  EvtComplex amp( 0., 0. ); // Phase space amplitude.
  amp += EvtComplex( 1.0          , 0.0        ) * a00_980 .evaluate( point );
  amp += EvtComplex( -.126314     ,  .188701   ) * phi     .evaluate( point );
  amp += EvtComplex( -.561428     ,  .0135338  ) * a0p_980 .evaluate( point );
  amp += EvtComplex(  .035        , -.00110488 ) * f0_1370 .evaluate( point );
  amp += EvtComplex( -.0872735    ,  .0791190  ) * a0m_980 .evaluate( point );
  amp += EvtComplex( 0.           , 0.         ) * f0_980  .evaluate( point );
  amp += EvtComplex(  .257341     , -.0408343  ) * f2_1270 .evaluate( point );
  amp += EvtComplex( -.0614342    , -.649930   ) * a00_1450.evaluate( point );
  amp += EvtComplex( -.104629     ,  .830120   ) * a0p_1450.evaluate( point );
  amp += EvtComplex( 0.           , 0.         ) * a0m_1450.evaluate( point );

  return 2.8 * amp; // Multiply by 2.8 in order to reuse the same probmax as Ks pi pi.
}


void EvtD0gammaDalitz::readPDGValues()
{
  // Define the EvtIds.
  _BP  = EvtPDL::getId( "B+"      );
  _BM  = EvtPDL::getId( "B-"      );
  _B0  = EvtPDL::getId( "B0"      );
  _B0B = EvtPDL::getId( "anti-B0" );
  _D0  = EvtPDL::getId( "D0"      );
  _D0B = EvtPDL::getId( "anti-D0" );
  _KM  = EvtPDL::getId( "K-"      );
  _KP  = EvtPDL::getId( "K+"      );
  _K0  = EvtPDL::getId( "K0"      );
  _K0B = EvtPDL::getId( "anti-K0" );
  _KL  = EvtPDL::getId( "K_L0"    );
  _KS  = EvtPDL::getId( "K_S0"    );
  _PIM = EvtPDL::getId( "pi-"     );
  _PIP = EvtPDL::getId( "pi+"     );

  // Read the relevant masses.
  _mD0 = EvtPDL::getMass( _D0  );
  _mKs = EvtPDL::getMass( _KS  );
  _mPi = EvtPDL::getMass( _PIP );
  _mK  = EvtPDL::getMass( _KP  );
}


void EvtD0gammaDalitz::reportInvalidAndExit() const
{
  report( Severity::Error, "EvtD0gammaDalitz" ) << "EvtD0gammaDalitz: Invalid mode." << std::endl;
  exit( 1 );
}
