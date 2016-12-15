#include "StageAndFilter.h"
#include "Riostream.h"
#include "TString.h"

int main(int argc, char** argv)
{
	if ( argc < 4 )
	{
		std::cout << "Usage : " << argv[0] << " --from alien://source --to destination --filter filterName [--verbose [level]]" << std::endl;
		return 1;
	}

	TString from;
	TString to;
	TString filter;
	Int_t verboseLevel;

	for ( int i = 0; i < argc; ++i )
	{
		if ( TString(argv[i]) == "--from" )
		{
			from = argv[i+1];
		}
		if ( TString(argv[i]) == "--to" )
		{
			to = argv[i+1];
		}
		if ( TString(argv[i]) == "--filter" )
		{
			filter = argv[i+1];
		}
		if ( TString(argv[i]) == "--verbose" )
		{
			verboseLevel = TString(argv[i+1]).Atoi();
		}
	}

	Int_t rv = AAF::StageAndFilter(from,to,filter,verboseLevel);

	return rv;
}
