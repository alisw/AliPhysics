////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>

using std::endl;
using std::cout;


#include "Tracking/IOInterface.hpp"
#include "Tracking/Tracker.hpp"
#include "Clustering/ClusterFinder.hpp"
#include "Decision/IOInterface.hpp"
#include "Utils.hpp"
#include "Debug/ClusterSource.hpp"
#include "System/Thread.hpp"
#include "System/Mutex.hpp"
#include "System/MutexCondition.hpp"
#include "Buffers/Queue.hpp"
#include "Tracking/IOHandler.hpp"
#include "Tracking/MansoTracker.hpp"
#include "Buffers/RegionTree.hpp"


using namespace dHLT;


int main()
{

	cout << "done." << endl;

	return 0;
};
