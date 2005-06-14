////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DEBUG_TRIGGER_SOURCE_HPP
#define dHLT_DEBUG_TRIGGER_SOURCE_HPP

#include <vector>
#include "DDL/TriggerInputInterface.hpp"

namespace dHLT
{
namespace Debug
{


class TriggerSource
{
public:

	TriggerSource();

	virtual ~TriggerSource();
	
	
	void NewEvent(const EventID id);
	
	void NewTriggerBlock();
	
	void AddTrigger(const TriggerRecord& t);
	
	void Clear();
	
	void Process();
	
	void Dump();
	
	void SetCallback(DDL::TriggerInputCallback* callback)
	{
		framework = callback;
	};
	
	DDL::TriggerInputCallback* GetCallback() const
	{
		return framework;
	};


protected:

	DDL::TriggerInputCallback* framework;

	struct DataBlock
	{
		std::vector<TriggerRecord> trigger;
	};
	
	struct EventBlock
	{
		EventID eventid;
		std::vector<DataBlock> block;
	};
	
	std::vector<EventBlock> event;
	
	EventBlock* currentevent;
	DataBlock* currentblock;

};


} // Debug
} // dHLT

#endif // dHLT_DEBUG_TRIGGER_SOURCE_HPP
