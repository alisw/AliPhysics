////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_BUFFERS_REGION_TREE_HPP
#define dHLT_BUFFERS_REGION_TREE_HPP

#include "Buffers/QuadTree.hpp"
#include "RegionOfInterest.hpp"
#include "new.hpp"

namespace dHLT
{
namespace Buffers
{


template <typename DataType>
class RegionTree
{
public:

	typedef QuadTree<DataType> Node;
	class Handler : public Node::Handler {};


	RegionTree()
	{
		for (UChar i = 0; i < 10; i++)
			root[i] = NULL;
	};


	~RegionTree()
	{
		for (UChar i = 0; i < 10; i++)
			if (root[i] != NULL) delete root[i];
	};


	bool Empty(const ChamberID chamber) const
	{
		Assert( chamber < 10 );
		return root[chamber] == NULL;
	};


	Node* Add(const ChamberID chamber)
	{
		Assert( chamber < 10 );
		Assert( root[chamber] == NULL );
		root[chamber] = new Node();
		return root[chamber];
	};


	void Remove(const ChamberID chamber)
	{
		Assert( chamber < 10 );
		Assert( root[chamber] != NULL );
		delete root[chamber];
		root[chamber] = NULL;
	};


	Node* Get(const ChamberID chamber)
	{
		Assert( chamber < 10 );
		return root[chamber];
	};


	const Node* Get(const ChamberID chamber) const
	{
		Assert( chamber < 10 );
		return root[chamber];
	};


	template <class Handler>
	void ForAllOverlapping(
			const ChamberID chamber, const UInt left, const UInt bottom, const UChar level,
			Handler& handler
		)
	{
		Node* node;
		if (Empty(chamber))
		{
			node = Add(chamber);
			handler.InitialiseFirst(node);
			if (level != 0)
				node->ForAllOverlapping(left, bottom, level, handler);
			else
				handler.AddTo(node);
		}
		else
		{
			node = Get(chamber);
			if (level != 0)
			{
				handler.FoundOverlapping(node);
				node->ForAllOverlapping(left, bottom, level, handler);
			}
			else
			{
				handler.FoundContained(node);
				handler.AddTo(node);
				node->Traverse(handler);
			};
		};
	};


	template <class Handler>
	void ForAllContainedIn(
			const ChamberID chamber, const UInt left, const UInt bottom, const UChar level,
			Handler handler
		)
	{
		if ( ! Empty(chamber) )
		{
			Node* treeroot = Get(chamber);
			if (level != 0)
			{
				Node* node = treeroot->Find(left, bottom, level);
				if (node != NULL)
				{
					handler.FoundContained(node);
					node->Traverse(handler);
				};

				node = treeroot->Find(left + 1, bottom, level);
				if (node != NULL)
				{
					handler.FoundContained(node);
					node->Traverse(handler);
				};

				node = treeroot->Find(left, bottom + 1, level);
				if (node != NULL)
				{
					handler.FoundContained(node);
					node->Traverse(handler);
				};

				node = treeroot->Find(left + 1, bottom + 1, level);
				if (node != NULL)
				{
					handler.FoundContained(node);
					node->Traverse(handler);
				};
			}
			else
				handler.FoundContained(treeroot);
		};
	};


private:

	Node* root[10];
};


} // Buffers
} // dHLT

#endif // dHLT_BUFFERS_REGION_TREE_HPP
