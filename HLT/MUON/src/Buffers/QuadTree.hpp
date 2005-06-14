////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_BUFFERS_QUAD_TREE_HPP
#define dHLT_BUFFERS_QUAD_TREE_HPP

#include "BasicTypes.hpp"
#include "Utils.hpp"
#include "new.hpp"

namespace dHLT
{
namespace Buffers
{


template <typename DataType>
class QuadTree
{
public:

	class Handler
	{
	public:

		/* Called when the specified node falls exactly within the region
		   passed to ForAllOverlapping.
		 */
		inline void AddTo(QuadTree* /*node*/) {};
		
		/* Called for the first node of a tree when the specified node has
		   just been created and needs to be initialised. This method 
		   corresponds to FoundOverlapping except it is called when the node
		   was not found and must be created.
		 */
		inline void InitialiseFirst(QuadTree* /*node*/) {};
		
		/* Called when the specified node has just been created and needs to be
		   initialised. The specified parent is the parent of the new node. This
		   method corresponds to FoundOverlapping except it is called when the
		   node was not found and must be created.
		 */
		inline void Initialise(QuadTree* /*node*/, QuadTree* /*parent*/) {};
		
		/* Called when a node was found that partialy overlaps the region given
		   to ForAllOverlapping. This method corresponds to Initialise(First)
		   except it is called when the node has been found.
		 */
		inline void FoundOverlapping(QuadTree* /*node*/) {};
		
		/* Called when a child node was found to be contained in the region
		   given to the ForAllOverlapping method.
		 */
		inline void FoundContained(QuadTree* /*node*/) {};
	};


	QuadTree()
	{
		child[0][0] = child[0][1] = child[1][0] = child[1][1] = NULL;
	};


	~QuadTree()
	{
		// Perform a recursive delete of all children.
		if (child[0][0] != NULL) delete child[0][0];
		if (child[0][1] != NULL) delete child[0][1];
		if (child[1][0] != NULL) delete child[1][0];
		if (child[1][1] != NULL) delete child[1][1];
	};


	/* Checks to see if the child node at grid position x, y on grid level 'level'
	   exists or not. True is retured if it exists otherwise false is returned.
	 */
	bool Empty(const UInt x, const UInt y, const UChar level) const
	{
		return child[(x >> level) & 1][(y >> level) & 1] == NULL;
	};


	/* Adds a new child node to this parent node at grid position x, y
	   on grid level 'level'.
	 */
	QuadTree* Add(const UInt x, const UInt y, const UChar level)
	{
		Assert( Empty(x, y, level) );
		QuadTree* newnode = new QuadTree();
		child[(x >> level) & 1][(y >> level) & 1] = newnode;
		return newnode;
	};


	/* Deletes the child subtree of this node at grid position x, y
	   on grid level 'level'.
	 */
	void Remove(const UInt x, const UInt y, const UChar level)
	{
		Assert( not Empty(x, y, level) );
		delete child[(x >> level) & 1][(y >> level) & 1];
		child[(x >> level) & 1][(y >> level) & 1] = NULL;
	};


	/* Gets the child of this node at grid position x, y on grid level 'level'.
	 */
	QuadTree* Child(const UInt x, const UInt y, const UChar level)
	{
		return child[(x >> level) & 1][(y >> level) & 1];
	};
	
	// Constant method version of Child above.
	const QuadTree* Child(const UInt x, const UInt y, const UChar level) const
	{
		return child[(x >> level) & 1][(y >> level) & 1];
	};


	/* Performs a depth first traversal to the specified level (tree depth),
	   calling FoundOverlapping or Initialise method of the given event handler.
	   Initialise is called if the node has to be created.
	   When the given level is reached then the handlers AddTo method is called
	   for all 4 nodes on that level, for grid positions: (left, bottom) ;
	   (left+1, bottom) ; (left, bottom+1) ; (left+1, bottom+1).
	   For all children nodes for levels i > 'level' the tree is traversed and
	   FoundContained method is called for every child node found.
	 */
	template <class Handler>
	void ForAllOverlapping(const UInt left, const UInt bottom, const UChar level, Handler& handler)
	{
		UInt right = left + 1;
		UInt top = bottom + 1;

		QuadTree* current = this;

		UChar i;
		for (i = level; i > 0; i--)
		{
			if ( OnSameBranch(left, right, i) )
			{
				if ( OnSameBranch(bottom, top, i) )
				{
					if ( current->Empty(left, bottom, i) )
					{
						QuadTree* parent = current;
						current = current->Add(left, bottom, i);
						handler.Initialise(current, parent);
					}
					else
					{
						current = current->Child(left, bottom, i);
						handler.FoundOverlapping(current);
					};
				}
				else
				{
					ForAllOverlappingHorizontal(current, left, right, bottom, i, handler);
					ForAllOverlappingHorizontal(current, left, right, top, i, handler);
					return;
				};
			}
			else
			{
				if ( OnSameBranch(bottom, top, i) )
				{
					ForAllOverlappingVertical(current, left, bottom, top, i, handler);
					ForAllOverlappingVertical(current, right, bottom, top, i, handler);
				}
				else
				{
					ForAllOverlapping(current, left, bottom, i, handler);
					ForAllOverlapping(current, left, top, i, handler);
					ForAllOverlapping(current, right, bottom, i, handler);
					ForAllOverlapping(current, right, top, i, handler);
				};
				return;
			};
		};

		QuadTree* leaf;
		if ( OnSameBranch(left, right, i) )
		{
			leaf = current->Fetch(left, bottom, i, handler);
			handler.AddTo(leaf);
			leaf->Traverse(handler);
			if ( not OnSameBranch(bottom, top, i) )
			{
				leaf = current->Fetch(left, top, i, handler);
				handler.AddTo(leaf);
				leaf->Traverse(handler);
			}
		}
		else
		{
			leaf = current->Fetch(left, bottom, i, handler);
			handler.AddTo(leaf);
			leaf->Traverse(handler);
			leaf = current->Fetch(right, bottom, i, handler);
			handler.AddTo(leaf);
			leaf->Traverse(handler);
			if ( not OnSameBranch(bottom, top, i) )
			{
				leaf = current->Fetch(left, top, i, handler);
				handler.AddTo(leaf);
				leaf->Traverse(handler);
				leaf = current->Fetch(right, top, i, handler);
				handler.AddTo(leaf);
				leaf->Traverse(handler);
			};
		};
	};


	/*
	template <class Function>
	void Traverse(Function process)
	{
		QuadTree* stack[14*4];
		QuadTree** bottom = &stack[14*4];
		QuadTree** top = &stack[14*4-1];
		*top = this;
		do
		{
			QuadTree* current = *top++;  // pop the stack;
			(*process)(current);
			// Push the children of the current node onto the stack.
			if (current->child[0][0] != NULL) *(--top) = current->child[0][0];
			if (current->child[0][1] != NULL) *(--top) = current->child[0][1];
			if (current->child[1][0] != NULL) *(--top) = current->child[1][0];
			if (current->child[1][1] != NULL) *(--top) = current->child[1][1];
		}
		while (top != bottom);
	};
	*/

	/* Traverses the whole sub tree stopping at every child and calling the
	   FoundContained method of the given event handler.
	 */
	template <class Handler>
	void Traverse(Handler& handler)
	{
		if (child[0][0] != NULL)
		{
			handler.FoundContained(child[0][0]);
			child[0][0]->Traverse(handler);
		};
		if (child[0][1] != NULL)
		{
			handler.FoundContained(child[0][1]);
			child[0][1]->Traverse(handler);
		};
		if (child[1][0] != NULL)
		{
			handler.FoundContained(child[1][0]);
			child[1][0]->Traverse(handler);
		};
		if (child[1][1] != NULL)
		{
			handler.FoundContained(child[1][1]);
			child[1][1]->Traverse(handler);
		};
	};


	/* Finds the child that corresponds to the x, y grid position on the
	   given grid level. (Level can also be though of as tree depth).
	 */
	QuadTree* Find(const UInt x, const UInt y, const UChar level)
	{
		QuadTree* current = this;
		for (Char i = level; i >= 0; i--)
		{
			current = current->Child(x, y, i);
			if (current == NULL) break;
		};
		return current;
	};
	
	// The constant version of Find above.
	const QuadTree* Find(const UInt x, const UInt y, const UChar level) const
	{
		QuadTree* current = this;
		for (Char i = level; i >= 0; i--)
		{
			current = current->Child(x, y, i);
			if (current == NULL) break;
		};
		return current;
	};

	
	DataType data;   // The data stored in the tree.

private:
	
	inline bool OnSameBranch(const UInt a, const UInt b, const UChar level) const
	{
		return a >> level == b >> level;
	};


	template <class Handler>
	inline QuadTree* Fetch(const UInt x, const UInt y, const UChar level, Handler& handler)
	{
		QuadTree* node;
		if ( Empty(x, y, level) )
		{
			node = Add(x, y, level);
			handler.Initialise(node, this);   // use this as the parent
		}
		else
		{
			node = Child(x, y, level);
			handler.FoundContained(node);
		};
		return node;
	};


	template <class Handler>
	inline void ForAllOverlapping(QuadTree* current, const UInt x, const UInt y, const UChar level, Handler& handler)
	{
		UChar i;
		for (i = level; i > 0; i--)
		{
			if ( current->Empty(x, y, i) )
			{
				QuadTree* parent = current;
				current = current->Add(x, y, i);
				handler.Initialise(current, parent);
			}
			else
			{
				current = current->Child(x, y, i);
				handler.FoundOverlapping(current);
			};
		};

		QuadTree* leaf = current->Fetch(x, y, i, handler);
		handler.AddTo(leaf);
		leaf->Traverse(handler);
	};


	template <class Handler>
	inline void ForAllOverlappingHorizontal(
			QuadTree* current, const UInt left, const UInt right, const UInt bottom, const UChar level, Handler& handler
		)
	{
		UChar i;
		for (i = level; i > 0; i--)
		{
			if ( OnSameBranch(left, right, i) )
			{
				if ( current->Empty(left, bottom, i) )
				{
					QuadTree* parent = current;
					current = current->Add(left, bottom, i);
					handler.Initialise(current, parent);
				}
				else
				{
					current = current->Child(left, bottom, i);
					handler.FoundOverlapping(current);
				};
			}
			else
			{
				ForAllOverlapping(current, left, bottom, i, handler);
				ForAllOverlapping(current, right, bottom, i, handler);
				return;
			};
		};

		QuadTree* leaf = current->Fetch(left, bottom, i, handler);
		handler.AddTo(leaf);
		leaf->Traverse(handler);
		if ( not OnSameBranch(left, right, i) )
		{
			leaf = current->Fetch(right, bottom, i, handler);
			handler.AddTo(leaf);
			leaf->Traverse(handler);
		};
	};


	template <class Handler>
	inline void ForAllOverlappingVertical(
			QuadTree* current, const UInt left, const UInt bottom, const UInt top, const UChar level, Handler& handler
		)
	{
		UChar i;
		for (i = level; i > 0; i--)
		{
			if ( OnSameBranch(bottom, top, i) )
			{
				if ( current->Empty(left, bottom, i) )
				{
					QuadTree* parent = current;
					current = current->Add(left, bottom, i);
					handler.Initialise(current, parent);
				}
				else
				{
					current = current->Child(left, bottom, i);
					handler.FoundOverlapping(current);
				};
			}
			else
			{
				ForAllOverlapping(current, left, bottom, i, handler);
				ForAllOverlapping(current, left, top, i, handler);
				return;
			};
		};

		QuadTree* leaf = current->Fetch(left, bottom, i, handler);
		handler.AddTo(leaf);
		leaf->Traverse(handler);
		if ( not OnSameBranch(bottom, top, i) )
		{
			leaf = current->Fetch(left, top, i, handler);
			handler.AddTo(leaf);
			leaf->Traverse(handler);
		};
	};

	
	QuadTree* child[2][2];
};


} // Buffers
} // dHLT

#endif // dHLT_BUFFERS_QUAD_TREE_HPP
