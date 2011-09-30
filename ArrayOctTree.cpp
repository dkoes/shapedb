/*
 * ArrayOctTree.cpp
 *
 *  Created on: Sep 29, 2011
 *      Author: dkoes
 */

#include "ArrayOctTree.h"
#include <cassert>

//recursive creation of oct tree
//push the oct tree defined by the overlay of cube and mol
//onto the back of tree and return the position of the created tree
ArrayOctTree::ChildNode ArrayOctTree::create(const Cube& cube,
		const vector<MolSphere>& mol)
{
	//does the mol overlap with this cube?
	vector<MolSphere> pruned;
	pruned.reserve(mol.size());
	for (unsigned i = 0, n = mol.size(); i < n; i++)
	{
		if (cube.intersectsSphere(mol[i]))
		{
			pruned.push_back(mol[i]);
		}
	}

	ChildNode ret;
	if (pruned.size() == 0)
	{
		//no overlap, all done
		ret.isLeaf = true;
		ret.isFull = false;
	}
	else if (cube.getDimension() <= resolution) //consider it full
	{
		ret.isLeaf = true;
		ret.isFull = true;
	}
	else //subdivide into children
	{
		//assume this is an interior node
		ret.isLeaf = false;
		unsigned fullcnt = 0;
		//allocate space
		unsigned pos = tree.size();
		ret.index = tree.size();
		tree.push_back(OctNode());
		for (unsigned i = 0; i < 8; i++)
		{
			Cube newc = cube.getOctant(i);
			ChildNode child = create(newc, pruned);
			tree[pos].children[i] = child;
			if (child.isLeaf && child.isFull)
				fullcnt++;
		}

		//are all the children full? then truncate and mark node as full
		if (fullcnt == 8)
		{
			tree.pop_back();
			ret.isLeaf = true;
			ret.isFull = true;
			ret.index = 0;
			assert(pos == tree.size());
		}
	}
	return ret;
}

//invert filled and unfilled
void ArrayOctTree::invert()
{
	for (unsigned i = 0, n = tree.size(); i < n; i++)
	{
		for (unsigned j = 0; j < 8; j++)
		{
			if (tree[i].children[j].isLeaf)
				tree[i].children[j].isFull = !tree[i].children[j].isFull;
		}
	}
}

//mogrifying intersection
bool ArrayOctTree::intersect(const OctTree* Rhs)
{
	const ArrayOctTree& rhs = dynamic_cast<const ArrayOctTree&> (*Rhs);
	vector<OctNode> newtree;
	bool changed = false;
	ChildNode newroot = root.intersect(tree,rhs.tree,rhs.root, newtree,changed);
	if(changed)
	{
		swap(root,newroot);
		swap(tree,newtree);
	}

	return changed;
}

//mogrifying union
bool ArrayOctTree::unionWith(const OctTree* Rhs)
{
	const ArrayOctTree& rhs = dynamic_cast<const ArrayOctTree&> (*Rhs);
	vector<OctNode> newtree;
	bool changed = false;
	ChildNode newroot = root.unionWith(tree,rhs.tree,rhs.root, newtree,changed);
	if(changed)
	{
		swap(root,newroot);
		swap(tree,newtree);
	}

	return changed;
}

//volume calculations that don't require creating a tmp tree
float ArrayOctTree::intersectVolume(const OctTree * Rhs) const
{
	const ArrayOctTree& rhs = dynamic_cast<const ArrayOctTree&> (*Rhs);

	return root.intersectVolume(tree, rhs.tree, rhs.root, dimension);
}

float ArrayOctTree::unionVolume(const OctTree *Rhs) const
{
	const ArrayOctTree& rhs = dynamic_cast<const ArrayOctTree&> (*Rhs);
	return root.unionVolume(tree, rhs.tree, rhs.root, dimension);
}

bool ArrayOctTree::containedIn(const OctTree *Rhs) const
{
	const ArrayOctTree& rhs = dynamic_cast<const ArrayOctTree&> (*Rhs);
	return root.containedIn(tree, rhs.tree, rhs.root);
}

//return total volume contained in octtree
float ArrayOctTree::volume() const
{
	return root.volume(tree, dimension);
}

//return number of leaves
unsigned ArrayOctTree::leaves() const
{
	unsigned cnt = 0;
	for (unsigned i = 0, n = tree.size(); i < n; i++)
	{
		for (unsigned j = 0; j < 8; j++)
		{
			if (tree[i].children[j].isLeaf)
				cnt++;
		}
	}
	return cnt;
}

void ArrayOctTree::clear()
{
	tree.clear();
	root = ChildNode(true,false);
}

void ArrayOctTree::fill()
{
	tree.clear();
	root = ChildNode(true,true);
}

void ArrayOctTree::write(ostream& out) const
{
	out.write((char*) &dimension, sizeof(dimension));
	out.write((char*) &resolution, sizeof(resolution));
	out.write((char*) &root, sizeof(root));
	unsigned sz = tree.size();
	out.write((char*) &sz, sizeof(sz));
	for (unsigned i = 0, n = tree.size(); i < n; i++)
	{
		out.write((char*) &tree[i], sizeof(tree[i]));
	}
}

void ArrayOctTree::read(istream& in)
{
	in.read((char*) &dimension, sizeof(dimension));
	in.read((char*) &resolution, sizeof(resolution));
	in.read((char*) &root, sizeof(root));
	unsigned sz = 0;
	in.read((char*) &sz, sizeof(sz));
	tree.resize(sz);
	for (unsigned i = 0; i < sz; i++)
	{
		in.read((char*) &tree[i], sizeof(tree[i]));
	}
}

unsigned ArrayOctTree::getOctantPattern(const vector<unsigned>& coord, bool MSV) const
{
	//find correct node
	ChildNode node = root; //root always has children
	for (unsigned i = 0, n = coord.size(); i < n; i++)
	{
		if (node.isLeaf)
			return node.getBitPattern(tree, MSV);
		else
			node = tree[node.index].children[coord[i]];
	}
	return node.getBitPattern(tree,MSV);
}

float ArrayOctTree::volumeDistance(const OctTree * Rhs) const
{
	const ArrayOctTree& rhs = dynamic_cast<const ArrayOctTree&> (*Rhs);
	float ival = 0, uval = 0;
	root.intersectUnionVolume(tree, rhs.tree, rhs.root, dimension, ival, uval);

	return 1 - ival/uval;
}


float ArrayOctTree::hausdorffDistance(const OctTree* B) const
{
	abort();
	return 0;
}

//copy the subtree of this childnode into new tree and return the new childnode
ArrayOctTree::ChildNode ArrayOctTree::ChildNode::copyTo(const vector<OctNode>& from, vector<OctNode>& to) const
{
	if(isLeaf)
		return *this;

	//create a node w/children in to
	unsigned pos = to.size();
	ChildNode ret(false, false, pos);
	to.push_back(OctNode());
	const OctNode& node = from[index];
	for(unsigned i = 0; i < 8; i++)
	{
		ChildNode child = node.children[i].copyTo(from, to);
		to[pos].children[i] = child;
	}
	return ret;
}

//compute the interesection and put it in newtree
//return the resulting childnode
ArrayOctTree::ChildNode ArrayOctTree::ChildNode::intersect(
		const vector<OctNode>& tree, const vector<OctNode>& rtree,
		const ChildNode& rhs, vector<OctNode>& newtree, bool& changed) const
{
	if(rhs.isLeaf && isLeaf && rhs.isFull == isFull)
	{
		return ChildNode(isLeaf, isFull);
	}
	else if(isLeaf && !isFull)
	{
		return ChildNode(true, false); //still empty
	}
	else if (rhs.isLeaf && !rhs.isFull)
	{
		changed = true;
		return ChildNode(true, false); //empty now
	}
	else if (rhs.isLeaf && rhs.isFull)
	{
		//copy this tree
		changed = true;
		return copyTo(tree, newtree);
	}
	else if (isLeaf && isFull)
	{
		//just copy from rhs
		changed = true;
		return rhs.copyTo(rtree, newtree);
	}
	else //both have children
	{
		unsigned numEmpty = 0;
		unsigned pos = newtree.size();
		newtree.push_back(OctNode());
		ChildNode ret(false, false);
		ret.index = pos;
		for (unsigned i = 0; i < 8; i++)
		{
			ChildNode nchild = tree[index].children[i].intersect(tree, rtree,
					rtree[rhs.index].children[i], newtree, changed);
			newtree[pos].children[i] = nchild;
			if(nchild.isLeaf && !nchild.isFull)
			{
				numEmpty++;
			}
		}

		if(numEmpty == 8)
		{
			newtree.pop_back();
			ret.isLeaf = true;
			ret.isFull = false;
			ret.index = 0;
		}
		return ret;
	}
}

ArrayOctTree::ChildNode ArrayOctTree::ChildNode::unionWith(
		const vector<OctNode>& tree, const vector<OctNode>& rtree,
		const ChildNode& rhs, vector<OctNode>& newtree, bool& changed) const
{
	if(rhs.isLeaf && isLeaf && rhs.isFull == isFull)
	{
		return ChildNode(isLeaf, isFull);
	}
	else if(isLeaf && isFull)
	{
		return ChildNode(true, true); //still full
	}
	else if (rhs.isLeaf && rhs.isFull)
	{
		changed = true;
		return ChildNode(true, true); //full now
	}
	else if (rhs.isLeaf && !rhs.isFull)
	{
		//copy this tree
		changed = true;
		return copyTo(tree, newtree);
	}
	else if (isLeaf && !isFull)
	{
		//just copy from rhs
		changed = true;
		return rhs.copyTo(rtree, newtree);
	}
	else //both have children
	{
		unsigned numFull = 0;
		unsigned pos = newtree.size();
		newtree.push_back(OctNode());
		ChildNode ret(false, false);
		ret.index = pos;
		for (unsigned i = 0; i < 8; i++)
		{
			ChildNode nchild = tree[index].children[i].unionWith(tree, rtree,
					rtree[rhs.index].children[i], newtree, changed);
			newtree[pos].children[i] = nchild;
			if(nchild.isLeaf && nchild.isFull)
			{
				numFull++;
			}
		}

		if(numFull == 8)
		{
			newtree.pop_back();
			ret.isLeaf = true;
			ret.isFull = true;
			ret.index = 0;
		}
		return ret;
	}
}


//compute the interesection volume
float ArrayOctTree::ChildNode::intersectVolume(
		const vector<OctNode>& tree, const vector<OctNode>& rtree,
		const ChildNode& rhs, float dim) const
{
	if(rhs.isLeaf && isLeaf && rhs.isFull == isFull)
	{
		if(isFull)
			return dim*dim*dim;
		else
			return 0;
	}
	else if(isLeaf && !isFull)
	{
		return 0;
	}
	else if (rhs.isLeaf && !rhs.isFull)
	{
		return 0;
	}
	else if (rhs.isLeaf && rhs.isFull)
	{
		//same as now
		return volume(tree, dim);
	}
	else if (isLeaf && isFull)
	{
		//same as rhs
		return rhs.volume(rtree, dim);
	}
	else //both have children
	{
		float ret = 0;
		for (unsigned i = 0; i < 8; i++)
		{
			ret += tree[index].children[i].intersectVolume(tree, rtree,
					rtree[rhs.index].children[i], dim/2);
		}

		return ret;
	}
}

float ArrayOctTree::ChildNode::unionVolume(
		const vector<OctNode>& tree, const vector<OctNode>& rtree,
		const ChildNode& rhs, float dim) const
{
	if(rhs.isLeaf && isLeaf && rhs.isFull == isFull)
	{
		if(isFull)
			return dim*dim*dim;
		else
			return 0;
	}
	else if(isLeaf && isFull)
	{
		return dim*dim*dim;
	}
	else if (rhs.isLeaf && rhs.isFull)
	{
		return dim*dim*dim;
	}
	else if (rhs.isLeaf && !rhs.isFull)
	{
		//vol of this tree
		return volume(tree, dim);
	}
	else if (isLeaf && !isFull)
	{
		// from rhs
		return rhs.volume(rtree, dim);
	}
	else //both have children
	{
		float ret = 0;
		for (unsigned i = 0; i < 8; i++)
		{
			ret += tree[index].children[i].unionVolume(tree, rtree,
					rtree[rhs.index].children[i], dim/2);
		}

		return ret;
	}
}

//compute both the intersection and union volume at once
void ArrayOctTree::ChildNode::intersectUnionVolume(
		const vector<OctNode>& tree, const vector<OctNode>& rtree,
		const ChildNode& rhs, float dim, float& intersectval, float& unionval) const
{
	if(rhs.isLeaf && isLeaf && rhs.isFull == isFull)
	{
		if(isFull)
		{
			intersectval += dim*dim*dim;
			unionval += dim*dim*dim;
		}
	}
	else if(isLeaf && isFull)
	{
		unionval += dim*dim*dim;
		intersectval += rhs.volume(rtree, dim);
	}
	else if (rhs.isLeaf && rhs.isFull)
	{
		unionval += dim*dim*dim;
		intersectval += volume(tree,dim);
	}
	else if (rhs.isLeaf && !rhs.isFull)
	{
		//vol of this tree
		unionval += volume(tree, dim);
	}
	else if (isLeaf && !isFull)
	{
		// from rhs
		unionval += rhs.volume(rtree, dim);
	}
	else //both have children
	{
		for (unsigned i = 0; i < 8; i++)
		{
			tree[index].children[i].intersectUnionVolume(tree, rtree,
					rtree[rhs.index].children[i], dim/2, intersectval,unionval);
		}
	}
}

//return true if this is contained in rhs - short circuit eval
bool ArrayOctTree::ChildNode::containedIn(
		const vector<OctNode>& tree, const vector<OctNode>& rtree,
		const ChildNode& rhs) const
{
	if(rhs.isLeaf && isLeaf && rhs.isFull == isFull)
	{
		return true;
	}
	else if(isLeaf && !isFull)
	{
		return true;
	}
	else if (rhs.isLeaf && !rhs.isFull)
	{
		return false; //somthing not in nothing (nothing handled above)
	}
	else if (rhs.isLeaf && rhs.isFull)
	{
		return true;
	}
	else if (isLeaf && isFull)
	{
		return false;
	}
	else //both have children
	{
		for (unsigned i = 0; i < 8; i++)
		{
			if(!tree[index].children[i].containedIn(tree, rtree,
					rtree[rhs.index].children[i]))
				return false;
		}
		return true;;
	}
}



unsigned ArrayOctTree::ChildNode::getBitPattern(const vector<OctNode>& tree, bool MSV) const
{
	if(isLeaf)
	{
		if(isFull)
			return 255;
		else
			return 0;
	}
	else
	{
		unsigned ret = 0;
		const OctNode& node = tree[index];
		for (unsigned i = 0; i < 8; i++)
		{
			if (node.children[i].isLeaf)
			{
				if (node.children[i].isFull)
					ret |= (1 << i);
			}
			else if (MSV)
			{
				ret |= (1 << i);
			}
		}
		return ret;
	}
}

float ArrayOctTree::ChildNode::volume(const vector<OctNode>& tree, float dim) const
{
	if(isLeaf)
	{
		if(isFull)
			return dim*dim*dim;
		else
			return 0;
	}
	else
	{
		if(volumeCache > 0)
			return volumeCache;
		float ret = 0;
		const OctNode& node = tree[index];
		for (unsigned i = 0; i < 8; i++)
		{
			ret += node.children[i].volume(tree, dim/2);
		}
		volumeCache = ret;
		return ret;
	}
}

