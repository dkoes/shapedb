/*
 * MappableOctTree.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: dkoes
 */

#include "MappableOctTree.h"
#include <cstring>
#include <cassert>
#include <boost/multi_array.hpp>
using namespace boost;

MappableOctTree* MappableOctTree::clone() const
{
	unsigned sz = bytes();
	void* mem = malloc(sz);
	memcpy(mem, this, sz);
	return (MappableOctTree*)mem;
}



MChildNode MappableOctTree::createFrom_r(unsigned N, MChildNode* nodes, const MappableOctTree** trees, vector<MOctNode>& newtree, bool isUnion)
{
	float bitvol = 0;
	unsigned numFilled = 0;
	unsigned andpat = 0xff, orpat = 0;
	for(unsigned i = 0; i < N; i++)
	{
		if(nodes[i].isLeaf)
		{
			numFilled++;
			andpat &= nodes[i].pattern;
			orpat |= nodes[i].pattern;

			if(nodes[i].pattern != 0 && bitvol == 0)
			{
				//compute volume of a single bit pattern
				bitvol = nodes[i].volume()/nodes[i].index;
			}
		}
	}

	if(isUnion && (numFilled == N || orpat == 0xff)) //just return a patterned node
	{
		unsigned nb = __builtin_popcount(orpat);
		MChildNode ret(true, orpat,nb, bitvol*nb);
		return ret;
	}
	else if(!isUnion && (numFilled == N || andpat == 0))
	{
		unsigned nb = __builtin_popcount(andpat);
		MChildNode ret(true, andpat,nb, bitvol*nb);
		return ret;
	}
	else //need to look at non-leaf children
	{
		unsigned pos = newtree.size();
		MChildNode ret(false, 0, pos, 0);
		newtree.push_back(MOctNode());

		unsigned numFullEmpty = 0;
		unsigned newPattern = 0;
		MChildNode nextnodes[N+1];
		const MappableOctTree *nexttrees[N+1];

		//setup nonleaf trees
		unsigned n = 0;
		for(unsigned i = 0; i < N; i++)
		{
			if(!nodes[i].isLeaf)
			{
				nexttrees[n] = trees[i];
				n++;
			}
		}
		if(numFilled > 0)
			nexttrees[n++] = NULL; //dummy for merge of leaves
		for (unsigned i = 0; i < 8; i++)
		{
			//create array of sub-octants from non-leafs
			unsigned cnt = 0;
			for(unsigned j = 0; j < N; j++)
			{
				if(!nodes[j].isLeaf)
				{
					nextnodes[cnt] = trees[j]->tree[nodes[j].index].children[i];
					cnt++;
				}
			}
			//summerize filled nodes
			if(numFilled > 0)
			{
				unsigned pat = isUnion ? orpat : andpat;
				if(pat & (1<<i))
					pat = 0xff;
				else
					pat = 0;
				nextnodes[cnt] = MChildNode(true, pat, pat ? 8 : 0, pat ? bitvol : 0);
				cnt++;
			}
			assert(cnt == n);

			MChildNode nchild = createFrom_r(n, nextnodes, nexttrees, newtree, isUnion);
			newtree[pos].children[i] = nchild;
			if(nchild.isLeaf && (nchild.pattern == 0 || nchild.pattern == 0xff))
			{
				numFullEmpty++;
				if(nchild.pattern == 0xff)
					newPattern |= 1<<i;
			}

			ret.vol += nchild.volume();
		}

		if(numFullEmpty == 8)
		{
			newtree.pop_back();
			ret.isLeaf = true;
			ret.pattern = newPattern;
			ret.index = __builtin_popcount(newPattern);
			//volume is correct
		}
		return ret;
	}
}

//invert tree, this can be done inplace
void MappableOctTree::invert(float dim)
{
	float expectedV = dim*dim*dim-root.volume();
	root.invert(tree, dim*dim*dim);
	assert(root.volume() == expectedV);
}

MappableOctTree* MappableOctTree::newFromVector(const vector<MOctNode>& newtree, const MChildNode& newroot)
{
	unsigned sz = sizeof(MappableOctTree)+newtree.size()*sizeof(MOctNode);
	void *mem = malloc(sz);
	return new (mem) MappableOctTree(newroot, newtree);
}

//union
MappableOctTree* MappableOctTree::createFromUnion(unsigned N, const MappableOctTree** in)
{
	MChildNode roots[N];
	for(unsigned i = 0; i < N; i++)
	{
		roots[i] = in[i]->root;
	}
	vector<MOctNode> newtree;
	MChildNode newroot = createFrom_r(N, roots, in,  newtree, true);

	return newFromVector(newtree, newroot);
}

//mutually recurse the N trees, round up/down nodes iff doing so will not
//change the local equality relationship between the N trees
bool MappableOctTree::createRoundedSet_r(unsigned N, MChildNode* nodes, const MappableOctTree** trees, bool roundUp,
			vector< vector<MOctNode> >& newtrees, vector<MChildNode>& newroots)
{
	unsigned numFilled = 0;
	unsigned numFull = 0;
	unsigned numEmpty = 0;
	int lastPat = -1;
	bool differentPats = false;
	float bitvol = 0;
	for(unsigned i = 0; i < N; i++)
	{
		if(nodes[i].isLeaf)
		{
			numFilled++;

			if(nodes[i].pattern == 0xff)
				numFull++;
			else if(nodes[i].pattern == 0)
				numEmpty++;
			else if(lastPat == -1)
			{
				lastPat = nodes[i].pattern;
				bitvol = nodes[i].vol/nodes[i].index;
			}
			else if(nodes[i].pattern != lastPat)
				differentPats = true;
		}
	}

	if (numFilled == N) //all leaves
	{
		if (roundUp)
		{
			//can only round up equivalent leaves if remaining leaves are empty
			if (!differentPats && lastPat > 0 && numFull == 0)
			{
				for(unsigned i = 0; i < N; i++)
				{
					if(nodes[i].pattern == 0)
						newroots[i] = nodes[i];
					else
						newroots[i] = MChildNode(true, 0xff, 8, bitvol*8);
				}
				return true;
			}
			else //just copy
			{
				for(unsigned i = 0; i < N; i++)
				{
					newroots[i] = nodes[i];
				}
				return false;
			}
		}
		else //round down
		{
			//leaves must be equivalent or full
			if(!differentPats && lastPat > 0 && numEmpty == 0)
			{
				for(unsigned i = 0; i < N; i++)
				{
					if(nodes[i].pattern == 0)
						newroots[i] = nodes[i];
					else
						newroots[i] = MChildNode(true, 0, 0, 0);
				}
				return true;
			}
			else //just duplicate
			{
				for(unsigned i = 0; i < N; i++)
				{
					newroots[i] = nodes[i];
				}
				return false;
			}
		}
	}
	else //must descend non-leaf tree parts
	{
		//create new octnodes for all non-leaf trees
		unsigned positions[N];
		for(unsigned i = 0; i < N; i++)
		{
			if(nodes[i].isLeaf)
			{
				//stays a leaf
				newroots[i] = nodes[i];
			}
			else
			{
				unsigned pos = newtrees[i].size();
				positions[i] = pos;
				newroots[i] = MChildNode(false, 0, pos, 0);
				newtrees[i].push_back(MOctNode());
			}
		}

		MChildNode subnodes[N];
		bool changed = false;
		unsigned leafCnts[N];
		unsigned leafPatterns[N];
		memset(leafCnts, 0, sizeof(leafCnts));
		memset(leafPatterns, 0, sizeof(leafPatterns));

		for(unsigned o = 0; o < 8; o++)
		{
			//setup subnodes, divide leaf nodes
			for(unsigned i = 0; i < N; i++)
			{
				if(nodes[i].isLeaf)
				{
					if(nodes[i].pattern & (1<<o))
						subnodes[i] = MChildNode(true, 0xff, 8, nodes[i].vol/nodes[i].index);
					else
						subnodes[i] = MChildNode(true, 0, 0, 0);
				}
				else
				{
					subnodes[i] = trees[i]->tree[nodes[i].index].children[o];
				}
			}

			vector<MChildNode> newnodes(N);
			changed |= createRoundedSet_r(N, subnodes, trees, roundUp, newtrees, newnodes);

			for(unsigned i = 0; i < N; i++)
			{
				if(!nodes[i].isLeaf)
				{
					newtrees[i][positions[i]].children[o] = newnodes[i];
					if(newnodes[i].isLeaf)
					{
						if(newnodes[i].pattern == 0xff)
						{
							leafCnts[i]++;
							leafPatterns[i] |= (1<<o);
						}
						else if(newnodes[i].pattern == 0)
							leafCnts[i]++;
					}
					newroots[i].vol += newnodes[i].vol;
				}
			}
		}

		//coalesce any all leaves
		for(unsigned i = 0; i < N; i++)
		{
			if(leafCnts[i] == 8)
			{
				newtrees[i].pop_back();
				newroots[i].isLeaf = true;
				newroots[i].pattern = leafPatterns[i];
				newroots[i].index = __builtin_popcount(leafPatterns[i]);
				//volume is correct
			}
		}

		return changed;
	}
}

bool MappableOctTree::createRoundedSet(unsigned N, const MappableOctTree**in, bool roundUp, MappableOctTree** out)
{
	MChildNode roots[N];
	for(unsigned i = 0; i < N; i++)
	{
		roots[i] = in[i]->root;
	}
	vector< vector<MOctNode> > newtrees(N);
	vector<MChildNode>  newroots(N);
	bool changed = createRoundedSet_r(N, roots, in,  roundUp, newtrees, newroots);

	for(unsigned i = 0; i < N; i++)
	{
		out[i] = newFromVector(newtrees[i], newroots[i]);
	}

	return changed;
}


MappableOctTree* MappableOctTree::createFromIntersection(unsigned N, const MappableOctTree** in)
{
	MChildNode roots[N];
	for(unsigned i = 0; i < N; i++)
	{
		roots[i] = in[i]->root;
	}
	vector<MOctNode> newtree;
	MChildNode newroot = createFrom_r(N, roots, in,  newtree, false);

	return newFromVector(newtree, newroot);
}

MChildNode MappableOctTree::createTruncated_r(const MChildNode& node, float dim, float res, bool fillIn, vector<MOctNode>& newtree) const
{
	if(dim <= res) //at level that needs to be truncated to a leaf
	{
		if(fillIn && volume() != 0)
		{
			return MChildNode(true, 0xff, 8, dim*dim*dim);
		}
		else
		{
			return MChildNode(true, 0, 0, 0);
		}
	}
	else if(node.isLeaf)
	{
		return node;
	}
	else //recursively descend
	{
		unsigned pos = newtree.size();
		MChildNode ret(false, 0, pos, 0);
		newtree.push_back(MOctNode());

		unsigned numFullEmpty = 0;
		unsigned newPattern = 0;
		for (unsigned i = 0; i < 8; i++)
		{
			MChildNode nchild = createTruncated_r(tree[node.index].children[i] , dim/2, res, fillIn, newtree);
			newtree[pos].children[i] = nchild;
			if(nchild.isLeaf && (nchild.pattern == 0 || nchild.pattern == 0xff))
			{
				numFullEmpty++;
				if(nchild.pattern == 0xff)
					newPattern |= 1<<i;
			}

			ret.vol += nchild.volume();
		}

		if(numFullEmpty == 8)
		{
			newtree.pop_back();
			ret.isLeaf = true;
			ret.pattern = newPattern;
			ret.index = __builtin_popcount(newPattern);
			//volume is correct
		}
		return ret;
	}
}

MappableOctTree* MappableOctTree::createTruncated(float dim, float res, bool fillIn) const
{
	vector<MOctNode> newtree;
	MChildNode newroot = createTruncated_r(root, dim, res, fillIn, newtree);

	return newFromVector(newtree, newroot);
}

MChildNode MappableOctTree::createRounded_r(const MChildNode& node, float dim, float volcutoff, bool fillIn, vector<MOctNode>& newtree) const
{
	float maxvol = dim*dim*dim;
	float vol = node.volume();
	float empty = maxvol-vol;

	if(fillIn && empty < volcutoff)
	{
		return MChildNode(true, 0xff, 8, maxvol);
	}
	else if(!fillIn && vol < volcutoff)
	{
		return MChildNode(true, 0, 0, 0);
	}
	else if(node.isLeaf)
	{
		return node;
	}
	else //recursively descend
	{
		unsigned pos = newtree.size();
		MChildNode ret(false, 0, pos, 0);
		newtree.push_back(MOctNode());

		unsigned numFullEmpty = 0;
		unsigned newPattern = 0;
		for (unsigned i = 0; i < 8; i++)
		{
			MChildNode nchild = createRounded_r(tree[node.index].children[i] , dim/2, volcutoff, fillIn, newtree);
			newtree[pos].children[i] = nchild;
			if(nchild.isLeaf && (nchild.pattern == 0 || nchild.pattern == 0xff))
			{
				numFullEmpty++;
				if(nchild.pattern == 0xff)
					newPattern |= 1<<i;
			}

			ret.vol += nchild.volume();
		}

		if(numFullEmpty == 8)
		{
			newtree.pop_back();
			ret.isLeaf = true;
			ret.pattern = newPattern;
			ret.index = __builtin_popcount(newPattern);
			//volume is correct
		}
		return ret;
	}
}

MappableOctTree* MappableOctTree::createRounded(float dim, float vol, bool fillIn) const
{
	vector<MOctNode> newtree;
	MChildNode newroot = createRounded_r(root, dim, vol, fillIn, newtree);

	return newFromVector(newtree, newroot);
}


//volume calculations that don't require creating a tmp tree
float MappableOctTree::intersectVolume(const MappableOctTree * rhs) const
{
	float ival = 0, uval = 0;
	root.intersectUnionVolume(tree,  rhs->root, rhs->tree, ival, uval);
	return ival;
}

float MappableOctTree::unionVolume(const MappableOctTree *rhs) const
{
	float ival = 0, uval = 0;
	root.intersectUnionVolume(tree,  rhs->root, rhs->tree, ival, uval);
	return uval;
}


bool MappableOctTree::containedIn(const MappableOctTree *rhs) const
{
	return root.containedIn(tree, rhs->root, rhs->tree);
}

//return total volume contained in octtree
float MappableOctTree::volume() const
{
	return root.volume();
}

//return number of leaves
unsigned MappableOctTree::leaves() const
{
	unsigned cnt = 0;
	for (unsigned i = 0; i < N; i++)
	{
		for (unsigned j = 0; j < 8; j++)
		{
			if (tree[i].children[j].isLeaf)
				cnt++;
		}
	}
	return cnt;
}


void MappableOctTree::write(ostream& out) const
{
	unsigned sz = sizeof(MappableOctTree)+N*sizeof(MOctNode);
	out.write((char*)this, sz);
}

float MappableOctTree::relativeVolumeDistance(const MappableOctTree * rhs) const
{
	float ival = 0, uval = 0;
	root.intersectUnionVolume(tree,  rhs->root, rhs->tree, ival, uval);

	return 1 - ival/uval;
}

float MappableOctTree::absoluteVolumeDistance(const MappableOctTree * rhs) const
{
	float ival = 0, uval = 0;
	root.intersectUnionVolume(tree, rhs->root,rhs->tree, ival, uval);

	return uval - ival;
}



void MChildNode::invert(MOctNode* tree, float maxvol)
{
	if(isLeaf)
	{
		pattern = ~pattern;
		index = 8 - index;
		vol = maxvol - vol;
	}
	else
	{
		vol = 0;
		for (unsigned i = 0; i < 8; i++)
		{
			tree[index].children[i].invert(tree, maxvol/8);
			vol += tree[index].children[i].vol;
		}
	}
}


//compute both the intersection and union volume at once
void MChildNode::intersectUnionVolume(
		const MOctNode* tree, const MChildNode& rhs, const MOctNode* rtree, float& intersectval, float& unionval) const

{
	if(rhs.isLeaf && isLeaf)
	{
		if(pattern == 0 || rhs.pattern == 0)
			intersectval += 0; //no intersect val
		else
		{
			float bitv = vol/index; //volume of each bit
			unsigned intpat = pattern & rhs.pattern;

			intersectval += __builtin_popcount(intpat)*bitv;
		}
		if(pattern == 0)
			unionval += rhs.volume();
		else if(rhs.pattern == 0)
			unionval += volume();
		else
		{
			float bitv = vol/index; //volume of each bit
			unsigned orpat = pattern | rhs.pattern;

			unionval += __builtin_popcount(orpat)*bitv;
		}
	}
	else if(isLeaf && pattern == 0xff)
	{
		unionval += volume();
		intersectval += rhs.volume();
	}
	else if (rhs.isLeaf && rhs.pattern == 0xff)
	{
		unionval += rhs.volume();
		intersectval += volume();
	}
	else if (rhs.isLeaf && rhs.pattern == 0)
	{
		//vol of this tree
		unionval += volume();
	}
	else if (isLeaf && pattern == 0)
	{
		// from rhs
		unionval += rhs.volume();
	}
	else if(isLeaf)
	{
		assert(!rhs.isLeaf);
		assert(index > 0);
		//rhs has children, have to compare bit by bit
		float bvol = vol/index;
		for(unsigned i = 0; i < 8; i++)
		{
			MChildNode rchild = rtree[rhs.index].children[i];
			if(pattern & (1<<i))
			{
				unionval += bvol;
				intersectval += rchild.volume();
			}
			else
			{
				unionval += rchild.volume();
			}
		}
	}
	else if(rhs.isLeaf)
	{
		assert(!isLeaf);
		assert(rhs.index > 0);
		//rhs has children, have to compare bit by bit
		float bvol = rhs.vol/rhs.index;
		for(unsigned i = 0; i < 8; i++)
		{
			MChildNode child = tree[index].children[i];
			if(rhs.pattern & (1<<i))
			{
				unionval += bvol;
				intersectval += child.volume();
			}
			else
			{
				unionval += child.volume();
			}
		}
	}
	else //both have children
	{
		for (unsigned i = 0; i < 8; i++)
		{
			tree[index].children[i].intersectUnionVolume(tree,
					rtree[rhs.index].children[i], rtree, intersectval,unionval);
		}
	}
}

//return true if this is contained in rhs - short circuit eval
bool MChildNode::containedIn(
		const MOctNode* tree, const MChildNode& rhs, const MOctNode* rtree) const
{
	if(rhs.isLeaf && isLeaf)
	{
		return (pattern & rhs.pattern) == pattern;
	}
	else if(isLeaf && pattern == 0)
	{
		return true;
	}
	else if (rhs.isLeaf && rhs.pattern == 0)
	{
		return false; //somthing not in nothing (nothing handled above)
	}
	else if (rhs.isLeaf && rhs.pattern == 0xff)
	{
		return true;
	}
	else if (isLeaf && pattern == 0xff)
	{
		return false;
	}
	else if(isLeaf)
	{
		assert(!rhs.isLeaf);
		assert(index > 0);
		//rhs has children, have to compare bit by bit
		for(unsigned i = 0; i < 8; i++)
		{
			if(pattern & (1<<i))
			{
				MChildNode rchild = rtree[rhs.index].children[i];
				//rchild must be a leaf and full
				if(!rchild.isLeaf || rchild.pattern != 0xff)
					return false;
			}
		}
		return true;
	}
	else if(rhs.isLeaf)
	{
		assert(!isLeaf);
		assert(rhs.index > 0);
		//rhs has children, have to compare bit by bit
		for(unsigned i = 0; i < 8; i++)
		{
			if(!(rhs.pattern & (1<<i)))
			{
				MChildNode child = tree[index].children[i];
				//if rhs's bit is empty, then our child must be empty
				if(!child.isLeaf || child.pattern != 0)
					return false;
			}
		}
		return true;
	}
	else //both have children
	{
		for (unsigned i = 0; i < 8; i++)
		{
			if(!tree[index].children[i].containedIn(tree,
					rtree[rhs.index].children[i], rtree))
				return false;
		}
		return true;;
	}
	return false;
}

static void addVertices(vector<Vertex>& vertices, Cube box)
{
	for(unsigned i = 0; i < 8; i++)
	{
		float x = box.x;
		float y = box.y;
		float z = box.z;

		if(i & 1)
			x += box.dim;
		if(i & 2)
			y += box.dim;
		if(i & 4)
			z += box.dim;

		vertices.push_back(Vertex(x,y,z));
	}
}

//append vertices of this child to vertices, this child's box is provided
void MChildNode::collectVertices(vector<Vertex>& vertices, Cube box, const MOctNode* tree) const
{
	if(isLeaf)
	{
		if(pattern == 0)
			return; //empty, no vertices
		else if(pattern == 0xff) //full, add this box's vertices
		{
			addVertices(vertices, box);
		}
		else
		{
			for(unsigned i = 0; i < 8; i++)
			{
				if(pattern&(1<<i))
				{
					addVertices(vertices, box.getOctant(i));
				}
			}
		}
	}
	else
	{
		for (unsigned i = 0; i < 8; i++)
		{
			tree[index].children[i].collectVertices(vertices, box.getOctant(i), tree);
		}
	}
}

//sort, then remove duplicates and remove vertices that appear 8 times altogether
static void uniqueRemove8ify(vector<Vertex>& V)
{
	sort(V.begin(), V.end());

	unsigned i = 1;
	unsigned n = V.size();
	unsigned cnt = 1;
	unsigned last = 0;
	while(i < n)
	{
		if(V[i] == V[last])
		{
			cnt++;
			i++;
		}
		else
		{
			if(cnt < 8)
				last++; //keep what we copied
			V[last] = V[i];
			cnt = 1;
		}
	}

	if(cnt < 8)
		V.resize(last+1);
	else
		V.resize(last);
}

//the Hausdorff distance is the maximum minimum distance between the
//two shapes; it is not symmetric, but this computes the symmetric
//(max H(A,B),H(B,A)) distance
//
//TODO: make this efficient if it turns out to be worthwhile, right now
//do the easiest implimentation: collect all leaf corner coordinates, sort,
//remove completely buried (8 versions), compute all pairs distances
float MappableOctTree::hausdorffDistance(const MappableOctTree* B, float dim) const
{
	vector<Vertex> Av, Bv;

	Cube box(dim, -dim / 2, -dim / 2, -dim / 2);
	root.collectVertices(Av, box, tree);
	B->root.collectVertices(Bv, box, B->tree);

	uniqueRemove8ify(Av);
	uniqueRemove8ify(Bv);

	//compute distances once
	vector<float> Bmin(Bv.size(), HUGE_VAL);
	float Amax = 0;
	for(unsigned i = 0, n = Av.size(); i < n; i++)
	{
		float Amin = HUGE_VAL;
		for(unsigned j = 0, m = Bv.size(); j < m; j++)
		{
			float dist = Av[i].distanceSq(Bv[j]);
			if(dist < Amin)
				Amin = dist;

			if(dist < Bmin[j])
				Bmin[j] = dist;
		}
		if(Amin > Amax)
			Amax = Amin;
	}

	float Bmax = *std::max_element(Bmin.begin(), Bmin.end());
//cout << Av.size() << " " << oldA << " " << Bv.size() << " " << oldB << " " << Amax << " " << Bmax << "\n";
	return sqrt(std::max(Amax, Bmax));
}

bool MChildNode::equals(const MOctNode* tree, const MChildNode& rhs, const MOctNode* rtree) const
{
	if(isLeaf)
	{
		if(rhs.isLeaf)
		{
			if(rhs.pattern == pattern)
				return true;
			else
				return false;
		}
		else
			return false;
	}
	else if(rhs.isLeaf)
	{
		return false;
	}
	else
	{
		for (unsigned i = 0; i < 8; i++)
		{
			if(!tree[index].children[i].equals(tree,
					rtree[rhs.index].children[i], rtree))
				return false;
		}
		return true;
	}
	return true;
}

bool MappableOctTree::equals(const MappableOctTree* rhs) const
{
	return root.equals(tree, rhs->root, rhs->tree);
}

