/*
 * GSSTree.cpp
 *
 *  Created on: Sep 12, 2011
 *      Author: dkoes
 */

#include "GSSTree.h"
#include "OctTree.h"
#include <math.h>
#include <iostream>


const unsigned GSSTree::MaxSplit = 8; //number of children in each node

GSSTree::GSSNode::~GSSNode() {}


void GSSTree::setBoundingBox(array<float,6>& box)
{
	float dims[3] = {0,0,0};
	for(unsigned i = 0; i < 3; i++)
	{
		min[i] = box[2*i];
		dims[i] = box[2*i+1]-min[i];
	}

	//to make things easier, cubify using the longest dimension
	dim = max(dims[2],max(dims[0],dims[1]));
	//round all dimensions up to a power of 2; use fp
	dim = pow(2, ceil(log(dim)/log(2)));
	if(dim < 1)
		dim = 1;
}

// perform any transformations necessary to conver the input coordinate system
//into the gss tree coordinate system
void GSSTree::transformMol(const vector<MolSphere>& m, vector<MolSphere>& ret)
{
	ret.clear();
	ret.reserve(m.size());

	for(unsigned i = 0, n = m.size(); i < n; i++)
	{
		ret.push_back(MolSphere(m[i].x-min[0], m[i].y-min[1], m[i].z-min[2],m[i].r));
	}
}

//add a single mol
void GSSTree::add(const vector<MolSphere>& m)
{
	//reposition mol
	vector<MolSphere> mol;
	transformMol(m, mol);
	OctTree oct(dim, maxres, mol);

	//insert into the tree - first find the leaf node to insert into, then
	//performe the insertion

	float dist = HUGE_VAL;
	GSSLeafNode *leaf = NULL;
	root->findInsertionPoint(oct, dist, leaf);
	assert(leaf != NULL);
	leaf->insert(*this, oct, LeafData(m));
}

//nearest neighbor search, return closest set of molspheres
void GSSTree::nn_search(const vector<MolSphere>& m, vector<MolSphere>& res)
{
	vector<MolSphere> mol;
	transformMol(m, mol);
	OctTree tree(dim, maxres, mol);
	LeafData data;
	float dist = HUGE_VAL;
	root->findNearest(tree, dist, data);
	swap(res, data.spheres);
}


void GSSTree::write(const filesystem::path& out)
{
	cerr << "Unimplemented: write\n";
	exit(-1);
}

//insert data into leaf node, if there isn't enough room, split
//call update on parents if need-be
void GSSTree::GSSLeafNode::insert(GSSTree& gTree, const OctTree& tree, const LeafData& m)
{
	if(data.size() < MaxSplit)
	{
		bool needUpdate = false;

		assert(data.size() == trees.size());
		//easy case, just add
		trees.push_back(new OctTree(tree));
		data.push_back(LeafData(m));

		//update MIV/MSV
		needUpdate |= MIV->intersect(tree);
		needUpdate |= MSV->unionWith(tree);

		if(needUpdate && parent != NULL)
		{
			parent->update(gTree, which, NULL);
		}
	}
	else //must split
	{
		//first add tree
		trees.push_back(new OctTree(tree));
		data.push_back(LeafData(m));

		//so that we can have one split function, have concept of MIV and MSV
		//even at leaf
		vector<unsigned> split1;
		vector<unsigned> split2;
		split(trees, trees, split1, split2);
		GSSLeafNode *newnode = new GSSLeafNode(gTree.dim, gTree.maxres);

		//put split2 in new node
		newnode->MIV->fill();
		newnode->data.resize(split2.size());
		for(unsigned i = 0, n = split2.size(); i < n; i++)
		{
			unsigned index = split2[i];
			newnode->trees.push_back(trees[index]);
			swap(newnode->data[i],data[index]);

			//update MIV/MSV
			newnode->MSV->unionWith(*trees[index]);
			newnode->MIV->intersect(*trees[index]);
		}
		//now mogirfy to split1
		vector<OctTree*> tmptrees; tmptrees.reserve(MaxSplit);
		vector<LeafData> tmpdata(split1.size());

		MIV->fill();
		MSV->clear();
		for(unsigned i = 0, n = split1.size(); i < n; i++)
		{
			unsigned index = split1[i];
			tmptrees.push_back(trees[index]);
			swap(tmpdata[i],data[index]);

			MIV->intersect(*trees[index]);
			MSV->unionWith(*trees[index]);
		}

		swap(tmpdata,data);
		swap(tmptrees,trees);

		if(parent != NULL) parent->update(gTree, which, newnode);
		else if(newnode != NULL) //need new root
		{
			gTree.createRoot(this, newnode);
		}
	}

}

//add a child without any checks, just get the pointers right
void GSSTree::GSSInternalNode::addChild(GSSNode *child)
{
	child->parent = this;
	child->which = children.size();
	children.push_back(child);
}

//recurse UP the tree updating MIV/MSV and splitting as necessary
void GSSTree::GSSInternalNode::update(GSSTree& tree, unsigned whichChild, GSSNode *newnode)
{
	//first update MIV/MSV
	assert(whichChild < children.size());

	bool needUpdate = false;
	needUpdate |= MIV->intersect(*children[whichChild]->MIV);
	needUpdate |= MSV->unionWith(*children[whichChild]->MSV);

	if(newnode != NULL)
	{
		//add node to this level
		if(children.size() < MaxSplit)
		{
			//easy case
			needUpdate |= MIV->intersect(*newnode->MIV);
			needUpdate |= MSV->unionWith(*newnode->MSV);
			addChild(newnode);
			newnode = NULL;
		}
		else //must split
		{
			//setup split structures
			vector<OctTree*> MIVs, MSVs;
			MIVs.reserve(MaxSplit);
			MSVs.reserve(MaxSplit);

			//add in new node briefly so it's part of the split
			children.push_back(newnode);

			for(unsigned i = 0, n = children.size(); i < n; i++)
			{
				MIVs.push_back(children[i]->MIV);
				MSVs.push_back(children[i]->MSV);
			}

			vector<unsigned> split1, split2;

			//do the split
			split(MIVs, MSVs, split1, split2);

			//save oldchildren
			vector<GSSNode*> oldchildren;
			oldchildren.reserve(MaxSplit);
			swap(oldchildren, children);

			//mogify this node to hold split1
			MIV->fill();
			MSV->clear();
			for(unsigned i = 0, n = split1.size(); i < n; i++)
			{
				GSSNode *child = oldchildren[split1[i]];
				addChild(child);
				MIV->intersect(*child->MIV);
				MSV->unionWith(*child->MSV);
			}

			//same thing for the new node and split2
			GSSInternalNode *newinode = new GSSInternalNode(tree.dim,tree.maxres);
			newinode->MIV->fill();
			newinode->MSV->clear();
			for(unsigned i = 0, n = split2.size(); i < n; i++)
			{
				GSSNode *child = oldchildren[split2[i]];
				newinode->addChild(child);
				newinode->MIV->intersect(*child->MIV);
				newinode->MSV->unionWith(*child->MSV);
			}
			newnode = newinode;
		}
	}

	if(parent != NULL)
	{
		parent->update(tree, which, newnode);
	}
	else if(newnode != NULL && parent == NULL) //need new root
	{
		tree.createRoot(this, newnode);
	}

}

//create a new root using left and right, which are assumed to
//be splits of the current root
void GSSTree::createRoot(GSSNode *left, GSSNode *right)
{
	assert(left == root);
	GSSInternalNode* newroot = new GSSInternalNode(dim, maxres);
	//update trees
	*newroot->MIV = *left->MIV;
	newroot->MIV->intersect(*right->MIV);

	*newroot->MSV = *left->MSV;
	newroot->MSV->unionWith(*right->MSV);

	//add two children
	newroot->addChild(left);
	newroot->addChild(right);

	root = newroot;
}


//return a distance to a single leaf, should be compatible with search Dist
float GSSTree::leafDist(const OctTree* obj, const OctTree *leaf)
{
	return 1 - obj->intersectVolume(*leaf)/obj->unionVolume(*leaf);
}

//return a "distance" between obj and MIV/MSV; the lower the distance
//the higher the more likely a node should be searched
//min and max should bookend the ultimate leaf distances
float GSSTree::searchDist(const OctTree* obj, const OctTree *MIV, const OctTree *MSV, float& min, float& max)
{
	min = 1 - obj->intersectVolume(*MSV)/obj->unionVolume(*MIV);
	max = 1 - obj->intersectVolume(*MIV)/obj->unionVolume(*MSV);

	return min+max;
}


/*
 * return a "distance" between approximations
 * for now, just the sum of the volume different of the MIVs and of the MSVs
 */
float GSSTree::splitDist(OctTree* leftMIV, OctTree* leftMSV, OctTree* rightMIV, OctTree* rightMSV)
{
	float mind = leftMIV->intersectVolume(*rightMIV)/leftMIV->unionVolume(*rightMIV);
	float maxd = leftMSV->intersectVolume(*rightMSV)/leftMSV->unionVolume(*rightMSV);

	return mind+maxd;
}


/* Decide how to split a node.  For now, use the canonical quadratic split
 * where you find the two most "distant" nodes, and then greedily partition
 * using those.
 */
void GSSTree::split(const vector<OctTree*>& MIV, const vector<OctTree*>& MSV,
		vector<unsigned>& s1, vector<unsigned>& s2)
{
	assert(MIV.size() == MSV.size());
	s1.clear();
	s2.clear();
	unsigned besti = 0;
	unsigned bestj = 0;
	float dist = 0;

	unsigned N = MIV.size();
	s1.reserve(N);
	s2.reserve(N);
	float distances[N][N];

	for(unsigned i = 0; i < N; i++)
	{
		for(unsigned j = i+1; j < N; j++)
		{
			float d = splitDist(MIV[i],MSV[i],MIV[j],MSV[j]);
			distances[i][j] = d;
			distances[j][i] = d;
			if(d > dist)
			{
				dist = d;
				besti = i;
				bestj = j;
			}
		}
		distances[i][i] = 0;
	}
	assert(besti != bestj);

	//seed with besti and bestj
	OctTree iMIV = *MIV[besti];
	OctTree iMSV = *MSV[besti];
	s1.push_back(besti);

	OctTree jMIV = *MIV[bestj];
	OctTree jMSV = *MSV[bestj];
	s2.push_back(bestj);

	for(unsigned i = 0; i < N; i++)
	{
		if(i != besti && i != bestj)
		{
			//choose based on distance to cumulative MIV/MSV
			//unless we've already filled one split
			if(s1.size() > N/2+1)
			{
				s2.push_back(i);
			}
			else if(s2.size() > N/2+1)
			{
				s1.push_back(i);
			}
			else
			{
				float di = splitDist(&iMIV,&iMSV,MIV[i],MSV[i]);
				float dj = splitDist(&jMIV, &jMSV, MIV[i],MSV[i]);

				if(di < dj)
				{
					s1.push_back(i);
					iMIV.intersect(*MIV[i]);
					iMSV.unionWith(*MSV[i]);
				}
				else if(dj < di)
				{
					s2.push_back(i);
					jMIV.intersect(*MIV[i]);
					jMSV.unionWith(*MSV[i]);
				}
				else //equal, seems unlikely
				{
					if(s1.size() < s2.size())
					{
						s1.push_back(i);
						iMIV.intersect(*MIV[i]);
						iMSV.unionWith(*MSV[i]);
					}
					else
					{
						s2.push_back(i);
						jMIV.intersect(*MIV[i]);
						jMSV.unionWith(*MSV[i]);
					}
				}

			}
		}
	}

}


//find the object with the best distance in this node, if it's better than
//the passed distance, update
void GSSTree::GSSLeafNode::findNearest(const OctTree& tree, float& distance, LeafData& d)
{
	for(unsigned i = 0, n = trees.size(); i < n; i++)
	{
		float dist = leafDist(&tree, trees[i]);
		if(dist < distance)
		{
			distance = dist;
			d = data[i];
		}
	}
}

//if this leaf is a more appropriate place for tree, return self
void GSSTree::GSSLeafNode::findInsertionPoint(const OctTree& tree, float& distance, GSSLeafNode*& leaf)
{
	if(trees.size() == 0)
	{
		distance = 0;
		leaf = this;
	}
	else
	{
		float min, max;
		float dist = searchDist(&tree, MIV, MSV, min, max);
		if(dist < distance)
		{
			distance = dist;
			leaf = this;
		}
	}
}


struct ScoreIndex
{
	unsigned index;
	float score;
	float min;

	ScoreIndex() {}
	ScoreIndex(unsigned i, float s, float m): index(i), score(s), min(m) {}

	bool operator<(const ScoreIndex& si) const
	{
		return score < si.score;
	}
};
//explore children to find closest value
void GSSTree::GSSInternalNode::findNearest(const OctTree& tree, float& distance, LeafData& data)
{
	vector<ScoreIndex> SIs;
	SIs.reserve(children.size());

	for(unsigned i = 0, n = children.size(); i < n; i++)
	{
		float min = 0, max = 0;
		float score = searchDist(&tree, children[i]->MIV, children[i]->MSV, min, max);
		if(min < distance) //there's hope
		{
			SIs.push_back(ScoreIndex(i,score,min));
		}
	}

	//explore in priority order
	sort(SIs.begin(), SIs.end());

	for(unsigned i = 0, n = SIs.size(); i < n; i++)
	{
		if(SIs[i].min < distance)
		{
			children[SIs[i].index]->findNearest(tree, distance, data);
		}
	}
}

//look for a good place for the data indexed by tree
void GSSTree::GSSInternalNode::findInsertionPoint(const OctTree& tree, float& distance, GSSLeafNode*& leaf)
{
	float min, max;
	float dist = searchDist(&tree, MIV, MSV, min, max);
	if(dist < distance)
	{
		//if this tree might contain something better, look at all children -> really should sort this
		for(unsigned i = 0, n = children.size(); i < n; i++)
		{
			children[i]->findInsertionPoint(tree, distance, leaf);
		}
	}
}


