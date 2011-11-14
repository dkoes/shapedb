#include "GSSTreeStructures.h"
#include "DataViewers.h"

//write out a leaf node to outNodes and the representative trees to outTrees
void GSSLeaf::writeLeaf(const DataViewer *data, const Cluster& cluster, ostream& outNodes, ostream& outTrees)
{
	GSSNodeCommon info;
	info.isLeaf = true;
	info.N = cluster.size();

	unsigned positions[info.N];

	outNodes.write((char*)&info, sizeof(info));

	//compute the position offsets
	unsigned curoffset = info.N*sizeof(unsigned); //offsets are from start of data section
	for(unsigned i = 0, n = cluster.size(); i < n; i++)
	{
		unsigned ind = cluster[i];
		positions[i] = curoffset;
		const MappableOctTree *tree = data->getMIV(ind);
		curoffset += sizeof(file_index);
		curoffset += tree->bytes();
	}

	//write positions
	outNodes.write((char*)positions, sizeof(unsigned)*info.N);

	//write tree info
	for(unsigned i = 0, n = cluster.size(); i < n; i++)
	{
		unsigned ind = cluster[i];
		//for leaf, write out object file index
		file_index index = data->getIndex(ind);
		outNodes.write((char*)&index, sizeof(file_index));

		//and single tree
		const MappableOctTree *tree = data->getMIV(ind);
		tree->write(outNodes);
	}

	//write out trees of cluster to a GSSDoubleTree
	unsigned offset = cluster.MIV->bytes();

	outTrees.write((char*)&offset, sizeof(offset));
	cluster.MIV->write(outTrees);
	cluster.MSV->write(outTrees);

}

//very similar to write leaf, but write out MIV/MSV of subnodes
void GSSInternalNode::writeNode(const DataViewer *data, const Cluster& cluster, ostream& outNodes, ostream& outTrees)
{
	GSSNodeCommon info;
	info.isLeaf = false;
	info.N = cluster.size();

	unsigned positions[info.N];
	outNodes.write((char*)&info, sizeof(info));

	//compute the position offsets
	unsigned curoffset = info.N*sizeof(unsigned); //offsets are from start of data section
	for(unsigned i = 0, n = cluster.size(); i < n; i++)
	{
		unsigned ind = cluster[i];
		positions[i] = curoffset;
		const MappableOctTree *miv = data->getMIV(ind);
		const MappableOctTree *msv = data->getMSV(ind);
		curoffset += sizeof(Child);
		curoffset += miv->bytes();
		curoffset += msv->bytes();
	}

	//write positions
	outNodes.write((char*)positions, sizeof(unsigned)*info.N);

	//write tree info
	for(unsigned i = 0, n = cluster.size(); i < n; i++)
	{
		unsigned ind = cluster[i];
		const MappableOctTree *miv = data->getMIV(ind);
		const MappableOctTree *msv = data->getMSV(ind);

		Child child;
		child.MSVindex = miv->bytes();
		child.node_pos = data->getIndex(ind);
		outNodes.write((char*)&child, sizeof(Child));

		miv->write(outNodes);
		msv->write(outNodes);
	}

	//write out trees of cluster to a GSSDoubleTree
	unsigned offset = cluster.MIV->bytes();

	outTrees.write((char*)&offset, sizeof(offset));
	cluster.MIV->write(outTrees);
	cluster.MSV->write(outTrees);

}

unsigned GSSLeaf::bytes() const
{
	unsigned ret = sizeof(GSSLeaf);
	if(info.N == 0) return ret;
	ret += child_positions[info.N-1]; //offset to last child
	const Child* child = (const Child*)&data[child_positions[info.N-1]];
	//add size of last child
	ret += sizeof(file_index);
	ret += child->tree.bytes();
	return ret;
}

unsigned GSSInternalNode::Child::bytes() const
{
	unsigned ret = sizeof(Child);
	ret += MSVindex;
	const MappableOctTree* msv = getMSV();
	ret += msv->bytes();

	return ret;
}

unsigned GSSInternalNode::bytes() const
{
	unsigned ret = sizeof(GSSInternalNode);
	if(info.N == 0) return ret;
	ret += child_positions[info.N-1]; //offset to last child
	const Child* child = (const Child*)&data[child_positions[info.N-1]];
	//add size of last child
	ret += child->bytes();
	return ret;
}

void GSSInternalNode::setChildPos(unsigned i, file_index newpos, bool isLeaf, file_index lstart, file_index lend)
{
	Child *child =  (Child*)&data[child_positions[i]];
	child->node_pos = newpos;
	child->isLeaf = isLeaf;
	child->leaves_start = lstart;
	child->leaves_end = lend;
}

//malloc a new node that has reduced sized trees
GSSInternalNode* GSSInternalNode::createTruncated(float dimension, float resolution)
{
	unsigned nc = size();
	MappableOctTree *MIVs[nc];
	MappableOctTree *MSVs[nc];
	unsigned positions[nc];

	unsigned curoffset = nc*sizeof(unsigned); //skip of positions
	for(unsigned c = 0; c < nc; c++)
	{
		const Child *child = getChild(c);
		MIVs[c] = child->getMIV()->createTruncated(dimension, resolution*2, false);
		MSVs[c] = child->getMSV()->createTruncated(dimension, resolution*2, true);

		positions[c] = curoffset;

		curoffset += sizeof(Child);
		curoffset += MIVs[c]->bytes();
		curoffset += MSVs[c]->bytes();
	}

	unsigned char* buffer = (unsigned char*)malloc(curoffset+sizeof(GSSInternalNode));

	memcpy(buffer, this, sizeof(GSSInternalNode));
	unsigned offset = sizeof(GSSInternalNode);
	memcpy(buffer+offset, &positions, sizeof(unsigned)*nc);
	offset += sizeof(unsigned)*nc;

	for(unsigned c = 0; c < nc; c++)
	{
		Child child = *getChild(c);
		child.MSVindex = MIVs[c]->bytes();

		memcpy(buffer+offset, &child, sizeof(Child));
		offset += sizeof(Child);
		memcpy(buffer+offset, MIVs[c], MIVs[c]->bytes());
		offset += MIVs[c]->bytes();
		memcpy(buffer+offset, MSVs[c], MSVs[c]->bytes());
		offset += MSVs[c]->bytes();
	}

	assert(offset == curoffset+sizeof(GSSInternalNode));
	return (GSSInternalNode*)buffer;
}

