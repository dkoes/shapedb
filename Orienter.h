/*
 * Orienter.h
 *
 *  Created on: Jul 29, 2014
 *      Author: dkoes
 *
 *  Maintains a translation vector and rotation matrix for reorienting coordinates.
 */

#ifndef ORIENTER_H_
#define ORIENTER_H_

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <rdkit/Geometry/point.h>
#include <iostream>

typedef Eigen::Matrix<double, Eigen::Dynamic, 3> ECoords;
typedef Eigen::Matrix3d EMatrix3; //3x3 double
typedef Eigen::Vector3d EVector3;

class Orienter
{
	EVector3 translate; //gets added to coords
	EMatrix3 rotate; //then multiply by this - actually transpose of rotation matrix for right multiply
public:
	Orienter() :
			translate(Eigen::Vector3d::Zero()), rotate(
					Eigen::Matrix3d::Identity())
	{

	}

	bool isIdentity() const
	{
		return translate == Eigen::Vector3d::Zero() && rotate == Eigen::Matrix3d::Identity();
	}

	//accumlate a translation vector
	void addTranslation(const EVector3& t)
	{
		translate += t;
	}

	void addRotation(const EMatrix3& r)
	{
		rotate *= r;
	}

	void reorient(ECoords& coords) const
	{
		coords.rowwise() += translate.transpose();
		coords *= rotate;
	}

	//reorient rd points
	void reorient(std::vector<RDGeom::Point3D>& coords) const
	{
		for (unsigned i = 0, n = coords.size(); i < n; i++)
		{
			EVector3 pt = EVector3(coords[i].x, coords[i].y, coords[i].z);
			pt = (pt + translate).transpose() * rotate;
			coords[i].x = pt(0);
			coords[i].y = pt(1);
			coords[i].z = pt(2);
		}
	}

	//apply inverse transformation
	void unorient(unsigned n, float *coords) const
	{
		using namespace Eigen;
		std::vector<EVector3> pnts(n);
		for (unsigned i = 0; i < n; i++)
		{
			pnts[i] = EVector3(coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]);
		}

		unorient(pnts);

		for (unsigned i = 0; i < n; i++)
		{
			coords[3 * i] = pnts[i].coeff(0);
			coords[3 * i + 1] = pnts[i].coeff(1);
			coords[3 * i + 2] = pnts[i].coeff(2);
		}
	}

	void unorient(ECoords& coords) const
	{
		coords = coords*rotate.inverse();
		coords.rowwise() -= translate.transpose();
	}

	void unorient(std::vector<RDGeom::Point3D>& coords) const
	{
		EMatrix3 inv = rotate.inverse();
		for (unsigned i = 0, n = coords.size(); i < n; i++)
		{
			EVector3 pt = EVector3(coords[i].x, coords[i].y, coords[i].z);
			pt = (pt.transpose()*inv).transpose() - translate;
			coords[i].x = pt(0);
			coords[i].y = pt(1);
			coords[i].z = pt(2);
		}
	}

	void unorient(std::vector<EVector3>& pnts) const
	{
		EMatrix3 inv = rotate.inverse();
		for (unsigned i = 0, n = pnts.size(); i < n; i++)
		{
			pnts[i] = (pnts[i].transpose()*inv).transpose() - translate;
		}
	}

	// for debugging
	void dump(std::ostream& out) const
	{
		out << translate << "\n";
		out << rotate << "\n";
	}
};

#endif /* ORIENTER_H_ */
