#ifndef INITIALISEMODEL_H_
#define INITIALISEMODEL_H_

#include "Model.h"

#include <vector>

using namespace std;


class InitialiseModel
{
public:
	InitialiseModel();
	virtual ~InitialiseModel();

	void CreateBone(const Coordinate &topLeftCorner, const Coordinate &bottomRightCorner,
	const LatticeType lattice, Model *model);
	
	void SetTumourComposition(const int patCellPopulation, const int pantCellPopulation,
		const int qatCellPopulation, const int qantCellPopulation, const int dCellPopulation, 
		const Coordinate &tumourCentre,	const int radius, const LatticeType lattice, Model *model);

	//Added by Marjan 1/3/12
	void SetTTheumourComposition(CellPopulation *cellPop);
	//Added by Marjan 1/3/12

	void SetTumourComposition(const char *patFileName, const char *pantFileName,
		const char *qatFileName, const char *qantFileName, const char *dFileName,
		const LatticeType lattice, Model *model);
	void SetRectangularTumourComposition(const int patCellPopulation, const int pantCellPopulation, 
		const int qatCellPopulation, const int qantCellPopulation, const int dCellPopulation, 
		const int tumourHeight, const int tumourWidth, const LatticeType lattice, Model *model);		
		
	void FindTumourGrowthRegion(const Coordinate &regionCentre, const int radius, Model *model);
	//Added by Marjan 8/3/12
	void FindTumourNonGrowthRegion(Model *model);
	//Added by Marjan 8/3/12
		
private:
	void FindDiskCoordinates(vector<Coordinate> *diskCoordinates, const Coordinate &tumourCentre, const int radius);
	void FindDiskCoordinates(vector<Coordinate> *diskCoordinates, 
		vector<Coordinate> *circleCoordinates, const Coordinate &tumourCentre, const int radius);
		
	void FindInnerCircularCoordinates(vector<Coordinate> *diskCoordinates, 
		vector<Coordinate> *circleCoordinates);
		 
	void FindOuterCircularCoordinates(vector<Coordinate> *coordinates, const Coordinate &centre, const int radius);
	
	void FindCircleCoordinates(vector<Coordinate> *coordinates, const Coordinate &centre, const int radius);
	void FindEdgeCoordinates(vector<Coordinate> *coordinates);
};

#endif /*INITIALISEMODEL_H_*/
