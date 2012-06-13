#include "InitialiseModel.h"

#include <cmath>
#include <math.h>
#include <algorithm>
#include <fstream>


/*
 * Default constructor
 */
InitialiseModel::InitialiseModel() {
}


/*
 * Clean up the initialised Model
 */
InitialiseModel::~InitialiseModel() {
}

/*
 * Create a slice of bone between the topLeftCorner coordinate and bottomRightCorner coordinate
 */
void InitialiseModel::CreateBone(const Coordinate &topLeftCorner, const Coordinate &bottomRightCorner, 
	const LatticeType lattice, Model *model) {

	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1);

	// [cells/site]
	int boneDensity = 0;
	switch (model->m_porosityEqu) {
		case HYP_POROSITY_106701:
			boneDensity = 106701;
			break;
		case HYP_POROSITY_1224:
			boneDensity = 1224;
			break;
		case HYP_POROSITY_691:
			boneDensity = 691;
			break;
		case HYP_POROSITY_300:
			boneDensity = 300;
			break;
		default:
			throw string("Error: porosity equation not recognised.");
	}

	for (int x = topLeftCorner.x; x <= bottomRightCorner.x; ++x) {
		for (int y = topLeftCorner.y; y <= bottomRightCorner.y; ++y) {
			// Set the bone's cellular composition
			(*l)[y][x]->SetCells(0, 0, 0, 0, 0, boneDensity);
			// Set the bone's chemical state
			(*l)[y][x]->SetGlucoseConcentration(0);
			(*l)[y][x]->SetOxygenConcentration(0);
		}
	}

	/* Create a ring of qat cells around the wall
	int offset = 2;
	for (int x = topLeftCorner.x-offset; x <= bottomRightCorner.x+offset; ++x) {
		int y = topLeftCorner.y-offset;
		(*l)[y][x]->SetCells(0, 0, 10, 0, 0);
	}	
	for (int x = topLeftCorner.x-offset; x <= bottomRightCorner.x+offset; ++x) {
		int y = bottomRightCorner.y+offset;
		(*l)[y][x]->SetCells(0, 0, 10, 0, 0);
	}
	for (int y = topLeftCorner.y-offset; y <= bottomRightCorner.y+offset; ++y) {
		int x = topLeftCorner.x-offset;
		(*l)[y][x]->SetCells(0, 0, 10, 0, 0);
	}
	for (int y = topLeftCorner.y-offset; y <= bottomRightCorner.y+offset; ++y) {
		int x = bottomRightCorner.x+offset;
		(*l)[y][x]->SetCells(0, 0, 10, 0, 0);
	} */
}

/*
 * Initialise the model from a file. The file must be in the same formats as those saved by the model: csv etc
 * 
 * ...FileName: the file name of the file containing the respective tumour cell metabolic population data.
 * lattice: the lattice to populate, either T or T1.
 * model: the model object.
 */
void InitialiseModel::SetTumourComposition(const char *patFileName, const char *pantFileName,
	const char *qatFileName, const char *qantFileName, const char *dFileName, 
	const LatticeType lattice, Model *model) {
		
	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1);		
		
	// Open the cellular description files		
	ifstream fPatIn(patFileName, ios_base::in); 
	// Check for errors
	if (fPatIn.fail()) throw string("Error: failed to open "+string(patFileName)+".");
	ifstream fPantIn(pantFileName, ios_base::in);
	if (fPantIn.fail()) throw string("Error: failed to open "+string(pantFileName)+".");
	ifstream fQatIn(qatFileName, ios_base::in);
	if (fQatIn.fail()) throw string("Error: failed to open "+string(qatFileName)+".");
	ifstream fQantIn(qantFileName, ios_base::in);
	if (fQantIn.fail()) throw string("Error: failed to open "+string(qantFileName)+".");
	ifstream fDIn(dFileName, ios_base::in);
	if (fDIn.fail()) throw string("Error: failed to open "+string(dFileName)+".");
	
	for (int i = 0; i < HEIGHT; ++i) {
		for (int j = 0; j < WIDTH; ++j) {					
			
			int patCellPopulation, pantCellPopulation, qatCellPopulation, qantCellPopulation, dCellPopulation;
			patCellPopulation = pantCellPopulation = qatCellPopulation = qantCellPopulation = dCellPopulation = 0;
			// Read the cell populations from file
			fPatIn >> patCellPopulation;
			fPantIn >> pantCellPopulation;
			fQatIn >> qatCellPopulation;
			fQantIn >> qantCellPopulation;
			fDIn >> dCellPopulation;
			
			// Set the lattice site cellular composition		
			(*l)[i][j]->SetCells(patCellPopulation, pantCellPopulation, qatCellPopulation, qantCellPopulation, dCellPopulation, 0);
			
			// Lines end with a non-comma character
			if (j < WIDTH-1) {
				for (int i = 0 ; i < 5; ++i) {	
					string comma("");
					switch (i) {
						case 0:
							fPatIn >> comma;
							if (comma != ",") throw string("Error: Contrary to expectations non ',' character read, "+comma+" read instead.");
							break;
						case 1:
							fPantIn >> comma;
							if (comma != ",") throw string("Error: Contrary to expectations non ',' charactor read, "+comma+" read instead.");
							break;
						case 2:
							fQatIn >> comma;
							if (comma != ",") throw string("Error: Contrary to expectation non ',' character read, "+comma+" read instead.");
							break;
						case 3:
							fQantIn >> comma;
							if (comma != ",") throw string("Error: Contrary to expectation non ',' character read, "+comma+" read instead.");
							break;
						case 4:
							fDIn >> comma;
							if (comma != ",") throw string("Error: Contrary to expectation non ',' character read, "+comma+" read instead.");
							break;
					}
				}
			}									
		}
	}	
	fPatIn.close();
	fPantIn.close();
	fQatIn.close();
	fQantIn.close();
	fDIn.close(); 
}

/*
 * Initialise a heterogeneous tumour with the given composition, location and radius.
 * 
 * ...CellPopulation: the respective metabolic population of a lattice site.
 * tumourCentre: the tumour's centre.
 * radius: the tumour's radius. 
 * lattice: the lattice to populate, either T or T1.
 * model: the model object.
 */
void InitialiseModel::SetTumourComposition(const int patCellPopulation, 
	const int pantCellPopulation, const int qatCellPopulation, const int qantCellPopulation, 
	const int dCellPopulation, 
	const Coordinate &tumourCentre, const int radius, const LatticeType lattice, Model *model) {
		
	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1);
		
	vector<Coordinate> *tumourSites = NULL;
	try {
		tumourSites = new vector<Coordinate>(); 
		// Find the lattice sites occupied by the tumour	
		FindDiskCoordinates(tumourSites, tumourCentre, radius);
		
		// Initiise the tumour lattice sites
		vector<Coordinate>::const_iterator it;
		for (it = tumourSites->begin(); it != tumourSites->end(); ++it) {
			int y = (*it).y;
			int x = (*it).x;
						
			// Set the lattice site cellular composition		
			(*l)[y][x]->SetCells(patCellPopulation, pantCellPopulation, qatCellPopulation, qantCellPopulation, dCellPopulation, 0);
		}
		delete tumourSites;
	}
	catch (string &errStr) {
		// Clean up
		delete tumourSites;
		throw errStr;	
	}
	catch (bad_alloc &e) {
		// Clean up
		delete tumourSites;
		throw e;
	}
	// Catch all errors
	catch (...) {
		// Clean up
		delete tumourSites;
	}		
}

//Added by Marjan 1/3/12
/*
 * Randomly assign positions to the cells
 */
void InitialiseModel::SetTTheumourComposition(CellPopulation *cell)
{
	// Assume that the lattice size is 121*121 and assign x and y to it
	srand((unsigned)time(0));
	map<long int, Cell>::iterator kk;

	for(kk = (cell->cells).begin(); kk != (cell->cells).end(); kk++)
	{
		double doubleYKeeper = ((double)rand() / (RAND_MAX)) * (99-25+1) + 1;//81 and 42
		double doubleXKeeper = ((double)rand() / (RAND_MAX)) * (106-16+1) + 16;
		(kk->second).cellPosition[0] = doubleYKeeper;
		(kk->second).cellPosition[1] = doubleXKeeper;
	}
	for(kk = (cell->cells).begin(); kk != (cell->cells).end(); kk++)
	{
		cout<< (kk->second).cellPosition[0]<< "\t"<< (kk->second).cellPosition[1]<< endl;
	}
}
//Added by Marjan 1/3/12

/*
 * Create a tumour in the centre of the lattice, located in a rectangle of lattice sites, with the given
 * composition.
 * 
 * ...CellPopulation: the respective metabolic population of a lattice site. 
 * tumourHeight: the height of the tumour.
 * tumourWidth: the width of the tumour. 
 * lattice: the lattice to populate, either T or T1.
 * model: the model object.
 */
void InitialiseModel::SetRectangularTumourComposition(const int patCellPopulation, 
	const int pantCellPopulation, const int qatCellPopulation, const int qantCellPopulation, 
	const int dCellPopulation, 
	const int tumourHeight, const int tumourWidth, const LatticeType lattice, Model *model) {
		
	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1);		
		
	int latticeCentreX = WIDTH / 2;
	int latticeCentreY = HEIGHT / 2;	
	for (int x = latticeCentreX - tumourWidth/2; x <= latticeCentreX + tumourWidth/2; ++x) {
		for (int y = latticeCentreY - tumourHeight/2; y <= latticeCentreY + tumourHeight/2; ++y) {						
			// Set the lattice site cellular composition		
			(*l)[y][x]->SetCells(patCellPopulation, pantCellPopulation, qatCellPopulation, qantCellPopulation, dCellPopulation, 0);
		}
	}
}

/*
 * Find the coordinates of lattice sites on the lattice perimeter
 */
void InitialiseModel::FindEdgeCoordinates(vector<Coordinate> *coordinates) {
	for (int i = 0; i < HEIGHT; ++i) {
		// The top edge
		Coordinate coordTop;
		coordTop.x = i;
		coordTop.y = 0;
		// The left edge
		Coordinate coordLeft;
		coordLeft.x = 0;
		coordLeft.y = i;
		// The right edge
		Coordinate coordRight;
		coordRight.x = HEIGHT - 1;
		coordRight.y = i;
		// The bottom edge
		Coordinate coordBottom;
		coordBottom.x = i;
		coordBottom.y = HEIGHT - 1;
		coordinates->push_back(coordTop);
		coordinates->push_back(coordLeft);
		coordinates->push_back(coordRight);
		coordinates->push_back(coordBottom);
	}	
}

/*
 * Find the lattice site that form a circle with given radius, which is centred on the centre of
 * the lattice.
 * 
 * Note: the circle occupies the last site of the radius.
 */
void InitialiseModel::FindCircleCoordinates(vector<Coordinate> *coordinates, const Coordinate &centre, const int radius) {
	// The circle is assumed to be centred on the lattice centre
	int radiusSquared = pow(radius, 2.0);
	
	// Find the circle coordinates horizontally with the given radius
	for (int i = centre.x - radius; i <= centre.x + radius; ++i) {
		// Check the tumour is going to be created within the lattice
		if (i < 0 || i >= WIDTH)
			throw string("Error: tumour is initialised outside of the lattice.");
		Coordinate coordBottom;
		Coordinate coordTop;
		coordBottom.x = i;
		coordTop.x = i;
		coordBottom.y = sqrt(radiusSquared - pow(i - centre.x, 2.0)) + centre.y;
		// Find the y coordinate reflected in the middle row
		coordTop.y = HEIGHT - 1 - coordBottom.y;
		// An alternatice circle
		//coordTop->y = -1 * sqrt(radiusSquared - pow(i - latticeCentreX, 2)) + latticeCentreY;
		//coordBottom->y = HEIGHT - 1 - coordTop->y;
		
		// If the coordinates have not already been added circle coordinates, add them
		if (find(coordinates->begin(), coordinates->end(), coordBottom) == coordinates->end())			
			coordinates->push_back(coordBottom);		
		if (find(coordinates->begin(), coordinates->end(), coordTop) == coordinates->end())			
			coordinates->push_back(coordTop);	
	}
	// Find the circle coordinates vertically with the given radius
	for (int i = centre.y - radius; i <= centre.y + radius; ++i) {
		// Check the tumour is going to be created within the lattice
		if (i < 0 || i >= HEIGHT)
			throw string("Error: tumour is initialised outside of the lattice.");
		Coordinate coordRight;
		Coordinate coordLeft;
		coordRight.y = i;
		coordLeft.y = i;
		coordRight.x = sqrt(radiusSquared - pow(i - centre.y, 2.0)) + centre.x;
		// Find the x coordinate reflected in the middle column
		coordLeft.x = WIDTH - 1 - coordRight.x;
		// An alternative circle
		//coordLeft->x = -1 * sqrt(radiusSquared - pow(i - latticeCentreY, 2)) + latticeCentreX;
		//coordRight->x = WIDTH - 1 - coordLeft->x;
		
		// If the coordinates have not already been added to the circle coordinats, add them
		if (find(coordinates->begin(), coordinates->end(), coordRight) == coordinates->end())			
			coordinates->push_back(coordRight);		
		if (find(coordinates->begin(), coordinates->end(), coordLeft) == coordinates->end())			
			coordinates->push_back(coordLeft);	
	}
} 

/*
 * Find the lattice sites between a circle and the lattice edge
 */
void InitialiseModel::FindOuterCircularCoordinates(vector<Coordinate> *coordinates, const Coordinate &centre, const int radius) {
	vector<Coordinate> *circleCoordinates = NULL;
	try {
		circleCoordinates = new vector<Coordinate>();
		FindCircleCoordinates(circleCoordinates, centre, radius);
		
		vector<Coordinate>::const_iterator it;
		for (it = circleCoordinates->begin(); it != circleCoordinates->end(); ++it) {
			int y = (*it).y;
			int x = (*it).x;
			// Fill in the lattice sites at the top of the lattice
			if (y < HEIGHT / 2) {
				// Ignore duplicated coordinates
				if (find(coordinates->begin(), coordinates->end(), *it) == coordinates->end()) {
					--y;
					while (y > 0) {
						Coordinate coord;
						coord.y = y;
						coord.x = x;
						coordinates->push_back(coord);
						--y;
					}
					coordinates->push_back(*it);
				}
			}
			// Fill in the lattice sites at the bottom of the lattice
			else if (y > HEIGHT / 2) {
				// Ignore duplicates
				if (find(coordinates->begin(), coordinates->end(), *it) == coordinates->end()) {
					++y;
					while (y < HEIGHT-1) {
						Coordinate coord;
						coord.y = y;
						coord.x = x;
						coordinates->push_back(coord);
						++y;
					}
					coordinates->push_back(*it);
				}
			}
			// Fill in the lattice sites above and below
			else {
				// Ignore duplicates
				if (find(coordinates->begin(), coordinates->end(), *it) == coordinates->end()) {
					--y;
					while (y > 0) {
						Coordinate coord;
						coord.y = y;
						coord.x = x;
						coordinates->push_back(coord);
						--y;
					}
					y = (*it).y;
					++y; 
					while (y < HEIGHT-1) {
						Coordinate coord;
						coord.y = y;
						coord.x = x;
						coordinates->push_back(coord);
						++y;
					}
					coordinates->push_back(*it);
				}
			}
		}
		// Clean up 
		delete circleCoordinates;
	}
 	catch (string &errStr) {
		// Clean up
		delete circleCoordinates;
		throw errStr;	
	}
	catch (bad_alloc &e) {
		// Clean up
		delete circleCoordinates;
		throw e;
	}
	// Catch all errors
	catch (...) {
		// Clean up
		delete circleCoordinates;
	}	
}

/*
 * Find the lattice sites within a circle and the circle coordinates
 */
void InitialiseModel::FindDiskCoordinates(vector<Coordinate> *diskCoordinates, 
	vector<Coordinate> *circleCoordinates, const Coordinate &tumourCentre, const int radius) {
		
	FindCircleCoordinates(circleCoordinates, tumourCentre, radius);
	
	FindInnerCircularCoordinates(diskCoordinates, circleCoordinates);
}

/*
 * Find the lattice sites within a circle
 */
void InitialiseModel::FindDiskCoordinates(vector<Coordinate> *diskCoordinates, const Coordinate &tumourCentre, const int radius) {
	vector<Coordinate> *circleCoordinates;
	try {
		circleCoordinates = new vector<Coordinate>();
		FindCircleCoordinates(circleCoordinates, tumourCentre, radius);
		
		FindInnerCircularCoordinates(diskCoordinates, circleCoordinates);
		delete circleCoordinates;
	}
 	catch (string &errStr) {
		// Clean up
		delete circleCoordinates;
		throw errStr;	
	}
	catch (bad_alloc &e) {
		// Clean up
		delete circleCoordinates;
		throw e;
	}
	// Catch all errors
	catch (...) {
		// Clean up
		delete circleCoordinates;
	}	
}
	
/*
 * Find the coordinates within a circle
 */	
void InitialiseModel::FindInnerCircularCoordinates(vector<Coordinate> *diskCoordinates, 
	vector<Coordinate> *circleCoordinates) {
	
	vector<Coordinate>::const_iterator it; 
	for (it = circleCoordinates->begin(); it != circleCoordinates->end(); ++it) {
		int y = (*it).y;
		int x = (*it).x; 	
		// Fill in the lattice sites in the top half of the circle
		if (y < HEIGHT / 2) {
			// Ignore duplicated coordinates
			if (find(diskCoordinates->begin(), diskCoordinates->end(), *it) == diskCoordinates->end()) {
				++y;
				while (y <= HEIGHT / 2) {
					Coordinate coord;
					coord.x = x;
					coord.y = y; 					
					diskCoordinates->push_back(coord);
					++y;
				} 				
				diskCoordinates->push_back(*it);
			}	
		}
		// Fill in the lattice site in the bottom half of the circle
		else if (y > HEIGHT / 2) {
			// Ingnore duplicated coordinates
			if (find(diskCoordinates->begin(), diskCoordinates->end(), *it) == diskCoordinates->end()) {
				--y;
				while (y > HEIGHT / 2) {
					Coordinate coord;
					coord.x = x;
					coord.y = y; 					
					diskCoordinates->push_back(coord);	
					--y;
				} 				
				diskCoordinates->push_back(*it);				
			}
		}
		else {			
			if (find(diskCoordinates->begin(), diskCoordinates->end(), *it) == diskCoordinates->end()) 				
				diskCoordinates->push_back(*it);
		}
	}
}

/*
 * Find a circular tumour growth region of the given radius, centred on the lattice centre
 * 
 * As chemical is only diffused within the growth region, including the source sites will
 * have the effect of maintaining a constant chemical concentration within these sites.* 
 */
 void InitialiseModel::FindTumourGrowthRegion(const Coordinate &regionCentre, const int radius, Model *model) {
 	
	// Find a filled circular region for the tumour to grow within
	FindDiskCoordinates(model->m_tumourGrowthSites, model->m_tumourGrowthEdgeSites, regionCentre, radius);	 	
 }

 //Added by Marjan 8/3/12
 /*
  * Find the tumour non growth region
  */
 void InitialiseModel::FindTumourNonGrowthRegion(Model *model)
 {
		Coordinate coord;
		vector<Coordinate>::const_iterator it,ij;

		bool flag = true;
		for(int yPos = 0; yPos < HEIGHT; yPos++){
			for(int xPos = 0; xPos < WIDTH; xPos++){
				for(it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); it++){
					if((xPos == it->x) && (yPos == it->y)){
						flag = false;
					}
				}
				if(flag){
					coord.x = xPos;
					coord.y = yPos;
					model->m_nonTumourGrowthRegion->push_back(coord);
				}
				flag = true;
			}
		}

		/*for(ij = model->m_tumourGrowthEdgeSites->begin(); ij != model->m_tumourGrowthEdgeSites->end(); ij++){
			cout << ij->x << '\t'<< ij->y<< endl;
		}*/
 }
