/*
 * CellPopulation.cpp
 *
 *  Created on: Oct 4, 2011
 *      Author: marjan
 */

#include "CellPopulation.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <map>


static long int nextKey;
static long int numberOfAliveCells;
static long int numberOfDeadCells;
static long int totalNumberOfCells = 0;
static int reverseSimulationTime;

using namespace std;


CellPopulation::CellPopulation(int initialNumberOfCells, int timeBetweenDivision, int lifeSpan)
{
	m_initialNumberOfCells = initialNumberOfCells;
	m_timeBetweenDivision = timeBetweenDivision;
	m_lifeSpan = lifeSpan;

	// Initialize the map container with as many cells as is decided in the main function
	Cell bCell[getInitialNumberOfCells()];
	totalNumberOfCells = getInitialNumberOfCells();

	for(int i = 0; i< getInitialNumberOfCells(); i++)
	{
		cells[i] = bCell[i];
		cells[i].setCellKey(i);
	}

	//generateRandomAge(&cells);

	nextKey = getInitialNumberOfCells() - 1;
	numberOfAliveCells = getInitialNumberOfCells();
	numberOfDeadCells = 0;
}

CellPopulation::~CellPopulation()
{
}

int CellPopulation::getInitialNumberOfCells(void)
{
	return m_initialNumberOfCells;
}

int CellPopulation::getLifeSpan(void)
{
	return m_lifeSpan;
}

int CellPopulation::getTimeBetweenDivision(void)
{
	return m_timeBetweenDivision;
}

void CellPopulation::generateRandomAge(map<long int, Cell> *cells)
{
	map<long int, Cell>::iterator jj;
	//Create random ages for the cells
	srand((unsigned)time(0));
	for(jj = cells->begin(); jj != cells->end(); jj++)
	{
		int random_integer = rand() % (getLifeSpan()+1);
		(jj->second).setRealAge(random_integer);
	}
}

void CellPopulation::evolve(CellPopulation *celPop, double copyLattice[121][121], double chemicalSource, double cellEvolutionParameter)
{
	//reverseSimulationTime = numberOfDays - 1;

	//while(numberOfDays-- > 0)
	//{
		//step(numberOfDays);
	//}
	step(celPop, copyLattice, chemicalSource, cellEvolutionParameter);
}

long int CellPopulation::generateNewCellKey(void)
{
	return ++nextKey;
}

void CellPopulation::step(CellPopulation *celPop, double copyLattice[121][121], double chemicalSource, double cellEvolutionParamter)
{
	map<long int, Cell>::iterator kk;
	map<long int, Cell>::iterator ll;


	//All the cells can contribute in the for loop within this time step
	for(ll = (celPop->cells).begin(); ll != (celPop->cells).end(); ll++)
	{
		(ll->second).setJustBorn(true);
		(ll->second).divideImmedialtely = false;
	}

	for(kk = (celPop->cells).begin(); kk != (celPop->cells).end(); kk++)
	{
		bool flagPregnancy = true;
		bool flagDeath = true;

		// This process accelerates the cell division
		int y = (kk->second).cellPosition[0];
		int x = (kk->second).cellPosition[1];
		double chem = copyLattice[y][x];

		if((chem >= cellEvolutionParamter * chemicalSource))
		{
			double condition = (kk->second).getJustBorn();
			if(condition == true)
			{
				(kk->second).divideImmedialtely = true;
			}
		}

		//cout<< (int)celPop->cellPop->size()<< endl;
		//cout<< (kk->second).m_cellKey<< endl;
		// if the it is not the time for the cell to divide yet
		if((kk->second).getPregnancyAge() < getTimeBetweenDivision() && (kk->second).getRealAge() < getLifeSpan() &&
				(kk->second).getLifeStatusOfCell() && (kk->second).getJustBorn())
		{
			//Increase the pregnancy age
			(kk->second). updatePregnancyAge();
			//cout<< "Cell key is: "<< (kk->second).m_cellKey <<" The pregnancy age is: " << (kk->second).getPregnancyAge() <<endl;
			//Don't execute the next if-statements
		    flagPregnancy = false;
		}

		// if it is the time for the cell to divide
		if(((kk->second).getPregnancyAge() == getTimeBetweenDivision() &&
				(kk->second).getRealAge() < getLifeSpan() && (kk->second).getLifeStatusOfCell() &&
				flagPregnancy && (kk->second).getJustBorn()) || ((kk->second).divideImmedialtely))
		{
			(kk->second).divideImmedialtely = false;
			long int key = generateNewCellKey();
			cells.insert(pair<long int,Cell>(key, (kk->second).divide(key,
					     kk->second.cellPosition[0], (kk->second).cellPosition[1])));
			//cout<< "The map size is: "<< (int)cells.size()<< endl;
			numberOfAliveCells++;
			totalNumberOfCells++;
			//cout<< "Alive cells are: " <<numberOfAliveCells <<endl;
			(kk->second).resetPregnancyAgeAfterGivingBirth();
		}

		// update the realAge before death
		if((kk->second).getRealAge() < getLifeSpan() && (kk->second).getJustBorn())
		{
			(kk->second).updateRealAge();
			//cout << "Cell key is: "<< (kk->second).m_cellKey<< " Real age is: " << (kk->second).getRealAge() <<endl;
			flagDeath = false;
		}

		// remove the deadCells from the map container
		if(((kk->second).getRealAge() == getLifeSpan()) && flagDeath && (kk->second).getJustBorn() && (kk->second).getLifeStatusOfCell())
		{
			cells.erase(kk);
			numberOfDeadCells++;
			numberOfAliveCells--;
			totalNumberOfCells--;
			//cout<< "Dead cells are: "<<numberOfDeadCells <<endl;
			(kk->second).setLifeStatusOfCell(false);
		}
	}//end if for loop

	//Log the number of alive and dead cells after iterating in the for loop
	//logFile(numberOfDays, numberOfAliveCells, numberOfDeadCells);

	// Commented out by Marjan for test 25/4/12
	//saveInFile(numberOfDays, totalNumberOfCells);
}

/*
void CellPopulation::logFile(int numberOfDays, long int numberOfAliveCells, long int numberOfDeadCells)
{
		ofstream cellPopulation("/Users/marjan/Desktop/cells03.txt", ios::app);

		if(!cellPopulation)
		{
			exit(0);
		}
		else
		{
			if(writeOnce)
			{
        		cellPopulation << "Starting with "<< getInitialNumberOfCells() <<" cells " <<" T(Division) = " << getTimeBetweenDivision() << " T(Death) = "<< getLifeSpan()<< endl;
        		writeOnce = false;
			}
			//cellPopulation << "Simulation time: " << abs(numberOfDays - reverseSimulationTime)<< " AliveCells: "<< numberOfAliveCells << " Dead cells: " << numberOfDeadCells << endl;
			cellPopulation << "Simulation time: " << abs(numberOfDays - reverseSimulationTime)<< " totalNumberOfCells: "<< totalNumberOfCells << endl;
		}
}
*/

void CellPopulation::saveInFile(int simulationTime, long int totalNumeberOfCells)
{
	ofstream cellPopulation("/Users/marjan/Desktop/cells08.txt", ios::app);

	if(!cellPopulation)
	{
		exit(0);
	}
	else
	{
		//cellPopulation << "Simulation time: " << abs(numberOfDays - reverseSimulationTime)<< " AliveCells: "<< numberOfAliveCells << " Dead cells: " << numberOfDeadCells << endl;
		float result = (float(abs(simulationTime - reverseSimulationTime))/(float)getLifeSpan());
		cellPopulation << result<< "\t\t"<<totalNumberOfCells << endl;
	}
}
