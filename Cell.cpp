/*
 * Cell.cpp
 *
 *  Created on: Oct 4, 2011
 *      Author: marjan
 */

#include "Cell.h"
#include <iostream>

using namespace std;

Cell::Cell()
{
	pregnancyAge = 0;
	alive = true;
	justBorn = true;
	realAge = 0;
	divideImmedialtely = false;
	pressureInduced = false;
	chosenForPressure = false;
}

Cell::~Cell()
{
}

// Give birth to a daughter cell
Cell Cell::divide(long int cellKey, double motherY, double motherX)
{
	Cell daughter;

	daughter.pregnancyAge = 1;
	daughter.realAge = 1;
	daughter.alive = true;
	daughter.justBorn = false;
	daughter.m_cellKey = cellKey;
	daughter.cellPosition[0] = motherY;
	daughter.cellPosition[1] = motherX;
	daughter.divideImmedialtely = false;
	daughter.pressureInduced = false;
	daughter.chosenForPressure = false;
	//cout<< "The new cellKey is: "<< daughter.m_cellKey<< endl;
	return daughter;
}

void Cell::setCellKey(long int key)
{
	m_cellKey = key;
}

long int Cell::getCellKey(void)
{
	return m_cellKey;
}

void Cell::updateRealAge(void)
{
	realAge++;
}

int Cell::getRealAge(void)
{
	return realAge;
}

void Cell::setRealAge(int rlAge)
{
	realAge = rlAge;
}

int Cell::getLifeSpan(void)
{
	return lifeSpan;
}

void Cell::setLifeSpan(int lfpSpan)
{
	lifeSpan = lfpSpan;
}

bool Cell::getLifeStatusOfCell(void)
{
	return alive;
}

void Cell::setLifeStatusOfCell(bool lfStatus)
{
	alive = lfStatus;
}

bool Cell::getJustBorn(void)
{
	return justBorn;
}

void Cell::setJustBorn(bool jstBorn)
{
	justBorn = jstBorn;
}

void Cell::updatePregnancyAge(void)
{
	pregnancyAge++;
}
int Cell::getPregnancyAge(void)
{
	return pregnancyAge;
}

void Cell::resetPregnancyAgeAfterGivingBirth(void)
{
	pregnancyAge = 1;
}
