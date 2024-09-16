#ifndef DAMPINGFORCE_H
#define DAMPINGFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class dampingForce
{
public:
	dampingForce(elasticPlate &m_plate, timeStepper &m_stepper, double m_viscosity);
	~dampingForce();
	void computeFd();
	void computeJd();

	void setFirstJacobian();

	VectorXd TotalForceVec;

private:
	elasticPlate *plate;
    timeStepper *stepper;

    double viscosity;

    Vector3d u1;
    Vector3d u2;
    Vector3d f1;
    Vector3d f2;

    Vector3d jac;

    double ind;

    int index1;
    int index2;

    double edgeLength;

    double sectionArea;

    double mass1, mass2;

};

#endif
