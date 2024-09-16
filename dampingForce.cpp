#include "dampingForce.h"

dampingForce::dampingForce(elasticPlate &m_plate, timeStepper &m_stepper, double m_viscosity)
{
	plate = &m_plate;
    stepper = &m_stepper;

	viscosity = m_viscosity;

	TotalForceVec = VectorXd::Zero(plate->ndof);

	// sectionArea = 2 * plate->radius;
	sectionArea = 1.0;
}

dampingForce::~dampingForce()
{
	;
}

void dampingForce::computeFd()
{

	/*
	TotalForceVec = VectorXd::Zero(plate->ndof);

	for (int i = 0; i < plate->ndof; i++)
	{
		stepper->addForce(i, viscosity * plate->massArray[i] * plate->u[i]);
	}
	*/

	for (int i = 0; i < plate->edgeNum; i++)
	{
		index1 = plate->v_edgeElement[i].nv_1;
		index2 = plate->v_edgeElement[i].nv_2;

		u1 = plate->getVelocity(index1);
		u2 = plate->getVelocity(index2);

		mass1 = plate->massArray(3 * index1 + 0);
		mass2 = plate->massArray(3 * index2 + 0);

		f1 = - 0.5 * u1 * mass1 * viscosity;
		for (int k = 0; k < 3; k++)
		{
			ind = 3 * index1 + k;
			
			stepper->addForce(ind, - f1[k]);
		}

		f2 = - 0.5 * u2 * mass1 * viscosity;
		for (int k = 0; k < 3; k++)
		{
			ind = 3 * index2 + k;
			
			stepper->addForce(ind, - f2[k]);
		}
	}
}

void dampingForce::computeJd()
{
	
	for (int i = 0; i < plate->nv; i++)
	{
		mass1 = plate->massArray(3 * index1 + 0);
		mass2 = plate->massArray(3 * index2 + 0);

		index1 = plate->v_edgeElement[i].nv_1;
		index2 = plate->v_edgeElement[i].nv_2;

		jac(0) = - 0.5 * mass1 * viscosity / plate->dt;
		jac(1) = - 0.5 * mass1 * viscosity / plate->dt;
		jac(2) = - 0.5 * mass1 * viscosity / plate->dt;
		
		for (int j = 0; j < 3; j++)
		{
			ind = 3 * index1 + j;
			stepper->addJacobian(ind, ind, - jac(j));
		}

		jac(0) = - 0.5 * mass2 * viscosity / plate->dt;
		jac(1) = - 0.5 * mass2 * viscosity / plate->dt;
		jac(2) = - 0.5 * mass2 * viscosity / plate->dt;

		for (int j = 0; j < 3; j++)
		{
			ind = 3 * index2 + j;
			stepper->addJacobian(ind, ind, - jac(j));
		}
	}
	

	/*

	for (int i = 0; i < plate->ndof; i++)
	{
		stepper->addJacobian(i, i, plate->massArray[i] * viscosity / plate->dt);
	}

	*/
}

void dampingForce::setFirstJacobian()
{
	for (int i = 0; i < plate->nv; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			ind = 3 * i + j;
			stepper->addJacobian(ind, ind, 1);
		}
	}
}
