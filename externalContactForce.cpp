#include "externalContactForce.h"

externalContactForce::externalContactForce(elasticPlate &m_plate, timeStepper &m_stepper, 
	double m_dBar, double m_stiffness, double m_mu, double m_epsilonV)
{
	plate = &m_plate;
    stepper = &m_stepper;

    dBar = m_dBar;

    stiffness = m_stiffness;

    mu = m_mu;

    epsilonV = m_epsilonV;

    f.setZero(3);

    dx.setZero(9);
    dxdx.setZero(9,9);
    jacobian.setZero(9,9);

    Id3<<1,0,0,
	     0,1,0,
	     0,0,1;
}

externalContactForce::~externalContactForce()
{
	;
}

void externalContactForce::computeFc()
{
	minDistance = 100000.00;
	contactIndex = 0;

	for (int i = 0; i < plate->contactNum; i++)
	{
		nv_0 = plate->v_ContactElement[i].nv_0;
		edgeLocal = plate->v_ContactElement[i].edgeIndex;

		nv_1 = plate->v_edgeElement[edgeLocal].nv_1;
		nv_2 = plate->v_edgeElement[edgeLocal].nv_2;

		x_0 = plate->getVertexOld(nv_0);
		x_1 = plate->getVertexOld(nv_1);
		x_2 = plate->getVertexOld(nv_2);

		e_1 = x_0 - x_1;
		e_2 = x_0 - x_2;
		e_3 = x_2 - x_1;

		e_4 = e_1.cross(e_2);

		d1 = ( x_0 - x_1 ).norm();
		d2 = ( x_0 - x_2 ).norm();
		d3 = e_4.norm() / e_3.norm();

		d = d1 + d2;

		if (d < minDistance)
		{
			minDistance = d;
			contactIndex = i;
		}
	}

	if (contactIndex != 0 && contactIndex != plate->nv - 2)
	{
		nv_0 = plate->v_ContactElement[contactIndex].nv_0;
		edgeLocal = plate->v_ContactElement[contactIndex].edgeIndex;

		nv_1 = plate->v_edgeElement[edgeLocal].nv_1;
		nv_2 = plate->v_edgeElement[edgeLocal].nv_2;

		x_0 = plate->getVertexOld(nv_0);
		x_1 = plate->getVertexOld(nv_1);
		x_2 = plate->getVertexOld(nv_2);

		e_1 = x_0 - x_1;
		e_2 = x_0 - x_2;
		e_3 = x_2 - x_1;

		e_4 = e_1.cross(e_2);

		normE1 = e_1.norm();
		normE2 = e_2.norm();
		normE3 = e_3.norm();
		normE4 = e_4.norm();

		t_1 = e_1 / normE1;
		t_2 = e_2 / normE2;
		t_3 = e_3 / normE3;
		t_4 = e_4 / normE4;

		d1 = ( x_0 - x_1 ).norm();
		d2 = ( x_0 - x_2 ).norm();
		d3 = e_4.norm() / e_3.norm();

		//	if (disType == 3)
		{
			arrayNum.setZero(9);

			arrayNum(0) = 3 * nv_0 + 0;
			arrayNum(1) = 3 * nv_0 + 1;
			arrayNum(2) = 3 * nv_0 + 2;
			arrayNum(3) = 3 * nv_1 + 0;
			arrayNum(4) = 3 * nv_1 + 1;
			arrayNum(5) = 3 * nv_1 + 2;
			arrayNum(6) = 3 * nv_2 + 0;
			arrayNum(7) = 3 * nv_2 + 1;
			arrayNum(8) = 3 * nv_2 + 2;

			fVec.setZero(9);

			d = e_4.norm() / e_3.norm();

			dEnergydD = - 2 * d * log( (-d+dBar) / dBar ) - d * d  / (d - dBar);

			de1 =   crossMat(e_2) * t_4 / normE3;
			de2 = - crossMat(e_1) * t_4 / normE3;
			de3 = - t_3 * normE4 / (normE3 * normE3);

			fVec.setZero(9);

			fVec.segment(0,3) =   de1 + de2;
			fVec.segment(3,3) = - de1 - de3;
			fVec.segment(6,3) = - de2 + de3;

			fVec = - fVec * stiffness * dEnergydD;

			for (int kk = 0; kk < 9; kk++)
			{
				ind = arrayNum(kk);
				stepper->addForce(ind, -fVec[kk]);
			}

		}

		// friction

		x_0 = plate->getVertexOld(nv_0);
		x_1 = plate->getVertexOld(nv_1);
		x_2 = plate->getVertexOld(nv_2);

		e_1 = x_0 - x_1;
		e_2 = x_0 - x_2;
		e_3 = x_2 - x_1;

		e_4 = e_1.cross(e_2);

		normE1 = e_1.norm();
		normE2 = e_2.norm();
		normE3 = e_3.norm();
		normE4 = e_4.norm();

		t_1 = e_1 / normE1;
		t_2 = e_2 / normE2;
		t_3 = e_3 / normE3;
		t_4 = e_4 / normE4;

		d1 = ( x_0 - x_1 ).norm();
		d2 = ( x_0 - x_2 ).norm();
		d3 = e_4.norm() / e_3.norm();

		v_0 = plate->getVelocityOld(nv_0);
		v_1 = plate->getVelocityOld(nv_1);
		v_2 = plate->getVelocityOld(nv_2);

		normDirection = t_3.cross(t_1.cross(t_3));

		v_r = (v_1 + v_2) / 2 - v_0;

		v_rt = v_r - v_r.dot(normDirection) * normDirection;

		if ( v_rt.norm() >= epsilonV )
		{
			fVelocity = 1.0;
				
			tK = v_rt / v_rt.norm();
		}
		else
		{
			fVelocity = - ( v_rt.norm() * v_rt.norm() ) / (epsilonV * epsilonV) + 2 * v_rt.norm() / epsilonV;

			tK = v_rt / (v_rt.norm() + 1e-15);
		}

		f = fVec.segment(0,3);

		friction = mu * f.norm() * fVelocity * tK;

		for (int k = 0; k < 3; k++)
		{
			ind = arrayNum(k);
			stepper->addForce(ind, -friction[k]);

			ind = arrayNum(k+3);
			stepper->addForce(ind, friction[k] / 2);

			ind = arrayNum(k+6);
			stepper->addForce(ind, friction[k] / 2);
		}

		/*

		// friction

		v_1 = plate->getVelocityOld(plate->v_ContactElement[contactIndex].nv_1);
		v_2 = plate->getVelocityOld(plate->v_ContactElement[contactIndex].nv_2);

		v_r = v_2 - v_1;

		v_rt = v_r - v_r.dot(dDdEdge) * dDdEdge;

		if ( v_rt.norm() >= epsilonV )
		{
			fVelocity = 1.0;
				
			tK = v_rt / v_rt.norm();
		}
		else
		{
			fVelocity = - ( v_rt.norm() * v_rt.norm() ) / (epsilonV * epsilonV) + 2 * v_rt.norm() / epsilonV;

			tK = v_rt / (v_rt.norm() + 1e-15);
		}

		friction = mu * f.norm() * fVelocity * tK;

		for (int k = 0; k < 3; k++)
		{
			ind = arrayNum(k);
			//stepper->addForce(ind, -friction[k]);

			ind = arrayNum(k+3);
			//stepper->addForce(ind, friction[k]);
		}

		//if (d > dBar)
		{
			;
		}
		*/
	}

	/*

	for (int i = 0; i < plate->contactNum; i++)
	{
		arrayNum = plate->v_ContactElement[i].arrayNum;	

		x_1 = plate->getVertex(plate->v_ContactElement[i].nv_1);
		x_2 = plate->getVertex(plate->v_ContactElement[i].nv_2);

		d = ( x_2 - x_1 ).norm();

		if (d < dBar)
		{
			cout << " contact! " << endl;

			// contact

			dDdEdge = ( x_2 - x_1 ) / d;

			dEnergydD = - 2 * (d - dBar) * log(d / dBar) - (d - dBar) * (d - dBar) / d;

			f = stiffness * dEnergydD * dDdEdge;

			for (int k = 0; k < 3; k++)
			{
				ind = arrayNum(k);
				stepper->addForce(ind, -f[k]);

				ind = arrayNum(k+3);
				stepper->addForce(ind, f[k]);
			}

			// friction

			v_1 = plate->getVelocityOld(plate->v_ContactElement[i].nv_1);
			v_2 = plate->getVelocityOld(plate->v_ContactElement[i].nv_2);

			v_r = v_2 - v_1;

			v_rt = v_r - v_r.dot(dDdEdge) * dDdEdge;

			if ( v_rt.norm() >= epsilonV )
			{
				fVelocity = 1.0;
				
				tK = v_rt / v_rt.norm();
			}
			else
			{
				fVelocity = - ( v_rt.norm() * v_rt.norm() ) / (epsilonV * epsilonV) + 2 * v_rt.norm() / epsilonV;

				tK = v_rt / (v_rt.norm() + 1e-15);
			}

			friction = mu * f.norm() * fVelocity * tK;

			for (int k = 0; k < 3; k++)
			{
				ind = arrayNum(k);
				stepper->addForce(ind, -friction[k]);

				ind = arrayNum(k+3);
				stepper->addForce(ind, friction[k]);
			}

		}
	}

	*/


	/*

	for (int i = 0; i < plate->contactNum; i++)
	{
		arrayNum = plate->v_ContactElement[i].arrayNum;	

		x_1 = plate->getVertexOld(plate->v_ContactElement[i].nv_1);
		x_2 = plate->getVertexOld(plate->v_ContactElement[i].nv_2);

		d = ( x_2 - x_1 ).norm();

		if (d < dBar)
		{
			rNorm = ( x_2 - x_1 ) / d;

			v_1 = plate->getVelocityOld(plate->v_ContactElement[i].nv_1);
			v_2 = plate->getVelocityOld(plate->v_ContactElement[i].nv_2);

			v_r = v_2 - v_1;

			v_rt = v_r - v_r.dot(rNorm) * rNorm;

			if ( v_rt.norm() >= epsilonV )
			{
				fVelocity = 1.0;

				tK = v_rt / v_rt.norm();
			}
			else
			{
				fVelocity = - ( v_rt.norm() * v_rt.norm() ) / (epsilonV * epsilonV) + 2 * v_rt.norm() / epsilonV;

				tK = v_rt / (v_rt.norm() + 1e-15);
			}

			friction = mu * f.norm() * fVelocity * tK;

			for (int k = 0; k < 3; k++)
			{
				ind = arrayNum(k);
				stepper->addForce(ind, -friction[k]);

				ind = arrayNum(k+3);
				stepper->addForce(ind, friction[k]);
			}
		}
	}

	*/
}

void externalContactForce::computeJc()
{

	minDistance = 100000.00;
	contactIndex = 0;

	for (int i = 0; i < plate->contactNum; i++)
	{
		nv_0 = plate->v_ContactElement[i].nv_0;
		edgeLocal = plate->v_ContactElement[i].edgeIndex;

		nv_1 = plate->v_edgeElement[edgeLocal].nv_1;
		nv_2 = plate->v_edgeElement[edgeLocal].nv_2;

		x_0 = plate->getVertex(nv_0);
		x_1 = plate->getVertex(nv_1);
		x_2 = plate->getVertex(nv_2);

		e_1 = x_0 - x_1;
		e_2 = x_0 - x_2;
		e_3 = x_2 - x_1;

		e_4 = e_1.cross(e_2);

		d1 = ( x_0 - x_1 ).norm();
		d2 = ( x_0 - x_2 ).norm();
		d3 = e_4.norm() / e_3.norm();

		d = d1 + d2;

		if (d < minDistance)
		{
			minDistance = d;
			contactIndex = i;
		}
	}

	if (contactIndex != 0 && contactIndex != plate->nv - 2)
	{
		nv_0 = plate->v_ContactElement[contactIndex].nv_0;
		edgeLocal = plate->v_ContactElement[contactIndex].edgeIndex;

		nv_1 = plate->v_edgeElement[edgeLocal].nv_1;
		nv_2 = plate->v_edgeElement[edgeLocal].nv_2;

		x_0 = plate->getVertex(nv_0);
		x_1 = plate->getVertex(nv_1);
		x_2 = plate->getVertex(nv_2);

		e_1 = x_0 - x_1;
		e_2 = x_0 - x_2;
		e_3 = x_2 - x_1;

		e_4 = e_1.cross(e_2);

		normE1 = e_1.norm();
		normE2 = e_2.norm();
		normE3 = e_3.norm();
		normE4 = e_4.norm();

		t_1 = e_1 / normE1;
		t_2 = e_2 / normE2;
		t_3 = e_3 / normE3;
		t_4 = e_4 / normE4;

		d1 = ( x_0 - x_1 ).norm();
		d2 = ( x_0 - x_2 ).norm();
		d3 = e_4.norm() / e_3.norm();

		d = e_4.norm() / e_3.norm();

		gradT3 = (Id3 - t_3 * t_3.transpose()) / normE3;
		gradT4 = (Id3 - t_4 * t_4.transpose()) / normE4;

		de1de1 = - crossMat(e_2) * gradT4 * crossMat(e_2) / normE3;
		de2de2 = - crossMat(e_1) * gradT4 * crossMat(e_1) / normE3;
		de3de3 = - (gradT3 * normE3 - 2 * t_3 * t_3.transpose()) * normE4 / (normE3 * normE3 * normE3);

		de1de2 = (   crossMat(t_4) + crossMat(e_2) * gradT4 * crossMat(e_1) ) / normE3;
		de1de3 = ( - crossMat(e_2) * t_4 * t_3.transpose() / ( normE3 * normE3) ).transpose();
		de2de3 = (   crossMat(e_1) * t_4 * t_3.transpose() / ( normE3 * normE3) ).transpose();

		de2de1 = de1de2.transpose();
		de3de1 = de1de3.transpose();
		de3de2 = de2de3.transpose();

		de1 =   crossMat(e_2) * t_4 / normE3;
		de2 = - crossMat(e_1) * t_4 / normE3;
		de3 = - t_3 * normE4 / (normE3 * normE3);

		dx.segment(0,3) =   de1 + de2;
		dx.segment(3,3) = - de1 - de3;
		dx.segment(6,3) = - de2 + de3;

		dxdx.block(0,0,3,3) = de1de1 + de2de2 + de1de2 + de2de1;
		dxdx.block(3,3,3,3) = de1de1 + de3de3 + de1de3 + de3de1;
		dxdx.block(6,6,3,3) = de2de2 + de3de3 - de2de3 - de3de2;

		dxdx.block(0,3,3,3) = - de1de1 - de1de3 - de2de1 - de2de3;
        dxdx.block(0,6,3,3) = - de2de2 + de1de3 + de2de3 - de1de2;
        dxdx.block(3,6,3,3) = - de3de3 - de1de3 + de1de2 + de3de2;

        dxdx.block(3,0,3,3) = dxdx.block(0,3,3,3).transpose();
        dxdx.block(6,0,3,3) = dxdx.block(0,6,3,3).transpose();
        dxdx.block(6,3,3,3) = dxdx.block(3,6,3,3).transpose();

        dEnergydD = - 2 * d * log( (-d+dBar) / dBar ) - d * d  / (d - dBar);
        d2EnergydD2 = d * d / ((-d+dBar) * (-d+dBar)) + 4 * d / (-d+dBar) - 2 * log( (-d+dBar) / dBar );

        jacobian = stiffness * (d2EnergydD2 * dx * dx.transpose() + dEnergydD * dxdx);

        for (int jj = 0; jj < 9; jj++)
        {
            for (int kk = 0; kk < 9; kk++)
            {
				ind1 = arrayNum(jj);
				ind2 = arrayNum(kk);

				stepper->addJacobian(ind1, ind2, jacobian(jj,kk));
            }
        }
	}


	/*
	minDistance = 100000.00;
	contactIndex = 0;

	for (int i = 0; i < plate->contactNum; i++)
	{
		arrayNum = plate->v_ContactElement[i].arrayNum;	

		x_1 = plate->getVertex(plate->v_ContactElement[i].nv_1);
		x_2 = plate->getVertex(plate->v_ContactElement[i].nv_2);

		d = ( x_2 - x_1 ).norm();

		if (d < minDistance)
		{
			minDistance = d;
			contactIndex = i;
		}
	}
	*/

	/*

	if (minDistance >= dBar / 100)
	{
		arrayNum = plate->v_ContactElement[contactIndex].arrayNum;	

		x_1 = plate->getVertex(plate->v_ContactElement[contactIndex].nv_1);
		x_2 = plate->getVertex(plate->v_ContactElement[contactIndex].nv_2);

		d = ( x_2 - x_1 ).norm();

		//if (d < dBar)
		{
			dDdEdge = ( x_2 - x_1 ) / d;

			d2DdEdge2 = ( Id3 - dDdEdge * dDdEdge.transpose() ) / d;

			dEnergydD = - 2 * d * log( (-d+dBar) / dBar ) - d * d  / (d - dBar);

			d2EnergydD2 = d * d / ((-d+dBar) * (-d+dBar)) + 4 * d / (-d+dBar) - 2 * log( (-d+dBar) / dBar );

			M0 = stiffness * d2EnergydD2 * dDdEdge * dDdEdge.transpose() + stiffness * dEnergydD * d2DdEdge2;

			Jss.block(0,0,3,3) =  - M0;
			Jss.block(3,3,3,3) =  - M0;
			Jss.block(3,0,3,3) =    M0;
			Jss.block(0,3,3,3) =    M0;

			for (int j = 0; j < 6; j++)
        	{
            	for (int k = 0; k < 6; k++)
            	{
					ind1 = arrayNum(j);
					ind2 = arrayNum(k);

					stepper->addJacobian(ind1, ind2, Jss(k,j));
            	}
        	}
		}
	}
	*/

	/*

	for (int i = 0; i < plate->contactNum; i++)
	{
		arrayNum = plate->v_ContactElement[i].arrayNum;	

		x_1 = plate->getVertex(plate->v_ContactElement[i].nv_1);
		x_2 = plate->getVertex(plate->v_ContactElement[i].nv_2);

		d = ( x_2 - x_1 ).norm();

		if (d < dBar)
		{
			dDdEdge = ( x_2 - x_1 ) / d;

			d2DdEdge2 = ( Id3 - dDdEdge * dDdEdge.transpose() ) / d;

			dEnergydD = - 2 * d * log( (-d+dBar) / dBar ) - d * d  / (d - dBar);

			d2EnergydD2 = d * d / ((-d+dBar) * (-d+dBar)) + 4 * d / (-d+dBar) - 2 * log( (-d+dBar) / dBar );

			M0 = stiffness * d2EnergydD2 * dDdEdge * dDdEdge.transpose() + stiffness * dEnergydD * d2DdEdge2;

			Jss.block(0,0,3,3) =  - M0;
			Jss.block(3,3,3,3) =  - M0;
			Jss.block(3,0,3,3) =    M0;
			Jss.block(0,3,3,3) =    M0;

			for (int j = 0; j < 6; j++)
        	{
            	for (int k = 0; k < 6; k++)
            	{
					ind1 = arrayNum(j);
					ind2 = arrayNum(k);

					stepper->addJacobian(ind1, ind2, Jss(k,j));
            	}
        	}
		}
	}

	*/
}

void externalContactForce::setFirstJacobian()
{
	for (int i = 0; i < plate->contactNum; i++)
	{
		nv_0 = plate->v_ContactElement[i].nv_0;
		edgeLocal = plate->v_ContactElement[i].edgeIndex;

		nv_1 = plate->v_edgeElement[edgeLocal].nv_1;
		nv_2 = plate->v_edgeElement[edgeLocal].nv_2;

		arrayNum.setZero(9);

		arrayNum(0) = 3 * nv_0 + 0;
		arrayNum(1) = 3 * nv_0 + 1;
		arrayNum(2) = 3 * nv_0 + 2;
		arrayNum(3) = 3 * nv_1 + 0;
		arrayNum(4) = 3 * nv_1 + 1;
		arrayNum(5) = 3 * nv_1 + 2;
		arrayNum(6) = 3 * nv_2 + 0;
		arrayNum(7) = 3 * nv_2 + 1;
		arrayNum(8) = 3 * nv_2 + 2;

		for (int j = 0; j < 9; j++)
        {
            for (int k = 0; k < 9; k++)
            {
				ind1 = arrayNum(j);
				ind2 = arrayNum(k);

				stepper->addJacobian(ind1, ind2, 1);
            }
        }
	}
}

Matrix3d externalContactForce::crossMat(Vector3d a)
{
	Matrix3d b;

	b<<0,-a(2),a(1),
	a(2),0,-a(0),
	-a(1),a(0),0;

	return b;
}