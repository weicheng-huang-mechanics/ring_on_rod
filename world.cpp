#include "world.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				
	saveData = m_inputData.GetBoolOpt("saveData");			
	deltaTime = m_inputData.GetScalarOpt("deltaTime");     
	totalTime = m_inputData.GetScalarOpt("totalTime");    
	YoungM = m_inputData.GetScalarOpt("YoungM");
	density = m_inputData.GetScalarOpt("density");
	rodRadius = m_inputData.GetScalarOpt("rodRadius");
	Possion = m_inputData.GetScalarOpt("Possion");
	stol = m_inputData.GetScalarOpt("stol");
	forceTol = m_inputData.GetScalarOpt("forceTol");
	scaleRendering = m_inputData.GetScalarOpt("scaleRendering");
	maxIter = m_inputData.GetIntOpt("maxIter");
	gVector = m_inputData.GetVecOpt("gVector");
	viscosityA = m_inputData.GetScalarOpt("viscosityA");
	viscosityB = m_inputData.GetScalarOpt("viscosityB");
	deltaLength = m_inputData.GetScalarOpt("deltaLength");
	pointMass = m_inputData.GetScalarOpt("pointMass");
	pointSize = m_inputData.GetScalarOpt("pointSize");
	dBar = m_inputData.GetScalarOpt("dBar");
	stiffness = m_inputData.GetScalarOpt("stiffness");
	mu = m_inputData.GetScalarOpt("mu");
	epsilonV = m_inputData.GetScalarOpt("epsilonV");
	nv = m_inputData.GetIntOpt("nv");
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) return;
	
	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	// Open an input file named after the current time
	ostringstream name;
	name.precision(6);
	name << fixed;

    name << "datafiles/simDiscreteNet";
    name << "_nv_" << nv;
    name << "_delta_" << deltaTime;
    name << "_dbar_" << dBar;
    name << "_epsilonV_" << epsilonV * 1e6;
    name << ".txt";

	outfile.open(name.str().c_str());
	outfile.precision(10);	
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}

	/*

	if ( timeStep % 10 != 0)
	{
		return;
	}

	for (int i = 0; i < plate->nv; i++)
	{
		Vector3d xCurrent = plate->getVertex(i);

		outfile << currentTime << " " << xCurrent(0) << " " << xCurrent(2) << endl;
	}

	*/

	Vector3d xEnd = plate->getVertex(plate->nv - 1);
	//Vector3d xVec = plate->getVelocity(plate->nv - 1);

	outfile << currentTime << " " << xEnd(0) << " " << xEnd(2) << endl;

	/*

	//if (timeStep == Nstep)
	{
		Vector3d xCenter = plate->getVertex(0);

		outfile << currentTime << " " << xCenter(2) << endl;
	}

	*/
}

void world::setPlateStepper()
{
	// Create the plate 
	plate = new elasticPlate(YoungM, density, rodRadius, Possion, deltaTime, pointMass, nv);

	plateBoundaryCondition();

	plate->setup();

	stepper = new timeStepper(*plate, viscosityB);

	// set up force
	m_inertialForce = new inertialForce(*plate, *stepper);
	m_gravityForce = new externalGravityForce(*plate, *stepper, gVector);
	m_stretchForce = new elasticStretchingForce(*plate, *stepper);
	m_bendingForce = new elasticBendingForce(*plate, *stepper);
	m_dampingForce = new dampingForce(*plate, *stepper, viscosityA);
	m_externalContactForce = new externalContactForce(*plate, *stepper, dBar, stiffness, mu, epsilonV);

	plate->updateTimeStep();

	// set up first jacobian
	m_inertialForce->setFirstJacobian();
	m_stretchForce->setFirstJacobian();
	m_bendingForce->setFirstJacobian();
	m_dampingForce->setFirstJacobian();
	m_externalContactForce->setFirstJacobian();

	stepper->first_time_PARDISO_setup();

	// time step 
	Nstep = totalTime / deltaTime;
	timeStep = 0;
	currentTime = 0.0;
}

void world::plateBoundaryCondition()
{
	plate->setVertexBoundaryCondition(plate->getVertex(0), 0);
	plate->setVertexBoundaryCondition(plate->getVertex(plate->nv - 2), plate->nv - 2);
}

void world::updateTimeStep()
{
	bool goodSolved = false;

	while (goodSolved == false)
	{
		// Start with a trial solution for our solution x
		plate->updateGuess(); // x = x0 + u * dt

		updateEachStep();

		goodSolved = true;
	}

	plate->updateTimeStep();

	if (render) 
	{
		cout << "time: " << currentTime << endl;
	}

	currentTime += deltaTime;
		
	timeStep++;

	stepper->clearElasticJacobian();
	m_stretchForce->computeJsViscus();
	m_bendingForce->computeJbViscus();
	//m_externalContactForce->computeJcViscus();
}

void world::updateEachStep()
{
	double normf = forceTol * 10.0;
	double normf0 = 0;
	
	bool solved = false;
	
	int iter = 0;
		
	while (solved == false)
	{
		plate->prepareForIteration();

		stepper->setZero();

		m_inertialForce->computeFi();
		m_gravityForce->computeFg();
		m_stretchForce->computeFs();
		m_bendingForce->computeFb();
		m_dampingForce->computeFd();
		m_externalContactForce->computeFc();
		stepper->addDampingForce();

		normf = stepper->GlobalForceVec.norm();

		if (iter == 0) 
		{
			normf0 = normf;
		}
		
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}

		normf = 0.0;
		
		if (solved == false)
		{
			m_inertialForce->computeJi();
			m_gravityForce->computeJg();
			m_stretchForce->computeJs();
			m_bendingForce->computeJb();
			m_dampingForce->computeJd();
			//m_externalContactForce->computeJc();
			stepper->addDampingJacobian();

			stepper->integrator(); // Solve equations of motion
			plate->updateNewtonMethod(stepper->GlobalMotionVec); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}

	if (render)
	{
		cout << "iter " << iter << endl;
	}
}

int world::simulationRunning()
{
	if (timeStep < Nstep) 
	{
		return 1;
	}
	else 
	{
		return -1;
	}
}

Vector3d world::getScaledCoordinate(int i, int j)
{
	Vector3d xCurrent;
	
	if (j == 0)
	{
		xCurrent = plate->v_edgeElement[i].x_1 * scaleRendering;
	}
	if (j == 1)
	{
		xCurrent = plate->v_edgeElement[i].x_2 * scaleRendering;
	}

	return xCurrent;
}

Vector3d world::getScaledPoint()
{
	return plate->getVertex(plate->nv - 1) * scaleRendering;
}

int world::numStretchingPair()
{
	return plate->edgeNum;
}