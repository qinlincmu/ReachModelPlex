#include "../flowstar/Continuous.h"

using namespace flowstar;
using namespace std;


int main(int argc, char *argv[])
{
	unsigned int numVars = 9;

	int x_id = stateVars.declareVar("x");
	int y_id = stateVars.declareVar("y");
        int vx_id = stateVars.declareVar("vx");
        int vy_id = stateVars.declareVar("vy");
	int psi_id = stateVars.declareVar("psi");
        int dpsi_id = stateVars.declareVar("dpsi");
        int a_id = stateVars.declareVar("a");
	int delta_id = stateVars.declareVar("delta");
	int t_id = stateVars.declareVar("t");

        double x_0_min = atof(argv[1]);
        double x_0_max = atof(argv[2]);
        double y_0_min = atof(argv[3]);
        double y_0_max = atof(argv[4]);
        double vx_0_min = atof(argv[5]);
        double vx_0_max = atof(argv[6]);
        double vy_0_min = atof(argv[7]);
        double vy_0_max = atof(argv[8]);
        double psi_0_min = atof(argv[9]);
        double psi_0_max = atof(argv[10]);
        double dpsi_0_min = atof(argv[11]);
        double dpsi_0_max = atof(argv[12]);
        double a_0_min = atof(argv[13]);
        double a_0_max = atof(argv[14]);
        double delta_0_min = atof(argv[15]);
        double delta_0_max = atof(argv[16]);
        std::string a22 = argv[17];
        std::string a24 = argv[18];
        std::string a42 = argv[19];
	std::string a44 = argv[20];
	std::string b2 = argv[21];
	std::string b4 = argv[22];
        double horizon = atof(argv[23]);
        double stepsize = atof(argv[24]);



	// define the dynamics
        // -0.01, 0.01, -0.01, 0.01, 0.6686432942638573, 0.6886432942638573, 0.0, 0.0025686002949781885, -0.010619283308751818, 0.009380716691248183, -0.07192833087518177, -0.05192833087518177, -0.006240490917944199, 0.013759509082055801, -121.36586207534677, -1.7343857500668522, -48.89978484066284, -95.40923313999672, +42.646502835538755, +276.9462365591398, 0.001, 0.001
        // ./RC_bicycle -0.01, 0.01, -0.01, 0.01, 0.0, 0.01, 0.0, 0.01, -0.01, 0.01, -0.01, 0.01, -0.033, -0.013, -125.67, -1.74, -50.63, -98.79, 42.64, 276.94, 0.01, 0.01
	Expression_AST<Real> ode_expression_x("vx * cos(psi) - vy * sin(psi)");//("vx * cos(psi) - vy * sin(psi)");
	Expression_AST<Real> ode_expression_y("vx * sin(psi) + vy * cos(psi)");//("vx * sin(psi) + vy * cos(psi)");
	Expression_AST<Real> ode_expression_vx("a");
       
	Expression_AST<Real> ode_expression_vy(a22+"*vy+"+a24+"*dpsi+"+b2+"*delta");
        //Expression_AST<Real> ode_expression_vy("-121.36586207534677*vy+(-1)*1.7343857500668522*dpsi+42.646502835538755*delta");
        //std::cout << (a22+"*vy"+a24+"*dpsi"+b2+"*delta");
        //Expression_AST<Real> ode_expression_vy("0");
	Expression_AST<Real> ode_expression_psi("dpsi");//("dpsi");
        Expression_AST<Real> ode_expression_dpsi(a42+"*vy+"+a44+"*dpsi+"+b4+"*delta");
        //Expression_AST<Real> ode_expression_dpsi("-48.89978484066284*vy+(-1)*95.40923313999672*dpsi+276.9462365591398*delta");
        //Expression_AST<Real> ode_expression_dpsi("276.946236559139*delta+(-1)*48.89978484066284*vy+(-1)*95.40923313999672*dpsi");
        //Expression_AST<Real> ode_expression_dpsi("0");
        //std::cout << (a42+"*vy"+a44+"*dpsi"+b4+"*delta");
        //Expression_AST<Real> ode_expression_dpsi("-50.63563063493295 * vy - 98.79607250176998 * dpsi + 276.9462365591398 * delta");
        Expression_AST<Real> ode_expression_a("0");
        Expression_AST<Real> ode_expression_delta("0");
	Expression_AST<Real> ode_expression_t("1");



	vector<Expression_AST<Real> > ode_rhs(numVars);
	ode_rhs[x_id] = ode_expression_x;
	ode_rhs[y_id] = ode_expression_y;
        ode_rhs[vx_id] = ode_expression_vx;
        ode_rhs[vy_id] = ode_expression_vy;
	ode_rhs[psi_id] = ode_expression_psi;
        ode_rhs[dpsi_id] = ode_expression_dpsi;
        ode_rhs[a_id] = ode_expression_a;
	ode_rhs[delta_id] = ode_expression_delta;
	ode_rhs[t_id] = ode_expression_t;



	Deterministic_Continuous_Dynamics dynamics(ode_rhs);



	// set the reachability parameters
	Computational_Setting setting;

	// set the stepsize and the order
	setting.setFixedStepsize(stepsize, 5);
//	setting.setFixedStepsize(0.04, 5, 8);
//	setting.setAdaptiveStepsize(0.01, 0.04, 5);

	// set the time horizon
	setting.setTime(horizon);

	// set the cutoff threshold
	setting.setCutoffThreshold(1e-8);

	// set the queue size for the symbolic remainder, it is 0 if symbolic remainder is not used
	setting.setQueueSize(100);

	// print out the computation steps
	setting.printOff();

	// set up the remainder estimation
	Interval I(-0.1, 0.1);
	vector<Interval> remainder_estimation(numVars, I);
	setting.setRemainderEstimation(remainder_estimation);

	// call this function when all of the parameters are defined
	setting.prepare();


	// define the initial set which is a box
	Interval init_x(x_0_min, x_0_max), init_y(y_0_min, y_0_max), init_vx(vx_0_min, vx_0_max), init_vy(vy_0_min, vy_0_max),
                init_psi(psi_0_min, psi_0_max), init_dpsi(dpsi_0_min, dpsi_0_max), init_a(a_0_min, a_0_max), init_delta(delta_0_min, delta_0_max), init_t;

	vector<Interval> initial_box(numVars);
	initial_box[x_id] = init_x;
	initial_box[y_id] = init_y;
        initial_box[vx_id] = init_vx;
        initial_box[vy_id] = init_vy;
	initial_box[psi_id] = init_psi;
        initial_box[dpsi_id] = init_dpsi;
        initial_box[a_id] = init_a;
	initial_box[delta_id] = init_delta;
	initial_box[t_id] = init_t;

	Flowpipe initialSet(initial_box);


	// empty unsafe set
	vector<Constraint> unsafeSet;
        //Constraint cons1(formula);
        //Constraint cons2("x_id-1");
        //unsafeSet.push_back(cons1);
        //cout << formula << endl;
        //cons2.push_back();

	/*
	 * The structure of the class Result_of_Reachability is defined as below:
	 * nonlinear_flowpipes: the list of computed flowpipes
	 * tmv_flowpipes: translation of the flowpipes, they will be used for further analysis
	 * fp_end_of_time: the flowpipe at the time T
	 */
	Result_of_Reachability result;

	// run the reachability computation
	clock_t begin, end;
	begin = clock();

	dynamics.reach(result, setting, initialSet, unsafeSet);
//	dynamics.reach(result, setting, result.fp_end_of_time, unsafeSet);
//	dynamics.reach(result, setting, result.fp_end_of_time, unsafeSet);

	end = clock();
	//printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	// flowpipes should be translated to single Taylor model vectors before plotting
	result.transformToTaylorModels(setting);
                

	//Plot_Setting plot_setting;
	//plot_setting.printOff();
	//plot_setting.setOutputDims(x_id, y_id);
        //int indicator = 1;
	//int indicator = plot_setting.plot_2D_interval_MATLAB("RC_bicycle", result);//Qin
        //plot_setting.plot_2D_interval_MATLAB("RC_bicycle", result);

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = result.tmv_flowpipes.begin();
	std::list<Flowpipe>::const_iterator fpIter = result.nonlinear_flowpipes.begin();
	std::vector<unsigned int> varIDs({0, 1, 2, 3, 4, 5});
	
	for(; tmvIter != result.tmv_flowpipes.end() ; ++tmvIter, ++fpIter)
	{
		std::vector<Interval> box;
		tmvIter->intEval(box, fpIter->domain, varIDs);

		Interval X = box[0], Y = box[1], VX = box[2], VY = box[3], PSI = box[4], DPSI = box[5]; 
		std::cout << X.inf() << "," << X.sup() << ",";
		std::cout << Y.inf() << "," << Y.sup() << ",";
                std::cout << VX.inf() << "," << VX.sup() << ",";
                std::cout << VY.inf() << "," << VY.sup() << ",";
		std::cout << PSI.inf() << "," << PSI.sup()<< ",";
                std::cout << DPSI.inf() << "," << DPSI.sup() << ";";
	}
	return 0;
}
