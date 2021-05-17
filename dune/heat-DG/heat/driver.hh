/*
 * driver.hh
 *
 */

#ifndef HEAT_DG_HH_
#define HEAT_DG_HH_


template <class GV, class PTree>
void driver(const GV &gv, // GridView
			const PTree& ptree, //ParameterTree 
			Dune::MPIHelper& helper //MPI-Helper 
		   )
{

	//	INCLUDE EXTERNAL CLASSES


	typedef Properties<GV, PTree> Properties;

	//	CHOOSE DOMAIN AND RANGE FIELD TYPE
	typedef typename GV::Grid::ctype Coord;
	typedef double Real;
	const int dim = GV::dimension;
	Real time = 0.0;
	Real dt = 0.0;
	
	Properties property(gv,ptree);



	/* Non-dimensionalize time prams */
	/*in input file dt, dt_initprob, t_END, t_END_initprob, dt_min, dt_max, t_OP are specified in kilo-annum*/
	//INITIAL PROBLEM time and dt
	auto Xc_t = property.characteristicValue.t_c;
	auto t_day_sec = 1.;
	// double dt_initprob  = ptree.get("initial_problem.dt_initial",(double)0.0001);
	// dt_initprob *= (1000.*364.*24.*60.*60.); /*convert to seconds*/
	// dt_initprob *= 1./Xc_t; /*ndim*/
	// double t_END_initprob  = ptree.get("initial_problem.time_end",(double)1);
	// t_END_initprob *= (1000.*364.*24.*60.*60.); /*convert to seconds*/
	// t_END_initprob 	 *= 1./Xc_t; /*ndim*/
	//MAIN PROBLEM time and dt
	dt  = ptree.get("time.dt_initial",(double)0.001);
	dt *= t_day_sec/Xc_t; /*convert to seconds then nondim*/
	double t_END  = ptree.get("time.time_end",(double)300);
	t_END *= t_day_sec/Xc_t; /*convert to seconds*/
	//t_END += t_END_initprob; //Total simulation time
	// output time interval
	double t_OP   = ptree.get("output.time_interval",(double)1000);
	t_OP *= t_day_sec/Xc_t; /*convert to seconds*/
	//adaptive time control
	bool adaptive_time_control = ptree.get("adaptive_time_control.flag",(bool)true);
	double dt_min = ptree.get("adaptive_time_control.dt_min",(double)1.e-6);
	dt_min *= t_day_sec/Xc_t; /*convert to seconds*/
	double dt_max = ptree.get("adaptive_time_control.dt_max",(double)10.);
	dt_max *= t_day_sec/Xc_t; /*convert to seconds*/


	
	// Real T_end = property.parameter.T_END();
	// double time_fraction = t_END / T_end;
	
	
	Real dtstart = dt;
	Real time_op = time;
	Real clock_time_elapsed = 0.;

	int maxAllowableIterations = ptree.get("adaptive_time_control.max_newton_steps",(int)10);
	int minAllowableIterations = ptree.get("adaptive_time_control.min_newton_steps",(int)4);
	int max_linear_iteration = ptree.get("newton.MaxLinearIteration",(int)10);
	
	const int degree_T = 1;
	const int degree_P = 0;

	//	GFS
#ifdef PARALLEL
	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON0;
#else
	typedef Dune::PDELab::ConformingDirichletConstraints CON0;	// pure Neumann: no constraints
#endif									
	using VBE0 = Dune::PDELab::ISTL::VectorBackend<> ;	// default block size: 1
	//typedef OPBLocalFiniteElementMap<Coord,Real,degree_P,dim,Dune::GeometryType::simplex > OPBSim;
	
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_T, dim, Dune::PDELab::QkDGBasisPolynomial::legendre> FEM_T;
	// typedef Dune::PDELab::PkLocalFiniteElementMap<GV, Coord, Real, degree_T  > FEM_T;
	FEM_T fem_T;
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_P, dim, Dune::PDELab::QkDGBasisPolynomial::lagrange> FEM_P; 
	// typedef Dune::PDELab::PkLocalFiniteElementMap<GV, Coord, Real, degree_P  > FEM_P;
	FEM_P fem_P;

	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_T, CON0, VBE0> GFS_T; // gfs
	GFS_T gfs_T(gv, fem_T);
	
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_P, CON0, VBE0> GFS_P; // gfs
	GFS_P gfs_P(gv, fem_P);

	using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;//  block size -> numOfPVs
	// gfs for composite system Pw,Sg,Sh,T,XCH4,YH2O,XC
	typedef Dune::PDELab::CompositeGridFunctionSpace<VBE, Dune::PDELab::EntityBlockedOrderingTag, GFS_T, GFS_P> GFS;
	GFS gfs(gfs_T,  gfs_P);
	
    typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
    CC cc;


	//	MAKE VECTOR CONTAINER FOR THE SOLUTION
	using U = Dune::PDELab::Backend::Vector<GFS, double>;
	U uold(gfs, 0.0);
	U unew(gfs, 0.0);

	auto T_Initiallamda = [&](const Dune::FieldVector<double,dim>& x){return 0.;};
	auto T = Dune::PDELab::makeGridFunctionFromCallable(gv, T_Initiallamda);
	auto P_Initiallamda = [&](const Dune::FieldVector<double,dim>& x){return 1.e0;};
	auto P = Dune::PDELab::makeGridFunctionFromCallable(gv, P_Initiallamda);
	using InitialTypelamda = Dune::PDELab::CompositeGridFunction<decltype(T), decltype(P)> ;
	InitialTypelamda initialamda(T,P);
	Dune::PDELab::interpolate(initialamda, gfs, uold); // Initialize the solution at t=0 (uold) with the given initial values
	
	// Dune::PDELab::interpolate(T, gfs_T, uold); // Initialize the solution at t=0 (uold) with the given initial values
	
	using PathP = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_P>>;
    using SUBGFS_P = Dune::PDELab::GridFunctionSubSpace<GFS,PathP>;
    SUBGFS_P    subgfs_P(gfs);
	using PathT = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_T>>;
    using SUBGFS_T = Dune::PDELab::GridFunctionSubSpace<GFS,PathT>;
    SUBGFS_T    subgfs_T(gfs);

    
	//	MAKE INSTATIONARY GRID OPERATOR SPACE
	double method_g = 0.;//ConvectionDiffusionDGMethod::SIPG;
	double method_w = 0.;//ConvectionDiffusionDGMethod::SIPG;
	double method_T = 0.;//ConvectionDiffusionDGMethod::NIPG;
	double method_x = 0.;//ConvectionDiffusionDGMethod::NIPG;
	double method_y = 0.;//ConvectionDiffusionDGMethod::NIPG;
	double alpha_g = 1.e0;//* property.characteristicValue.X_source_mass;
	double alpha_w = 1.e0;//* property.characteristicValue.P_c;
	double alpha_s = 1.e0;
	double alpha_T = 1.e0;//* property.characteristicValue.T_c;
	double alpha_x = 1.e0 ;//* property.characteristicValue.X_source_mass;
	double alpha_y = 1.e0;
	double intorder=2; 

	typedef ProblemBoundaryConditions<GV,Properties> BoundaryConditions ;
	BoundaryConditions bc( gv,property ) ;

	typedef LocalOperator<GV, Properties, U, GFS, FEM_T, FEM_P> LOP; // spatial part
	LOP lop(gv, property, &unew, gfs, &time, &dt, intorder, method_g, method_w, method_T, method_x, method_y, alpha_g, alpha_w, alpha_s, alpha_T, alpha_x, alpha_y);

	typedef TimeOperator<GV, Properties> TLOP; // temporal part
	TLOP tlop(gv, property, intorder);

	typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
	MBE mbe(100);

	typedef Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, Real, Real, Real, CC, CC> GOLOP;
	GOLOP goLOP(gfs, cc, gfs, cc, lop, mbe);

	// How well did we estimate the number of entries per matrix row?
	// => print Jacobian pattern statistics0
	typename GOLOP::Traits::Jacobian jac(goLOP);
	std::cout << " LOP DONE ! " << std::endl;

	typedef Dune::PDELab::GridOperator<GFS, GFS, TLOP, MBE, Real, Real, Real, CC, CC> GOTLOP;
	GOTLOP goTLOP(gfs, cc, gfs, cc, tlop, mbe);

	typedef Dune::PDELab::OneStepGridOperator<GOLOP, GOTLOP> IGO;
	IGO igo(goLOP, goTLOP);
	std::cout << " IGO DONE ! " << std::endl;

	// SELECT A LINEAR SOLVER BACKEND
#ifdef PARALLEL

	// typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS, CC> LS; //doesn't work
	// LS ls(gfs, cc, 100, 2);

	// typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC> LS;
	// LS ls(gfs, cc, max_linear_iteration, 5,1);

	typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO> LS; //works
	LS ls(gfs, max_linear_iteration, 1, true, true);
	/* 	NOTES:
		LINEAR SOLVER STATISTICS
		res.iterations = i;
		res.reduction = def/def0;
		res.conv_rate  = pow(res.reduction,1.0/i);
		res.elapsed = watch.elapsed();
		// final print
		if (_verbose>0)
		{
			std::cout << "=== rate=" << res.conv_rate
					<< ", T=" << res.elapsed
					<< ", TIT=" << res.elapsed/i
					<< ", IT=" << i << std::endl;
		}
	*/

	

	// using LS = Dune::PDELab::ISTLBackend_OVLP_AMG_4_DG< IGO, CC, GFS, CC,
	// Dune::PDELab::CG2DGProlongation,
	// Dune::SeqSSOR,
	// Dune::BiCGSTABSolver>;
	// LS ls(igo, cc, gfs, cc, max_linear_iteration,2,true,true);
	// // set parameters for AMG in CG-subspace
	// Dune::Amg::Parameters params = ls.parameters();
	// params.setCoarsenTarget(10000);
	// params.setMaxLevel(20);
	// params.setProlongationDampingFactor(0.1);
	// params.setNoPreSmoothSteps(3);
	// params.setNoPostSmoothSteps(3);
	// params.setGamma(1);
	// params.setAdditive(false);
	// ls.setParameters(params);
	// typedef Dune::PDELab::ISTLBackend_BCGS_AMG_ILU0<IGO> LS;
	// LS ls(gfs,max_linear_iteration,1,true,true);

	// typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILUn<GFS, CC> LS; //works
	// LS ls(gfs, cc);

	// typedef Dune::PDELab::ISTLBackend_CG_AMG_SSOR<IGO> LS; // should be checked
	// int verbose = 0;
	// if (gfs.gridView().comm().rank() == 0)
	// 	verbose = 1;
	// LS ls(gfs, 100, verbose);
	auto param = ls.parameters();
	//param.setMaxLevel(3); // max number of coarsening levels
	param.setCoarsenTarget(10000000); // max DoF at coarsest level
	ls.setParameters(param);

	std::cout << " PARALLEL LS DONE ! " << std::endl;

#else
	// typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
	// LS ls(100, true);

	// using LS = Dune::PDELab::ISTLBackend_SEQ_AMG_4_DG<IGO,GFS,Dune::PDELab::CG2DGProlongation,Dune::SeqSSOR,Dune::BiCGSTABSolver>; //works
	// LS ls(igo,gfs,max_linear_iteration,2,true,true);
	// // // set parameters for AMG in CG-subspace
	// Dune::Amg::Parameters params = ls.parameters();
	// params.setCoarsenTarget(10000);
	// params.setMaxLevel(20);
	// params.setProlongationDampingFactor(0.1);
	// params.setNoPreSmoothSteps(3);
	// params.setNoPostSmoothSteps(3);
	// params.setGamma(1);
	// params.setAdditive(false);
	// ls.setParameters(params);


	typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
	LS ls;
	std::cout << " LS DONE ! " << std::endl;
#endif

    //    SELECT SOLVER FOR NON-LINEAR PROBLEM
    // using PDESOLVER = Dune::PDELab::Newton< IGO, LS, U >;
    // PDESOLVER pdesolver( igo, ls );
    // // select control parameters for non-linear PDE-solver
	// pdesolver.setLineSearchStrategy(ptree.get("newton.line_search_strategy",(std::string)"noLineSearch"));//Strategy {  hackbuschReusken, hackbuschReuskenAcceptBest }
    // //pdesolver.setLineSearchStrategy(PDESOLVER::Strategy::hackbuschReuskenAcceptBest);
	// pdesolver.setReassembleThreshold(0.0);
    // pdesolver.setVerbosityLevel(2);
    // pdesolver.setReduction(ptree.get("newton.reduction",(double)1e-5));
    // pdesolver.setMinLinearReduction(ptree.get("newton.min_linear_reduction",(double)1.e-9));
	// pdesolver.setMaxIterations(ptree.get("newton.max_iterations",(int)15));
    // pdesolver.setForceIteration(ptree.get("newton.force_iterations",(bool)true));
	// pdesolver.setAbsoluteLimit(ptree.get("newton.abs_error",(double)1.e-4)); 

	//TODO: CHECK NEW NEWTON PARAMS
	//	SELECT SOLVER FOR NON-LINEAR PROBLEM
	typedef Dune::PDELab::NewtonMethod<IGO, LS> PDESOLVER;
	PDESOLVER pdesolver(igo, ls);

	// 	select control parameters for non-linear PDE-solver
	// typedef Dune::PDELab::LineSearchNone<PDESOLVER> lineSearchStrategy;
	// typedef Dune::PDELab::LineSearchHackbuschReusken<PDESOLVER> lineSearchStrategy;
	// lineSearchStrategy linesearchstrategy(pdesolver);
	pdesolver.setParameters(ptree.sub("newton"));
	//linesearchstrategy.setParameters(ptree.sub("newton.line_search"));
	// pdesolver.setVerbosityLevel(ptree.get("newton.VerbosityLevel",(int)3));
	// pdesolver.setReduction(ptree.get("newton.Reduction",(double)1e-5));
    // // pdesolver.setMinLinearReduction(ptree.get("newton.min_linear_reduction",(double)1.e-9));
	// pdesolver.setAbsoluteLimit(ptree.get("newton.AbsoluteLimit",(double)1.e-4)); 
	
	std::cout << " PDESOLVER DONE ! " << std::endl;

	// SELECT TIME-STEPPER
	Dune::PDELab::ImplicitEulerParameter<Real> method1;
	Dune::PDELab::OneStepThetaParameter<Real> method2(0.5); //Crank-Nicholson -> 0.5, Implicit Euler -> 1.0, Explicit Euler -> 0.0
	Dune::PDELab::Alexander2Parameter<Real> method3;
	Dune::PDELab::Alexander3Parameter<Real> method4;
	Dune::PDELab::FractionalStepParameter<Real> method5;
	Dune::PDELab::HeunParameter<Real> method6;
	Dune::PDELab::Shu3Parameter<Real> method7;
	Dune::PDELab::RK4Parameter<Real> method8; // didnot work

	Dune::PDELab::TimeSteppingParameterInterface<Real> *pmethod = &method1;

	Dune::PDELab::OneStepMethod<Real, IGO, PDESOLVER, U, U> osm(*pmethod, igo, pdesolver);
	osm.setVerbosityLevel(2);
	std::cout << " OSM DONE ! " << std::endl;

	//	GRAPHICS FOR INITIAL GUESS
	// primary variables
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_P, U> DGF_P;
	DGF_P dgf_P(subgfs_P, uold);
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_T, U> DGF_T;
	DGF_T dgf_T(subgfs_T, uold);

	//	VTK
	std::string fileName = ptree.get("output.file_name",(std::string)"test");
	std::string pathName = ptree.get("output.path_name",(std::string)"test");
	pathName += "outputs/";
	pathName += fileName ;
	std::string fileNameDefects = fileName;
	std::string pathNameDefects = pathName+"/"+fileName;
	std::time_t now = std::time(0);
	struct tm *tstruct;
	char buf [80];
	std::time(&now);
	tstruct = std::localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d-%H-%M-%S", tstruct );
	std::string timeStr(buf); 
	pathNameDefects.append("_"+timeStr);
    std::string jacPath = pathNameDefects;
		jacPath +="/";
		jacPath +=fileName;

	if(helper.rank()==0){
		
		std::filesystem::create_directory(pathNameDefects);
	}
	// pathName += "/" ;
	// const std::string str = "";
	//Dune::PDELab::FilenameHelper fn(pathName + fileName);

	int subsampling = 1;
	Dune::RefinementIntervals RefInt(subsampling);

	using VTKWRITER = Dune::SubsamplingVTKWriter<GV> ;
	VTKWRITER vtkwriter(gv, RefInt, false, Dune::VTK::Precision::float32); // vtk nonconforming
	using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV> ;
	VTKSEQUENCEWRITER vtkSequenceWriter(std::make_shared<VTKWRITER>(vtkwriter), fileName, pathName, "");

	// add data field for all components of the space to the VTK writer
	// primary variables
	
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_T>>(dgf_T, "T"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_P>>(dgf_P, "P"));

	vtkSequenceWriter.write(time, Dune::VTK::appendedraw);
	vtkSequenceWriter.clear();
	
	std::string dgmethod_g = std::__cxx11::to_string(method_g);
	std::string dgmethod_w = std::__cxx11::to_string(method_w);
	std::string dgmethod_T = std::__cxx11::to_string(method_T);
	std::string dgmethod_x = std::__cxx11::to_string(method_x);
	std::string dgmethod_y = std::__cxx11::to_string(method_y);
	double dissCoeff = property.parameter.HydrateDissociationRateConstant();
	double formCoeff = property.parameter.HydrateFormationRateConstant();
	if(helper.rank()==0){
		std::string parameters_file = pathNameDefects;
		parameters_file +="/";
		parameters_file +=fileName;
		parameters_file +="_parameters";
		parameters_file += ".txt";
		property.ReportParameters( 	parameters_file,
									dgmethod_g, dgmethod_w, dgmethod_T, dgmethod_x, dgmethod_y,
									alpha_g, alpha_w, alpha_s, alpha_T, alpha_x, alpha_y, dissCoeff, formCoeff,
									property.characteristicValue.permeability_c,
									property.characteristicValue.X_gravity,
									property.characteristicValue.X_source_mass,
									property.characteristicValue.X_source_heat,
									property.characteristicValue.dispersivity_c,
									property.characteristicValue.specificheat_c);
	}

	//	INITIALIZE
	unew = uold;
	int opcount = 1;
	double timecount = time;
	double dtLast = dtstart;
	int dtFlag = 0;

	bool exceptionCaught = false;

	int newton_iterations = 0;
	double newton_first_defect = 0.;
	double newton_defect = 0.;
	//Dune::BlockVector<Dune::FieldVector<Real, 1>> newton_defects;
	
	
	//	BEGIN TIME LOOP
	while ( time < (t_END - 1e-3/Xc_t))
	{
		if( exceptionCaught==false ){
				dt = std::max(dt,dt_min);
		}

		if(helper.rank()==0){
			std::cout<< "_____________________________________________________" <<std::endl;
			std::cout<< " current opcount = " << opcount - 1 << std::endl;
		}

		clock_t start = clock();
		try{
			if(helper.rank()==0){
			std::cout<<"****************************" << std::endl;
			std::cout<<"  CALLING osm.apply() !"	  << std::endl;
			std::cout<<"****************************" << std::endl;
			}

			osm.apply( time, dt, uold, unew );
			
			newton_iterations = osm.getPDESolver().result().iterations;
			newton_first_defect = osm.getPDESolver().result().first_defect;
			newton_defect = osm.getPDESolver().result().defect;
            auto jacobian = osm.getPDESolver().getJacobian();
			if(helper.rank()==0 &&  (opcount%10==0)){
			Dune::writeMatrixToMatlab ( Dune::PDELab::Backend::native(jacobian), jacPath+"jacobian");
			Dune::writeVectorToMatlab(Dune::PDELab::Backend::native(unew),jacPath+"solution");
			}
			auto newton_defects = osm.getPDESolver().result().defects;
			auto u_norm_two = unew.two_norm();
			auto u_norm_one = unew.one_norm();
			auto u_norm_infinity = unew.infinity_norm();
			if(helper.rank()==0 &&  (newton_iterations>1)){//((time+dt )/(t_OP * opcount) > (1.-1.e-6)) and ((time+dt ) / (t_OP * opcount)< (1. + 1.e-6))
				std::string s_OP = std::__cxx11::to_string(time);
				std::string parameters_file = pathNameDefects;
				parameters_file +="/";
				parameters_file +=fileNameDefects;
				parameters_file +="_";
				parameters_file +=s_OP;
				parameters_file += ".txt";
				property.ReportNewton( parameters_file,
							time /*s*/,
							dt /*s*/,
							newton_iterations, newton_defects,
							u_norm_two, u_norm_one,  u_norm_infinity);
			}
			//std::cout << "  ================////////=============== " <<  newton_defects[newton_iterations] << std::endl;
			/*
				std::cout << "  Newton iteration " << std::setw(2)(Sets the field width to be used on output operations) << this->res_.iterations
                           << ".  New defect: "
                           << std::setw(12) << std::setprecision(4) << std::scientific (Sets the floatfield format flag for the str stream to scientific)
                           << this->res_.defect
                           << ".  Reduction (this): "
                           << std::setw(12) << std::setprecision(4) << std::scientific
                           << this->res_.defect/this->prev_defect_
                           << ".  Reduction (total): "
                           << std::setw(12) << std::setprecision(4) << std::scientific
                           << this->res_.reduction << std::endl;
			*/

			exceptionCaught = false;

		}catch ( Dune::Exception &e ) {
			exceptionCaught = true;
			if( (dt) > (1e-3/Xc_t) ){

				if(helper.rank()==0){
					std::cout << "Catched Error, Dune reported error: " << e << std::endl;
				}

				unew = uold;

				newton_iterations = 0;

				dt *= 0.5;
				dtLast = dt;
					continue;
			}
			else
			{
				if(helper.rank()==0){
					std::cout << "ABORTING, due to DUNE error: " << e << std::endl;
				}
				exit(0);
			}
		}
		clock_t end = clock();
		double clock_time_this_step = (double) (end-start) / CLOCKS_PER_SEC;
		clock_time_elapsed += clock_time_this_step;

		if(helper.rank()==0){
			std::cout<<"DONE"<<std::endl;
			std::cout<<"_____________________________________________________"<<std::endl;
		}

		/*********************************************************************************************
		 * OUTPUT
		 *********************************************************************************************/
		/* At each time step: **Statistics**
		 * t_new,
		 * dt,
		 * newton iterations per fixed point iteration,
		 * total newton iterations
		 */
		
		if(helper.rank()==0){
			std::string statistics_file = pathNameDefects;
			statistics_file +="/";
			statistics_file +=fileName;
			statistics_file +="_statistics";
			statistics_file += ".txt";
			property.ReportStatistics( 	statistics_file,
										time,
										dt,
										newton_iterations,
										newton_first_defect,
										newton_defect,
										clock_time_elapsed );
		}
		// GRAPHICS FOR NEW OUTPUT
		// primary variables
		
		DGF_T dgf_T(subgfs_T, unew);
		DGF_P dgf_P(subgfs_P, unew);

		/*********************************************************************************************
			 * OUTPUT
			 *********************************************************************************************/
		if (((time + dt)/(t_OP * opcount) > (1.-1.e-3)) and ((time + dt) / (t_OP * opcount)< (1. + 1.e-3)))
		{
			// primary variables
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_T>>(dgf_T, "T"));

			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_P>>(dgf_P, "P"));

			vtkSequenceWriter.write(time, Dune::VTK::appendedraw);
			vtkSequenceWriter.clear();
			if(helper.rank()==0){
				std::cout<< " ******************************************************************* " << std::endl;
				std::cout<< " OUTPUT WRITTEN " << opcount << " ----processor: " << helper.rank() << std::endl;
				std::cout<< " ******************************************************************* " << std::endl;
				std::cout<< std::flush;
			}
			timecount = time;
			opcount = opcount + 1;
		}

		//		PREPARE FOR NEXT TIME INTEGRATION
		//		1. ASSIGN THE 'NEXT' VALUE TO 'OLD' VARIABLE
		uold = unew;
		//		2. ADVANCE TIME:
		time += dt;
		if(helper.rank()==0){
			std::cout<<" "<< std::endl;
			std::cout<< " time = " << time*Xc_t ;
			std::cout<< std::flush;
		}
		if (adaptive_time_control)
		{
			if (newton_iterations > maxAllowableIterations)
			{
				dt = std::max(dt*0.9 , dt_min);
			}
			else if (newton_iterations <= minAllowableIterations)
			{
				dt = std::min(dt * 1.2, dt_max);
			}
			if (dtFlag == -1)
			{
				dt = dtLast;//std::max(dt, dtLast);
			}
			dtFlag = 0;
		}
		else
		{
			dt = dtstart;
		}
		if(helper.rank()==0){
			std::cout << " , time+dt = " << (time + dt)*Xc_t
					  << " , opTime = "  << t_OP * opcount * Xc_t ;
			std::cout<< std::flush;
		}
		dtLast = dt;
		if ((time + dt) > (t_OP * opcount + 1.e-6) and time < (t_OP * opcount - 1.e-6))
		{
			
			dt = t_OP * opcount - time;

			if(helper.rank()==0){
				std::cout<< " , because timeNext > opNext , dt set to : " << dt*Xc_t << std::endl;
				std::cout<< std::flush;
			}
			dtFlag = -1;
		}

		if(helper.rank()==0){
			std::cout<< " , dt  : " << dt*Xc_t << std::endl;
			std::cout<<" "<< std::endl;
			std::cout << " READY FOR NEXT ITERATION. " << std::endl;
			std::cout<< std::flush;
		}
	}
};

#endif /* HEAT_DG_HH_ */