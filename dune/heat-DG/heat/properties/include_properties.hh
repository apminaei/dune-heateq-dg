// CONSTITUTIVE AND MATERIAL PROPERTIES
// #include"salt.hh"
// #include"H2O.hh"
// #include"CH4.hh"
// #include"eosCH4.hh"
// #include"hydrate.hh"
// #include"mixture.hh"
// #include"soil.hh"
// #include"hydraulic_properties.hh"
// #include"hydrate_phase_change.hh"

template<typename GV, typename PTree>
class Properties
{
private:
	  const GV& gv;
	  const PTree& ptree;
	//   double *time;
	//   double *dt;

	  const static int dim = GV::dimension;
	  constexpr static double eps = 1.e-6;
	  
public:

	//PARAMETERS AND PROPERTIES
  	Indices index;
  	CharacteristicValues characteristicValue;
  	MeshParameters<PTree> mesh;
  	Parameters<PTree> parameter;
//   	Methane<PTree> gas;
// // #ifdef STATEINDEPENDENTPROPERTIES
// //   	BaseEoS<PTree> eos;
// // #elif defined(PENG_ROBINSON_EOS)
//   	PengRobinson<PTree> eos;
// // #endif
//   	Water<PTree> water;
//   	Salt salt;
//   	Mixture<PTree> mixture;
//   	Hydrate<GV, PTree> hydrate;
//   	Soil<GV,Parameters<PTree>> soil;
//   	HydraulicProperties<GV,Parameters<PTree>> hydraulicProperty;
//   	HydratePhaseChangeKinetics<GV,PTree> kinetics;
  	
  	//! construct from grid view
  	Properties ( const GV& gv_ , 
  				 const PTree& ptree_)
	: gv( gv_ ),
	  ptree(ptree_),
	  mesh(ptree_),
	  parameter(ptree_)
	  
  	{}
	
	double dt_initial = ptree.get("time.dt_initial",(double)1.);
	int time_interval = ptree.get("output.time_interval",(double)1);
	/******************************************************************************/

  	void ReportStatistics( std::string file_name,
  						   double time /*s*/,
						   double dt /*s*/,
						   int total_newton_iterations,
						   double newton_first_defect,
						   double newton_reduction,
						   double clock_time_elapsed /*s*/) {

  		std::fstream result;

  		if(time == 0. ){
  			result.open(file_name, std::fstream::out | std::fstream::trunc);
  			result	<< "time [s]" << '\t'
  					<< "dt [s]"	<< '\t'
					<< "total no. of newton iterations" << '\t'
					<< "  newton first defect" << '\t'
					<< "  newton reduaction" << '\t'
					<< "clock time [s]" << '\t'
  					<< mesh.X_cells << "*" << mesh.Z_cells << '\t'
					<< std::endl;
  			result.close();
  		}

  		result.open(file_name, std::fstream::app);
  		double t_new = time+dt;

		result	<< time	<< '\t'
				<< dt	<< '\t'
				<< total_newton_iterations << '\t'
				<< newton_first_defect << '\t'
				<< newton_reduction << '\t'
				<< clock_time_elapsed
				<< std::endl;
		result.close();
  	}

	/******************************************************************************/
  	void ReportParameters( std::string file_name,
  						   std::string method_g,
							std::string method_w,
							std::string method_T,
							std::string method_x,
							std::string method_y,
							double alpha_g , double alpha_w, double alpha_s, double alpha_T, double alpha_x, double alpha_y,
							double dissCoeff, double formCoeff,
							double X_P, double X_G, double X_M, double X_H, double X_D, double X_S) {

  		std::fstream result;

  		
  			result.open(file_name, std::fstream::out | std::fstream::trunc);
  			result	<< "penalty coeff.  " << alpha_g << '\t' << alpha_w << '\t'<< alpha_s << '\t'<< alpha_T << '\t'<< alpha_x << '\t'<< alpha_y << '\n'
					<< " S=-1, N=0, I=1,  " << method_g << '\t'<< method_w << '\t'<< method_T << '\t'<< method_x << '\t'<< method_y << '\n'
					<<  " dissCoeff " << dissCoeff << '\t' << " formCoeff " << formCoeff<< '\n'
					<< "Xc_Permeability= " << X_P << '\t' << "Xc_gravity= "<< X_G << '\t'<< "X_source_mass= " << X_M << '\t'
					<< "X_source_heat= " <<X_H << '\t'<<"X_dispersivity= " << X_D << '\t'<< "X_specific_heat= " << X_S  
					<< std::endl;
  			result.close();
  		

  	}


  	void ReportNewton( std::string file_name,
  						   double time /*s*/,
						   double dt /*s*/,
						   int newton_iterations,
						   Dune::BlockVector<Dune::FieldVector<double, 1>> newton_defects,
						   double u_normtwo, double u_normone, double u_infnorm) {

  		std::fstream result;

  		
  			result.open(file_name, std::fstream::out | std::fstream::trunc);
  			result	<< " time [s] = " << time << '\t'
  					<< " dt [s] = "	<< dt << '\t'
					<< " total no. of newton iterations = " << newton_iterations << '\t'
					<< " Norm two= " << u_normtwo << '\t'
  					<< " Norm one  = "	<< u_normone << '\t'
					<< " Norm inf = " << u_infnorm << '\t'
  					<< mesh.X_cells << "*" << mesh.Z_cells << '\n';
					
  			//result.close();
  		

  		//result.open(file_name, std::fstream::app);
  		double t_new = time+dt;
		for (int i=0; i < newton_iterations+1; i++){
			result	<< i	<< '\t'
					<< newton_defects[i] << '\n';	
					
		}
		// result << std::endl;
		
		result.close();
  	}


	/******************************************************************************/

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

};
