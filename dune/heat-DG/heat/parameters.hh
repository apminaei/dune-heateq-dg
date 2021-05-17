/*
 * parameters.hh
 *
 * 
 */

template<typename PTree>
class Parameters
{
private:
	const PTree& ptree;
	const double pi = 3.14159265358979323846;
	const static int dim = MeshParameters<PTree>::dimension;
	constexpr static double eps = 1.0e-6;
	CharacteristicValues X_c;
	double *time;
	double *dt;

	double T_t0;
	double Sg_t0;
	double Pw_t0;
	double Pw_Well;
	double Pg_t0;
	double Sh_t0;
	double XCH4_t0;
	double YH2O_t0;
	double XC_t0;
	double Sg_x0;
	double Pw_x0;
	double Sgin_x0;
	double gradPx;
	double gradPz;
	double gradTx;
	double gradTz;

	double t_END;
	double T_end;
	double kd;
	double kf;
	int numMaterials;
	int numProps;
	std::vector<std::vector<double> > prop;

	double ref_salinity;
	double ref_saltconcentration;
	double ref_temperature;
	double ref_pressure;

	bool gravity_flag;
	double g_magnitude;

public:

	MeshParameters<PTree> mesh;
	//! constructor
	Parameters (const PTree& ptree_) 
	:ptree(ptree_),
    mesh(ptree_)
	{
			Sg_t0 = ptree.get("initial.Sg",(double)0.0);
			Pw_t0 = ptree.get("initial.Pw",(double)2.e6);
			Pg_t0 = ptree.get("initial.Pg",(double)2.0848e6);
			T_t0 = ptree.get("initial.T",(double)4.) ; // in Celcius
			Sh_t0 = ptree.get("initial.Sh",(double)0.3);
			YH2O_t0 = ptree.get("initial.YH2O",(double)0.0005);
			XCH4_t0 = ptree.get("initial.XCH4",(double)0.);
			XC_t0 = ptree.get("initial.XC",(double)5.5e-3);
			t_END = ptree.get("time.time_end",(double)2.16e6);
			T_end = ptree.get("time.T_end",(double)2.16e6);
			Pw_Well = ptree.get("boundary.Pw_at_Well",(double)8.e6);
			Pw_x0 = ptree.get("boundary.Pw_at_left",(double)2.e6);
			Sgin_x0 = ptree.get("boundary.Sg_at_inlet",(double)0.0);
			gradPx = ptree.get("grad.gradPx",(double)0.0);
			gradPz = ptree.get("grad.gradPz",(double)0.0);
			gradTx = ptree.get("grad.gradTx",(double)0.0);
			gradTz = ptree.get("grad.gradTz",(double)0.035);

			numMaterials = ptree.get("sediment.number_of_materials",(int)1);
			numProps = 8;
			prop = std::vector<std::vector<double> > (numMaterials,std::vector<double>(numProps, 0.));
			for(int n_mat=0; n_mat<numMaterials; n_mat++ ){
				std::string name = "sediment.material"+std::to_string(n_mat);
				prop[n_mat][0] = ptree.get(name+".por",	(double)0.5);
				prop[n_mat][1] = ptree.get(name+".K",	(double)1.e-15);
				prop[n_mat][2] = ptree.get(name+".pentry",(double)5.e4);
				prop[n_mat][3] = ptree.get(name+".lambda",(double)1.2);
				prop[n_mat][4] = ptree.get(name+".swr",	(double)0.);
				prop[n_mat][5] = ptree.get(name+".sgr",	(double)0.);
				prop[n_mat][6] = ptree.get(name+".m",	(double)1.);
				prop[n_mat][7] = ptree.get(name+".beta",(double)1.);
			}
			

			//reference state
			ref_salinity = ptree.get("reference_state.salinity",(double)0.);
			ref_saltconcentration = ref_salinity * (18.0/58.4); /*MolarMass_H2O/MolarMass_salt*/
			ref_temperature = (ptree.get("reference_state.temperature",(double)0.));/*K*/
			ref_pressure = ptree.get("reference_state.pressure",(double)1.01e5);


			kd = ptree.get("hydrate_phase_change.dissociation_rate",(double)1.e-14);/*mol/m².Pa.s*/
			kf = ptree.get("hydrate_phase_change.formation_rate",(double)1.e-13);/*mol/m².Pa.s*/
			
			

			//gravity
			gravity_flag = ptree.get("gravity.flag",(bool)true);
			g_magnitude = ptree.get("gravity.magnitude",(double)9.81);
	}

	/**********************************************************************
	 * INPUTS
	 **********
	 * z_domain : height of the computational domain [m]
	 * z_cells	: no. of cells along Z-axis
	 * 
	 *
	 *
	 *
	 **********************************************************************/

	//2. Initial Values

	double InitialSg(Dune::FieldVector< double,dim > xglobal) const {
		double Sg = Sg_t0;
		return Sg;
	}

	double InitialPw(Dune::FieldVector< double,dim > xglobal) const {
		double Pw = Pw_t0 + 1030.21 * 9.81 * (0.-xglobal[1])*X_c.x_c;
		return Pw; /* Pa */
	}

	double Pw_at_Well(Dune::FieldVector< double,dim > xglobal) const {
		double Pw = Pw_Well;
		return Pw; /* Pa */
	}

	double InitialPg(Dune::FieldVector< double,dim > xglobal) const {
		double Pg = Pg_t0;
		return Pg;/* Pa */
	}

	double InitialT(Dune::FieldVector< double,dim > xglobal) const {
		double T = T_t0 + ref_temperature  ; /*K*/
		return T; /* K */
	}
	

	double InitialSh(Dune::FieldVector< double,dim > xglobal) const {
		double Sh = 0.0 ;
		double GHSZ_width = mesh.Z_GHSZ_top - mesh.Z_GHSZ_bottom;
		
		if( mesh.isGHSZ(xglobal)){
			Sh = 1.2 * (xglobal[1]-mesh.Z_GHSZ_bottom)/GHSZ_width * (xglobal[1]-mesh.Z_GHSZ_top)/(-GHSZ_width);//* (rand()%2) + 0.001;//
		}
		return Sh;
	}

	double InitialXCH4(Dune::FieldVector< double,dim > xglobal) const {
		double XCH4 = XCH4_t0;
		return XCH4;
	}
	double InitialYH2O(Dune::FieldVector< double,dim > xglobal) const {
		double YH2O = YH2O_t0;
		return YH2O;
	}

	double InitialXC(Dune::FieldVector< double,dim > xglobal) const {
		
		double XC = XC_t0;
		return XC;
	}

	//3. Boundary values
	double InletSg(Dune::FieldVector< double,dim > xglobal) const {
		double Sg = Sgin_x0;
		return Sg;
	}

	double LeftPw(Dune::FieldVector< double,dim > xglobal) const {
		double Pw = Pw_x0;
		return Pw; /* Pa */
	}


	//4. Material properties
	std::vector< std::vector<double> > layer_properties() const {
		return prop;
	}
	//numMaterials
	int num_materials() const {
		return numMaterials;
	}
	/**********************************************************************/
	/* REFERENCE STATE VALUES */
	double ReferenceSalinity() const {
		double refSal = ref_salinity;
		return refSal;
	}
	double ReferenceTemperature() const {
		double ref = ref_temperature;
		return ref; /*K*/
	}
	double ReferencePressure() const {
		double ref = ref_pressure;
		return ref; /*Pa*/
	}
	double HydrateDissociationRateConstant() const {
		return kd ;//* (t_END/T_end) ;
	}
	double HydrateFormationRateConstant() const {
		return kf ;//* (t_END/T_end);
	}
	double time_end() const {
		return t_END;
	}
	double T_END() const {
		return T_end;
	}

#ifdef STATEINDEPENDENTPROPERTIES
	// (P,T,sal) reference state for CASE1
	double RefP() const {
		return 10.* ReferencePressure(); /*Pa*/
	}

	double RefT() const {
		return ReferenceTemperature() + 9.0; /*K*/
	}

	double RefSal() const {
		return 0.035; /*kg/kg*/;
	}
#endif

    /**********************************************************************/
	Dune::FieldVector<double,dim>
	SedimentVelocity ( double time, double dt ) const {

		Dune::FieldVector<double,dim> vs( 0. );
		vs[dim-1] = 0.;
		vs[0] = 0.;

		return vs; /*m/s*/
	}

	/* GRAVITY VECTOR */
	Dune::FieldVector<double,dim>
	g( ) const {
		Dune::FieldVector<double,dim> gravity( 0. );
		double g = 0.;
		if(gravity_flag) g = -g_magnitude;
		gravity[dim-1] = g;
		gravity[0] = 0.;
		return gravity; /*N/kg*/
	}

};

