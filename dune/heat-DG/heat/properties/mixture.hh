/* ALL PARAMETERS ARE NONDIMENSIONAL */
template<typename PTree>
class Mixture
{
private:
	const PTree& ptree;
	Parameters<PTree> parameter;
	Methane<PTree> methane;
	Water<PTree> water;
	Salt salt;
	CharacteristicValues X_c;

public:
	  //! construct from grid view
	Mixture (const PTree& ptree_  )
	: ptree(ptree_),
	parameter(ptree_),
		water(ptree_),
		methane(ptree_)
	{}

	/* MOLE FRACTIONS ( X -> liquid ; Y -> gas ) */

	std::vector<double> EquilibriumMoleFractions( double T/*K*/, double Pg/*Pa*/, double Xc, double z )const{

		// double S = Xc * (salt.MolarMass()/methane.MolarMass());
		// double f_CH4 = z*Pg/(methane.SolubilityCoefficient(T,S)*X_c.P_c);
		// double f_H2O = Pg/(water.SaturatedVaporPressure( T,S )*X_c.P_c);
		
		

		double S = Xc * (salt.MolarMass()/water.MolarMass());
		double f_CH4 = z*Pg/methane.SolubilityCoefficient(T,S)/X_c.P_c;
		double f_H2O = Pg/water.SaturatedVaporPressure( T,S )/X_c.P_c;
		double c = z*water.SaturatedVaporPressure( T,S )*X_c.P_c-methane.SolubilityCoefficient(T,S)*X_c.P_c ;
		double c_yh2o =	z*water.SaturatedVaporPressure( T,S )*X_c.P_c-methane.SolubilityCoefficient(T,S)*X_c.P_c*water.SaturatedVaporPressure( T,S )*X_c.P_c/Pg*(1. - Xc);
		double c_xch4 =	-z*Pg+z*water.SaturatedVaporPressure( T,S )*X_c.P_c*(1. - Xc);

		double Y_H2O = c_yh2o/c;
		double X_CH4 = c_xch4/c;
		double Y_CH4 = X_CH4*methane.SolubilityCoefficient(T,S)*X_c.P_c/(z*Pg);
		double X_H2O = Y_H2O * f_H2O;

		// double Y_H2O = ((1.-Xc)-f_CH4)/(f_H2O-f_CH4);//std::max(0., std::min(1., ));
		// double Y_CH4 = (1.-Y_H2O);//std::max(0., std::min(1., ));
		// double X_H2O = Y_H2O * f_H2O;
		// double X_CH4 = (1. - Xc - X_H2O);//std::max(0., std::min(1., ));

		std::vector<double> X(Indices::numOfComps,0.);
		X[Indices::compId_XCH4] = X_CH4;
		X[Indices::compId_XH2O] = X_H2O;
		X[Indices::compId_YCH4] = Y_CH4;
		X[Indices::compId_YH2O] = Y_H2O;

		return X;
	}

	double YCH4( double X_CH4, double T, double Pg, double Xc, double z )const{

		// NOTE: it is not necessary to check case1,2 for fncs f_CH4 and f_H2O because the cases are already determined within classes CH4 and H2O.
		double S = Xc * (salt.MolarMass()/water.MolarMass());
		double Y_CH4 = X_CH4 * (methane.SolubilityCoefficient(T,S)*X_c.P_c) / ( z * Pg ) ;
		return Y_CH4;
	}

	double XH2O( double Y_H2O, double T, double Pg, double Xc )const{

		// NOTE: it is not necessary to check case1,2 for fncs f_CH4 and f_H2O because the cases are already determined within classes CH4 and H2O.
		double S = Xc * (salt.MolarMass()/water.MolarMass());
		double X_H2O = Y_H2O * Pg / (water.SaturatedVaporPressure( T,S )*X_c.P_c);
		return X_H2O;
	}

	/* MASS TRANSFER COEFFICIENTS */

	double DiffCoeffH2OInGas( double T, double Pg ) const {

		double D; /* m^2/s */

#ifdef STATEINDEPENDENTPROPERTIES
		double T_ref = parameter.RefT();
		double P_ref = parameter.RefP();
		T = T_ref;
		Pg = P_ref;
#endif

		double a0 = 0. ;
		double a1 = 2.26e-9;
		double a2 = 0.002554;

		D = 1.e0*( a0 + a1*T + a2/Pg ) ;

		return D/X_c.dispersivity_c ;

	}

	double DiffCoeffCH4InLiquid( double T, double Pw ) const {

		double D; /* m^2/s */

#ifdef STATEINDEPENDENTPROPERTIES
		double T_ref = parameter.RefT();
		double P_ref = parameter.RefP();
		T = T_ref;
		Pw = P_ref;
#endif

		double A = 0.003475 ; 	/* K */
		double B = 1.e0*1.57e-5;		/* cm^2/s */

		D = B * exp(-A/T) * 1.0e-6;//pow((Pw/1.e5),2.) * B * exp(-A/T) * 1.0e-6;

		return D/X_c.dispersivity_c;
	}

};
