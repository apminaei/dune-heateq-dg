/* ALL PARAMETERS ARE NONDIMENSIONAL */
template<typename PTree>
class Water
{
private:
	CharacteristicValues characteristicValue;
	const PTree& ptree;
	Parameters<PTree> parameter;
public:

	Water (const PTree& ptree_  )
	 : ptree(ptree_),
	 parameter(ptree_)
  	{}

	double CriticalTemperature( ) const {
		return 647.096 ; /* [K] */
	}

	double CriticalPressure( ) const {
		return 22.064 * 1.0e6 ; /* [Pa] */
	}

	double MolarMass( ) const {
		/* unit -> kg/mol */
		return 18.0/1000;
	}

	double Density( double T, double Pw, double S ) const {

		double rho;
		/* rho: unit -> kg/m^3 */

		/*
		 * averages values & expansion coefficients: ρ0=1027 kg/m^3,  T0=10°C,  S_0=35 g/kg
		 * Thermal expansion: \alpha_T=0.15 kg/(m^3 °C)
		 * Salinity contraction: \alpha_S=0.78 kg/(m^3 g/kg)
		 * Pressure compressibility: \alpha_P=0.0045 kg/(m^3 dbar)
		 * UNESCO EOS-80 : Equation of state for seawater
		 * We use a linear EOS (web ref:http://www.ccpo.odu.edu/~atkinson/OEAS405/Chapter2_Stratified_Ocean/Lec_04_DensityEOS.pdf)
		 */

		double rho_0 = 1027.0;
		double T_0 = 10.;
		double S_0 = 0.03115;
		double alpha_T = -0.15;
		double alpha_S = 0.78*1e3;//108.524;//0.96706917808*1e2;//
		double alpha_P = 0.0045;

// #ifdef STATEINDEPENDENTPROPERTIES
// 		double T_ref = parameter.RefT();
// 		double P_ref = parameter.RefP();
// 		double S_ref = parameter.RefSal();
// 		T  = T_ref;
// 		Pw = P_ref;
// 		S  = S_ref;
// #endif

		rho = rho_0
			+ (   alpha_P*(Pw*1.e-4)
				+ alpha_T*((T-273.15)-T_0)
				+ alpha_S*(S-S_0)
			  );

		return rho/characteristicValue.density_c;
	}

	double MolarDensity(double T, double Pw, double S)const{
		return Density( T,Pw,S)*characteristicValue.density_c/MolarMass();
	}

	double DynamicViscosity( double T, double Pw, double S ) const {
		double mu;
		/* mu: unit -> Pa.s */

		// REFERENCE:
		double mu_0 = 0.001792 ; // kg/m/s
		double a = - 1.94 ;
		double b = - 4.80 ;
		double c =  6.74 ;
		double T0 = 273.15 ; // K

#ifdef STATEINDEPENDENTPROPERTIES
		double T_ref = parameter.RefT();
		double Tr = T0/T_ref ;
#else
		double Tr = T0/T ;
#endif

		mu = mu_0 * exp( a + b * Tr + c * Tr*Tr );

		return mu/characteristicValue.viscosity_c;
	}

	double ThermalConductivity( double T, double Pw, double S ) const {

		double kth;
		/* kth: unit -> W.m^-1 K^-1 */

#ifdef STATEINDEPENDENTPROPERTIES
		double T_ref = parameter.RefT();
		double P_ref = parameter.RefP();
		//double S_ref = parameter.RefSal();
		T  = T_ref;
		Pw = P_ref;
		//S  = S_ref;
#endif

		kth = 0.57153*( 1 + 0.003*(T-273.15) - 1.025e-5*(T-273.15)*(T-273.15) + 6.53e-10*Pw - 0.29*S );// 0.0245650
		// std::cout << kth << std::endl;
		// exit(0);
		return kth/characteristicValue.thermalconductivity_c ;
	}

	double Cp( double T, double Pw, double S ) const {
		double Cp;
		/* Cp: unit -> J*kg^-1*K^-1 */

		Cp = 3945.0 ;

		return Cp/characteristicValue.specificheat_c;

	}

	double Cv( double T, double Pw, double S ) const {
		double Cv;
		/* mu: unit -> J*kg^-1*K^-1 */

		Cv = Cp( T, Pw, S )*characteristicValue.specificheat_c;

		return Cv/characteristicValue.volumetricheat_c;

	}

	double SaturatedVaporPressure( double T /*K*/,double S ) const {

		double psat;   /* [Pa] */

// #ifdef STATEINDEPENDENTPROPERTIES
// 		double T_ref = parameter.RefT();
// 		T = T_ref;
// #endif

		// REF: SUGAR TOOLBOX

		double Pc = CriticalPressure(); // in Pa
		double Tc = CriticalTemperature(); // in K
		double Tr = T/Tc;

		double c1 = -7.85951783;
		double c2 = 1.84408259;
		double c3 = -11.7866497;
		double c4 = 22.6807411;
		double c5 = -15.9618719;
		double c6 = 1.80122502;

		double lnppc = 1./Tr * (  c1*(1-Tr)
								+ c2*pow((1-Tr),1.5)
								+ c3*pow((1-Tr),3)
								+ c4*pow((1-Tr),3.5)
								+ c5*pow((1-Tr),4)
								+ c6*pow((1-Tr),7.5) );

		psat = Pc * exp(lnppc);   /* [Pa] */

		return psat/characteristicValue.P_c;
	}

};
