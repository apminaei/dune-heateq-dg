/* ALL PARAMETERS ARE NONDIMENSIONAL */
template<typename GV, typename PTree>
class Hydrate {

private:
	CharacteristicValues characteristicValue;
	const PTree& ptree;
	const GV& gv;
	Parameters<PTree> parameter;
public:

	Hydrate (const GV& gv_, const PTree& ptree_  )
	 : ptree(ptree_), gv(gv_),
	 parameter(ptree_)
  	{}

	double Density( ) const {
		/* unit -> kg/m^3 */

		double rho_h = 920.0 ;

		return rho_h/characteristicValue.density_c;

	}

	double MolarMass() const
	{
		/* unit -> kg/mol */
		return 119.5/1000;
	}

	double MolarDensity() const {
		return Density()*characteristicValue.density_c/MolarMass();
	}

	double HydrationNumber() const
	{
		return 5.90;
	}

	double ThermalConductivity( double T, double P ) const
	{
		double kth;
		/* kth: unit -> W.m^-1 K^-1 */

		kth = 0.5 ;

		return kth/characteristicValue.thermalconductivity_c;
	}

	double Cp( double T, double P ) const
	{
		double Cp;
		/* Cp: unit -> J/kg.K */

#ifdef STATEINDEPENDENTPROPERTIES
		double T_ref = parameter.RefT();
		T = T_ref;
#endif

		Cp = ( 1.9370547e-05*T*T*T - 1.5151760e-02*T*T + 3.9553876*T - 342.70565 )*1.0e3;


	
		return Cp/characteristicValue.specificheat_c;

	}

	double Cv( double T, double P ) const
	{
		double Cv;
		/* Cv: unit -> J/kg.K */

		Cv = Cp( T, P )*characteristicValue.specificheat_c ;

		return Cv/characteristicValue.volumetricheat_c;

	}
};