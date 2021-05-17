class Salt
{
private:
	CharacteristicValues characteristicValue;

public:

	double MolarMass( ) const {
		/* unit -> kg/mol */
		return 58.4/1000;
	}

	double DiffCoeff( double T, double Pw ) const {

		double D = 1.e0*1.0e-9;	/* m^2/s */

		return D/characteristicValue.dispersivity_c ;
	}

	double Source( ) const {
		return 0.;///characteristicValue.X_source_mass;
	}

};
