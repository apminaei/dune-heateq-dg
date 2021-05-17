/* All Values are dimensional and transfer to nondim in Initial.hh */
template<typename GV, typename Properties>
class ProblemInitialConditions
{
private:
	const GV& gv;
	const Properties& property;
	const static int dim = GV::dimension;
	constexpr static double eps = 1.e-6;


public:

		//! construct from grid view
		ProblemInitialConditions(const GV& gv_, const Properties& property_)
		: 	gv( gv_ ),
			property(property_)
		{}

		/* Initial Conditions */
		double
		evaluate (const typename GV::Traits::template Codim<0>::Entity& element,
			  	const Dune::FieldVector<double,dim>& xlocal) const {

			auto xglobal = element.geometry().global(xlocal);

			/******************************************************************************/
			
			double T = property.parameter.InitialT(xglobal)/  property.characteristicValue.T_c; /*K*/
			// std::cout << " T = " << T << std::endl;
			//exit(0);
			/******************************************************************************/
			// MOLE FRACTIONS
			
			
		  /******************************************************************************/

		  return T; /*ndim*/
	  	}

	  //! get a reference to the grid view
	  inline const GV& getGridView () {return gv;}
};

