/* ALL VALUES ARE NONDIMENSIONAL */
template<typename GV,typename Properties>
class ProblemBoundaryConditions
{
private :
	const GV& gv ;
	const Properties& property;
	const static int dim = GV::dimension;
	Indices indices;
	CharacteristicValues characteristicValues;
	ProblemInitialConditions<GV,Properties> icvalue;
	// double time_fraction = property.parameter.time_end() / 31.536e6; 
	double Xc_time = 1. / (36.*24.*36. * 1.e3);
	double press_rate = 1.e-3;

public :

	// ! construct from gridview
	ProblemBoundaryConditions (const GV& gv_,const Properties& property_)
	: gv ( gv_ ),
	  property(property_),
	  icvalue(gv_, property_)
	{}
	
	/* NOTE:
	 * dirichlet: BC_w -> Pw, BC_g -> Sg, 
	 * neumann: total mass flux (convective+diffusive) BC_w -> f_w, BC_g -> f_g
	 * */
	
	/* boundary types */
	template<typename I> 
	int
	type( I& intersection,/*const typename GV::Traits::template Codim<0>::Entity& element,*/
			const Dune::FieldVector<double,dim-1>& xlocal,
		  double time /*ndims*/,
		  double dt /*ndims*/ ) const {
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);
		
		int bctype = indices.BCId_neumann;
		// if(property.mesh.isLeftBoundary(globalPos))// || property.mesh.isTopBoundary(globalPos))
			bctype = indices.BCId_dirichlet;
		
		
		return bctype;
	}

	/* boundary values */
	template<typename I> 
	double
	value ( I& intersection,/*const typename GV::Traits::template Codim<0>::Entity& element,*/
			const Dune::FieldVector<double,dim-1>& xlocal,
			double time/*s*/,
			double dt/*s*/ ) const {
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);
		
		// References to inside and outside cells
		const auto &cell_inside = intersection.inside();
	   	auto icv /*ndim*/ = icvalue.evaluate(cell_inside,iplocal);
		auto bcvalue = 0.;
		// if(property.mesh.isLeftBoundary(globalPos))// || property.mesh.isTopBoundary(globalPos))
			bcvalue = icv;//+2 * std::sin(3.14 * (globalPos[0] * globalPos[1]) * (time+dt));
		
		
		
		return bcvalue;
	}

	// BC regarding phase Velocities
	template<typename I> 
	int
	velType( I& intersection,/*const typename GV::Traits::template Codim<0>::Entity& element,*/
			const Dune::FieldVector<double,dim-1>& xlocal,
		  double time /*ndims*/,
		  double dt /*ndims*/ ) const {
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);

		
		int bctype = indices.BCId_neumann;

		// if(property.mesh.isLeftBoundary(globalPos) );//|| property.mesh.isTopBoundary(globalPos))
			bctype = indices.BCId_dirichlet;
		
		return bctype;
	}

	template<typename I> 
	double
	velValue ( I& intersection,/*const typename GV::Traits::template Codim<0>::Entity& element,*/
			const Dune::FieldVector<double,dim-1>& xlocal,
			double time/*s*/,
			double dt/*s*/ ) const {
		auto iplocal = intersection.geometryInInside().global(xlocal);
		auto globalPos = intersection.inside().geometry().global(iplocal);
		
		// References to inside and outside cells
		const auto &cell_inside = intersection.inside();
		auto icv /*ndim*/ = icvalue.evaluate(cell_inside,iplocal);
		auto bcvalue = 0.;//icv;//2 * std::sin(3.14 * (globalPos[0] * globalPos[1]) * (time+dt));;
		
		// if(property.mesh.isLeftBoundary(globalPos) )//|| property.mesh.isTopBoundary(globalPos))
			bcvalue = icv ;//+ 2 * std::sin(3.14 * (globalPos[0] * globalPos[1]) * (time+dt));
		
		return bcvalue;
	}

	// ! get a reference to the gridview
	inline const GV& getGridView () { return gv ; }
};


	

