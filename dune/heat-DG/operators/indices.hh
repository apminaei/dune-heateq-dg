/*
 * indices.hh
 *
 * 			DONE
 */

class Indices{
private:

public:



	/* PRIMARY VARIABLES */  // NOTE: The order of the indices must be the same as the order in the GFS, initial and boundary conditions
	static const int numOfPVs 	 = 2 ;
	static const int PVId_T  	= 0;
	static const int PVId_P  	= 1;
	

	/* BOUNDARY CONDITION */
	static const unsigned int numOfVelBCs 	= 1; // V_g , V_w
	static const unsigned int BCId_heat		= 0;
	static const unsigned int BCId_dirichlet		= 1;
	static const unsigned int BCId_neumann 	 		= 0;
	static const unsigned int BCId_depressurization = 2;



};


