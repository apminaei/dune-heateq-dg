using namespace Dune::PDELab;

template <class GV, typename Params>
class TimeOperator:
		public Dune::PDELab::NumericalJacobianApplyVolume<TimeOperator<GV, Params>>,
		public Dune::PDELab::NumericalJacobianVolume<TimeOperator<GV, Params>>,
		public Dune::PDELab::FullVolumePattern,
		public Dune::PDELab::LocalOperatorDefaultFlags,
		public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
	const GV &gv;
	const Params&	  property;

	constexpr static double eps = 1.0e-6;
	constexpr static double pi = atan(1.) * 4;
	unsigned int intorder;
	double Xc_K;
	double Xc_mu;
	double Xc_rho;
	double Xc_kth;
	double Xc_C;
	double Xc_P;
	double Xc_T;
	double Xc_X;
	double Xc_Y;
	double T_ref;


public:
	// pattern assembly flags
	enum{doPatternVolume = true	};
	// residual assembly flags
	enum{doAlphaVolume = true	};

  	typedef typename GV::IndexSet IndexSet;

	// constructor remembers parameters
	TimeOperator(const GV &gv_, const Params&	 property_, unsigned int intorder_ = 4)
		:gv(gv_), property( property_ ),intorder(intorder_)
	{
		Xc_K = property.characteristicValue.permeability_c;
		Xc_mu = property.characteristicValue.viscosity_c;
		Xc_rho = property.characteristicValue.density_c;
		Xc_kth = property.characteristicValue.thermalconductivity_c;
		Xc_C = property.characteristicValue.specificheat_c;
		Xc_P = property.characteristicValue.P_c;
		Xc_T = property.characteristicValue.T_c;
		Xc_X = property.characteristicValue.x_c;
		Xc_Y = property.characteristicValue.x_c;
		T_ref = property.parameter.ReferenceTemperature()/Xc_T;
	}
	// volume integral depending on test and ansatz functions
	template <typename EG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume(const EG &eg, const LFSU &lfsu, const X &x, const LFSV &lfsv, R &r) const
	{

		// subspaces
		//Temperature
		//Temperature
		const auto &lfsv_T = lfsv.template child<Indices::PVId_T>();
		const auto &lfsu_T = lfsu.template child<Indices::PVId_T>();
		
		// define types
		
    	using RF = double;
		
		// Reference to cell
	  	const auto& cell = eg.entity();
		const IndexSet &indexSet = gv.indexSet();
		int cell_number = indexSet.index(cell);

		// Get geometry
		auto geo = eg.geometry();

		// Transformation matrix
		typename EG::Geometry::JacobianInverseTransposed jac;

		// loop over quadrature points
		for (const auto &ip : quadratureRule(geo, intorder))
		{
			// evaluate test shape functions
			std::vector<Dune::FieldVector<double,1>> phi_T(lfsu_T.size());
			lfsu_T.finiteElement().localBasis().evaluateFunction(ip.position(), phi_T);
			std::vector<Dune::FieldVector<double,1>> psi_T(lfsv_T.size());
			lfsv_T.finiteElement().localBasis().evaluateFunction(ip.position(), psi_T);

			
			// evaluate T
			RF T = 0.0;
			for (int i = 0; i < lfsu_T.size(); i++)
				T += x(lfsu_T, i) * phi_T[i];

			auto factor = ip.weight() * geo.integrationElement(ip.position());
			for (int i = 0; i < lfsv_T.size(); i++)
			{
				r.accumulate(lfsv_T, i, (T-T_ref) * psi_T[i] * factor);
			}

		} 	//End Quadrature Rule
	}	// End of alpha volume
	
 	// jacobian contribution of volume term
	// template <typename EG, typename LFSU, typename X, typename LFSV, typename M>
	// void jacobian_volume(const EG &eg, const LFSU &lfsu, const X &x, const LFSV &lfsv, M& mat) const
	// {

	// 	// subspaces
	// 	//Gas pressure
	// 	const auto &lfsv_Pw = lfsv.template child<Indices::PVId_Pw>();
	// 	const auto &lfsu_Pw = lfsu.template child<Indices::PVId_Pw>();
	// 	//Water Saturation
	// 	const auto &lfsv_Sg = lfsv.template child<Indices::PVId_Sg>();
	// 	const auto &lfsu_Sg = lfsu.template child<Indices::PVId_Sg>();

	// 	//Hydrate Saturation
	// 	const auto &lfsv_Sh = lfsv.template child<Indices::PVId_Sh>();
	// 	const auto &lfsu_Sh = lfsu.template child<Indices::PVId_Sh>();

	// 	//Temperature
	// 	const auto &lfsv_T = lfsv.template child<Indices::PVId_T>();
	// 	const auto &lfsu_T = lfsu.template child<Indices::PVId_T>();

	// 	//Methane mole fraction
    // 	const auto &lfsv_XCH4 = lfsv.template child<Indices::PVId_XCH4>();
    // 	const auto &lfsu_XCH4 = lfsu.template child<Indices::PVId_XCH4>();

    // 	//Water mole fraction
    // 	const auto &lfsv_YH2O = lfsv.template child<Indices::PVId_YH2O>();
    // 	const auto &lfsu_YH2O = lfsu.template child<Indices::PVId_YH2O>();
		
	// 	//Salt mole fraction
    // 	const auto &lfsv_XC = lfsv.template child<Indices::PVId_C>();
    // 	const auto &lfsu_XC = lfsu.template child<Indices::PVId_C>();

	// 	// define types
	// 	using RF = typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::FiniteElementType::
	// 		Traits::LocalBasisType::Traits::RangeFieldType;
	// 	using RangeType = typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::FiniteElementType::
	// 		Traits::LocalBasisType::Traits::RangeType;
	// 	using size_type = typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::SizeType;
	
	// 	// Reference to cell
	//   	const auto& cell = eg.entity();
	// 	const IndexSet &indexSet = gv.indexSet();
	// 	int cell_number = indexSet.index(cell);
		
	// 	// Get geometry
	// 	auto geo = eg.geometry();

	// 	// Transformation matrix
	// 	typename EG::Geometry::JacobianInverseTransposed jac;

	// 	// loop over quadrature points
	// 	for (const auto &ip : quadratureRule(geo, intorder))
	// 	{
	// 		// evaluate test shape functions
	// 		std::vector<RangeType> phi_Pw(lfsu_Pw.size());
	// 		lfsu_Pw.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Pw);
	// 		std::vector<RangeType> psi_Pw(lfsv_Pw.size());
	// 		lfsv_Pw.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Pw);

	// 		std::vector<RangeType> phi_Sg(lfsu_Sg.size());
	// 		lfsu_Sg.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Sg);
	// 		std::vector<RangeType> psi_Sg(lfsv_Sg.size());
	// 		lfsv_Sg.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Sg);

	// 		std::vector<RangeType> phi_Sh(lfsu_Sh.size());
	// 		lfsu_Sh.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Sh);
	// 		std::vector<RangeType> psi_Sh(lfsv_Sh.size());
	// 		lfsv_Sh.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Sh);

	// 		std::vector<RangeType> phi_T(lfsu_T.size());
	// 		lfsu_T.finiteElement().localBasis().evaluateFunction(ip.position(), phi_T);
	// 		std::vector<RangeType> psi_T(lfsv_T.size());
	// 		lfsv_T.finiteElement().localBasis().evaluateFunction(ip.position(), psi_T);

	// 		std::vector<RangeType> phi_XCH4(lfsu_XCH4.size());
	// 		lfsu_XCH4.finiteElement().localBasis().evaluateFunction(ip.position(), phi_XCH4);
	// 		std::vector<RangeType> psi_XCH4(lfsv_XCH4.size());
	// 		lfsv_XCH4.finiteElement().localBasis().evaluateFunction(ip.position(), psi_XCH4);

	// 		std::vector<RangeType> phi_YH2O(lfsu_YH2O.size());
	// 		lfsu_YH2O.finiteElement().localBasis().evaluateFunction(ip.position(), phi_YH2O);
	// 		std::vector<RangeType> psi_YH2O(lfsv_YH2O.size());
	// 		lfsv_YH2O.finiteElement().localBasis().evaluateFunction(ip.position(), psi_YH2O);

	// 		std::vector<RangeType> phi_XC(lfsu_XC.size());
	// 		lfsu_XC.finiteElement().localBasis().evaluateFunction(ip.position(), phi_XC);
	// 		std::vector<RangeType> psi_XC(lfsv_XC.size());
	// 		lfsv_XC.finiteElement().localBasis().evaluateFunction(ip.position(), psi_XC);

	// 		auto ip_global = geo.global(ip.position());
	// 		auto ip_local = geo.local(ip_global);
	// 		// evaluate Pw
	// 		RF Pw = 0.0;
	// 		for (size_type i = 0; i < lfsu_Pw.size(); i++)
	// 			Pw += x(lfsu_Pw, i) * phi_Pw[i];

	// 		// evaluate Sg
	// 		RF Sg = 0.0;
	// 		for (size_type i = 0; i < lfsu_Sg.size(); i++)
	// 			Sg += x(lfsu_Sg, i) * phi_Sg[i];

	// 		// evaluate Sh
	// 		RF Sh = 0.0;
	// 		for (size_type i = 0; i < lfsu_Sh.size(); i++)
	// 			Sh += x(lfsu_Sh, i) * phi_Sh[i];

	// 		// evaluate T
	// 		RF T = 0.0;
	// 		for (size_type i = 0; i < lfsu_T.size(); i++)
	// 			T += x(lfsu_T, i) * phi_T[i];

	// 		// evaluate XCH4
    //   		RF XCH4 = 0.0;
    //   		for (size_type i = 0; i < lfsu_XCH4.size(); i++)
    //     		XCH4 += x(lfsu_XCH4, i) * phi_XCH4[i];

    //   		// evaluate YH2O
    //   		RF YH2O = 0.0;
    //   		for (size_type i = 0; i < lfsu_YH2O.size(); i++)
    //     		YH2O += x(lfsu_YH2O, i) * phi_YH2O[i];

	// 		// evaluate XC
    //   		RF XC = 0.0;
    //   		for (size_type i = 0; i < lfsu_XC.size(); i++)
    //     		XC += x(lfsu_XC, i) * phi_XC[i];

	// 		RF Sw = 1. - Sg - Sh;
	// 		// evaluate Pw
	// 		auto BrooksCParams = property.hydraulicProperty.BrooksCoreyParameters(cell, ip_local);/*BrooksCParams[0] gives Pentry in Pa*/
    //   	auto por = property.soil.SedimentPorosity(cell, ip_local);
    //   	auto Pc = property.hydraulicProperty.CapillaryPressure(cell, ip_local, Sw, Sh, por) ; /* ndim */
	// 		RF Pg = Pw + Pc;
	// 		RF Peff = (Pg * Sg + Pw * Sw) / (1. - Sh);
			
	// 		//auto por = property.soil.SedimentPorosity(cell, ip_local);
	// 		double S = XC * (property.salt.MolarMass()/property.gas.MolarMass());
    //   	auto zCH4 = property.eos.EvaluateCompressibilityFactor(T_dim, Pg_dim);
	// 		auto rho_g = property.gas.Density(T_dim, Pg_dim, zCH4);
	// 		auto rho_w = property.water.Density(T_dim, Pw * Xc_P, S);
	// 		auto rho_h = property.hydrate.Density() ;
	// 		auto rho_s = property.soil.Density() ;
	// 		auto Cv_g = property.gas.Cv(T_dim, Pg_dim, zCH4) ;
	// 		auto Cv_w = property.water.Cv(T_dim, Pw * Xc_P, S) ;
	// 		auto Cv_h = property.hydrate.Cv(T_dim, Peff * Xc_P) ;
	// 		auto Cv_s = property.soil.Cv();
	// 		auto Cv_eff = (1. - por) * rho_s * Cv_s + por * (rho_g * (1. - Sw - Sh) * Cv_g + rho_w * Sw * Cv_w + rho_h * Sh * Cv_h);

	// 		// integrate (A grad u - bu)*grad phi_i + a*u*phi_i
	// 		RF factor = ip.weight() * geo.integrationElement(ip.position());
	// 		for (size_type i = 0; i < lfsv_Sg.size(); i++)
	// 		{
	// 			for (size_type j = 0; j < lfsu_XCH4.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_Sg, i, lfsu_XCH4, j , por * rho_w * Sw * phi_XCH4[j] * psi_Sg[i] * factor);
	// 			}
	// 			for (size_type j = 0; j < lfsu_Sg.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_Sg, i, lfsu_Sg, j , por * (-rho_w * XCH4 + rho_g * (1. - YH2O)) * phi_Sg[j] * psi_Sg[i] * factor);
	// 			}
	// 			for (size_type j = 0; j < lfsu_Sh.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_Sg, i, lfsu_Sh, j , por * rho_w * XCH4 * -phi_Sh[j] * psi_Sg[i] * factor);
	// 			}
	// 			for (size_type j = 0; j < lfsu_YH2O.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_Sg, i, lfsu_YH2O, j , por * rho_g * Sg * -phi_YH2O[j] * psi_Sg[i] * factor);
	// 			}
	// 		}

	// 		for (size_type i = 0; i < lfsv_XC.size(); i++)
	// 		{
	// 			for (size_type j = 0; j < lfsu_XC.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_XC, i, lfsu_XC, j , por * rho_w * Sw * phi_XC[j] * psi_XC[i] * factor);
	// 			}
	// 			for (size_type j = 0; j < lfsu_Sg.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_XC, i, lfsu_Sg, j , por * rho_w * XC  * -phi_Sg[j] * psi_XC[i] * factor);
	// 			}
	// 			for (size_type j = 0; j < lfsu_Sh.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_XC, i, lfsu_Sh, j , por * rho_w * XC * -phi_Sh[j] * psi_XC[i] * factor);
	// 			}
	// 		}
	// 		for (size_type i = 0; i < lfsv_Pw.size(); i++)
	// 		{
	// 			for (size_type j = 0; j < lfsu_XCH4.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_Pw, i, lfsu_XCH4, j , por * rho_w * Sw * -phi_XCH4[j] * psi_Pw[i] * factor);
	// 			}
	// 			for (size_type j = 0; j < lfsu_XC.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_Pw, i, lfsu_XC, j , por * rho_w * Sw * -phi_XC[j] * psi_Pw[i] * factor);
	// 			}
	// 			for (size_type j = 0; j < lfsu_Sg.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_Pw, i, lfsu_Sg, j , por * (-rho_w * (1. - XC - XCH4) + rho_g * YH2O) * phi_Sg[j] * psi_Pw[i] * factor);
	// 			}
	// 			for (size_type j = 0; j < lfsu_Sh.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_Pw, i, lfsu_Sh, j , por * rho_w * (1. - XC - XCH4) * -phi_Sh[j] * psi_Pw[i] * factor);
	// 			}
	// 			for (size_type j = 0; j < lfsu_YH2O.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_Pw, i, lfsu_YH2O, j , por * rho_g * Sg * phi_YH2O[j] * psi_Pw[i] * factor);
	// 			}
	// 		}
	// 		for (size_type i = 0; i < lfsv_Sh.size(); i++)
	// 		{
	// 			for (size_type j = 0; j < lfsu_Sh.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_Sh, i, lfsu_Sh, j , por * rho_h * phi_Sh[j] * psi_Sh[i] * factor);
	// 			}
	// 		}
	// 		for (size_type i = 0; i < lfsv_T.size(); i++)
	// 		{
	// 			for (size_type j = 0; j < lfsu_T.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_T, i, lfsu_T, j , Cv_eff * phi_T[j] * psi_T[i] * factor);
	// 			}
	// 			for (size_type j = 0; j < lfsu_Sg.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_T, i, lfsu_Sg, j , por * (-rho_w * Cv_w + rho_g * Cv_g) * (T-T_ref) * phi_Sg[j] * psi_T[i] * factor);
	// 			}
	// 			for (size_type j = 0; j < lfsu_Sh.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_T, i, lfsu_Sh, j , por * (-rho_w * Cv_w + rho_h * Cv_h) * (T-T_ref) * phi_Sh[j] * psi_T[i] * factor);
	// 			}
	// 		}

	// 	} 	//End Quadrature Rule
	// }	// End of jacobian volume

};