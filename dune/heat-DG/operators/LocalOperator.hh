/*
 * LocalOperator.hh
 *
 *
 */



using namespace Dune::PDELab;

// struct ConvectionDiffusionDGMethod
// {
//   enum Type { NIPG, SIPG, IIPG };
// };

template <class GV, typename Params, typename U, class GFS, class FEM_T, class FEM_P>
class LocalOperator :
                      public Dune::PDELab::NumericalJacobianApplyVolume<LocalOperator<GV, Params, U, GFS, FEM_T, FEM_P>>,
                      public Dune::PDELab::NumericalJacobianVolume<LocalOperator<GV, Params, U, GFS, FEM_T, FEM_P>>,
                      public Dune::PDELab::NumericalJacobianApplySkeleton<LocalOperator<GV, Params, U, GFS, FEM_T, FEM_P>>,
                      public Dune::PDELab::NumericalJacobianSkeleton<LocalOperator<GV, Params, U, GFS, FEM_T, FEM_P>>,
                      public Dune::PDELab::NumericalJacobianApplyBoundary<LocalOperator<GV, Params, U, GFS, FEM_T, FEM_P>>,
                      public Dune::PDELab::NumericalJacobianBoundary<LocalOperator<GV, Params, U, GFS, FEM_T, FEM_P>>,
                      public Dune::PDELab::FullSkeletonPattern, // matrix entries skeleton
                      public Dune::PDELab::FullVolumePattern,
                      public Dune::PDELab::LocalOperatorDefaultFlags,
                      public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:


  const GV &gv;
  const Params&	  property;
  typedef ProblemBoundaryConditions<GV,Params> BC ;

  U *unew;
  GFS gfs;
  double *time;
  double *dt;
  double alpha_g;
  double alpha_w;
  double alpha_s;
  double alpha_T;
  double alpha_x;
  double alpha_y;
  double method_g;
  double method_w;
  double method_T;
  double method_x;
  double method_y;
  double theta_g;
  double theta_w;
  double theta_T;
  double theta_x;
  double theta_y;
  constexpr static double eps = 1.0e-6;
  constexpr static double eps_ap	= 0.;
  constexpr static double pi = atan(1.) * 4;
  unsigned int intorder;
  double Xc_conv_m;
  double Xc_conv_h;
  double Xc_source_m;
  double Xc_source_h;
  double Xc_diff_m;
  double Xc_diff_h;
  double Xc_grav;
  double Xc_K;
  double Xc_mu;
  double Xc_rho;
  double Xc_kth;
  double Xc_C;
  double Xc_P;
  double Xc_T;
  double Xc_t;
  double Xc_x;
  double T_ref;

public:
  // pattern assembly flags
	  enum { doPatternVolume	= true };
	  enum { doPatternSkeleton	= true };

	  // residual assembly flags
	  enum { doAlphaVolume  	= true };
	  enum { doAlphaSkeleton	= true };
	  enum { doAlphaBoundary	= true };


  typedef typename GV::IndexSet IndexSet;

  typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
  typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
  typedef typename U::template LocalView<LFSCache> VectorView;


  using LocalBasisType_T = typename FEM_T::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_T = Dune::PDELab::LocalBasisCache<LocalBasisType_T>;
  using LocalBasisType_P = typename FEM_P::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_P = Dune::PDELab::LocalBasisCache<LocalBasisType_P>;

  
  std::vector<Cache_T> cache_T;
  std::vector<Cache_P> cache_P;
  // constructor stores parameters
  LocalOperator(const GV &gv_, const Params&	 property_,
                     U *unew_,
                     GFS gfs_,
                     double *time_,
                     double *dt_,
                     unsigned int intorder_ = 6,
                     double method_g_ = 0.,
                     double method_w_ = 0.,
                     double method_T_ = 0.,
                     double method_x_ = 0.,
                     double method_y_ = 0.,
                     double alpha_g_ = 1., double alpha_w_ = 1., double alpha_s_ = 1., double alpha_T_ = 1., double alpha_x_ = 1., double alpha_y_ = 1.)
      : gv(gv_), property( property_ ),
        unew(unew_),
        gfs(gfs_),
        time(time_),
        dt(dt_),
        intorder(intorder_),
        method_g(method_g_), method_w(method_w_), method_T(method_T_), method_x(method_x_), method_y(method_y_),
        alpha_g(alpha_g_), alpha_w(alpha_w_), alpha_s(alpha_s_),  alpha_T(alpha_T_), alpha_x(alpha_x_), alpha_y(alpha_y_),
        cache_T(20),cache_P(20)
  {
    theta_g = method_g;
    // if (method_g == -1.)
    //   theta_g = -1.0;
    // if (method_g == 1.)
    //   theta_g = 1.0;

    theta_w = method_w;
    // if (method_w == -1.)
    //   theta_w = -1.0;
    // if (method_w == 1.)
    //   theta_w = 1.0;

    theta_T = method_T;
    // if (method_T == -1.)
    //   theta_T = -1.0;
    // if (method_T == 1.)
    //   theta_T = 1.0;

    theta_x = method_x;
    // if (method_x == ConvectionDiffusionDGMethod::SIPG)
    //   theta_x = -1.0;
    // if (method_w == ConvectionDiffusionDGMethod::NIPG)
    //   theta_x = 1.0;

    theta_y = method_y;
    // if (method_y == ConvectionDiffusionDGMethod::SIPG)
    //   theta_y = -1.0;
    // if (method_y == ConvectionDiffusionDGMethod::NIPG)
    //   theta_y = 1.0;

    Xc_conv_m = property.characteristicValue.X_convective_mass;
    Xc_conv_h = property.characteristicValue.X_convective_heat;
    Xc_source_m = property.characteristicValue.X_source_mass;
    Xc_source_h = property.characteristicValue.X_source_heat;
    Xc_diff_m = property.characteristicValue.X_diffusive_mass;
    Xc_diff_h = property.characteristicValue.X_diffusive_heat;
    Xc_grav = property.characteristicValue.X_gravity;

    Xc_K = property.characteristicValue.permeability_c;
    Xc_mu = property.characteristicValue.viscosity_c;
    Xc_rho = property.characteristicValue.density_c;
    Xc_kth = property.characteristicValue.thermalconductivity_c;
    Xc_C = property.characteristicValue.specificheat_c;
    Xc_P = property.characteristicValue.P_c;
    Xc_T = property.characteristicValue.T_c;
    Xc_t = property.characteristicValue.t_c;
    Xc_x = property.characteristicValue.x_c;
    T_ref = property.parameter.ReferenceTemperature()/Xc_T;

  }

  // volume integral depending on test and ansatz functions
  template <typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG &eg, const LFSU &lfsu, const X &x, const LFSV &lfsv, R &r) const
  {
    // subspaces
    //Temperature
    const auto &lfsv_P = lfsv.template child<Indices::PVId_P>();
    const auto &lfsu_P = lfsu.template child<Indices::PVId_P>();

    const auto &lfsv_T = lfsv.template child<Indices::PVId_T>();
    const auto &lfsu_T = lfsu.template child<Indices::PVId_T>();

    // define types
    // using RF = double;
    using RF = typename LFSU::template Child<Indices::PVId_T>::Type::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeFieldType;
    typedef typename LFSU::template Child<Indices::PVId_T>::Type::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::JacobianType JacobianType;
    using size_type = typename LFSU::template Child<Indices::PVId_T>::Type::Traits::SizeType;

    // dimensions
    const int dim = EG::Entity::dimension;
    
    const int order_t = std::max(lfsu_T.finiteElement().localBasis().order(),
                               lfsv_T.finiteElement().localBasis().order());
    const int order_p = std::max(lfsu_P.finiteElement().localBasis().order(),
                               lfsv_P.finiteElement().localBasis().order());/* If different degrees are used for different functions ? */

    // std::cout << order_p << "  " << order_x << "  " << order_s << "  " << order_t << std::endl;
    // exit(0);
    // Reference to cell
	  const auto& cell = eg.entity();
		const IndexSet &indexSet = gv.indexSet();
		int cell_number = indexSet.index(cell);

    // Get geometry
    auto geo = eg.geometry();

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T(lfsu_T.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T(lfsv_T.size());

    Dune::FieldVector<RF, dim> gradu_T(0.0);

    // Transformation matrix
    typename EG::Geometry::JacobianInverseTransposed jac;

    // loop over quadrature points
    //      auto intorder = intorderadd + quadrature_factor * order;
    for (const auto &ip : quadratureRule(geo, intorder))
    {
      // evaluate basis functions
      auto &phi_T  = cache_T[order_t].evaluateFunction(ip.position(), lfsu_T.finiteElement().localBasis());
      auto &psi_T  = cache_T[order_t].evaluateFunction(ip.position(), lfsv_T.finiteElement().localBasis());
      auto &phi_P = cache_P[order_p].evaluateFunction(ip.position(), lfsu_P.finiteElement().localBasis());
      auto &psi_P = cache_P[order_p].evaluateFunction(ip.position(), lfsv_P.finiteElement().localBasis());
      
      auto ip_global = geo.global(ip.position());
      auto ip_local = geo.local(ip_global);
      
      // evaluate T
      RF T = 0.0;
      for (int i = 0; i < lfsu_T.size(); i++)
        T += x(lfsu_T, i) * phi_T[i];

      // evaluate T
      RF P = 0.0;
      for (int i = 0; i < lfsu_P.size(); i++)
        P += x(lfsu_P, i) * phi_P[i];

      
      auto T_dim = T * Xc_T;

      // evaluate gradient of basis functions
      auto &js_T = cache_T[order_t].evaluateJacobian(ip.position(), lfsu_T.finiteElement().localBasis());
      auto &js_v_T = cache_T[order_t].evaluateJacobian(ip.position(), lfsv_T.finiteElement().localBasis());
      
      // transform gradients of shape functions to real element
      jac = geo.jacobianInverseTransposed(ip.position());

      for (int i = 0; i < lfsu_T.size(); i++)
        jac.mv(js_T[i][0], gradphi_T[i]);
      for (int i = 0; i < lfsv_T.size(); i++)
        jac.mv(js_v_T[i][0], gradpsi_T[i]);

      // compute gradient of T
      gradu_T = 0.0;
      for (int i = 0; i < lfsu_T.size(); i++)
        gradu_T.axpy(x(lfsu_T, i), gradphi_T[i]);

      auto diffusiveflux_Heat =  gradu_T;


      
      auto f =   2*std::sin(3.14 *  (ip_global[0] *  ((*time)+(*dt)))); // (ip.position()[0] * ip.position()[1])* ((*time)+(*dt))*
      
      if (ip_global[0] > 0.5)
        diffusiveflux_Heat *= 0.5;
      // f += -2*((*time))*((*time))* ip.position()[1] * ip.position()[1]
      //   *3.14*3.14 *std::sin(3.14 * (ip.position()[0] * ip.position()[1]) * ((*time)));
      // f += -2*(*time)*(*time)*ip.position()[0]*ip.position()[0] 
      //   *3.14*3.14 *std::sin(3.14 * (ip.position()[0] * ip.position()[1]) * (*time));
      auto factor = ip.weight() * geo.integrationElement(ip.position());
      
      for (int i = 0; i < lfsv_T.size(); i++)
      {
        r.accumulate(lfsv_T, i, ( diffusiveflux_Heat * gradpsi_T[i] - P * psi_T[i] - f * psi_T[i]) * factor);//- 25.e0  * psi_T[i]
      }
      // if (T <= 0)  
      for (int i = 0; i < lfsv_P.size(); i++)
      {
        r.accumulate(lfsv_P, i, (  (P-std::max(0., (P-(T))))* psi_P[i]) * factor);// 
      }
      
     
    } //End Quadrature Rule
    
    
  }  // End of alpha_volume



  // skeleton integral depending on test and ansatz functions
  // each face is only visited ONCE!
  template <typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_skeleton(const IG &ig,
                      const LFSU &lfsu_s, const X &x_s, const LFSV &lfsv_s,
                      const LFSU &lfsu_n, const X &x_n, const LFSV &lfsv_n,
                      R &r_s, R &r_n) const
  {
    // subspaces
     //Temperature
    const auto &lfsv_T_s = lfsv_s.template child<Indices::PVId_T>();
    const auto &lfsu_T_s = lfsu_s.template child<Indices::PVId_T>();
    const auto &lfsv_T_n = lfsv_n.template child<Indices::PVId_T>();
    const auto &lfsu_T_n = lfsu_n.template child<Indices::PVId_T>();

     //Temperature
    const auto &lfsv_P_s = lfsv_s.template child<Indices::PVId_P>();
    const auto &lfsu_P_s = lfsu_s.template child<Indices::PVId_P>();
    const auto &lfsv_P_n = lfsv_n.template child<Indices::PVId_P>();
    const auto &lfsu_P_n = lfsu_n.template child<Indices::PVId_P>();


  
    // define types
    using RF = typename LFSU::template Child<Indices::PVId_T>::Type::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeFieldType;
    using size_type = typename LFSU::template Child<Indices::PVId_T>::Type::Traits::SizeType;

    // using RF = double;

    // dimensions
    const int dim= IG::Entity::dimension;
    const int dimension = GV::dimension;
    const int order_t = std::max(lfsu_T_s.finiteElement().localBasis().order(),
                               lfsv_T_s.finiteElement().localBasis().order());
    const int order_p = std::max(lfsu_P_s.finiteElement().localBasis().order(),
                                lfsv_P_s.finiteElement().localBasis().order());/* If different degrees are used for different functions ? */
    

    // References to inside and outside cells
    const auto &cell_inside = ig.inside();
    const auto &cell_outside = ig.outside();

    // Get geometries
    auto geo = ig.geometry();
    //const auto dimension = geo.mydimension;
    auto geo_inside = cell_inside.geometry();
    auto geo_outside = cell_outside.geometry();

    // Get geometry of intersection in local coordinates of cell_inside and cell_outside
    auto geo_in_inside = ig.geometryInInside();
    auto geo_in_outside = ig.geometryInOutside();
    // cell geometries
    auto ref_el_inside 	= referenceElement(geo_inside);
	  auto ref_el_outside = referenceElement(geo_outside);
	  auto inside_cell_center_local 	= ref_el_inside.position(0,0);
	  auto outside_cell_center_local 	= ref_el_outside.position(0,0);
    // face diameter; this should be revised for anisotropic meshes?
    auto h_F = std::min(geo_inside.volume(), geo_outside.volume()) / geo.volume(); // Houston!

    // compute weights
    RF omega_s;
    RF omega_n;
    RF harmonic_average(0.0);
    omega_s = 0.5;
    omega_n = 0.5;
    harmonic_average = 1.0;

        //  if (true)
        //    {
        //      RF delta_s = (An_F_s*n_F);
        //      RF delta_n = (An_F_n*n_F);
        //      omega_s = delta_n/(delta_s+delta_n+1e-20);
        //      omega_n = delta_s/(delta_s+delta_n+1e-20);
        //      harmonic_average = 2.0*delta_s*delta_n/(delta_s+delta_n+1e-20);
        //    }
        //  else
        //    {
        //      omega_s = omega_n = 0.5;
        //      harmonic_average = 1.0;
        //    }

    // get polynomial degree
    auto order_i = lfsu_T_s.finiteElement().localBasis().order();
    auto order_o = lfsv_T_n.finiteElement().localBasis().order();
    auto degree = std::max(order_i, order_o);

    // penalty factor
    auto penalty_factor_g = (alpha_g / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_w = (alpha_w / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_s = (alpha_s / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_T = (alpha_T / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_x = (alpha_x / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_y = (alpha_y / h_F) * harmonic_average * degree * (degree + dim - 1);

    // std::cout << geo_inside.volume() << "  " <<  geo_outside.volume()<< "  " << geo.volume() << "  "<< penalty_factor_g << std::endl;
    // exit(0);

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T_s(lfsu_T_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T_s(lfsv_T_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T_n(lfsu_T_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T_n(lfsv_T_n.size());
    Dune::FieldVector<RF, dim> gradu_T_s(0.0);
    Dune::FieldVector<RF, dim> gradu_T_n(0.0);
    
    // Transformation matrix
    typename IG::Entity::Geometry::JacobianInverseTransposed jac;

    // loop over quadrature points
    auto intorder1 = intorder;
    for (const auto &ip : quadratureRule(geo, intorder1))
    {
      // exact normal
      auto n_F_local = ig.unitOuterNormal(ip.position());

      // position of quadrature point in local coordinates of elements
      auto iplocal_s = geo_in_inside.global(ip.position());
      auto iplocal_n = geo_in_outside.global(ip.position());

      auto ip_global_s = geo_inside.global(iplocal_s);
      auto ip_global_n = geo_outside.global(iplocal_n);

      // evaluate basis functions
      auto &phi_T_s = cache_T[order_t].evaluateFunction(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &psi_T_s = cache_T[order_t].evaluateFunction(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &phi_T_n = cache_T[order_t].evaluateFunction(iplocal_n, lfsu_T_n.finiteElement().localBasis());
      auto &psi_T_n = cache_T[order_t].evaluateFunction(iplocal_n, lfsv_T_n.finiteElement().localBasis());
      // evaluate basis functions
      auto &phi_P_s = cache_P[order_p].evaluateFunction(iplocal_s, lfsu_P_s.finiteElement().localBasis());
      auto &psi_P_s = cache_P[order_p].evaluateFunction(iplocal_s, lfsv_P_s.finiteElement().localBasis());
      auto &phi_P_n = cache_P[order_p].evaluateFunction(iplocal_n, lfsu_P_n.finiteElement().localBasis());
      auto &psi_P_n = cache_P[order_p].evaluateFunction(iplocal_n, lfsv_P_n.finiteElement().localBasis());
      
      // evaluate T
      RF T_s = 0.0;
      for (int i = 0; i < lfsu_T_s.size(); i++)
        T_s += x_s(lfsu_T_s, i) * phi_T_s[i];
      RF T_n = 0.0;
      for (int i = 0; i < lfsu_T_n.size(); i++)
        T_n += x_n(lfsu_T_n, i) * phi_T_n[i];

      // evaluate T
      RF P_s = 0.0;
      for (int i = 0; i < lfsu_P_s.size(); i++)
        P_s += x_s(lfsu_P_s, i) * phi_P_s[i];
      RF P_n = 0.0;
      for (int i = 0; i < lfsu_P_n.size(); i++)
        P_n += x_n(lfsu_P_n, i) * phi_P_n[i];

      
      auto T_s_dim = T_s * Xc_T;
      auto T_n_dim = T_n * Xc_T;



      // evaluate gradient of basis functions
      auto &js_T_s = cache_T[order_t].evaluateJacobian(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &js_v_T_s = cache_T[order_t].evaluateJacobian(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &js_T_n = cache_T[order_t].evaluateJacobian(iplocal_n, lfsu_T_n.finiteElement().localBasis());
      auto &js_v_T_n = cache_T[order_t].evaluateJacobian(iplocal_n, lfsv_T_n.finiteElement().localBasis());
      
      // transform gradients of shape functions to real element
      jac = geo_inside.jacobianInverseTransposed(iplocal_s);
      for (int i = 0; i < lfsu_T_s.size(); i++)
        jac.mv(js_T_s[i][0], gradphi_T_s[i]);
      for (int i = 0; i < lfsv_T_s.size(); i++)
        jac.mv(js_v_T_s[i][0], gradpsi_T_s[i]);

      
      jac = geo_outside.jacobianInverseTransposed(iplocal_n);
      for (int i = 0; i < lfsu_T_n.size(); i++)
        jac.mv(js_T_n[i][0], gradphi_T_n[i]);
      for (int i = 0; i < lfsv_T_n.size(); i++)
        jac.mv(js_v_T_n[i][0], gradpsi_T_n[i]);

      // compute gradient of T
      gradu_T_s = 0.0;
      for (int i = 0; i < lfsu_T_s.size(); i++)
        gradu_T_s.axpy(x_s(lfsu_T_s, i), gradphi_T_s[i]);
      gradu_T_n = 0.0;
      for (int i = 0; i < lfsu_T_n.size(); i++)
        gradu_T_n.axpy(x_n(lfsu_T_n, i), gradphi_T_n[i]);

      
      double normalflux_T = (omega_s * gradu_T_s + omega_n * gradu_T_n) * n_F_local;
      
      double normalflux_T_diff = ( gradu_T_s - gradu_T_n) * n_F_local;
     
      RF omegaup_T_s, omegaup_T_n;
      if (normalflux_T>=0.0)
      {
        omegaup_T_s = 0.5;
        omegaup_T_n = 0.5;
      }
      else
      {
        omegaup_T_s = 0.5;
        omegaup_T_n = 0.5;
      }

      
      auto diffusiveflux_Heat_s = omegaup_T_s * gradu_T_s;
      // *******************   //
      
      if (ip_global_s[0] > 0.5)
        diffusiveflux_Heat_s *= 0.2;
      auto diffusiveflux_Heat_n = omegaup_T_n * gradu_T_n;

      if (ip_global_n[0] > 0.5)
        diffusiveflux_Heat_n *= 0.2;


      
      auto diffusiveflux_Heat = -(diffusiveflux_Heat_s + diffusiveflux_Heat_n) * n_F_local;

      /*ACCCUMULATE RESIDUALS*/
			double tmp=0.;
      
      auto factor = ip.weight() * geo.integrationElement(ip.position());
      // ENERGY balance
      tmp =  diffusiveflux_Heat;
      double term_nipg_T = theta_T * (T_s - T_n);
      double term_penalty_T = penalty_factor_T * (T_s - T_n) ;
      // diffusion term
      for (int i = 0; i < lfsv_T_s.size(); i++)
      {
        r_s.accumulate(lfsv_T_s, i, tmp * psi_T_s[i] * factor);
      }
      for (int i = 0; i < lfsv_T_n.size(); i++)
      {
        r_n.accumulate(lfsv_T_n, i, tmp * -psi_T_n[i] * factor);
      }
      // (non-)symmetric IP term
      for (int i = 0; i < lfsv_T_s.size(); i++)
      {
        r_s.accumulate(lfsv_T_s, i, -omegaup_T_s *  term_nipg_T * n_F_local * gradpsi_T_s[i] * factor);
      }
      for (int i = 0; i < lfsv_T_n.size(); i++)
      {
        r_n.accumulate(lfsv_T_n, i, -omegaup_T_n *  term_nipg_T * n_F_local * gradpsi_T_n[i] * factor);
      }
      // standard IP term integral
      for (int i = 0; i < lfsv_T_s.size(); i++)
      {
        r_s.accumulate(lfsv_T_s, i, term_penalty_T * psi_T_s[i] * factor);
      }
      for (int i = 0; i < lfsv_T_n.size(); i++)
      {
        r_n.accumulate(lfsv_T_n, i, term_penalty_T * -psi_T_n[i] * factor);
      }

      // for (int i = 0; i < lfsv_P_s.size(); i++)
      // {
      //   r_s.accumulate(lfsv_P_s, i, term_penalty_T * psi_P_s[i] * factor);
      // }
      // for (int i = 0; i < lfsv_P_n.size(); i++)
      // {
      //   r_n.accumulate(lfsv_P_n, i, term_penalty_T * -psi_P_n[i] * factor);
      // }
     
    } //End Quadrature Rule
  }//End of alpha_skeleton





  template <typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary(const IG &ig,
                      const LFSU &lfsu,
                      const X &x,
                      const LFSV &lfsv,
                      R &r) const
  {
    // subspaces
    //Temperature
     //Temperature
    const auto &lfsv_T = lfsv.template child<Indices::PVId_T>();
    const auto &lfsu_T = lfsu.template child<Indices::PVId_T>();
    // define types
    

    using RF = double;
    // dimensions
    const int dimension = GV::dimension;
    const int dim = IG::Entity::dimension;
    const int order_t = std::max(lfsu_T.finiteElement().localBasis().order(),
                               lfsv_T.finiteElement().localBasis().order());

    // References to inside and outside cells
    const auto &cell_inside = ig.inside();

    // Get geometries
    auto geo = ig.geometry();
    //const auto dimension = geo.mydimension;
    auto geo_inside = cell_inside.geometry();

    // Get geometry of intersection in local coordinates of cell_inside and cell_outside
    auto geo_in_inside = ig.geometryInInside();

    // cell geometries
    auto ref_el_inside 	= referenceElement(geo_inside);
    auto inside_cell_center_local 	= ref_el_inside.position(0,0);
    auto inside_cell_center_global 	= geo_inside.center();

    // face geometry
    auto ref_el = referenceElement(geo);
    auto face_center_local = ref_el.position(0,0);
    auto face_center_global = geo.center();

    // face diameter; this should be revised for anisotropic meshes?
    auto h_F = geo_inside.volume() / geo.volume(); // Houston!

    // compute weights
    RF omega_s;
    RF omega_n;
    RF harmonic_average(0.0);
    harmonic_average = 1.0;

    // get polynomial degree
    auto degree = lfsv_T.finiteElement().localBasis().order();

    // penalty factor
    auto penalty_factor_g = (alpha_g / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_w = (alpha_w / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_s = (alpha_s / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_T = (alpha_T / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_x = (alpha_x / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_y = (alpha_y / h_F) * harmonic_average * degree * (degree + dim - 1);

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T_s(lfsu_T.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T_s(lfsv_T.size());
    Dune::FieldVector<RF, dim> gradu_T_s(0.0);
    
    // Transformation matrix
    typename IG::Entity::Geometry::JacobianInverseTransposed jac;

    auto intorder1 = intorder;
    // loop over quadrature points
    for (const auto &ip : quadratureRule(geo, intorder1))
    {
      // integration factor
      auto factor = ip.weight() * geo.integrationElement(ip.position());

      // exact normal
      auto n_F_local = ig.unitOuterNormal(ip.position());

      // position of quadrature point in local coordinates of elements
      auto iplocal_s = geo_in_inside.global(ip.position());
      auto ip_global_s = geo_inside.global(iplocal_s);

      BC bc( gv,property ) ;

      // evaluate boundary condition types for {Pw,Sg} or {Fw,Fg}
			auto bctype = bc.type(ig, ip.position(), (*time)*Xc_t, (*dt)*Xc_t) ;
      auto veltype = bc.velType(ig, ip.position(), (*time)*Xc_t, (*dt)*Xc_t) ;
			// evaluate boundary condition values for {Pw,Sg} or {Fw,Fg}
			auto bcvalue = bc.value(ig, ip.position(), (*time)*Xc_t, (*dt) *Xc_t) ;
      auto velvalue = bc.velValue(ig, ip.position(), (*time)*Xc_t, (*dt) *Xc_t) ;

      // evaluate basis functions at local quadrature points
      auto &psi_T_s = cache_T[order_t].evaluateFunction(iplocal_s, lfsv_T.finiteElement().localBasis());
      auto &phi_T_s = cache_T[order_t].evaluateFunction(iplocal_s, lfsu_T.finiteElement().localBasis());
      
      // evaluate T
      RF T_s = 0.0;
      for (int i = 0; i < lfsu_T.size(); i++)
        T_s += x(lfsu_T, i) * phi_T_s[i];

      RF T_n = T_s;
      if (bctype == Indices::BCId_dirichlet)
      {
        T_n = bcvalue ;
      }


      omega_s = 0.5;
      omega_n = 0.5;

      // evaluate gradient of basis functions
      auto &js_T_s = cache_T[order_t].evaluateJacobian(iplocal_s, lfsu_T.finiteElement().localBasis());
      auto &js_v_T_s = cache_T[order_t].evaluateJacobian(iplocal_s, lfsv_T.finiteElement().localBasis());
      
      // transform gradients of shape functions to real element
      jac = geo_inside.jacobianInverseTransposed(iplocal_s);
      for (int i = 0; i < lfsu_T.size(); i++)
        jac.mv(js_T_s[i][0], gradphi_T_s[i]);
      for (int i = 0; i < lfsv_T.size(); i++)
        jac.mv(js_v_T_s[i][0], gradpsi_T_s[i]);

      
      // compute gradient of T
      gradu_T_s = 0.0;
      for (int i = 0; i < lfsu_T.size(); i++)
        gradu_T_s.axpy(x(lfsu_T, i), gradphi_T_s[i]);

      
      // evaluate normal flux of T
      RF grad_T_s = gradu_T_s * n_F_local;
      RF grad_T_n = grad_T_s;
      if (veltype == Indices::BCId_neumann)
      {
        grad_T_n = velvalue;
      }

      
      
      double normalflux_T = ( omega_s * grad_T_s +omega_n * grad_T_n);//
     
      
      RF omegaup_T_s, omegaup_T_n;
      if (normalflux_T>=0.0)
      {
        omegaup_T_s = 0.5;
        omegaup_T_n = 0.5;
      }
      else
      {
        omegaup_T_s = 0.5;
        omegaup_T_n = 0.5;
      }


      //   fluxes and diff. flux
      
      auto diffusiveflux_Heat_s = omegaup_T_s * grad_T_s; // k_eff will be harmonic_average of k_eff_s and k_eff_n

      //   *******************   //
      
      auto diffusiveflux_Heat_n = omegaup_T_n * grad_T_n; // k_eff will be harmonic_average of k_eff_s and k_eff_n


      auto diffusiveflux_Heat = (diffusiveflux_Heat_s + diffusiveflux_Heat_n);
      if (veltype == Indices::BCId_neumann){
        diffusiveflux_Heat = ( grad_T_n);
      }


      
      //  ACCCUMULATE RESIDUALS  //
			RF tmp=0.;
      // ENERGY balance
      tmp =  -diffusiveflux_Heat;
      double term_nipg_T = method_T * (T_s - T_n);
      double term_penalty_T = penalty_factor_T * (T_s - T_n) ;

      for (int i = 0; i < lfsv_T.size(); i++)
      {
        r.accumulate(lfsv_T, i, tmp * psi_T_s[i] * factor);
      }

      // (non-)symmetric IP term
      for (int i = 0; i < lfsv_T.size(); i++)
      {
        r.accumulate(lfsv_T, i,    term_nipg_T * n_F_local * gradpsi_T_s[i] * factor); // in the run testAveragingXC-T there is no upwinding for sym terms
      }

      // standard IP term integral
      for (int i = 0; i < lfsv_T.size(); i++)
      {
        r.accumulate(lfsv_T, i, term_penalty_T * psi_T_s[i] * factor);
      }

      
    } // end of quadrature rule
  } // end of alpha_boundary



};
