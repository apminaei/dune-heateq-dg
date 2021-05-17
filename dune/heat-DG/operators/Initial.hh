/*
 * Initial.hh
 * 
 *  class of initial primary variables that depends on 
 *  user defined initial values and properties of the problem
 *   
 *  All Values from initial_conditions are dimensional and must be transfered to nondim here 
 *    
 */

#ifndef INITIAL_HH_
#define INITIAL_HH_



/** \brief A function for initial values of T
 */
template<typename GV, typename Properties,  typename RF>
class T_Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, T_Initial<GV,Properties,RF> >
{
private:
	  const GV& gv;
    const Properties& property;
    ProblemInitialConditions<GV,Properties> icvalue;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  T_Initial ( const GV& gv_, const Properties& property_ )
  : gv( gv_ ), property(property_), icvalue(gv_, property_)  {}
  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    y= icvalue.evaluate(e,xlocal);///CharacteristicValues::T_c ; //initial temperature
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};

#endif /* INITIAL_HH_ */
