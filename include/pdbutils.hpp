#include "global.hpp"
#include "exceptions.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>
#include <set>


namespace tpp {


  class BondMatrix {
    private:
      std::set<int> curset;
      boost::numeric::ublas::matrix<bool> mtx;
      void rec_check(int i);
      bool check();
    public:
      BondMatrix(const OBMol&);
      ~BondMatrix() {;}
      const boost::numeric::ublas::matrix<bool> &operator *() {
        return mtx;
      }
      void print (std::ostream &os);
  };

  extern std::pair<boost::numeric::ublas::matrix<int>, 
         boost::numeric::ublas::matrix<int> > top_neig
   (OBMol &, boost::numeric::ublas::matrix<int>, int, int, int, int);
  extern boost::numeric::ublas::vector<int> generate_long_tail(OBMol &);
  extern boost::numeric::ublas::vector<int> generate_long_tail1(OBMol &);
  extern t_atom_array mol_renum(OBMol &, t_atom_array &, boost::numeric::ublas::vector<int>);
  extern t_atom_array mol_renum1(OBMol &, t_atom_array &, boost::numeric::ublas::vector<int>);
}

