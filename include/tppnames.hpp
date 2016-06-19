#ifndef TPP_TOPIO_H
#define TPP_TOPIO_H

#include "global.hpp"


namespace tpp {

  /*
   * Base abstract class for naming TPP instances.
   */
  template<typename T>
  class name_generator {
    protected:
     T instance;
    public:
     name_generator(T &i): instance(i) {;}
     virtual string getName() = 0;
  };

  /*
   * Class for naming force field parameters.
   */
  class ttc_name_generator: name_generator<t_top_coord> {

    private:
      array<string,4> btypes;

      string getBondName();
      string getAngleName();
      string getDihedralName();

    public:

      ttc_name_generator(t_top_coord &i): name_generator<t_top_coord>(i) {;}

      /*
       * Add information about bonding types
       */
      ttc_name_generator &set_btypes(array<string,4> );

      /*
       * getName calls private subfunctions in depend on i.f
       */
      virtual string getName();
  };

  /*
   * Common feature to be able the following with name_generator instance ng:
   *  cout << ng;
   * or the following:
   *  s = lexical_cast<int>(ng)
   */
  template<typename T>
  ostream &operator << (ostream &out, name_generator<T> &ng) {
    out << ng.getName();
    return out;
  }

}

#endif
