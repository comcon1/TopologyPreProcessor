#ifndef TPP_TOPIO_H
#define TPP_TOPIO_H

#include "core.hpp"


namespace tpp {

  /**
   *  \brief Base abstract class for naming TPP instances.
   */
  template<typename T>
  class name_generator {
    protected:
     T instance;
    public:
     name_generator(T &i): instance(i) {;}
     virtual std::string getName() = 0;
  };

  /**
   *  \brief Class for naming force field parameters.
   */
  class ttc_name_generator: name_generator<TopCoord> {

    private:
      boost::array<std::string,4> btypes;

      std::string getBondName();
      std::string getAngleName();
      std::string getDihedralName();

    public:

      ttc_name_generator(TopCoord &i): name_generator<TopCoord>(i) {;}

      /**
       * Add information about bonding types
       */
      ttc_name_generator &set_btypes(boost::array<std::string,4> );

      /**
       * getName calls private subfunctions in depend on i.f
       */
      virtual std::string getName();
  };

  /**
   * \brief Common feature that extends the stream usage.
   *
   * Example: with name_generator instance ng:
   *  cout << ng;
   * or:
   *  s = lexical_cast<int>(ng)
   */
  template<typename T>
  std::ostream &operator << (std::ostream &out, name_generator<T> &ng) {
    out << ng.getName();
    return out;
  }

}

#endif
