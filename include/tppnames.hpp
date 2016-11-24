#ifndef TPP_TOPIO_H
#define TPP_TOPIO_H

#include "core.hpp"


namespace tpp {

  /**
   *  \brief Base abstract class for naming TPP instances.
   */
  template<typename T>
  class NameGenerator {
    protected:
     T instance;
    public:
     NameGenerator(T &i): instance(i) {;}
     virtual ~NameGenerator() {};
     virtual std::string getName() = 0;
  };

  /**
   *  \brief Class for naming force field parameters.
   */
  class TTCNameGenerator: NameGenerator<TopCoord> {

    private:
      boost::array<std::string,4> btypes;

      std::string getBondName();
      std::string getAngleName();
      std::string getDihedralName();

    public:

      TTCNameGenerator(TopCoord &i): NameGenerator<TopCoord>(i) {;}

      /**
       * Add information about bonding types
       */
      TTCNameGenerator &set_btypes(boost::array<std::string,4> );

      /**
       * getName calls private subfunctions in depend on i.f
       */
      virtual std::string getName();
  };

  /**
   * \brief Common feature that extends the stream usage.
   *
   * Example: with NameGenerator instance ng:
   *  cout << ng;
   * or:
   *  s = lexical_cast<int>(ng)
   */
  template<typename T>
  std::ostream &operator << (std::ostream &out, NameGenerator<T> &ng) {
    out << ng.getName();
    return out;
  }

}

#endif
