#ifndef TPP_TPPNAMES_H
#define TPP_TPPNAMES_H

#include "core.hpp"

namespace tpp {

  /**
   *  \brief Base abstract class for naming TPP instances.
   */
  template<typename T>
  class NameGenerator {
    protected:
     const T &instance;
     T instance2;
    public:
     NameGenerator(T &i): instance(i) {;}
     NameGenerator(const T &i): instance2(i), instance(instance2) {;}
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

      TTCNameGenerator(const TopCoord &i): NameGenerator<TopCoord>(i) {;}

      /**
       * Add information about bonding types
       */
      TTCNameGenerator &set_btypes(boost::array<std::string,4> );

      /**
       * getName calls private subfunctions in depend on i.f
       */
      virtual std::string getName();
  };

  /** \brief Class for naming atoms.
   *
   */
  class AtomNameGenerator: NameGenerator<Atom> {
    protected:
      unsigned heavyNum = 0;
      unsigned lightNum = 1;
      bool b36Flag = false;
      const char* defaultAtomName = "X";
      std::string an2str(int);
    public:
      AtomNameGenerator(tpp::Atom &i): NameGenerator<tpp::Atom>(i) {;}
      AtomNameGenerator &setNums(unsigned h, unsigned l, bool f);

      virtual std::string getName();
  };

  /** \brief Class for naming residue according to input filename.
   */
  class ResidueNameGenerator: NameGenerator<std::string> {
    public:
      ResidueNameGenerator(const std::string &i): NameGenerator<std::string>(i) {;}

      /** \brief Guess residue name from filename
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
