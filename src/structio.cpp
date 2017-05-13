#include "structio.hpp"
#include "exceptions.hpp"

#include "strutil.hpp"
#include "logger.hpp"

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/algorithm/string.hpp>

#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>

#include <cctype>

using std::string;
using std::to_string;
using std::pair;
using std::ios;
using std::fstream;

using boost::numeric_cast;
using boost::lexical_cast;
using boost::format;

using OpenBabel::OBConversion;
using OpenBabel::OBMol;
using OpenBabel::OBMolAtomIter;
using OpenBabel::OBAtomAtomIter;

namespace tpp {


  StructureIO::StructureIO(bool ign, bool rtp) : ignoreIndexFlag(ign), rtpoutput_file(rtp){
    // TODO: some interface stuff
    ;
  }


  void StructureIO::loadFromFile(Topology &tp, InputFormat ifm,
                                 const char *fname)  {
    // test if file exists
    TPPI << format("Trying to read structure from '%s' ..") % fname;
    fstream inf(fname, ios::in);
    if (!inf.is_open()) {
      TPPE << "Fail to open file for read.";
      Exception e("Can't open specified file for read.");
      e.add("procname", "tpp::load_struct")
       .add("error", "invalid filename")
       .add("filename", fname);
      throw e;
    }
    // loading from stream

    loadFromStream(tp, ifm, &inf);
    inf.close();
    TPPI << format("Structure file '%s' has been processed.") % fname;
  } // end StructureIO::loadFromFile


  void StructureIO::loadFromStream(Topology &tp, InputFormat ifm, std::istream *inf) {
    try {
      // test if tp variable contains structure
      if (!tp.atoms.empty()) {
        TPPI << "Loading will replacing your current structure.";
        tp.atoms.clear();
      }

      // Reading molecule in OpenBabel format
      // -----------------------------------
      OBConversion conv(inf);
      switch (ifm) {
        case TPP_IF_PDB: conv.SetInFormat("PDB"); break;
        case TPP_IF_GRO: conv.SetInFormat("GRO"); break;
        case TPP_IF_G96: conv.SetInFormat("G96"); break;
        case TPP_IF_GAMOPT: ;
        case TPP_IF_GAMHESS: ;
        case TPP_IF_GAMSP: conv.SetInFormat("GAMOUT"); break;
        default:
          Exception e("Invalid input format");
          e.add("procname", "tpp::loadFromStream");
          e.add("error", "OpenBabel: parsing error");
          e.add("format", to_string(ifm));
          throw e;
      }
      // reading from file with OpenBabel function
      TPPI << "Reading molecule by OpenBabel..";
      OBMol mol;
      if ( (!conv.Read(&mol)) || (!mol.NumAtoms()) ) {
        Exception e("Can't read file format.");
        e.add("procname", "tpp::load_struct");
        e.add("error", "OpenBabel: parsing error");
        throw e;
      }
      TPPI << "                               ..DONE.";

      tp.mol = mol;

      TPPI << " ----> Reading molecule by internal procedure..";
      // Reading molecule in AtomArray format
      // ---------------------------------------
      switch(ifm) {
        case TPP_IF_PDB: {
                  char *s0  = new char[300];
                  string curString, field, _aName, _rName, _qName;
                  int  fieldCounter, strc, incrementalIndex = 0;
                  std::pair<AtomArray::iterator, bool> at_it;
                  Atom cur0;
                  float __x, __y, __z;
                  strc = 0;
                  inf->clear();
                  inf->seekg(0);
                  while (! inf->eof() ) {
                    inf->getline(s0, 300);
                    curString = string(s0);
                    strc++;

                    fieldCounter = 0;
                    _qName = "un";

                    // ATOM | HETATM field
                    field = curString.substr(0,6);
                    boost::trim( field );
                    if ( (field != "ATOM") && (field != "HETATM") )
                      continue;
                    fieldCounter++;
                    try {
                      // Atom index
                      field = curString.substr(6,5);
                      boost::trim(field);
                      cur0.oldindex = lexical_cast<unsigned>(field);
                      cur0.index = ignoreIndexFlag ? 65535 : cur0.oldindex;
                      fieldCounter += (field.size() > 0);
                      // Atom name
                      field = curString.substr(11,5);
                      boost::trim( field );
                      fieldCounter += (field.size() > 0);
                      _aName = field;
                      // Residue name
                      field = curString.substr(16,4);
                      boost::trim( field );
                      fieldCounter += (field.size() > 0);
                      _rName = field;
                      // pass chain letter 21 -> 23
                      // Molecule number
                      field = curString.substr(23,4);
                      boost::trim(field);
                      cur0.mol_id = numeric_cast<unsigned char>( lexical_cast<unsigned>(field) );
                      fieldCounter += (field.size() > 0);
                      // Coordinates
                      field = curString.substr(30,8);
                      boost::trim(field);
                      __x = lexical_cast<float>(field);
                      fieldCounter += (field.size() > 0);
                      field = curString.substr(38,8);
                      boost::trim(field);
                      __y = lexical_cast<float>(field);
                      fieldCounter += (field.size() > 0);
                      field = curString.substr(46,8);
                      boost::trim(field);
                      __z = lexical_cast<float>(field);
                      fieldCounter += (field.size() > 0);
                    } catch (const boost::bad_lexical_cast& e) {
                      TPPE << format("** Caught bad lexical cast with error: %s") % e.what();
                      TPPE << format("** Last field read: [%1$s]") % field;
                      TPPE << curString;
                      Exception tppe("Failed to extract numbers from the ATOM string.");
                      tppe.add("classname", "StructIO");
                      tppe.add("procname", "loadFromStream");
                      tppe.add("error", "PDB parsing error");
                      tppe.add("line", lexical_cast<string>(strc).c_str());
                      throw tppe;
                    } catch (const boost::bad_numeric_cast &e) {
                      TPPE << format("** Caught bad numeric cast with error: ") % e.what();
                      TPPE << format("** Last field read: [%1$s]") % field;
                      TPPE << curString;
                      Exception tppe("Some of ATOM indexes is out of bonds.");
                      tppe.add("classname", "StructIO");
                      tppe.add("procname", "loadFromStream");
                      tppe.add("error", "PDB numeric conversion error");
                      tppe.add("line", lexical_cast<string>(strc).c_str());
                      throw tppe;
                    }
                    // Chemical atom
                    if (curString.size() > 76) {
                      field = curString.substr(76,2);
                      boost::trim( field );
                      fieldCounter += (field.size() > 0);
                      _qName = field;
                    }

                    // internal check
                    if ( fieldCounter < 8 ) {
                      Exception e("Invalid ATOM string in your PDB file.");
                      e.add("classname", "StructIO");
                      e.add("procname", "loadFromStream");
                      e.add("error", "PDB parsing error");
                      e.add("line", lexical_cast<string>(strc).c_str());
                      throw e;
                    }

                    cur0.old_aname = string(_aName);
                    cur0.atom_name = string(_aName);
                    if (rtpoutput_file) {
                      // need to replace 1H2 to H21
                      if (isdigit(_aName[0]))
                        cur0.atom_name = cur0.old_aname.substr(1) + cur0.old_aname.substr(0,1);
                    }
                    cur0.res_name  = string(_rName);
                    cur0.qmname    = string(_qName);
                    cur0.coord(0)  = numeric_cast<double>(__x);
                    cur0.coord(1)  = numeric_cast<double>(__y);
                    cur0.coord(2)  = numeric_cast<double>(__z);
                    cur0.comment   = string("QMname: ") + _qName;

                    TPPD << format("  %4s - %4s: %3d(%3d) [%8.3f,%8.3f,%8.3f] %s") % cur0.res_name % cur0.atom_name
                          % cur0.oldindex % cur0.index % cur0.coord(0) % cur0.coord(1) % cur0.coord(2) % cur0.comment;
                    at_it = tp.atoms.insert( cur0 );

                    if (! at_it.second) {
                      TPPE << "Failed to insert just readed atom.";
                      Exception e("Repeat index or something else.");
                      e.add("classname", "StructIO");
                      e.add("procname", "loadFromStream");
                      e.add("error", "PDB parsing error");
                      e.add("line", lexical_cast<string>(strc).c_str());
                      throw e;
                    }

                    if (ignoreIndexFlag) {
                      incrementalIndex++;
                      cur0.index = incrementalIndex;
                      tp.atoms.replace(at_it.first, cur0);
                    }

                  } // end-while
                  delete[] s0;
                  /* FINISH PARSING PDB */
                } // end-case PDB
                         break;
        case TPP_IF_G96:   molToAtoms(tp); // #TODO 1
                           break;
        case TPP_IF_GRO:   molToAtoms(tp); // #TODO 2
                           break;
        case TPP_IF_GAMOPT: ;
        case TPP_IF_GAMHESS: ;
        case TPP_IF_GAMSP:   molToAtoms(tp); // #TODO 3
                             break;
        default:
          Exception e("Invalid input format");
          e.add("classname", "StructIO");
          e.add("procname", "loadFromStream");
          e.add("error", "PDB parsing error");
          throw e;
      };

      // setting up residue name
      if ( strutil::trim(tp.atoms.begin()->res_name) != "") {
        tp.res_name = tp.atoms.begin()->res_name;
        // checking res_name correctness
        for (int i=0; i<tp.res_name.size(); ++i) {
          if (!isalnum(tp.res_name[i])) {
            TPPD << "Incorrect residue name in PDB file, using 'RES' instead.";
            tp.res_name = "RES";
            for (AtomArray::iterator it = tp.atoms.begin(); it != tp.atoms.end(); ++it) {
              Atom at(*it);
              at.res_name = "RES";
              tp.atoms.replace(it, at);
            }
          }
        }
      } else {
        tp.res_name = "VMO";
      }

      if ( tp.atoms.empty() ) {
        Exception e("No atoms in structure file or bad format.");
        e.add("classname", "StructIO");
        e.add("procname", "loadFromStream");
        e.add("error", "no atoms");
        throw e;
      } else {
        TPPD << "  ========";
        TPPI << format("  Successfully readed %d atoms!") % tp.atoms.size();
      }
//              ----> Reading molecule by internal procedure..
      TPPI << "                                               ..DONE. <----";
    } catch(const Exception &e) {
      //e.fix_log();
      throw e;
    }

  } // end loadFromStream

  void StructureIO::saveToFile(Topology &tp, OutputFormat ofm, const char *fname) {
    try {
    // test if file exists
    TPPI << format("Starting to write structure into '%s'") % fname;
    fstream out(fname, ios::out);
    if (!out.is_open()) {
      TPPE << "Fail to open file for write.";
      Exception e("Can't open specified file for write.");
      e.add("classname", "StructIO");
      e.add("procname", "saveToFile");
      e.add("error", "invalid filename");
      e.add("filename", fname);
      throw e;
    }
    switch (ofm) {
      case TPP_OF_PDB:
        out << format("TITLE  Written by TPP: %1$s (topoplogy for residue %2$-4s)\n") % tp.name % tp.res_name;
        for (AtomArray::iterator it = tp.atoms.begin(); it != tp.atoms.end(); ++it) {
          out << format("%1$-6s%2$5d %3$4s%4$4s %5$c%6$4d    %7$8.3f%8$8.3f%9$8.3f  0.00  0.00          %10$2s\n")
            % "ATOM" % (int)it->index % it->atom_name % tp.res_name % 'A'
            % (int)it->mol_id % it->coord(0) % it->coord(1) % it->coord(2) % it->qmname;
        }
        out << "END\n";
        break;
      case TPP_OF_GRO: /* #TODO 4*/
  /*       os << format("%1$5d%2$-3s    %3$-3s%4$5d%5$8.3f%6$8.3f%7$8.3f\n");*/
        break;
      case TPP_OF_G96: /* #TODO 5 */
  /*       os << format("%1$5d %2$-4s  %3$-4s %4$7d%5$15.9f%6$15.9f%7$15.9f\n")*/
        break;
      case TPP_OF_GAMIN:
        out << format("\n\
  $DATA\n\
  %1$s\n\
  C1\n") % tp.name;
        for (AtomArray::iterator it = tp.atoms.begin(); it != tp.atoms.end(); ++it) {
          out << format("%1$-4s %2$2.1f  %3$15.10f %4$15.10f %5$15.10f\n")
            % it->atom_name % numeric_cast<float>(it->ncharge) % it->coord(0) % it->coord(1) % it->coord(2);
        }
        out << " $END\n";
        break;
    };
    out.close();
    // header
    } catch (const Exception &e) {
      //e.fix_log();
      throw e;
    }

  } // end StructureIO::saveToFile

  void StructureIO::molToAtoms(Topology &tp) {
    tp.atoms.clear();
    pair<AtomArray::iterator, bool> res;
    FOR_ATOMS_OF_MOL(it,tp.mol) {
      Atom cur0;
      cur0.index = it->GetIdx();
      cur0.coord(0) = it->GetX();
      cur0.coord(1) = it->GetY();
      cur0.coord(2) = it->GetZ();
      cur0.mol_id   = 1;
      cur0.res_name = it->GetResidue()->GetName();
      cur0.atom_name = string("") + it->GetType();
      cur0.atom_type = string("") + it->GetType();
      cur0.charge = it->GetPartialCharge();
      cur0.mass = it->GetAtomicMass();
      cur0.ncharge = it->GetAtomicNum();
      cur0.num_connects = it->GetValence();
      int iii=0;
      FOR_NBORS_OF_ATOM(ait, &*it) {
        cur0.connects[iii] = ait->GetIdx();
        iii++;
      }
      res = tp.atoms.insert(cur0);
      if (!res.second) {
        Exception e("..something to log..");
        e.add("classname", "StructIO");
        e.add("procname", "molToAtoms");
        e.add("error", "PDB parsing error. Repeated index.");
        throw e;
      }
    }

  } // end StructureIO::molToAtoms

} // end namespace tpp
