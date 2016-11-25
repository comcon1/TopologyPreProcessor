#include "structio.hpp"
#include "exceptions.hpp"
#include "runtime.hpp"

#include "strutil.h"

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>

#include <cctype>


namespace tpp {

  using std::string;
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

  StructureIO::StructureIO() {
    // TODO: some interface stuff
    ;
  } 


  void StructureIO::loadFromFile(Topology &tp, InputFormat ifm, const char *fname) {
    // test if file exists
    runtime.log_write(string("Trying to read structure from '")+fname+"'.\n");
    fstream inf(fname, ios::in);
    if (!inf.is_open()) { 
      BOOST_CHECK(0);
      runtime.log_write("Fail to open file for read.\n");
      Parameters params;
      params.add( "procname", "tpp::load_struct");
      params.add("error", "invalid filename");
      params.add("filename", fname);
      Exception e("Can't open specified file for read.", params);
      e.fix_log();
      throw e;
    }
    // loading from stream

    loadFromStream(tp, ifm, &inf);
    inf.close();
  } // end StructureIO::loadFromFile


  void StructureIO::loadFromStream(Topology &tp, InputFormat ifm, std::istream *inf) {
    try {
      // test if tp variable contains structure
      if (!tp.atoms.empty()) {
        runtime.log_write("Replacing your current structure.\n");
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
        default: BOOST_FAIL(0);
      }
      // reading from file with OpenBabel function
      runtime.log_write("Reading by OpenBabel..");
      OBMol mol;
      if ( (!conv.Read(&mol)) || (!mol.NumAtoms()) ) {
        Parameters params;
        params.add("procname", "tpp::load_struct");
        params.add("error", "OpenBabel: parsing error");
        throw Exception("Can't read file format.", params);
      }
      runtime.log_write("OK.\n");

      tp.mol = mol;

      // Reading molecule in AtomArray format
      // ---------------------------------------
      switch(ifm) {
        case TPP_IF_PDB: {
                  char *s0  = new char[300], *SEL = new char[7], *NAM = new char[5], *RES = new char[5], *qat = new char[2]; 
                  int  res, strc, incrementalIndex = 0;
                  std::pair<AtomArray::iterator, bool> at_it;
                  Atom cur0;
                  float __x,__y,__z;
                  strc = 0;
                  bool ignoreIndexFlag = (cmdline.read("ignore_index") == "on");
                  inf->clear();
                  inf->seekg(0);
                  while (! inf->eof() ) {
                    inf->getline(s0, 300);
                    strc++;
                    //         cerr << strc << "|";

                    res = 0;
                    res += sscanf(s0, "%6s", SEL);
                    //         cerr << SEL << ".";
                    if ( strcmp(SEL, "ATOM") && strcmp(SEL, "HETATM") ) continue;

                    res += sscanf(s0 + 6, "%5u", &(cur0.oldindex));
                    cur0.index = ignoreIndexFlag ? 65535 : cur0.oldindex;
                    //         cur0.atom_id = lexical_cast<unsigned>(strcpy(s0, 
                    res += sscanf(s0 + 11, "%5s",  NAM);
                    res += sscanf(s0 + 16, "%4s",  RES);
                    // pass chain letter 21 -> 23
                    res += sscanf(s0 + 23,"%4u", &(cur0.mol_id) );
                    res += sscanf(s0 + 30, "%f", &__x );
                    res += sscanf(s0 + 38, "%f", &__y );
                    res += sscanf(s0 + 46, "%f", &__z );

                    strcpy(qat, "un");
                    if (strlen(s0) > 76)
                      res += sscanf(s0 + 76, "%2s", qat);

                    if ( res < 8 ) {
                      BOOST_CHECK(0);
                      Parameters params;
                      params.add("procname", "tpp::load_struct");
                      params.add("error", "PDB parsing error");
                      params.add("line", lexical_cast<string>(strc).c_str());
                      throw Exception("Invalid ATOM string in your PDB file.", params);
                    }

                    //         cerr << strc << ")";
                    cur0.old_aname = string(NAM);
                    cur0.atom_name = string(NAM);
                    if (cmdline.exists("rtpoutput_file")) {
                      // need to replace 1H2 to H21
                      if ((NAM[0] >= 48) && (NAM[0] <= 57))
                        cur0.atom_name = cur0.old_aname.substr(1) + cur0.old_aname.substr(0,1);
                    }
                    cur0.res_name  = string(RES);
                    cur0.qmname    = string(qat);
                    cur0.coord(0)  = numeric_cast<double>(__x);
                    cur0.coord(1)  = numeric_cast<double>(__y);
                    cur0.coord(2)  = numeric_cast<double>(__z);
                    cur0.comment   = string("QMname: ") + qat;

                    runtime.log_write( (format("%s - %s: %d(%d) [%8.3f,%8.3f,%8.3f] %s \n") % cur0.res_name % cur0.atom_name 
                          % cur0.oldindex % cur0.index % cur0.coord(0) % cur0.coord(1) % cur0.coord(2) % cur0.comment).str() );
                    at_it = tp.atoms.insert( cur0 );
                    //         cerr << strc << "@";

                    if (! at_it.second) {
                      runtime.log_write("ERROR: bad insertingo..\n");
                      Parameters params;
                      params.add("procname", "tpp::load_struct");
                      params.add("error", "PDB parsing error");
                      params.add("line", lexical_cast<string>(strc).c_str());
                      throw Exception("Repeat index or something else.", params);
                    }

                    if (ignoreIndexFlag) {
                      incrementalIndex++;
                      cur0.index = incrementalIndex;
                      tp.atoms.replace(at_it.first, cur0);
                    }
                    //         cerr << strc << "$";

                  } // end-while
                  delete[] s0; delete[] SEL; delete[] NAM; delete[] RES;
                  /* FINISH PARSIGN PDB */
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
        default: BOOST_FAIL(0);
                 break;
      };

      // setting up residue name
      if ( strutil::trim(tp.atoms.begin()->res_name) != "") {
        tp.res_name = tp.atoms.begin()->res_name;
        // checking res_name correctness    
        for (int i=0; i<tp.res_name.size(); ++i) {
          if (!isalnum(tp.res_name[i])) {
            runtime.log_write("Incorrect residue name in PDB file, using 'RES' instead.\n");
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
        Parameters params;
        params.add("procname", "tpp::load_struct");
        params.add("error", "no atoms");
        throw Exception("No atoms in structure file or bad format.", params);
      } else {
        runtime.log_write(string("Successfully readed ")+lexical_cast<string>(tp.atoms.size())+ " atoms!\n");
      }
    } catch(Exception e) { e.fix_log(); throw e;  }
  } // end loadFromStream

  void StructureIO::saveToFile(Topology &tp, OutputFormat ofm, const char *fname) {
    try {
    // test if file exists
    runtime.log_write(string("Trying to write structure into '")+fname+"'.\n");
    fstream out(fname, ios::out);
    if (!out.is_open()) { 
      runtime.log_write("Fail to open file for write.\n");
      BOOST_CHECK(0);
      Parameters params;
      params.add("procname", "tpp::load_struct");
      params.add("error", "invalid filename");
      params.add("filename", fname);
      throw Exception("Can't open specified file for write.", params);
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
    } catch (Exception e) {
      e.fix_log();
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
                Parameters params;
                params.add("procname", "tpp::StructureIO::molToAtoms");
                params.add("error", "PDB parsing error. Repeated index.");
                throw Exception("..something to log..", params);
      }
    }

  } // end StructureIO::molToAtoms

} // end namespace tpp
