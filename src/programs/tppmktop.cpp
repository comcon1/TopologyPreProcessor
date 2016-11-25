/*! \file tppmktop.cpp
 *
 *	\briefThis file provides executable for TPPMKTOP utility.
 *
 */

#include "global.hpp"
#include "core.hpp"
#include "exceptions.hpp"
#include "runtime.hpp"
#include "topio.hpp"
#include "structio.hpp"
#include "db_base.hpp"
#include "atom_definer.hpp"
#include "bond_definer.hpp"

#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>

#include <boost/format.hpp>

namespace p_o = boost::program_options;
using tpp::cmdline;
using tpp::Parameters;

using boost::format;

using std::cout;
using std::cerr;
using std::endl;

using std::string;
void helpscreen();
double sumcharge(const tpp::Topology &);

int main(int argc, char * argv[]) {
	string progname("Execution rules for TPPMKTOP ");
	progname = progname + VERSION;
	p_o::options_description desc(progname);
	p_o::variables_map vars;
	desc.add_options()("input,i", p_o::value<std::string>(),
			"Input filename (any format)")("output,o",
			p_o::value<std::string>(), "Output filename (itp format)")(
			"rtp-output,r", p_o::value<std::string>(),
			"Output filename (rtp format)")("forcefield,f",
			p_o::value<std::string>(), "Forcefield name")

	("lack-file,l", p_o::value<std::string>(),
			"Topology lack filename (default 'lack.itp')")("sqlserver,s",
			p_o::value<std::string>(),
			"Mysql-server adress (default 'localhost')")("sqlport,t",
			p_o::value<unsigned>(), "Mysql-server port (default '3306')")(
			"sqluser,u", p_o::value<std::string>(),
			"Mysql-user (default 'tppuser')")("sqlpassword,p",
			p_o::value<std::string>(), "Mysql-password (default 'estatic')")(
			"nocalculate,n", "Create final topology (don't create lack-file)")(
			"max-bonds,m",
			"Maximize amount of bonds, angles and dihedrals by selecting other atom-types.")(
			"verbose,v", "Verbose mode")("help,h", "Print this message");
	try {
		try { // parsing boost::program_options
			p_o::store(p_o::parse_command_line(argc, argv, desc), vars);
			p_o::notify(vars);

			// boolean options
			if ((vars.count("verbose") > 1) || (vars.count("nocalculate") > 1)
					|| (vars.count("max-bonds") > 1))
				throw 1;
			cmdline.add("verbose_flag",
					vars.count("verbose") ? "on" : "off");
			cmdline.add("nocalculate_flag",
					vars.count("nocalculate") ? "on" : "off");
			if (vars.count("max-bonds"))
				cmdline.add("max-bonds", "on");
			if (vars.count("help") == 1)
				helpscreen();

			// main string options
			if (vars.count("input") == 1) {
				cmdline.add("input_file",
						vars["input"].as<std::string>());
			} else
				throw 1;
			if (vars.count("output") == 1) {
				cmdline.add("output_file",
						vars["output"].as<std::string>());
			} else
				throw 1;
			if (vars.count("forcefield") == 1) {
				cmdline.add("forcefield",
						vars["forcefield"].as<std::string>());
			} else
				throw 1;
			// optional parameters
			if (vars.count("lack-file") == 1) {
				cmdline.add("lack_file",
						vars["lack-file"].as<std::string>());
			} else if (vars.count("lack-file") == 0) {
				cmdline.add("lack_file", "lack.itp");
			} else
				throw 1;
			if (vars.count("rtp-output") == 1) {
				cmdline.add("rtpoutput_file",
						vars["rtp-output"].as<std::string>());
			} else if (vars.count("rtp-output") > 0) {
				throw 1;
			}
			// SQL parameters
			if (vars.count("sqlserver") == 1) {
				cmdline.add("sqlserver",
						vars["sqlserver"].as<std::string>());
			} else if (vars.count("sqlserver") == 0) {
				cmdline.add("sqlserver", "localhost");
			} else
				throw 1;
			if (vars.count("sqluser") == 1) {
				cmdline.add("sqluser",
						vars["sqluser"].as<std::string>());
			} else if (vars.count("sqluser") == 0) {
				cmdline.add("sqluser", "tppuser");
			} else
				throw 1;
			if (vars.count("sqlport") == 1) {
				unsigned o = vars["sqlport"].as<unsigned>();
				cmdline.add("sqlport",
						boost::lexical_cast<string>(
								vars["sqlport"].as<unsigned>()));
			} else if (vars.count("sqlport") == 0) {
				cmdline.add("sqlport", "3306");
			} else
				throw 1;
			if (vars.count("sqlpassword") == 1) {
				cmdline.add("sqlpassword",
						vars["sqlpassword"].as<std::string>());
			} else if (vars.count("sqlpassword") == 0) {
				cmdline.add("sqlpassword", "estatic");
			} else
				throw 1;
		} catch (boost::program_options::error &e) {
			throw 1;
		}
	} catch (int ExC) {
		if (ExC) {
			cerr
					<< format(
							"\nTPPMKTOP %1% : Error in input parameters.\n\n") % VERSION;
			cerr << desc;
		}
		return (ExC);
	}
	// finish analysing
	// starting work with input and output files

	if (cmdline.read("verbose_flag") == "on") {
		cout
				<< format(
						"\
**********************************************************************\n\
*   Biology faculty, Department of biophysics, Erg Research Group    *\n\
*   Moscow, Lomonosov's Moscow State University                      *\n\
*   for more info, see homepage  http://erg.biophys.msu.ru/          *\n\
*                                                                    *\n\
*   Authors:       comcon1, dr.zoidberg, piton, leela                *\n\
*                                                                    *\n\
*   Product:       program  TPPMKTOP-%1$-6s                          *\n\
*                                                                    *\n\
*    Utilite for generating final topology from PDB file according   *\n\
* to SMARTS-patterns in database. Also you can maximize correctness  *\n\
* of topology by using special algorythmes.                          *\n\
*                                                                    *\n\
* Configured:     %2$-19s                                *\n\
**********************************************************************\n\
\n\n") % PACKAGE_VERSION % CONFIGURE_CDATE;
	} else {
		cout << format("Starting %1$s program.\n") % "TPPMKTOP";
	}

	// INPUT analysing
	tpp::InputFormat iform;
	string::size_type ind = cmdline.read("input_file").find(".", 0);
	if (ind == string::npos) {
		cerr << "ERROR:\n";
		cerr
				<< "Couldn't determine format of input file. Please specify extension.\n";
		return 1;
	}
	string subs = cmdline.read("input_file").substr(ind + 1);
	if (subs == "pdb")
		iform = tpp::TPP_IF_PDB;
	else if (subs == "gro")
		iform = tpp::TPP_IF_GRO;
	else if (subs == "g96")
		iform = tpp::TPP_IF_G96;
	else if ((subs == "log") || (subs == "out"))
		iform = tpp::TPP_IF_GAMSP;
	else {
		cerr << "ERROR:\n";
		cerr
				<< "Couldn't determine format of input file. Please specify other extension.\n";
		return 1;
	}

	if (cmdline.read("verbose_flag") == "on") {
		switch (iform) {
		case tpp::TPP_IF_PDB:
			cout << "Input file format: Protein Data Bank." << endl;
			break;
		case tpp::TPP_IF_GRO:
			cout << "Input file format: GROmacs structure." << endl;
			break;
		case tpp::TPP_IF_G96:
			cout << "Input file format: Gromos 96 structure." << endl;
			break;
		case tpp::TPP_IF_GAMSP:
			cout << "Input file format: GAMess output." << endl;
			break;
		};
	}

	if (cmdline.exists("nocalculate")) {
		cout << "TPPMKTOP will try to make full-determined topology!" << endl;
	}

	// program body, using modules
	try {
		tpp::Topology TOP;
                tpp::StructureIO sio;
		// setting up common topology parameters
		TOP.res_name = cmdline.read("input_file").substr(0, 3);
		TOP.nrexcl = 3;
		// ;-)
		sio.loadFromFile(TOP, iform, cmdline.read("input_file").c_str());
		// customization of 2-nd level parameters
		tpp::Parameters par0;
		par0.add("host", cmdline.read("sqlserver"));
		par0.add("dbname", "tppforcefield");
		par0.add("user", cmdline.read("sqluser"));
		par0.add("password", cmdline.read("sqlpassword"));
		par0.add("port", cmdline.read("sqlport"));
		par0.add("ffname", cmdline.read("forcefield"));
		if (cmdline.exists("max-bonds")) {
			par0.add("maxbonds", "on");
			par0.add("maxangles", "on");
			par0.add("maxdihedrals", "on");
		}
		// initial DB queries
		tpp::DbInfo DI(par0);
		par0.add("ffid", boost::lexical_cast<string>(DI.get_ffid()));
		TOP.ffinclude = DI.get_ffinclude().c_str();
		TOP.ffinfo = par0.read("ffname") + " revision " + DI.get_ffrev();
		//
		if (cmdline.read("verbose_flag") == "on") {
			cout << DI.get_statistics();
		}
		// starting program body
		tpp::AtomDefiner AD(par0, TOP);
		AD.proceed();
		AD.log_scores();
		AD.atom_align();
		tpp::BondDefiner BD(par0, TOP);
		BD.bond_align();
		tpp::save_topology(TOP, cmdline.read("output_file").c_str());
		tpp::save_lack(TOP, cmdline.read("lack_file").c_str());
		if (cmdline.exists("rtpoutput_file")) {
			tpp::save_topology_rtp(TOP,
					cmdline.read("rtpoutput_file").c_str());
		}
		cout
				<< format(
						"Please, correct your charges according to sum: %1$8.3f.\n")
						% sumcharge(TOP);
	} catch (tpp::SqlException &e) {
		e.fix_log();
		cerr << "TPP_SQL_EXCEPTION FROM: " << e["procname"] << endl;
		cerr << "more info see in log-file." << endl;
		return 3;
	} catch (tpp::DbException &e) {
		e.fix_log();
		cerr << "TPP_DB_EXCEPTION FROM: " << e["procname"] << endl;
		cerr << "more info see in log-file." << endl;
		return 2;
	} catch (tpp::Exception &e) {
		e.fix_log();
		cerr << "TPP_EXCEPTION FROM: " << e["procname"] << endl;
		cerr << "more info see in log-file." << endl;
		return 1;
	}

	cout << "TPPMKTOP finished normally!" << endl;
}

void helpscreen() {
	cout
			<< format(
					"\n\
--------------------------------*****---------------------------------\n\
                                          ERG Research Group,         \n\
                                          Department of biophysics,   \n\
                                          Biology faculty, MSU, Russia\n\
\n\
           THE PART OF TOPOLOGY PREPROCESSOR PROJECT                  \n\
        ---      (comcon1, zoidberg, piton, leela)      ---           \n\
  TPP version: %1$-3s, compiled at %2$-8s on GCC %3$s.\n\
  BOOST version:  %4$-8s \n\
  OpenBabel version: %5$-8s \n\
  OpenBabel share data: %6$-8s \n\
\n\
                                TPPMKTOP\n\
   Utilite for checking your structure file to be suite for  next-step\n\
programs.  Also it adapts names and  position of atoms in file to make\n\
following topology file more obvious.\n\
\n\
 USAGE: \n\
 tppmktop -i <inp> -o <out> -f <f.field> [-l <lack>] [other opt-s]\n\
\n\
      -i  the name of (I)nput-file, in PDB or GRO/G96 format.           \n\
      -o  the name of (O)utput-file, contained prepared structure.      \n\
      -f  the (F)orcefield name (f.i. OPLS-AA)                          \n\
      -v  (V)erbose mode, typing more information during the execution\n\
 [ special topopolgy generation settings ]\n\
      -n  do (N)ot calculate force parameters. Write final ITP.\n\
      -l  specify topology (L)ACK-file definition.\n\
      -m  (M)aximize amount of bonded interactions.\n\
 [ database options ] \n\
      -s  MySQL (S)erver host name or IP\n\
      -t  MySQL server (P)ort number\n\
      -u  MySQL (U)ser\n\
      -p  MySQL (P)assword\n\
      -h  print this message.                                         \n\
\n\
--------------------------------*****---------------------------------\n\
")
					% PACKAGE_VERSION % CONFIGURE_CDATE % __VERSION__
					% BOOST_LIB_VERSION % BABEL_VERSION % BABEL_DATADIR << endl;
	throw 0;
}

double sumcharge(const tpp::Topology &tp) {
	double sum = 0.0;
	for (tpp::AtomArray::iterator it = tp.atoms.begin();
			it != tp.atoms.end(); ++it) {
		sum += it->charge;
	}
	return sum;
}
