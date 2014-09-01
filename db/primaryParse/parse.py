#!/usr/bin/python
from numpy import *
from matplotlib.mlab import *

print 'INSERT INTO `atoms` (`id`, `ffield`, `uname`, `name`, `znuc`, `mass`, `charge`, `ptype`, `c6`, `c12`) VALUES '
ff = loadtxt('/usr/share/gromacs/top/parmbsc0.ff/ffnonbonded.itp',skiprows=1,dtype=str, comments=';')
i = 0
for atom in ff:
    i += 1
    print "(NULL, 8, '%s','%s', '%d', '%8.3f', '%8.3f', '%1s', '%10.6g', '%10.6g' )"\
            % (atom[0], atom[0], int(atom[1]), float(atom[2]), float(atom[3]),\
            atom[4], float(atom[5]), float(atom[6]) ),
    if (i < len(ff)):
        print ','
    else:
        print ';\n'

ff = loadtxt('/usr/share/gromacs/top/parmbsc0.ff/atomtypes.atp',dtype=str, delimiter=';')
for atom in ff:
    print "UPDATE `atoms` SET comment = '%s' WHERE `ffield` = 8 and uname = '%s';" \
            % ( atom[1].replace("'","\\'"), atom[0].split()[0]  )

ff1 = genfromtxt('/usr/share/gromacs/top/parmbsc0.ff/ffbonded.itp', skip_header=2, skip_footer=479, comments=';',dtype=str)
ff2 = genfromtxt('/usr/share/gromacs/top/parmbsc0.ff/ffbonded.itp', skip_header=2, skip_footer=479, delimiter=' ; ',dtype=str)
print "\n\nINSERT INTO `bonds` (`id`, `ffield`, `i`, `j`, `f`, `c1`, `c2`, `comment`) VALUES"
for i in range(len(ff1)):
    opts = ff1[i,:]
    com = ff2[i,1]
    print "(NULL, 8, '%s', '%s', %d, %10.5f, %10.5f, '%s')" \
            % (opts[0], opts[1], int(opts[2]), float(opts[3]), float(opts[4]), com.replace("'", "\\'")),
    if (i < len(ff1)-1):
        print ','
    else:
        print ';\n'

ff1 = genfromtxt('/usr/share/gromacs/top/parmbsc0.ff/ffbonded.itp', skip_header=127, skip_footer=212, comments=';',dtype=str)
ff2 = genfromtxt('/usr/share/gromacs/top/parmbsc0.ff/ffbonded.itp', skip_header=127, skip_footer=212, delimiter=' ;',dtype=str)
print "\n\nINSERT INTO `angles` (`id`, `ffield`, `i`, `j`,`k`, `f`, `c1`, `c2`, `comment`) VALUES"
for i in range(len(ff1)):
    opts = ff1[i,:]
    com = ff2[i,1]
    print "(NULL, 8, '%s', '%s', '%s', %d, %10.5f, %10.5f, '%s')" \
            % (opts[0], opts[1], opts[2], int(opts[3]), float(opts[4]), float(opts[5]), com.replace("'", "\\'")),
    if (i < len(ff1)-1):
        print ','
    else:
        print ';\n'


ff1 = genfromtxt('/usr/share/gromacs/top/parmbsc0.ff/ffbonded.itp', skip_header=372, skip_footer=3, comments=';',dtype=str)
ff2 = genfromtxt('/usr/share/gromacs/top/parmbsc0.ff/ffbonded.itp', skip_header=372, skip_footer=3, delimiter=' ;',dtype=str)
print "\n\nINSERT INTO `angles` (`id`, `ffield`, `i`, `j`,`k`,`l`, `f`, `c1`, `c2`, `c3`, `comment`) VALUES"
for i in range(len(ff1)):
    opts = ff1[i,:]
    com = ff2[i,1]
    print "(NULL, 8, '%s', '%s', '%s', '%s', %d, %10.5f, %10.5f, %d, '%s')" \
            % (opts[0], opts[1], opts[2], opts[3], int(opts[4]), float(opts[5]), float(opts[6]), int(opts[7]), com.replace("'", "\\'")),
    if (i < len(ff1)-1):
        print ','
    else:
        print ';\n'
