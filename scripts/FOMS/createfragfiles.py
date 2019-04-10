#!/usr/local/bin/python -u

# given a director and fragment suffix, generate conformer and database files

import sys, subprocess, os
from shutil import rmtree

def checkf(fname):
    if not os.path.exists(fname):
        print "Missing %s" % fname
        sys.exit(-1)

if len(sys.argv) < 3:
    print "Need directory and fragment suffix"
    sys.exit(-1)

d = sys.argv[1]
frag = sys.argv[2]

print "Running %s %s" % (d,frag)

if not os.path.isdir(d):
        print "%s is not a valid directory" % d
        sys.exit(-1)

os.chdir(d)

checkf('actives.sdf')

#align to frag
f = frag.lstrip('f')
fragfile = 'frag%s.pdb' % f
checkf(fragfile)

smarts = 'smarts%s' % f
checkf(smarts)

factives = "actives_%s.sdf" % frag
infile = "aligned_actives_%s.sdf" % frag
if os.path.exists(infile) and len(sys.argv) <= 3:
    print "%s %s already computed, add -f to override" % (d,frag)
    sys.exit(0)

print "%s %s actives" % (d, frag)

os.system("obabel -s `cat %s` -isdf actives.sdf -osdf -O %s" % (smarts,factives))
os.system("obfitall `cat %s` %s %s > %s 2> /dev/null" % (smarts, fragfile, factives, infile))
#create db
dbfile = "aligned_actives_%s.db" % frag
if os.path.isdir(dbfile):
        rmtree(dbfile)
os.system("ShapeDB -Create -in %s -db %s" % (infile,dbfile)) 

#selfalign this subset
print "%s %s selfalign actives" % (d, frag)
sinfile = "selfaligned_actives_%s.sdf" % frag
if os.path.exists(sinfile):
    os.remove(sinfile)
os.system("aligncanonical.py %s %s" % (factives, sinfile))

#create dbs
dbfile = "selfaligned_actives_%s.db" % frag
if os.path.isdir(dbfile):
    rmtree(dbfile)

os.system("ShapeDB -Create -in %s -db %s" % (sinfile,dbfile))



print "%s %s decoys " % (d, frag)
fdecoys = "decoys_%s.sdf" % frag
infile = "aligned_decoys_%s.sdf" % frag
os.system("obabel -s `cat %s` -isdf decoys.sdf -osdf -O %s" % (smarts,fdecoys))
os.system("obfitall `cat %s` %s %s > %s 2> /dev/null" % (smarts, fragfile, fdecoys, infile))
#create db
dbfile = "aligned_decoys_%s.db" % frag
if os.path.isdir(dbfile):
    rmtree(dbfile)
os.system("ShapeDB -Create -in %s -db %s" % (infile,dbfile)) 

#selfalign this subset
print "%s %s selfalign decoys " % (d, frag)
sinfile = "selfaligned_decoys_%s.sdf" % frag
if os.path.exists(sinfile):
    os.remove(sinfile)
os.system("aligncanonical.py %s %s" % (fdecoys, sinfile))

#create dbs
dbfile = "selfaligned_decoys_%s.db" % frag
if os.path.isdir(dbfile):
    rmtree(dbfile)

os.system("ShapeDB -Create -in %s -db %s" % (sinfile,dbfile))

