#!/usr/local/bin/python

#take a benchmark and a fragment suffix and run all the test (if necessary)
#to get similarity and shape constraint searches

import argparse, sys, os, time, re, subprocess

def clean_verbose_helper(lines, fname):
    '''Write out none verbose lines to fname and return time'''
    f = open(fname,'w')
    t = 0
    for line in lines:
        if line.startswith('pseudo') or line.startswith(' level'):
            continue
        m = re.search(r'Found.* (\S+) s ', line)
        if m:
            t = float(m.group(1))
        else:
            f.write(line)
    f.close()
    return t
            
        
def clean_all_verbose_output(allactives, alldecoys, timename):
    '''split searchallpoints files into individual files and delete original'''
    afiles = []
    dfiles = []
    tfiles = []
    for aline in open(allactives):
        if len(aline) == 0 or aline == '\n':
            continue
        m = re.search('^# DCPTS (\d+)',aline)
        if m:
            name = allactives.replace('.allout','_%s.out' % m.group(1))
            afiles.append(name)
            tfiles.append('%s_%s.time' % (timename, m.group(1)))
            out = open(name,'w')
        else:
            out.write(aline)
            
    for dline in open(alldecoys):
        m = re.search('^# DCPTS (\d+)',dline)
        if m:
            name = alldecoys.replace('.allout','_%s.out' % m.group(1))
            dfiles.append(name)
            out = open(name,'w')
        else:
            out.write(dline)
    out.close()
                        
    for (actives,decoys,time) in zip(afiles,dfiles,tfiles):
        alines = open(actives).readlines()
        dlines = open(decoys).readlines()
        t = clean_verbose_helper(alines, actives) + clean_verbose_helper(dlines, decoys)
        with open(time,'w') as f:
            f.write('Time %f\n' % t)

        
        
def clean_verbose_output(actives, decoys, time):
    '''strip out debug information from header of file and put time in time file'''
    alines = open(actives).readlines()
    dlines = open(decoys).readlines()
    t = clean_verbose_helper(alines, actives) + clean_verbose_helper(dlines, decoys)
    with open(time,'w') as f:
        f.write('Time %f\n' % t)        

parser = argparse.ArgumentParser(description='Collect FOMS data')
parser.add_argument('dir',metavar='target directory')
parser.add_argument('fragment',metavar='fragment suffix')
parser.add_argument('-f','--force',action='store_true',help='Regenerate data even if it already exists')
args = parser.parse_args()

if not os.path.isdir(args.dir):
    print "%s is not a valid directory" % args.dir
    sys.exit(-1)

os.chdir(args.dir)

print "Running %s %s" % (args.dir, args.fragment)
#to processes a fragment suffix we need
# aligned_actives_f.db  - for FOMS
# aligned_decoys_f.db   FOMS
# selfaligned_actives_f.db VAMS
# selfaligned_actives_f.sdf rdkit
# selfaligned_decoys_f.db VAMS
# selfaligned_decoys_f.sdf rdkit

frag = args.fragment

neededfiles = ['aligned_actives_%s.db',
                'aligned_decoys_%s.db',
                'selfaligned_actives_%s.db',
                'selfaligned_actives_%s.sdf',
                'selfaligned_decoys_%s.db',
                'selfaligned_decoys_%s.sdf' ]

for f in neededfiles:
    if not os.path.exists(f % frag):
        print "Missing file %s/%s" % (args.dir, f % frag)
        sys.exit(-1)
    
#also need the lig and rec file
fragnum = frag.lstrip('f')

otherneeded = ['lig%s.pdb' % fragnum, 'rec%s.pdb' % fragnum, 'salig.pdb', 'salig.sdf','sarec.pdb']
for f in otherneeded:
    if not os.path.exists(f):
        print "Missing file %s/%s" % (args.dir, f )
        sys.exit(-1)

#get and verify counts
if not os.path.exists('active.%s.cnt' % frag) or not os.path.exists('decoy.%s.cnt' % frag) or args.force:
    a1 = subprocess.check_output('sdsorter -printCnt -reduceconfs 1 aligned_actives_%s.sdf' % frag, shell=True)
    a2 = subprocess.check_output('sdsorter -printCnt -reduceconfs 1 selfaligned_actives_%s.sdf' % frag, shell=True)
    if a1 != a2:
        print "Inconsistent active counts"
        sys.exit(-1)
    m = re.search(r'Count: (\d+)', a1)
    if not m:
        print "Unexpected count %s" % a1
        sys.exit(-1)
    with open('active.%s.cnt' % frag,'w') as f:
        f.write('%s\n' % m.group(1))

    d1 = subprocess.check_output('sdsorter -printCnt -reduceconfs 1 aligned_decoys_%s.sdf' % frag, shell=True)
    d2 = subprocess.check_output('sdsorter -printCnt -reduceconfs 1 selfaligned_decoys_%s.sdf' % frag, shell=True)
    if d1 != d2:
        print "Inconsistent decoy counts"
        sys.exit(-1)
    m = re.search(r'Count: (\d+)', d1)
    if not m:
        print "Unexpected count %s" % d1
        sys.exit(-1)
    with open('decoy.%s.cnt' %frag, 'w') as f:
        f.write('%s\n' % m.group(1))        
        
#do rdkit similarity with ligand
if not os.path.exists('rdkit.%s.time' % frag) or not os.path.exists('rdkit.decoys.%s.out' % frag) or not os.path.exists('rdkit.actives.%s.out' % frag) or args.force:
    start = time.time()
    os.system('shape.py salig.sdf selfaligned_actives_{0}.sdf -c -t > rdkit.actives.{0}.out'.format(frag))
    os.system('shape.py salig.sdf selfaligned_decoys_{0}.sdf -c -t > rdkit.decoys.{0}.out'.format(frag))
    t = time.time() - start
    with open('rdkit.%s.time' % frag,'w') as f:
        f.write('Time %f\n' % t)

#do VAMS similarity with ligand
if not os.path.exists('vams.%s.time' % frag) or not os.path.exists('vams.decoys.%s.out' % frag) or not os.path.exists('vams.actives.%s.out' % frag) or args.force:
    start = time.time()
    os.system('ShapeDB -NNSearch -k 0 -less 0 -more 0 -ligand salig.pdb -db selfaligned_actives_{0}.db/ -print -single-conformer  > vams.actives.{0}.out'.format(frag))
    os.system('ShapeDB -NNSearch -k 0 -less 0 -more 0 -ligand salig.pdb -db selfaligned_decoys_{0}.db/ -print -single-conformer  > vams.decoys.{0}.out'.format(frag))
    t = time.time() - start
    with open('vams.%s.time' % frag,'w') as f:
        f.write('Time %f\n' % t)
        
#do FOMS similarity with ligand
if not os.path.exists('foms.%s.time' % frag) or not os.path.exists('foms.decoys.%s.out' % frag) or not os.path.exists('foms.actives.%s.out' % frag) or args.force:
    start = time.time()
    os.system('ShapeDB -NNSearch -k 0 -less 0 -more 0 -ligand lig{1}.pdb -db aligned_actives_{0}.db/ -print -single-conformer  > foms.actives.{0}.out'.format(frag, fragnum))
    os.system('ShapeDB -NNSearch -k 0 -less 0 -more 0 -ligand lig{1}.pdb -db aligned_decoys_{0}.db/ -print -single-conformer  > foms.decoys.{0}.out'.format(frag, fragnum))
    t = time.time() - start
    with open('foms.%s.time' % frag,'w') as f:
        f.write('Time %f\n' % t)
        
#do FOMS similarity with ligand/receptor
if not os.path.exists('foms.rec.%s.time' % frag) or not os.path.exists('foms.rec.decoys.%s.out' % frag) or not os.path.exists('foms.rec.actives.%s.out' % frag) or args.force:
    start = time.time()
    os.system('ShapeDB -NNSearch -k 0 -less 0 -more 0 -ligand lig{1}.pdb -receptor rec{1}.pdb -db aligned_actives_{0}.db/ -print -single-conformer  > foms.rec.actives.{0}.out'.format(frag, fragnum))
    os.system('ShapeDB -NNSearch -k 0 -less 0 -more 0 -ligand lig{1}.pdb -receptor rec{1}.pdb -db aligned_decoys_{0}.db/ -print -single-conformer  > foms.rec.decoys.{0}.out'.format(frag, fragnum))
    t = time.time() - start
    with open('foms.rec.%s.time' % frag,'w') as f:
        f.write('Time %f\n' % t)
        
#shape constraint search
more = [0, .5, 1, 1.5, 2]
less = [0, .5, 1, 1.5, 2]
irad = [0, 0.5, 1, 2]

for m in more:
    for l in less:
        #ligand receptor
        name = 'foms.sc.%.2f_%.2f.%s' % (m,l,frag)
        if not os.path.exists('%s.time' % name) or not os.path.exists('%s.decoys.out' % name) or not os.path.exists('%s.actives.out' % name) or args.force:
            os.system('ShapeDB -DCSearch -less %.2f -more %.2f -ligand lig%s.pdb -receptor rec%s.pdb -db aligned_actives_%s.db/ -print -single-conformer -v > %s.actives.out' % (l, m, fragnum, fragnum, frag, name))
            os.system('ShapeDB -DCSearch -less %.2f -more %.2f  -ligand lig%s.pdb -receptor rec%s.pdb -db aligned_decoys_%s.db/ -print -single-conformer -v > %s.decoys.out' % (l, m, fragnum, fragnum, frag, name))
            clean_verbose_output('%s.actives.out'%name,'%s.decoys.out'%name,'%s.time'%name)
        
        for r in irad:
            #do FOMS shape constraint with ligand/receptor and interaction points
            name = 'foms.sci.%.2f_%.2f_%.2f.%s' % (m,l,r,frag)
            if not os.path.exists('%s.time' % name) or not os.path.exists('%s.decoys.out' % name) or not os.path.exists('%s.actives.out' % name) or args.force:
                os.system('ShapeDB -DCSearch -less %.2f -more %.2f -use-interaction-points -interaction-point-radius=%.2f -ligand lig%s.pdb -receptor rec%s.pdb -db aligned_actives_%s.db/ -print -single-conformer -v > %s.actives.out' % (l, m, r, fragnum, fragnum, frag, name))
                os.system('ShapeDB -DCSearch -less %.2f -more %.2f -use-interaction-points -interaction-point-radius=%.2f -ligand lig%s.pdb -receptor rec%s.pdb -db aligned_decoys_%s.db/ -print -single-conformer -v > %s.decoys.out' % (l, m, r, fragnum, fragnum, frag, name))
                clean_verbose_output('%s.actives.out'%name,'%s.decoys.out'%name,'%s.time'%name)

# sample all interaction pts
for m in more:
    name = 'foms.sciall.%.2f.%s' % (m,frag)
    if not os.path.exists('%s.time' % name) or not os.path.exists('%s.decoys.allout' % name) or not os.path.exists('%s.actives.allout' % name) or args.force:
        os.system('ShapeDB -SearchAllPointCombos -more %.2f -use-interaction-points -ligand lig%s.pdb -receptor rec%s.pdb -db aligned_actives_%s.db/ -print -single-conformer -v > %s.actives.allout' % (m, fragnum, fragnum, frag, name))
        os.system('ShapeDB -SearchAllPointCombos -more %.2f -use-interaction-points -ligand lig%s.pdb -receptor rec%s.pdb -db aligned_decoys_%s.db/ -print -single-conformer -v > %s.decoys.allout' % (m, fragnum, fragnum, frag, name))
        clean_all_verbose_output('%s.actives.allout'%name,'%s.decoys.allout'%name,name)
        
        
