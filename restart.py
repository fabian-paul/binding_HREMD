# Truncate all xtc files after the frame that was written at the last restart step.

import struct
import os
f = open('0.restart.bin')
# TODO: check if restart files are complete (simulation might have crashed during save)
# restart files are: [0-14].restart.(coor|vel|idx|xsc). Should be present and have size > 0.
t = struct.unpack('Q',f.read(8))[0]
f.close()
print 'Number of last integration step that was backed up:', t

N_rep = 14
# TODO: get this number automatically.
last_frame = t/2500 # 2500 = xtcfreq in ACEMD
# We have to calculate the frame number, because ACEMD doesn't 
# store the time information in the xtc file.
# TODO: get this number from input file.
print 'Numer of last frame (corresponding to this last integration step):', last_frame

try:
  os.mkdir('restart-backup')
except OSError, e:
  if e.errno = 17: # directory exists, just go on and overwrite content
    pass
  else:
    raise e

for i in xrange(N_rep):
  # make backup of xtc
  name = '%d.md.xtc'%i
  backup = 'restart-backup/%d.md.xtc'%i
  os.rename(name, backup)
  print 'Processing', name, '...',
  os.system(os.path.dirname(__file__)+os.sep+'truncate-xtc %d %s %s'%(last_frame,backup,name))
  print 'You may now delete the backup of the original trajectories in the subdirectory restart-backup/'



