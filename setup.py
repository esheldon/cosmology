#!/usr/bin/env python

import sys
import os
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup

from optparse import OptionParser
parser=OptionParser(__doc__)
parser.add_option("--prefix",default=None, help="the prefix")

def create_ups(prefix):
    import distutils.config
    main_libdir=distutils.sysconfig.get_python_lib()
    pylib_install_subdir = main_libdir.replace(distutils.sysconfig.PREFIX+os.sep,'')

    pylib_install_subdir = pylib_install_subdir.replace('dist-packages','site-packages')


    upstext="""
setupOptional("numpy")
envPrepend(PYTHONPATH,${PRODUCT_DIR}/%s)
""" % pylib_install_subdir 


    upsdir=os.path.join(prefix,'ups')
    upsdir = os.path.expandvars(upsdir)
    upsdir = os.path.expanduser(upsdir)
    if not os.path.exists(upsdir):
        sys.stdout.write("Creating ups dir: %s\n" % upsdir)
        os.makedirs(upsdir)

    upsname=os.path.join(upsdir,'cosmology.table')
    sys.stdout.write('Writing ups table file: %s\n' % upsname)
    upsfile=open(upsname, 'w')
    upsfile.write(upstext)
    upsfile.close()
    del upsfile

def configuration(parent_package='',top_path=None):
    config = Configuration('cosmology',parent_package,top_path)
    config.add_extension('_cosmolib', ['cosmolib.f90'])
    return config

if __name__ == "__main__":
    with_ups=False
    try:
        ind = sys.argv.index('with_ups')
        del sys.argv[ind]
        with_ups=True
    except:
        pass
    setup(configuration=configuration)

    # only install ups if not default prefix
    if with_ups:
        options, args = parser.parse_args(sys.argv[1:])
        prefix=options.prefix
        if prefix is not None:
            create_ups(prefix)
