#!/usr/bin/env python

import sys
import os
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup
import distutils.config

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

    upsname=os.path.join(upsdir,'fimage.table')
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
    if sys.argv.count("install_ups") != 0:
        # this is special, just install ups
        # I did it this way because a command class can *only*
        # install under the main package dir, not prefix/ups for example.
        # and setting config_fc --fcompiler=gnu95 disallows the with_ups
        # option for some reason
        options, args = parser.parse_args(sys.argv[1:])
        prefix=options.prefix
        if prefix is not None:
            create_ups(prefix)
        else:
            raise ValueError("can only install ups when prefix is set")

        sys.exit(0)

    setup(configuration=configuration)


