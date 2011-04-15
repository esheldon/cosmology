import distutils
from distutils.core import setup, Extension, Command
import os
import numpy

data_files=[]

class AddUPS(Command):
    _data_files = data_files
    user_options=[]
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        # create the ups table

        main_libdir=distutils.sysconfig.get_python_lib()
        pylib_install_subdir = main_libdir.replace(distutils.sysconfig.PREFIX+os.sep,'')
        pylib_install_subdir = pylib_install_subdir.replace('dist-packages','site-packages')


        if not os.path.exists('ups'):
            os.mkdir('ups')
        tablefile=open('ups/cosmology.table','w')
        tab="""
        setupOptional("python")
        setupOptional("numpy")
        envPrepend(PYTHONPATH,${PRODUCT_DIR}/%s)
        """ % pylib_install_subdir
        tablefile.write(tab)
        tablefile.close()

        AddUPS._data_files.append(('ups',['ups/cosmology.table']))

ext=Extension("cosmology._cosmolib", 
              ["cosmology/cosmolib_pywrap.c","cosmology/cosmolib.c"])
setup(name="cosmolib", 
      packages=['cosmology'],
      version="1.0",
      cmdclass={"with_ups": AddUPS},
      data_files=data_files,
      ext_modules=[ext],
      include_dirs=numpy.get_include())




