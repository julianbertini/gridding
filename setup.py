from distutils.core import setup, Extension

module1 = Extension('grid',
                    sources = ['gridroutines_py.c'])

setup (name = 'GridPackage',
       version = '1.0',
       description = 'This is a gridding package.',
       ext_modules = [module1])
