from distutils.core import setup, Extension

module1 = Extension('spam',
                    sources = ['gridroutines_py.c'])

setup (name = 'SpamPackage',
       version = '1.0',
       description = 'This is a demo package',
       ext_modules = [module1])
