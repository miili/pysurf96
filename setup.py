#!/usr/bin/env python
try:
    from numpy.distutils.core import Extension
    from numpy.distutils.core import setup
except ImportError:
    class NoNumpy(Exception):
        pass

    raise NoNumpy('Numpy Needs to be installed '
                  'for extensions to be compiled.')


setup(
    name='pysurf96',
    version='0.1',
    description='Surface Wave Dispersion Python Wrapper for surf96',
    author='MDM',
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    package_dir={
        'pysurf96': 'src'
    },
    packages=[
        'pysurf96',
    ],
    ext_modules=[
        Extension(
            name='pysurf96.surfdisp96_ext',
            sources=['src/surfdisp96.f'],
            extra_f77_compile_args='-g3 -ffixed-line-length-none -ffloat-store -W -fbounds-check -m64 -mcmodel=medium'.split(), # noqa
            f2py_options=['only:', 'surfdisp96', ':'],
            language='f77')
        ]
)
