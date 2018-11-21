#!/usr/bin/env python
try:
    from numpy.distutils.core import Extension
    from numpy.distutils.core import setup
except ImportError:
    raise ImportError('Numpy needs to be installed!')


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
            extra_f77_compile_args='-O3 -ffixed-line-length-none -fbounds-check -m64'.split(), # noqa
            f2py_options=['only:', 'surfdisp96', ':'],
            language='f77')
        ]
)
