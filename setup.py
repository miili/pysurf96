from numpy.distutils.core import setup
from numpy.distutils.extension import Extension

setup(
    ext_modules=[
        Extension(
            name="pysurf96.surfdisp96_ext",
            sources=["pysurf96/surfdisp96.f"],
            extra_f77_compile_args=[
                "-O3",
                "-ffixed-line-length-none",
                "-fbounds-check",
                "-m64",
            ],
            f2py_options=["only:", "surfdisp96", ":"],
            language="f77",
        )
    ]
)
