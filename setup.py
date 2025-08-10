from mypyc.build import mypycify
from setuptools import setup

setup(
    name="twosquares",

    # mypyc docs say to just set packages simply like this:
    #   packages=['twosquares'],
    #
    # However: When I do that, twosquares/__init__.py *itself* is included in the wheel which we don't want,
    #   because then the python version will be used instead of the mypyc-compiled pyd version.
    packages=["twosquares-stubs"],
    include_package_data=True,
    package_data={'twosquares-stubs': ["*.pyi"]},

    ext_modules=mypycify([
        "twosquares/__init__.py",
    ]),

    license="MIT",
)
