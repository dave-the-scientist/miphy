# VERY IMPORTANT that the files in PROGRAM_NAME.egg-info/ (created while running
# python setup.py sdist) have read permissions for any user. Otherwise running
# 'pip' to do anything on a user's system will break.

# When testing is implemented, update MANIFEST.in to specifically include only what we're looking for

from miphy import __version__ as mi_v
from miphy_resources import __version__ as mir_v
if mi_v != mir_v:
    print('\nError: miphy and miphy_resources have different version numbers.\n')
    exit()

from setuptools import setup
setup(
    name="MIPhy",
    version=mi_v,

    author="Dave Curran",
    author_email="daves.binf.tools@gmail.com",
    url="https://github.com/dave-the-scientist/miphy",

    packages=["miphy_resources"],
    # Scripts are isntalled to /usr/local/bin
    scripts=["miphy.py", "miphy-tools.py"],
    # Include additional files into the package
    include_package_data=True,

    license="LICENSE.txt",
    description="Useful towel-related stuff.",

    long_description=open("README.txt").read(),

    # Dependent packages (distributions)
    install_requires=[
        "flask", "numpy"
    ],
)

import os
egg_dir = 'MIPhy.egg-info'
if not os.path.isdir(egg_dir):
    print('\nAttempting to check file permissions in %s/, but the folder was not found.' % egg_dir)
else:
    for fname in os.listdir(egg_dir):
        fpath = os.path.join(egg_dir, fname)
        prmsn = oct(os.stat(fpath).st_mode & 0o777)
        if int(prmsn[-1]) < 4:
            print('\nWarning: file %s with permission %s does not allow global reading; change its permissions (666 should work) and re-run this setup.py to avoid installation issues for users.' % (fpath, prmsn))
