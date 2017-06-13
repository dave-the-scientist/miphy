from setuptools import setup
# VERY IMPORTANT that the files in PROGRAM_NAME.egg-info/ (created while running python setup.py sdist) have read permissions for any user. Otherwise running 'pip' to do anything on a user's system will break.

setup(
    name="MIPhy",
    version="0.8.0",

    author="Dave Curran",
    author_email="dmcurran@ucalgary.ca",
    url="https://github.com/dave-the-scientist/miphy",

    packages=["miphy_resources"],
    # Scripts are isntalled to /usr/local/bin
    scripts=["miphy.py", "miphy-tools.py"],
    # Include additional files into the package
    include_package_data=True,

    # license="LICENSE.txt",
    description="Useful towel-related stuff.",

    # long_description=open("README.txt").read(),

    # Dependent packages (distributions)
    install_requires=[
        "flask", "numpy"
    ],
)
