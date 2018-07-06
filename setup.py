from setuptools import setup

setup(
    name='covfuzze',
    version='0.1.2',
    py_modules=['covfuzze'],
    author='al-mcintyre',
    author_email='abm237@cornell.edu',
    url='https://github.com/al-mcintyre/CovFuzze',
    install_requires=["numpy >= 1.9.1","pysam >= 0.6", "pandas >= 0.17.0", "seaborn >= 0.7.1", "matplotlib >= 2.0.2", \
            "pybedtools >= 0.7.10"],
    license='LICENSE',
    package_dir = {'': 'covfuzze'},
    entry_points={
       'console_scripts': [
           'covfuzze = covfuzze:main',
       ],
    },
    long_description=open('README.md').read(),
)
