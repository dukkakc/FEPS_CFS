from setuptools import find_packages, setup
setup(
    name='FEPS',
    packages=find_packages(),
    version='0.1.0',
    description='FEPS library',
    author='Clarence White',
    install_requires=['biopython','propy3','cogent3','numpy<1.21,>=1.17'],
    setup_requires=['pytest-runner'],
    tests_require=['pytest==4.4.1'],
    test_suite='tests',
)
