import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

with open('LICENSE', 'r', encoding='utf-8') as fh:
    license = fh.read()

setuptools.setup(
    name='crackling',
    version='2.0.1',
    author='Jake Bradford, Carl Schmitz, Timothy Chappell, Dimitri Perrin',
    author_email='jake.bradford, dimitri.perrin (add @.qut.edu.au)',
    description='Faster and better CRISPR guide RNA design with the Crackling method',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/bmds-lab/Crackling',
    project_urls = {
        'Bug Tracker': 'https://github.com/bmds-lab/Crackling/issues',
        'Lab website': 'http://biomedicaldatascience.com/'
    },
    package_dir={'': 'src'},
    license=license,
    install_requires=[],
    python_requires='>=3.8',
    entry_points = {
        'console_scripts': [
            'Crackling=crackling.utils.Crackling_cli:main',
            'countHitTranscripts=crackling.utils.countHitTranscripts:main',
            'extractOfftargets=crackling.utils.extractOfftargets:main',
            'trainModel=crackling.utils.trainModel:main'
        ],
    },
    include_package_data=True,
    package_data={'': ['data/*']},
)
