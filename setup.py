import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

with open('LICENSE', 'r', encoding='utf-8') as fh:
    license = fh.read()

setuptools.setup(
    name='Crackling',
    version='2.0.0',
    author='Jake Bradford, Timothy Chappell, Dimitri Perrin',
    author_email='jake.bradford, dimitri.perrin (add @.qut.edu.au)',
    description='Faster and better CRISPR guide RNA design with the Crackling method',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/bmds-lab/Crackling',
    project_urls = {
        'Bug Tracker': 'https://github.com/bmds-lab/Crackling/issues',
        'Lab website': 'http://biomedicaldatascience.com/'
    },
    package_dir={'': 'src/Crackling'},
    license=license,
    install_requires=[],
    python_requires='>=3.6',
    entry_points = {
        'console_scripts': [
            'Crackling=utils.Crackling_cli:main',
            'countHitTranscripts=utils.countHitTranscripts:main',
            'extractOfftargets=utils.extractOfftargets:main',
        ],
    }
)
