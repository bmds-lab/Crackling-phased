import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

with open('LICENSE', 'r', encoding='utf-8') as fh:
    license = fh.read()

setuptools.setup(
    name='haplocrackling',
    version='1.0.0',
    author='Jake Bradford, Timothy Chappell, Dimitri Perrin',
    author_email='jake.bradford, dimitri.perrin (add @.qut.edu.au)',
    description='Faster and better CRISPR guide RNA design with the Crackling method',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/bmds-lab/Crackling-phased',
    project_urls = {
        'Bug Tracker': 'https://github.com/bmds-lab/Crackling-phased/issues',
        'Lab website': 'http://biomedicaldatascience.com/'
    },
    package_dir={'': 'src'},
    packages=setuptools.find_packages(where='src'),
    license=license,
    install_requires=[],
    python_requires='>=3.6',
    entry_points = {
        'console_scripts': [
            'Haplocrackling=haplocrackling.utils.HaploCrackling_cli:main',
            'countHitTranscripts=haplocrackling.utils.countHitTranscripts:main',
            'extractOfftargets=haplocrackling.utils.extractOfftargets:main',
            'trainModel=haplocrackling.utils.trainModel:main'
        ],
    },
    include_package_data=True,
    package_data={'': ['data/*']},
)
