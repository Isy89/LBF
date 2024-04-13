from setuptools import find_namespace_packages, setup

entry_points = [
    "lbfextract = lbfextract.__main__:cli",
]

with open("requirements.txt", "r", encoding="utf-8") as f:
    requires = []
    for line in f:
        req = line.split("#", 1)[0].strip()
        requires.append(req)


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name="LBFextract",
    version="0.1.0a1",
    author="Isaac Lazzeri",
    description="Cli to extract different liquid biopsy related genomics signals from bam files",
    long_description=readme(),
    url="https://github.com/isy89/LBF",
    packages=find_namespace_packages(where='src'),
    package_dir={
        '': 'src',
    },
    classifiers=[
        'Programming Language :: Python :: 3.9',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Development Status :: 3 - Alpha'
    ],
    entry_points={"console_scripts": entry_points,
                  "lbfextract": [
                      "entropy = lbfextract.fextract_entropy.plugin:hook",
                      "coverage_in_batch = lbfextract.fextract_batch_coverage.plugin:hook",
                      "entropy_in_batch = lbfextract.fextract_entropy_in_batch.plugin:hook",
                      "fragment_length_distribution = lbfextract.fextract_fragment_length_distribution.plugin:hook",
                      "fragment_length_distribution_in_batch = lbfextract.fextract_fragment_length_distribution_in_batch.plugin:hook",
                      "relative_entropy_to_flanking = lbfextract.fextract_relative_entropy_to_flanking.plugin:hook",
                      "relative_entropy_to_flanking_in_batch = lbfextract.fextract_relative_entropy_to_flanking_in_batch.plugin:hook",
                      "fragment_length_ratios = lbfextract.fextract_fragment_length_ratios.plugin:hook",
                  ],
                  "lbfextract_cli": [
                      "extract_entropy = lbfextract.fextract_entropy.plugin:hook_cli",
                      "extract_coverage_in_batch = lbfextract.fextract_batch_coverage.plugin:hook_cli",
                      "extract_entropy_in_batch = lbfextract.fextract_entropy_in_batch.plugin:hook_cli",
                      "extract_fragment_length_distribution = lbfextract.fextract_fragment_length_distribution.plugin:hook_cli",
                      "extract_fragment_length_distribution_in_batch = lbfextract.fextract_fragment_length_distribution_in_batch.plugin:hook_cli",
                      "extract_relative_entropy_to_flanking = lbfextract.fextract_relative_entropy_to_flanking.plugin:hook_cli",
                      "extract_relative_entropy_to_flanking_in_batch = lbfextract.fextract_relative_entropy_to_flanking_in_batch.plugin:hook_cli",
                      "extract_fragment_length_ratios = lbfextract.fextract_fragment_length_ratios.plugin:hook_cli",

                  ]},
    project_urls={"Homepage": "https://github.com/Isy89/LBF",
                  "Bug Tracker": "https://github.com/Isy89/LBF/issues"},
    include_package_data=True,
    install_requires=requires,
    python_requires=">=3.10",
    extras_require={},
)
